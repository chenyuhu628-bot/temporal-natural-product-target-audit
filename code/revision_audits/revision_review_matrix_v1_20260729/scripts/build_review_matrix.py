#!/usr/bin/env python
"""Build an aggregate-only audit of the major-revision reviewer report.

The script reads identifier-bearing locked inputs locally, but writes only
aggregate counts, distribution summaries, hashes, and a reviewer-action
matrix. Temporary MMseqs2 FASTA and hit files are deleted before exit.
"""

from __future__ import annotations

import argparse
import calendar
import csv
import gzip
import hashlib
import json
import subprocess
import tempfile
from collections import defaultdict
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator
from sklearn import __version__ as sklearn_version
from sklearn.feature_extraction.text import TfidfVectorizer


EXPECTED_PROTOCOL_SHA256 = "bd4902476160cc7c5cbacaf0cfd0f1a28c5300bde22232b3a2cc6c1f3c143dc2"
EXPECTED_REVIEW_SHA256 = "c40be1c0ee7e3d077da6d0bb89476a94c8b894508c239633826e7a4e10e24888"
CUTOFF = date(2022, 8, 31)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("rt", encoding="utf-8", newline="")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, Any]]) -> None:
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_bool(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes", "y"}:
        return True
    if normalized in {"0", "false", "no", "n"}:
        return False
    raise ValueError(f"Unrecognized Boolean value: {value!r}")


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []
    with path.open("rt", encoding="utf-8", newline="") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    sequences[current] = "".join(chunks)
                current = line[1:].split()[0]
                if current in sequences:
                    raise ValueError("Duplicate FASTA identifier")
                chunks = []
            else:
                if current is None:
                    raise ValueError("Sequence encountered before FASTA header")
                chunks.append(line)
    if current is not None:
        sequences[current] = "".join(chunks)
    return sequences


def distribution(values: Iterable[float]) -> dict[str, Any]:
    array = np.asarray(list(values), dtype=np.float64)
    if array.size == 0:
        return {
            "n_observed": 0,
            "n_at_upper_bound": 0,
            "min": "",
            "q1": "",
            "median": "",
            "q3": "",
            "max": "",
            "mean": "",
        }
    array = np.clip(array, 0.0, 1.0)
    return {
        "n_observed": int(array.size),
        "n_at_upper_bound": int(np.isclose(array, 1.0, atol=1e-7, rtol=0.0).sum()),
        "min": f"{array.min():.10f}",
        "q1": f"{np.quantile(array, 0.25):.10f}",
        "median": f"{np.quantile(array, 0.50):.10f}",
        "q3": f"{np.quantile(array, 0.75):.10f}",
        "max": f"{array.max():.10f}",
        "mean": f"{array.mean():.10f}",
    }


def summary_row(
    *,
    family: str,
    unit: str,
    scope: str,
    values: Iterable[float],
    n_total: int,
    n_no_detected_alignment: int = 0,
    n_partially_censored: int = 0,
    n_scope_units_with_all_relations_detected: int | str = "",
    interpretation: str,
) -> dict[str, Any]:
    stats = distribution(values)
    return {
        "similarity_family": family,
        "analysis_unit": unit,
        "scope": scope,
        "n_total": n_total,
        **stats,
        "n_no_detected_alignment": n_no_detected_alignment,
        "n_partially_censored": n_partially_censored,
        "n_scope_units_with_all_relations_detected": n_scope_units_with_all_relations_detected,
        "interpretation": interpretation,
    }


def build_scopes(
    evaluation_rows: list[dict[str, str]],
    scaffold_rows: list[dict[str, str]],
    homology_rows: dict[str, list[dict[str, str]]],
) -> tuple[dict[str, list[dict[str, str]]], dict[str, Any]]:
    scaffold_pairs = {
        row["canonical_pair_key"]
        for row in scaffold_rows
        if parse_bool(row["audit_scaffold_cold_under_selected_policy"])
    }
    cold_targets = {
        threshold: {
            row["uniprot_canonical_accession"]
            for row in rows
            if parse_bool(row["is_future_target_homology_cold_candidate"])
        }
        for threshold, rows in homology_rows.items()
    }
    scopes: dict[str, list[dict[str, str]]] = {
        "temporal_strict_ab": list(evaluation_rows),
        "scaffold_cold": [
            row for row in evaluation_rows if row["canonical_pair_key"] in scaffold_pairs
        ],
    }
    for threshold in ("0_30", "0_50", "0_70"):
        scopes[f"joint_scaffold_homology_{threshold}"] = [
            row
            for row in evaluation_rows
            if row["canonical_pair_key"] in scaffold_pairs
            and row["uniprot_canonical_accession"] in cold_targets[threshold]
        ]
    keys_050 = {
        row["canonical_pair_key"] for row in scopes["joint_scaffold_homology_0_50"]
    }
    keys_070 = {
        row["canonical_pair_key"] for row in scopes["joint_scaffold_homology_0_70"]
    }
    if keys_050 != keys_070:
        raise AssertionError("The locked 0.50 and 0.70 relation masks are not identical")
    scopes["joint_scaffold_homology_0_50_0_70_identical"] = scopes[
        "joint_scaffold_homology_0_50"
    ]
    del scopes["joint_scaffold_homology_0_50"]
    del scopes["joint_scaffold_homology_0_70"]
    counts: dict[str, Any] = {}
    for scope, rows in scopes.items():
        counts[scope] = {
            "relations": len(rows),
            "queries": len({row["query_id"] for row in rows}),
            "targets": len({row["uniprot_canonical_accession"] for row in rows}),
        }
    counts["identity_0_50_0_70_relation_mask_equal"] = True
    counts["identity_0_50_0_70_relation_set_sha256"] = hashlib.sha256(
        "\n".join(sorted(keys_050)).encode("utf-8")
    ).hexdigest()
    return scopes, counts


def morgan_maxima(
    historical_compounds: list[dict[str, str]],
    query_compounds: list[dict[str, str]],
) -> dict[str, float]:
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

    def fingerprint(smiles: str):
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError("RDKit could not parse a locked compound structure")
        return generator.GetFingerprint(molecule)

    historical_fingerprints = [
        fingerprint(row["representative_smiles"]) for row in historical_compounds
    ]
    maxima: dict[str, float] = {}
    for row in query_compounds:
        query_fp = fingerprint(row["representative_smiles"])
        similarities = DataStructs.BulkTanimotoSimilarity(
            query_fp, historical_fingerprints
        )
        maxima[row["inchikey_full"]] = float(max(similarities))
    return maxima


def sequence_3mer_maxima(
    candidate_targets: list[str],
    sequences: dict[str, str],
    historical_targets: set[str],
) -> dict[str, float]:
    if set(candidate_targets).difference(sequences):
        raise ValueError("Candidate target is absent from the locked sequence FASTA")
    index = {target: position for position, target in enumerate(candidate_targets)}
    if historical_targets.difference(index):
        raise ValueError("Historical target is absent from candidate universe")
    vectorizer = TfidfVectorizer(
        analyzer="char",
        ngram_range=(3, 3),
        lowercase=False,
        norm="l2",
        dtype=np.float32,
    )
    matrix = vectorizer.fit_transform([sequences[target] for target in candidate_targets])
    historical_indices = np.asarray(
        [index[target] for target in sorted(historical_targets)], dtype=np.int64
    )
    historical_matrix = matrix[historical_indices]
    maxima = np.zeros(len(candidate_targets), dtype=np.float32)
    block_size = 256
    for start in range(0, len(candidate_targets), block_size):
        stop = min(start + block_size, len(candidate_targets))
        similarities = (matrix[start:stop] @ historical_matrix.T).toarray()
        maxima[start:stop] = np.asarray(similarities.max(axis=1), dtype=np.float32)
    maxima = np.clip(maxima, 0.0, 1.0)
    return {
        target: float(maxima[position]) for target, position in index.items()
    }


def summarize_complete_similarity(
    *,
    family: str,
    value_by_entity: dict[str, float],
    entity_field: str,
    scopes: dict[str, list[dict[str, str]]],
    query_reducer: str,
    interpretation: str,
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for scope, rows in scopes.items():
        relation_values = [value_by_entity[row[entity_field]] for row in rows]
        by_query: dict[str, list[float]] = defaultdict(list)
        for row in rows:
            by_query[row["query_id"]].append(value_by_entity[row[entity_field]])
        if query_reducer == "constant":
            query_values = [values[0] for values in by_query.values()]
        elif query_reducer == "max":
            query_values = [max(values) for values in by_query.values()]
        else:
            raise ValueError(f"Unsupported query reducer: {query_reducer}")
        output.append(
            summary_row(
                family=family,
                unit="relation_weighted",
                scope=scope,
                values=relation_values,
                n_total=len(rows),
                interpretation=interpretation,
            )
        )
        output.append(
            summary_row(
                family=family,
                unit="query",
                scope=scope,
                values=query_values,
                n_total=len(by_query),
                interpretation=interpretation,
            )
        )
    return output


def run_mmseqs_detected_identity(
    *,
    mmseqs_exe: Path,
    sequences: dict[str, str],
    future_targets: set[str],
    historical_targets: set[str],
    threads: int,
) -> tuple[dict[str, float], dict[str, Any]]:
    with tempfile.TemporaryDirectory(prefix="npass_aggregate_mmseqs_") as temp_name:
        temp_dir = Path(temp_name)
        future_fasta = temp_dir / "future.fasta"
        historical_fasta = temp_dir / "historical.fasta"
        hits_path = temp_dir / "hits.tsv"
        search_temp = temp_dir / "search_tmp"

        def write_fasta(path: Path, targets: set[str]) -> None:
            with path.open("wt", encoding="utf-8", newline="\n") as handle:
                for target in sorted(targets):
                    handle.write(f">{target}\n{sequences[target]}\n")

        write_fasta(future_fasta, future_targets)
        write_fasta(historical_fasta, historical_targets)
        command = [
            str(mmseqs_exe),
            "easy-search",
            str(future_fasta),
            str(historical_fasta),
            str(hits_path),
            str(search_temp),
            "--min-seq-id",
            "0.000000",
            "-c",
            "0.800000",
            "--cov-mode",
            "0",
            "--seq-id-mode",
            "0",
            "--filter-hits",
            "1",
            "-s",
            "7.500000",
            "--max-seqs",
            "1000000",
            "--format-output",
            "query,target,pident,alnlen,qcov,tcov,evalue,bits",
            "--threads",
            str(threads),
        ]
        result = subprocess.run(
            command,
            check=False,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        if result.returncode != 0 or not hits_path.is_file():
            raise RuntimeError(
                "Pinned MMseqs2 aggregate audit failed; "
                f"return code {result.returncode}, stderr tail {result.stderr[-500:]!r}"
            )
        maxima: dict[str, float] = {}
        raw_text = hits_path.read_text(encoding="utf-8", errors="strict").replace("\r", "")
        for line_number, line in enumerate(raw_text.splitlines(), start=1):
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) != 8:
                raise ValueError(
                    f"Malformed MMseqs2 aggregate-audit row at line {line_number}"
                )
            query, target, raw_pident = parts[0], parts[1], float(parts[2])
            if query not in future_targets or target not in historical_targets:
                raise ValueError("MMseqs2 returned an identifier outside the locked subsets")
            if not 0.0 <= raw_pident <= 100.0:
                raise ValueError(
                    "MMseqs2 pident is outside its documented percentage scale"
                )
            pident = raw_pident / 100.0
            maxima[query] = max(maxima.get(query, 0.0), pident)
        raw_hash = sha256(hits_path)
        version_result = subprocess.run(
            [str(mmseqs_exe), "version"],
            check=True,
            capture_output=True,
            text=True,
            encoding="utf-8",
            errors="replace",
        )
        audit = {
            "mmseqs_version": version_result.stdout.strip(),
            "parameters": {
                "minimum_sequence_identity": 0.0,
                "minimum_bidirectional_coverage": 0.8,
                "coverage_mode": 0,
                "sequence_identity_mode": 0,
                "filter_hits": 1,
                "sensitivity": 7.5,
                "maximum_sequences": 1000000,
                "threads": threads,
            },
            "raw_pident_unit": "percent",
            "reported_identity_unit": "fraction",
            "temporary_identifier_bearing_output_deleted": True,
            "temporary_hit_table_sha256": raw_hash,
            "future_target_count": len(future_targets),
            "historical_target_count": len(historical_targets),
            "future_target_count_with_detected_alignment": len(maxima),
            "future_target_count_without_detected_alignment": len(future_targets) - len(maxima),
            "interpretation": (
                "Detected-alignment conditional maxima under the pinned MMseqs2 "
                "heuristic; non-hits are censored and are not assigned zero identity."
            ),
        }
    return maxima, audit


def summarize_mmseqs(
    maxima: dict[str, float],
    scopes: dict[str, list[dict[str, str]]],
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    interpretation = (
        "Detected-alignment conditional MMseqs2 identity; missing alignments are "
        "reported as censored, not as identity zero."
    )
    for scope, rows in scopes.items():
        relation_observed = [
            maxima[row["uniprot_canonical_accession"]]
            for row in rows
            if row["uniprot_canonical_accession"] in maxima
        ]
        output.append(
            summary_row(
                family="mmseqs2_detected_alignment_identity",
                unit="relation_weighted",
                scope=scope,
                values=relation_observed,
                n_total=len(rows),
                n_no_detected_alignment=len(rows) - len(relation_observed),
                interpretation=interpretation,
            )
        )
        by_query: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            by_query[row["query_id"]].append(row)
        query_observed: list[float] = []
        no_detection = 0
        partial = 0
        all_detected = 0
        for query_rows in by_query.values():
            values = [
                maxima[row["uniprot_canonical_accession"]]
                for row in query_rows
                if row["uniprot_canonical_accession"] in maxima
            ]
            if not values:
                no_detection += 1
                continue
            query_observed.append(max(values))
            if len(values) == len(query_rows):
                all_detected += 1
            else:
                partial += 1
        output.append(
            summary_row(
                family="mmseqs2_detected_alignment_identity",
                unit="query_maximum",
                scope=scope,
                values=query_observed,
                n_total=len(by_query),
                n_no_detected_alignment=no_detection,
                n_partially_censored=partial,
                n_scope_units_with_all_relations_detected=all_detected,
                interpretation=interpretation,
            )
        )
    return output


def date_precision_audit(
    row_ledger: Path,
    pair_ledger: Path,
) -> dict[str, Any]:
    rows = read_tsv(row_ledger)
    non_day = [
        row for row in rows if row["row_eligibility_status"] == "excluded_non_day_precision"
    ]
    definitely_before = 0
    crossing = 0
    definitely_after = 0

    def bounds(value: str, precision: str) -> tuple[date, date]:
        if precision == "year":
            year = int(value)
            return date(year, 1, 1), date(year, 12, 31)
        if precision == "month":
            year, month = (int(item) for item in value.split("-"))
            return date(year, month, 1), date(
                year, month, calendar.monthrange(year, month)[1]
            )
        if precision == "day":
            parsed = date.fromisoformat(value)
            return parsed, parsed
        raise ValueError(f"Unsupported date precision: {precision}")

    for row in non_day:
        lower, upper = bounds(row["publication_date"], row["date_precision"])
        if upper <= CUTOFF:
            definitely_before += 1
        elif lower > CUTOFF:
            definitely_after += 1
        else:
            crossing += 1
    interval_eligible_rows = [
        row
        for row in rows
        if row["row_eligibility_status"] == "eligible_pre_cutoff"
        or (
            row["row_eligibility_status"] == "excluded_non_day_precision"
            and bounds(row["publication_date"], row["date_precision"])[1] <= CUTOFF
        )
    ]
    by_pair: dict[str, list[str]] = defaultdict(list)
    for row in interval_eligible_rows:
        by_pair[row["canonical_pair_key"]].append(row["evidence_tier_v1_1"])

    def tier_rank(value: str) -> int:
        if value.startswith("A_"):
            return 0
        if value.startswith("B_"):
            return 1
        raise ValueError(f"Unexpected strict evidence tier: {value!r}")

    interval_tier = {
        key: min(values, key=tier_rank) for key, values in by_pair.items()
    }
    pair_rows = read_tsv(pair_ledger)
    interval_vs_corrected = 0
    interval_vs_old = 0
    corrected_vs_old = 0
    corrected_changes_reverted = 0
    corrected_changes_persist = 0
    for row in pair_rows:
        key = row["canonical_pair_key"]
        candidate = interval_tier[key]
        corrected = row["corrected_best_strict_evidence_tier"]
        old = row["old_best_strict_evidence_tier"]
        interval_vs_corrected += candidate != corrected
        interval_vs_old += candidate != old
        corrected_vs_old += corrected != old
        corrected_changes_reverted += corrected != old and candidate == old
        corrected_changes_persist += corrected != old and candidate != old
    return {
        "locked_row_count": len(rows),
        "day_precision_eligible_row_count": sum(
            row["row_eligibility_status"] == "eligible_pre_cutoff" for row in rows
        ),
        "non_day_precision_row_count": len(non_day),
        "non_day_interval_definitely_before_or_on_cutoff": definitely_before,
        "non_day_interval_crosses_cutoff": crossing,
        "non_day_interval_definitely_after_cutoff": definitely_after,
        "interval_certain_eligible_row_count": len(interval_eligible_rows),
        "non_numeric_or_missing_reference_row_count": sum(
            row["row_eligibility_status"] == "excluded_not_numeric_pmid"
            for row in rows
        ),
        "historical_relation_count_under_interval_certain_policy": len(interval_tier),
        "corrected_vs_old_tier_change_count": corrected_vs_old,
        "interval_certain_vs_day_only_tier_change_count": interval_vs_corrected,
        "interval_certain_vs_old_tier_change_count": interval_vs_old,
        "day_only_changes_reverted_under_interval_certain_policy": corrected_changes_reverted,
        "day_only_changes_persisting_vs_old_under_interval_certain_policy": corrected_changes_persist,
        "interpretation": (
            "All non-day intervals are certainly no later than the cutoff. The "
            "remaining differences from the legacy tier cannot be attributed "
            "solely to after-cutoff removal."
        ),
    }


def reviewer_rows() -> list[dict[str, str]]:
    return [
        {
            "comment_id": "C1",
            "severity": "critical",
            "review_request": "Establish an executable external verification route, code and environment package, controlled row-level access, independent rerun, and rights matrix.",
            "current_evidence": "Local package machinery passes static checks, but formal packages are blocked by author metadata and every external release gate remains false.",
            "artifact_mapping": "scripts/reviewer_package_readiness_v3.json; plan/release_gate_state_v1.json; Table_S6_reproducibility_and_release.tsv",
            "action_class": "direct_text;external_human_gate",
            "v4_action": "Retain truthful access limits now; complete author metadata, rights clearance, software licensing, recipient/channel/encryption, and second-person rerun before submission.",
            "external_human_gate": "yes",
            "protocol_coverage": "external gate, not solvable by analysis protocol",
            "verification": "Reviewer is correct that no external route currently exists; journal policy requires a data-availability statement and access conditions, while public deposition is strongly encouraged rather than categorically mandatory.",
        },
        {
            "comment_id": "C2",
            "severity": "critical",
            "review_request": "Propagate exact-score tie uncertainty and separate operational from all-zero sequence queries.",
            "current_evidence": "The frozen realization is reproducible but many top-50 memberships are boundary-tie dependent; most sequence-query score vectors are structurally all zero.",
            "artifact_mapping": "Table_5_score_degeneracy_and_ties.tsv; corrective_prediction_ranks.tsv.gz",
            "action_class": "recompute_existing_data",
            "v4_action": "Report exact expected, best, and worst tie-aware Recall and NDCG; classify score-identifiable versus tie-dependent; stratify sequence queries.",
            "external_human_gate": "no",
            "protocol_coverage": "covered by Analysis A",
            "verification": "Numerical premise verified.",
        },
        {
            "comment_id": "C3",
            "severity": "critical",
            "review_request": "Separate temporal repair from date-precision exclusion using interval-censored scenarios.",
            "current_evidence": "All non-day date intervals end on or before the cutoff; admitting interval-certain rows preserves all historical relations and reverses most day-only tier changes.",
            "artifact_mapping": "restricted_ledger/historical_row_eligibility.tsv.gz; restricted_ledger/historical_pair_before_after.tsv.gz",
            "action_class": "recompute_existing_data",
            "v4_action": "Make interval-certain, day-only conservative, and explicitly bounded scenarios central; report tier, representation, rank, scope, and metric effects.",
            "external_human_gate": "no",
            "protocol_coverage": "covered by Analysis C",
            "verification": "Reviewer identified the right issue, but the audit strengthens it: no non-day interval crosses the cutoff.",
        },
        {
            "comment_id": "I1",
            "severity": "important",
            "review_request": "Keep the work retrospective and do not retrofit a new model on the outcome-visible endpoint.",
            "current_evidence": "The manuscript already labels the repair outcome-visible, author-run, non-independent, and post hoc.",
            "artifact_mapping": "manuscript_full_v3.md; frozen revision protocol",
            "action_class": "direct_text;external_future_design",
            "v4_action": "Preserve methodological-audit positioning; reserve validation or superiority claims for a preregistered unseen snapshot with an independent executor.",
            "external_human_gate": "future validation only",
            "protocol_coverage": "covered by protocol scientific status",
            "verification": "Reviewer recommendation is methodologically sound.",
        },
        {
            "comment_id": "I2",
            "severity": "important",
            "review_request": "Elevate source-document-component uncertainty and add finite-sample upper bounds for zero-hit scopes.",
            "current_evidence": "The endpoint has a concentrated document graph and component resampling materially widens some intervals.",
            "artifact_mapping": "Table_S4_pmid_document_dependence.tsv; source_concentration_aggregate_summary.json",
            "action_class": "recompute_existing_data;direct_text",
            "v4_action": "Co-report query and document-component intervals with separate estimands and one-sided finite-sample zero-hit bounds.",
            "external_human_gate": "no",
            "protocol_coverage": "covered by Analysis B",
            "verification": "Document, component, and largest-component counts verified.",
        },
        {
            "comment_id": "I3",
            "severity": "important",
            "review_request": "Test chemical normalization, Tier B weights, and unresolved-entity bounds.",
            "current_evidence": "Current structures are role-separated but not parent or tautomer normalized; Tier B uses a policy weight; unresolved candidates are a substantial fraction of the initial pool.",
            "artifact_mapping": "role_separated_compound_structure_audit.tsv.gz; Table_S5_frozen_unresolved_exclusions.tsv; unresolved aggregate audit",
            "action_class": "recompute_existing_data;external_human_mapping_limit",
            "v4_action": "Run fixed normalization and weight scenarios, compare available resolved versus unresolved aggregates, and provide best/worst endpoint bounds without inventing missing mappings.",
            "external_human_gate": "only if recovering missing mappings",
            "protocol_coverage": "covered by Analyses D, E, and F",
            "verification": "Reviewer numerical exclusion fraction verified; individual unresolved mappings may remain unidentified.",
        },
        {
            "comment_id": "I4",
            "severity": "important",
            "review_request": "Rename double-cold scope, merge identical masks, and report maximum structure and sequence similarity distributions.",
            "current_evidence": "The scope is project-defined exact-scaffold plus homology stress testing; the identity 0.50 and 0.70 masks are exactly identical.",
            "artifact_mapping": "scaffold_audit.tsv.gz; homology threshold ledgers; max_similarity_summary.tsv",
            "action_class": "direct_text;recompute_existing_data",
            "v4_action": "Use project-defined joint scaffold-homology cold scope, merge identical masks, and add aggregate maximum-similarity distributions.",
            "external_human_gate": "no",
            "protocol_coverage": "rename and mask merge covered; maximum-similarity output added by this audit",
            "verification": "Exact mask identity and aggregate similarity distributions verified.",
        },
        {
            "comment_id": "I5",
            "severity": "important",
            "review_request": "Show exhaustive-neighbour sensitivity alongside the locked approximation and sharpen positioning versus current methods.",
            "current_evidence": "Exhaustive scoring changes many full ranks and some top-50 memberships but minimally changes temporal Recall at 50.",
            "artifact_mapping": "Table_S2_top100_exhaustive_fidelity.tsv; DataSAIL reference; introduction and discussion",
            "action_class": "direct_text;reuse_existing_analysis",
            "v4_action": "Co-display exhaustive sensitivity, explain why the frozen approximation remains primary, and position against leakage benchmarks and modern prediction methods without adding a new endpoint-tuned model.",
            "external_human_gate": "no",
            "protocol_coverage": "existing analysis available; explicit reporting elevation and comparator text were omitted from the frozen analysis list",
            "verification": "Rank-change, membership-change, and headline-metric difference verified.",
        },
        {
            "comment_id": "S1",
            "severity": "suggestion",
            "review_request": "Add a claim-evidence-unsupported-interpretation table.",
            "current_evidence": "Claim boundaries are distributed across the manuscript rather than presented as a reader-facing matrix.",
            "artifact_mapping": "discussion; data availability; new reader-facing table",
            "action_class": "direct_text",
            "v4_action": "Add a compact claim-boundary table.",
            "external_human_gate": "no",
            "protocol_coverage": "not an analysis item",
            "verification": "Useful and low-risk.",
        },
        {
            "comment_id": "S2",
            "severity": "suggestion",
            "review_request": "Visually combine query, component, and tie-aware uncertainty and label zero bootstrap intervals as empirically degenerate.",
            "current_evidence": "Current figure and table separate some dependence and tie diagnostics.",
            "artifact_mapping": "Figure_3; Table_4",
            "action_class": "direct_text;figure_table_revision",
            "v4_action": "Update the figure and table after Analyses A and B.",
            "external_human_gate": "no",
            "protocol_coverage": "dependent on Analyses A and B",
            "verification": "Recommended presentation change.",
        },
        {
            "comment_id": "S3",
            "severity": "suggestion",
            "review_request": "State recommended and prohibited uses.",
            "current_evidence": "Limitations are present, but operational use boundaries are not consolidated.",
            "artifact_mapping": "discussion; conclusion; benchmark card",
            "action_class": "direct_text",
            "v4_action": "Permit database-quality and retrieval-audit use; prohibit direct experimental confirmation, activity-probability, or clinical-priority interpretation.",
            "external_human_gate": "no",
            "protocol_coverage": "not an analysis item",
            "verification": "Recommended.",
        },
        {
            "comment_id": "S4",
            "severity": "suggestion",
            "review_request": "Compress abstract and keep the article within journal word limits.",
            "current_evidence": "Locked validation reports an abstract within the limit and a pre-reference manuscript close to the article limit.",
            "artifact_mapping": "manuscript_validation_v3.json; abstract; full manuscript",
            "action_class": "direct_text",
            "v4_action": "Do not repair a nonexistent abstract violation; compress proactively because new analyses will add text.",
            "external_human_gate": "no",
            "protocol_coverage": "not an analysis item",
            "verification": "Reviewer word counts are not reproducible; audited counts are 234 and 7,871.",
        },
        {
            "comment_id": "S5",
            "severity": "suggestion",
            "review_request": "Unify manuscript and supplementary version dates.",
            "current_evidence": "The current release artifacts carry adjacent dates.",
            "artifact_mapping": "manuscript front matter; supplementary workbook metadata",
            "action_class": "direct_text;build_validation",
            "v4_action": "Use one release date across all v4 outputs and assert it in validation.",
            "external_human_gate": "no",
            "protocol_coverage": "not an analysis item",
            "verification": "Reviewer observation verified.",
        },
    ]


def markdown_matrix(rows: list[dict[str, str]], receipt: dict[str, Any]) -> str:
    lines = [
        "# Reviewer-comment disposition and minimal v4 scope",
        "",
        "This package is an aggregate-only audit. It contains no compound, target, pair, query, or source-document identifiers.",
        "",
        "> Correction notice: the initial MMseqs2 summary interpreted percentage-scale pident as a fraction and clipped detected values to 1.0. That build is superseded. The current table divides pident by 100 and is guarded by unit-specific validation.",
        "",
        "## Verified corrections to the review report",
        "",
        f"- The locked manuscript validator reports **{receipt['word_counts']['abstract']} abstract words**, not 258.",
        f"- It reports **{receipt['word_counts']['pre_reference_manuscript']} words before references**, not approximately 7,979. The margin is nevertheless too small for the planned revision.",
        f"- The identity 0.50 and 0.70 joint masks are exactly identical: **{receipt['scope_counts']['joint_scaffold_homology_0_50_0_70_identical']['relations']} relations, {receipt['scope_counts']['joint_scaffold_homology_0_50_0_70_identical']['queries']} queries, and {receipt['scope_counts']['joint_scaffold_homology_0_50_0_70_identical']['targets']} targets**.",
        f"- All **{receipt['date_precision_audit']['non_day_precision_row_count']} non-day-precision rows** have interval upper bounds on or before the cutoff; none crosses it.",
        "",
        "## Comment matrix",
        "",
        "| ID | Severity | Action class | Protocol coverage | External or human gate | v4 disposition |",
        "|---|---|---|---|---|---|",
    ]
    for row in rows:
        values = [
            row["comment_id"],
            row["severity"],
            row["action_class"],
            row["protocol_coverage"],
            row["external_human_gate"],
            row["v4_action"],
        ]
        values = [value.replace("|", "/").replace("\n", " ") for value in values]
        lines.append("| " + " | ".join(values) + " |")
    lines.extend(
        [
            "",
            "## Minimal publishable v4 scope",
            "",
            "The smallest scientifically defensible revision retains the frozen endpoint and baselines and does **not** add an outcome-tuned model.",
            "",
            "1. Complete exact tie-aware expected and best-to-worst retrieval estimates, including the operational versus all-zero sequence-query split.",
            "2. Make interval-certain date handling the central temporal sensitivity, with day-only handling clearly labeled as a conservative date-precision policy.",
            "3. Co-report query and source-document-component uncertainty plus finite-sample zero-hit upper bounds.",
            "4. Complete fixed chemistry-normalization, Tier B weighting, and unresolved-endpoint bound analyses.",
            "5. Rename the joint scaffold-homology scope, merge the identical 0.50 and 0.70 masks, and report the aggregate maximum-similarity distributions supplied here.",
            "6. Elevate the already-computed exhaustive-neighbour comparison to a co-displayed sensitivity rather than claiming ranking equivalence.",
            "7. Add the claim-boundary table, intended/prohibited-use statements, modern-method positioning, unified dates, and proactive word-count compression.",
            "",
            "## Submission gates that computation cannot close",
            "",
            "- Supply author metadata and choose a software license.",
            "- Complete a source-by-source rights matrix and obtain any required provider confirmation.",
            "- Configure the authorized recipient, controlled channel, and encryption for restricted reviewer access.",
            "- Have a genuinely independent person rerun the locked workflow and sign the consistency report.",
            "",
            "These are submission-governance gates. They should not be represented as completed merely because a local package builder exists.",
            "",
            "## Maximum-similarity interpretation",
            "",
            "Morgan and native 3-mer cosine summaries are complete under the frozen feature definitions. MMseqs2 rows are conditional on detected alignments under the pinned heuristic. Missing alignments are explicitly censored and are never converted to identity zero.",
            "",
        ]
    )
    return "\n".join(lines)


def build(args: argparse.Namespace) -> None:
    script_path = Path(__file__).resolve()
    output_dir = script_path.parents[1]
    project_root = script_path.parents[2]
    reviewer_report = Path(args.review_report).resolve()
    if not reviewer_report.is_file():
        raise FileNotFoundError("Reviewer report attachment was not found")

    paths = {
        "protocol": project_root
        / "manuscript_molecular_diversity_v3_20260728/plan/revision_analysis_protocol_v4_20260729.md",
        "manuscript": project_root
        / "manuscript_molecular_diversity_v3_20260728/manuscript_full_v3.md",
        "manuscript_validation": project_root
        / "manuscript_molecular_diversity_v3_20260728/manuscript_validation_v3.json",
        "package_readiness": project_root
        / "manuscript_molecular_diversity_v3_20260728/scripts/reviewer_package_readiness_v3.json",
        "release_gate_state": project_root
        / "manuscript_molecular_diversity_v3_20260728/plan/release_gate_state_v1.json",
        "table_ties": project_root
        / "manuscript_molecular_diversity_v3_20260728/tables/Table_5_score_degeneracy_and_ties.tsv",
        "table_exhaustive": project_root
        / "manuscript_molecular_diversity_v3_20260728/tables/Table_S2_top100_exhaustive_fidelity.tsv",
        "table_dependence": project_root
        / "manuscript_molecular_diversity_v3_20260728/tables/Table_S4_pmid_document_dependence.tsv",
        "table_unresolved": project_root
        / "manuscript_molecular_diversity_v3_20260728/tables/Table_S5_frozen_unresolved_exclusions.tsv",
        "table_release": project_root
        / "manuscript_molecular_diversity_v3_20260728/tables/Table_S6_reproducibility_and_release.tsv",
        "source_concentration": project_root
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/audit/source_concentration_v1/source_concentration_aggregate_summary.json",
        "rebuild_summary": project_root
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/audit/asof_rebuild_summary.json",
        "unresolved_summary": project_root
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/audit/entity_unresolved_v1/unresolved_reason_aggregate_summary.json",
        "historical_row_eligibility": project_root
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/restricted_ledger/historical_row_eligibility.tsv.gz",
        "historical_pair_before_after": project_root
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/restricted_ledger/historical_pair_before_after.tsv.gz",
        "role_separated_structure_audit": project_root
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/restricted_ledger/role_separated_compound_structure_audit.tsv.gz",
        "historical_compounds": project_root
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/scoring_inputs/historical_compounds.tsv.gz",
        "query_compounds": project_root
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/scoring_inputs/query_compounds.tsv.gz",
        "historical_pairs": project_root
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/scoring_inputs/historical_pairs.tsv.gz",
        "candidate_targets": project_root
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/scoring_inputs/candidate_targets.tsv.gz",
        "candidate_sequences": project_root
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/scoring_inputs/candidate_sequences.fasta",
        "evaluation_pairs": project_root
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/evaluation_inputs/evaluation_pairs.tsv.gz",
        "scaffold_audit": project_root
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/evaluation_inputs/scaffold_audit.tsv.gz",
        "homology_0_30": project_root
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/evaluation_inputs/homology_0_30.tsv.gz",
        "homology_0_50": project_root
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/evaluation_inputs/homology_0_50.tsv.gz",
        "homology_0_70": project_root
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/evaluation_inputs/homology_0_70.tsv.gz",
        "score_feature_definition": project_root / "scripts/pu_retrieval_scores.py",
        "mmseqs_launcher": project_root
        / "tools/mmseqs2/18-8cc5c/mmseqs/mmseqs.bat",
        "mmseqs_executable": project_root
        / "tools/mmseqs2/18-8cc5c/mmseqs/bin/mmseqs.exe",
    }
    missing = [label for label, path in paths.items() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing locked inputs: {', '.join(missing)}")
    if sha256(paths["protocol"]) != EXPECTED_PROTOCOL_SHA256:
        raise AssertionError("Frozen v4 protocol hash mismatch")
    if sha256(reviewer_report) != EXPECTED_REVIEW_SHA256:
        raise AssertionError("Reviewer report hash mismatch")

    generated_names = [
        "reviewer_comment_matrix.tsv",
        "reviewer_comment_matrix.md",
        "max_similarity_summary.tsv",
        "correction_record.json",
        "input_hashes.json",
        "execution_receipt.json",
        "manifest.json",
    ]
    existing = [name for name in generated_names if (output_dir / name).exists()]
    if existing and not args.force:
        raise FileExistsError(
            "Refusing to overwrite existing outputs without --force: "
            + ", ".join(existing)
        )

    evaluation_rows = read_tsv(paths["evaluation_pairs"])
    scaffold_rows = read_tsv(paths["scaffold_audit"])
    homology_rows = {
        "0_30": read_tsv(paths["homology_0_30"]),
        "0_50": read_tsv(paths["homology_0_50"]),
        "0_70": read_tsv(paths["homology_0_70"]),
    }
    scopes, scope_counts = build_scopes(
        evaluation_rows, scaffold_rows, homology_rows
    )

    historical_compounds = read_tsv(paths["historical_compounds"])
    query_compounds = read_tsv(paths["query_compounds"])
    historical_pairs = read_tsv(paths["historical_pairs"])
    candidate_targets = [
        row["uniprot_canonical_accession"]
        for row in read_tsv(paths["candidate_targets"])
    ]
    sequences = read_fasta(paths["candidate_sequences"])
    historical_targets = {
        row["uniprot_canonical_accession"] for row in historical_pairs
    }
    future_targets = {
        row["uniprot_canonical_accession"] for row in evaluation_rows
    }

    similarity_rows: list[dict[str, Any]] = []
    morgan = morgan_maxima(historical_compounds, query_compounds)
    similarity_rows.extend(
        summarize_complete_similarity(
            family="morgan_radius2_2048_tanimoto",
            value_by_entity=morgan,
            entity_field="inchikey_full",
            scopes=scopes,
            query_reducer="constant",
            interpretation=(
                "Complete maximum query-to-historical-compound similarity under "
                "the frozen Morgan radius-2, 2048-bit representation."
            ),
        )
    )
    kmer = sequence_3mer_maxima(
        candidate_targets, sequences, historical_targets
    )
    similarity_rows.extend(
        summarize_complete_similarity(
            family="native_sequence_3mer_tfidf_cosine",
            value_by_entity=kmer,
            entity_field="uniprot_canonical_accession",
            scopes=scopes,
            query_reducer="max",
            interpretation=(
                "Complete maximum endpoint-target-to-historical-target cosine "
                "similarity under the frozen native 3-mer feature definition."
            ),
        )
    )
    mmseqs_maxima, mmseqs_audit = run_mmseqs_detected_identity(
        mmseqs_exe=paths["mmseqs_executable"],
        sequences=sequences,
        future_targets=future_targets,
        historical_targets=historical_targets,
        threads=args.threads,
    )
    similarity_rows.extend(summarize_mmseqs(mmseqs_maxima, scopes))

    similarity_fields = [
        "similarity_family",
        "analysis_unit",
        "scope",
        "n_total",
        "n_observed",
        "n_no_detected_alignment",
        "n_partially_censored",
        "n_scope_units_with_all_relations_detected",
        "n_at_upper_bound",
        "min",
        "q1",
        "median",
        "q3",
        "max",
        "mean",
        "interpretation",
    ]
    write_tsv(
        output_dir / "max_similarity_summary.tsv",
        similarity_fields,
        similarity_rows,
    )

    rows = reviewer_rows()
    matrix_fields = [
        "comment_id",
        "severity",
        "review_request",
        "current_evidence",
        "artifact_mapping",
        "action_class",
        "v4_action",
        "external_human_gate",
        "protocol_coverage",
        "verification",
    ]
    write_tsv(output_dir / "reviewer_comment_matrix.tsv", matrix_fields, rows)

    validation = read_json(paths["manuscript_validation"])
    abstract_check = next(
        item["check"]
        for item in validation["checks"]
        if item["check"].startswith("abstract word count is journal-compliant:")
    )
    manuscript_check = next(
        item["check"]
        for item in validation["checks"]
        if item["check"].startswith(
            "pre-reference manuscript word count is within Research Article limit:"
        )
    )
    word_counts = {
        "abstract": int(abstract_check.rsplit(":", 1)[1].strip()),
        "pre_reference_manuscript": int(
            manuscript_check.rsplit(":", 1)[1].strip()
        ),
    }
    source_summary = read_json(paths["source_concentration"])
    endpoint_source = next(
        item for item in source_summary["cohorts"] if item["cohort"] == "endpoint"
    )
    date_audit = date_precision_audit(
        paths["historical_row_eligibility"], paths["historical_pair_before_after"]
    )
    correction_record = {
        "schema_version": "revision_review_matrix_correction_record_v1",
        "status": "CORRECTED_AND_SUPERSEDED",
        "recorded_at_utc": utc_now(),
        "superseded_manifest_sha256": (
            "ad5a9692e93049004e9d0436e70ee1122ef53463e5667bdbe4f54ff9461897f1"
        ),
        "affected_artifact": "max_similarity_summary.tsv",
        "affected_rows": (
            "The eight MMseqs2 detected-alignment summary rows in the initial build."
        ),
        "error": (
            "MMseqs2 pident was emitted on a 0-to-100 percentage scale. The initial "
            "build treated it as a 0-to-1 fraction and the generic distribution "
            "guard clipped every detected value above one to one."
        ),
        "correction": (
            "Convert raw pident to fractional identity by dividing by 100 before "
            "aggregation. Assert the raw percentage range, the final fractional "
            "range, the detected 0.30-scope value, and the merged 0.50/0.70 range."
        ),
        "unaffected_artifacts": [
            "reviewer_comment_matrix.tsv",
            "reviewer_comment_matrix.md",
            "Morgan similarity rows",
            "native 3-mer similarity rows",
            "date-precision audit",
            "scope and source-dependence counts",
        ],
        "identifier_bearing_outputs_retained": 0,
    }
    write_json(output_dir / "correction_record.json", correction_record)
    receipt: dict[str, Any] = {
        "schema_version": "revision_review_matrix_execution_receipt_v1",
        "status": "PASS",
        "created_at_utc": utc_now(),
        "claim_boundary": (
            "Aggregate audit only; no external validation, biological validation, "
            "or evidence that an undetected alignment has identity zero."
        ),
        "reviewer_comment_count": len(rows),
        "word_counts": word_counts,
        "scope_counts": scope_counts,
        "source_dependence": {
            "source_document_count": endpoint_source["unique_source_document_count"],
            "query_source_component_count": endpoint_source[
                "query_source_component_summary"
            ]["component_count"],
            "largest_query_component_count": endpoint_source[
                "query_source_component_summary"
            ]["largest_component_left_node_count"],
            "largest_query_component_fraction": endpoint_source[
                "query_source_component_summary"
            ]["largest_component_left_node_fraction"],
        },
        "date_precision_audit": date_audit,
        "mmseqs_detected_alignment_audit": mmseqs_audit,
        "software": {
            "python": __import__("sys").version.split()[0],
            "numpy": np.__version__,
            "rdkit": rdBase.rdkitVersion,
            "scikit_learn": sklearn_version,
        },
        "aggregate_similarity_row_count": len(similarity_rows),
        "identifier_bearing_outputs_retained": 0,
        "correction_history": correction_record,
    }
    write_json(output_dir / "execution_receipt.json", receipt)
    (output_dir / "reviewer_comment_matrix.md").write_text(
        markdown_matrix(rows, receipt), encoding="utf-8", newline="\n"
    )

    input_hash_entries: dict[str, Any] = {
        "reviewer_report_attachment": {
            "location_class": "external_attachment",
            "sha256": sha256(reviewer_report),
        }
    }
    for label, path in paths.items():
        input_hash_entries[label] = {
            "location_class": "project_relative",
            "path": path.relative_to(project_root).as_posix(),
            "sha256": sha256(path),
        }
    input_hashes = {
        "schema_version": "revision_review_matrix_input_hashes_v1",
        "frozen_protocol_sha256": EXPECTED_PROTOCOL_SHA256,
        "inputs": input_hash_entries,
    }
    write_json(output_dir / "input_hashes.json", input_hashes)

    output_files = [
        "reviewer_comment_matrix.tsv",
        "reviewer_comment_matrix.md",
        "max_similarity_summary.tsv",
        "correction_record.json",
        "input_hashes.json",
        "execution_receipt.json",
    ]
    validator_path = output_dir / "scripts/validate_review_matrix.py"
    manifest = {
        "schema_version": "revision_review_matrix_manifest_v1",
        "package_id": "revision_review_matrix_v1_20260729",
        "status": "BUILD_COMPLETE_CORRECTED",
        "aggregate_only": True,
        "identifier_bearing_outputs_retained": 0,
        "frozen_protocol_sha256": EXPECTED_PROTOCOL_SHA256,
        "scripts": {
            "scripts/build_review_matrix.py": sha256(script_path),
            "scripts/validate_review_matrix.py": sha256(validator_path),
        },
        "outputs": {
            name: {"sha256": sha256(output_dir / name)}
            for name in output_files
        },
        "validation_report": (
            "validation_report.json is generated after this immutable manifest "
            "is written and therefore is not self-listed."
        ),
    }
    write_json(output_dir / "manifest.json", manifest)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--review-report", required=True)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    if args.threads < 1:
        parser.error("--threads must be positive")
    build(args)


if __name__ == "__main__":
    main()
