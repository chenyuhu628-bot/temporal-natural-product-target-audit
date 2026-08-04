"""Build a restricted, author-run strict A/B successor input bundle.

This builder reads only the explicitly named C31-primary interim/processed
artifacts, never legacy results, outer ledgers, raw data, or BindingDB payloads.
It writes a new isolated bundle with the smallest schemas required by the
author-run score/evaluation runner. The bundle is restricted to internal use;
it is not a release package, does not establish third-party data rights, and
does not claim human gate completion, endpoint access control, or independence.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


PROTOCOL_ID = "npass_strict_ab_doublecold_successor_v1_20260719"
RUN_ID = "npass_strict_ab_doublecold_author_run_v1_20260719"
ROOT = Path(__file__).resolve().parents[2]
WORKSPACE = ROOT

SOURCE_PATHS = {
    "historical_pairs": WORKSPACE
    / "data"
    / "interim"
    / "cold_start_inputs_v1_1_pmid_verified_C31_primary"
    / "strict_temporal_training_pairs.tsv.gz",
    "future_pairs": WORKSPACE
    / "data"
    / "interim"
    / "cold_start_inputs_v1_1_pmid_verified_C31_primary"
    / "strict_temporal_future_pairs.tsv.gz",
    "compounds": WORKSPACE
    / "data"
    / "interim"
    / "cold_start_inputs_v1_1_pmid_verified_C31_primary"
    / "compounds_for_rdkit_with_validated_inchi_repairs.tsv.gz",
    "scaffold_audit": WORKSPACE
    / "data"
    / "interim"
    / "scaffold_coldness_audit_v1_1_pmid_verified_C31_primary_with_validated_inchi_repair"
    / "future_pair_scaffold_coldness_audit.tsv.gz",
    "homology_0_30": WORKSPACE
    / "data"
    / "interim"
    / "mmseqs2_target_coldness_v1_1_pmid_verified_C31_primary"
    / "identity_0_30"
    / "future_target_coldness_audit.tsv.gz",
    "homology_0_50": WORKSPACE
    / "data"
    / "interim"
    / "mmseqs2_target_coldness_v1_1_pmid_verified_C31_primary"
    / "identity_0_50"
    / "future_target_coldness_audit.tsv.gz",
    "homology_0_70": WORKSPACE
    / "data"
    / "interim"
    / "mmseqs2_target_coldness_v1_1_pmid_verified_C31_primary"
    / "identity_0_70"
    / "future_target_coldness_audit.tsv.gz",
    "candidate_targets": WORKSPACE
    / "data"
    / "processed"
    / "chembl31_human_target_catalogue_v1"
    / "chembl31_human_single_protein_targets_v1.tsv.gz",
    "candidate_sequences": WORKSPACE
    / "data"
    / "processed"
    / "chembl31_human_target_catalogue_v1"
    / "chembl31_human_single_protein_targets_v1.fasta",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv_gz(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing TSV header: {path}")
        return list(reader.fieldnames), list(reader)


def write_tsv_gz(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite output: {path}")
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(text, fieldnames=fields, delimiter="\t", extrasaction="raise")
                writer.writeheader()
                writer.writerows(rows)


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    accession = ""
    pieces: list[str] = []
    for raw in path.read_text(encoding="ascii").splitlines():
        if raw.startswith(">"):
            if accession:
                if accession in records:
                    raise ValueError(f"Duplicate FASTA accession: {accession}")
                records[accession] = "".join(pieces)
            accession = raw[1:].split("|")[0].strip()
            pieces = []
        else:
            pieces.append(raw.strip())
    if accession:
        if accession in records:
            raise ValueError(f"Duplicate FASTA accession: {accession}")
        records[accession] = "".join(pieces)
    return records


def write_fasta(path: Path, records: dict[str, str]) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite output: {path}")
    with path.open("wt", encoding="ascii", newline="\n") as handle:
        for accession in sorted(records):
            sequence = records[accession]
            if not sequence:
                raise ValueError(f"Empty sequence: {accession}")
            handle.write(f">{accession}\n")
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")


def require_fields(fields: list[str], required: set[str], label: str) -> None:
    missing = sorted(required.difference(fields))
    if missing:
        raise ValueError(f"{label} lacks fields: {missing}")


def assert_unique(rows: list[dict[str, str]], key_fields: tuple[str, ...], label: str) -> None:
    observed: set[tuple[str, ...]] = set()
    for row in rows:
        key = tuple(row[field] for field in key_fields)
        if not all(key):
            raise ValueError(f"{label} has an empty key: {key_fields}")
        if key in observed:
            raise ValueError(f"{label} has a duplicate key: {key_fields}")
        observed.add(key)


def best_tier(value: str, label: str) -> str:
    tiers = {item.strip() for item in value.split(";") if item.strip()}
    allowed = {"A_affinity_candidate", "B_quantitative_functional_candidate"}
    if not tiers or not tiers.issubset(allowed):
        raise ValueError(f"{label} contains a non-strict tier set: {value!r}")
    return "A_affinity_candidate" if "A_affinity_candidate" in tiers else "B_quantitative_functional_candidate"


def assert_new_isolated_output(path: Path) -> None:
    resolved = path.resolve()
    if resolved.exists():
        raise FileExistsError(f"Output root already exists: {resolved}")
    if ROOT in (resolved, *resolved.parents):
        raise ValueError("Restricted successor input bundle must be outside the protocol folder")
    for blocked in (WORKSPACE / "data", WORKSPACE / "results", WORKSPACE / "manifests"):
        try:
            resolved.relative_to(blocked.resolve())
        except ValueError:
            continue
        raise ValueError(f"Output root may not be inside legacy tree: {blocked}")


def write_json(path: Path, payload: dict[str, Any]) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite output: {path}")
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> int:
    global WORKSPACE, SOURCE_PATHS
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=ROOT,
        help="isolated reconstruction workspace containing data/",
    )
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()
    original_workspace = WORKSPACE
    WORKSPACE = args.project_root.resolve()
    SOURCE_PATHS = {
        name: WORKSPACE / path.relative_to(original_workspace)
        for name, path in SOURCE_PATHS.items()
    }
    output_root = args.output_root.resolve()
    assert_new_isolated_output(output_root)
    for label, path in SOURCE_PATHS.items():
        if not path.is_file():
            raise FileNotFoundError(f"Required authorized source is absent: {label} ({path})")

    history_fields, history_source = read_tsv_gz(SOURCE_PATHS["historical_pairs"])
    future_fields, future_source = read_tsv_gz(SOURCE_PATHS["future_pairs"])
    compound_fields, compound_source = read_tsv_gz(SOURCE_PATHS["compounds"])
    scaffold_fields, scaffold_source = read_tsv_gz(SOURCE_PATHS["scaffold_audit"])
    target_fields, target_source = read_tsv_gz(SOURCE_PATHS["candidate_targets"])
    homology_sources = {
        threshold: read_tsv_gz(SOURCE_PATHS[f"homology_{threshold}"])
        for threshold in ("0_30", "0_50", "0_70")
    }
    sequences = read_fasta(SOURCE_PATHS["candidate_sequences"])

    require_fields(
        history_fields,
        {
            "canonical_pair_key",
            "inchikey_full",
            "uniprot_canonical_accession",
            "primary_evidence_tiers",
            "decision",
            "unrecorded_pair_policy",
        },
        "historical source",
    )
    require_fields(
        future_fields,
        {
            "canonical_pair_key",
            "inchikey_full",
            "uniprot_canonical_accession",
            "future_primary_evidence_tiers",
            "future_decision",
            "leakage_gate_decision",
            "negative_label_emitted",
        },
        "future source",
    )
    require_fields(compound_fields, {"inchikey_full", "representative_smiles"}, "compound source")
    require_fields(
        scaffold_fields,
        {
            "canonical_pair_key",
            "audit_scaffold_cold_under_selected_policy",
            "audit_outcome",
            "audit_eligibility_or_exclusion_reason",
        },
        "scaffold audit source",
    )
    require_fields(target_fields, {"uniprot_accession"}, "candidate-target source")
    for threshold, (fields, _) in homology_sources.items():
        require_fields(
            fields,
            {
                "uniprot_canonical_accession",
                "is_future_target_homology_cold_candidate",
                "future_target_coldness_status",
            },
            f"homology {threshold} source",
        )

    if len(history_source) != 4990 or len(future_source) != 358:
        raise ValueError("The locked strict A/B historical/future cardinalities are not 4,990/358")
    assert_unique(history_source, ("canonical_pair_key",), "historical source")
    assert_unique(history_source, ("inchikey_full", "uniprot_canonical_accession"), "historical source")
    assert_unique(future_source, ("canonical_pair_key",), "future source")
    assert_unique(future_source, ("inchikey_full", "uniprot_canonical_accession"), "future source")
    history_keys = {row["canonical_pair_key"] for row in history_source}
    future_keys = {row["canonical_pair_key"] for row in future_source}
    history_tuples = {(row["inchikey_full"], row["uniprot_canonical_accession"]) for row in history_source}
    future_tuples = {(row["inchikey_full"], row["uniprot_canonical_accession"]) for row in future_source}
    if history_keys.intersection(future_keys) or history_tuples.intersection(future_tuples):
        raise ValueError("Historical and future strict A/B source pairs overlap")
    if {row["decision"] for row in history_source} != {"strict_pre_cutoff_training_candidate"}:
        raise ValueError("Historical source decision differs from strict pre-cutoff contract")
    if {row["future_decision"] for row in future_source} != {"strict_post_cutoff_future_candidate"}:
        raise ValueError("Future source decision differs from strict post-cutoff contract")
    if {row["leakage_gate_decision"] for row in future_source} != {
        "include_primary_future_candidate_after_C31_leakage_screen"
    }:
        raise ValueError("Future source fails the C31 primary leakage-gate contract")
    if {row["negative_label_emitted"].casefold() for row in future_source} != {"false"}:
        raise ValueError("Future source includes an emitted negative label")

    history = [
        {
            "canonical_pair_key": row["canonical_pair_key"],
            "inchikey_full": row["inchikey_full"],
            "uniprot_canonical_accession": row["uniprot_canonical_accession"],
            "best_strict_evidence_tier": best_tier(row["primary_evidence_tiers"], "historical source"),
            "decision": "strict_pre_cutoff_training_candidate",
            "unrecorded_pair_policy": "unlabeled_not_negative",
        }
        for row in history_source
    ]
    future_tier_by_pair = {
        row["canonical_pair_key"]: best_tier(row["future_primary_evidence_tiers"], "future source")
        for row in future_source
    }
    ordered_future_compounds = sorted({row["inchikey_full"] for row in future_source})
    if len(ordered_future_compounds) != 222:
        raise ValueError(f"The locked strict A/B future query count is not 222: {len(ordered_future_compounds)}")
    query_id_by_compound = {
        compound: f"query_{index:04d}" for index, compound in enumerate(ordered_future_compounds, start=1)
    }
    scoring_queries = [
        {"query_id": query_id_by_compound[compound], "inchikey_full": compound}
        for compound in ordered_future_compounds
    ]
    endpoint = [
        {
            "canonical_pair_key": row["canonical_pair_key"],
            "query_id": query_id_by_compound[row["inchikey_full"]],
            "inchikey_full": row["inchikey_full"],
            "uniprot_canonical_accession": row["uniprot_canonical_accession"],
            "best_strict_evidence_tier": future_tier_by_pair[row["canonical_pair_key"]],
            "decision": "strict_post_cutoff_future_candidate",
            "c31_leakage_gate_status": "pass_no_historical_activity",
        }
        for row in future_source
    ]

    assert_unique(compound_source, ("inchikey_full",), "compound source")
    required_compounds = {row["inchikey_full"] for row in history}.union(ordered_future_compounds)
    compound_by_key = {row["inchikey_full"]: row for row in compound_source}
    missing_compounds = sorted(required_compounds.difference(compound_by_key))
    if missing_compounds:
        raise ValueError(f"Required compound structures are missing: {missing_compounds[:5]}")
    compounds = [
        {
            "inchikey_full": key,
            "representative_smiles": compound_by_key[key]["representative_smiles"],
        }
        for key in sorted(required_compounds)
    ]
    if any(not row["representative_smiles"].strip() for row in compounds):
        raise ValueError("A required successor compound has no representative SMILES")

    assert_unique(target_source, ("uniprot_accession",), "candidate-target source")
    candidate_target_ids = sorted(row["uniprot_accession"] for row in target_source)
    if len(candidate_target_ids) != 4123 or set(candidate_target_ids) != set(sequences):
        raise ValueError("C31 candidate target TSV/FASTA contract is not exactly 4,123 matching accessions")
    missing_history_targets = sorted(
        {row["uniprot_canonical_accession"] for row in history}.difference(candidate_target_ids)
    )
    if missing_history_targets:
        raise ValueError(f"Historical targets absent from C31 candidate universe: {missing_history_targets[:5]}")
    candidate_targets = [{"uniprot_canonical_accession": item} for item in candidate_target_ids]
    candidate_sequences = {item: sequences[item] for item in candidate_target_ids}

    scaffold_by_pair = {row["canonical_pair_key"]: row for row in scaffold_source}
    if set(scaffold_by_pair) != future_keys or len(scaffold_by_pair) != len(scaffold_source):
        raise ValueError("Scaffold audit does not have the exact future endpoint keyset")
    scaffold = [
        {
            "canonical_pair_key": pair_key,
            "audit_scaffold_cold_under_selected_policy": scaffold_by_pair[pair_key][
                "audit_scaffold_cold_under_selected_policy"
            ],
            "audit_outcome": scaffold_by_pair[pair_key]["audit_outcome"],
            "audit_eligibility_or_exclusion_reason": scaffold_by_pair[pair_key][
                "audit_eligibility_or_exclusion_reason"
            ],
        }
        for pair_key in sorted(future_keys)
    ]

    endpoint_targets = {row["uniprot_canonical_accession"] for row in endpoint}
    homology: dict[str, list[dict[str, str]]] = {}
    for threshold, (_, rows) in homology_sources.items():
        assert_unique(rows, ("uniprot_canonical_accession",), f"homology {threshold} source")
        source_by_target = {row["uniprot_canonical_accession"]: row for row in rows}
        if set(source_by_target) != endpoint_targets:
            raise ValueError(f"Homology {threshold} source does not have the exact future-target keyset")
        homology[threshold] = [
            {
                "uniprot_canonical_accession": target,
                "is_future_target_homology_cold_candidate": source_by_target[target][
                    "is_future_target_homology_cold_candidate"
                ],
                "future_target_coldness_status": source_by_target[target]["future_target_coldness_status"],
            }
            for target in sorted(endpoint_targets)
        ]

    scoring_dir = output_root / "scoring_inputs"
    evaluation_dir = output_root / "evaluation_inputs"
    metadata_dir = output_root / "metadata"
    scoring_dir.mkdir(parents=True, exist_ok=False)
    evaluation_dir.mkdir(parents=True, exist_ok=False)
    metadata_dir.mkdir(parents=True, exist_ok=False)

    scoring_paths = {
        "historical_pairs": scoring_dir / "historical_pairs.tsv.gz",
        "scoring_queries": scoring_dir / "scoring_queries.tsv.gz",
        "compounds": scoring_dir / "compounds.tsv.gz",
        "candidate_targets": scoring_dir / "candidate_targets.tsv.gz",
        "candidate_sequences": scoring_dir / "candidate_sequences.fasta",
    }
    write_tsv_gz(scoring_paths["historical_pairs"], list(history[0]), history)
    write_tsv_gz(scoring_paths["scoring_queries"], list(scoring_queries[0]), scoring_queries)
    write_tsv_gz(scoring_paths["compounds"], list(compounds[0]), compounds)
    write_tsv_gz(scoring_paths["candidate_targets"], list(candidate_targets[0]), candidate_targets)
    write_fasta(scoring_paths["candidate_sequences"], candidate_sequences)

    evaluation_paths = {
        "endpoint": evaluation_dir / "evaluation_pairs.tsv.gz",
        "scaffold_audit": evaluation_dir / "scaffold_audit.tsv.gz",
        "homology_0_30": evaluation_dir / "homology_0_30.tsv.gz",
        "homology_0_50": evaluation_dir / "homology_0_50.tsv.gz",
        "homology_0_70": evaluation_dir / "homology_0_70.tsv.gz",
    }
    write_tsv_gz(evaluation_paths["endpoint"], list(endpoint[0]), endpoint)
    write_tsv_gz(evaluation_paths["scaffold_audit"], list(scaffold[0]), scaffold)
    for threshold, path_key in (("0_30", "homology_0_30"), ("0_50", "homology_0_50"), ("0_70", "homology_0_70")):
        write_tsv_gz(evaluation_paths[path_key], list(homology[threshold][0]), homology[threshold])

    scoring_manifest = {
        "protocol_id": PROTOCOL_ID,
        "run_id": RUN_ID,
        "execution_mode": "author_run_non_independent",
        "project_lead_authorized_internal_use": True,
        "input_kind": "scoring_without_endpoint_file",
        "legacy_outer_or_result_input": False,
        "access_level": "restricted_author_run_internal",
        "authorization_basis": "Project-lead statement recorded in strict_ab_doublecold_successor_v1_20260719/governance/receipt_records/D0_project_lead_authorization_note_20260719.md; no public-release, human-gate, endpoint-access-control, or independent-evaluation claim.",
        "file_sha256": {name: sha256(path) for name, path in scoring_paths.items()},
    }
    evaluation_manifest = {
        "protocol_id": PROTOCOL_ID,
        "run_id": RUN_ID,
        "execution_mode": "author_run_non_independent",
        "project_lead_authorized_internal_use": True,
        "input_kind": "author_run_evaluation_endpoint",
        "legacy_outer_or_result_input": False,
        "access_level": "restricted_author_run_internal",
        "authorization_basis": "Project-lead statement recorded in strict_ab_doublecold_successor_v1_20260719/governance/receipt_records/D0_project_lead_authorization_note_20260719.md; no public-release, human-gate, endpoint-access-control, or independent-evaluation claim.",
        "file_sha256": {name: sha256(path) for name, path in evaluation_paths.items()},
    }
    scoring_manifest_path = scoring_dir / "author_run_input_manifest.json"
    evaluation_manifest_path = evaluation_dir / "author_run_input_manifest.json"
    write_json(scoring_manifest_path, scoring_manifest)
    write_json(evaluation_manifest_path, evaluation_manifest)

    author_run_receipt = {
        "protocol_id": PROTOCOL_ID,
        "run_id": RUN_ID,
        "execution_mode": "author_run_non_independent",
        "project_lead_authorized_internal_use": True,
        "authorization_basis": "Assistant-recorded project-lead statement; see the D0 authorization note in the successor governance receipt directory.",
        "technical_lock": {
            "protocol_spec_sha256": sha256(ROOT / "configs" / "successor_evaluation_spec_v1.json"),
            "author_run_runner_sha256": sha256(Path(__file__).parent / "run_author_run_strict_ab_successor.py"),
            "score_core_sha256": sha256(Path(__file__).parent / "score_successor_blind.py"),
            "evaluation_core_sha256": sha256(Path(__file__).parent / "evaluate_successor_sealed.py"),
        },
        "endpoint_handling": "The scoring command receives no endpoint or cold-scope file. No personnel or operating-system endpoint access separation is claimed; this is an author-run non-independent calculation.",
        "human_gate_status": "not_claimed",
        "release_status": "restricted_internal_no_public_release",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    author_run_receipt_path = metadata_dir / "author_run_protocol_receipt.json"
    write_json(author_run_receipt_path, author_run_receipt)
    provenance = {
        "protocol_id": PROTOCOL_ID,
        "run_id": RUN_ID,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "access_level": "restricted_author_run_internal",
        "source_files": {
            name: {"path": str(path), "sha256": sha256(path)} for name, path in SOURCE_PATHS.items()
        },
        "derived_file_hashes": {
            **{f"scoring/{name}": sha256(path) for name, path in scoring_paths.items()},
            **{f"evaluation/{name}": sha256(path) for name, path in evaluation_paths.items()},
            "scoring/author_run_input_manifest": sha256(scoring_manifest_path),
            "evaluation/author_run_input_manifest": sha256(evaluation_manifest_path),
            "metadata/author_run_protocol_receipt": sha256(author_run_receipt_path),
        },
        "counts": {
            "historical_pairs": len(history),
            "future_endpoint_pairs": len(endpoint),
            "scoring_queries": len(scoring_queries),
            "scoring_compounds": len(compounds),
            "candidate_targets": len(candidate_targets),
            "future_endpoint_targets": len(endpoint_targets),
            "scaffold_audit_pairs": len(scaffold),
            "homology_targets_each_threshold": {threshold: len(rows) for threshold, rows in homology.items()},
        },
        "normalizations": {
            "best_strict_evidence_tier": "When a source pair has both A and B records, retain A as the best strict tier; otherwise retain B.",
            "unrecorded_pair_policy": "Source wording normalized to unlabeled_not_negative without creating a negative label.",
            "c31_leakage_gate_status": "All 358 source rows passed include_primary_future_candidate_after_C31_leakage_screen; output uses the fixed normalized pass_no_historical_activity value required by the evaluator.",
        },
        "claim_boundary": "Author-run, non-independent successor evaluation only; no independent external validation, human-gate completion, endpoint access control, public release, direct-binding claim, or human P3 claim.",
    }
    provenance_path = metadata_dir / "bundle_provenance.json"
    write_json(provenance_path, provenance)
    print(
        json.dumps(
            {
                "output_root": str(output_root),
                "historical_pairs": len(history),
                "future_endpoint_pairs": len(endpoint),
                "scoring_queries": len(scoring_queries),
                "candidate_targets": len(candidate_targets),
                "figures_generated": False,
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
