"""Compute aggregate-only bounds for the frozen 65 entity-unresolved relations."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import platform
import statistics
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
OUTPUT = ROOT / "outputs"
ANALYSIS_ID = "revision_unresolved_bounds_v1_20260729"
PARENT_PROTOCOL_SHA256 = "bd4902476160cc7c5cbacaf0cfd0f1a28c5300bde22232b3a2cc6c1f3c143dc2"
BASELINES = (
    "weighted_target_popularity",
    "sequence_3mer_transfer",
    "weighted_morgan_transfer",
    "structure_sequence_pair_neighbor",
)
PRIMARY_RELATIONS = 358
HISTORICAL_ACTIVITY_EXCLUSIONS = 19
UNRESOLVED_RELATIONS = 65
INITIAL_RELATIONS = 442
MAX_ENDPOINT_RELATIONS = PRIMARY_RELATIONS + UNRESOLVED_RELATIONS
LOCKED = {
    "initial_decision_ledger": (
        "results/strict_temporal_future_v1_1_pmid_verified_chembl31_leakage_decision_ledger.csv.gz",
        "f22327c5def9cb44cc8ea5077e09f7d6bd50164cc3ef574fb61be54637259dcb",
    ),
    "unresolved_ledger": (
        "data/processed/strict_temporal_future_v1_1_pmid_verified_chembl31_C31_entity_unresolved.csv.gz",
        "1d9ef5960a63b8955cd9a24513c05cf409e23119628707c597d8440f4ad56b9c",
    ),
    "preliminary_mapping": (
        "data/interim/chembl_31_future_candidate_entity_mapping.csv.gz",
        "e5175f8bfd18a2e8abc0dffd72eb32bb25fec3a05ab6c1aa79c63be5f5ddc46e",
    ),
    "sqlite_validation": (
        "data/interim/chembl_31_future_candidate_sqlite_entity_validation.csv.gz",
        "3e99ad8faea406f527a7a1a2f942d4d2a63faca68f16febedc575947b092178e",
    ),
    "endpoint": (
        "author_run_strict_ab_asof_cutoff_execution_v1_20260728/evaluation_inputs/evaluation_pairs.tsv.gz",
        "09296b066a23197a7c178f00514f2b3d9ed7e6f3c459ea92a55e01a6010d1132",
    ),
    "ranks": (
        "author_run_strict_ab_asof_cutoff_execution_v1_20260728/score/corrective_prediction_ranks.tsv.gz",
        "87739aa818744c7084088d13c386444aa41bbef38c257083325298003181479e",
    ),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def open_rows(path: Path, delimiter: str) -> list[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        require(bool(reader.fieldnames), f"Missing header: {path.name}")
        return list(reader)


def parse_bool(value: str) -> bool:
    normalized = value.strip().casefold()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no"}:
        return False
    raise ValueError(f"Invalid boolean {value!r}")


def split_count(value: str) -> int:
    return len({part.strip() for part in value.split(";") if part.strip()})


def distribution(values: Iterable[int]) -> dict[str, float | int]:
    array = np.asarray(list(values), dtype=float)
    require(array.size > 0, "Cannot summarize an empty distribution")
    return {
        "min": int(array.min()),
        "q25": float(np.quantile(array, 0.25)),
        "median": float(np.median(array)),
        "mean": float(array.mean()),
        "q75": float(np.quantile(array, 0.75)),
        "max": int(array.max()),
    }


def verify_sources() -> tuple[dict[str, Path], dict[str, str]]:
    paths: dict[str, Path] = {}
    hashes: dict[str, str] = {}
    for role, (relative, expected) in LOCKED.items():
        path = WORKSPACE / relative
        require(path.is_file(), f"Missing locked input: {relative}")
        actual = sha256(path)
        require(actual == expected, f"Locked input hash mismatch: {role}")
        paths[role] = path
        hashes[relative] = actual
    return paths, dict(sorted(hashes.items()))


def cohort_summary(label: str, rows: list[dict[str, str]]) -> dict[str, object]:
    numeric_fields = (
        "future_v3_all_record_count",
        "future_v3_primary_A_B_P1_record_count",
        "future_v3_strict_entity_primary_record_count",
        "audit_npass_record_count",
    )
    result: dict[str, object] = {
        "cohort": label,
        "relation_count": len(rows),
        "distinct_compound_count": len({row["inchikey_full"] for row in rows}),
        "distinct_target_count": len({row["uniprot_canonical_accession"] for row in rows}),
        "tier_A_present_count": sum("A_affinity_candidate" in row["future_primary_evidence_tiers"] for row in rows),
        "tier_A_present_fraction": f"{sum('A_affinity_candidate' in row['future_primary_evidence_tiers'] for row in rows) / len(rows):.17g}",
        "tier_B_present_count": sum("B_quantitative_functional_candidate" in row["future_primary_evidence_tiers"] for row in rows),
        "tier_B_present_fraction": f"{sum('B_quantitative_functional_candidate' in row['future_primary_evidence_tiers'] for row in rows) / len(rows):.17g}",
        "preliminary_both_entities_exact_count": sum(parse_bool(row["audit_preliminary_both_entities_exactly_mapped"]) for row in rows),
        "preliminary_both_entities_exact_fraction": f"{sum(parse_bool(row['audit_preliminary_both_entities_exactly_mapped']) for row in rows) / len(rows):.17g}",
        "sqlite_validated_mapping_positive_count": sum(int(row["audit_sqlite_validated_entity_mapping_count"]) > 0 for row in rows),
        "sqlite_validated_mapping_positive_fraction": f"{sum(int(row['audit_sqlite_validated_entity_mapping_count']) > 0 for row in rows) / len(rows):.17g}",
    }
    for field in numeric_fields:
        summary = distribution(int(row[field]) for row in rows)
        for statistic, value in summary.items():
            result[f"{field}_{statistic}"] = value
    for derived, values in (
        ("primary_reference_count", [split_count(row["future_primary_references"]) for row in rows]),
        ("primary_activity_type_count", [split_count(row["future_primary_activity_types"]) for row in rows]),
    ):
        summary = distribution(values)
        for statistic, value in summary.items():
            result[f"{derived}_{statistic}"] = value
    return result


def failure_reason(
    pair: str,
    preliminary: dict[str, dict[str, str]],
    sqlite_rows: dict[str, list[dict[str, str]]],
) -> str:
    mapping = preliminary[pair]
    compound_exact = mapping["chembl_compound_match_status"] == "full_inchikey_exact"
    target_exact = mapping["chembl_target_match_status"] == "source_uniprot_exact"
    declared = parse_bool(mapping["both_entities_exactly_mapped"])
    require(declared == (compound_exact and target_exact), "Preliminary joint mapping flag is inconsistent")
    if not compound_exact and not target_exact:
        return "preliminary_compound_and_target_unmatched"
    if not compound_exact:
        return "preliminary_compound_unmatched"
    if not target_exact:
        return "preliminary_target_unmatched"
    candidates = sqlite_rows.get(pair, [])
    require(not any(parse_bool(row["sqlite_entity_pair_validated"]) for row in candidates), "Unresolved relation has a validated SQLite pair")
    return "sqlite_no_validated_joint_mapping"


def extract_top50_hits(
    rank_path: Path,
    endpoint: list[dict[str, str]],
) -> tuple[dict[str, int], int]:
    relevant: dict[str, set[str]] = defaultdict(set)
    for row in endpoint:
        relevant[row["query_id"]].add(row["uniprot_canonical_accession"])
    hit_counts = Counter({baseline: 0 for baseline in BASELINES})
    found = 0
    with gzip.open(rank_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["target_uniprot_accession"] not in relevant.get(row["query_id"], set()):
                continue
            baseline = row["baseline"]
            require(baseline in BASELINES, "Unknown baseline in rank ledger")
            found += 1
            hit_counts[baseline] += int(row["rank"]) <= 50
    require(found == len(BASELINES) * len(endpoint), "Endpoint rank extraction is incomplete")
    return dict(hit_counts), found


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    require(bool(rows), f"Empty output: {path.name}")
    require(not path.exists(), f"Refusing to overwrite: {path.name}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    require(not path.exists(), f"Refusing to overwrite: {path.name}")
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def summary_markdown(
    failures: list[dict[str, object]],
    comparability: list[dict[str, object]],
    bounds: list[dict[str, object]],
) -> str:
    lines = [
        "# Entity-unresolved endpoint bounds",
        "",
        "Author-run, outcome-visible, post hoc descriptive sensitivity; no mappings were invented.",
        "",
        "## Mapping-failure strata",
        "",
        "| Stratum | Relations | Fraction |",
        "|---|---:|---:|",
    ]
    for row in failures:
        lines.append(f"| {row['failure_stratum']} | {row['relation_count']} | {float(row['relation_fraction']):.4f} |")
    lines += [
        "",
        "## Observable source-record comparison",
        "",
        "| Cohort | Relations | Compounds | Targets | Mean v3 records | A present | B present |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in comparability:
        lines.append(
            f"| {row['cohort']} | {row['relation_count']} | {row['distinct_compound_count']} | "
            f"{row['distinct_target_count']} | {float(row['future_v3_all_record_count_mean']):.3f} | "
            f"{float(row['tier_A_present_fraction']):.3f} | {float(row['tier_B_present_fraction']):.3f} |"
        )
    lines += [
        "",
        "These observable summaries describe the recorded candidate rows; they do not establish missing-at-random mapping failure.",
        "",
        "## Endpoint and relation-level top-50 bounds",
        "",
        "The identified endpoint cardinality is 358–423 relations. The 19 historical-activity exclusions remain excluded.",
        "",
        "| Baseline | Observed hits / 358 | All 65 fail | All 65 succeed |",
        "|---|---:|---:|---:|",
    ]
    for row in bounds:
        lines.append(
            f"| {row['baseline']} | {row['observed_top50_hit_relation_count']} / 358 "
            f"({float(row['observed_primary_relation_hit_fraction']):.6f}) | "
            f"{float(row['all_unresolved_fail_lower_fraction']):.6f} | "
            f"{float(row['all_unresolved_succeed_upper_fraction']):.6f} |"
        )
    lines += [
        "",
        "The bounds are relation-weighted temporal top-50 bounds only. Query-macro Recall/NDCG/MRR and scaffold/homology-scope bounds are not identifiable without the missing mappings and ranks.",
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    started = time.perf_counter()
    require(OUTPUT.is_dir() and not any(OUTPUT.iterdir()), "Create-once output directory must be empty")
    paths, source_hashes = verify_sources()
    initial = open_rows(paths["initial_decision_ledger"], ",")
    unresolved = open_rows(paths["unresolved_ledger"], ",")
    preliminary_rows = open_rows(paths["preliminary_mapping"], ",")
    sqlite_validation = open_rows(paths["sqlite_validation"], ",")
    endpoint = open_rows(paths["endpoint"], "\t")

    require(len(initial) == INITIAL_RELATIONS, "Initial decision ledger count mismatch")
    require(len(unresolved) == UNRESOLVED_RELATIONS, "Unresolved ledger count mismatch")
    require(len(endpoint) == PRIMARY_RELATIONS, "Primary endpoint count mismatch")
    require(INITIAL_RELATIONS == PRIMARY_RELATIONS + HISTORICAL_ACTIVITY_EXCLUSIONS + UNRESOLVED_RELATIONS, "Endpoint flow arithmetic failed")

    by_stratum: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in initial:
        by_stratum[row["leakage_gate_stratum"]].append(row)
    retained = by_stratum["primary_C31_validated_no_historical_activity"]
    historical_overlap = by_stratum["historical_C31_activity_hit"]
    unresolved_initial = by_stratum["C31_entity_unresolved"]
    require(len(retained) == PRIMARY_RELATIONS, "Retained stratum count mismatch")
    require(len(historical_overlap) == HISTORICAL_ACTIVITY_EXCLUSIONS, "Historical-overlap count mismatch")
    require(len(unresolved_initial) == UNRESOLVED_RELATIONS, "Unresolved stratum count mismatch")
    unresolved_keys = {row["canonical_pair_key"] for row in unresolved}
    require(unresolved_keys == {row["canonical_pair_key"] for row in unresolved_initial}, "Frozen unresolved keyset differs from initial ledger")
    require(all(not parse_bool(row["negative_label_emitted"]) for row in unresolved), "Unresolved relation emitted a negative label")
    require(all(int(row["audit_sqlite_validated_entity_mapping_count"]) == 0 for row in unresolved), "Unresolved relation has a validated mapping")

    preliminary = {row["pair_key"]: row for row in preliminary_rows}
    require(len(preliminary) == len(preliminary_rows), "Preliminary mapping ledger has duplicate keys")
    sqlite_rows: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in sqlite_validation:
        sqlite_rows[row["pair_key"]].append(row)
    reason_counts = Counter(
        failure_reason(row["canonical_pair_key"], preliminary, sqlite_rows) for row in unresolved
    )
    require(sum(reason_counts.values()) == UNRESOLVED_RELATIONS, "Mapping failure reasons are incomplete")
    failure_rows = [
        {
            "failure_stratum": reason,
            "relation_count": count,
            "relation_fraction": f"{count / UNRESOLVED_RELATIONS:.17g}",
        }
        for reason, count in sorted(reason_counts.items())
    ]

    comparability_rows = [
        cohort_summary("entity_unresolved", unresolved_initial),
        cohort_summary("primary_resolved_no_historical_activity", retained),
        cohort_summary("entity_resolved_all", retained + historical_overlap),
    ]
    indexed = {row["cohort"]: row for row in comparability_rows}
    unresolved_summary = indexed["entity_unresolved"]
    retained_summary = indexed["primary_resolved_no_historical_activity"]
    difference_rows: list[dict[str, object]] = []
    for measure in (
        "future_v3_all_record_count_mean",
        "future_v3_primary_A_B_P1_record_count_mean",
        "future_v3_strict_entity_primary_record_count_mean",
        "audit_npass_record_count_mean",
        "primary_reference_count_mean",
        "primary_activity_type_count_mean",
        "tier_A_present_fraction",
        "tier_B_present_fraction",
    ):
        unresolved_value = float(unresolved_summary[measure])
        retained_value = float(retained_summary[measure])
        difference_rows.append(
            {
                "measure": measure,
                "entity_unresolved_value": f"{unresolved_value:.17g}",
                "primary_resolved_value": f"{retained_value:.17g}",
                "unresolved_minus_resolved": f"{unresolved_value - retained_value:.17g}",
                "interpretation": "observable_descriptive_difference_only_not_a_missing_at_random_test",
            }
        )

    hit_counts, extracted_rank_cells = extract_top50_hits(paths["ranks"], endpoint)
    bound_rows: list[dict[str, object]] = []
    for baseline in BASELINES:
        observed_hits = hit_counts[baseline]
        bound_rows.append(
            {
                "baseline": baseline,
                "observed_endpoint_relation_count": PRIMARY_RELATIONS,
                "unresolved_relation_count": UNRESOLVED_RELATIONS,
                "maximum_potential_endpoint_relation_count": MAX_ENDPOINT_RELATIONS,
                "observed_top50_hit_relation_count": observed_hits,
                "observed_primary_relation_hit_fraction": f"{observed_hits / PRIMARY_RELATIONS:.17g}",
                "all_unresolved_fail_top50_hit_count": observed_hits,
                "all_unresolved_fail_lower_fraction": f"{observed_hits / MAX_ENDPOINT_RELATIONS:.17g}",
                "all_unresolved_succeed_top50_hit_count": observed_hits + UNRESOLVED_RELATIONS,
                "all_unresolved_succeed_upper_fraction": f"{(observed_hits + UNRESOLVED_RELATIONS) / MAX_ENDPOINT_RELATIONS:.17g}",
                "estimand": "relation_weighted_temporal_top50_hit_fraction",
                "status": "sharp_assumption_bounds_no_mapping_or_rank_imputation",
            }
        )

    cardinality_rows = [
        {
            "initial_candidate_relations": INITIAL_RELATIONS,
            "definitive_historical_activity_exclusions": HISTORICAL_ACTIVITY_EXCLUSIONS,
            "frozen_primary_endpoint_relations": PRIMARY_RELATIONS,
            "entity_unresolved_relations": UNRESOLVED_RELATIONS,
            "identified_endpoint_cardinality_lower": PRIMARY_RELATIONS,
            "identified_endpoint_cardinality_upper": MAX_ENDPOINT_RELATIONS,
            "unresolved_negative_labels": 0,
            "unresolved_readmissions": 0,
        }
    ]
    nonidentified_rows = [
        {
            "estimand_family": "query_macro_Recall_NDCG_MRR",
            "status": "not_identifiable",
            "exact_blocker": "Missing entity resolution can change query membership, per-query relevant counts, target identity, candidate masks, and ranks; relation-level success assumptions do not bound equal-query-weighted metrics.",
            "invented_values_used": "false",
        },
        {
            "estimand_family": "scaffold_cold_scope_metrics",
            "status": "not_identifiable",
            "exact_blocker": "Unresolved compound identity or validated representation prevents assigning a defensible scaffold-scope membership and the corresponding target rank.",
            "invented_values_used": "false",
        },
        {
            "estimand_family": "joint_scaffold_homology_scope_metrics",
            "status": "not_identifiable",
            "exact_blocker": "Unresolved compound or target mapping prevents assigning both scaffold and homology membership and prevents reading a frozen target rank.",
            "invented_values_used": "false",
        },
    ]

    tables = {
        "mapping_failure_strata.tsv": failure_rows,
        "cohort_comparability.tsv": comparability_rows,
        "comparability_differences.tsv": difference_rows,
        "endpoint_cardinality_bounds.tsv": cardinality_rows,
        "relation_level_top50_bounds.tsv": bound_rows,
        "non_identifiable_estimands.tsv": nonidentified_rows,
    }
    for filename, rows in tables.items():
        write_tsv(OUTPUT / filename, rows)
    summary_path = OUTPUT / "SUMMARY.md"
    summary_path.write_text(summary_markdown(failure_rows, comparability_rows, bound_rows), encoding="utf-8")

    output_hashes = {path.name: sha256(path) for path in sorted(OUTPUT.iterdir()) if path.is_file()}
    receipt = {
        "analysis_id": ANALYSIS_ID,
        "status": "PASS",
        "execution_mode": "author_run_outcome_visible_post_hoc_descriptive_missing_endpoint_bounds",
        "claim_boundary": "No remapping, rank imputation, query-macro bound, scope-specific bound, negative label, readmission, external validation, or biological confirmation claim.",
        "initial_candidate_partition": {
            "retained_primary": len(retained),
            "historical_activity_excluded": len(historical_overlap),
            "entity_unresolved": len(unresolved_initial),
        },
        "endpoint_cardinality_bounds": [PRIMARY_RELATIONS, MAX_ENDPOINT_RELATIONS],
        "endpoint_rank_cells_extracted": extracted_rank_cells,
        "relation_level_bounds_computed": len(bound_rows),
        "query_macro_bounds_computed": 0,
        "scope_specific_bounds_computed": 0,
        "identifier_bearing_output": False,
        "absolute_paths_emitted": False,
        "main_manuscript_modified": False,
        "runtime_seconds": time.perf_counter() - started,
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "platform": platform.platform(),
        },
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_json(OUTPUT / "EXECUTION_RECEIPT.json", receipt)
    output_hashes["EXECUTION_RECEIPT.json"] = sha256(OUTPUT / "EXECUTION_RECEIPT.json")

    code_files = [
        ROOT / "PROTOCOL.md",
        ROOT / "IMPLEMENTATION_LOCK.json",
        ROOT / "CODE_LOCK.json",
        Path(__file__),
        ROOT / "scripts" / "test_unresolved_bounds.py",
        ROOT / "scripts" / "validate_unresolved_bounds.py",
    ]
    manifest = {
        "schema_version": "1.0",
        "analysis_id": ANALYSIS_ID,
        "parent_protocol_sha256": PARENT_PROTOCOL_SHA256,
        "source_files": source_hashes,
        "code_files": {
            str(path.relative_to(WORKSPACE)).replace("\\", "/"): sha256(path) for path in code_files
        },
        "outputs_before_manifest": dict(sorted(output_hashes.items())),
        "row_counts": {name: len(rows) for name, rows in tables.items()},
        "aggregate_only": True,
        "identified_endpoint_cardinality_bounds": [PRIMARY_RELATIONS, MAX_ENDPOINT_RELATIONS],
        "identified_retrieval_bound_family": "relation_weighted_temporal_top50",
        "nonidentified_bound_families": [
            "query_macro_metrics",
            "scaffold_cold_metrics",
            "joint_scaffold_homology_metrics",
        ],
    }
    write_json(OUTPUT / "MANIFEST.json", manifest)
    print(json.dumps({"status": "PASS", "bounds": len(bound_rows), "runtime_seconds": receipt["runtime_seconds"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
