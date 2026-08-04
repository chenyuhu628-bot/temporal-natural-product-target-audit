#!/usr/bin/env python3
"""Fail-closed validation for the aggregate-only v4 table set."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "tables"
REPORT = TABLES / "major_revision_table_validation_v4.json"

EXPECTED_ROWS = {
    "Table_1_temporal_repair_flow.tsv": 68,
    "Table_2_historical_before_after_audit.tsv": 145,
    "Table_3_corrected_aggregate_performance.tsv": 20,
    "Table_4_corrected_bootstrap_summaries.tsv": 18,
    "Table_5_score_degeneracy_and_ties.tsv": 4,
    "Table_6_claim_evidence_use_boundaries.tsv": 7,
    "Table_S1_scope_mask_integrity.tsv": 44,
    "Table_S2_top100_exhaustive_fidelity.tsv": 55,
    "Table_S3_zero_and_failure_accounting.tsv": 20,
    "Table_S4_pmid_document_dependence.tsv": 262,
    "Table_S5_frozen_unresolved_exclusions.tsv": 2,
    "Table_S6_reproducibility_and_release.tsv": 45,
    "Table_S7_tie_aware_retrieval.tsv": 1816,
    "Table_S8_date_precision_policy.tsv": 537,
    "Table_S10_unresolved_entity_bounds.tsv": 251,
    "Table_S11_maximum_similarity_distributions.tsv": 24,
    "Table_S12_rights_and_controlled_access.tsv": 10,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(name: str) -> list[dict[str, str]]:
    path = TABLES / name
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"Empty table: {name}")
    return rows


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def close(observed: str, expected: float, tolerance: float = 1e-12) -> bool:
    return math.isclose(float(observed), expected, rel_tol=0, abs_tol=tolerance)


def main() -> int:
    checks: list[dict[str, str]] = []

    def check(name: str, condition: bool, detail: str) -> None:
        require(condition, f"{name}: {detail}")
        checks.append({"check": name, "status": "PASS", "detail": detail})

    tables = {name: read_tsv(name) for name in EXPECTED_ROWS}
    check(
        "fixed_row_counts",
        all(len(tables[name]) == expected for name, expected in EXPECTED_ROWS.items()),
        str({name: len(rows) for name, rows in tables.items()}),
    )

    table4 = tables["Table_4_corrected_bootstrap_summaries.tsv"]
    check(
        "table4_scope_grid",
        sum(row["query_subset"] == "all_queries" for row in table4) == 16
        and sum(row["query_subset"] != "all_queries" for row in table4) == 2,
        "16 all-query rows plus two temporal 3-mer strata",
    )
    morgan = next(
        row
        for row in table4
        if row["scope"] == "temporal_strict_ab"
        and row["baseline"] == "weighted_morgan_transfer"
        and row["query_subset"] == "all_queries"
    )
    check(
        "temporal_morgan_tie_values",
        close(morgan["salted_recall_at_50"], 0.24422994422994423)
        and close(morgan["tie_expected_recall_at_50"], 0.24343164162004741)
        and close(morgan["tie_worst_recall_at_50"], 0.23972543972543975)
        and close(morgan["tie_best_recall_at_50"], 0.26029601029601029),
        "locked salt, expectation, and exact bounds",
    )
    joint = [
        row
        for row in table4
        if row["scope"].startswith("project_defined_joint_scaffold_homology")
    ]
    check(
        "joint_scope_zero_hit_bounds",
        len(joint) == 8
        and {row["empirical_zero_hit_status"] for row in joint}
        == {"empirical_degenerate_zero_hits"}
        and {
            round(float(row["one_sided_cp95_upper"]), 10)
            for row in joint
        }
        == {round(0.14586850331224344, 10), round(0.12730543165483876, 10)},
        "eight merged-scope rows with empirical-degenerate labels and exact upper bounds",
    )

    table8 = tables["Table_S8_date_precision_policy.tsv"]
    check(
        "date_policy_core",
        any(
            row["section"] == "history"
            and row["scenario_or_policy"] == "interval_certain_pre_cutoff"
            and row["item"] == "selected_source_row_count"
            and row["value"] == "20455"
            for row in table8
        )
        and any(
            row["section"] == "interval_status"
            and row["subgroup"] == "crossing_cutoff"
            and row["item"] == "source_row_count"
            and row["value"] == "0"
            for row in table8
        ),
        "20,455 interval-certain rows and zero crossing rows",
    )

    table9 = read_tsv("Table_S9_weight_and_structure_policy.tsv")
    check(
        "table9_weight_grid",
        any(
            row["section"] == "weight_metric"
            and row["scenario_or_policy"] == "A1_B1_0_all_equal"
            and row["scope"] == "temporal_strict_ab"
            and row["baseline"] == "weighted_morgan_transfer"
            and row["item"] == "Recall@50"
            and close(row["value"], 0.2758472758472758)
            for row in table9
        ),
        "all-equal temporal Morgan Recall@50 retained",
    )
    check(
        "table9_structure_outputs",
        any(row["section"].startswith("structure_") for row in table9),
        "validated structure-policy TSV records incorporated",
    )

    table10 = tables["Table_S10_unresolved_entity_bounds.tsv"]
    check(
        "unresolved_bounds",
        any(
            row["section"] == "relation_top50_bound"
            and row["baseline"] == "weighted_morgan_transfer"
            and row["item"] == "all_unresolved_fail_lower_fraction"
            and close(row["value"], 0.18439716312056736)
            for row in table10
        )
        and any(
            row["section"] == "relation_top50_bound"
            and row["baseline"] == "weighted_morgan_transfer"
            and row["item"] == "all_unresolved_succeed_upper_fraction"
            and close(row["value"], 0.3380614657210402)
            for row in table10
        ),
        "relation-level Morgan all-fail/all-success bounds",
    )

    table11 = tables["Table_S11_maximum_similarity_distributions.tsv"]
    corrected_mmseq = next(
        row
        for row in table11
        if row["similarity_family"] == "mmseqs2_detected_alignment_identity"
        and row["analysis_unit"] == "query_maximum"
        and row["scope"] == "joint_scaffold_homology_0_50_0_70_identical"
    )
    check(
        "corrected_mmseq_units",
        corrected_mmseq["n_observed"] == "5"
        and corrected_mmseq["n_no_detected_alignment"] == "17"
        and close(corrected_mmseq["min"], 0.252)
        and close(corrected_mmseq["max"], 0.465),
        "MMseqs2 pident is fractional, not clipped percentage",
    )

    manifest_path = TABLES / "major_revision_table_assembly_manifest_v4.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    check(
        "manifest_inventory",
        manifest.get("aggregate_only") is True
        and len(manifest.get("outputs", [])) == 18,
        "18 aggregate-only main/supplementary tables",
    )
    check(
        "manifest_hashes",
        all(
            (ROOT / item["path"]).is_file()
            and (ROOT / item["path"]).stat().st_size == item["bytes"]
            and sha256(ROOT / item["path"]) == item["sha256"]
            for item in manifest["outputs"]
        ),
        "all table hashes and byte counts match",
    )

    absolute_path = re.compile(r"(?i)(?:[A-Z]:[\\/]|\\\\[^\\\s]+\\[^\\\s]+)")
    inchikey = re.compile(r"\b[A-Z]{14}-[A-Z]{10}-[A-Z]\b")
    uniprot = re.compile(
        r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9])\b"
    )
    inspected = [ROOT / item["path"] for item in manifest["outputs"]]
    check(
        "aggregate_safety_scan",
        all(
            not absolute_path.search(path.read_text(encoding="utf-8"))
            and not inchikey.search(path.read_text(encoding="utf-8"))
            and not uniprot.search(path.read_text(encoding="utf-8"))
            for path in inspected
        ),
        "no absolute path or recognizable entity identifier",
    )

    payload = {
        "schema_version": "major_revision_table_validation_v4",
        "status": "PASS",
        "check_count": len(checks),
        "checks": checks,
        "validated_manifest_sha256": sha256(manifest_path),
    }
    REPORT.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"status": "PASS", "checks": len(checks)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
