"""Validate aggregate-only contracts for the tie-aware revision outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


EXPECTED_FILES = {
    "tie_aware_metrics.tsv",
    "three_mer_operability.tsv",
    "double_cold_query_hit_upper_bounds.tsv",
    "uncertainty_estimands.tsv",
    "tie_aware_summary.json",
    "execution_receipt.json",
    "run_manifest.json",
}
DISPLAY_SCOPES = {
    "temporal_strict_ab",
    "scaffold_cold_strict_ab",
    "project_defined_joint_scaffold_homology_cold_0_30",
    "project_defined_joint_scaffold_homology_cold_0_50_0_70_identical_mask",
}
BASELINES = {
    "weighted_target_popularity",
    "sequence_3mer_transfer",
    "weighted_morgan_transfer",
    "structure_sequence_pair_neighbor",
}
FORBIDDEN_IDENTIFIER_HEADERS = {
    "query_id",
    "canonical_pair_key",
    "query_compound_inchikey_full",
    "inchikey_full",
    "target_uniprot_accession",
    "uniprot_canonical_accession",
    "ref_id",
    "pmid",
    "component_membership",
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or [])
        return fields, list(reader)


def float_value(row: dict[str, str], field: str) -> float:
    value = float(row[field])
    require(math.isfinite(value), f"Nonfinite value in {field}")
    return value


def validate_no_identifier_leakage(output_dir: Path) -> None:
    inchikey = re.compile(r"\b[A-Z]{14}-[A-Z]{10}-[A-Z]\b")
    query_id = re.compile(r"\bquery_[0-9]{4,}\b", flags=re.IGNORECASE)
    absolute_path = re.compile(r"(?:\b[A-Za-z]:[\\/]|(?:^|[\s\"'])/[A-Za-z0-9_.-]+/)")
    for path in sorted(output_dir.iterdir()):
        if path.suffix not in {".tsv", ".json"}:
            continue
        text = path.read_text(encoding="utf-8")
        require(not inchikey.search(text), f"InChIKey-like identifier leaked in {path.name}")
        require(not query_id.search(text), f"Query identifier leaked in {path.name}")
        require(not absolute_path.search(text), f"Absolute path leaked in {path.name}")
        if path.suffix == ".tsv":
            fields, _ = read_tsv(path)
            require(
                not FORBIDDEN_IDENTIFIER_HEADERS.intersection(fields),
                f"Identifier-bearing column leaked in {path.name}",
            )


def validate_metric_rows(rows: list[dict[str, str]]) -> None:
    require(len(rows) == 48, "Tie-aware metric row count is not 48")
    require({row["scope"] for row in rows} == DISPLAY_SCOPES, "Display scopes changed")
    require({row["baseline"] for row in rows} == BASELINES, "Baselines changed")
    require({int(row["k"]) for row in rows} == {10, 50}, "Cutoffs changed")
    seen: set[tuple[str, str, str, int]] = set()
    for row in rows:
        key = (row["scope"], row["baseline"], row["query_subset"], int(row["k"]))
        require(key not in seen, "Tie-aware metric cell is duplicated")
        seen.add(key)
        query_count = int(row["query_count"])
        relation_count = int(row["relevant_relation_count"])
        require(query_count >= 0 and relation_count >= query_count, "Invalid denominators")
        if query_count == 0:
            require(
                row["baseline"] == "sequence_3mer_transfer"
                and row["query_subset"] in {
                    "score_operational",
                    "structural_all_zero_non_operational",
                },
                "An unplanned cell is empty",
            )
            require(relation_count == 0, "Empty query subset has relevant relations")
            for field in (
                "legacy_salted_recall",
                "tie_expected_fractional_recall",
                "tie_worst_recall",
                "tie_best_recall",
                "legacy_salted_ndcg",
                "tie_expected_ndcg",
                "tie_worst_ndcg",
                "tie_best_ndcg",
                "legacy_salted_query_any_hit_rate",
                "tie_expected_query_any_hit_probability",
            ):
                require(row[field] == "", f"Empty query subset has a value in {field}")
            require(
                row["query_bootstrap_status"] == "not_estimable_no_queries",
                "Empty query subset has the wrong query-bootstrap status",
            )
            require(
                row["pmid_component_bootstrap_status"]
                == "not_estimable_component_count_lt_2",
                "Empty query subset has the wrong component-bootstrap status",
            )
            require(
                row["tie_interpretation"] == "not_estimable_empty_prespecified_subset",
                "Empty query subset lacks an explicit interpretation",
            )
        else:
            for metric in ("recall", "ndcg"):
                expected_field = (
                    "tie_expected_fractional_recall"
                    if metric == "recall"
                    else "tie_expected_ndcg"
                )
                lower = float_value(row, f"tie_worst_{metric}")
                expected = float_value(row, expected_field)
                upper = float_value(row, f"tie_best_{metric}")
                salted = float_value(row, f"legacy_salted_{metric}")
                require(
                    -1e-15 <= lower <= expected + 1e-15,
                    "Expected metric below lower bound",
                )
                require(
                    expected <= upper + 1e-15 <= 1.0 + 2e-15,
                    "Expected metric above upper bound",
                )
                require(0.0 <= salted <= 1.0, "Salted metric outside [0,1]")
        membership_sum = sum(
            int(row[field])
            for field in (
                "membership_score_identifiable_relation_count",
                "membership_boundary_tie_dependent_relation_count",
                "membership_not_retrieved_relation_count",
            )
        )
        require(membership_sum == relation_count, "Membership classes do not sum")
        require(
            int(row["recall_score_identifiable_query_count"])
            + int(row["recall_tie_dependent_query_count"])
            == query_count,
            "Recall query classes do not sum",
        )
        require(
            int(row["ndcg_score_identifiable_query_count"])
            + int(row["ndcg_tie_dependent_query_count"])
            == query_count,
            "NDCG query classes do not sum",
        )
        require(
            row["query_bootstrap_status"]
            in {
                "estimable_descriptive_query_bootstrap",
                "n=1_descriptive_point_only",
                "not_estimable_no_queries",
            },
            "Unknown query-bootstrap status",
        )
        require(
            row["pmid_component_bootstrap_status"]
            in {
                "estimable_descriptive_pmid_component_sensitivity",
                "not_estimable_component_count_lt_2",
            },
            "Unknown component-bootstrap status",
        )
        if row["query_subset"] == "structural_all_zero_non_operational":
            require(
                row["tie_interpretation"] == "non_operational_uniform_tie_allocation",
                "Structural all-zero rows are not marked non-operational",
            )
    expected_cell_count = 4 * (4 + 2) * 2
    require(len(seen) == expected_cell_count, "Metric matrix is incomplete")


def validate_operability(rows: list[dict[str, str]]) -> None:
    require(len(rows) == 4, "3-mer operability row count is not four")
    temporal = next(row for row in rows if row["scope"] == "temporal_strict_ab")
    require(int(temporal["all_query_count"]) == 222, "Temporal query count changed")
    require(
        int(temporal["score_operational_query_count"]) == 60,
        "Operational 3-mer query count changed",
    )
    require(
        int(temporal["structural_all_zero_query_count"]) == 162,
        "Structural all-zero 3-mer query count changed",
    )
    for row in rows:
        require(
            int(row["score_operational_query_count"])
            + int(row["structural_all_zero_query_count"])
            == int(row["all_query_count"]),
            "3-mer operability query counts do not sum",
        )


def validate_zero_hit(rows: list[dict[str, str]]) -> None:
    require(len(rows) == 8, "Zero-hit row count is not eight")
    require({row["baseline"] for row in rows} == BASELINES, "Zero-hit baselines changed")
    for row in rows:
        n_queries = int(row["query_count"])
        require(
            int(row["legacy_salted_query_hit_count_at_50"]) == 0,
            "A zero-hit cell contains a hit",
        )
        require(
            row["empirical_query_bootstrap_status"] == "empirical_degenerate_zero_hits",
            "Zero-width empirical interval is not labelled degenerate",
        )
        observed = float_value(row, "one_sided_clopper_pearson_upper")
        expected = 1.0 - 0.05 ** (1.0 / n_queries)
        require(
            math.isclose(observed, expected, rel_tol=1e-14, abs_tol=1e-14),
            "Clopper-Pearson upper bound is incorrect",
        )


def validate_manifest(output_dir: Path, payload: dict[str, Any]) -> None:
    require(payload.get("aggregate_only") is True, "Manifest is not aggregate-only")
    require(
        payload.get("identifier_bearing_output") is False,
        "Manifest asserts identifier-bearing output",
    )
    require(
        payload.get("absolute_paths_recorded") is False,
        "Manifest asserts absolute paths",
    )
    outputs = payload.get("outputs")
    require(isinstance(outputs, list) and len(outputs) == 6, "Manifest output list changed")
    for item in outputs:
        path = output_dir / item["basename"]
        require(path.is_file(), f"Manifest output is absent: {item['basename']}")
        require(path.stat().st_size == item["bytes"], "Manifest byte count mismatch")
        require(sha256(path) == item["sha256"], "Manifest output hash mismatch")
    contract = payload.get("output_contract", {})
    require(contract.get("tie_aware_metrics_rows") == 48, "Manifest metric contract changed")
    require(
        contract.get("homology_0_50_0_70_displayed_once") is True,
        "Identical-mask display contract changed",
    )


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--output-dir", required=True, type=Path)
    result.add_argument("--receipt", required=True, type=Path)
    return result


def main() -> int:
    args = parser().parse_args()
    output_dir = args.output_dir.resolve()
    receipt_path = args.receipt.resolve()
    require(output_dir.is_dir(), "Output directory is absent")
    require(not receipt_path.exists(), "Validation receipt already exists")
    require(
        {path.name for path in output_dir.iterdir()} == EXPECTED_FILES,
        "Output file inventory changed before validation",
    )
    metric_fields, metric_rows = read_tsv(output_dir / "tie_aware_metrics.tsv")
    _, operability_rows = read_tsv(output_dir / "three_mer_operability.tsv")
    _, zero_hit_rows = read_tsv(output_dir / "double_cold_query_hit_upper_bounds.tsv")
    _, estimand_rows = read_tsv(output_dir / "uncertainty_estimands.tsv")
    require(metric_fields, "Metric table has no header")
    validate_metric_rows(metric_rows)
    validate_operability(operability_rows)
    validate_zero_hit(zero_hit_rows)
    require(len(estimand_rows) == 4, "Uncertainty-estimand row count is not four")
    manifest_path = output_dir / "run_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    validate_manifest(output_dir, manifest)
    validate_no_identifier_leakage(output_dir)
    summary = json.loads((output_dir / "tie_aware_summary.json").read_text(encoding="utf-8"))
    require(
        summary.get("output_boundary", {}).get("identifier_bearing_rows") is False,
        "Summary output boundary changed",
    )
    require(
        summary.get("three_mer_temporal_operability", {}).get(
            "score_operational_query_count"
        )
        == 60,
        "Summary 3-mer operational count changed",
    )
    checks = [
        "exact output inventory",
        "48-cell merged-scope tie-aware metric matrix",
        "metric expectations lie within exact best-worst bounds",
        "relation membership classifications are exhaustive",
        "query identifiability classifications are exhaustive",
        "60 operational and 162 structural-all-zero temporal 3-mer queries",
        "eight zero-hit cells labelled empirical degenerate",
        "one-sided exact Clopper-Pearson bounds",
        "query and PMID-component estimands separately documented",
        "manifest output hashes and byte counts",
        "no identifier-bearing columns or recognizable query/InChIKey values",
        "no absolute paths",
    ]
    receipt = {
        "schema_version": "1.0",
        "analysis_id": "revision_tie_aware_v1_20260729",
        "status": "PASS",
        "check_count": len(checks),
        "checks": checks,
        "validated_manifest": {
            "basename": manifest_path.name,
            "sha256": sha256(manifest_path),
        },
        "validated_output_file_count": len(EXPECTED_FILES),
        "identifier_bearing_output_detected": False,
        "absolute_path_detected": False,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    receipt_path.parent.mkdir(parents=True, exist_ok=True)
    receipt_path.write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"status": "PASS", "check_count": len(checks)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
