"""Validate aggregate-only Tier B weight-policy sensitivity outputs."""

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
    "weight_variants.tsv",
    "aggregate_metrics.tsv",
    "metric_deltas_vs_0_7.tsv",
    "complete_rank_top50_changes_vs_0_7.tsv",
    "scope_cardinality_invariance.tsv",
    "weight_policy_summary.json",
    "execution_receipt.json",
    "run_manifest.json",
}
VARIANTS = {"A1_B0_5", "A1_B0_7_primary", "A1_B1_0_all_equal"}
BASELINES = {
    "weighted_target_popularity",
    "sequence_3mer_transfer",
    "weighted_morgan_transfer",
    "structure_sequence_pair_neighbor",
}
SCOPES = {
    "temporal_strict_ab",
    "scaffold_cold_strict_ab",
    "project_defined_joint_scaffold_homology_cold_0_30",
    "project_defined_joint_scaffold_homology_cold_0_50_0_70_identical_mask",
}
METRICS = {"Recall@10", "Recall@50", "NDCG@10", "NDCG@50", "MRR"}
FORBIDDEN_HEADERS = {
    "query_id",
    "canonical_pair_key",
    "inchikey_full",
    "query_compound_inchikey_full",
    "uniprot_canonical_accession",
    "target_uniprot_accession",
    "rank",
    "score",
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
        return list(reader.fieldnames or []), list(reader)


def finite(value: str, label: str) -> float:
    parsed = float(value)
    require(math.isfinite(parsed), f"Nonfinite value in {label}")
    return parsed


def validate_release_boundary(output_dir: Path) -> None:
    inchikey = re.compile(r"\b[A-Z]{14}-[A-Z]{10}-[A-Z]\b")
    query_id = re.compile(r"\bquery_[0-9]{4,}\b", flags=re.IGNORECASE)
    absolute_path = re.compile(r"(?:\b[A-Za-z]:[\\/]|(?:^|[\s\"'])/[A-Za-z0-9_.-]+/)")
    for path in sorted(output_dir.iterdir()):
        text = path.read_text(encoding="utf-8")
        require(not inchikey.search(text), f"InChIKey-like identifier in {path.name}")
        require(not query_id.search(text), f"Query identifier in {path.name}")
        require(not absolute_path.search(text), f"Absolute path in {path.name}")
        if path.suffix == ".tsv":
            fields, _ = read_tsv(path)
            require(
                not FORBIDDEN_HEADERS.intersection(fields),
                f"Identifier-bearing rank/score field in {path.name}",
            )


def validate_weight_rows(rows: list[dict[str, str]]) -> None:
    require(len(rows) == 3, "Weight-variant table is not three rows")
    require({row["weight_variant"] for row in rows} == VARIANTS, "Variants changed")
    indexed = {row["weight_variant"]: row for row in rows}
    require(
        finite(indexed["A1_B0_5"]["tier_B_weight"], "B weight") == 0.5,
        "0.5 variant changed",
    )
    require(
        finite(indexed["A1_B0_7_primary"]["tier_B_weight"], "B weight") == 0.7,
        "0.7 variant changed",
    )
    all_equal = indexed["A1_B1_0_all_equal"]
    require(
        finite(all_equal["tier_A_weight"], "A weight") == 1.0
        and finite(all_equal["tier_B_weight"], "B weight") == 1.0
        and all_equal["all_equal_alias"] == "all_equal_A1_B1",
        "All-equal alias is incorrect",
    )


def validate_aggregate(rows: list[dict[str, str]]) -> dict[tuple[str, str, str], dict[str, str]]:
    require(len(rows) == 48, "Aggregate metric table is not 48 rows")
    require({row["weight_variant"] for row in rows} == VARIANTS, "Aggregate variants changed")
    require({row["baseline"] for row in rows} == BASELINES, "Aggregate baselines changed")
    require({row["scope"] for row in rows} == SCOPES, "Aggregate scopes changed")
    index: dict[tuple[str, str, str], dict[str, str]] = {}
    for row in rows:
        key = (row["weight_variant"], row["scope"], row["baseline"])
        require(key not in index, "Aggregate cell duplicated")
        index[key] = row
        require(int(row["query_count"]) > 0, "Aggregate query denominator is empty")
        for metric in METRICS:
            value = finite(row[metric], metric)
            require(0.0 <= value <= 1.0, f"{metric} lies outside [0,1]")
    return index


def validate_deltas(
    rows: list[dict[str, str]],
    aggregate: dict[tuple[str, str, str], dict[str, str]],
) -> None:
    require(len(rows) == 240, "Metric-delta table is not 240 rows")
    seen: set[tuple[str, str, str, str]] = set()
    for row in rows:
        key = (
            row["weight_variant"],
            row["scope"],
            row["baseline"],
            row["metric"],
        )
        require(key not in seen, "Metric-delta cell duplicated")
        seen.add(key)
        require(row["reference_variant"] == "A1_B0_7_primary", "Reference variant changed")
        observed = float(aggregate[key[:3]][row["metric"]])
        reference = float(
            aggregate[("A1_B0_7_primary", row["scope"], row["baseline"])][
                row["metric"]
            ]
        )
        require(
            math.isclose(finite(row["value"], "delta value"), observed, abs_tol=1e-15),
            "Delta observed value mismatch",
        )
        require(
            math.isclose(
                finite(row["reference_value"], "reference value"),
                reference,
                abs_tol=1e-15,
            ),
            "Delta reference mismatch",
        )
        require(
            math.isclose(
                finite(row["absolute_difference"], "absolute difference"),
                observed - reference,
                abs_tol=1e-15,
            ),
            "Metric difference is inconsistent",
        )
        if row["weight_variant"] == "A1_B0_7_primary":
            require(
                finite(row["absolute_difference"], "primary difference") == 0.0,
                "Primary metric difference is nonzero",
            )
    require(len(seen) == 240, "Metric-delta matrix is incomplete")


def validate_rank_changes(rows: list[dict[str, str]]) -> None:
    require(len(rows) == 12, "Rank-change table is not 12 rows")
    require({row["weight_variant"] for row in rows} == VARIANTS, "Rank variants changed")
    require({row["baseline"] for row in rows} == BASELINES, "Rank baselines changed")
    for row in rows:
        require(int(row["eligible_candidate_rows"]) == 914_532, "Eligible row count changed")
        require(int(row["endpoint_relation_count"]) == 358, "Endpoint relation count changed")
        require(
            int(row["rank_permutation_blocks_checked"]) == 222,
            "Rank block count changed",
        )
        eligible = int(row["eligible_candidate_rows"])
        changed = int(row["rank_changed_candidate_count"])
        require(
            math.isclose(
                finite(row["rank_changed_candidate_fraction"], "rank fraction"),
                changed / eligible,
                abs_tol=1e-15,
            ),
            "Rank-change fraction is inconsistent",
        )
        if row["weight_variant"] == "A1_B0_7_primary":
            for field in (
                "score_changed_candidate_count",
                "rank_changed_candidate_count",
                "absolute_rank_change_sum" if "absolute_rank_change_sum" in row else "",
                "maximum_absolute_rank_change",
                "query_count_with_any_rank_change",
                "top50_symmetric_difference_membership_count",
                "query_count_with_any_top50_membership_change",
                "endpoint_relation_rank_changed_count",
                "endpoint_relation_top50_membership_changed_count",
            ):
                if field:
                    require(int(row[field]) == 0, f"Primary change is nonzero: {field}")


def validate_scope(rows: list[dict[str, str]]) -> None:
    require(len(rows) == 12, "Scope-invariance table is not 12 rows")
    for row in rows:
        require(
            int(row["relation_count_change_vs_0_7"]) == 0
            and int(row["query_count_change_vs_0_7"]) == 0
            and int(row["target_count_change_vs_0_7"]) == 0,
            "A scope cardinality changed with evidence weight",
        )


def validate_manifest(output_dir: Path, manifest: dict[str, Any]) -> None:
    require(manifest.get("aggregate_only") is True, "Manifest is not aggregate-only")
    require(
        manifest.get("identifier_bearing_output") is False,
        "Manifest asserts identifier-bearing output",
    )
    require(
        manifest.get("absolute_paths_recorded") is False,
        "Manifest asserts absolute paths",
    )
    outputs = manifest.get("outputs")
    require(isinstance(outputs, list) and len(outputs) == 7, "Manifest outputs changed")
    for item in outputs:
        path = output_dir / item["basename"]
        require(path.is_file(), f"Manifest output absent: {item['basename']}")
        require(path.stat().st_size == item["bytes"], "Manifest byte count mismatch")
        require(sha256(path) == item["sha256"], "Manifest output hash mismatch")


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
        "Output inventory changed",
    )
    _, weight_rows = read_tsv(output_dir / "weight_variants.tsv")
    _, aggregate_rows = read_tsv(output_dir / "aggregate_metrics.tsv")
    _, delta_rows = read_tsv(output_dir / "metric_deltas_vs_0_7.tsv")
    _, rank_rows = read_tsv(output_dir / "complete_rank_top50_changes_vs_0_7.tsv")
    _, scope_rows = read_tsv(output_dir / "scope_cardinality_invariance.tsv")
    validate_weight_rows(weight_rows)
    aggregate = validate_aggregate(aggregate_rows)
    validate_deltas(delta_rows, aggregate)
    validate_rank_changes(rank_rows)
    validate_scope(scope_rows)
    summary = json.loads((output_dir / "weight_policy_summary.json").read_text(encoding="utf-8"))
    require(
        summary.get("complete_score_rank_contract", {}).get("computed_score_rank_rows")
        == 10_974_384,
        "Complete computed row count changed",
    )
    require(
        summary.get("complete_score_rank_contract", {}).get(
            "primary_0_7_rows_reproduced_exactly"
        )
        == 3_658_128,
        "Primary exact reproduction count changed",
    )
    require(
        summary.get("complete_score_rank_contract", {}).get(
            "complete_rank_ledgers_written"
        )
        is False,
        "Complete identifier-bearing ranks were written",
    )
    manifest_path = output_dir / "run_manifest.json"
    validate_manifest(
        output_dir, json.loads(manifest_path.read_text(encoding="utf-8"))
    )
    validate_release_boundary(output_dir)
    checks = [
        "exact output inventory",
        "three frozen weight variants and all-equal alias",
        "48-cell aggregate metric matrix",
        "240-cell metric-delta matrix with exact arithmetic consistency",
        "12-cell complete-rank/top-50 change matrix",
        "primary 0.7 aggregate deltas are zero",
        "scope cardinality invariant across weights",
        "10,974,384 complete score/rank rows computed",
        "3,658,128 primary score/rank rows exactly reproduced",
        "manifest output hashes and byte counts",
        "no identifier-bearing score/rank columns or recognizable identifiers",
        "no absolute paths",
    ]
    receipt = {
        "schema_version": "1.0",
        "analysis_id": "revision_weight_policy_v1_20260729",
        "status": "PASS",
        "check_count": len(checks),
        "checks": checks,
        "validated_manifest": {
            "basename": manifest_path.name,
            "sha256": sha256(manifest_path),
        },
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
