#!/usr/bin/env python3
"""Assemble aggregate-only major-revision tables from validated analysis outputs."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Iterable


ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parent
TABLES = ROOT / "tables"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"Empty TSV: {path}")
    return rows


def write_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, object]]) -> None:
    materialized = list(rows)
    if not materialized:
        raise ValueError(f"Refusing to write empty TSV: {path}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        for row in materialized:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalized_records(
    source_path: Path,
    section: str,
    dimension_fields: tuple[str, ...],
) -> list[dict[str, object]]:
    rows = read_tsv(source_path)
    output: list[dict[str, object]] = []
    for source_row_index, row in enumerate(rows, start=1):
        dimensions = {field: row.get(field, "") for field in dimension_fields}
        for field, value in row.items():
            if field in dimension_fields:
                continue
            output.append(
                {
                    "section": section,
                    "source_row": source_row_index,
                    "scenario_or_policy": dimensions.get("scenario", "")
                    or dimensions.get("weight_variant", "")
                    or dimensions.get("policy", "")
                    or dimensions.get("cohort", ""),
                    "scope": dimensions.get("scope", "")
                    or dimensions.get("display_scope", "")
                    or dimensions.get("provenance_scope", ""),
                    "baseline": dimensions.get("baseline", ""),
                    "subgroup": dimensions.get("query_subset", "")
                    or dimensions.get("analysis_unit", "")
                    or dimensions.get("reason_category", "")
                    or dimensions.get("estimand", "")
                    or dimensions.get("interval_status", "")
                    or dimensions.get("item", "")
                    or dimensions.get("reference_variant", ""),
                    "item": field,
                    "value": value,
                    "status": row.get("status", ""),
                    "source_artifact": source_path.name,
                }
            )
    return output


def update_table_1() -> Path:
    path = TABLES / "Table_1_temporal_repair_flow.tsv"
    rows = [row for row in read_tsv(path) if row["section"] != "date_policy_sensitivity"]
    date_rows = [
        ("day_only_selected_rows", "13885"),
        ("interval_certain_selected_rows", "20455"),
        ("date_resolved_rows_definitely_before_or_on_cutoff", "20455"),
        ("date_interval_crossing_cutoff_rows", "0"),
        ("date_interval_definitely_after_cutoff_rows", "0"),
        ("date_unresolved_or_non_numeric_reference_rows", "192"),
        ("tier_B_to_A_changes_interval_vs_day_only", "141"),
        ("representative_structure_changes_interval_vs_day_only", "0"),
        ("scope_membership_changes_interval_vs_day_only", "0"),
    ]
    rows.extend(
        {
            "section": "date_policy_sensitivity",
            "item": item,
            "value": value,
            "status": "verified",
        }
        for item, value in date_rows
    )
    write_tsv(path, ["section", "item", "value", "status"], rows)
    return path


def build_table_4() -> Path:
    tie_path = PROJECT / "revision_tie_aware_v1_20260729" / "outputs" / "tie_aware_metrics.tsv"
    bound_path = PROJECT / "revision_tie_aware_v1_20260729" / "outputs" / "double_cold_query_hit_upper_bounds.tsv"
    tie_rows = read_tsv(tie_path)
    bounds = {
        (row["scope"], row["baseline"]): row
        for row in read_tsv(bound_path)
    }
    output: list[dict[str, object]] = []
    for row in tie_rows:
        if row["k"] != "50":
            continue
        if row["query_subset"] != "all_queries" and not (
            row["scope"] == "temporal_strict_ab"
            and row["baseline"] == "sequence_3mer_transfer"
        ):
            continue
        key = (row["scope"], row["baseline"])
        bound = bounds.get(key, {})
        tie_queries = int(row["recall_tie_dependent_query_count"])
        if row["query_subset"] == "structural_all_zero_non_operational":
            result_class = "non_operational_uniform_tie_allocation"
        elif tie_queries:
            result_class = "tie_dependent"
        else:
            result_class = "score_identifiable"
        output.append(
            {
                "scope": row["scope"],
                "baseline": row["baseline"],
                "query_subset": row["query_subset"],
                "query_count": row["query_count"],
                "relevant_relation_count": row["relevant_relation_count"],
                "salted_recall_at_50": row["legacy_salted_recall"],
                "tie_expected_recall_at_50": row["tie_expected_fractional_recall"],
                "tie_worst_recall_at_50": row["tie_worst_recall"],
                "tie_best_recall_at_50": row["tie_best_recall"],
                "query_bootstrap_ci95_low": row["query_bootstrap_expected_recall_ci95_low"],
                "query_bootstrap_ci95_high": row["query_bootstrap_expected_recall_ci95_high"],
                "pmid_component_ci95_low": row["pmid_component_expected_recall_ci95_low"],
                "pmid_component_ci95_high": row["pmid_component_expected_recall_ci95_high"],
                "score_identifiable_query_count": row["recall_score_identifiable_query_count"],
                "tie_dependent_query_count": row["recall_tie_dependent_query_count"],
                "empirical_zero_hit_status": bound.get("empirical_query_bootstrap_status", ""),
                "one_sided_cp95_upper": bound.get("one_sided_clopper_pearson_upper", ""),
                "result_class": result_class,
            }
        )
    path = TABLES / "Table_4_corrected_bootstrap_summaries.tsv"
    fields = [
        "scope",
        "baseline",
        "query_subset",
        "query_count",
        "relevant_relation_count",
        "salted_recall_at_50",
        "tie_expected_recall_at_50",
        "tie_worst_recall_at_50",
        "tie_best_recall_at_50",
        "query_bootstrap_ci95_low",
        "query_bootstrap_ci95_high",
        "pmid_component_ci95_low",
        "pmid_component_ci95_high",
        "score_identifiable_query_count",
        "tie_dependent_query_count",
        "empirical_zero_hit_status",
        "one_sided_cp95_upper",
        "result_class",
    ]
    write_tsv(path, fields, output)
    return path


def build_s7() -> Path:
    tie_path = PROJECT / "revision_tie_aware_v1_20260729" / "outputs" / "tie_aware_metrics.tsv"
    bound_path = PROJECT / "revision_tie_aware_v1_20260729" / "outputs" / "double_cold_query_hit_upper_bounds.tsv"
    operability_path = PROJECT / "revision_tie_aware_v1_20260729" / "outputs" / "three_mer_operability.tsv"
    rows: list[dict[str, object]] = []
    rows.extend(
        normalized_records(
            tie_path,
            "tie_aware_metric",
            ("scope", "baseline", "query_subset", "k", "status"),
        )
    )
    rows.extend(
        normalized_records(
            bound_path,
            "joint_cold_zero_hit_bound",
            ("scope", "baseline"),
        )
    )
    rows.extend(
        normalized_records(
            operability_path,
            "sequence_operability",
            ("scope",),
        )
    )
    path = TABLES / "Table_S7_tie_aware_retrieval.tsv"
    write_tsv(
        path,
        [
            "section",
            "source_row",
            "scenario_or_policy",
            "scope",
            "baseline",
            "subgroup",
            "item",
            "value",
            "status",
            "source_artifact",
        ],
        rows,
    )
    return path


def build_s8() -> Path:
    base = PROJECT / "revision_date_policy_v1_20260729" / "outputs"
    specifications = [
        ("interval_status", "interval_status_counts.tsv", ("interval_status",)),
        ("history", "scenario_history_summary.tsv", ("scenario",)),
        ("structure", "scenario_structure_summary.tsv", ("scenario",)),
        ("scope", "scope_denominators.tsv", ("scenario", "provenance_scope", "display_scope")),
        ("recall_at_50", "recall_at_50.tsv", ("scenario", "provenance_scope", "display_scope", "baseline")),
        ("score_rank_change", "score_rank_change_summary.tsv", ("scenario", "baseline")),
        ("scenario_equivalence", "scenario_equivalence.tsv", ("scenario",)),
    ]
    rows: list[dict[str, object]] = []
    for section, name, dimensions in specifications:
        rows.extend(normalized_records(base / name, section, dimensions))
    path = TABLES / "Table_S8_date_precision_policy.tsv"
    write_tsv(
        path,
        [
            "section",
            "source_row",
            "scenario_or_policy",
            "scope",
            "baseline",
            "subgroup",
            "item",
            "value",
            "status",
            "source_artifact",
        ],
        rows,
    )
    return path


def build_s9() -> Path:
    base = PROJECT / "revision_weight_policy_v1_20260729" / "outputs"
    specifications = [
        ("weight_variant", "weight_variants.tsv", ("weight_variant",)),
        ("weight_metric", "aggregate_metrics.tsv", ("weight_variant", "scope", "baseline")),
        (
            "weight_rank_change",
            "complete_rank_top50_changes_vs_0_7.tsv",
            ("weight_variant", "reference_variant", "baseline"),
        ),
        (
            "weight_scope_invariance",
            "scope_cardinality_invariance.tsv",
            ("weight_variant", "scope"),
        ),
    ]
    rows: list[dict[str, object]] = []
    for section, name, dimensions in specifications:
        rows.extend(normalized_records(base / name, section, dimensions))

    structure_base = PROJECT / "revision_structure_policy_v1_20260729"
    if structure_base.is_dir():
        for source_name in (
            "structure_policy_summary.tsv",
            "scaffold_scope_changes.tsv",
            "rank_change_summary.tsv",
            "scope_recall_at_50.tsv",
        ):
            source = structure_base / source_name
            rows.extend(
                normalized_records(
                    source,
                    f"structure_{source.stem}",
                    ("policy", "scope", "baseline", "analysis_unit", "status"),
                )
            )

    path = TABLES / "Table_S9_weight_and_structure_policy.tsv"
    write_tsv(
        path,
        [
            "section",
            "source_row",
            "scenario_or_policy",
            "scope",
            "baseline",
            "subgroup",
            "item",
            "value",
            "status",
            "source_artifact",
        ],
        rows,
    )
    return path


def build_s10() -> Path:
    base = PROJECT / "revision_unresolved_bounds_v1_20260729" / "outputs"
    specifications = [
        ("mapping_failure", "mapping_failure_strata.tsv", ("reason_category",)),
        ("cohort_comparability", "cohort_comparability.tsv", ("cohort",)),
        ("comparability_difference", "comparability_differences.tsv", ("item",)),
        ("endpoint_cardinality", "endpoint_cardinality_bounds.tsv", ()),
        ("relation_top50_bound", "relation_level_top50_bounds.tsv", ("baseline",)),
        ("non_identifiable", "non_identifiable_estimands.tsv", ("estimand",)),
    ]
    rows: list[dict[str, object]] = []
    for section, name, dimensions in specifications:
        rows.extend(normalized_records(base / name, section, dimensions))
    path = TABLES / "Table_S10_unresolved_entity_bounds.tsv"
    write_tsv(
        path,
        [
            "section",
            "source_row",
            "scenario_or_policy",
            "scope",
            "baseline",
            "subgroup",
            "item",
            "value",
            "status",
            "source_artifact",
        ],
        rows,
    )
    return path


def build_s11() -> Path:
    source = PROJECT / "revision_review_matrix_v1_20260729" / "max_similarity_summary.tsv"
    rows = read_tsv(source)
    path = TABLES / "Table_S11_maximum_similarity_distributions.tsv"
    write_tsv(path, list(rows[0]), rows)
    return path


def main() -> int:
    update_table_1()
    build_table_4()
    build_s7()
    build_s8()
    build_s9()
    build_s10()
    build_s11()
    outputs = [
        TABLES / "Table_1_temporal_repair_flow.tsv",
        TABLES / "Table_2_historical_before_after_audit.tsv",
        TABLES / "Table_3_corrected_aggregate_performance.tsv",
        TABLES / "Table_4_corrected_bootstrap_summaries.tsv",
        TABLES / "Table_5_score_degeneracy_and_ties.tsv",
        TABLES / "Table_6_claim_evidence_use_boundaries.tsv",
        TABLES / "Table_S1_scope_mask_integrity.tsv",
        TABLES / "Table_S2_top100_exhaustive_fidelity.tsv",
        TABLES / "Table_S3_zero_and_failure_accounting.tsv",
        TABLES / "Table_S4_pmid_document_dependence.tsv",
        TABLES / "Table_S5_frozen_unresolved_exclusions.tsv",
        TABLES / "Table_S6_reproducibility_and_release.tsv",
        TABLES / "Table_S7_tie_aware_retrieval.tsv",
        TABLES / "Table_S8_date_precision_policy.tsv",
        TABLES / "Table_S9_weight_and_structure_policy.tsv",
        TABLES / "Table_S10_unresolved_entity_bounds.tsv",
        TABLES / "Table_S11_maximum_similarity_distributions.tsv",
        TABLES / "Table_S12_rights_and_controlled_access.tsv",
    ]
    for path in outputs:
        if not path.is_file():
            raise FileNotFoundError(path)
    manifest = {
        "schema_version": "major_revision_table_assembly_v4",
        "aggregate_only": True,
        "outputs": [
            {
                "path": path.relative_to(ROOT).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
            for path in outputs
        ],
    }
    manifest_path = TABLES / "major_revision_table_assembly_manifest_v4.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"status": "PASS", "outputs": len(outputs)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
