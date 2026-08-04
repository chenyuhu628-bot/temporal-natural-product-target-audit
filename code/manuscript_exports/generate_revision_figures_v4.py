#!/usr/bin/env python3
"""Generate aggregate-only PNG figures for the major-revision manuscript v4."""

from __future__ import annotations

import csv
import hashlib
import json
from datetime import date
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Patch, Rectangle
from matplotlib.ticker import PercentFormatter
from PIL import Image


BASELINES = [
    "weighted_target_popularity",
    "sequence_3mer_transfer",
    "weighted_morgan_transfer",
    "structure_sequence_pair_neighbor",
]
BASELINE_LABELS = {
    "weighted_target_popularity": "Target popularity",
    "sequence_3mer_transfer": "Sequence 3-mer",
    "weighted_morgan_transfer": "Morgan transfer",
    "structure_sequence_pair_neighbor": "Pair neighbour",
}
BASELINE_SHORT = {
    "weighted_target_popularity": "Popularity",
    "sequence_3mer_transfer": "3-mer",
    "weighted_morgan_transfer": "Morgan",
    "structure_sequence_pair_neighbor": "Pair",
}
FIGURE_PALETTE = {
    "ink": "#24333D",
    "blue": "#277DA1",
    "pale_blue": "#DCEBF2",
    "slate": "#5B7083",
    "teal": "#2A9D8F",
    "purple": "#B279A2",
    "coral": "#E76F51",
    "amber": "#F4A261",
    "mid_gray": "#B9BFC3",
    "pale_gray": "#E9EDF0",
    "grid": "#DCE3E8",
}
_ARIAL_NARROW_PATH = Path("C:/Windows/Fonts/arialn.ttf")
MATRIX_CELL_FONT = (
    FontProperties(fname=str(_ARIAL_NARROW_PATH), size=7.3, weight="bold")
    if _ARIAL_NARROW_PATH.is_file()
    else FontProperties(family=["Arial", "DejaVu Sans"], stretch="condensed", size=7.3, weight="bold")
)
SCOPE_CELL_FONT = (
    FontProperties(fname=str(_ARIAL_NARROW_PATH), size=7.4)
    if _ARIAL_NARROW_PATH.is_file()
    else FontProperties(family=["Arial", "DejaVu Sans"], stretch="condensed", size=7.4)
)
COLORS = {
    "weighted_target_popularity": FIGURE_PALETTE["slate"],
    "sequence_3mer_transfer": FIGURE_PALETTE["purple"],
    "weighted_morgan_transfer": FIGURE_PALETTE["teal"],
    "structure_sequence_pair_neighbor": FIGURE_PALETTE["coral"],
}
DISPLAY_SCOPES = [
    "temporal_strict_ab",
    "scaffold_cold_strict_ab",
    "project_defined_joint_scaffold_homology_cold_0_30",
    "project_defined_joint_scaffold_homology_cold_0_50_0_70_identical_mask",
]
SCOPE_LABELS = {
    "temporal_strict_ab": "Temporal strict A/B",
    "scaffold_cold_strict_ab": "Scaffold cold",
    "project_defined_joint_scaffold_homology_cold_0_30": "Joint scaffold–homology 0.30",
    "project_defined_joint_scaffold_homology_cold_0_50_0_70_identical_mask": "Joint scaffold–homology 0.50/0.70",
}
METRICS = ["Recall@10", "Recall@50", "NDCG@10", "NDCG@50", "MRR"]


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"Empty table: {path}")
    return rows


def require_unique_row(rows: list[dict[str, str]], label: str, **criteria: str) -> dict[str, str]:
    matches = [row for row in rows if all(row.get(key, "") == value for key, value in criteria.items())]
    if len(matches) != 1:
        rendered = ", ".join(f"{key}={value!r}" for key, value in criteria.items())
        raise ValueError(f"{label}: expected one row for {rendered}; found {len(matches)}")
    return matches[0]


def read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a JSON object: {path}")
    return payload


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def save_rgb_png(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, dpi=600, bbox_inches="tight", facecolor="white", format="png")
    plt.close(fig)
    with Image.open(path) as image:
        image.convert("RGB").save(path, format="PNG", dpi=(600, 600), optimize=True)


def save_rgb_png_fixed_canvas(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, dpi=600, facecolor="white", format="png")
    plt.close(fig)
    with Image.open(path) as image:
        image.convert("RGB").save(path, format="PNG", dpi=(600, 600), optimize=True)


def save_rgb_png_and_pdf(fig: plt.Figure, png_path: Path, pdf_path: Path) -> None:
    fig.savefig(pdf_path, bbox_inches="tight", facecolor="white", format="pdf")
    save_rgb_png(fig, png_path)
    if not pdf_path.is_file() or pdf_path.stat().st_size == 0 or not pdf_path.read_bytes().startswith(b"%PDF"):
        raise ValueError(f"PDF output contract failed for {pdf_path.name}")


def set_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": 8.0,
            "axes.titlesize": 9.2,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "text.color": FIGURE_PALETTE["ink"],
            "axes.labelcolor": FIGURE_PALETTE["ink"],
            "xtick.color": FIGURE_PALETTE["ink"],
            "ytick.color": FIGURE_PALETTE["ink"],
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(0.0, 1.02, label.lower(), transform=ax.transAxes, fontsize=10.5, fontweight="bold", va="bottom", color=FIGURE_PALETTE["ink"])


def figure1(project: Path, figures: Path) -> tuple[Path, Path, list[dict[str, str]]]:
    """Draw the current corrective lineage as a linear, aggregate-only workflow."""

    rebuild_path = project / "author_run_strict_ab_asof_cutoff_execution_v1_20260728" / "audit" / "asof_rebuild_summary.json"
    future_path = project / "results" / "strict_temporal_future_v1_1_pmid_verified_chembl31_leakage_gate_summary.json"
    source_path = project / "author_run_strict_ab_asof_cutoff_execution_v1_20260728" / "audit" / "source_concentration_v1" / "source_concentration_aggregate_summary.json"
    score_path = project / "author_run_strict_ab_asof_cutoff_execution_v1_20260728" / "score" / "corrective_score_manifest.json"
    legacy_c37_summary_path = project / "results" / "chembl37_p2_source_overlap_audit_v1_summary.json"
    legacy_c37_config_path = project / "configs" / "chembl37_p2_manual_review_queue_v1_1.json"
    date_status_path = project / "revision_date_policy_v1_20260729" / "outputs" / "interval_status_counts.tsv"
    date_history_path = project / "revision_date_policy_v1_20260729" / "outputs" / "scenario_history_summary.tsv"

    rebuild = read_json(rebuild_path)
    future = read_json(future_path)
    source = read_json(source_path)
    score = read_json(score_path)
    legacy_c37_summary = read_json(legacy_c37_summary_path)
    legacy_c37_config = read_json(legacy_c37_config_path)
    date_status = {row["interval_status"]: row for row in read_tsv(date_status_path)}
    date_history = {row["scenario"]: row for row in read_tsv(date_history_path)}

    counts = rebuild["counts"]
    rows = rebuild["row_eligibility_counts"]
    endpoint_status = future["status_counts_in_frozen_future_table"]
    endpoint_cohort = next(item for item in source["cohorts"] if item["cohort"] == "endpoint")
    historical_cohort = next(item for item in source["cohorts"] if item["cohort"] == "historical")
    overlap = source["cross_cohort_source_overlap"][0]
    day_history = date_history["day_only_conservative"]
    interval_history = date_history["interval_certain_pre_cutoff"]
    before_status = date_status["definitely_before_or_on"]
    crossing_status = date_status["crossing_cutoff"]
    after_status = date_status["definitely_after"]
    unresolved_status = date_status["unresolved_interval"]
    day_rows = int(day_history["selected_source_row_count"])
    interval_rows = int(interval_history["selected_source_row_count"])
    non_day_before_rows = int(before_status["month_precision_count"]) + int(before_status["year_precision_count"])
    tier_a_upgrades = int(interval_history["tier_A_pair_count"]) - int(day_history["tier_A_pair_count"])
    historical_source_rows = int(counts["historical_strict_v2_rows"])
    unresolved_rows = int(unresolved_status["source_row_count"])
    initial_endpoint_candidates = int(future["inputs"]["frozen_future_pair_count"])
    excluded_historical_overlap = int(endpoint_status["historical_activity_recorded_in_chembl31"])
    excluded_entity_unresolved = int(endpoint_status["entity_pair_not_sqlite_validated"])
    after_historical_overlap = initial_endpoint_candidates - excluded_historical_overlap
    endpoint_relations = int(counts["endpoint_relations"])
    endpoint_queries = int(counts["endpoint_queries"])
    endpoint_pmids = int(endpoint_cohort["unique_source_document_count"])
    component_summary = endpoint_cohort["query_source_component_summary"]
    largest_component_queries = int(component_summary["largest_component_left_node_count"])
    largest_component_pmids = int(component_summary["largest_component_source_document_count"])
    remaining_component_queries = endpoint_queries - largest_component_queries
    remaining_component_pmids = endpoint_pmids - largest_component_pmids
    largest_query_fraction = largest_component_queries / endpoint_queries
    largest_pmid_fraction = largest_component_pmids / endpoint_pmids

    if sum(int(value) for value in rows.values()) != int(counts["historical_strict_v2_rows"]):
        raise ValueError("Figure 1 row-eligibility counts do not sum to the historical source-row total")
    if int(future["inputs"]["frozen_future_pair_count"]) != 442:
        raise ValueError("Unexpected frozen future pair count")
    if 442 - int(endpoint_status["historical_activity_recorded_in_chembl31"]) - int(endpoint_status["entity_pair_not_sqlite_validated"]) != int(counts["endpoint_relations"]):
        raise ValueError("Figure 1 endpoint flow does not close arithmetically")
    if int(endpoint_status["no_activity_found_in_validated_chembl31_entity_pair"]) != int(counts["endpoint_relations"]):
        raise ValueError("Figure 1 endpoint count disagrees across frozen sources")
    if score["baselines"] != BASELINES or int(score["target_count"]) != int(counts["candidate_targets"]):
        raise ValueError("Figure 1 score-manifest contract mismatch")
    if int(overlap["shared_source_document_count"]) != 0:
        raise ValueError("Figure 1 current-cohort PMID overlap is no longer zero")
    if int(legacy_c37_config["expected_frozen_pair_counts"]["all"]) != 846 or int(legacy_c37_summary["pairs_with_at_least_one_shared_PMID"]) != 835:
        raise ValueError("Legacy C37 lineage receipt changed; review exclusion provenance")
    if day_rows != 13885 or interval_rows != 20455 or interval_rows - day_rows != 6570:
        raise ValueError("Figure 1 date-policy row-count contract changed")
    if non_day_before_rows != 6570 or int(crossing_status["source_row_count"]) != 0 or int(after_status["source_row_count"]) != 0:
        raise ValueError("Figure 1 interval classification contract changed")
    if int(before_status["source_row_count"]) != interval_rows or int(unresolved_status["source_row_count"]) != 192:
        raise ValueError("Figure 1 date-interval totals no longer close")
    if tier_a_upgrades != 141 or int(interval_history["tier_or_weight_changed_pair_count_vs_day_only"]) != 141:
        raise ValueError("Figure 1 B-to-A upgrade contract changed")
    if int(day_history["historical_pair_count"]) != int(interval_history["historical_pair_count"]) or int(day_history["historical_pair_count"]) != int(counts["historical_pairs"]):
        raise ValueError("Figure 1 history membership is no longer invariant across date policies")
    if day_rows + non_day_before_rows + unresolved_rows != historical_source_rows:
        raise ValueError("Figure 1 displayed date-policy partition does not sum to the historical source-row total")
    if after_historical_overlap - excluded_entity_unresolved != endpoint_relations:
        raise ValueError("Figure 1 displayed endpoint-filtering stages do not close arithmetically")
    if not np.isclose(largest_query_fraction, float(component_summary["largest_component_left_node_fraction"])):
        raise ValueError("Figure 1 largest query-component fraction disagrees with its source summary")
    if not np.isclose(largest_pmid_fraction, float(component_summary["largest_component_source_document_fraction"])):
        raise ValueError("Figure 1 largest PMID-component fraction disagrees with its source summary")

    def rel(path: Path) -> str:
        return path.relative_to(project).as_posix()

    source_rows: list[dict[str, Any]] = [
        {"panel": "A", "item": "cutoff", "value": rebuild["cutoff"], "status": "current_corrective_branch", "upstream_artifact": rel(rebuild_path)},
        {"panel": "A", "item": "historical_strict_v2_source_rows", "value": counts["historical_strict_v2_rows"], "status": "verified", "upstream_artifact": rel(rebuild_path)},
        {"panel": "A", "item": "day_only_selected_rows", "value": day_rows, "status": "primary_date_policy", "upstream_artifact": rel(date_history_path)},
        {"panel": "A", "item": "interval_certain_selected_rows", "value": interval_rows, "status": "post_hoc_date_policy_sensitivity", "upstream_artifact": rel(date_history_path)},
        {"panel": "A", "item": "non_day_precision_definitely_pre_cutoff_rows", "value": non_day_before_rows, "status": "all_definitely_before_or_on_cutoff", "upstream_artifact": rel(date_status_path)},
        {"panel": "A", "item": "crossing_cutoff_rows", "value": crossing_status["source_row_count"], "status": "interval_audit", "upstream_artifact": rel(date_status_path)},
        {"panel": "A", "item": "definitely_after_cutoff_rows", "value": after_status["source_row_count"], "status": "interval_audit", "upstream_artifact": rel(date_status_path)},
        {"panel": "A", "item": "unresolved_interval_rows", "value": unresolved_status["source_row_count"], "status": "excluded_under_both_policies", "upstream_artifact": rel(date_status_path)},
        {"panel": "A", "item": "relation_keys_in_both_date_policies", "value": counts["historical_pairs"], "status": "membership_invariant", "upstream_artifact": rel(date_history_path)},
        {"panel": "A", "item": "historical_targets_in_both_date_policies", "value": counts["historical_targets"], "status": "membership_invariant", "upstream_artifact": rel(date_history_path)},
        {"panel": "A", "item": "tier_B_to_A_upgrades", "value": tier_a_upgrades, "status": "interval_certain_vs_day_only", "upstream_artifact": rel(date_history_path)},
        {"panel": "A", "item": "day_only_source_row_percent", "value": f"{100 * day_rows / historical_source_rows:.1f}", "status": "derived_for_display", "upstream_artifact": rel(date_history_path)},
        {"panel": "A", "item": "additional_interval_certain_source_row_percent", "value": f"{100 * non_day_before_rows / historical_source_rows:.1f}", "status": "derived_for_display", "upstream_artifact": rel(date_status_path)},
        {"panel": "A", "item": "unresolved_source_row_percent", "value": f"{100 * unresolved_rows / historical_source_rows:.1f}", "status": "derived_for_display", "upstream_artifact": rel(date_status_path)},
        {"panel": "B", "item": "initial_later_candidate_pairs", "value": future["inputs"]["frozen_future_pair_count"], "status": "frozen_prior_endpoint_lineage", "upstream_artifact": rel(future_path)},
        {"panel": "B", "item": "excluded_C31_historical_activity", "value": endpoint_status["historical_activity_recorded_in_chembl31"], "status": "excluded", "upstream_artifact": rel(future_path)},
        {"panel": "B", "item": "remaining_after_historical_activity_filter", "value": after_historical_overlap, "status": "derived_for_display", "upstream_artifact": rel(future_path)},
        {"panel": "B", "item": "excluded_entity_unresolved", "value": endpoint_status["entity_pair_not_sqlite_validated"], "status": "excluded", "upstream_artifact": rel(future_path)},
        {"panel": "B", "item": "final_endpoint_relations", "value": counts["endpoint_relations"], "status": "frozen", "upstream_artifact": rel(rebuild_path)},
        {"panel": "B", "item": "final_endpoint_queries", "value": counts["endpoint_queries"], "status": "frozen", "upstream_artifact": rel(rebuild_path)},
        {"panel": "B", "item": "final_endpoint_targets", "value": counts["endpoint_targets"], "status": "frozen", "upstream_artifact": rel(rebuild_path)},
        {"panel": "C", "item": "candidate_target_universe", "value": counts["candidate_targets"], "status": "fixed", "upstream_artifact": rel(score_path)},
        {"panel": "C", "item": "baselines", "value": "; ".join(score["baselines"]), "status": "fixed", "upstream_artifact": rel(score_path)},
        {"panel": "C", "item": "same_query_historical_target_mask", "value": score["parameters"]["mask_historical_targets_for_same_query"], "status": "fixed", "upstream_artifact": rel(score_path)},
        {"panel": "C", "item": "endpoint_read_by_score_engine", "value": score["endpoint_read_by_score_engine"], "status": "evaluation_only", "upstream_artifact": rel(score_path)},
        {"panel": "D", "item": "endpoint_evidence_rows", "value": endpoint_cohort["evidence_row_count"], "status": "current_corrective_branch", "upstream_artifact": rel(source_path)},
        {"panel": "D", "item": "endpoint_numeric_pmids", "value": endpoint_cohort["unique_source_document_count"], "status": "current_corrective_branch", "upstream_artifact": rel(source_path)},
        {"panel": "D", "item": "endpoint_query_pmid_components", "value": endpoint_cohort["query_source_component_summary"]["component_count"], "status": "current_corrective_branch", "upstream_artifact": rel(source_path)},
        {"panel": "D", "item": "largest_component_queries", "value": endpoint_cohort["query_source_component_summary"]["largest_component_left_node_count"], "status": "current_corrective_branch", "upstream_artifact": rel(source_path)},
        {"panel": "D", "item": "remaining_queries_outside_largest_component", "value": remaining_component_queries, "status": "derived_for_display", "upstream_artifact": rel(source_path)},
        {"panel": "D", "item": "largest_component_query_percent", "value": f"{100 * largest_query_fraction:.1f}", "status": "derived_for_display", "upstream_artifact": rel(source_path)},
        {"panel": "D", "item": "largest_component_endpoint_pmids", "value": largest_component_pmids, "status": "current_corrective_branch", "upstream_artifact": rel(source_path)},
        {"panel": "D", "item": "remaining_endpoint_pmids_outside_largest_component", "value": remaining_component_pmids, "status": "derived_for_display", "upstream_artifact": rel(source_path)},
        {"panel": "D", "item": "largest_component_endpoint_pmid_percent", "value": f"{100 * largest_pmid_fraction:.1f}", "status": "derived_for_display", "upstream_artifact": rel(source_path)},
        {"panel": "D", "item": "historical_numeric_pmids", "value": historical_cohort["unique_source_document_count"], "status": "current_corrective_branch", "upstream_artifact": rel(source_path)},
        {"panel": "D", "item": "shared_pmids_current_endpoint_vs_history", "value": overlap["shared_source_document_count"], "status": "provenance_only_not_external_validation", "upstream_artifact": rel(source_path)},
        {"panel": "not_displayed", "item": "legacy_C37_shared_PMID_candidates", "value": legacy_c37_summary["pairs_with_at_least_one_shared_PMID"], "status": "removed_from_figure_broader_846_pair_lineage_not_current_358_endpoint", "upstream_artifact": rel(legacy_c37_summary_path)},
    ]
    source_data_path = figures / "Figure_1_source_data.tsv"
    write_tsv(source_data_path, ["panel", "item", "value", "status", "upstream_artifact"], source_rows)

    set_style()
    with plt.rc_context(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": 8.0,
            "axes.titlesize": 9.2,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 8.0,
            "ytick.labelsize": 8.0,
            "axes.linewidth": 0.8,
            "hatch.linewidth": 0.55,
        }
    ):
        fig = plt.figure(figsize=(7.1, 4.6), facecolor="white")
        outer = fig.add_gridspec(2, 1, left=0.075, right=0.985, top=0.935, bottom=0.075, hspace=0.58, height_ratios=[1.0, 1.03])
        top = outer[0].subgridspec(1, 2, wspace=0.32)
        bottom = outer[1].subgridspec(1, 2, width_ratios=[1.70, 1.0], wspace=0.27)
        ax_a = fig.add_subplot(top[0, 0])
        ax_b = fig.add_subplot(top[0, 1])
        ax_c = fig.add_subplot(bottom[0, 0])
        ax_d = fig.add_subplot(bottom[0, 1])

        ink = FIGURE_PALETTE["ink"]
        blue = FIGURE_PALETTE["blue"]
        pale_blue = FIGURE_PALETTE["pale_blue"]
        coral = FIGURE_PALETTE["coral"]
        amber = FIGURE_PALETTE["amber"]
        pale_gray = FIGURE_PALETTE["pale_gray"]
        mid_gray = FIGURE_PALETTE["mid_gray"]
        line_gray = FIGURE_PALETTE["mid_gray"]
        grid_gray = FIGURE_PALETTE["grid"]
        muted = FIGURE_PALETTE["slate"]
        lw = 0.8

        def panel_heading(ax: plt.Axes, letter: str, title: str) -> None:
            ax.text(0.0, 1.07, letter.lower(), transform=ax.transAxes, ha="left", va="bottom", fontsize=10.5, fontweight="bold", color=ink)
            ax.text(0.075, 1.07, title, transform=ax.transAxes, ha="left", va="bottom", fontsize=9.2, fontweight="bold", color=ink)

        # a | Mutually exclusive source-row partition under the date policy.
        panel_heading(ax_a, "a", "Date-policy partition")
        ax_a.set_xlim(0, historical_source_rows)
        ax_a.set_ylim(-0.82, 1.03)
        bar_y = 0.28
        bar_height = 0.29
        cutoff_value = date.fromisoformat(str(rebuild["cutoff"]))
        cutoff_label = f"{cutoff_value.day} {cutoff_value.strftime('%B %Y')}"
        ax_a.barh(bar_y, day_rows, left=0, height=bar_height, color=blue, edgecolor="white", linewidth=lw)
        ax_a.barh(bar_y, non_day_before_rows, left=day_rows, height=bar_height, color=pale_blue, edgecolor=blue, linewidth=0.55, hatch="//")
        ax_a.barh(bar_y, unresolved_rows, left=interval_rows, height=bar_height, color=pale_gray, edgecolor=ink, linewidth=0.55, hatch="..")
        ax_a.text(0, 0.92, f"NPASS v2 strict A/B source rows (n = {historical_source_rows:,})", ha="left", va="center", fontsize=8.2, fontweight="bold", color=ink)
        ax_a.text(0, 0.73, f"Cutoff: {cutoff_label}", ha="left", va="center", fontsize=8.0, color=ink)
        ax_a.text(day_rows / 2, bar_y, f"Day-only eligible\n{day_rows:,} ({100 * day_rows / historical_source_rows:.1f}%)", ha="center", va="center", fontsize=8.0, fontweight="bold", color="white", linespacing=1.15)
        ax_a.text(day_rows + non_day_before_rows / 2, bar_y, f"{non_day_before_rows:,} ({100 * non_day_before_rows / historical_source_rows:.1f}%)", ha="center", va="center", fontsize=8.0, fontweight="bold", color=ink)
        ax_a.annotate(
            "Additional interval-certain\nmonth/year rows",
            xy=(day_rows + non_day_before_rows / 2, bar_y - bar_height / 2),
            xytext=(interval_rows - historical_source_rows * 0.01, 0.07),
            ha="right",
            va="top",
            fontsize=8.0,
            color=ink,
            linespacing=1.0,
            arrowprops={"arrowstyle": "-", "color": line_gray, "linewidth": lw, "shrinkA": 2, "shrinkB": 1},
        )
        ax_a.annotate(
            f"Unresolved\n{unresolved_rows:,} ({100 * unresolved_rows / historical_source_rows:.1f}%)",
            xy=(interval_rows + unresolved_rows / 2, bar_y + bar_height / 2),
            xytext=(historical_source_rows, 0.64),
            ha="right",
            va="center",
            fontsize=8.0,
            color=ink,
            arrowprops={"arrowstyle": "-", "color": line_gray, "linewidth": lw, "shrinkA": 2, "shrinkB": 1},
        )
        ax_a.text(0, -0.43, f"Interval-certain sensitivity = day-only primary +\nadditional rows: {interval_rows:,}; crossing/after cutoff = {int(crossing_status['source_row_count']) + int(after_status['source_row_count'])}", ha="left", va="center", fontsize=8.0, fontweight="bold", color=ink, linespacing=1.1)
        ax_a.text(0, -0.73, f"Same {int(counts['historical_pairs']):,} relations and {int(counts['historical_targets']):,} targets; {tier_a_upgrades:,} Tier B→A upgrades", ha="left", va="center", fontsize=8.0, color=ink)
        ax_a.axis("off")

        # b | Remaining relations after each endpoint eligibility filter.
        panel_heading(ax_b, "b", "Endpoint filtering")
        stage_x = np.arange(3)
        stage_values = [initial_endpoint_candidates, after_historical_overlap, endpoint_relations]
        ax_b.bar(stage_x[0], stage_values[0], width=0.54, color=mid_gray, edgecolor="white", linewidth=lw, zorder=3)
        ax_b.bar(stage_x[1], stage_values[1], width=0.54, color=pale_gray, edgecolor="white", linewidth=lw, zorder=3)
        ax_b.bar(stage_x[1], excluded_historical_overlap, bottom=stage_values[1], width=0.54, color=mid_gray, edgecolor="white", linewidth=lw, hatch="//", zorder=3)
        ax_b.bar(stage_x[2], stage_values[2], width=0.54, color=coral, edgecolor="white", linewidth=lw, zorder=3)
        ax_b.bar(stage_x[2], excluded_entity_unresolved, bottom=stage_values[2], width=0.54, color=mid_gray, edgecolor="white", linewidth=lw, hatch="//", zorder=3)
        ax_b.plot([stage_x[0] + 0.27, stage_x[1] - 0.27], [stage_values[0], stage_values[0]], color=line_gray, linewidth=lw, linestyle=(0, (2, 2)), zorder=2)
        ax_b.plot([stage_x[1] + 0.27, stage_x[2] - 0.27], [stage_values[1], stage_values[1]], color=line_gray, linewidth=lw, linestyle=(0, (2, 2)), zorder=2)
        ax_b.text(stage_x[0], stage_values[0] - 17, f"{stage_values[0]:,}", ha="center", va="top", fontsize=8.5, fontweight="bold", color=ink)
        ax_b.text(stage_x[1], stage_values[1] - 16, f"{stage_values[1]:,}", ha="center", va="top", fontsize=8.5, fontweight="bold", color=ink)
        ax_b.text(stage_x[2], stage_values[2] - 16, f"{stage_values[2]:,}", ha="center", va="top", fontsize=8.5, fontweight="bold", color="white")
        ax_b.text(0.5, 520, f"−{excluded_historical_overlap:,} historical-\nactivity overlaps", ha="center", va="top", fontsize=8.0, color=muted, linespacing=1.05)
        ax_b.annotate(
            f"−{excluded_entity_unresolved:,} entity-\nunresolved pairs",
            xy=(stage_x[2] + 0.12, stage_values[2] + excluded_entity_unresolved / 2),
            xytext=(2.92, 525),
            ha="right",
            va="top",
            fontsize=8.0,
            color=muted,
            linespacing=1.05,
            arrowprops={"arrowstyle": "-", "color": line_gray, "linewidth": lw, "shrinkA": 2, "shrinkB": 2},
        )
        ax_b.plot([2.27, 2.34], [63, 63], color=line_gray, linewidth=lw)
        ax_b.text(2.38, 63, f"{int(counts['endpoint_queries']):,} queries\n{int(counts['endpoint_targets']):,} targets", ha="left", va="center", fontsize=8.0, fontweight="bold", color=ink, linespacing=1.1)
        ax_b.set_xlim(-0.5, 3.0)
        ax_b.set_ylim(0, 540)
        ax_b.set_ylabel("Candidate pairs", labelpad=3, color=ink)
        ax_b.set_xticks(stage_x, ["Initial\ncandidates", "After overlap\nscreen", "Frozen\nendpoint"])
        ax_b.set_yticks([0, 100, 200, 300, 400, 500])
        ax_b.tick_params(axis="x", length=0, pad=3, colors=ink)
        ax_b.tick_params(axis="y", width=lw, length=3, colors=ink)
        ax_b.yaxis.grid(True, color=grid_gray, linewidth=0.65, zorder=0)
        ax_b.spines["left"].set_color(line_gray)
        ax_b.spines["bottom"].set_color(line_gray)

        # c | Fixed target-universe retrieval and evaluation sequence.
        panel_heading(ax_c, "c", "Fixed retrieval protocol")
        ax_c.set_xlim(0, 1)
        ax_c.set_ylim(0, 1)
        node_x = [0.03, 0.23, 0.49, 0.72, 0.97]
        process_y = 0.58
        for start, end in zip(node_x[:-1], node_x[1:]):
            ax_c.add_patch(FancyArrowPatch((start + 0.018, process_y), (end - 0.018, process_y), arrowstyle="-|>", mutation_scale=7.5, linewidth=lw, color=muted))
        ax_c.scatter(node_x, [process_y] * len(node_x), s=18, facecolors="white", edgecolors=ink, linewidths=lw, zorder=3)
        ax_c.text(0.005, 0.72, "Query q", ha="left", va="bottom", fontsize=8.0, fontweight="bold", color=ink)
        ax_c.text(node_x[1], 0.48, "Mask corrected\nhistorical targets", ha="center", va="top", fontsize=8.0, color=ink, linespacing=1.05)
        ax_c.text(node_x[2], 0.72, f"Rank fixed universe of\n{int(counts['candidate_targets']):,} human single-protein targets", ha="center", va="bottom", fontsize=8.0, fontweight="bold", color=ink, linespacing=1.0)
        ax_c.text(node_x[3], 0.48, "Complete\nranked list", ha="center", va="top", fontsize=8.0, color=ink, linespacing=1.05)
        ax_c.text(0.995, 0.77, "Evaluate\nlater-recorded\nstrict A/B targets", ha="right", va="bottom", fontsize=8.0, color=ink, linespacing=1.0)
        ax_c.text(0.995, 0.68, "evaluation only", ha="right", va="bottom", fontsize=8.0, fontweight="bold", color=coral)
        ax_c.plot([node_x[2], node_x[2]], [process_y - 0.02, 0.33], color=line_gray, linewidth=lw)
        ax_c.plot([0.28, 0.70], [0.33, 0.33], color=line_gray, linewidth=lw)
        ax_c.text(node_x[2], 0.255, "Four fixed baselines", ha="center", va="center", fontsize=8.0, fontweight="bold", color=ink)
        baseline_labels = [BASELINE_LABELS[name] for name in score["baselines"]]
        for label, (x_pos, y_pos) in zip(baseline_labels, [(0.35, 0.15), (0.64, 0.15), (0.35, 0.05), (0.64, 0.05)]):
            ax_c.text(x_pos, y_pos, label, ha="center", va="center", fontsize=8.0, color=ink)
        ax_c.axis("off")

        # d | Concentration of query and PMID nodes in the largest component.
        panel_heading(ax_d, "d", "Evidence dependence")
        y_positions = [1.72, 0.95]
        largest_fractions = [largest_query_fraction, largest_pmid_fraction]
        largest_counts = [largest_component_queries, largest_component_pmids]
        remaining_counts = [remaining_component_queries, remaining_component_pmids]
        totals = [endpoint_queries, endpoint_pmids]
        for y_pos, largest_fraction, largest_count, remaining_count, total in zip(y_positions, largest_fractions, largest_counts, remaining_counts, totals):
            ax_d.barh(y_pos, largest_fraction, left=0, height=0.32, color=amber, edgecolor="white", linewidth=lw)
            ax_d.barh(y_pos, 1 - largest_fraction, left=largest_fraction, height=0.32, color=pale_gray, edgecolor="white", linewidth=lw, hatch="//")
            ax_d.text(largest_fraction / 2, y_pos, f"{100 * largest_fraction:.1f}%", ha="center", va="center", fontsize=8.0, fontweight="bold", color=ink)
            ax_d.text(largest_fraction + (1 - largest_fraction) / 2, y_pos, f"{remaining_count}/{total}\n{100 * (1 - largest_fraction):.1f}%", ha="center", va="center", fontsize=8.0, color=ink, linespacing=1.0)
        ax_d.text(0.0, 2.45, "Share of the endpoint query–PMID graph", ha="left", va="center", fontsize=8.0, color=muted)
        ax_d.text(0.0, 2.12, "Largest component", ha="left", va="center", fontsize=8.0, fontweight="bold", color=ink)
        ax_d.text(1.0, 2.12, "Remaining", ha="right", va="center", fontsize=8.0, color=muted)
        ax_d.set_xlim(0, 1)
        ax_d.set_ylim(-0.38, 2.58)
        ax_d.set_yticks(y_positions, [f"{largest_component_queries}/{endpoint_queries}\nqueries", f"{largest_component_pmids}/{endpoint_pmids}\nendpoint PMIDs"])
        ax_d.set_xticks([])
        ax_d.tick_params(axis="y", length=0, pad=4, colors=ink)
        for spine in ax_d.spines.values():
            spine.set_visible(False)
        ax_d.text(0, 0.46, f"{int(component_summary['component_count']):,} query–PMID components | {int(endpoint_cohort['evidence_row_count']):,} evidence rows", ha="left", va="center", fontsize=8.0, color=ink)
        ax_d.text(0, 0.15, f"{endpoint_pmids:,} endpoint PMIDs | {int(historical_cohort['unique_source_document_count']):,} historical PMIDs", ha="left", va="center", fontsize=8.0, color=ink)
        ax_d.text(0, -0.16, f"PMID intersection = {int(overlap['shared_source_document_count'])}", ha="left", va="center", fontsize=8.0, fontweight="bold", color=ink)

        output = figures / "Figure_1_source_aware_temporal_endpoint.png"
        save_rgb_png(fig, output)

    upstream = [
        {"path": rel(rebuild_path), "sha256": sha256_file(rebuild_path), "role": "current_row_level_rebuild"},
        {"path": rel(future_path), "sha256": sha256_file(future_path), "role": "frozen_endpoint_flow"},
        {"path": rel(source_path), "sha256": sha256_file(source_path), "role": "current_strict_ab_provenance_audit"},
        {"path": rel(score_path), "sha256": sha256_file(score_path), "role": "fixed_retrieval_task"},
        {"path": rel(date_status_path), "sha256": sha256_file(date_status_path), "role": "date_interval_classification"},
        {"path": rel(date_history_path), "sha256": sha256_file(date_history_path), "role": "day_only_and_interval_certain_history_scenarios"},
        {"path": rel(legacy_c37_summary_path), "sha256": sha256_file(legacy_c37_summary_path), "role": "excluded_legacy_broader_lineage_not_displayed"},
        {"path": rel(legacy_c37_config_path), "sha256": sha256_file(legacy_c37_config_path), "role": "excluded_legacy_846_pair_lineage_definition"},
    ]
    return output, source_data_path, upstream


def figure2(project: Path, figures: Path) -> tuple[Path, Path, list[dict[str, str]]]:
    """Compare date policies and expose exact-tie operability at Recall@50."""

    date_recall_path = project / "revision_date_policy_v1_20260729" / "outputs" / "recall_at_50.tsv"
    tie_metrics_path = project / "revision_tie_aware_v1_20260729" / "outputs" / "tie_aware_metrics.tsv"
    operability_path = project / "revision_tie_aware_v1_20260729" / "outputs" / "three_mer_operability.tsv"
    date_rows = read_tsv(date_recall_path)
    tie_rows = read_tsv(tie_metrics_path)
    operability_rows = read_tsv(operability_path)

    date_lookup = {
        (row["scenario"], row["provenance_scope"], row["baseline"]): row
        for row in date_rows
        if row["scenario"] in {"day_only_conservative", "interval_certain_pre_cutoff"}
    }
    tie_lookup = {
        (row["scope"], row["baseline"]): row
        for row in tie_rows
        if row["query_subset"] == "all_queries" and row["k"] == "50"
    }
    operability_lookup = {row["scope"]: row for row in operability_rows}
    primary_scopes = ("temporal_strict_ab", "scaffold_cold_strict_ab")
    joint_scopes = (
        "project_defined_joint_scaffold_homology_cold_0_30",
        "project_defined_joint_scaffold_homology_cold_0_50_0_70_identical_mask",
    )

    for scope in primary_scopes:
        for baseline in BASELINES:
            day = date_lookup[("day_only_conservative", scope, baseline)]
            interval = date_lookup[("interval_certain_pre_cutoff", scope, baseline)]
            tie = tie_lookup[(scope, baseline)]
            if not np.isclose(float(day["Recall@50"]), float(tie["legacy_salted_recall"]), atol=1e-15):
                raise ValueError(f"Figure 2 day-only/tie-aware mismatch for {scope}/{baseline}")
            if int(day["evaluable_query_count"]) != int(interval["evaluable_query_count"]) or int(day["evaluable_query_count"]) != int(tie["query_count"]):
                raise ValueError(f"Figure 2 denominator mismatch for {scope}/{baseline}")
    for scope in joint_scopes:
        for baseline in BASELINES:
            tie = tie_lookup[(scope, baseline)]
            if float(tie["legacy_salted_recall"]) != 0:
                raise ValueError(f"Figure 2 joint-scope salted Recall@50 changed for {scope}/{baseline}")
            if baseline == "sequence_3mer_transfer":
                if float(tie["tie_worst_recall"]) != 0 or float(tie["tie_best_recall"]) != 1:
                    raise ValueError(f"Figure 2 non-operational 3-mer bounds changed for {scope}")
            elif float(tie["tie_expected_fractional_recall"]) != 0 or int(tie["recall_tie_dependent_query_count"]) != 0:
                raise ValueError(f"Figure 2 score-identifiable zero contract changed for {scope}/{baseline}")

    def rel(path: Path) -> str:
        return path.relative_to(project).as_posix()

    source_rows: list[dict[str, Any]] = []
    for panel, scope in zip(("A", "B"), primary_scopes):
        for baseline in BASELINES:
            for scenario in ("day_only_conservative", "interval_certain_pre_cutoff"):
                row = date_lookup[(scenario, scope, baseline)]
                source_rows.append(
                    {
                        "panel": panel,
                        "display_scope": SCOPE_LABELS[scope],
                        "provenance_scope": scope,
                        "baseline": baseline,
                        "scenario_or_subset": scenario,
                        "metric": "Recall@50",
                        "estimate": row["Recall@50"],
                        "ci95_low": "",
                        "ci95_high": "",
                        "tie_worst": "",
                        "tie_best": "",
                        "query_count": row["evaluable_query_count"],
                        "relevant_relation_count": "",
                        "status": row["analysis_label"],
                        "upstream_artifact": rel(date_recall_path),
                    }
                )
            tie = tie_lookup[(scope, baseline)]
            source_rows.append(
                {
                    "panel": panel,
                    "display_scope": SCOPE_LABELS[scope],
                    "provenance_scope": scope,
                    "baseline": baseline,
                    "scenario_or_subset": "day_only_exact_tie",
                    "metric": "Recall@50",
                    "estimate": tie["tie_expected_fractional_recall"],
                    "ci95_low": "",
                    "ci95_high": "",
                    "tie_worst": tie["tie_worst_recall"],
                    "tie_best": tie["tie_best_recall"],
                    "query_count": tie["query_count"],
                    "relevant_relation_count": tie["relevant_relation_count"],
                    "status": tie["tie_interpretation"],
                    "upstream_artifact": rel(tie_metrics_path),
                }
            )
    for scope in joint_scopes:
        for baseline in BASELINES:
            tie = tie_lookup[(scope, baseline)]
            status = "non_operational_uniform_tie_0_to_1" if baseline == "sequence_3mer_transfer" else "score_identifiable_zero"
            source_rows.append(
                {
                    "panel": "C",
                    "display_scope": SCOPE_LABELS[scope],
                    "provenance_scope": scope,
                    "baseline": baseline,
                    "scenario_or_subset": "all_queries",
                    "metric": "Recall@50",
                    "estimate": tie["tie_expected_fractional_recall"],
                    "ci95_low": "",
                    "ci95_high": "",
                    "tie_worst": tie["tie_worst_recall"],
                    "tie_best": tie["tie_best_recall"],
                    "query_count": tie["query_count"],
                    "relevant_relation_count": tie["relevant_relation_count"],
                    "status": status,
                    "upstream_artifact": rel(tie_metrics_path),
                }
            )
    for scope in DISPLAY_SCOPES:
        row = operability_lookup[scope]
        source_rows.append(
            {
                "panel": "D",
                "display_scope": SCOPE_LABELS[scope],
                "provenance_scope": scope,
                "baseline": "sequence_3mer_transfer",
                "scenario_or_subset": "operability_partition",
                "metric": "query_count",
                "estimate": row["score_operational_query_count"],
                "ci95_low": "",
                "ci95_high": "",
                "tie_worst": "",
                "tie_best": row["structural_all_zero_query_count"],
                "query_count": row["all_query_count"],
                "relevant_relation_count": row["all_relevant_relation_count"],
                "status": row["structural_all_zero_interpretation"],
                "upstream_artifact": rel(operability_path),
            }
        )

    source_path = figures / "Figure_2_source_data.tsv"
    source_fields = [
        "panel",
        "display_scope",
        "provenance_scope",
        "baseline",
        "scenario_or_subset",
        "metric",
        "estimate",
        "ci95_low",
        "ci95_high",
        "tie_worst",
        "tie_best",
        "query_count",
        "relevant_relation_count",
        "status",
        "upstream_artifact",
    ]
    write_tsv(source_path, source_fields, source_rows)

    set_style()
    fig = plt.figure(figsize=(11.8, 8.2), facecolor="white")
    grid = fig.add_gridspec(2, 6, left=0.075, right=0.985, top=0.90, bottom=0.105, hspace=0.42, wspace=0.90, height_ratios=[1.0, 1.05])
    x_limit = 0.43

    for column, scope in enumerate(primary_scopes):
        ax = fig.add_subplot(grid[0, :3] if column == 0 else grid[0, 3:])
        for yi, baseline in enumerate(BASELINES):
            day = float(date_lookup[("day_only_conservative", scope, baseline)]["Recall@50"])
            interval = float(date_lookup[("interval_certain_pre_cutoff", scope, baseline)]["Recall@50"])
            tie = tie_lookup[(scope, baseline)]
            expected = float(tie["tie_expected_fractional_recall"])
            worst = float(tie["tie_worst_recall"])
            best = float(tie["tie_best_recall"])
            shown_best = min(best, x_limit - 0.006)
            ax.hlines(yi, worst, shown_best, color=FIGURE_PALETTE["purple"], lw=2.2, alpha=0.72, zorder=1)
            ax.plot([worst, shown_best], [yi, yi], "|", color=FIGURE_PALETTE["purple"], ms=8, zorder=2)
            if best > x_limit:
                ax.plot(x_limit - 0.006, yi, marker=">", color=FIGURE_PALETTE["purple"], ms=5, clip_on=False, zorder=3)
            ax.scatter(day, yi - 0.13, marker="o", s=34, color=COLORS[baseline], edgecolor="white", linewidth=0.6, zorder=4)
            ax.scatter(interval, yi + 0.13, marker="D", s=32, facecolor="white", edgecolor=COLORS[baseline], linewidth=1.2, zorder=4)
            ax.scatter(expected, yi, marker="X", s=38, color=FIGURE_PALETTE["ink"], edgecolor="white", linewidth=0.5, zorder=5)
            ax.text(min(day + 0.008, x_limit - 0.04), yi - 0.13, f"{day:.3f}", va="center", fontsize=8.0, color=FIGURE_PALETTE["slate"])
            ax.text(min(interval + 0.008, x_limit - 0.04), yi + 0.13, f"{interval:.3f}", va="center", fontsize=8.0, color=FIGURE_PALETTE["slate"])
            if best - worst > 1e-12:
                ax.text(x_limit - 0.004, yi + 0.31, f"tie E={expected:.3f} [{worst:.3f}–{best:.3f}]", ha="right", va="center", fontsize=8.0, color=FIGURE_PALETTE["purple"])
        ax.set_yticks(np.arange(len(BASELINES)), [BASELINE_LABELS[item] for item in BASELINES])
        ax.set_ylim(len(BASELINES) - 0.50, -0.50)
        ax.set_xlim(0, x_limit)
        ax.set_xlabel("Macro Recall@50")
        count = "358 relations / 222 queries" if scope == "temporal_strict_ab" else "123 relations / 88 queries"
        ax.set_title(f"{SCOPE_LABELS[scope]}\n{count}", loc="left", x=0.075, fontweight="bold")
        ax.grid(axis="x", color=FIGURE_PALETTE["grid"], lw=0.7)
        panel_label(ax, "A" if column == 0 else "B")

    ax = fig.add_subplot(grid[1, :4])
    for row_index, scope in enumerate(joint_scopes):
        for column_index, baseline in enumerate(BASELINES):
            tie = tie_lookup[(scope, baseline)]
            non_operational = baseline == "sequence_3mer_transfer"
            face = FIGURE_PALETTE["pale_gray"] if non_operational else FIGURE_PALETTE["pale_blue"]
            edge = FIGURE_PALETTE["purple"] if non_operational else FIGURE_PALETTE["teal"]
            ax.add_patch(FancyBboxPatch((column_index - 0.46, row_index - 0.39), 0.92, 0.78, boxstyle="round,pad=0.010,rounding_size=0.04", facecolor=face, edgecolor=edge, linewidth=1.0))
            if non_operational:
                text_value = f"E={float(tie['tie_expected_fractional_recall']):.3f}\nbounds 0–1\nnon-operational"
                text_color = FIGURE_PALETTE["purple"]
            else:
                text_value = "R@50 = 0\nscore-identifiable"
                text_color = FIGURE_PALETTE["teal"]
            ax.text(column_index, row_index, text_value, ha="center", va="center", fontsize=8.0, color=text_color, fontweight="bold")
    ax.set_xlim(-0.52, len(BASELINES) - 0.48)
    ax.set_ylim(len(joint_scopes) - 0.52, -0.52)
    ax.set_xticks(np.arange(len(BASELINES)), [BASELINE_SHORT[item] for item in BASELINES])
    joint_ylabels = ["0.30\n24 rel / 19 q", "0.50/0.70 identical\n29 rel / 22 q"]
    ax.set_yticks(np.arange(len(joint_scopes)), joint_ylabels)
    ax.tick_params(length=0)
    ax.set_title("Project-defined joint scaffold–homology scopes\nlegacy salted Recall@50 = 0 in every cell", loc="left", x=0.075, fontweight="bold")
    for spine in ax.spines.values():
        spine.set_visible(False)
    panel_label(ax, "C")

    ax = fig.add_subplot(grid[1, 4:])
    operability_labels = ["Temporal", "Scaffold", "Joint 0.30", "Joint 0.50/0.70"]
    operational = np.array([int(operability_lookup[scope]["score_operational_query_count"]) for scope in DISPLAY_SCOPES])
    all_zero = np.array([int(operability_lookup[scope]["structural_all_zero_query_count"]) for scope in DISPLAY_SCOPES])
    y = np.arange(len(DISPLAY_SCOPES))
    ax.barh(
        y,
        operational,
        color=FIGURE_PALETTE["teal"],
        edgecolor=FIGURE_PALETTE["ink"],
        linewidth=0.45,
        hatch="//",
        label="Score-operational",
    )
    ax.barh(
        y,
        all_zero,
        left=operational,
        color=FIGURE_PALETTE["purple"],
        edgecolor=FIGURE_PALETTE["ink"],
        linewidth=0.45,
        hatch="..",
        label="All-zero / non-operational",
    )
    for yi, (operational_count, all_zero_count) in enumerate(zip(operational, all_zero)):
        ax.text(operational_count + all_zero_count + 3, yi, f"{operational_count} / {all_zero_count}", va="center", fontsize=8.0)
    ax.set_yticks(y, operability_labels)
    ax.invert_yaxis()
    ax.set_xlim(0, 245)
    ax.set_xlabel("Queries (operational / all-zero)")
    ax.set_title("Sequence 3-mer operability", loc="left", x=0.075, fontweight="bold")
    ax.grid(axis="x", color=FIGURE_PALETTE["grid"], lw=0.7)
    ax.legend(frameon=False, fontsize=8.0, loc="lower right")
    panel_label(ax, "D")

    legend_handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=FIGURE_PALETTE["slate"], markeredgecolor="white", markersize=6, label="Day-only salted"),
        Line2D([0], [0], marker="D", color="none", markerfacecolor="white", markeredgecolor=FIGURE_PALETTE["slate"], markersize=5.5, label="Interval-certain salted"),
        Line2D([0], [0], marker="X", color="none", markerfacecolor=FIGURE_PALETTE["ink"], markeredgecolor="white", markersize=6, label="Day-only exact-tie expectation"),
        Line2D([0], [0], color=FIGURE_PALETTE["purple"], lw=2.2, marker="|", markersize=8, label="Exact-tie best–worst range"),
    ]
    fig.legend(handles=legend_handles, frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.53, 0.975), fontsize=8.0)
    fig.text(0.53, 0.028, "Date-policy and exact-tie analyses are author-run, outcome-visible, post hoc descriptive sensitivities; tie bounds are conditional ranking uncertainty, not sampling intervals.", ha="center", va="bottom", fontsize=8.0, color=FIGURE_PALETTE["slate"])

    output = figures / "Figure_2_corrected_retrieval_performance.png"
    save_rgb_png(fig, output)
    upstream = [
        {"path": rel(date_recall_path), "sha256": sha256_file(date_recall_path), "role": "day_only_and_interval_certain_recall50"},
        {"path": rel(tie_metrics_path), "sha256": sha256_file(tie_metrics_path), "role": "exact_tie_expectations_bounds_and_status"},
        {"path": rel(operability_path), "sha256": sha256_file(operability_path), "role": "sequence_3mer_operability_partition"},
    ]
    return output, source_path, upstream


def figure3(project: Path, tables: Path, figures: Path) -> tuple[Path, Path, list[dict[str, str]]]:
    """Display tie-aware intervals and finite-sample zero-hit bounds."""

    tie_metrics_path = project / "revision_tie_aware_v1_20260729" / "outputs" / "tie_aware_metrics.tsv"
    upper_bound_path = project / "revision_tie_aware_v1_20260729" / "outputs" / "double_cold_query_hit_upper_bounds.tsv"
    tie_rows = read_tsv(tie_metrics_path)
    upper_rows = read_tsv(upper_bound_path)
    tie_lookup = {
        (row["scope"], row["baseline"]): row
        for row in tie_rows
        if row["query_subset"] == "all_queries" and row["k"] == "50"
    }
    upper_lookup = {(row["scope"], row["baseline"]): row for row in upper_rows}
    primary_scopes = ("temporal_strict_ab", "scaffold_cold_strict_ab")
    joint_scopes = (
        "project_defined_joint_scaffold_homology_cold_0_30",
        "project_defined_joint_scaffold_homology_cold_0_50_0_70_identical_mask",
    )

    for scope in primary_scopes:
        for baseline in BASELINES:
            row = tie_lookup[(scope, baseline)]
            if row["query_bootstrap_status"] != "estimable_descriptive_query_bootstrap" or row["pmid_component_bootstrap_status"] != "estimable_descriptive_pmid_component_sensitivity":
                raise ValueError(f"Figure 3 interval status changed for {scope}/{baseline}")
    for scope in joint_scopes:
        bounds = {float(upper_lookup[(scope, baseline)]["one_sided_clopper_pearson_upper"]) for baseline in BASELINES}
        if len(bounds) != 1 or any(int(upper_lookup[(scope, baseline)]["legacy_salted_query_hit_count_at_50"]) != 0 for baseline in BASELINES):
            raise ValueError(f"Figure 3 joint-scope upper-bound contract changed for {scope}")
        for baseline in BASELINES:
            tie = tie_lookup[(scope, baseline)]
            worst = float(tie["tie_worst_recall"])
            best = float(tie["tie_best_recall"])
            if baseline == "sequence_3mer_transfer":
                if not np.isclose(worst, 0.0, atol=1e-15) or not np.isclose(best, 1.0, atol=1e-15):
                    raise ValueError(f"Figure 3 3-mer operability contract changed for {scope}")
            elif not np.isclose(worst, 0.0, atol=1e-15) or not np.isclose(best, 0.0, atol=1e-15):
                raise ValueError(f"Figure 3 score-identifiable zero contract changed for {scope}/{baseline}")

    def rel(path: Path) -> str:
        return path.relative_to(project).as_posix()

    selected_rows: list[dict[str, Any]] = []
    for panel, scope in zip(("A", "B"), primary_scopes):
        for baseline in BASELINES:
            row = tie_lookup[(scope, baseline)]
            selected_rows.append(
                {
                    "panel": panel,
                    "record_type": "tie_expected_with_two_interval_sensitivities",
                    "display_scope": SCOPE_LABELS[scope],
                    "provenance_scope": scope,
                    "baseline": baseline,
                    "metric_or_item": "Recall@50",
                    "point_estimate": row["tie_expected_fractional_recall"],
                    "query_ci95_low": row["query_bootstrap_expected_recall_ci95_low"],
                    "query_ci95_high": row["query_bootstrap_expected_recall_ci95_high"],
                    "component_ci95_low": row["pmid_component_expected_recall_ci95_low"],
                    "component_ci95_high": row["pmid_component_expected_recall_ci95_high"],
                    "tie_worst": row["tie_worst_recall"],
                    "tie_best": row["tie_best_recall"],
                    "status": row["tie_interpretation"],
                    "upstream_artifact": rel(tie_metrics_path),
                }
            )
    for scope in joint_scopes:
        for baseline in BASELINES:
            upper = upper_lookup[(scope, baseline)]
            tie = tie_lookup[(scope, baseline)]
            tie_status = "non_operational_uniform_tie_0_to_1" if baseline == "sequence_3mer_transfer" else "score_identifiable_zero"
            selected_rows.append(
                {
                    "panel": "C",
                    "record_type": "one_sided_query_hit_upper_bound",
                    "display_scope": SCOPE_LABELS[scope],
                    "provenance_scope": scope,
                    "baseline": baseline,
                    "metric_or_item": "query_any_hit_at_50",
                    "point_estimate": upper["one_sided_clopper_pearson_upper"],
                    "query_ci95_low": "",
                    "query_ci95_high": "",
                    "component_ci95_low": "",
                    "component_ci95_high": "",
                    "tie_worst": tie["tie_worst_recall"],
                    "tie_best": tie["tie_best_recall"],
                    "status": tie_status,
                    "upstream_artifact": rel(upper_bound_path),
                }
            )

    source_path = figures / "Figure_3_source_data.tsv"
    fields = [
        "panel",
        "record_type",
        "display_scope",
        "provenance_scope",
        "baseline",
        "metric_or_item",
        "point_estimate",
        "query_ci95_low",
        "query_ci95_high",
        "component_ci95_low",
        "component_ci95_high",
        "tie_worst",
        "tie_best",
        "status",
        "upstream_artifact",
    ]
    write_tsv(source_path, fields, selected_rows)

    set_style()
    fig = plt.figure(figsize=(7.5, 5.2), facecolor="white")
    outer_grid = fig.add_gridspec(2, 1, left=0.13, right=0.96, top=0.865, bottom=0.12, hspace=0.76, height_ratios=[1.0, 0.65])
    top_grid = outer_grid[0, 0].subgridspec(1, 2, wspace=0.40)
    bottom_grid = outer_grid[1, 0].subgridspec(1, 2, width_ratios=[0.64, 0.36], wspace=0.20)

    for column, scope in enumerate(primary_scopes):
        ax = fig.add_subplot(top_grid[0, column])
        y = np.arange(len(BASELINES))
        for yi, baseline in enumerate(BASELINES):
            row = tie_lookup[(scope, baseline)]
            point = float(row["tie_expected_fractional_recall"])
            query_low = float(row["query_bootstrap_expected_recall_ci95_low"])
            query_high = float(row["query_bootstrap_expected_recall_ci95_high"])
            component_low = float(row["pmid_component_expected_recall_ci95_low"])
            component_high = float(row["pmid_component_expected_recall_ci95_high"])
            query_error = np.array([[max(point - query_low, 0.0)], [max(query_high - point, 0.0)]])
            component_error = np.array([[max(point - component_low, 0.0)], [max(component_high - point, 0.0)]])
            ax.errorbar(point, yi - 0.11, xerr=query_error, fmt="o", ms=4.2, capsize=2.5, lw=1.2, color=FIGURE_PALETTE["blue"], markerfacecolor="white", zorder=3)
            ax.errorbar(point, yi + 0.11, xerr=component_error, fmt="s", ms=4.0, capsize=2.5, lw=1.2, color=FIGURE_PALETTE["amber"], markerfacecolor="white", zorder=3)
            ax.scatter(point, yi, marker="D", s=28, color=FIGURE_PALETTE["ink"], edgecolor="white", linewidth=0.5, zorder=4)
            label_anchor = max(query_high, component_high)
            ax.annotate(
                f"E={point:.3f}",
                xy=(label_anchor, yi),
                xytext=(5, 0),
                textcoords="offset points",
                ha="left",
                va="center",
                fontsize=8.0,
                color=FIGURE_PALETTE["ink"],
            )
        ax.set_yticks(y, [BASELINE_LABELS[item] for item in BASELINES])
        ax.set_ylim(len(BASELINES) - 0.48, -0.48)
        ax.set_xlim(0, 0.50)
        ax.set_xlabel("Tie-expected macro Recall@50")
        first = tie_lookup[(scope, BASELINES[0])]
        ax.set_title(
            f"{SCOPE_LABELS[scope]}\n{first['pmid_source_document_count']} PMIDs · {first['pmid_component_count']} components · {first['query_count']} queries",
            loc="left",
            x=0.075,
            fontweight="bold",
        )
        ax.grid(axis="x", color=FIGURE_PALETTE["grid"], lw=0.7)
        panel_label(ax, "A" if column == 0 else "B")

    ax_bounds = fig.add_subplot(bottom_grid[0, 0])
    y_positions = np.array([0.0, 1.0])
    upper_values: list[float] = []
    query_counts: list[int] = []
    for scope in joint_scopes:
        row = upper_lookup[(scope, BASELINES[0])]
        upper_values.append(float(row["one_sided_clopper_pearson_upper"]))
        query_counts.append(int(row["query_count"]))
    bars = ax_bounds.barh(y_positions, upper_values, height=0.36, color=(FIGURE_PALETTE["blue"], FIGURE_PALETTE["teal"]), alpha=0.92)
    for yi, (bar, upper_value, query_count) in enumerate(zip(bars, upper_values, query_counts)):
        ax_bounds.text(upper_value + 0.004, yi, f"0/{query_count}   U95={upper_value:.3f}", va="center", fontsize=8.0, fontweight="bold")
    ax_bounds.set_yticks(y_positions, ["Joint 0.30", "Joint 0.50/0.70\nidentical mask"])
    ax_bounds.set_ylim(1.48, -0.48)
    ax_bounds.set_xlim(0, 0.18)
    ax_bounds.set_xlabel("One-sided exact 95% upper bound")
    ax_bounds.set_title("Zero-hit query-level upper bounds", loc="left", x=0.075, fontweight="bold")
    ax_bounds.grid(axis="x", color=FIGURE_PALETTE["grid"], lw=0.7)
    panel_label(ax_bounds, "C")

    ax_status = fig.add_subplot(bottom_grid[0, 1])
    ax_status.set_xlim(0, 1)
    ax_status.set_ylim(0, 1)
    ax_status.axis("off")
    ax_status.plot([0.02, 0.02], [0.02, 0.98], color=FIGURE_PALETTE["grid"], lw=0.8, transform=ax_status.transAxes, clip_on=False)
    ax_status.text(0.10, 0.72, "Score-identifiable:", fontsize=8.2, fontweight="bold", va="top")
    ax_status.text(0.10, 0.57, "Popularity, Morgan, Pair\n0 observed top-50 hits", fontsize=8.0, va="top", linespacing=1.15)
    ax_status.text(0.10, 0.29, "Sequence 3-mer:", fontsize=8.2, fontweight="bold", va="top")
    ax_status.text(0.10, 0.14, "non-operational under exact ties", fontsize=8.0, color=FIGURE_PALETTE["slate"], va="top")

    legend_handles = [
        Line2D([0], [0], marker="D", color="none", markerfacecolor=FIGURE_PALETTE["ink"], markeredgecolor="white", markersize=5.5, label="Exact-tie expectation"),
        Line2D([0], [0], marker="o", color=FIGURE_PALETTE["blue"], lw=1.4, markerfacecolor="white", markersize=5, label="Query-bootstrap 95% interval"),
        Line2D([0], [0], marker="s", color=FIGURE_PALETTE["amber"], lw=1.4, markerfacecolor="white", markersize=5, label="PMID-component 95% sensitivity interval"),
    ]
    fig.legend(handles=legend_handles, frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.545, 0.982), fontsize=7.7)

    output = figures / "Figure_3_dependence_and_fidelity_audits.png"
    save_rgb_png_fixed_canvas(fig, output)
    upstream = [
        {"path": rel(tie_metrics_path), "sha256": sha256_file(tie_metrics_path), "role": "tie_expectation_query_bootstrap_and_pmid_component_intervals"},
        {"path": rel(upper_bound_path), "sha256": sha256_file(upper_bound_path), "role": "one_sided_exact_query_hit_upper_bounds"},
    ]
    return output, source_path, upstream


def figure4(project: Path, tables: Path, figures: Path) -> tuple[Path, Path, list[dict[str, str]]]:
    """Display policy robustness independently of the interval diagnostics."""

    date_policy_path = tables / "Table_S8_date_precision_policy.tsv"
    policy_path = tables / "Table_S9_weight_and_structure_policy.tsv"
    top100_path = tables / "Table_S2_top100_exhaustive_fidelity.tsv"
    date_policy_rows = read_tsv(date_policy_path)
    policy_rows = read_tsv(policy_path)
    top100_rows = read_tsv(top100_path)
    difference_lookup = {
        (row["scope"], row["metric_or_item"]): row
        for row in top100_rows
        if row["record_type"] == "metric_difference"
    }
    for metric in METRICS:
        left = difference_lookup[("double_cold_0_50", metric)]
        right = difference_lookup[("double_cold_0_70", metric)]
        if left["relation_count"] != right["relation_count"] or left["query_count"] != right["query_count"] or not np.isclose(float(left["exhaustive_minus_top100"]), float(right["exhaustive_minus_top100"]), atol=1e-15):
            raise ValueError("Figure 4 0.50/0.70 top-100 masks are no longer identical")

    def rel(path: Path) -> str:
        return path.relative_to(project).as_posix()

    date_reference = "day_only_conservative"
    date_sensitivity = "interval_certain_pre_cutoff"
    weight_reference = "A1_B0_7_primary"
    weight_sensitivities = ("A1_B0_5", "A1_B1_0_all_equal")
    structure_reference = "raw_primary"
    structure_sensitivity = "cleanup_parent_charge_tautomer"
    matrix_scopes = ("temporal_strict_ab", "scaffold_cold_strict_ab")
    structure_scope_map = {
        "temporal_strict_ab": "temporal_strict_ab",
        "scaffold_cold_strict_ab": "scaffold_cold",
    }

    def require_status(row: dict[str, str], expected: str, label: str) -> None:
        if row.get("status", "") != expected:
            raise ValueError(f"{label}: expected status {expected!r}; found {row.get('status', '')!r}")

    def numeric_value(row: dict[str, str], label: str) -> float:
        try:
            return float(row["value"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"{label}: non-numeric value {row.get('value')!r}") from exc

    require_unique_row(date_policy_rows, "Figure 4 date-policy reference", section="history", scenario_or_policy=date_reference, item="selected_source_row_count")
    require_unique_row(date_policy_rows, "Figure 4 date-policy sensitivity", section="history", scenario_or_policy=date_sensitivity, item="selected_source_row_count")
    expected_weights = {"A1_B0_5": 0.5, weight_reference: 0.7, "A1_B1_0_all_equal": 1.0}
    for policy, expected_weight in expected_weights.items():
        weight_row = require_unique_row(policy_rows, "Figure 4 Tier-B policy", section="weight_variant", scenario_or_policy=policy, item="tier_B_weight")
        if not np.isclose(numeric_value(weight_row, f"Tier-B weight {policy}"), expected_weight, atol=1e-12):
            raise ValueError(f"Figure 4 Tier-B weight changed for {policy}")
    for policy in (structure_reference, structure_sensitivity):
        policy_summary = [
            row
            for row in policy_rows
            if row["section"] == "structure_structure_policy_summary"
            and row["scenario_or_policy"] == policy
            and row["item"] == "role"
        ]
        if len(policy_summary) != 2 or {row["value"] for row in policy_summary} != {"historical", "query"}:
            raise ValueError(f"Figure 4 structure policy roles missing or duplicated for {policy}")
        if any(row["status"] != "complete" for row in policy_summary):
            raise ValueError(f"Figure 4 structure policy status changed for {policy}")

    def date_value(policy: str, scope: str, baseline: str) -> float:
        row = require_unique_row(date_policy_rows, "Figure 4 date-policy Recall@50", section="recall_at_50", scenario_or_policy=policy, scope=scope, baseline=baseline, item="Recall@50")
        require_status(row, "", f"Figure 4 date-policy Recall@50 {policy}/{scope}/{baseline}")
        return numeric_value(row, f"Figure 4 date-policy Recall@50 {policy}/{scope}/{baseline}")

    def date_counts(policy: str, scope: str) -> tuple[int, int]:
        relation_row = require_unique_row(date_policy_rows, "Figure 4 date-policy relation count", section="scope", scenario_or_policy=policy, scope=scope, baseline="", item="candidate_relation_count")
        query_row = require_unique_row(date_policy_rows, "Figure 4 date-policy query count", section="scope", scenario_or_policy=policy, scope=scope, baseline="", item="query_count")
        require_status(relation_row, "", f"Figure 4 date-policy relation count {policy}/{scope}")
        require_status(query_row, "", f"Figure 4 date-policy query count {policy}/{scope}")
        return int(relation_row["value"]), int(query_row["value"])

    def weight_value(policy: str, scope: str, baseline: str) -> float:
        row = require_unique_row(policy_rows, "Figure 4 weight-policy Recall@50", section="weight_metric", scenario_or_policy=policy, scope=scope, baseline=baseline, item="Recall@50")
        require_status(row, "descriptive_fixed_salt", f"Figure 4 weight-policy Recall@50 {policy}/{scope}/{baseline}")
        return numeric_value(row, f"Figure 4 weight-policy Recall@50 {policy}/{scope}/{baseline}")

    def weight_counts(policy: str, scope: str) -> tuple[int, int]:
        counts: set[tuple[int, int]] = set()
        for baseline in BASELINES:
            relation_row = require_unique_row(policy_rows, "Figure 4 weight-policy relation count", section="weight_metric", scenario_or_policy=policy, scope=scope, baseline=baseline, item="relation_count")
            query_row = require_unique_row(policy_rows, "Figure 4 weight-policy query count", section="weight_metric", scenario_or_policy=policy, scope=scope, baseline=baseline, item="query_count")
            require_status(relation_row, "descriptive_fixed_salt", f"Figure 4 weight-policy relation count {policy}/{scope}/{baseline}")
            require_status(query_row, "descriptive_fixed_salt", f"Figure 4 weight-policy query count {policy}/{scope}/{baseline}")
            counts.add((int(relation_row["value"]), int(query_row["value"])))
        if len(counts) != 1:
            raise ValueError(f"Figure 4 weight-policy scope counts differ across baselines for {policy}/{scope}")
        return next(iter(counts))

    def structure_value(policy: str, scope: str, baseline: str) -> float:
        row = require_unique_row(policy_rows, "Figure 4 structure-policy Recall@50", section="structure_scope_recall_at_50", scenario_or_policy=policy, scope=scope, baseline=baseline, item="recall_at_50")
        require_status(row, "estimable", f"Figure 4 structure-policy Recall@50 {policy}/{scope}/{baseline}")
        return numeric_value(row, f"Figure 4 structure-policy Recall@50 {policy}/{scope}/{baseline}")

    def structure_counts(policy: str, scope: str) -> tuple[int, int]:
        counts: set[tuple[int, int]] = set()
        for baseline in BASELINES:
            relation_row = require_unique_row(policy_rows, "Figure 4 structure-policy relation count", section="structure_scope_recall_at_50", scenario_or_policy=policy, scope=scope, baseline=baseline, item="relation_count")
            query_row = require_unique_row(policy_rows, "Figure 4 structure-policy query count", section="structure_scope_recall_at_50", scenario_or_policy=policy, scope=scope, baseline=baseline, item="query_count")
            require_status(relation_row, "estimable", f"Figure 4 structure-policy relation count {policy}/{scope}/{baseline}")
            require_status(query_row, "estimable", f"Figure 4 structure-policy query count {policy}/{scope}/{baseline}")
            counts.add((int(relation_row["value"]), int(query_row["value"])))
        if len(counts) != 1:
            raise ValueError(f"Figure 4 structure-policy scope counts differ across baselines for {policy}/{scope}")
        return next(iter(counts))

    def exhaustive_row(scope: str, metric: str = "Recall@50") -> dict[str, str]:
        row = require_unique_row(top100_rows, "Figure 4 exhaustive comparison", record_type="metric_difference", scope=scope, metric_or_item=metric)
        if row["status"] != "post_hoc_sensitivity_only":
            raise ValueError(f"Figure 4 exhaustive-comparison status changed for {scope}/{metric}")
        return row

    date_joint_scopes = {"joint_scaffold_homology_0.30": 1, "joint_scaffold_homology_0.50/0.70": 2}
    for policy in (date_reference, date_sensitivity):
        for scope, expected_count in date_joint_scopes.items():
            for baseline in BASELINES:
                rows = [
                    row
                    for row in date_policy_rows
                    if row["section"] == "recall_at_50"
                    and row["scenario_or_policy"] == policy
                    and row["scope"] == scope
                    and row["baseline"] == baseline
                    and row["item"] == "Recall@50"
                ]
                if len(rows) != expected_count or any(row["status"] != "" or not np.isclose(float(row["value"]), 0.0, atol=1e-15) for row in rows):
                    raise ValueError(f"Figure 4 date-policy joint-zero contract changed for {policy}/{scope}/{baseline}")
    joint_scopes = (
        "project_defined_joint_scaffold_homology_cold_0_30",
        "project_defined_joint_scaffold_homology_cold_0_50_0_70_identical_mask",
    )
    for policy in (weight_reference, *weight_sensitivities):
        for scope in joint_scopes:
            for baseline in BASELINES:
                if not np.isclose(weight_value(policy, scope, baseline), 0.0, atol=1e-15):
                    raise ValueError(f"Figure 4 weight-policy joint-zero contract changed for {policy}/{scope}/{baseline}")
    structure_joint_scopes = ("joint_scaffold_homology_0_30", "joint_scaffold_homology_0_50_0_70_identical")
    for policy in (structure_reference, structure_sensitivity):
        for scope in structure_joint_scopes:
            for baseline in BASELINES:
                if not np.isclose(structure_value(policy, scope, baseline), 0.0, atol=1e-15):
                    raise ValueError(f"Figure 4 structure-policy joint-zero contract changed for {policy}/{scope}/{baseline}")
    exhaustive_joint_values: set[tuple[float, float]] = set()
    for scope in ("double_cold_0_30", "double_cold_0_50", "double_cold_0_70"):
        row = exhaustive_row(scope)
        joint_values = (float(row["top100_value"]), float(row["exhaustive_value"]))
        exhaustive_joint_values.add(joint_values)
        if not np.isclose(joint_values[0], 0.0, atol=1e-15) or not np.isclose(joint_values[1], 0.0, atol=1e-15):
            raise ValueError(f"Figure 4 exhaustive joint-zero contract changed for {scope}")
    if len(exhaustive_joint_values) != 1:
        raise ValueError("Figure 4 exhaustive joint scopes no longer share one zero transition")
    exhaustive_joint_primary, exhaustive_joint_sensitivity = next(iter(exhaustive_joint_values))

    matrix = np.full((5, 8), np.nan, dtype=float)
    matrix_rows: list[dict[str, Any]] = []
    matrix_scope_text = ["Unchanged", "Unchanged", "Unchanged", "", "Unchanged"]
    structure_scope_order = [
        ("T", "temporal_strict_ab"),
        ("S", "scaffold_cold"),
        ("J0.30", "joint_scaffold_homology_0_30"),
        ("J0.50/0.70", "joint_scaffold_homology_0_50_0_70_identical"),
    ]
    expected_structure_counts = {
        "temporal_strict_ab": ((358, 222), (358, 222)),
        "scaffold_cold": ((123, 88), (117, 81)),
        "joint_scaffold_homology_0_30": ((24, 19), (26, 21)),
        "joint_scaffold_homology_0_50_0_70_identical": ((29, 22), (31, 24)),
    }
    structure_summary_parts: list[str] = []
    structure_counts_by_scope: dict[str, tuple[tuple[int, int], tuple[int, int]]] = {}
    for short_label, scope in structure_scope_order:
        primary_counts = structure_counts(structure_reference, scope)
        sensitivity_counts = structure_counts(structure_sensitivity, scope)
        if (primary_counts, sensitivity_counts) != expected_structure_counts[scope]:
            raise ValueError(f"Figure 4 structure scope counts changed for {scope}: {primary_counts} -> {sensitivity_counts}")
        structure_counts_by_scope[scope] = (primary_counts, sensitivity_counts)
        primary_text = f"{primary_counts[0]}/{primary_counts[1]}"
        sensitivity_text = f"{sensitivity_counts[0]}/{sensitivity_counts[1]}"
        if primary_counts == sensitivity_counts:
            structure_summary_parts.append(f"{short_label}: {primary_text} unchanged")
        else:
            structure_summary_parts.append(f"{short_label}: {primary_text} → {sensitivity_text}")
    structure_count_summary = "; ".join(structure_summary_parts)
    matrix_scope_text[3] = structure_count_summary
    temporal_primary, temporal_sensitivity = structure_counts_by_scope["temporal_strict_ab"]
    scaffold_primary, scaffold_sensitivity = structure_counts_by_scope["scaffold_cold"]
    if temporal_primary != temporal_sensitivity:
        raise ValueError("Figure 4 temporal structure-policy scope is no longer unchanged")
    scope_exception_display = (
        f"Scaffold: {scaffold_primary[0]}/{scaffold_primary[1]} → "
        f"{scaffold_sensitivity[0]}/{scaffold_sensitivity[1]}\n"
        f"Temporal: {temporal_primary[0]}/{temporal_primary[1]} unchanged"
    )

    def add_policy_matrix_row(
        row_index: int,
        perturbation: str,
        display_row: str,
        reference_policy: str,
        sensitivity_policy: str,
        value_reader,
        count_reader,
        source_path: Path,
        source_scope_map: dict[str, str] | None = None,
    ) -> None:
        mapping = source_scope_map or {scope: scope for scope in matrix_scopes}
        for scope_index, display_scope in enumerate(matrix_scopes):
            source_scope = mapping[display_scope]
            primary_counts = count_reader(reference_policy, source_scope)
            sensitivity_counts = count_reader(sensitivity_policy, source_scope)
            if row_index != 3 and primary_counts != sensitivity_counts:
                raise ValueError(f"Figure 4 {perturbation} unexpectedly changed scope membership for {display_scope}")
            for baseline_index, baseline in enumerate(BASELINES):
                primary_value = value_reader(reference_policy, source_scope, baseline)
                sensitivity_value = value_reader(sensitivity_policy, source_scope, baseline)
                delta = sensitivity_value - primary_value
                matrix[row_index, scope_index * len(BASELINES) + baseline_index] = delta
                matrix_rows.append(
                    {
                        "panel": "A",
                        "record_type": "policy_robustness_delta_recall_at_50",
                        "display_scope": "Temporal strict A/B" if display_scope == "temporal_strict_ab" else "Scaffold-cold",
                        "provenance_scope": source_scope,
                        "baseline": baseline,
                        "metric_or_item": "Recall@50",
                        "point_estimate": delta,
                        "perturbation": perturbation,
                        "display_row": display_row,
                        "reference_policy": reference_policy,
                        "sensitivity_policy": sensitivity_policy,
                        "scope": display_scope,
                        "primary_value": primary_value,
                        "sensitivity_value": sensitivity_value,
                        "delta": delta,
                        "scope_relation_count_primary": primary_counts[0],
                        "scope_query_count_primary": primary_counts[1],
                        "scope_relation_count_sensitivity": sensitivity_counts[0],
                        "scope_query_count_sensitivity": sensitivity_counts[1],
                        "scope_change": matrix_scope_text[row_index],
                        "joint_zero_status": "zero_retained_caption_qualified",
                        "status": "descriptive_policy_sensitivity",
                        "upstream_artifact": rel(source_path),
                    }
                )

    add_policy_matrix_row(0, "date_policy", "Date: interval-certain", date_reference, date_sensitivity, date_value, date_counts, date_policy_path)
    add_policy_matrix_row(1, "tier_B_weight", "Tier B = 0.5", weight_reference, weight_sensitivities[0], weight_value, weight_counts, policy_path)
    add_policy_matrix_row(2, "tier_B_weight", "Tier B = 1.0", weight_reference, weight_sensitivities[1], weight_value, weight_counts, policy_path)
    add_policy_matrix_row(3, "structure_representation", "Structure policy", structure_reference, structure_sensitivity, structure_value, structure_counts, policy_path, structure_scope_map)

    for scope_index, scope in enumerate(matrix_scopes):
        row = exhaustive_row(scope)
        primary_value = float(row["top100_value"])
        sensitivity_value = float(row["exhaustive_value"])
        delta = float(row["exhaustive_minus_top100"])
        if not np.isclose(sensitivity_value - primary_value, delta, atol=1e-15):
            raise ValueError(f"Figure 4 exhaustive delta arithmetic changed for {scope}")
        counts = (int(row["relation_count"]), int(row["query_count"]))
        for baseline_index, baseline in enumerate(BASELINES):
            evaluated = baseline == "structure_sequence_pair_neighbor"
            if evaluated:
                matrix[4, scope_index * len(BASELINES) + baseline_index] = delta
            matrix_rows.append(
                {
                    "panel": "A",
                    "record_type": "policy_robustness_delta_recall_at_50" if evaluated else "not_evaluated",
                    "display_scope": "Temporal strict A/B" if scope == "temporal_strict_ab" else "Scaffold-cold",
                    "provenance_scope": scope,
                    "baseline": baseline,
                    "metric_or_item": "Recall@50",
                    "point_estimate": delta if evaluated else "",
                    "perturbation": "exhaustive_ranking",
                    "display_row": "Exhaustive ranking",
                    "reference_policy": "frozen_top_100",
                    "sensitivity_policy": "exhaustive_all_1131_historical_targets",
                    "scope": scope,
                    "primary_value": primary_value if evaluated else "",
                    "sensitivity_value": sensitivity_value if evaluated else "",
                    "delta": delta if evaluated else "",
                    "scope_relation_count_primary": counts[0],
                    "scope_query_count_primary": counts[1],
                    "scope_relation_count_sensitivity": counts[0],
                    "scope_query_count_sensitivity": counts[1],
                    "scope_change": "Unchanged",
                    "joint_zero_status": "pair_zero_retained" if evaluated else "not_evaluated",
                    "status": "post_hoc_sensitivity_only" if evaluated else "not_evaluated",
                    "upstream_artifact": rel(top100_path),
                }
            )

    expected_matrix = np.array(
        [
            [0.000, 0.000, 0.009, 0.003, 0.000, 0.000, 0.021, 0.009],
            [0.002, 0.000, -0.021, -0.051, 0.031, 0.000, -0.055, -0.055],
            [0.067, -0.007, 0.032, 0.025, 0.160, 0.000, 0.076, 0.087],
            [0.000, 0.000, -0.006, -0.010, 0.015, 0.001, -0.037, -0.025],
            [np.nan, np.nan, np.nan, 0.0015015015015014677, np.nan, np.nan, np.nan, 0.000],
        ],
        dtype=float,
    )
    for row_index in range(matrix.shape[0]):
        for column_index in range(matrix.shape[1]):
            expected = expected_matrix[row_index, column_index]
            actual = matrix[row_index, column_index]
            if np.isnan(expected):
                if not np.isnan(actual):
                    raise ValueError(f"Figure 4 not-evaluated cell became numeric at {row_index}/{column_index}")
            elif row_index == 4:
                if not np.isclose(actual, expected, atol=1e-12):
                    raise ValueError(f"Figure 4 exhaustive expected delta changed at {row_index}/{column_index}: {actual}")
            elif not np.isclose(actual, expected, atol=5.1e-4):
                raise ValueError(f"Figure 4 policy expected delta changed at {row_index}/{column_index}: {actual}")

    all_policy_upstream = ";".join((rel(date_policy_path), rel(policy_path), rel(top100_path)))
    summary_rows = [
        {
            "panel": "B",
            "record_type": "scope_membership_exception",
            "display_row": "Only structure policy changed scope",
            "figure_display_text": scope_exception_display.replace("\n", "; "),
            "scope_change": structure_count_summary,
            "status": "descriptive_policy_sensitivity",
            "upstream_artifact": rel(policy_path),
        },
        {
            "panel": "B",
            "record_type": "scope_membership_general_result",
            "display_row": "All other policies",
            "figure_display_text": "All other policies: unchanged",
            "scope_change": "Unchanged",
            "status": "descriptive_policy_sensitivity",
            "upstream_artifact": all_policy_upstream,
        },
        {
            "panel": "B",
            "record_type": "joint_cold_general_result",
            "display_row": "All policies",
            "figure_display_text": "Zero-hit conclusion retained under all policies†",
            "joint_zero_status": "zero_retained_caption_qualified",
            "status": "descriptive_policy_sensitivity",
            "upstream_artifact": all_policy_upstream,
        },
        {
            "panel": "B",
            "record_type": "joint_cold_exhaustive_pair_result",
            "display_row": "Exhaustive Pair",
            "figure_display_text": f"Exhaustive Pair: {exhaustive_joint_primary:g} → {exhaustive_joint_sensitivity:g}",
            "joint_zero_status": "pair_zero_retained",
            "status": "post_hoc_sensitivity_only",
            "upstream_artifact": rel(top100_path),
        },
    ]
    source_rows = matrix_rows + summary_rows
    source_path = figures / "Figure_4_source_data.tsv"
    fields = [
        "panel",
        "record_type",
        "display_scope",
        "provenance_scope",
        "baseline",
        "metric_or_item",
        "point_estimate",
        "perturbation",
        "display_row",
        "figure_display_text",
        "reference_policy",
        "sensitivity_policy",
        "scope",
        "primary_value",
        "sensitivity_value",
        "delta",
        "scope_relation_count_primary",
        "scope_query_count_primary",
        "scope_relation_count_sensitivity",
        "scope_query_count_sensitivity",
        "scope_change",
        "joint_zero_status",
        "status",
        "upstream_artifact",
    ]
    write_tsv(source_path, fields, source_rows)

    set_style()
    fig = plt.figure(figsize=(7.3, 4.5), facecolor="white")
    outer_grid = fig.add_gridspec(2, 2, left=0.055, right=0.985, top=0.81, bottom=0.15, width_ratios=[0.73, 0.27], height_ratios=[1.0, 0.10], wspace=0.11, hspace=0.28)
    heat_grid = outer_grid[0, 0].subgridspec(1, 2, width_ratios=[2.25, 8.0], wspace=0.02)
    ax_rows = fig.add_subplot(heat_grid[0, 0])
    ax_heat = fig.add_subplot(heat_grid[0, 1])
    ax_summary = fig.add_subplot(outer_grid[0, 1])
    colorbar_grid = outer_grid[1, 0].subgridspec(1, 2, width_ratios=[2.25, 8.0], wspace=0.02)
    ax_colorbar_label = fig.add_subplot(colorbar_grid[0, 0])
    cax = fig.add_subplot(colorbar_grid[0, 1])
    ax_colorbar_label.axis("off")

    row_labels = ["Date: interval-certain", "Tier B = 0.5", "Tier B = 1.0", "Structure policy", "Exhaustive ranking"]
    ax_rows.set_xlim(0, 1)
    ax_rows.set_ylim(4.5, -0.5)
    ax_rows.axis("off")
    for row_index, label in enumerate(row_labels):
        ax_rows.text(1.0, row_index, label, ha="right", va="center", fontsize=8.0, color=FIGURE_PALETTE["ink"])
    ax_rows.text(-0.28, 1.28, "a", transform=ax_rows.transAxes, ha="left", va="bottom", fontsize=10.5, fontweight="bold", color=FIGURE_PALETTE["ink"], clip_on=False)
    ax_rows.text(0.03, 1.28, "Recall@50 sensitivity", transform=ax_rows.transAxes, ha="left", va="bottom", fontsize=9.2, fontweight="bold", color=FIGURE_PALETTE["ink"], clip_on=False)

    cmap = plt.get_cmap("PuOr").copy()
    cmap.set_bad(FIGURE_PALETTE["pale_gray"])
    norm = TwoSlopeNorm(vmin=-0.16, vcenter=0.0, vmax=0.16)
    image = ax_heat.imshow(np.ma.masked_invalid(matrix), aspect="auto", cmap=cmap, norm=norm, interpolation="nearest")
    for row_index in range(matrix.shape[0]):
        for column_index in range(matrix.shape[1]):
            value = matrix[row_index, column_index]
            if np.isnan(value):
                label = "—"
                text_color = FIGURE_PALETTE["slate"]
            else:
                label = f"{value:+.3f}".replace("+0.", "+.").replace("-0.", "-.")
                red, green, blue, _ = cmap(norm(value))
                luminance = 0.2126 * red + 0.7152 * green + 0.0722 * blue
                text_color = "white" if luminance < 0.46 else FIGURE_PALETTE["ink"]
            ax_heat.text(column_index, row_index, label, ha="center", va="center", fontproperties=MATRIX_CELL_FONT, color=text_color)
    ax_heat.set_xticks(np.arange(8), ["Pop", "3-mer", "Morgan", "Pair"] * 2, fontsize=7.5)
    ax_heat.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False, length=0, pad=3)
    ax_heat.set_yticks([])
    ax_heat.tick_params(axis="y", length=0)
    ax_heat.text(0.25, 1.16, "Temporal strict A/B", transform=ax_heat.transAxes, ha="center", va="bottom", fontsize=8.0, fontweight="bold", color=FIGURE_PALETTE["ink"])
    ax_heat.text(0.75, 1.16, "Scaffold-cold", transform=ax_heat.transAxes, ha="center", va="bottom", fontsize=8.0, fontweight="bold", color=FIGURE_PALETTE["ink"])
    ax_heat.axvline(3.5, color="white", linewidth=1.5)
    ax_heat.set_xticks(np.arange(-0.5, 8, 1), minor=True)
    ax_heat.set_yticks(np.arange(-0.5, 5, 1), minor=True)
    ax_heat.grid(which="minor", color="white", linewidth=1.3)
    ax_heat.tick_params(which="minor", length=0)
    for spine in ax_heat.spines.values():
        spine.set_visible(False)

    ax_summary.set_xlim(0, 1)
    ax_summary.set_ylim(0, 1)
    ax_summary.axis("off")
    ax_summary.text(-0.20, 1.28, "b", transform=ax_summary.transAxes, ha="left", va="bottom", fontsize=10.5, fontweight="bold", color=FIGURE_PALETTE["ink"], clip_on=False)
    ax_summary.text(0.03, 1.28, "Scope and conclusion stability", transform=ax_summary.transAxes, ha="left", va="bottom", fontsize=9.2, fontweight="bold", color=FIGURE_PALETTE["ink"], clip_on=False)

    ax_summary.text(0.02, 0.96, "Scope membership", ha="left", va="top", fontsize=8.6, fontweight="bold", color=FIGURE_PALETTE["ink"])
    ax_summary.text(0.02, 0.885, "relations / queries", ha="left", va="top", fontsize=7.3, color=FIGURE_PALETTE["slate"])
    ax_summary.text(0.02, 0.78, "Only structure policy\nchanged scope*", ha="left", va="top", fontsize=8.3, fontweight="bold", color=FIGURE_PALETTE["ink"], linespacing=1.05)
    ax_summary.add_patch(Rectangle((0.02, 0.44), 0.96, 0.20, facecolor=FIGURE_PALETTE["amber"], edgecolor="none", alpha=0.18))
    ax_summary.add_patch(Rectangle((0.02, 0.44), 0.96, 0.20, fill=False, edgecolor=FIGURE_PALETTE["amber"], linewidth=0.9))
    ax_summary.text(0.07, 0.54, scope_exception_display, ha="left", va="center", fontsize=7.6, color=FIGURE_PALETTE["ink"], linespacing=1.22)
    ax_summary.text(0.02, 0.39, "All other policies: unchanged", ha="left", va="top", fontsize=7.6, color=FIGURE_PALETTE["slate"])

    ax_summary.text(0.02, 0.26, "Joint-cold conclusion", ha="left", va="top", fontsize=8.6, fontweight="bold", color=FIGURE_PALETTE["ink"])
    ax_summary.text(0.02, 0.15, "Zero-hit conclusion retained\nunder all policies†", ha="left", va="top", fontsize=8.2, fontweight="bold", color=FIGURE_PALETTE["ink"], linespacing=1.05)
    ax_summary.text(0.02, 0.005, f"Exhaustive Pair: {exhaustive_joint_primary:g} → {exhaustive_joint_sensitivity:g}", ha="left", va="bottom", fontsize=7.6, color=FIGURE_PALETTE["slate"])

    colorbar = fig.colorbar(image, cax=cax, orientation="horizontal")
    colorbar.set_label("Δ Recall@50 vs primary", fontsize=8.0, labelpad=3)
    colorbar.set_ticks([-0.16, -0.08, 0.0, 0.08, 0.16])
    colorbar.ax.tick_params(labelsize=7.5, length=2, width=0.6)
    ax_unused = fig.add_subplot(outer_grid[1, 1])
    ax_unused.axis("off")

    png_output = figures / "Figure_4_policy_robustness.png"
    save_rgb_png(fig, png_output)
    upstream = [
        {"path": rel(date_policy_path), "sha256": sha256_file(date_policy_path), "role": "date_policy_robustness"},
        {"path": rel(policy_path), "sha256": sha256_file(policy_path), "role": "tier_B_weight_and_structure_policy_robustness"},
        {"path": rel(top100_path), "sha256": sha256_file(top100_path), "role": "exhaustive_pair_neighbour_robustness"},
    ]
    return png_output, source_path, upstream


def figure5(project: Path, tables: Path, figures: Path) -> tuple[Path, Path, list[dict[str, str]]]:
    """Display representation-dependent similarity without inventing individual observations."""

    similarity_path = tables / "Table_S11_maximum_similarity_distributions.tsv"
    similarity_rows = read_tsv(similarity_path)

    def rel(path: Path) -> str:
        return path.relative_to(project).as_posix()

    scope_order = (
        "temporal_strict_ab",
        "scaffold_cold",
        "joint_scaffold_homology_0_30",
        "joint_scaffold_homology_0_50_0_70_identical",
    )
    scope_titles = {
        "temporal_strict_ab": "Temporal strict A/B",
        "scaffold_cold": "Scaffold-cold",
        "joint_scaffold_homology_0_30": "Joint 0.30",
        "joint_scaffold_homology_0_50_0_70_identical": "Joint 0.50/0.70",
    }
    expected_query_counts = {
        "temporal_strict_ab": 222,
        "scaffold_cold": 88,
        "joint_scaffold_homology_0_30": 19,
        "joint_scaffold_homology_0_50_0_70_identical": 22,
    }
    morgan_family = "morgan_radius2_2048_tanimoto"
    native_family = "native_sequence_3mer_tfidf_cosine"
    mmseqs_family = "mmseqs2_detected_alignment_identity"

    def selected_row(family: str, unit: str, scope: str) -> dict[str, str]:
        return require_unique_row(
            similarity_rows,
            "Figure 5 maximum-similarity summary",
            similarity_family=family,
            analysis_unit=unit,
            scope=scope,
        )

    def numeric(row: dict[str, str], field: str) -> float:
        try:
            return float(row[field])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"Figure 5 invalid {field} for {row.get('similarity_family')}/{row.get('scope')}") from exc

    query_rows = {
        morgan_family: [selected_row(morgan_family, "query", scope) for scope in scope_order],
        native_family: [selected_row(native_family, "query", scope) for scope in scope_order],
    }
    for family, rows in query_rows.items():
        for scope, row in zip(scope_order, rows):
            n_total = int(row["n_total"])
            if n_total != expected_query_counts[scope] or int(row["n_observed"]) != n_total:
                raise ValueError(f"Figure 5 query-count contract changed for {family}/{scope}")
            if int(row["n_no_detected_alignment"]) != 0:
                raise ValueError(f"Figure 5 complete similarity unexpectedly censored for {family}/{scope}")
            five_number = [numeric(row, field) for field in ("min", "q1", "median", "q3", "max")]
            if five_number != sorted(five_number) or five_number[0] < 0 or five_number[-1] > 1:
                raise ValueError(f"Figure 5 invalid five-number summary for {family}/{scope}: {five_number}")

    joint_030 = "joint_scaffold_homology_0_30"
    joint_identical = "joint_scaffold_homology_0_50_0_70_identical"
    morgan_lookup = dict(zip(scope_order, query_rows[morgan_family]))
    native_lookup = dict(zip(scope_order, query_rows[native_family]))
    if not np.isclose(numeric(morgan_lookup[joint_030], "median"), 0.421875, atol=5e-10) or not np.isclose(numeric(morgan_lookup[joint_030], "max"), 0.72, atol=5e-10):
        raise ValueError("Figure 5 Joint 0.30 Morgan summary changed")
    if not np.isclose(numeric(morgan_lookup[joint_identical], "median"), 0.4436961207, atol=5e-10) or not np.isclose(numeric(morgan_lookup[joint_identical], "max"), 1.0, atol=5e-10):
        raise ValueError("Figure 5 identical-mask Morgan summary changed")
    if not np.isclose(numeric(native_lookup[joint_030], "median"), 0.2044200450, atol=5e-10):
        raise ValueError("Figure 5 Joint 0.30 native 3-mer median changed")
    if not np.isclose(numeric(native_lookup[joint_identical], "median"), 0.2054470778, atol=5e-10):
        raise ValueError("Figure 5 identical-mask native 3-mer median changed")
    for scope in (joint_030, joint_identical):
        if not np.isclose(numeric(native_lookup[scope], "max"), 0.4691654444, atol=5e-10):
            raise ValueError(f"Figure 5 native 3-mer maximum changed for {scope}")

    mmseqs_scopes = (joint_030, joint_identical)
    mmseqs_rows = [selected_row(mmseqs_family, "query_maximum", scope) for scope in mmseqs_scopes]
    expected_detection_counts = {
        joint_030: (1, 18),
        joint_identical: (5, 17),
    }
    for scope, row in zip(mmseqs_scopes, mmseqs_rows):
        detected = int(row["n_observed"])
        censored = int(row["n_no_detected_alignment"])
        if (detected, censored) != expected_detection_counts[scope] or detected + censored != int(row["n_total"]):
            raise ValueError(f"Figure 5 MMseqs2 detection counts changed for {scope}")

    individual_mask_scopes = {
        "joint_scaffold_homology_0_50",
        "joint_scaffold_homology_0_70",
    }
    if any(
        row["scope"] in individual_mask_scopes
        and row["similarity_family"] in {morgan_family, native_family, mmseqs_family}
        for row in similarity_rows
    ):
        raise ValueError("Figure 5 found separate 0.50/0.70 rows; the identical mask must be displayed once")

    source_rows: list[dict[str, Any]] = []
    for panel, family in (("A", morgan_family), ("B", native_family)):
        for display_order, (scope, row) in enumerate(zip(scope_order, query_rows[family]), start=1):
            figure_annotation = f"median {numeric(row, 'median'):.3f}"
            if family == morgan_family and scope == joint_identical:
                figure_annotation += "; maximum = 1.000 despite scaffold absence"
            output_row = dict(row)
            output_row.update(
                {
                    "panel": panel,
                    "display_order": display_order,
                    "display_scope": f"{scope_titles[scope]} (n = {int(row['n_total'])})",
                    "figure_annotation": figure_annotation,
                    "upstream_artifact": rel(similarity_path),
                }
            )
            source_rows.append(output_row)
    for display_order, (scope, row) in enumerate(zip(mmseqs_scopes, mmseqs_rows), start=1):
        detected = int(row["n_observed"])
        total = int(row["n_total"])
        identity_min = numeric(row, "min")
        identity_max = numeric(row, "max")
        annotation = (
            f"{detected}/{total} detected; identity = {identity_min:.3f}"
            if np.isclose(identity_min, identity_max, atol=1e-15)
            else f"{detected}/{total} detected; identity range = {identity_min:.3f}–{identity_max:.3f}"
        )
        output_row = dict(row)
        output_row.update(
            {
                "panel": "C",
                "display_order": display_order,
                "display_scope": scope_titles[scope],
                "figure_annotation": annotation,
                "upstream_artifact": rel(similarity_path),
            }
        )
        source_rows.append(output_row)

    source_path = figures / "Figure_5_source_data.tsv"
    source_fields = ["panel", "display_order", "display_scope", *similarity_rows[0].keys(), "figure_annotation", "upstream_artifact"]
    write_tsv(source_path, source_fields, source_rows)

    set_style()
    fig = plt.figure(figsize=(7.3, 5.2), facecolor="white")
    outer_grid = fig.add_gridspec(
        2,
        2,
        left=0.17,
        right=0.965,
        top=0.92,
        bottom=0.15,
        height_ratios=[1.0, 0.72],
        hspace=0.64,
        wspace=0.34,
    )
    ax_morgan = fig.add_subplot(outer_grid[0, 0])
    ax_native = fig.add_subplot(outer_grid[0, 1])
    bottom_grid = outer_grid[1, :].subgridspec(1, 2, width_ratios=[0.72, 0.28], wspace=0.06)
    ax_mmseqs = fig.add_subplot(bottom_grid[0, 0])
    ax_mmseqs_labels = fig.add_subplot(bottom_grid[0, 1])

    scope_styles = {
        "temporal_strict_ab": (FIGURE_PALETTE["pale_blue"], FIGURE_PALETTE["blue"]),
        "scaffold_cold": (FIGURE_PALETTE["pale_gray"], FIGURE_PALETTE["slate"]),
        joint_030: (FIGURE_PALETTE["amber"], FIGURE_PALETTE["amber"]),
        joint_identical: (FIGURE_PALETTE["amber"], FIGURE_PALETTE["amber"]),
    }

    def draw_five_number_panel(
        ax: plt.Axes,
        rows: list[dict[str, str]],
        title: str,
        xlabel: str,
        panel: str,
        show_scope_labels: bool,
    ) -> None:
        y_positions = np.arange(len(scope_order), dtype=float)
        for y_pos, scope, row in zip(y_positions, scope_order, rows):
            minimum, q1, median, q3, maximum = [numeric(row, field) for field in ("min", "q1", "median", "q3", "max")]
            facecolor, edgecolor = scope_styles[scope]
            face_alpha = 1.0 if scope in {"temporal_strict_ab", "scaffold_cold"} else 0.30
            ax.hlines(y_pos, minimum, maximum, color=edgecolor, linewidth=0.9, zorder=2)
            ax.vlines([minimum, maximum], y_pos - 0.11, y_pos + 0.11, color=edgecolor, linewidth=0.9, zorder=2)
            ax.add_patch(
                Rectangle(
                    (q1, y_pos - 0.20),
                    q3 - q1,
                    0.40,
                    facecolor=facecolor,
                    edgecolor=edgecolor,
                    linewidth=0.9,
                    alpha=face_alpha,
                    zorder=3,
                )
            )
            ax.vlines(median, y_pos - 0.20, y_pos + 0.20, color=FIGURE_PALETTE["ink"], linewidth=1.15, zorder=4)
            ax.plot(
                median,
                y_pos,
                marker="o",
                markersize=3.6,
                markerfacecolor=FIGURE_PALETTE["ink"],
                markeredgecolor=FIGURE_PALETTE["ink"],
                markeredgewidth=0.0,
                linestyle="none",
                zorder=5,
            )
        ax.set_xlim(0, 1.02)
        ax.set_ylim(3.62, -0.50)
        ax.set_xticks(np.linspace(0, 1, 6))
        ax.set_xlabel(xlabel)
        ax.set_yticks(y_positions)
        if show_scope_labels:
            ax.set_yticklabels([f"{scope_titles[scope]} (n = {expected_query_counts[scope]})" for scope in scope_order])
        else:
            ax.set_yticklabels([])
        ax.tick_params(axis="y", length=0, pad=5)
        ax.grid(axis="x", color=FIGURE_PALETTE["grid"], linewidth=0.7)
        ax.set_title(title, loc="left", x=0.075, fontweight="bold")
        panel_label(ax, panel)

    draw_five_number_panel(
        ax_morgan,
        query_rows[morgan_family],
        "Maximum Morgan similarity to history",
        "Maximum Morgan Tanimoto similarity",
        "a",
        True,
    )
    draw_five_number_panel(
        ax_native,
        query_rows[native_family],
        "Maximum native 3-mer similarity to history",
        "Maximum native 3-mer cosine similarity",
        "b",
        False,
    )
    ax_morgan.annotate(
        "Max = 1.000 despite\nscaffold absence",
        xy=(numeric(morgan_lookup[joint_identical], "max"), 3.0),
        xytext=(0.66, 3.42),
        ha="left",
        va="center",
        fontsize=7.1,
        color=FIGURE_PALETTE["slate"],
        arrowprops={"arrowstyle": "-", "color": FIGURE_PALETTE["slate"], "linewidth": 0.7},
        annotation_clip=False,
    )

    y_positions = np.arange(len(mmseqs_scopes), dtype=float)
    detected_percentages: list[float] = []
    censored_percentages: list[float] = []
    mmseqs_annotations: list[str] = []
    for row in mmseqs_rows:
        total = int(row["n_total"])
        detected = int(row["n_observed"])
        detected_percentages.append(100.0 * detected / total)
        censored_percentages.append(100.0 * int(row["n_no_detected_alignment"]) / total)
        identity_min = numeric(row, "min")
        identity_max = numeric(row, "max")
        mmseqs_annotations.append(
            f"{detected}/{total} detected (ID {identity_min:.3f})"
            if np.isclose(identity_min, identity_max, atol=1e-15)
            else f"{detected}/{total} detected (ID {identity_min:.3f}–{identity_max:.3f})"
        )
    ax_mmseqs.barh(y_positions, detected_percentages, height=0.34, color=FIGURE_PALETTE["teal"], edgecolor="white", linewidth=0.8, label="Alignment detected")
    ax_mmseqs.barh(
        y_positions,
        censored_percentages,
        left=detected_percentages,
        height=0.34,
        color=FIGURE_PALETTE["pale_gray"],
        edgecolor=FIGURE_PALETTE["slate"],
        linewidth=0.7,
        hatch="///",
        label="No alignment detected (censored)",
    )
    ax_mmseqs.set_xlim(0, 100)
    ax_mmseqs.set_ylim(1.52, -0.52)
    ax_mmseqs.set_xticks(np.arange(0, 101, 20))
    ax_mmseqs.xaxis.set_major_formatter(PercentFormatter(xmax=100, decimals=0))
    ax_mmseqs.set_xlabel("Queries (%)")
    ax_mmseqs.set_yticks(y_positions, [scope_titles[scope] for scope in mmseqs_scopes])
    ax_mmseqs.tick_params(axis="y", length=0, pad=5)
    ax_mmseqs.grid(axis="x", color=FIGURE_PALETTE["grid"], linewidth=0.7)
    ax_mmseqs.set_title("MMseqs2 alignment detection in joint scopes", loc="left", x=0.075, fontweight="bold")
    panel_label(ax_mmseqs, "c")
    legend_handles = [
        Patch(facecolor=FIGURE_PALETTE["teal"], edgecolor="none", label="Alignment detected"),
        Patch(facecolor=FIGURE_PALETTE["pale_gray"], edgecolor=FIGURE_PALETTE["slate"], hatch="///", label="No alignment detected (censored)"),
    ]
    ax_mmseqs.legend(handles=legend_handles, frameon=False, ncol=2, loc="upper right", bbox_to_anchor=(1.0, 0.98), fontsize=7.5, borderaxespad=0)

    ax_mmseqs_labels.set_xlim(0, 1)
    ax_mmseqs_labels.set_ylim(1.52, -0.52)
    ax_mmseqs_labels.axis("off")
    for y_pos, annotation in zip(y_positions, mmseqs_annotations):
        ax_mmseqs_labels.text(0.0, y_pos, annotation, ha="left", va="center", fontsize=7.8, color=FIGURE_PALETTE["ink"])
    fig.text(0.17, 0.045, "Non-detections are censored and are not assigned identity zero.", ha="left", va="bottom", fontsize=7.4, color=FIGURE_PALETTE["slate"])

    png_output = figures / "Figure_5_cold_scope_similarity.png"
    pdf_output = figures / "Figure_5_cold_scope_similarity.pdf"
    svg_output = figures / "Figure_5_cold_scope_similarity.svg"
    fig.savefig(pdf_output, bbox_inches="tight", facecolor="white", format="pdf")
    fig.savefig(svg_output, bbox_inches="tight", facecolor="white", format="svg")
    save_rgb_png(fig, png_output)
    if not pdf_output.is_file() or pdf_output.stat().st_size == 0 or not pdf_output.read_bytes().startswith(b"%PDF"):
        raise ValueError(f"PDF output contract failed for {pdf_output.name}")
    if not svg_output.is_file() or svg_output.stat().st_size == 0 or b"<svg" not in svg_output.read_bytes()[:2048]:
        raise ValueError(f"SVG output contract failed for {svg_output.name}")
    upstream = [
        {
            "path": rel(similarity_path),
            "sha256": sha256_file(similarity_path),
            "role": "maximum_similarity_five_number_summaries_and_mmseqs2_censoring",
        }
    ]
    return png_output, source_path, upstream


def supplementary_figure_s1(project: Path, tables: Path, figures: Path) -> tuple[Path, Path, list[dict[str, str]]]:
    """Move the frozen top-100/exhaustive fidelity heatmap to Supplementary Figure S1."""

    top100_path = tables / "Table_S2_top100_exhaustive_fidelity.tsv"
    top100_rows = read_tsv(top100_path)
    top100_scope_inputs = {
        "temporal_strict_ab": ("temporal_strict_ab",),
        "scaffold_cold_strict_ab": ("scaffold_cold_strict_ab",),
        "project_defined_joint_scaffold_homology_cold_0_30": ("double_cold_0_30",),
        "project_defined_joint_scaffold_homology_cold_0_50_0_70_identical_mask": ("double_cold_0_50", "double_cold_0_70"),
    }

    def rel(path: Path) -> str:
        return path.relative_to(project).as_posix()

    selected_rows: list[dict[str, Any]] = []
    fidelity_matrix = np.zeros((len(DISPLAY_SCOPES), len(METRICS)), dtype=float)
    for row_index, display_scope in enumerate(DISPLAY_SCOPES):
        input_scopes = top100_scope_inputs[display_scope]
        for column_index, metric in enumerate(METRICS):
            rows = [
                require_unique_row(
                    top100_rows,
                    "Supplementary Figure S1 metric difference",
                    record_type="metric_difference",
                    scope=input_scope,
                    metric_or_item=metric,
                )
                for input_scope in input_scopes
            ]
            if any(row["status"] != "post_hoc_sensitivity_only" for row in rows):
                raise ValueError(f"Supplementary Figure S1 status changed for {display_scope}/{metric}")
            first = rows[0]
            first_values = (
                int(first["relation_count"]),
                int(first["query_count"]),
                float(first["top100_value"]),
                float(first["exhaustive_value"]),
                float(first["exhaustive_minus_top100"]),
            )
            for other in rows[1:]:
                other_values = (
                    int(other["relation_count"]),
                    int(other["query_count"]),
                    float(other["top100_value"]),
                    float(other["exhaustive_value"]),
                    float(other["exhaustive_minus_top100"]),
                )
                if any(not np.isclose(left, right, atol=1e-15) for left, right in zip(first_values, other_values)):
                    raise ValueError(f"Supplementary Figure S1 duplicated 0.50/0.70 masks diverged for {metric}")
            if not np.isclose(first_values[3] - first_values[2], first_values[4], atol=1e-15):
                raise ValueError(f"Supplementary Figure S1 arithmetic changed for {display_scope}/{metric}")
            fidelity_matrix[row_index, column_index] = first_values[4]
            selected_rows.append(
                {
                    "record_type": "exhaustive_minus_frozen_top100",
                    "display_scope": SCOPE_LABELS[display_scope],
                    "provenance_scope": display_scope if len(input_scopes) == 1 else "double_cold_0_50_and_0_70_identical",
                    "baseline": "structure_sequence_pair_neighbor",
                    "metric": metric,
                    "relation_count": first_values[0],
                    "query_count": first_values[1],
                    "top100_value": first_values[2],
                    "exhaustive_value": first_values[3],
                    "exhaustive_minus_top100": first_values[4],
                    "status": first["status"] + ("; identical_0_50_0_70_mask_merged" if len(rows) == 2 else ""),
                    "upstream_artifact": rel(top100_path),
                }
            )

    source_path = figures / "Figure_S1_source_data.tsv"
    write_tsv(
        source_path,
        [
            "record_type",
            "display_scope",
            "provenance_scope",
            "baseline",
            "metric",
            "relation_count",
            "query_count",
            "top100_value",
            "exhaustive_value",
            "exhaustive_minus_top100",
            "status",
            "upstream_artifact",
        ],
        selected_rows,
    )

    set_style()
    fig, ax = plt.subplots(figsize=(7.1, 3.7), facecolor="white")
    fig.subplots_adjust(left=0.19, right=0.90, top=0.82, bottom=0.22)
    fidelity_limit = 0.0016
    image = ax.imshow(fidelity_matrix, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-fidelity_limit, vcenter=0.0, vmax=fidelity_limit))
    for row_index, values in enumerate(fidelity_matrix):
        for column_index, value in enumerate(values):
            label = "0" if np.isclose(value, 0.0, atol=1e-15) else f"{value:+.1e}"
            text_color = "white" if abs(value) > 0.58 * fidelity_limit else FIGURE_PALETTE["ink"]
            ax.text(column_index, row_index, label, ha="center", va="center", fontsize=8.0, color=text_color)
    ax.set_xticks(np.arange(len(METRICS)), METRICS, rotation=22, ha="right")
    ax.set_yticks(
        np.arange(len(DISPLAY_SCOPES)),
        ["Temporal", "Scaffold", "Joint 0.30", "Joint 0.50/0.70\nidentical mask"],
    )
    ax.tick_params(length=0)
    ax.set_title("Pair-neighbour exhaustive − frozen top-100 fidelity", loc="left", fontsize=9.2, fontweight="bold", pad=8)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.042, pad=0.025)
    colorbar.set_label("Metric difference", fontsize=8.0)
    colorbar.ax.tick_params(labelsize=8.0)
    fig.text(
        0.50,
        0.055,
        "This is a post hoc fidelity sensitivity and does not replace the locked top-100 primary representation.",
        ha="center",
        va="bottom",
        fontsize=8.0,
        color=FIGURE_PALETTE["slate"],
    )

    png_output = figures / "Figure_S1_top100_exhaustive_fidelity.png"
    save_rgb_png(fig, png_output)
    upstream = [
        {"path": rel(top100_path), "sha256": sha256_file(top100_path), "role": "exhaustive_minus_frozen_top100_pair_neighbour_fidelity"},
    ]
    return png_output, source_path, upstream


def main() -> int:
    project = Path(__file__).resolve().parents[2]
    manuscript = project / "manuscript_molecular_diversity_v4_20260729"
    tables = manuscript / "tables"
    figures = manuscript / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    figure1_path, figure1_source, figure1_upstream = figure1(project, figures)
    figure2_path, figure2_source, figure2_upstream = figure2(project, figures)
    figure3_path, figure3_source, figure3_upstream = figure3(project, tables, figures)
    figure4_path, figure4_source, figure4_upstream = figure4(project, tables, figures)
    figure5_path, figure5_source, figure5_upstream = figure5(project, tables, figures)
    figure_s1_path, figure_s1_source, figure_s1_upstream = supplementary_figure_s1(project, tables, figures)
    manifest_path = figures / "revision_figure_generation_manifest_v4.json"
    script_hash = sha256_file(Path(__file__))

    def image_receipt(path: Path) -> dict[str, Any]:
        with Image.open(path) as image:
            dpi = image.info.get("dpi", (0.0, 0.0))
            if image.mode != "RGB" or min(dpi) < 599.0:
                raise ValueError(f"Figure output contract failed for {path.name}: mode={image.mode}, dpi={dpi}")
            return {
                "color_mode": image.mode,
                "dpi": [round(float(dpi[0]), 3), round(float(dpi[1]), 3)],
                "pixel_dimensions": [int(image.width), int(image.height)],
            }

    payload = {
        "schema_version": "1.5",
        "format": "PNG_only",
        "identifier_bearing_outputs": False,
        "narrative_term_double_cold_used_in_figures": False,
        "figures": [
            {
                "figure": figure1_path.name,
                "sha256": sha256_file(figure1_path),
                "bytes": figure1_path.stat().st_size,
                **image_receipt(figure1_path),
                "source_data": figure1_source.name,
                "source_data_sha256": sha256_file(figure1_source),
                "generation_script": Path(__file__).name,
                "generation_script_sha256": script_hash,
                "provenance": "redrawn_from_current_corrective_aggregate_sources; legacy copied preview and its broader-lineage C37 835/835 value removed",
                "upstream_artifacts": figure1_upstream,
            },
            {
                "figure": figure2_path.name,
                "sha256": sha256_file(figure2_path),
                "bytes": figure2_path.stat().st_size,
                **image_receipt(figure2_path),
                "source_data": figure2_source.name,
                "source_data_sha256": sha256_file(figure2_source),
                "generation_script": Path(__file__).name,
                "generation_script_sha256": script_hash,
                "upstream_artifacts": figure2_upstream,
            },
            {
                "figure": figure3_path.name,
                "sha256": sha256_file(figure3_path),
                "bytes": figure3_path.stat().st_size,
                **image_receipt(figure3_path),
                "source_data": figure3_source.name,
                "source_data_sha256": sha256_file(figure3_source),
                "generation_script": Path(__file__).name,
                "generation_script_sha256": script_hash,
                "upstream_artifacts": figure3_upstream,
            },
            {
                "figure": figure4_path.name,
                "figure_label": "Figure 4",
                "sha256": sha256_file(figure4_path),
                "bytes": figure4_path.stat().st_size,
                **image_receipt(figure4_path),
                "source_data": figure4_source.name,
                "source_data_sha256": sha256_file(figure4_source),
                "generation_script": Path(__file__).name,
                "generation_script_sha256": script_hash,
                "provenance": "date_weight_structure_and_exhaustive_policy_robustness_split_from_Figure_3",
                "upstream_artifacts": figure4_upstream,
            },
            {
                "figure": figure5_path.name,
                "figure_label": "Figure 5",
                "sha256": sha256_file(figure5_path),
                "bytes": figure5_path.stat().st_size,
                **image_receipt(figure5_path),
                "source_data": figure5_source.name,
                "source_data_sha256": sha256_file(figure5_source),
                "generation_script": Path(__file__).name,
                "generation_script_sha256": script_hash,
                "provenance": "aggregate_only_representation_dependent_similarity_and_mmseqs2_censoring",
                "upstream_artifacts": figure5_upstream,
            },
            {
                "figure": figure_s1_path.name,
                "figure_label": "Supplementary Figure S1",
                "sha256": sha256_file(figure_s1_path),
                "bytes": figure_s1_path.stat().st_size,
                **image_receipt(figure_s1_path),
                "source_data": figure_s1_source.name,
                "source_data_sha256": sha256_file(figure_s1_source),
                "generation_script": Path(__file__).name,
                "generation_script_sha256": script_hash,
                "provenance": "post_hoc_pair_neighbour_fidelity_sensitivity; frozen_top100_remains_primary",
                "upstream_artifacts": figure_s1_upstream,
            },
        ],
    }
    with manifest_path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    print(json.dumps(payload, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
