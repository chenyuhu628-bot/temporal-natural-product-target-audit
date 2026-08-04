#!/usr/bin/env python3
"""Generate aggregate-only PNG figures for the corrective manuscript v3."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
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
COLORS = {
    "weighted_target_popularity": "#5B7083",
    "sequence_3mer_transfer": "#B279A2",
    "weighted_morgan_transfer": "#2A9D8F",
    "structure_sequence_pair_neighbor": "#E76F51",
}
SCOPES = [
    "temporal_strict_ab",
    "scaffold_cold_strict_ab",
    "double_cold_0_30",
    "double_cold_0_50",
    "double_cold_0_70",
]
SCOPE_LABELS = {
    "temporal_strict_ab": "Temporal strict A/B",
    "scaffold_cold_strict_ab": "Scaffold cold",
    "double_cold_0_30": "Double cold 0.30",
    "double_cold_0_50": "Double cold 0.50",
    "double_cold_0_70": "Double cold 0.70",
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


def set_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 8.5,
            "axes.titlesize": 10.5,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.8,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(-0.12, 1.06, label, transform=ax.transAxes, fontsize=12, fontweight="bold", va="top")


def figure1(project: Path, figures: Path) -> tuple[Path, Path, list[dict[str, str]]]:
    """Draw the current corrective lineage as a linear, aggregate-only workflow."""

    rebuild_path = project / "author_run_strict_ab_asof_cutoff_execution_v1_20260728" / "audit" / "asof_rebuild_summary.json"
    future_path = project / "results" / "strict_temporal_future_v1_1_pmid_verified_chembl31_leakage_gate_summary.json"
    source_path = project / "author_run_strict_ab_asof_cutoff_execution_v1_20260728" / "audit" / "source_concentration_v1" / "source_concentration_aggregate_summary.json"
    score_path = project / "author_run_strict_ab_asof_cutoff_execution_v1_20260728" / "score" / "corrective_score_manifest.json"
    legacy_c37_summary_path = project / "results" / "chembl37_p2_source_overlap_audit_v1_summary.json"
    legacy_c37_config_path = project / "configs" / "chembl37_p2_manual_review_queue_v1_1.json"

    rebuild = read_json(rebuild_path)
    future = read_json(future_path)
    source = read_json(source_path)
    score = read_json(score_path)
    legacy_c37_summary = read_json(legacy_c37_summary_path)
    legacy_c37_config = read_json(legacy_c37_config_path)

    counts = rebuild["counts"]
    rows = rebuild["row_eligibility_counts"]
    endpoint_status = future["status_counts_in_frozen_future_table"]
    endpoint_cohort = next(item for item in source["cohorts"] if item["cohort"] == "endpoint")
    historical_cohort = next(item for item in source["cohorts"] if item["cohort"] == "historical")
    overlap = source["cross_cohort_source_overlap"][0]

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

    def rel(path: Path) -> str:
        return path.relative_to(project).as_posix()

    source_rows: list[dict[str, Any]] = [
        {"panel": "A", "item": "cutoff", "value": rebuild["cutoff"], "status": "current_corrective_branch", "upstream_artifact": rel(rebuild_path)},
        {"panel": "A", "item": "historical_strict_v2_source_rows", "value": counts["historical_strict_v2_rows"], "status": "verified", "upstream_artifact": rel(rebuild_path)},
        {"panel": "A", "item": "eligible_pre_cutoff", "value": rows["eligible_pre_cutoff"], "status": "verified", "upstream_artifact": rel(rebuild_path)},
        {"panel": "A", "item": "excluded_non_day_precision", "value": rows["excluded_non_day_precision"], "status": "verified", "upstream_artifact": rel(rebuild_path)},
        {"panel": "A", "item": "excluded_not_numeric_pmid", "value": rows["excluded_not_numeric_pmid"], "status": "verified", "upstream_artifact": rel(rebuild_path)},
        {"panel": "A", "item": "excluded_pubmed_not_found", "value": rows["excluded_pubmed_not_found"], "status": "verified", "upstream_artifact": rel(rebuild_path)},
        {"panel": "A", "item": "excluded_after_cutoff", "value": rows["excluded_after_cutoff"], "status": "verified", "upstream_artifact": rel(rebuild_path)},
        {"panel": "A", "item": "historical_pairs_with_any_excluded_row", "value": counts["historical_pairs_with_any_excluded_row"], "status": "verified", "upstream_artifact": rel(rebuild_path)},
        {"panel": "A", "item": "historical_tier_weight_changes", "value": counts["historical_tier_weight_changes"], "status": "verified", "upstream_artifact": rel(rebuild_path)},
        {"panel": "A", "item": "historical_relation_keys", "value": counts["historical_pairs"], "status": "membership_frozen", "upstream_artifact": rel(rebuild_path)},
        {"panel": "A", "item": "historical_targets", "value": counts["historical_targets"], "status": "membership_frozen", "upstream_artifact": rel(rebuild_path)},
        {"panel": "B", "item": "initial_later_candidate_pairs", "value": future["inputs"]["frozen_future_pair_count"], "status": "frozen_prior_endpoint_lineage", "upstream_artifact": rel(future_path)},
        {"panel": "B", "item": "excluded_C31_historical_activity", "value": endpoint_status["historical_activity_recorded_in_chembl31"], "status": "excluded", "upstream_artifact": rel(future_path)},
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
        {"panel": "D", "item": "historical_numeric_pmids", "value": historical_cohort["unique_source_document_count"], "status": "current_corrective_branch", "upstream_artifact": rel(source_path)},
        {"panel": "D", "item": "shared_pmids_current_endpoint_vs_history", "value": overlap["shared_source_document_count"], "status": "provenance_only_not_external_validation", "upstream_artifact": rel(source_path)},
        {"panel": "not_displayed", "item": "legacy_C37_shared_PMID_candidates", "value": legacy_c37_summary["pairs_with_at_least_one_shared_PMID"], "status": "removed_from_figure_broader_846_pair_lineage_not_current_358_endpoint", "upstream_artifact": rel(legacy_c37_summary_path)},
    ]
    source_data_path = figures / "Figure_1_source_data.tsv"
    write_tsv(source_data_path, ["panel", "item", "value", "status", "upstream_artifact"], source_rows)

    set_style()
    fig = plt.figure(figsize=(12.0, 6.2), facecolor="white")
    grid = fig.add_gridspec(1, 4, left=0.025, right=0.985, top=0.975, bottom=0.08, wspace=0.065)

    navy = "#19324A"
    teal = "#2A9D8F"
    blue = "#277DA1"
    orange = "#E76F51"
    amber = "#F4A261"
    purple = "#6C5B7B"
    pale = "#F5F7F8"
    line = "#CBD5DB"
    ink = "#24333D"

    def make_panel(index: int, title: str, accent: str) -> plt.Axes:
        ax = fig.add_subplot(grid[0, index])
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")
        ax.add_patch(FancyBboxPatch((0.015, 0.015), 0.97, 0.97, boxstyle="round,pad=0.012,rounding_size=0.025", facecolor="white", edgecolor=line, linewidth=1.2))
        ax.add_patch(FancyBboxPatch((0.035, 0.905), 0.11, 0.065, boxstyle="round,pad=0.008,rounding_size=0.025", facecolor=accent, edgecolor=accent))
        ax.text(0.09, 0.938, ("A", "B", "C", "D")[index], ha="center", va="center", color="white", fontsize=10.5, fontweight="bold")
        ax.text(0.17, 0.938, title, ha="left", va="center", color=navy, fontsize=9.2, fontweight="bold")
        return ax

    def box(ax: plt.Axes, xy: tuple[float, float], wh: tuple[float, float], text_value: str, *, face: str = pale, edge: str = line, color: str = ink, size: float = 7.4, weight: str = "normal") -> None:
        x, y = xy
        w, h = wh
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.010,rounding_size=0.018", facecolor=face, edgecolor=edge, linewidth=1.0))
        ax.text(x + w / 2, y + h / 2, text_value, ha="center", va="center", color=color, fontsize=size, fontweight=weight, linespacing=1.25)

    def down(ax: plt.Axes, y0: float, y1: float, *, color: str = "#8697A0") -> None:
        ax.add_patch(FancyArrowPatch((0.50, y0), (0.50, y1), arrowstyle="-|>", mutation_scale=10, linewidth=1.0, color=color))

    ax = make_panel(0, "Historical source repair", teal)
    box(ax, (0.09, 0.785), (0.82, 0.10), "NPASS v2 strict A/B candidate rows\n20,647", face="#EAF6F3", edge=teal, weight="bold")
    down(ax, 0.782, 0.720, color=teal)
    box(ax, (0.11, 0.625), (0.78, 0.09), "Row-level PMID + date gate\ncutoff: 31 August 2022", face="#EFF5F8", edge=blue, weight="bold")
    down(ax, 0.622, 0.565, color=blue)
    box(ax, (0.07, 0.465), (0.86, 0.095), "13,885 eligible pre-cutoff rows", face="#EAF6F3", edge=teal, color="#17695F", size=8.1, weight="bold")
    box(ax, (0.07, 0.315), (0.86, 0.125), "Excluded from historical features\n6,570 non-day precision · 192 non-PMID/\nnonnumeric · 0 PubMed-not-found · 0 after cutoff", face="#FFF1EC", edge=orange, color="#8F3F2E", size=6.5)
    box(ax, (0.07, 0.190), (0.86, 0.095), "1,222 keys had ≥1 excluded row\n166 tier/weight assignments changed", face="#FFF7E8", edge=amber, color="#7A521C", size=7.0)
    down(ax, 0.185, 0.142, color=teal)
    box(ax, (0.07, 0.055), (0.86, 0.085), "4,990 frozen relation keys\n1,131 targets", face=teal, edge=teal, color="white", size=7.3, weight="bold")

    ax = make_panel(1, "Endpoint eligibility", blue)
    box(ax, (0.07, 0.775), (0.86, 0.105), "NPASS v3 later-recorded candidates\n442 canonical pairs", face="#EAF2F8", edge=blue, weight="bold")
    box(ax, (0.07, 0.635), (0.86, 0.105), "Frozen temporal rule\nv3-only · strict A/B · all PMIDs day-precise after cutoff\nno non-PMID reference on the pair", face=pale, edge=line, size=6.6)
    down(ax, 0.630, 0.575, color=orange)
    box(ax, (0.10, 0.485), (0.80, 0.085), "−19 exact-entity C31\nhistorical-activity overlaps", face="#FFF1EC", edge=orange, color="#8F3F2E", size=6.8)
    down(ax, 0.480, 0.425, color=orange)
    box(ax, (0.10, 0.335), (0.80, 0.085), "−65 entity-unresolved pairs", face="#FFF1EC", edge=orange, color="#8F3F2E", size=7.2)
    down(ax, 0.330, 0.265, color=blue)
    box(ax, (0.07, 0.125), (0.86, 0.135), "FINAL FROZEN ENDPOINT\n358 relations\n222 queries · 156 targets", face=blue, edge=blue, color="white", size=8.0, weight="bold")
    ax.text(0.50, 0.065, "No excluded pair was relabelled as a negative", ha="center", va="center", fontsize=6.8, color="#5B6870")

    ax = make_panel(2, "Fixed retrieval task", purple)
    box(ax, (0.15, 0.795), (0.70, 0.085), "Query compound q", face="#F1EDF6", edge=purple, color="#4E3F62", weight="bold")
    down(ax, 0.790, 0.735, color=purple)
    box(ax, (0.09, 0.650), (0.82, 0.080), "Mask corrected historical targets of q", face="#FFF7E8", edge=amber, color="#7A521C", size=7.0)
    down(ax, 0.645, 0.590, color=purple)
    box(ax, (0.09, 0.505), (0.82, 0.080), "Fixed universe\n4,123 human single-protein targets", face="#EFF5F8", edge=blue, size=6.8, weight="bold")
    down(ax, 0.500, 0.445, color=purple)
    box(ax, (0.09, 0.320), (0.82, 0.120), "Four unchanged baselines\nPopularity · sequence 3-mer\nMorgan transfer · pair neighbour", face="#F1EDF6", edge=purple, color="#4E3F62", size=7.1)
    down(ax, 0.315, 0.260, color=purple)
    box(ax, (0.15, 0.185), (0.70, 0.070), "Complete ranked target list", face=purple, edge=purple, color="white", size=7.3, weight="bold")
    down(ax, 0.180, 0.135, color=blue)
    box(ax, (0.09, 0.055), (0.82, 0.075), "Later strict A/B targets\nloaded only for evaluation", face="#EAF2F8", edge=blue, color="#164F6A", size=6.6, weight="bold")

    ax = make_panel(3, "Interpretation + provenance", amber)
    box(ax, (0.07, 0.775), (0.86, 0.105), "CURRENT STRICT A/B AUDIT\n486 endpoint evidence rows\n124 numeric PMIDs", face="#FFF7E8", edge=amber, color="#6E4A19", size=6.9, weight="bold")
    box(ax, (0.07, 0.620), (0.86, 0.120), "95 query–PMID components\nlargest component: 51 of 222 queries", face="#F1EDF6", edge=purple, color="#4E3F62", size=7.2)
    box(ax, (0.07, 0.475), (0.86, 0.105), "0 PMIDs shared with eligible history\n(5,108 historical PMIDs)", face="#EAF6F3", edge=teal, color="#17695F", size=7.3, weight="bold")
    box(ax, (0.07, 0.310), (0.86, 0.125), "0 shared PMIDs ≠ external validation\nsame database lineage · outcome-visible\nauthor-run", face="#FFF1EC", edge=orange, color="#8F3F2E", size=6.7, weight="bold")
    box(ax, (0.07, 0.175), (0.86, 0.095), "Unrecorded ≠ negative", face=pale, edge=line, size=7.5, weight="bold")
    box(ax, (0.07, 0.055), (0.86, 0.095), "Later recorded ≠ first biological discovery", face=pale, edge=line, size=6.9, weight="bold")

    for x0, x1 in ((0.247, 0.263), (0.493, 0.509), (0.739, 0.755)):
        fig.add_artist(FancyArrowPatch((x0, 0.46), (x1, 0.46), transform=fig.transFigure, arrowstyle="-|>", mutation_scale=12, linewidth=1.2, color="#8A9AA3"))

    fig.text(0.50, 0.025, "Green: eligible historical evidence   Blue: frozen endpoint/evaluation   Orange: exclusions and claim boundaries   Purple: scoring and audit", ha="center", va="bottom", fontsize=7.2, color="#53636C")
    output = figures / "Figure_1_source_aware_temporal_endpoint.png"
    save_rgb_png(fig, output)

    upstream = [
        {"path": rel(rebuild_path), "sha256": sha256_file(rebuild_path), "role": "current_row_level_rebuild"},
        {"path": rel(future_path), "sha256": sha256_file(future_path), "role": "frozen_endpoint_flow"},
        {"path": rel(source_path), "sha256": sha256_file(source_path), "role": "current_strict_ab_provenance_audit"},
        {"path": rel(score_path), "sha256": sha256_file(score_path), "role": "fixed_retrieval_task"},
        {"path": rel(legacy_c37_summary_path), "sha256": sha256_file(legacy_c37_summary_path), "role": "excluded_legacy_broader_lineage_not_displayed"},
        {"path": rel(legacy_c37_config_path), "sha256": sha256_file(legacy_c37_config_path), "role": "excluded_legacy_846_pair_lineage_definition"},
    ]
    return output, source_data_path, upstream


def figure2(tables: Path, figures: Path) -> tuple[Path, Path]:
    bootstrap = read_tsv(tables / "Table_4_corrected_bootstrap_summaries.tsv")
    aggregate = read_tsv(tables / "Table_3_corrected_aggregate_performance.tsv")
    attribution = read_tsv(tables / "Table_S3_zero_and_failure_accounting.tsv")
    lookup_boot = {
        (row["scope"], row["baseline_or_left"], row["metric"]): row
        for row in bootstrap
        if row["record_type"] == "baseline"
    }
    lookup_agg = {(row["scope"], row["baseline"]): row for row in aggregate}
    attribution_lookup = {(row["scope"], row["baseline"]): row for row in attribution}

    source_rows: list[dict[str, Any]] = []
    for scope in ("temporal_strict_ab", "scaffold_cold_strict_ab"):
        for baseline in BASELINES:
            row = lookup_boot[(scope, baseline, "Recall@50")]
            source_rows.append(
                {
                    "panel": "primary_recall50",
                    "scope": scope,
                    "baseline": baseline,
                    "metric": "Recall@50",
                    "estimate": row["estimate"],
                    "ci95_low": row["ci95_low"],
                    "ci95_high": row["ci95_high"],
                    "attribution_class": "salt_derived_all_zero_scope" if scope == "scaffold_cold_strict_ab" and baseline == "sequence_3mer_transfer" else "mixed_or_score_separated",
                    "status": row["status"],
                }
            )
    for scope in ("double_cold_0_30", "double_cold_0_50", "double_cold_0_70"):
        for baseline in BASELINES:
            row = lookup_agg[(scope, baseline)]
            attribution_row = attribution_lookup[(scope, baseline)]
            if int(attribution_row["endpoint_zero_score_relation_count"]) == int(attribution_row["candidate_relation_count"]):
                attribution_class = "zero_score_tie_salt_derived"
            elif int(attribution_row["endpoint_positive_unique_score_relation_count"]) == int(attribution_row["candidate_relation_count"]) and int(attribution_row["endpoint_rank_le_50_relation_count"]) == 0:
                attribution_class = "positive_unique_scores_all_below_50"
            else:
                attribution_class = "mixed"
            for metric in METRICS:
                source_rows.append(
                    {
                        "panel": "double_cold_all_metrics",
                        "scope": scope,
                        "baseline": baseline,
                        "metric": metric,
                        "estimate": row[metric],
                        "ci95_low": "",
                        "ci95_high": "",
                        "attribution_class": attribution_class,
                        "status": row["status"],
                    }
                )

    source_path = figures / "Figure_2_source_data.tsv"
    write_tsv(source_path, ["panel", "scope", "baseline", "metric", "estimate", "ci95_low", "ci95_high", "attribution_class", "status"], source_rows)

    set_style()
    fig = plt.figure(figsize=(10.8, 7.6), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.55])

    for column, scope in enumerate(("temporal_strict_ab", "scaffold_cold_strict_ab")):
        ax = fig.add_subplot(grid[0, column])
        y = np.arange(len(BASELINES))
        for yi, baseline in enumerate(BASELINES):
            row = lookup_boot[(scope, baseline, "Recall@50")]
            estimate = float(row["estimate"])
            low = float(row["ci95_low"])
            high = float(row["ci95_high"])
            ax.errorbar(
                estimate,
                yi,
                xerr=np.array([[estimate - low], [high - estimate]]),
                fmt="o",
                ms=6,
                capsize=3,
                lw=1.4,
                color=COLORS[baseline],
                ecolor=COLORS[baseline],
            )
            suffix = "†" if scope == "scaffold_cold_strict_ab" and baseline == "sequence_3mer_transfer" else ""
            ax.text(high + 0.008, yi, f"{estimate:.3f}{suffix}", va="center", fontsize=7.5)
        ax.set_yticks(y, [BASELINE_LABELS[item] for item in BASELINES])
        ax.invert_yaxis()
        ax.set_xlim(0, 0.36)
        ax.set_xlabel("Macro Recall@50 (primary query-bootstrap 95% CI)")
        count = "358 relations / 222 queries" if column == 0 else "123 relations / 88 queries"
        ax.set_title(f"{SCOPE_LABELS[scope]}\n{count}", loc="left", fontweight="bold")
        ax.grid(axis="x", color="#D9E0E5", lw=0.7)
        panel_label(ax, "A" if column == 0 else "B")

    ax = fig.add_subplot(grid[1, :])
    heat_rows = [(scope, baseline) for scope in SCOPES[2:] for baseline in BASELINES]
    matrix = np.array([[float(lookup_agg[(scope, baseline)][metric]) for metric in METRICS] for scope, baseline in heat_rows])
    vmax = float(np.max(matrix)) if np.max(matrix) > 0 else 1.0
    image = ax.imshow(matrix, aspect="auto", cmap="Blues", norm=Normalize(vmin=0, vmax=vmax))
    for row_index, values in enumerate(matrix):
        for column_index, value in enumerate(values):
            if value == 0:
                label = "0"
            else:
                baseline = heat_rows[row_index][1]
                suffix = "‡" if baseline == "structure_sequence_pair_neighbor" else "†"
                label = f"{value:.1e}{suffix}"
            color = "white" if value > 0.58 * vmax else "#183247"
            ax.text(column_index, row_index, label, ha="center", va="center", fontsize=7.2, color=color)
    ax.set_xticks(np.arange(len(METRICS)), METRICS)
    ylabels = [f"{SCOPE_LABELS[scope]} — {BASELINE_SHORT[baseline]}" for scope, baseline in heat_rows]
    ax.set_yticks(np.arange(len(heat_rows)), ylabels)
    ax.tick_params(length=0)
    ax.set_title("Strict double-cold stress scopes: all planned metrics retained", loc="left", fontweight="bold")
    panel_label(ax, "C")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.018, pad=0.015)
    colorbar.set_label("Raw macro metric value")
    fig.text(
        0.5,
        -0.035,
        "Recall and NDCG are exactly zero in every double-cold cell; the 0.50 and 0.70 masks are identical.\n† Post hoc, outcome-visible attribution: zero-score endpoint ties, ranks set by salt (including scaffold 3-mer in B).  ‡ Pair: positive unique endpoint scores, but all observed relations rank below 50.",
        ha="center",
        va="top",
        fontsize=7.4,
        color="#4A5964",
    )

    output = figures / "Figure_2_corrected_retrieval_performance.png"
    save_rgb_png(fig, output)
    return output, source_path


def figure3(tables: Path, figures: Path) -> tuple[Path, Path]:
    source_dependence = read_tsv(tables / "Table_S4_pmid_document_dependence.tsv")
    top100 = read_tsv(tables / "Table_S2_top100_exhaustive_fidelity.tsv")

    focus_rows = [row for row in source_dependence if row["record_type"] == "document_component_focus_contrast"]
    baseline_component_rows = [row for row in source_dependence if row["record_type"] == "document_component_baseline_metric"]
    scope_rows = [row for row in source_dependence if row["record_type"] == "document_component_scope"]
    difference_rows = [row for row in top100 if row["record_type"] == "metric_difference"]

    focus_lookup = {(row["scope_or_cohort"], row["metric_or_item"]): row for row in focus_rows}
    component_lookup = {
        (row["scope_or_cohort"], row["baseline_or_left"], row["metric_or_item"]): row
        for row in baseline_component_rows
    }
    scope_lookup = {(row["scope_or_cohort"], row["metric_or_item"]): row["point_estimate"] for row in scope_rows}
    difference_lookup = {(row["scope"], row["metric_or_item"]): row for row in difference_rows}

    selected_rows: list[dict[str, Any]] = []
    for scope in ("temporal_strict_ab", "scaffold_cold_strict_ab"):
        for metric in METRICS:
            row = focus_lookup[(scope, metric)]
            selected_rows.append({"panel": "focus_contrast", **row})
        for baseline in BASELINES:
            row = component_lookup[(scope, baseline, "Recall@50")]
            selected_rows.append({"panel": "interval_width_ratio", **row})
        for item in ("query_count", "source_document_count", "component_count", "component_query_size_max", "largest_component_query_fraction"):
            selected_rows.append(
                {
                    "panel": "scope_provenance",
                    "record_type": "document_component_scope",
                    "scope_or_cohort": scope,
                    "baseline_or_left": "",
                    "right_baseline": "",
                    "metric_or_item": item,
                    "point_estimate": scope_lookup[(scope, item)],
                    "primary_ci95_low": "",
                    "primary_ci95_high": "",
                    "component_ci95_low": "",
                    "component_ci95_high": "",
                    "status": "descriptive_provenance_audit",
                }
            )
    for row in difference_rows:
        selected_rows.append(
            {
                "panel": "top100_exhaustive",
                "record_type": row["record_type"],
                "scope_or_cohort": row["scope"],
                "baseline_or_left": "structure_sequence_pair_neighbor",
                "right_baseline": "",
                "metric_or_item": row["metric_or_item"],
                "point_estimate": row["exhaustive_minus_top100"],
                "primary_ci95_low": "",
                "primary_ci95_high": "",
                "component_ci95_low": "",
                "component_ci95_high": "",
                "status": row["status"],
            }
        )
    source_path = figures / "Figure_3_source_data.tsv"
    fields = ["panel", "record_type", "scope_or_cohort", "baseline_or_left", "right_baseline", "metric_or_item", "point_estimate", "primary_ci95_low", "primary_ci95_high", "component_ci95_low", "component_ci95_high", "status"]
    write_tsv(source_path, fields, selected_rows)

    set_style()
    fig = plt.figure(figsize=(11.4, 7.8), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0])

    for column, scope in enumerate(("temporal_strict_ab", "scaffold_cold_strict_ab")):
        ax = fig.add_subplot(grid[0, column])
        y = np.arange(len(METRICS))
        offset = 0.10
        for yi, metric in enumerate(METRICS):
            row = focus_lookup[(scope, metric)]
            point = float(row["point_estimate"])
            q_low, q_high = float(row["primary_ci95_low"]), float(row["primary_ci95_high"])
            c_low, c_high = float(row["component_ci95_low"]), float(row["component_ci95_high"])
            ax.errorbar(point, yi - offset, xerr=np.array([[point - q_low], [q_high - point]]), fmt="o", ms=4.5, capsize=2.5, color="#277DA1", label="Query bootstrap" if yi == 0 else "")
            ax.errorbar(point, yi + offset, xerr=np.array([[point - c_low], [c_high - point]]), fmt="s", ms=4.2, capsize=2.5, color="#F4A261", label="PMID-component sensitivity" if yi == 0 else "")
        ax.axvline(0, color="#404040", lw=0.9, ls="--")
        ax.set_yticks(y, METRICS)
        ax.invert_yaxis()
        ax.set_xlabel("Pair neighbour minus Morgan transfer")
        largest = int(float(scope_lookup[(scope, "component_query_size_max")]))
        queries = int(float(scope_lookup[(scope, "query_count")]))
        pmids = int(float(scope_lookup[(scope, "source_document_count")]))
        components = int(float(scope_lookup[(scope, "component_count")]))
        ax.set_title(f"{SCOPE_LABELS[scope]}\n{pmids} PMIDs; {components} components; largest {largest}/{queries} queries", loc="left", fontweight="bold")
        ax.grid(axis="x", color="#E1E6EA", lw=0.7)
        panel_label(ax, "A" if column == 0 else "B")

    ax = fig.add_subplot(grid[1, 0])
    x = np.arange(len(BASELINES))
    width = 0.34
    for offset_index, scope in enumerate(("temporal_strict_ab", "scaffold_cold_strict_ab")):
        ratios = []
        for baseline in BASELINES:
            row = component_lookup[(scope, baseline, "Recall@50")]
            primary_width = float(row["primary_ci95_high"]) - float(row["primary_ci95_low"])
            component_width = float(row["component_ci95_high"]) - float(row["component_ci95_low"])
            ratios.append(component_width / primary_width if primary_width > 0 else np.nan)
        positions = x + (offset_index - 0.5) * width
        ax.bar(positions, ratios, width=width, color=("#277DA1" if offset_index == 0 else "#F4A261"), label=SCOPE_LABELS[scope])
    ax.axhline(1, color="#404040", lw=0.9, ls="--")
    ax.set_xticks(x, [BASELINE_SHORT[item] for item in BASELINES], rotation=18, ha="right")
    ax.set_ylabel("Component-CI width / query-CI width")
    ax.set_title("Recall@50 interval widening under PMID components", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=7.5)
    ax.grid(axis="y", color="#E1E6EA", lw=0.7)
    panel_label(ax, "C")

    ax = fig.add_subplot(grid[1, 1])
    delta_matrix = np.array([[float(difference_lookup[(scope, metric)]["exhaustive_minus_top100"]) for metric in METRICS] for scope in SCOPES])
    max_abs = max(float(np.max(np.abs(delta_matrix))), 1e-6)
    image = ax.imshow(delta_matrix, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-max_abs, vcenter=0, vmax=max_abs))
    for row_index, values in enumerate(delta_matrix):
        for column_index, value in enumerate(values):
            label = "0" if value == 0 else f"{value:+.1e}"
            text_color = "white" if abs(value) > 0.58 * max_abs else "#1D2933"
            ax.text(column_index, row_index, label, ha="center", va="center", fontsize=7.0, color=text_color)
    ax.set_xticks(np.arange(len(METRICS)), METRICS, rotation=22, ha="right")
    ax.set_yticks(np.arange(len(SCOPES)), [SCOPE_LABELS[item] for item in SCOPES])
    ax.tick_params(length=0)
    ax.set_title("Exhaustive minus frozen top-100 pair-neighbour metrics", loc="left", fontweight="bold")
    colorbar = fig.colorbar(image, ax=ax, fraction=0.045, pad=0.02)
    colorbar.set_label("Metric difference")
    panel_label(ax, "D")

    legend_handles = [
        Line2D([0], [0], marker="o", color="#277DA1", lw=1.4, markersize=5, label="Query bootstrap"),
        Line2D([0], [0], marker="s", color="#F4A261", lw=1.4, markersize=5, label="PMID-component sensitivity"),
    ]
    fig.legend(handles=legend_handles, frameon=False, ncol=2, loc="lower center", bbox_to_anchor=(0.5, -0.015), fontsize=8)
    output = figures / "Figure_3_dependence_and_fidelity_audits.png"
    save_rgb_png(fig, output)
    return output, source_path


def main() -> int:
    project = Path(__file__).resolve().parents[2]
    manuscript = project / "manuscript_molecular_diversity_v3_20260728"
    tables = manuscript / "tables"
    figures = manuscript / "figures"
    figures.mkdir(parents=True, exist_ok=True)

    figure1_path, figure1_source, figure1_upstream = figure1(project, figures)
    figure2_path, figure2_source = figure2(tables, figures)
    figure3_path, figure3_source = figure3(tables, figures)
    manifest_path = figures / "corrective_figure_generation_manifest_v3.json"
    script_hash = sha256_file(Path(__file__))

    def table_source(path: Path, role: str) -> dict[str, str]:
        return {"path": path.relative_to(project).as_posix(), "sha256": sha256_file(path), "role": role}

    payload = {
        "schema_version": "1.1",
        "format": "PNG_only",
        "identifier_bearing_outputs": False,
        "figures": [
            {
                "figure": figure1_path.name,
                "sha256": sha256_file(figure1_path),
                "bytes": figure1_path.stat().st_size,
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
                "source_data": figure2_source.name,
                "source_data_sha256": sha256_file(figure2_source),
                "generation_script": Path(__file__).name,
                "generation_script_sha256": script_hash,
                "upstream_artifacts": [
                    table_source(tables / "Table_3_corrected_aggregate_performance.tsv", "corrected_metrics"),
                    table_source(tables / "Table_4_corrected_bootstrap_summaries.tsv", "primary_query_bootstrap"),
                    table_source(tables / "Table_S3_zero_and_failure_accounting.tsv", "post_hoc_endpoint_score_attribution"),
                ],
            },
            {
                "figure": figure3_path.name,
                "sha256": sha256_file(figure3_path),
                "bytes": figure3_path.stat().st_size,
                "source_data": figure3_source.name,
                "source_data_sha256": sha256_file(figure3_source),
                "generation_script": Path(__file__).name,
                "generation_script_sha256": script_hash,
                "upstream_artifacts": [
                    table_source(tables / "Table_S2_top100_exhaustive_fidelity.tsv", "approximation_fidelity"),
                    table_source(tables / "Table_S4_pmid_document_dependence.tsv", "document_dependence"),
                ],
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
