"""Evaluate corrective ranks against the frozen strict A/B endpoint and scopes."""

from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import json
import math
import sys
import time
import tracemalloc
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
sys.path.insert(0, str(WORKSPACE / "scripts"))

from pu_retrieval_metrics import macro_average, query_metrics

from asof_successor_common import (
    BASELINES,
    BOOTSTRAP_REPLICATES,
    BOOTSTRAP_SEED,
    ENDPOINT_DECISION,
    EXPECTED_CANDIDATE_TARGETS,
    EXPECTED_COMPLETE_RANK_ROWS,
    EXPECTED_ENDPOINT_RELATIONS,
    EXPECTED_ENDPOINT_TARGETS,
    EXPECTED_QUERIES,
    FOCUS_LEFT,
    FOCUS_RIGHT,
    LEGACY_TIE_SALT,
    METRICS,
    PROTOCOL_ID,
    RUN_ID,
    RUN_MODE,
    SCOPES,
    STRICT_TIERS,
    assert_isolated_input,
    assert_new_output_dir,
    code_hashes,
    environment_receipt,
    load_input_manifest,
    load_receipt,
    parse_bool,
    peak_rss_bytes,
    read_tsv_gz,
    require,
    require_fields,
    require_unique,
    sha256,
    write_json,
    write_tsv_gz,
)


EVALUATION_INPUT_KIND = "corrective_evaluation_endpoint"
FROZEN_ENDPOINT_SHA256 = "09296b066a23197a7c178f00514f2b3d9ed7e6f3c459ea92a55e01a6010d1132"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    return parser


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--receipt", required=True, type=Path)
    parser.add_argument("--input-manifest", required=True, type=Path)
    parser.add_argument("--score-dir", required=True, type=Path)
    parser.add_argument("--endpoint", required=True, type=Path)
    parser.add_argument("--scaffold-audit", required=True, type=Path)
    parser.add_argument("--homology-0-30", required=True, type=Path)
    parser.add_argument("--homology-0-50", required=True, type=Path)
    parser.add_argument("--homology-0-70", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)


def read_bool_map(path: Path, key: str, flag: str, status: str, label: str) -> dict[str, bool]:
    fields, rows = read_tsv_gz(path)
    require_fields(fields, {key, flag, status}, label)
    require_unique(rows, (key,), label)
    result: dict[str, bool] = {}
    for row in rows:
        result[row[key]] = parse_bool(row[flag], f"{label} {row[key]}")
        require(bool(row[status].strip()), f"{label} has an empty status")
    return result


def load_endpoint(path: Path) -> list[dict[str, str]]:
    require(sha256(path) == FROZEN_ENDPOINT_SHA256, "Frozen endpoint byte hash changed")
    fields, rows = read_tsv_gz(path)
    require_fields(
        fields,
        {
            "canonical_pair_key",
            "query_id",
            "inchikey_full",
            "uniprot_canonical_accession",
            "best_strict_evidence_tier",
            "decision",
            "c31_leakage_gate_status",
        },
        "frozen corrective endpoint",
    )
    require_unique(rows, ("canonical_pair_key",), "frozen corrective endpoint")
    require_unique(rows, ("inchikey_full", "uniprot_canonical_accession"), "frozen corrective endpoint")
    require(len(rows) == EXPECTED_ENDPOINT_RELATIONS, "Frozen endpoint relation count changed")
    require(
        len({row["query_id"] for row in rows}) == EXPECTED_QUERIES,
        "Frozen endpoint query count changed",
    )
    require(
        len({row["uniprot_canonical_accession"] for row in rows}) == EXPECTED_ENDPOINT_TARGETS,
        "Frozen endpoint target count changed",
    )
    query_compounds: dict[str, str] = {}
    for row in rows:
        require(row["best_strict_evidence_tier"] in STRICT_TIERS, "Endpoint has a non-strict tier")
        require(row["decision"] == ENDPOINT_DECISION, "Endpoint has an invalid decision")
        require(
            row["c31_leakage_gate_status"] == "pass_no_historical_activity",
            "Endpoint fails the normalized C31 leakage gate",
        )
        query_id = row["query_id"]
        compound = row["inchikey_full"]
        if query_id in query_compounds:
            require(query_compounds[query_id] == compound, "Endpoint query ID maps to multiple compounds")
        query_compounds[query_id] = compound
    return rows


def verify_recorded_code_hashes(score_manifest: dict[str, Any]) -> None:
    code = score_manifest.get("code")
    require(isinstance(code, dict) and code, "Score manifest lacks executable code hashes")
    for label, item in code.items():
        require(isinstance(item, dict), f"Invalid code receipt for {label}")
        path = Path(str(item.get("path", "")))
        require(path.is_file(), f"Recorded score code is absent: {path}")
        require(item.get("sha256") == sha256(path), f"Recorded score code hash changed: {label}")


def load_prediction_ranks(
    path: Path,
    expected_row_count: int,
    needed_targets: dict[str, set[str]],
) -> tuple[
    dict[str, dict[str, dict[str, int]]],
    dict[tuple[str, str], dict[str, int]],
    dict[str, str],
]:
    selected: dict[str, dict[str, dict[str, int]]] = {
        baseline: defaultdict(dict) for baseline in BASELINES
    }
    audit: dict[tuple[str, str], dict[str, int]] = defaultdict(
        lambda: {
            "count": 0,
            "rank_sum": 0,
            "rank_square_sum": 0,
            "rank_min": 2**31 - 1,
            "rank_max": 0,
            "candidate_count": -1,
        }
    )
    query_compounds: dict[str, str] = {}
    seen_rows = 0
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Corrective prediction ranks lack a header")
        require_fields(
            reader.fieldnames,
            {
                "protocol_id",
                "baseline",
                "query_id",
                "query_compound_inchikey_full",
                "target_uniprot_accession",
                "rank",
                "score",
                "eligible_candidate_target_count",
            },
            "corrective prediction ranks",
        )
        for row in reader:
            seen_rows += 1
            require(row["protocol_id"] == PROTOCOL_ID, "Prediction rank protocol ID mismatch")
            require(row["baseline"] in BASELINES, "Prediction rank has an unknown baseline")
            query_id = row["query_id"]
            compound = row["query_compound_inchikey_full"]
            if query_id in query_compounds:
                require(query_compounds[query_id] == compound, "Rank query ID maps to multiple compounds")
            query_compounds[query_id] = compound
            rank = int(row["rank"])
            score = float(row["score"])
            require(math.isfinite(score), "Corrective rank contains a non-finite score")
            candidate_count = int(row["eligible_candidate_target_count"])
            require(1 <= rank <= candidate_count, "Prediction rank is outside its candidate range")
            key = (row["baseline"], query_id)
            item = audit[key]
            require(
                item["candidate_count"] in {-1, candidate_count},
                "Candidate count changes within a baseline/query block",
            )
            item["candidate_count"] = candidate_count
            item["count"] += 1
            item["rank_sum"] += rank
            item["rank_square_sum"] += rank * rank
            item["rank_min"] = min(item["rank_min"], rank)
            item["rank_max"] = max(item["rank_max"], rank)
            target = row["target_uniprot_accession"]
            if target in needed_targets.get(query_id, set()):
                require(
                    target not in selected[row["baseline"]][query_id],
                    "Endpoint target occurs more than once in a prediction block",
                )
                selected[row["baseline"]][query_id][target] = rank
    require(seen_rows == expected_row_count, "Corrective prediction row count mismatch")
    require(len(query_compounds) == EXPECTED_QUERIES, "Corrective prediction query count mismatch")
    expected_blocks = {
        (baseline, query_id) for baseline in BASELINES for query_id in query_compounds
    }
    require(set(audit) == expected_blocks, "Corrective prediction baseline/query blocks are incomplete")
    for key, item in audit.items():
        count = item["count"]
        expected_sum = count * (count + 1) // 2
        expected_square_sum = count * (count + 1) * (2 * count + 1) // 6
        require(
            count == item["candidate_count"]
            and item["rank_min"] == 1
            and item["rank_max"] == count
            and item["rank_sum"] == expected_sum
            and item["rank_square_sum"] == expected_square_sum,
            f"Corrective rank permutation failed for {key}",
        )
    return selected, audit, query_compounds


def build_scope_relevance(
    endpoint: list[dict[str, str]],
    scaffold: dict[str, bool],
    homology: dict[str, dict[str, bool]],
) -> dict[str, dict[str, list[dict[str, str]]]]:
    endpoint_keys = {row["canonical_pair_key"] for row in endpoint}
    require(set(scaffold) == endpoint_keys, "Scaffold audit keyset differs from endpoint")
    endpoint_targets = {row["uniprot_canonical_accession"] for row in endpoint}
    for threshold, flags in homology.items():
        require(set(flags) == endpoint_targets, f"Homology keyset differs from endpoint at {threshold}")
    result: dict[str, dict[str, list[dict[str, str]]]] = {
        scope: defaultdict(list) for scope in SCOPES
    }
    for row in endpoint:
        query_id = row["query_id"]
        target = row["uniprot_canonical_accession"]
        result["temporal_strict_ab"][query_id].append(row)
        if scaffold[row["canonical_pair_key"]]:
            result["scaffold_cold_strict_ab"][query_id].append(row)
            for threshold in ("0_30", "0_50", "0_70"):
                if homology[threshold][target]:
                    result[f"double_cold_{threshold}"][query_id].append(row)
    return result


def bootstrap_baselines(
    arrays: dict[str, np.ndarray],
    scope: str,
    indices: np.ndarray | None,
) -> list[dict[str, object]]:
    n_queries = next(iter(arrays.values())).shape[0] if arrays else 0
    rows: list[dict[str, object]] = []
    for baseline in BASELINES:
        for metric_index, metric in enumerate(METRICS):
            if n_queries == 0:
                rows.append(
                    {
                        "scope": scope,
                        "baseline": baseline,
                        "metric": metric,
                        "query_count": 0,
                        "mean": "",
                        "ci95_low": "",
                        "ci95_high": "",
                        "status": "not_estimable",
                    }
                )
                continue
            observed = float(arrays[baseline][:, metric_index].mean())
            if n_queries == 1:
                low = high = observed
                status = "n=1_descriptive_only"
            else:
                require(indices is not None, "Bootstrap indices are absent for an estimable scope")
                values = arrays[baseline][indices, metric_index].mean(axis=1)
                low, high = (float(item) for item in np.percentile(values, [2.5, 97.5]))
                status = "estimable_descriptive"
            rows.append(
                {
                    "scope": scope,
                    "baseline": baseline,
                    "metric": metric,
                    "query_count": n_queries,
                    "mean": f"{observed:.17g}",
                    "ci95_low": f"{low:.17g}",
                    "ci95_high": f"{high:.17g}",
                    "status": status,
                }
            )
    return rows


def ordered_baseline_pairs() -> list[tuple[str, str, str]]:
    pairs: list[tuple[str, str, str]] = []
    for first, second in itertools.combinations(BASELINES, 2):
        if {first, second} == {FOCUS_LEFT, FOCUS_RIGHT}:
            pairs.append((FOCUS_LEFT, FOCUS_RIGHT, "prespecified_focus"))
        else:
            pairs.append((first, second, "supplementary_unordered_pair"))
    return pairs


def bootstrap_pairs(
    arrays: dict[str, np.ndarray],
    scope: str,
    indices: np.ndarray | None,
) -> list[dict[str, object]]:
    n_queries = next(iter(arrays.values())).shape[0] if arrays else 0
    rows: list[dict[str, object]] = []
    for left, right, role in ordered_baseline_pairs():
        for metric_index, metric in enumerate(METRICS):
            if n_queries == 0:
                rows.append(
                    {
                        "scope": scope,
                        "left_baseline": left,
                        "right_baseline": right,
                        "comparison_role": role,
                        "metric": metric,
                        "query_count": 0,
                        "mean_difference_left_minus_right": "",
                        "ci95_low": "",
                        "ci95_high": "",
                        "status": "not_estimable",
                    }
                )
                continue
            difference = arrays[left][:, metric_index] - arrays[right][:, metric_index]
            observed = float(difference.mean())
            if n_queries == 1:
                low = high = observed
                status = "n=1_descriptive_only"
            else:
                require(indices is not None, "Paired-bootstrap indices are absent")
                values = difference[indices].mean(axis=1)
                low, high = (float(item) for item in np.percentile(values, [2.5, 97.5]))
                status = "estimable_descriptive"
            rows.append(
                {
                    "scope": scope,
                    "left_baseline": left,
                    "right_baseline": right,
                    "comparison_role": role,
                    "metric": metric,
                    "query_count": n_queries,
                    "mean_difference_left_minus_right": f"{observed:.17g}",
                    "ci95_low": f"{low:.17g}",
                    "ci95_high": f"{high:.17g}",
                    "status": status,
                }
            )
    return rows


def run(args: argparse.Namespace) -> int:
    total_started = time.perf_counter()
    tracemalloc.start()
    phase_seconds: dict[str, float] = {}

    validation_started = time.perf_counter()
    receipt_path = assert_isolated_input(args.receipt)
    receipt = load_receipt(receipt_path)
    input_paths = {
        "endpoint": assert_isolated_input(args.endpoint),
        "scaffold_audit": assert_isolated_input(args.scaffold_audit),
        "homology_0_30": assert_isolated_input(args.homology_0_30),
        "homology_0_50": assert_isolated_input(args.homology_0_50),
        "homology_0_70": assert_isolated_input(args.homology_0_70),
    }
    input_manifest_path = assert_isolated_input(args.input_manifest)
    evaluation_manifest = load_input_manifest(
        input_manifest_path, input_paths, EVALUATION_INPUT_KIND
    )
    score_dir = assert_isolated_input(args.score_dir)
    output_dir = assert_new_output_dir(args.output_dir)
    score_manifest_path = score_dir / "corrective_score_manifest.json"
    rank_path = score_dir / "corrective_prediction_ranks.tsv.gz"
    require(score_manifest_path.is_file() and rank_path.is_file(), "Corrective score directory is incomplete")
    score_manifest = json.loads(score_manifest_path.read_text(encoding="utf-8"))
    require(score_manifest.get("protocol_id") == PROTOCOL_ID, "Score manifest protocol ID mismatch")
    require(score_manifest.get("run_id") == RUN_ID, "Score manifest run ID mismatch")
    require(score_manifest.get("stage") == "corrective_score", "Score manifest stage mismatch")
    require(score_manifest.get("execution_mode") == RUN_MODE, "Score manifest execution mode mismatch")
    require(score_manifest.get("endpoint_file_supplied_to_score_command") is False, "Score received endpoint")
    require(score_manifest.get("endpoint_read_by_score_engine") is False, "Score engine read endpoint")
    require(score_manifest.get("legacy_outer_or_result_read") is False, "Score read a legacy result")
    require(score_manifest.get("target_count") == EXPECTED_CANDIDATE_TARGETS, "Score target count changed")
    require(score_manifest.get("query_count") == EXPECTED_QUERIES, "Score query count changed")
    require(score_manifest.get("row_count") == EXPECTED_COMPLETE_RANK_ROWS, "Score row count changed")
    require(score_manifest.get("baselines") == BASELINES, "Score baseline order or membership changed")
    require(
        score_manifest.get("parameters", {}).get("tie_salt") == LEGACY_TIE_SALT,
        "Score tie salt changed",
    )
    require(
        score_manifest.get("parameters", {}).get("pair_neighbor_sequence_top_k") == 100,
        "Pair-neighbor top-k changed",
    )
    require(
        score_manifest.get("parameters", {}).get("morgan_radius") == 2
        and score_manifest.get("parameters", {}).get("morgan_bits") == 2048,
        "Morgan fingerprint parameters changed",
    )
    require(
        score_manifest.get("parameters", {}).get("mask_historical_targets_for_same_query") is True,
        "Historical same-query target masking changed",
    )
    role = score_manifest.get("role_separation", {})
    require(role.get("historical_and_query_files_distinct") is True, "Score structures were not role-separated")
    require(role.get("pooled_structure_map_used") is False, "Score used a pooled structure map")
    require(score_manifest.get("rank_output", {}).get("sha256") == sha256(rank_path), "Rank hash mismatch")
    verify_recorded_code_hashes(score_manifest)

    endpoint = load_endpoint(input_paths["endpoint"])
    scaffold = read_bool_map(
        input_paths["scaffold_audit"],
        "canonical_pair_key",
        "audit_scaffold_cold_under_selected_policy",
        "audit_outcome",
        "corrective scaffold audit",
    )
    homology = {
        "0_30": read_bool_map(
            input_paths["homology_0_30"],
            "uniprot_canonical_accession",
            "is_future_target_homology_cold_candidate",
            "future_target_coldness_status",
            "corrective homology 0.30",
        ),
        "0_50": read_bool_map(
            input_paths["homology_0_50"],
            "uniprot_canonical_accession",
            "is_future_target_homology_cold_candidate",
            "future_target_coldness_status",
            "corrective homology 0.50",
        ),
        "0_70": read_bool_map(
            input_paths["homology_0_70"],
            "uniprot_canonical_accession",
            "is_future_target_homology_cold_candidate",
            "future_target_coldness_status",
            "corrective homology 0.70",
        ),
    }
    scoped = build_scope_relevance(endpoint, scaffold, homology)
    needed_targets: dict[str, set[str]] = defaultdict(set)
    for row in endpoint:
        needed_targets[row["query_id"]].add(row["uniprot_canonical_accession"])
    ranks, rank_audit, rank_query_compounds = load_prediction_ranks(
        rank_path, EXPECTED_COMPLETE_RANK_ROWS, needed_targets
    )
    for row in endpoint:
        require(
            rank_query_compounds.get(row["query_id"]) == row["inchikey_full"],
            f"Endpoint/query compound mismatch for {row['query_id']}",
        )
    phase_seconds["input_load_validation_and_rank_integrity"] = time.perf_counter() - validation_started

    scope_audit_rows: list[dict[str, object]] = []
    for scope in SCOPES:
        scope_rows = [row for rows in scoped[scope].values() for row in rows]
        scope_audit_rows.append(
            {
                "scope": scope,
                "candidate_relation_count": len(scope_rows),
                "query_compound_count": len(scoped[scope]),
                "target_count": len({row["uniprot_canonical_accession"] for row in scope_rows}),
                "A_affinity_candidate_count": sum(
                    row["best_strict_evidence_tier"] == "A_affinity_candidate" for row in scope_rows
                ),
                "B_quantitative_functional_candidate_count": sum(
                    row["best_strict_evidence_tier"]
                    == "B_quantitative_functional_candidate"
                    for row in scope_rows
                ),
                "status": "estimable" if scoped[scope] else "not_estimable",
            }
        )

    metric_rows: list[dict[str, object]] = []
    baseline_bootstrap_rows: list[dict[str, object]] = []
    contrast_rows: list[dict[str, object]] = []
    generator = np.random.Generator(np.random.PCG64(BOOTSTRAP_SEED))
    scope_seconds: dict[str, float] = {}
    for scope in SCOPES:
        scope_started = time.perf_counter()
        relevant = scoped[scope]
        queries = sorted(relevant)
        arrays: dict[str, np.ndarray] = {}
        for baseline in BASELINES:
            per_query: list[dict[str, float]] = []
            candidate_counts: list[int] = []
            for query_id in queries:
                positive_ranks: list[int] = []
                for endpoint_row in relevant[query_id]:
                    target = endpoint_row["uniprot_canonical_accession"]
                    try:
                        positive_ranks.append(ranks[baseline][query_id][target])
                    except KeyError as exc:
                        raise ValueError(
                            f"Endpoint target lacks a corrective rank for {baseline}, {query_id}, {target}"
                        ) from exc
                per_query.append(query_metrics(positive_ranks, (10, 50)))
                candidate_counts.append(rank_audit[(baseline, query_id)]["candidate_count"])
            arrays[baseline] = (
                np.asarray([[row[metric] for metric in METRICS] for row in per_query], dtype=float)
                if per_query
                else np.empty((0, len(METRICS)), dtype=float)
            )
            if per_query:
                summary = macro_average(per_query)
                candidate_values = np.asarray(candidate_counts, dtype=float)
                metric_rows.append(
                    {
                        "scope": scope,
                        "baseline": baseline,
                        "evaluable_query_count": len(per_query),
                        "candidate_relation_count": sum(len(relevant[item]) for item in queries),
                        "candidate_target_median": f"{np.median(candidate_values):.17g}",
                        "candidate_target_iqr": f"{(np.percentile(candidate_values, 75) - np.percentile(candidate_values, 25)):.17g}",
                        "candidate_target_min": int(candidate_values.min()),
                        "candidate_target_max": int(candidate_values.max()),
                        **{metric: f"{summary[metric]:.17g}" for metric in METRICS},
                        "zero_recall_at_10_queries": sum(row["Recall@10"] == 0.0 for row in per_query),
                        "zero_recall_at_50_queries": sum(row["Recall@50"] == 0.0 for row in per_query),
                        "zero_mrr_queries": sum(row["MRR"] == 0.0 for row in per_query),
                        "status": "estimable" if len(per_query) > 1 else "n=1_descriptive_only",
                    }
                )
            else:
                metric_rows.append(
                    {
                        "scope": scope,
                        "baseline": baseline,
                        "evaluable_query_count": 0,
                        "candidate_relation_count": 0,
                        "candidate_target_median": "",
                        "candidate_target_iqr": "",
                        "candidate_target_min": "",
                        "candidate_target_max": "",
                        **{metric: "" for metric in METRICS},
                        "zero_recall_at_10_queries": 0,
                        "zero_recall_at_50_queries": 0,
                        "zero_mrr_queries": 0,
                        "status": "not_estimable",
                    }
                )
        n_queries = len(queries)
        indices = (
            generator.integers(
                0,
                n_queries,
                size=(BOOTSTRAP_REPLICATES, n_queries),
                endpoint=False,
            )
            if n_queries > 1
            else None
        )
        baseline_bootstrap_rows.extend(bootstrap_baselines(arrays, scope, indices))
        contrast_rows.extend(bootstrap_pairs(arrays, scope, indices))
        scope_seconds[scope] = time.perf_counter() - scope_started

    focus_rows = [
        row for row in contrast_rows if row["comparison_role"] == "prespecified_focus"
    ]
    require(len(metric_rows) == len(SCOPES) * len(BASELINES), "Aggregate metric matrix incomplete")
    require(
        len(baseline_bootstrap_rows) == len(SCOPES) * len(BASELINES) * len(METRICS),
        "Baseline-bootstrap matrix incomplete",
    )
    require(
        len(contrast_rows) == len(SCOPES) * 6 * len(METRICS),
        "Paired-contrast matrix incomplete",
    )
    require(len(focus_rows) == len(SCOPES) * len(METRICS), "Focus-contrast matrix incomplete")
    require(
        all(
            row["left_baseline"] == FOCUS_LEFT and row["right_baseline"] == FOCUS_RIGHT
            for row in focus_rows
        ),
        "Focus contrast direction is not pair-neighbor minus Morgan",
    )
    phase_seconds["metric_and_bootstrap_evaluation"] = sum(scope_seconds.values())

    output_dir.mkdir(parents=True, exist_ok=False)
    metric_path = output_dir / "corrective_aggregate_metrics.tsv.gz"
    scope_path = output_dir / "corrective_scope_denominator_audit.tsv.gz"
    baseline_ci_path = output_dir / "corrective_baseline_bootstrap_metrics.tsv.gz"
    contrast_path = output_dir / "corrective_paired_bootstrap_contrasts.tsv.gz"
    focus_path = output_dir / "corrective_focus_pair_neighbor_minus_morgan.tsv.gz"
    write_tsv_gz(metric_path, list(metric_rows[0]), metric_rows)
    write_tsv_gz(scope_path, list(scope_audit_rows[0]), scope_audit_rows)
    write_tsv_gz(baseline_ci_path, list(baseline_bootstrap_rows[0]), baseline_bootstrap_rows)
    write_tsv_gz(contrast_path, list(contrast_rows[0]), contrast_rows)
    write_tsv_gz(focus_path, list(focus_rows[0]), focus_rows)

    python_current, python_peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    peak_rss, peak_rss_method = peak_rss_bytes()
    phase_seconds["total_before_manifest_write"] = time.perf_counter() - total_started
    manifest_path = output_dir / "corrective_evaluation_manifest.json"
    shared_metric_path = WORKSPACE / "scripts" / "pu_retrieval_metrics.py"
    executed_entrypoint = Path(sys.argv[0]).resolve()
    executed_code = [Path(__file__), ROOT / "scripts" / "asof_successor_common.py", shared_metric_path]
    if executed_entrypoint.is_file() and executed_entrypoint not in {path.resolve() for path in executed_code}:
        executed_code.append(executed_entrypoint)
    outputs = {
        "aggregate_metrics": metric_path,
        "scope_denominator_audit": scope_path,
        "baseline_bootstrap": baseline_ci_path,
        "paired_bootstrap": contrast_path,
        "focus_pair_neighbor_minus_morgan": focus_path,
    }
    write_json(
        manifest_path,
        {
            "schema_version": "1.0",
            "protocol_id": PROTOCOL_ID,
            "run_id": evaluation_manifest.get("run_id"),
            "stage": "corrective_evaluation",
            "execution_mode": RUN_MODE,
            "claim_boundary": "Author-run non-independent post hoc temporal-purity correction; no blinded or external-validation claim.",
            "command_argv": list(sys.argv),
            "receipt": {"path": str(receipt_path), "sha256": sha256(receipt_path)},
            "evaluation_input_manifest": {
                "path": str(input_manifest_path),
                "sha256": sha256(input_manifest_path),
            },
            "score_manifest": {"path": str(score_manifest_path), "sha256": sha256(score_manifest_path)},
            "score_rank": {"path": str(rank_path), "sha256": sha256(rank_path)},
            "evaluation_inputs": {
                name: {"path": str(path), "sha256": sha256(path)}
                for name, path in input_paths.items()
            },
            "outputs": {
                name: {"path": str(path), "sha256": sha256(path)} for name, path in outputs.items()
            },
            "all_scopes_reported": SCOPES,
            "all_baselines_reported": BASELINES,
            "all_metrics_reported": METRICS,
            "focus_contrast": {
                "left_baseline": FOCUS_LEFT,
                "right_baseline": FOCUS_RIGHT,
                "estimand": "left_minus_right",
                "row_count": len(focus_rows),
            },
            "bootstrap": {
                "unit": "query_compound",
                "replicates": BOOTSTRAP_REPLICATES,
                "prng": "PCG64",
                "seed": BOOTSTRAP_SEED,
                "interval": "95% percentile",
                "n_zero_rule": "emit_not_estimable_for_every_planned_cell",
            },
            "runtime": {
                "wall_seconds_by_phase": phase_seconds,
                "wall_seconds_by_scope": scope_seconds,
                "python_tracemalloc_current_bytes": python_current,
                "python_tracemalloc_peak_bytes": python_peak,
                "process_peak_rss_bytes": peak_rss,
                "process_peak_rss_method": peak_rss_method,
            },
            "environment": environment_receipt(),
            "code": code_hashes(executed_code),
            "endpoint_hash_frozen": True,
            "endpoint_access_control_claim": "none",
            "figures_generated": False,
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
        },
    )
    print(
        json.dumps(
            {
                "output_dir": str(output_dir),
                "aggregate_rows": len(metric_rows),
                "paired_contrast_rows": len(contrast_rows),
                "focus_contrast_rows": len(focus_rows),
                "figures_generated": False,
            }
        )
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
