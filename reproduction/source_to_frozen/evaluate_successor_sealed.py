"""Score sealed strict A/B successor endpoints after blind predictions are locked."""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import itertools
import json
import platform
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
sys.path.insert(0, str(WORKSPACE / "scripts"))

from pu_baseline_io import read_tsv_gz, sha256 as legacy_sha256
from pu_retrieval_metrics import macro_average, query_metrics
from successor_common import (
    ENDPOINT_DECISION,
    PROTOCOL_ID,
    STRICT_TIERS,
    assert_isolated_input,
    input_manifest,
    load_authorized_manifest,
    load_receipt,
    parse_bool,
    require_fields,
    require_unique,
)


BASELINES = [
    "weighted_target_popularity",
    "sequence_3mer_transfer",
    "weighted_morgan_transfer",
    "structure_sequence_pair_neighbor",
]
SCOPES = [
    "temporal_strict_ab",
    "scaffold_cold_strict_ab",
    "double_cold_0_30",
    "double_cold_0_50",
    "double_cold_0_70",
]
METRICS = ["Recall@10", "Recall@50", "NDCG@10", "NDCG@50", "MRR"]


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--receipt", required=True, type=Path)
    result.add_argument("--sealed-input-manifest", required=True, type=Path)
    result.add_argument("--blind-prediction-dir", required=True, type=Path)
    result.add_argument("--sealed-endpoint", required=True, type=Path)
    result.add_argument("--sealed-scaffold-audit", required=True, type=Path)
    result.add_argument("--sealed-homology-0-30", required=True, type=Path)
    result.add_argument("--sealed-homology-0-50", required=True, type=Path)
    result.add_argument("--sealed-homology-0-70", required=True, type=Path)
    result.add_argument("--output-dir", required=True, type=Path)
    return result


def write_tsv_gz(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite output: {path}")
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(text, fieldnames=fields, delimiter="\t", extrasaction="raise")
                writer.writeheader()
                writer.writerows(rows)


def read_bool_map(path: Path, key: str, flag: str, status: str, label: str) -> dict[str, bool]:
    fields, rows = read_tsv_gz(path)
    require_fields(fields, {key, flag, status}, label)
    require_unique(rows, (key,), label)
    result: dict[str, bool] = {}
    for row in rows:
        value = parse_bool(row[flag], f"{label} {row[key]}")
        if not row[status].strip():
            raise ValueError(f"{label} has an empty status")
        result[row[key]] = value
    return result


def load_endpoint(path: Path) -> list[dict[str, str]]:
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
        "sealed endpoint",
    )
    require_unique(rows, ("canonical_pair_key",), "sealed endpoint")
    require_unique(rows, ("inchikey_full", "uniprot_canonical_accession"), "sealed endpoint")
    for row in rows:
        if row["best_strict_evidence_tier"] not in STRICT_TIERS:
            raise ValueError("Sealed endpoint contains a non-strict evidence tier")
        if row["decision"] != ENDPOINT_DECISION:
            raise ValueError("Sealed endpoint contains an invalid decision")
        if row["c31_leakage_gate_status"] != "pass_no_historical_activity":
            raise ValueError("Sealed endpoint fails the normalized C31 leakage gate")
    return rows


def load_prediction_ranks(
    path: Path,
    expected_row_count: int,
    needed_targets: dict[str, set[str]],
) -> tuple[
    dict[str, dict[str, dict[str, int]]],
    dict[tuple[str, str], dict[str, int]],
    dict[str, str],
]:
    selected: dict[str, dict[str, dict[str, int]]] = {baseline: defaultdict(dict) for baseline in BASELINES}
    audit: dict[tuple[str, str], dict[str, int]] = defaultdict(lambda: {"count": 0, "rank_sum": 0, "rank_min": 2**31 - 1, "rank_max": 0, "candidate_count": -1})
    query_compounds: dict[str, str] = {}
    seen_rows = 0
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError("Blind prediction ranks lack a header")
        require_fields(
            list(reader.fieldnames),
            {
                "protocol_id",
                "baseline",
                "query_id",
                "query_compound_inchikey_full",
                "target_uniprot_accession",
                "rank",
                "eligible_candidate_target_count",
            },
            "blind prediction ranks",
        )
        for row in reader:
            seen_rows += 1
            if row["protocol_id"] != PROTOCOL_ID or row["baseline"] not in BASELINES:
                raise ValueError("Blind prediction rank has an invalid protocol ID or baseline")
            query_id = row["query_id"]
            query_compound = row["query_compound_inchikey_full"]
            if query_id in query_compounds and query_compounds[query_id] != query_compound:
                raise ValueError("Blind prediction query_id maps to multiple compounds")
            query_compounds[query_id] = query_compound
            rank = int(row["rank"])
            candidate_count = int(row["eligible_candidate_target_count"])
            if rank < 1 or rank > candidate_count:
                raise ValueError("Blind prediction rank is outside its candidate range")
            key = (row["baseline"], row["query_id"])
            item = audit[key]
            if item["candidate_count"] not in {-1, candidate_count}:
                raise ValueError("Candidate count changes within a baseline/query rank block")
            item["candidate_count"] = candidate_count
            item["count"] += 1
            item["rank_sum"] += rank
            item["rank_min"] = min(item["rank_min"], rank)
            item["rank_max"] = max(item["rank_max"], rank)
            if row["target_uniprot_accession"] in needed_targets.get(query_id, set()):
                selected[row["baseline"]][query_id][row["target_uniprot_accession"]] = rank
    if seen_rows != expected_row_count:
        raise ValueError(f"Blind prediction row count mismatch: {seen_rows} != {expected_row_count}")
    for key, item in audit.items():
        count = item["count"]
        expected_sum = count * (count + 1) // 2
        if count != item["candidate_count"] or item["rank_min"] != 1 or item["rank_max"] != count or item["rank_sum"] != expected_sum:
            raise ValueError(f"Blind prediction rank permutation failed for {key}")
    return selected, audit, query_compounds


def build_scope_relevance(
    endpoint: list[dict[str, str]],
    scaffold: dict[str, bool],
    homology: dict[str, dict[str, bool]],
) -> dict[str, dict[str, list[dict[str, str]]]]:
    endpoint_keys = {row["canonical_pair_key"] for row in endpoint}
    if set(scaffold) != endpoint_keys:
        raise ValueError("Scaffold audit keyset differs from sealed endpoint keyset")
    endpoint_targets = {row["uniprot_canonical_accession"] for row in endpoint}
    for threshold, flags in homology.items():
        if set(flags) != endpoint_targets:
            raise ValueError(f"Homology target keyset differs from sealed endpoint at {threshold}")

    result: dict[str, dict[str, list[dict[str, str]]]] = {scope: defaultdict(list) for scope in SCOPES}
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
    result: list[dict[str, object]] = []
    n_queries = next(iter(arrays.values())).shape[0] if arrays else 0
    for baseline in BASELINES:
        for metric_index, metric in enumerate(METRICS):
            if n_queries == 0:
                result.append(
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
                if indices is None:
                    raise ValueError("Bootstrap indices are missing for an estimable scope")
                values = arrays[baseline][indices, metric_index].mean(axis=1)
                low, high = (float(item) for item in np.percentile(values, [2.5, 97.5]))
                status = "estimable_descriptive"
            result.append(
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
    return result


def bootstrap_pairs(
    arrays: dict[str, np.ndarray],
    scope: str,
    indices: np.ndarray | None,
) -> list[dict[str, object]]:
    result: list[dict[str, object]] = []
    n_queries = next(iter(arrays.values())).shape[0] if arrays else 0
    if n_queries == 0:
        return result
    for left, right in itertools.combinations(BASELINES, 2):
        difference = arrays[left] - arrays[right]
        for metric_index, metric in enumerate(METRICS):
            observed = float(difference[:, metric_index].mean())
            if n_queries == 1:
                low = high = observed
                status = "n=1_descriptive_only"
            else:
                if indices is None:
                    raise ValueError("Bootstrap indices are missing for an estimable scope")
                values = difference[indices, metric_index].mean(axis=1)
                low, high = (float(item) for item in np.percentile(values, [2.5, 97.5]))
                status = "estimable_descriptive"
            result.append(
                {
                    "scope": scope,
                    "left_baseline": left,
                    "right_baseline": right,
                    "metric": metric,
                    "query_count": n_queries,
                    "mean_difference_left_minus_right": f"{observed:.17g}",
                    "ci95_low": f"{low:.17g}",
                    "ci95_high": f"{high:.17g}",
                    "status": status,
                }
            )
    return result


def main() -> int:
    args = parser().parse_args()
    receipt = load_receipt(args.receipt, "evaluate")
    if receipt.get("actor_role") not in {"independent_evaluator", "author_run_evaluator"}:
        raise ValueError("Evaluation receipt actor_role must be independent_evaluator or author_run_evaluator")
    independence = receipt.get("evaluation_independence")
    if independence not in {"independent", "author_run_non_independent"}:
        raise ValueError("Evaluation receipt must state independent or author_run_non_independent")

    prediction_dir = assert_isolated_input(args.blind_prediction_dir)
    input_paths = {
        "sealed_endpoint": assert_isolated_input(args.sealed_endpoint),
        "sealed_scaffold_audit": assert_isolated_input(args.sealed_scaffold_audit),
        "sealed_homology_0_30": assert_isolated_input(args.sealed_homology_0_30),
        "sealed_homology_0_50": assert_isolated_input(args.sealed_homology_0_50),
        "sealed_homology_0_70": assert_isolated_input(args.sealed_homology_0_70),
    }
    sealed_authorization = load_authorized_manifest(
        assert_isolated_input(args.sealed_input_manifest),
        input_paths,
        {
            "authorized_internal_successor_use": True,
            "sealed_endpoint_included": True,
            "legacy_outer_or_result_input": False,
        },
    )
    output_dir = assert_isolated_input(args.output_dir)
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite successor evaluation directory: {output_dir}")
    prediction_manifest_path = prediction_dir / "blind_prediction_manifest.json"
    rank_path = prediction_dir / "blind_prediction_ranks.tsv.gz"
    if not prediction_manifest_path.is_file() or not rank_path.is_file():
        raise FileNotFoundError("Blind prediction directory is incomplete")
    prediction_manifest = json.loads(prediction_manifest_path.read_text(encoding="utf-8"))
    if prediction_manifest.get("protocol_id") != PROTOCOL_ID or prediction_manifest.get("sealed_endpoint_read") is not False:
        raise ValueError("Blind prediction manifest does not demonstrate endpoint-blind scoring")

    endpoint = load_endpoint(input_paths["sealed_endpoint"])
    scaffold = read_bool_map(
        input_paths["sealed_scaffold_audit"],
        "canonical_pair_key",
        "audit_scaffold_cold_under_selected_policy",
        "audit_outcome",
        "sealed scaffold audit",
    )
    homology = {
        "0_30": read_bool_map(
            input_paths["sealed_homology_0_30"],
            "uniprot_canonical_accession",
            "is_future_target_homology_cold_candidate",
            "future_target_coldness_status",
            "sealed homology 0.30",
        ),
        "0_50": read_bool_map(
            input_paths["sealed_homology_0_50"],
            "uniprot_canonical_accession",
            "is_future_target_homology_cold_candidate",
            "future_target_coldness_status",
            "sealed homology 0.50",
        ),
        "0_70": read_bool_map(
            input_paths["sealed_homology_0_70"],
            "uniprot_canonical_accession",
            "is_future_target_homology_cold_candidate",
            "future_target_coldness_status",
            "sealed homology 0.70",
        ),
    }
    scoped = build_scope_relevance(endpoint, scaffold, homology)
    needed_targets: dict[str, set[str]] = defaultdict(set)
    for endpoint_row in endpoint:
        needed_targets[endpoint_row["query_id"]].add(endpoint_row["uniprot_canonical_accession"])
    ranks, rank_audit, query_compounds = load_prediction_ranks(
        rank_path,
        int(prediction_manifest["row_count"]),
        needed_targets,
    )
    for endpoint_row in endpoint:
        query_id = endpoint_row["query_id"]
        if query_compounds.get(query_id) != endpoint_row["inchikey_full"]:
            raise ValueError(f"Endpoint/query compound mismatch for {query_id}")

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
                    row["best_strict_evidence_tier"] == "B_quantitative_functional_candidate" for row in scope_rows
                ),
                "status": "estimable" if scoped[scope] else "not_estimable",
            }
        )

    metric_rows: list[dict[str, object]] = []
    baseline_bootstrap_rows: list[dict[str, object]] = []
    bootstrap_rows: list[dict[str, object]] = []
    generator = np.random.Generator(np.random.PCG64(20260719))
    for scope in SCOPES:
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
                        rank = ranks[baseline][query_id][target]
                    except KeyError as exc:
                        raise ValueError(f"Endpoint target lacks a blind rank for {baseline}, {query_id}, {target}") from exc
                    positive_ranks.append(rank)
                per_query.append(query_metrics(positive_ranks, (10, 50)))
                candidate_counts.append(rank_audit[(baseline, query_id)]["candidate_count"])
            if per_query:
                summary = macro_average(per_query)
                arrays[baseline] = np.asarray([[row[metric] for metric in METRICS] for row in per_query], dtype=float)
                zero_at_10 = sum(row["Recall@10"] == 0.0 for row in per_query)
                zero_at_50 = sum(row["Recall@50"] == 0.0 for row in per_query)
                zero_mrr = sum(row["MRR"] == 0.0 for row in per_query)
                candidate_values = np.asarray(candidate_counts, dtype=float)
                status = "estimable" if len(per_query) > 1 else "n=1_descriptive_only"
                values = {metric: f"{summary[metric]:.17g}" for metric in METRICS}
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
                        **values,
                        "zero_recall_at_10_queries": zero_at_10,
                        "zero_recall_at_50_queries": zero_at_50,
                        "zero_mrr_queries": zero_mrr,
                        "status": status,
                    }
                )
            else:
                arrays[baseline] = np.empty((0, len(METRICS)), dtype=float)
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
        n_scope_queries = len(queries)
        indices = (
            generator.integers(0, n_scope_queries, size=(10000, n_scope_queries), endpoint=False)
            if n_scope_queries > 1
            else None
        )
        baseline_bootstrap_rows.extend(bootstrap_baselines(arrays, scope, indices))
        bootstrap_rows.extend(bootstrap_pairs(arrays, scope, indices))

    output_dir.mkdir(parents=True, exist_ok=False)
    metric_fields = list(metric_rows[0]) if metric_rows else []
    bootstrap_fields = list(bootstrap_rows[0]) if bootstrap_rows else [
        "scope", "left_baseline", "right_baseline", "metric", "query_count",
        "mean_difference_left_minus_right", "ci95_low", "ci95_high", "status",
    ]
    baseline_bootstrap_fields = list(baseline_bootstrap_rows[0]) if baseline_bootstrap_rows else [
        "scope", "baseline", "metric", "query_count", "mean", "ci95_low", "ci95_high", "status",
    ]
    metric_path = output_dir / "aggregate_metrics.tsv.gz"
    scope_audit_path = output_dir / "scope_denominator_audit.tsv.gz"
    baseline_bootstrap_path = output_dir / "baseline_bootstrap_metrics.tsv.gz"
    bootstrap_path = output_dir / "paired_bootstrap_contrasts.tsv.gz"
    write_tsv_gz(metric_path, metric_fields, metric_rows)
    write_tsv_gz(scope_audit_path, list(scope_audit_rows[0]), scope_audit_rows)
    write_tsv_gz(baseline_bootstrap_path, baseline_bootstrap_fields, baseline_bootstrap_rows)
    write_tsv_gz(bootstrap_path, bootstrap_fields, bootstrap_rows)
    input_paths["blind_prediction_manifest"] = prediction_manifest_path
    input_paths["blind_prediction_ranks"] = rank_path
    evaluation_manifest = input_manifest(input_paths, receipt)
    evaluation_manifest.update(
        {
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "stage": "sealed_evaluation",
            "protocol_id": PROTOCOL_ID,
            "script": str(Path(__file__).resolve()),
            "script_sha256": legacy_sha256(Path(__file__)),
            "sealed_input_manifest": {
                "path": sealed_authorization["_path"],
                "sha256": legacy_sha256(Path(sealed_authorization["_path"])),
            },
            "prediction_manifest_sha256": legacy_sha256(prediction_manifest_path),
            "prediction_rank_sha256": legacy_sha256(rank_path),
            "outputs": {
                "aggregate_metrics": {"path": str(metric_path), "sha256": legacy_sha256(metric_path)},
                "scope_denominator_audit": {"path": str(scope_audit_path), "sha256": legacy_sha256(scope_audit_path)},
                "baseline_bootstrap": {"path": str(baseline_bootstrap_path), "sha256": legacy_sha256(baseline_bootstrap_path)},
                "paired_bootstrap": {"path": str(bootstrap_path), "sha256": legacy_sha256(bootstrap_path)},
            },
            "evaluation_independence": independence,
            "all_scopes_reported": SCOPES,
            "all_baselines_reported": BASELINES,
            "all_metrics_reported": METRICS,
            "bootstrap": {"replicates": 10000, "prng": "PCG64", "seed": 20260719, "interval": "95% percentile"},
            "environment": {"python": sys.version, "platform": platform.platform()},
        }
    )
    manifest_path = output_dir / "sealed_evaluation_manifest.json"
    manifest_path.write_text(json.dumps(evaluation_manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output_dir": str(output_dir), "evaluation_independence": independence, "figures_generated": False}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
