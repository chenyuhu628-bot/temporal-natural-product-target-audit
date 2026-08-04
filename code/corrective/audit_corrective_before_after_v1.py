"""Stream old and corrected ranks and emit aggregate-only difference summaries."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import statistics
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator


PROTOCOL_ID = "npass_strict_ab_asof_cutoff_corrective_successor_v1_20260728"
PROTOCOL_LOCK_SHA256 = "96befee13ae1d41ad433c8697fac92ccd30fb25e24c3cf1279c6b4b7e040abd9"
BASELINES = [
    "weighted_target_popularity",
    "sequence_3mer_transfer",
    "weighted_morgan_transfer",
    "structure_sequence_pair_neighbor",
]
METRICS = ["Recall@10", "Recall@50", "NDCG@10", "NDCG@50", "MRR"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rank_groups(path: Path) -> Iterator[tuple[tuple[str, str], dict[str, tuple[int, float]]]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"baseline", "query_id", "rank", "score", "eligible_candidate_target_count"}
        target_field = (
            "target_uniprot_accession"
            if "target_uniprot_accession" in (reader.fieldnames or [])
            else "target_id"
        )
        required.add(target_field)
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Rank file lacks fields: {sorted(missing)}")
        current: tuple[str, str] | None = None
        group: dict[str, tuple[int, float]] = {}
        expected = 0
        for row in reader:
            key = (row["baseline"], row["query_id"])
            if current is not None and key != current:
                if len(group) != expected:
                    raise ValueError("Rank group cardinality differs from declared candidate count")
                yield current, group
                group = {}
            if key != current:
                current = key
                expected = int(row["eligible_candidate_target_count"])
            target = row[target_field]
            if target in group:
                raise ValueError("Duplicate target in rank group")
            rank = int(row["rank"])
            score = float(row["score"])
            if not math.isfinite(score):
                raise ValueError("Nonfinite score")
            group[target] = (rank, score)
        if current is not None:
            if len(group) != expected:
                raise ValueError("Final rank group cardinality differs from declared candidate count")
            yield current, group


def read_metrics(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"scope", "baseline", *METRICS}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Metric table lacks fields: {sorted(missing)}")
        result: dict[tuple[str, str], dict[str, str]] = {}
        for row in reader:
            key = (row["scope"], row["baseline"])
            if key in result:
                raise ValueError("Duplicate scope/baseline metric cell")
            result[key] = row
        return result


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protocol-lock", required=True, type=Path)
    parser.add_argument("--old-ranks", required=True, type=Path)
    parser.add_argument("--corrected-ranks", required=True, type=Path)
    parser.add_argument("--old-metrics", required=True, type=Path)
    parser.add_argument("--corrected-metrics", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    args = parser.parse_args()
    if sha256(args.protocol_lock) != PROTOCOL_LOCK_SHA256:
        raise ValueError("Protocol lock hash mismatch")
    if args.output_dir.exists():
        raise FileExistsError(args.output_dir)

    accum = {
        baseline: {
            "row_count": 0,
            "score_changed_rows": 0,
            "rank_changed_rows": 0,
            "absolute_rank_change_sum": 0,
            "maximum_absolute_rank_change": 0,
            "query_blocks": 0,
            "query_blocks_with_any_rank_change": 0,
            "query_blocks_with_top10_membership_change": 0,
            "query_blocks_with_top50_membership_change": 0,
            "top10_symmetric_difference_sum": 0,
            "top50_symmetric_difference_sum": 0,
        }
        for baseline in BASELINES
    }
    old_iter = rank_groups(args.old_ranks)
    new_iter = rank_groups(args.corrected_ranks)
    group_count = 0
    while True:
        try:
            old_key, old_group = next(old_iter)
        except StopIteration:
            try:
                next(new_iter)
            except StopIteration:
                break
            raise ValueError("Corrected rank file has additional groups")
        try:
            new_key, new_group = next(new_iter)
        except StopIteration as exc:
            raise ValueError("Old rank file has additional groups") from exc
        if old_key != new_key or set(old_group) != set(new_group):
            raise ValueError("Old and corrected rank group keysets differ")
        baseline, _ = old_key
        if baseline not in accum:
            raise ValueError(f"Unexpected baseline: {baseline}")
        item = accum[baseline]
        item["query_blocks"] += 1
        group_count += 1
        any_rank_change = False
        old_top10: set[str] = set()
        new_top10: set[str] = set()
        old_top50: set[str] = set()
        new_top50: set[str] = set()
        for target in old_group:
            old_rank, old_score = old_group[target]
            new_rank, new_score = new_group[target]
            item["row_count"] += 1
            if old_score != new_score:
                item["score_changed_rows"] += 1
            delta = abs(old_rank - new_rank)
            if delta:
                any_rank_change = True
                item["rank_changed_rows"] += 1
                item["absolute_rank_change_sum"] += delta
                item["maximum_absolute_rank_change"] = max(item["maximum_absolute_rank_change"], delta)
            if old_rank <= 10:
                old_top10.add(target)
            if new_rank <= 10:
                new_top10.add(target)
            if old_rank <= 50:
                old_top50.add(target)
            if new_rank <= 50:
                new_top50.add(target)
        if any_rank_change:
            item["query_blocks_with_any_rank_change"] += 1
        top10_diff = len(old_top10.symmetric_difference(new_top10))
        top50_diff = len(old_top50.symmetric_difference(new_top50))
        item["top10_symmetric_difference_sum"] += top10_diff
        item["top50_symmetric_difference_sum"] += top50_diff
        item["query_blocks_with_top10_membership_change"] += top10_diff > 0
        item["query_blocks_with_top50_membership_change"] += top50_diff > 0
    if group_count != 222 * 4:
        raise ValueError(f"Rank block count is not 888: {group_count}")

    rank_rows: list[dict[str, object]] = []
    for baseline in BASELINES:
        item = accum[baseline]
        row_count = int(item["row_count"])
        rank_rows.append(
            {
                "baseline": baseline,
                **item,
                "score_changed_row_fraction": item["score_changed_rows"] / row_count,
                "rank_changed_row_fraction": item["rank_changed_rows"] / row_count,
                "mean_absolute_rank_change": item["absolute_rank_change_sum"] / row_count,
            }
        )

    old_metrics = read_metrics(args.old_metrics)
    corrected_metrics = read_metrics(args.corrected_metrics)
    if set(old_metrics) != set(corrected_metrics) or len(old_metrics) != 20:
        raise ValueError("Old and corrected metric matrices differ in keyset or cardinality")
    metric_rows: list[dict[str, object]] = []
    for scope, baseline in sorted(old_metrics):
        for metric in METRICS:
            old_value = float(old_metrics[(scope, baseline)][metric])
            new_value = float(corrected_metrics[(scope, baseline)][metric])
            metric_rows.append(
                {
                    "scope": scope,
                    "baseline": baseline,
                    "metric": metric,
                    "old_value": f"{old_value:.17g}",
                    "corrected_value": f"{new_value:.17g}",
                    "corrected_minus_old": f"{(new_value - old_value):.17g}",
                }
            )

    args.output_dir.mkdir(parents=True, exist_ok=False)
    rank_path = args.output_dir / "rank_change_aggregate_by_baseline.tsv"
    metric_path = args.output_dir / "metric_change_aggregate.tsv"
    write_tsv(rank_path, list(rank_rows[0]), rank_rows)
    write_tsv(metric_path, list(metric_rows[0]), metric_rows)
    metric_deltas = [abs(float(row["corrected_minus_old"])) for row in metric_rows]
    summary = {
        "schema_version": "1.0",
        "protocol_id": PROTOCOL_ID,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "old_ranks": {"path": str(args.old_ranks), "sha256": sha256(args.old_ranks)},
            "corrected_ranks": {"path": str(args.corrected_ranks), "sha256": sha256(args.corrected_ranks)},
            "old_metrics": {"path": str(args.old_metrics), "sha256": sha256(args.old_metrics)},
            "corrected_metrics": {"path": str(args.corrected_metrics), "sha256": sha256(args.corrected_metrics)},
        },
        "rank_aggregate": rank_rows,
        "metric_delta_cells": len(metric_rows),
        "metric_cells_changed": sum(value != 0.0 for value in metric_deltas),
        "maximum_absolute_metric_delta": max(metric_deltas),
        "median_absolute_metric_delta": statistics.median(metric_deltas),
        "outputs": {
            "rank_change": {"path": str(rank_path), "sha256": sha256(rank_path)},
            "metric_change": {"path": str(metric_path), "sha256": sha256(metric_path)},
        },
        "identifiers_emitted": False,
        "claim_boundary": "Aggregate-only audit of a prespecified temporal-purity correction; no model selection or external-validation claim.",
    }
    summary_path = args.output_dir / "before_after_aggregate_summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(summary_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

