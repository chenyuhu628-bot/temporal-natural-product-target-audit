"""Aggregate-only score degeneracy and rank-boundary tie audit."""

from __future__ import annotations

import argparse
import math
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from audit_common import (
    PROTOCOL_ID,
    distribution,
    finalize_manifest,
    input_descriptor,
    open_dict_reader,
    require_new_output_dir,
    require_protocol_lock,
    write_json_new,
    write_tsv_new,
)


AUDIT_ID = "score_degeneracy_and_boundary_ties_v1"
LOCKED_BASELINES = (
    "weighted_target_popularity",
    "sequence_3mer_transfer",
    "weighted_morgan_transfer",
    "structure_sequence_pair_neighbor",
)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--protocol-lock", required=True, type=Path)
    result.add_argument("--ranks", required=True, type=Path)
    result.add_argument("--output-dir", required=True, type=Path)
    result.add_argument("--expected-query-count", type=int, default=222)
    result.add_argument("--boundary-k", type=int, action="append", default=None)
    return result


def boundary_metrics(sorted_rows: list[dict[str, Any]], k: int) -> dict[str, Any]:
    n = len(sorted_rows)
    if n < k:
        return {
            "defined": False,
            "crosses_boundary": False,
            "tie_block_size": 0,
            "selected_from_boundary_tie": 0,
            "excluded_from_boundary_tie": 0,
            "boundary_score_is_zero": False,
        }
    boundary_score = sorted_rows[k - 1]["score"]
    equal_ranks = [row["rank"] for row in sorted_rows if row["score"] == boundary_score]
    start = min(equal_ranks)
    end = max(equal_ranks)
    crosses = start <= k < end
    return {
        "defined": True,
        "crosses_boundary": crosses,
        "tie_block_size": len(equal_ranks),
        "selected_from_boundary_tie": k - start + 1,
        "excluded_from_boundary_tie": end - k,
        "boundary_score_is_zero": boundary_score == 0.0,
    }


def summarize_group(rows: list[dict[str, str]], boundaries: tuple[int, ...]) -> dict[str, Any]:
    if not rows:
        raise ValueError("Cannot summarize an empty query-baseline group")
    baseline = rows[0]["baseline"]
    query = rows[0]["query_id"]
    compounds = {row["query_compound_inchikey_full"] for row in rows}
    if len(compounds) != 1:
        raise ValueError("A query-baseline group maps to multiple compounds")
    if any(row["baseline"] != baseline or row["query_id"] != query for row in rows):
        raise ValueError("Mixed group supplied to group summarizer")

    parsed = []
    targets: set[str] = set()
    eligible_counts: set[int] = set()
    for row in rows:
        if row["protocol_id"] != PROTOCOL_ID:
            raise ValueError("Rank file protocol_id does not match corrective protocol")
        target = row["target_uniprot_accession"]
        if target in targets:
            raise ValueError("Duplicate target in a query-baseline group")
        targets.add(target)
        rank = int(row["rank"])
        score = float(row["score"])
        if not math.isfinite(score):
            raise ValueError("Nonfinite score encountered")
        if score < 0.0:
            raise ValueError("Locked baselines must emit nonnegative scores")
        eligible_counts.add(int(row["eligible_candidate_target_count"]))
        parsed.append({"rank": rank, "score": score})
    if len(eligible_counts) != 1:
        raise ValueError("Eligible candidate count changes within a group")
    eligible = next(iter(eligible_counts))
    if eligible != len(parsed):
        raise ValueError("Rank row count differs from eligible candidate count")
    parsed.sort(key=lambda row: row["rank"])
    if [row["rank"] for row in parsed] != list(range(1, eligible + 1)):
        raise ValueError("Ranks are not a complete 1..N permutation")
    scores = [row["score"] for row in parsed]
    if any(left < right for left, right in zip(scores, scores[1:])):
        raise ValueError("Scores increase with worsening rank")

    tie_sizes = list(Counter(scores).values())
    unique_scores = len(tie_sizes)
    positive = sum(score > 0.0 for score in scores)
    result: dict[str, Any] = {
        "baseline": baseline,
        "eligible_candidate_count": eligible,
        "positive_target_count": positive,
        "zero_target_count": eligible - positive,
        "unique_score_count": unique_scores,
        "largest_tie_block_size": max(tie_sizes),
        "non_singleton_tie_block_count": sum(size > 1 for size in tie_sizes),
        "has_any_tie": unique_scores < eligible,
        "all_zero": positive == 0,
        "constant_score": unique_scores == 1,
        "boundaries": {},
    }
    for k in boundaries:
        result["boundaries"][str(k)] = boundary_metrics(parsed, k)
    return result


def aggregate_baseline(groups: list[dict[str, Any]], boundaries: tuple[int, ...]) -> dict[str, Any]:
    query_count = len(groups)
    total_rows = sum(group["eligible_candidate_count"] for group in groups)
    positive_rows = sum(group["positive_target_count"] for group in groups)
    output: dict[str, Any] = {
        "baseline": groups[0]["baseline"],
        "query_count": query_count,
        "rank_row_count": total_rows,
        "all_zero_query_count": sum(group["all_zero"] for group in groups),
        "all_zero_query_fraction": sum(group["all_zero"] for group in groups) / query_count,
        "constant_score_query_count": sum(group["constant_score"] for group in groups),
        "query_count_with_any_tie": sum(group["has_any_tie"] for group in groups),
        "positive_rank_row_count": positive_rows,
        "positive_rank_row_fraction": positive_rows / total_rows,
        "eligible_candidate_count_distribution": distribution(
            [group["eligible_candidate_count"] for group in groups]
        ),
        "positive_target_count_distribution": distribution(
            [group["positive_target_count"] for group in groups]
        ),
        "unique_score_count_distribution": distribution(
            [group["unique_score_count"] for group in groups]
        ),
        "largest_tie_block_size_distribution": distribution(
            [group["largest_tie_block_size"] for group in groups]
        ),
        "boundary_audits": {},
    }
    for k in boundaries:
        records = [group["boundaries"][str(k)] for group in groups]
        defined = [record for record in records if record["defined"]]
        output["boundary_audits"][str(k)] = {
            "defined_query_count": len(defined),
            "boundary_tie_query_count": sum(record["crosses_boundary"] for record in defined),
            "boundary_tie_query_fraction": (
                sum(record["crosses_boundary"] for record in defined) / len(defined) if defined else None
            ),
            "zero_score_at_boundary_query_count": sum(
                record["boundary_score_is_zero"] for record in defined
            ),
            "boundary_tie_block_size_distribution": distribution(
                [record["tie_block_size"] for record in defined if record["crosses_boundary"]]
            ),
            "tie_members_selected_by_salt_total": sum(
                record["selected_from_boundary_tie"]
                for record in defined
                if record["crosses_boundary"]
            ),
            "tie_members_excluded_by_salt_total": sum(
                record["excluded_from_boundary_tie"]
                for record in defined
                if record["crosses_boundary"]
            ),
        }
    return output


def main() -> int:
    args = parser().parse_args()
    require_protocol_lock(args.protocol_lock)
    output_dir = require_new_output_dir(args.output_dir)
    boundaries = tuple(sorted(set(args.boundary_k or [10, 50])))
    if not boundaries or any(k < 1 for k in boundaries):
        raise ValueError("Every boundary K must be positive")
    rank_input = input_descriptor("corrective_full_rank_file", args.ranks)

    required = {
        "protocol_id",
        "baseline",
        "query_id",
        "query_compound_inchikey_full",
        "target_uniprot_accession",
        "rank",
        "score",
        "eligible_candidate_target_count",
    }
    groups_by_baseline: dict[str, list[dict[str, Any]]] = defaultdict(list)
    seen_groups: set[tuple[str, str]] = set()
    current_key: tuple[str, str] | None = None
    current_rows: list[dict[str, str]] = []

    def flush() -> None:
        nonlocal current_rows, current_key
        if current_key is None:
            return
        if current_key in seen_groups:
            raise ValueError("A query-baseline group is noncontiguous or repeated")
        seen_groups.add(current_key)
        summary = summarize_group(current_rows, boundaries)
        groups_by_baseline[summary["baseline"]].append(summary)
        current_rows = []

    with open_dict_reader(args.ranks) as reader:
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Rank file lacks required fields: {sorted(missing)}")
        for row in reader:
            key = (row["baseline"], row["query_id"])
            if current_key is None:
                current_key = key
            elif key != current_key:
                flush()
                current_key = key
            current_rows.append(row)
    flush()

    if set(groups_by_baseline) != set(LOCKED_BASELINES):
        raise ValueError(
            f"Baseline set differs from lock: {sorted(groups_by_baseline)}"
        )
    aggregate = []
    for baseline in LOCKED_BASELINES:
        groups = groups_by_baseline[baseline]
        if len(groups) != args.expected_query_count:
            raise ValueError(
                f"{baseline} has {len(groups)} queries, expected {args.expected_query_count}"
            )
        aggregate.append(aggregate_baseline(groups, boundaries))

    summary = {
        "audit_id": AUDIT_ID,
        "protocol_id": PROTOCOL_ID,
        "score_tie_definition": "exact equality after parsing the emitted decimal as IEEE-754 float",
        "boundary_tie_definition": "the equal-score block containing rank K extends below K",
        "deterministic_tie_interpretation": (
            "membership within a boundary-crossing tie is assigned by the locked salted hash; "
            "it is not evidence of score separation"
        ),
        "boundary_k": list(boundaries),
        "baselines": aggregate,
    }
    flat_rows = []
    for item in aggregate:
        flat: dict[str, Any] = {
            "baseline": item["baseline"],
            "query_count": item["query_count"],
            "rank_row_count": item["rank_row_count"],
            "all_zero_query_count": item["all_zero_query_count"],
            "all_zero_query_fraction": item["all_zero_query_fraction"],
            "constant_score_query_count": item["constant_score_query_count"],
            "query_count_with_any_tie": item["query_count_with_any_tie"],
            "positive_rank_row_count": item["positive_rank_row_count"],
            "positive_rank_row_fraction": item["positive_rank_row_fraction"],
        }
        for k in boundaries:
            boundary = item["boundary_audits"][str(k)]
            flat[f"top{k}_boundary_tie_query_count"] = boundary["boundary_tie_query_count"]
            flat[f"top{k}_zero_score_boundary_query_count"] = boundary[
                "zero_score_at_boundary_query_count"
            ]
        flat_rows.append(flat)

    output_dir.mkdir(parents=True, exist_ok=False)
    summary_name = "rank_score_aggregate_summary.json"
    table_name = "rank_score_aggregate_by_baseline.tsv"
    write_json_new(output_dir / summary_name, summary)
    table_fields = list(flat_rows[0])
    write_tsv_new(output_dir / table_name, table_fields, flat_rows)
    manifest = finalize_manifest(
        output_dir=output_dir,
        audit_id=AUDIT_ID,
        script_path=Path(__file__),
        inputs=[rank_input],
        output_names=[summary_name, table_name],
        extra={"query_baseline_group_count": sum(len(value) for value in groups_by_baseline.values())},
    )
    write_json_new(output_dir / "run_manifest.json", manifest)
    print(
        f"{AUDIT_ID}: wrote aggregate-only summaries for "
        f"{sum(len(value) for value in groups_by_baseline.values())} query-baseline groups"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
