#!/usr/bin/env python3
"""Aggregate post hoc attribution of endpoint reciprocal ranks to score ties."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import shutil
import tempfile
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


PROTOCOL_ID = "npass_strict_ab_asof_cutoff_corrective_successor_v1_20260728"
PROTOCOL_LOCK_SHA256 = "96befee13ae1d41ad433c8697fac92ccd30fb25e24c3cf1279c6b4b7e040abd9"
EXPECTED_INPUT_SHA256 = {
    "endpoint": "09296b066a23197a7c178f00514f2b3d9ed7e6f3c459ea92a55e01a6010d1132",
    "ranks": "87739aa818744c7084088d13c386444aa41bbef38c257083325298003181479e",
    "scaffold": "fa0029ef5b7822ad5ca93f7bd93ac808f85f1e0c02e827fa91be375031b2d7af",
    "homology_0_30": "3a8247ed8f683fe6fce5fb345f56e3ec73a872b065eca922e92e494f084a1793",
    "homology_0_50": "ec3bbd435f73bc1c724efdfd94ac10c32d6f9a55bd6c9a9349233a01e5dc7f5b",
    "homology_0_70": "ec3bbd435f73bc1c724efdfd94ac10c32d6f9a55bd6c9a9349233a01e5dc7f5b",
}
SCOPES = ["temporal_strict_ab", "scaffold_cold_strict_ab", "double_cold_0_30", "double_cold_0_50", "double_cold_0_70"]
BASELINES = ["weighted_target_popularity", "sequence_3mer_transfer", "weighted_morgan_transfer", "structure_sequence_pair_neighbor"]
EXPECTED_SCOPE_COUNTS = {
    "temporal_strict_ab": (358, 222, 156),
    "scaffold_cold_strict_ab": (123, 88, 70),
    "double_cold_0_30": (24, 19, 17),
    "double_cold_0_50": (29, 22, 21),
    "double_cold_0_70": (29, 22, 21),
}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def read_gzip_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        fields = list(reader.fieldnames or [])
    return fields, rows


def require_columns(observed: Iterable[str], expected: set[str], label: str) -> None:
    if set(observed) != expected:
        raise ValueError(f"{label} columns mismatch")


def parse_true(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no", ""}:
        return False
    raise ValueError(f"Invalid boolean value: {value!r}")


def load_masks(args: argparse.Namespace) -> dict[str, set[tuple[str, str]]]:
    endpoint_fields, endpoint_rows = read_gzip_tsv(args.endpoint)
    require_columns(endpoint_fields, {"canonical_pair_key", "query_id", "inchikey_full", "uniprot_canonical_accession", "best_strict_evidence_tier", "decision", "c31_leakage_gate_status"}, "endpoint")
    if len(endpoint_rows) != 358:
        raise ValueError("Endpoint must contain exactly 358 relations")
    endpoint_by_pair: dict[str, tuple[str, str]] = {}
    for row in endpoint_rows:
        pair = row["canonical_pair_key"]
        value = (row["query_id"], row["uniprot_canonical_accession"])
        if pair in endpoint_by_pair or not all(value):
            raise ValueError("Duplicate or incomplete endpoint relation")
        endpoint_by_pair[pair] = value

    scaffold_fields, scaffold_rows = read_gzip_tsv(args.scaffold_audit)
    require_columns(scaffold_fields, {"canonical_pair_key", "audit_scaffold_cold_under_selected_policy", "audit_outcome", "audit_eligibility_or_exclusion_reason"}, "scaffold audit")
    scaffold_pairs = {row["canonical_pair_key"] for row in scaffold_rows if parse_true(row["audit_scaffold_cold_under_selected_policy"])}
    if not scaffold_pairs <= set(endpoint_by_pair):
        raise ValueError("Scaffold audit contains non-endpoint pairs")

    homology_targets: dict[str, set[str]] = {}
    for label, path in (("0_30", args.homology_0_30), ("0_50", args.homology_0_50), ("0_70", args.homology_0_70)):
        fields, rows = read_gzip_tsv(path)
        require_columns(fields, {"uniprot_canonical_accession", "is_future_target_homology_cold_candidate", "future_target_coldness_status"}, f"homology {label}")
        targets: set[str] = set()
        for row in rows:
            if parse_true(row["is_future_target_homology_cold_candidate"]):
                targets.add(row["uniprot_canonical_accession"])
        homology_targets[label] = targets

    masks: dict[str, set[tuple[str, str]]] = {
        "temporal_strict_ab": set(endpoint_by_pair.values()),
        "scaffold_cold_strict_ab": {endpoint_by_pair[pair] for pair in scaffold_pairs},
    }
    for label in ("0_30", "0_50", "0_70"):
        masks[f"double_cold_{label}"] = {
            endpoint_by_pair[pair]
            for pair in scaffold_pairs
            if endpoint_by_pair[pair][1] in homology_targets[label]
        }
    for scope, relations in masks.items():
        expected_relations, expected_queries, expected_targets = EXPECTED_SCOPE_COUNTS[scope]
        counts = (len(relations), len({query for query, _ in relations}), len({target for _, target in relations}))
        if counts != (expected_relations, expected_queries, expected_targets):
            raise ValueError(f"Scope count mismatch for {scope}: {counts}")
    return masks


def empty_aggregate(scope: str, baseline: str, relations: set[tuple[str, str]]) -> dict[str, Any]:
    return {
        "scope": scope,
        "baseline": baseline,
        "endpoint_relation_count": len(relations),
        "endpoint_query_count": len({query for query, _ in relations}),
        "endpoint_target_count": len({target for _, target in relations}),
        "endpoint_zero_score_relation_count": 0,
        "endpoint_positive_score_relation_count": 0,
        "endpoint_exact_tied_relation_count": 0,
        "endpoint_zero_score_tied_relation_count": 0,
        "endpoint_positive_score_tied_relation_count": 0,
        "endpoint_positive_unique_score_relation_count": 0,
        "endpoint_rank_le_50_relation_count": 0,
        "endpoint_rank_gt_50_relation_count": 0,
        "scope_all_zero_vector_query_count": 0,
        "mrr_first_hit_zero_score_tied_query_count": 0,
        "mrr_first_hit_positive_score_tied_query_count": 0,
        "mrr_first_hit_positive_unique_score_query_count": 0,
        "mrr_first_hit_other_query_count": 0,
        "status": "post_hoc_outcome_visible_attribution_only",
    }


def audit_group(
    baseline: str,
    query_id: str,
    rows: list[tuple[str, int, float, int]],
    relevant_by_scope_query: dict[str, dict[str, set[str]]],
    aggregates: dict[tuple[str, str], dict[str, Any]],
) -> None:
    if baseline not in BASELINES:
        raise ValueError(f"Unexpected baseline: {baseline}")
    if not rows:
        raise ValueError("Empty rank block")
    eligible_counts = {row[3] for row in rows}
    if len(eligible_counts) != 1 or next(iter(eligible_counts)) != len(rows):
        raise ValueError("Eligible-candidate count mismatch within rank block")
    ranks = [row[1] for row in rows]
    if sorted(ranks) != list(range(1, len(rows) + 1)):
        raise ValueError("Rank block is not an exact permutation")
    targets = [row[0] for row in rows]
    if len(targets) != len(set(targets)):
        raise ValueError("Duplicate target in rank block")
    score_counts = Counter(row[2] for row in rows)
    all_zero = all(score == 0.0 for _, _, score, _ in rows)
    rank_by_target = {target: (rank, score) for target, rank, score, _ in rows}

    for scope in SCOPES:
        relevant = relevant_by_scope_query[scope].get(query_id)
        if not relevant:
            continue
        aggregate = aggregates[(scope, baseline)]
        if all_zero:
            aggregate["scope_all_zero_vector_query_count"] += 1
        first_rank = math.inf
        first_score = math.nan
        first_tie_size = 0
        for target in relevant:
            if target not in rank_by_target:
                raise ValueError(f"Endpoint target missing from eligible ranks in {scope}")
            rank, score = rank_by_target[target]
            tie_size = score_counts[score]
            if score == 0.0:
                aggregate["endpoint_zero_score_relation_count"] += 1
            elif score > 0.0:
                aggregate["endpoint_positive_score_relation_count"] += 1
            else:
                raise ValueError("Negative score is outside the frozen baseline contract")
            if tie_size > 1:
                aggregate["endpoint_exact_tied_relation_count"] += 1
                if score == 0.0:
                    aggregate["endpoint_zero_score_tied_relation_count"] += 1
                else:
                    aggregate["endpoint_positive_score_tied_relation_count"] += 1
            elif score > 0.0:
                aggregate["endpoint_positive_unique_score_relation_count"] += 1
            if rank <= 50:
                aggregate["endpoint_rank_le_50_relation_count"] += 1
            else:
                aggregate["endpoint_rank_gt_50_relation_count"] += 1
            if rank < first_rank:
                first_rank, first_score, first_tie_size = rank, score, tie_size
        if first_score == 0.0 and first_tie_size > 1:
            aggregate["mrr_first_hit_zero_score_tied_query_count"] += 1
        elif first_score > 0.0 and first_tie_size > 1:
            aggregate["mrr_first_hit_positive_score_tied_query_count"] += 1
        elif first_score > 0.0 and first_tie_size == 1:
            aggregate["mrr_first_hit_positive_unique_score_query_count"] += 1
        else:
            aggregate["mrr_first_hit_other_query_count"] += 1


def run(args: argparse.Namespace) -> None:
    if args.output_dir.exists():
        raise FileExistsError(f"Create-once output already exists: {args.output_dir}")
    if sha256_file(args.protocol_lock) != PROTOCOL_LOCK_SHA256:
        raise ValueError("Protocol-lock hash mismatch")
    implementation_lock = load_json(args.implementation_lock)
    if implementation_lock.get("lock_state") != "LOCKED_BEFORE_REAL_ENDPOINT_SCORE_ATTRIBUTION_EXECUTION":
        raise ValueError("Endpoint-score attribution implementation lock is not executable")
    if implementation_lock.get("protocol_id") != PROTOCOL_ID:
        raise ValueError("Implementation-lock protocol mismatch")
    if implementation_lock.get("real_audit_executed_before_lock") is not False:
        raise ValueError("Implementation lock does not preserve pre-execution custody")
    if implementation_lock.get("outcome_visible_before_lock") is not True:
        raise ValueError("Implementation lock must disclose outcome visibility")
    if implementation_lock.get("script_sha256") != sha256_file(Path(__file__)):
        raise ValueError("Implementation-lock script hash mismatch")
    if implementation_lock.get("governance_amendment_sha256") != sha256_file(args.governance_amendment):
        raise ValueError("Implementation-lock governance hash mismatch")
    input_paths = {
        "endpoint": args.endpoint,
        "ranks": args.ranks,
        "scaffold": args.scaffold_audit,
        "homology_0_30": args.homology_0_30,
        "homology_0_50": args.homology_0_50,
        "homology_0_70": args.homology_0_70,
    }
    for role, path in input_paths.items():
        if sha256_file(path) != EXPECTED_INPUT_SHA256[role]:
            raise ValueError(f"Frozen input hash mismatch: {role}")

    masks = load_masks(args)
    relevant_by_scope_query: dict[str, dict[str, set[str]]] = {}
    for scope, relations in masks.items():
        mapping: dict[str, set[str]] = defaultdict(set)
        for query, target in relations:
            mapping[query].add(target)
        relevant_by_scope_query[scope] = dict(mapping)
    aggregates = {(scope, baseline): empty_aggregate(scope, baseline, masks[scope]) for scope in SCOPES for baseline in BASELINES}

    expected_rank_fields = {"protocol_id", "baseline", "query_id", "query_compound_inchikey_full", "target_uniprot_accession", "rank", "score", "eligible_candidate_target_count"}
    current_key: tuple[str, str] | None = None
    current_rows: list[tuple[str, int, float, int]] = []
    observed_blocks: set[tuple[str, str]] = set()
    rank_row_count = 0
    with gzip.open(args.ranks, "rt", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_columns(reader.fieldnames or [], expected_rank_fields, "ranks")
        for raw in reader:
            if raw["protocol_id"] != PROTOCOL_ID:
                raise ValueError("Rank-row protocol mismatch")
            key = (raw["baseline"], raw["query_id"])
            if current_key is None:
                current_key = key
            if key != current_key:
                if key in observed_blocks:
                    raise ValueError("Non-contiguous recurring rank block")
                observed_blocks.add(current_key)
                audit_group(current_key[0], current_key[1], current_rows, relevant_by_scope_query, aggregates)
                current_key = key
                current_rows = []
            score = float(raw["score"])
            if not math.isfinite(score):
                raise ValueError("Nonfinite score")
            current_rows.append((raw["target_uniprot_accession"], int(raw["rank"]), score, int(raw["eligible_candidate_target_count"])))
            rank_row_count += 1
    if current_key is None:
        raise ValueError("Empty rank ledger")
    if current_key in observed_blocks:
        raise ValueError("Recurring final rank block")
    observed_blocks.add(current_key)
    audit_group(current_key[0], current_key[1], current_rows, relevant_by_scope_query, aggregates)
    if rank_row_count != 3_658_128 or len(observed_blocks) != 888:
        raise ValueError("Rank-ledger dimensions changed")

    rows = [aggregates[(scope, baseline)] for scope in SCOPES for baseline in BASELINES]
    for row in rows:
        if row["endpoint_zero_score_relation_count"] + row["endpoint_positive_score_relation_count"] != row["endpoint_relation_count"]:
            raise ValueError("Endpoint score attribution does not sum to the relation denominator")
        if row["endpoint_rank_le_50_relation_count"] + row["endpoint_rank_gt_50_relation_count"] != row["endpoint_relation_count"]:
            raise ValueError("Rank-boundary attribution does not sum to the relation denominator")
        mrr_query_sum = sum(row[key] for key in (
            "mrr_first_hit_zero_score_tied_query_count",
            "mrr_first_hit_positive_score_tied_query_count",
            "mrr_first_hit_positive_unique_score_query_count",
            "mrr_first_hit_other_query_count",
        ))
        if mrr_query_sum != row["endpoint_query_count"]:
            raise ValueError("MRR attribution does not sum to the query denominator")

    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    temporary = Path(tempfile.mkdtemp(prefix=f".{args.output_dir.name}.", dir=args.output_dir.parent))
    try:
        output_tsv = temporary / "endpoint_score_attribution_aggregate.tsv"
        fields = list(rows[0])
        with output_tsv.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        summary = {
            "audit_id": "endpoint_score_attribution_posthoc_v1",
            "protocol_id": PROTOCOL_ID,
            "timing": "outcome_visible_post_hoc_explanatory_audit",
            "row_count": len(rows),
            "aggregate_only": True,
            "identifier_bearing_output": False,
            "rank_row_count_verified": rank_row_count,
            "rank_block_count_verified": len(observed_blocks),
            "interpretation_boundary": "A reciprocal rank inside a zero-score tie block is salt-derived; positive unique score does not imply top-50 retrieval.",
        }
        summary_path = temporary / "endpoint_score_attribution_summary.json"
        summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        governance_sha = sha256_file(args.governance_amendment)
        manifest = {
            "audit_id": summary["audit_id"],
            "protocol_id": PROTOCOL_ID,
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "aggregate_only": True,
            "identifier_bearing_output": False,
            "governance": {"basename": args.governance_amendment.name, "sha256": governance_sha},
            "implementation_lock": {"basename": args.implementation_lock.name, "sha256": sha256_file(args.implementation_lock)},
            "script": {"basename": Path(__file__).name, "sha256": sha256_file(Path(__file__))},
            "inputs": [
                {"role": role, "basename": path.name, "bytes": path.stat().st_size, "sha256": sha256_file(path)}
                for role, path in input_paths.items()
            ],
            "outputs": [
                {"basename": output_tsv.name, "bytes": output_tsv.stat().st_size, "sha256": sha256_file(output_tsv)},
                {"basename": summary_path.name, "bytes": summary_path.stat().st_size, "sha256": sha256_file(summary_path)},
            ],
        }
        (temporary / "run_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        os.replace(temporary, args.output_dir)
    except BaseException:
        shutil.rmtree(temporary, ignore_errors=True)
        raise
    print(f"Wrote 20 aggregate endpoint-score attribution rows to {args.output_dir}")


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--protocol-lock", required=True, type=Path)
    result.add_argument("--governance-amendment", required=True, type=Path)
    result.add_argument("--implementation-lock", required=True, type=Path)
    result.add_argument("--endpoint", required=True, type=Path)
    result.add_argument("--ranks", required=True, type=Path)
    result.add_argument("--scaffold-audit", required=True, type=Path)
    result.add_argument("--homology-0-30", required=True, type=Path)
    result.add_argument("--homology-0-50", required=True, type=Path)
    result.add_argument("--homology-0-70", required=True, type=Path)
    result.add_argument("--output-dir", required=True, type=Path)
    return result


if __name__ == "__main__":
    run(parser().parse_args())
