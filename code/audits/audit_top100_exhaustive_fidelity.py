"""Aggregate-only top-100 versus exhaustive pair-neighbour fidelity sensitivity."""

from __future__ import annotations

import argparse
import math
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np
from rdkit import DataStructs

from audit_common import (
    PROTOCOL_ID,
    distribution,
    finalize_manifest,
    input_descriptor,
    open_dict_reader,
    parse_bool,
    require_new_output_dir,
    require_protocol_lock,
    write_json_new,
    write_tsv_new,
)


SUITE_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = SUITE_ROOT.parents[1]
sys.path.insert(0, str(WORKSPACE / "scripts"))

from pu_baseline_io import read_fasta, read_tsv_gz  # noqa: E402
from pu_retrieval_metrics import macro_average, query_metrics, rank_scores  # noqa: E402
from pu_retrieval_scores import build_sequence_kmer_matrix, build_train_maps, morgan_fingerprints  # noqa: E402


AUDIT_ID = "pair_neighbor_top100_vs_exhaustive_fidelity_v1"
LOCKED_TOP_K = 100
LOCKED_SCORE_ATOL = 1e-7
LOCKED_TIE_SALT = "npass_strict_ab_doublecold_successor_v1_20260719"
LOCKED_WEIGHTS = {
    "A_affinity_candidate": 1.0,
    "B_quantitative_functional_candidate": 0.7,
}
DEFAULT_BOUNDARIES = (10, 50)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--protocol-lock", required=True, type=Path)
    result.add_argument("--historical-pairs", required=True, type=Path)
    result.add_argument("--scoring-queries", required=True, type=Path)
    result.add_argument("--historical-compounds", required=True, type=Path)
    result.add_argument("--query-compounds", required=True, type=Path)
    result.add_argument("--candidate-targets", required=True, type=Path)
    result.add_argument("--candidate-sequences", required=True, type=Path)
    result.add_argument("--evaluation-pairs", required=True, type=Path)
    result.add_argument("--scaffold-audit", required=True, type=Path)
    result.add_argument("--homology-0-30", required=True, type=Path)
    result.add_argument("--homology-0-50", required=True, type=Path)
    result.add_argument("--homology-0-70", required=True, type=Path)
    result.add_argument("--output-dir", required=True, type=Path)
    result.add_argument("--expected-historical-pair-count", type=int, default=4990)
    result.add_argument("--expected-query-count", type=int, default=222)
    result.add_argument("--expected-target-count", type=int, default=4123)
    result.add_argument("--exhaustive-chunk-size", type=int, default=512)
    return result


def load_role_compounds(path: Path, role: str) -> dict[str, str]:
    fields, rows = read_tsv_gz(path)
    required = {"inchikey_full", "representative_smiles"}
    if not required.issubset(fields):
        raise ValueError(f"{role} compound table lacks {sorted(required.difference(fields))}")
    output: dict[str, str] = {}
    for row in rows:
        compound = row["inchikey_full"]
        smiles = row["representative_smiles"]
        if not compound or not smiles:
            raise ValueError(f"{role} compound table contains an empty key or structure")
        if compound in output and output[compound] != smiles:
            raise ValueError(f"{role} compound table has conflicting structures")
        output[compound] = smiles
    if not output:
        raise ValueError(f"{role} compound table is empty")
    return output


def build_activation(
    *,
    query_compound: str,
    train_by_compound: dict[str, list[tuple[int, float]]],
    historical_fingerprints: dict[str, Any],
    query_fingerprints: dict[str, Any],
    target_to_historical_column: dict[int, int],
) -> np.ndarray:
    historical_compounds = list(train_by_compound)
    activation = np.zeros(len(target_to_historical_column), dtype=np.float32)
    if not historical_compounds:
        return activation
    similarities = DataStructs.BulkTanimotoSimilarity(
        query_fingerprints[query_compound],
        [historical_fingerprints[compound] for compound in historical_compounds],
    )
    for compound, similarity in zip(historical_compounds, similarities):
        for target_index, evidence_weight in train_by_compound[compound]:
            column = target_to_historical_column[target_index]
            candidate = float(similarity) * float(evidence_weight)
            if candidate > activation[column]:
                activation[column] = candidate
    return activation


def exhaustive_scores(
    sequence_similarity: np.ndarray, activation: np.ndarray, chunk_size: int
) -> np.ndarray:
    if chunk_size < 1:
        raise ValueError("Exhaustive chunk size must be positive")
    output = np.empty(sequence_similarity.shape[0], dtype=np.float32)
    for start in range(0, sequence_similarity.shape[0], chunk_size):
        stop = min(start + chunk_size, sequence_similarity.shape[0])
        output[start:stop] = np.max(
            sequence_similarity[start:stop] * activation[np.newaxis, :], axis=1
        ).astype(np.float32, copy=False)
    return output


def spearman_complete_ranks(left: np.ndarray, right: np.ndarray) -> float:
    if left.shape != right.shape or left.ndim != 1:
        raise ValueError("Rank arrays must be aligned one-dimensional vectors")
    n = left.size
    if n < 2:
        return 1.0
    differences = left.astype(np.float64) - right.astype(np.float64)
    return float(1.0 - 6.0 * np.dot(differences, differences) / (n * (n * n - 1.0)))


def top_set_metrics(
    top_order: np.ndarray, exhaustive_order: np.ndarray, k: int
) -> dict[str, Any]:
    size = min(k, len(top_order), len(exhaustive_order))
    left = set(int(value) for value in top_order[:size])
    right = set(int(value) for value in exhaustive_order[:size])
    intersection = len(left.intersection(right))
    union = len(left.union(right))
    return {
        "identical": left == right,
        "intersection_count": intersection,
        "jaccard": intersection / union if union else 1.0,
    }


def load_bool_map(path: Path, key_field: str, value_field: str, role: str) -> dict[str, bool]:
    output: dict[str, bool] = {}
    with open_dict_reader(path) as reader:
        fields = reader.fieldnames or []
        missing = {key_field, value_field}.difference(fields)
        if missing:
            raise ValueError(f"{role} lacks required fields: {sorted(missing)}")
        for row in reader:
            key = row[key_field].strip()
            if not key or key in output:
                raise ValueError(f"{role} has an empty or duplicate key")
            output[key] = parse_bool(row[value_field])
    if not output:
        raise ValueError(f"{role} is empty")
    return output


def build_locked_scope_pairs(
    *,
    evaluation_pairs: Path,
    scaffold_audit: Path,
    homology_paths: dict[str, Path],
    query_id_by_compound: dict[str, str],
    target_index: dict[str, int],
) -> dict[str, dict[str, list[int]]]:
    endpoint_rows: list[dict[str, str]] = []
    with open_dict_reader(evaluation_pairs) as reader:
        fields = reader.fieldnames or []
        required = {
            "canonical_pair_key",
            "query_id",
            "inchikey_full",
            "uniprot_canonical_accession",
        }
        missing = required.difference(fields)
        if missing:
            raise ValueError(f"Evaluation pair ledger lacks {sorted(missing)}")
        seen_pairs: set[str] = set()
        for row in reader:
            pair = row["canonical_pair_key"].strip()
            query_id = row["query_id"].strip()
            compound = row["inchikey_full"].strip()
            target = row["uniprot_canonical_accession"].strip()
            if pair in seen_pairs:
                raise ValueError("Evaluation pair ledger contains a duplicate relation")
            seen_pairs.add(pair)
            if query_id_by_compound.get(compound) != query_id:
                raise ValueError("Evaluation pair query mapping differs from scoring queries")
            if target not in target_index:
                raise ValueError("Evaluation pair target is outside the candidate universe")
            endpoint_rows.append(row)
    if len(endpoint_rows) != 358:
        raise ValueError(f"Frozen endpoint relation count is {len(endpoint_rows)}, expected 358")

    scaffold = load_bool_map(
        scaffold_audit,
        "canonical_pair_key",
        "audit_scaffold_cold_under_selected_policy",
        "scaffold audit",
    )
    endpoint_keys = {row["canonical_pair_key"] for row in endpoint_rows}
    if set(scaffold) != endpoint_keys:
        raise ValueError("Scaffold audit keyset differs from frozen endpoint")
    endpoint_targets = {row["uniprot_canonical_accession"] for row in endpoint_rows}
    homology = {
        threshold: load_bool_map(
            path,
            "uniprot_canonical_accession",
            "is_future_target_homology_cold_candidate",
            f"homology {threshold}",
        )
        for threshold, path in homology_paths.items()
    }
    for threshold, flags in homology.items():
        if set(flags) != endpoint_targets:
            raise ValueError(f"Homology {threshold} keyset differs from frozen endpoint targets")

    output: dict[str, dict[str, list[int]]] = {
        "temporal_strict_ab": defaultdict(list),
        "scaffold_cold_strict_ab": defaultdict(list),
        "double_cold_0_30": defaultdict(list),
        "double_cold_0_50": defaultdict(list),
        "double_cold_0_70": defaultdict(list),
    }
    for row in endpoint_rows:
        query_id = row["query_id"]
        target = row["uniprot_canonical_accession"]
        target_idx = target_index[target]
        output["temporal_strict_ab"][query_id].append(target_idx)
        if scaffold[row["canonical_pair_key"]]:
            output["scaffold_cold_strict_ab"][query_id].append(target_idx)
            for threshold in ("0_30", "0_50", "0_70"):
                if homology[threshold][target]:
                    output[f"double_cold_{threshold}"][query_id].append(target_idx)
    return {scope: dict(rows) for scope, rows in output.items()}


def aggregate_fidelity(query_rows: list[dict[str, Any]]) -> dict[str, Any]:
    total_eligible = sum(row["eligible_candidate_count"] for row in query_rows)
    exact_changed = sum(row["exact_score_changed_target_count"] for row in query_rows)
    tolerance_changed = sum(
        row["tolerance_score_changed_target_count"] for row in query_rows
    )
    exact_rank_changed = sum(row["rank_changed_target_count"] for row in query_rows)
    output: dict[str, Any] = {
        "query_count": len(query_rows),
        "eligible_query_target_count": total_eligible,
        "query_count_with_any_exact_score_change": sum(
            row["exact_score_changed_target_count"] > 0 for row in query_rows
        ),
        "query_count_with_any_tolerance_score_change": sum(
            row["tolerance_score_changed_target_count"] > 0 for row in query_rows
        ),
        "exact_score_changed_target_count": exact_changed,
        "exact_score_changed_target_fraction": exact_changed / total_eligible,
        "tolerance_score_changed_target_count": tolerance_changed,
        "tolerance_score_changed_target_fraction": tolerance_changed / total_eligible,
        "exhaustive_strictly_higher_target_count": sum(
            row["exhaustive_strictly_higher_target_count"] for row in query_rows
        ),
        "query_count_with_any_rank_change": sum(
            row["rank_changed_target_count"] > 0 for row in query_rows
        ),
        "rank_changed_target_count": exact_rank_changed,
        "rank_changed_target_fraction": exact_rank_changed / total_eligible,
        "changed_target_count_per_query_distribution": distribution(
            [row["tolerance_score_changed_target_count"] for row in query_rows]
        ),
        "maximum_absolute_score_error_distribution": distribution(
            [row["maximum_absolute_score_error"] for row in query_rows]
        ),
        "mean_absolute_score_error_distribution": distribution(
            [row["mean_absolute_score_error"] for row in query_rows]
        ),
        "rank_spearman_distribution": distribution(
            [row["rank_spearman"] for row in query_rows]
        ),
        "mean_absolute_rank_shift_distribution": distribution(
            [row["mean_absolute_rank_shift"] for row in query_rows]
        ),
        "maximum_absolute_rank_shift_distribution": distribution(
            [row["maximum_absolute_rank_shift"] for row in query_rows]
        ),
        "top_k": {},
    }
    for k in DEFAULT_BOUNDARIES:
        output["top_k"][str(k)] = {
            "query_count_with_changed_membership": sum(
                not row["top_k"][str(k)]["identical"] for row in query_rows
            ),
            "jaccard_distribution": distribution(
                [row["top_k"][str(k)]["jaccard"] for row in query_rows]
            ),
        }
    return output


def scope_metric_summary(
    *,
    scopes: dict[str, dict[str, list[int]]],
    ranks_top100: dict[str, np.ndarray],
    ranks_exhaustive: dict[str, np.ndarray],
) -> list[dict[str, Any]]:
    rows = []
    for scope, positives_by_query in sorted(scopes.items()):
        top_query_metrics = []
        exhaustive_query_metrics = []
        relation_count = 0
        for query_id, target_indices in positives_by_query.items():
            if query_id not in ranks_top100:
                raise ValueError("Scope contains a query absent from scored rankings")
            top_ranks = [int(ranks_top100[query_id][index]) for index in target_indices]
            exhaustive_ranks = [int(ranks_exhaustive[query_id][index]) for index in target_indices]
            if any(rank < 1 for rank in top_ranks + exhaustive_ranks):
                raise ValueError("A scope positive is masked or unranked")
            top_query_metrics.append(query_metrics(top_ranks, DEFAULT_BOUNDARIES))
            exhaustive_query_metrics.append(query_metrics(exhaustive_ranks, DEFAULT_BOUNDARIES))
            relation_count += len(target_indices)
        top_aggregate = macro_average(top_query_metrics)
        exhaustive_aggregate = macro_average(exhaustive_query_metrics)
        for metric in sorted(top_aggregate):
            rows.append(
                {
                    "scope": scope,
                    "metric": metric,
                    "relation_count": relation_count,
                    "query_count": len(positives_by_query),
                    "top100_value": top_aggregate[metric],
                    "exhaustive_value": exhaustive_aggregate[metric],
                    "exhaustive_minus_top100": exhaustive_aggregate[metric]
                    - top_aggregate[metric],
                }
            )
    return rows


def main() -> int:
    args = parser().parse_args()
    started = time.perf_counter()
    require_protocol_lock(args.protocol_lock)
    output_dir = require_new_output_dir(args.output_dir)
    if args.exhaustive_chunk_size < 1:
        raise ValueError("Exhaustive chunk size must be positive")
    input_paths = {
        "corrected_historical_pairs": args.historical_pairs,
        "frozen_scoring_queries": args.scoring_queries,
        "corrected_historical_role_compounds": args.historical_compounds,
        "frozen_query_role_compounds": args.query_compounds,
        "frozen_candidate_targets": args.candidate_targets,
        "frozen_candidate_sequences": args.candidate_sequences,
        "frozen_evaluation_pairs": args.evaluation_pairs,
        "frozen_scaffold_audit": args.scaffold_audit,
        "frozen_homology_0_30": args.homology_0_30,
        "frozen_homology_0_50": args.homology_0_50,
        "frozen_homology_0_70": args.homology_0_70,
    }
    inputs = [input_descriptor(role, path) for role, path in input_paths.items()]

    history_fields, history = read_tsv_gz(args.historical_pairs)
    query_fields, queries = read_tsv_gz(args.scoring_queries)
    target_fields, target_rows = read_tsv_gz(args.candidate_targets)
    sequences = read_fasta(args.candidate_sequences)
    if len(history) != args.expected_historical_pair_count:
        raise ValueError(
            f"Historical relation count is {len(history)}, expected {args.expected_historical_pair_count}"
        )
    if len(queries) != args.expected_query_count:
        raise ValueError(f"Query count is {len(queries)}, expected {args.expected_query_count}")
    if "uniprot_canonical_accession" not in target_fields:
        raise ValueError("Candidate target table lacks canonical accession")
    target_ids = [row["uniprot_canonical_accession"] for row in target_rows]
    if len(target_ids) != args.expected_target_count or len(set(target_ids)) != len(target_ids):
        raise ValueError("Candidate target count or uniqueness differs from lock")
    if set(target_ids) != set(sequences):
        raise ValueError("Candidate target and sequence sets differ")
    if not {"query_id", "inchikey_full"}.issubset(query_fields):
        raise ValueError("Scoring query table lacks query_id or compound key")
    query_id_by_compound: dict[str, str] = {}
    for row in queries:
        compound = row["inchikey_full"]
        query_id = row["query_id"]
        if compound in query_id_by_compound and query_id_by_compound[compound] != query_id:
            raise ValueError("A scoring compound maps to multiple query IDs")
        query_id_by_compound[compound] = query_id

    required_history = {
        "inchikey_full",
        "uniprot_canonical_accession",
        "best_strict_evidence_tier",
    }
    if not required_history.issubset(history_fields):
        raise ValueError(f"Historical pair table lacks {sorted(required_history.difference(history_fields))}")
    helper_history = []
    for row in history:
        tier = row["best_strict_evidence_tier"]
        if tier not in LOCKED_WEIGHTS:
            raise ValueError(f"Historical pair uses a tier outside strict A/B: {tier!r}")
        helper_history.append(
            {
                "inchikey_full": row["inchikey_full"],
                "uniprot_canonical_accession": row["uniprot_canonical_accession"],
                "best_strict_evidence_tier_v1_1": tier,
            }
        )

    historical_smiles = load_role_compounds(args.historical_compounds, "historical")
    query_smiles = load_role_compounds(args.query_compounds, "query")
    required_historical_compounds = {row["inchikey_full"] for row in history}
    required_query_compounds = {row["inchikey_full"] for row in queries}
    if required_historical_compounds.difference(historical_smiles):
        raise ValueError("Historical role structure table is incomplete")
    if required_query_compounds.difference(query_smiles):
        raise ValueError("Query role structure table is incomplete")
    historical_fingerprints = morgan_fingerprints(
        {key: historical_smiles[key] for key in required_historical_compounds}, 2, 2048
    )
    query_fingerprints = morgan_fingerprints(
        {key: query_smiles[key] for key in required_query_compounds}, 2, 2048
    )

    target_index = {target: index for index, target in enumerate(target_ids)}
    train_by_compound, _ = build_train_maps(helper_history, target_index, LOCKED_WEIGHTS)
    _, sequence_matrix = build_sequence_kmer_matrix(target_ids, sequences)
    historical_target_indices = np.asarray(
        sorted({target_index[row["uniprot_canonical_accession"]] for row in history}),
        dtype=np.int32,
    )
    if historical_target_indices.size < LOCKED_TOP_K:
        raise ValueError("Historical target set is smaller than the locked top-100")
    sequence_similarity = (
        sequence_matrix @ sequence_matrix[historical_target_indices].T
    ).toarray().astype(np.float32, copy=False)
    top_columns = np.argpartition(sequence_similarity, -LOCKED_TOP_K, axis=1)[
        :, -LOCKED_TOP_K:
    ].astype(np.int32, copy=False)
    top_similarities = np.take_along_axis(sequence_similarity, top_columns, axis=1)
    target_to_historical_column = {
        int(target_idx): column for column, target_idx in enumerate(historical_target_indices)
    }

    sequence_boundary_tie_count = 0
    sequence_boundary_tie_sizes = []
    for values in sequence_similarity:
        threshold = np.partition(values, -LOCKED_TOP_K)[-LOCKED_TOP_K]
        above = int(np.count_nonzero(values > threshold))
        equal = int(np.count_nonzero(values == threshold))
        if above < LOCKED_TOP_K < above + equal:
            sequence_boundary_tie_count += 1
            sequence_boundary_tie_sizes.append(equal)

    scope_inputs = build_locked_scope_pairs(
        evaluation_pairs=args.evaluation_pairs,
        scaffold_audit=args.scaffold_audit,
        homology_paths={
            "0_30": args.homology_0_30,
            "0_50": args.homology_0_50,
            "0_70": args.homology_0_70,
        },
        query_id_by_compound=query_id_by_compound,
        target_index=target_index,
    )
    query_rows = []
    top_ranks_by_query: dict[str, np.ndarray] = {}
    exhaustive_ranks_by_query: dict[str, np.ndarray] = {}
    for query_row in queries:
        query_id = query_row["query_id"]
        query_compound = query_row["inchikey_full"]
        activation = build_activation(
            query_compound=query_compound,
            train_by_compound=train_by_compound,
            historical_fingerprints=historical_fingerprints,
            query_fingerprints=query_fingerprints,
            target_to_historical_column=target_to_historical_column,
        )
        top_scores = np.max(
            top_similarities * activation[top_columns], axis=1
        ).astype(np.float32, copy=False)
        full_scores = exhaustive_scores(
            sequence_similarity, activation, args.exhaustive_chunk_size
        )
        if np.any(top_scores > full_scores + LOCKED_SCORE_ATOL):
            raise ValueError("Top-100 score exceeds exhaustive score beyond numerical tolerance")
        allowed = np.ones(len(target_ids), dtype=bool)
        for target_idx, _ in train_by_compound.get(query_compound, []):
            allowed[target_idx] = False
        top_order, top_ranks = rank_scores(
            top_scores, allowed, query_id, target_ids, LOCKED_TIE_SALT
        )
        full_order, full_ranks = rank_scores(
            full_scores, allowed, query_id, target_ids, LOCKED_TIE_SALT
        )
        top_ranks_by_query[query_id] = top_ranks
        exhaustive_ranks_by_query[query_id] = full_ranks
        differences = full_scores[allowed].astype(np.float64) - top_scores[allowed].astype(np.float64)
        rank_differences = (
            full_ranks[allowed].astype(np.int64) - top_ranks[allowed].astype(np.int64)
        )
        query_rows.append(
            {
                "eligible_candidate_count": int(np.count_nonzero(allowed)),
                "exact_score_changed_target_count": int(np.count_nonzero(differences != 0.0)),
                "tolerance_score_changed_target_count": int(
                    np.count_nonzero(np.abs(differences) > LOCKED_SCORE_ATOL)
                ),
                "exhaustive_strictly_higher_target_count": int(
                    np.count_nonzero(differences > LOCKED_SCORE_ATOL)
                ),
                "maximum_absolute_score_error": float(np.max(np.abs(differences))),
                "mean_absolute_score_error": float(np.mean(np.abs(differences))),
                "rank_changed_target_count": int(np.count_nonzero(rank_differences != 0)),
                "rank_spearman": spearman_complete_ranks(
                    top_ranks[allowed], full_ranks[allowed]
                ),
                "mean_absolute_rank_shift": float(np.mean(np.abs(rank_differences))),
                "maximum_absolute_rank_shift": int(np.max(np.abs(rank_differences))),
                "top_k": {
                    str(k): top_set_metrics(top_order, full_order, k)
                    for k in DEFAULT_BOUNDARIES
                },
            }
        )

    fidelity = aggregate_fidelity(query_rows)
    scope_metric_rows = scope_metric_summary(
        scopes=scope_inputs,
        ranks_top100=top_ranks_by_query,
        ranks_exhaustive=exhaustive_ranks_by_query,
    )
    summary = {
        "audit_id": AUDIT_ID,
        "protocol_id": PROTOCOL_ID,
        "sensitivity_status": "post_hoc_sensitivity_only_may_not_replace_primary",
        "locked_primary_top_k": LOCKED_TOP_K,
        "locked_score_absolute_tolerance": LOCKED_SCORE_ATOL,
        "locked_tie_salt": LOCKED_TIE_SALT,
        "candidate_target_count": len(target_ids),
        "historical_target_count": int(historical_target_indices.size),
        "sequence_top100_boundary_tie_target_count": sequence_boundary_tie_count,
        "sequence_top100_boundary_tie_target_fraction": sequence_boundary_tie_count
        / len(target_ids),
        "sequence_top100_boundary_tie_size_distribution": distribution(
            sequence_boundary_tie_sizes
        ),
        "rank_and_score_fidelity": fidelity,
        "scope_metric_difference_available": bool(scope_metric_rows),
        "scope_metric_difference_row_count": len(scope_metric_rows),
        "runtime_seconds_internal_wall_clock": time.perf_counter() - started,
        "alternative_rank_ledger_written": False,
        "interpretation_boundary": (
            "The exhaustive calculation audits truncation fidelity only. It cannot replace the "
            "locked top-100 primary method regardless of performance direction."
        ),
    }
    flat = {
        "query_count": fidelity["query_count"],
        "eligible_query_target_count": fidelity["eligible_query_target_count"],
        "query_count_with_any_tolerance_score_change": fidelity[
            "query_count_with_any_tolerance_score_change"
        ],
        "tolerance_score_changed_target_count": fidelity[
            "tolerance_score_changed_target_count"
        ],
        "tolerance_score_changed_target_fraction": fidelity[
            "tolerance_score_changed_target_fraction"
        ],
        "query_count_with_any_rank_change": fidelity["query_count_with_any_rank_change"],
        "rank_changed_target_count": fidelity["rank_changed_target_count"],
        "rank_changed_target_fraction": fidelity["rank_changed_target_fraction"],
        "top10_changed_membership_query_count": fidelity["top_k"]["10"][
            "query_count_with_changed_membership"
        ],
        "top50_changed_membership_query_count": fidelity["top_k"]["50"][
            "query_count_with_changed_membership"
        ],
        "sequence_top100_boundary_tie_target_count": sequence_boundary_tie_count,
    }

    output_dir.mkdir(parents=True, exist_ok=False)
    summary_name = "top100_exhaustive_aggregate_summary.json"
    fidelity_table_name = "top100_exhaustive_fidelity.tsv"
    write_json_new(output_dir / summary_name, summary)
    write_tsv_new(output_dir / fidelity_table_name, list(flat), [flat])
    output_names = [summary_name, fidelity_table_name]
    if scope_metric_rows:
        metric_name = "top100_exhaustive_metric_differences.tsv"
        write_tsv_new(output_dir / metric_name, list(scope_metric_rows[0]), scope_metric_rows)
        output_names.append(metric_name)
    manifest = finalize_manifest(
        output_dir=output_dir,
        audit_id=AUDIT_ID,
        script_path=Path(__file__),
        inputs=inputs,
        output_names=output_names,
        extra={
            "top_k": LOCKED_TOP_K,
            "score_absolute_tolerance": LOCKED_SCORE_ATOL,
            "tie_salt": LOCKED_TIE_SALT,
            "alternative_rank_ledger_written": False,
        },
    )
    write_json_new(output_dir / "run_manifest.json", manifest)
    print(
        f"{AUDIT_ID}: wrote aggregate fidelity for {len(query_rows)} queries; "
        "no alternative rank ledger was written"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
