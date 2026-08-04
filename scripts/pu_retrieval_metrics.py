"""Tie-safe observed-candidate retrieval metrics for PU baselines."""

from __future__ import annotations

import hashlib
import math

import numpy as np


def deterministic_ties(query_id: str, target_ids: list[str], salt: str) -> np.ndarray:
    return np.fromiter(
        (int.from_bytes(hashlib.sha256(f"{salt}|{query_id}|{target}".encode("utf-8")).digest()[:8], "big") for target in target_ids),
        dtype=np.uint64,
        count=len(target_ids),
    )


def rank_scores(scores: np.ndarray, allowed: np.ndarray, query_id: str, target_ids: list[str], salt: str) -> tuple[np.ndarray, np.ndarray]:
    if scores.shape != allowed.shape or scores.size != len(target_ids):
        raise ValueError("Score, candidate-mask, and target-universe dimensions differ")
    if not np.any(allowed):
        raise ValueError("Query has no eligible candidate target")
    masked = np.where(allowed, scores, -np.inf)
    order = np.lexsort((deterministic_ties(query_id, target_ids, salt), -masked))
    order = order[allowed[order]]
    ranks = np.full(scores.shape, -1, dtype=np.int32)
    ranks[order] = np.arange(1, len(order) + 1, dtype=np.int32)
    return order, ranks


def query_metrics(positive_ranks: list[int], k_values: tuple[int, ...]) -> dict[str, float]:
    if not positive_ranks or any(rank < 1 for rank in positive_ranks):
        raise ValueError("Every observed evaluation target must have an eligible candidate rank")
    output: dict[str, float] = {f"Recall@{k}": sum(rank <= k for rank in positive_ranks) / len(positive_ranks) for k in k_values}
    for k in k_values:
        dcg = sum(1.0 / math.log2(rank + 1) for rank in positive_ranks if rank <= k)
        ideal = sum(1.0 / math.log2(index + 2) for index in range(min(k, len(positive_ranks))))
        output[f"NDCG@{k}"] = dcg / ideal if ideal else 0.0
    output["MRR"] = 1.0 / min(positive_ranks)
    return output


def macro_average(values: list[dict[str, float]]) -> dict[str, float]:
    if not values:
        raise ValueError("Cannot average zero query metric rows")
    keys = values[0]
    return {key: float(np.mean([row[key] for row in values])) for key in keys}

