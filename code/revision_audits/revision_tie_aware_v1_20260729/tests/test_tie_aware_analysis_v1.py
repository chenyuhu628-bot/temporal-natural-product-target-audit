from __future__ import annotations

import itertools
import math
import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import run_tie_aware_analysis_v1 as analysis


def realized_metrics(relevance: tuple[int, ...], k: int) -> tuple[float, float, float]:
    relevant_count = sum(relevance)
    recall = sum(relevance[:k]) / relevant_count
    dcg = sum(
        value / math.log2(rank + 1)
        for rank, value in enumerate(relevance[:k], start=1)
    )
    ideal = sum(
        1.0 / math.log2(rank + 1)
        for rank in range(1, min(k, relevant_count) + 1)
    )
    hit = float(any(relevance[:k]))
    return recall, dcg / ideal, hit


class TieMetricTests(unittest.TestCase):
    def test_exact_expectation_and_bounds_match_enumeration(self) -> None:
        first_block = set(itertools.permutations((1, 0, 0)))
        second_block = set(itertools.permutations((1, 0)))
        realizations = [
            realized_metrics(first + second, 4)
            for first in first_block
            for second in second_block
        ]
        result = analysis.tie_query_metrics(
            [(1, 3, 3, 1), (4, 5, 2, 1)],
            relevant_count=2,
            salted_relevant_ranks=[1, 5],
            k=4,
        )
        recalls = [item[0] for item in realizations]
        ndcgs = [item[1] for item in realizations]
        hits = [item[2] for item in realizations]
        self.assertAlmostEqual(result["tie_expected_recall"], np.mean(recalls), places=15)
        self.assertAlmostEqual(result["tie_worst_recall"], min(recalls), places=15)
        self.assertAlmostEqual(result["tie_best_recall"], max(recalls), places=15)
        self.assertAlmostEqual(result["tie_expected_ndcg"], np.mean(ndcgs), places=15)
        self.assertAlmostEqual(result["tie_worst_ndcg"], min(ndcgs), places=15)
        self.assertAlmostEqual(result["tie_best_ndcg"], max(ndcgs), places=15)
        self.assertAlmostEqual(
            result["tie_expected_any_hit_probability"], np.mean(hits), places=15
        )
        self.assertEqual(result["membership_score_identifiable"], 1)
        self.assertEqual(result["membership_boundary_tie_dependent"], 1)
        self.assertEqual(result["membership_not_retrieved"], 0)

    def test_unique_scores_are_score_identifiable(self) -> None:
        result = analysis.tie_query_metrics(
            [(1, 1, 1, 1), (3, 3, 1, 1)],
            relevant_count=2,
            salted_relevant_ranks=[1, 3],
            k=10,
        )
        self.assertTrue(result["recall_score_identifiable"])
        self.assertTrue(result["ndcg_score_identifiable"])
        self.assertEqual(result["tie_expected_recall"], 1.0)
        self.assertEqual(result["tie_worst_recall"], 1.0)
        self.assertEqual(result["tie_best_recall"], 1.0)

    def test_all_zero_block_fractional_recall(self) -> None:
        result = analysis.tie_query_metrics(
            [(1, 100, 100, 2)],
            relevant_count=2,
            salted_relevant_ranks=[20, 80],
            k=10,
        )
        self.assertAlmostEqual(result["tie_expected_recall"], 0.1, places=15)
        self.assertEqual(result["tie_worst_recall"], 0.0)
        self.assertEqual(result["tie_best_recall"], 1.0)
        expected_hit = 1.0 - math.comb(98, 10) / math.comb(100, 10)
        self.assertAlmostEqual(
            result["tie_expected_any_hit_probability"], expected_hit, places=15
        )
        self.assertEqual(result["membership_boundary_tie_dependent"], 2)

    def test_nondiscounted_membership_can_be_fixed_while_ndcg_is_tied(self) -> None:
        result = analysis.tie_query_metrics(
            [(1, 3, 3, 1)],
            relevant_count=1,
            salted_relevant_ranks=[2],
            k=3,
        )
        self.assertTrue(result["recall_score_identifiable"])
        self.assertFalse(result["ndcg_score_identifiable"])
        self.assertEqual(result["membership_score_identifiable"], 1)


class StatisticalUtilityTests(unittest.TestCase):
    def test_zero_success_clopper_pearson_formula(self) -> None:
        for n in (19, 22):
            expected = 1.0 - 0.05 ** (1.0 / n)
            self.assertAlmostEqual(
                analysis.clopper_pearson_upper(0, n), expected, places=15
            )

    def test_derived_seed_is_stable_and_cell_specific(self) -> None:
        first = analysis.derived_seed(123, "scope", "baseline", "subset")
        second = analysis.derived_seed(123, "scope", "baseline", "subset")
        changed = analysis.derived_seed(123, "scope", "other", "subset")
        self.assertEqual(first, second)
        self.assertNotEqual(first, changed)
        self.assertGreaterEqual(first, 0)
        self.assertLess(first, 2**64)

    def test_component_partition(self) -> None:
        selected = ["synthetic_a", "synthetic_b", "synthetic_c"]
        query_pairs = {
            "synthetic_a": {"pair_a"},
            "synthetic_b": {"pair_b"},
            "synthetic_c": {"pair_c"},
        }
        pair_documents = {
            "pair_a": {"doc_shared"},
            "pair_b": {"doc_shared"},
            "pair_c": {"doc_single"},
        }
        components, document_count, edge_count = analysis.build_components(
            selected, query_pairs, pair_documents
        )
        self.assertEqual(sorted(len(component) for component in components), [1, 2])
        self.assertEqual(document_count, 2)
        self.assertEqual(edge_count, 3)

    def test_component_bootstrap_preserves_constant_estimand(self) -> None:
        selected = ["synthetic_a", "synthetic_b", "synthetic_c"]
        values = np.full((3, 2), 0.25)
        low, high, status = analysis.component_bootstrap_interval(
            values,
            selected,
            [["synthetic_a", "synthetic_b"], ["synthetic_c"]],
            seed=42,
        )
        self.assertEqual(status, "estimable_descriptive_pmid_component_sensitivity")
        np.testing.assert_allclose(low, [0.25, 0.25])
        np.testing.assert_allclose(high, [0.25, 0.25])

    def test_empty_prespecified_subset_is_not_estimable(self) -> None:
        values = np.empty((0, 2), dtype=np.float64)
        low, high, status = analysis.query_bootstrap_interval(values, seed=42)
        self.assertIsNone(low)
        self.assertIsNone(high)
        self.assertEqual(status, "not_estimable_no_queries")
        component_low, component_high, component_status = (
            analysis.component_bootstrap_interval(values, [], [], seed=43)
        )
        self.assertIsNone(component_low)
        self.assertIsNone(component_high)
        self.assertEqual(component_status, "not_estimable_component_count_lt_2")


if __name__ == "__main__":
    unittest.main()
