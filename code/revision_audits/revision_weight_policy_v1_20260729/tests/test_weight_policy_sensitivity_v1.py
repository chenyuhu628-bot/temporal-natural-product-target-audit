from __future__ import annotations

import sys
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
import run_weight_policy_sensitivity_v1 as analysis


class WeightPolicyTests(unittest.TestCase):
    def test_variant_contract_and_all_equal_alias(self) -> None:
        self.assertEqual(
            [variant["tier_B_weight"] for variant in analysis.VARIANTS],
            [0.5, 0.7, 1.0],
        )
        all_equal = analysis.VARIANTS[2]
        self.assertEqual(all_equal["tier_A_weight"], 1.0)
        self.assertEqual(all_equal["tier_B_weight"], 1.0)
        self.assertEqual(all_equal["alias"], "all_equal_A1_B1")

    def test_popularity_uses_only_frozen_tier_weights(self) -> None:
        history = [
            {
                "inchikey_full": "synthetic_compound_a",
                "uniprot_canonical_accession": "synthetic_target_a",
                "best_strict_evidence_tier": analysis.TIER_A,
            },
            {
                "inchikey_full": "synthetic_compound_b",
                "uniprot_canonical_accession": "synthetic_target_a",
                "best_strict_evidence_tier": analysis.TIER_B,
            },
            {
                "inchikey_full": "synthetic_compound_b",
                "uniprot_canonical_accession": "synthetic_target_b",
                "best_strict_evidence_tier": analysis.TIER_B,
            },
        ]
        _, _, popularity, historical_indices, _ = analysis.build_history_maps(
            history, {"synthetic_target_a": 0, "synthetic_target_b": 1}
        )
        np.testing.assert_allclose(
            popularity,
            np.asarray(
                [[1.5, 0.5], [1.7, 0.7], [2.0, 1.0]], dtype=np.float32
            ),
        )
        np.testing.assert_array_equal(historical_indices, [0, 1])

    def test_chemical_activation_reuses_similarities_across_weights(self) -> None:
        similarities = [0.4, 0.8]
        compounds = ["synthetic_a", "synthetic_b"]
        relations = {
            "synthetic_a": [(0, analysis.TIER_A)],
            "synthetic_b": [(0, analysis.TIER_B), (1, analysis.TIER_B)],
        }
        activation = analysis.compute_weighted_activation(
            similarities, compounds, relations, {0: 0, 1: 1}
        )
        expected = np.asarray(
            [
                [0.4, 0.4],
                [0.56, 0.56],
                [0.8, 0.8],
            ],
            dtype=np.float32,
        )
        np.testing.assert_allclose(activation, expected)

    def test_rank_change_stats_are_exact(self) -> None:
        stats = {
            "eligible_candidate_rows": 0,
            "score_changed_candidate_count": 0,
            "rank_changed_candidate_count": 0,
            "absolute_rank_change_sum": 0,
            "maximum_absolute_rank_change": 0,
            "query_count_with_any_rank_change": 0,
            "top50_symmetric_difference_membership_count": 0,
            "query_count_with_any_top50_membership_change": 0,
            "endpoint_relation_count": 0,
            "endpoint_relation_rank_changed_count": 0,
            "endpoint_relation_top50_membership_changed_count": 0,
            "rank_permutation_blocks_checked": 0,
        }
        primary_scores = np.asarray([3.0, 2.0, 1.0], dtype=np.float32)
        candidate_scores = np.asarray([2.0, 3.0, 1.0], dtype=np.float32)
        primary_ranks = np.asarray([1, 2, 51], dtype=np.int32)
        candidate_ranks = np.asarray([2, 1, 49], dtype=np.int32)
        analysis.update_rank_change_stats(
            stats,
            candidate_scores,
            candidate_ranks,
            primary_scores,
            primary_ranks,
            np.asarray([True, True, True]),
            np.asarray([0, 2], dtype=np.int32),
        )
        self.assertEqual(stats["score_changed_candidate_count"], 2)
        self.assertEqual(stats["rank_changed_candidate_count"], 3)
        self.assertEqual(stats["absolute_rank_change_sum"], 4)
        self.assertEqual(stats["maximum_absolute_rank_change"], 2)
        self.assertEqual(stats["top50_symmetric_difference_membership_count"], 1)
        self.assertEqual(stats["endpoint_relation_rank_changed_count"], 2)
        self.assertEqual(
            stats["endpoint_relation_top50_membership_changed_count"], 1
        )

    def test_primary_variant_has_zero_change(self) -> None:
        stats = analysis.initialize_rank_change_stats()[
            ("A1_B0_7_primary", "weighted_target_popularity")
        ]
        scores = np.asarray([2.0, 1.0], dtype=np.float32)
        ranks = np.asarray([1, 2], dtype=np.int32)
        analysis.update_rank_change_stats(
            stats,
            scores,
            ranks,
            scores,
            ranks,
            np.asarray([True, True]),
            np.asarray([0], dtype=np.int32),
        )
        self.assertEqual(stats["score_changed_candidate_count"], 0)
        self.assertEqual(stats["rank_changed_candidate_count"], 0)
        self.assertEqual(stats["top50_symmetric_difference_membership_count"], 0)

    def test_scope_cardinality_is_weight_invariant(self) -> None:
        rows = {
            scope: [
                {
                    "query_id": "synthetic_query",
                    "uniprot_canonical_accession": "synthetic_target",
                }
            ]
            for scope in analysis.DISPLAY_SCOPES
        }
        output = analysis.scope_cardinality_rows(rows)
        self.assertEqual(len(output), 12)
        self.assertTrue(
            all(
                row["relation_count_change_vs_0_7"] == 0
                and row["query_count_change_vs_0_7"] == 0
                and row["target_count_change_vs_0_7"] == 0
                for row in output
            )
        )

    def test_complete_row_contract(self) -> None:
        self.assertEqual(
            analysis.EXPECTED_PRIMARY_RANK_ROWS * len(analysis.VARIANTS),
            analysis.EXPECTED_COMPUTED_SCORE_RANK_ROWS,
        )


if __name__ == "__main__":
    unittest.main()
