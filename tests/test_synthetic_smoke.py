"""Synthetic smoke tests requiring no third-party source data."""

from __future__ import annotations

import csv
import hashlib
import json
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code/audits"))
sys.path.insert(0, str(ROOT / "scripts"))

from audit_rank_score_structure import boundary_metrics
from pu_retrieval_metrics import query_metrics, rank_scores


class SyntheticSmokeTest(unittest.TestCase):
    def test_fixture_contract(self) -> None:
        with (ROOT / "tests/fixtures/ranking_fixture.tsv").open("r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        expected = json.loads((ROOT / "tests/fixtures/expected_metrics.json").read_text(encoding="utf-8"))
        self.assertEqual(len(rows), expected["expected_row_count"])
        self.assertEqual(len({row["query_label"] for row in rows}), expected["expected_query_count"])
        recalls = []
        for query in sorted({row["query_label"] for row in rows}):
            block = [row for row in rows if row["query_label"] == query]
            ordered = sorted(block, key=lambda row: (-float(row["score"]), row["target_label"]))
            relevant = sum(int(row["observed_candidate"]) for row in ordered)
            recalls.append(sum(int(row["observed_candidate"]) for row in ordered[:2]) / relevant)
        self.assertEqual(float(np.mean(recalls)), expected["expected_macro_recall_at_2"])

    def test_locked_helpers_on_synthetic_values(self) -> None:
        _, ranks = rank_scores(
            np.asarray([0.9, 0.5, 0.5]),
            np.asarray([True, True, True]),
            "synthetic_query",
            ["one", "two", "three"],
            "synthetic_salt",
        )
        self.assertEqual(sorted(ranks.tolist()), [1, 2, 3])
        metrics = query_metrics([1, 3], (10,))
        self.assertAlmostEqual(metrics["Recall@10"], 1.0)
        tie = boundary_metrics([{"rank": 1, "score": 1.0}, {"rank": 2, "score": 0.5}, {"rank": 3, "score": 0.5}], 2)
        self.assertTrue(tie["crosses_boundary"])

    def test_authoritative_chain_hashes(self) -> None:
        chain = json.loads((ROOT / "manifests/CURRENT_EXECUTION_CHAIN.json").read_text(encoding="utf-8"))
        for step in chain["authoritative_steps"]:
            observed = hashlib.sha256((ROOT / step["package_path"]).read_bytes()).hexdigest()
            self.assertEqual(observed, step["sha256"], step["step_id"])


if __name__ == "__main__":
    unittest.main()
