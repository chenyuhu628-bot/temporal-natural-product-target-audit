"""Unit tests for unresolved endpoint arithmetic and parsers."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from run_unresolved_bounds import (  # noqa: E402
    HISTORICAL_ACTIVITY_EXCLUSIONS,
    INITIAL_RELATIONS,
    MAX_ENDPOINT_RELATIONS,
    PRIMARY_RELATIONS,
    UNRESOLVED_RELATIONS,
    parse_bool,
    split_count,
)


class UnresolvedBoundsTests(unittest.TestCase):
    def test_endpoint_partition(self) -> None:
        self.assertEqual(
            INITIAL_RELATIONS,
            PRIMARY_RELATIONS + HISTORICAL_ACTIVITY_EXCLUSIONS + UNRESOLVED_RELATIONS,
        )
        self.assertEqual(MAX_ENDPOINT_RELATIONS, 423)

    def test_relation_bound_formula(self) -> None:
        hits = 54
        self.assertAlmostEqual(hits / MAX_ENDPOINT_RELATIONS, 54 / 423)
        self.assertAlmostEqual((hits + UNRESOLVED_RELATIONS) / MAX_ENDPOINT_RELATIONS, 119 / 423)

    def test_multivalue_count(self) -> None:
        self.assertEqual(split_count("PMID:1;PMID:2;PMID:1"), 2)
        self.assertEqual(split_count(""), 0)

    def test_boolean_parser(self) -> None:
        self.assertTrue(parse_bool("true"))
        self.assertFalse(parse_bool("False"))
        with self.assertRaises(ValueError):
            parse_bool("unknown")


if __name__ == "__main__":
    unittest.main()

