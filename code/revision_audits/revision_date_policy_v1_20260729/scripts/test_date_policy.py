"""Unit tests for calendar interval and policy semantics."""

from __future__ import annotations

import sys
import unittest
from datetime import date
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))

from run_date_policy_sensitivity import (  # noqa: E402
    CUTOFF,
    eligible_under,
    interval_bounds,
    interval_status,
)


def row(publication_date: str, precision: str) -> dict[str, str]:
    return {
        "publication_date": publication_date,
        "date_precision": precision,
        "ref_id_type": "PMID",
        "ref_id": "123",
        "found_in_pubmed": "True",
    }


class DatePolicyTests(unittest.TestCase):
    def test_calendar_bounds(self) -> None:
        self.assertEqual(interval_bounds("2020-02", "month"), (date(2020, 2, 1), date(2020, 2, 29)))
        self.assertEqual(interval_bounds("2022", "year"), (date(2022, 1, 1), date(2022, 12, 31)))
        self.assertEqual(interval_bounds("2022-08-31", "day"), (CUTOFF, CUTOFF))

    def test_interval_states(self) -> None:
        self.assertEqual(interval_status(row("2021", "year")), "definitely_before_or_on")
        self.assertEqual(interval_status(row("2022", "year")), "crossing_cutoff")
        self.assertEqual(interval_status(row("2022-09", "month")), "definitely_after")
        self.assertEqual(interval_status(row("", "missing")), "unresolved_interval")

    def test_policy_boundary(self) -> None:
        crossing = row("2022", "year")
        self.assertFalse(eligible_under(crossing, "day_only_conservative"))
        self.assertFalse(eligible_under(crossing, "interval_certain_pre_cutoff"))
        self.assertTrue(eligible_under(crossing, "interval_earliest_bound"))

    def test_non_pmid_never_eligible(self) -> None:
        item = row("2020", "year")
        item["ref_id_type"] = "DOI"
        item["ref_id"] = "10.1/example"
        for scenario in (
            "day_only_conservative",
            "interval_certain_pre_cutoff",
            "interval_earliest_bound",
        ):
            self.assertFalse(eligible_under(item, scenario))


if __name__ == "__main__":
    unittest.main()

