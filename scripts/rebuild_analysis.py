"""Fail-closed input gate for the full corrective analysis."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    inputs = config.get("inputs", {})
    with (ROOT / "reproduction" / "INPUT_RECONSTRUCTION_MATRIX.tsv").open(
        "r", encoding="utf-8", newline=""
    ) as handle:
        expected = {
            row["input_name"]: row["reference_sha256"]
            for row in csv.DictReader(handle, delimiter="\t")
        }
    missing = []
    mismatched = []
    for name, value in sorted(inputs.items()):
        path = Path(value)
        if not path.is_absolute():
            path = ROOT / path
        if not path.is_file():
            missing.append({"input": name, "path": value})
            continue
        observed = sha256(path)
        if name not in expected or observed != expected[name]:
            mismatched.append(
                {
                    "input": name,
                    "expected_sha256": expected.get(name),
                    "observed_sha256": observed,
                }
            )
    undeclared = sorted(set(expected) - set(inputs))
    report = {
        "status": "PASS" if not missing and not mismatched and not undeclared else "BLOCKED",
        "missing": missing,
        "mismatched": mismatched,
        "undeclared_inputs": undeclared,
        "author_side_clean_environment_reproduction": "COMPLETED",
        "independent_third_party_reproduction": "NOT YET PERFORMED",
        "note": "PASS verifies the 16 locked prerequisite hashes only; it does not grant redistribution rights or claim an independent reproduction.",
    }
    print(json.dumps(report, indent=2))
    if missing or mismatched or undeclared:
        return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
