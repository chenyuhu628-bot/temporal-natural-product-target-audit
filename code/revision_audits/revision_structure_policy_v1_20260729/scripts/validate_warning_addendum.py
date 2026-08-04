#!/usr/bin/env python
"""Validate the aggregate RDKit warning addendum."""

from __future__ import annotations

import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path


OUTPUT_DIR = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    summary_path = OUTPUT_DIR / "runtime_warning_summary.json"
    manifest_path = OUTPUT_DIR / "warning_manifest_addendum.json"
    summary = load(summary_path)
    manifest = load(manifest_path)
    checks = []

    def check(value: bool, label: str) -> None:
        checks.append({"check": label, "status": "PASS" if value else "FAIL"})

    check(manifest["status"] == "COMPLETE", "warning addendum manifest complete")
    check(manifest["aggregate_only"] is True, "warning addendum is aggregate-only")
    check(
        manifest["primary_manifest_sha256"] == sha256(OUTPUT_DIR / "manifest.json"),
        "warning addendum binds primary manifest",
    )
    check(
        manifest["outputs"]["runtime_warning_summary.json"] == sha256(summary_path),
        "warning summary hash matches",
    )
    for name, expected in manifest["scripts"].items():
        check(sha256(OUTPUT_DIR / name) == expected, f"warning script hash matches: {name}")
    check(summary["status"] == "AGGREGATED", "warning summary status aggregated")
    check(summary["rdkit_version"] == "2026.03.4", "warning summary RDKit version locked")
    check(
        summary["warning_event_count"]
        == summary["max_tautomers_warning_event_count"]
        + summary["max_transforms_warning_event_count"]
        + summary["kekulization_warning_event_count"]
        + summary["uncategorized_warning_event_count"],
        "warning category accounting closes",
    )
    check(summary["max_tautomers_warning_event_count"] == 56, "max-tautomer warning count matches")
    check(summary["max_transforms_warning_event_count"] == 280, "max-transform warning count matches")
    check(summary["kekulization_warning_event_count"] == 10, "kekulization warning count matches")
    check(summary["uncategorized_warning_event_count"] == 0, "no uncategorized warning events")
    check(summary["source_stdout_byte_count"] == 0, "transient stdout was empty")
    check(summary["temporary_logs_deleted_after_aggregation"] is True, "transient logs declared deleted")
    check(not (OUTPUT_DIR / "build.stderr.log").exists() and not (OUTPUT_DIR / "build.stdout.log").exists(), "transient logs actually deleted")
    text = summary_path.read_text(encoding="utf-8") + manifest_path.read_text(encoding="utf-8")
    findings = []
    for label, pattern in {
        "full InChIKey": re.compile(r"\b[A-Z]{14}-[A-Z]{10}-[A-Z]\b"),
        "absolute path": re.compile(r"\b[A-Za-z]:\\"),
        "query identifier": re.compile(r"\bquery[_-][0-9]{2,}\b", re.IGNORECASE),
    }.items():
        match = pattern.search(text)
        if match:
            findings.append({"type": label, "match": match.group(0)})
    check(not findings, "warning addendum contains no entity identifiers or absolute paths")
    status = "PASS" if all(row["status"] == "PASS" for row in checks) else "FAIL"
    report = {
        "schema_version": "structure_policy_warning_validation_v1",
        "status": status,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "check_count": len(checks),
        "checks": checks,
        "identifier_scan_findings": findings,
        "validated_warning_manifest_sha256": sha256(manifest_path),
    }
    path = OUTPUT_DIR / "warning_validation_report.json"
    if path.exists():
        raise FileExistsError("Create-once warning validation report already exists")
    path.write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    if status != "PASS":
        raise SystemExit("Warning addendum validation failed")
    print(f"PASS: {len(checks)} checks")


if __name__ == "__main__":
    main()
