#!/usr/bin/env python
"""Validate the aggregate-only structure-policy sensitivity package."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


OUTPUT_DIR = Path(__file__).resolve().parents[1]
PARENT_PROTOCOL_SHA256 = "bd4902476160cc7c5cbacaf0cfd0f1a28c5300bde22232b3a2cc6c1f3c143dc2"
POLICIES = {
    "raw_primary",
    "cleanup_fragment_parent",
    "cleanup_charge_normalized",
    "cleanup_canonical_tautomer",
    "cleanup_parent_charge_tautomer",
}
BASELINES = {
    "weighted_target_popularity",
    "sequence_3mer_transfer",
    "weighted_morgan_transfer",
    "structure_sequence_pair_neighbor",
}
SCOPES = {
    "temporal_strict_ab",
    "scaffold_cold",
    "joint_scaffold_homology_0_30",
    "joint_scaffold_homology_0_50_0_70_identical",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("rt", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    checks: list[dict[str, str]] = []

    def check(condition: bool, label: str) -> None:
        checks.append({"check": label, "status": "PASS" if condition else "FAIL"})

    required = [
        "PROTOCOL.md",
        "structure_policy_summary.tsv",
        "scaffold_scope_changes.tsv",
        "rank_change_summary.tsv",
        "scope_recall_at_50.tsv",
        "blockers.tsv",
        "calibration_summary.json",
        "input_hashes.json",
        "execution_receipt.json",
        "test_results.json",
        "manifest.json",
        "scripts/build_structure_policy.py",
        "scripts/validate_structure_policy.py",
    ]
    check(all((OUTPUT_DIR / name).is_file() for name in required), "all required files exist")
    manifest = read_json(OUTPUT_DIR / "manifest.json")
    check(manifest["status"] == "BUILD_COMPLETE", "manifest build status complete")
    check(manifest["aggregate_only"] is True, "manifest declares aggregate-only outputs")
    check(manifest["identifier_bearing_outputs_retained"] == 0, "manifest retains zero identifier rows")
    check(manifest["parent_protocol_sha256"] == PARENT_PROTOCOL_SHA256, "parent protocol hash matches")
    check(manifest["protocol"]["PROTOCOL.md"] == sha256(OUTPUT_DIR / "PROTOCOL.md"), "child protocol hash matches")
    for name, record in manifest["outputs"].items():
        check(sha256(OUTPUT_DIR / name) == record["sha256"], f"manifest output hash matches: {name}")
    for name, expected in manifest["scripts"].items():
        check(sha256(OUTPUT_DIR / name) == expected, f"manifest script hash matches: {name}")

    input_hashes = read_json(OUTPUT_DIR / "input_hashes.json")
    check(input_hashes["parent_protocol_sha256"] == PARENT_PROTOCOL_SHA256, "input lock records parent protocol")
    check(input_hashes["child_protocol"]["sha256"] == sha256(OUTPUT_DIR / "PROTOCOL.md"), "input lock records child protocol")
    check(
        input_hashes["inputs"]["homology_0_50"]["sha256"]
        == input_hashes["inputs"]["homology_0_70"]["sha256"],
        "0.50 and 0.70 locked mask files have identical hashes",
    )

    structures = read_tsv(OUTPUT_DIR / "structure_policy_summary.tsv")
    check(len(structures) == 10, "structure summary has five policies by two roles")
    check({row["policy"] for row in structures} == POLICIES, "all structure policies retained")
    check({row["role"] for row in structures} == {"historical", "query"}, "roles remain separated")
    check(all(row["status"] == "complete" for row in structures), "all role-policy transformations complete")
    for row in structures:
        expected = 1726 if row["role"] == "historical" else 222
        check(int(row["input_record_count"]) == expected, f"role count matches: {row['policy']} {row['role']}")
        check(int(row["parse_failure_count"]) == 0 and int(row["transform_failure_count"]) == 0, f"no transformation failure: {row['policy']} {row['role']}")
        check(int(row["nonempty_scaffold_count"]) + int(row["empty_or_acyclic_scaffold_count"]) == expected, f"scaffold accounting closes: {row['policy']} {row['role']}")

    scope_changes = read_tsv(OUTPUT_DIR / "scaffold_scope_changes.tsv")
    check(len(scope_changes) == 5 and {row["policy"] for row in scope_changes} == POLICIES, "scope-change table covers five policies")
    raw_scope = next(row for row in scope_changes if row["policy"] == "raw_primary")
    check(
        int(raw_scope["scaffold_cold_relation_count"]) == 123
        and int(raw_scope["scaffold_cold_query_count"]) == 88
        and int(raw_scope["joint_0_30_relation_count"]) == 24
        and int(raw_scope["joint_0_30_query_count"]) == 19
        and int(raw_scope["joint_0_50_0_70_relation_count"]) == 29
        and int(raw_scope["joint_0_50_0_70_query_count"]) == 22,
        "raw scaffold and joint scopes reproduce frozen counts",
    )
    check(
        all(int(raw_scope[field]) == 0 for field in [
            "relation_entered_vs_raw_count",
            "relation_exited_vs_raw_count",
            "query_entered_vs_raw_count",
            "query_exited_vs_raw_count",
        ]),
        "raw scope has zero self-change",
    )

    ranks = read_tsv(OUTPUT_DIR / "rank_change_summary.tsv")
    check(len(ranks) == 20, "rank-change table has five policies by four baselines")
    check({row["baseline"] for row in ranks} == BASELINES, "all four baselines retained")
    check({row["policy"] for row in ranks} == POLICIES, "rank table retains all policies")
    check(all(int(row["eligible_rank_cell_count"]) == 914532 for row in ranks), "rank denominator fixed at 914,532 per baseline")
    invariant = [row for row in ranks if row["baseline"] in {"weighted_target_popularity", "sequence_3mer_transfer"}]
    check(all(int(row["score_changed_cell_count"]) == 0 and int(row["rank_changed_cell_count"]) == 0 for row in invariant), "structure-independent baselines are invariant")
    raw_rank = [row for row in ranks if row["policy"] == "raw_primary"]
    check(all(int(row["score_changed_cell_count"]) == 0 and int(row["rank_changed_cell_count"]) == 0 for row in raw_rank), "raw rank comparison has zero self-change")

    recalls = read_tsv(OUTPUT_DIR / "scope_recall_at_50.tsv")
    check(len(recalls) == 80, "Recall table has five policies by four baselines by four scopes")
    check({row["policy"] for row in recalls} == POLICIES, "Recall table retains all policies")
    check({row["baseline"] for row in recalls} == BASELINES, "Recall table retains all baselines")
    check({row["scope"] for row in recalls} == SCOPES, "Recall table uses nonduplicated scopes")
    check(all(row["status"] == "estimable" for row in recalls), "all policy-scope Recall values estimable")
    check(all(0.0 <= float(row["recall_at_50"]) <= 1.0 for row in recalls), "all Recall values lie in [0,1]")
    raw_recalls = [row for row in recalls if row["policy"] == "raw_primary"]
    check(all(abs(float(row["delta_total_vs_frozen_raw_scope"])) <= 1e-15 for row in raw_recalls), "raw Recall exactly reproduces frozen scope values")

    calibration = read_json(OUTPUT_DIR / "calibration_summary.json")
    check(calibration["status"] == "PASS", "calibration status PASS")
    check(calibration["rdkit_version"] == "2026.03.4", "RDKit version is 2026.03.4")
    check(calibration["raw_scaffold_endpoint_cells_checked"] == 358 and calibration["raw_scaffold_membership_mismatch_count"] == 0, "all raw scaffold endpoint cells match")
    rank_cal = calibration["raw_complete_rank_calibration"]
    check(rank_cal["cells_checked"] == 3658128, "all 3,658,128 raw rank cells checked")
    check(all(value == 0 for key, value in rank_cal.items() if key != "cells_checked"), "raw score and rank cells match exactly")
    check(calibration["raw_metric_cells_checked"] == 20 and calibration["raw_metric_mismatch_count"] == 0, "all 20 raw metric cells match")
    check(calibration["identity_0_50_0_70_input_hash_equal"] is True and calibration["identity_0_50_0_70_target_mask_equal"] is True, "identical homology masks verified twice")

    blockers = read_tsv(OUTPUT_DIR / "blockers.tsv")
    check(len(blockers) == 0, "no policy blocker or imputation")
    tests = read_json(OUTPUT_DIR / "test_results.json")
    check(tests["status"] == "PASS" and all(row["status"] == "PASS" for row in tests["checks"]), "all internal tests PASS")
    receipt = read_json(OUTPUT_DIR / "execution_receipt.json")
    check(receipt["status"] == "PASS", "execution receipt PASS")
    check(receipt["rdkit_version"] == "2026.03.4", "receipt records locked RDKit")
    check(receipt["identifier_bearing_outputs_retained"] == 0 and receipt["standardized_structure_rows_retained"] == 0, "receipt retains no identifiers or structures")
    check(receipt["role_separation"]["historical_and_query_maps_distinct"] is True and receipt["role_separation"]["pooled_structure_map_used"] is False, "receipt confirms role separation")

    scan_files = [OUTPUT_DIR / name for name in required if not name.startswith("scripts/")]
    patterns = {
        "full InChIKey": re.compile(r"\b[A-Z]{14}-[A-Z]{10}-[A-Z]\b"),
        "UniProt-like accession": re.compile(r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9])\b"),
        "explicit source-document identifier": re.compile(r"\b(?:PMID|pmid)\s*[:=]\s*[0-9]{4,}\b"),
        "query identifier": re.compile(r"\bquery[_-][0-9]{2,}\b", re.IGNORECASE),
        "Windows absolute path": re.compile(r"\b[A-Za-z]:\\"),
    }
    findings: list[dict[str, str]] = []
    for path in scan_files:
        text = path.read_text(encoding="utf-8")
        for label, pattern in patterns.items():
            match = pattern.search(text)
            if match:
                findings.append({"file": path.name, "finding_type": label, "match": match.group(0)})
    check(not findings, "aggregate package contains no entity identifiers or absolute paths")

    status = "PASS" if all(row["status"] == "PASS" for row in checks) else "FAIL"
    report = {
        "schema_version": "structure_policy_validation_report_v1",
        "status": status,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "check_count": len(checks),
        "checks": checks,
        "identifier_scan_findings": findings,
        "validated_manifest_sha256": sha256(OUTPUT_DIR / "manifest.json"),
    }
    (OUTPUT_DIR / "validation_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    if status != "PASS":
        failed = [row["check"] for row in checks if row["status"] == "FAIL"]
        raise SystemExit("Validation failed: " + "; ".join(failed))
    print(f"PASS: {len(checks)} checks")


if __name__ == "__main__":
    main()
