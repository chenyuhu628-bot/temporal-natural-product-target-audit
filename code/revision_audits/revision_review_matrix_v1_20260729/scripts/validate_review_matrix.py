#!/usr/bin/env python
"""Validate the aggregate-only reviewer matrix and similarity audit."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


EXPECTED_PROTOCOL_SHA256 = "bd4902476160cc7c5cbacaf0cfd0f1a28c5300bde22232b3a2cc6c1f3c143dc2"


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


def close(actual: str, expected: float, tolerance: float = 2e-6) -> bool:
    return abs(float(actual) - expected) <= tolerance


def main() -> None:
    script_path = Path(__file__).resolve()
    output_dir = script_path.parents[1]
    checks: list[dict[str, str]] = []

    def check(condition: bool, label: str) -> None:
        checks.append({"check": label, "status": "PASS" if condition else "FAIL"})

    required = [
        "reviewer_comment_matrix.tsv",
        "reviewer_comment_matrix.md",
        "max_similarity_summary.tsv",
        "correction_record.json",
        "input_hashes.json",
        "execution_receipt.json",
        "manifest.json",
    ]
    check(all((output_dir / name).is_file() for name in required), "all required outputs exist")

    manifest = read_json(output_dir / "manifest.json")
    check(manifest["aggregate_only"] is True, "manifest declares aggregate-only package")
    check(
        manifest["status"] == "BUILD_COMPLETE_CORRECTED",
        "manifest build status records completed correction",
    )
    check(
        manifest["identifier_bearing_outputs_retained"] == 0,
        "manifest declares zero retained identifier-bearing outputs",
    )
    check(
        manifest["frozen_protocol_sha256"] == EXPECTED_PROTOCOL_SHA256,
        "manifest protocol hash matches frozen protocol",
    )
    for name, record in manifest["outputs"].items():
        check(
            sha256(output_dir / name) == record["sha256"],
            f"manifest output hash matches: {name}",
        )
    for name, expected_hash in manifest["scripts"].items():
        check(
            sha256(output_dir / name) == expected_hash,
            f"manifest script hash matches: {name}",
        )

    input_hashes = read_json(output_dir / "input_hashes.json")
    check(
        input_hashes["frozen_protocol_sha256"] == EXPECTED_PROTOCOL_SHA256,
        "input lock records frozen protocol hash",
    )
    check(
        input_hashes["inputs"]["protocol"]["sha256"] == EXPECTED_PROTOCOL_SHA256,
        "protocol input itself has frozen hash",
    )
    check(
        input_hashes["inputs"]["reviewer_report_attachment"]["sha256"]
        == "c40be1c0ee7e3d077da6d0bb89476a94c8b894508c239633826e7a4e10e24888",
        "review report attachment hash matches audited report",
    )
    check(
        input_hashes["inputs"]["homology_0_50"]["sha256"]
        == input_hashes["inputs"]["homology_0_70"]["sha256"],
        "locked 0.50 and 0.70 homology inputs have identical hashes",
    )

    matrix = read_tsv(output_dir / "reviewer_comment_matrix.tsv")
    check(len(matrix) == 13, "review matrix contains 13 classified comments")
    check(
        {row["comment_id"] for row in matrix}
        == {f"C{i}" for i in range(1, 4)}
        | {f"I{i}" for i in range(1, 6)}
        | {f"S{i}" for i in range(1, 6)},
        "review matrix comment identifiers are complete",
    )
    check(
        all(row["action_class"] and row["v4_action"] for row in matrix),
        "every review comment has an action class and v4 disposition",
    )
    check(
        sum(row["external_human_gate"] == "yes" for row in matrix) == 1,
        "categorical external-human submission gate is isolated",
    )

    receipt = read_json(output_dir / "execution_receipt.json")
    check(receipt["status"] == "PASS", "execution receipt status is PASS")
    check(receipt["reviewer_comment_count"] == 13, "receipt reviewer-comment count matches")
    correction = read_json(output_dir / "correction_record.json")
    check(
        correction["status"] == "CORRECTED_AND_SUPERSEDED"
        and correction["superseded_manifest_sha256"]
        == "ad5a9692e93049004e9d0436e70ee1122ef53463e5667bdbe4f54ff9461897f1",
        "correction record transparently supersedes the initial erroneous manifest",
    )
    check(
        "dividing by 100" in correction["correction"]
        and correction["identifier_bearing_outputs_retained"] == 0,
        "correction record states the pident unit conversion",
    )
    check(
        receipt["correction_history"]["superseded_manifest_sha256"]
        == correction["superseded_manifest_sha256"],
        "execution receipt carries the same correction history",
    )
    check(receipt["word_counts"]["abstract"] == 234, "abstract word count is 234")
    check(
        receipt["word_counts"]["pre_reference_manuscript"] == 7871,
        "pre-reference manuscript word count is 7,871",
    )
    date_audit = receipt["date_precision_audit"]
    check(date_audit["locked_row_count"] == 20647, "date ledger has 20,647 rows")
    check(
        date_audit["non_day_precision_row_count"] == 6570,
        "date ledger has 6,570 non-day rows",
    )
    check(
        date_audit["non_day_interval_definitely_before_or_on_cutoff"] == 6570
        and date_audit["non_day_interval_crosses_cutoff"] == 0
        and date_audit["non_day_interval_definitely_after_cutoff"] == 0,
        "every non-day interval is certainly no later than cutoff",
    )
    check(
        date_audit["interval_certain_eligible_row_count"] == 20455,
        "interval-certain policy admits 20,455 rows",
    )
    check(
        date_audit["historical_relation_count_under_interval_certain_policy"] == 4990,
        "interval-certain policy preserves 4,990 historical relations",
    )
    check(
        date_audit["day_only_changes_reverted_under_interval_certain_policy"] == 141
        and date_audit[
            "day_only_changes_persisting_vs_old_under_interval_certain_policy"
        ]
        == 25,
        "interval-certain policy reverses 141 of 166 day-only tier changes",
    )
    scopes = receipt["scope_counts"]
    check(
        scopes["temporal_strict_ab"]
        == {"relations": 358, "queries": 222, "targets": 156},
        "temporal scope counts match",
    )
    check(
        scopes["scaffold_cold"]
        == {"relations": 123, "queries": 88, "targets": 70},
        "scaffold-cold scope counts match",
    )
    check(
        scopes["joint_scaffold_homology_0_30"]
        == {"relations": 24, "queries": 19, "targets": 17},
        "joint 0.30 scope counts match",
    )
    check(
        scopes["joint_scaffold_homology_0_50_0_70_identical"]
        == {"relations": 29, "queries": 22, "targets": 21},
        "merged 0.50 and 0.70 scope counts match",
    )
    check(
        scopes["identity_0_50_0_70_relation_mask_equal"] is True,
        "receipt records exact 0.50 and 0.70 mask equality",
    )
    source = receipt["source_dependence"]
    check(
        source["source_document_count"] == 124
        and source["query_source_component_count"] == 95
        and source["largest_query_component_count"] == 51,
        "source-document dependence counts match",
    )

    summary = read_tsv(output_dir / "max_similarity_summary.tsv")
    check(len(summary) == 24, "similarity summary contains 24 aggregate rows")
    check(
        {row["scope"] for row in summary}
        == {
            "temporal_strict_ab",
            "scaffold_cold",
            "joint_scaffold_homology_0_30",
            "joint_scaffold_homology_0_50_0_70_identical",
        },
        "similarity summary uses four nonduplicated scopes",
    )
    check(
        not any(row["scope"] in {"joint_scaffold_homology_0_50", "joint_scaffold_homology_0_70"} for row in summary),
        "identical 0.50 and 0.70 scopes are not duplicated",
    )
    by_key = {
        (row["similarity_family"], row["analysis_unit"], row["scope"]): row
        for row in summary
    }

    morgan_temporal_query = by_key[
        ("morgan_radius2_2048_tanimoto", "query", "temporal_strict_ab")
    ]
    check(
        morgan_temporal_query["n_total"] == "222"
        and close(morgan_temporal_query["median"], 0.6523199023)
        and close(morgan_temporal_query["mean"], 0.6702981, 3e-6)
        and morgan_temporal_query["n_at_upper_bound"] == "65",
        "temporal query-level Morgan distribution matches",
    )
    morgan_joint = by_key[
        (
            "morgan_radius2_2048_tanimoto",
            "query",
            "joint_scaffold_homology_0_50_0_70_identical",
        )
    ]
    check(
        morgan_joint["n_total"] == "22"
        and close(morgan_joint["median"], 0.4436961, 3e-6)
        and close(morgan_joint["mean"], 0.4748473, 3e-6),
        "merged joint-scope Morgan distribution matches",
    )
    kmer_temporal_relation = by_key[
        (
            "native_sequence_3mer_tfidf_cosine",
            "relation_weighted",
            "temporal_strict_ab",
        )
    ]
    check(
        kmer_temporal_relation["n_total"] == "358"
        and close(kmer_temporal_relation["min"], 0.1075003, 3e-6)
        and close(kmer_temporal_relation["mean"], 0.7748950, 3e-6),
        "temporal relation-level native 3-mer distribution matches",
    )
    kmer_joint_query = by_key[
        (
            "native_sequence_3mer_tfidf_cosine",
            "query",
            "joint_scaffold_homology_0_30",
        )
    ]
    check(
        kmer_joint_query["n_total"] == "19"
        and close(kmer_joint_query["median"], 0.2044200, 3e-6)
        and close(kmer_joint_query["mean"], 0.2412888, 3e-6),
        "joint 0.30 query-level native 3-mer distribution matches",
    )
    mmseq_temporal_relation = by_key[
        (
            "mmseqs2_detected_alignment_identity",
            "relation_weighted",
            "temporal_strict_ab",
        )
    ]
    check(
        mmseq_temporal_relation["n_total"] == "358"
        and mmseq_temporal_relation["n_observed"] == "265"
        and mmseq_temporal_relation["n_no_detected_alignment"] == "93",
        "temporal relation-level MMseqs2 detected/censored counts match",
    )
    check(
        close(mmseq_temporal_relation["min"], 0.224, 1e-9)
        and close(mmseq_temporal_relation["mean"], 0.9687396226, 3e-9)
        and close(mmseq_temporal_relation["max"], 1.0, 1e-9),
        "MMseqs2 pident percentage is converted to fractional identity",
    )
    mmseq_joint_030 = by_key[
        (
            "mmseqs2_detected_alignment_identity",
            "relation_weighted",
            "joint_scaffold_homology_0_30",
        )
    ]
    check(
        mmseq_joint_030["n_observed"] == "1"
        and close(mmseq_joint_030["min"], 0.252, 1e-9)
        and close(mmseq_joint_030["max"], 0.252, 1e-9),
        "joint 0.30 detected identity is 0.252 rather than an upper-bound artifact",
    )
    mmseq_joint_rel = by_key[
        (
            "mmseqs2_detected_alignment_identity",
            "relation_weighted",
            "joint_scaffold_homology_0_50_0_70_identical",
        )
    ]
    check(
        mmseq_joint_rel["n_observed"] == "6"
        and close(mmseq_joint_rel["min"], 0.252, 1e-9)
        and close(mmseq_joint_rel["q1"], 0.345, 1e-9)
        and close(mmseq_joint_rel["median"], 0.366, 1e-9)
        and close(mmseq_joint_rel["q3"], 0.378, 1e-9)
        and close(mmseq_joint_rel["max"], 0.465, 1e-9)
        and close(mmseq_joint_rel["mean"], 0.3615, 1e-9),
        "merged joint-scope detected relation identities match audited range",
    )
    mmseq_joint_query = by_key[
        (
            "mmseqs2_detected_alignment_identity",
            "query_maximum",
            "joint_scaffold_homology_0_50_0_70_identical",
        )
    ]
    check(
        mmseq_joint_query["n_total"] == "22"
        and mmseq_joint_query["n_observed"] == "5"
        and mmseq_joint_query["n_no_detected_alignment"] == "17"
        and mmseq_joint_query[
            "n_scope_units_with_all_relations_detected"
        ]
        == "4",
        "merged joint-scope MMseqs2 detected/censored query counts match",
    )
    check(
        close(mmseq_joint_query["min"], 0.252, 1e-9)
        and close(mmseq_joint_query["q1"], 0.342, 1e-9)
        and close(mmseq_joint_query["median"], 0.378, 1e-9)
        and close(mmseq_joint_query["q3"], 0.378, 1e-9)
        and close(mmseq_joint_query["max"], 0.465, 1e-9)
        and close(mmseq_joint_query["mean"], 0.363, 1e-9),
        "merged joint-scope detected query identities match audited range",
    )
    check(
        all(
            0.0 <= float(row["min"]) <= float(row["max"]) <= 1.0
            for row in summary
            if row["similarity_family"]
            == "mmseqs2_detected_alignment_identity"
            and row["n_observed"] != "0"
        ),
        "all reported MMseqs2 identities use fractional scale",
    )
    check(
        all(
            "censored" in row["interpretation"].lower()
            for row in summary
            if row["similarity_family"]
            == "mmseqs2_detected_alignment_identity"
        ),
        "every MMseqs2 row labels non-detections as censored",
    )

    scan_files = [
        output_dir / "reviewer_comment_matrix.tsv",
        output_dir / "reviewer_comment_matrix.md",
        output_dir / "max_similarity_summary.tsv",
        output_dir / "correction_record.json",
        output_dir / "input_hashes.json",
        output_dir / "execution_receipt.json",
        output_dir / "manifest.json",
    ]
    patterns = {
        "full InChIKey": re.compile(r"\b[A-Z]{14}-[A-Z]{10}-[A-Z]\b"),
        "UniProt-like accession": re.compile(
            r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z0-9]{3}[0-9])\b"
        ),
        "explicit source-document identifier": re.compile(
            r"\b(?:PMID|pmid)\s*[:=]\s*[0-9]{4,}\b"
        ),
        "query identifier": re.compile(r"\bquery[_-][0-9]{2,}\b", re.IGNORECASE),
        "Windows absolute path": re.compile(r"\b[A-Za-z]:\\"),
    }
    findings: list[dict[str, str]] = []
    for path in scan_files:
        text = path.read_text(encoding="utf-8")
        for label, pattern in patterns.items():
            match = pattern.search(text)
            if match:
                findings.append(
                    {"file": path.name, "finding_type": label, "match": match.group(0)}
                )
    check(not findings, "aggregate outputs contain no entity identifiers or absolute paths")

    status = "PASS" if all(item["status"] == "PASS" for item in checks) else "FAIL"
    report = {
        "schema_version": "revision_review_matrix_validation_v1",
        "status": status,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "check_count": len(checks),
        "checks": checks,
        "identifier_scan_findings": findings,
        "validated_manifest_sha256": sha256(output_dir / "manifest.json"),
    }
    (output_dir / "validation_report.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    if status != "PASS":
        failed = [item["check"] for item in checks if item["status"] == "FAIL"]
        raise SystemExit("Validation failed: " + "; ".join(failed))
    print(f"PASS: {len(checks)} checks")


if __name__ == "__main__":
    main()
