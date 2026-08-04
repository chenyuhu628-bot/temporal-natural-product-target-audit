#!/usr/bin/env python3
"""Create an auditable ChEMBL 31 leakage gate for frozen future candidates.

The script is deliberately stdlib-only and never edits its two inputs or raw
data.  A C31 no-hit is a leakage-screen result, not a biological negative or
final positive label.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import os
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


FUTURE_REQUIRED = {
    "canonical_pair_key",
    "inchikey_full",
    "uniprot_canonical_accession",
    "decision",
    "label_status",
    "unrecorded_pair_policy",
}
AUDIT_REQUIRED = {
    "pair_key",
    "historical_activity_hit_status",
    "temporal_audit_action",
    "historical_activity_row_count",
    "sqlite_validated_entity_mapping_count",
    "direct_binding_asserted",
}

STATUS_RULES = {
    "no_activity_found_in_validated_chembl31_entity_pair": {
        "decision": "include_primary_future_candidate_after_C31_leakage_screen",
        "stratum": "primary_C31_validated_no_historical_activity",
        "rationale": (
            "Exact full-InChIKey and human SINGLE PROTEIN target mapping were "
            "validated in ChEMBL 31 and no activity row was found. This is a "
            "historical leakage screen only, not evidence of inactivity."
        ),
    },
    "historical_activity_recorded_in_chembl31": {
        "decision": "exclude_from_primary_future_pool_historical_C31_activity",
        "stratum": "historical_C31_activity_hit",
        "rationale": (
            "At least one exact-entity ChEMBL 31 activity row was recorded "
            "before the frozen cutoff snapshot; retain the observed row in the "
            "ledger and exclude this pair from the primary future pool."
        ),
    },
    "entity_pair_not_sqlite_validated": {
        "decision": "holdout_until_C31_entity_mapping_is_resolved",
        "stratum": "C31_entity_unresolved",
        "rationale": (
            "The required exact ChEMBL 31 compound-target entity validation was "
            "not completed; do not infer historical absence and do not include "
            "the pair in the primary future pool."
        ),
    },
}

EXPECTED_AUDIT_ACTIONS = {
    "no_activity_found_in_validated_chembl31_entity_pair":
        "retain_pending_next_evidence_audit_no_chembl31_activity_found",
    "historical_activity_recorded_in_chembl31":
        "exclude_from_future_candidate_pool_historical_chembl31_activity",
    "entity_pair_not_sqlite_validated":
        "exclude_until_entity_mapping_is_resolved",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_gzip_csv(path: Path, required: set[str], key: str) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise ValueError(f"No header in {path}")
        fields = list(reader.fieldnames)
        missing = sorted(required.difference(fields))
        if missing:
            raise ValueError(f"{path} is missing required columns: {missing}")
        rows = list(reader)
    keys = [row.get(key, "") for row in rows]
    if any(not value for value in keys):
        raise ValueError(f"Blank {key} in {path}")
    duplicate_count = len(keys) - len(set(keys))
    if duplicate_count:
        raise ValueError(f"{path} has {duplicate_count} duplicate {key} values")
    return fields, rows


def write_gzip_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as handle:
                writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="raise")
                writer.writeheader()
                writer.writerows(rows)
    os.replace(temporary, path)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def distinct_count(rows: list[dict[str, str]], column: str) -> int:
    return len({row[column] for row in rows})


def validate_audit_row(row: dict[str, str]) -> None:
    status = row["historical_activity_hit_status"]
    if status not in STATUS_RULES:
        raise ValueError(f"Unexpected historical_activity_hit_status: {status}")
    if row["temporal_audit_action"] != EXPECTED_AUDIT_ACTIONS[status]:
        raise ValueError(f"Unexpected action for {row['pair_key']}: {row['temporal_audit_action']}")
    activity_rows = int(row["historical_activity_row_count"])
    validated_entities = int(row["sqlite_validated_entity_mapping_count"])
    if status == "no_activity_found_in_validated_chembl31_entity_pair":
        valid = activity_rows == 0 and validated_entities > 0
    elif status == "historical_activity_recorded_in_chembl31":
        valid = activity_rows > 0 and validated_entities > 0
    else:
        valid = activity_rows == 0 and validated_entities == 0
    if not valid:
        raise ValueError(f"Inconsistent C31 audit counts for {row['pair_key']} ({status})")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--future-table",
        type=Path,
        default=Path("data/processed/strict_temporal_future_candidates_v1_1_pmid_verified.csv.gz"),
    )
    parser.add_argument(
        "--pair-audit",
        type=Path,
        default=Path("data/interim/chembl_31_future_candidate_historical_pair_audit.csv.gz"),
    )
    parser.add_argument("--tag", default="v1_1_pmid_verified_chembl31")
    args = parser.parse_args()
    root = args.project_root.resolve()
    future_path = args.future_table if args.future_table.is_absolute() else root / args.future_table
    audit_path = args.pair_audit if args.pair_audit.is_absolute() else root / args.pair_audit

    future_fields, future_rows = read_gzip_csv(future_path, FUTURE_REQUIRED, "canonical_pair_key")
    audit_fields, audit_rows = read_gzip_csv(audit_path, AUDIT_REQUIRED, "pair_key")
    future_by_key = {row["canonical_pair_key"]: row for row in future_rows}
    audit_by_key = {row["pair_key"]: row for row in audit_rows}
    missing_audit = sorted(set(future_by_key).difference(audit_by_key))
    if missing_audit:
        raise ValueError(f"Frozen future pairs absent from C31 audit ({len(missing_audit)}): {missing_audit[:5]}")
    invalid_future = [row["canonical_pair_key"] for row in future_rows if row["decision"] != "strict_post_cutoff_future_candidate"]
    if invalid_future:
        raise ValueError(f"Frozen future table contains non-future decisions: {invalid_future[:5]}")

    fixed_future = ["canonical_pair_key", "inchikey_full", "uniprot_canonical_accession"]
    future_rest = [field for field in future_fields if field not in fixed_future]
    audit_output = ["audit_" + field for field in audit_fields]
    gate_fields = [
        "leakage_gate_decision",
        "leakage_gate_stratum",
        "leakage_gate_rationale",
        "label_status_after_C31",
        "negative_label_emitted",
        "direct_binding_asserted_by_gate",
    ]
    output_fields = fixed_future + ["future_" + field for field in future_rest] + audit_output + gate_fields
    ledger: list[dict[str, str]] = []
    for pair_key in sorted(future_by_key):
        future = future_by_key[pair_key]
        audit = audit_by_key[pair_key]
        validate_audit_row(audit)
        if audit.get("inchikey_full") and audit["inchikey_full"] != future["inchikey_full"]:
            raise ValueError(f"InChIKey mismatch on {pair_key}")
        if audit.get("npass_uniprot_source") and audit["npass_uniprot_source"] != future["uniprot_canonical_accession"]:
            raise ValueError(f"UniProt mismatch on {pair_key}")
        rule = STATUS_RULES[audit["historical_activity_hit_status"]]
        entry = {
            "canonical_pair_key": pair_key,
            "inchikey_full": future["inchikey_full"],
            "uniprot_canonical_accession": future["uniprot_canonical_accession"],
            **{"future_" + field: future[field] for field in future_rest},
            **{"audit_" + field: audit[field] for field in audit_fields},
            "leakage_gate_decision": rule["decision"],
            "leakage_gate_stratum": rule["stratum"],
            "leakage_gate_rationale": rule["rationale"],
            "label_status_after_C31": "P1_candidate_only__P2_assay_or_paper_review_required",
            "negative_label_emitted": "false",
            "direct_binding_asserted_by_gate": "false",
        }
        ledger.append(entry)

    strata = {rule["stratum"]: [] for rule in STATUS_RULES.values()}
    for row in ledger:
        strata[row["leakage_gate_stratum"]].append(row)
    if sum(map(len, strata.values())) != len(ledger):
        raise AssertionError("Strata do not partition the frozen future table")
    if any(row["negative_label_emitted"] != "false" for row in ledger):
        raise AssertionError("This gate must never emit a negative label")

    base = f"strict_temporal_future_{args.tag}"
    paths = {
        "ledger": root / "results" / f"{base}_leakage_decision_ledger.csv.gz",
        "primary": root / "data" / "processed" / f"{base}_primary_C31_validated_no_historical_activity.csv.gz",
        "historical_hits": root / "data" / "processed" / f"{base}_historical_C31_activity_hits.csv.gz",
        "entity_unresolved": root / "data" / "processed" / f"{base}_C31_entity_unresolved.csv.gz",
        "summary": root / "results" / f"{base}_leakage_gate_summary.json",
        "manifest": root / "manifests" / f"{base}_leakage_gate_manifest.json",
    }
    write_gzip_csv(paths["ledger"], output_fields, ledger)
    write_gzip_csv(paths["primary"], output_fields, strata["primary_C31_validated_no_historical_activity"])
    write_gzip_csv(paths["historical_hits"], output_fields, strata["historical_C31_activity_hit"])
    write_gzip_csv(paths["entity_unresolved"], output_fields, strata["C31_entity_unresolved"])

    counts = {
        name: {
            "pairs": len(rows),
            "unique_compounds_full_InChIKey": distinct_count(rows, "inchikey_full") if rows else 0,
            "unique_targets_UniProt": distinct_count(rows, "uniprot_canonical_accession") if rows else 0,
        }
        for name, rows in {"full_ledger": ledger, **strata}.items()
    }
    summary = {
        "audit_name": "C31 leakage gate applied to frozen strict temporal future table",
        "run_at_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "future_table": str(future_path),
            "future_table_sha256": sha256(future_path),
            "pair_audit": str(audit_path),
            "pair_audit_sha256": sha256(audit_path),
            "frozen_future_pair_count": len(future_rows),
            "C31_audit_pair_count": len(audit_rows),
            "future_pairs_found_once_in_C31_audit": len(future_rows),
            "C31_audit_pairs_outside_frozen_future_table": len(audit_rows) - len(future_rows),
        },
        "status_counts_in_frozen_future_table": dict(sorted(Counter(row["audit_historical_activity_hit_status"] for row in ledger).items())),
        "stratum_counts": counts,
        "output_files": {name: str(path) for name, path in paths.items() if name != "manifest"},
        "warnings": [
            "All retained primary rows remain P1 candidates only; P2 assay/source-paper review is required before final positives.",
            "A C31 no-hit is not a negative label, evidence of inactivity, or proof of no pre-cutoff literature.",
            "The gate makes no direct-binding assertion and emits no unrecorded compound-target pairs as negatives.",
        ],
    }
    write_json(paths["summary"], summary)
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "script": str(Path(__file__).resolve()),
        "script_sha256": sha256(Path(__file__).resolve()),
        "input_sha256": {"future_table": sha256(future_path), "pair_audit": sha256(audit_path)},
        "output_sha256": {name: sha256(path) for name, path in paths.items() if name != "manifest"},
        "validation": {
            "future_table_unique_key": "canonical_pair_key",
            "audit_unique_key": "pair_key",
            "all_frozen_future_keys_matched_once": True,
            "all_rows_partitioned_once": True,
            "negative_labels_emitted": False,
            "direct_binding_asserted_by_gate": False,
        },
    }
    write_json(paths["manifest"], manifest)
    print(json.dumps({"stratum_counts": counts, "summary": str(paths["summary"]), "manifest": str(paths["manifest"])}, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise

