#!/usr/bin/env python3
"""Build conservative date-verified NPASS temporal pair candidate tables.

The script deliberately writes no unrecorded molecular--target pair and never
uses an absence from a database as a negative label.  Its primary candidates
are exact P1 NPASS records with A/B evidence tiers and strict current UniProt
mapping.  A date-verified historical table is produced only when a separately
archived v2 PMID metadata table is available.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Iterator

from reproducible_io import deterministic_gzip_text


STRICT_MAPPING = "strict_one_to_one_reviewed_human"
PRIMARY_TIERS = {"A_affinity_candidate", "B_quantitative_functional_candidate"}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_rows(path: Path, delimiter: str) -> Iterator[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter=delimiter)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def joined(values: set[str]) -> str:
    return ";".join(sorted(value for value in values if value))


def new_state(inchikey: str, accession: str) -> dict:
    return {
        "inchikey_full": inchikey,
        "uniprot_canonical_accession": accession,
        "source_versions": set(), "source_pair_keys": set(), "source_np_ids": set(),
        "source_target_ids": set(), "primary_evidence_tiers": set(),
        "primary_activity_types": set(), "primary_ref_ids": set(),
        "v2_all_record_count": 0, "v3_all_record_count": 0,
        "v2_primary_record_count": 0, "v3_primary_record_count": 0,
        "v2_strict_primary_record_count": 0, "v3_strict_primary_record_count": 0,
        "v2_pre_cutoff_record_count": 0, "v2_date_reasons": set(),
        "v3_cross_version_statuses": set(), "v3_temporal_screen_statuses": set(),
    }


def canonical_pair(row: dict[str, str]) -> tuple[str, str, str]:
    inchikey = row.get("inchikey_full", "").strip()
    accession = row.get("uniprot_canonical_accession", "").strip()
    if not accession:
        accession = row.get("uniprot_raw", "").strip()
    return f"{inchikey}|{accession}", inchikey, accession


def is_primary(row: dict[str, str]) -> bool:
    return (
        row.get("automatic_verification_level", "") == "P1_npass_raw_exact_candidate"
        and row.get("evidence_tier_v1_1", "") in PRIMARY_TIERS
    )


def has_strict_entity(row: dict[str, str]) -> bool:
    return (
        bool(row.get("inchikey_full", "").strip())
        and bool(row.get("uniprot_canonical_accession", "").strip())
        and row.get("mapping_status", "") == STRICT_MAPPING
        and row.get("sequence_found", "").strip().casefold() == "true"
    )


def load_pmid_dates(path: Path) -> dict[str, tuple[str, str]]:
    if not path.exists():
        return {}
    dates: dict[str, tuple[str, str]] = {}
    for row in read_rows(path, ","):
        dates[row.get("pmid", "").strip()] = (
            row.get("publication_date", "").strip(), row.get("date_precision", "").strip()
        )
    return dates


def pre_cutoff_status(row: dict[str, str], dates: dict[str, tuple[str, str]], cutoff: date) -> str:
    if row.get("ref_id_type", "").strip().upper() != "PMID" or not row.get("ref_id", "").strip().isdigit():
        return "v2_primary_reference_not_pmid"
    publication_date, precision = dates.get(row["ref_id"].strip(), ("", "missing"))
    if precision != "day" or not publication_date:
        return "v2_primary_pmid_missing_or_not_day_precise"
    return "pre_cutoff" if date.fromisoformat(publication_date) <= cutoff else "v2_primary_pmid_after_cutoff"


def ingest(path: Path, version: str, states: dict[str, dict], dates: dict[str, tuple[str, str]], cutoff: date, v2_dates_available: bool) -> None:
    for row in read_rows(path, "\t"):
        key, inchikey, accession = canonical_pair(row)
        state = states.setdefault(key, new_state(inchikey, accession))
        state["source_versions"].add(version)
        state["source_pair_keys"].add(row.get("pair_key", ""))
        state["source_np_ids"].add(row.get("source_np_id", ""))
        state["source_target_ids"].add(row.get("source_target_id", ""))
        state[f"{version}_all_record_count"] += 1
        if version == "v3":
            state["v3_cross_version_statuses"].add(row.get("cross_version_status", ""))
            state["v3_temporal_screen_statuses"].add(row.get("temporal_screen_status", ""))
        if not is_primary(row):
            continue
        state[f"{version}_primary_record_count"] += 1
        state["primary_evidence_tiers"].add(row.get("evidence_tier_v1_1", ""))
        state["primary_activity_types"].add(row.get("activity_type", ""))
        state["primary_ref_ids"].add(f"{row.get('ref_id_type', '')}:{row.get('ref_id', '')}")
        if not has_strict_entity(row):
            continue
        state[f"{version}_strict_primary_record_count"] += 1
        if version == "v2":
            if not v2_dates_available:
                state["v2_date_reasons"].add("v2_pmid_metadata_not_available")
            else:
                status = pre_cutoff_status(row, dates, cutoff)
                if status == "pre_cutoff":
                    state["v2_pre_cutoff_record_count"] += 1
                else:
                    state["v2_date_reasons"].add(status)


def materialise(key: str, state: dict, v2_dates_available: bool) -> dict[str, str]:
    reasons: list[str] = []
    has_primary = state["v2_primary_record_count"] + state["v3_primary_record_count"] > 0
    has_strict_primary = state["v2_strict_primary_record_count"] + state["v3_strict_primary_record_count"] > 0
    has_v2 = state["v2_all_record_count"] > 0
    has_v3 = state["v3_all_record_count"] > 0
    training = state["v2_pre_cutoff_record_count"] > 0
    future = (
        not has_v2 and state["v3_strict_primary_record_count"] > 0 and has_v3
        and state["v3_cross_version_statuses"] == {"v3_only"}
        and state["v3_temporal_screen_statuses"] == {"future_candidate_pmid_only"}
    )
    if training:
        decision = "strict_pre_cutoff_training_candidate"
        rationale = "at_least_one_primary_A_or_B_P1_record_has_day_precise_PMID_on_or_before_cutoff"
    elif future:
        decision = "strict_post_cutoff_future_candidate"
        rationale = "absent_from_v2_and_all_v3_records_are_v3_only_PMID_screened_future_candidates"
    elif state["v2_strict_primary_record_count"]:
        decision = "historical_snapshot_candidate_not_date_verified_for_training"
        rationale = "v2_snapshot_primary_candidate_retained_separately_until_record_level_PubMed_dates_are_archived"
        reasons.extend(sorted(state["v2_date_reasons"]))
    else:
        decision = "excluded_from_strict_temporal_pair_tables"
        rationale = "does_not_meet_all_primary_evidence_entity_and_temporal_conditions"
    if not has_primary:
        reasons.append("no_primary_A_or_B_P1_exact_record")
    elif not has_strict_primary:
        reasons.append("no_primary_record_with_strict_reviewed_human_sequence_mapping")
    if has_v3 and not future and not training:
        if has_v2:
            reasons.append("historical_pair_present_in_v2")
        if state["v3_cross_version_statuses"] != {"v3_only"}:
            reasons.append("v3_pair_not_exclusively_v3_only")
        reasons.extend(f"v3_temporal_screen_{value or 'missing'}" for value in sorted(state["v3_temporal_screen_statuses"]) if value != "future_candidate_pmid_only")
    if has_v2 and not v2_dates_available and state["v2_strict_primary_record_count"]:
        reasons.append("v2_pmid_metadata_not_available")
    return {
        "canonical_pair_key": key,
        "inchikey_full": state["inchikey_full"],
        "uniprot_canonical_accession": state["uniprot_canonical_accession"],
        "source_versions": joined(state["source_versions"]),
        "source_pair_keys": joined(state["source_pair_keys"]),
        "source_np_ids": joined(state["source_np_ids"]),
        "source_target_ids": joined(state["source_target_ids"]),
        "v2_all_record_count": str(state["v2_all_record_count"]),
        "v3_all_record_count": str(state["v3_all_record_count"]),
        "v2_primary_A_B_P1_record_count": str(state["v2_primary_record_count"]),
        "v3_primary_A_B_P1_record_count": str(state["v3_primary_record_count"]),
        "v2_strict_entity_primary_record_count": str(state["v2_strict_primary_record_count"]),
        "v3_strict_entity_primary_record_count": str(state["v3_strict_primary_record_count"]),
        "v2_day_precise_pre_cutoff_primary_record_count": str(state["v2_pre_cutoff_record_count"]),
        "primary_evidence_tiers": joined(state["primary_evidence_tiers"]),
        "primary_activity_types": joined(state["primary_activity_types"]),
        "primary_references": joined(state["primary_ref_ids"]),
        "v3_cross_version_statuses": joined(state["v3_cross_version_statuses"]),
        "v3_temporal_screen_statuses": joined(state["v3_temporal_screen_statuses"]),
        "decision": decision,
        "decision_rationale": rationale,
        "exclusion_or_holdout_reasons": ";".join(sorted(set(reasons))),
        "label_status": "P1_candidate_only__P2_assay_or_paper_review_required",
        "unrecorded_pair_policy": "not_represented_and_never_emitted_as_negative",
    }


def write_table(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else [
        "canonical_pair_key", "inchikey_full", "uniprot_canonical_accession", "source_versions", "source_pair_keys",
        "source_np_ids", "source_target_ids", "v2_all_record_count", "v3_all_record_count",
        "v2_primary_A_B_P1_record_count", "v3_primary_A_B_P1_record_count", "v2_strict_entity_primary_record_count",
        "v3_strict_entity_primary_record_count", "v2_day_precise_pre_cutoff_primary_record_count", "primary_evidence_tiers",
        "primary_activity_types", "primary_references", "v3_cross_version_statuses", "v3_temporal_screen_statuses",
        "decision", "decision_rationale", "exclusion_or_holdout_reasons", "label_status", "unrecorded_pair_policy",
    ]
    with deterministic_gzip_text(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--cutoff", default="2022-08-31")
    parser.add_argument("--tag", default="v1")
    parser.add_argument("--v2-pubmed-metadata", type=Path, default=None,
                        help="Optional CSV/GZ with pmid, publication_date, and date_precision columns.")
    args = parser.parse_args()
    root = args.project_root.resolve()
    cutoff = date.fromisoformat(args.cutoff)
    v2_path = root / "data/processed/npass_v2_evidence_records_v1_1_uniprot_mapped.tsv.gz"
    v3_path = root / "data/processed/npass_v3_evidence_records_v1_1_uniprot_mapped.tsv.gz"
    metadata_path = args.v2_pubmed_metadata or root / "data/interim/pubmed_v2_pmid_metadata.csv.gz"
    metadata_path = metadata_path.resolve()
    v2_dates_available = metadata_path.exists()
    states: dict[str, dict] = {}
    pmid_dates = load_pmid_dates(metadata_path)
    ingest(v2_path, "v2", states, pmid_dates, cutoff, v2_dates_available)
    ingest(v3_path, "v3", states, pmid_dates, cutoff, v2_dates_available)
    rows = [materialise(key, states[key], v2_dates_available) for key in sorted(states)]
    training = [row for row in rows if row["decision"] == "strict_pre_cutoff_training_candidate"]
    future = [row for row in rows if row["decision"] == "strict_post_cutoff_future_candidate"]
    snapshot = [row for row in rows if int(row["v2_strict_entity_primary_record_count"]) > 0]
    tag = args.tag
    paths = {
        "training": root / f"data/processed/strict_temporal_training_candidates_{tag}.csv.gz",
        "future": root / f"data/processed/strict_temporal_future_candidates_{tag}.csv.gz",
        "historical_snapshot": root / f"data/processed/historical_snapshot_primary_candidates_{tag}.csv.gz",
        "ledger": root / f"results/strict_temporal_pair_decision_ledger_{tag}.csv.gz",
        "summary": root / f"results/strict_temporal_pair_tables_{tag}_summary.json",
        "manifest": root / f"manifests/strict_temporal_pair_tables_{tag}_manifest.json",
    }
    for name, subset in (("training", training), ("future", future), ("historical_snapshot", snapshot), ("ledger", rows)):
        write_table(paths[name], subset)
    decisions = Counter(row["decision"] for row in rows)
    reason_counts = Counter(reason for row in rows for reason in row["exclusion_or_holdout_reasons"].split(";") if reason)
    summary = {
        "created_at": utc_now(), "cutoff": cutoff.isoformat(), "primary_evidence_tiers": sorted(PRIMARY_TIERS),
        "strict_entity_rule": "full InChIKey + strict one-to-one reviewed human UniProt mapping + retrieved sequence",
        "training_rule": "v2 exact A/B P1 record with a day-precise PMID date on or before cutoff",
        "future_rule": "no v2 record, exact A/B P1 strict-entity record, and all v3 records labelled v3_only + future_candidate_pmid_only",
        "v2_pmid_metadata": {"path": str(metadata_path), "available": v2_dates_available, "rows": len(pmid_dates)},
        "unrecorded_pair_policy": "Unrecorded compound-target pairs are outside every output table and are never emitted as negative labels.",
        "input_sha256": {str(v2_path): sha256(v2_path), str(v3_path): sha256(v3_path)},
        "counts": {"observed_canonical_pairs": len(rows), "strict_pre_cutoff_training_candidates": len(training),
                   "strict_post_cutoff_future_candidates": len(future), "historical_snapshot_primary_candidates": len(snapshot),
                   "decision_counts": dict(sorted(decisions.items())), "reason_counts": dict(sorted(reason_counts.items()))},
        "outputs": {name: str(path) for name, path in paths.items()},
        "warning": "P1 candidates are not final positive labels; P2 assay/source-paper review remains required.",
    }
    for target in (paths["summary"], paths["manifest"]):
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(paths["summary"])


if __name__ == "__main__":
    main()
