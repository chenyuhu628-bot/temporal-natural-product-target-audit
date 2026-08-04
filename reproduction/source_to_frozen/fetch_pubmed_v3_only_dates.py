#!/usr/bin/env python3
"""Resolve PubMed dates for NPASS v3-only human-protein candidate records.

The script preserves each E-utilities response under data/raw/pubmed/ and
creates a derived date table and a conservative pair-level temporal audit. It
does not assert that a quantitative endpoint is direct binding.
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import time
import urllib.parse
import urllib.request
from datetime import date, datetime, timezone
from pathlib import Path

import pandas as pd

from reproducible_io import PANDAS_GZIP


EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
MONTHS = {name: index for index, name in enumerate(("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"), start=1)}


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def parse_pubmed_date(value: str) -> tuple[str, str]:
    """Return ISO date plus precision, never inventing a day that was absent."""
    text = str(value or "").strip()
    if not text:
        return "", "missing"
    match = re.search(r"\b(19|20)\d{2}\b", text)
    if not match:
        return "", "unparseable"
    year = int(match.group(0))
    rest = text[match.end():].strip(" .;,-")
    month_match = re.search(r"\b([A-Za-z]{3,9})\b", rest)
    if not month_match:
        return f"{year:04d}", "year"
    month = MONTHS.get(month_match.group(1)[:3].casefold())
    if not month:
        return f"{year:04d}", "year"
    after_month = rest[month_match.end():]
    day_match = re.search(r"\b([0-2]?\d|3[01])\b", after_month)
    if not day_match:
        return f"{year:04d}-{month:02d}", "month"
    try:
        return date(year, month, int(day_match.group(1))).isoformat(), "day"
    except ValueError:
        return f"{year:04d}-{month:02d}", "month"


def choose_date(record: dict) -> tuple[str, str, str]:
    """Prefer an electronic publication date, then issue publication date."""
    for source, value in (("epubdate", record.get("epubdate", "")), ("pubdate", record.get("pubdate", "")), ("sortpubdate", record.get("sortpubdate", ""))):
        parsed, precision = parse_pubmed_date(value)
        if precision in {"day", "month", "year"}:
            return parsed, precision, source
    return "", "missing", ""


def fetch_batch(pmids: list[str], raw_path: Path, retries: int) -> dict:
    parameters = {"db": "pubmed", "retmode": "json", "id": ",".join(pmids), "tool": "npass_temporal_audit"}
    url = EUTILS_URL + "?" + urllib.parse.urlencode(parameters)
    request = urllib.request.Request(url, headers={"User-Agent": "NPASS-Temporal-DoubleCold-PU/0.1"})
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=60) as response:
                raw = response.read()
            raw_path.write_bytes(raw)
            return json.loads(raw.decode("utf-8"))
        except Exception as exc:  # NCBI can briefly throttle anonymous clients.
            last_error = exc
            time.sleep(min(20, attempt * 2))
    raise RuntimeError(f"E-utilities batch failed: {last_error}")


def collect_v3_only_pmids(root: Path) -> tuple[pd.DataFrame, list[str]]:
    membership = pd.read_csv(root / "results/npass_cross_version_pair_membership.csv.gz", dtype=str, keep_default_na=False)
    v3_only = set(membership.loc[membership["cross_version_status"].eq("v3_only"), "pair_key"])
    records = pd.read_csv(root / "data/interim/npass_v3_human_single_protein_records.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
    records = records.loc[records["pair_key"].isin(v3_only)].copy()
    pmids = sorted({value.strip() for value in records.loc[records["ref_id_type"].str.strip().str.upper().eq("PMID"), "ref_id"] if value.strip().isdigit()})
    return records, pmids


def metadata_table(batch_payloads: list[dict], requested_pmids: list[str]) -> pd.DataFrame:
    rows = []
    documents: dict[str, dict] = {}
    for payload in batch_payloads:
        documents.update({key: value for key, value in payload.get("result", {}).items() if key != "uids"})
    for pmid in requested_pmids:
        record = documents.get(pmid, {})
        publication_date, date_precision, date_source = choose_date(record)
        article_ids = {entry.get("idtype"): entry.get("value") for entry in record.get("articleids", []) if isinstance(entry, dict)}
        rows.append({
            "pmid": pmid,
            "found_in_pubmed": bool(record),
            "pubdate_raw": record.get("pubdate", ""),
            "epubdate_raw": record.get("epubdate", ""),
            "sortpubdate_raw": record.get("sortpubdate", ""),
            "publication_date": publication_date,
            "date_precision": date_precision,
            "date_source": date_source,
            "doi": article_ids.get("doi", ""),
            "journal": record.get("fulljournalname", ""),
            "title": record.get("title", ""),
        })
    return pd.DataFrame(rows)


def classify_pairs(records: pd.DataFrame, metadata: pd.DataFrame, cutoff: date) -> pd.DataFrame:
    evidence = records.copy()
    evidence["is_pmid"] = evidence["ref_id_type"].str.strip().str.upper().eq("PMID") & evidence["ref_id"].str.strip().str.isdigit()
    evidence = evidence.merge(metadata[["pmid", "publication_date", "date_precision"]], left_on="ref_id", right_on="pmid", how="left")

    rows = []
    for pair_key, group in evidence.groupby("pair_key", sort=False):
        pmid_rows = group.loc[group["is_pmid"]]
        non_pmid = int((~group["is_pmid"]).sum())
        parsed_dates = []
        unresolved = False
        for _, record in pmid_rows.iterrows():
            value, precision = record["publication_date"], record["date_precision"]
            if precision == "day" and value:
                parsed_dates.append(date.fromisoformat(value))
            else:
                unresolved = True
        if pmid_rows.empty:
            status = "no_pmid_reference"
        elif non_pmid:
            status = "non_pmid_reference_present"
        elif unresolved:
            status = "pmid_date_not_day_precise"
        elif any(value <= cutoff for value in parsed_dates):
            status = "pre_cutoff_or_same_day"
        else:
            status = "future_candidate_pmid_only"
        rows.append({
            "pair_key": pair_key,
            "v3_records": int(len(group)),
            "pmid_records": int(len(pmid_rows)),
            "unique_pmids": int(pmid_rows["ref_id"].nunique()),
            "non_pmid_records": non_pmid,
            "earliest_day_precise_pmid_date": min(parsed_dates).isoformat() if parsed_dates else "",
            "latest_day_precise_pmid_date": max(parsed_dates).isoformat() if parsed_dates else "",
            "temporal_screen_status": status,
        })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--cutoff", default="2022-08-31", help="strict temporal cutoff, YYYY-MM-DD")
    parser.add_argument("--batch-size", type=int, default=200)
    parser.add_argument("--retries", type=int, default=5)
    args = parser.parse_args()
    cutoff = date.fromisoformat(args.cutoff)
    root = args.project_root.resolve()
    records, pmids = collect_v3_only_pmids(root)
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    raw_dir = root / "data/raw/pubmed" / f"v3_only_{run_id}"
    raw_dir.mkdir(parents=True, exist_ok=False)

    payloads = []
    for start in range(0, len(pmids), args.batch_size):
        batch = pmids[start:start + args.batch_size]
        payloads.append(fetch_batch(batch, raw_dir / f"esummary_{start // args.batch_size:04d}.json", args.retries))
        if start + args.batch_size < len(pmids):
            time.sleep(0.35)  # E-utilities' unauthenticated request-rate guidance.

    metadata = metadata_table(payloads, pmids)
    metadata_path = root / "data/interim/pubmed_v3_only_pmid_metadata.csv.gz"
    metadata.to_csv(metadata_path, index=False, compression=PANDAS_GZIP)
    pair_audit = classify_pairs(records, metadata, cutoff)
    pair_path = root / "results/temporal_v3_only_pmid_screen.csv.gz"
    pair_audit.to_csv(pair_path, index=False, compression=PANDAS_GZIP)
    summary = {
        "created_at": utc_now(),
        "cutoff": cutoff.isoformat(),
        "scope": "NPASS v3-only human-single-protein candidate pairs; PMID-date screen only",
        "v3_only_records": int(len(records)),
        "v3_only_pairs": int(records["pair_key"].nunique()),
        "unique_pmids_requested": len(pmids),
        "raw_responses": str(raw_dir),
        "metadata_file": str(metadata_path),
        "pair_screen_file": str(pair_path),
        "date_precision_counts": metadata["date_precision"].value_counts().to_dict(),
        "pair_status_counts": pair_audit["temporal_screen_status"].value_counts().to_dict(),
    }
    write_json(root / "results/temporal_v3_only_pmid_screen_summary.json", summary)
    write_json(root / "manifests/pubmed_v3_only_query_manifest.json", {
        "queried_at": utc_now(), "endpoint": EUTILS_URL, "database": "pubmed", "batch_size": args.batch_size,
        "pmid_count": len(pmids), "raw_response_directory": str(raw_dir), "cutoff": cutoff.isoformat(),
    })
    print(root / "results/temporal_v3_only_pmid_screen_summary.json")


if __name__ == "__main__":
    main()
