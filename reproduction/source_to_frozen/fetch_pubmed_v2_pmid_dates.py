#!/usr/bin/env python3
"""Archive PubMed dates for strict NPASS v2 A/B P1 training candidates.

The v2 snapshot is not treated as a sufficient temporal label.  This script
retrieves and retains the PubMed ESummary response for every PMID attached to
an exact A/B P1 record with a strict reviewed-human UniProt mapping.  It never
creates negative examples and it preserves date precision rather than making
up a publication day.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import time
import urllib.parse
import urllib.request
from datetime import date, datetime, timezone
from pathlib import Path

import pandas as pd

from reproducible_io import PANDAS_GZIP


EUTILS_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
MONTHS = {
    name: index
    for index, name in enumerate(
        ("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"),
        start=1,
    )
}
PRIMARY_TIERS = {"A_affinity_candidate", "B_quantitative_functional_candidate"}
STRICT_MAPPING = "strict_one_to_one_reviewed_human"


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_pubmed_date(value: str) -> tuple[str, str]:
    """Return an ISO date and its real precision; never fabricate a day."""
    text = str(value or "").strip()
    if not text:
        return "", "missing"
    match = re.search(r"\b(19|20)\d{2}\b", text)
    if not match:
        return "", "unparseable"
    year = int(match.group(0))
    rest = text[match.end() :].strip(" .;,-")
    month_match = re.search(r"\b([A-Za-z]{3,9})\b", rest)
    if not month_match:
        return f"{year:04d}", "year"
    month = MONTHS.get(month_match.group(1)[:3].casefold())
    if not month:
        return f"{year:04d}", "year"
    day_match = re.search(r"\b([0-2]?\d|3[01])\b", rest[month_match.end() :])
    if not day_match:
        return f"{year:04d}-{month:02d}", "month"
    try:
        return date(year, month, int(day_match.group(1))).isoformat(), "day"
    except ValueError:
        return f"{year:04d}-{month:02d}", "month"


def choose_date(record: dict) -> tuple[str, str, str]:
    for source, value in (
        ("epubdate", record.get("epubdate", "")),
        ("pubdate", record.get("pubdate", "")),
        ("sortpubdate", record.get("sortpubdate", "")),
    ):
        parsed, precision = parse_pubmed_date(value)
        if precision in {"day", "month", "year"}:
            return parsed, precision, source
    return "", "missing", ""


def collect_strict_v2_pmids(root: Path) -> tuple[pd.DataFrame, list[str], Path]:
    source = root / "data/processed/npass_v2_evidence_records_v1_1_uniprot_mapped.tsv.gz"
    if not source.exists():
        raise FileNotFoundError(f"Missing v2 evidence table: {source}")
    records = pd.read_csv(source, sep="\t", dtype=str, keep_default_na=False)
    eligible = (
        records["automatic_verification_level"].eq("P1_npass_raw_exact_candidate")
        & records["evidence_tier_v1_1"].isin(PRIMARY_TIERS)
        & records["mapping_status"].eq(STRICT_MAPPING)
        & records["sequence_found"].str.casefold().eq("true")
    )
    records = records.loc[eligible].copy()
    pmids = sorted(
        {
            value.strip()
            for value in records.loc[
                records["ref_id_type"].str.strip().str.upper().eq("PMID"), "ref_id"
            ]
            if value.strip().isdigit()
        }
    )
    return records, pmids, source


def fetch_batch(pmids: list[str], raw_path: Path, retries: int) -> dict:
    params = {
        "db": "pubmed",
        "retmode": "json",
        "id": ",".join(pmids),
        "tool": "npass_temporal_audit",
    }
    request = urllib.request.Request(
        EUTILS_URL + "?" + urllib.parse.urlencode(params),
        headers={"User-Agent": "NPASS-Temporal-DoubleCold-PU/0.1"},
    )
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            with urllib.request.urlopen(request, timeout=90) as response:
                raw = response.read()
            payload = json.loads(raw.decode("utf-8"))
            raw_path.write_bytes(raw)
            return payload
        except Exception as exc:  # transient NCBI throttling / transport failures
            last_error = exc
            if attempt < retries:
                time.sleep(min(30, attempt * 3))
    raise RuntimeError(f"PubMed ESummary batch failed after {retries} attempts: {last_error}")


def metadata_table(payloads: list[dict], requested_pmids: list[str]) -> pd.DataFrame:
    documents: dict[str, dict] = {}
    for payload in payloads:
        documents.update(
            {key: value for key, value in payload.get("result", {}).items() if key != "uids"}
        )
    rows = []
    for pmid in requested_pmids:
        record = documents.get(pmid, {})
        publication_date, date_precision, date_source = choose_date(record)
        article_ids = {
            entry.get("idtype"): entry.get("value")
            for entry in record.get("articleids", [])
            if isinstance(entry, dict)
        }
        rows.append(
            {
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
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--batch-size", type=int, default=200)
    parser.add_argument("--retries", type=int, default=6)
    parser.add_argument("--request-delay", type=float, default=0.38)
    args = parser.parse_args()
    if args.batch_size < 1 or args.batch_size > 200:
        parser.error("--batch-size must be between 1 and 200")
    if args.request_delay < 0.34:
        parser.error("--request-delay must be at least 0.34 seconds without an NCBI API key")

    root = args.project_root.resolve()
    records, pmids, source = collect_strict_v2_pmids(root)
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    raw_dir = root / "data/raw/pubmed" / f"v2_strict_primary_{run_id}"
    raw_dir.mkdir(parents=True, exist_ok=False)

    payloads: list[dict] = []
    for index, start in enumerate(range(0, len(pmids), args.batch_size)):
        batch = pmids[start : start + args.batch_size]
        payloads.append(fetch_batch(batch, raw_dir / f"esummary_{index:04d}.json", args.retries))
        if start + args.batch_size < len(pmids):
            time.sleep(args.request_delay)

    metadata = metadata_table(payloads, pmids)
    output = root / "data/interim/pubmed_v2_pmid_metadata.csv.gz"
    temporary = output.with_suffix(output.suffix + ".tmp")
    metadata.to_csv(temporary, index=False, compression=PANDAS_GZIP)
    os.replace(temporary, output)
    summary = {
        "created_at": utc_now(),
        "scope": "NPASS v2 exact A/B P1 records with strict reviewed-human UniProt mapping and sequence",
        "source_evidence_table": str(source),
        "source_evidence_sha256": sha256(source),
        "strict_candidate_records": int(len(records)),
        "strict_candidate_pairs": int(records["pair_key"].nunique()),
        "non_pmid_records_excluded_from_PubMed_query": int(
            (~records["ref_id_type"].str.strip().str.upper().eq("PMID")).sum()
        ),
        "unique_pmids_requested": len(pmids),
        "raw_responses": str(raw_dir),
        "metadata_file": str(output),
        "date_precision_counts": metadata["date_precision"].value_counts().to_dict(),
        "pubmed_found_count": int(metadata["found_in_pubmed"].sum()),
        "missing_from_pubmed_count": int((~metadata["found_in_pubmed"]).sum()),
    }
    write_json(root / "results/temporal_v2_pmid_screen_summary.json", summary)
    write_json(
        root / "manifests/pubmed_v2_strict_primary_query_manifest.json",
        {
            "queried_at": utc_now(),
            "endpoint": EUTILS_URL,
            "database": "pubmed",
            "batch_size": args.batch_size,
            "request_delay_seconds": args.request_delay,
            "pmid_count": len(pmids),
            "raw_response_directory": str(raw_dir),
            "scope": summary["scope"],
        },
    )
    print(root / "results/temporal_v2_pmid_screen_summary.json")


if __name__ == "__main__":
    main()
