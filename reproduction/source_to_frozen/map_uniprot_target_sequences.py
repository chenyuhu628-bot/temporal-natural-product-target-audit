#!/usr/bin/env python3
"""Map NPASS UniProt accessions, retain source responses, and fetch strict FASTA.

The mapping is deliberately conservative: one-to-many, deleted, unreviewed, or
non-human source accessions are retained in audit output but excluded from the
strict sequence/model set.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from reproducible_io import PANDAS_GZIP


API_ROOT = "https://rest.uniprot.org"


def now() -> str:
    return datetime.now(timezone.utc).isoformat()


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def request(url: str, *, data: bytes | None = None, headers: dict[str, str] | None = None, retries: int = 5) -> tuple[bytes, dict[str, str], int]:
    last_error: Exception | None = None
    for attempt in range(1, retries + 1):
        try:
            req = urllib.request.Request(url, data=data, headers={"User-Agent": "NPASS-Temporal-DoubleCold-PU/0.1", **(headers or {})})
            with urllib.request.urlopen(req, timeout=120) as response:
                return response.read(), dict(response.headers.items()), response.status
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as exc:
            last_error = exc
            if isinstance(exc, urllib.error.HTTPError) and exc.code not in {429, 500, 502, 503, 504}:
                raise
            time.sleep(min(30, 2 ** attempt))
    raise RuntimeError(f"UniProt request failed after {retries} attempts: {last_error}")


class NoRedirect(urllib.request.HTTPRedirectHandler):
    def redirect_request(self, req, fp, code, msg, headers, newurl):  # type: ignore[override]
        return None


def submit_mapping(accessions: list[str], raw_dir: Path) -> tuple[str, dict[str, str]]:
    payload = urllib.parse.urlencode({"ids": ",".join(accessions), "from": "UniProtKB_AC-ID", "to": "UniProtKB"}).encode()
    raw, headers, _ = request(f"{API_ROOT}/idmapping/run", data=payload, headers={"Content-Type": "application/x-www-form-urlencoded"})
    (raw_dir / "mapping_submission_response.json").write_bytes(raw)
    job_id = json.loads(raw.decode("utf-8"))["jobId"]
    return job_id, headers


def wait_for_job(job_id: str, raw_dir: Path) -> None:
    opener = urllib.request.build_opener(NoRedirect())
    url = f"{API_ROOT}/idmapping/status/{job_id}"
    for attempt in range(1, 61):
        req = urllib.request.Request(url, headers={"User-Agent": "NPASS-Temporal-DoubleCold-PU/0.1"})
        try:
            with opener.open(req, timeout=60) as response:
                body = response.read()
                (raw_dir / f"mapping_status_{attempt:03d}.json").write_bytes(body)
                status = json.loads(body.decode("utf-8")).get("jobStatus", "")
                if status not in {"RUNNING", "NEW"}:
                    return
        except urllib.error.HTTPError as exc:
            if exc.code == 303:
                (raw_dir / f"mapping_status_{attempt:03d}.headers.json").write_text(json.dumps(dict(exc.headers.items()), indent=2) + "\n", encoding="utf-8")
                return
            if exc.code not in {429, 500, 502, 503, 504}:
                raise
        time.sleep(2)
    raise TimeoutError(f"UniProt ID-mapping job {job_id} did not complete in two minutes")


def next_link(header: str) -> str | None:
    match = re.search(r"<([^>]+)>;\s*rel=\"next\"", header or "")
    return match.group(1) if match else None


def fetch_mapping_pages(job_id: str, raw_dir: Path) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    fields = "accession,reviewed,organism_id,length,protein_name,gene_primary,sequence_version,date_sequence_modified,version"
    url = f"{API_ROOT}/idmapping/uniprotkb/results/{job_id}?" + urllib.parse.urlencode({"format": "tsv", "fields": fields, "size": "500"})
    frames, header_records, index = [], [], 0
    while url:
        raw, headers, status = request(url)
        path = raw_dir / f"mapping_results_{index:04d}.tsv"
        path.write_bytes(raw)
        frames.append(pd.read_csv(io.BytesIO(raw), sep="\t", dtype=str, keep_default_na=False))
        header_records.append({"url": url, "status": status, "sha256": digest(raw), "headers": headers})
        url = next_link(headers.get("Link", ""))
        index += 1
    return pd.concat(frames, ignore_index=True), header_records


def fetch_strict_fasta(job_id: str, raw_dir: Path) -> tuple[dict[str, str], dict[str, str], bytes]:
    query = "(reviewed:true) AND (organism_id:9606)"
    url = f"{API_ROOT}/idmapping/uniprotkb/results/stream/{job_id}?" + urllib.parse.urlencode({"format": "fasta", "compressed": "true", "query": query})
    raw, headers, _ = request(url)
    path = raw_dir / "strict_reviewed_human_sequences.fasta.gz"
    path.write_bytes(raw)
    plain = gzip.decompress(raw) if raw.startswith(b"\x1f\x8b") else raw
    if not plain.startswith(b">"):
        raise RuntimeError("UniProt FASTA response did not begin with a FASTA header")
    sequences, accession, chunks = {}, "", []
    for line in plain.decode("utf-8").splitlines():
        if line.startswith(">"):
            if accession:
                sequences[accession] = "".join(chunks)
            fields = line[1:].split("|", 2)
            accession, chunks = (fields[1] if len(fields) > 1 else line[1:].split()[0]), []
        else:
            chunks.append(line.strip())
    if accession:
        sequences[accession] = "".join(chunks)
    return sequences, headers, raw


def column(frame: pd.DataFrame, choices: list[str]) -> str:
    for choice in choices:
        if choice in frame.columns:
            return choice
    raise KeyError(f"Missing expected mapping column; available={list(frame.columns)}")


def audit_mapping(accessions: list[str], mapping: pd.DataFrame, sequences: dict[str, str]) -> pd.DataFrame:
    source_col = column(mapping, ["From"])
    entry_col = column(mapping, ["Entry"])
    reviewed_col = column(mapping, ["Reviewed"])
    organism_col = column(mapping, ["Organism (ID)", "Organism ID"])
    rows = []
    for source in accessions:
        group = mapping.loc[mapping[source_col].eq(source)].copy()
        entries = sorted(set(group[entry_col]) - {""})
        if not entries:
            status, canonical = "unmapped_or_deleted", ""
        elif len(entries) != 1:
            status, canonical = "one_to_many_mapping", ""
        else:
            canonical = entries[0]
            target = group.loc[group[entry_col].eq(canonical)].iloc[0]
            reviewed = str(target[reviewed_col]).strip().casefold() == "reviewed"
            human = str(target[organism_col]).strip() == "9606"
            status = "strict_one_to_one_reviewed_human" if reviewed and human else "not_reviewed_human"
        sequence = sequences.get(canonical, "") if status == "strict_one_to_one_reviewed_human" else ""
        rows.append({
            "uniprot_source_accession": source, "uniprot_canonical_accession": canonical, "mapping_status": status,
            "mapping_row_count": len(group), "mapping_entry_count": len(entries), "sequence_found": bool(sequence),
            "sequence_length": len(sequence), "sequence_md5": hashlib.md5(sequence.encode("ascii")).hexdigest() if sequence else "",
        })
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.project_root.resolve()
    accessions = set()
    for version in ("v2", "v3"):
        values = pd.read_csv(root / "data/interim" / f"npass_{version}_human_single_protein_records.tsv.gz", sep="\t", usecols=["uniprot_raw"], dtype=str, keep_default_na=False)
        accessions.update(values["uniprot_raw"].str.strip())
    accessions = sorted(accessions - {""})
    run_id = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    raw_dir = root / "data/raw/uniprot" / f"id_mapping_{run_id}"
    raw_dir.mkdir(parents=True, exist_ok=False)
    (raw_dir / "input_accessions.txt").write_text("\n".join(accessions) + "\n", encoding="utf-8")

    job_id, submission_headers = submit_mapping(accessions, raw_dir)
    wait_for_job(job_id, raw_dir)
    mapping, page_headers = fetch_mapping_pages(job_id, raw_dir)
    sequences, fasta_headers, fasta_raw = fetch_strict_fasta(job_id, raw_dir)
    mapping_audit = audit_mapping(accessions, mapping, sequences)
    audit_path = root / "data/interim/uniprot_npass_target_mapping_audit.csv.gz"
    mapping_audit.to_csv(audit_path, index=False, compression=PANDAS_GZIP)

    for version in ("v2", "v3"):
        source = root / "data/processed" / f"npass_{version}_evidence_records_v1_1.tsv.gz"
        records = pd.read_csv(source, sep="\t", dtype=str, keep_default_na=False)
        records = records.merge(mapping_audit, left_on="uniprot_raw", right_on="uniprot_source_accession", how="left", validate="many_to_one")
        records.to_csv(root / "data/processed" / f"npass_{version}_evidence_records_v1_1_uniprot_mapped.tsv.gz", sep="\t", index=False, compression=PANDAS_GZIP)

    summary = {
        "queried_at": now(), "id_mapping_job_id": job_id, "input_accessions": len(accessions),
        "raw_response_directory": str(raw_dir), "mapping_audit": str(audit_path),
        "mapping_status_counts": mapping_audit["mapping_status"].value_counts().to_dict(),
        "strict_sequence_count": int(mapping_audit["sequence_found"].sum()), "submission_headers": submission_headers,
        "fasta_headers": fasta_headers, "fasta_sha256": digest(fasta_raw), "mapping_page_headers": page_headers,
        "warning": "Current UniProt mapping/sequence retrieval. Final temporal work should also archive the relevant historical release."
    }
    (root / "results/uniprot_target_mapping_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(root / "results/uniprot_target_mapping_summary.json")


if __name__ == "__main__":
    main()
