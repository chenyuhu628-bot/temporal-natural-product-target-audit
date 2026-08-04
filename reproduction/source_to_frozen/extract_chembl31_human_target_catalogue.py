#!/usr/bin/env python3
"""Extract a historical C31 human single-protein candidate target catalogue.

The C31 target catalogue is a historical *candidate universe*, not an activity
label source. It is derived read-only from the verified ChEMBL 31 SQLite
snapshot so that future P1 target identities do not define the ranking space.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import sqlite3
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_tsv_gz(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(text, fieldnames=fields, delimiter="\t", extrasaction="raise")
                writer.writeheader()
                writer.writerows(rows)


def write_fasta(path: Path, rows: list[dict[str, str]]) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="ascii", newline="\n") as handle:
        for row in rows:
            handle.write(f">{row['uniprot_accession']}|{row['chembl_target_id']}|{row['chembl_tid']}\n")
            sequence = row["sequence"]
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index:index + 80] + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--tag", default="v1")
    args = parser.parse_args()
    root = args.project_root.resolve()
    config_path = root / "configs/chembl31_human_target_catalogue_v1.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if args.tag != config["version"]:
        raise ValueError("Tag must match the locked target-catalogue configuration")
    extraction_manifest_path = root / config["extraction_manifest_relative_path"]
    extraction_manifest = json.loads(extraction_manifest_path.read_text(encoding="utf-8"))
    if extraction_manifest.get("archive_sha256") != config["archive_sha256"]:
        raise ValueError("C31 extraction manifest archive hash does not match locked configuration")
    db_path = root / config["database_relative_path"]
    if not db_path.is_file() or db_path.stat().st_size != config["database_expected_bytes"]:
        raise ValueError("C31 database path or byte count does not match locked configuration")
    relative_db = db_path.relative_to(Path(extraction_manifest["destination"]))
    listed = {Path(item["path"]) for item in extraction_manifest.get("files", [])}
    if relative_db not in listed:
        raise ValueError("C31 database is not listed in the verified extraction manifest")
    uri = f"file:{db_path.as_posix()}?mode=ro"
    connection = sqlite3.connect(uri, uri=True)
    connection.execute("PRAGMA query_only=ON")
    query = """
        SELECT
            cs.accession AS uniprot_accession,
            td.tid AS chembl_tid,
            td.chembl_id AS chembl_target_id,
            td.pref_name AS target_pref_name,
            td.tax_id AS target_tax_id,
            td.organism AS target_organism,
            cs.sequence AS sequence,
            cs.sequence_md5sum AS sequence_md5,
            cs.db_source AS sequence_db_source,
            cs.db_version AS sequence_db_version
        FROM target_dictionary AS td
        JOIN target_components AS tc ON td.tid = tc.tid
        JOIN component_sequences AS cs ON tc.component_id = cs.component_id
        WHERE td.target_type = ?
          AND cs.organism = ?
          AND cs.accession IS NOT NULL AND cs.accession != ''
          AND cs.sequence IS NOT NULL AND cs.sequence != ''
        ORDER BY cs.accession
    """
    fetched = connection.execute(query, (config["target_type"], config["organism"])).fetchall()
    connection.close()
    if len(fetched) != config["expected_unique_accession_count"]:
        raise ValueError(f"Unexpected C31 target row count: {len(fetched)}")
    rows: list[dict[str, str]] = []
    seen: set[str] = set()
    for accession, tid, chembl_id, name, tax_id, organism, sequence, sequence_md5, db_source, db_version in fetched:
        accession = str(accession).strip()
        clean_sequence = "".join(str(sequence).split()).upper()
        if accession in seen or not clean_sequence:
            raise ValueError("C31 target catalogue must have one nonempty sequence per accession")
        seen.add(accession)
        rows.append({
            "uniprot_accession": accession,
            "chembl_tid": str(tid),
            "chembl_target_id": str(chembl_id),
            "target_pref_name": str(name or ""),
            "target_tax_id": str(tax_id or ""),
            "target_organism": str(organism or ""),
            "sequence": clean_sequence,
            "sequence_length": str(len(clean_sequence)),
            "sequence_md5": str(sequence_md5 or ""),
            "sequence_db_source": str(db_source or ""),
            "sequence_db_version": str(db_version or ""),
            "catalogue_role": "historical_C31_candidate_target_universe_only__not_activity_labels",
        })
    out_dir = root / "data/processed/chembl31_human_target_catalogue_v1"
    outputs = {
        "targets": out_dir / f"chembl31_human_single_protein_targets_{args.tag}.tsv.gz",
        "sequences": out_dir / f"chembl31_human_single_protein_targets_{args.tag}.fasta",
        "summary": root / "results" / f"chembl31_human_target_catalogue_{args.tag}_summary.json",
        "manifest": root / "manifests" / f"chembl31_human_target_catalogue_{args.tag}_manifest.json",
    }
    if any(path.exists() for path in outputs.values()):
        raise FileExistsError("Refusing to overwrite C31 target-catalogue outputs")
    target_fields = [key for key in rows[0] if key != "sequence"]
    write_tsv_gz(outputs["targets"], target_fields, [{key: value for key, value in row.items() if key != "sequence"} for row in rows])
    write_fasta(outputs["sequences"], rows)
    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "Historical C31 candidate target universe; no ChEMBL activity labels were extracted or used.",
        "target_type": config["target_type"], "organism": config["organism"], "unique_accession_count": len(rows),
        "sequence_length_min": min(int(row["sequence_length"]) for row in rows), "sequence_length_max": max(int(row["sequence_length"]) for row in rows),
        "sequence_db_versions": sorted({row["sequence_db_version"] for row in rows}),
        "outputs": {name: str(path) for name, path in outputs.items() if name != "manifest"},
        "validation": {"database_opened_read_only": True, "activity_labels_extracted": False, "future_P1_labels_used_to_construct_catalogue": False},
    }
    outputs["summary"].write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(), "script": str(Path(__file__).resolve()), "script_sha256": sha256(Path(__file__).resolve()),
        "config": {"path": str(config_path), "sha256": sha256(config_path)}, "extraction_manifest": {"path": str(extraction_manifest_path), "sha256": sha256(extraction_manifest_path)},
        "database": {"path": str(db_path), "bytes": db_path.stat().st_size},
        "outputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in outputs.items() if name != "manifest"},
        "validation": summary["validation"],
    }
    outputs["manifest"].write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

