#!/usr/bin/env python3
"""Read-only historical ChEMBL 31 activity audit for strict NPASS future candidates.

The script never downloads, extracts, writes to, or opens the ChEMBL database
in read-write mode.  It validates the precomputed *full* InChIKey / source
UniProt mappings against the frozen SQLite database, then extracts only
historical activities for validated molecule--target combinations.

It deliberately does not call any activity a direct-binding observation.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import sqlite3
import sys
import time
from collections import defaultdict
from datetime import UTC, datetime
from pathlib import Path
from typing import Any, Iterable

from reproducible_io import deterministic_gzip_text


SQLITE_HEADER = b"SQLite format 3\x00"
DEFAULT_MAPPING = (
    "data/interim/chembl_31_future_candidate_entity_mapping.csv.gz"
)
EXPECTED_ARCHIVE = "data/raw/chembl/chembl_31/chembl_31_sqlite.tar.gz"


def utc_now() -> str:
    return datetime.now(UTC).replace(microsecond=0).isoformat()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("r", encoding="utf-8", newline="")


def split_ids(value: str | None) -> list[str]:
    return sorted(
        {
            part.strip().upper()
            for part in (value or "").split(";")
            if part.strip()
        }
    )


def true_value(value: str | bool | None) -> bool:
    return str(value or "").strip().casefold() in {"true", "1", "yes"}


def qident(name: str) -> str:
    """Quote a SQLite identifier selected from database metadata."""
    return '"' + name.replace('"', '""') + '"'


def chunks(values: list[Any], size: int) -> Iterable[list[Any]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def sqlite_file_ok(path: Path) -> bool:
    try:
        with path.open("rb") as handle:
            return handle.read(len(SQLITE_HEADER)) == SQLITE_HEADER
    except OSError:
        return False


def candidate_database_files(base: Path) -> list[Path]:
    matches: list[Path] = []
    if not base.exists():
        return matches
    for suffix in ("*.db", "*.sqlite", "*.sqlite3"):
        for path in base.rglob(suffix):
            if ".parts" not in path.parts and sqlite_file_ok(path):
                matches.append(path.resolve())
    return sorted(set(matches))


def database_readiness(
    base: Path, explicit: Path | None
) -> tuple[Path | None, dict[str, Any]]:
    archive = base / Path(EXPECTED_ARCHIVE).name
    archive_size = archive.stat().st_size if archive.exists() else None
    candidates = candidate_database_files(base)
    state: dict[str, Any] = {
        "checked_at_utc": utc_now(),
        "database_search_root": str(base),
        "archive_path": str(archive),
        "archive_exists": archive.exists(),
        "archive_bytes": archive_size,
        "database_candidates": [str(item) for item in candidates],
    }
    if explicit is not None:
        db = explicit.expanduser().resolve()
        state["explicit_database_path"] = str(db)
        if sqlite_file_ok(db):
            state["status"] = "ready"
            return db, state
        state["status"] = "explicit_database_missing_or_not_sqlite"
        return None, state
    if len(candidates) == 1:
        state["status"] = "ready"
        return candidates[0], state
    if len(candidates) == 0:
        state["status"] = "no_unpacked_sqlite_database_found"
    else:
        state["status"] = "multiple_sqlite_databases_found"
    return None, state


def wait_for_database(
    base: Path, explicit: Path | None, wait_seconds: int
) -> tuple[Path | None, dict[str, Any]]:
    deadline = time.monotonic() + wait_seconds
    while True:
        db, readiness = database_readiness(base, explicit)
        if db is not None or time.monotonic() >= deadline:
            readiness["wait_seconds_requested"] = wait_seconds
            return db, readiness
        time.sleep(min(10, max(0.1, deadline - time.monotonic())))


def read_candidates(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    if not path.exists():
        raise FileNotFoundError(f"Candidate entity-mapping file is missing: {path}")
    required = {
        "pair_key",
        "inchikey_full",
        "npass_uniprot_source",
        "best_evidence_tier_v1_1",
        "npass_record_count",
        "chembl_compound_ids",
        "chembl_target_ids",
        "both_entities_exactly_mapped",
    }
    with open_text(path) as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        missing = sorted(required.difference(fields))
        if missing:
            raise ValueError(f"Candidate mapping lacks required columns: {', '.join(missing)}")
        rows = []
        for row in reader:
            clean = {key: (value or "").strip() for key, value in row.items()}
            clean["inchikey_full"] = clean["inchikey_full"].upper()
            clean["npass_uniprot_source"] = clean["npass_uniprot_source"].upper()
            rows.append(clean)
    if not rows:
        raise ValueError(f"Candidate mapping has no rows: {path}")
    return rows, fields


def expand_preliminary_mappings(candidates: list[dict[str, str]]) -> list[dict[str, str]]:
    """Expand semicolon-separated ChEMBL identifiers without transferring labels."""
    expanded: list[dict[str, str]] = []
    seen: set[tuple[str, str, str]] = set()
    for row in candidates:
        if not true_value(row.get("both_entities_exactly_mapped")):
            continue
        compound_ids = split_ids(row.get("chembl_compound_ids"))
        target_ids = split_ids(row.get("chembl_target_ids"))
        for compound_id in compound_ids:
            for target_id in target_ids:
                key = (row["pair_key"], compound_id, target_id)
                if key in seen:
                    continue
                seen.add(key)
                expanded.append(
                    {
                        "candidate_mapping_id": f"map_{len(expanded) + 1:06d}",
                        "pair_key": row["pair_key"],
                        "inchikey_full": row["inchikey_full"],
                        "npass_uniprot_source": row["npass_uniprot_source"],
                        "npass_uniprot_canonical": row.get("npass_uniprot_canonical", ""),
                        "best_evidence_tier_v1_1": row["best_evidence_tier_v1_1"],
                        "npass_record_count": row["npass_record_count"],
                        "chembl_compound_id": compound_id,
                        "chembl_target_id": target_id,
                    }
                )
    return expanded


def connect_read_only(database: Path) -> sqlite3.Connection:
    uri = database.resolve().as_uri() + "?mode=ro"
    connection = sqlite3.connect(uri, uri=True)
    connection.row_factory = sqlite3.Row
    connection.execute("PRAGMA query_only = ON")
    connection.execute("PRAGMA temp_store = MEMORY")
    return connection


def table_names(connection: sqlite3.Connection) -> set[str]:
    return {
        row["name"]
        for row in connection.execute(
            "SELECT name FROM sqlite_master WHERE type = 'table'"
        )
    }


def column_names(connection: sqlite3.Connection, table: str) -> set[str]:
    return {row["name"] for row in connection.execute(f"PRAGMA table_info({qident(table)})")}


def ensure_schema(connection: sqlite3.Connection) -> dict[str, set[str]]:
    required_tables = {
        "activities",
        "assays",
        "target_dictionary",
        "molecule_dictionary",
        "compound_structures",
        "target_components",
        "component_sequences",
    }
    present = table_names(connection)
    missing_tables = sorted(required_tables.difference(present))
    if missing_tables:
        raise RuntimeError(
            "The SQLite file does not expose the expected ChEMBL activity schema; "
            f"missing tables: {', '.join(missing_tables)}"
        )
    schema = {table: column_names(connection, table) for table in present}
    required_columns = {
        "activities": {"activity_id", "assay_id", "molregno"},
        "assays": {"assay_id", "tid"},
        "target_dictionary": {"tid", "chembl_id", "target_type"},
        "molecule_dictionary": {"molregno", "chembl_id"},
        "compound_structures": {"molregno", "standard_inchi_key"},
        "target_components": {"tid", "component_id"},
        "component_sequences": {"component_id", "accession"},
    }
    problems = []
    for table, columns in required_columns.items():
        absent = sorted(columns.difference(schema[table]))
        if absent:
            problems.append(f"{table}: {', '.join(absent)}")
    if problems:
        raise RuntimeError(
            "The SQLite schema is not compatible with this audit: " + "; ".join(problems)
        )
    if not ({"tax_id", "organism"} & schema["target_dictionary"]):
        raise RuntimeError(
            "target_dictionary has neither tax_id nor organism; human-target validation is impossible"
        )
    return schema


def select_optional(alias: str, table_columns: set[str], column: str, output: str) -> str:
    if column in table_columns:
        return f"{alias}.{qident(column)} AS {qident(output)}"
    return f"NULL AS {qident(output)}"


def query_compound_keys(
    connection: sqlite3.Connection, compound_ids: list[str]
) -> dict[str, set[str]]:
    observed: dict[str, set[str]] = defaultdict(set)
    for batch in chunks(compound_ids, 400):
        placeholders = ",".join("?" for _ in batch)
        query = f"""
            SELECT md.chembl_id AS chembl_compound_id,
                   UPPER(cs.standard_inchi_key) AS standard_inchi_key
            FROM molecule_dictionary AS md
            JOIN compound_structures AS cs ON cs.molregno = md.molregno
            WHERE md.chembl_id IN ({placeholders})
        """
        for row in connection.execute(query, batch):
            if row["standard_inchi_key"]:
                observed[row["chembl_compound_id"]].add(row["standard_inchi_key"])
    return observed


def query_target_metadata(
    connection: sqlite3.Connection, target_ids: list[str], schema: dict[str, set[str]]
) -> dict[str, dict[str, Any]]:
    metadata: dict[str, dict[str, Any]] = {}
    organism_expr = select_optional("td", schema["target_dictionary"], "organism", "target_organism")
    tax_expr = select_optional("td", schema["target_dictionary"], "tax_id", "target_tax_id")
    name_expr = select_optional("td", schema["target_dictionary"], "pref_name", "target_pref_name")
    for batch in chunks(target_ids, 400):
        placeholders = ",".join("?" for _ in batch)
        query = f"""
            SELECT td.chembl_id AS chembl_target_id,
                   td.tid AS tid,
                   td.target_type AS target_type,
                   {name_expr},
                   {organism_expr},
                   {tax_expr},
                   UPPER(cs.accession) AS component_accession
            FROM target_dictionary AS td
            LEFT JOIN target_components AS tc ON tc.tid = td.tid
            LEFT JOIN component_sequences AS cs ON cs.component_id = tc.component_id
            WHERE td.chembl_id IN ({placeholders})
        """
        for row in connection.execute(query, batch):
            target_id = row["chembl_target_id"]
            current = metadata.setdefault(
                target_id,
                {
                    "target_type": row["target_type"],
                    "target_pref_name": row["target_pref_name"],
                    "target_organism": row["target_organism"],
                    "target_tax_id": row["target_tax_id"],
                    "component_accessions": set(),
                },
            )
            if row["component_accession"]:
                current["component_accessions"].add(row["component_accession"])
    return metadata


def is_human_target(metadata: dict[str, Any]) -> bool:
    tax_id = str(metadata.get("target_tax_id") or "").strip()
    organism = str(metadata.get("target_organism") or "").strip().casefold()
    return tax_id == "9606" or organism == "homo sapiens"


def validate_entities(
    mappings: list[dict[str, str]],
    compound_keys: dict[str, set[str]],
    target_metadata: dict[str, dict[str, Any]],
) -> tuple[list[dict[str, Any]], dict[str, dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    mapping_lookup: dict[str, dict[str, Any]] = {}
    for mapping in mappings:
        compound_id = mapping["chembl_compound_id"]
        target_id = mapping["chembl_target_id"]
        metadata = target_metadata.get(target_id, {})
        observed_keys = compound_keys.get(compound_id, set())
        component_accessions = sorted(metadata.get("component_accessions", set()))
        compound_ok = mapping["inchikey_full"] in observed_keys
        target_type_ok = str(metadata.get("target_type") or "").strip().upper() == "SINGLE PROTEIN"
        human_ok = is_human_target(metadata)
        accession_ok = mapping["npass_uniprot_source"] in set(component_accessions)
        target_ok = target_type_ok and human_ok and accession_ok
        entity_ok = compound_ok and target_ok
        if entity_ok:
            status = "validated_full_inchikey_human_single_protein_source_uniprot"
        elif not compound_ok and not target_ok:
            status = "compound_and_target_validation_failed"
        elif not compound_ok:
            status = "compound_full_inchikey_validation_failed"
        else:
            status = "human_single_protein_source_uniprot_validation_failed"
        row: dict[str, Any] = {
            **mapping,
            "sqlite_compound_observed_inchikeys": ";".join(sorted(observed_keys)),
            "sqlite_compound_full_inchikey_exact": compound_ok,
            "sqlite_target_type": metadata.get("target_type", ""),
            "sqlite_target_pref_name": metadata.get("target_pref_name", ""),
            "sqlite_target_organism": metadata.get("target_organism", ""),
            "sqlite_target_tax_id": metadata.get("target_tax_id", ""),
            "sqlite_target_component_accessions": ";".join(component_accessions),
            "sqlite_target_single_protein": target_type_ok,
            "sqlite_target_human": human_ok,
            "sqlite_target_source_uniprot_exact": accession_ok,
            "sqlite_entity_pair_validated": entity_ok,
            "sqlite_entity_validation_status": status,
        }
        rows.append(row)
        mapping_lookup[mapping["candidate_mapping_id"]] = row
    return rows, mapping_lookup


def document_join_and_selects(schema: dict[str, set[str]]) -> tuple[str, list[str]]:
    """Return a safe optional document join and a stable set of select fields."""
    document_table = schema.get("docs")
    compound_record_table = schema.get("compound_records")
    join = ""
    doc_alias: str | None = None
    if document_table and compound_record_table and {"record_id", "doc_id"}.issubset(compound_record_table) and "record_id" in schema["activities"] and "doc_id" in document_table:
        join = (
            "LEFT JOIN compound_records AS cr ON cr.record_id = act.record_id\n"
            "LEFT JOIN docs AS doc ON doc.doc_id = cr.doc_id"
        )
        doc_alias = "doc"
    elif document_table and "doc_id" in schema["activities"] and "doc_id" in document_table:
        join = "LEFT JOIN docs AS doc ON doc.doc_id = act.doc_id"
        doc_alias = "doc"
    elif document_table and "doc_id" in schema["assays"] and "doc_id" in document_table:
        join = "LEFT JOIN docs AS doc ON doc.doc_id = assay.doc_id"
        doc_alias = "doc"
    fields: list[str] = []
    for column, output in (
        ("chembl_id", "document_chembl_id"),
        ("pubmed_id", "document_pubmed_id"),
        ("doi", "document_doi"),
        ("year", "document_year"),
        ("journal", "document_journal"),
        ("title", "document_title"),
    ):
        if doc_alias is not None:
            fields.append(select_optional(doc_alias, document_table or set(), column, output))
        else:
            fields.append(f"NULL AS {qident(output)}")
    return join, fields


def extract_activities(
    connection: sqlite3.Connection,
    validated_mappings: list[dict[str, Any]],
    schema: dict[str, set[str]],
) -> list[dict[str, Any]]:
    if not validated_mappings:
        return []
    document_join, document_selects = document_join_and_selects(schema)
    activity_selects = [
        select_optional("act", schema["activities"], "activity_id", "activity_id"),
        select_optional("act", schema["activities"], "record_id", "record_id"),
        select_optional("act", schema["activities"], "doc_id", "activity_doc_id"),
        select_optional("act", schema["activities"], "standard_type", "standard_type"),
        select_optional("act", schema["activities"], "standard_relation", "standard_relation"),
        select_optional("act", schema["activities"], "standard_value", "standard_value"),
        select_optional("act", schema["activities"], "standard_units", "standard_units"),
        select_optional("act", schema["activities"], "standard_flag", "standard_flag"),
        select_optional("act", schema["activities"], "pchembl_value", "pchembl_value"),
        select_optional("act", schema["activities"], "activity_comment", "activity_comment"),
        select_optional("act", schema["activities"], "data_validity_comment", "data_validity_comment"),
        select_optional("act", schema["activities"], "potential_duplicate", "potential_duplicate"),
        select_optional("act", schema["activities"], "bao_endpoint", "bao_endpoint"),
        select_optional("act", schema["activities"], "type", "reported_type"),
        select_optional("act", schema["activities"], "relation", "reported_relation"),
        select_optional("act", schema["activities"], "value", "reported_value"),
        select_optional("act", schema["activities"], "units", "reported_units"),
    ]
    assay_selects = [
        select_optional("assay", schema["assays"], "assay_id", "assay_id"),
        select_optional("assay", schema["assays"], "assay_type", "assay_type"),
        select_optional("assay", schema["assays"], "description", "assay_description"),
        select_optional("assay", schema["assays"], "assay_organism", "assay_organism"),
        select_optional("assay", schema["assays"], "assay_tax_id", "assay_tax_id"),
        select_optional("assay", schema["assays"], "assay_strain", "assay_strain"),
        select_optional("assay", schema["assays"], "assay_tissue", "assay_tissue"),
        select_optional("assay", schema["assays"], "assay_cell_type", "assay_cell_type"),
        select_optional("assay", schema["assays"], "confidence_score", "assay_confidence_score"),
        select_optional("assay", schema["assays"], "relationship_type", "assay_relationship_type"),
        select_optional("assay", schema["assays"], "bao_format", "assay_bao_format"),
        select_optional("assay", schema["assays"], "assay_category", "assay_category"),
        select_optional("assay", schema["assays"], "src_id", "assay_source_id"),
    ]
    target_selects = [
        "md.chembl_id AS chembl_compound_id",
        "md.molregno AS chembl_molregno",
        "td.chembl_id AS chembl_target_id",
        "td.tid AS chembl_tid",
        select_optional("td", schema["target_dictionary"], "pref_name", "activity_target_pref_name"),
        select_optional("td", schema["target_dictionary"], "target_type", "activity_target_type"),
        select_optional("td", schema["target_dictionary"], "organism", "activity_target_organism"),
        select_optional("td", schema["target_dictionary"], "tax_id", "activity_target_tax_id"),
    ]
    selects = ["wanted.candidate_mapping_id"] + activity_selects + assay_selects + target_selects + document_selects
    lookup = {row["candidate_mapping_id"]: row for row in validated_mappings}
    extracted: list[dict[str, Any]] = []
    seen: set[tuple[str, str]] = set()
    for batch in chunks(validated_mappings, 200):
        values_sql = ", ".join("(?, ?, ?)" for _ in batch)
        parameters: list[str] = []
        for item in batch:
            parameters.extend(
                [
                    item["candidate_mapping_id"],
                    item["chembl_compound_id"],
                    item["chembl_target_id"],
                ]
            )
        query = f"""
            WITH wanted(candidate_mapping_id, compound_chembl_id, target_chembl_id) AS (
                VALUES {values_sql}
            )
            SELECT {', '.join(selects)}
            FROM wanted
            JOIN molecule_dictionary AS md
              ON md.chembl_id = wanted.compound_chembl_id
            JOIN activities AS act ON act.molregno = md.molregno
            JOIN assays AS assay ON assay.assay_id = act.assay_id
            JOIN target_dictionary AS td
              ON td.tid = assay.tid
             AND td.chembl_id = wanted.target_chembl_id
            {document_join}
        """
        for result in connection.execute(query, parameters):
            row = dict(result)
            map_id = str(row["candidate_mapping_id"])
            activity_id = str(row.get("activity_id") or "")
            dedupe_key = (map_id, activity_id)
            if activity_id and dedupe_key in seen:
                continue
            if activity_id:
                seen.add(dedupe_key)
            source = lookup[map_id]
            assay_type = str(row.get("assay_type") or "").strip().upper()
            if assay_type == "B":
                context = "binding_assay_candidate_manual_review_required"
            elif assay_type == "F":
                context = "functional_assay_candidate_manual_review_required"
            else:
                context = "other_or_unclassified_assay_manual_review_required"
            extracted.append(
                {
                    "candidate_mapping_id": map_id,
                    "pair_key": source["pair_key"],
                    "inchikey_full": source["inchikey_full"],
                    "npass_uniprot_source": source["npass_uniprot_source"],
                    "best_evidence_tier_v1_1": source["best_evidence_tier_v1_1"],
                    "npass_record_count": source["npass_record_count"],
                    "sqlite_entity_pair_validated": True,
                    **row,
                    "historical_activity_context": context,
                    "direct_binding_asserted": False,
                    "manual_assay_and_source_review_required": True,
                }
            )
    return extracted


def distinct_nonempty(values: Iterable[Any]) -> str:
    return ";".join(sorted({str(value).strip() for value in values if str(value or "").strip()}))


def pair_audit_rows(
    candidates: list[dict[str, str]],
    entity_rows: list[dict[str, Any]],
    activity_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    mappings_by_pair: dict[str, list[dict[str, Any]]] = defaultdict(list)
    activities_by_pair: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for row in entity_rows:
        mappings_by_pair[row["pair_key"]].append(row)
    for row in activity_rows:
        activities_by_pair[row["pair_key"]].append(row)
    output: list[dict[str, Any]] = []
    for candidate in candidates:
        pair_key = candidate["pair_key"]
        mappings = mappings_by_pair[pair_key]
        activities = activities_by_pair[pair_key]
        validated = [item for item in mappings if item["sqlite_entity_pair_validated"]]
        activity_ids = {str(item.get("activity_id") or "") for item in activities}
        activity_ids.discard("")
        if activities:
            temporal_action = "exclude_from_future_candidate_pool_historical_chembl31_activity"
            hit_status = "historical_activity_recorded_in_chembl31"
        elif validated:
            temporal_action = "retain_pending_next_evidence_audit_no_chembl31_activity_found"
            hit_status = "no_activity_found_in_validated_chembl31_entity_pair"
        else:
            temporal_action = "exclude_until_entity_mapping_is_resolved"
            hit_status = "entity_pair_not_sqlite_validated"
        output.append(
            {
                "pair_key": pair_key,
                "inchikey_full": candidate["inchikey_full"],
                "npass_uniprot_source": candidate["npass_uniprot_source"],
                "best_evidence_tier_v1_1": candidate["best_evidence_tier_v1_1"],
                "npass_record_count": candidate["npass_record_count"],
                "preliminary_both_entities_exactly_mapped": true_value(candidate.get("both_entities_exactly_mapped")),
                "expanded_entity_mapping_count": len(mappings),
                "sqlite_validated_entity_mapping_count": len(validated),
                "historical_activity_row_count": len(activities),
                "historical_unique_activity_id_count": len(activity_ids),
                "historical_activity_hit_status": hit_status,
                "temporal_audit_action": temporal_action,
                "observed_assay_types": distinct_nonempty(item.get("assay_type") for item in activities),
                "observed_standard_types": distinct_nonempty(item.get("standard_type") for item in activities),
                "observed_assay_confidence_scores": distinct_nonempty(item.get("assay_confidence_score") for item in activities),
                "observed_document_pubmed_ids": distinct_nonempty(item.get("document_pubmed_id") for item in activities),
                "observed_document_years": distinct_nonempty(item.get("document_year") for item in activities),
                "binding_assay_candidate_row_count": sum(
                    str(item.get("assay_type") or "").strip().upper() == "B" for item in activities
                ),
                "functional_assay_candidate_row_count": sum(
                    str(item.get("assay_type") or "").strip().upper() == "F" for item in activities
                ),
                "direct_binding_asserted": False,
                "interpretation": (
                    "A ChEMBL31 activity hit is historical-record evidence for leakage control, "
                    "not a direct-binding claim."
                ),
            }
        )
    return output


def serialise(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    return value


def write_csv_gz(path: Path, rows: list[dict[str, Any]], preferred: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(preferred)
    observed = {key for row in rows for key in row}
    fieldnames.extend(sorted(observed.difference(fieldnames)))
    with deterministic_gzip_text(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: serialise(row.get(key, "")) for key in fieldnames})


def schema_summary(connection: sqlite3.Connection, schema: dict[str, set[str]]) -> dict[str, Any]:
    return {
        "sqlite_version": connection.execute("SELECT sqlite_version()").fetchone()[0],
        "user_version": connection.execute("PRAGMA user_version").fetchone()[0],
        "schema_tables_used": {
            table: sorted(schema[table])
            for table in sorted(
                {
                    "activities",
                    "assays",
                    "target_dictionary",
                    "molecule_dictionary",
                    "compound_structures",
                    "target_components",
                    "component_sequences",
                    "docs",
                    "compound_records",
                }.intersection(schema)
            )
        },
    }


def run_audit(root: Path, database: Path, mapping_path: Path) -> dict[str, Any]:
    candidates, mapping_columns = read_candidates(mapping_path)
    preliminary = expand_preliminary_mappings(candidates)
    with connect_read_only(database) as connection:
        schema = ensure_schema(connection)
        compound_ids = sorted({item["chembl_compound_id"] for item in preliminary})
        target_ids = sorted({item["chembl_target_id"] for item in preliminary})
        compound_keys = query_compound_keys(connection, compound_ids)
        target_metadata = query_target_metadata(connection, target_ids, schema)
        entity_rows, mapping_lookup = validate_entities(
            preliminary, compound_keys, target_metadata
        )
        validated = [
            mapping_lookup[item["candidate_mapping_id"]]
            for item in preliminary
            if mapping_lookup[item["candidate_mapping_id"]]["sqlite_entity_pair_validated"]
        ]
        activities = extract_activities(connection, validated, schema)
        schema_info = schema_summary(connection, schema)
    pair_rows = pair_audit_rows(candidates, entity_rows, activities)
    interim = root / "data/interim"
    results = root / "results"
    entity_path = interim / "chembl_31_future_candidate_sqlite_entity_validation.csv.gz"
    activity_path = interim / "chembl_31_future_candidate_historical_activity_rows.csv.gz"
    pair_path = interim / "chembl_31_future_candidate_historical_pair_audit.csv.gz"
    summary_path = results / "chembl_31_historical_activity_audit_summary.json"
    write_csv_gz(
        entity_path,
        entity_rows,
        [
            "candidate_mapping_id", "pair_key", "inchikey_full", "npass_uniprot_source",
            "best_evidence_tier_v1_1", "chembl_compound_id", "chembl_target_id",
            "sqlite_compound_full_inchikey_exact", "sqlite_target_single_protein",
            "sqlite_target_human", "sqlite_target_source_uniprot_exact",
            "sqlite_entity_pair_validated", "sqlite_entity_validation_status",
        ],
    )
    write_csv_gz(
        activity_path,
        activities,
        [
            "candidate_mapping_id", "pair_key", "inchikey_full", "npass_uniprot_source",
            "best_evidence_tier_v1_1", "chembl_compound_id", "chembl_target_id",
            "activity_id", "assay_id", "assay_type", "assay_confidence_score",
            "standard_type", "standard_relation", "standard_value", "standard_units",
            "pchembl_value", "document_pubmed_id", "document_doi", "document_year",
            "historical_activity_context", "direct_binding_asserted",
        ],
    )
    write_csv_gz(
        pair_path,
        pair_rows,
        [
            "pair_key", "inchikey_full", "npass_uniprot_source",
            "best_evidence_tier_v1_1", "historical_activity_hit_status",
            "temporal_audit_action", "historical_activity_row_count",
            "historical_unique_activity_id_count", "observed_assay_types",
            "observed_standard_types", "direct_binding_asserted",
        ],
    )
    hit_pairs = [
        item
        for item in pair_rows
        if item["historical_activity_hit_status"] == "historical_activity_recorded_in_chembl31"
    ]
    summary: dict[str, Any] = {
        "audit_name": "historical ChEMBL 31 exact-entity activity audit",
        "run_at_utc": utc_now(),
        "database_path": str(database),
        "database_bytes": database.stat().st_size,
        "database_open_mode": "SQLite URI mode=ro with PRAGMA query_only=ON",
        "candidate_mapping_path": str(mapping_path),
        "candidate_mapping_columns": mapping_columns,
        "candidate_pair_count": len(candidates),
        "preliminary_exact_entity_mapped_pair_count": sum(
            true_value(item.get("both_entities_exactly_mapped")) for item in candidates
        ),
        "expanded_entity_mapping_count": len(preliminary),
        "sqlite_validated_entity_mapping_count": sum(
            item["sqlite_entity_pair_validated"] for item in entity_rows
        ),
        "candidate_pairs_with_sqlite_validated_entity_mapping": sum(
            item["sqlite_validated_entity_mapping_count"] > 0 for item in pair_rows
        ),
        "candidate_pairs_with_historical_activity_hits": len(hit_pairs),
        "historical_activity_row_count": len(activities),
        "output_files": {
            "entity_validation": str(entity_path),
            "historical_activity_rows": str(activity_path),
            "pair_audit": str(pair_path),
        },
        "interpretation": [
            "A hit means ChEMBL 31 contains at least one activity row for a full-InChIKey and exact human SINGLE PROTEIN UniProt-validated mapping.",
            "For conservative temporal-leakage control, hit pairs are excluded from the future-candidate pool.",
            "No hit means no such row was found in ChEMBL 31; it is not evidence of biological inactivity or absence from all pre-cutoff literature.",
            "Assay type B and quantitative fields are metadata for manual evidence review; this audit never asserts direct binding.",
        ],
        **schema_info,
    }
    results.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, default=str) + "\n",
        encoding="utf-8",
    )
    return {"summary_path": str(summary_path), **summary}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
        help="Project root containing data/, results/, and scripts/.",
    )
    parser.add_argument(
        "--sqlite-db",
        type=Path,
        help="Explicit, already-unpacked ChEMBL 31 SQLite database. The script does not unpack archives.",
    )
    parser.add_argument(
        "--candidate-mapping",
        type=Path,
        help=f"Precomputed exact entity map (default: {DEFAULT_MAPPING}).",
    )
    parser.add_argument(
        "--wait-seconds",
        type=int,
        default=0,
        help="Optionally poll for an unpacked SQLite database; never downloads or extracts anything.",
    )
    parser.add_argument(
        "--readiness-only",
        action="store_true",
        help="Print database readiness JSON and exit without querying SQLite or writing outputs.",
    )
    args = parser.parse_args()
    if args.wait_seconds < 0:
        parser.error("--wait-seconds must be non-negative")
    root = args.project_root.expanduser().resolve()
    mapping_path = (
        args.candidate_mapping.expanduser().resolve()
        if args.candidate_mapping
        else root / DEFAULT_MAPPING
    )
    chembl_root = root / "data/raw/chembl/chembl_31"
    database, readiness = wait_for_database(chembl_root, args.sqlite_db, args.wait_seconds)
    if args.readiness_only:
        print(json.dumps(readiness, ensure_ascii=False, indent=2))
        return 0
    if database is None:
        print(json.dumps(readiness, ensure_ascii=False, indent=2), file=sys.stderr)
        print(
            "No unpacked SQLite database is ready. This script does not download or extract the "
            "archive; verify the archive checksum and unpack it, then rerun with --sqlite-db.",
            file=sys.stderr,
        )
        return 2
    result = run_audit(root, database, mapping_path)
    print(result["summary_path"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
