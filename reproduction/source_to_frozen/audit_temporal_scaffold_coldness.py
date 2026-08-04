#!/usr/bin/env python3
"""Audit scaffold coldness of observed strict temporal future pairs.

This program is deliberately an *audit*, not a split generator or a negative
sampler.  It reads date-verified, pair-level historical and future tables plus
an RDKit-derived Bemis-Murcko scaffold table.  It writes a ledger for every
observed future pair and never enumerates unrecorded compound-target pairs.

Only rows with ``scaffold_status == 'ok'`` and a non-empty ``scaffold_key``
are allowed to form scaffold groups.  Invalid structures and acyclic/empty
Bemis-Murcko results are recorded explicitly but never share an empty group.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import shutil
import sys
import uuid
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, TextIO

from reproducible_io import deterministic_gzip_text


SCHEMA_VERSION = "1.0"
HISTORICAL_DECISION = "strict_pre_cutoff_training_candidate"
FUTURE_DECISION = "strict_post_cutoff_future_candidate"


class InputContractError(ValueError):
    """Raised when an input cannot support a reproducible coldness audit."""


@dataclass(frozen=True)
class PairRow:
    pair_key: str
    inchikey_full: str
    uniprot_accession: str
    source_row_count: int
    values: dict[str, str]


@dataclass(frozen=True)
class ScaffoldRow:
    inchikey_full: str
    scaffold_status: str
    scaffold_key: str
    bemis_murcko_smiles: str
    representative_smiles: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def open_text(path: Path, mode: str) -> TextIO:
    if path.suffix.lower() == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def detect_delimiter(path: Path) -> str:
    with open_text(path, "rt") as handle:
        for line in handle:
            if line.strip():
                # These project tables are CSV or TSV.  Counting is more
                # predictable than csv.Sniffer for fields containing semicolon
                # lists and long rationale strings.
                return "\t" if line.count("\t") > line.count(",") else ","
    raise InputContractError(f"Empty input table: {path}")


def normalized(value: object) -> str:
    return str(value or "").strip()


def normalized_inchikey(value: object) -> str:
    return normalized(value).upper()


def normalized_accession(value: object) -> str:
    return normalized(value).upper()


def first_present(fieldnames: Iterable[str], aliases: tuple[str, ...], label: str, path: Path) -> str:
    available = set(fieldnames)
    for alias in aliases:
        if alias in available:
            return alias
    raise InputContractError(
        f"{path} lacks the required {label} column. Expected one of {list(aliases)}; "
        f"found {sorted(available)}"
    )


def assert_regular_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise FileNotFoundError(f"{label} does not exist or is not a file: {path}")


def pair_key_from_values(inchikey_full: str, uniprot_accession: str) -> str:
    return f"{inchikey_full}|{uniprot_accession}"


def load_pair_table(
    path: Path,
    *,
    expected_decision: str | None,
    require_strict_decision: bool,
) -> tuple[list[PairRow], list[str], int]:
    """Read a pair-level strict temporal table and collapse exact duplicates.

    Every output row remains an observed input pair.  Duplicate source rows for
    the same chemical-protein pair are collapsed only after confirming the
    identity is consistent; their multiplicity is retained in
    ``source_row_count``.
    """

    assert_regular_file(path, "Pair table")
    delimiter = detect_delimiter(path)
    entries: dict[str, PairRow] = {}
    source_rows = 0
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise InputContractError(f"Missing header in pair table: {path}")
        fieldnames = list(reader.fieldnames)
        inchikey_column = first_present(
            fieldnames,
            ("inchikey_full", "full_inchikey"),
            "full InChIKey",
            path,
        )
        accession_column = first_present(
            fieldnames,
            ("uniprot_canonical_accession", "uniprot_accession", "target_accession"),
            "canonical UniProt accession",
            path,
        )
        supplied_pair_column = next(
            (column for column in ("canonical_pair_key", "pair_key") if column in fieldnames),
            None,
        )
        if require_strict_decision and "decision" not in fieldnames:
            raise InputContractError(
                f"{path} has no 'decision' column, so it cannot be verified as a strict temporal table. "
                "Use --allow-unverified-decision only for an explicitly documented non-primary audit."
            )

        for line_number, raw in enumerate(reader, start=2):
            source_rows += 1
            inchikey_full = normalized_inchikey(raw.get(inchikey_column))
            accession = normalized_accession(raw.get(accession_column))
            if not inchikey_full or not accession:
                raise InputContractError(
                    f"{path}:{line_number} has a blank full InChIKey or canonical UniProt accession"
                )
            if "|" in inchikey_full or "|" in accession:
                raise InputContractError(
                    f"{path}:{line_number} contains '|' inside an entity identifier, which is unsupported"
                )
            calculated_pair_key = pair_key_from_values(inchikey_full, accession)
            if supplied_pair_column:
                supplied_pair_key = normalized(raw.get(supplied_pair_column))
                if supplied_pair_key and supplied_pair_key.upper() != calculated_pair_key:
                    raise InputContractError(
                        f"{path}:{line_number} has inconsistent {supplied_pair_column}={supplied_pair_key!r}; "
                        f"expected {calculated_pair_key!r} from its entity columns"
                    )
            if require_strict_decision:
                observed_decision = normalized(raw.get("decision"))
                if observed_decision != expected_decision:
                    raise InputContractError(
                        f"{path}:{line_number} decision={observed_decision!r}; expected {expected_decision!r}. "
                        "Do not mix a non-strict candidate pool into this audit."
                    )
            cleaned = {column: normalized(raw.get(column)) for column in fieldnames}
            # Preserve canonical, normalized identity in any common source field.
            cleaned[inchikey_column] = inchikey_full
            cleaned[accession_column] = accession
            if supplied_pair_column:
                cleaned[supplied_pair_column] = calculated_pair_key

            existing = entries.get(calculated_pair_key)
            if existing is None:
                entries[calculated_pair_key] = PairRow(
                    pair_key=calculated_pair_key,
                    inchikey_full=inchikey_full,
                    uniprot_accession=accession,
                    source_row_count=1,
                    values=cleaned,
                )
            else:
                # The pair identity has already been checked.  Non-identity
                # metadata may legitimately be repeated by a pair-level export;
                # retain the first row and report the multiplicity rather than
                # inventing an aggregate label.
                entries[calculated_pair_key] = PairRow(
                    pair_key=existing.pair_key,
                    inchikey_full=existing.inchikey_full,
                    uniprot_accession=existing.uniprot_accession,
                    source_row_count=existing.source_row_count + 1,
                    values=existing.values,
                )
    if not entries:
        raise InputContractError(f"No pair rows were read from {path}")
    return [entries[key] for key in sorted(entries)], fieldnames, source_rows


def load_scaffold_table(path: Path) -> tuple[dict[str, ScaffoldRow], int, list[str]]:
    """Load one scaffold assignment per full InChIKey without empty grouping."""

    assert_regular_file(path, "Scaffold table")
    delimiter = detect_delimiter(path)
    rows: dict[str, ScaffoldRow] = {}
    source_rows = 0
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise InputContractError(f"Missing header in scaffold table: {path}")
        fieldnames = list(reader.fieldnames)
        inchikey_column = first_present(
            fieldnames,
            ("inchikey_full", "full_inchikey", "molecule_id"),
            "full InChIKey/molecule identifier",
            path,
        )
        status_column = first_present(fieldnames, ("scaffold_status",), "scaffold status", path)
        key_column = first_present(
            fieldnames,
            ("scaffold_key", "bemis_murcko_smiles"),
            "non-empty scaffold key field",
            path,
        )
        murcko_column = "bemis_murcko_smiles" if "bemis_murcko_smiles" in fieldnames else key_column
        smiles_column = "representative_smiles" if "representative_smiles" in fieldnames else ""

        for line_number, raw in enumerate(reader, start=2):
            source_rows += 1
            inchikey_full = normalized_inchikey(raw.get(inchikey_column))
            if not inchikey_full:
                raise InputContractError(f"{path}:{line_number} has a blank compound identifier")
            scaffold = ScaffoldRow(
                inchikey_full=inchikey_full,
                scaffold_status=normalized(raw.get(status_column)),
                scaffold_key=normalized(raw.get(key_column)),
                bemis_murcko_smiles=normalized(raw.get(murcko_column)),
                representative_smiles=normalized(raw.get(smiles_column)) if smiles_column else "",
            )
            existing = rows.get(inchikey_full)
            if existing is not None and existing != scaffold:
                raise InputContractError(
                    f"{path}:{line_number} conflicts with an earlier scaffold assignment for {inchikey_full}. "
                    "Resolve the compound identity/scaffold conflict before auditing coldness."
                )
            rows[inchikey_full] = scaffold
    if not rows:
        raise InputContractError(f"No scaffold rows were read from {path}")
    return rows, source_rows, fieldnames


def scaffold_mapping_category(scaffold: ScaffoldRow | None) -> str:
    """Classify a molecule without ever treating an empty key as a group."""

    if scaffold is None:
        return "missing_scaffold_row"
    status = scaffold.scaffold_status.casefold()
    if status == "ok" and scaffold.scaffold_key:
        return "eligible_nonempty_scaffold"
    if status == "invalid_smiles":
        return "invalid_smiles"
    if status == "acyclic_or_empty_bemis_murcko":
        return "acyclic_or_empty_bemis_murcko"
    if status == "ok" and not scaffold.scaffold_key:
        return "ok_status_but_empty_scaffold_key"
    if not status:
        return "missing_scaffold_status"
    return "unsupported_scaffold_status"


def write_tsv_gz(path: Path, fieldnames: list[str], rows: Iterable[dict[str, object]]) -> None:
    with deterministic_gzip_text(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def safe_rmtree(path: Path, project_root: Path) -> None:
    """Remove only a private temporary output below the declared project root."""

    resolved_path = path.resolve()
    resolved_root = project_root.resolve()
    try:
        resolved_path.relative_to(resolved_root)
    except ValueError as exc:
        raise RuntimeError(f"Refusing to remove path outside project root: {resolved_path}") from exc
    if resolved_path.exists():
        shutil.rmtree(resolved_path)


def compound_pair_counts(pairs: Iterable[PairRow]) -> dict[str, int]:
    counts: Counter[str] = Counter()
    for pair in pairs:
        counts[pair.inchikey_full] += 1
    return dict(counts)


def audit(
    *,
    project_root: Path,
    historical_pairs_path: Path,
    future_pairs_path: Path,
    scaffold_table_path: Path,
    output_dir: Path,
    coverage_policy: str,
    require_strict_decision: bool,
) -> dict[str, object]:
    historical_pairs, historical_fields, historical_source_rows = load_pair_table(
        historical_pairs_path,
        expected_decision=HISTORICAL_DECISION,
        require_strict_decision=require_strict_decision,
    )
    future_pairs, future_fields, future_source_rows = load_pair_table(
        future_pairs_path,
        expected_decision=FUTURE_DECISION,
        require_strict_decision=require_strict_decision,
    )
    scaffold_rows, scaffold_source_rows, scaffold_fields = load_scaffold_table(scaffold_table_path)

    historical_pair_keys = {pair.pair_key for pair in historical_pairs}
    historical_compounds = {pair.inchikey_full for pair in historical_pairs}
    future_compounds = {pair.inchikey_full for pair in future_pairs}
    historical_pair_counts = compound_pair_counts(historical_pairs)
    future_pair_counts = compound_pair_counts(future_pairs)

    historical_categories = {
        compound: scaffold_mapping_category(scaffold_rows.get(compound))
        for compound in historical_compounds
    }
    future_categories = {
        compound: scaffold_mapping_category(scaffold_rows.get(compound))
        for compound in future_compounds
    }
    historical_mapping_counts = Counter(historical_categories.values())
    future_mapping_counts = Counter(future_categories.values())

    # A non-empty, valid scaffold is the *only* entity allowed into a group.
    # In particular, no call to setdefault() occurs for an empty scaffold key.
    historical_group_compounds: dict[str, set[str]] = defaultdict(set)
    historical_group_pairs: Counter[str] = Counter()
    for pair in historical_pairs:
        scaffold = scaffold_rows.get(pair.inchikey_full)
        if scaffold_mapping_category(scaffold) == "eligible_nonempty_scaffold":
            assert scaffold is not None and scaffold.scaffold_key
            historical_group_compounds[scaffold.scaffold_key].add(pair.inchikey_full)
            historical_group_pairs[scaffold.scaffold_key] += 1
    historical_group_keys = set(historical_group_compounds)

    # Acyclic/empty Bemis-Murcko results are known non-groupable molecules; they
    # cannot match a non-empty future scaffold.  Missing/invalid/unsupported
    # rows are unresolved and therefore block an absolute new-scaffold claim in
    # the default strict-coverage mode.
    historical_coverage_blocker_categories = {
        "missing_scaffold_row",
        "invalid_smiles",
        "ok_status_but_empty_scaffold_key",
        "missing_scaffold_status",
        "unsupported_scaffold_status",
    }
    historical_coverage_blocker_compounds = sorted(
        compound
        for compound, category in historical_categories.items()
        if category in historical_coverage_blocker_categories
    )
    historical_coverage_complete = not historical_coverage_blocker_compounds

    # All future pair rows are preserved.  The first outcome is deliberately
    # checked before compound coldness so a pair that leaked into history is
    # never presented as a valid future cold-start evaluation item.
    audit_rows: list[dict[str, object]] = []
    outcome_counts: Counter[str] = Counter()
    cold_compounds_under_policy: set[str] = set()
    cold_pairs_under_policy: set[str] = set()
    candidate_cold_pairs_against_mapped_history: set[str] = set()
    future_group_compounds: dict[str, set[str]] = defaultdict(set)
    future_group_pairs: Counter[str] = Counter()

    audit_fields = [
        "audit_pair_key",
        "audit_inchikey_full",
        "audit_uniprot_canonical_accession",
        "audit_source_row_count",
        "audit_future_pair_also_in_historical_table",
        "audit_future_compound_seen_in_historical_table",
        "audit_scaffold_mapping_category",
        "audit_scaffold_status",
        "audit_scaffold_key",
        "audit_historical_same_scaffold_compound_count",
        "audit_historical_same_scaffold_pair_count",
        "audit_historical_scaffold_coverage_complete",
        "audit_historical_scaffold_coverage_blocker_compound_count",
        "audit_candidate_cold_against_mapped_eligible_history",
        "audit_scaffold_cold_under_selected_policy",
        "audit_coldness_scope",
        "audit_outcome",
        "audit_eligibility_or_exclusion_reason",
        "audit_unrecorded_pair_policy",
    ]
    for pair in future_pairs:
        scaffold = scaffold_rows.get(pair.inchikey_full)
        category = future_categories[pair.inchikey_full]
        scaffold_key = scaffold.scaffold_key if scaffold is not None else ""
        same_scaffold_compound_count = (
            len(historical_group_compounds.get(scaffold_key, set()))
            if category == "eligible_nonempty_scaffold"
            else 0
        )
        same_scaffold_pair_count = (
            historical_group_pairs.get(scaffold_key, 0)
            if category == "eligible_nonempty_scaffold"
            else 0
        )
        pair_also_historical = pair.pair_key in historical_pair_keys
        compound_seen_historically = pair.inchikey_full in historical_compounds
        candidate_cold = category == "eligible_nonempty_scaffold" and scaffold_key not in historical_group_keys
        if candidate_cold:
            candidate_cold_pairs_against_mapped_history.add(pair.pair_key)

        cold_under_policy = False
        coldness_scope = "not_applicable"
        if pair_also_historical:
            outcome = "excluded_future_pair_also_observed_in_historical_table"
            reason = "The same observed compound-target pair is present in the historical strict table."
        elif category != "eligible_nonempty_scaffold":
            outcome = f"excluded_future_{category}"
            reason = (
                "This future compound cannot enter a Bemis-Murcko scaffold group under the recorded scaffold "
                "assignment; it is not converted into an empty or shared scaffold group."
            )
        elif scaffold_key in historical_group_keys:
            outcome = "not_scaffold_cold__scaffold_seen_in_eligible_historical_compounds"
            reason = "A non-empty, valid scaffold key matching this future compound was observed among historical compounds."
            coldness_scope = "eligible_historical_nonempty_scaffold_groups"
        elif coverage_policy == "strict" and not historical_coverage_complete:
            outcome = "coldness_unresolved__historical_scaffold_coverage_incomplete"
            reason = (
                "The scaffold key is absent from mapped historical non-empty scaffold groups, but at least one "
                "historical compound has a missing, invalid, empty-status, or unsupported scaffold assignment. "
                "Under strict coverage, this pair is not eligible for a primary new-scaffold claim."
            )
            coldness_scope = "mapped_eligible_history_only__not_primary_due_to_incomplete_historical_coverage"
        else:
            outcome = "eligible_scaffold_cold__scaffold_absent_from_historical_nonempty_groups"
            reason = "The non-empty, valid scaffold key is absent from all eligible historical compound scaffold groups."
            cold_under_policy = True
            coldness_scope = (
                "all_historical_compounds_have_resolved_or_known_non_groupable_scaffold_assignments"
                if historical_coverage_complete
                else "mapped_eligible_history_only__non_strict_policy"
            )
            cold_compounds_under_policy.add(pair.inchikey_full)
            cold_pairs_under_policy.add(pair.pair_key)

        outcome_counts[outcome] += 1
        row = dict(pair.values)
        row.update(
            {
                "audit_pair_key": pair.pair_key,
                "audit_inchikey_full": pair.inchikey_full,
                "audit_uniprot_canonical_accession": pair.uniprot_accession,
                "audit_source_row_count": pair.source_row_count,
                "audit_future_pair_also_in_historical_table": pair_also_historical,
                "audit_future_compound_seen_in_historical_table": compound_seen_historically,
                "audit_scaffold_mapping_category": category,
                "audit_scaffold_status": scaffold.scaffold_status if scaffold is not None else "",
                "audit_scaffold_key": scaffold_key,
                "audit_historical_same_scaffold_compound_count": same_scaffold_compound_count,
                "audit_historical_same_scaffold_pair_count": same_scaffold_pair_count,
                "audit_historical_scaffold_coverage_complete": historical_coverage_complete,
                "audit_historical_scaffold_coverage_blocker_compound_count": len(historical_coverage_blocker_compounds),
                "audit_candidate_cold_against_mapped_eligible_history": candidate_cold,
                "audit_scaffold_cold_under_selected_policy": cold_under_policy,
                "audit_coldness_scope": coldness_scope,
                "audit_outcome": outcome,
                "audit_eligibility_or_exclusion_reason": reason,
                "audit_unrecorded_pair_policy": "Only observed strict future pairs are emitted; unrecorded pairs are neither emitted nor labeled negative.",
            }
        )
        audit_rows.append(row)
        if category == "eligible_nonempty_scaffold":
            future_group_compounds[scaffold_key].add(pair.inchikey_full)
            future_group_pairs[scaffold_key] += 1

    compound_rows: list[dict[str, object]] = []
    for compound in sorted(historical_compounds | future_compounds):
        scaffold = scaffold_rows.get(compound)
        category = scaffold_mapping_category(scaffold)
        key = scaffold.scaffold_key if scaffold is not None else ""
        compound_rows.append(
            {
                "inchikey_full": compound,
                "in_historical_strict_pair_table": compound in historical_compounds,
                "in_future_strict_pair_table": compound in future_compounds,
                "historical_unique_pair_count": historical_pair_counts.get(compound, 0),
                "future_unique_pair_count": future_pair_counts.get(compound, 0),
                "scaffold_mapping_category": category,
                "scaffold_status": scaffold.scaffold_status if scaffold is not None else "",
                "scaffold_key": key,
                "bemis_murcko_smiles": scaffold.bemis_murcko_smiles if scaffold is not None else "",
                "representative_smiles": scaffold.representative_smiles if scaffold is not None else "",
                "is_eligible_nonempty_scaffold_group_member": category == "eligible_nonempty_scaffold",
                "historical_same_scaffold_compound_count": len(historical_group_compounds.get(key, set())) if key else 0,
                "historical_same_scaffold_pair_count": historical_group_pairs.get(key, 0) if key else 0,
                "historical_scaffold_coverage_complete": historical_coverage_complete,
                "compound_cold_against_mapped_eligible_history": (
                    compound in future_compounds
                    and category == "eligible_nonempty_scaffold"
                    and key not in historical_group_keys
                ),
                "compound_cold_under_selected_policy": compound in cold_compounds_under_policy,
                "compound_note": (
                    "Invalid/acyclic/empty scaffolds are excluded from all scaffold groups; they are never assigned a shared blank group."
                    if category != "eligible_nonempty_scaffold"
                    else "Non-empty scaffold membership is evaluated only against eligible historical compounds."
                ),
            }
        )

    group_rows: list[dict[str, object]] = []
    for scaffold_key in sorted(historical_group_keys | set(future_group_compounds)):
        historical_compound_count = len(historical_group_compounds.get(scaffold_key, set()))
        future_compound_count = len(future_group_compounds.get(scaffold_key, set()))
        group_rows.append(
            {
                "scaffold_key": scaffold_key,
                "historical_eligible_compound_count": historical_compound_count,
                "historical_eligible_pair_count": historical_group_pairs.get(scaffold_key, 0),
                "future_eligible_compound_count": future_compound_count,
                "future_eligible_pair_count": future_group_pairs.get(scaffold_key, 0),
                "seen_in_historical_eligible_compounds": historical_compound_count > 0,
                "future_scaffold_cold_under_selected_policy": (
                    historical_compound_count == 0
                    and any(
                        row["audit_scaffold_key"] == scaffold_key
                        and row["audit_scaffold_cold_under_selected_policy"]
                        for row in audit_rows
                    )
                ),
            }
        )

    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite existing output directory: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    temp_dir = output_dir.parent / f".{output_dir.name}.tmp-{uuid.uuid4().hex}"
    temp_dir.mkdir()
    try:
        pair_output = temp_dir / "future_pair_scaffold_coldness_audit.tsv.gz"
        compound_output = temp_dir / "compound_scaffold_coldness_audit.tsv.gz"
        group_output = temp_dir / "nonempty_scaffold_group_summary.tsv.gz"
        pair_output_fields = list(future_fields)
        for field in audit_fields:
            if field not in pair_output_fields:
                pair_output_fields.append(field)
        write_tsv_gz(pair_output, pair_output_fields, audit_rows)
        write_tsv_gz(
            compound_output,
            [
                "inchikey_full",
                "in_historical_strict_pair_table",
                "in_future_strict_pair_table",
                "historical_unique_pair_count",
                "future_unique_pair_count",
                "scaffold_mapping_category",
                "scaffold_status",
                "scaffold_key",
                "bemis_murcko_smiles",
                "representative_smiles",
                "is_eligible_nonempty_scaffold_group_member",
                "historical_same_scaffold_compound_count",
                "historical_same_scaffold_pair_count",
                "historical_scaffold_coverage_complete",
                "compound_cold_against_mapped_eligible_history",
                "compound_cold_under_selected_policy",
                "compound_note",
            ],
            compound_rows,
        )
        write_tsv_gz(
            group_output,
            [
                "scaffold_key",
                "historical_eligible_compound_count",
                "historical_eligible_pair_count",
                "future_eligible_compound_count",
                "future_eligible_pair_count",
                "seen_in_historical_eligible_compounds",
                "future_scaffold_cold_under_selected_policy",
            ],
            group_rows,
        )
        output_hashes = {path.name: sha256(path) for path in sorted(temp_dir.iterdir()) if path.is_file()}
        summary = {
            "schema_version": SCHEMA_VERSION,
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "purpose": "Scaffold-coldness audit of observed strict temporal future pairs only; not a negative-sampling or model-training table.",
            "inputs": {
                "historical_strict_pairs": {"path": str(historical_pairs_path), "sha256": sha256(historical_pairs_path)},
                "future_strict_pairs": {"path": str(future_pairs_path), "sha256": sha256(future_pairs_path)},
                "bemis_murcko_scaffolds": {"path": str(scaffold_table_path), "sha256": sha256(scaffold_table_path)},
            },
            "input_columns": {
                "historical_pair_table": historical_fields,
                "future_pair_table": future_fields,
                "scaffold_table": scaffold_fields,
            },
            "policies": {
                "strict_decision_validation": require_strict_decision,
                "historical_expected_decision": HISTORICAL_DECISION if require_strict_decision else None,
                "future_expected_decision": FUTURE_DECISION if require_strict_decision else None,
                "historical_coverage_policy": coverage_policy,
                "group_membership_rule": "Only scaffold_status='ok' with non-empty scaffold_key is groupable.",
                "empty_scaffold_rule": "invalid_smiles and acyclic_or_empty_bemis_murcko are excluded from groups and never share a blank scaffold group.",
                "unrecorded_pair_policy": "No unrecorded compound-target pairs are created, scored, or labeled as negative by this audit.",
            },
            "counts": {
                "historical": {
                    "source_rows": historical_source_rows,
                    "unique_observed_pairs": len(historical_pairs),
                    "unique_compounds": len(historical_compounds),
                    "scaffold_mapping_categories": dict(sorted(historical_mapping_counts.items())),
                    "eligible_nonempty_scaffold_groups": len(historical_group_keys),
                },
                "future": {
                    "source_rows": future_source_rows,
                    "unique_observed_pairs": len(future_pairs),
                    "unique_compounds": len(future_compounds),
                    "scaffold_mapping_categories": dict(sorted(future_mapping_counts.items())),
                    "audit_outcomes_by_pair": dict(sorted(outcome_counts.items())),
                    "candidate_cold_pairs_against_mapped_eligible_history": len(candidate_cold_pairs_against_mapped_history),
                    "scaffold_cold_pairs_under_selected_policy": len(cold_pairs_under_policy),
                    "scaffold_cold_compounds_under_selected_policy": len(cold_compounds_under_policy),
                },
                "cross_table": {
                    "future_pairs_also_observed_historically": len({pair.pair_key for pair in future_pairs} & historical_pair_keys),
                    "future_compounds_also_observed_historically": len(future_compounds & historical_compounds),
                    "scaffold_rows_source": scaffold_source_rows,
                    "scaffold_rows_not_referenced_by_strict_pair_tables": len(set(scaffold_rows) - historical_compounds - future_compounds),
                },
                "historical_scaffold_coverage": {
                    "complete_for_strict_new_scaffold_claim": historical_coverage_complete,
                    "blocker_compound_count": len(historical_coverage_blocker_compounds),
                    "blocker_compounds": historical_coverage_blocker_compounds,
                    "known_non_groupable_acyclic_or_empty_compound_count": historical_mapping_counts.get("acyclic_or_empty_bemis_murcko", 0),
                },
            },
            "outputs_sha256": output_hashes,
            "limitations": [
                "Scaffold coldness is chemical-group isolation only; it does not establish protein homology coldness or double-cold status.",
                "Input labels remain whatever evidence status they had in the strict temporal tables; this audit does not promote P1 candidates to P2/P3.",
                "A non-empty Bemis-Murcko scaffold is a structural grouping rule, not a claim of pharmacological novelty.",
            ],
        }
        summary_path = temp_dir / "summary.json"
        summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        # Rename only after every output has been written.  The destination did
        # not exist at the start, so this cannot overwrite a prior audit.
        temp_dir.rename(output_dir)
    except Exception:
        safe_rmtree(temp_dir, project_root)
        raise

    return json.loads((output_dir / "summary.json").read_text(encoding="utf-8"))


def default_path(project_root: Path, relative: str) -> Path:
    return project_root / Path(relative)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument(
        "--historical-pairs",
        type=Path,
        default=None,
        help="Date-verified strict pre-cutoff pair table (default: v1_1 PMID-verified table).",
    )
    parser.add_argument(
        "--future-pairs",
        type=Path,
        default=None,
        help="Date-verified strict post-cutoff pair table (default: v1_1 PMID-verified table).",
    )
    parser.add_argument(
        "--scaffold-table",
        required=True,
        type=Path,
        help="RDKit-derived table with inchikey_full/molecule_id, scaffold_status, and scaffold_key.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="New output directory; the script refuses to overwrite an existing directory.",
    )
    parser.add_argument(
        "--historical-coverage-policy",
        choices=("strict", "mapped-history-only"),
        default="strict",
        help=(
            "'strict' (default) leaves new-scaffold candidates unresolved if any historical compound has an "
            "unresolved scaffold assignment. 'mapped-history-only' reports the narrower, non-primary comparison."
        ),
    )
    parser.add_argument(
        "--allow-unverified-decision",
        action="store_true",
        help="Disable strict decision-column validation. Only use for a separately documented exploratory audit.",
    )
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    historical_pairs = (args.historical_pairs or default_path(
        project_root, "data/processed/strict_temporal_training_candidates_v1_1_pmid_verified.csv.gz"
    )).resolve()
    future_pairs = (args.future_pairs or default_path(
        project_root, "data/processed/strict_temporal_future_candidates_v1_1_pmid_verified.csv.gz"
    )).resolve()
    scaffold_table = args.scaffold_table.resolve()
    output_dir = (args.output_dir or default_path(
        project_root, "data/interim/scaffold_coldness_audit_v1"
    )).resolve()

    # A caller may choose another data volume, but output must remain below the
    # declared project root so temporary cleanup cannot affect unrelated paths.
    try:
        output_dir.relative_to(project_root)
    except ValueError as exc:
        raise SystemExit(f"--output-dir must be inside --project-root: {output_dir}") from exc

    summary = audit(
        project_root=project_root,
        historical_pairs_path=historical_pairs,
        future_pairs_path=future_pairs,
        scaffold_table_path=scaffold_table,
        output_dir=output_dir,
        coverage_policy=args.historical_coverage_policy,
        require_strict_decision=not args.allow_unverified_decision,
    )
    print(output_dir / "summary.json")
    print(
        "eligible scaffold-cold future pairs under selected policy: "
        f"{summary['counts']['future']['scaffold_cold_pairs_under_selected_policy']}"
    )


if __name__ == "__main__":
    main()
