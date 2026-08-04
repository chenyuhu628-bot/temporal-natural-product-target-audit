#!/usr/bin/env python3
"""Create an auditable, scaffold-only repair table from raw NPASS InChI.

This utility is deliberately narrow.  It reads a prepared compound table and
the corresponding RDKit scaffold table, then considers *only* compounds whose
``scaffold_status`` is ``invalid_smiles``.  For each such full InChIKey it
looks up exact-key rows in explicitly supplied raw NPASS structure tables.

The raw NPASS SMILES is retained as provenance, but is never used to derive a
replacement.  A replacement is eligible only when all of the following hold:

1. RDKit parses the raw NPASS InChI;
2. the InChI-derived molecule regenerates the requested full InChIKey exactly;
3. RDKit serializes that molecule to canonical isomeric SMILES;
4. that replacement SMILES parses and regenerates the same full InChIKey.

All source rows and validation fields are written to a repair-candidate table.
The separate ``--repaired-compounds`` output is a *new* table: it preserves
``original_representative_smiles`` and explicitly identifies whether
``representative_smiles`` is an InChI-validated replacement.  It never alters
the supplied input or any raw NPASS file.  It is appropriate only as a
chemical-structure input for a later scaffold derivation; it makes no claim
about activity, binding, or other pharmacology.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


INVALID_SCAFFOLD_STATUS = "invalid_smiles"
VALIDATED_STATUS = "validated_inchi_to_canonical_isomeric_smiles"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path, mode: str):
    return gzip.open(path, mode, encoding="utf-8", newline="") if path.suffix == ".gz" else path.open(
        mode, encoding="utf-8", newline=""
    )


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing header: {path}")
        return list(reader.fieldnames), list(reader)


def write_tsv_gz(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> None:
    """Write a deterministic gzip TSV, leaving inputs untouched."""
    with path.open("wb") as binary:
        with gzip.GzipFile(filename="", mode="wb", fileobj=binary, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(text, fieldnames=fields, delimiter="\t", extrasaction="ignore")
                writer.writeheader()
                writer.writerows(rows)


def require_columns(path: Path, fields: list[str], required: set[str]) -> None:
    missing = sorted(required - set(fields))
    if missing:
        raise ValueError(f"{path} is missing required columns: {missing}")


def compact_json(values: Iterable[str]) -> str:
    return json.dumps(sorted(set(value for value in values if value)), ensure_ascii=False, separators=(",", ":"))


@dataclass(frozen=True)
class RawStructureRecord:
    source_version: str
    source_row_number: int
    source_np_id: str
    raw_inchi: str
    raw_inchikey_full: str
    raw_smiles: str


def read_matching_raw_structures(path: Path, source_version: str, requested_keys: set[str]) -> dict[str, list[RawStructureRecord]]:
    """Read only exact-full-InChIKey matches from one raw NPASS structure file."""
    matches: dict[str, list[RawStructureRecord]] = defaultdict(list)
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing header: {path}")
        require_columns(path, list(reader.fieldnames), {"np_id", "InChI", "InChIKey", "SMILES"})
        for line_number, row in enumerate(reader, start=2):
            key = (row.get("InChIKey") or "").strip()
            if key not in requested_keys:
                continue
            matches[key].append(
                RawStructureRecord(
                    source_version=source_version,
                    source_row_number=line_number,
                    source_np_id=(row.get("np_id") or "").strip(),
                    raw_inchi=(row.get("InChI") or "").strip(),
                    raw_inchikey_full=key,
                    raw_smiles=(row.get("SMILES") or "").strip(),
                )
            )
    return matches


def validate_candidate(
    *,
    record: RawStructureRecord,
    requested_key: str,
    input_smiles: str,
    scaffold_status: str,
    Chem,
    rdBase,
) -> dict[str, object]:
    """Validate a candidate derived solely from raw InChI, retaining raw SMILES."""
    result: dict[str, object] = {
        "inchikey_full": requested_key,
        "source_version": record.source_version,
        "source_row_number": record.source_row_number,
        "source_np_id": record.source_np_id,
        "input_representative_smiles": input_smiles,
        "input_scaffold_status": scaffold_status,
        "raw_inchi": record.raw_inchi,
        "raw_inchikey_full": record.raw_inchikey_full,
        "raw_smiles": record.raw_smiles,
        "raw_smiles_rdkit_parse_status": "not_checked",
        "inchi_parsed": False,
        "regenerated_inchikey_from_raw_inchi": "",
        "repaired_canonical_isomeric_smiles": "",
        "repaired_smiles_parsed": False,
        "regenerated_inchikey_from_repaired_smiles": "",
        "validation_status": "",
        "validated_for_scaffold_repair": False,
    }
    # This status is diagnostic only; raw SMILES is never used as a repair input.
    with rdBase.BlockLogs():
        raw_smiles_molecule = Chem.MolFromSmiles(record.raw_smiles) if record.raw_smiles else None
    result["raw_smiles_rdkit_parse_status"] = "parseable" if raw_smiles_molecule is not None else "invalid_or_missing"

    if record.raw_inchikey_full != requested_key:
        result["validation_status"] = "raw_inchikey_does_not_match_requested_full_key"
        return result
    if not record.raw_inchi:
        result["validation_status"] = "raw_inchi_missing"
        return result
    with rdBase.BlockLogs():
        molecule = Chem.MolFromInchi(record.raw_inchi)
    if molecule is None:
        result["validation_status"] = "raw_inchi_not_parsed_by_rdkit"
        return result
    result["inchi_parsed"] = True
    regenerated_from_inchi = Chem.MolToInchiKey(molecule)
    result["regenerated_inchikey_from_raw_inchi"] = regenerated_from_inchi
    if regenerated_from_inchi != requested_key:
        result["validation_status"] = "raw_inchi_regenerated_key_mismatch"
        return result
    repaired_smiles = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    result["repaired_canonical_isomeric_smiles"] = repaired_smiles
    with rdBase.BlockLogs():
        repaired_molecule = Chem.MolFromSmiles(repaired_smiles)
    if repaired_molecule is None:
        result["validation_status"] = "inchi_derived_smiles_not_parsed_by_rdkit"
        return result
    result["repaired_smiles_parsed"] = True
    regenerated_from_repaired = Chem.MolToInchiKey(repaired_molecule)
    result["regenerated_inchikey_from_repaired_smiles"] = regenerated_from_repaired
    if regenerated_from_repaired != requested_key:
        result["validation_status"] = "repaired_smiles_regenerated_key_mismatch"
        return result
    result["validation_status"] = VALIDATED_STATUS
    result["validated_for_scaffold_repair"] = True
    return result


def missing_candidate_row(inchikey: str, input_smiles: str, scaffold_status: str) -> dict[str, object]:
    return {
        "inchikey_full": inchikey,
        "source_version": "",
        "source_row_number": "",
        "source_np_id": "",
        "input_representative_smiles": input_smiles,
        "input_scaffold_status": scaffold_status,
        "raw_inchi": "",
        "raw_inchikey_full": "",
        "raw_smiles": "",
        "raw_smiles_rdkit_parse_status": "not_checked",
        "inchi_parsed": False,
        "regenerated_inchikey_from_raw_inchi": "",
        "repaired_canonical_isomeric_smiles": "",
        "repaired_smiles_parsed": False,
        "regenerated_inchikey_from_repaired_smiles": "",
        "validation_status": "no_exact_raw_structure_match",
        "validated_for_scaffold_repair": False,
    }


def parser() -> argparse.ArgumentParser:
    description = "Derive exact-key, InChI-validated scaffold-only repair candidates from raw NPASS structure records."
    cli = argparse.ArgumentParser(description=description)
    cli.add_argument("--compounds", required=True, type=Path, help="Prepared strict cold-start compounds TSV/TSV.GZ.")
    cli.add_argument("--scaffolds", required=True, type=Path, help="Existing RDKit scaffold TSV/TSV.GZ for the same compound set.")
    cli.add_argument("--npass-v2-structures", required=True, type=Path, help="Raw NPASS v2 structure table.")
    cli.add_argument("--npass-v3-structures", required=True, type=Path, help="Raw NPASS v3 structure table.")
    cli.add_argument("--repair-table", required=True, type=Path, help="New provenance-rich repair-candidate TSV.GZ.")
    cli.add_argument("--repaired-compounds", required=True, type=Path, help="New compounds TSV.GZ for a later RDKit scaffold run.")
    cli.add_argument("--summary", required=True, type=Path, help="New JSON audit summary.")
    return cli


def main() -> None:
    args = parser().parse_args()
    try:
        from rdkit import Chem, rdBase
    except ImportError as exc:
        raise SystemExit("RDKit is required. Run with the project RDKit environment, then rerun unchanged.") from exc

    for input_path in (args.compounds, args.scaffolds, args.npass_v2_structures, args.npass_v3_structures):
        if not input_path.is_file():
            raise FileNotFoundError(input_path)
    outputs = (args.repair_table, args.repaired_compounds, args.summary)
    for output_path in outputs:
        if output_path.exists():
            raise FileExistsError(f"Refusing to overwrite an existing output: {output_path}")
        if not output_path.parent.is_dir():
            raise FileNotFoundError(f"Output directory does not exist: {output_path.parent}")

    compound_fields, compound_rows = read_tsv(args.compounds)
    require_columns(args.compounds, compound_fields, {"inchikey_full", "representative_smiles"})
    compounds: dict[str, dict[str, str]] = {}
    for row in compound_rows:
        key = row["inchikey_full"].strip()
        if not key:
            raise ValueError(f"Blank full InChIKey in {args.compounds}")
        if key in compounds:
            raise ValueError(f"Duplicate full InChIKey in {args.compounds}: {key}")
        compounds[key] = row

    scaffold_fields, scaffold_rows = read_tsv(args.scaffolds)
    require_columns(args.scaffolds, scaffold_fields, {"inchikey_full", "scaffold_status"})
    scaffold_status: dict[str, str] = {}
    for row in scaffold_rows:
        key = row["inchikey_full"].strip()
        status = row["scaffold_status"].strip()
        if not key:
            raise ValueError(f"Blank full InChIKey in {args.scaffolds}")
        if key in scaffold_status and scaffold_status[key] != status:
            raise ValueError(f"Conflicting scaffold statuses for {key}")
        scaffold_status[key] = status
    missing_scaffold = sorted(set(compounds) - set(scaffold_status))
    if missing_scaffold:
        raise ValueError(f"Scaffold table lacks {len(missing_scaffold)} compounds; first: {missing_scaffold[:5]}")

    invalid_keys = sorted(key for key in compounds if scaffold_status[key] == INVALID_SCAFFOLD_STATUS)
    v2_matches = read_matching_raw_structures(args.npass_v2_structures, "v2", set(invalid_keys))
    v3_matches = read_matching_raw_structures(args.npass_v3_structures, "v3", set(invalid_keys))

    candidate_rows: list[dict[str, object]] = []
    candidates_by_key: dict[str, list[dict[str, object]]] = defaultdict(list)
    for key in invalid_keys:
        raw_records = [*v2_matches.get(key, []), *v3_matches.get(key, [])]
        if not raw_records:
            candidate = missing_candidate_row(key, compounds[key]["representative_smiles"], scaffold_status[key])
            candidate_rows.append(candidate)
            candidates_by_key[key].append(candidate)
            continue
        for record in raw_records:
            candidate = validate_candidate(
                record=record,
                requested_key=key,
                input_smiles=compounds[key]["representative_smiles"],
                scaffold_status=scaffold_status[key],
                Chem=Chem,
                rdBase=rdBase,
            )
            candidate_rows.append(candidate)
            candidates_by_key[key].append(candidate)

    candidate_fields = [
        "inchikey_full", "source_version", "source_row_number", "source_np_id", "input_representative_smiles",
        "input_scaffold_status", "raw_inchi", "raw_inchikey_full", "raw_smiles", "raw_smiles_rdkit_parse_status",
        "inchi_parsed", "regenerated_inchikey_from_raw_inchi", "repaired_canonical_isomeric_smiles",
        "repaired_smiles_parsed", "regenerated_inchikey_from_repaired_smiles", "validation_status",
        "validated_for_scaffold_repair",
    ]

    selection: dict[str, dict[str, object]] = {}
    for key in invalid_keys:
        candidates = candidates_by_key[key]
        valid = [row for row in candidates if row["validation_status"] == VALIDATED_STATUS]
        candidate_smiles = sorted({str(row["repaired_canonical_isomeric_smiles"]) for row in valid})
        if not valid:
            selection[key] = {
                "repair_status": "no_validated_inchi_repair",
                "repair_applied": False,
                "repaired_smiles": "",
                "candidate_rows": candidates,
            }
        elif len(candidate_smiles) != 1:
            selection[key] = {
                "repair_status": "ambiguous_validated_inchi_repairs_not_applied",
                "repair_applied": False,
                "repaired_smiles": "",
                "candidate_rows": candidates,
            }
        else:
            selection[key] = {
                "repair_status": "validated_inchi_repair_applied_to_new_output_only",
                "repair_applied": True,
                "repaired_smiles": candidate_smiles[0],
                "candidate_rows": candidates,
            }

    repair_columns = [
        "original_representative_smiles", "structure_for_scaffold_source", "structure_repair_status",
        "structure_repair_applied", "repair_validation_rule", "repair_source_versions", "repair_source_np_ids",
        "repair_source_row_numbers", "repair_source_raw_inchis", "repair_source_raw_smiles",
        "repair_regenerated_inchikey", "validated_repaired_smiles",
    ]
    repaired_rows: list[dict[str, object]] = []
    for key, original in compounds.items():
        item: dict[str, object] = dict(original)
        original_smiles = original["representative_smiles"]
        choice = selection.get(key)
        if choice is None:
            item.update(
                {
                    "original_representative_smiles": original_smiles,
                    "structure_for_scaffold_source": "original_compounds_input",
                    "structure_repair_status": "not_selected__scaffold_status_not_invalid_smiles",
                    "structure_repair_applied": False,
                    "repair_validation_rule": "not_applicable",
                    "repair_source_versions": "[]",
                    "repair_source_np_ids": "[]",
                    "repair_source_row_numbers": "[]",
                    "repair_source_raw_inchis": "[]",
                    "repair_source_raw_smiles": "[]",
                    "repair_regenerated_inchikey": "",
                    "validated_repaired_smiles": "",
                }
            )
        else:
            candidate_rows_for_key = choice["candidate_rows"]
            valid = [row for row in candidate_rows_for_key if row["validation_status"] == VALIDATED_STATUS]
            repaired_smiles = str(choice["repaired_smiles"])
            item["representative_smiles"] = repaired_smiles if choice["repair_applied"] else original_smiles
            item.update(
                {
                    "original_representative_smiles": original_smiles,
                    "structure_for_scaffold_source": "validated_npass_inchi_repair" if choice["repair_applied"] else "original_compounds_input",
                    "structure_repair_status": choice["repair_status"],
                    "structure_repair_applied": choice["repair_applied"],
                    "repair_validation_rule": (
                        "raw_npass_inchi_parses__both_regenerated_full_inchikeys_equal_requested_key__all_validated_candidates_agree"
                        if choice["repair_applied"]
                        else "no_replacement_without_unambiguous_full_inchikey_validation"
                    ),
                    "repair_source_versions": compact_json(str(row["source_version"]) for row in valid),
                    "repair_source_np_ids": compact_json(str(row["source_np_id"]) for row in valid),
                    "repair_source_row_numbers": compact_json(str(row["source_row_number"]) for row in valid),
                    "repair_source_raw_inchis": compact_json(str(row["raw_inchi"]) for row in valid),
                    "repair_source_raw_smiles": compact_json(str(row["raw_smiles"]) for row in valid),
                    "repair_regenerated_inchikey": compact_json(
                        str(row["regenerated_inchikey_from_repaired_smiles"]) for row in valid
                    ),
                    "validated_repaired_smiles": repaired_smiles,
                }
            )
        repaired_rows.append(item)

    repaired_fields = [*compound_fields, *(field for field in repair_columns if field not in compound_fields)]
    write_tsv_gz(args.repair_table, candidate_fields, candidate_rows)
    write_tsv_gz(args.repaired_compounds, repaired_fields, repaired_rows)

    candidate_statuses = Counter(str(row["validation_status"]) for row in candidate_rows)
    selection_statuses = Counter(str(value["repair_status"]) for value in selection.values())
    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "script": Path(__file__).name,
        "rdkit_version": rdBase.rdkitVersion,
        "scope": {
            "selected_only_when_existing_scaffold_status_equals": INVALID_SCAFFOLD_STATUS,
            "repair_source": "raw NPASS InChI only; raw NPASS SMILES is preserved for provenance and diagnostic parsing only",
            "identity_rule": "both InChI-derived molecule and its canonical isomeric SMILES must regenerate the requested full InChIKey exactly",
            "replacement_rule": "apply only if every validated raw-source candidate agrees on one replacement SMILES",
            "non_claim": "The output is a scaffold-input repair audit and makes no pharmacological or interaction claim.",
        },
        "inputs": {
            "compounds": {"path": str(args.compounds), "sha256": sha256(args.compounds)},
            "scaffolds": {"path": str(args.scaffolds), "sha256": sha256(args.scaffolds)},
            "npass_v2_structures": {"path": str(args.npass_v2_structures), "sha256": sha256(args.npass_v2_structures)},
            "npass_v3_structures": {"path": str(args.npass_v3_structures), "sha256": sha256(args.npass_v3_structures)},
        },
        "outputs": {
            "repair_table": {"path": str(args.repair_table), "sha256": sha256(args.repair_table)},
            "repaired_compounds": {"path": str(args.repaired_compounds), "sha256": sha256(args.repaired_compounds)},
        },
        "counts": {
            "input_compounds": len(compounds),
            "input_scaffold_rows": len(scaffold_status),
            "selected_invalid_smiles_compounds": len(invalid_keys),
            "raw_exact_matches_v2": sum(len(rows) for rows in v2_matches.values()),
            "raw_exact_matches_v3": sum(len(rows) for rows in v3_matches.values()),
            "repair_candidate_rows": len(candidate_rows),
            "candidate_validation_statuses": dict(sorted(candidate_statuses.items())),
            "repair_selection_statuses": dict(sorted(selection_statuses.items())),
            "compounds_with_validated_repair_applied_to_new_output": sum(
                bool(value["repair_applied"]) for value in selection.values()
            ),
            "compounds_left_with_original_smiles": len(compounds) - sum(
                bool(value["repair_applied"]) for value in selection.values()
            ),
        },
        "selected_invalid_smiles_full_inchikeys": invalid_keys,
    }
    args.summary.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(args.summary)


if __name__ == "__main__":
    main()

