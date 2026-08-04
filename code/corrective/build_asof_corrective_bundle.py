"""Build the restricted row-level as-of-cutoff corrective input bundle.

The builder keeps the frozen 4,990 historical pair keys and 358-relation
endpoint, but reconstructs every historical tier, weight, structure, fingerprint
and scaffold from day-precise v2 rows proven to be on or before 2022-08-31.
Historical and query compound representations remain role separated.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import sys
from collections import Counter, defaultdict
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any

from asof_common import (
    assert_unique,
    membership_sha256,
    parse_bool,
    read_rows,
    read_table,
    require_fields,
    sha256,
    write_json,
    write_tsv_gz,
)


PROTOCOL_ID = "npass_strict_ab_asof_cutoff_corrective_successor_v1_20260728"
RUN_ID = "npass_strict_ab_asof_cutoff_author_run_v1_20260728"
CUTOFF = date(2022, 8, 31)
P1 = "P1_npass_raw_exact_candidate"
STRICT_MAPPING = "strict_one_to_one_reviewed_human"
TIERS = {"A_affinity_candidate", "B_quantitative_functional_candidate"}
WEIGHTS = {"A_affinity_candidate": 1.0, "B_quantitative_functional_candidate": 0.7}
PROTOCOL_LOCK_SHA256 = "96befee13ae1d41ad433c8697fac92ccd30fb25e24c3cf1279c6b4b7e040abd9"

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
PROTOCOL_LOCK = ROOT / "manifests" / "protocol_lock_manifest_v1.json"
SPEC = ROOT / "configs" / "asof_cutoff_rebuild_spec_v1.json"

SOURCES = {
    "v2_evidence": WORKSPACE / "data/processed/npass_v2_evidence_records_v1_1_uniprot_mapped.tsv.gz",
    "v3_evidence": WORKSPACE / "data/processed/npass_v3_evidence_records_v1_1_uniprot_mapped.tsv.gz",
    "v2_pubmed": WORKSPACE / "data/interim/pubmed_v2_pmid_metadata.csv.gz",
    "v2_raw_structures": WORKSPACE / "data/raw/npass2/NPASSv2.0_download_naturalProducts_structureInfo.txt",
    "old_history": WORKSPACE / "author_run_strict_ab_doublecold_bundle_v2_20260720/scoring_inputs/historical_pairs.tsv.gz",
    "old_queries": WORKSPACE / "author_run_strict_ab_doublecold_bundle_v2_20260720/scoring_inputs/scoring_queries.tsv.gz",
    "old_compounds": WORKSPACE / "author_run_strict_ab_doublecold_bundle_v2_20260720/scoring_inputs/compounds.tsv.gz",
    "candidate_targets": WORKSPACE / "author_run_strict_ab_doublecold_bundle_v2_20260720/scoring_inputs/candidate_targets.tsv.gz",
    "candidate_sequences": WORKSPACE / "author_run_strict_ab_doublecold_bundle_v2_20260720/scoring_inputs/candidate_sequences.fasta",
    "endpoint": WORKSPACE / "author_run_strict_ab_doublecold_bundle_v2_20260720/evaluation_inputs/evaluation_pairs.tsv.gz",
    "old_scaffold_audit": WORKSPACE / "author_run_strict_ab_doublecold_bundle_v2_20260720/evaluation_inputs/scaffold_audit.tsv.gz",
    "homology_0_30": WORKSPACE / "author_run_strict_ab_doublecold_bundle_v2_20260720/evaluation_inputs/homology_0_30.tsv.gz",
    "homology_0_50": WORKSPACE / "author_run_strict_ab_doublecold_bundle_v2_20260720/evaluation_inputs/homology_0_50.tsv.gz",
    "homology_0_70": WORKSPACE / "author_run_strict_ab_doublecold_bundle_v2_20260720/evaluation_inputs/homology_0_70.tsv.gz",
    "catalogue_source": WORKSPACE / "data/processed/chembl31_human_target_catalogue_v1/chembl31_human_single_protein_targets_v1.tsv.gz",
    "sequence_source": WORKSPACE / "data/processed/chembl31_human_target_catalogue_v1/chembl31_human_single_protein_targets_v1.fasta",
}

EXPECTED_HASHES = {
    "v2_evidence": "69802f62470930405d2658cc37e99f668aeb7db7795696ff78c76563d552d582",
    "v3_evidence": "cc80a1cfe11c4bfd007cf0a9b110026429bffc99b23abff500d2cfddeb7fadbe",
    "v2_pubmed": "f94a892fd9ce106a62e6306b16a0a2f0b5482987c1610d33e6cba21e090029b0",
    "v2_raw_structures": "cbe688e9b6fdd0960c78d1a93d7f487ca4a2bbad3017275515c3d90d2a0f72fa",
    "old_history": "75a01dc27fbffd677154865f6817e51f00d8f0dacdbbe4cf8ef7dbcf31a2e959",
    "old_queries": "0e6068d2e25cb3ea325656fb3517563788cd496e88cfaa3de761890fec9e9318",
    "old_compounds": "80c4c654d8726202e7221d561c5f9f5a9e94b21a1b7009128b1d59cb32aa2674",
    "candidate_targets": "0ee86746b306fb388a1f74a6b88ce4d1eba01b7a4eb473315f6b3def57145cdc",
    "candidate_sequences": "a83421dba2482f236fe18340dd592cc7d5ed22c98c4fc39435c40f04f289b442",
    "endpoint": "09296b066a23197a7c178f00514f2b3d9ed7e6f3c459ea92a55e01a6010d1132",
    "old_scaffold_audit": "fa0029ef5b7822ad5ca93f7bd93ac808f85f1e0c02e827fa91be375031b2d7af",
    "homology_0_30": "3a8247ed8f683fe6fce5fb345f56e3ec73a872b065eca922e92e494f084a1793",
    "homology_0_50": "ec3bbd435f73bc1c724efdfd94ac10c32d6f9a55bd6c9a9349233a01e5dc7f5b",
    "homology_0_70": "ec3bbd435f73bc1c724efdfd94ac10c32d6f9a55bd6c9a9349233a01e5dc7f5b",
    "catalogue_source": "9d50f497b5304028585b12506f52a5377301db8c6eb14d9da16ac8c24e3cdc3a",
    "sequence_source": "e8c3f0b17a3231853b1e86d4a63c232a716692f62d38c238ce0eb90fffb676bc",
}


def canonical_pair(row: dict[str, str]) -> str:
    return f"{row.get('inchikey_full', '').strip()}|{row.get('uniprot_canonical_accession', '').strip()}"


def is_strict_primary(row: dict[str, str]) -> bool:
    return (
        row.get("automatic_verification_level", "") == P1
        and row.get("evidence_tier_v1_1", "") in TIERS
        and bool(row.get("inchikey_full", "").strip())
        and bool(row.get("uniprot_canonical_accession", "").strip())
        and row.get("mapping_status", "") == STRICT_MAPPING
        and row.get("sequence_found", "").strip().casefold() == "true"
    )


def best_tier(values: set[str]) -> str:
    if not values or not values.issubset(TIERS):
        raise ValueError(f"Invalid eligible tier set: {sorted(values)}")
    return "A_affinity_candidate" if "A_affinity_candidate" in values else "B_quantitative_functional_candidate"


def modal_smiles(counts: Counter[str]) -> str:
    if not counts:
        return ""
    top = max(counts.values())
    return sorted(value for value, count in counts.items() if count == top)[0]


def validate_sources() -> None:
    if sha256(PROTOCOL_LOCK) != PROTOCOL_LOCK_SHA256:
        raise ValueError("Protocol lock hash differs from the pre-result lock")
    protocol = json.loads(PROTOCOL_LOCK.read_text(encoding="utf-8"))
    if protocol.get("protocol_id") != PROTOCOL_ID:
        raise ValueError("Protocol lock ID mismatch")
    for name, path in SOURCES.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        actual = sha256(path)
        if actual != EXPECTED_HASHES[name]:
            raise ValueError(f"Frozen source hash mismatch for {name}: {actual}")


def load_pubmed() -> dict[str, dict[str, str]]:
    fields, rows = read_table(SOURCES["v2_pubmed"], ",")
    require_fields(
        fields,
        {"pmid", "found_in_pubmed", "publication_date", "date_precision", "date_source"},
        "v2 PubMed metadata",
    )
    assert_unique(rows, ("pmid",), "v2 PubMed metadata")
    return {row["pmid"].strip(): row for row in rows}


def classify_row(row: dict[str, str], pubmed: dict[str, dict[str, str]]) -> tuple[str, dict[str, str]]:
    ref_type = row.get("ref_id_type", "").strip().upper()
    ref_id = row.get("ref_id", "").strip()
    meta = pubmed.get(ref_id, {})
    details = {
        "publication_date": meta.get("publication_date", "").strip(),
        "date_precision": meta.get("date_precision", "missing").strip() or "missing",
        "date_source": meta.get("date_source", "").strip(),
        "found_in_pubmed": meta.get("found_in_pubmed", "False").strip() or "False",
    }
    if ref_type != "PMID" or not ref_id.isdigit():
        return "excluded_not_numeric_pmid", details
    if not meta or not parse_bool(details["found_in_pubmed"]):
        return "excluded_pubmed_not_found", details
    if details["date_precision"] != "day" or not details["publication_date"]:
        return "excluded_non_day_precision", details
    try:
        published = date.fromisoformat(details["publication_date"])
    except ValueError as exc:
        raise ValueError(f"Invalid archived day-precise publication date for PMID {ref_id}") from exc
    return ("eligible_pre_cutoff" if published <= CUTOFF else "excluded_after_cutoff"), details


def load_raw_v2_records(requested: set[str]) -> dict[str, list[dict[str, str]]]:
    matches: dict[str, list[dict[str, str]]] = defaultdict(list)
    with SOURCES["v2_raw_structures"].open("rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require_fields(reader.fieldnames or [], {"np_id", "InChI", "InChIKey", "SMILES"}, "raw NPASS v2 structures")
        for row_number, row in enumerate(reader, start=2):
            key = row.get("InChIKey", "").strip()
            if key not in requested:
                continue
            matches[key].append(
                {
                    "source_row_number": str(row_number),
                    "source_np_id": row.get("np_id", "").strip(),
                    "raw_inchi": row.get("InChI", "").strip(),
                    "raw_smiles": row.get("SMILES", "").strip(),
                }
            )
    return matches


def scaffold_record(smiles: str, Chem: Any, MurckoScaffold: Any, rdBase: Any) -> dict[str, str]:
    with rdBase.BlockLogs():
        molecule = Chem.MolFromSmiles(smiles) if smiles else None
    if molecule is None:
        return {"scaffold_status": "invalid_smiles", "bemis_murcko_smiles": "", "scaffold_key": ""}
    scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
    scaffold_smiles = Chem.MolToSmiles(scaffold, canonical=True) if scaffold.GetNumAtoms() else ""
    if not scaffold_smiles:
        return {
            "scaffold_status": "acyclic_or_empty_bemis_murcko",
            "bemis_murcko_smiles": "",
            "scaffold_key": "",
        }
    return {"scaffold_status": "ok", "bemis_murcko_smiles": scaffold_smiles, "scaffold_key": scaffold_smiles}


def validated_v2_repair(
    inchikey: str,
    allowed_np_ids: set[str],
    raw_records: list[dict[str, str]],
    Chem: Any,
    rdBase: Any,
) -> tuple[str, list[dict[str, str]]]:
    audit: list[dict[str, str]] = []
    validated: list[tuple[str, dict[str, str]]] = []
    for record in raw_records:
        status = "source_np_id_not_eligible"
        repaired = ""
        if record["source_np_id"] in allowed_np_ids:
            status = "raw_inchi_missing"
            if record["raw_inchi"]:
                with rdBase.BlockLogs():
                    molecule = Chem.MolFromInchi(record["raw_inchi"])
                if molecule is None:
                    status = "raw_inchi_not_parsed"
                elif Chem.MolToInchiKey(molecule) != inchikey:
                    status = "raw_inchi_key_mismatch"
                else:
                    repaired = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
                    with rdBase.BlockLogs():
                        roundtrip = Chem.MolFromSmiles(repaired)
                    if roundtrip is None or Chem.MolToInchiKey(roundtrip) != inchikey:
                        status = "repaired_smiles_key_mismatch"
                        repaired = ""
                    else:
                        status = "validated_v2_eligible_inchi_repair"
                        validated.append((repaired, record))
        audit.append(
            {
                "inchikey_full": inchikey,
                "source_version": "v2",
                "source_np_id": record["source_np_id"],
                "source_row_number": record["source_row_number"],
                "eligible_source_np_id": str(record["source_np_id"] in allowed_np_ids),
                "validation_status": status,
                "repaired_smiles": repaired,
            }
        )
    choices = sorted({value for value, _ in validated})
    if len(choices) != 1:
        raise ValueError(
            f"Historical structure repair for {inchikey} is not a unique eligible-v2 exact-key structure: {len(choices)} choices"
        )
    return choices[0], audit


def copy_exact(source: Path, target: Path) -> None:
    if target.exists():
        raise FileExistsError(f"Refusing to overwrite output: {target}")
    shutil.copyfile(source, target)
    if sha256(source) != sha256(target):
        raise IOError(f"Byte-exact copy validation failed: {source} -> {target}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--code-lock-manifest", required=True, type=Path)
    args = parser.parse_args()

    try:
        from rdkit import Chem, DataStructs, rdBase, rdFingerprintGenerator
        from rdkit.Chem.Scaffolds import MurckoScaffold
    except ImportError as exc:
        raise SystemExit("The locked RDKit environment is required") from exc

    output_root = args.output_root.resolve()
    if output_root.exists():
        raise FileExistsError(f"Corrective output root already exists: {output_root}")
    if output_root.parent != WORKSPACE.resolve():
        raise ValueError(f"Corrective output root must be a new direct child of {WORKSPACE}")
    if not args.code_lock_manifest.is_file():
        raise FileNotFoundError(args.code_lock_manifest)
    code_lock = json.loads(args.code_lock_manifest.read_text(encoding="utf-8"))
    if code_lock.get("protocol_id") != PROTOCOL_ID or code_lock.get("lock_state") != "LOCKED_BEFORE_EXECUTION":
        raise ValueError("Implementation code lock is absent or not executable")

    validate_sources()
    pubmed = load_pubmed()

    old_history_fields, old_history_rows = read_table(SOURCES["old_history"])
    require_fields(
        old_history_fields,
        {"canonical_pair_key", "inchikey_full", "uniprot_canonical_accession", "best_strict_evidence_tier"},
        "frozen historical pair table",
    )
    if len(old_history_rows) != 4990:
        raise ValueError("Frozen historical pair cardinality is not 4,990")
    assert_unique(old_history_rows, ("canonical_pair_key",), "frozen historical pair table")
    old_history = {row["canonical_pair_key"]: row for row in old_history_rows}
    historical_keys = set(old_history)

    endpoint_fields, endpoint_rows = read_table(SOURCES["endpoint"])
    require_fields(
        endpoint_fields,
        {"canonical_pair_key", "query_id", "inchikey_full", "uniprot_canonical_accession", "best_strict_evidence_tier"},
        "frozen endpoint",
    )
    if len(endpoint_rows) != 358:
        raise ValueError("Frozen endpoint relation count is not 358")
    assert_unique(endpoint_rows, ("canonical_pair_key",), "frozen endpoint")
    endpoint_keys = {row["canonical_pair_key"] for row in endpoint_rows}
    if historical_keys.intersection(endpoint_keys):
        raise ValueError("Frozen historical and endpoint keysets overlap")
    query_ids = {row["query_id"] for row in endpoint_rows}
    query_compound_ids = {row["inchikey_full"] for row in endpoint_rows}
    endpoint_targets = {row["uniprot_canonical_accession"] for row in endpoint_rows}
    if (len(query_ids), len(query_compound_ids), len(endpoint_targets)) != (222, 222, 156):
        raise ValueError("Frozen endpoint 222-query/156-target contract failed")

    query_fields, query_rows = read_table(SOURCES["old_queries"])
    require_fields(query_fields, {"query_id", "inchikey_full"}, "frozen query map")
    assert_unique(query_rows, ("query_id",), "frozen query map")
    assert_unique(query_rows, ("inchikey_full",), "frozen query map")
    if {(row["query_id"], row["inchikey_full"]) for row in query_rows} != {
        (row["query_id"], row["inchikey_full"]) for row in endpoint_rows
    }:
        raise ValueError("Frozen query map differs from endpoint query mapping")

    old_compound_fields, old_compound_rows = read_table(SOURCES["old_compounds"])
    require_fields(old_compound_fields, {"inchikey_full", "representative_smiles"}, "old compound map")
    assert_unique(old_compound_rows, ("inchikey_full",), "old compound map")
    old_smiles = {row["inchikey_full"]: row["representative_smiles"] for row in old_compound_rows}

    row_ledger: list[dict[str, object]] = []
    row_status_counts: Counter[str] = Counter()
    eligible_rows_by_pair: dict[str, list[dict[str, str]]] = defaultdict(list)
    excluded_rows_by_pair: Counter[str] = Counter()
    history_smiles_counts: dict[str, Counter[str]] = defaultdict(Counter)
    eligible_np_ids_by_compound: dict[str, set[str]] = defaultdict(set)
    historical_source_row_number = 1

    for row in read_rows(SOURCES["v2_evidence"]):
        historical_source_row_number += 1
        pair_key = canonical_pair(row)
        if pair_key not in historical_keys or not is_strict_primary(row):
            continue
        status, date_details = classify_row(row, pubmed)
        row_status_counts[status] += 1
        ledger_row = {
            "canonical_pair_key": pair_key,
            "source_version": "v2",
            "source_row_number": historical_source_row_number,
            "source_np_id": row.get("source_np_id", ""),
            "source_target_id": row.get("source_target_id", ""),
            "inchikey_full": row.get("inchikey_full", ""),
            "uniprot_canonical_accession": row.get("uniprot_canonical_accession", ""),
            "evidence_tier_v1_1": row.get("evidence_tier_v1_1", ""),
            "activity_type": row.get("activity_type", ""),
            "ref_id_type": row.get("ref_id_type", ""),
            "ref_id": row.get("ref_id", ""),
            "found_in_pubmed": date_details["found_in_pubmed"],
            "publication_date": date_details["publication_date"],
            "date_precision": date_details["date_precision"],
            "date_source": date_details["date_source"],
            "row_eligibility_status": status,
            "smiles": row.get("smiles", ""),
        }
        row_ledger.append(ledger_row)
        if status == "eligible_pre_cutoff":
            eligible_rows_by_pair[pair_key].append(row)
            compound = row["inchikey_full"].strip()
            eligible_np_ids_by_compound[compound].add(row.get("source_np_id", "").strip())
            smiles = row.get("smiles", "").strip()
            if smiles:
                history_smiles_counts[compound][smiles] += 1
        else:
            excluded_rows_by_pair[pair_key] += 1

    expected_row_counts = {
        "eligible_pre_cutoff": 13885,
        "excluded_non_day_precision": 6570,
        "excluded_not_numeric_pmid": 192,
        "excluded_pubmed_not_found": 0,
        "excluded_after_cutoff": 0,
    }
    if {key: row_status_counts.get(key, 0) for key in expected_row_counts} != expected_row_counts:
        raise ValueError(f"Row eligibility counts differ from the locked precheck: {dict(row_status_counts)}")
    normalized_row_status_counts = {
        key: row_status_counts.get(key, 0) for key in expected_row_counts
    }
    if set(eligible_rows_by_pair) != historical_keys:
        missing = sorted(historical_keys.difference(eligible_rows_by_pair))
        extra = sorted(set(eligible_rows_by_pair).difference(historical_keys))
        raise ValueError(f"Corrected history does not retain exactly 4,990 keys; missing={len(missing)}, extra={len(extra)}")

    corrected_history: list[dict[str, object]] = []
    pair_diff_rows: list[dict[str, object]] = []
    tier_change_count = 0
    for pair_key in sorted(historical_keys):
        old = old_history[pair_key]
        eligible = eligible_rows_by_pair[pair_key]
        tier = best_tier({row["evidence_tier_v1_1"] for row in eligible})
        changed = tier != old["best_strict_evidence_tier"]
        tier_change_count += changed
        corrected_history.append(
            {
                "canonical_pair_key": pair_key,
                "inchikey_full": old["inchikey_full"],
                "uniprot_canonical_accession": old["uniprot_canonical_accession"],
                "best_strict_evidence_tier": tier,
                "eligible_pre_cutoff_v2_row_count": len(eligible),
                "decision": "strict_pre_cutoff_training_candidate",
                "unrecorded_pair_policy": "unlabeled_not_negative",
            }
        )
        pair_diff_rows.append(
            {
                "canonical_pair_key": pair_key,
                "inchikey_full": old["inchikey_full"],
                "uniprot_canonical_accession": old["uniprot_canonical_accession"],
                "old_best_strict_evidence_tier": old["best_strict_evidence_tier"],
                "corrected_best_strict_evidence_tier": tier,
                "old_weight": WEIGHTS[old["best_strict_evidence_tier"]],
                "corrected_weight": WEIGHTS[tier],
                "tier_or_weight_changed": changed,
                "eligible_pre_cutoff_row_count": len(eligible),
                "excluded_row_count": excluded_rows_by_pair[pair_key],
            }
        )
    if tier_change_count != 166:
        raise ValueError(f"Corrected tier changes differ from locked precheck: {tier_change_count}")

    endpoint_strict_rows: list[dict[str, str]] = []
    query_smiles_counts: dict[str, Counter[str]] = defaultdict(Counter)
    for row in read_rows(SOURCES["v3_evidence"]):
        pair_key = canonical_pair(row)
        if pair_key not in endpoint_keys or not is_strict_primary(row):
            continue
        if row.get("cross_version_status", "") != "v3_only" or row.get("temporal_screen_status", "") != "future_candidate_pmid_only":
            raise ValueError("A frozen endpoint source row fails the v3-only future temporal contract")
        endpoint_strict_rows.append(row)
        smiles = row.get("smiles", "").strip()
        if smiles:
            query_smiles_counts[row["inchikey_full"].strip()][smiles] += 1
    if set(query_smiles_counts) != query_compound_ids:
        raise ValueError("Every frozen query must have at least one nonempty strict v3 endpoint-side SMILES")

    history_compound_ids = {row["inchikey_full"] for row in corrected_history}
    if set(history_smiles_counts) != history_compound_ids:
        raise ValueError("Every historical compound must have at least one nonempty eligible-v2 SMILES")

    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    historical_representation: dict[str, str] = {
        compound: modal_smiles(counts) for compound, counts in history_smiles_counts.items()
    }
    query_representation: dict[str, str] = {
        compound: modal_smiles(counts) for compound, counts in query_smiles_counts.items()
    }

    invalid_history: set[str] = set()
    with rdBase.BlockLogs():
        for compound, smiles in historical_representation.items():
            if Chem.MolFromSmiles(smiles) is None:
                invalid_history.add(compound)
        invalid_query = {
            compound for compound, smiles in query_representation.items() if Chem.MolFromSmiles(smiles) is None
        }
    if invalid_query:
        raise ValueError(f"Protocol does not permit an unplanned query-side repair; invalid query compounds={len(invalid_query)}")

    raw_v2 = load_raw_v2_records(invalid_history)
    repair_audit_rows: list[dict[str, str]] = []
    repaired_history: set[str] = set()
    for compound in sorted(invalid_history):
        repaired, audit = validated_v2_repair(
            compound,
            eligible_np_ids_by_compound[compound],
            raw_v2.get(compound, []),
            Chem,
            rdBase,
        )
        historical_representation[compound] = repaired
        repaired_history.add(compound)
        repair_audit_rows.extend(audit)

    historical_scaffolds: dict[str, dict[str, str]] = {}
    query_scaffolds: dict[str, dict[str, str]] = {}
    historical_fingerprints: dict[str, Any] = {}
    query_fingerprints: dict[str, Any] = {}
    for compound, smiles in historical_representation.items():
        with rdBase.BlockLogs():
            molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError("A corrected historical structure remains unparsable")
        historical_fingerprints[compound] = generator.GetFingerprint(molecule)
        historical_scaffolds[compound] = scaffold_record(smiles, Chem, MurckoScaffold, rdBase)
    for compound, smiles in query_representation.items():
        with rdBase.BlockLogs():
            molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError("A corrected query structure remains unparsable")
        query_fingerprints[compound] = generator.GetFingerprint(molecule)
        query_scaffolds[compound] = scaffold_record(smiles, Chem, MurckoScaffold, rdBase)

    old_fingerprints: dict[str, Any] = {}
    old_scaffolds: dict[str, dict[str, str]] = {}
    for compound in history_compound_ids.union(query_compound_ids):
        old_value = old_smiles.get(compound, "")
        with rdBase.BlockLogs():
            molecule = Chem.MolFromSmiles(old_value) if old_value else None
        if molecule is None:
            raise ValueError("Old frozen compound map contains an unparsable structure after its recorded repair")
        old_fingerprints[compound] = generator.GetFingerprint(molecule)
        old_scaffolds[compound] = scaffold_record(old_value, Chem, MurckoScaffold, rdBase)

    historical_compound_rows: list[dict[str, object]] = []
    query_compound_rows: list[dict[str, object]] = []
    structure_audit_rows: list[dict[str, object]] = []
    for compound in sorted(history_compound_ids):
        corrected = historical_representation[compound]
        old_value = old_smiles[compound]
        fp_equal = bool(DataStructs.TanimotoSimilarity(historical_fingerprints[compound], old_fingerprints[compound]) == 1.0)
        scaffold_equal = historical_scaffolds[compound] == old_scaffolds[compound]
        source = "validated_v2_eligible_inchi_repair" if compound in repaired_history else "mode_of_eligible_pre_cutoff_v2_strict_rows"
        row = {
            "inchikey_full": compound,
            "representative_smiles": corrected,
            "structure_role": "historical",
            "selection_rule": source,
        }
        historical_compound_rows.append(row)
        structure_audit_rows.append(
            {
                **row,
                "old_combined_representative_smiles": old_value,
                "representative_smiles_changed": corrected != old_value,
                "morgan_2048_radius2_equal_to_old": fp_equal,
                "scaffold_status": historical_scaffolds[compound]["scaffold_status"],
                "scaffold_key": historical_scaffolds[compound]["scaffold_key"],
                "scaffold_assignment_equal_to_old": scaffold_equal,
                "nonempty_source_smiles_row_count": sum(history_smiles_counts[compound].values()),
                "distinct_source_smiles_count": len(history_smiles_counts[compound]),
            }
        )
    for compound in sorted(query_compound_ids):
        corrected = query_representation[compound]
        old_value = old_smiles[compound]
        fp_equal = bool(DataStructs.TanimotoSimilarity(query_fingerprints[compound], old_fingerprints[compound]) == 1.0)
        scaffold_equal = query_scaffolds[compound] == old_scaffolds[compound]
        row = {
            "inchikey_full": compound,
            "representative_smiles": corrected,
            "structure_role": "query",
            "selection_rule": "mode_of_frozen_endpoint_v3_strict_rows",
        }
        query_compound_rows.append(row)
        structure_audit_rows.append(
            {
                **row,
                "old_combined_representative_smiles": old_value,
                "representative_smiles_changed": corrected != old_value,
                "morgan_2048_radius2_equal_to_old": fp_equal,
                "scaffold_status": query_scaffolds[compound]["scaffold_status"],
                "scaffold_key": query_scaffolds[compound]["scaffold_key"],
                "scaffold_assignment_equal_to_old": scaffold_equal,
                "nonempty_source_smiles_row_count": sum(query_smiles_counts[compound].values()),
                "distinct_source_smiles_count": len(query_smiles_counts[compound]),
            }
        )

    valid_history_scaffold_keys = {
        item["scaffold_key"] for item in historical_scaffolds.values() if item["scaffold_status"] == "ok"
    }
    scaffold_rows: list[dict[str, object]] = []
    for endpoint in sorted(endpoint_rows, key=lambda row: row["canonical_pair_key"]):
        scaffold = query_scaffolds[endpoint["inchikey_full"]]
        if scaffold["scaffold_status"] != "ok":
            cold = False
            outcome = "excluded_future_acyclic_or_empty_bemis_murcko"
            reason = (
                "This future compound cannot enter a Bemis-Murcko scaffold group under the recorded scaffold assignment; "
                "it is not converted into an empty or shared scaffold group."
            )
        elif scaffold["scaffold_key"] in valid_history_scaffold_keys:
            cold = False
            outcome = "not_scaffold_cold__scaffold_seen_in_eligible_historical_compounds"
            reason = "A non-empty, valid scaffold key matching this future compound was observed among historical compounds."
        else:
            cold = True
            outcome = "eligible_scaffold_cold__scaffold_absent_from_historical_nonempty_groups"
            reason = "The non-empty, valid scaffold key is absent from all eligible historical compound scaffold groups."
        scaffold_rows.append(
            {
                "canonical_pair_key": endpoint["canonical_pair_key"],
                "audit_scaffold_cold_under_selected_policy": cold,
                "audit_outcome": outcome,
                "audit_eligibility_or_exclusion_reason": reason,
            }
        )

    old_scaffold_fields, old_scaffold_rows = read_table(SOURCES["old_scaffold_audit"])
    require_fields(
        old_scaffold_fields,
        {"canonical_pair_key", "audit_scaffold_cold_under_selected_policy"},
        "old scaffold audit",
    )
    old_scaffold_flags = {
        row["canonical_pair_key"]: parse_bool(row["audit_scaffold_cold_under_selected_policy"])
        for row in old_scaffold_rows
    }
    if set(old_scaffold_flags) != endpoint_keys:
        raise ValueError("Old scaffold audit keyset differs from the frozen endpoint")
    new_scaffold_flags = {
        row["canonical_pair_key"]: bool(row["audit_scaffold_cold_under_selected_policy"])
        for row in scaffold_rows
    }
    scaffold_flag_changes = sum(old_scaffold_flags[key] != new_scaffold_flags[key] for key in endpoint_keys)

    candidate_fields, candidate_rows = read_table(SOURCES["candidate_targets"])
    require_fields(candidate_fields, {"uniprot_canonical_accession"}, "candidate targets")
    assert_unique(candidate_rows, ("uniprot_canonical_accession",), "candidate targets")
    candidate_ids = {row["uniprot_canonical_accession"] for row in candidate_rows}
    if len(candidate_ids) != 4123 or not endpoint_targets.issubset(candidate_ids):
        raise ValueError("Endpoint targets are not an exact subset of the fixed 4,123-target universe")
    historical_target_ids = {row["uniprot_canonical_accession"] for row in corrected_history}
    old_historical_target_ids = {row["uniprot_canonical_accession"] for row in old_history_rows}
    if historical_target_ids != old_historical_target_ids or len(historical_target_ids) != 1131:
        raise ValueError("Historical target membership changed; homology masks may not be reused")

    homology_flags: dict[str, dict[str, bool]] = {}
    for threshold in ("0_30", "0_50", "0_70"):
        fields, rows = read_table(SOURCES[f"homology_{threshold}"])
        require_fields(
            fields,
            {"uniprot_canonical_accession", "is_future_target_homology_cold_candidate"},
            f"homology {threshold}",
        )
        assert_unique(rows, ("uniprot_canonical_accession",), f"homology {threshold}")
        flags = {
            row["uniprot_canonical_accession"]: parse_bool(row["is_future_target_homology_cold_candidate"])
            for row in rows
        }
        if set(flags) != endpoint_targets:
            raise ValueError(f"Homology {threshold} target keyset differs from frozen endpoint targets")
        homology_flags[threshold] = flags

    scope_summary: dict[str, dict[str, int]] = {}
    scope_memberships: dict[str, list[dict[str, str]]] = {
        "temporal_strict_ab": endpoint_rows,
        "scaffold_cold_strict_ab": [
            row for row in endpoint_rows if new_scaffold_flags[row["canonical_pair_key"]]
        ],
    }
    for threshold in ("0_30", "0_50", "0_70"):
        scope_memberships[f"double_cold_{threshold}"] = [
            row
            for row in endpoint_rows
            if new_scaffold_flags[row["canonical_pair_key"]]
            and homology_flags[threshold][row["uniprot_canonical_accession"]]
        ]
    for scope, rows in scope_memberships.items():
        scope_summary[scope] = {
            "relations": len(rows),
            "queries": len({row["query_id"] for row in rows}),
            "targets": len({row["uniprot_canonical_accession"] for row in rows}),
        }

    output_root.mkdir(parents=False, exist_ok=False)
    restricted_dir = output_root / "restricted_ledger"
    scoring_dir = output_root / "scoring_inputs"
    evaluation_dir = output_root / "evaluation_inputs"
    audit_dir = output_root / "audit"
    metadata_dir = output_root / "metadata"
    for directory in (restricted_dir, scoring_dir, evaluation_dir, audit_dir, metadata_dir):
        directory.mkdir()

    ledger_path = restricted_dir / "historical_row_eligibility.tsv.gz"
    pair_diff_path = restricted_dir / "historical_pair_before_after.tsv.gz"
    structure_audit_path = restricted_dir / "role_separated_compound_structure_audit.tsv.gz"
    repair_audit_path = restricted_dir / "historical_v2_inchi_repair_audit.tsv.gz"
    write_tsv_gz(ledger_path, list(row_ledger[0]), row_ledger)
    write_tsv_gz(pair_diff_path, list(pair_diff_rows[0]), pair_diff_rows)
    write_tsv_gz(structure_audit_path, list(structure_audit_rows[0]), structure_audit_rows)
    repair_fields = [
        "inchikey_full",
        "source_version",
        "source_np_id",
        "source_row_number",
        "eligible_source_np_id",
        "validation_status",
        "repaired_smiles",
    ]
    write_tsv_gz(repair_audit_path, repair_fields, repair_audit_rows)

    scoring_paths = {
        "historical_pairs": scoring_dir / "historical_pairs.tsv.gz",
        "scoring_queries": scoring_dir / "scoring_queries.tsv.gz",
        "historical_compounds": scoring_dir / "historical_compounds.tsv.gz",
        "query_compounds": scoring_dir / "query_compounds.tsv.gz",
        "candidate_targets": scoring_dir / "candidate_targets.tsv.gz",
        "candidate_sequences": scoring_dir / "candidate_sequences.fasta",
    }
    write_tsv_gz(scoring_paths["historical_pairs"], list(corrected_history[0]), corrected_history)
    copy_exact(SOURCES["old_queries"], scoring_paths["scoring_queries"])
    write_tsv_gz(scoring_paths["historical_compounds"], list(historical_compound_rows[0]), historical_compound_rows)
    write_tsv_gz(scoring_paths["query_compounds"], list(query_compound_rows[0]), query_compound_rows)
    copy_exact(SOURCES["candidate_targets"], scoring_paths["candidate_targets"])
    copy_exact(SOURCES["candidate_sequences"], scoring_paths["candidate_sequences"])

    evaluation_paths = {
        "endpoint": evaluation_dir / "evaluation_pairs.tsv.gz",
        "scaffold_audit": evaluation_dir / "scaffold_audit.tsv.gz",
        "homology_0_30": evaluation_dir / "homology_0_30.tsv.gz",
        "homology_0_50": evaluation_dir / "homology_0_50.tsv.gz",
        "homology_0_70": evaluation_dir / "homology_0_70.tsv.gz",
    }
    copy_exact(SOURCES["endpoint"], evaluation_paths["endpoint"])
    write_tsv_gz(evaluation_paths["scaffold_audit"], list(scaffold_rows[0]), scaffold_rows)
    for threshold in ("0_30", "0_50", "0_70"):
        copy_exact(SOURCES[f"homology_{threshold}"], evaluation_paths[f"homology_{threshold}"])

    created_at = datetime.now(timezone.utc).isoformat()
    authorization = (
        "Project-lead local repair authorization recorded in "
        "strict_ab_asof_cutoff_successor_v1_20260728/governance/project_lead_authorization_20260728.md; "
        "no public redistribution or independent-evaluation claim."
    )
    scoring_manifest = {
        "protocol_id": PROTOCOL_ID,
        "run_id": RUN_ID,
        "execution_mode": "author_run_non_independent_corrective_successor",
        "input_kind": "corrective_role_separated_scoring_without_endpoint_file",
        "project_lead_authorized_internal_use": True,
        "legacy_outer_or_result_input": False,
        "endpoint_file_included": False,
        "access_level": "restricted_author_run_internal",
        "authorization_basis": authorization,
        "file_sha256": {name: sha256(path) for name, path in scoring_paths.items()},
    }
    evaluation_manifest = {
        "protocol_id": PROTOCOL_ID,
        "run_id": RUN_ID,
        "execution_mode": "author_run_non_independent_corrective_successor",
        "input_kind": "corrective_evaluation_endpoint",
        "project_lead_authorized_internal_use": True,
        "legacy_outer_or_result_input": False,
        "endpoint_file_included": True,
        "access_level": "restricted_author_run_internal",
        "authorization_basis": authorization,
        "file_sha256": {name: sha256(path) for name, path in evaluation_paths.items()},
    }
    scoring_manifest_path = scoring_dir / "author_run_input_manifest.json"
    evaluation_manifest_path = evaluation_dir / "author_run_input_manifest.json"
    write_json(scoring_manifest_path, scoring_manifest)
    write_json(evaluation_manifest_path, evaluation_manifest)

    receipt = {
        "protocol_id": PROTOCOL_ID,
        "run_id": RUN_ID,
        "created_at_utc": created_at,
        "execution_mode": "author_run_non_independent_corrective_successor",
        "project_lead_authorized_internal_use": True,
        "public_release_authorized": False,
        "protocol_lock": {"path": str(PROTOCOL_LOCK), "sha256": sha256(PROTOCOL_LOCK)},
        "code_lock": {"path": str(args.code_lock_manifest.resolve()), "sha256": sha256(args.code_lock_manifest)},
        "spec": {"path": str(SPEC), "sha256": sha256(SPEC)},
        "endpoint_handling": "Frozen endpoint was used to define query IDs and evaluation only; the scorer receives no endpoint file.",
        "claim_boundary": "Author-run post hoc correction; no blind, independent, external-validation, direct-binding, or public-release claim.",
    }
    receipt_path = metadata_dir / "author_run_protocol_receipt.json"
    write_json(receipt_path, receipt)

    audit_summary = {
        "protocol_id": PROTOCOL_ID,
        "run_id": RUN_ID,
        "created_at_utc": created_at,
        "rdkit_version": rdBase.rdkitVersion,
        "cutoff": CUTOFF.isoformat(),
        "row_eligibility_counts": normalized_row_status_counts,
        "counts": {
            "historical_strict_v2_rows": len(row_ledger),
            "historical_pairs": len(corrected_history),
            "historical_pairs_with_any_excluded_row": sum(value > 0 for value in excluded_rows_by_pair.values()),
            "historical_tier_weight_changes": tier_change_count,
            "historical_compounds": len(historical_compound_rows),
            "query_compounds": len(query_compound_rows),
            "shared_role_compounds": len(history_compound_ids.intersection(query_compound_ids)),
            "historical_inchi_repairs": len(repaired_history),
            "historical_smiles_changed_vs_old_combined": sum(
                row["structure_role"] == "historical" and row["representative_smiles_changed"]
                for row in structure_audit_rows
            ),
            "query_smiles_changed_vs_old_combined": sum(
                row["structure_role"] == "query" and row["representative_smiles_changed"]
                for row in structure_audit_rows
            ),
            "historical_morgan_changed_vs_old_combined": sum(
                row["structure_role"] == "historical" and not row["morgan_2048_radius2_equal_to_old"]
                for row in structure_audit_rows
            ),
            "query_morgan_changed_vs_old_combined": sum(
                row["structure_role"] == "query" and not row["morgan_2048_radius2_equal_to_old"]
                for row in structure_audit_rows
            ),
            "historical_scaffold_assignment_changes": sum(
                row["structure_role"] == "historical" and not row["scaffold_assignment_equal_to_old"]
                for row in structure_audit_rows
            ),
            "query_scaffold_assignment_changes": sum(
                row["structure_role"] == "query" and not row["scaffold_assignment_equal_to_old"]
                for row in structure_audit_rows
            ),
            "endpoint_scaffold_flag_changes": scaffold_flag_changes,
            "endpoint_relations": len(endpoint_rows),
            "endpoint_queries": len(query_ids),
            "endpoint_targets": len(endpoint_targets),
            "candidate_targets": len(candidate_ids),
            "historical_targets": len(historical_target_ids),
            "endpoint_v3_strict_source_rows": len(endpoint_strict_rows),
        },
        "scope_counts": scope_summary,
        "membership_sha256": {
            "historical_pairs": membership_sha256(historical_keys),
            "endpoint_pairs": membership_sha256(endpoint_keys),
            "historical_targets": membership_sha256(historical_target_ids),
            "candidate_targets": membership_sha256(candidate_ids),
        },
        "hard_gate_checks": {
            "protocol_lock_verified": True,
            "all_source_hashes_verified": True,
            "history_exactly_4990": True,
            "every_history_pair_has_eligible_v2_row": True,
            "no_v3_history_contribution": True,
            "historical_query_structure_maps_separate": True,
            "endpoint_hash_unchanged": sha256(evaluation_paths["endpoint"]) == EXPECTED_HASHES["endpoint"],
            "candidate_universe_4123": len(candidate_ids) == 4123,
            "historical_target_membership_unchanged": historical_target_ids == old_historical_target_ids,
            "homology_masks_reused_only_after_target_keyset_check": True,
        },
        "claim_boundary": "Restricted aggregate audit of an author-run corrective successor; no external-validation or biological claim.",
    }
    audit_summary_path = audit_dir / "asof_rebuild_summary.json"
    write_json(audit_summary_path, audit_summary)

    provenance = {
        "protocol_id": PROTOCOL_ID,
        "run_id": RUN_ID,
        "created_at_utc": created_at,
        "access_level": "restricted_author_run_internal",
        "source_files": {
            name: {"path": str(path), "sha256": sha256(path)} for name, path in SOURCES.items()
        },
        "restricted_outputs": {
            "historical_row_eligibility": {"path": str(ledger_path), "sha256": sha256(ledger_path)},
            "historical_pair_before_after": {"path": str(pair_diff_path), "sha256": sha256(pair_diff_path)},
            "role_separated_structure_audit": {"path": str(structure_audit_path), "sha256": sha256(structure_audit_path)},
            "historical_v2_inchi_repair_audit": {"path": str(repair_audit_path), "sha256": sha256(repair_audit_path)},
        },
        "scoring_outputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in scoring_paths.items()},
        "evaluation_outputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in evaluation_paths.items()},
        "manifests": {
            "scoring": {"path": str(scoring_manifest_path), "sha256": sha256(scoring_manifest_path)},
            "evaluation": {"path": str(evaluation_manifest_path), "sha256": sha256(evaluation_manifest_path)},
            "receipt": {"path": str(receipt_path), "sha256": sha256(receipt_path)},
            "audit_summary": {"path": str(audit_summary_path), "sha256": sha256(audit_summary_path)},
        },
        "counts": audit_summary["counts"],
        "scope_counts": scope_summary,
        "release_status": "restricted_internal_no_public_release",
    }
    provenance_path = metadata_dir / "bundle_provenance.json"
    write_json(provenance_path, provenance)

    print(
        json.dumps(
            {
                "output_root": str(output_root),
                "row_eligibility_counts": normalized_row_status_counts,
                "historical_tier_changes": tier_change_count,
                "scaffold_flag_changes": scaffold_flag_changes,
                "scope_counts": scope_summary,
                "figures_generated": False,
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

