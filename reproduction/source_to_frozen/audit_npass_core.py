#!/usr/bin/env python3
"""Create a source-faithful first-pass NPASS entity-alignment audit.

This script deliberately does not call any activity a verified direct binding
event. It only identifies quantitative records whose structure, human
single-protein target, UniProt accession, and source reference fields can be
traced to the downloaded NPASS files.
"""

from __future__ import annotations

import argparse
import json
import re
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from reproducible_io import PANDAS_GZIP


SOURCES = {
    "v2": {
        "raw": "data/raw/npass2",
        "structure": "NPASSv2.0_download_naturalProducts_structureInfo.txt",
        "target": "NPASSv2.0_download_naturalProducts_targetInfo.txt",
        "activity": "NPASSv2.0_download_naturalProducts_activities.txt",
    },
    "v3": {
        "raw": "data/raw/npass3",
        "structure": "NPASS3.0_naturalproducts_structure.txt",
        "target": "NPASS3.0_target.txt",
        "activity": "NPASS3.0_activities.txt",
    },
}

STRUCTURE_COLUMNS = ["np_id", "InChI", "InChIKey", "SMILES"]
TARGET_COLUMNS = ["target_id", "target_type", "target_name", "target_organism_tax_id", "target_organism", "uniprot_id"]
ACTIVITY_COLUMNS = [
    "np_id", "target_id", "activity_type_grouped", "activity_relation", "activity_type", "activity_value",
    "activity_units", "assay_organism", "assay_tax_id", "assay_strain", "assay_tissue", "assay_cell_type",
    "ref_id", "ref_id_type",
]
CANONICAL_COLUMNS = [
    "source_version", "source_np_id", "source_target_id", "inchikey_full", "inchikey_connectivity", "smiles",
    "uniprot_raw", "target_type", "target_name", "target_tax_id", "activity_type_grouped", "activity_relation",
    "activity_type", "activity_value", "activity_units", "assay_organism", "assay_tax_id", "assay_strain",
    "assay_tissue", "assay_cell_type", "ref_id", "ref_id_type", "pair_key",
]


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def read_table(path: Path, expected: list[str]) -> pd.DataFrame:
    table = pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, low_memory=False)
    missing = sorted(set(expected).difference(table.columns))
    if missing:
        raise ValueError(f"{path.name} is missing expected columns: {missing}")
    return table


def normalize_uniprot(value: str) -> str:
    """Retain the first accession-like token but preserve raw text separately."""
    value = str(value or "").strip()
    if not value or value.lower() in {"n.a.", "na", "nan"}:
        return ""
    return re.split(r"[;,|\s]+", value, maxsplit=1)[0].upper()


def assert_files(root: Path, source: dict[str, str]) -> dict[str, Path]:
    raw = root / source["raw"]
    paths = {key: raw / source[key] for key in ("structure", "target", "activity")}
    missing = [str(path) for path in paths.values() if not path.exists()]
    if missing:
        raise FileNotFoundError("Missing raw input(s):\n" + "\n".join(missing))
    return paths


def prepare_lookup_tables(paths: dict[str, Path]) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    structures = read_table(paths["structure"], STRUCTURE_COLUMNS)[STRUCTURE_COLUMNS].copy()
    structures["InChIKey"] = structures["InChIKey"].str.strip().str.upper()
    structures["SMILES"] = structures["SMILES"].str.strip()
    structures = structures.drop_duplicates(subset="np_id", keep="first")

    targets = read_table(paths["target"], TARGET_COLUMNS)[TARGET_COLUMNS].copy()
    targets["uniprot_norm"] = targets["uniprot_id"].map(normalize_uniprot)
    targets["target_tax_norm"] = targets["target_organism_tax_id"].str.strip().str.replace(".0", "", regex=False)
    type_norm = targets["target_type"].str.strip().str.casefold()
    human_single = (targets["target_tax_norm"] == "9606") & type_norm.isin({"individual protein", "single protein"})
    eligible_targets = targets.loc[human_single & targets["uniprot_norm"].ne("")].drop_duplicates(subset="target_id", keep="first")
    metrics = {
        "structures": int(len(structures)),
        "structures_with_inchikey": int(structures["InChIKey"].ne("").sum()),
        "targets": int(len(targets)),
        "human_single_protein_targets": int(human_single.sum()),
        "human_single_protein_targets_with_uniprot": int(len(eligible_targets)),
    }
    return structures, eligible_targets, metrics


def write_candidates(version: str, paths: dict[str, Path], interim_dir: Path) -> dict[str, object]:
    structures, targets, metrics = prepare_lookup_tables(paths)
    output = interim_dir / f"npass_{version}_human_single_protein_records.tsv.gz"
    if output.exists():
        output.unlink()

    total_records = 0
    candidate_records = 0
    pmid_records = 0
    first = True
    for activity in pd.read_csv(paths["activity"], sep="\t", dtype=str, keep_default_na=False, chunksize=200_000, low_memory=False):
        missing = sorted(set(ACTIVITY_COLUMNS).difference(activity.columns))
        if missing:
            raise ValueError(f"{paths['activity'].name} is missing expected columns: {missing}")
        total_records += len(activity)
        joined = activity.merge(targets, on="target_id", how="inner", validate="many_to_one")
        joined = joined.merge(structures[["np_id", "InChIKey", "SMILES"]], on="np_id", how="left", validate="many_to_one")
        joined["InChIKey"] = joined["InChIKey"].fillna("").str.strip().str.upper()
        joined["SMILES"] = joined["SMILES"].fillna("").str.strip()
        joined = joined.loc[joined["InChIKey"].ne("")].copy()
        if joined.empty:
            continue
        joined["source_version"] = version
        joined["source_np_id"] = joined["np_id"]
        joined["source_target_id"] = joined["target_id"]
        joined["inchikey_full"] = joined["InChIKey"]
        joined["inchikey_connectivity"] = joined["InChIKey"].str.split("-", n=1).str[0]
        joined["smiles"] = joined["SMILES"]
        joined["uniprot_raw"] = joined["uniprot_norm"]
        joined["target_tax_id"] = joined["target_tax_norm"]
        joined["pair_key"] = joined["inchikey_full"] + "|" + joined["uniprot_raw"]
        candidate = joined[CANONICAL_COLUMNS]
        candidate_records += len(candidate)
        pmid_records += int(candidate["ref_id_type"].str.strip().str.upper().eq("PMID").sum())
        candidate.to_csv(output, sep="\t", index=False, mode="wt" if first else "at", header=first, compression=PANDAS_GZIP)
        first = False

    if first:
        pd.DataFrame(columns=CANONICAL_COLUMNS).to_csv(output, sep="\t", index=False, compression=PANDAS_GZIP)
    metrics.update({
        "activity_records": int(total_records),
        "candidate_human_single_protein_records": int(candidate_records),
        "candidate_records_with_pmid": int(pmid_records),
        "candidate_file": str(output),
    })
    return metrics


def make_cross_version_pair_audit(interim_dir: Path, results_dir: Path) -> dict[str, int]:
    frames = []
    for version in ("v2", "v3"):
        path = interim_dir / f"npass_{version}_human_single_protein_records.tsv.gz"
        if path.exists():
            frame = pd.read_csv(path, sep="\t", usecols=["source_version", "pair_key"], dtype=str, keep_default_na=False)
            frames.append(frame)
    if len(frames) != 2:
        return {}
    pairs = pd.concat(frames, ignore_index=True)
    counts = pairs.groupby(["pair_key", "source_version"], as_index=False).size().pivot(index="pair_key", columns="source_version", values="size").fillna(0).astype(int).reset_index()
    counts["in_v2"] = counts.get("v2", 0).gt(0)
    counts["in_v3"] = counts.get("v3", 0).gt(0)
    counts["cross_version_status"] = "v3_only"
    counts.loc[counts["in_v2"] & counts["in_v3"], "cross_version_status"] = "shared"
    counts.loc[counts["in_v2"] & ~counts["in_v3"], "cross_version_status"] = "v2_only"
    output = results_dir / "npass_cross_version_pair_membership.csv.gz"
    counts.to_csv(output, index=False, compression=PANDAS_GZIP)
    return {
        "unique_pairs": int(len(counts)),
        "v2_only_pairs": int((counts["cross_version_status"] == "v2_only").sum()),
        "v3_only_pairs": int((counts["cross_version_status"] == "v3_only").sum()),
        "shared_pairs": int((counts["cross_version_status"] == "shared").sum()),
        "pair_membership_file": str(output),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--versions", nargs="+", choices=("v2", "v3"), default=("v2", "v3"))
    args = parser.parse_args()
    root = args.project_root.resolve()
    interim_dir = root / "data/interim"
    results_dir = root / "results"
    interim_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    summary: dict[str, object] = {"created_at": utc_now(), "scope": "entity alignment only", "versions": {}}
    for version in args.versions:
        summary["versions"][version] = write_candidates(version, assert_files(root, SOURCES[version]), interim_dir)
    if set(args.versions) == {"v2", "v3"}:
        summary["cross_version_pair_membership"] = make_cross_version_pair_audit(interim_dir, results_dir)
    output = results_dir / "initial_entity_alignment_audit.json"
    output.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
