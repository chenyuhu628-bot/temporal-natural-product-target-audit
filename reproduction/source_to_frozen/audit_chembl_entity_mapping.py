#!/usr/bin/env python3
"""Audit exact NPASS future-candidate entity mapping against frozen ChEMBL releases.

This is only a molecule/target identity audit. It does not infer ChEMBL
activities and deliberately avoids connectivity-only structure matches.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from reproducible_io import PANDAS_GZIP


TIER_RANK = {
    "A_affinity_candidate": 3,
    "B_quantitative_functional_candidate": 2,
    "C_contextual_functional_exact": 1,
}


def source_paths(root: Path, version: str) -> tuple[Path, Path]:
    base = root / "data/raw/chembl" / f"chembl_{version}"
    chemreps = base / f"chembl_{version}_chemreps.txt.gz"
    mapping = base / "chembl_uniprot_mapping.txt"
    if not chemreps.exists() or not mapping.exists():
        raise FileNotFoundError(f"Missing minimal ChEMBL {version} files under {base}")
    return chemreps, mapping


def future_candidates(root: Path) -> pd.DataFrame:
    source = root / "data/processed/npass_v3_evidence_records_v1_1_uniprot_mapped.tsv.gz"
    records = pd.read_csv(source, sep="\t", dtype=str, keep_default_na=False)
    records = records.loc[
        records["temporal_screen_status"].eq("future_candidate_pmid_only")
        & records["automatic_verification_level"].eq("P1_npass_raw_exact_candidate")
        & records["mapping_status"].eq("strict_one_to_one_reviewed_human")
    ].copy()
    records["tier_rank"] = records["evidence_tier_v1_1"].map(TIER_RANK).fillna(0).astype(int)
    rows = []
    for pair_key, group in records.groupby("pair_key", sort=False):
        representative = group.sort_values("tier_rank", ascending=False).iloc[0]
        rows.append({
            "pair_key": pair_key,
            "inchikey_full": representative["inchikey_full"].upper(),
            "npass_uniprot_source": representative["uniprot_raw"],
            "npass_uniprot_canonical": representative["uniprot_canonical_accession"],
            "best_evidence_tier_v1_1": representative["evidence_tier_v1_1"],
            "npass_record_count": len(group),
        })
    return pd.DataFrame(rows)


def read_chemrep_matches(chemreps: Path, candidate_keys: set[str]) -> dict[str, list[str]]:
    matches: dict[str, list[str]] = {}
    for chunk in pd.read_csv(chemreps, sep="\t", compression="gzip", dtype=str, keep_default_na=False, chunksize=250_000, usecols=["chembl_id", "standard_inchi_key"]):
        chunk["standard_inchi_key"] = chunk["standard_inchi_key"].str.strip().str.upper()
        hit = chunk.loc[chunk["standard_inchi_key"].isin(candidate_keys)]
        for key, group in hit.groupby("standard_inchi_key"):
            matches.setdefault(key, []).extend(group["chembl_id"].tolist())
    return {key: sorted(set(values)) for key, values in matches.items()}


def read_target_mapping(mapping: Path) -> dict[str, list[str]]:
    table = pd.read_csv(mapping, sep="\t", comment="#", header=None, names=["uniprot", "chembl_target_id", "target_name", "target_type"], dtype=str, keep_default_na=False)
    table["uniprot"] = table["uniprot"].str.strip().str.upper()
    table = table.loc[table["target_type"].str.strip().str.upper().eq("SINGLE PROTEIN")]
    return table.groupby("uniprot")["chembl_target_id"].agg(lambda values: sorted(set(values))).to_dict()


def audit_release(root: Path, version: str, candidates: pd.DataFrame) -> tuple[pd.DataFrame, dict[str, object]]:
    chemreps, mapping_path = source_paths(root, version)
    compound_map = read_chemrep_matches(chemreps, set(candidates["inchikey_full"]))
    target_map = read_target_mapping(mapping_path)
    output = candidates.copy()
    output["chembl_release"] = version
    output["chembl_compound_ids"] = output["inchikey_full"].map(lambda key: ";".join(compound_map.get(key, [])))
    output["chembl_compound_match_status"] = output["chembl_compound_ids"].ne("").map({True: "full_inchikey_exact", False: "unmatched"})
    output["chembl_target_ids"] = output["npass_uniprot_source"].str.upper().map(lambda key: ";".join(target_map.get(key, [])))
    output["chembl_target_match_status"] = output["chembl_target_ids"].ne("").map({True: "source_uniprot_exact", False: "unmatched"})
    output["both_entities_exactly_mapped"] = output["chembl_compound_ids"].ne("") & output["chembl_target_ids"].ne("")
    output_path = root / "data/interim" / f"chembl_{version}_future_candidate_entity_mapping.csv.gz"
    output.to_csv(output_path, index=False, compression=PANDAS_GZIP)
    summary = {
        "chembl_release": version,
        "candidate_pairs": int(len(output)),
        "compound_full_inchikey_matched_pairs": int(output["chembl_compound_ids"].ne("").sum()),
        "target_source_uniprot_matched_pairs": int(output["chembl_target_ids"].ne("").sum()),
        "both_entities_exactly_mapped_pairs": int(output["both_entities_exactly_mapped"].sum()),
        "mapping_file": str(output_path),
        "warning": "Both-entity mapping is not evidence that ChEMBL records this pair; SQLite activity extraction is required."
    }
    return output, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--releases", nargs="+", default=["31"], choices=["31", "37"])
    args = parser.parse_args()
    root = args.project_root.resolve()
    candidates = future_candidates(root)
    summaries = []
    for release in args.releases:
        _, summary = audit_release(root, release, candidates)
        summaries.append(summary)
    result = {"scope": "strict UniProt-mapped v3 future P1 candidates; full InChIKey and source-UniProt only", "candidate_pair_count": int(len(candidates)), "releases": summaries}
    output = root / "results/chembl_future_candidate_entity_mapping_summary.json"
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
