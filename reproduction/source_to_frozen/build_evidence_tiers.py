#!/usr/bin/env python3
"""Apply the frozen v1 evidence policy to aligned NPASS records.

The output is a conservative, record-level screening table. It does not replace
assay metadata review or source-paper verification.
"""

from __future__ import annotations

import argparse
import json
import math
import re
from pathlib import Path

import pandas as pd

from reproducible_io import PANDAS_GZIP


MOLAR_FACTORS = {
    "M": 1.0, "MM": 1e-3, "UM": 1e-6, "NM": 1e-9, "PM": 1e-12, "FM": 1e-15,
}


def normalise_unit(value: str) -> str:
    value = str(value or "").strip().upper().replace("μ", "U").replace("µ", "U")
    value = value.replace("MOLAR", "M")
    return value


def numeric_value(value: str) -> float | None:
    try:
        parsed = float(str(value).strip())
        return parsed if math.isfinite(parsed) else None
    except (TypeError, ValueError):
        return None


def classify(record: pd.Series, policy: dict) -> tuple[str, float, str, float | None, float | None]:
    relation = str(record["activity_relation"] or "").strip()
    endpoint = str(record["activity_type"] or "").strip().upper()
    value = numeric_value(record["activity_value"])
    unit = normalise_unit(record["activity_units"])
    molar = value * MOLAR_FACTORS[unit] if value is not None and unit in MOLAR_FACTORS else None
    p_activity = -math.log10(molar) if molar is not None and molar > 0 else None

    if relation != "=":
        return "D_censored_or_unresolved", 0.0, "non_exact_relation", molar, p_activity
    if value is None:
        return "D_censored_or_unresolved", 0.0, "missing_or_non_numeric_value", molar, p_activity
    affinity_tier = "A_affinity_candidate"
    potency_tier = "B_quantitative_functional_candidate" if "B_quantitative_functional_candidate" in policy["tiers"] else "B_quantitative_potency_candidate"
    contextual_tier = "C_contextual_functional_exact" if "C_contextual_functional_exact" in policy["tiers"] else "C_functional_or_other_exact"
    quarantine_tier = "Q_unverified_screen_or_other_exact" if "Q_unverified_screen_or_other_exact" in policy["tiers"] else contextual_tier
    potency_endpoints = policy.get("functional_potency_endpoints", policy.get("potency_endpoints", []))
    contextual_endpoints = policy.get("contextual_functional_endpoints", policy.get("functional_endpoints", []))
    quarantine_endpoints = policy.get("quarantine_endpoints", [])
    if endpoint in set(policy["affinity_endpoints"]):
        tier = affinity_tier if molar is not None else quarantine_tier
    elif endpoint in set(potency_endpoints):
        tier = potency_tier if molar is not None else quarantine_tier
    elif endpoint in set(contextual_endpoints):
        tier = contextual_tier
    elif endpoint in set(quarantine_endpoints):
        tier = quarantine_tier
    else:
        tier = quarantine_tier
    return tier, float(policy["tiers"][tier]["weight"]), "exact_record", molar, p_activity


def measurement_class(record: pd.Series) -> str:
    endpoint = str(record["activity_type"] or "").strip().upper()
    if endpoint in {"KI", "KD"}:
        return "affinity_or_inhibition_constant_candidate"
    if endpoint in {"IC50", "EC50"}:
        return "quantitative_functional_potency_candidate"
    if endpoint in {"AC50", "POTENCY"}:
        return "screen_or_potency_endpoint_needs_assay_context"
    if endpoint in {"INHIBITION", "ACTIVITY", "EMAX", "EFFICACY", "%INHIB (MEAN)", "%MAX (MEAN)", "RESIDUAL ACTIVITY"}:
        return "functional_effect_needs_assay_context"
    return "other_or_unclassified_endpoint"


def verification_fields(tier: str) -> tuple[str, bool, str, bool]:
    if tier in {"A_affinity_candidate", "B_quantitative_functional_candidate", "C_contextual_functional_exact", "B_quantitative_potency_candidate", "C_functional_or_other_exact"}:
        return "P1_npass_raw_exact_candidate", False, "weak_candidate_only", False
    return "P0_censored_or_unresolved", False, "exclude_from_positive_and_unlabeled", False


def load_temporal_maps(root: Path) -> tuple[dict[str, str], dict[str, str]]:
    membership = pd.read_csv(root / "results/npass_cross_version_pair_membership.csv.gz", dtype=str, keep_default_na=False)
    temporal = pd.read_csv(root / "results/temporal_v3_only_pmid_screen.csv.gz", dtype=str, keep_default_na=False)
    return (
        dict(zip(membership["pair_key"], membership["cross_version_status"])),
        dict(zip(temporal["pair_key"], temporal["temporal_screen_status"])),
    )


def apply_version(root: Path, version: str, policy: dict, policy_tag: str, membership: dict[str, str], temporal: dict[str, str]) -> pd.DataFrame:
    input_path = root / "data/interim" / f"npass_{version}_human_single_protein_records.tsv.gz"
    records = pd.read_csv(input_path, sep="\t", dtype=str, keep_default_na=False)
    classified = records.apply(lambda row: classify(row, policy), axis=1, result_type="expand")
    classified.columns = [f"evidence_tier_{policy_tag}", f"evidence_weight_{policy_tag}", "evidence_screen_reason", "activity_value_molar", "p_activity"]
    records = pd.concat([records, classified], axis=1)
    records["measurement_class"] = records.apply(measurement_class, axis=1)
    verification = records[f"evidence_tier_{policy_tag}"].map(verification_fields)
    records[["automatic_verification_level", "training_eligible_before_assay_review", "pu_role_before_assay_review", "direct_binding_permitted"]] = pd.DataFrame(verification.tolist(), index=records.index)
    records["cross_version_status"] = records["pair_key"].map(membership).fillna("not_in_membership_audit")
    records["temporal_screen_status"] = records["pair_key"].map(temporal).fillna("not_v3_only_or_not_screened")
    output_path = root / "data/processed" / f"npass_{version}_evidence_records_{policy_tag}.tsv.gz"
    records.to_csv(output_path, sep="\t", index=False, compression=PANDAS_GZIP)
    return records


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--policy", type=Path, default=None, help="policy JSON; defaults to evidence_policy_v1.json")
    parser.add_argument("--tag", default="v1", help="short tag used in output filenames and columns")
    args = parser.parse_args()
    root = args.project_root.resolve()
    policy_path = args.policy.resolve() if args.policy else root / "configs/evidence_policy_v1.json"
    policy = json.loads(policy_path.read_text(encoding="utf-8"))
    membership, temporal = load_temporal_maps(root)
    summaries = []
    for version in ("v2", "v3"):
        records = apply_version(root, version, policy, args.tag, membership, temporal)
        tier_column = f"evidence_tier_{args.tag}"
        grouped = records.groupby(["source_version", tier_column], dropna=False).agg(
            records=("pair_key", "size"), pairs=("pair_key", "nunique"), compounds=("inchikey_full", "nunique"), targets=("uniprot_raw", "nunique")
        ).reset_index()
        summaries.append(grouped)
        if version == "v3":
            future = records.loc[records["temporal_screen_status"].eq("future_candidate_pmid_only")]
            summaries.append(future.groupby(["source_version", tier_column], dropna=False).agg(
                records=("pair_key", "size"), pairs=("pair_key", "nunique"), compounds=("inchikey_full", "nunique"), targets=("uniprot_raw", "nunique")
            ).reset_index().assign(scope="v3_only_future_pmid_screen"))
    summary = pd.concat(summaries, ignore_index=True)
    summary["scope"] = summary.get("scope", pd.Series(index=summary.index, dtype=str)).fillna("all_aligned_records")
    summary_path = root / "results" / f"evidence_tier_{args.tag}_summary.csv"
    summary.to_csv(summary_path, index=False)
    (root / "results" / f"evidence_tier_{args.tag}_metadata.json").write_text(json.dumps({
        "policy": str(policy_path),
        "summary": str(summary_path),
        "warning": "Tier A is an affinity candidate, not an asserted direct-binding label."
    }, indent=2) + "\n", encoding="utf-8")
    print(summary_path)


if __name__ == "__main__":
    main()
