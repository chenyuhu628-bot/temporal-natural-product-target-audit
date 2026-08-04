"""Validate aggregate-only unresolved endpoint bounds."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "outputs"
EXPECTED_ROWS = {
    "mapping_failure_strata.tsv": 2,
    "cohort_comparability.tsv": 3,
    "comparability_differences.tsv": 8,
    "endpoint_cardinality_bounds.tsv": 1,
    "relation_level_top50_bounds.tsv": 4,
    "non_identifiable_estimands.tsv": 3,
}
EXPECTED_HITS = {
    "weighted_target_popularity": 54,
    "sequence_3mer_transfer": 9,
    "weighted_morgan_transfer": 78,
    "structure_sequence_pair_neighbor": 81,
}
FORBIDDEN_HEADERS = {
    "canonical_pair_key",
    "pair_key",
    "inchikey_full",
    "uniprot_canonical_accession",
    "query_id",
    "target_uniprot_accession",
    "smiles",
    "ref_id",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def read_tsv(name: str) -> list[dict[str, str]]:
    with (OUTPUT / name).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(bool(reader.fieldnames), f"Missing header: {name}")
        require(not FORBIDDEN_HEADERS.intersection(reader.fieldnames or []), f"Identifier-bearing output: {name}")
        return list(reader)


def main() -> int:
    checks: list[dict[str, object]] = []

    def check(name: str, condition: bool, detail: str) -> None:
        require(condition, f"{name}: {detail}")
        checks.append({"check": name, "status": "PASS", "detail": detail})

    tables = {name: read_tsv(name) for name in EXPECTED_ROWS}
    check(
        "row_counts",
        all(len(tables[name]) == expected for name, expected in EXPECTED_ROWS.items()),
        str({name: len(rows) for name, rows in tables.items()}),
    )
    failures = {row["failure_stratum"]: int(row["relation_count"]) for row in tables["mapping_failure_strata.tsv"]}
    check(
        "mapping_failure_partition",
        failures == {"preliminary_compound_unmatched": 51, "preliminary_target_unmatched": 14},
        str(failures),
    )
    cardinality = tables["endpoint_cardinality_bounds.tsv"][0]
    check(
        "endpoint_partition",
        int(cardinality["initial_candidate_relations"]) == 442
        and int(cardinality["definitive_historical_activity_exclusions"]) == 19
        and int(cardinality["frozen_primary_endpoint_relations"]) == 358
        and int(cardinality["entity_unresolved_relations"]) == 65,
        "442 = 358 + 19 + 65",
    )
    check(
        "endpoint_cardinality_bounds",
        int(cardinality["identified_endpoint_cardinality_lower"]) == 358
        and int(cardinality["identified_endpoint_cardinality_upper"]) == 423,
        "358–423",
    )
    check(
        "no_negative_or_readmission",
        int(cardinality["unresolved_negative_labels"]) == 0
        and int(cardinality["unresolved_readmissions"]) == 0,
        "zero negative labels and readmissions",
    )

    bounds = {row["baseline"]: row for row in tables["relation_level_top50_bounds.tsv"]}
    check("baseline_completeness", set(bounds) == set(EXPECTED_HITS), str(sorted(bounds)))
    for baseline, hits in EXPECTED_HITS.items():
        row = bounds[baseline]
        check(
            f"observed_hits_{baseline}",
            int(row["observed_top50_hit_relation_count"]) == hits,
            str(hits),
        )
        check(
            f"bound_formula_{baseline}",
            math.isclose(float(row["all_unresolved_fail_lower_fraction"]), hits / 423, abs_tol=1e-15)
            and math.isclose(float(row["all_unresolved_succeed_upper_fraction"]), (hits + 65) / 423, abs_tol=1e-15),
            f"[{hits}/423, ({hits}+65)/423]",
        )
        check(
            f"bound_order_{baseline}",
            0.0
            <= float(row["all_unresolved_fail_lower_fraction"])
            <= float(row["observed_primary_relation_hit_fraction"])
            <= float(row["all_unresolved_succeed_upper_fraction"])
            <= 1.0,
            "lower ≤ observed ≤ upper",
        )

    nonidentified = tables["non_identifiable_estimands.tsv"]
    check(
        "nonidentified_status",
        all(row["status"] == "not_identifiable" and row["invented_values_used"] == "false" for row in nonidentified),
        "three estimand families withheld without imputation",
    )
    check(
        "no_query_macro_bound_output",
        all(row["estimand"] == "relation_weighted_temporal_top50_hit_fraction" for row in bounds.values()),
        "relation-level temporal bounds only",
    )

    receipt = json.loads((OUTPUT / "EXECUTION_RECEIPT.json").read_text(encoding="utf-8"))
    check("execution_receipt", receipt.get("status") == "PASS", "runner PASS")
    check(
        "rank_cells",
        int(receipt.get("endpoint_rank_cells_extracted", 0)) == 1432,
        "358 relations × 4 baselines",
    )
    check(
        "receipt_claim_boundary",
        receipt.get("query_macro_bounds_computed") == 0
        and receipt.get("scope_specific_bounds_computed") == 0
        and receipt.get("identifier_bearing_output") is False
        and receipt.get("absolute_paths_emitted") is False,
        "no nonidentified bounds or identifying output",
    )

    manifest = json.loads((OUTPUT / "MANIFEST.json").read_text(encoding="utf-8"))
    check(
        "parent_protocol",
        manifest.get("parent_protocol_sha256")
        == "bd4902476160cc7c5cbacaf0cfd0f1a28c5300bde22232b3a2cc6c1f3c143dc2",
        "frozen v4 protocol hash",
    )
    check(
        "manifest_hashes",
        all((OUTPUT / name).is_file() and sha256(OUTPUT / name) == value for name, value in manifest["outputs_before_manifest"].items()),
        f"{len(manifest['outputs_before_manifest'])} outputs verified",
    )

    absolute_pattern = re.compile(r"(?:[A-Za-z]:[\\/]|E:\\\\NPASS|C:\\\\Users)")
    artifacts = [path for path in OUTPUT.iterdir() if path.is_file()]
    check(
        "no_absolute_paths",
        all(not absolute_pattern.search(path.read_text(encoding="utf-8")) for path in artifacts),
        f"{len(artifacts)} text artifacts inspected",
    )

    validation = {
        "analysis_id": "revision_unresolved_bounds_v1_20260729",
        "status": "PASS",
        "check_count": len(checks),
        "checks": checks,
        "validated_output_hashes": {path.name: sha256(path) for path in sorted(artifacts)},
    }
    target = ROOT / "VALIDATION.json"
    require(not target.exists(), "Refusing to overwrite VALIDATION.json")
    target.write_text(json.dumps(validation, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "checks": len(checks)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

