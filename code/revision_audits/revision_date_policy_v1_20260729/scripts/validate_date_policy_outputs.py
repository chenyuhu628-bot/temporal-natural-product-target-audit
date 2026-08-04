"""Validate aggregate-only date-policy outputs and write a validation receipt."""

from __future__ import annotations

import csv
import hashlib
import json
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "outputs"
SCENARIOS = {
    "day_only_conservative",
    "interval_certain_pre_cutoff",
    "interval_earliest_bound",
}
EXPECTED_ROWS = {
    "interval_status_counts.tsv": 4,
    "scenario_history_summary.tsv": 3,
    "scenario_structure_summary.tsv": 3,
    "scope_denominators.tsv": 15,
    "recall_at_50.tsv": 60,
    "score_rank_change_summary.tsv": 12,
    "scenario_equivalence.tsv": 3,
}
FORBIDDEN_HEADERS = {
    "canonical_pair_key",
    "inchikey_full",
    "uniprot_canonical_accession",
    "target_uniprot_accession",
    "query_id",
    "representative_smiles",
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
        require(not FORBIDDEN_HEADERS.intersection(reader.fieldnames or []), f"Identifier-bearing header: {name}")
        return list(reader)


def main() -> int:
    checks: list[dict[str, object]] = []

    def check(name: str, condition: bool, detail: str) -> None:
        require(condition, f"{name}: {detail}")
        checks.append({"check": name, "status": "PASS", "detail": detail})

    tables = {name: read_tsv(name) for name in EXPECTED_ROWS}
    check(
        "aggregate_row_counts",
        all(len(tables[name]) == count for name, count in EXPECTED_ROWS.items()),
        str({name: len(rows) for name, rows in tables.items()}),
    )
    observed_scenarios = {
        row["scenario"]
        for name, rows in tables.items()
        if name != "interval_status_counts.tsv"
        for row in rows
    }
    check("scenario_completeness", observed_scenarios == SCENARIOS, str(sorted(observed_scenarios)))

    interval = {row["interval_status"]: int(row["source_row_count"]) for row in tables["interval_status_counts.tsv"]}
    check(
        "interval_partition",
        interval
        == {
            "definitely_before_or_on": 20455,
            "crossing_cutoff": 0,
            "definitely_after": 0,
            "unresolved_interval": 192,
        },
        str(interval),
    )

    history = {row["scenario"]: row for row in tables["scenario_history_summary.tsv"]}
    check("day_only_rows", int(history["day_only_conservative"]["selected_source_row_count"]) == 13885, "13,885")
    check(
        "interval_rows",
        int(history["interval_certain_pre_cutoff"]["selected_source_row_count"]) == 20455
        and int(history["interval_earliest_bound"]["selected_source_row_count"]) == 20455,
        "20,455 in both interval policies",
    )
    check(
        "history_membership_fixed",
        all(int(row["historical_pair_count"]) == 4990 and int(row["historical_membership_change_count_vs_day_only"]) == 0 for row in history.values()),
        "4,990 pairs; zero membership changes",
    )

    equivalence = {row["scenario"]: row for row in tables["scenario_equivalence.tsv"]}
    check(
        "interval_policy_state_equivalence",
        equivalence["interval_certain_pre_cutoff"]["scoring_state_sha256"]
        == equivalence["interval_earliest_bound"]["scoring_state_sha256"],
        "certain-pre and earliest-bound collapse to one scoring state",
    )

    scopes = tables["scope_denominators.tsv"]
    for scenario in SCENARIOS:
        per = [row for row in scopes if row["scenario"] == scenario]
        half = next(row for row in per if row["provenance_scope"] == "joint_scaffold_homology_0.50")
        seven = next(row for row in per if row["provenance_scope"] == "joint_scaffold_homology_0.70")
        check(
            f"homology_mask_identity_{scenario}",
            half["membership_sha256"] == seven["membership_sha256"]
            and half["candidate_relation_count"] == seven["candidate_relation_count"]
            and half["query_count"] == seven["query_count"],
            "0.50 and 0.70 joint masks are identical",
        )

    recall = tables["recall_at_50.tsv"]
    for scenario in SCENARIOS:
        half = {
            row["baseline"]: row["Recall@50"]
            for row in recall
            if row["scenario"] == scenario and row["provenance_scope"] == "joint_scaffold_homology_0.50"
        }
        seven = {
            row["baseline"]: row["Recall@50"]
            for row in recall
            if row["scenario"] == scenario and row["provenance_scope"] == "joint_scaffold_homology_0.70"
        }
        check(f"homology_metric_identity_{scenario}", half == seven, "Recall@50 rows are identical")

    receipt = json.loads((OUTPUT / "EXECUTION_RECEIPT.json").read_text(encoding="utf-8"))
    check("execution_receipt", receipt.get("status") == "PASS", "runner status PASS")
    check(
        "frozen_endpoint_reproduction",
        int(receipt.get("frozen_endpoint_score_rank_cells_reproduced", 0)) == 1432,
        "1,432 endpoint baseline cells",
    )
    check(
        "aggregate_only_receipt",
        receipt.get("identifier_bearing_rows_emitted") is False
        and receipt.get("absolute_paths_emitted") is False
        and receipt.get("main_manuscript_modified") is False,
        "no identifier rows, absolute paths, or manuscript writes",
    )

    manifest = json.loads((OUTPUT / "MANIFEST.json").read_text(encoding="utf-8"))
    check(
        "parent_protocol_lock",
        manifest.get("parent_protocol_sha256")
        == "bd4902476160cc7c5cbacaf0cfd0f1a28c5300bde22232b3a2cc6c1f3c143dc2",
        "frozen v4 protocol hash",
    )
    output_hashes = manifest["outputs_before_manifest"]
    check(
        "manifest_output_hashes",
        all((OUTPUT / name).is_file() and sha256(OUTPUT / name) == value for name, value in output_hashes.items()),
        f"{len(output_hashes)} pre-manifest outputs",
    )

    absolute_pattern = re.compile(r"(?:[A-Za-z]:[\\/]|E:\\\\NPASS|C:\\\\Users)")
    inspected = [path for path in OUTPUT.iterdir() if path.is_file()]
    check(
        "no_absolute_paths",
        all(not absolute_pattern.search(path.read_text(encoding="utf-8")) for path in inspected),
        f"{len(inspected)} text artifacts inspected",
    )

    validation_path = ROOT / "VALIDATION.json"
    require(not validation_path.exists(), "Refusing to overwrite VALIDATION.json")
    payload = {
        "analysis_id": "revision_date_policy_v1_20260729",
        "status": "PASS",
        "check_count": len(checks),
        "checks": checks,
        "validated_output_hashes": {
            path.name: sha256(path) for path in sorted(OUTPUT.iterdir()) if path.is_file()
        },
    }
    validation_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": "PASS", "checks": len(checks)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

