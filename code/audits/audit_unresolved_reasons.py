"""Aggregate-only reason decomposition for the frozen C31 entity-unresolved relations."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from audit_common import (
    PROTOCOL_ID,
    choose_field,
    distribution,
    finalize_manifest,
    input_descriptor,
    open_dict_reader,
    parse_bool,
    require_new_output_dir,
    require_protocol_lock,
    write_json_new,
    write_tsv_new,
)


AUDIT_ID = "frozen_entity_unresolved_reason_decomposition_v1"


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--protocol-lock", required=True, type=Path)
    result.add_argument("--unresolved-ledger", required=True, type=Path)
    result.add_argument("--preliminary-mapping-ledger", required=True, type=Path)
    result.add_argument("--sqlite-validation-ledger", required=True, type=Path)
    result.add_argument("--output-dir", required=True, type=Path)
    result.add_argument("--expected-relation-count", type=int, default=65)
    return result


def load_preliminary(path: Path) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    with open_dict_reader(path) as reader:
        fields = reader.fieldnames or []
        pair_field = choose_field(fields, ["pair_key", "canonical_pair_key"], "preliminary mapping")
        required = {
            "chembl_compound_match_status",
            "chembl_target_match_status",
            "both_entities_exactly_mapped",
        }
        missing = required.difference(fields)
        if missing:
            raise ValueError(f"Preliminary mapping lacks required fields: {sorted(missing)}")
        for row in reader:
            pair = row[pair_field].strip()
            if pair in result:
                raise ValueError("Preliminary mapping has duplicate relation keys")
            result[pair] = row
    return result


def load_sqlite_rows(path: Path) -> dict[str, list[dict[str, str]]]:
    result: dict[str, list[dict[str, str]]] = defaultdict(list)
    with open_dict_reader(path) as reader:
        fields = reader.fieldnames or []
        pair_field = choose_field(fields, ["pair_key", "canonical_pair_key"], "SQLite validation")
        required = {
            "sqlite_compound_full_inchikey_exact",
            "sqlite_target_single_protein",
            "sqlite_target_human",
            "sqlite_target_source_uniprot_exact",
            "sqlite_entity_pair_validated",
            "sqlite_entity_validation_status",
        }
        missing = required.difference(fields)
        if missing:
            raise ValueError(f"SQLite validation lacks required fields: {sorted(missing)}")
        for row in reader:
            result[row[pair_field].strip()].append(row)
    return result


def sqlite_failure_reason(rows: list[dict[str, str]]) -> str:
    if not rows:
        return "preliminary_both_exact_no_sqlite_candidate_rows"
    if any(parse_bool(row["sqlite_entity_pair_validated"]) for row in rows):
        raise ValueError("Frozen unresolved relation unexpectedly has a validated SQLite entity pair")
    compound_ok = [parse_bool(row["sqlite_compound_full_inchikey_exact"]) for row in rows]
    target_ok = [
        parse_bool(row["sqlite_target_single_protein"])
        and parse_bool(row["sqlite_target_human"])
        and parse_bool(row["sqlite_target_source_uniprot_exact"])
        for row in rows
    ]
    joint_ok = [left and right for left, right in zip(compound_ok, target_ok)]
    if any(joint_ok):
        raise ValueError("SQLite flags imply a valid joint entity pair but validated flag is false")
    if not any(compound_ok) and not any(target_ok):
        return "sqlite_compound_and_target_validation_failed"
    if not any(compound_ok):
        return "sqlite_compound_validation_failed"
    if not any(target_ok):
        return "sqlite_target_validation_failed"
    return "sqlite_no_joint_valid_mapping"


def main() -> int:
    args = parser().parse_args()
    require_protocol_lock(args.protocol_lock)
    output_dir = require_new_output_dir(args.output_dir)
    inputs = [
        input_descriptor("frozen_entity_unresolved_ledger", args.unresolved_ledger),
        input_descriptor("preliminary_entity_mapping_ledger", args.preliminary_mapping_ledger),
        input_descriptor("sqlite_entity_validation_ledger", args.sqlite_validation_ledger),
    ]
    preliminary = load_preliminary(args.preliminary_mapping_ledger)
    sqlite_rows = load_sqlite_rows(args.sqlite_validation_ledger)

    reason_counts: Counter[str] = Counter()
    compounds: Counter[str] = Counter()
    targets: Counter[str] = Counter()
    seen_pairs: set[str] = set()
    with open_dict_reader(args.unresolved_ledger) as reader:
        fields = reader.fieldnames or []
        pair_field = choose_field(
            fields,
            ["audit_pair_key", "canonical_pair_key", "pair_key"],
            "unresolved ledger",
        )
        compound_field = choose_field(
            fields,
            ["inchikey_full", "audit_inchikey_full"],
            "unresolved ledger",
        )
        target_field = choose_field(
            fields,
            ["uniprot_canonical_accession", "audit_npass_uniprot_source"],
            "unresolved ledger",
        )
        for row in reader:
            pair = row[pair_field].strip()
            if not pair:
                raise ValueError("Unresolved ledger contains an empty relation key")
            if pair in seen_pairs:
                raise ValueError("Unresolved ledger contains duplicate relations")
            seen_pairs.add(pair)
            compound = row[compound_field].strip()
            target = row[target_field].strip()
            if not compound or not target:
                raise ValueError("Unresolved ledger contains an empty entity identifier")
            compounds[compound] += 1
            targets[target] += 1

            if "leakage_gate_stratum" in row and row["leakage_gate_stratum"].strip() != "C31_entity_unresolved":
                raise ValueError("Frozen unresolved ledger contains another leakage stratum")
            if "audit_sqlite_validated_entity_mapping_count" in row:
                if int(row["audit_sqlite_validated_entity_mapping_count"]) != 0:
                    raise ValueError("Unresolved relation has a positive validated mapping count")
            if "negative_label_emitted" in row and parse_bool(row["negative_label_emitted"]):
                raise ValueError("Unresolved relation must not emit a negative label")

            mapping = preliminary.get(pair)
            if mapping is None:
                raise ValueError("Frozen unresolved relation is absent from preliminary mapping ledger")
            compound_exact = mapping["chembl_compound_match_status"] == "full_inchikey_exact"
            target_exact = mapping["chembl_target_match_status"] == "source_uniprot_exact"
            declared_both = parse_bool(mapping["both_entities_exactly_mapped"])
            if declared_both != (compound_exact and target_exact):
                raise ValueError("Preliminary both-entity flag contradicts component match statuses")
            if not compound_exact and not target_exact:
                reason = "preliminary_compound_and_target_unmatched"
            elif not compound_exact:
                reason = "preliminary_compound_unmatched"
            elif not target_exact:
                reason = "preliminary_target_unmatched"
            else:
                reason = sqlite_failure_reason(sqlite_rows.get(pair, []))
            reason_counts[reason] += 1

    if len(seen_pairs) != args.expected_relation_count:
        raise ValueError(
            f"Frozen unresolved relation count is {len(seen_pairs)}, "
            f"expected {args.expected_relation_count}"
        )
    if sum(reason_counts.values()) != len(seen_pairs):
        raise ValueError("Reason decomposition is not exhaustive")

    reason_rows = [
        {
            "reason_category": reason,
            "relation_count": count,
            "relation_fraction": count / len(seen_pairs),
        }
        for reason, count in sorted(reason_counts.items())
    ]
    summary = {
        "audit_id": AUDIT_ID,
        "protocol_id": PROTOCOL_ID,
        "frozen_unresolved_relation_count": len(seen_pairs),
        "distinct_compound_count": len(compounds),
        "distinct_target_count": len(targets),
        "relations_per_compound_distribution": distribution(list(compounds.values())),
        "relations_per_target_distribution": distribution(list(targets.values())),
        "reason_decomposition": reason_rows,
        "negative_labels_emitted": 0,
        "readmission_performed": False,
        "interpretation_boundary": (
            "This is an aggregate decomposition of frozen exclusions; no relation is silently "
            "reintroduced and unresolved absence is not a negative label"
        ),
    }
    output_dir.mkdir(parents=True, exist_ok=False)
    summary_name = "unresolved_reason_aggregate_summary.json"
    table_name = "unresolved_reason_counts.tsv"
    write_json_new(output_dir / summary_name, summary)
    write_tsv_new(
        output_dir / table_name,
        ["reason_category", "relation_count", "relation_fraction"],
        reason_rows,
    )
    manifest = finalize_manifest(
        output_dir=output_dir,
        audit_id=AUDIT_ID,
        script_path=Path(__file__),
        inputs=inputs,
        output_names=[summary_name, table_name],
        extra={"expected_relation_count": args.expected_relation_count},
    )
    write_json_new(output_dir / "run_manifest.json", manifest)
    print(f"{AUDIT_ID}: decomposed {len(seen_pairs)} relations without identifier output")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
