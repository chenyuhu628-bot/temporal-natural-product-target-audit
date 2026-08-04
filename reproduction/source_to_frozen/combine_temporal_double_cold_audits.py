#!/usr/bin/env python3
"""Combine observed temporal future pairs with scaffold and target-coldness audits.

This is an audit join, not a label-generation or negative-sampling program.
It emits only rows already observed in the supplied strict future-pair table.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from reproducible_io import deterministic_gzip_text


TRUE_VALUES = {"1", "true", "t", "yes", "y"}
FALSE_VALUES = {"0", "false", "f", "no", "n", ""}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return deterministic_gzip_text(path) if mode.startswith("w") else gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def delimiter_for(path: Path, declared: str | None = None) -> str:
    if declared == "csv" or (declared is None and path.name.endswith(".csv.gz")) or path.suffix == ".csv":
        return ","
    if declared == "tsv" or (declared is None and path.name.endswith(".tsv.gz")) or path.suffix == ".tsv":
        return "\t"
    raise ValueError(f"Cannot infer delimiter for {path}; set a csv/tsv declaration in the configuration.")


def read_table(path: Path, label: str, declared_delimiter: str | None = None) -> tuple[list[str], list[dict[str, str]]]:
    if not path.is_file():
        raise FileNotFoundError(f"{label} does not exist: {path}")
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter_for(path, declared_delimiter))
        fields = reader.fieldnames
        if not fields:
            raise ValueError(f"{label} has no header: {path}")
        if len(fields) != len(set(fields)):
            raise ValueError(f"{label} has duplicate header names: {path}")
        rows = list(reader)
    return fields, rows


def require_columns(fields: list[str], required: list[str], label: str) -> None:
    missing = [column for column in required if column not in fields]
    if missing:
        raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")


def index_unique(rows: list[dict[str, str]], column: str, label: str) -> dict[str, dict[str, str]]:
    index: dict[str, dict[str, str]] = {}
    duplicates: list[str] = []
    blank_count = 0
    for row in rows:
        key = (row.get(column) or "").strip()
        if not key:
            blank_count += 1
            continue
        if key in index:
            duplicates.append(key)
        else:
            index[key] = row
    if blank_count:
        raise ValueError(f"{label} contains {blank_count} blank values in required key column {column!r}.")
    if duplicates:
        examples = ", ".join(sorted(set(duplicates))[:5])
        raise ValueError(f"{label} must have one row per {column!r}; duplicate values include: {examples}")
    return index


def parse_bool(value: str | None, label: str) -> bool:
    normalized = (value or "").strip().lower()
    if normalized in TRUE_VALUES:
        return True
    if normalized in FALSE_VALUES:
        return False
    raise ValueError(f"Expected a boolean-like value for {label}, got {value!r}.")


def normalized_threshold(value: Any) -> str:
    try:
        threshold = float(str(value))
    except ValueError as exc:
        raise ValueError(f"Invalid MMseqs identity threshold: {value!r}") from exc
    if not 0.0 < threshold <= 1.0:
        raise ValueError(f"MMseqs identity threshold must be in (0, 1], got {value!r}")
    return f"{threshold:.2f}"


def resolve_path(project_root: Path, value: str) -> Path:
    path = Path(value)
    return path if path.is_absolute() else project_root / path


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Configuration must be a JSON object: {path}")
    return data


def parse_threshold_overrides(values: list[str] | None) -> dict[str, str]:
    overrides: dict[str, str] = {}
    for value in values or []:
        if "=" not in value:
            raise ValueError("Each --mmseqs-ledger override must use THRESHOLD=PATH.")
        threshold, path = value.split("=", 1)
        normalized = normalized_threshold(threshold)
        if not path.strip():
            raise ValueError(f"No path supplied for MMseqs threshold {normalized}.")
        overrides[normalized] = path.strip()
    return overrides


def sample(values: set[str], limit: int = 5) -> list[str]:
    return sorted(values)[:limit]


def prefixed_fields(prefix: str, fields: list[str]) -> list[str]:
    return [f"{prefix}__{field}" for field in fields]


def source_values(row: dict[str, str] | None, fields: list[str], prefix: str) -> dict[str, str]:
    return {f"{prefix}__{field}": "" if row is None else row.get(field, "") for field in fields}


def stable_counts(rows: list[dict[str, str]], column: str) -> dict[str, int]:
    counter = Counter((row.get(column) or "").strip() or "<blank>" for row in rows)
    return dict(sorted(counter.items()))


def make_argument_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", required=True, help="Project root used to resolve relative configuration paths.")
    parser.add_argument("--config", required=True, help="JSON configuration path, relative to --project-root or absolute.")
    parser.add_argument("--future-pairs", help="Override inputs.future_pair_table.")
    parser.add_argument("--scaffold-ledger", help="Override inputs.scaffold_future_pair_ledger.")
    parser.add_argument("--mmseqs-ledger", action="append", metavar="THRESHOLD=PATH", help="Override one MMseqs target ledger; repeat for multiple thresholds.")
    parser.add_argument("--output-dir", help="Override output_dir. The directory must not already exist.")
    parser.add_argument("--check-only", action="store_true", help="Validate all inputs and print a compact report without creating outputs.")
    return parser


def main() -> int:
    args = make_argument_parser().parse_args()
    project_root = Path(args.project_root).resolve()
    config_path = resolve_path(project_root, args.config).resolve()
    config = load_json(config_path)

    inputs = config.get("inputs", {})
    columns = config.get("columns", {})
    policies = config.get("policies", {})
    if not isinstance(inputs, dict) or not isinstance(columns, dict) or not isinstance(policies, dict):
        raise ValueError("Configuration inputs, columns, and policies must each be JSON objects.")

    future_path_value = args.future_pairs or inputs.get("future_pair_table")
    scaffold_path_value = args.scaffold_ledger or inputs.get("scaffold_future_pair_ledger")
    output_path_value = args.output_dir or config.get("output_dir")
    if not future_path_value or not scaffold_path_value or not output_path_value:
        raise ValueError("Configuration requires future_pair_table, scaffold_future_pair_ledger, and output_dir.")
    future_path = resolve_path(project_root, str(future_path_value)).resolve()
    scaffold_path = resolve_path(project_root, str(scaffold_path_value)).resolve()
    output_dir = resolve_path(project_root, str(output_path_value)).resolve()

    configured_ledgers = inputs.get("mmseqs_future_target_ledgers", {})
    if not isinstance(configured_ledgers, dict):
        raise ValueError("inputs.mmseqs_future_target_ledgers must be an object keyed by identity threshold.")
    mmseqs_paths = {normalized_threshold(key): str(value) for key, value in configured_ledgers.items()}
    mmseqs_paths.update(parse_threshold_overrides(args.mmseqs_ledger))
    thresholds = [normalized_threshold(value) for value in config.get("identity_thresholds", sorted(mmseqs_paths))]
    if not thresholds or len(thresholds) != len(set(thresholds)):
        raise ValueError("identity_thresholds must contain one or more unique values.")
    missing_threshold_paths = [threshold for threshold in thresholds if threshold not in mmseqs_paths]
    if missing_threshold_paths:
        raise ValueError(f"No MMseqs target ledger configured for: {', '.join(missing_threshold_paths)}")

    future_pair_key = columns.get("future_pair_key", "canonical_pair_key")
    future_target = columns.get("future_target", "uniprot_canonical_accession")
    future_decision = columns.get("future_decision", "decision")
    future_label_status = columns.get("future_label_status", "label_status")
    future_unrecorded_policy = columns.get("future_unrecorded_pair_policy", "unrecorded_pair_policy")
    scaffold_pair_key = columns.get("scaffold_pair_key", "canonical_pair_key")
    scaffold_cold_flag = columns.get("scaffold_cold_flag", "audit_scaffold_cold_under_selected_policy")
    scaffold_outcome = columns.get("scaffold_outcome", "audit_outcome")
    scaffold_reason = columns.get("scaffold_reason", "audit_eligibility_or_exclusion_reason")
    target_accession = columns.get("target_accession", "uniprot_canonical_accession")
    target_cold_flag = columns.get("target_cold_flag", "is_future_target_homology_cold_candidate")
    target_status = columns.get("target_status", "future_target_coldness_status")

    # An explicit CLI path may be a different (but supported) delimited
    # artifact than the configured default, for example the post-C31 TSV
    # future table. In that case infer its delimiter from its own suffix.
    future_declared_delimiter = None if args.future_pairs else inputs.get("future_pair_table_format")
    scaffold_declared_delimiter = None if args.scaffold_ledger else inputs.get("scaffold_future_pair_ledger_format")
    future_fields, future_rows = read_table(future_path, "Strict temporal future-pair table", future_declared_delimiter)
    scaffold_fields, scaffold_rows = read_table(scaffold_path, "Scaffold-coldness audit ledger", scaffold_declared_delimiter)
    require_columns(future_fields, [future_pair_key, future_target, future_decision, future_label_status, future_unrecorded_policy], "Strict temporal future-pair table")
    require_columns(scaffold_fields, [scaffold_pair_key, scaffold_cold_flag, scaffold_outcome, scaffold_reason], "Scaffold-coldness audit ledger")
    future_by_pair = index_unique(future_rows, future_pair_key, "Strict temporal future-pair table")
    scaffold_by_pair = index_unique(scaffold_rows, scaffold_pair_key, "Scaffold-coldness audit ledger")
    future_targets = {row[future_target].strip() for row in future_rows}

    expected_decision = policies.get("expected_future_decision", "strict_post_cutoff_future_candidate")
    temporal_invalid_keys = {
        key for key, row in future_by_pair.items() if (row.get(future_decision) or "").strip() != expected_decision
    }
    p1_marker = str(policies.get("p1_candidate_marker", "P1_candidate_only"))
    p1_invalid_keys = {
        key for key, row in future_by_pair.items() if p1_marker not in (row.get(future_label_status) or "")
    }
    if policies.get("require_expected_future_decision", True) and temporal_invalid_keys:
        raise ValueError(f"Future table contains {len(temporal_invalid_keys)} row(s) without expected decision {expected_decision!r}; examples: {sample(temporal_invalid_keys)}")
    if policies.get("require_p1_candidate_marker", True) and p1_invalid_keys:
        raise ValueError(f"Future table contains {len(p1_invalid_keys)} row(s) without P1 marker {p1_marker!r}; examples: {sample(p1_invalid_keys)}")

    missing_scaffold_keys = set(future_by_pair) - set(scaffold_by_pair)
    extra_scaffold_keys = set(scaffold_by_pair) - set(future_by_pair)
    if policies.get("require_exact_scaffold_pair_keyset", True) and (missing_scaffold_keys or extra_scaffold_keys):
        raise ValueError(
            "Strict future/scaffold pair-key alignment failed: "
            f"missing ledger rows={len(missing_scaffold_keys)} examples={sample(missing_scaffold_keys)}; "
            f"extra ledger rows={len(extra_scaffold_keys)} examples={sample(extra_scaffold_keys)}"
        )

    target_ledgers: dict[str, tuple[list[str], dict[str, dict[str, str]], Path]] = {}
    target_alignment: dict[str, dict[str, Any]] = {}
    for threshold in thresholds:
        target_path = resolve_path(project_root, mmseqs_paths[threshold]).resolve()
        fields, rows = read_table(target_path, f"MMseqs target audit at {threshold}", inputs.get("mmseqs_future_target_ledger_format"))
        require_columns(fields, [target_accession, target_cold_flag, target_status], f"MMseqs target audit at {threshold}")
        by_target = index_unique(rows, target_accession, f"MMseqs target audit at {threshold}")
        missing_targets = future_targets - set(by_target)
        extra_targets = set(by_target) - future_targets
        if policies.get("require_exact_mmseqs_target_keyset", True) and (missing_targets or extra_targets):
            raise ValueError(
                f"Strict future/MMseqs target alignment failed at {threshold}: "
                f"missing={len(missing_targets)} examples={sample(missing_targets)}; "
                f"extra={len(extra_targets)} examples={sample(extra_targets)}"
            )
        target_ledgers[threshold] = (fields, by_target, target_path)
        target_alignment[threshold] = {
            "rows": len(rows),
            "missing_future_targets": len(missing_targets),
            "extra_audit_targets": len(extra_targets),
            "missing_target_examples": sample(missing_targets),
            "extra_target_examples": sample(extra_targets),
        }

    validation = {
        "future_pairs": len(future_rows),
        "future_targets": len(future_targets),
        "temporal_decision_mismatch_count": len(temporal_invalid_keys),
        "p1_marker_mismatch_count": len(p1_invalid_keys),
        "scaffold_missing_future_pair_count": len(missing_scaffold_keys),
        "scaffold_extra_pair_count": len(extra_scaffold_keys),
        "mmseqs_target_alignment": target_alignment,
    }
    if args.check_only:
        print(json.dumps({"status": "validation_passed", "config": str(config_path), "validation": validation}, ensure_ascii=False, indent=2))
        return 0
    if output_dir.exists():
        raise FileExistsError(f"Output directory must be new and must not be overwritten: {output_dir}")

    scaffold_prefix = "scaffold_ledger"
    target_prefixes = {threshold: f"mmseqs_identity_{threshold.replace('.', '_')}" for threshold in thresholds}
    combined_fields = [
        "combined_observed_future_pair_only",
        "combined_temporal_future_eligible",
        "combined_temporal_reason",
        "combined_scaffold_cold_eligible",
        "combined_scaffold_reason",
    ]
    for threshold in thresholds:
        combined_fields.extend([
            f"combined_homology_cold_eligible_identity_{threshold.replace('.', '_')}",
            f"combined_homology_reason_identity_{threshold.replace('.', '_')}",
            f"combined_temporal_scaffold_homology_eligible_identity_{threshold.replace('.', '_')}",
            f"combined_temporal_scaffold_homology_reason_identity_{threshold.replace('.', '_')}",
        ])
    combined_fields.extend([
        "combined_label_status_preserved",
        "combined_label_policy",
        "combined_unrecorded_pair_policy",
        "combined_negative_label_count_for_row",
        "combined_label_promotion_count_for_row",
    ])
    output_fields = list(future_fields) + prefixed_fields(scaffold_prefix, scaffold_fields)
    for threshold in thresholds:
        output_fields += prefixed_fields(target_prefixes[threshold], target_ledgers[threshold][0])
    output_fields += combined_fields

    output_dir.mkdir(parents=True, exist_ok=False)
    ledger_path = output_dir / "temporal_scaffold_homology_pair_eligibility.tsv.gz"
    counts_path = output_dir / "threshold_pair_eligibility_counts.csv"
    manifest_path = output_dir / "summary.json"
    config_snapshot_path = output_dir / "config_snapshot.json"

    counts: dict[str, Any] = {
        "observed_future_pair_count": 0,
        "temporal_future_eligible_pair_count": 0,
        "p1_candidate_only_pair_count": 0,
        "scaffold_cold_eligible_pair_count": 0,
        "scaffold_outcomes": Counter(),
        "thresholds": {threshold: Counter() for threshold in thresholds},
    }
    with open_text(ledger_path, "wt") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_fields, delimiter="\t", lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        for pair_key, future_row in future_by_pair.items():
            row = dict(future_row)
            scaffold_row = scaffold_by_pair.get(pair_key)
            row.update(source_values(scaffold_row, scaffold_fields, scaffold_prefix))
            temporal_eligible = (future_row.get(future_decision) or "").strip() == expected_decision
            temporal_reason = "expected_strict_post_cutoff_future_decision" if temporal_eligible else f"unexpected_future_decision:{future_row.get(future_decision, '')}"
            scaffold_eligible = False
            scaffold_reason_value = "missing_scaffold_ledger_row"
            if scaffold_row is not None:
                scaffold_eligible = parse_bool(scaffold_row.get(scaffold_cold_flag), f"{pair_key} scaffold coldness")
                scaffold_reason_value = (scaffold_row.get(scaffold_reason) or scaffold_row.get(scaffold_outcome) or "scaffold_audit_reason_blank").strip()
                counts["scaffold_outcomes"][(scaffold_row.get(scaffold_outcome) or "<blank>").strip() or "<blank>"] += 1
            label_status = future_row.get(future_label_status, "")
            row.update({
                "combined_observed_future_pair_only": "True",
                "combined_temporal_future_eligible": str(temporal_eligible),
                "combined_temporal_reason": temporal_reason,
                "combined_scaffold_cold_eligible": str(scaffold_eligible),
                "combined_scaffold_reason": scaffold_reason_value,
                "combined_label_status_preserved": label_status,
                "combined_label_policy": "Source label status is preserved. This audit does not promote P1 candidates to P2/P3 or create a supervised positive label.",
                "combined_unrecorded_pair_policy": "Only rows observed in the supplied strict temporal future-pair table are emitted; no unrecorded pair is created or labeled negative.",
                "combined_negative_label_count_for_row": "0",
                "combined_label_promotion_count_for_row": "0",
            })
            for threshold in thresholds:
                fields, by_target, _ = target_ledgers[threshold]
                prefix = target_prefixes[threshold]
                target_row = by_target.get((future_row.get(future_target) or "").strip())
                row.update(source_values(target_row, fields, prefix))
                homology_eligible = False
                homology_reason = "missing_mmseqs_target_audit_row"
                if target_row is not None:
                    homology_eligible = parse_bool(target_row.get(target_cold_flag), f"{pair_key} MMseqs {threshold} coldness")
                    homology_reason = (target_row.get(target_status) or "target_audit_status_blank").strip()
                combined_eligible = temporal_eligible and scaffold_eligible and homology_eligible
                if combined_eligible:
                    combined_reason = "eligible_observed_P1_candidate_only__temporal_future_plus_scaffold_cold_plus_homology_cold"
                else:
                    failed = []
                    if not temporal_eligible:
                        failed.append(f"temporal:{temporal_reason}")
                    if not scaffold_eligible:
                        failed.append(f"scaffold:{scaffold_reason_value}")
                    if not homology_eligible:
                        failed.append(f"homology_{threshold}:{homology_reason}")
                    combined_reason = "; ".join(failed)
                suffix = threshold.replace(".", "_")
                row[f"combined_homology_cold_eligible_identity_{suffix}"] = str(homology_eligible)
                row[f"combined_homology_reason_identity_{suffix}"] = homology_reason
                row[f"combined_temporal_scaffold_homology_eligible_identity_{suffix}"] = str(combined_eligible)
                row[f"combined_temporal_scaffold_homology_reason_identity_{suffix}"] = combined_reason
                threshold_counts = counts["thresholds"][threshold]
                threshold_counts["future_pairs_with_homology_cold_target"] += int(homology_eligible)
                threshold_counts["temporal_scaffold_homology_eligible_pairs"] += int(combined_eligible)
            writer.writerow(row)
            counts["observed_future_pair_count"] += 1
            counts["temporal_future_eligible_pair_count"] += int(temporal_eligible)
            counts["p1_candidate_only_pair_count"] += int(p1_marker in label_status)
            counts["scaffold_cold_eligible_pair_count"] += int(scaffold_eligible)

    count_rows = []
    for threshold in thresholds:
        by_target = target_ledgers[threshold][1]
        count_rows.append({
            "identity_threshold": threshold,
            "observed_future_pair_count": counts["observed_future_pair_count"],
            "temporal_future_eligible_pair_count": counts["temporal_future_eligible_pair_count"],
            "scaffold_cold_eligible_pair_count": counts["scaffold_cold_eligible_pair_count"],
            "homology_cold_future_target_count": sum(parse_bool(row.get(target_cold_flag), f"MMseqs {threshold} summary") for row in by_target.values()),
            "future_pairs_with_homology_cold_target": counts["thresholds"][threshold]["future_pairs_with_homology_cold_target"],
            "temporal_scaffold_homology_eligible_pairs": counts["thresholds"][threshold]["temporal_scaffold_homology_eligible_pairs"],
            "negative_labels_created": 0,
            "label_promotions_created": 0,
        })
    with counts_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(count_rows[0]))
        writer.writeheader()
        writer.writerows(count_rows)

    config_snapshot_path.write_text(json.dumps(config, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    manifest = {
        "schema_version": "1.0",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "purpose": "Pair-level join of observed strict temporal future pairs with precomputed scaffold and MMseqs homology-coldness audits. No interaction labels or unrecorded pairs are created.",
        "config": {"path": str(config_path), "sha256": sha256_file(config_path)},
        "inputs": {
            "future_pair_table": {"path": str(future_path), "sha256": sha256_file(future_path)},
            "scaffold_future_pair_ledger": {"path": str(scaffold_path), "sha256": sha256_file(scaffold_path)},
            "mmseqs_future_target_ledgers": {threshold: {"path": str(path), "sha256": sha256_file(path)} for threshold, (_, _, path) in target_ledgers.items()},
        },
        "policies": policies,
        "validation": validation,
        "counts": {
            "observed_future_pair_count": counts["observed_future_pair_count"],
            "temporal_future_eligible_pair_count": counts["temporal_future_eligible_pair_count"],
            "p1_candidate_only_pair_count": counts["p1_candidate_only_pair_count"],
            "scaffold_cold_eligible_pair_count": counts["scaffold_cold_eligible_pair_count"],
            "scaffold_outcomes": dict(sorted(counts["scaffold_outcomes"].items())),
            "thresholds": {threshold: dict(counts["thresholds"][threshold]) for threshold in thresholds},
            "unrecorded_pairs_emitted": 0,
            "negative_labels_created": 0,
            "label_promotions_created": 0,
        },
        "outputs": {
            "pair_eligibility_ledger": {"path": str(ledger_path), "sha256": sha256_file(ledger_path)},
            "threshold_pair_eligibility_counts": {"path": str(counts_path), "sha256": sha256_file(counts_path)},
            "config_snapshot": {"path": str(config_snapshot_path), "sha256": sha256_file(config_snapshot_path)},
        },
        "limitations": [
            "Every output row was already observed in the input strict temporal future-pair table; absence is never converted into a negative label.",
            "Source label_status values are retained verbatim. In the configured P1 mode, a row is a candidate only and not a final positive, direct-binding claim, or assay-verified label.",
            "Scaffold and homology flags are split-eligibility audit fields under their own frozen protocols, not guarantees of chemical or biological novelty.",
        ],
    }
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "output_dir": str(output_dir), "counts": manifest["counts"]}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, FileExistsError, ValueError, OSError, json.JSONDecodeError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
