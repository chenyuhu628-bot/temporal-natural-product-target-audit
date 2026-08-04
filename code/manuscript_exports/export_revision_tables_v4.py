#!/usr/bin/env python3
"""Export manuscript-v4 tables mechanically from locked corrective and revision runs.

The exporter emits aggregate-only TSV files.  It does not parse or expose
identifier-rich endpoint or rank ledgers; the complete rank ledger is read only
as bytes for custody hashing.  Every input and output is recorded by SHA-256 in
a relative-path manifest so that the public/review-safe package can scan it.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import itertools
import json
import os
import re
import tempfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable


RUN_DIRNAME = "author_run_strict_ab_asof_cutoff_execution_v1_20260728"
MANUSCRIPT_DIRNAME = "manuscript_molecular_diversity_v4_20260729"
STRICT_DIRNAME = "strict_ab_asof_cutoff_successor_v1_20260728"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"Expected JSON object: {path}")
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    kwargs: dict[str, Any] = {"mode": "rt", "encoding": "utf-8-sig", "newline": ""}
    with opener(path, **kwargs) as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows:
        raise ValueError(f"Empty TSV: {path}")
    return rows


def atomic_write_tsv(path: Path, fieldnames: list[str], rows: Iterable[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    os.close(fd)
    temporary_path = Path(temporary)
    try:
        with temporary_path.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fieldnames,
                delimiter="\t",
                lineterminator="\n",
                extrasaction="raise",
            )
            writer.writeheader()
            for row in rows:
                writer.writerow({field: row.get(field, "") for field in fieldnames})
        os.replace(temporary_path, path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def atomic_write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", suffix=".tmp", dir=path.parent)
    os.close(fd)
    temporary_path = Path(temporary)
    try:
        with temporary_path.open("w", encoding="utf-8", newline="\n") as handle:
            json.dump(payload, handle, indent=2, ensure_ascii=False, sort_keys=True)
            handle.write("\n")
        os.replace(temporary_path, path)
    finally:
        if temporary_path.exists():
            temporary_path.unlink()


def rel(path: Path, root: Path) -> str:
    return path.resolve().relative_to(root.resolve()).as_posix()


def record_input(path: Path, project_root: Path, registry: dict[str, dict[str, Any]]) -> None:
    resolved = path.resolve()
    key = rel(resolved, project_root)
    registry[key] = {"bytes": resolved.stat().st_size, "sha256": sha256_file(resolved)}


def key_value_rows(section: str, mapping: dict[str, Any], status: str = "verified") -> list[dict[str, Any]]:
    return [
        {"section": section, "item": key, "value": value, "status": status}
        for key, value in sorted(mapping.items())
    ]


def require_exact_set(observed: set[str], expected: set[str], label: str) -> None:
    if observed != expected:
        raise ValueError(f"{label} mismatch: observed={sorted(observed)} expected={sorted(expected)}")


def require_columns(rows: list[dict[str, str]], expected: set[str], label: str) -> None:
    observed = set(rows[0])
    if observed != expected:
        raise ValueError(f"{label} columns mismatch: observed={sorted(observed)} expected={sorted(expected)}")


def require_unique_keys(rows: list[dict[str, str]], fields: tuple[str, ...], expected: set[tuple[str, ...]], label: str) -> None:
    keys = [tuple(row[field] for field in fields) for row in rows]
    if len(keys) != len(set(keys)):
        raise ValueError(f"{label} contains duplicate keys")
    if set(keys) != expected:
        missing = expected - set(keys)
        extra = set(keys) - expected
        raise ValueError(f"{label} grid mismatch: missing={sorted(missing)} extra={sorted(extra)}")


def declared_output_sha256(manifest: dict[str, Any], basename: str) -> str:
    outputs = manifest.get("outputs")
    if isinstance(outputs, list):
        matches = [item for item in outputs if isinstance(item, dict) and item.get("basename") == basename]
    elif isinstance(outputs, dict):
        matches = []
        for item in outputs.values():
            if not isinstance(item, dict):
                continue
            candidate = item.get("basename") or Path(str(item.get("path", ""))).name
            if candidate == basename:
                matches.append(item)
    else:
        matches = []
    if len(matches) != 1 or not isinstance(matches[0].get("sha256"), str):
        raise ValueError(f"Exactly one declared output hash required for {basename}")
    return str(matches[0]["sha256"])


def verify_manifest_outputs(manifest_path: Path, output_paths: Iterable[Path], protocol_id: str, protocol_lock_sha256: str | None = None) -> dict[str, Any]:
    manifest = load_json(manifest_path)
    if manifest.get("protocol_id") != protocol_id:
        raise ValueError(f"Protocol mismatch in {manifest_path.name}")
    declared_protocol_sha = manifest.get("protocol_lock_sha256")
    if protocol_lock_sha256 is not None and declared_protocol_sha is not None and declared_protocol_sha != protocol_lock_sha256:
        raise ValueError(f"Protocol-lock hash mismatch in {manifest_path.name}")
    for path in output_paths:
        declared = declared_output_sha256(manifest, path.name)
        actual = sha256_file(path)
        if actual != declared:
            raise ValueError(f"Upstream output hash mismatch for {path}")
    return manifest


def scan_aggregate_outputs(paths: Iterable[Path]) -> None:
    absolute_windows = re.compile(r"(?i)(?:^|[\s\t\"'])[A-Z]:\\")
    unc_path = re.compile(r"\\\\[^\s\\]+\\[^\s\\]+")
    absolute_posix = re.compile(r"(?:^|[\s\t\"'])/(?:[^/\s\t\"']+/)+[^/\s\t\"']+")
    inchi_key = re.compile(r"\b[A-Z]{14}-[A-Z]{10}-[A-Z]\b")
    uniprot = re.compile(r"\b(?:[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9][A-Z][A-Z0-9]{2}[0-9])\b")
    forbidden_header_tokens = {
        "query_id",
        "target_id",
        "relation_id",
        "pair_id",
        "pmid",
        "inchikey",
        "smiles",
        "sequence",
        "uniprot_accession",
    }
    for path in paths:
        text = path.read_text(encoding="utf-8")
        header = text.splitlines()[0].lower().split("\t")
        for field in header:
            if field in forbidden_header_tokens or field.endswith("_id"):
                raise ValueError(f"Forbidden identifier field in aggregate output {path.name}: {field}")
        if absolute_windows.search(text) or unc_path.search(text) or absolute_posix.search(text):
            raise ValueError(f"Absolute path detected in aggregate output {path.name}")
        if inchi_key.search(text) or uniprot.search(text):
            raise ValueError(f"Entity identifier detected in aggregate output {path.name}")


def build_input_paths(project_root: Path) -> dict[str, Path]:
    run = project_root / RUN_DIRNAME
    manuscript = project_root / MANUSCRIPT_DIRNAME
    strict = project_root / STRICT_DIRNAME
    return {
        "rebuild": run / "audit/asof_rebuild_summary.json",
        "endpoint_gate_summary": project_root / "results/strict_temporal_future_v1_1_pmid_verified_chembl31_leakage_gate_summary.json",
        "before_after": run / "audit/before_after_v1/before_after_aggregate_summary.json",
        "metric_change": run / "audit/before_after_v1/metric_change_aggregate.tsv",
        "rank_change": run / "audit/before_after_v1/rank_change_aggregate_by_baseline.tsv",
        "metrics": run / "evaluation/corrective_aggregate_metrics.tsv.gz",
        "baseline_bootstrap": run / "evaluation/corrective_baseline_bootstrap_metrics.tsv.gz",
        "contrasts": run / "evaluation/corrective_paired_bootstrap_contrasts.tsv.gz",
        "denominators": run / "evaluation/corrective_scope_denominator_audit.tsv.gz",
        "evaluation_manifest": run / "evaluation/corrective_evaluation_manifest.json",
        "evaluation_input_manifest": run / "evaluation_inputs/author_run_input_manifest.json",
        "protocol_receipt": run / "metadata/author_run_protocol_receipt.json",
        "scoring_input_manifest": run / "scoring_inputs/author_run_input_manifest.json",
        "score_manifest": run / "score/corrective_score_manifest.json",
        "complete_rank": run / "score/corrective_prediction_ranks.tsv.gz",
        "ties": run / "audit/score_degeneracy_ties_v1/rank_score_aggregate_by_baseline.tsv",
        "ties_summary": run / "audit/score_degeneracy_ties_v1/rank_score_aggregate_summary.json",
        "ties_manifest": run / "audit/score_degeneracy_ties_v1/run_manifest.json",
        "tie_runtime": run / "runtime/score_degeneracy_ties.json",
        "top100_fidelity": run / "audit/top100_exhaustive_v1/top100_exhaustive_fidelity.tsv",
        "top100_metrics": run / "audit/top100_exhaustive_v1/top100_exhaustive_metric_differences.tsv",
        "top100_summary": run / "audit/top100_exhaustive_v1/top100_exhaustive_aggregate_summary.json",
        "top100_manifest": run / "audit/top100_exhaustive_v1/run_manifest.json",
        "top100_failed_runtime": run / "runtime/top100_exhaustive.json",
        "top100_retry_runtime": run / "runtime/top100_exhaustive_retry1.json",
        "source_cohort": run / "audit/source_concentration_v1/source_concentration_by_cohort.tsv",
        "source_overlap": run / "audit/source_concentration_v1/cross_cohort_source_overlap.tsv",
        "source_date_attrition": run / "audit/source_concentration_v1/row_date_precision_attrition.tsv",
        "source_summary": run / "audit/source_concentration_v1/source_concentration_aggregate_summary.json",
        "source_manifest": run / "audit/source_concentration_v1/run_manifest.json",
        "component_scopes": run / "audit/document_component_bootstrap_v1/document_component_scope_summary.tsv",
        "component_metrics": run / "audit/document_component_bootstrap_v1/document_component_bootstrap_metrics.tsv",
        "component_focus": run / "audit/document_component_bootstrap_v1/document_component_bootstrap_focus_contrast.tsv",
        "component_summary": run / "audit/document_component_bootstrap_v1/document_component_bootstrap_summary.json",
        "component_manifest": run / "audit/document_component_bootstrap_v1/run_manifest.json",
        "component_runtime": run / "runtime/document_component_bootstrap_v1_4.json",
        "unresolved": run / "audit/entity_unresolved_v1/unresolved_reason_counts.tsv",
        "unresolved_summary": run / "audit/entity_unresolved_v1/unresolved_reason_aggregate_summary.json",
        "unresolved_manifest": run / "audit/entity_unresolved_v1/run_manifest.json",
        "endpoint_attribution": run / "audit/endpoint_score_attribution_posthoc_v1/endpoint_score_attribution_aggregate.tsv",
        "endpoint_attribution_summary": run / "audit/endpoint_score_attribution_posthoc_v1/endpoint_score_attribution_summary.json",
        "endpoint_attribution_manifest": run / "audit/endpoint_score_attribution_posthoc_v1/run_manifest.json",
        "endpoint_attribution_runtime": run / "runtime/endpoint_score_attribution_posthoc_v1.json",
        "protocol_lock": strict / "manifests/protocol_lock_manifest_v1.json",
        "implementation_base_lock": strict / "manifests/implementation_code_lock_v1.json",
        "implementation_lock": strict / "manifests/implementation_code_lock_v1_1.json",
        "release_local_implementation_lock": project_root / "manifests/implementation_code_lock_v1_1.json",
        "audit_implementation_lock": strict / "manifests/audit_implementation_code_lock_v1_2.json",
        "document_component_base_lock": strict / "manifests/document_component_audit_implementation_code_lock_v1_3.json",
        "document_component_lock": strict / "manifests/document_component_audit_implementation_code_lock_v1_4.json",
        "endpoint_attribution_implementation_lock": strict / "manifests/endpoint_score_attribution_implementation_lock_v1.json",
        "endpoint_attribution_governance": strict / "governance/endpoint_score_attribution_posthoc_amendment_20260728_08.md",
        "asof_successor_common_script": strict / "scripts/asof_successor_common.py",
        "evaluate_script": strict / "scripts/evaluate_asof_cutoff_successor.py",
        "run_author_script": strict / "scripts/run_author_run_asof_cutoff_successor.py",
        "score_script": strict / "scripts/score_asof_cutoff_successor.py",
        "pu_metrics_script": project_root / "scripts/pu_retrieval_metrics.py",
        "audit_rank_script": strict / "audit_suite_v1_20260728/scripts/audit_rank_score_structure.py",
        "audit_top100_script": strict / "audit_suite_v1_20260728/scripts/audit_top100_exhaustive_fidelity.py",
        "audit_source_script": strict / "audit_suite_v1_20260728/scripts/audit_source_concentration.py",
        "audit_unresolved_script": strict / "audit_suite_v1_20260728/scripts/audit_unresolved_reasons.py",
        "audit_runtime_wrapper": strict / "audit_suite_v1_20260728/scripts/run_with_runtime_receipt.py",
        "endpoint_attribution_script": strict / "audit_suite_v1_20260728/scripts/audit_endpoint_score_attribution_posthoc_v1.py",
        "release_gate": manuscript / "plan/release_gate_state_v1.json",
    }


def validate_lock_and_execution_chain(
    paths: dict[str, Path],
    protocol_id: str,
    protocol_lock_sha256: str,
) -> dict[str, Any]:
    """Fail closed on the locked code, audit, runtime, and evaluation custody chain."""

    checks: dict[str, bool] = {}

    def require_link(label: str, observed: Any, expected: Any) -> None:
        if observed != expected:
            raise ValueError(f"Custody link mismatch ({label}): observed={observed!r} expected={expected!r}")
        checks[label] = True

    protocol_lock = load_json(paths["protocol_lock"])
    implementation_base = load_json(paths["implementation_base_lock"])
    implementation_lock = load_json(paths["implementation_lock"])
    audit_lock = load_json(paths["audit_implementation_lock"])
    component_base_lock = load_json(paths["document_component_base_lock"])
    component_lock = load_json(paths["document_component_lock"])
    attribution_lock = load_json(paths["endpoint_attribution_implementation_lock"])
    component_manifest = load_json(paths["component_manifest"])
    attribution_manifest = load_json(paths["endpoint_attribution_manifest"])
    evaluation_manifest = load_json(paths["evaluation_manifest"])
    score_manifest = load_json(paths["score_manifest"])
    protocol_receipt = load_json(paths["protocol_receipt"])
    evaluation_input_manifest = load_json(paths["evaluation_input_manifest"])
    scoring_input_manifest = load_json(paths["scoring_input_manifest"])

    require_link("protocol_lock_protocol_id", protocol_lock.get("protocol_id"), protocol_id)
    require_link("v1_to_protocol", implementation_base["protocol_lock"]["sha256"], protocol_lock_sha256)
    require_link("v1_1_to_v1", implementation_lock["base_implementation_lock"]["sha256"], sha256_file(paths["implementation_base_lock"]))
    require_link("audit_v1_2_to_v1_1", audit_lock["base_lock"]["sha256"], sha256_file(paths["implementation_lock"]))
    require_link("document_v1_3_to_audit_v1_2", component_base_lock["base_lock"]["sha256"], sha256_file(paths["audit_implementation_lock"]))
    require_link("document_v1_4_to_v1_3", component_lock["base_lock"]["sha256"], sha256_file(paths["document_component_base_lock"]))
    require_link(
        "component_manifest_to_release_local_v1_3",
        component_manifest["implementation_lock"]["sha256"],
        sha256_file(paths["document_component_base_lock"]),
    )
    require_link("component_manifest_to_protocol", component_manifest["protocol_lock_sha256"], protocol_lock_sha256)

    require_link("attribution_lock_protocol_id", attribution_lock.get("protocol_id"), protocol_id)
    require_link("attribution_lock_to_protocol", attribution_lock["protocol_lock_sha256"], protocol_lock_sha256)
    require_link("attribution_lock_to_governance", attribution_lock["governance_amendment_sha256"], sha256_file(paths["endpoint_attribution_governance"]))
    require_link("attribution_lock_to_script", attribution_lock["script_sha256"], sha256_file(paths["endpoint_attribution_script"]))
    require_link("attribution_manifest_to_lock", attribution_manifest["implementation_lock"]["sha256"], sha256_file(paths["endpoint_attribution_implementation_lock"]))
    require_link("attribution_manifest_to_governance", attribution_manifest["governance"]["sha256"], sha256_file(paths["endpoint_attribution_governance"]))
    require_link("attribution_manifest_to_script", attribution_manifest["script"]["sha256"], sha256_file(paths["endpoint_attribution_script"]))

    audit_script_links = (
        ("score_tie_script_to_base_lock", "ties_manifest", "audit_rank_script", "audit_suite_v1_20260728/scripts/audit_rank_score_structure.py"),
        ("top100_script_to_base_lock", "top100_manifest", "audit_top100_script", "audit_suite_v1_20260728/scripts/audit_top100_exhaustive_fidelity.py"),
        ("source_script_to_base_lock", "source_manifest", "audit_source_script", "audit_suite_v1_20260728/scripts/audit_source_concentration.py"),
        ("unresolved_script_to_base_lock", "unresolved_manifest", "audit_unresolved_script", "audit_suite_v1_20260728/scripts/audit_unresolved_reasons.py"),
    )
    for label, manifest_name, script_name, lock_key in audit_script_links:
        audit_manifest = load_json(paths[manifest_name])
        actual_script_sha = sha256_file(paths[script_name])
        require_link(f"{label}_manifest", audit_manifest["script"]["sha256"], actual_script_sha)
        require_link(label, implementation_base["source_sha256"][lock_key], actual_script_sha)

    for lock_name, lock in (
        ("implementation_base", implementation_base),
        ("implementation_v1_1", implementation_lock),
        ("audit_v1_2", audit_lock),
        ("document_v1_3", component_base_lock),
        ("document_v1_4", component_lock),
    ):
        require_link(f"{lock_name}_protocol_id", lock.get("protocol_id"), protocol_id)

    evaluation_manifest = verify_manifest_outputs(
        paths["evaluation_manifest"],
        [paths["metrics"], paths["baseline_bootstrap"], paths["contrasts"], paths["denominators"]],
        protocol_id,
        protocol_lock_sha256,
    )
    verify_manifest_outputs(paths["before_after"], [paths["metric_change"], paths["rank_change"]], protocol_id, protocol_lock_sha256)
    verify_manifest_outputs(paths["ties_manifest"], [paths["ties"], paths["ties_summary"]], protocol_id, protocol_lock_sha256)
    verify_manifest_outputs(paths["top100_manifest"], [paths["top100_fidelity"], paths["top100_metrics"], paths["top100_summary"]], protocol_id, protocol_lock_sha256)
    verify_manifest_outputs(paths["source_manifest"], [paths["source_cohort"], paths["source_overlap"], paths["source_date_attrition"], paths["source_summary"]], protocol_id, protocol_lock_sha256)
    verify_manifest_outputs(paths["unresolved_manifest"], [paths["unresolved"], paths["unresolved_summary"]], protocol_id, protocol_lock_sha256)
    verify_manifest_outputs(paths["endpoint_attribution_manifest"], [paths["endpoint_attribution"], paths["endpoint_attribution_summary"]], protocol_id, protocol_lock_sha256)
    component_manifest = verify_manifest_outputs(paths["component_manifest"], [paths["component_scopes"], paths["component_metrics"], paths["component_focus"], paths["component_summary"]], protocol_id, protocol_lock_sha256)

    require_link("protocol_receipt_protocol_id", protocol_receipt.get("protocol_id"), protocol_id)
    require_link("protocol_receipt_to_protocol_lock", protocol_receipt["protocol_lock"]["sha256"], protocol_lock_sha256)
    require_link(
        "protocol_receipt_to_release_local_code_lock",
        protocol_receipt["code_lock"]["sha256"],
        sha256_file(paths["release_local_implementation_lock"]),
    )
    require_link("evaluation_manifest_to_receipt", evaluation_manifest["receipt"]["sha256"], sha256_file(paths["protocol_receipt"]))
    require_link("evaluation_manifest_to_input_manifest", evaluation_manifest["evaluation_input_manifest"]["sha256"], sha256_file(paths["evaluation_input_manifest"]))
    require_link("evaluation_manifest_to_score_manifest", evaluation_manifest["score_manifest"]["sha256"], sha256_file(paths["score_manifest"]))
    require_link("evaluation_input_manifest_protocol_id", evaluation_input_manifest.get("protocol_id"), protocol_id)
    require_link("score_manifest_protocol_id", score_manifest.get("protocol_id"), protocol_id)
    require_link("score_manifest_to_receipt", score_manifest["receipt"]["sha256"], sha256_file(paths["protocol_receipt"]))
    require_link("score_manifest_to_input_manifest", score_manifest["input_manifest"]["sha256"], sha256_file(paths["scoring_input_manifest"]))
    require_link("scoring_input_manifest_protocol_id", scoring_input_manifest.get("protocol_id"), protocol_id)
    complete_rank_sha256 = sha256_file(paths["complete_rank"])
    require_link("complete_rank_to_score_manifest", complete_rank_sha256, score_manifest["rank_output"]["sha256"])
    require_link("evaluation_rank_to_score_rank", evaluation_manifest["score_rank"]["sha256"], score_manifest["rank_output"]["sha256"])

    score_runtime = score_manifest.get("runtime", {})
    score_wall_seconds = score_runtime.get("wall_seconds_by_phase", {}).get("total_before_manifest_write")
    score_peak_rss_bytes = score_runtime.get("process_peak_rss_bytes")
    if type(score_wall_seconds) not in (int, float) or score_wall_seconds <= 0:
        raise ValueError("Score manifest lacks a positive total scoring wall time")
    checks["score_runtime_wall_seconds_positive"] = True
    if type(score_peak_rss_bytes) is not int or score_peak_rss_bytes <= 0:
        raise ValueError("Score manifest lacks a positive measured peak RSS")
    checks["score_runtime_peak_rss_positive"] = True
    if not isinstance(score_runtime.get("process_peak_rss_method"), str) or not score_runtime["process_peak_rss_method"].strip():
        raise ValueError("Score manifest lacks a peak-RSS measurement method")
    checks["score_runtime_peak_rss_method_recorded"] = True

    ties_manifest = load_json(paths["ties_manifest"])
    tie_rank_inputs = [
        item
        for item in ties_manifest.get("inputs", [])
        if isinstance(item, dict)
        and item.get("role") == "corrective_full_rank_file"
        and item.get("basename") == paths["complete_rank"].name
    ]
    if len(tie_rank_inputs) != 1:
        raise ValueError("Tie-audit manifest must declare exactly one complete-rank input")
    require_link("tie_manifest_to_complete_rank", tie_rank_inputs[0].get("sha256"), complete_rank_sha256)

    evaluation_code_links = (
        ("asof_successor_common.py", "asof_successor_common_script", "scripts/asof_successor_common.py"),
        ("evaluate_asof_cutoff_successor.py", "evaluate_script", "scripts/evaluate_asof_cutoff_successor.py"),
        ("pu_retrieval_metrics.py", "pu_metrics_script", "workspace/scripts/pu_retrieval_metrics.py"),
    )
    for basename, script_name, lock_key in evaluation_code_links:
        actual_script_sha = sha256_file(paths[script_name])
        require_link(f"evaluation_{basename}_to_actual", evaluation_manifest["code"][basename]["sha256"], actual_script_sha)

    scoring_code_links = (
        ("asof_successor_common.py", "asof_successor_common_script", "scripts/asof_successor_common.py"),
        ("score_asof_cutoff_successor.py", "score_script", "scripts/score_asof_cutoff_successor.py"),
        ("pu_retrieval_metrics.py", "pu_metrics_script", "workspace/scripts/pu_retrieval_metrics.py"),
    )
    for basename, script_name, lock_key in scoring_code_links:
        actual_script_sha = sha256_file(paths[script_name])
        require_link(f"score_{basename}_to_actual", score_manifest["code"][basename]["sha256"], actual_script_sha)

    for input_name, input_item in evaluation_manifest["evaluation_inputs"].items():
        require_link(
            f"evaluation_input_{input_name}_to_input_manifest",
            input_item["sha256"],
            evaluation_input_manifest["file_sha256"]["endpoint" if input_name == "endpoint" else input_name],
        )

    wrapper_sha = sha256_file(paths["audit_runtime_wrapper"])
    require_link("runtime_wrapper_to_base_lock", implementation_base["source_sha256"]["audit_suite_v1_20260728/scripts/run_with_runtime_receipt.py"], wrapper_sha)
    component_runtime = load_json(paths["component_runtime"])
    endpoint_runtime = load_json(paths["endpoint_attribution_runtime"])
    tie_runtime = load_json(paths["tie_runtime"])
    top100_failed_runtime = load_json(paths["top100_failed_runtime"])
    top100_retry_runtime = load_json(paths["top100_retry_runtime"])

    def validate_runtime_receipt(
        label: str,
        runtime: dict[str, Any],
        expected_audit: str,
        expected_success: bool,
    ) -> None:
        exit_code = runtime.get("exit_code")
        if expected_success:
            require_link(f"{label}_exit_code", exit_code, 0)
        else:
            if type(exit_code) is not int or exit_code == 0:
                raise ValueError(f"Custody link mismatch ({label}_nonzero_exit_code): observed={exit_code!r}")
            checks[f"{label}_nonzero_exit_code"] = True
        require_link(f"{label}_audit_label", runtime.get("audit_label"), expected_audit)
        require_link(f"{label}_protocol_id", runtime.get("protocol_id"), protocol_id)
        require_link(f"{label}_protocol_lock", runtime.get("protocol_lock", {}).get("sha256"), protocol_lock_sha256)
        require_link(f"{label}_wrapper", runtime.get("wrapper_script", {}).get("sha256"), wrapper_sha)
        wall_seconds = runtime.get("wall_seconds")
        if type(wall_seconds) not in (int, float) or wall_seconds <= 0:
            raise ValueError(f"Runtime receipt {label} lacks a positive wall time")
        checks[f"{label}_wall_seconds_positive"] = True
        argv_sha256 = runtime.get("command", {}).get("argv_sha256")
        if not isinstance(argv_sha256, str) or not re.fullmatch(r"[0-9a-f]{64}", argv_sha256):
            raise ValueError(f"Runtime receipt {label} lacks a valid command fingerprint")
        checks[f"{label}_command_fingerprint_recorded"] = True

    for label, runtime, expected_audit, expected_success in (
        (
            "component_runtime",
            component_runtime,
            "document_component_bootstrap_v1_release_local_v1_3",
            True,
        ),
        ("endpoint_attribution_runtime", endpoint_runtime, "endpoint_score_attribution_posthoc_v1", True),
        ("tie_runtime", tie_runtime, "score_degeneracy_ties", True),
        ("top100_failed_runtime", top100_failed_runtime, "top100_exhaustive", False),
        ("top100_retry_runtime", top100_retry_runtime, "top100_exhaustive_retry1", True),
    ):
        validate_runtime_receipt(label, runtime, expected_audit, expected_success)

    require_link(
        "top100_failed_retry_same_command",
        top100_failed_runtime["command"]["argv_sha256"],
        top100_retry_runtime["command"]["argv_sha256"],
    )
    try:
        failed_end = datetime.fromisoformat(top100_failed_runtime["ended_at_utc"])
        retry_start = datetime.fromisoformat(top100_retry_runtime["started_at_utc"])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("Top-100 runtime receipts lack parseable attempt timestamps") from exc
    if failed_end >= retry_start:
        raise ValueError("Top-100 retry did not start after the failed attempt ended")
    checks["top100_failed_before_retry"] = True

    def require_manifest_created_during_runtime(
        label: str,
        manifest: dict[str, Any],
        runtime: dict[str, Any],
    ) -> None:
        try:
            created = datetime.fromisoformat(manifest["created_at_utc"])
            started = datetime.fromisoformat(runtime["started_at_utc"])
            ended = datetime.fromisoformat(runtime["ended_at_utc"])
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"{label} lacks parseable custody timestamps") from exc
        if not started <= created <= ended:
            raise ValueError(f"{label} manifest timestamp falls outside its runtime receipt")
        checks[f"{label}_manifest_created_during_runtime"] = True

    require_manifest_created_during_runtime("tie_runtime", ties_manifest, tie_runtime)
    require_manifest_created_during_runtime(
        "top100_retry_runtime",
        load_json(paths["top100_manifest"]),
        top100_retry_runtime,
    )

    require_link("attribution_lock_state", attribution_lock.get("lock_state"), "LOCKED_BEFORE_REAL_ENDPOINT_SCORE_ATTRIBUTION_EXECUTION")
    require_link("attribution_outcome_visible_disclosure", attribution_lock.get("outcome_visible_before_lock"), True)

    return {
        "checks": checks,
        "evaluation_manifest": evaluation_manifest,
        "score_manifest": score_manifest,
        "component_manifest": component_manifest,
        "component_runtime": component_runtime,
        "endpoint_attribution_runtime": endpoint_runtime,
        "tie_runtime": tie_runtime,
        "top100_failed_runtime": top100_failed_runtime,
        "top100_retry_runtime": top100_retry_runtime,
    }


LOCAL_STAGING_GATE_FIELD = "local_restricted_staging_authorized"
EXTERNAL_TRANSFER_PREREQUISITE_FIELDS = (
    "external_transfer_authorized",
    "rights_clearance_complete",
    "recipient_and_channel_configured",
    "encryption_configured",
)
EXTERNAL_TRANSFER_READY_FIELD = "ready_for_external_transfer"
RELEASE_GATE_FIELDS = (
    LOCAL_STAGING_GATE_FIELD,
    *EXTERNAL_TRANSFER_PREREQUISITE_FIELDS,
    EXTERNAL_TRANSFER_READY_FIELD,
)


def validate_release_gate(release_gate: dict[str, Any]) -> tuple[str, ...]:
    for field in RELEASE_GATE_FIELDS:
        if type(release_gate.get(field)) is not bool:
            raise ValueError(f"Release gate {field} must be a JSON boolean")
    computed_ready = all(release_gate[field] for field in EXTERNAL_TRANSFER_PREREQUISITE_FIELDS)
    if release_gate[EXTERNAL_TRANSFER_READY_FIELD] != computed_ready:
        raise ValueError("Release gate ready flag does not equal the required conjunction")
    return RELEASE_GATE_FIELDS


def build_tables(project_root: Path) -> dict[str, Any]:
    run = project_root / RUN_DIRNAME
    manuscript = project_root / MANUSCRIPT_DIRNAME
    strict = project_root / STRICT_DIRNAME
    output_dir = manuscript / "tables"
    output_dir.mkdir(parents=True, exist_ok=True)

    inputs: dict[str, dict[str, Any]] = {}
    outputs: dict[str, dict[str, Any]] = {}

    paths = build_input_paths(project_root)
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing locked table inputs:\n" + "\n".join(missing))
    for path in paths.values():
        record_input(path, project_root, inputs)

    rebuild = load_json(paths["rebuild"])
    endpoint_gate = load_json(paths["endpoint_gate_summary"])
    before_after = load_json(paths["before_after"])
    unresolved_summary = load_json(paths["unresolved_summary"])
    protocol_lock_sha256 = sha256_file(paths["protocol_lock"])
    if protocol_lock_sha256 != "96befee13ae1d41ad433c8697fac92ccd30fb25e24c3cf1279c6b4b7e040abd9":
        raise ValueError("Unexpected corrective protocol-lock hash")

    if not all(rebuild.get("hard_gate_checks", {}).values()):
        raise ValueError("At least one rebuild hard gate is false")

    protocol_id = str(rebuild["protocol_id"])
    custody = validate_lock_and_execution_chain(paths, protocol_id, protocol_lock_sha256)
    evaluation_manifest = custody["evaluation_manifest"]
    score_manifest = custody["score_manifest"]
    component_manifest = custody["component_manifest"]
    component_runtime = custody["component_runtime"]
    endpoint_attribution_runtime = custody["endpoint_attribution_runtime"]
    tie_runtime = custody["tie_runtime"]
    top100_failed_runtime = custody["top100_failed_runtime"]
    top100_retry_runtime = custody["top100_retry_runtime"]
    release_gate = load_json(paths["release_gate"])
    release_fields = validate_release_gate(release_gate)
    custody["checks"]["release_gate_fields_strict_json_boolean"] = True
    custody["checks"]["release_gate_ready_rule_consistent"] = True
    custody["checks"]["local_restricted_staging_gate_recorded"] = True

    expected_scopes = {
        "temporal_strict_ab",
        "scaffold_cold_strict_ab",
        "double_cold_0_30",
        "double_cold_0_50",
        "double_cold_0_70",
    }
    expected_baselines = {
        "weighted_target_popularity",
        "sequence_3mer_transfer",
        "weighted_morgan_transfer",
        "structure_sequence_pair_neighbor",
    }
    expected_metrics = {"Recall@10", "Recall@50", "NDCG@10", "NDCG@50", "MRR"}

    # T1: row-level cutoff repair, scope preservation, and frozen key-set hashes.
    t1_rows: list[dict[str, Any]] = []
    t1_rows.extend(key_value_rows("row_eligibility", rebuild["row_eligibility_counts"]))
    endpoint_status = endpoint_gate["status_counts_in_frozen_future_table"]
    endpoint_flow = {
        "initial_later_recorded_candidate_pairs": endpoint_gate["inputs"]["frozen_future_pair_count"],
        "excluded_chembl31_historical_activity_pairs": endpoint_status["historical_activity_recorded_in_chembl31"],
        "excluded_entity_unresolved_pairs": endpoint_status["entity_pair_not_sqlite_validated"],
        "final_frozen_endpoint_relations": endpoint_status["no_activity_found_in_validated_chembl31_entity_pair"],
    }
    if (
        endpoint_flow["initial_later_recorded_candidate_pairs"]
        - endpoint_flow["excluded_chembl31_historical_activity_pairs"]
        - endpoint_flow["excluded_entity_unresolved_pairs"]
        != endpoint_flow["final_frozen_endpoint_relations"]
    ):
        raise ValueError("Endpoint construction flow does not close arithmetically")
    if endpoint_flow["final_frozen_endpoint_relations"] != rebuild["counts"]["endpoint_relations"]:
        raise ValueError("Endpoint construction flow disagrees with corrective rebuild")
    t1_rows.extend(key_value_rows("endpoint_construction", endpoint_flow))
    t1_rows.extend(key_value_rows("frozen_counts", rebuild["counts"]))
    for scope, counts in sorted(rebuild["scope_counts"].items()):
        for item, value in sorted(counts.items()):
            t1_rows.append({"section": f"scope:{scope}", "item": item, "value": value, "status": "verified"})
    t1_rows.extend(key_value_rows("hard_gate", rebuild["hard_gate_checks"]))
    t1_rows.extend(key_value_rows("membership_sha256", rebuild["membership_sha256"]))
    t1_path = output_dir / "Table_1_temporal_repair_flow.tsv"
    atomic_write_tsv(t1_path, ["section", "item", "value", "status"], t1_rows)

    # T2: historical before/after audit, including provenance concentration.
    source_rows = read_tsv(paths["source_cohort"])
    require_columns(
        source_rows,
        {
            "cohort", "relation_count", "query_count", "unique_source_document_count",
            "relation_source_document_edge_count", "relation_count_without_source_document",
            "largest_source_document_relation_count", "largest_source_document_relation_fraction",
            "top10_source_document_edge_share", "source_document_edge_hhi",
            "source_document_relation_degree_gini", "relation_source_component_count",
            "largest_relation_component_fraction", "query_source_component_count",
            "largest_query_component_fraction",
        },
        "source concentration",
    )
    require_exact_set({row["cohort"] for row in source_rows}, {"historical", "endpoint"}, "source cohorts")
    if len(source_rows) != 2:
        raise ValueError("Source concentration must contain exactly two cohort rows")
    historical_sources = next(row for row in source_rows if row["cohort"] == "historical")
    rank_change_rows = read_tsv(paths["rank_change"])
    require_columns(
        rank_change_rows,
        {
            "baseline", "row_count", "score_changed_rows", "rank_changed_rows",
            "absolute_rank_change_sum", "maximum_absolute_rank_change", "query_blocks",
            "query_blocks_with_any_rank_change", "query_blocks_with_top10_membership_change",
            "query_blocks_with_top50_membership_change", "top10_symmetric_difference_sum",
            "top50_symmetric_difference_sum", "score_changed_row_fraction",
            "rank_changed_row_fraction", "mean_absolute_rank_change",
        },
        "before/after rank changes",
    )
    require_exact_set({row["baseline"] for row in rank_change_rows}, expected_baselines, "before/after baselines")
    metric_change_rows = read_tsv(paths["metric_change"])
    require_columns(metric_change_rows, {"scope", "baseline", "metric", "old_value", "corrected_value", "corrected_minus_old"}, "before/after metric changes")
    require_unique_keys(
        metric_change_rows,
        ("scope", "baseline", "metric"),
        set(itertools.product(expected_scopes, expected_baselines, expected_metrics)),
        "before/after metric grid",
    )
    date_attrition_rows = read_tsv(paths["source_date_attrition"])
    require_columns(date_attrition_rows, {"cohort", "row_status", "date_source", "date_precision", "source_row_count", "distinct_relation_count"}, "source-date attrition")
    t2_rows: list[dict[str, Any]] = []
    historical_items = {
        "historical_strict_v2_rows": rebuild["counts"]["historical_strict_v2_rows"],
        "eligible_pre_cutoff_rows": rebuild["row_eligibility_counts"]["eligible_pre_cutoff"],
        "excluded_non_day_precision_rows": rebuild["row_eligibility_counts"]["excluded_non_day_precision"],
        "excluded_not_numeric_pmid_rows": rebuild["row_eligibility_counts"]["excluded_not_numeric_pmid"],
        "historical_pairs_with_any_excluded_row": rebuild["counts"]["historical_pairs_with_any_excluded_row"],
        "historical_tier_weight_changes": rebuild["counts"]["historical_tier_weight_changes"],
        "historical_inchi_repairs": rebuild["counts"]["historical_inchi_repairs"],
        "historical_smiles_changes": rebuild["counts"]["historical_smiles_changed_vs_old_combined"],
        "historical_morgan_changes": rebuild["counts"]["historical_morgan_changed_vs_old_combined"],
        "historical_scaffold_assignment_changes": rebuild["counts"]["historical_scaffold_assignment_changes"],
        "metric_cells_changed": before_after["metric_cells_changed"],
        "metric_delta_cells": before_after["metric_delta_cells"],
        "maximum_absolute_metric_delta": before_after["maximum_absolute_metric_delta"],
        "historical_unique_pmids": historical_sources["unique_source_document_count"],
        "historical_largest_pmid_relation_fraction": historical_sources["largest_source_document_relation_fraction"],
        "historical_largest_query_pmid_component_fraction": historical_sources["largest_query_component_fraction"],
    }
    for item, value in historical_items.items():
        t2_rows.append({"category": "repair_summary", "scope": "", "baseline": "", "metric": "", "item": item, "value": value, "old_value": "", "corrected_value": "", "corrected_minus_old": "", "status": "verified"})
    for row in rank_change_rows:
        for item in (
            "score_changed_row_fraction",
            "rank_changed_row_fraction",
            "mean_absolute_rank_change",
            "maximum_absolute_rank_change",
            "query_blocks_with_top10_membership_change",
            "query_blocks_with_top50_membership_change",
        ):
            t2_rows.append(
                {
                    "category": "rank_change_by_baseline",
                    "scope": "",
                    "item": item,
                    "baseline": row["baseline"],
                    "metric": "",
                    "value": row[item],
                    "old_value": "",
                    "corrected_value": "",
                    "corrected_minus_old": "",
                    "status": "verified",
                }
            )
    for row in metric_change_rows:
        t2_rows.append(
            {
                "category": "metric_before_after",
                "scope": row["scope"],
                "baseline": row["baseline"],
                "metric": row["metric"],
                "item": "",
                "value": "",
                "old_value": row["old_value"],
                "corrected_value": row["corrected_value"],
                "corrected_minus_old": row["corrected_minus_old"],
                "status": "verified",
            }
        )
    for row in date_attrition_rows:
        t2_rows.append(
            {
                "category": "source_date_precision_attrition",
                "scope": row["cohort"],
                "baseline": "",
                "metric": "",
                "item": f"{row['row_status']}|{row['date_source']}|{row['date_precision']}",
                "value": row["source_row_count"],
                "old_value": "",
                "corrected_value": "",
                "corrected_minus_old": "",
                "status": f"distinct_relations={row['distinct_relation_count']}",
            }
        )
    t2_path = output_dir / "Table_2_historical_before_after_audit.tsv"
    atomic_write_tsv(t2_path, ["category", "scope", "baseline", "metric", "item", "value", "old_value", "corrected_value", "corrected_minus_old", "status"], t2_rows)

    # T3: exact corrected aggregate metrics; no zero cell may be removed.
    metric_rows = read_tsv(paths["metrics"])
    metric_columns = {
        "scope", "baseline", "evaluable_query_count", "candidate_relation_count",
        "candidate_target_median", "candidate_target_iqr", "candidate_target_min",
        "candidate_target_max", "Recall@10", "Recall@50", "NDCG@10", "NDCG@50",
        "MRR", "zero_recall_at_10_queries", "zero_recall_at_50_queries",
        "zero_mrr_queries", "status",
    }
    require_columns(metric_rows, metric_columns, "T3")
    require_unique_keys(metric_rows, ("scope", "baseline"), set(itertools.product(expected_scopes, expected_baselines)), "T3")
    t3_path = output_dir / "Table_3_corrected_aggregate_performance.tsv"
    t3_fields = [
        "scope", "baseline", "evaluable_query_count", "candidate_relation_count",
        "candidate_target_median", "candidate_target_iqr", "candidate_target_min",
        "candidate_target_max", "Recall@10", "Recall@50", "NDCG@10", "NDCG@50",
        "MRR", "zero_recall_at_10_queries", "zero_recall_at_50_queries",
        "zero_mrr_queries", "status",
    ]
    atomic_write_tsv(t3_path, t3_fields, metric_rows)

    # T4: normalized query-bootstrap baseline cells and all paired contrasts.
    bootstrap = evaluation_manifest["bootstrap"]
    baseline_rows = read_tsv(paths["baseline_bootstrap"])
    contrast_rows = read_tsv(paths["contrasts"])
    require_columns(baseline_rows, {"scope", "baseline", "metric", "query_count", "mean", "ci95_low", "ci95_high", "status"}, "T4 baseline bootstrap")
    require_columns(contrast_rows, {"scope", "left_baseline", "right_baseline", "comparison_role", "metric", "query_count", "mean_difference_left_minus_right", "ci95_low", "ci95_high", "status"}, "T4 paired contrasts")
    require_unique_keys(
        baseline_rows,
        ("scope", "baseline", "metric"),
        set(itertools.product(expected_scopes, expected_baselines, expected_metrics)),
        "T4 baseline bootstrap",
    )
    contrast_keys = [(row["scope"], frozenset((row["left_baseline"], row["right_baseline"])), row["metric"]) for row in contrast_rows]
    expected_pairs = {frozenset(pair) for pair in itertools.combinations(expected_baselines, 2)}
    expected_contrast_keys = set(itertools.product(expected_scopes, expected_pairs, expected_metrics))
    if len(contrast_keys) != len(set(contrast_keys)) or set(contrast_keys) != expected_contrast_keys:
        raise ValueError("T4 paired-contrast grid is not the complete unique 5 x 6 x 5 design")
    focus_rows = [row for row in contrast_rows if row["comparison_role"] == "prespecified_focus"]
    expected_focus = set(itertools.product(expected_scopes, expected_metrics))
    if {(row["scope"], row["metric"]) for row in focus_rows} != expected_focus or len(focus_rows) != 25:
        raise ValueError("T4 focus contrast is incomplete")
    if any(row["left_baseline"] != "structure_sequence_pair_neighbor" or row["right_baseline"] != "weighted_morgan_transfer" for row in focus_rows):
        raise ValueError("T4 focus contrast direction changed")
    t4_rows: list[dict[str, Any]] = []
    for row in baseline_rows:
        t4_rows.append(
            {
                "record_type": "baseline",
                "scope": row["scope"],
                "baseline_or_left": row["baseline"],
                "right_baseline": "",
                "comparison_role": "",
                "metric": row["metric"],
                "query_count": row["query_count"],
                "estimate": row["mean"],
                "ci95_low": row["ci95_low"],
                "ci95_high": row["ci95_high"],
                "status": row["status"],
                "replicates": bootstrap["replicates"],
                "seed": bootstrap["seed"],
                "resampling_unit": bootstrap["unit"],
            }
        )
    for row in contrast_rows:
        t4_rows.append(
            {
                "record_type": "paired_contrast_left_minus_right",
                "scope": row["scope"],
                "baseline_or_left": row["left_baseline"],
                "right_baseline": row["right_baseline"],
                "comparison_role": row["comparison_role"],
                "metric": row["metric"],
                "query_count": row["query_count"],
                "estimate": row["mean_difference_left_minus_right"],
                "ci95_low": row["ci95_low"],
                "ci95_high": row["ci95_high"],
                "status": row["status"],
                "replicates": bootstrap["replicates"],
                "seed": bootstrap["seed"],
                "resampling_unit": bootstrap["unit"],
            }
        )
    t4_path = output_dir / "Table_4_corrected_bootstrap_summaries.tsv"
    atomic_write_tsv(
        t4_path,
        [
            "record_type",
            "scope",
            "baseline_or_left",
            "right_baseline",
            "comparison_role",
            "metric",
            "query_count",
            "estimate",
            "ci95_low",
            "ci95_high",
            "status",
            "replicates",
            "seed",
            "resampling_unit",
        ],
        t4_rows,
    )

    # T5: score degeneracy and top-k boundary ties for the fixed query set.
    tie_rows = read_tsv(paths["ties"])
    require_columns(
        tie_rows,
        {
            "baseline", "query_count", "rank_row_count", "all_zero_query_count",
            "all_zero_query_fraction", "constant_score_query_count", "query_count_with_any_tie",
            "positive_rank_row_count", "positive_rank_row_fraction",
            "top10_boundary_tie_query_count", "top10_zero_score_boundary_query_count",
            "top50_boundary_tie_query_count", "top50_zero_score_boundary_query_count",
        },
        "T5",
    )
    require_exact_set({row["baseline"] for row in tie_rows}, expected_baselines, "T5 baselines")
    if len(tie_rows) != 4:
        raise ValueError(f"T5 row count must be 4, observed {len(tie_rows)}")
    tie_summary = load_json(paths["ties_summary"])
    summary_by_baseline = {row["baseline"]: row for row in tie_summary["baselines"]}
    require_exact_set(set(summary_by_baseline), expected_baselines, "T5 summary baselines")
    t5_rows: list[dict[str, Any]] = []
    for row in tie_rows:
        summary = summary_by_baseline[row["baseline"]]
        enriched: dict[str, Any] = dict(row)
        for prefix, distribution in (
            ("unique_score_count", summary["unique_score_count_distribution"]),
            ("largest_tie_block_size", summary["largest_tie_block_size_distribution"]),
        ):
            for statistic in ("min", "q25", "median", "q75", "max"):
                enriched[f"{prefix}_{statistic}"] = distribution[statistic]
        for boundary in ("10", "50"):
            audit = summary["boundary_audits"][boundary]
            distribution = audit["boundary_tie_block_size_distribution"]
            for statistic in ("min", "median", "max"):
                enriched[f"top{boundary}_boundary_block_{statistic}"] = distribution[statistic]
            enriched[f"top{boundary}_tie_members_selected_by_salt"] = audit["tie_members_selected_by_salt_total"]
            enriched[f"top{boundary}_tie_members_excluded_by_salt"] = audit["tie_members_excluded_by_salt_total"]
        t5_rows.append(enriched)
    t5_path = output_dir / "Table_5_score_degeneracy_and_ties.tsv"
    t5_fields = list(tie_rows[0].keys()) + [
        "unique_score_count_min", "unique_score_count_q25", "unique_score_count_median",
        "unique_score_count_q75", "unique_score_count_max", "largest_tie_block_size_min",
        "largest_tie_block_size_q25", "largest_tie_block_size_median",
        "largest_tie_block_size_q75", "largest_tie_block_size_max",
        "top10_boundary_block_min", "top10_boundary_block_median", "top10_boundary_block_max",
        "top10_tie_members_selected_by_salt", "top10_tie_members_excluded_by_salt",
        "top50_boundary_block_min", "top50_boundary_block_median", "top50_boundary_block_max",
        "top50_tie_members_selected_by_salt", "top50_tie_members_excluded_by_salt",
    ]
    atomic_write_tsv(t5_path, t5_fields, t5_rows)

    # S1: scope denominators plus execution-integrity gates and immutable hashes.
    denominator_rows = read_tsv(paths["denominators"])
    require_columns(denominator_rows, {"scope", "candidate_relation_count", "query_compound_count", "target_count", "A_affinity_candidate_count", "B_quantitative_functional_candidate_count", "status"}, "S1 denominators")
    require_unique_keys(denominator_rows, ("scope",), {(scope,) for scope in expected_scopes}, "S1 denominators")
    s1_rows: list[dict[str, Any]] = []
    for row in denominator_rows:
        for key, value in row.items():
            if key == "scope":
                continue
            s1_rows.append(
                {
                    "record_type": "scope_denominator",
                    "scope_or_group": row["scope"],
                    "item": key,
                    "value": value,
                    "status": row["status"],
                }
            )
    for key, value in sorted(rebuild["hard_gate_checks"].items()):
        s1_rows.append({"record_type": "hard_gate", "scope_or_group": "global", "item": key, "value": value, "status": "verified"})
    for key, value in sorted(rebuild["membership_sha256"].items()):
        s1_rows.append({"record_type": "membership_sha256", "scope_or_group": "global", "item": key, "value": value, "status": "verified"})
    s1_path = output_dir / "Table_S1_scope_mask_integrity.tsv"
    atomic_write_tsv(s1_path, ["record_type", "scope_or_group", "item", "value", "status"], s1_rows)

    # S2: fixed top-100 primary representation versus exhaustive sensitivity.
    fidelity_rows = read_tsv(paths["top100_fidelity"])
    difference_rows = read_tsv(paths["top100_metrics"])
    require_columns(
        fidelity_rows,
        {
            "query_count", "eligible_query_target_count", "query_count_with_any_tolerance_score_change",
            "tolerance_score_changed_target_count", "tolerance_score_changed_target_fraction",
            "query_count_with_any_rank_change", "rank_changed_target_count", "rank_changed_target_fraction",
            "top10_changed_membership_query_count", "top50_changed_membership_query_count",
            "sequence_top100_boundary_tie_target_count",
        },
        "S2 fidelity",
    )
    require_columns(difference_rows, {"scope", "metric", "relation_count", "query_count", "top100_value", "exhaustive_value", "exhaustive_minus_top100"}, "S2 metric differences")
    if len(fidelity_rows) != 1:
        raise ValueError("S2 fidelity summary must contain one row")
    require_unique_keys(difference_rows, ("scope", "metric"), set(itertools.product(expected_scopes, expected_metrics)), "S2 metric differences")
    s2_rows: list[dict[str, Any]] = []
    for key, value in fidelity_rows[0].items():
        s2_rows.append(
            {
                "record_type": "fidelity_summary",
                "scope": "all_queries",
                "metric_or_item": key,
                "relation_count": "",
                "query_count": fidelity_rows[0].get("query_count", ""),
                "fidelity_value": value,
                "top100_value": "",
                "exhaustive_value": "",
                "exhaustive_minus_top100": "",
                "status": "post_hoc_sensitivity_only",
            }
        )
    top100_summary = load_json(paths["top100_summary"])
    fidelity_summary = top100_summary["rank_and_score_fidelity"]
    supplemental_fidelity_values = {
        "exact_score_changed_target_count": fidelity_summary["exact_score_changed_target_count"],
        "exact_score_changed_target_fraction": fidelity_summary["exact_score_changed_target_fraction"],
        "rank_spearman_min": fidelity_summary["rank_spearman_distribution"]["min"],
        "rank_spearman_q25": fidelity_summary["rank_spearman_distribution"]["q25"],
        "rank_spearman_median": fidelity_summary["rank_spearman_distribution"]["median"],
        "rank_spearman_q75": fidelity_summary["rank_spearman_distribution"]["q75"],
        "rank_spearman_max": fidelity_summary["rank_spearman_distribution"]["max"],
        "rank_spearman_mean": fidelity_summary["rank_spearman_distribution"]["mean"],
        "mean_absolute_rank_shift_mean": fidelity_summary["mean_absolute_rank_shift_distribution"]["mean"],
        "maximum_absolute_rank_shift_max": fidelity_summary["maximum_absolute_rank_shift_distribution"]["max"],
        "mean_absolute_score_error_mean": fidelity_summary["mean_absolute_score_error_distribution"]["mean"],
        "maximum_absolute_score_error_max": fidelity_summary["maximum_absolute_score_error_distribution"]["max"],
        "top10_jaccard_mean": fidelity_summary["top_k"]["10"]["jaccard_distribution"]["mean"],
        "top50_jaccard_mean": fidelity_summary["top_k"]["50"]["jaccard_distribution"]["mean"],
        "locked_primary_top_k": top100_summary["locked_primary_top_k"],
        "locked_score_absolute_tolerance": top100_summary["locked_score_absolute_tolerance"],
        "historical_target_count": top100_summary["historical_target_count"],
        "candidate_target_count": top100_summary["candidate_target_count"],
        "alternative_rank_ledger_written": top100_summary["alternative_rank_ledger_written"],
    }
    for key, value in supplemental_fidelity_values.items():
        s2_rows.append(
            {
                "record_type": "fidelity_summary",
                "scope": "all_queries",
                "metric_or_item": key,
                "relation_count": "",
                "query_count": fidelity_rows[0]["query_count"],
                "fidelity_value": value,
                "top100_value": "",
                "exhaustive_value": "",
                "exhaustive_minus_top100": "",
                "status": "post_hoc_sensitivity_only",
            }
        )
    for row in difference_rows:
        s2_rows.append(
            {
                "record_type": "metric_difference",
                "scope": row["scope"],
                "metric_or_item": row["metric"],
                "relation_count": row["relation_count"],
                "query_count": row["query_count"],
                "fidelity_value": "",
                "top100_value": row["top100_value"],
                "exhaustive_value": row["exhaustive_value"],
                "exhaustive_minus_top100": row["exhaustive_minus_top100"],
                "status": "post_hoc_sensitivity_only",
            }
        )
    s2_path = output_dir / "Table_S2_top100_exhaustive_fidelity.tsv"
    atomic_write_tsv(
        s2_path,
        ["record_type", "scope", "metric_or_item", "relation_count", "query_count", "fidelity_value", "top100_value", "exhaustive_value", "exhaustive_minus_top100", "status"],
        s2_rows,
    )

    # S3: complete zero/failure accounting derived from the unfiltered T3 rows.
    attribution_rows = read_tsv(paths["endpoint_attribution"])
    attribution_columns = {
        "scope", "baseline", "endpoint_relation_count", "endpoint_query_count", "endpoint_target_count",
        "endpoint_zero_score_relation_count", "endpoint_positive_score_relation_count",
        "endpoint_exact_tied_relation_count", "endpoint_zero_score_tied_relation_count",
        "endpoint_positive_score_tied_relation_count", "endpoint_positive_unique_score_relation_count",
        "endpoint_rank_le_50_relation_count", "endpoint_rank_gt_50_relation_count",
        "scope_all_zero_vector_query_count", "mrr_first_hit_zero_score_tied_query_count",
        "mrr_first_hit_positive_score_tied_query_count", "mrr_first_hit_positive_unique_score_query_count",
        "mrr_first_hit_other_query_count", "status",
    }
    require_columns(attribution_rows, attribution_columns, "S3 endpoint-score attribution")
    require_unique_keys(attribution_rows, ("scope", "baseline"), set(itertools.product(expected_scopes, expected_baselines)), "S3 endpoint-score attribution")
    attribution_lookup = {(row["scope"], row["baseline"]): row for row in attribution_rows}
    s3_rows = []
    for row in metric_rows:
        attribution = attribution_lookup[(row["scope"], row["baseline"])]
        if attribution["endpoint_relation_count"] != row["candidate_relation_count"] or attribution["endpoint_query_count"] != row["evaluable_query_count"]:
            raise ValueError("S3 metric and endpoint-score attribution denominators disagree")
        relation_count = int(attribution["endpoint_relation_count"])
        query_count = int(attribution["endpoint_query_count"])
        relation_attribution_total = sum(
            int(attribution[field])
            for field in (
                "endpoint_zero_score_tied_relation_count",
                "endpoint_positive_score_tied_relation_count",
                "endpoint_positive_unique_score_relation_count",
            )
        )
        if relation_attribution_total != relation_count:
            raise ValueError("S3 endpoint-score three-way attribution does not equal the relation denominator")
        mrr_attribution_total = sum(
            int(attribution[field])
            for field in (
                "mrr_first_hit_zero_score_tied_query_count",
                "mrr_first_hit_positive_score_tied_query_count",
                "mrr_first_hit_positive_unique_score_query_count",
                "mrr_first_hit_other_query_count",
            )
        )
        if mrr_attribution_total != query_count:
            raise ValueError("S3 MRR first-hit four-way attribution does not equal the query denominator")
        s3_rows.append(
            {
                "scope": row["scope"],
                "baseline": row["baseline"],
                "evaluable_query_count": row["evaluable_query_count"],
                "candidate_relation_count": row["candidate_relation_count"],
                "endpoint_target_count": attribution["endpoint_target_count"],
                "zero_recall_at_10_queries": row["zero_recall_at_10_queries"],
                "zero_recall_at_50_queries": row["zero_recall_at_50_queries"],
                "zero_mrr_queries": row["zero_mrr_queries"],
                "endpoint_zero_score_relation_count": attribution["endpoint_zero_score_relation_count"],
                "endpoint_positive_score_relation_count": attribution["endpoint_positive_score_relation_count"],
                "endpoint_exact_tied_relation_count": attribution["endpoint_exact_tied_relation_count"],
                "endpoint_zero_score_tied_relation_count": attribution["endpoint_zero_score_tied_relation_count"],
                "endpoint_positive_score_tied_relation_count": attribution["endpoint_positive_score_tied_relation_count"],
                "endpoint_positive_unique_score_relation_count": attribution["endpoint_positive_unique_score_relation_count"],
                "endpoint_rank_le_50_relation_count": attribution["endpoint_rank_le_50_relation_count"],
                "endpoint_rank_gt_50_relation_count": attribution["endpoint_rank_gt_50_relation_count"],
                "scope_all_zero_vector_query_count": attribution["scope_all_zero_vector_query_count"],
                "mrr_first_hit_zero_score_tied_query_count": attribution["mrr_first_hit_zero_score_tied_query_count"],
                "mrr_first_hit_positive_score_tied_query_count": attribution["mrr_first_hit_positive_score_tied_query_count"],
                "mrr_first_hit_positive_unique_score_query_count": attribution["mrr_first_hit_positive_unique_score_query_count"],
                "mrr_first_hit_other_query_count": attribution["mrr_first_hit_other_query_count"],
                "metric_status": row["status"],
                "attribution_status": attribution["status"],
            }
        )
    s3_path = output_dir / "Table_S3_zero_and_failure_accounting.tsv"
    atomic_write_tsv(s3_path, list(s3_rows[0].keys()), s3_rows)

    # S4: aggregate PMID concentration and document-component sensitivity.
    overlap_rows = read_tsv(paths["source_overlap"])
    component_scope_rows = read_tsv(paths["component_scopes"])
    component_metric_rows = read_tsv(paths["component_metrics"])
    component_focus_rows = read_tsv(paths["component_focus"])
    require_columns(overlap_rows, {"left_cohort", "right_cohort", "shared_source_document_count", "left_source_document_overlap_fraction", "right_source_document_overlap_fraction", "left_relation_count_with_cross_cohort_source", "left_relation_fraction_with_cross_cohort_source", "right_relation_count_with_cross_cohort_source", "right_relation_fraction_with_cross_cohort_source"}, "S4 cross-cohort overlap")
    if len(overlap_rows) != 1 or {overlap_rows[0]["left_cohort"], overlap_rows[0]["right_cohort"]} != {"endpoint", "historical"}:
        raise ValueError("S4 requires exactly one endpoint/historical overlap row")
    require_columns(component_scope_rows, {"scope", "relation_count", "query_count", "source_document_count", "query_source_document_edge_count", "component_count", "component_query_size_min", "component_query_size_q25", "component_query_size_median", "component_query_size_q75", "component_query_size_max", "component_query_size_mean", "component_source_document_size_min", "component_source_document_size_median", "component_source_document_size_max", "largest_component_query_fraction", "largest_component_source_document_fraction", "component_bootstrap_status"}, "S4 component scopes")
    require_unique_keys(component_scope_rows, ("scope",), {(scope,) for scope in expected_scopes}, "S4 component scopes")
    require_columns(component_metric_rows, {"scope", "baseline", "metric", "relation_count", "query_count", "source_document_count", "component_count", "point_estimate", "primary_query_bootstrap_ci95_low", "primary_query_bootstrap_ci95_high", "primary_query_bootstrap_status", "component_bootstrap_ci95_low", "component_bootstrap_ci95_high", "component_bootstrap_status", "component_bootstrap_replicates", "primary_query_interval_width", "component_interval_width", "component_to_primary_interval_width_ratio", "interval_width_ratio_status"}, "S4 component metrics")
    require_unique_keys(component_metric_rows, ("scope", "baseline", "metric"), set(itertools.product(expected_scopes, expected_baselines, expected_metrics)), "S4 component metrics")
    require_columns(component_focus_rows, {"scope", "left_baseline", "right_baseline", "metric", "relation_count", "query_count", "source_document_count", "component_count", "point_difference_left_minus_right", "primary_query_bootstrap_ci95_low", "primary_query_bootstrap_ci95_high", "primary_query_bootstrap_status", "component_bootstrap_ci95_low", "component_bootstrap_ci95_high", "component_bootstrap_status", "component_bootstrap_replicates", "primary_query_interval_width", "component_interval_width", "component_to_primary_interval_width_ratio", "interval_width_ratio_status"}, "S4 component focus")
    require_unique_keys(component_focus_rows, ("scope", "metric"), set(itertools.product(expected_scopes, expected_metrics)), "S4 component focus")
    if any(row["left_baseline"] != "structure_sequence_pair_neighbor" or row["right_baseline"] != "weighted_morgan_transfer" for row in component_focus_rows):
        raise ValueError("S4 component focus direction changed")
    s4_rows: list[dict[str, Any]] = []
    for row in source_rows:
        for key, value in row.items():
            if key == "cohort":
                continue
            s4_rows.append(
                {
                    "record_type": "source_concentration",
                    "scope_or_cohort": row["cohort"],
                    "baseline_or_left": "",
                    "right_baseline": "",
                    "metric_or_item": key,
                    "point_estimate": value,
                    "primary_ci95_low": "",
                    "primary_ci95_high": "",
                    "component_ci95_low": "",
                    "component_ci95_high": "",
                    "status": "descriptive_provenance_audit",
                }
            )
    source_summary = load_json(paths["source_summary"])
    source_summary_by_cohort = {row["cohort"]: row for row in source_summary["cohorts"]}
    require_exact_set(set(source_summary_by_cohort), {"historical", "endpoint"}, "S4 source-summary cohorts")
    for cohort, row in sorted(source_summary_by_cohort.items()):
        extra_values = {
            "source_document_effective_number": row["source_document_effective_number"],
            "top1pct_source_document_edge_share": row["top1pct_source_document_edge_share"],
            "top5_source_document_edge_share": row["top5_source_document_edge_share"],
            "source_documents_per_relation_median": row["source_documents_per_relation_distribution"]["median"],
            "source_documents_per_relation_max": row["source_documents_per_relation_distribution"]["max"],
            "source_documents_per_query_median": row["source_documents_per_query_distribution"]["median"],
            "source_documents_per_query_max": row["source_documents_per_query_distribution"]["max"],
            "relations_per_source_document_median": row["relations_per_source_document_distribution"]["median"],
            "relations_per_source_document_max": row["relations_per_source_document_distribution"]["max"],
            "queries_per_source_document_median": row["queries_per_source_document_distribution"]["median"],
            "queries_per_source_document_max": row["queries_per_source_document_distribution"]["max"],
        }
        for key, value in extra_values.items():
            s4_rows.append(
                {
                    "record_type": "source_concentration",
                    "scope_or_cohort": cohort,
                    "baseline_or_left": "",
                    "right_baseline": "",
                    "metric_or_item": key,
                    "point_estimate": value,
                    "primary_ci95_low": "",
                    "primary_ci95_high": "",
                    "component_ci95_low": "",
                    "component_ci95_high": "",
                    "status": "descriptive_provenance_audit",
                }
            )
    for row in overlap_rows:
        label = f"{row['left_cohort']}_vs_{row['right_cohort']}"
        for key, value in row.items():
            if key in {"left_cohort", "right_cohort"}:
                continue
            s4_rows.append(
                {
                    "record_type": "cross_cohort_pmid_overlap",
                    "scope_or_cohort": label,
                    "baseline_or_left": "",
                    "right_baseline": "",
                    "metric_or_item": key,
                    "point_estimate": value,
                    "primary_ci95_low": "",
                    "primary_ci95_high": "",
                    "component_ci95_low": "",
                    "component_ci95_high": "",
                    "status": "provenance_only_not_external_validation",
                }
            )
    for row in component_scope_rows:
        for key, value in row.items():
            if key in {"scope", "component_bootstrap_status"}:
                continue
            s4_rows.append(
                {
                    "record_type": "document_component_scope",
                    "scope_or_cohort": row["scope"],
                    "baseline_or_left": "",
                    "right_baseline": "",
                    "metric_or_item": key,
                    "point_estimate": value,
                    "primary_ci95_low": "",
                    "primary_ci95_high": "",
                    "component_ci95_low": "",
                    "component_ci95_high": "",
                    "status": row["component_bootstrap_status"],
                }
            )
    for row in component_metric_rows:
        s4_rows.append(
            {
                "record_type": "document_component_baseline_metric",
                "scope_or_cohort": row["scope"],
                "baseline_or_left": row["baseline"],
                "right_baseline": "",
                "metric_or_item": row["metric"],
                "point_estimate": row["point_estimate"],
                "primary_ci95_low": row["primary_query_bootstrap_ci95_low"],
                "primary_ci95_high": row["primary_query_bootstrap_ci95_high"],
                "component_ci95_low": row["component_bootstrap_ci95_low"],
                "component_ci95_high": row["component_bootstrap_ci95_high"],
                "status": row["component_bootstrap_status"],
            }
        )
    for row in component_focus_rows:
        s4_rows.append(
            {
                "record_type": "document_component_focus_contrast",
                "scope_or_cohort": row["scope"],
                "baseline_or_left": row["left_baseline"],
                "right_baseline": row["right_baseline"],
                "metric_or_item": row["metric"],
                "point_estimate": row["point_difference_left_minus_right"],
                "primary_ci95_low": row["primary_query_bootstrap_ci95_low"],
                "primary_ci95_high": row["primary_query_bootstrap_ci95_high"],
                "component_ci95_low": row["component_bootstrap_ci95_low"],
                "component_ci95_high": row["component_bootstrap_ci95_high"],
                "status": row["component_bootstrap_status"],
            }
        )
    s4_path = output_dir / "Table_S4_pmid_document_dependence.tsv"
    atomic_write_tsv(
        s4_path,
        ["record_type", "scope_or_cohort", "baseline_or_left", "right_baseline", "metric_or_item", "point_estimate", "primary_ci95_low", "primary_ci95_high", "component_ci95_low", "component_ci95_high", "status"],
        s4_rows,
    )

    # S5: frozen unresolved reasons, without relation identifiers or readmission.
    unresolved_rows = read_tsv(paths["unresolved"])
    require_columns(unresolved_rows, {"reason_category", "relation_count", "relation_fraction"}, "S5 unresolved")
    require_exact_set({row["reason_category"] for row in unresolved_rows}, {"preliminary_compound_unmatched", "preliminary_target_unmatched"}, "S5 unresolved reasons")
    if len(unresolved_rows) != 2:
        raise ValueError("S5 must contain exactly two frozen unresolved reason rows")
    s5_rows = [
        {
            **row,
            "frozen_unresolved_relation_count": unresolved_summary["frozen_unresolved_relation_count"],
            "distinct_compound_count": unresolved_summary["distinct_compound_count"],
            "distinct_target_count": unresolved_summary["distinct_target_count"],
            "negative_labels_emitted": unresolved_summary["negative_labels_emitted"],
            "readmission_performed": unresolved_summary["readmission_performed"],
            "status": "frozen_exclusion_no_negative_label",
        }
        for row in unresolved_rows
    ]
    s5_path = output_dir / "Table_S5_frozen_unresolved_exclusions.tsv"
    atomic_write_tsv(s5_path, list(s5_rows[0].keys()), s5_rows)

    # S6: reproducibility ledger and release state.  This table records only
    # aggregate/manuscript artifacts and immutable locks; package-specific
    # scanner reports are appended by the package builders after local staging.
    s6_rows: list[dict[str, Any]] = []
    shareability_by_name = {
        "protocol_lock": "non_identifying_local_review_candidate",
        "implementation_base_lock": "non_identifying_local_review_candidate",
        "implementation_lock": "non_identifying_local_review_candidate",
        "audit_implementation_lock": "non_identifying_local_review_candidate",
        "document_component_base_lock": "non_identifying_local_review_candidate",
        "document_component_lock": "non_identifying_local_review_candidate",
        "endpoint_attribution_implementation_lock": "non_identifying_local_review_candidate",
        "endpoint_attribution_governance": "non_identifying_local_review_candidate",
        "evaluation_manifest": "restricted_paths_present",
        "scoring_input_manifest": "non_identifying_local_review_candidate",
        "score_manifest": "restricted_paths_present",
        "complete_rank": "restricted_identifier",
        "tie_runtime": "non_identifying_local_review_candidate",
        "top100_failed_runtime": "non_identifying_local_review_candidate",
        "top100_retry_runtime": "non_identifying_local_review_candidate",
        "release_gate": "non_identifying_local_review_candidate",
        "component_manifest": "non_identifying_local_review_candidate",
        "endpoint_attribution_manifest": "non_identifying_local_review_candidate",
    }
    status_by_name = {
        "scoring_input_manifest": "verified_score_input_custody",
        "score_manifest": "verified_corrected_scoring_manifest",
        "complete_rank": "verified_against_score_manifest",
        "tie_runtime": "verified_success_receipt",
        "top100_failed_runtime": "verified_failed_attempt_receipt",
        "top100_retry_runtime": "verified_successful_retry_receipt",
        "release_gate": "verified_release_gate_source",
    }
    for name in (
        "protocol_lock",
        "implementation_base_lock",
        "implementation_lock",
        "audit_implementation_lock",
        "document_component_base_lock",
        "document_component_lock",
        "endpoint_attribution_implementation_lock",
        "endpoint_attribution_governance",
        "evaluation_manifest",
        "scoring_input_manifest",
        "score_manifest",
        "complete_rank",
        "tie_runtime",
        "top100_failed_runtime",
        "top100_retry_runtime",
        "release_gate",
        "component_manifest",
        "endpoint_attribution_manifest",
    ):
        path = paths[name]
        s6_rows.append(
            {
                "artifact_or_check": name,
                "relative_path": rel(path, project_root),
                "sha256": sha256_file(path),
                "value": path.stat().st_size,
                "status": status_by_name.get(name, "verified"),
                "shareability_tier": shareability_by_name[name],
            }
        )
    score_runtime = score_manifest["runtime"]
    score_peak_rss_bytes = score_runtime["process_peak_rss_bytes"]
    evidence_rows = [
        ("corrected_scoring_row_count", score_manifest["row_count"], "verified_from_score_manifest"),
        ("corrected_scoring_query_count", score_manifest["query_count"], "verified_from_score_manifest"),
        ("corrected_scoring_target_count", score_manifest["target_count"], "verified_from_score_manifest"),
        ("scoring_wall_seconds", score_runtime["wall_seconds_by_phase"]["total_before_manifest_write"], "measured"),
        ("scoring_peak_rss_bytes", score_peak_rss_bytes, "measured"),
        ("scoring_peak_rss_megabytes_decimal", score_peak_rss_bytes / 1_000_000, "derived_from_measured_bytes"),
        ("scoring_peak_rss_method", score_runtime["process_peak_rss_method"], "verified_from_score_manifest"),
        ("score_endpoint_file_supplied_to_command", score_manifest["endpoint_file_supplied_to_score_command"], "verified_from_score_manifest"),
        ("score_endpoint_read_by_engine", score_manifest["endpoint_read_by_score_engine"], "verified_from_score_manifest"),
        ("evaluation_wall_seconds", evaluation_manifest["runtime"]["wall_seconds_by_phase"]["total_before_manifest_write"], "measured"),
        ("evaluation_peak_rss_bytes", evaluation_manifest["runtime"]["process_peak_rss_bytes"], "measured"),
        ("tie_audit_wall_seconds", tie_runtime["wall_seconds"], "measured"),
        ("tie_audit_exit_code", tie_runtime["exit_code"], "successful_execution"),
        ("top100_first_attempt_wall_seconds", top100_failed_runtime["wall_seconds"], "measured_failed_attempt"),
        ("top100_first_attempt_exit_code", top100_failed_runtime["exit_code"], "failed_attempt"),
        ("top100_retry_wall_seconds", top100_retry_runtime["wall_seconds"], "measured_successful_retry"),
        ("top100_retry_exit_code", top100_retry_runtime["exit_code"], "successful_retry"),
        (
            "top100_failed_retry_same_command_sha256",
            top100_failed_runtime["command"]["argv_sha256"],
            "verified_same_command_fingerprint",
        ),
        ("document_component_wall_seconds", component_runtime["wall_seconds"], "measured"),
        ("document_component_peak_rss_bytes_internal", load_json(paths["component_summary"])["runtime"]["process_peak_rss_bytes"], "measured"),
        ("endpoint_score_attribution_wall_seconds", endpoint_attribution_runtime["wall_seconds"], "measured"),
    ]
    evidence_rows.extend(
        (field, release_gate[field], "gate_true" if release_gate[field] else "gate_false")
        for field in release_fields
    )
    for item, value, status in evidence_rows:
        s6_rows.append(
            {
                "artifact_or_check": item,
                "relative_path": "",
                "sha256": "",
                "value": value,
                "status": status,
                "shareability_tier": "non_identifying_local_review_candidate",
            }
        )
    s6_path = output_dir / "Table_S6_reproducibility_and_release.tsv"
    atomic_write_tsv(s6_path, ["artifact_or_check", "relative_path", "sha256", "value", "status", "shareability_tier"], s6_rows)

    table_paths = [t1_path, t2_path, t3_path, t4_path, t5_path, s1_path, s2_path, s3_path, s4_path, s5_path, s6_path]
    scan_aggregate_outputs(table_paths)
    for path in table_paths:
        outputs[rel(path, project_root)] = {
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
            "row_count_excluding_header": sum(1 for _ in path.open("r", encoding="utf-8")) - 1,
            "identifier_bearing": False,
        }

    manifest_path = output_dir / "revision_table_export_manifest_v4.json"
    manifest = {
        "schema_version": "1.0",
        "protocol_id": rebuild["protocol_id"],
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "claim_boundary": "Aggregate-only corrective tables; no external-validation or biological-discovery claim.",
        "identifier_bearing_outputs": False,
        "exporter": {
            "relative_path": rel(Path(__file__), project_root),
            "sha256": sha256_file(Path(__file__)),
        },
        "custody_checks": dict(sorted(custody["checks"].items())),
        "inputs": dict(sorted(inputs.items())),
        "outputs": dict(sorted(outputs.items())),
        "table_contract": {
            "all_five_scopes_retained": True,
            "all_four_baselines_retained": True,
            "all_five_metrics_retained": True,
            "zero_and_not_estimable_cells_not_filtered": True,
            "top100_remains_primary": True,
            "document_component_bootstrap_is_post_hoc_sensitivity": True,
            "table_content_is_deterministic_for_fixed_inputs": True,
            "manifest_timestamp_is_execution_metadata": True,
        },
    }
    atomic_write_json(manifest_path, manifest)
    return {"manifest": manifest_path, "outputs": outputs}


def parse_args() -> argparse.Namespace:
    default_root = Path(__file__).resolve().parents[2]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=default_root)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = build_tables(args.project_root.resolve())
    print(f"Wrote {len(result['outputs'])} aggregate-only tables")
    print(f"Manifest: {result['manifest']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
