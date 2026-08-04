"""Fail-closed v1.4 entry point for the locked document-component audit."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

import audit_document_component_bootstrap_v1 as base


V1_3_IMPLEMENTATION_LOCK_SHA256 = (
    "fed2a1cbe7840411fc383e9fc6e132ed2dd20caea4f5bbad8666e3fd5044af2b"
)
PROTOCOL_LOCK_SHA256 = (
    "96befee13ae1d41ad433c8697fac92ccd30fb25e24c3cf1279c6b4b7e040abd9"
)
SOURCE_AUTHORIZATION_SHA256 = (
    "7bbb511d034a96edeaaadc9e86d28176b17007f9fe7aece6bc0573037b2da42f"
)
AUDIT_SPEC_LOCK_SHA256 = (
    "c36dc3a6efbe304b7fbfef1a59d4acddc751e1281351b84e3b1341f4978b2d22"
)
FROZEN_EXECUTION_INPUT_SHA256 = {
    "corrective_score_manifest": "64df609ee4e83c6dd7efcb2285c38a0baec2240c99252d9849e25f1577379731",
    "corrective_evaluation_manifest": "938d243eb8064e949e750d246855e27c88c73e2b20c91551760084994afa5148",
    "primary_query_bootstrap_baselines": "07e1adc7c14b073de4a06e576cb7435da1d280488547411b9ab07ba5197c4c8f",
    "primary_query_bootstrap_focus_contrast": "151377f24ed98a82830cd3ce7f17617d4b7f42e863204153317958b6eed28590",
}
V1_4_IMPLEMENTATION_LOCK_NAME = "document_component_audit_implementation_code_lock_v1_4.json"
SOURCE_AUTHORIZATION_NAME = "document_component_source_evidence_authorization_v1.json"
AUDIT_SPEC_LOCK_NAME = "document_component_audit_spec_lock_v1.json"
SCORE_MANIFEST_NAME = "corrective_score_manifest.json"
EVALUATION_MANIFEST_NAME = "corrective_evaluation_manifest.json"
V1_3_IMPLEMENTATION_LOCK_RELATIVE = (
    "manifests/document_component_audit_implementation_code_lock_v1_3.json"
)
V1_4_AMENDMENT_ID = "document_cluster_audit_fail_closed_amendment_20260728_07"
V1_4_ADDED_FILES = frozenset(
    {
        "governance/document_cluster_audit_fail_closed_amendment_20260728_07.md",
        "audit_suite_v1_20260728/scripts/audit_document_component_bootstrap_v1_4.py",
        "audit_suite_v1_20260728/tests/test_document_component_bootstrap_v1_4.py",
        "manifests/document_component_source_evidence_authorization_v1.json",
        "manifests/document_component_audit_spec_lock_v1.json",
        "scripts/lock_document_component_audit_implementation_v1_4.py",
    }
)
PRIMARY_STATUS_VALUES = frozenset(
    {"not_estimable", "n=1_descriptive_only", "estimable_descriptive"}
)
SUCCESSOR_ROOT = Path(__file__).resolve().parents[2]


def parser() -> argparse.ArgumentParser:
    result = base.parser()
    result.description = __doc__
    result.add_argument("--source-evidence-manifest", required=True, type=Path)
    result.add_argument("--audit-spec-lock", required=True, type=Path)
    return result


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def verify_locked_manifest_file(path: Path, expected_name: str, expected_hash: str) -> dict[str, Any]:
    resolved = path.resolve()
    require(resolved.name == expected_name, f"Unexpected locked manifest name for {expected_name}")
    require(resolved.is_file(), f"Locked manifest is absent: {expected_name}")
    require(base.sha256(resolved) == expected_hash, f"Locked manifest hash changed: {expected_name}")
    return json.loads(resolved.read_text(encoding="utf-8"))


def verify_frozen_input_file(
    path: Path, expected_hash: str, label: str, expected_name: str | None = None
) -> Path:
    resolved = path.resolve()
    require(resolved.is_file(), f"Frozen {label} is absent")
    if expected_name is not None:
        require(resolved.name == expected_name, f"Frozen {label} basename changed")
    require(base.sha256(resolved) == expected_hash, f"Frozen {label} hash changed")
    return resolved


def verify_posthoc_prelock_disclosure(payload: dict[str, Any]) -> None:
    require(
        payload.get("source_evidence_rows_had_been_processed_before_lock") is True,
        "v1.4 lock understates prior source-evidence processing",
    )
    require(
        payload.get("corrected_rank_rows_had_been_processed_before_lock") is True,
        "v1.4 lock understates prior rank processing",
    )
    require(
        payload.get("primary_metric_rows_had_been_inspected_before_lock") is True,
        "v1.4 lock understates prior primary-result inspection",
    )
    require(
        payload.get("used_to_select_component_algorithm") is False,
        "Prior results were used to select the component algorithm",
    )
    require(
        payload.get("source_concentration_result_read_before_lock") is False,
        "Source-concentration results were read before the v1.4 lock",
    )
    require(
        payload.get("real_document_component_audit_executed_before_lock") is False,
        "A real component audit was executed before the v1.4 lock",
    )
    require(
        payload.get("document_component_implementation_testing_before_lock")
        == "synthetic_only",
        "Document-component implementation testing was not synthetic-only",
    )


def verify_source_authorization_payload(
    payload: dict[str, Any], source_evidence: Path
) -> dict[str, Any]:
    require(payload.get("schema_version") == "1.0", "Source authorization schema changed")
    require(payload.get("protocol_id") == base.PROTOCOL_ID, "Source authorization protocol mismatch")
    require(
        payload.get("authorization_scope") == "local_read_only_document_dependence_audit",
        "Source authorization scope changed",
    )
    require(payload.get("endpoint_pair_filter_required") is True, "Endpoint pair filter is not required")
    require(payload.get("numeric_pmid_rows_only") is True, "Numeric-PMID-only rule is not authorized")
    require(payload.get("external_redistribution_authorized") is False, "External redistribution was asserted")
    require(payload.get("identifier_release_authorized") is False, "Identifier release was asserted")
    require(
        payload.get("output_boundary")
        == "aggregate_only_no_query_pair_target_or_pmid_identifiers",
        "Source authorization output boundary changed",
    )
    source = payload.get("source_evidence")
    require(isinstance(source, dict), "Source authorization lacks a source-evidence receipt")
    resolved_source = source_evidence.resolve()
    recorded_path = Path(str(source.get("path", ""))).resolve()
    require(recorded_path == resolved_source, "Provided source evidence differs from authorized local path")
    require(resolved_source.is_file(), "Authorized source evidence is absent")
    require(source.get("basename") == resolved_source.name, "Authorized source basename changed")
    require(source.get("bytes") == resolved_source.stat().st_size, "Authorized source byte count changed")
    require(source.get("sha256") == base.sha256(resolved_source), "Authorized source SHA-256 changed")
    return source


def verify_audit_spec_payload(
    payload: dict[str, Any], source_authorization_path: Path, protocol_lock_path: Path
) -> None:
    require(payload.get("schema_version") == "1.0", "Audit-spec lock schema changed")
    require(payload.get("protocol_id") == base.PROTOCOL_ID, "Audit-spec lock protocol mismatch")
    require(
        payload.get("lock_state")
        == "LOCKED_BEFORE_SOURCE_CONCENTRATION_RESULT_INSPECTION_AND_BEFORE_COMPONENT_BOOTSTRAP_EXECUTION",
        "Audit spec was not locked at the required time",
    )
    require(payload.get("scientific_protocol_changed") is False, "Audit spec changes the scientific protocol")
    require(payload.get("identifier_release_authorized") is False, "Audit spec authorizes identifiers")
    require(payload.get("external_redistribution_authorized") is False, "Audit spec authorizes redistribution")
    protocol_receipt = payload.get("base_protocol_lock")
    require(isinstance(protocol_receipt, dict), "Audit spec lacks protocol-lock receipt")
    require(
        protocol_receipt.get("path") == "manifests/protocol_lock_manifest_v1.json",
        "Audit spec protocol-lock path changed",
    )
    require(
        protocol_receipt.get("sha256") == PROTOCOL_LOCK_SHA256
        and base.sha256(protocol_lock_path) == PROTOCOL_LOCK_SHA256,
        "Audit spec protocol-lock hash mismatch",
    )
    require(
        protocol_lock_path.resolve()
        == (SUCCESSOR_ROOT / "manifests/protocol_lock_manifest_v1.json").resolve(),
        "Provided protocol lock differs from the audit-spec path",
    )
    bootstrap = payload.get("bootstrap")
    require(isinstance(bootstrap, dict), "Audit spec lacks bootstrap rules")
    require(bootstrap.get("replicates") == base.BOOTSTRAP_REPLICATES, "Audit replicate count changed")
    require(bootstrap.get("prng") == "PCG64", "Audit PRNG changed")
    require(bootstrap.get("seed") == base.BOOTSTRAP_SEED, "Audit seed changed")
    require(
        bootstrap.get("resampling_unit")
        == "scope_specific_query_pmid_bipartite_connected_component",
        "Audit resampling unit changed",
    )
    locked_files = payload.get("locked_files")
    require(isinstance(locked_files, dict) and locked_files, "Audit spec has no locked-file inventory")
    expected_source_relative = "manifests/document_component_source_evidence_authorization_v1.json"
    require(
        locked_files.get(expected_source_relative) == SOURCE_AUTHORIZATION_SHA256,
        "Audit spec does not pin the source authorization",
    )
    require(
        source_authorization_path.resolve()
        == (SUCCESSOR_ROOT / expected_source_relative).resolve(),
        "Source authorization path differs from audit spec",
    )
    for relative_name, expected_hash in sorted(locked_files.items()):
        candidate = (SUCCESSOR_ROOT / relative_name).resolve()
        try:
            candidate.relative_to(SUCCESSOR_ROOT.resolve())
        except ValueError as error:
            raise ValueError("Audit-spec locked path escapes successor root") from error
        require(candidate.is_file(), f"Audit-spec locked file is absent: {relative_name}")
        require(base.sha256(candidate) == expected_hash, f"Audit-spec locked file drifted: {relative_name}")


def verify_v1_4_implementation_lock(path: Path) -> dict[str, Any]:
    resolved = path.resolve()
    require(resolved.name == V1_4_IMPLEMENTATION_LOCK_NAME, "Unexpected v1.4 implementation lock name")
    require(resolved.is_file(), "v1.4 implementation lock is absent")
    payload = json.loads(resolved.read_text(encoding="utf-8"))
    require(payload.get("schema_version") == "1.4", "Implementation lock schema is not 1.4")
    require(payload.get("protocol_id") == base.PROTOCOL_ID, "Implementation lock protocol mismatch")
    require(
        payload.get("lock_state") == "LOCKED_BEFORE_DOCUMENT_COMPONENT_AUDIT_V1_4_EXECUTION",
        "v1.4 implementation was not locked before execution",
    )
    base_lock = payload.get("base_lock")
    require(isinstance(base_lock, dict), "v1.4 implementation lock lacks its base")
    require(
        base_lock.get("path") == V1_3_IMPLEMENTATION_LOCK_RELATIVE
        and base_lock.get("sha256") == V1_3_IMPLEMENTATION_LOCK_SHA256,
        "v1.3 custody anchor differs from the hard-coded hash",
    )
    v1_3_path = SUCCESSOR_ROOT / V1_3_IMPLEMENTATION_LOCK_RELATIVE
    require(
        v1_3_path.is_file() and base.sha256(v1_3_path) == V1_3_IMPLEMENTATION_LOCK_SHA256,
        "v1.3 implementation lock drifted",
    )
    require(payload.get("amendment_id") == V1_4_AMENDMENT_ID, "v1.4 amendment ID changed")
    require(
        payload.get("frozen_execution_input_sha256") == FROZEN_EXECUTION_INPUT_SHA256,
        "v1.4 execution-input custody anchors changed",
    )
    added = payload.get("added_file_sha256")
    require(
        isinstance(added, dict) and set(added) == V1_4_ADDED_FILES,
        "v1.4 added-file inventory changed",
    )
    for relative_name, expected_hash in sorted(added.items()):
        candidate = (SUCCESSOR_ROOT / relative_name).resolve()
        try:
            candidate.relative_to(SUCCESSOR_ROOT.resolve())
        except ValueError as error:
            raise ValueError("v1.4 locked path escapes successor root") from error
        require(candidate.is_file(), f"v1.4 locked file is absent: {relative_name}")
        require(base.sha256(candidate) == expected_hash, f"v1.4 locked file drifted: {relative_name}")
    script_relative = "audit_suite_v1_20260728/scripts/audit_document_component_bootstrap_v1_4.py"
    require(added.get(script_relative) == base.sha256(Path(__file__)), "Executing v1.4 script is not locked")
    require(payload.get("scientific_protocol_changed") is False, "v1.4 lock changes the protocol")
    require(
        payload.get("endpoint_score_rank_or_primary_metric_changed") is False,
        "v1.4 lock changes a frozen scientific input",
    )
    verify_posthoc_prelock_disclosure(payload)
    require(
        payload.get("external_redistribution_authorized") is False,
        "v1.4 lock authorizes external redistribution",
    )
    require(payload.get("identifier_release_authorized") is False, "v1.4 lock authorizes identifiers")
    return payload


def load_prediction_ranks_strict(
    path: Path,
    needed_targets: dict[str, set[str]],
    endpoint_query_compounds: dict[str, str],
    expected_row_count: int = base.EXPECTED_RANK_ROWS,
) -> dict[str, dict[str, dict[str, int]]]:
    selected: dict[str, dict[str, dict[str, int]]] = {
        baseline: defaultdict(dict) for baseline in base.BASELINES
    }
    seen_groups: set[tuple[str, str]] = set()
    query_compounds: dict[str, str] = {}
    current_key: tuple[str, str] | None = None
    current_candidate_count = -1
    current_rank_seen: np.ndarray | None = None
    current_targets: set[str] = set()
    current_row_count = 0
    total_rows = 0
    candidate_count_by_query: dict[str, int] = {}

    def finish_group() -> None:
        nonlocal current_key, current_candidate_count, current_rank_seen
        nonlocal current_targets, current_row_count
        if current_key is None:
            return
        require(current_key not in seen_groups, "A rank baseline/query block recurs noncontiguously")
        seen_groups.add(current_key)
        require(current_rank_seen is not None, "Rank occupancy vector is absent")
        require(current_row_count == current_candidate_count, "Rank block row count differs from candidate count")
        require(len(current_targets) == current_candidate_count, "Rank block contains duplicate target rows")
        require(bool(np.all(current_rank_seen[1:])), "Rank block is not an exact 1..N permutation")
        current_key = None
        current_candidate_count = -1
        current_rank_seen = None
        current_targets = set()
        current_row_count = 0

    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = set(reader.fieldnames or [])
        required = {
            "protocol_id",
            "baseline",
            "query_id",
            "query_compound_inchikey_full",
            "target_uniprot_accession",
            "rank",
            "score",
            "eligible_candidate_target_count",
        }
        require(required.issubset(fields), "Corrective rank file lacks required fields")
        for row in reader:
            total_rows += 1
            require(row["protocol_id"] == base.PROTOCOL_ID, "Rank protocol ID mismatch")
            baseline = row["baseline"]
            query = row["query_id"]
            require(baseline in base.BASELINES, "Rank file contains an unknown baseline")
            require(query in endpoint_query_compounds, "Rank file contains an unknown query")
            key = (baseline, query)
            candidate_count = int(row["eligible_candidate_target_count"])
            require(
                1 <= candidate_count <= base.EXPECTED_CANDIDATE_TARGETS,
                "Rank block candidate count lies outside the frozen universe",
            )
            require(
                query not in candidate_count_by_query
                or candidate_count_by_query[query] == candidate_count,
                "Candidate count changes across baselines for a query",
            )
            candidate_count_by_query[query] = candidate_count
            if current_key != key:
                finish_group()
                require(key not in seen_groups, "A rank baseline/query block recurs noncontiguously")
                current_key = key
                current_candidate_count = candidate_count
                current_rank_seen = np.zeros(candidate_count + 1, dtype=np.bool_)
            require(candidate_count == current_candidate_count, "Candidate count changes within a rank block")
            rank = int(row["rank"])
            require(1 <= rank <= candidate_count, "Rank lies outside its candidate range")
            require(current_rank_seen is not None and not bool(current_rank_seen[rank]), "Duplicate rank in a block")
            current_rank_seen[rank] = True
            target = row["target_uniprot_accession"]
            require(target not in current_targets, "Duplicate target in a rank block")
            current_targets.add(target)
            current_row_count += 1
            score = float(row["score"])
            require(math.isfinite(score), "Rank file contains a nonfinite score")
            compound = row["query_compound_inchikey_full"]
            require(
                query not in query_compounds or query_compounds[query] == compound,
                "Rank query maps to multiple compounds",
            )
            query_compounds[query] = compound
            if target in needed_targets.get(query, set()):
                require(target not in selected[baseline][query], "Endpoint target repeats within a rank block")
                selected[baseline][query][target] = rank
    finish_group()
    require(total_rows == expected_row_count, "Corrective rank row count differs from the locked contract")
    require(query_compounds == endpoint_query_compounds, "Rank query/compound map differs from endpoint")
    expected_groups = {
        (baseline, query)
        for baseline in base.BASELINES
        for query in endpoint_query_compounds
    }
    require(seen_groups == expected_groups, "Rank baseline/query blocks are incomplete")
    for baseline in base.BASELINES:
        for query, targets in needed_targets.items():
            require(set(selected[baseline][query]) == targets, "Endpoint target ranks are incomplete")
    return {baseline: dict(rows) for baseline, rows in selected.items()}


def validate_primary_statuses(rows: dict[Any, dict[str, str]], label: str) -> None:
    for row in rows.values():
        require(row.get("status") in PRIMARY_STATUS_VALUES, f"{label} contains an unknown status")


def base_argv(args: argparse.Namespace) -> list[str]:
    return [
        str(Path(__file__).resolve()),
        "--protocol-lock",
        str(args.protocol_lock),
        "--implementation-lock",
        str(args.implementation_lock),
        "--score-manifest",
        str(args.score_manifest),
        "--evaluation-manifest",
        str(args.evaluation_manifest),
        "--endpoint",
        str(args.endpoint),
        "--source-evidence",
        str(args.source_evidence),
        "--ranks",
        str(args.ranks),
        "--scaffold-audit",
        str(args.scaffold_audit),
        "--homology-0-30",
        str(args.homology_0_30),
        "--homology-0-50",
        str(args.homology_0_50),
        "--homology-0-70",
        str(args.homology_0_70),
        "--primary-baseline-bootstrap",
        str(args.primary_baseline_bootstrap),
        "--primary-focus-contrast",
        str(args.primary_focus_contrast),
        "--output-dir",
        str(args.output_dir),
    ]


def main() -> int:
    args = parser().parse_args()
    base.require_protocol_lock(args.protocol_lock)
    verify_v1_4_implementation_lock(args.implementation_lock)
    source_authorization = verify_locked_manifest_file(
        args.source_evidence_manifest, SOURCE_AUTHORIZATION_NAME, SOURCE_AUTHORIZATION_SHA256
    )
    audit_spec = verify_locked_manifest_file(
        args.audit_spec_lock, AUDIT_SPEC_LOCK_NAME, AUDIT_SPEC_LOCK_SHA256
    )
    verify_audit_spec_payload(
        audit_spec, args.source_evidence_manifest, args.protocol_lock
    )
    verify_frozen_input_file(
        args.score_manifest,
        FROZEN_EXECUTION_INPUT_SHA256["corrective_score_manifest"],
        "corrective score manifest",
        SCORE_MANIFEST_NAME,
    )
    verify_frozen_input_file(
        args.evaluation_manifest,
        FROZEN_EXECUTION_INPUT_SHA256["corrective_evaluation_manifest"],
        "corrective evaluation manifest",
        EVALUATION_MANIFEST_NAME,
    )
    verify_frozen_input_file(
        args.primary_baseline_bootstrap,
        FROZEN_EXECUTION_INPUT_SHA256["primary_query_bootstrap_baselines"],
        "primary query-bootstrap baseline table",
    )
    verify_frozen_input_file(
        args.primary_focus_contrast,
        FROZEN_EXECUTION_INPUT_SHA256["primary_query_bootstrap_focus_contrast"],
        "primary query-bootstrap focus table",
    )
    verify_source_authorization_payload(source_authorization, args.source_evidence)

    original_base_file = base.__file__
    original_verify_implementation_lock = base.verify_implementation_lock
    original_load_prediction_ranks = base.load_prediction_ranks
    original_finalize = base.finalize_manifest
    original_load_primary_baseline = base.load_primary_baseline
    original_load_primary_focus = base.load_primary_focus

    def load_primary_baseline_strict(path: Path) -> dict[tuple[str, str, str], dict[str, str]]:
        rows = original_load_primary_baseline(path)
        validate_primary_statuses(rows, "Primary baseline-bootstrap table")
        return rows

    def load_primary_focus_strict(path: Path) -> dict[tuple[str, str], dict[str, str]]:
        rows = original_load_primary_focus(path)
        validate_primary_statuses(rows, "Primary focus-contrast table")
        return rows

    def finalize_with_v1_4_receipts(*, inputs: Any, extra: Any = None, **kwargs: Any) -> dict[str, Any]:
        augmented_inputs = list(inputs) + [
            base.input_descriptor(
                "source_evidence_authorization_manifest", args.source_evidence_manifest
            ),
            base.input_descriptor("document_component_audit_spec_lock", args.audit_spec_lock),
        ]
        augmented_extra = dict(extra or {})
        augmented_extra.update(
            {
                "fail_closed_implementation_version": "1.4",
                "v1_3_custody_anchor_sha256": V1_3_IMPLEMENTATION_LOCK_SHA256,
                "source_evidence_authorization": {
                    "basename": args.source_evidence_manifest.name,
                    "sha256": SOURCE_AUTHORIZATION_SHA256,
                },
                "audit_spec_lock": {
                    "basename": args.audit_spec_lock.name,
                    "sha256": AUDIT_SPEC_LOCK_SHA256,
                },
                "rank_permutation_validation": "per_block_boolean_occupancy_and_target_uniqueness",
                "frozen_execution_input_sha256": dict(FROZEN_EXECUTION_INPUT_SHA256),
            }
        )
        return original_finalize(inputs=augmented_inputs, extra=augmented_extra, **kwargs)

    base.verify_implementation_lock = verify_v1_4_implementation_lock
    base.load_prediction_ranks = load_prediction_ranks_strict
    base.load_primary_baseline = load_primary_baseline_strict
    base.load_primary_focus = load_primary_focus_strict
    base.finalize_manifest = finalize_with_v1_4_receipts
    base.__file__ = str(Path(__file__).resolve())
    saved_argv = sys.argv
    try:
        sys.argv = base_argv(args)
        return base.main()
    finally:
        sys.argv = saved_argv
        base.__file__ = original_base_file
        base.verify_implementation_lock = original_verify_implementation_lock
        base.load_prediction_ranks = original_load_prediction_ranks
        base.load_primary_baseline = original_load_primary_baseline
        base.load_primary_focus = original_load_primary_focus
        base.finalize_manifest = original_finalize


if __name__ == "__main__":
    raise SystemExit(main())
