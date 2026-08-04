"""Aggregate-only query-PMID connected-component bootstrap sensitivity."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import sys
import time
import tracemalloc
from collections import defaultdict
from pathlib import Path
from typing import Any

import numpy as np

from audit_common import (
    PROTOCOL_ID,
    SUITE_ID,
    distribution,
    finalize_manifest,
    input_descriptor,
    open_dict_reader,
    parse_bool,
    require_new_output_dir,
    require_protocol_lock,
    sha256,
    utc_now,
    write_json_new,
    write_tsv_new,
)


AUDIT_ID = "source_document_connected_component_bootstrap_v1"
RUN_ID = "npass_strict_ab_asof_cutoff_author_run_v1_20260728"
RUN_MODE = "author_run_non_independent_corrective_successor"
FROZEN_ENDPOINT_SHA256 = "09296b066a23197a7c178f00514f2b3d9ed7e6f3c459ea92a55e01a6010d1132"
BASE_IMPLEMENTATION_LOCK_SHA256 = "5cca17a44abbec99e2f978c956a8a49eb4f104f9089d966adc9cf035e8921ad4"
IMPLEMENTATION_LOCK_NAME = "document_component_audit_implementation_code_lock_v1_3.json"

BASELINES = [
    "weighted_target_popularity",
    "sequence_3mer_transfer",
    "weighted_morgan_transfer",
    "structure_sequence_pair_neighbor",
]
SCOPES = [
    "temporal_strict_ab",
    "scaffold_cold_strict_ab",
    "double_cold_0_30",
    "double_cold_0_50",
    "double_cold_0_70",
]
METRICS = ["Recall@10", "Recall@50", "NDCG@10", "NDCG@50", "MRR"]
FOCUS_LEFT = "structure_sequence_pair_neighbor"
FOCUS_RIGHT = "weighted_morgan_transfer"
EXPECTED_ENDPOINT_RELATIONS = 358
EXPECTED_ENDPOINT_QUERIES = 222
EXPECTED_ENDPOINT_TARGETS = 156
EXPECTED_RANK_ROWS = 3_658_128
EXPECTED_CANDIDATE_TARGETS = 4_123
BOOTSTRAP_REPLICATES = 10_000
BOOTSTRAP_SEED = 20_260_719
BOOTSTRAP_BATCH_SIZE = 256

SUITE_ROOT = Path(__file__).resolve().parents[1]
SUCCESSOR_ROOT = SUITE_ROOT.parent


class UnionFind:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}
        self.size: dict[str, int] = {}

    def add(self, value: str) -> None:
        if value not in self.parent:
            self.parent[value] = value
            self.size[value] = 1

    def find(self, value: str) -> str:
        parent = self.parent[value]
        if parent != value:
            self.parent[value] = self.find(parent)
        return self.parent[value]

    def union(self, left: str, right: str) -> None:
        self.add(left)
        self.add(right)
        root_left = self.find(left)
        root_right = self.find(right)
        if root_left == root_right:
            return
        if self.size[root_left] < self.size[root_right]:
            root_left, root_right = root_right, root_left
        self.parent[root_right] = root_left
        self.size[root_left] += self.size[root_right]


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--protocol-lock", required=True, type=Path)
    result.add_argument("--implementation-lock", required=True, type=Path)
    result.add_argument("--score-manifest", required=True, type=Path)
    result.add_argument("--evaluation-manifest", required=True, type=Path)
    result.add_argument("--endpoint", required=True, type=Path)
    result.add_argument("--source-evidence", required=True, type=Path)
    result.add_argument("--ranks", required=True, type=Path)
    result.add_argument("--scaffold-audit", required=True, type=Path)
    result.add_argument("--homology-0-30", required=True, type=Path)
    result.add_argument("--homology-0-50", required=True, type=Path)
    result.add_argument("--homology-0-70", required=True, type=Path)
    result.add_argument("--primary-baseline-bootstrap", required=True, type=Path)
    result.add_argument("--primary-focus-contrast", required=True, type=Path)
    result.add_argument("--output-dir", required=True, type=Path)
    return result


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def verify_implementation_lock(path: Path) -> dict[str, Any]:
    path = path.resolve()
    require(path.name == IMPLEMENTATION_LOCK_NAME, "Unexpected document-component implementation lock name")
    require(path.is_file(), "Document-component implementation lock is absent")
    payload = json.loads(path.read_text(encoding="utf-8"))
    require(payload.get("schema_version") == "1.3", "Implementation lock schema is not 1.3")
    require(payload.get("protocol_id") == PROTOCOL_ID, "Implementation lock protocol mismatch")
    require(
        payload.get("lock_state") == "LOCKED_BEFORE_DOCUMENT_COMPONENT_AUDIT_EXECUTION",
        "Document-component code was not locked before execution",
    )
    base = payload.get("base_lock")
    require(isinstance(base, dict), "Implementation lock lacks its base lock")
    require(base.get("sha256") == BASE_IMPLEMENTATION_LOCK_SHA256, "Base implementation lock hash mismatch")
    base_path = SUCCESSOR_ROOT / str(base.get("path", ""))
    require(base_path.is_file() and sha256(base_path) == BASE_IMPLEMENTATION_LOCK_SHA256, "Base lock drifted")
    added = payload.get("added_file_sha256")
    require(isinstance(added, dict) and added, "Implementation lock has no added-file inventory")
    for relative_name, expected_hash in sorted(added.items()):
        candidate = (SUCCESSOR_ROOT / relative_name).resolve()
        try:
            candidate.relative_to(SUCCESSOR_ROOT.resolve())
        except ValueError as error:
            raise ValueError("Implementation lock path escapes successor root") from error
        require(candidate.is_file(), f"Locked document-component file is absent: {relative_name}")
        require(sha256(candidate) == expected_hash, f"Locked document-component file drifted: {relative_name}")
    script_relative = "audit_suite_v1_20260728/scripts/audit_document_component_bootstrap_v1.py"
    require(added.get(script_relative) == sha256(Path(__file__)), "Executing audit script is not the locked version")
    require(payload.get("identifier_release_authorized") is False, "Implementation lock cannot authorize identifiers")
    return payload


def verify_file_receipt(item: Any, label: str) -> Path:
    require(isinstance(item, dict), f"Manifest lacks {label} file receipt")
    path = Path(str(item.get("path", "")))
    require(path.is_file(), f"Recorded {label} file is absent")
    require(item.get("sha256") == sha256(path), f"Recorded {label} file hash changed")
    return path.resolve()


def verify_receipt_map(items: Any, label: str) -> dict[str, Path]:
    require(isinstance(items, dict) and items, f"Manifest lacks {label} receipts")
    return {name: verify_file_receipt(item, f"{label}.{name}") for name, item in items.items()}


def verify_corrected_manifests(
    *,
    score_manifest_path: Path,
    evaluation_manifest_path: Path,
    ranks: Path,
    endpoint: Path,
    scaffold: Path,
    homology: dict[str, Path],
    primary_baseline: Path,
    primary_focus: Path,
) -> tuple[dict[str, Any], dict[str, Any]]:
    score_manifest = json.loads(score_manifest_path.read_text(encoding="utf-8"))
    require(score_manifest.get("protocol_id") == PROTOCOL_ID, "Score manifest protocol mismatch")
    require(score_manifest.get("run_id") == RUN_ID, "Score manifest run mismatch")
    require(score_manifest.get("stage") == "corrective_score", "Score manifest stage mismatch")
    require(score_manifest.get("execution_mode") == RUN_MODE, "Score manifest mode mismatch")
    require(score_manifest.get("row_count") == EXPECTED_RANK_ROWS, "Score manifest rank-row count changed")
    require(score_manifest.get("query_count") == EXPECTED_ENDPOINT_QUERIES, "Score manifest query count changed")
    require(score_manifest.get("target_count") == EXPECTED_CANDIDATE_TARGETS, "Score target universe changed")
    require(score_manifest.get("baselines") == BASELINES, "Score baseline order or membership changed")
    require(score_manifest.get("endpoint_file_supplied_to_score_command") is False, "Score received endpoint")
    require(score_manifest.get("endpoint_read_by_score_engine") is False, "Score engine read endpoint")
    require(score_manifest.get("legacy_outer_or_result_read") is False, "Score read legacy results")
    recorded_rank = verify_file_receipt(score_manifest.get("rank_output"), "score rank output")
    require(recorded_rank == ranks.resolve(), "Provided rank path differs from score manifest")
    verify_file_receipt(score_manifest.get("receipt"), "score execution receipt")
    verify_receipt_map(score_manifest.get("inputs"), "score inputs")
    verify_receipt_map(score_manifest.get("code"), "score code")

    evaluation_manifest = json.loads(evaluation_manifest_path.read_text(encoding="utf-8"))
    require(evaluation_manifest.get("protocol_id") == PROTOCOL_ID, "Evaluation manifest protocol mismatch")
    require(evaluation_manifest.get("run_id") == RUN_ID, "Evaluation manifest run mismatch")
    require(evaluation_manifest.get("stage") == "corrective_evaluation", "Evaluation manifest stage mismatch")
    require(evaluation_manifest.get("execution_mode") == RUN_MODE, "Evaluation manifest mode mismatch")
    require(evaluation_manifest.get("all_scopes_reported") == SCOPES, "Evaluation scope contract changed")
    require(evaluation_manifest.get("all_baselines_reported") == BASELINES, "Evaluation baseline contract changed")
    require(evaluation_manifest.get("all_metrics_reported") == METRICS, "Evaluation metric contract changed")
    bootstrap = evaluation_manifest.get("bootstrap", {})
    require(bootstrap.get("unit") == "query_compound", "Primary bootstrap is not query-level")
    require(bootstrap.get("replicates") == BOOTSTRAP_REPLICATES, "Primary bootstrap replicate count changed")
    require(bootstrap.get("prng") == "PCG64", "Primary bootstrap PRNG changed")
    require(bootstrap.get("seed") == BOOTSTRAP_SEED, "Primary bootstrap seed changed")
    focus = evaluation_manifest.get("focus_contrast", {})
    require(focus.get("left_baseline") == FOCUS_LEFT, "Focus contrast left baseline changed")
    require(focus.get("right_baseline") == FOCUS_RIGHT, "Focus contrast right baseline changed")
    require(focus.get("row_count") == 25, "Primary focus contrast row count changed")
    recorded_score = verify_file_receipt(evaluation_manifest.get("score_manifest"), "evaluation score manifest")
    require(recorded_score == score_manifest_path.resolve(), "Evaluation records another score manifest")
    recorded_eval_rank = verify_file_receipt(evaluation_manifest.get("score_rank"), "evaluation score rank")
    require(recorded_eval_rank == ranks.resolve(), "Evaluation records another rank file")
    verify_file_receipt(evaluation_manifest.get("receipt"), "evaluation execution receipt")
    verify_file_receipt(evaluation_manifest.get("evaluation_input_manifest"), "evaluation input manifest")
    evaluation_inputs = verify_receipt_map(evaluation_manifest.get("evaluation_inputs"), "evaluation inputs")
    expected_inputs = {
        "endpoint": endpoint.resolve(),
        "scaffold_audit": scaffold.resolve(),
        "homology_0_30": homology["0_30"].resolve(),
        "homology_0_50": homology["0_50"].resolve(),
        "homology_0_70": homology["0_70"].resolve(),
    }
    require(evaluation_inputs == expected_inputs, "Provided evaluation input paths differ from manifest")
    outputs = verify_receipt_map(evaluation_manifest.get("outputs"), "evaluation outputs")
    require(outputs.get("baseline_bootstrap") == primary_baseline.resolve(), "Primary baseline table mismatch")
    require(
        outputs.get("focus_pair_neighbor_minus_morgan") == primary_focus.resolve(),
        "Primary focus table mismatch",
    )
    verify_receipt_map(evaluation_manifest.get("code"), "evaluation code")
    return score_manifest, evaluation_manifest


def load_endpoint(path: Path) -> list[dict[str, str]]:
    require(sha256(path) == FROZEN_ENDPOINT_SHA256, "Frozen endpoint byte hash changed")
    rows: list[dict[str, str]] = []
    with open_dict_reader(path) as reader:
        fields = set(reader.fieldnames or [])
        required = {
            "canonical_pair_key",
            "query_id",
            "inchikey_full",
            "uniprot_canonical_accession",
            "best_strict_evidence_tier",
            "decision",
            "c31_leakage_gate_status",
        }
        require(required.issubset(fields), "Frozen endpoint lacks required fields")
        rows.extend(reader)
    require(len(rows) == EXPECTED_ENDPOINT_RELATIONS, "Endpoint relation count is not 358")
    pair_keys = [row["canonical_pair_key"] for row in rows]
    relation_keys = [(row["inchikey_full"], row["uniprot_canonical_accession"]) for row in rows]
    require(len(set(pair_keys)) == len(pair_keys), "Endpoint relation keys are not unique")
    require(len(set(relation_keys)) == len(relation_keys), "Endpoint compound-target pairs are not unique")
    require(len({row["query_id"] for row in rows}) == EXPECTED_ENDPOINT_QUERIES, "Endpoint query count is not 222")
    require(
        len({row["uniprot_canonical_accession"] for row in rows}) == EXPECTED_ENDPOINT_TARGETS,
        "Endpoint target count is not 156",
    )
    query_compounds: dict[str, str] = {}
    for row in rows:
        require(row["best_strict_evidence_tier"] in {"A_affinity_candidate", "B_quantitative_functional_candidate"}, "Endpoint tier is not strict A/B")
        require(row["decision"] == "strict_post_cutoff_future_candidate", "Endpoint decision changed")
        require(row["c31_leakage_gate_status"] == "pass_no_historical_activity", "Endpoint leakage gate changed")
        query = row["query_id"]
        compound = row["inchikey_full"]
        require(query not in query_compounds or query_compounds[query] == compound, "Query maps to multiple compounds")
        query_compounds[query] = compound
    require(
        len(set(query_compounds.values())) == EXPECTED_ENDPOINT_QUERIES,
        "Endpoint query IDs are not one-to-one with query compounds",
    )
    return rows


def load_source_documents(
    path: Path, endpoint_pair_keys: set[str]
) -> tuple[dict[str, set[str]], dict[str, int]]:
    pair_documents: dict[str, set[str]] = defaultdict(set)
    total_rows = 0
    non_endpoint_rows = 0
    endpoint_ineligible_reference_rows = 0
    endpoint_eligible_reference_rows = 0
    with open_dict_reader(path) as reader:
        fields = set(reader.fieldnames or [])
        required = {"canonical_pair_key", "ref_id_type", "ref_id"}
        require(required.issubset(fields), "Source evidence lacks required fields")
        for row in reader:
            total_rows += 1
            pair = row["canonical_pair_key"].strip()
            if pair not in endpoint_pair_keys:
                non_endpoint_rows += 1
                continue
            reference_type = row["ref_id_type"].strip().upper()
            reference = row["ref_id"].strip()
            if reference_type != "PMID" or not reference.isdigit():
                endpoint_ineligible_reference_rows += 1
                continue
            pair_documents[pair].add(reference)
            endpoint_eligible_reference_rows += 1
    missing_pair_count = len(endpoint_pair_keys.difference(pair_documents))
    require(missing_pair_count == 0, "At least one endpoint relation lacks a numeric PMID")
    return dict(pair_documents), {
        "total_source_evidence_row_count": total_rows,
        "non_endpoint_source_evidence_row_count": non_endpoint_rows,
        "endpoint_ineligible_reference_row_count": endpoint_ineligible_reference_rows,
        "endpoint_eligible_reference_row_count": endpoint_eligible_reference_rows,
        "endpoint_relation_count_with_numeric_pmid": len(pair_documents),
    }


def load_bool_map(path: Path, key_field: str, flag_field: str, status_field: str) -> dict[str, bool]:
    output: dict[str, bool] = {}
    with open_dict_reader(path) as reader:
        fields = set(reader.fieldnames or [])
        require({key_field, flag_field, status_field}.issubset(fields), "Scope mask lacks required fields")
        for row in reader:
            key = row[key_field].strip()
            require(bool(key) and key not in output, "Scope mask has an empty or duplicate key")
            require(bool(row[status_field].strip()), "Scope mask has an empty status")
            output[key] = parse_bool(row[flag_field])
    require(bool(output), "Scope mask is empty")
    return output


def build_scope_rows(
    endpoint: list[dict[str, str]], scaffold: dict[str, bool], homology: dict[str, dict[str, bool]]
) -> dict[str, list[dict[str, str]]]:
    endpoint_keys = {row["canonical_pair_key"] for row in endpoint}
    endpoint_targets = {row["uniprot_canonical_accession"] for row in endpoint}
    require(set(scaffold) == endpoint_keys, "Scaffold mask keyset differs from endpoint")
    for threshold, flags in homology.items():
        require(set(flags) == endpoint_targets, f"Homology {threshold} keyset differs from endpoint")
    output = {scope: [] for scope in SCOPES}
    for row in endpoint:
        output["temporal_strict_ab"].append(row)
        if scaffold[row["canonical_pair_key"]]:
            output["scaffold_cold_strict_ab"].append(row)
            for threshold in ("0_30", "0_50", "0_70"):
                if homology[threshold][row["uniprot_canonical_accession"]]:
                    output[f"double_cold_{threshold}"].append(row)
    return output


def build_query_document_components(
    scope_rows: list[dict[str, str]], pair_documents: dict[str, set[str]]
) -> tuple[list[dict[str, list[str]]], int, int]:
    query_documents: dict[str, set[str]] = defaultdict(set)
    for row in scope_rows:
        documents = pair_documents[row["canonical_pair_key"]]
        require(bool(documents), "Scope relation lacks an eligible source document")
        query_documents[row["query_id"]].update(documents)
    if not query_documents:
        return [], 0, 0

    union = UnionFind()
    for query, documents in query_documents.items():
        query_node = f"query:{query}"
        union.add(query_node)
        for document in documents:
            document_node = f"document:{document}"
            union.union(query_node, document_node)
    component_queries: dict[str, list[str]] = defaultdict(list)
    component_documents: dict[str, list[str]] = defaultdict(list)
    for query in query_documents:
        component_queries[union.find(f"query:{query}")].append(query)
    all_documents = set().union(*query_documents.values())
    for document in all_documents:
        component_documents[union.find(f"document:{document}")].append(document)
    roots = set(component_queries) | set(component_documents)
    components = [
        {
            "queries": sorted(component_queries[root]),
            "documents": sorted(component_documents[root]),
        }
        for root in roots
    ]
    components.sort(key=lambda item: tuple(item["queries"]))
    require(sum(len(item["queries"]) for item in components) == len(query_documents), "Component query partition failed")
    edge_count = sum(len(documents) for documents in query_documents.values())
    return components, len(all_documents), edge_count


def load_prediction_ranks(
    path: Path, needed_targets: dict[str, set[str]], endpoint_query_compounds: dict[str, str]
) -> dict[str, dict[str, dict[str, int]]]:
    selected: dict[str, dict[str, dict[str, int]]] = {
        baseline: defaultdict(dict) for baseline in BASELINES
    }
    audit: dict[tuple[str, str], dict[str, int]] = defaultdict(
        lambda: {
            "count": 0,
            "rank_sum": 0,
            "rank_square_sum": 0,
            "rank_min": 2**31 - 1,
            "rank_max": 0,
            "candidate_count": -1,
        }
    )
    query_compounds: dict[str, str] = {}
    seen_rows = 0
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
            seen_rows += 1
            require(row["protocol_id"] == PROTOCOL_ID, "Rank protocol ID mismatch")
            baseline = row["baseline"]
            require(baseline in BASELINES, "Rank file contains an unknown baseline")
            query = row["query_id"]
            compound = row["query_compound_inchikey_full"]
            require(query not in query_compounds or query_compounds[query] == compound, "Rank query maps to multiple compounds")
            query_compounds[query] = compound
            rank = int(row["rank"])
            score = float(row["score"])
            candidate_count = int(row["eligible_candidate_target_count"])
            require(math.isfinite(score), "Rank file contains a nonfinite score")
            require(1 <= rank <= candidate_count, "Rank lies outside candidate range")
            key = (baseline, query)
            item = audit[key]
            require(item["candidate_count"] in {-1, candidate_count}, "Candidate count changes within a rank block")
            item["candidate_count"] = candidate_count
            item["count"] += 1
            item["rank_sum"] += rank
            item["rank_square_sum"] += rank * rank
            item["rank_min"] = min(item["rank_min"], rank)
            item["rank_max"] = max(item["rank_max"], rank)
            target = row["target_uniprot_accession"]
            if target in needed_targets.get(query, set()):
                require(target not in selected[baseline][query], "Endpoint target repeats within a rank block")
                selected[baseline][query][target] = rank
    require(seen_rows == EXPECTED_RANK_ROWS, "Corrective rank row count is not 3,658,128")
    require(query_compounds == endpoint_query_compounds, "Rank query/compound map differs from endpoint")
    expected_blocks = {(baseline, query) for baseline in BASELINES for query in endpoint_query_compounds}
    require(set(audit) == expected_blocks, "Rank baseline/query blocks are incomplete")
    for item in audit.values():
        count = item["count"]
        expected_sum = count * (count + 1) // 2
        expected_square_sum = count * (count + 1) * (2 * count + 1) // 6
        require(
            count == item["candidate_count"]
            and item["rank_min"] == 1
            and item["rank_max"] == count
            and item["rank_sum"] == expected_sum
            and item["rank_square_sum"] == expected_square_sum,
            "Corrective rank permutation check failed",
        )
    for baseline in BASELINES:
        for query, targets in needed_targets.items():
            require(set(selected[baseline][query]) == targets, "Endpoint target ranks are incomplete")
    return {baseline: dict(rows) for baseline, rows in selected.items()}


def query_metric_vector(positive_ranks: list[int]) -> np.ndarray:
    require(bool(positive_ranks) and all(rank >= 1 for rank in positive_ranks), "Invalid positive ranks")
    values: list[float] = []
    for k in (10, 50):
        values.append(sum(rank <= k for rank in positive_ranks) / len(positive_ranks))
    for k in (10, 50):
        dcg = sum(1.0 / math.log2(rank + 1) for rank in positive_ranks if rank <= k)
        ideal = sum(1.0 / math.log2(index + 2) for index in range(min(k, len(positive_ranks))))
        values.append(dcg / ideal if ideal else 0.0)
    values.append(1.0 / min(positive_ranks))
    return np.asarray(values, dtype=np.float64)


def build_scope_arrays(
    scope_rows: list[dict[str, str]], selected: dict[str, dict[str, dict[str, int]]]
) -> tuple[list[str], dict[str, np.ndarray]]:
    positives: dict[str, list[str]] = defaultdict(list)
    for row in scope_rows:
        positives[row["query_id"]].append(row["uniprot_canonical_accession"])
    queries = sorted(positives)
    arrays: dict[str, np.ndarray] = {}
    for baseline in BASELINES:
        arrays[baseline] = (
            np.vstack(
                [
                    query_metric_vector(
                        [selected[baseline][query][target] for target in positives[query]]
                    )
                    for query in queries
                ]
            )
            if queries
            else np.empty((0, len(METRICS)), dtype=np.float64)
        )
    return queries, arrays


def parse_optional_float(value: str, label: str) -> float | None:
    if not value.strip():
        return None
    parsed = float(value)
    require(math.isfinite(parsed), f"{label} is nonfinite")
    return parsed


def load_primary_baseline(path: Path) -> dict[tuple[str, str, str], dict[str, str]]:
    output: dict[tuple[str, str, str], dict[str, str]] = {}
    with open_dict_reader(path) as reader:
        fields = set(reader.fieldnames or [])
        required = {"scope", "baseline", "metric", "query_count", "mean", "ci95_low", "ci95_high", "status"}
        require(required.issubset(fields), "Primary baseline-bootstrap table lacks required fields")
        for row in reader:
            key = (row["scope"], row["baseline"], row["metric"])
            require(key not in output, "Primary baseline-bootstrap table contains duplicate cells")
            output[key] = row
    expected = {(scope, baseline, metric) for scope in SCOPES for baseline in BASELINES for metric in METRICS}
    require(set(output) == expected and len(output) == 100, "Primary baseline-bootstrap matrix is not 100 planned cells")
    return output


def load_primary_focus(path: Path) -> dict[tuple[str, str], dict[str, str]]:
    output: dict[tuple[str, str], dict[str, str]] = {}
    with open_dict_reader(path) as reader:
        fields = set(reader.fieldnames or [])
        required = {
            "scope",
            "left_baseline",
            "right_baseline",
            "metric",
            "query_count",
            "mean_difference_left_minus_right",
            "ci95_low",
            "ci95_high",
            "status",
        }
        require(required.issubset(fields), "Primary focus-contrast table lacks required fields")
        for row in reader:
            require(row["left_baseline"] == FOCUS_LEFT and row["right_baseline"] == FOCUS_RIGHT, "Primary focus direction changed")
            key = (row["scope"], row["metric"])
            require(key not in output, "Primary focus table contains duplicate cells")
            output[key] = row
    expected = {(scope, metric) for scope in SCOPES for metric in METRICS}
    require(set(output) == expected and len(output) == 25, "Primary focus matrix is not 25 planned cells")
    return output


def component_bootstrap_means(
    *,
    queries: list[str],
    arrays: dict[str, np.ndarray],
    components: list[dict[str, list[str]]],
    rng: np.random.Generator,
) -> np.ndarray | None:
    if len(components) < 2:
        return None
    query_index = {query: index for index, query in enumerate(queries)}
    stacked = np.stack([arrays[baseline] for baseline in BASELINES], axis=0)
    component_sizes = np.asarray([len(component["queries"]) for component in components], dtype=np.int64)
    component_sums = np.stack(
        [
            stacked[:, [query_index[query] for query in component["queries"]], :].sum(axis=1)
            for component in components
        ],
        axis=0,
    )
    output = np.empty((BOOTSTRAP_REPLICATES, len(BASELINES), len(METRICS)), dtype=np.float64)
    component_count = len(components)
    for start in range(0, BOOTSTRAP_REPLICATES, BOOTSTRAP_BATCH_SIZE):
        stop = min(start + BOOTSTRAP_BATCH_SIZE, BOOTSTRAP_REPLICATES)
        draws = rng.integers(0, component_count, size=(stop - start, component_count), dtype=np.int32)
        denominators = component_sizes[draws].sum(axis=1)
        require(np.all(denominators > 0), "Component bootstrap produced an empty replicate")
        numerators = component_sums[draws].sum(axis=1)
        output[start:stop] = numerators / denominators[:, np.newaxis, np.newaxis]
    return output


def format_float(value: float | None) -> str:
    return "" if value is None else f"{float(value):.17g}"


def interval_comparison(
    primary_low: float | None,
    primary_high: float | None,
    component_low: float | None,
    component_high: float | None,
) -> dict[str, str]:
    primary_width = None
    component_width = None
    ratio = None
    status = "not_defined_missing_interval"
    if primary_low is not None and primary_high is not None:
        require(primary_high >= primary_low, "Primary interval endpoints are reversed")
        primary_width = primary_high - primary_low
    if component_low is not None and component_high is not None:
        require(component_high >= component_low, "Component interval endpoints are reversed")
        component_width = component_high - component_low
    if primary_width is not None and component_width is not None:
        if primary_width > 0.0:
            ratio = component_width / primary_width
            status = "defined"
        else:
            status = "not_defined_primary_interval_width_zero"
    return {
        "primary_query_interval_width": format_float(primary_width),
        "component_interval_width": format_float(component_width),
        "component_to_primary_interval_width_ratio": format_float(ratio),
        "interval_width_ratio_status": status,
    }


def peak_rss_bytes() -> tuple[int | None, str]:
    if os.name == "nt":
        try:
            import ctypes
            from ctypes import wintypes

            class Counters(ctypes.Structure):
                _fields_ = [
                    ("cb", wintypes.DWORD),
                    ("PageFaultCount", wintypes.DWORD),
                    ("PeakWorkingSetSize", ctypes.c_size_t),
                    ("WorkingSetSize", ctypes.c_size_t),
                    ("QuotaPeakPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                    ("QuotaNonPagedPoolUsage", ctypes.c_size_t),
                    ("PagefileUsage", ctypes.c_size_t),
                    ("PeakPagefileUsage", ctypes.c_size_t),
                ]

            kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
            psapi = ctypes.WinDLL("psapi", use_last_error=True)
            kernel32.GetCurrentProcess.restype = wintypes.HANDLE
            function = psapi.GetProcessMemoryInfo
            function.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
            function.restype = wintypes.BOOL
            counters = Counters()
            counters.cb = ctypes.sizeof(counters)
            if function(kernel32.GetCurrentProcess(), ctypes.byref(counters), counters.cb):
                return int(counters.PeakWorkingSetSize), "windows_GetProcessMemoryInfo_PeakWorkingSetSize"
        except (AttributeError, OSError, TypeError):
            pass
    return None, "not_available"


def main() -> int:
    args = parser().parse_args()
    total_started = time.perf_counter()
    tracemalloc.start()
    require_protocol_lock(args.protocol_lock)
    implementation_lock = verify_implementation_lock(args.implementation_lock)
    output_dir = require_new_output_dir(args.output_dir)

    input_paths = {
        "protocol_lock": args.protocol_lock,
        "implementation_lock": args.implementation_lock,
        "score_manifest": args.score_manifest,
        "evaluation_manifest": args.evaluation_manifest,
        "frozen_endpoint": args.endpoint,
        "future_source_evidence": args.source_evidence,
        "corrective_full_ranks": args.ranks,
        "scaffold_audit": args.scaffold_audit,
        "homology_0_30": args.homology_0_30,
        "homology_0_50": args.homology_0_50,
        "homology_0_70": args.homology_0_70,
        "primary_query_bootstrap_baselines": args.primary_baseline_bootstrap,
        "primary_query_bootstrap_focus_contrast": args.primary_focus_contrast,
    }
    inputs = [input_descriptor(role, path) for role, path in input_paths.items()]
    verify_corrected_manifests(
        score_manifest_path=args.score_manifest,
        evaluation_manifest_path=args.evaluation_manifest,
        ranks=args.ranks,
        endpoint=args.endpoint,
        scaffold=args.scaffold_audit,
        homology={"0_30": args.homology_0_30, "0_50": args.homology_0_50, "0_70": args.homology_0_70},
        primary_baseline=args.primary_baseline_bootstrap,
        primary_focus=args.primary_focus_contrast,
    )

    endpoint = load_endpoint(args.endpoint)
    endpoint_pair_keys = {row["canonical_pair_key"] for row in endpoint}
    pair_documents, source_input_audit = load_source_documents(args.source_evidence, endpoint_pair_keys)
    scaffold = load_bool_map(
        args.scaffold_audit,
        "canonical_pair_key",
        "audit_scaffold_cold_under_selected_policy",
        "audit_outcome",
    )
    homology = {
        "0_30": load_bool_map(
            args.homology_0_30,
            "uniprot_canonical_accession",
            "is_future_target_homology_cold_candidate",
            "future_target_coldness_status",
        ),
        "0_50": load_bool_map(
            args.homology_0_50,
            "uniprot_canonical_accession",
            "is_future_target_homology_cold_candidate",
            "future_target_coldness_status",
        ),
        "0_70": load_bool_map(
            args.homology_0_70,
            "uniprot_canonical_accession",
            "is_future_target_homology_cold_candidate",
            "future_target_coldness_status",
        ),
    }
    scoped_rows = build_scope_rows(endpoint, scaffold, homology)
    needed_targets: dict[str, set[str]] = defaultdict(set)
    endpoint_query_compounds: dict[str, str] = {}
    for row in endpoint:
        needed_targets[row["query_id"]].add(row["uniprot_canonical_accession"])
        endpoint_query_compounds[row["query_id"]] = row["inchikey_full"]
    selected = load_prediction_ranks(args.ranks, dict(needed_targets), endpoint_query_compounds)
    primary_baseline = load_primary_baseline(args.primary_baseline_bootstrap)
    primary_focus = load_primary_focus(args.primary_focus_contrast)

    rng = np.random.Generator(np.random.PCG64(BOOTSTRAP_SEED))
    scope_summary_rows: list[dict[str, Any]] = []
    baseline_output_rows: list[dict[str, Any]] = []
    focus_output_rows: list[dict[str, Any]] = []
    for scope in SCOPES:
        rows = scoped_rows[scope]
        queries, arrays = build_scope_arrays(rows, selected)
        components, document_count, edge_count = build_query_document_components(rows, pair_documents)
        query_sizes = [len(component["queries"]) for component in components]
        document_sizes = [len(component["documents"]) for component in components]
        query_distribution = distribution(query_sizes)
        document_distribution = distribution(document_sizes)
        component_status = (
            "estimable_descriptive_document_component_sensitivity"
            if len(components) >= 2
            else "not_estimable_component_count_lt_2"
        )
        scope_summary_rows.append(
            {
                "scope": scope,
                "relation_count": len(rows),
                "query_count": len(queries),
                "source_document_count": document_count,
                "query_source_document_edge_count": edge_count,
                "component_count": len(components),
                "component_query_size_min": query_distribution["min"],
                "component_query_size_q25": query_distribution["q25"],
                "component_query_size_median": query_distribution["median"],
                "component_query_size_q75": query_distribution["q75"],
                "component_query_size_max": query_distribution["max"],
                "component_query_size_mean": query_distribution["mean"],
                "component_source_document_size_min": document_distribution["min"],
                "component_source_document_size_median": document_distribution["median"],
                "component_source_document_size_max": document_distribution["max"],
                "largest_component_query_fraction": (
                    max(query_sizes) / len(queries) if queries else None
                ),
                "largest_component_source_document_fraction": (
                    max(document_sizes) / document_count if document_count else None
                ),
                "component_bootstrap_status": component_status,
            }
        )
        replicate_means = component_bootstrap_means(
            queries=queries, arrays=arrays, components=components, rng=rng
        )
        point_stack = (
            np.stack([arrays[baseline].mean(axis=0) for baseline in BASELINES], axis=0)
            if queries
            else None
        )
        if replicate_means is None:
            component_low = component_high = None
        else:
            component_low, component_high = np.percentile(replicate_means, [2.5, 97.5], axis=0)

        for baseline_index, baseline in enumerate(BASELINES):
            for metric_index, metric in enumerate(METRICS):
                primary = primary_baseline[(scope, baseline, metric)]
                require(int(primary["query_count"]) == len(queries), "Primary baseline query denominator changed")
                primary_mean = parse_optional_float(primary["mean"], "primary baseline mean")
                point = (
                    None
                    if point_stack is None
                    else float(point_stack[baseline_index, metric_index])
                )
                require(
                    (primary_mean is None and point is None)
                    or (
                        primary_mean is not None
                        and point is not None
                        and math.isclose(primary_mean, point, rel_tol=1e-12, abs_tol=1e-12)
                    ),
                    "Primary baseline point estimate was not preserved",
                )
                primary_low = parse_optional_float(primary["ci95_low"], "primary baseline CI low")
                primary_high = parse_optional_float(primary["ci95_high"], "primary baseline CI high")
                low = None if component_low is None else float(component_low[baseline_index, metric_index])
                high = None if component_high is None else float(component_high[baseline_index, metric_index])
                baseline_output_rows.append(
                    {
                        "scope": scope,
                        "baseline": baseline,
                        "metric": metric,
                        "relation_count": len(rows),
                        "query_count": len(queries),
                        "source_document_count": document_count,
                        "component_count": len(components),
                        "point_estimate": format_float(point),
                        "primary_query_bootstrap_ci95_low": format_float(primary_low),
                        "primary_query_bootstrap_ci95_high": format_float(primary_high),
                        "primary_query_bootstrap_status": primary["status"],
                        "component_bootstrap_ci95_low": format_float(low),
                        "component_bootstrap_ci95_high": format_float(high),
                        "component_bootstrap_status": component_status,
                        "component_bootstrap_replicates": (
                            BOOTSTRAP_REPLICATES if replicate_means is not None else 0
                        ),
                        **interval_comparison(primary_low, primary_high, low, high),
                    }
                )

        focus_difference = (
            point_stack[BASELINES.index(FOCUS_LEFT)]
            - point_stack[BASELINES.index(FOCUS_RIGHT)]
            if point_stack is not None
            else None
        )
        replicate_focus = None
        if replicate_means is not None:
            replicate_focus = (
                replicate_means[:, BASELINES.index(FOCUS_LEFT), :]
                - replicate_means[:, BASELINES.index(FOCUS_RIGHT), :]
            )
            focus_low, focus_high = np.percentile(replicate_focus, [2.5, 97.5], axis=0)
        else:
            focus_low = focus_high = None
        for metric_index, metric in enumerate(METRICS):
            primary = primary_focus[(scope, metric)]
            require(int(primary["query_count"]) == len(queries), "Primary focus query denominator changed")
            primary_mean = parse_optional_float(
                primary["mean_difference_left_minus_right"], "primary focus mean"
            )
            point = (
                None if focus_difference is None else float(focus_difference[metric_index])
            )
            require(
                (primary_mean is None and point is None)
                or (
                    primary_mean is not None
                    and point is not None
                    and math.isclose(primary_mean, point, rel_tol=1e-12, abs_tol=1e-12)
                ),
                "Primary focus point estimate was not preserved",
            )
            primary_low = parse_optional_float(primary["ci95_low"], "primary focus CI low")
            primary_high = parse_optional_float(primary["ci95_high"], "primary focus CI high")
            low = None if focus_low is None else float(focus_low[metric_index])
            high = None if focus_high is None else float(focus_high[metric_index])
            focus_output_rows.append(
                {
                    "scope": scope,
                    "left_baseline": FOCUS_LEFT,
                    "right_baseline": FOCUS_RIGHT,
                    "metric": metric,
                    "relation_count": len(rows),
                    "query_count": len(queries),
                    "source_document_count": document_count,
                    "component_count": len(components),
                    "point_difference_left_minus_right": format_float(point),
                    "primary_query_bootstrap_ci95_low": format_float(primary_low),
                    "primary_query_bootstrap_ci95_high": format_float(primary_high),
                    "primary_query_bootstrap_status": primary["status"],
                    "component_bootstrap_ci95_low": format_float(low),
                    "component_bootstrap_ci95_high": format_float(high),
                    "component_bootstrap_status": component_status,
                    "component_bootstrap_replicates": (
                        BOOTSTRAP_REPLICATES if replicate_means is not None else 0
                    ),
                    **interval_comparison(primary_low, primary_high, low, high),
                }
            )

    require(len(scope_summary_rows) == 5, "Scope summary row count is not five")
    require(len(baseline_output_rows) == 100, "Baseline component-bootstrap output is not 100 cells")
    require(len(focus_output_rows) == 25, "Focus component-bootstrap output is not 25 cells")

    current_bytes, peak_bytes = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    peak_rss, peak_rss_method = peak_rss_bytes()
    summary = {
        "audit_id": AUDIT_ID,
        "audit_suite_id": SUITE_ID,
        "protocol_id": PROTOCOL_ID,
        "created_at_utc": utc_now(),
        "dependence_graph": "scope-specific undirected query_compound--PMID bipartite connected components",
        "source_document_type": "PMID",
        "bootstrap": {
            "estimand": "equally_query_weighted_macro_mean",
            "unit": "scope_specific_query_source_document_connected_component",
            "replicates": BOOTSTRAP_REPLICATES,
            "prng": "PCG64",
            "seed": BOOTSTRAP_SEED,
            "scope_order": SCOPES,
            "draw_count_per_replicate": "observed_component_count",
            "interval": "95_percentile_2.5_97.5",
            "component_count_lt_2_rule": "not_estimable_component_count_lt_2",
        },
        "scope_summaries": scope_summary_rows,
        "baseline_metric_cell_count": len(baseline_output_rows),
        "focus_contrast_cell_count": len(focus_output_rows),
        "focus_contrast": {
            "left_baseline": FOCUS_LEFT,
            "right_baseline": FOCUS_RIGHT,
            "estimand": "left_minus_right",
        },
        "primary_query_bootstrap_point_estimates_verified_unchanged": True,
        "endpoint_contract_verified": {
            "relation_count": EXPECTED_ENDPOINT_RELATIONS,
            "query_count": EXPECTED_ENDPOINT_QUERIES,
            "target_count": EXPECTED_ENDPOINT_TARGETS,
        },
        "rank_row_contract_verified": EXPECTED_RANK_ROWS,
        "every_endpoint_relation_has_numeric_pmid": True,
        "aggregate_only": True,
        "identifier_bearing_output": False,
        "component_membership_output": False,
        "runtime": {
            "wall_seconds": time.perf_counter() - total_started,
            "python_tracemalloc_current_bytes": current_bytes,
            "python_tracemalloc_peak_bytes": peak_bytes,
            "process_peak_rss_bytes": peak_rss,
            "process_peak_rss_method": peak_rss_method,
        },
        "interpretation_boundary": (
            "This post hoc sensitivity diagnoses source-document dependence. It does not create "
            "document-disjoint validation, restore blindness, or replace the primary query bootstrap."
        ),
    }

    output_dir.mkdir(parents=True, exist_ok=False)
    scope_name = "document_component_scope_summary.tsv"
    metric_name = "document_component_bootstrap_metrics.tsv"
    focus_name = "document_component_bootstrap_focus_contrast.tsv"
    summary_name = "document_component_bootstrap_summary.json"
    write_tsv_new(output_dir / scope_name, list(scope_summary_rows[0]), scope_summary_rows)
    write_tsv_new(output_dir / metric_name, list(baseline_output_rows[0]), baseline_output_rows)
    write_tsv_new(output_dir / focus_name, list(focus_output_rows[0]), focus_output_rows)
    write_json_new(output_dir / summary_name, summary)
    manifest = finalize_manifest(
        output_dir=output_dir,
        audit_id=AUDIT_ID,
        script_path=Path(__file__),
        inputs=inputs,
        output_names=[scope_name, metric_name, focus_name, summary_name],
        extra={
            "implementation_lock": {
                "basename": args.implementation_lock.name,
                "sha256": sha256(args.implementation_lock),
                "amendment_id": implementation_lock.get("amendment_id"),
            },
            "source_evidence_input_audit": source_input_audit,
            "output_contract": {
                "scope_summary_rows": 5,
                "baseline_metric_cells": 100,
                "focus_contrast_cells": 25,
            },
            "component_membership_written": False,
            "primary_query_bootstrap_carried_through": True,
        },
    )
    write_json_new(output_dir / "run_manifest.json", manifest)
    print(
        f"{AUDIT_ID}: wrote 5 scope rows, 100 baseline cells, and 25 focus cells; "
        "no identifiers or component memberships were written"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

