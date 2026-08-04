"""Run the frozen aggregate-only exact tie-aware retrieval sensitivity."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import platform
import sys
import time
import tracemalloc
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
PARENT_PROTOCOL_SHA256 = "bd4902476160cc7c5cbacaf0cfd0f1a28c5300bde22232b3a2cc6c1f3c143dc2"
ANALYSIS_ID = "revision_tie_aware_v1_20260729"
PROTOCOL_ID = "npass_strict_ab_major_revision_v4_20260729"
EXPECTED_RANK_ROWS = 3_658_128
EXPECTED_ENDPOINT_RELATIONS = 358
EXPECTED_ENDPOINT_QUERIES = 222
EXPECTED_ENDPOINT_TARGETS = 156
BOOTSTRAP_REPLICATES = 10_000
QUERY_BOOTSTRAP_BASE_SEED = 2_026_072_901
COMPONENT_BOOTSTRAP_BASE_SEED = 2_026_072_902
IDENTIFIABLE_TOLERANCE = 1e-15
ALPHA = 0.05

BASELINES = [
    "weighted_target_popularity",
    "sequence_3mer_transfer",
    "weighted_morgan_transfer",
    "structure_sequence_pair_neighbor",
]
INTERNAL_SCOPES = [
    "temporal_strict_ab",
    "scaffold_cold_strict_ab",
    "double_cold_0_30",
    "double_cold_0_50",
    "double_cold_0_70",
]
DISPLAY_SCOPE_TO_INTERNAL = {
    "temporal_strict_ab": "temporal_strict_ab",
    "scaffold_cold_strict_ab": "scaffold_cold_strict_ab",
    "project_defined_joint_scaffold_homology_cold_0_30": "double_cold_0_30",
    "project_defined_joint_scaffold_homology_cold_0_50_0_70_identical_mask": "double_cold_0_50",
}
DISPLAY_SCOPES = list(DISPLAY_SCOPE_TO_INTERNAL)
K_VALUES = (10, 50)

EXPECTED_INPUTS = {
    "ranks": (
        "corrective_prediction_ranks.tsv.gz",
        "87739aa818744c7084088d13c386444aa41bbef38c257083325298003181479e",
    ),
    "endpoint": (
        "evaluation_pairs.tsv.gz",
        "09296b066a23197a7c178f00514f2b3d9ed7e6f3c459ea92a55e01a6010d1132",
    ),
    "scaffold": (
        "scaffold_audit.tsv.gz",
        "fa0029ef5b7822ad5ca93f7bd93ac808f85f1e0c02e827fa91be375031b2d7af",
    ),
    "homology_0_30": (
        "homology_0_30.tsv.gz",
        "3a8247ed8f683fe6fce5fb345f56e3ec73a872b065eca922e92e494f084a1793",
    ),
    "homology_0_50": (
        "homology_0_50.tsv.gz",
        "ec3bbd435f73bc1c724efdfd94ac10c32d6f9a55bd6c9a9349233a01e5dc7f5b",
    ),
    "homology_0_70": (
        "homology_0_70.tsv.gz",
        "ec3bbd435f73bc1c724efdfd94ac10c32d6f9a55bd6c9a9349233a01e5dc7f5b",
    ),
    "source_evidence": (
        "future_active_source_strict_primary_evidence.tsv.gz",
        "e800a181976de63a7ee027ed929576ae16fd925b0f458bff414d41075cefa21f",
    ),
}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json_new(path: Path, payload: Any) -> None:
    require(not path.exists(), f"Refusing to overwrite {path.name}")
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


def write_tsv_new(path: Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    require(not path.exists(), f"Refusing to overwrite {path.name}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def open_dict_reader(path: Path) -> csv.DictReader:
    if path.suffix == ".gz":
        handle = gzip.open(path, "rt", encoding="utf-8", newline="")
    else:
        handle = path.open("r", encoding="utf-8", newline="")
    return csv.DictReader(handle, delimiter="\t")


def parse_bool(value: str) -> bool:
    normalized = value.strip().lower()
    require(normalized in {"true", "false"}, f"Invalid boolean value: {value!r}")
    return normalized == "true"


def input_descriptor(role: str, path: Path) -> dict[str, Any]:
    return {
        "role": role,
        "basename": path.name,
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def output_descriptor(path: Path) -> dict[str, Any]:
    return {
        "basename": path.name,
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def format_float(value: float | None) -> str:
    return "" if value is None else f"{float(value):.17g}"


def derived_seed(base_seed: int, *labels: str) -> int:
    material = "|".join([str(base_seed), *labels]).encode("utf-8")
    return int.from_bytes(hashlib.sha256(material).digest()[:8], "big", signed=False)


def verify_locks(
    plan_lock: Path, parent_protocol: Path, implementation_lock_path: Path
) -> dict[str, Any]:
    require(plan_lock.is_file(), "Local analysis-plan lock is absent")
    require(parent_protocol.is_file(), "Parent major-revision protocol is absent")
    require(implementation_lock_path.is_file(), "Implementation lock is absent")
    require(sha256(parent_protocol) == PARENT_PROTOCOL_SHA256, "Parent protocol hash changed")
    plan = json.loads(plan_lock.read_text(encoding="utf-8"))
    require(plan.get("analysis_id") == ANALYSIS_ID, "Local plan analysis ID changed")
    require(
        plan.get("lock_state") == "LOCKED_BEFORE_TIE_AWARE_RESULT_COMPUTATION",
        "Local plan was not frozen before result computation",
    )
    implementation = json.loads(implementation_lock_path.read_text(encoding="utf-8"))
    require(
        implementation.get("lock_state") == "LOCKED_BEFORE_REAL_INPUT_EXECUTION",
        "Implementation was not locked before real-input execution",
    )
    require(
        implementation.get("parent_protocol_sha256") == PARENT_PROTOCOL_SHA256,
        "Implementation lock parent-protocol anchor changed",
    )
    require(
        implementation.get("local_plan_sha256") == sha256(plan_lock),
        "Implementation lock local-plan anchor changed",
    )
    files = implementation.get("implementation_files")
    require(isinstance(files, dict) and files, "Implementation lock lacks its file inventory")
    for relative_name, expected_hash in sorted(files.items()):
        candidate = (ROOT / relative_name).resolve()
        try:
            candidate.relative_to(ROOT.resolve())
        except ValueError as exc:
            raise ValueError("Implementation file escapes analysis root") from exc
        require(candidate.is_file(), f"Locked implementation file is absent: {relative_name}")
        require(sha256(candidate) == expected_hash, f"Locked implementation drifted: {relative_name}")
    require(
        files.get("scripts/run_tie_aware_analysis_v1.py") == sha256(Path(__file__).resolve()),
        "Executing analysis script differs from the implementation lock",
    )
    return implementation


def verify_scientific_inputs(paths: dict[str, Path]) -> list[dict[str, Any]]:
    descriptors: list[dict[str, Any]] = []
    for role, path in paths.items():
        require(path.is_file(), f"Input is absent: {role}")
        expected_name, expected_hash = EXPECTED_INPUTS[role]
        require(path.name == expected_name, f"Unexpected basename for {role}")
        actual_hash = sha256(path)
        require(actual_hash == expected_hash, f"Frozen input hash changed: {role}")
        descriptors.append(input_descriptor(role, path))
    return descriptors


def load_endpoint(path: Path) -> tuple[list[dict[str, str]], dict[str, str]]:
    rows: list[dict[str, str]] = []
    query_compounds: dict[str, str] = {}
    seen_pairs: set[str] = set()
    seen_compound_target: set[tuple[str, str]] = set()
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "canonical_pair_key",
            "query_id",
            "inchikey_full",
            "uniprot_canonical_accession",
            "best_strict_evidence_tier",
            "decision",
            "c31_leakage_gate_status",
        }
        require(required.issubset(set(reader.fieldnames or [])), "Endpoint fields are incomplete")
        for row in reader:
            pair = row["canonical_pair_key"]
            compound_target = (row["inchikey_full"], row["uniprot_canonical_accession"])
            require(pair and pair not in seen_pairs, "Endpoint pair key is empty or duplicated")
            require(
                compound_target not in seen_compound_target,
                "Endpoint compound-target relation is duplicated",
            )
            seen_pairs.add(pair)
            seen_compound_target.add(compound_target)
            require(
                row["decision"] == "strict_post_cutoff_future_candidate",
                "Endpoint decision changed",
            )
            require(
                row["c31_leakage_gate_status"] == "pass_no_historical_activity",
                "Endpoint leakage-gate state changed",
            )
            require(
                row["best_strict_evidence_tier"]
                in {"A_affinity_candidate", "B_quantitative_functional_candidate"},
                "Endpoint contains a non-strict tier",
            )
            query = row["query_id"]
            compound = row["inchikey_full"]
            require(
                query not in query_compounds or query_compounds[query] == compound,
                "Endpoint query maps to multiple compounds",
            )
            query_compounds[query] = compound
            rows.append(row)
    require(len(rows) == EXPECTED_ENDPOINT_RELATIONS, "Endpoint relation count changed")
    require(len(query_compounds) == EXPECTED_ENDPOINT_QUERIES, "Endpoint query count changed")
    require(
        len({row["uniprot_canonical_accession"] for row in rows})
        == EXPECTED_ENDPOINT_TARGETS,
        "Endpoint target count changed",
    )
    return rows, query_compounds


def load_bool_map(path: Path, key_field: str, flag_field: str, status_field: str) -> dict[str, bool]:
    output: dict[str, bool] = {}
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(
            {key_field, flag_field, status_field}.issubset(set(reader.fieldnames or [])),
            f"Scope input {path.name} lacks required fields",
        )
        for row in reader:
            key = row[key_field].strip()
            require(key and key not in output, f"Scope input {path.name} has a duplicate key")
            require(row[status_field].strip(), f"Scope input {path.name} has an empty status")
            output[key] = parse_bool(row[flag_field])
    require(output, f"Scope input {path.name} is empty")
    return output


def build_internal_scopes(
    endpoint: list[dict[str, str]],
    scaffold: dict[str, bool],
    homology: dict[str, dict[str, bool]],
) -> dict[str, list[dict[str, str]]]:
    pair_keys = {row["canonical_pair_key"] for row in endpoint}
    targets = {row["uniprot_canonical_accession"] for row in endpoint}
    require(set(scaffold) == pair_keys, "Scaffold keyset differs from endpoint")
    for threshold, flags in homology.items():
        require(set(flags) == targets, f"Homology {threshold} keyset differs from endpoint")
    output = {scope: [] for scope in INTERNAL_SCOPES}
    for row in endpoint:
        output["temporal_strict_ab"].append(row)
        if scaffold[row["canonical_pair_key"]]:
            output["scaffold_cold_strict_ab"].append(row)
            target = row["uniprot_canonical_accession"]
            for threshold in ("0_30", "0_50", "0_70"):
                if homology[threshold][target]:
                    output[f"double_cold_{threshold}"].append(row)
    mask_050 = {row["canonical_pair_key"] for row in output["double_cold_0_50"]}
    mask_070 = {row["canonical_pair_key"] for row in output["double_cold_0_70"]}
    require(mask_050 == mask_070, "0.50 and 0.70 scope masks are not identical")
    return output


def load_pair_documents(
    path: Path, endpoint_pairs: set[str]
) -> tuple[dict[str, set[str]], dict[str, int]]:
    pair_documents: dict[str, set[str]] = defaultdict(set)
    total_rows = 0
    nonendpoint_rows = 0
    ineligible_rows = 0
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(
            {"canonical_pair_key", "ref_id_type", "ref_id"}.issubset(
                set(reader.fieldnames or [])
            ),
            "Source-evidence fields are incomplete",
        )
        for row in reader:
            total_rows += 1
            pair = row["canonical_pair_key"].strip()
            if pair not in endpoint_pairs:
                nonendpoint_rows += 1
                continue
            reference = row["ref_id"].strip()
            if row["ref_id_type"].strip().upper() != "PMID" or not reference.isdigit():
                ineligible_rows += 1
                continue
            pair_documents[pair].add(reference)
    require(set(pair_documents) == endpoint_pairs, "An endpoint relation lacks a numeric PMID")
    require(total_rows == 486, "Source-evidence row count changed")
    require(nonendpoint_rows == 0 and ineligible_rows == 0, "Source-evidence eligibility changed")
    return dict(pair_documents), {
        "source_evidence_row_count": total_rows,
        "endpoint_relation_count_with_numeric_pmid": len(pair_documents),
        "nonendpoint_row_count": nonendpoint_rows,
        "ineligible_reference_row_count": ineligible_rows,
    }


def make_relevance(
    display_rows: dict[str, list[dict[str, str]]],
) -> tuple[dict[str, dict[str, set[str]]], dict[str, dict[str, set[str]]]]:
    targets: dict[str, dict[str, set[str]]] = {}
    pairs: dict[str, dict[str, set[str]]] = {}
    for scope, rows in display_rows.items():
        by_query_targets: dict[str, set[str]] = defaultdict(set)
        by_query_pairs: dict[str, set[str]] = defaultdict(set)
        for row in rows:
            query = row["query_id"]
            by_query_targets[query].add(row["uniprot_canonical_accession"])
            by_query_pairs[query].add(row["canonical_pair_key"])
        targets[scope] = dict(by_query_targets)
        pairs[scope] = dict(by_query_pairs)
    return targets, pairs


def discounts(start: int, stop: int) -> list[float]:
    return [1.0 / math.log2(rank + 1.0) for rank in range(start, stop + 1)]


def probability_no_relevant(m: int, r: int, s: int) -> float:
    require(0 <= r <= m and 0 <= s <= m, "Invalid hypergeometric arguments")
    if s == 0:
        return 1.0
    if m - r < s:
        return 0.0
    probability = 1.0
    for index in range(s):
        probability *= (m - r - index) / (m - index)
    return probability


def tie_query_metrics(
    blocks: list[tuple[int, int, int, int]],
    relevant_count: int,
    salted_relevant_ranks: list[int],
    k: int,
) -> dict[str, float | int | bool]:
    """Return exact uniform-tie expectations and attainable bounds for one query."""
    require(relevant_count > 0, "A query metric requires at least one relevant target")
    require(len(salted_relevant_ranks) == relevant_count, "Salted relevant ranks are incomplete")
    expected_relevant = 0.0
    worst_relevant = 0
    best_relevant = 0
    expected_dcg = 0.0
    worst_dcg = 0.0
    best_dcg = 0.0
    membership_identifiable = 0
    membership_boundary = 0
    membership_not_retrieved = 0
    guaranteed_hit = False
    possible_hit = False
    expected_hit_probability = 0.0
    boundary_hit_probability: float | None = None

    for start, end, block_size, relevant_in_block in blocks:
        require(end - start + 1 == block_size, "Tie-block positions are inconsistent")
        require(relevant_in_block > 0, "Stored tie block has no relevant member")
        evaluated_stop = min(end, k)
        s = max(0, evaluated_stop - start + 1)
        if s:
            expected_relevant += relevant_in_block * s / block_size
            forced = max(0, relevant_in_block - (block_size - s))
            possible = min(relevant_in_block, s)
            worst_relevant += forced
            best_relevant += possible
            block_discounts = discounts(start, evaluated_stop)
            expected_dcg += (relevant_in_block / block_size) * sum(block_discounts)
            if forced:
                worst_dcg += sum(block_discounts[-forced:])
            if possible:
                best_dcg += sum(block_discounts[:possible])
        if end <= k:
            membership_identifiable += relevant_in_block
            guaranteed_hit = True
            possible_hit = True
        elif start <= k < end:
            membership_boundary += relevant_in_block
            possible_hit = True
            boundary_hit_probability = 1.0 - probability_no_relevant(
                block_size, relevant_in_block, k - start + 1
            )
            if max(0, relevant_in_block - (block_size - (k - start + 1))) > 0:
                guaranteed_hit = True
        else:
            membership_not_retrieved += relevant_in_block

    require(
        membership_identifiable + membership_boundary + membership_not_retrieved
        == relevant_count,
        "Membership classification does not cover every relevant relation",
    )
    if guaranteed_hit:
        expected_hit_probability = 1.0
    elif boundary_hit_probability is not None:
        expected_hit_probability = boundary_hit_probability
    else:
        expected_hit_probability = 0.0
    ideal_dcg = sum(discounts(1, min(k, relevant_count)))
    salted_hits = sum(rank <= k for rank in salted_relevant_ranks)
    salted_dcg = sum(
        1.0 / math.log2(rank + 1.0) for rank in salted_relevant_ranks if rank <= k
    )
    result: dict[str, float | int | bool] = {
        "relevant_count": relevant_count,
        "legacy_salted_recall": salted_hits / relevant_count,
        "tie_expected_recall": expected_relevant / relevant_count,
        "tie_worst_recall": worst_relevant / relevant_count,
        "tie_best_recall": best_relevant / relevant_count,
        "legacy_salted_ndcg": salted_dcg / ideal_dcg,
        "tie_expected_ndcg": expected_dcg / ideal_dcg,
        "tie_worst_ndcg": worst_dcg / ideal_dcg,
        "tie_best_ndcg": best_dcg / ideal_dcg,
        "legacy_salted_any_hit": bool(salted_hits),
        "tie_expected_any_hit_probability": expected_hit_probability,
        "score_guaranteed_any_hit": guaranteed_hit,
        "score_possible_any_hit": possible_hit,
        "membership_score_identifiable": membership_identifiable,
        "membership_boundary_tie_dependent": membership_boundary,
        "membership_not_retrieved": membership_not_retrieved,
    }
    for metric in ("recall", "ndcg"):
        lower = float(result[f"tie_worst_{metric}"])
        expected = float(result[f"tie_expected_{metric}"])
        upper = float(result[f"tie_best_{metric}"])
        require(
            lower - IDENTIFIABLE_TOLERANCE
            <= expected
            <= upper + IDENTIFIABLE_TOLERANCE,
            f"Tie-aware {metric} expectation lies outside attainable bounds",
        )
        result[f"{metric}_score_identifiable"] = (
            abs(upper - lower) <= IDENTIFIABLE_TOLERANCE
        )
    return result


class UnionFind:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

    def add(self, item: str) -> None:
        if item not in self.parent:
            self.parent[item] = item

    def find(self, item: str) -> str:
        parent = self.parent[item]
        if parent != item:
            self.parent[item] = self.find(parent)
        return self.parent[item]

    def union(self, left: str, right: str) -> None:
        self.add(left)
        self.add(right)
        root_left = self.find(left)
        root_right = self.find(right)
        if root_left != root_right:
            self.parent[root_right] = root_left


def build_components(
    selected_queries: list[str],
    query_pairs: dict[str, set[str]],
    pair_documents: dict[str, set[str]],
) -> tuple[list[list[str]], int, int]:
    query_documents: dict[str, set[str]] = {}
    for query in selected_queries:
        documents: set[str] = set()
        for pair in query_pairs[query]:
            documents.update(pair_documents[pair])
        require(documents, "A selected query has no numeric PMID")
        query_documents[query] = documents
    union = UnionFind()
    for query, documents in query_documents.items():
        query_node = f"q:{query}"
        union.add(query_node)
        for document in documents:
            union.union(query_node, f"d:{document}")
    component_queries: dict[str, list[str]] = defaultdict(list)
    for query in selected_queries:
        component_queries[union.find(f"q:{query}")].append(query)
    components = [sorted(items) for items in component_queries.values()]
    components.sort(key=lambda items: tuple(items))
    require(
        sum(len(component) for component in components) == len(selected_queries),
        "Component construction did not partition selected queries",
    )
    all_documents = set().union(*query_documents.values()) if query_documents else set()
    edge_count = sum(len(documents) for documents in query_documents.values())
    return components, len(all_documents), edge_count


def query_bootstrap_interval(
    values: np.ndarray, seed: int
) -> tuple[np.ndarray | None, np.ndarray | None, str]:
    require(values.ndim == 2 and values.shape[1] == 2, "Query-bootstrap array is invalid")
    n_queries = values.shape[0]
    if n_queries == 0:
        return None, None, "not_estimable_no_queries"
    point = values.mean(axis=0)
    if n_queries == 1:
        return point, point, "n=1_descriptive_point_only"
    rng = np.random.Generator(np.random.PCG64(seed))
    draws = rng.integers(
        0, n_queries, size=(BOOTSTRAP_REPLICATES, n_queries), dtype=np.int32
    )
    replicate_means = values[draws].mean(axis=1)
    low, high = np.percentile(replicate_means, [2.5, 97.5], axis=0, method="linear")
    return low, high, "estimable_descriptive_query_bootstrap"


def component_bootstrap_interval(
    values: np.ndarray,
    selected_queries: list[str],
    components: list[list[str]],
    seed: int,
) -> tuple[np.ndarray | None, np.ndarray | None, str]:
    require(values.ndim == 2 and values.shape[1] == 2, "Component array is invalid")
    if len(components) < 2:
        return None, None, "not_estimable_component_count_lt_2"
    query_index = {query: index for index, query in enumerate(selected_queries)}
    component_sizes = np.asarray([len(component) for component in components], dtype=np.int64)
    component_sums = np.vstack(
        [
            values[[query_index[query] for query in component]].sum(axis=0)
            for component in components
        ]
    )
    rng = np.random.Generator(np.random.PCG64(seed))
    output = np.empty((BOOTSTRAP_REPLICATES, 2), dtype=np.float64)
    component_count = len(components)
    batch_size = 1_000
    for start in range(0, BOOTSTRAP_REPLICATES, batch_size):
        stop = min(start + batch_size, BOOTSTRAP_REPLICATES)
        draws = rng.integers(
            0,
            component_count,
            size=(stop - start, component_count),
            dtype=np.int32,
        )
        denominators = component_sizes[draws].sum(axis=1)
        numerators = component_sums[draws].sum(axis=1)
        output[start:stop] = numerators / denominators[:, np.newaxis]
    low, high = np.percentile(output, [2.5, 97.5], axis=0, method="linear")
    return low, high, "estimable_descriptive_pmid_component_sensitivity"


def binomial_cdf(x: int, n: int, probability: float) -> float:
    return sum(
        math.comb(n, index)
        * probability**index
        * (1.0 - probability) ** (n - index)
        for index in range(x + 1)
    )


def clopper_pearson_upper(x: int, n: int, alpha: float = ALPHA) -> float:
    require(n > 0 and 0 <= x <= n and 0 < alpha < 1, "Invalid binomial-bound inputs")
    if x == n:
        return 1.0
    if x == 0:
        return 1.0 - alpha ** (1.0 / n)
    low = x / n
    high = 1.0
    for _ in range(200):
        midpoint = (low + high) / 2.0
        if binomial_cdf(x, n, midpoint) > alpha:
            low = midpoint
        else:
            high = midpoint
    return (low + high) / 2.0


def load_rank_metrics(
    path: Path,
    endpoint_query_compounds: dict[str, str],
    relevance: dict[str, dict[str, set[str]]],
) -> tuple[
    dict[tuple[str, str, str, int], dict[str, float | int | bool]],
    dict[str, bool],
    dict[str, int],
]:
    records: dict[tuple[str, str, str, int], dict[str, float | int | bool]] = {}
    sequence_operational: dict[str, bool] = {}
    candidate_count_by_query: dict[str, int] = {}
    seen_groups: set[tuple[str, str]] = set()
    total_rows = 0
    current_key: tuple[str, str] | None = None
    current_compound = ""
    current_candidate_count = -1
    current_rows: list[tuple[str, int, float]] = []

    def finish_group() -> None:
        nonlocal current_key, current_compound, current_candidate_count, current_rows
        if current_key is None:
            return
        query, baseline = current_key
        require(current_key not in seen_groups, "A rank group recurs non-contiguously")
        seen_groups.add(current_key)
        require(
            len(current_rows) == current_candidate_count,
            "Rank-group row count differs from candidate count",
        )
        targets = [row[0] for row in current_rows]
        ranks = [row[1] for row in current_rows]
        require(len(set(targets)) == current_candidate_count, "Rank group has duplicate targets")
        require(set(ranks) == set(range(1, current_candidate_count + 1)), "Ranks are not 1..N")
        ordered_by_rank = sorted(current_rows, key=lambda row: row[1])
        require(
            all(
                ordered_by_rank[index - 1][2] >= ordered_by_rank[index][2]
                for index in range(1, len(ordered_by_rank))
            ),
            "Salted ranks violate descending score order",
        )
        score_groups: dict[float, list[tuple[str, int]]] = defaultdict(list)
        for target, rank, score in current_rows:
            score_groups[score].append((target, rank))
        block_layout: list[tuple[int, int, list[tuple[str, int]]]] = []
        start = 1
        for score in sorted(score_groups, reverse=True):
            members = score_groups[score]
            end = start + len(members) - 1
            member_ranks = [rank for _, rank in members]
            require(
                min(member_ranks) == start
                and max(member_ranks) == end
                and len(set(member_ranks)) == len(members),
                "Exact-score block is not a contiguous salted-rank interval",
            )
            block_layout.append((start, end, members))
            start = end + 1
        require(start - 1 == current_candidate_count, "Score blocks do not cover the rank group")
        if baseline == "sequence_3mer_transfer":
            sequence_operational[query] = any(score > 0.0 for _, _, score in current_rows)
        rank_by_target = {target: rank for target, rank, _ in current_rows}
        for display_scope in DISPLAY_SCOPES:
            relevant_targets = relevance[display_scope].get(query)
            if not relevant_targets:
                continue
            require(
                relevant_targets.issubset(rank_by_target),
                "A relevant target is absent from the eligible rank group",
            )
            blocks: list[tuple[int, int, int, int]] = []
            for block_start, block_end, members in block_layout:
                relevant_in_block = sum(target in relevant_targets for target, _ in members)
                if relevant_in_block:
                    blocks.append(
                        (
                            block_start,
                            block_end,
                            len(members),
                            relevant_in_block,
                        )
                    )
            salted_relevant_ranks = [rank_by_target[target] for target in relevant_targets]
            for k in K_VALUES:
                key = (display_scope, baseline, query, k)
                require(key not in records, "A query metric record is duplicated")
                records[key] = tie_query_metrics(
                    blocks, len(relevant_targets), salted_relevant_ranks, k
                )
        current_key = None
        current_compound = ""
        current_candidate_count = -1
        current_rows = []

    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
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
        require(required.issubset(set(reader.fieldnames or [])), "Rank fields are incomplete")
        for row in reader:
            total_rows += 1
            baseline = row["baseline"]
            query = row["query_id"]
            require(baseline in BASELINES, "Rank file contains an unknown baseline")
            require(query in endpoint_query_compounds, "Rank file contains an unknown query")
            key = (query, baseline)
            candidate_count = int(row["eligible_candidate_target_count"])
            if key != current_key:
                finish_group()
                current_key = key
                current_compound = row["query_compound_inchikey_full"]
                current_candidate_count = candidate_count
            require(
                row["query_compound_inchikey_full"] == current_compound
                == endpoint_query_compounds[query],
                "Rank query/compound mapping differs from endpoint",
            )
            require(
                candidate_count == current_candidate_count,
                "Candidate count changes within a rank group",
            )
            require(
                query not in candidate_count_by_query
                or candidate_count_by_query[query] == candidate_count,
                "Candidate count changes across baselines",
            )
            candidate_count_by_query[query] = candidate_count
            rank = int(row["rank"])
            score = float(row["score"])
            require(math.isfinite(score), "Rank file contains a nonfinite score")
            require(1 <= rank <= candidate_count, "Rank is outside the candidate range")
            current_rows.append((row["target_uniprot_accession"], rank, score))
    finish_group()
    require(total_rows == EXPECTED_RANK_ROWS, "Complete-rank row count changed")
    expected_groups = {
        (query, baseline)
        for query in endpoint_query_compounds
        for baseline in BASELINES
    }
    require(seen_groups == expected_groups, "Rank baseline/query groups are incomplete")
    require(
        len(sequence_operational) == EXPECTED_ENDPOINT_QUERIES,
        "3-mer operability is incomplete",
    )
    require(
        sum(sequence_operational.values()) == 60,
        "Frozen 3-mer operational-query count is not 60",
    )
    require(
        len(sequence_operational) - sum(sequence_operational.values()) == 162,
        "Frozen 3-mer structural-all-zero count is not 162",
    )
    return records, sequence_operational, candidate_count_by_query


def selected_queries_for_subset(
    scope_queries: Iterable[str],
    baseline: str,
    subset: str,
    sequence_operational: dict[str, bool],
) -> list[str]:
    queries = sorted(scope_queries)
    if subset == "all_queries":
        return queries
    require(baseline == "sequence_3mer_transfer", "3-mer subset used for another baseline")
    if subset == "score_operational":
        return [query for query in queries if sequence_operational[query]]
    require(
        subset == "structural_all_zero_non_operational",
        "Unknown query subset",
    )
    return [query for query in queries if not sequence_operational[query]]


def aggregate_metrics(
    records: dict[tuple[str, str, str, int], dict[str, float | int | bool]],
    relevance: dict[str, dict[str, set[str]]],
    query_pairs: dict[str, dict[str, set[str]]],
    pair_documents: dict[str, set[str]],
    sequence_operational: dict[str, bool],
) -> tuple[list[dict[str, Any]], dict[tuple[str, str, str, int], dict[str, Any]]]:
    output_rows: list[dict[str, Any]] = []
    index: dict[tuple[str, str, str, int], dict[str, Any]] = {}
    for scope in DISPLAY_SCOPES:
        scope_queries = relevance[scope]
        for baseline in BASELINES:
            subsets = ["all_queries"]
            if baseline == "sequence_3mer_transfer":
                subsets.extend(
                    ["score_operational", "structural_all_zero_non_operational"]
                )
            for subset in subsets:
                selected_queries = selected_queries_for_subset(
                    scope_queries, baseline, subset, sequence_operational
                )
                components, document_count, edge_count = build_components(
                    selected_queries, query_pairs[scope], pair_documents
                )
                for k in K_VALUES:
                    query_records = [
                        records[(scope, baseline, query, k)] for query in selected_queries
                    ]
                    values = (
                        np.asarray(
                            [
                                [
                                    float(row["tie_expected_recall"]),
                                    float(row["tie_expected_ndcg"]),
                                ]
                                for row in query_records
                            ],
                            dtype=np.float64,
                        )
                        if query_records
                        else np.empty((0, 2), dtype=np.float64)
                    )

                    def record_mean(field: str) -> float | None:
                        if not query_records:
                            return None
                        return float(np.mean([float(item[field]) for item in query_records]))

                    query_low, query_high, query_status = query_bootstrap_interval(
                        values,
                        derived_seed(
                            QUERY_BOOTSTRAP_BASE_SEED, scope, baseline, subset
                        ),
                    )
                    component_low, component_high, component_status = (
                        component_bootstrap_interval(
                            values,
                            selected_queries,
                            components,
                            derived_seed(
                                COMPONENT_BOOTSTRAP_BASE_SEED,
                                scope,
                                baseline,
                                subset,
                            ),
                        )
                    )
                    row: dict[str, Any] = {
                        "scope": scope,
                        "baseline": baseline,
                        "query_subset": subset,
                        "k": k,
                        "query_count": len(selected_queries),
                        "relevant_relation_count": sum(
                            int(item["relevant_count"]) for item in query_records
                        ),
                        "legacy_salted_recall": format_float(
                            record_mean("legacy_salted_recall")
                        ),
                        "tie_expected_fractional_recall": format_float(
                            None if not query_records else float(values[:, 0].mean())
                        ),
                        "tie_worst_recall": format_float(record_mean("tie_worst_recall")),
                        "tie_best_recall": format_float(record_mean("tie_best_recall")),
                        "query_bootstrap_expected_recall_ci95_low": format_float(
                            None if query_low is None else float(query_low[0])
                        ),
                        "query_bootstrap_expected_recall_ci95_high": format_float(
                            None if query_high is None else float(query_high[0])
                        ),
                        "query_bootstrap_status": query_status,
                        "pmid_component_expected_recall_ci95_low": format_float(
                            None if component_low is None else float(component_low[0])
                        ),
                        "pmid_component_expected_recall_ci95_high": format_float(
                            None if component_high is None else float(component_high[0])
                        ),
                        "legacy_salted_ndcg": format_float(record_mean("legacy_salted_ndcg")),
                        "tie_expected_ndcg": format_float(
                            None if not query_records else float(values[:, 1].mean())
                        ),
                        "tie_worst_ndcg": format_float(record_mean("tie_worst_ndcg")),
                        "tie_best_ndcg": format_float(record_mean("tie_best_ndcg")),
                        "query_bootstrap_expected_ndcg_ci95_low": format_float(
                            None if query_low is None else float(query_low[1])
                        ),
                        "query_bootstrap_expected_ndcg_ci95_high": format_float(
                            None if query_high is None else float(query_high[1])
                        ),
                        "pmid_component_expected_ndcg_ci95_low": format_float(
                            None if component_low is None else float(component_low[1])
                        ),
                        "pmid_component_expected_ndcg_ci95_high": format_float(
                            None if component_high is None else float(component_high[1])
                        ),
                        "pmid_component_bootstrap_status": component_status,
                        "membership_score_identifiable_relation_count": sum(
                            int(item["membership_score_identifiable"])
                            for item in query_records
                        ),
                        "membership_boundary_tie_dependent_relation_count": sum(
                            int(item["membership_boundary_tie_dependent"])
                            for item in query_records
                        ),
                        "membership_not_retrieved_relation_count": sum(
                            int(item["membership_not_retrieved"])
                            for item in query_records
                        ),
                        "recall_score_identifiable_query_count": sum(
                            bool(item["recall_score_identifiable"])
                            for item in query_records
                        ),
                        "recall_tie_dependent_query_count": sum(
                            not bool(item["recall_score_identifiable"])
                            for item in query_records
                        ),
                        "ndcg_score_identifiable_query_count": sum(
                            bool(item["ndcg_score_identifiable"])
                            for item in query_records
                        ),
                        "ndcg_tie_dependent_query_count": sum(
                            not bool(item["ndcg_score_identifiable"])
                            for item in query_records
                        ),
                        "legacy_salted_query_any_hit_rate": format_float(
                            record_mean("legacy_salted_any_hit")
                        ),
                        "tie_expected_query_any_hit_probability": format_float(
                            record_mean("tie_expected_any_hit_probability")
                        ),
                        "score_guaranteed_any_hit_query_count": sum(
                            bool(item["score_guaranteed_any_hit"])
                            for item in query_records
                        ),
                        "score_possible_any_hit_query_count": sum(
                            bool(item["score_possible_any_hit"])
                            for item in query_records
                        ),
                        "pmid_source_document_count": document_count,
                        "query_pmid_edge_count": edge_count,
                        "pmid_component_count": len(components),
                        "tie_interpretation": (
                            "not_estimable_empty_prespecified_subset"
                            if not query_records
                            else (
                                "non_operational_uniform_tie_allocation"
                                if subset == "structural_all_zero_non_operational"
                                else "conditional_score_ranking_sensitivity"
                            )
                        ),
                    }
                    require(
                        row["membership_score_identifiable_relation_count"]
                        + row["membership_boundary_tie_dependent_relation_count"]
                        + row["membership_not_retrieved_relation_count"]
                        == row["relevant_relation_count"],
                        "Aggregate membership classifications do not sum to relation count",
                    )
                    require(
                        row["recall_score_identifiable_query_count"]
                        + row["recall_tie_dependent_query_count"]
                        == row["query_count"],
                        "Recall query classifications do not sum to query count",
                    )
                    require(
                        row["ndcg_score_identifiable_query_count"]
                        + row["ndcg_tie_dependent_query_count"]
                        == row["query_count"],
                        "NDCG query classifications do not sum to query count",
                    )
                    key = (scope, baseline, subset, k)
                    require(key not in index, "Aggregate metric cell is duplicated")
                    index[key] = row
                    output_rows.append(row)
    require(len(output_rows) == 48, "Merged-scope tie-aware metric table is not 48 rows")
    return output_rows, index


def build_three_mer_operability_rows(
    relevance: dict[str, dict[str, set[str]]],
    sequence_operational: dict[str, bool],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for scope in DISPLAY_SCOPES:
        queries = sorted(relevance[scope])
        operational = [query for query in queries if sequence_operational[query]]
        structural_zero = [query for query in queries if not sequence_operational[query]]
        rows.append(
            {
                "scope": scope,
                "all_query_count": len(queries),
                "all_relevant_relation_count": sum(
                    len(relevance[scope][query]) for query in queries
                ),
                "score_operational_query_count": len(operational),
                "score_operational_relevant_relation_count": sum(
                    len(relevance[scope][query]) for query in operational
                ),
                "structural_all_zero_query_count": len(structural_zero),
                "structural_all_zero_relevant_relation_count": sum(
                    len(relevance[scope][query]) for query in structural_zero
                ),
                "structural_all_zero_query_fraction": format_float(
                    len(structural_zero) / len(queries)
                ),
                "structural_all_zero_interpretation": (
                    "non_operational; any top-k allocation is a tie-rule realization"
                ),
            }
        )
    require(len(rows) == 4, "Merged-scope 3-mer operability table is not four rows")
    temporal = next(row for row in rows if row["scope"] == "temporal_strict_ab")
    require(
        temporal["score_operational_query_count"] == 60
        and temporal["structural_all_zero_query_count"] == 162,
        "Temporal 3-mer operability split changed",
    )
    return rows


def build_zero_hit_rows(
    metric_index: dict[tuple[str, str, str, int], dict[str, Any]]
) -> list[dict[str, Any]]:
    cold_scopes = [
        "project_defined_joint_scaffold_homology_cold_0_30",
        "project_defined_joint_scaffold_homology_cold_0_50_0_70_identical_mask",
    ]
    rows: list[dict[str, Any]] = []
    for scope in cold_scopes:
        for baseline in BASELINES:
            metric = metric_index[(scope, baseline, "all_queries", 50)]
            n_queries = int(metric["query_count"])
            hit_rate = float(metric["legacy_salted_query_any_hit_rate"])
            hit_count = int(round(hit_rate * n_queries))
            require(
                math.isclose(hit_rate * n_queries, hit_count, abs_tol=1e-12),
                "Legacy hit rate does not correspond to an integer count",
            )
            require(hit_count == 0, "A frozen joint-cold top-50 query hit is nonzero")
            rows.append(
                {
                    "scope": scope,
                    "baseline": baseline,
                    "query_count": n_queries,
                    "legacy_salted_query_hit_count_at_50": hit_count,
                    "legacy_salted_query_hit_rate_at_50": format_float(hit_rate),
                    "empirical_query_bootstrap_ci95_low": "0",
                    "empirical_query_bootstrap_ci95_high": "0",
                    "empirical_query_bootstrap_status": "empirical_degenerate_zero_hits",
                    "one_sided_clopper_pearson_confidence": "0.95",
                    "one_sided_clopper_pearson_upper": format_float(
                        clopper_pearson_upper(hit_count, n_queries, ALPHA)
                    ),
                    "uniform_tie_expected_query_hit_probability_at_50": metric[
                        "tie_expected_query_any_hit_probability"
                    ],
                    "score_guaranteed_hit_query_count_at_50": metric[
                        "score_guaranteed_any_hit_query_count"
                    ],
                    "score_possible_hit_query_count_at_50": metric[
                        "score_possible_any_hit_query_count"
                    ],
                    "pmid_component_count": metric["pmid_component_count"],
                    "bound_assumption": (
                        "exchangeable independent Bernoulli queries; source-document "
                        "dependence is reported separately"
                    ),
                }
            )
    require(len(rows) == 8, "Merged-scope zero-hit table is not eight rows")
    return rows


def uncertainty_estimand_rows() -> list[dict[str, str]]:
    return [
        {
            "method": "uniform_exact_tie_expectation",
            "point_estimand": (
                "equally query-weighted macro mean conditional on frozen scores, "
                "relevance, scopes, and uniform permutations within exact-score blocks"
            ),
            "sampling_or_randomization_unit": "within-score tie permutation",
            "interval_or_bound": "exact expected value plus attainable best-worst bounds",
            "interpretation_boundary": (
                "conditional ranking uncertainty; not biological or dataset sampling uncertainty"
            ),
        },
        {
            "method": "query_bootstrap",
            "point_estimand": (
                "equally query-weighted macro mean of within-tie expected query metrics"
            ),
            "sampling_or_randomization_unit": (
                "observed endpoint query, resampled independently with replacement"
            ),
            "interval_or_bound": "10000-replicate 95% percentile interval",
            "interpretation_boundary": (
                "query-to-query empirical sampling variability under an independence approximation"
            ),
        },
        {
            "method": "pmid_connected_component_bootstrap",
            "point_estimand": (
                "same equally query-weighted macro mean of within-tie expected query metrics"
            ),
            "sampling_or_randomization_unit": (
                "scope/subset-specific non-overlapping query-PMID connected component; "
                "all member queries retained"
            ),
            "interval_or_bound": "10000-replicate 95% percentile sensitivity interval",
            "interpretation_boundary": (
                "source-document dependence sensitivity; not document-disjoint external validation"
            ),
        },
        {
            "method": "one_sided_exact_clopper_pearson",
            "point_estimand": (
                "probability that a query has at least one observed relevant target "
                "in the legacy fixed-salt top 50"
            ),
            "sampling_or_randomization_unit": "binary query-hit event",
            "interval_or_bound": "one-sided exact 95% upper confidence bound",
            "interpretation_boundary": (
                "finite-sample descriptive bound assuming exchangeable independent queries; "
                "does not absorb PMID clustering"
            ),
        },
    ]


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--plan-lock", required=True, type=Path)
    result.add_argument("--parent-protocol", required=True, type=Path)
    result.add_argument("--implementation-lock", required=True, type=Path)
    result.add_argument("--ranks", required=True, type=Path)
    result.add_argument("--endpoint", required=True, type=Path)
    result.add_argument("--scaffold", required=True, type=Path)
    result.add_argument("--homology-0-30", required=True, type=Path)
    result.add_argument("--homology-0-50", required=True, type=Path)
    result.add_argument("--homology-0-70", required=True, type=Path)
    result.add_argument("--source-evidence", required=True, type=Path)
    result.add_argument("--output-dir", required=True, type=Path)
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    started_at = utc_now()
    started = time.perf_counter()
    tracemalloc.start()
    implementation = verify_locks(
        args.plan_lock, args.parent_protocol, args.implementation_lock
    )
    paths = {
        "ranks": args.ranks.resolve(),
        "endpoint": args.endpoint.resolve(),
        "scaffold": args.scaffold.resolve(),
        "homology_0_30": args.homology_0_30.resolve(),
        "homology_0_50": args.homology_0_50.resolve(),
        "homology_0_70": args.homology_0_70.resolve(),
        "source_evidence": args.source_evidence.resolve(),
    }
    inputs = verify_scientific_inputs(paths)
    output_dir = args.output_dir.resolve()
    require(not output_dir.exists(), "Output directory already exists")

    endpoint, query_compounds = load_endpoint(paths["endpoint"])
    scaffold = load_bool_map(
        paths["scaffold"],
        "canonical_pair_key",
        "audit_scaffold_cold_under_selected_policy",
        "audit_outcome",
    )
    homology = {
        "0_30": load_bool_map(
            paths["homology_0_30"],
            "uniprot_canonical_accession",
            "is_future_target_homology_cold_candidate",
            "future_target_coldness_status",
        ),
        "0_50": load_bool_map(
            paths["homology_0_50"],
            "uniprot_canonical_accession",
            "is_future_target_homology_cold_candidate",
            "future_target_coldness_status",
        ),
        "0_70": load_bool_map(
            paths["homology_0_70"],
            "uniprot_canonical_accession",
            "is_future_target_homology_cold_candidate",
            "future_target_coldness_status",
        ),
    }
    internal_scopes = build_internal_scopes(endpoint, scaffold, homology)
    display_rows = {
        display_scope: internal_scopes[internal_scope]
        for display_scope, internal_scope in DISPLAY_SCOPE_TO_INTERNAL.items()
    }
    relevance, query_pairs = make_relevance(display_rows)
    pair_documents, source_audit = load_pair_documents(
        paths["source_evidence"],
        {row["canonical_pair_key"] for row in endpoint},
    )
    metric_records, sequence_operational, candidate_counts = load_rank_metrics(
        paths["ranks"], query_compounds, relevance
    )
    metric_rows, metric_index = aggregate_metrics(
        metric_records,
        relevance,
        query_pairs,
        pair_documents,
        sequence_operational,
    )
    operability_rows = build_three_mer_operability_rows(relevance, sequence_operational)
    zero_hit_rows = build_zero_hit_rows(metric_index)
    estimand_rows = uncertainty_estimand_rows()

    output_dir.mkdir(parents=True, exist_ok=False)
    metric_path = output_dir / "tie_aware_metrics.tsv"
    operability_path = output_dir / "three_mer_operability.tsv"
    zero_hit_path = output_dir / "double_cold_query_hit_upper_bounds.tsv"
    estimand_path = output_dir / "uncertainty_estimands.tsv"
    summary_path = output_dir / "tie_aware_summary.json"
    receipt_path = output_dir / "execution_receipt.json"
    manifest_path = output_dir / "run_manifest.json"

    write_tsv_new(metric_path, list(metric_rows[0]), metric_rows)
    write_tsv_new(operability_path, list(operability_rows[0]), operability_rows)
    write_tsv_new(zero_hit_path, list(zero_hit_rows[0]), zero_hit_rows)
    write_tsv_new(estimand_path, list(estimand_rows[0]), estimand_rows)

    temporal_operability = next(
        row for row in operability_rows if row["scope"] == "temporal_strict_ab"
    )
    headline_rows = [
        {
            "scope": row["scope"],
            "baseline": row["baseline"],
            "legacy_salted_recall_at_50": row["legacy_salted_recall"],
            "tie_expected_fractional_recall_at_50": row[
                "tie_expected_fractional_recall"
            ],
            "tie_worst_recall_at_50": row["tie_worst_recall"],
            "tie_best_recall_at_50": row["tie_best_recall"],
            "legacy_salted_ndcg_at_50": row["legacy_salted_ndcg"],
            "tie_expected_ndcg_at_50": row["tie_expected_ndcg"],
            "tie_worst_ndcg_at_50": row["tie_worst_ndcg"],
            "tie_best_ndcg_at_50": row["tie_best_ndcg"],
        }
        for row in metric_rows
        if row["query_subset"] == "all_queries" and row["k"] == 50
    ]
    summary = {
        "schema_version": "1.0",
        "analysis_id": ANALYSIS_ID,
        "protocol_id": PROTOCOL_ID,
        "analysis_role": "reviewer_requested_post_hoc_sensitivity",
        "claim_boundary": (
            "Author-run, outcome-visible, non-independent retrospective methodological "
            "audit; no external-validation or biological-discovery claim."
        ),
        "endpoint_contract_verified": {
            "relation_count": len(endpoint),
            "query_count": len(query_compounds),
            "target_count": len(
                {row["uniprot_canonical_accession"] for row in endpoint}
            ),
        },
        "rank_contract_verified": {
            "complete_rank_row_count": EXPECTED_RANK_ROWS,
            "baseline_query_group_count": len(query_compounds) * len(BASELINES),
            "candidate_target_count_min": min(candidate_counts.values()),
            "candidate_target_count_max": max(candidate_counts.values()),
        },
        "display_scope_policy": {
            "display_scope_count": len(DISPLAY_SCOPES),
            "homology_0_50_and_0_70_identical": True,
            "identical_masks_displayed_once": True,
            "separate_input_receipts_retained": True,
        },
        "tie_analysis": {
            "within_block_distribution": "uniform_over_all_exact_score_permutations",
            "metric_rows": len(metric_rows),
            "membership_classes": [
                "score_identifiable",
                "boundary_tie_dependent",
                "not_retrieved",
            ],
            "multi_salt_repetition_run": False,
        },
        "three_mer_temporal_operability": temporal_operability,
        "bootstrap": {
            "replicates": BOOTSTRAP_REPLICATES,
            "prng": "PCG64",
            "query_bootstrap_base_seed": QUERY_BOOTSTRAP_BASE_SEED,
            "pmid_component_bootstrap_base_seed": COMPONENT_BOOTSTRAP_BASE_SEED,
            "interval": "95_percentile_numpy_linear",
            "estimands_are_distinctly_labelled": True,
        },
        "source_document_audit": source_audit,
        "zero_hit_bound": {
            "rows": len(zero_hit_rows),
            "confidence": 0.95,
            "method": "one_sided_exact_Clopper_Pearson",
            "empirical_zero_width_intervals_labelled_degenerate": True,
        },
        "headline_all_query_metrics_at_50": headline_rows,
        "output_boundary": {
            "aggregate_only": True,
            "identifier_bearing_rows": False,
            "component_membership_written": False,
            "alternative_rank_ledger_written": False,
            "absolute_paths_written": False,
        },
        "interpretation_boundary": (
            "Tie-aware expectations and bounds are conditional ranking sensitivities. "
            "The PMID-component bootstrap diagnoses source-document dependence and "
            "does not create document-disjoint validation."
        ),
        "created_at_utc": utc_now(),
    }
    write_json_new(summary_path, summary)

    current_memory, peak_memory = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    receipt = {
        "schema_version": "1.0",
        "analysis_id": ANALYSIS_ID,
        "status": "completed",
        "started_at_utc": started_at,
        "ended_at_utc": utc_now(),
        "wall_seconds": time.perf_counter() - started,
        "environment": {
            "python": sys.version,
            "platform": platform.platform(),
            "numpy": np.__version__,
            "python_tracemalloc_current_bytes": current_memory,
            "python_tracemalloc_peak_bytes": peak_memory,
        },
        "command_receipt": {
            "raw_argv_recorded": False,
            "argument_switches": [
                "--plan-lock",
                "--parent-protocol",
                "--implementation-lock",
                "--ranks",
                "--endpoint",
                "--scaffold",
                "--homology-0-30",
                "--homology-0-50",
                "--homology-0-70",
                "--source-evidence",
                "--output-dir",
            ],
        },
        "locks": {
            "parent_protocol": {
                "basename": args.parent_protocol.name,
                "sha256": sha256(args.parent_protocol),
            },
            "local_plan": {
                "basename": args.plan_lock.name,
                "sha256": sha256(args.plan_lock),
            },
            "implementation": {
                "basename": args.implementation_lock.name,
                "sha256": sha256(args.implementation_lock),
            },
        },
        "preexecution_synthetic_test_receipt": implementation.get(
            "synthetic_test_receipt"
        ),
        "deviations_from_local_plan": [
            {
                "item": "display_scope_row_contract",
                "local_plan": (
                    "five separate scopes; 60 tie rows, five operability rows, "
                    "and 12 zero-hit rows"
                ),
                "implemented": (
                    "the byte-identical 0.50 and 0.70 masks are displayed once; "
                    "48 tie rows, four operability rows, and eight zero-hit rows"
                ),
                "reason": (
                    "Required by the later parent frozen major-revision protocol v4; "
                    "both mask inputs remain separately hash-receipted."
                ),
                "scientific_calculation_changed": False,
            }
        ],
        "external_transfer_performed": False,
        "identifier_bearing_output_written": False,
        "absolute_paths_written": False,
        "outputs_before_receipt": [
            output_descriptor(path)
            for path in (
                metric_path,
                operability_path,
                zero_hit_path,
                estimand_path,
                summary_path,
            )
        ],
    }
    write_json_new(receipt_path, receipt)

    manifest = {
        "schema_version": "1.0",
        "analysis_id": ANALYSIS_ID,
        "protocol_id": PROTOCOL_ID,
        "aggregate_only": True,
        "identifier_bearing_output": False,
        "absolute_paths_recorded": False,
        "inputs": inputs,
        "locks": receipt["locks"],
        "implementation_files": [
            {
                "relative_name": relative_name,
                "sha256": expected_hash,
            }
            for relative_name, expected_hash in sorted(
                implementation["implementation_files"].items()
            )
        ],
        "outputs": [
            output_descriptor(path)
            for path in (
                metric_path,
                operability_path,
                zero_hit_path,
                estimand_path,
                summary_path,
                receipt_path,
            )
        ],
        "output_contract": {
            "tie_aware_metrics_rows": len(metric_rows),
            "three_mer_operability_rows": len(operability_rows),
            "double_cold_query_hit_upper_bound_rows": len(zero_hit_rows),
            "uncertainty_estimand_rows": len(estimand_rows),
            "homology_0_50_0_70_displayed_once": True,
        },
        "created_at_utc": utc_now(),
    }
    write_json_new(manifest_path, manifest)
    print(
        json.dumps(
            {
                "analysis_id": ANALYSIS_ID,
                "status": "completed",
                "metric_rows": len(metric_rows),
                "operability_rows": len(operability_rows),
                "zero_hit_rows": len(zero_hit_rows),
                "identifier_bearing_output": False,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
