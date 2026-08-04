"""Run the frozen aggregate-only Tier B weight-policy sensitivity."""

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
from typing import Any

import numpy as np
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator
from sklearn import __version__ as sklearn_version
from sklearn.feature_extraction.text import TfidfVectorizer


ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
sys.path.insert(0, str(WORKSPACE / "scripts"))
from pu_retrieval_metrics import macro_average, query_metrics, rank_scores


ANALYSIS_ID = "revision_weight_policy_v1_20260729"
PROTOCOL_ID = "npass_strict_ab_major_revision_v4_20260729_analysis_D_weight_policy"
PARENT_PROTOCOL_SHA256 = "bd4902476160cc7c5cbacaf0cfd0f1a28c5300bde22232b3a2cc6c1f3c143dc2"
RANK_IMPLEMENTATION_SHA256 = "7e40e80c3203a4a6cf95ca675cc9c333168f8feb234449856f26b5e132ec3165"
EXPECTED_PRIMARY_RANK_ROWS = 3_658_128
EXPECTED_COMPUTED_SCORE_RANK_ROWS = 10_974_384
EXPECTED_HISTORY_RELATIONS = 4_990
EXPECTED_QUERIES = 222
EXPECTED_TARGETS = 4_123
EXPECTED_ENDPOINT_RELATIONS = 358
EXPECTED_ENDPOINT_TARGETS = 156
TIER_A = "A_affinity_candidate"
TIER_B = "B_quantitative_functional_candidate"
SALT = "npass_strict_ab_doublecold_successor_v1_20260719"
PAIR_NEIGHBOR_TOP_K = 100
MORGAN_RADIUS = 2
MORGAN_BITS = 2048

VARIANTS = [
    {
        "key": "A1_B0_5",
        "tier_A_weight": 1.0,
        "tier_B_weight": 0.5,
        "role": "lower_B_policy_sensitivity",
        "alias": "",
    },
    {
        "key": "A1_B0_7_primary",
        "tier_A_weight": 1.0,
        "tier_B_weight": 0.7,
        "role": "frozen_primary_reference",
        "alias": "",
    },
    {
        "key": "A1_B1_0_all_equal",
        "tier_A_weight": 1.0,
        "tier_B_weight": 1.0,
        "role": "all_equal_policy_sensitivity",
        "alias": "all_equal_A1_B1",
    },
]
PRIMARY_VARIANT_INDEX = 1
BASELINES = [
    "weighted_target_popularity",
    "sequence_3mer_transfer",
    "weighted_morgan_transfer",
    "structure_sequence_pair_neighbor",
]
METRICS = ["Recall@10", "Recall@50", "NDCG@10", "NDCG@50", "MRR"]
DISPLAY_SCOPE_TO_INTERNAL = {
    "temporal_strict_ab": "temporal_strict_ab",
    "scaffold_cold_strict_ab": "scaffold_cold_strict_ab",
    "project_defined_joint_scaffold_homology_cold_0_30": "double_cold_0_30",
    "project_defined_joint_scaffold_homology_cold_0_50_0_70_identical_mask": "double_cold_0_50",
}
DISPLAY_SCOPES = list(DISPLAY_SCOPE_TO_INTERNAL)

EXPECTED_INPUTS = {
    "historical_pairs": (
        "historical_pairs.tsv.gz",
        "cef748ae8ac277e49784d7e1fbf08e085beb14a84fc6c651a0fb8d99e88710d7",
    ),
    "scoring_queries": (
        "scoring_queries.tsv.gz",
        "0e6068d2e25cb3ea325656fb3517563788cd496e88cfaa3de761890fec9e9318",
    ),
    "historical_compounds": (
        "historical_compounds.tsv.gz",
        "f1f82793b5c652007a042699c19cd5640a8b68e8b1f0d4f94e4dc4f54045060c",
    ),
    "query_compounds": (
        "query_compounds.tsv.gz",
        "f51670fffd21d2e9109b4376dd53aab55bbacb09d2e4795dfa515fcaca98b113",
    ),
    "candidate_targets": (
        "candidate_targets.tsv.gz",
        "0ee86746b306fb388a1f74a6b88ce4d1eba01b7a4eb473315f6b3def57145cdc",
    ),
    "candidate_sequences": (
        "candidate_sequences.fasta",
        "a83421dba2482f236fe18340dd592cc7d5ed22c98c4fc39435c40f04f289b442",
    ),
    "primary_complete_ranks": (
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


def format_float(value: float | None) -> str:
    return "" if value is None else f"{float(value):.17g}"


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


def read_tsv_gz(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or []), list(reader)


def parse_bool(value: str) -> bool:
    normalized = value.strip().lower()
    require(normalized in {"true", "false"}, f"Invalid boolean value: {value!r}")
    return normalized == "true"


def load_fasta(path: Path) -> dict[str, str]:
    output: dict[str, str] = {}
    current: str | None = None
    parts: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    require(current not in output, "FASTA identifier is duplicated")
                    output[current] = "".join(parts)
                current = line[1:].split()[0]
                parts = []
            else:
                require(current is not None, "FASTA sequence appears before a header")
                parts.append(line)
    if current is not None:
        require(current not in output, "FASTA identifier is duplicated")
        output[current] = "".join(parts)
    return output


def verify_locks(
    plan_lock: Path, parent_protocol: Path, implementation_lock_path: Path
) -> dict[str, Any]:
    require(plan_lock.is_file(), "Analysis-plan lock is absent")
    require(parent_protocol.is_file(), "Parent protocol is absent")
    require(implementation_lock_path.is_file(), "Implementation lock is absent")
    require(sha256(parent_protocol) == PARENT_PROTOCOL_SHA256, "Parent protocol hash changed")
    plan = json.loads(plan_lock.read_text(encoding="utf-8"))
    require(plan.get("analysis_id") == ANALYSIS_ID, "Plan analysis ID changed")
    require(
        plan.get("lock_state")
        == "LOCKED_BEFORE_CURRENT_CORRECTED_WEIGHT_SENSITIVITY_EXECUTION",
        "Plan was not frozen before current execution",
    )
    implementation = json.loads(implementation_lock_path.read_text(encoding="utf-8"))
    require(
        implementation.get("lock_state") == "LOCKED_BEFORE_REAL_INPUT_EXECUTION",
        "Implementation was not locked before real-input execution",
    )
    require(
        implementation.get("plan_sha256") == sha256(plan_lock),
        "Implementation lock plan anchor changed",
    )
    require(
        implementation.get("parent_protocol_sha256") == PARENT_PROTOCOL_SHA256,
        "Implementation lock parent anchor changed",
    )
    files = implementation.get("implementation_files")
    require(isinstance(files, dict) and files, "Implementation file inventory is absent")
    for relative_name, expected_hash in sorted(files.items()):
        candidate = (ROOT / relative_name).resolve()
        try:
            candidate.relative_to(ROOT.resolve())
        except ValueError as exc:
            raise ValueError("Implementation path escapes analysis root") from exc
        require(candidate.is_file(), f"Implementation file is absent: {relative_name}")
        require(sha256(candidate) == expected_hash, f"Implementation drifted: {relative_name}")
    require(
        files.get("scripts/run_weight_policy_sensitivity_v1.py")
        == sha256(Path(__file__).resolve()),
        "Executing script differs from implementation lock",
    )
    ranking_path = WORKSPACE / "scripts" / "pu_retrieval_metrics.py"
    require(
        ranking_path.is_file() and sha256(ranking_path) == RANK_IMPLEMENTATION_SHA256,
        "Frozen ranking implementation changed",
    )
    return implementation


def verify_inputs(paths: dict[str, Path]) -> list[dict[str, Any]]:
    descriptors: list[dict[str, Any]] = []
    for role, path in paths.items():
        require(path.is_file(), f"Input is absent: {role}")
        expected_name, expected_hash = EXPECTED_INPUTS[role]
        require(path.name == expected_name, f"Unexpected input basename: {role}")
        require(sha256(path) == expected_hash, f"Frozen input hash changed: {role}")
        descriptors.append(input_descriptor(role, path))
    return descriptors


def load_bool_map(path: Path, key_field: str, flag_field: str, status_field: str) -> dict[str, bool]:
    fields, rows = read_tsv_gz(path)
    require({key_field, flag_field, status_field}.issubset(fields), "Scope fields missing")
    output: dict[str, bool] = {}
    for row in rows:
        key = row[key_field].strip()
        require(key and key not in output, "Scope key is empty or duplicated")
        require(row[status_field].strip(), "Scope status is empty")
        output[key] = parse_bool(row[flag_field])
    return output


def load_and_validate_inputs(paths: dict[str, Path]) -> dict[str, Any]:
    history_fields, history = read_tsv_gz(paths["historical_pairs"])
    require(
        {
            "canonical_pair_key",
            "inchikey_full",
            "uniprot_canonical_accession",
            "best_strict_evidence_tier",
        }.issubset(history_fields),
        "Historical-pair fields are incomplete",
    )
    require(len(history) == EXPECTED_HISTORY_RELATIONS, "Historical relation count changed")
    require(
        len({row["canonical_pair_key"] for row in history}) == len(history),
        "Historical pairs are duplicated",
    )
    require(
        all(row["best_strict_evidence_tier"] in {TIER_A, TIER_B} for row in history),
        "Historical pairs contain a non-strict tier",
    )

    query_fields, queries = read_tsv_gz(paths["scoring_queries"])
    require({"query_id", "inchikey_full"}.issubset(query_fields), "Query fields missing")
    require(len(queries) == EXPECTED_QUERIES, "Scoring-query count changed")
    require(
        len({row["query_id"] for row in queries}) == len(queries)
        and len({row["inchikey_full"] for row in queries}) == len(queries),
        "Scoring queries are not one-to-one",
    )

    historical_structure_fields, historical_structure_rows = read_tsv_gz(
        paths["historical_compounds"]
    )
    query_structure_fields, query_structure_rows = read_tsv_gz(paths["query_compounds"])
    for fields, rows, role in (
        (historical_structure_fields, historical_structure_rows, "historical"),
        (query_structure_fields, query_structure_rows, "query"),
    ):
        require(
            {"inchikey_full", "representative_smiles", "structure_role"}.issubset(fields),
            f"{role} structure fields missing",
        )
        require(
            len({row["inchikey_full"] for row in rows}) == len(rows),
            f"{role} structures are duplicated",
        )
        require(
            all(row["structure_role"] == role for row in rows),
            f"{role} structure role changed",
        )
    historical_smiles = {
        row["inchikey_full"]: row["representative_smiles"]
        for row in historical_structure_rows
    }
    query_smiles = {
        row["inchikey_full"]: row["representative_smiles"] for row in query_structure_rows
    }
    require(
        set(historical_smiles) == {row["inchikey_full"] for row in history},
        "Historical structure keyset changed",
    )
    require(
        set(query_smiles) == {row["inchikey_full"] for row in queries},
        "Query structure keyset changed",
    )

    target_fields, target_rows = read_tsv_gz(paths["candidate_targets"])
    require(
        "uniprot_canonical_accession" in target_fields,
        "Candidate-target field missing",
    )
    target_ids = [row["uniprot_canonical_accession"] for row in target_rows]
    require(
        len(target_ids) == EXPECTED_TARGETS and len(set(target_ids)) == EXPECTED_TARGETS,
        "Candidate-target universe changed",
    )
    sequences = load_fasta(paths["candidate_sequences"])
    require(set(sequences) == set(target_ids), "Candidate FASTA keyset changed")

    endpoint_fields, endpoint = read_tsv_gz(paths["endpoint"])
    require(
        {
            "canonical_pair_key",
            "query_id",
            "inchikey_full",
            "uniprot_canonical_accession",
        }.issubset(endpoint_fields),
        "Endpoint fields missing",
    )
    require(len(endpoint) == EXPECTED_ENDPOINT_RELATIONS, "Endpoint relation count changed")
    require(
        len({row["query_id"] for row in endpoint}) == EXPECTED_QUERIES,
        "Endpoint query count changed",
    )
    require(
        len({row["uniprot_canonical_accession"] for row in endpoint})
        == EXPECTED_ENDPOINT_TARGETS,
        "Endpoint target count changed",
    )
    query_compounds = {row["query_id"]: row["inchikey_full"] for row in queries}
    for row in endpoint:
        require(
            query_compounds[row["query_id"]] == row["inchikey_full"],
            "Endpoint query/compound mapping changed",
        )

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
    endpoint_pairs = {row["canonical_pair_key"] for row in endpoint}
    endpoint_targets = {row["uniprot_canonical_accession"] for row in endpoint}
    require(set(scaffold) == endpoint_pairs, "Scaffold keyset differs from endpoint")
    for threshold, flags in homology.items():
        require(set(flags) == endpoint_targets, f"Homology {threshold} keyset changed")

    internal_scopes: dict[str, list[dict[str, str]]] = {
        "temporal_strict_ab": [],
        "scaffold_cold_strict_ab": [],
        "double_cold_0_30": [],
        "double_cold_0_50": [],
        "double_cold_0_70": [],
    }
    for row in endpoint:
        internal_scopes["temporal_strict_ab"].append(row)
        if scaffold[row["canonical_pair_key"]]:
            internal_scopes["scaffold_cold_strict_ab"].append(row)
            target = row["uniprot_canonical_accession"]
            for threshold in ("0_30", "0_50", "0_70"):
                if homology[threshold][target]:
                    internal_scopes[f"double_cold_{threshold}"].append(row)
    require(
        {row["canonical_pair_key"] for row in internal_scopes["double_cold_0_50"]}
        == {row["canonical_pair_key"] for row in internal_scopes["double_cold_0_70"]},
        "0.50 and 0.70 masks are not identical",
    )
    display_scope_rows = {
        display: internal_scopes[internal]
        for display, internal in DISPLAY_SCOPE_TO_INTERNAL.items()
    }
    relevance: dict[str, dict[str, list[str]]] = {}
    for scope, rows in display_scope_rows.items():
        by_query: dict[str, list[str]] = defaultdict(list)
        for row in rows:
            by_query[row["query_id"]].append(row["uniprot_canonical_accession"])
        relevance[scope] = dict(by_query)

    return {
        "history": history,
        "queries": queries,
        "historical_smiles": historical_smiles,
        "query_smiles": query_smiles,
        "target_ids": target_ids,
        "sequences": sequences,
        "endpoint": endpoint,
        "display_scope_rows": display_scope_rows,
        "relevance": relevance,
    }


def build_fingerprints(smiles_by_compound: dict[str, str]) -> dict[str, Any]:
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=MORGAN_RADIUS, fpSize=MORGAN_BITS
    )
    output: dict[str, Any] = {}
    for compound, smiles in smiles_by_compound.items():
        molecule = Chem.MolFromSmiles(smiles)
        require(molecule is not None, "RDKit could not parse a frozen structure")
        output[compound] = generator.GetFingerprint(molecule)
    return output


def build_sequence_matrix(target_ids: list[str], sequences: dict[str, str]) -> Any:
    vectorizer = TfidfVectorizer(
        analyzer="char",
        ngram_range=(3, 3),
        lowercase=False,
        norm="l2",
        dtype=np.float32,
    )
    return vectorizer.fit_transform([sequences[target] for target in target_ids])


def variant_weight(variant_index: int, tier: str) -> float:
    variant = VARIANTS[variant_index]
    if tier == TIER_A:
        return float(variant["tier_A_weight"])
    require(tier == TIER_B, "Unknown evidence tier")
    return float(variant["tier_B_weight"])


def build_history_maps(
    history: list[dict[str, str]], target_index: dict[str, int]
) -> tuple[
    dict[str, list[tuple[int, str]]],
    list[str],
    np.ndarray,
    np.ndarray,
    dict[int, int],
]:
    by_compound: dict[str, list[tuple[int, str]]] = defaultdict(list)
    popularity = np.zeros((len(VARIANTS), len(target_index)), dtype=np.float32)
    for row in history:
        target_idx = target_index[row["uniprot_canonical_accession"]]
        tier = row["best_strict_evidence_tier"]
        by_compound[row["inchikey_full"]].append((target_idx, tier))
        for variant_index in range(len(VARIANTS)):
            popularity[variant_index, target_idx] += variant_weight(variant_index, tier)
    historical_compounds = list(by_compound)
    historical_target_indices = np.asarray(
        sorted(
            {
                target_idx
                for relations in by_compound.values()
                for target_idx, _ in relations
            }
        ),
        dtype=np.int32,
    )
    target_to_historical_column = {
        int(target_idx): column
        for column, target_idx in enumerate(historical_target_indices)
    }
    return (
        dict(by_compound),
        historical_compounds,
        popularity,
        historical_target_indices,
        target_to_historical_column,
    )


def compute_weighted_activation(
    chemical_similarities: list[float],
    historical_compounds: list[str],
    train_by_compound: dict[str, list[tuple[int, str]]],
    target_to_historical_column: dict[int, int],
) -> np.ndarray:
    output = np.zeros(
        (len(VARIANTS), len(target_to_historical_column)), dtype=np.float32
    )
    for compound, similarity in zip(historical_compounds, chemical_similarities):
        for target_idx, tier in train_by_compound[compound]:
            column = target_to_historical_column[target_idx]
            for variant_index in range(len(VARIANTS)):
                value = float(similarity) * variant_weight(variant_index, tier)
                if value > float(output[variant_index, column]):
                    output[variant_index, column] = value
    return output


def compute_sequence_scores(
    query_compound: str,
    train_by_compound: dict[str, list[tuple[int, str]]],
    sequence_matrix: Any,
    target_count: int,
) -> np.ndarray:
    output = np.zeros((len(VARIANTS), target_count), dtype=np.float32)
    known = train_by_compound.get(query_compound, [])
    if not known:
        return output
    indices = [target_idx for target_idx, _ in known]
    similarities = (sequence_matrix @ sequence_matrix[indices].T).toarray().astype(
        np.float32, copy=False
    )
    for variant_index in range(len(VARIANTS)):
        weights = np.asarray(
            [variant_weight(variant_index, tier) for _, tier in known],
            dtype=np.float32,
        )
        output[variant_index] = np.max(
            similarities * weights[np.newaxis, :], axis=1
        ).astype(np.float32, copy=False)
    return output


def precompute_pair_topk(
    sequence_matrix: Any, historical_target_indices: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    similarities = (
        sequence_matrix @ sequence_matrix[historical_target_indices].T
    ).toarray().astype(np.float32, copy=False)
    k = min(PAIR_NEIGHBOR_TOP_K, int(historical_target_indices.size))
    columns = np.argpartition(similarities, -k, axis=1)[:, -k:]
    values = np.take_along_axis(similarities, columns, axis=1)
    return columns.astype(np.int32, copy=False), values.astype(np.float32, copy=False)


def complete_scores_for_query(
    query_compound: str,
    query_fingerprint: Any,
    historical_fingerprints: dict[str, Any],
    historical_compounds: list[str],
    train_by_compound: dict[str, list[tuple[int, str]]],
    target_to_historical_column: dict[int, int],
    historical_target_indices: np.ndarray,
    popularity: np.ndarray,
    sequence_matrix: Any,
    pair_top_columns: np.ndarray,
    pair_top_similarities: np.ndarray,
    target_count: int,
) -> dict[str, np.ndarray]:
    chemical_similarities = DataStructs.BulkTanimotoSimilarity(
        query_fingerprint,
        [historical_fingerprints[compound] for compound in historical_compounds],
    )
    activation = compute_weighted_activation(
        chemical_similarities,
        historical_compounds,
        train_by_compound,
        target_to_historical_column,
    )
    morgan = np.zeros((len(VARIANTS), target_count), dtype=np.float32)
    morgan[:, historical_target_indices] = activation
    pair = np.empty((len(VARIANTS), target_count), dtype=np.float32)
    for variant_index in range(len(VARIANTS)):
        pair[variant_index] = np.max(
            pair_top_similarities
            * activation[variant_index, pair_top_columns],
            axis=1,
        ).astype(np.float32, copy=False)
    sequence = compute_sequence_scores(
        query_compound, train_by_compound, sequence_matrix, target_count
    )
    return {
        "weighted_target_popularity": popularity,
        "sequence_3mer_transfer": sequence,
        "weighted_morgan_transfer": morgan,
        "structure_sequence_pair_neighbor": pair,
    }


def initialize_rank_change_stats() -> dict[tuple[str, str], dict[str, int]]:
    output: dict[tuple[str, str], dict[str, int]] = {}
    for variant in VARIANTS:
        for baseline in BASELINES:
            output[(str(variant["key"]), baseline)] = {
                "eligible_candidate_rows": 0,
                "score_changed_candidate_count": 0,
                "rank_changed_candidate_count": 0,
                "absolute_rank_change_sum": 0,
                "maximum_absolute_rank_change": 0,
                "query_count_with_any_rank_change": 0,
                "top50_symmetric_difference_membership_count": 0,
                "query_count_with_any_top50_membership_change": 0,
                "endpoint_relation_count": 0,
                "endpoint_relation_rank_changed_count": 0,
                "endpoint_relation_top50_membership_changed_count": 0,
                "rank_permutation_blocks_checked": 0,
            }
    return output


def update_rank_change_stats(
    stats: dict[str, int],
    candidate_scores: np.ndarray,
    candidate_ranks: np.ndarray,
    primary_scores: np.ndarray,
    primary_ranks: np.ndarray,
    allowed: np.ndarray,
    endpoint_indices: np.ndarray,
) -> None:
    score_changed = candidate_scores[allowed] != primary_scores[allowed]
    rank_deltas = (
        candidate_ranks[allowed].astype(np.int64)
        - primary_ranks[allowed].astype(np.int64)
    )
    rank_changed = rank_deltas != 0
    absolute_deltas = np.abs(rank_deltas)
    top50_changed = (candidate_ranks[allowed] <= 50) != (primary_ranks[allowed] <= 50)
    endpoint_rank_changed = (
        candidate_ranks[endpoint_indices] != primary_ranks[endpoint_indices]
    )
    endpoint_top50_changed = (
        (candidate_ranks[endpoint_indices] <= 50)
        != (primary_ranks[endpoint_indices] <= 50)
    )
    stats["eligible_candidate_rows"] += int(allowed.sum())
    stats["score_changed_candidate_count"] += int(np.count_nonzero(score_changed))
    stats["rank_changed_candidate_count"] += int(np.count_nonzero(rank_changed))
    stats["absolute_rank_change_sum"] += int(absolute_deltas.sum())
    if absolute_deltas.size:
        stats["maximum_absolute_rank_change"] = max(
            stats["maximum_absolute_rank_change"], int(absolute_deltas.max())
        )
    stats["query_count_with_any_rank_change"] += int(bool(np.any(rank_changed)))
    stats["top50_symmetric_difference_membership_count"] += int(
        np.count_nonzero(top50_changed)
    )
    stats["query_count_with_any_top50_membership_change"] += int(
        bool(np.any(top50_changed))
    )
    stats["endpoint_relation_count"] += int(endpoint_indices.size)
    stats["endpoint_relation_rank_changed_count"] += int(
        np.count_nonzero(endpoint_rank_changed)
    )
    stats["endpoint_relation_top50_membership_changed_count"] += int(
        np.count_nonzero(endpoint_top50_changed)
    )
    stats["rank_permutation_blocks_checked"] += 1


def validate_primary_block(
    reader: csv.DictReader,
    query_id: str,
    query_compound: str,
    baseline: str,
    target_ids: list[str],
    allowed: np.ndarray,
    scores: np.ndarray,
    ranks: np.ndarray,
) -> int:
    verified = 0
    for target_index, target in enumerate(target_ids):
        if not allowed[target_index]:
            continue
        try:
            row = next(reader)
        except StopIteration as exc:
            raise ValueError("Frozen primary rank ledger ended early") from exc
        require(row["baseline"] == baseline, "Frozen primary baseline order changed")
        require(row["query_id"] == query_id, "Frozen primary query order changed")
        require(
            row["query_compound_inchikey_full"] == query_compound,
            "Frozen primary query/compound mapping changed",
        )
        require(
            row["target_uniprot_accession"] == target,
            "Frozen primary target order changed",
        )
        require(
            int(row["eligible_candidate_target_count"]) == int(allowed.sum()),
            "Frozen primary candidate count changed",
        )
        require(int(row["rank"]) == int(ranks[target_index]), "Primary rank not reproduced")
        require(
            float(row["score"]) == float(scores[target_index]),
            "Primary score not reproduced exactly",
        )
        verified += 1
    return verified


def scope_cardinality_rows(
    display_scope_rows: dict[str, list[dict[str, str]]]
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    primary_counts = {
        scope: (
            len(scope_rows),
            len({row["query_id"] for row in scope_rows}),
            len({row["uniprot_canonical_accession"] for row in scope_rows}),
        )
        for scope, scope_rows in display_scope_rows.items()
    }
    for variant in VARIANTS:
        for scope in DISPLAY_SCOPES:
            relation_count, query_count, target_count = primary_counts[scope]
            rows.append(
                {
                    "weight_variant": variant["key"],
                    "scope": scope,
                    "relation_count": relation_count,
                    "query_count": query_count,
                    "target_count": target_count,
                    "relation_count_change_vs_0_7": 0,
                    "query_count_change_vs_0_7": 0,
                    "target_count_change_vs_0_7": 0,
                    "status": "invariant_by_frozen_endpoint_and_scope_masks",
                }
            )
    return rows


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--plan-lock", required=True, type=Path)
    result.add_argument("--parent-protocol", required=True, type=Path)
    result.add_argument("--implementation-lock", required=True, type=Path)
    result.add_argument("--historical-pairs", required=True, type=Path)
    result.add_argument("--scoring-queries", required=True, type=Path)
    result.add_argument("--historical-compounds", required=True, type=Path)
    result.add_argument("--query-compounds", required=True, type=Path)
    result.add_argument("--candidate-targets", required=True, type=Path)
    result.add_argument("--candidate-sequences", required=True, type=Path)
    result.add_argument("--primary-complete-ranks", required=True, type=Path)
    result.add_argument("--endpoint", required=True, type=Path)
    result.add_argument("--scaffold", required=True, type=Path)
    result.add_argument("--homology-0-30", required=True, type=Path)
    result.add_argument("--homology-0-50", required=True, type=Path)
    result.add_argument("--homology-0-70", required=True, type=Path)
    result.add_argument("--output-dir", required=True, type=Path)
    return result


def main(argv: list[str] | None = None) -> int:
    args = parser().parse_args(argv)
    started_at = utc_now()
    started = time.perf_counter()
    tracemalloc.start()
    implementation = verify_locks(
        args.plan_lock.resolve(),
        args.parent_protocol.resolve(),
        args.implementation_lock.resolve(),
    )
    paths = {
        "historical_pairs": args.historical_pairs.resolve(),
        "scoring_queries": args.scoring_queries.resolve(),
        "historical_compounds": args.historical_compounds.resolve(),
        "query_compounds": args.query_compounds.resolve(),
        "candidate_targets": args.candidate_targets.resolve(),
        "candidate_sequences": args.candidate_sequences.resolve(),
        "primary_complete_ranks": args.primary_complete_ranks.resolve(),
        "endpoint": args.endpoint.resolve(),
        "scaffold": args.scaffold.resolve(),
        "homology_0_30": args.homology_0_30.resolve(),
        "homology_0_50": args.homology_0_50.resolve(),
        "homology_0_70": args.homology_0_70.resolve(),
    }
    inputs = verify_inputs(paths)
    output_dir = args.output_dir.resolve()
    require(not output_dir.exists(), "Output directory already exists")
    loaded = load_and_validate_inputs(paths)

    target_ids: list[str] = loaded["target_ids"]
    target_index = {target: index for index, target in enumerate(target_ids)}
    (
        train_by_compound,
        historical_compounds,
        popularity,
        historical_target_indices,
        target_to_historical_column,
    ) = build_history_maps(loaded["history"], target_index)
    historical_fingerprints = build_fingerprints(loaded["historical_smiles"])
    query_fingerprints = build_fingerprints(loaded["query_smiles"])
    sequence_matrix = build_sequence_matrix(target_ids, loaded["sequences"])
    pair_top_columns, pair_top_similarities = precompute_pair_topk(
        sequence_matrix, historical_target_indices
    )

    rank_stats = initialize_rank_change_stats()
    per_query_metrics: dict[
        tuple[str, str, str], list[dict[str, float]]
    ] = defaultdict(list)
    primary_verified_rows = 0
    computed_score_rank_rows = 0

    with gzip.open(
        paths["primary_complete_ranks"], "rt", encoding="utf-8", newline=""
    ) as primary_handle:
        primary_reader = csv.DictReader(primary_handle, delimiter="\t")
        required_rank_fields = {
            "baseline",
            "query_id",
            "query_compound_inchikey_full",
            "target_uniprot_accession",
            "rank",
            "score",
            "eligible_candidate_target_count",
        }
        require(
            required_rank_fields.issubset(set(primary_reader.fieldnames or [])),
            "Frozen primary rank fields are incomplete",
        )
        for query_row in loaded["queries"]:
            query_id = query_row["query_id"]
            query_compound = query_row["inchikey_full"]
            allowed = np.ones(len(target_ids), dtype=bool)
            for target_idx, _ in train_by_compound.get(query_compound, []):
                allowed[target_idx] = False
            require(bool(np.any(allowed)), "A query has no eligible candidate target")
            query_scores = complete_scores_for_query(
                query_compound=query_compound,
                query_fingerprint=query_fingerprints[query_compound],
                historical_fingerprints=historical_fingerprints,
                historical_compounds=historical_compounds,
                train_by_compound=train_by_compound,
                target_to_historical_column=target_to_historical_column,
                historical_target_indices=historical_target_indices,
                popularity=popularity,
                sequence_matrix=sequence_matrix,
                pair_top_columns=pair_top_columns,
                pair_top_similarities=pair_top_similarities,
                target_count=len(target_ids),
            )
            temporal_endpoint_indices = np.asarray(
                [
                    target_index[target]
                    for target in loaded["relevance"]["temporal_strict_ab"][query_id]
                ],
                dtype=np.int32,
            )
            for baseline in BASELINES:
                scores_by_variant = query_scores[baseline]
                require(
                    scores_by_variant.shape == (len(VARIANTS), len(target_ids)),
                    "Score matrix dimensions changed",
                )
                require(np.all(np.isfinite(scores_by_variant)), "A score is nonfinite")
                ranks_by_variant: list[np.ndarray] = []
                for variant_index in range(len(VARIANTS)):
                    _, ranks = rank_scores(
                        scores_by_variant[variant_index],
                        allowed,
                        query_id,
                        target_ids,
                        SALT,
                    )
                    require(
                        np.array_equal(
                            np.sort(ranks[allowed]),
                            np.arange(1, int(allowed.sum()) + 1, dtype=np.int32),
                        )
                        and np.all(ranks[~allowed] == -1),
                        "A complete rank block is not a 1..N permutation",
                    )
                    ranks_by_variant.append(ranks)
                    computed_score_rank_rows += int(allowed.sum())
                primary_scores = scores_by_variant[PRIMARY_VARIANT_INDEX]
                primary_ranks = ranks_by_variant[PRIMARY_VARIANT_INDEX]
                primary_verified_rows += validate_primary_block(
                    primary_reader,
                    query_id,
                    query_compound,
                    baseline,
                    target_ids,
                    allowed,
                    primary_scores,
                    primary_ranks,
                )
                for variant_index, variant in enumerate(VARIANTS):
                    stats = rank_stats[(str(variant["key"]), baseline)]
                    update_rank_change_stats(
                        stats,
                        scores_by_variant[variant_index],
                        ranks_by_variant[variant_index],
                        primary_scores,
                        primary_ranks,
                        allowed,
                        temporal_endpoint_indices,
                    )
                    for scope in DISPLAY_SCOPES:
                        relevant_targets = loaded["relevance"][scope].get(query_id)
                        if not relevant_targets:
                            continue
                        relevant_ranks = [
                            int(ranks_by_variant[variant_index][target_index[target]])
                            for target in relevant_targets
                        ]
                        require(
                            all(rank >= 1 for rank in relevant_ranks),
                            "An endpoint target is masked from its query",
                        )
                        per_query_metrics[
                            (str(variant["key"]), scope, baseline)
                        ].append(query_metrics(relevant_ranks, (10, 50)))
        require(next(primary_reader, None) is None, "Frozen primary rank ledger has extra rows")

    require(
        primary_verified_rows == EXPECTED_PRIMARY_RANK_ROWS,
        "Primary complete-rank row verification count changed",
    )
    require(
        computed_score_rank_rows == EXPECTED_COMPUTED_SCORE_RANK_ROWS,
        "Computed complete score/rank row count changed",
    )

    weight_rows = [
        {
            "weight_variant": variant["key"],
            "tier_A_weight": format_float(float(variant["tier_A_weight"])),
            "tier_B_weight": format_float(float(variant["tier_B_weight"])),
            "role": variant["role"],
            "all_equal_alias": variant["alias"],
            "historical_relation_count": len(loaded["history"]),
            "historical_compound_count": len(train_by_compound),
            "historical_target_count": len(historical_target_indices),
            "endpoint_or_scope_changed": "False",
        }
        for variant in VARIANTS
    ]
    aggregate_rows: list[dict[str, Any]] = []
    aggregate_index: dict[tuple[str, str, str], dict[str, Any]] = {}
    for variant in VARIANTS:
        variant_key = str(variant["key"])
        for scope in DISPLAY_SCOPES:
            scope_rows = loaded["display_scope_rows"][scope]
            for baseline in BASELINES:
                values = per_query_metrics[(variant_key, scope, baseline)]
                expected_queries = len(loaded["relevance"][scope])
                require(
                    len(values) == expected_queries,
                    "Per-query metric denominator changed",
                )
                means = macro_average(values)
                row = {
                    "weight_variant": variant_key,
                    "tier_A_weight": format_float(float(variant["tier_A_weight"])),
                    "tier_B_weight": format_float(float(variant["tier_B_weight"])),
                    "scope": scope,
                    "baseline": baseline,
                    "query_count": len(values),
                    "relation_count": len(scope_rows),
                    **{metric: format_float(means[metric]) for metric in METRICS},
                    "zero_recall_at_50_query_count": sum(
                        item["Recall@50"] == 0.0 for item in values
                    ),
                    "status": "descriptive_fixed_salt",
                }
                aggregate_rows.append(row)
                aggregate_index[(variant_key, scope, baseline)] = row
    require(len(aggregate_rows) == 48, "Aggregate metric matrix is not 48 rows")

    delta_rows: list[dict[str, Any]] = []
    primary_key = str(VARIANTS[PRIMARY_VARIANT_INDEX]["key"])
    for variant in VARIANTS:
        variant_key = str(variant["key"])
        for scope in DISPLAY_SCOPES:
            for baseline in BASELINES:
                observed = aggregate_index[(variant_key, scope, baseline)]
                reference = aggregate_index[(primary_key, scope, baseline)]
                for metric in METRICS:
                    observed_value = float(observed[metric])
                    reference_value = float(reference[metric])
                    delta = observed_value - reference_value
                    if reference_value == 0.0:
                        relative_change = None
                        relative_status = "undefined_reference_zero"
                    else:
                        relative_change = delta / reference_value
                        relative_status = "defined"
                    delta_rows.append(
                        {
                            "weight_variant": variant_key,
                            "reference_variant": primary_key,
                            "scope": scope,
                            "baseline": baseline,
                            "metric": metric,
                            "value": format_float(observed_value),
                            "reference_value": format_float(reference_value),
                            "absolute_difference": format_float(delta),
                            "relative_difference": format_float(relative_change),
                            "relative_difference_status": relative_status,
                            "direction_descriptive": (
                                "increase"
                                if delta > 0
                                else ("decrease" if delta < 0 else "no_change")
                            ),
                        }
                    )
    require(len(delta_rows) == 240, "Metric-delta matrix is not 240 rows")

    rank_change_rows: list[dict[str, Any]] = []
    for variant in VARIANTS:
        variant_key = str(variant["key"])
        for baseline in BASELINES:
            stats = rank_stats[(variant_key, baseline)]
            eligible = stats["eligible_candidate_rows"]
            changed = stats["rank_changed_candidate_count"]
            rank_change_rows.append(
                {
                    "weight_variant": variant_key,
                    "reference_variant": primary_key,
                    "baseline": baseline,
                    "eligible_candidate_rows": eligible,
                    "score_changed_candidate_count": stats[
                        "score_changed_candidate_count"
                    ],
                    "score_changed_candidate_fraction": format_float(
                        stats["score_changed_candidate_count"] / eligible
                    ),
                    "rank_changed_candidate_count": changed,
                    "rank_changed_candidate_fraction": format_float(changed / eligible),
                    "mean_absolute_rank_change": format_float(
                        stats["absolute_rank_change_sum"] / eligible
                    ),
                    "maximum_absolute_rank_change": stats[
                        "maximum_absolute_rank_change"
                    ],
                    "query_count_with_any_rank_change": stats[
                        "query_count_with_any_rank_change"
                    ],
                    "top50_symmetric_difference_membership_count": stats[
                        "top50_symmetric_difference_membership_count"
                    ],
                    "query_count_with_any_top50_membership_change": stats[
                        "query_count_with_any_top50_membership_change"
                    ],
                    "endpoint_relation_count": stats["endpoint_relation_count"],
                    "endpoint_relation_rank_changed_count": stats[
                        "endpoint_relation_rank_changed_count"
                    ],
                    "endpoint_relation_top50_membership_changed_count": stats[
                        "endpoint_relation_top50_membership_changed_count"
                    ],
                    "rank_permutation_blocks_checked": stats[
                        "rank_permutation_blocks_checked"
                    ],
                }
            )
    require(len(rank_change_rows) == 12, "Rank-change matrix is not 12 rows")
    scope_rows = scope_cardinality_rows(loaded["display_scope_rows"])
    require(len(scope_rows) == 12, "Scope-invariance matrix is not 12 rows")

    output_dir.mkdir(parents=True, exist_ok=False)
    weight_path = output_dir / "weight_variants.tsv"
    aggregate_path = output_dir / "aggregate_metrics.tsv"
    delta_path = output_dir / "metric_deltas_vs_0_7.tsv"
    rank_change_path = output_dir / "complete_rank_top50_changes_vs_0_7.tsv"
    scope_path = output_dir / "scope_cardinality_invariance.tsv"
    summary_path = output_dir / "weight_policy_summary.json"
    receipt_path = output_dir / "execution_receipt.json"
    manifest_path = output_dir / "run_manifest.json"
    write_tsv_new(weight_path, list(weight_rows[0]), weight_rows)
    write_tsv_new(aggregate_path, list(aggregate_rows[0]), aggregate_rows)
    write_tsv_new(delta_path, list(delta_rows[0]), delta_rows)
    write_tsv_new(rank_change_path, list(rank_change_rows[0]), rank_change_rows)
    write_tsv_new(scope_path, list(scope_rows[0]), scope_rows)

    recall50_rows = [
        {
            "weight_variant": row["weight_variant"],
            "scope": row["scope"],
            "baseline": row["baseline"],
            "Recall@50": row["Recall@50"],
            "difference_vs_0_7": next(
                delta["absolute_difference"]
                for delta in delta_rows
                if delta["weight_variant"] == row["weight_variant"]
                and delta["scope"] == row["scope"]
                and delta["baseline"] == row["baseline"]
                and delta["metric"] == "Recall@50"
            ),
        }
        for row in aggregate_rows
    ]
    summary = {
        "schema_version": "1.0",
        "analysis_id": ANALYSIS_ID,
        "protocol_id": PROTOCOL_ID,
        "analysis_role": "reviewer_requested_post_hoc_descriptive_policy_sensitivity",
        "claim_boundary": (
            "Author-run, outcome-visible, non-independent retrospective sensitivity; "
            "no tuning, superiority, validation, or biological claim."
        ),
        "weight_variants": VARIANTS,
        "all_equal_alias_verified": {
            "variant": "A1_B1_0_all_equal",
            "tier_A_weight": 1.0,
            "tier_B_weight": 1.0,
            "additional_variant": False,
        },
        "complete_score_rank_contract": {
            "computed_score_rank_rows": computed_score_rank_rows,
            "primary_0_7_rows_reproduced_exactly": primary_verified_rows,
            "complete_rank_ledgers_written": False,
            "rank_permutation_blocks_checked": sum(
                row["rank_permutation_blocks_checked"] for row in rank_change_rows
            ),
        },
        "reuse_contract": {
            "chemical_similarity_computations_per_query": 1,
            "sequence_matrix_computations": 1,
            "pair_neighbor_top100_computations": 1,
            "query_sequence_similarity_computations_per_query": 1,
            "weights_evaluated_in_same_process": True,
        },
        "scope_invariance": {
            "all_cardinality_differences_zero": True,
            "homology_0_50_and_0_70_inputs_identical": True,
            "identical_masks_displayed_once": True,
        },
        "recall_at_50": recall50_rows,
        "row_counts": {
            "weight_variants": len(weight_rows),
            "aggregate_metrics": len(aggregate_rows),
            "metric_deltas": len(delta_rows),
            "complete_rank_top50_changes": len(rank_change_rows),
            "scope_cardinality": len(scope_rows),
        },
        "output_boundary": {
            "aggregate_only": True,
            "identifier_bearing_rows": False,
            "complete_rank_rows_written": False,
            "absolute_paths_written": False,
        },
        "uncertainty_boundary": (
            "No bootstrap or p-value is added. Fixed-salt policy perturbations are "
            "descriptive; exact-tie uncertainty is reported separately."
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
            "rdkit": rdBase.rdkitVersion,
            "sklearn": sklearn_version,
            "python_tracemalloc_current_bytes": current_memory,
            "python_tracemalloc_peak_bytes": peak_memory,
        },
        "locks": {
            "plan": {
                "basename": args.plan_lock.name,
                "sha256": sha256(args.plan_lock),
            },
            "parent_protocol": {
                "basename": args.parent_protocol.name,
                "sha256": sha256(args.parent_protocol),
            },
            "implementation": {
                "basename": args.implementation_lock.name,
                "sha256": sha256(args.implementation_lock),
            },
        },
        "preexecution_synthetic_test_receipt": implementation.get(
            "synthetic_test_receipt"
        ),
        "complete_computation": {
            "score_rank_rows_computed": computed_score_rank_rows,
            "primary_score_rank_rows_exactly_verified": primary_verified_rows,
            "complete_rank_rows_written": 0,
        },
        "reuse_verified": {
            "single_process": True,
            "chemical_similarities_reused_across_weights": True,
            "sequence_features_reused_across_weights": True,
            "pair_neighbor_top100_reused_across_weights": True,
        },
        "deviations_from_plan": [],
        "identifier_bearing_output_written": False,
        "absolute_paths_written": False,
        "external_transfer_performed": False,
        "outputs_before_receipt": [
            output_descriptor(path)
            for path in (
                weight_path,
                aggregate_path,
                delta_path,
                rank_change_path,
                scope_path,
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
            {"relative_name": name, "sha256": digest}
            for name, digest in sorted(
                implementation["implementation_files"].items()
            )
        ],
        "outputs": [
            output_descriptor(path)
            for path in (
                weight_path,
                aggregate_path,
                delta_path,
                rank_change_path,
                scope_path,
                summary_path,
                receipt_path,
            )
        ],
        "output_contract": {
            "weight_variant_rows": len(weight_rows),
            "aggregate_metric_rows": len(aggregate_rows),
            "metric_delta_rows": len(delta_rows),
            "complete_rank_change_rows": len(rank_change_rows),
            "scope_cardinality_rows": len(scope_rows),
        },
        "created_at_utc": utc_now(),
    }
    write_json_new(manifest_path, manifest)
    print(
        json.dumps(
            {
                "analysis_id": ANALYSIS_ID,
                "status": "completed",
                "computed_score_rank_rows": computed_score_rank_rows,
                "primary_rows_reproduced": primary_verified_rows,
                "identifier_bearing_output": False,
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
