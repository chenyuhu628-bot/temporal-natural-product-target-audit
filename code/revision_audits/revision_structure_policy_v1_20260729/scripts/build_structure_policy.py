#!/usr/bin/env python
"""Run aggregate-only RDKit structure-policy sensitivity analysis."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.Scaffolds import MurckoScaffold


OUTPUT_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCORER_DIR = PROJECT_ROOT / "strict_ab_asof_cutoff_successor_v1_20260728/scripts"
SHARED_SCRIPT_DIR = PROJECT_ROOT / "scripts"
sys.path.insert(0, str(SCORER_DIR))
sys.path.insert(0, str(SHARED_SCRIPT_DIR))

from asof_successor_common import (  # noqa: E402
    BASELINES,
    LEGACY_TIE_SALT,
    MORGAN_BITS,
    MORGAN_RADIUS,
    STRICT_TIERS,
)
from pu_retrieval_metrics import rank_scores  # noqa: E402
from score_asof_cutoff_successor import (  # noqa: E402
    build_sequence_matrix,
    build_train_maps,
    morgan_fingerprints,
    pair_neighbor_scores,
    precompute_sequence_topk,
    sequence_transfer_scores,
    weighted_morgan_transfer_scores,
)


PARENT_PROTOCOL_SHA256 = "bd4902476160cc7c5cbacaf0cfd0f1a28c5300bde22232b3a2cc6c1f3c143dc2"
EXPECTED_RDKIT = "2026.03.4"
POLICIES = [
    "raw_primary",
    "cleanup_fragment_parent",
    "cleanup_charge_normalized",
    "cleanup_canonical_tautomer",
    "cleanup_parent_charge_tautomer",
]
INVARIANT_BASELINES = {"weighted_target_popularity", "sequence_3mer_transfer"}
SCOPES = [
    "temporal_strict_ab",
    "scaffold_cold",
    "joint_scaffold_homology_0_30",
    "joint_scaffold_homology_0_50_0_70_identical",
]


EXPECTED_INPUT_HASHES = {
    "parent_protocol": PARENT_PROTOCOL_SHA256,
    "historical_pairs": "cef748ae8ac277e49784d7e1fbf08e085beb14a84fc6c651a0fb8d99e88710d7",
    "scoring_queries": "0e6068d2e25cb3ea325656fb3517563788cd496e88cfaa3de761890fec9e9318",
    "historical_compounds": "f1f82793b5c652007a042699c19cd5640a8b68e8b1f0d4f94e4dc4f54045060c",
    "query_compounds": "f51670fffd21d2e9109b4376dd53aab55bbacb09d2e4795dfa515fcaca98b113",
    "candidate_targets": "0ee86746b306fb388a1f74a6b88ce4d1eba01b7a4eb473315f6b3def57145cdc",
    "candidate_sequences": "a83421dba2482f236fe18340dd592cc7d5ed22c98c4fc39435c40f04f289b442",
    "evaluation_pairs": "09296b066a23197a7c178f00514f2b3d9ed7e6f3c459ea92a55e01a6010d1132",
    "scaffold_audit": "fa0029ef5b7822ad5ca93f7bd93ac808f85f1e0c02e827fa91be375031b2d7af",
    "homology_0_30": "3a8247ed8f683fe6fce5fb345f56e3ec73a872b065eca922e92e494f084a1793",
    "homology_0_50": "ec3bbd435f73bc1c724efdfd94ac10c32d6f9a55bd6c9a9349233a01e5dc7f5b",
    "homology_0_70": "ec3bbd435f73bc1c724efdfd94ac10c32d6f9a55bd6c9a9349233a01e5dc7f5b",
    "frozen_ranks": "87739aa818744c7084088d13c386444aa41bbef38c257083325298003181479e",
    "frozen_metrics": "fac75f28185cc4e2d320ab45ddba5d07d4857c85ce477e62ec4ae960ff656ea8",
    "frozen_scorer": "7b8263828ab6eaf4756307f315df4502b1615a40bffa8219284612585ad9bdc8",
    "frozen_metrics_code": "7e40e80c320ab45e2d320ab45ddba5d07d4857c85ce477e62ec4ae960ff656ea8",
}
# Correct the code-lock value separately so the dictionary is explicit and easy to audit.
EXPECTED_INPUT_HASHES["frozen_metrics_code"] = (
    "7e40e80c3203a4a6cf95ca675cc9c333168f8feb234449856f26b5e132ec3165"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return path.open("rt", encoding="utf-8", newline="")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    current: str | None = None
    chunks: list[str] = []
    with path.open("rt", encoding="utf-8", newline="") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current is not None:
                    sequences[current] = "".join(chunks)
                current = line[1:].split()[0]
                chunks = []
            else:
                if current is None:
                    raise ValueError("Sequence before FASTA header")
                chunks.append(line)
    if current is not None:
        sequences[current] = "".join(chunks)
    return sequences


def write_tsv(path: Path, fields: list[str], rows: Iterable[dict[str, Any]]) -> None:
    with path.open("wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def parse_bool(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no"}:
        return False
    raise ValueError(f"Invalid Boolean: {value!r}")


def transform_molecule(molecule: Chem.Mol, policy: str) -> Chem.Mol:
    result = Chem.Mol(molecule)
    if policy == "raw_primary":
        pass
    elif policy == "cleanup_fragment_parent":
        result = rdMolStandardize.Cleanup(result)
        result = rdMolStandardize.FragmentParent(result, skipStandardize=True)
    elif policy == "cleanup_charge_normalized":
        result = rdMolStandardize.Cleanup(result)
        result = rdMolStandardize.Uncharger().uncharge(result)
    elif policy == "cleanup_canonical_tautomer":
        result = rdMolStandardize.Cleanup(result)
        result = rdMolStandardize.TautomerEnumerator().Canonicalize(result)
    elif policy == "cleanup_parent_charge_tautomer":
        result = rdMolStandardize.Cleanup(result)
        result = rdMolStandardize.FragmentParent(result, skipStandardize=True)
        result = rdMolStandardize.Uncharger().uncharge(result)
        result = rdMolStandardize.TautomerEnumerator().Canonicalize(result)
    else:
        raise ValueError(f"Unknown policy: {policy}")
    Chem.SanitizeMol(result)
    if result.GetNumAtoms() == 0:
        raise ValueError("Standardization produced an empty molecule")
    return result


def molecule_record(molecule: Chem.Mol, generator: Any) -> dict[str, Any]:
    canonical = Chem.MolToSmiles(molecule, canonical=True, isomericSmiles=True)
    if not canonical:
        raise ValueError("Canonical SMILES generation failed")
    inchi_key = Chem.MolToInchiKey(molecule)
    if not inchi_key or len(inchi_key) < 14:
        raise ValueError("InChIKey generation failed")
    scaffold_molecule = MurckoScaffold.GetScaffoldForMol(molecule)
    scaffold = (
        Chem.MolToSmiles(scaffold_molecule, canonical=True, isomericSmiles=True)
        if scaffold_molecule.GetNumAtoms()
        else ""
    )
    return {
        "mol": molecule,
        "canonical": canonical,
        "inchi_key": inchi_key,
        "connectivity": inchi_key[:14],
        "scaffold": scaffold,
        "fingerprint": generator.GetFingerprint(molecule),
    }


def collision_summary(values: list[str]) -> tuple[int, int, int]:
    counts = Counter(values)
    return len(counts), sum(count > 1 for count in counts.values()), sum(
        count - 1 for count in counts.values() if count > 1
    )


def transform_role(
    rows: list[dict[str, str]],
    role: str,
    generator: Any,
) -> tuple[dict[str, dict[str, dict[str, Any]]], list[dict[str, Any]], list[dict[str, str]]]:
    raw_parsed: dict[str, Chem.Mol] = {}
    parse_failures = 0
    for row in rows:
        key = row["inchikey_full"]
        molecule = Chem.MolFromSmiles(row["representative_smiles"])
        if molecule is None:
            parse_failures += 1
        else:
            raw_parsed[key] = molecule
    transformed: dict[str, dict[str, dict[str, Any]]] = {}
    summaries: list[dict[str, Any]] = []
    blockers: list[dict[str, str]] = []
    raw_records: dict[str, dict[str, Any]] = {}
    for policy in POLICIES:
        records: dict[str, dict[str, Any]] = {}
        transform_failures = 0
        for row in rows:
            key = row["inchikey_full"]
            molecule = raw_parsed.get(key)
            if molecule is None:
                continue
            try:
                records[key] = molecule_record(
                    transform_molecule(molecule, policy), generator
                )
            except Exception:
                transform_failures += 1
        transformed[policy] = records
        if policy == "raw_primary":
            raw_records = records
        canonical_values = [record["canonical"] for record in records.values()]
        inchi_values = [record["inchi_key"] for record in records.values()]
        connectivity_values = [record["connectivity"] for record in records.values()]
        scaffold_values = [record["scaffold"] for record in records.values() if record["scaffold"]]
        canonical_distinct, canonical_collision_groups, canonical_collision_excess = collision_summary(
            canonical_values
        )
        inchi_distinct, inchi_collision_groups, inchi_collision_excess = collision_summary(
            inchi_values
        )
        connectivity_distinct, connectivity_collision_groups, connectivity_collision_excess = collision_summary(
            connectivity_values
        )
        comparable = set(records).intersection(raw_records) if raw_records else set()
        canonical_changed = sum(
            records[key]["canonical"] != raw_records[key]["canonical"] for key in comparable
        )
        fingerprint_changed = sum(
            records[key]["fingerprint"].ToBitString()
            != raw_records[key]["fingerprint"].ToBitString()
            for key in comparable
        )
        scaffold_changed = sum(
            records[key]["scaffold"] != raw_records[key]["scaffold"] for key in comparable
        )
        derived_full_changed = sum(
            records[key]["inchi_key"] != raw_records[key]["inchi_key"] for key in comparable
        )
        derived_connectivity_changed = sum(
            records[key]["connectivity"] != raw_records[key]["connectivity"]
            for key in comparable
        )
        input_full_changed = sum(
            key in records and records[key]["inchi_key"] != key for key in raw_parsed
        )
        input_connectivity_changed = sum(
            key in records and records[key]["connectivity"] != key[:14]
            for key in raw_parsed
        )
        status = (
            "complete"
            if parse_failures == 0 and transform_failures == 0 and len(records) == len(rows)
            else "blocked_no_imputation"
        )
        if status != "complete":
            blockers.append(
                {
                    "policy": policy,
                    "role": role,
                    "blocker": (
                        "Role-specific structure parsing, transformation, sanitization, "
                        "or InChIKey generation failed; downstream policy scoring is not estimable."
                    ),
                    "failed_record_count": str(parse_failures + transform_failures),
                    "imputation_performed": "false",
                }
            )
        summaries.append(
            {
                "policy": policy,
                "role": role,
                "status": status,
                "input_record_count": len(rows),
                "parse_success_count": len(raw_parsed),
                "parse_failure_count": parse_failures,
                "transform_success_count": len(records),
                "transform_failure_count": transform_failures,
                "canonical_structure_changed_vs_raw": canonical_changed,
                "morgan_fingerprint_changed_vs_raw": fingerprint_changed,
                "full_inchikey_changed_vs_raw_derived": derived_full_changed,
                "connectivity_layer_changed_vs_raw_derived": derived_connectivity_changed,
                "full_inchikey_changed_vs_input_identifier": input_full_changed,
                "connectivity_layer_changed_vs_input_identifier": input_connectivity_changed,
                "scaffold_changed_vs_raw": scaffold_changed,
                "nonempty_scaffold_count": len(scaffold_values),
                "empty_or_acyclic_scaffold_count": len(records) - len(scaffold_values),
                "distinct_canonical_structure_count": canonical_distinct,
                "canonical_collision_group_count": canonical_collision_groups,
                "canonical_collision_excess_count": canonical_collision_excess,
                "distinct_full_inchikey_count": inchi_distinct,
                "full_inchikey_collision_group_count": inchi_collision_groups,
                "full_inchikey_collision_excess_count": inchi_collision_excess,
                "distinct_connectivity_layer_count": connectivity_distinct,
                "connectivity_collision_group_count": connectivity_collision_groups,
                "connectivity_collision_excess_count": connectivity_collision_excess,
            }
        )
    return transformed, summaries, blockers


def make_scopes(
    policy: str,
    evaluation_rows: list[dict[str, str]],
    historical_records: dict[str, dict[str, Any]],
    query_records: dict[str, dict[str, Any]],
    homology_cold: dict[str, set[str]],
) -> tuple[dict[str, list[dict[str, str]]], dict[str, Any]]:
    historical_scaffolds = {
        record["scaffold"] for record in historical_records.values() if record["scaffold"]
    }
    scaffold_rows: list[dict[str, str]] = []
    empty_rows: list[dict[str, str]] = []
    for row in evaluation_rows:
        scaffold = query_records[row["inchikey_full"]]["scaffold"]
        if not scaffold:
            empty_rows.append(row)
        elif scaffold not in historical_scaffolds:
            scaffold_rows.append(row)
    scopes = {
        "temporal_strict_ab": list(evaluation_rows),
        "scaffold_cold": scaffold_rows,
        "joint_scaffold_homology_0_30": [
            row
            for row in scaffold_rows
            if row["uniprot_canonical_accession"] in homology_cold["0_30"]
        ],
        "joint_scaffold_homology_0_50_0_70_identical": [
            row
            for row in scaffold_rows
            if row["uniprot_canonical_accession"] in homology_cold["0_50"]
        ],
    }
    joined_070 = {
        row["canonical_pair_key"]
        for row in scaffold_rows
        if row["uniprot_canonical_accession"] in homology_cold["0_70"]
    }
    joined_050 = {
        row["canonical_pair_key"]
        for row in scopes["joint_scaffold_homology_0_50_0_70_identical"]
    }
    if joined_050 != joined_070:
        raise AssertionError(f"Policy {policy} produces different 0.50 and 0.70 joins")
    details = {
        "historical_nonempty_scaffold_count": len(historical_scaffolds),
        "query_nonempty_scaffold_count": sum(
            bool(record["scaffold"]) for record in query_records.values()
        ),
        "query_empty_scaffold_count": sum(
            not bool(record["scaffold"]) for record in query_records.values()
        ),
        "endpoint_empty_scaffold_relation_count": len(empty_rows),
        "endpoint_empty_scaffold_query_count": len({row["query_id"] for row in empty_rows}),
    }
    return scopes, details


def recall_at_50(
    ranks: np.ndarray,
    scope_rows: list[dict[str, str]],
    query_index: dict[str, int],
    target_index: dict[str, int],
) -> tuple[float | None, int, int, int]:
    by_query: dict[str, list[int]] = defaultdict(list)
    for row in scope_rows:
        rank = int(
            ranks[
                query_index[row["query_id"]],
                target_index[row["uniprot_canonical_accession"]],
            ]
        )
        if rank < 1:
            raise AssertionError("Endpoint target is masked or unranked")
        by_query[row["query_id"]].append(rank)
    if not by_query:
        return None, 0, 0, 0
    query_values = [
        sum(rank <= 50 for rank in positive_ranks) / len(positive_ranks)
        for positive_ranks in by_query.values()
    ]
    zero_queries = sum(value == 0.0 for value in query_values)
    relation_hits = sum(rank <= 50 for values in by_query.values() for rank in values)
    return float(np.mean(query_values)), len(by_query), zero_queries, relation_hits


def fmt_optional(value: float | None) -> str:
    return "" if value is None else f"{value:.17g}"


def build_raw_scores(
    queries: list[dict[str, str]],
    target_ids: list[str],
    train_by_compound: dict[str, list[tuple[int, float]]],
    popularity: np.ndarray,
    raw_historical_fps: dict[str, Any],
    raw_query_fps: dict[str, Any],
    sequence_matrix: Any,
    historical_compounds_order: list[str],
    historical_column: dict[int, int],
    top_columns: np.ndarray,
    top_similarities: np.ndarray,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], np.ndarray]:
    query_count, target_count = len(queries), len(target_ids)
    scores = {
        baseline: np.empty((query_count, target_count), dtype=np.float32)
        for baseline in BASELINES
    }
    ranks = {
        baseline: np.empty((query_count, target_count), dtype=np.int32)
        for baseline in BASELINES
    }
    allowed_matrix = np.ones((query_count, target_count), dtype=bool)
    for query_position, query_row in enumerate(queries):
        compound = query_row["inchikey_full"]
        allowed = allowed_matrix[query_position]
        for target_position, _ in train_by_compound.get(compound, []):
            allowed[target_position] = False
        query_fp = raw_query_fps[compound]
        current = {
            "weighted_target_popularity": popularity.copy(),
            "sequence_3mer_transfer": sequence_transfer_scores(
                compound, train_by_compound, sequence_matrix, target_count
            ),
            "weighted_morgan_transfer": weighted_morgan_transfer_scores(
                query_fp,
                historical_compounds_order,
                train_by_compound,
                raw_historical_fps,
                target_count,
            ),
            "structure_sequence_pair_neighbor": pair_neighbor_scores(
                query_fp,
                historical_compounds_order,
                train_by_compound,
                raw_historical_fps,
                historical_column,
                top_columns,
                top_similarities,
            ),
        }
        for baseline in BASELINES:
            values = np.asarray(current[baseline], dtype=np.float32)
            _, rank_values = rank_scores(
                values, allowed, query_row["query_id"], target_ids, LEGACY_TIE_SALT
            )
            scores[baseline][query_position] = values
            ranks[baseline][query_position] = rank_values
    return scores, ranks, allowed_matrix


def calibrate_rank_ledger(
    rank_path: Path,
    queries: list[dict[str, str]],
    target_ids: list[str],
    raw_scores: dict[str, np.ndarray],
    raw_ranks: dict[str, np.ndarray],
    allowed_matrix: np.ndarray,
) -> dict[str, int]:
    counters = Counter()
    with gzip.open(rank_path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for query_position, query_row in enumerate(queries):
            for baseline in BASELINES:
                for target_position, target_id in enumerate(target_ids):
                    if not allowed_matrix[query_position, target_position]:
                        continue
                    row = next(reader, None)
                    if row is None:
                        counters["missing_frozen_rows"] += 1
                        continue
                    counters["cells_checked"] += 1
                    if (
                        row["query_id"] != query_row["query_id"]
                        or row["baseline"] != baseline
                        or row["target_uniprot_accession"] != target_id
                    ):
                        counters["ordering_or_key_mismatch"] += 1
                    expected_score = f"{float(raw_scores[baseline][query_position, target_position]):.17g}"
                    if row["score"] != expected_score:
                        counters["score_string_mismatch"] += 1
                    if int(row["rank"]) != int(
                        raw_ranks[baseline][query_position, target_position]
                    ):
                        counters["rank_mismatch"] += 1
                    if int(row["eligible_candidate_target_count"]) != int(
                        allowed_matrix[query_position].sum()
                    ):
                        counters["candidate_count_mismatch"] += 1
        if next(reader, None) is not None:
            counters["extra_frozen_rows"] += 1
    for key in (
        "missing_frozen_rows",
        "ordering_or_key_mismatch",
        "score_string_mismatch",
        "rank_mismatch",
        "candidate_count_mismatch",
        "extra_frozen_rows",
    ):
        counters.setdefault(key, 0)
    if any(counters[key] for key in counters if key != "cells_checked"):
        raise AssertionError("Raw score/rank ledger calibration failed closed")
    return dict(counters)


def score_policy(
    queries: list[dict[str, str]],
    target_ids: list[str],
    train_by_compound: dict[str, list[tuple[int, float]]],
    historical_fps: dict[str, Any],
    query_fps: dict[str, Any],
    historical_compounds_order: list[str],
    historical_column: dict[int, int],
    top_columns: np.ndarray,
    top_similarities: np.ndarray,
    allowed_matrix: np.ndarray,
    raw_scores: dict[str, np.ndarray],
    raw_ranks: dict[str, np.ndarray],
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    scores = {baseline: raw_scores[baseline] for baseline in INVARIANT_BASELINES}
    ranks = {baseline: raw_ranks[baseline] for baseline in INVARIANT_BASELINES}
    target_count = len(target_ids)
    for baseline in set(BASELINES).difference(INVARIANT_BASELINES):
        scores[baseline] = np.empty((len(queries), target_count), dtype=np.float32)
        ranks[baseline] = np.empty((len(queries), target_count), dtype=np.int32)
    for query_position, query_row in enumerate(queries):
        compound = query_row["inchikey_full"]
        query_fp = query_fps[compound]
        current = {
            "weighted_morgan_transfer": weighted_morgan_transfer_scores(
                query_fp,
                historical_compounds_order,
                train_by_compound,
                historical_fps,
                target_count,
            ),
            "structure_sequence_pair_neighbor": pair_neighbor_scores(
                query_fp,
                historical_compounds_order,
                train_by_compound,
                historical_fps,
                historical_column,
                top_columns,
                top_similarities,
            ),
        }
        for baseline, values in current.items():
            values = np.asarray(values, dtype=np.float32)
            _, rank_values = rank_scores(
                values,
                allowed_matrix[query_position],
                query_row["query_id"],
                target_ids,
                LEGACY_TIE_SALT,
            )
            scores[baseline][query_position] = values
            ranks[baseline][query_position] = rank_values
    return scores, ranks


def rank_change_row(
    policy: str,
    baseline: str,
    scores: np.ndarray,
    ranks: np.ndarray,
    raw_scores: np.ndarray,
    raw_ranks: np.ndarray,
    allowed_matrix: np.ndarray,
) -> dict[str, Any]:
    eligible_scores = scores[allowed_matrix]
    raw_eligible_scores = raw_scores[allowed_matrix]
    eligible_ranks = ranks[allowed_matrix]
    raw_eligible_ranks = raw_ranks[allowed_matrix]
    score_changed = eligible_scores != raw_eligible_scores
    rank_changed = eligible_ranks != raw_eligible_ranks
    top50_changed = (eligible_ranks <= 50) != (raw_eligible_ranks <= 50)
    score_delta = np.abs(
        eligible_scores.astype(np.float64) - raw_eligible_scores.astype(np.float64)
    )
    rank_delta = np.abs(
        eligible_ranks.astype(np.int64) - raw_eligible_ranks.astype(np.int64)
    )
    query_score_changed = 0
    query_rank_changed = 0
    query_top50_changed = 0
    for index in range(allowed_matrix.shape[0]):
        allowed = allowed_matrix[index]
        query_score_changed += bool(np.any(scores[index, allowed] != raw_scores[index, allowed]))
        query_rank_changed += bool(np.any(ranks[index, allowed] != raw_ranks[index, allowed]))
        query_top50_changed += bool(
            np.any((ranks[index, allowed] <= 50) != (raw_ranks[index, allowed] <= 50))
        )
    total = int(allowed_matrix.sum())
    return {
        "policy": policy,
        "baseline": baseline,
        "status": "complete",
        "eligible_rank_cell_count": total,
        "score_changed_cell_count": int(score_changed.sum()),
        "score_changed_cell_fraction": f"{float(score_changed.mean()):.17g}",
        "rank_changed_cell_count": int(rank_changed.sum()),
        "rank_changed_cell_fraction": f"{float(rank_changed.mean()):.17g}",
        "top50_membership_changed_cell_count": int(top50_changed.sum()),
        "query_count_with_any_score_change": query_score_changed,
        "query_count_with_any_rank_change": query_rank_changed,
        "query_count_with_any_top50_membership_change": query_top50_changed,
        "mean_absolute_rank_change_all_cells": f"{float(rank_delta.mean()):.17g}",
        "mean_absolute_rank_change_changed_cells": (
            f"{float(rank_delta[rank_changed].mean()):.17g}" if np.any(rank_changed) else "0"
        ),
        "maximum_absolute_rank_change": int(rank_delta.max()),
        "maximum_absolute_score_change": f"{float(score_delta.max()):.17g}",
        "invariant_by_design": str(baseline in INVARIANT_BASELINES).lower(),
    }


def endpoint_change_counts(
    current_scores: np.ndarray,
    current_ranks: np.ndarray,
    raw_scores: np.ndarray,
    raw_ranks: np.ndarray,
    rows: list[dict[str, str]],
    query_index: dict[str, int],
    target_index: dict[str, int],
) -> dict[str, int]:
    score_changed = rank_changed = top50_changed = 0
    changed_queries: set[str] = set()
    for row in rows:
        q = query_index[row["query_id"]]
        t = target_index[row["uniprot_canonical_accession"]]
        score_changed += current_scores[q, t] != raw_scores[q, t]
        rank_changed += current_ranks[q, t] != raw_ranks[q, t]
        changed = (current_ranks[q, t] <= 50) != (raw_ranks[q, t] <= 50)
        top50_changed += changed
        if changed:
            changed_queries.add(row["query_id"])
    return {
        "endpoint_relation_score_changed_count": score_changed,
        "endpoint_relation_rank_changed_count": rank_changed,
        "endpoint_relation_top50_membership_changed_count": top50_changed,
        "query_count_with_endpoint_top50_membership_change": len(changed_queries),
    }


def main() -> None:
    started = time.perf_counter()
    if rdBase.rdkitVersion != EXPECTED_RDKIT:
        raise RuntimeError(
            f"RDKit version mismatch: required {EXPECTED_RDKIT}, found {rdBase.rdkitVersion}"
        )
    paths = {
        "parent_protocol": PROJECT_ROOT
        / "manuscript_molecular_diversity_v3_20260728/plan/revision_analysis_protocol_v4_20260729.md",
        "historical_pairs": PROJECT_ROOT
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/scoring_inputs/historical_pairs.tsv.gz",
        "scoring_queries": PROJECT_ROOT
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/scoring_inputs/scoring_queries.tsv.gz",
        "historical_compounds": PROJECT_ROOT
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/scoring_inputs/historical_compounds.tsv.gz",
        "query_compounds": PROJECT_ROOT
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/scoring_inputs/query_compounds.tsv.gz",
        "candidate_targets": PROJECT_ROOT
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/scoring_inputs/candidate_targets.tsv.gz",
        "candidate_sequences": PROJECT_ROOT
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/scoring_inputs/candidate_sequences.fasta",
        "evaluation_pairs": PROJECT_ROOT
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/evaluation_inputs/evaluation_pairs.tsv.gz",
        "scaffold_audit": PROJECT_ROOT
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/evaluation_inputs/scaffold_audit.tsv.gz",
        "homology_0_30": PROJECT_ROOT
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/evaluation_inputs/homology_0_30.tsv.gz",
        "homology_0_50": PROJECT_ROOT
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/evaluation_inputs/homology_0_50.tsv.gz",
        "homology_0_70": PROJECT_ROOT
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/evaluation_inputs/homology_0_70.tsv.gz",
        "frozen_ranks": PROJECT_ROOT
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/score/corrective_prediction_ranks.tsv.gz",
        "frozen_metrics": PROJECT_ROOT
        / "author_run_strict_ab_asof_cutoff_execution_v1_20260728/evaluation/corrective_aggregate_metrics.tsv.gz",
        "frozen_scorer": PROJECT_ROOT
        / "strict_ab_asof_cutoff_successor_v1_20260728/scripts/score_asof_cutoff_successor.py",
        "frozen_metrics_code": PROJECT_ROOT / "scripts/pu_retrieval_metrics.py",
    }
    protocol_path = OUTPUT_DIR / "PROTOCOL.md"
    script_path = Path(__file__).resolve()
    validator_path = OUTPUT_DIR / "scripts/validate_structure_policy.py"
    generated = [
        "structure_policy_summary.tsv",
        "scaffold_scope_changes.tsv",
        "rank_change_summary.tsv",
        "scope_recall_at_50.tsv",
        "blockers.tsv",
        "calibration_summary.json",
        "input_hashes.json",
        "execution_receipt.json",
        "test_results.json",
        "manifest.json",
        "validation_report.json",
    ]
    existing = [name for name in generated if (OUTPUT_DIR / name).exists()]
    if existing:
        raise FileExistsError(
            "Create-once output already exists: " + ", ".join(existing)
        )
    for label, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(f"Missing locked input: {label}")
        actual = sha256(path)
        if actual != EXPECTED_INPUT_HASHES[label]:
            raise AssertionError(f"Locked input hash mismatch: {label}")

    history = read_tsv(paths["historical_pairs"])
    queries = read_tsv(paths["scoring_queries"])
    historical_rows = read_tsv(paths["historical_compounds"])
    query_rows = read_tsv(paths["query_compounds"])
    target_ids = [
        row["uniprot_canonical_accession"]
        for row in read_tsv(paths["candidate_targets"])
    ]
    sequences = read_fasta(paths["candidate_sequences"])
    evaluation_rows = read_tsv(paths["evaluation_pairs"])
    frozen_scaffold_rows = read_tsv(paths["scaffold_audit"])
    frozen_metrics_rows = read_tsv(paths["frozen_metrics"])
    if len(history) != 4990 or len(queries) != 222 or len(target_ids) != 4123:
        raise AssertionError("Locked cardinalities changed")
    if len(historical_rows) != 1726 or len(query_rows) != 222 or len(evaluation_rows) != 358:
        raise AssertionError("Locked role or endpoint cardinalities changed")
    if set(target_ids) != set(sequences):
        raise AssertionError("Target and sequence keysets differ")
    if {row["inchikey_full"] for row in historical_rows} != {
        row["inchikey_full"] for row in history
    }:
        raise AssertionError("Historical role structure keyset differs")
    if {row["inchikey_full"] for row in query_rows} != {
        row["inchikey_full"] for row in queries
    }:
        raise AssertionError("Query role structure keyset differs")

    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=MORGAN_RADIUS, fpSize=MORGAN_BITS
    )
    historical_transformed, historical_summary, historical_blockers = transform_role(
        historical_rows, "historical", generator
    )
    query_transformed, query_summary, query_blockers = transform_role(
        query_rows, "query", generator
    )
    structure_summary = historical_summary + query_summary
    blockers = historical_blockers + query_blockers
    policy_status = {
        policy: (
            "complete"
            if len(historical_transformed[policy]) == len(historical_rows)
            and len(query_transformed[policy]) == len(query_rows)
            else "blocked_no_imputation"
        )
        for policy in POLICIES
    }
    if policy_status["raw_primary"] != "complete":
        raise AssertionError("Raw role structures failed; calibration cannot proceed")

    homology_cold: dict[str, set[str]] = {}
    for threshold in ("0_30", "0_50", "0_70"):
        homology_cold[threshold] = {
            row["uniprot_canonical_accession"]
            for row in read_tsv(paths[f"homology_{threshold}"])
            if parse_bool(row["is_future_target_homology_cold_candidate"])
        }
    if homology_cold["0_50"] != homology_cold["0_70"]:
        raise AssertionError("Locked 0.50 and 0.70 target masks differ")

    policy_scopes: dict[str, dict[str, list[dict[str, str]]]] = {}
    scope_details: dict[str, dict[str, Any]] = {}
    for policy in POLICIES:
        if policy_status[policy] != "complete":
            continue
        policy_scopes[policy], scope_details[policy] = make_scopes(
            policy,
            evaluation_rows,
            historical_transformed[policy],
            query_transformed[policy],
            homology_cold,
        )
    raw_scaffold_flags = {
        row["canonical_pair_key"]
        for row in policy_scopes["raw_primary"]["scaffold_cold"]
    }
    frozen_scaffold_flags = {
        row["canonical_pair_key"]
        for row in frozen_scaffold_rows
        if parse_bool(row["audit_scaffold_cold_under_selected_policy"])
    }
    raw_scaffold_mismatch_count = len(raw_scaffold_flags.symmetric_difference(frozen_scaffold_flags))
    if raw_scaffold_mismatch_count:
        raise AssertionError("Raw scaffold calibration failed closed")

    target_index = {target: index for index, target in enumerate(target_ids)}
    query_index = {row["query_id"]: index for index, row in enumerate(queries)}
    train_by_compound, popularity = build_train_maps(history, target_index)
    sequence_matrix = build_sequence_matrix(target_ids, sequences)
    historical_target_indices = np.asarray(
        sorted(
            {
                target_position
                for pairs in train_by_compound.values()
                for target_position, _ in pairs
            }
        ),
        dtype=np.int32,
    )
    top_columns, top_similarities = precompute_sequence_topk(
        sequence_matrix, historical_target_indices
    )
    historical_column = {
        target_position: column
        for column, target_position in enumerate(historical_target_indices)
    }
    historical_compounds_order = list(train_by_compound)
    raw_historical_smiles = {
        key: record["canonical"]
        for key, record in historical_transformed["raw_primary"].items()
    }
    raw_query_smiles = {
        key: record["canonical"]
        for key, record in query_transformed["raw_primary"].items()
    }
    raw_historical_fps = morgan_fingerprints(raw_historical_smiles, "historical")
    raw_query_fps = morgan_fingerprints(raw_query_smiles, "query")
    raw_scores, raw_ranks, allowed_matrix = build_raw_scores(
        queries,
        target_ids,
        train_by_compound,
        popularity,
        raw_historical_fps,
        raw_query_fps,
        sequence_matrix,
        historical_compounds_order,
        historical_column,
        top_columns,
        top_similarities,
    )
    rank_calibration = calibrate_rank_ledger(
        paths["frozen_ranks"],
        queries,
        target_ids,
        raw_scores,
        raw_ranks,
        allowed_matrix,
    )

    frozen_metric_map = {
        (row["scope"], row["baseline"]): float(row["Recall@50"])
        for row in frozen_metrics_rows
    }
    raw_scope_to_frozen = {
        "temporal_strict_ab": ["temporal_strict_ab"],
        "scaffold_cold": ["scaffold_cold_strict_ab"],
        "joint_scaffold_homology_0_30": ["double_cold_0_30"],
        "joint_scaffold_homology_0_50_0_70_identical": [
            "double_cold_0_50",
            "double_cold_0_70",
        ],
    }
    raw_scope_recall: dict[tuple[str, str], float] = {}
    raw_metric_cells_checked = 0
    raw_metric_mismatch_count = 0
    for baseline in BASELINES:
        for scope in SCOPES:
            observed, _, _, _ = recall_at_50(
                raw_ranks[baseline],
                policy_scopes["raw_primary"][scope],
                query_index,
                target_index,
            )
            if observed is None:
                raise AssertionError("Frozen raw scope unexpectedly empty")
            raw_scope_recall[(baseline, scope)] = observed
            for frozen_scope in raw_scope_to_frozen[scope]:
                raw_metric_cells_checked += 1
                if abs(observed - frozen_metric_map[(frozen_scope, baseline)]) > 1e-15:
                    raw_metric_mismatch_count += 1
    if raw_metric_mismatch_count:
        raise AssertionError("Raw Recall@50 calibration failed closed")

    raw_relation_keys = {
        row["canonical_pair_key"]
        for row in policy_scopes["raw_primary"]["scaffold_cold"]
    }
    raw_query_keys = {
        row["query_id"] for row in policy_scopes["raw_primary"]["scaffold_cold"]
    }
    scope_change_rows: list[dict[str, Any]] = []
    rank_change_rows: list[dict[str, Any]] = []
    recall_rows: list[dict[str, Any]] = []
    invariant_failures = 0
    policy_runtime: dict[str, float] = {}
    for policy in POLICIES:
        if policy_status[policy] != "complete":
            continue
        policy_started = time.perf_counter()
        current_relation_keys = {
            row["canonical_pair_key"]
            for row in policy_scopes[policy]["scaffold_cold"]
        }
        current_query_keys = {
            row["query_id"] for row in policy_scopes[policy]["scaffold_cold"]
        }
        detail = scope_details[policy]
        scope_change_rows.append(
            {
                "policy": policy,
                "status": "complete",
                **detail,
                "scaffold_cold_relation_count": len(current_relation_keys),
                "scaffold_cold_query_count": len(current_query_keys),
                "relation_entered_vs_raw_count": len(current_relation_keys - raw_relation_keys),
                "relation_exited_vs_raw_count": len(raw_relation_keys - current_relation_keys),
                "relation_membership_symmetric_difference_count": len(
                    current_relation_keys.symmetric_difference(raw_relation_keys)
                ),
                "query_entered_vs_raw_count": len(current_query_keys - raw_query_keys),
                "query_exited_vs_raw_count": len(raw_query_keys - current_query_keys),
                "query_membership_symmetric_difference_count": len(
                    current_query_keys.symmetric_difference(raw_query_keys)
                ),
                "joint_0_30_relation_count": len(
                    policy_scopes[policy]["joint_scaffold_homology_0_30"]
                ),
                "joint_0_30_query_count": len(
                    {
                        row["query_id"]
                        for row in policy_scopes[policy]["joint_scaffold_homology_0_30"]
                    }
                ),
                "joint_0_50_0_70_relation_count": len(
                    policy_scopes[policy][
                        "joint_scaffold_homology_0_50_0_70_identical"
                    ]
                ),
                "joint_0_50_0_70_query_count": len(
                    {
                        row["query_id"]
                        for row in policy_scopes[policy][
                            "joint_scaffold_homology_0_50_0_70_identical"
                        ]
                    }
                ),
            }
        )
        if policy == "raw_primary":
            current_scores, current_ranks = raw_scores, raw_ranks
        else:
            historical_fps = {
                key: record["fingerprint"]
                for key, record in historical_transformed[policy].items()
            }
            query_fps = {
                key: record["fingerprint"]
                for key, record in query_transformed[policy].items()
            }
            current_scores, current_ranks = score_policy(
                queries,
                target_ids,
                train_by_compound,
                historical_fps,
                query_fps,
                historical_compounds_order,
                historical_column,
                top_columns,
                top_similarities,
                allowed_matrix,
                raw_scores,
                raw_ranks,
            )
        for baseline in BASELINES:
            rank_row = rank_change_row(
                policy,
                baseline,
                current_scores[baseline],
                current_ranks[baseline],
                raw_scores[baseline],
                raw_ranks[baseline],
                allowed_matrix,
            )
            rank_change_rows.append(rank_row)
            if baseline in INVARIANT_BASELINES and (
                int(rank_row["score_changed_cell_count"]) != 0
                or int(rank_row["rank_changed_cell_count"]) != 0
            ):
                invariant_failures += 1
            for scope in SCOPES:
                rows = policy_scopes[policy][scope]
                current_recall, query_count, zero_queries, relation_hits = recall_at_50(
                    current_ranks[baseline], rows, query_index, target_index
                )
                raw_same_recall, _, _, raw_same_hits = recall_at_50(
                    raw_ranks[baseline], rows, query_index, target_index
                )
                frozen_raw = raw_scope_recall[(baseline, scope)]
                changes = endpoint_change_counts(
                    current_scores[baseline],
                    current_ranks[baseline],
                    raw_scores[baseline],
                    raw_ranks[baseline],
                    rows,
                    query_index,
                    target_index,
                )
                recall_rows.append(
                    {
                        "policy": policy,
                        "baseline": baseline,
                        "scope": scope,
                        "status": "estimable" if current_recall is not None else "not_estimable_empty_scope",
                        "relation_count": len(rows),
                        "query_count": query_count,
                        "zero_recall_at_50_query_count": zero_queries,
                        "relation_hit_at_50_count": relation_hits,
                        "recall_at_50": fmt_optional(current_recall),
                        "raw_scores_same_policy_scope_relation_hit_at_50_count": raw_same_hits,
                        "raw_scores_same_policy_scope_recall_at_50": fmt_optional(raw_same_recall),
                        "delta_ranking_only_same_membership": fmt_optional(
                            None
                            if current_recall is None or raw_same_recall is None
                            else current_recall - raw_same_recall
                        ),
                        "frozen_raw_scope_recall_at_50": f"{frozen_raw:.17g}",
                        "delta_total_vs_frozen_raw_scope": fmt_optional(
                            None if current_recall is None else current_recall - frozen_raw
                        ),
                        **changes,
                    }
                )
        policy_runtime[policy] = time.perf_counter() - policy_started
    if invariant_failures:
        raise AssertionError("Invariant baselines changed under a structure policy")

    structure_fields = list(structure_summary[0])
    write_tsv(OUTPUT_DIR / "structure_policy_summary.tsv", structure_fields, structure_summary)
    scope_fields = list(scope_change_rows[0])
    write_tsv(OUTPUT_DIR / "scaffold_scope_changes.tsv", scope_fields, scope_change_rows)
    rank_fields = list(rank_change_rows[0])
    write_tsv(OUTPUT_DIR / "rank_change_summary.tsv", rank_fields, rank_change_rows)
    recall_fields = list(recall_rows[0])
    write_tsv(OUTPUT_DIR / "scope_recall_at_50.tsv", recall_fields, recall_rows)
    blocker_fields = [
        "policy",
        "role",
        "blocker",
        "failed_record_count",
        "imputation_performed",
    ]
    write_tsv(OUTPUT_DIR / "blockers.tsv", blocker_fields, blockers)

    calibration = {
        "status": "PASS",
        "rdkit_version": rdBase.rdkitVersion,
        "raw_scaffold_endpoint_cells_checked": len(frozen_scaffold_rows),
        "raw_scaffold_membership_mismatch_count": raw_scaffold_mismatch_count,
        "raw_complete_rank_calibration": rank_calibration,
        "raw_metric_cells_checked": raw_metric_cells_checked,
        "raw_metric_mismatch_count": raw_metric_mismatch_count,
        "identity_0_50_0_70_input_hash_equal": sha256(paths["homology_0_50"])
        == sha256(paths["homology_0_70"]),
        "identity_0_50_0_70_target_mask_equal": homology_cold["0_50"]
        == homology_cold["0_70"],
        "raw_scope_counts": {
            scope: {
                "relations": len(rows),
                "queries": len({row["query_id"] for row in rows}),
            }
            for scope, rows in policy_scopes["raw_primary"].items()
        },
    }
    write_json(OUTPUT_DIR / "calibration_summary.json", calibration)

    input_hashes = {
        "schema_version": "structure_policy_input_hashes_v1",
        "parent_protocol_sha256": PARENT_PROTOCOL_SHA256,
        "child_protocol": {
            "path": "PROTOCOL.md",
            "sha256": sha256(protocol_path),
        },
        "inputs": {
            label: {
                "path": path.relative_to(PROJECT_ROOT).as_posix(),
                "sha256": sha256(path),
            }
            for label, path in paths.items()
        },
    }
    write_json(OUTPUT_DIR / "input_hashes.json", input_hashes)

    test_checks = [
        (rdBase.rdkitVersion == EXPECTED_RDKIT, "RDKit version locked at 2026.03.4"),
        (raw_scaffold_mismatch_count == 0, "raw scaffold flags reproduce all endpoint cells"),
        (rank_calibration["cells_checked"] == 3658128, "all frozen complete-rank cells checked"),
        (
            all(value == 0 for key, value in rank_calibration.items() if key != "cells_checked"),
            "raw scores and ranks reproduce every frozen cell",
        ),
        (raw_metric_cells_checked == 20, "all frozen baseline-by-scope metric cells checked"),
        (raw_metric_mismatch_count == 0, "raw Recall@50 reproduces all frozen cells"),
        (homology_cold["0_50"] == homology_cold["0_70"], "0.50 and 0.70 masks are identical"),
        (invariant_failures == 0, "popularity and sequence baselines remain invariant"),
        (
            all(status == "complete" for status in policy_status.values()),
            "all structure policies transform both roles without imputation",
        ),
        (not blockers, "no downstream policy is blocked"),
    ]
    test_results = {
        "schema_version": "structure_policy_test_results_v1",
        "status": "PASS" if all(value for value, _ in test_checks) else "FAIL",
        "check_count": len(test_checks),
        "checks": [
            {"check": label, "status": "PASS" if value else "FAIL"}
            for value, label in test_checks
        ],
    }
    write_json(OUTPUT_DIR / "test_results.json", test_results)
    if test_results["status"] != "PASS":
        raise AssertionError("Internal structure-policy tests failed")

    receipt = {
        "schema_version": "structure_policy_execution_receipt_v1",
        "status": "PASS",
        "protocol_id": "npass_structure_policy_sensitivity_v1_20260729",
        "parent_protocol_sha256": PARENT_PROTOCOL_SHA256,
        "child_protocol_sha256": sha256(protocol_path),
        "created_at_utc": utc_now(),
        "claim_boundary": (
            "Outcome-visible author-run descriptive structure sensitivity; no "
            "identity adjudication, biological validation, or independent validation."
        ),
        "role_separation": {
            "historical_and_query_maps_distinct": True,
            "historical_compound_count": len(historical_rows),
            "query_compound_count": len(query_rows),
            "pooled_structure_map_used": False,
        },
        "policy_status": policy_status,
        "policy_count": len(POLICIES),
        "baseline_count": len(BASELINES),
        "nonduplicated_scope_count": len(SCOPES),
        "aggregate_output_rows": {
            "structure_policy_summary": len(structure_summary),
            "scaffold_scope_changes": len(scope_change_rows),
            "rank_change_summary": len(rank_change_rows),
            "scope_recall_at_50": len(recall_rows),
            "blockers": len(blockers),
        },
        "calibration": calibration,
        "policy_runtime_seconds": policy_runtime,
        "total_runtime_seconds": time.perf_counter() - started,
        "rdkit_version": rdBase.rdkitVersion,
        "numpy_version": np.__version__,
        "identifier_bearing_outputs_retained": 0,
        "standardized_structure_rows_retained": 0,
    }
    write_json(OUTPUT_DIR / "execution_receipt.json", receipt)

    output_names = [
        "structure_policy_summary.tsv",
        "scaffold_scope_changes.tsv",
        "rank_change_summary.tsv",
        "scope_recall_at_50.tsv",
        "blockers.tsv",
        "calibration_summary.json",
        "input_hashes.json",
        "execution_receipt.json",
        "test_results.json",
    ]
    manifest = {
        "schema_version": "structure_policy_manifest_v1",
        "package_id": "revision_structure_policy_v1_20260729",
        "status": "BUILD_COMPLETE",
        "aggregate_only": True,
        "identifier_bearing_outputs_retained": 0,
        "parent_protocol_sha256": PARENT_PROTOCOL_SHA256,
        "child_protocol_sha256": sha256(protocol_path),
        "scripts": {
            "scripts/build_structure_policy.py": sha256(script_path),
            "scripts/validate_structure_policy.py": sha256(validator_path),
        },
        "outputs": {
            name: {"sha256": sha256(OUTPUT_DIR / name)} for name in output_names
        },
        "protocol": {"PROTOCOL.md": sha256(protocol_path)},
        "validation_report": (
            "validation_report.json is generated after this immutable manifest "
            "and is attested by its own validated_manifest_sha256 field."
        ),
    }
    write_json(OUTPUT_DIR / "manifest.json", manifest)


if __name__ == "__main__":
    main()
