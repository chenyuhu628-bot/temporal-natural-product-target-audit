"""Run aggregate-only date-precision policy sensitivity analyses.

The implementation follows Analysis C of the frozen major-revision protocol.
It reads the isolated corrective ledgers, reconstructs scenario-specific
historical evidence and chemical representatives, recomputes the four frozen
baselines, and emits no row-level identifiers.
"""

from __future__ import annotations

import calendar
import csv
import gzip
import hashlib
import json
import math
import platform
import sys
import time
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import date, datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold


ANALYSIS_ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ANALYSIS_ROOT.parent
EXECUTION = WORKSPACE / "author_run_strict_ab_asof_cutoff_execution_v1_20260728"
SUCCESSOR_SCRIPTS = WORKSPACE / "strict_ab_asof_cutoff_successor_v1_20260728" / "scripts"
WORKSPACE_SCRIPTS = WORKSPACE / "scripts"
sys.path.insert(0, str(SUCCESSOR_SCRIPTS))
sys.path.insert(0, str(WORKSPACE_SCRIPTS))

from score_asof_cutoff_successor import (  # noqa: E402
    build_sequence_matrix,
    morgan_fingerprints,
    pair_neighbor_scores,
    precompute_sequence_topk,
    sequence_transfer_scores,
    weighted_morgan_transfer_scores,
)
from pu_retrieval_metrics import query_metrics, rank_scores  # noqa: E402


ANALYSIS_ID = "revision_date_policy_v1_20260729"
PARENT_PROTOCOL_ID = "npass_strict_ab_major_revision_v4_20260729"
PARENT_PROTOCOL_SHA256 = "bd4902476160cc7c5cbacaf0cfd0f1a28c5300bde22232b3a2cc6c1f3c143dc2"
CUTOFF = date(2022, 8, 31)
TIE_SALT = "npass_strict_ab_doublecold_successor_v1_20260719"
TIER_A = "A_affinity_candidate"
TIER_B = "B_quantitative_functional_candidate"
WEIGHTS = {TIER_A: 1.0, TIER_B: 0.7}
BASELINES = (
    "weighted_target_popularity",
    "sequence_3mer_transfer",
    "weighted_morgan_transfer",
    "structure_sequence_pair_neighbor",
)
SCENARIOS = (
    "day_only_conservative",
    "interval_certain_pre_cutoff",
    "interval_earliest_bound",
)
SCOPES = (
    "temporal_strict_ab",
    "scaffold_cold_strict_ab",
    "joint_scaffold_homology_0.30",
    "joint_scaffold_homology_0.50",
    "joint_scaffold_homology_0.70",
)
EXPECTED = {
    "ledger_rows": 20647,
    "history_pairs": 4990,
    "history_compounds": 1726,
    "history_targets": 1131,
    "queries": 222,
    "endpoint_relations": 358,
    "endpoint_targets": 156,
    "candidate_targets": 4123,
}
LOCKED_PARENT_HASHES = {
    "score/corrective_prediction_ranks.tsv.gz": "87739aa818744c7084088d13c386444aa41bbef38c257083325298003181479e",
    "evaluation_inputs/evaluation_pairs.tsv.gz": "09296b066a23197a7c178f00514f2b3d9ed7e6f3c459ea92a55e01a6010d1132",
    "evaluation_inputs/scaffold_audit.tsv.gz": "fa0029ef5b7822ad5ca93f7bd93ac808f85f1e0c02e827fa91be375031b2d7af",
    "evaluation_inputs/homology_0_30.tsv.gz": "3a8247ed8f683fe6fce5fb345f56e3ec73a872b065eca922e92e494f084a1793",
    "evaluation_inputs/homology_0_50.tsv.gz": "ec3bbd435f73bc1c724efdfd94ac10c32d6f9a55bd6c9a9349233a01e5dc7f5b",
    "evaluation_inputs/homology_0_70.tsv.gz": "ec3bbd435f73bc1c724efdfd94ac10c32d6f9a55bd6c9a9349233a01e5dc7f5b",
    "restricted_ledger/historical_row_eligibility.tsv.gz": "af9cea32b35b1a2e8b437294e1dc244d2b05170c9c78aed1f71b950bf8243fa5",
    "restricted_ledger/historical_pair_before_after.tsv.gz": "a30fbefc68bc53b7a903ac3199a49acd3598e4966cfe8a51a1b669477764ab12",
    "restricted_ledger/role_separated_compound_structure_audit.tsv.gz": "58a04027c7f6d3e4cf7ddb692af61ae2a07011f6165aa76986e33b3323ad169a",
}


@dataclass
class ScenarioState:
    name: str
    selected_rows: list[dict[str, str]]
    history: list[dict[str, str]]
    structures: dict[str, str]
    row_set_sha256: str
    state_sha256: str
    repaired_compound_count: int


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def hash_lines(lines: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for line in lines:
        digest.update(line.encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def read_tsv_gz(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        require(bool(reader.fieldnames), f"Missing header: {path.name}")
        return list(reader.fieldnames or []), list(reader)


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    key = ""
    pieces: list[str] = []
    for raw in path.read_text(encoding="ascii").splitlines():
        if raw.startswith(">"):
            if key:
                records[key] = "".join(pieces)
            key = raw[1:].split("|")[0].strip()
            pieces = []
        else:
            pieces.append(raw.strip())
    if key:
        records[key] = "".join(pieces)
    require(bool(records) and all(records.values()), "Candidate FASTA is empty or malformed")
    return records


def parse_bool(value: str) -> bool:
    normalized = value.strip().casefold()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no"}:
        return False
    raise ValueError(f"Invalid boolean {value!r}")


def interval_bounds(value: str, precision: str) -> tuple[date, date] | None:
    if precision == "day":
        parsed = date.fromisoformat(value)
        return parsed, parsed
    if precision == "month":
        year, month = (int(item) for item in value.split("-"))
        return date(year, month, 1), date(year, month, calendar.monthrange(year, month)[1])
    if precision == "year":
        year = int(value)
        return date(year, 1, 1), date(year, 12, 31)
    return None


def interval_status(row: dict[str, str]) -> str:
    bounds = interval_bounds(row["publication_date"], row["date_precision"])
    if bounds is None:
        return "unresolved_interval"
    lower, upper = bounds
    if upper <= CUTOFF:
        return "definitely_before_or_on"
    if lower > CUTOFF:
        return "definitely_after"
    return "crossing_cutoff"


def eligible_under(row: dict[str, str], scenario: str) -> bool:
    if row["ref_id_type"].strip().upper() != "PMID" or not row["ref_id"].strip().isdigit():
        return False
    if not parse_bool(row["found_in_pubmed"]):
        return False
    bounds = interval_bounds(row["publication_date"], row["date_precision"])
    if bounds is None:
        return False
    lower, upper = bounds
    if scenario == "day_only_conservative":
        return row["date_precision"] == "day" and upper <= CUTOFF
    if scenario == "interval_certain_pre_cutoff":
        return upper <= CUTOFF
    if scenario == "interval_earliest_bound":
        return lower <= CUTOFF
    raise ValueError(f"Unknown scenario {scenario}")


def best_tier(rows: list[dict[str, str]]) -> str:
    tiers = {row["evidence_tier_v1_1"] for row in rows}
    require(tiers and tiers.issubset(WEIGHTS), f"Invalid strict tier set: {sorted(tiers)}")
    return TIER_A if TIER_A in tiers else TIER_B


def modal_smiles(rows: list[dict[str, str]]) -> str:
    counts = Counter(row["smiles"].strip() for row in rows if row["smiles"].strip())
    require(bool(counts), "Scenario-specific historical compound has no nonempty SMILES")
    maximum = max(counts.values())
    return sorted(value for value, count in counts.items() if count == maximum)[0]


def scaffold_key(smiles: str) -> str:
    with rdBase.BlockLogs():
        molecule = Chem.MolFromSmiles(smiles)
    require(molecule is not None, "Scaffold input is not parseable")
    scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
    return Chem.MolToSmiles(scaffold, canonical=True) if scaffold.GetNumAtoms() else ""


def verify_inputs() -> dict[str, str]:
    observed: dict[str, str] = {}
    for relative, expected_hash in LOCKED_PARENT_HASHES.items():
        path = EXECUTION / relative
        require(path.is_file(), f"Locked input is missing: {relative}")
        actual = sha256(path)
        require(actual == expected_hash, f"Locked parent hash mismatch: {relative}")
        observed[relative] = actual

    scoring_manifest_path = EXECUTION / "scoring_inputs" / "author_run_input_manifest.json"
    evaluation_manifest_path = EXECUTION / "evaluation_inputs" / "author_run_input_manifest.json"
    for manifest_path, directory in (
        (scoring_manifest_path, EXECUTION / "scoring_inputs"),
        (evaluation_manifest_path, EXECUTION / "evaluation_inputs"),
    ):
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        declared = manifest["file_sha256"]
        name_map = {
            "historical_pairs": "historical_pairs.tsv.gz",
            "scoring_queries": "scoring_queries.tsv.gz",
            "historical_compounds": "historical_compounds.tsv.gz",
            "query_compounds": "query_compounds.tsv.gz",
            "candidate_targets": "candidate_targets.tsv.gz",
            "candidate_sequences": "candidate_sequences.fasta",
            "endpoint": "evaluation_pairs.tsv.gz",
            "scaffold_audit": "scaffold_audit.tsv.gz",
            "homology_0_30": "homology_0_30.tsv.gz",
            "homology_0_50": "homology_0_50.tsv.gz",
            "homology_0_70": "homology_0_70.tsv.gz",
        }
        for key, expected_hash in declared.items():
            filename = name_map[key]
            path = directory / filename
            actual = sha256(path)
            require(actual == expected_hash, f"Isolated manifest hash mismatch: {filename}")
            observed[str(path.relative_to(EXECUTION)).replace("\\", "/")] = actual
        observed[str(manifest_path.relative_to(EXECUTION)).replace("\\", "/")] = sha256(manifest_path)

    repair_path = EXECUTION / "restricted_ledger" / "historical_v2_inchi_repair_audit.tsv.gz"
    observed[str(repair_path.relative_to(EXECUTION)).replace("\\", "/")] = sha256(repair_path)
    return dict(sorted(observed.items()))


def build_scenario(
    name: str,
    ledger: list[dict[str, str]],
    frozen_history: list[dict[str, str]],
    repair_map: dict[str, str],
) -> ScenarioState:
    selected = [row for row in ledger if eligible_under(row, name)]
    by_pair: dict[str, list[dict[str, str]]] = defaultdict(list)
    by_compound: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in selected:
        by_pair[row["canonical_pair_key"]].append(row)
        by_compound[row["inchikey_full"]].append(row)
    frozen_keys = {row["canonical_pair_key"] for row in frozen_history}
    require(set(by_pair) == frozen_keys, f"Scenario {name} changes the 4,990-pair keyset")

    history: list[dict[str, str]] = []
    frozen_by_key = {row["canonical_pair_key"]: row for row in frozen_history}
    for key in sorted(by_pair):
        frozen = frozen_by_key[key]
        history.append(
            {
                "canonical_pair_key": key,
                "inchikey_full": frozen["inchikey_full"],
                "uniprot_canonical_accession": frozen["uniprot_canonical_accession"],
                "best_strict_evidence_tier": best_tier(by_pair[key]),
                "selected_row_count": str(len(by_pair[key])),
            }
        )

    structures: dict[str, str] = {}
    repaired = 0
    for compound in sorted(by_compound):
        chosen = modal_smiles(by_compound[compound])
        with rdBase.BlockLogs():
            molecule = Chem.MolFromSmiles(chosen)
        if molecule is None:
            require(compound in repair_map, "Unparseable modal SMILES lacks a locked validated repair")
            chosen = repair_map[compound]
            repaired += 1
            with rdBase.BlockLogs():
                molecule = Chem.MolFromSmiles(chosen)
            require(molecule is not None, "Locked validated repair is not parseable")
        structures[compound] = chosen

    require(len(history) == EXPECTED["history_pairs"], "Historical-pair count changed")
    require(len(structures) == EXPECTED["history_compounds"], "Historical-compound count changed")
    require(
        len({row["uniprot_canonical_accession"] for row in history}) == EXPECTED["history_targets"],
        "Historical-target count changed",
    )
    row_set_hash = hash_lines(
        sorted(f"{row['canonical_pair_key']}|{row['source_row_number']}" for row in selected)
    )
    state_hash = hash_lines(
        [
            *(f"PAIR|{row['canonical_pair_key']}|{row['best_strict_evidence_tier']}" for row in history),
            *(f"STRUCTURE|{compound}|{structures[compound]}" for compound in sorted(structures)),
        ]
    )
    return ScenarioState(name, selected, history, structures, row_set_hash, state_hash, repaired)


def build_train_maps(
    history: list[dict[str, str]], target_index: dict[str, int]
) -> tuple[dict[str, list[tuple[int, float]]], np.ndarray]:
    by_compound: dict[str, list[tuple[int, float]]] = defaultdict(list)
    popularity = np.zeros(len(target_index), dtype=np.float32)
    for row in history:
        target_idx = target_index[row["uniprot_canonical_accession"]]
        weight = WEIGHTS[row["best_strict_evidence_tier"]]
        by_compound[row["inchikey_full"]].append((target_idx, weight))
        popularity[target_idx] += weight
    return dict(by_compound), popularity


def load_homology(endpoint_targets: set[str]) -> tuple[dict[str, dict[str, bool]], dict[str, str]]:
    output: dict[str, dict[str, bool]] = {}
    source_hashes: dict[str, str] = {}
    for threshold in ("0_30", "0_50", "0_70"):
        path = EXECUTION / "evaluation_inputs" / f"homology_{threshold}.tsv.gz"
        _, rows = read_tsv_gz(path)
        flags = {
            row["uniprot_canonical_accession"]: parse_bool(row["is_future_target_homology_cold_candidate"])
            for row in rows
        }
        require(set(flags) == endpoint_targets, f"Homology {threshold} target set mismatch")
        output[threshold] = flags
        source_hashes[threshold] = sha256(path)
    require(output["0_50"] == output["0_70"], "0.50 and 0.70 homology masks are not identical")
    return output, source_hashes


def mask_hash(pair_keys: Iterable[str]) -> str:
    return hash_lines(sorted(pair_keys))


def build_scopes(
    state: ScenarioState,
    endpoint: list[dict[str, str]],
    query_structures: dict[str, str],
    homology: dict[str, dict[str, bool]],
) -> tuple[dict[str, list[dict[str, str]]], dict[str, str]]:
    historical_scaffolds = {scaffold_key(smiles) for smiles in state.structures.values()}
    historical_scaffolds.discard("")
    query_scaffolds = {compound: scaffold_key(smiles) for compound, smiles in query_structures.items()}
    scaffold_cold: dict[str, bool] = {}
    for row in endpoint:
        key = query_scaffolds[row["inchikey_full"]]
        scaffold_cold[row["canonical_pair_key"]] = bool(key) and key not in historical_scaffolds
    scopes = {
        "temporal_strict_ab": list(endpoint),
        "scaffold_cold_strict_ab": [row for row in endpoint if scaffold_cold[row["canonical_pair_key"]]],
    }
    for threshold in ("0_30", "0_50", "0_70"):
        scopes[f"joint_scaffold_homology_{threshold.replace('_', '.')}"] = [
            row
            for row in endpoint
            if scaffold_cold[row["canonical_pair_key"]]
            and homology[threshold][row["uniprot_canonical_accession"]]
        ]
    hashes = {scope: mask_hash(row["canonical_pair_key"] for row in rows) for scope, rows in scopes.items()}
    return scopes, hashes


def load_frozen_endpoint_ranks(
    endpoint: list[dict[str, str]],
) -> dict[tuple[str, str, str], tuple[int, float]]:
    relevant = defaultdict(set)
    for row in endpoint:
        relevant[row["query_id"]].add(row["uniprot_canonical_accession"])
    output: dict[tuple[str, str, str], tuple[int, float]] = {}
    path = EXECUTION / "score" / "corrective_prediction_ranks.tsv.gz"
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            query_id = row["query_id"]
            target = row["target_uniprot_accession"]
            if target in relevant.get(query_id, set()):
                output[(row["baseline"], query_id, target)] = (int(row["rank"]), float(row["score"]))
    require(len(output) == len(BASELINES) * len(endpoint), "Frozen endpoint rank extraction is incomplete")
    return output


def score_unique_states(
    states: dict[str, ScenarioState],
    queries: list[dict[str, str]],
    query_structures: dict[str, str],
    target_ids: list[str],
    sequences: dict[str, str],
    endpoint: list[dict[str, str]],
    frozen_endpoint_ranks: dict[tuple[str, str, str], tuple[int, float]],
) -> tuple[
    dict[str, dict[str, dict[str, dict[str, int]]]],
    dict[tuple[str, str], dict[str, Any]],
    dict[str, int],
]:
    target_index = {target: idx for idx, target in enumerate(target_ids)}
    endpoint_targets_by_query: dict[str, set[str]] = defaultdict(set)
    for row in endpoint:
        endpoint_targets_by_query[row["query_id"]].add(row["uniprot_canonical_accession"])

    sequence_matrix = build_sequence_matrix(target_ids, sequences)
    query_fps = morgan_fingerprints(query_structures, "query")
    prepared: dict[str, dict[str, Any]] = {}
    for key, state in states.items():
        train, popularity = build_train_maps(state.history, target_index)
        history_fps = morgan_fingerprints(state.structures, "historical")
        historical_target_indices = np.asarray(
            sorted({target for values in train.values() for target, _ in values}), dtype=np.int32
        )
        top_columns, top_similarities = precompute_sequence_topk(sequence_matrix, historical_target_indices)
        prepared[key] = {
            "train": train,
            "popularity": popularity,
            "history_fps": history_fps,
            "historical_compounds": list(train),
            "historical_column": {target: column for column, target in enumerate(historical_target_indices)},
            "top_columns": top_columns,
            "top_similarities": top_similarities,
        }

    primary_matches = [key for key, state in states.items() if state.name == "day_only_conservative"]
    require(len(primary_matches) == 1, "Unique-state map lacks exactly one day-only primary state")
    primary_key = primary_matches[0]
    endpoint_ranks: dict[str, dict[str, dict[str, dict[str, int]]]] = {
        key: {baseline: defaultdict(dict) for baseline in BASELINES} for key in prepared
    }
    changes: dict[tuple[str, str], dict[str, Any]] = {}
    for key in prepared:
        for baseline in BASELINES:
            changes[(key, baseline)] = {
                "eligible_candidate_count": 0,
                "score_changed_candidate_count": 0,
                "rank_changed_candidate_count": 0,
                "top50_membership_changed_query_count": 0,
                "endpoint_rank_changed_relation_count": 0,
                "endpoint_absolute_rank_changes": [],
                "maximum_absolute_candidate_rank_change": 0,
            }

    frozen_rank_matches = 0
    query_candidate_counts: dict[str, int] = {}
    for query in queries:
        query_id = query["query_id"]
        compound = query["inchikey_full"]
        scores_by_state: dict[str, dict[str, np.ndarray]] = defaultdict(dict)
        ranks_by_state: dict[str, dict[str, np.ndarray]] = defaultdict(dict)
        allowed: np.ndarray | None = None
        for key, item in prepared.items():
            train = item["train"]
            current_allowed = np.ones(len(target_ids), dtype=bool)
            for target, _ in train.get(compound, []):
                current_allowed[target] = False
            if allowed is None:
                allowed = current_allowed
            else:
                require(np.array_equal(allowed, current_allowed), "Scenario changed a query candidate mask")
            query_fp = query_fps[compound]
            scorers = {
                "weighted_target_popularity": lambda item=item: item["popularity"].copy(),
                "sequence_3mer_transfer": lambda item=item: sequence_transfer_scores(
                    compound, item["train"], sequence_matrix, len(target_ids)
                ),
                "weighted_morgan_transfer": lambda item=item: weighted_morgan_transfer_scores(
                    query_fp,
                    item["historical_compounds"],
                    item["train"],
                    item["history_fps"],
                    len(target_ids),
                ),
                "structure_sequence_pair_neighbor": lambda item=item: pair_neighbor_scores(
                    query_fp,
                    item["historical_compounds"],
                    item["train"],
                    item["history_fps"],
                    item["historical_column"],
                    item["top_columns"],
                    item["top_similarities"],
                ),
            }
            for baseline, scorer in scorers.items():
                scores = scorer()
                require(np.all(np.isfinite(scores)), "Non-finite scenario score")
                _, ranks = rank_scores(scores, current_allowed, query_id, target_ids, TIE_SALT)
                candidate_ranks = ranks[current_allowed]
                require(
                    candidate_ranks.min() == 1
                    and candidate_ranks.max() == int(current_allowed.sum())
                    and np.unique(candidate_ranks).size == int(current_allowed.sum()),
                    "Candidate ranks are not a permutation",
                )
                scores_by_state[key][baseline] = scores
                ranks_by_state[key][baseline] = ranks
                for target in endpoint_targets_by_query[query_id]:
                    endpoint_ranks[key][baseline][query_id][target] = int(ranks[target_index[target]])
        require(allowed is not None, "Missing query candidate mask")
        query_candidate_counts[query_id] = int(allowed.sum())

        for baseline in BASELINES:
            primary_scores = scores_by_state[primary_key][baseline]
            primary_ranks = ranks_by_state[primary_key][baseline]
            for target in endpoint_targets_by_query[query_id]:
                idx = target_index[target]
                frozen_rank, frozen_score = frozen_endpoint_ranks[(baseline, query_id, target)]
                require(int(primary_ranks[idx]) == frozen_rank, "Primary endpoint rank does not reproduce frozen score ledger")
                require(float(primary_scores[idx]) == frozen_score, "Primary endpoint score does not reproduce frozen score ledger")
                frozen_rank_matches += 1
            for key in prepared:
                item = changes[(key, baseline)]
                scores = scores_by_state[key][baseline]
                ranks = ranks_by_state[key][baseline]
                item["eligible_candidate_count"] += int(allowed.sum())
                item["score_changed_candidate_count"] += int(np.count_nonzero(scores[allowed] != primary_scores[allowed]))
                rank_diff = np.abs(ranks[allowed].astype(np.int64) - primary_ranks[allowed].astype(np.int64))
                item["rank_changed_candidate_count"] += int(np.count_nonzero(rank_diff))
                item["maximum_absolute_candidate_rank_change"] = max(
                    item["maximum_absolute_candidate_rank_change"], int(rank_diff.max(initial=0))
                )
                primary_top50 = set(np.flatnonzero(allowed & (primary_ranks <= 50)))
                current_top50 = set(np.flatnonzero(allowed & (ranks <= 50)))
                item["top50_membership_changed_query_count"] += primary_top50 != current_top50
                for target in endpoint_targets_by_query[query_id]:
                    idx = target_index[target]
                    delta = abs(int(ranks[idx]) - int(primary_ranks[idx]))
                    item["endpoint_absolute_rank_changes"].append(delta)
                    item["endpoint_rank_changed_relation_count"] += delta != 0

    require(frozen_rank_matches == len(BASELINES) * len(endpoint), "Frozen endpoint verification count mismatch")
    return endpoint_ranks, changes, query_candidate_counts


def scope_summary_rows(
    scenario_to_state: dict[str, ScenarioState],
    scopes_by_state: dict[str, dict[str, list[dict[str, str]]]],
    scope_hashes_by_state: dict[str, dict[str, str]],
    homology_source_hashes: dict[str, str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for scenario in SCENARIOS:
        state_key = scenario_to_state[scenario].state_sha256
        for scope in SCOPES:
            members = scopes_by_state[state_key][scope]
            provenance = ""
            if scope.endswith("0.30"):
                provenance = homology_source_hashes["0_30"]
            elif scope.endswith("0.50"):
                provenance = homology_source_hashes["0_50"]
            elif scope.endswith("0.70"):
                provenance = homology_source_hashes["0_70"]
            display_scope = (
                "joint_scaffold_homology_0.50/0.70"
                if scope in {"joint_scaffold_homology_0.50", "joint_scaffold_homology_0.70"}
                else scope
            )
            rows.append(
                {
                    "scenario": scenario,
                    "provenance_scope": scope,
                    "display_scope": display_scope,
                    "candidate_relation_count": len(members),
                    "query_count": len({row["query_id"] for row in members}),
                    "target_count": len({row["uniprot_canonical_accession"] for row in members}),
                    "membership_sha256": scope_hashes_by_state[state_key][scope],
                    "homology_source_sha256": provenance,
                }
            )
    return rows


def evaluate_recall_rows(
    scenario_to_state: dict[str, ScenarioState],
    scopes_by_state: dict[str, dict[str, list[dict[str, str]]]],
    endpoint_ranks: dict[str, dict[str, dict[str, dict[str, int]]]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for scenario in SCENARIOS:
        state_key = scenario_to_state[scenario].state_sha256
        for scope in SCOPES:
            grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
            for relation in scopes_by_state[state_key][scope]:
                grouped[relation["query_id"]].append(relation)
            for baseline in BASELINES:
                per_query: list[float] = []
                hit_queries = 0
                for query_id in sorted(grouped):
                    ranks = [
                        endpoint_ranks[state_key][baseline][query_id][row["uniprot_canonical_accession"]]
                        for row in grouped[query_id]
                    ]
                    recall = query_metrics(ranks, (50,))["Recall@50"]
                    per_query.append(recall)
                    hit_queries += recall > 0.0
                require(bool(per_query), f"Empty evaluation scope: {scope}")
                rows.append(
                    {
                        "analysis_label": "author_run_outcome_visible_post_hoc_descriptive_sensitivity",
                        "scenario": scenario,
                        "provenance_scope": scope,
                        "display_scope": (
                            "joint_scaffold_homology_0.50/0.70"
                            if scope in {"joint_scaffold_homology_0.50", "joint_scaffold_homology_0.70"}
                            else scope
                        ),
                        "baseline": baseline,
                        "evaluable_query_count": len(per_query),
                        "queries_with_at_least_one_top50_hit": hit_queries,
                        "Recall@50": f"{float(np.mean(per_query)):.17g}",
                    }
                )
    return rows


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    require(bool(rows), f"Refusing to write empty output: {path.name}")
    require(not path.exists(), f"Refusing to overwrite output: {path.name}")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    require(not path.exists(), f"Refusing to overwrite output: {path.name}")
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def build_summary(
    interval_rows: list[dict[str, object]],
    history_rows: list[dict[str, object]],
    structure_rows: list[dict[str, object]],
    recall_rows: list[dict[str, object]],
) -> str:
    lines = [
        "# Date-precision policy sensitivity",
        "",
        "Author-run, outcome-visible, post hoc descriptive sensitivity; not independent validation.",
        "",
        "## Interval-status audit",
        "",
        "| Status | Rows | Historical pairs |",
        "|---|---:|---:|",
    ]
    for row in interval_rows:
        lines.append(f"| {row['interval_status']} | {row['source_row_count']} | {row['historical_pair_count']} |")
    lines += [
        "",
        "All date-resolved month/year rows were definitely before the cutoff; no cutoff-crossing interval was observed.",
        "",
        "## Scenario reconstruction",
        "",
        "| Scenario | Selected rows | Pairs | Tier changes vs day-only | A pairs | B pairs |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in history_rows:
        lines.append(
            f"| {row['scenario']} | {row['selected_source_row_count']} | {row['historical_pair_count']} | "
            f"{row['tier_or_weight_changed_pair_count_vs_day_only']} | {row['tier_A_pair_count']} | {row['tier_B_pair_count']} |"
        )
    lines += [
        "",
        "| Scenario | SMILES changes | Morgan changes | Scaffold changes |",
        "|---|---:|---:|---:|",
    ]
    for row in structure_rows:
        lines.append(
            f"| {row['scenario']} | {row['representative_smiles_changed_count_vs_day_only']} | "
            f"{row['morgan_fingerprint_changed_count_vs_day_only']} | {row['scaffold_assignment_changed_count_vs_day_only']} |"
        )
    lines += [
        "",
        "## Recall@50",
        "",
        "The 0.50 and 0.70 masks are identical and are displayed once below; both provenance rows remain in the TSV.",
        "",
        "| Scenario | Scope | Baseline | Recall@50 |",
        "|---|---|---|---:|",
    ]
    seen: set[tuple[str, str, str]] = set()
    for row in recall_rows:
        key = (str(row["scenario"]), str(row["display_scope"]), str(row["baseline"]))
        if key in seen:
            continue
        seen.add(key)
        lines.append(f"| {key[0]} | {key[1]} | {key[2]} | {float(row['Recall@50']):.6f} |")
    lines += [
        "",
        "`interval_certain_pre_cutoff` and `interval_earliest_bound` are retained as separate policy labels but collapse to the same selected row set in this ledger because crossing and definitely-after intervals are absent.",
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    started = time.perf_counter()
    output_dir = ANALYSIS_ROOT / "outputs"
    require(output_dir.is_dir(), "Output directory is missing")
    require(not any(output_dir.iterdir()), "Output directory must be empty for a create-once run")

    source_hashes = verify_inputs()
    _, ledger = read_tsv_gz(EXECUTION / "restricted_ledger" / "historical_row_eligibility.tsv.gz")
    _, frozen_history = read_tsv_gz(EXECUTION / "scoring_inputs" / "historical_pairs.tsv.gz")
    _, frozen_history_structures_rows = read_tsv_gz(
        EXECUTION / "scoring_inputs" / "historical_compounds.tsv.gz"
    )
    _, query_rows = read_tsv_gz(EXECUTION / "scoring_inputs" / "scoring_queries.tsv.gz")
    _, query_structure_rows = read_tsv_gz(EXECUTION / "scoring_inputs" / "query_compounds.tsv.gz")
    _, target_rows = read_tsv_gz(EXECUTION / "scoring_inputs" / "candidate_targets.tsv.gz")
    sequences = read_fasta(EXECUTION / "scoring_inputs" / "candidate_sequences.fasta")
    _, endpoint = read_tsv_gz(EXECUTION / "evaluation_inputs" / "evaluation_pairs.tsv.gz")
    _, repair_rows = read_tsv_gz(
        EXECUTION / "restricted_ledger" / "historical_v2_inchi_repair_audit.tsv.gz"
    )

    require(len(ledger) == EXPECTED["ledger_rows"], "Ledger row count mismatch")
    require(len(frozen_history) == EXPECTED["history_pairs"], "Frozen history count mismatch")
    require(len(query_rows) == EXPECTED["queries"], "Query count mismatch")
    require(len(endpoint) == EXPECTED["endpoint_relations"], "Endpoint relation count mismatch")
    require(len(target_rows) == EXPECTED["candidate_targets"], "Candidate target count mismatch")
    target_ids = [row["uniprot_canonical_accession"] for row in target_rows]
    require(target_ids == sorted(target_ids), "Candidate targets are not deterministically sorted")
    require(set(target_ids) == set(sequences), "Candidate target and FASTA sets differ")

    frozen_structures = {row["inchikey_full"]: row["representative_smiles"] for row in frozen_history_structures_rows}
    query_structures = {row["inchikey_full"]: row["representative_smiles"] for row in query_structure_rows}
    repair_map = {row["inchikey_full"]: row["repaired_smiles"] for row in repair_rows}

    interval_rows: list[dict[str, object]] = []
    for status in ("definitely_before_or_on", "crossing_cutoff", "definitely_after", "unresolved_interval"):
        members = [row for row in ledger if interval_status(row) == status]
        interval_rows.append(
            {
                "interval_status": status,
                "source_row_count": len(members),
                "historical_pair_count": len({row["canonical_pair_key"] for row in members}),
                "day_precision_count": sum(row["date_precision"] == "day" for row in members),
                "month_precision_count": sum(row["date_precision"] == "month" for row in members),
                "year_precision_count": sum(row["date_precision"] == "year" for row in members),
                "unresolved_precision_count": sum(row["date_precision"] not in {"day", "month", "year"} for row in members),
            }
        )

    scenario_to_state = {
        scenario: build_scenario(scenario, ledger, frozen_history, repair_map) for scenario in SCENARIOS
    }
    primary = scenario_to_state["day_only_conservative"]
    require(primary.structures == frozen_structures, "Day-only reconstructed structures do not match frozen history")
    require(
        [(row["canonical_pair_key"], row["best_strict_evidence_tier"]) for row in primary.history]
        == [(row["canonical_pair_key"], row["best_strict_evidence_tier"]) for row in sorted(frozen_history, key=lambda item: item["canonical_pair_key"])],
        "Day-only reconstructed tiers do not match frozen history",
    )

    unique_states = {state.state_sha256: state for state in scenario_to_state.values()}
    endpoint_targets = {row["uniprot_canonical_accession"] for row in endpoint}
    require(len(endpoint_targets) == EXPECTED["endpoint_targets"], "Endpoint target count mismatch")
    homology, homology_source_hashes = load_homology(endpoint_targets)
    scopes_by_state: dict[str, dict[str, list[dict[str, str]]]] = {}
    scope_hashes_by_state: dict[str, dict[str, str]] = {}
    for state_key, state in unique_states.items():
        scopes_by_state[state_key], scope_hashes_by_state[state_key] = build_scopes(
            state, endpoint, query_structures, homology
        )

    _, frozen_scaffold_rows = read_tsv_gz(EXECUTION / "evaluation_inputs" / "scaffold_audit.tsv.gz")
    frozen_scaffold_keys = {
        row["canonical_pair_key"]
        for row in frozen_scaffold_rows
        if parse_bool(row["audit_scaffold_cold_under_selected_policy"])
    }
    require(
        scope_hashes_by_state[primary.state_sha256]["scaffold_cold_strict_ab"]
        == mask_hash(frozen_scaffold_keys),
        "Primary recomputed scaffold mask does not match frozen scaffold audit",
    )
    for state_key in unique_states:
        require(
            scope_hashes_by_state[state_key]["joint_scaffold_homology_0.50"]
            == scope_hashes_by_state[state_key]["joint_scaffold_homology_0.70"],
            "0.50 and 0.70 joint masks differ",
        )

    frozen_endpoint_ranks = load_frozen_endpoint_ranks(endpoint)
    endpoint_ranks, score_changes, query_candidate_counts = score_unique_states(
        unique_states,
        query_rows,
        query_structures,
        target_ids,
        sequences,
        endpoint,
        frozen_endpoint_ranks,
    )

    primary_tiers = {row["canonical_pair_key"]: row["best_strict_evidence_tier"] for row in primary.history}
    primary_selected_counts = {
        row["canonical_pair_key"]: row["selected_row_count"] for row in primary.history
    }
    primary_fps = morgan_fingerprints(primary.structures, "historical")
    primary_scaffolds = {compound: scaffold_key(smiles) for compound, smiles in primary.structures.items()}
    history_rows: list[dict[str, object]] = []
    structure_rows: list[dict[str, object]] = []
    for scenario in SCENARIOS:
        state = scenario_to_state[scenario]
        tiers = {row["canonical_pair_key"]: row["best_strict_evidence_tier"] for row in state.history}
        history_rows.append(
            {
                "scenario": scenario,
                "selected_source_row_count": len(state.selected_rows),
                "selected_row_set_sha256": state.row_set_sha256,
                "state_sha256": state.state_sha256,
                "historical_pair_count": len(state.history),
                "historical_compound_count": len(state.structures),
                "historical_target_count": len({row["uniprot_canonical_accession"] for row in state.history}),
                "historical_membership_change_count_vs_day_only": len(set(tiers).symmetric_difference(primary_tiers)),
                "tier_or_weight_changed_pair_count_vs_day_only": sum(tiers[key] != primary_tiers[key] for key in tiers),
                "tier_A_pair_count": sum(value == TIER_A for value in tiers.values()),
                "tier_B_pair_count": sum(value == TIER_B for value in tiers.values()),
                "pairs_with_selected_row_count_change_vs_day_only": sum(
                    row["selected_row_count"] != primary_selected_counts[row["canonical_pair_key"]]
                    for row in state.history
                ),
                "validated_repair_compound_count": state.repaired_compound_count,
            }
        )
        fps = morgan_fingerprints(state.structures, "historical")
        scaffolds = {compound: scaffold_key(smiles) for compound, smiles in state.structures.items()}
        structure_rows.append(
            {
                "scenario": scenario,
                "representative_smiles_changed_count_vs_day_only": sum(
                    state.structures[key] != primary.structures[key] for key in state.structures
                ),
                "morgan_fingerprint_changed_count_vs_day_only": sum(
                    DataStructs.TanimotoSimilarity(fps[key], primary_fps[key]) < 1.0 for key in fps
                ),
                "scaffold_assignment_changed_count_vs_day_only": sum(
                    scaffolds[key] != primary_scaffolds[key] for key in scaffolds
                ),
                "valid_structure_count": len(state.structures),
                "invalid_structure_count": 0,
            }
        )

    scope_rows = scope_summary_rows(
        scenario_to_state, scopes_by_state, scope_hashes_by_state, homology_source_hashes
    )
    recall_rows = evaluate_recall_rows(scenario_to_state, scopes_by_state, endpoint_ranks)
    change_rows: list[dict[str, object]] = []
    for scenario in SCENARIOS:
        state_key = scenario_to_state[scenario].state_sha256
        for baseline in BASELINES:
            item = score_changes[(state_key, baseline)]
            abs_endpoint = np.asarray(item["endpoint_absolute_rank_changes"], dtype=float)
            change_rows.append(
                {
                    "scenario": scenario,
                    "baseline": baseline,
                    "eligible_candidate_count": item["eligible_candidate_count"],
                    "score_changed_candidate_count_vs_day_only": item["score_changed_candidate_count"],
                    "score_changed_candidate_fraction_vs_day_only": f"{item['score_changed_candidate_count'] / item['eligible_candidate_count']:.17g}",
                    "rank_changed_candidate_count_vs_day_only": item["rank_changed_candidate_count"],
                    "rank_changed_candidate_fraction_vs_day_only": f"{item['rank_changed_candidate_count'] / item['eligible_candidate_count']:.17g}",
                    "top50_membership_changed_query_count_vs_day_only": item["top50_membership_changed_query_count"],
                    "endpoint_rank_changed_relation_count_vs_day_only": item["endpoint_rank_changed_relation_count"],
                    "endpoint_absolute_rank_change_median": f"{float(np.median(abs_endpoint)):.17g}",
                    "endpoint_absolute_rank_change_mean": f"{float(np.mean(abs_endpoint)):.17g}",
                    "endpoint_absolute_rank_change_max": int(abs_endpoint.max(initial=0)),
                    "candidate_absolute_rank_change_max": item["maximum_absolute_candidate_rank_change"],
                }
            )

    equivalence_rows = [
        {
            "scenario": scenario,
            "selected_row_set_sha256": scenario_to_state[scenario].row_set_sha256,
            "scoring_state_sha256": scenario_to_state[scenario].state_sha256,
            "equivalent_to_day_only": scenario_to_state[scenario].state_sha256 == primary.state_sha256,
            "equivalent_to_interval_certain": scenario_to_state[scenario].state_sha256
            == scenario_to_state["interval_certain_pre_cutoff"].state_sha256,
            "unique_state_scored_once_and_reused_for_equivalent_policy_label": sum(
                item.state_sha256 == scenario_to_state[scenario].state_sha256
                for item in scenario_to_state.values()
            )
            > 1,
        }
        for scenario in SCENARIOS
    ]

    outputs = {
        "interval_status_counts.tsv": interval_rows,
        "scenario_history_summary.tsv": history_rows,
        "scenario_structure_summary.tsv": structure_rows,
        "scope_denominators.tsv": scope_rows,
        "recall_at_50.tsv": recall_rows,
        "score_rank_change_summary.tsv": change_rows,
        "scenario_equivalence.tsv": equivalence_rows,
    }
    for filename, rows in outputs.items():
        write_tsv(output_dir / filename, rows)
    summary_path = output_dir / "SUMMARY.md"
    summary_path.write_text(
        build_summary(interval_rows, history_rows, structure_rows, recall_rows), encoding="utf-8"
    )

    output_hashes = {
        path.name: sha256(path)
        for path in sorted(output_dir.iterdir())
        if path.is_file()
    }
    runtime_seconds = time.perf_counter() - started
    receipt = {
        "analysis_id": ANALYSIS_ID,
        "parent_protocol_id": PARENT_PROTOCOL_ID,
        "parent_protocol_sha256": PARENT_PROTOCOL_SHA256,
        "status": "PASS",
        "execution_mode": "author_run_outcome_visible_post_hoc_descriptive_sensitivity",
        "claim_boundary": "No independent validation, blinded scoring, policy superiority, external validation, or biological confirmation claim.",
        "scenario_count": len(SCENARIOS),
        "unique_scoring_state_count": len(unique_states),
        "scenario_state_equivalence": {
            scenario: scenario_to_state[scenario].state_sha256 for scenario in SCENARIOS
        },
        "frozen_endpoint_score_rank_cells_reproduced": len(BASELINES) * len(endpoint),
        "complete_candidate_rank_blocks_checked": len(unique_states) * len(BASELINES) * len(query_rows),
        "absolute_paths_emitted": False,
        "identifier_bearing_rows_emitted": False,
        "main_manuscript_modified": False,
        "deviations_from_frozen_analysis_c": [],
        "runtime_seconds": runtime_seconds,
        "environment": {
            "python": platform.python_version(),
            "python_implementation": platform.python_implementation(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "rdkit": rdBase.rdkitVersion,
        },
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
    }
    write_json(output_dir / "EXECUTION_RECEIPT.json", receipt)
    output_hashes["EXECUTION_RECEIPT.json"] = sha256(output_dir / "EXECUTION_RECEIPT.json")

    code_paths = {
        "runner": Path(__file__),
        "protocol": ANALYSIS_ROOT / "PROTOCOL.md",
        "implementation_lock": ANALYSIS_ROOT / "IMPLEMENTATION_LOCK.json",
        "code_lock": ANALYSIS_ROOT / "CODE_LOCK.json",
        "score_reference": SUCCESSOR_SCRIPTS / "score_asof_cutoff_successor.py",
        "rank_metrics_reference": WORKSPACE_SCRIPTS / "pu_retrieval_metrics.py",
    }
    manifest = {
        "schema_version": "1.0",
        "analysis_id": ANALYSIS_ID,
        "parent_protocol_id": PARENT_PROTOCOL_ID,
        "parent_protocol_sha256": PARENT_PROTOCOL_SHA256,
        "claim_boundary": receipt["claim_boundary"],
        "source_files": source_hashes,
        "code_files": {
            str(path.relative_to(WORKSPACE)).replace("\\", "/"): sha256(path)
            for path in code_paths.values()
        },
        "outputs_before_manifest": dict(sorted(output_hashes.items())),
        "row_counts": {filename: len(rows) for filename, rows in outputs.items()},
        "scenario_names": list(SCENARIOS),
        "unique_scoring_state_count": len(unique_states),
        "homology_mask_policy": {
            "0.50_and_0.70_identical": True,
            "0.30_source_sha256": homology_source_hashes["0_30"],
            "0.50_source_sha256": homology_source_hashes["0_50"],
            "0.70_source_sha256": homology_source_hashes["0_70"],
        },
        "aggregate_only": True,
    }
    write_json(output_dir / "MANIFEST.json", manifest)
    print(
        json.dumps(
            {
                "status": "PASS",
                "analysis_id": ANALYSIS_ID,
                "scenarios": len(SCENARIOS),
                "unique_states": len(unique_states),
                "runtime_seconds": runtime_seconds,
            },
            ensure_ascii=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
