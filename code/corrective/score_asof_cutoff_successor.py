"""Score the four frozen baselines with role-separated corrective structures.

The command deliberately has no endpoint or cold-scope argument. Historical
and query structures are separate inputs, including when a full InChIKey is
present in both roles.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import sys
import time
import tracemalloc
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

import numpy as np
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator
from sklearn.feature_extraction.text import TfidfVectorizer

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
sys.path.insert(0, str(WORKSPACE / "scripts"))

from pu_retrieval_metrics import rank_scores

from asof_successor_common import (
    BASELINES,
    EXPECTED_CANDIDATE_TARGETS,
    EXPECTED_COMPLETE_RANK_ROWS,
    LEGACY_TIE_SALT,
    MORGAN_BITS,
    MORGAN_RADIUS,
    PAIR_NEIGHBOR_TOP_K,
    PROTOCOL_ID,
    RUN_MODE,
    STRICT_TIERS,
    assert_isolated_input,
    assert_new_output_dir,
    code_hashes,
    environment_receipt,
    finite_scores,
    load_input_manifest,
    load_receipt,
    peak_rss_bytes,
    read_fasta,
    read_tsv_gz,
    require,
    sha256,
    validate_historical_pairs,
    validate_queries,
    validate_role_structures,
    validate_targets,
    write_json,
)


SCORING_INPUT_KIND = "corrective_role_separated_scoring_without_endpoint_file"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    add_arguments(parser)
    return parser


def add_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--receipt", required=True, type=Path)
    parser.add_argument("--input-manifest", required=True, type=Path)
    parser.add_argument("--historical-pairs", required=True, type=Path)
    parser.add_argument("--scoring-queries", required=True, type=Path)
    parser.add_argument("--historical-compounds", required=True, type=Path)
    parser.add_argument("--query-compounds", required=True, type=Path)
    parser.add_argument("--candidate-targets", required=True, type=Path)
    parser.add_argument("--candidate-sequences", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)


def morgan_fingerprints(smiles_by_compound: dict[str, str], role: str) -> dict[str, object]:
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=MORGAN_RADIUS, fpSize=MORGAN_BITS)
    fingerprints: dict[str, object] = {}
    for compound, smiles in smiles_by_compound.items():
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(f"RDKit cannot parse locked {role} structure: {compound}")
        fingerprints[compound] = generator.GetFingerprint(molecule)
    return fingerprints


def build_train_maps(
    history: list[dict[str, str]], target_index: dict[str, int]
) -> tuple[dict[str, list[tuple[int, float]]], np.ndarray]:
    by_compound: dict[str, list[tuple[int, float]]] = defaultdict(list)
    popularity = np.zeros(len(target_index), dtype=np.float32)
    for row in history:
        target = row["uniprot_canonical_accession"]
        if target not in target_index:
            raise ValueError(f"Historical target absent from candidate universe: {target}")
        target_idx = target_index[target]
        weight = STRICT_TIERS[row["best_strict_evidence_tier"]]
        by_compound[row["inchikey_full"]].append((target_idx, weight))
        popularity[target_idx] += weight
    return dict(by_compound), popularity


def build_sequence_matrix(target_ids: list[str], sequences: dict[str, str]) -> object:
    vectorizer = TfidfVectorizer(
        analyzer="char",
        ngram_range=(3, 3),
        lowercase=False,
        norm="l2",
        dtype=np.float32,
    )
    return vectorizer.fit_transform([sequences[target] for target in target_ids])


def sequence_transfer_scores(
    query_compound: str,
    train_by_compound: dict[str, list[tuple[int, float]]],
    sequence_matrix: object,
    target_count: int,
) -> np.ndarray:
    known = train_by_compound.get(query_compound, [])
    if not known:
        return np.zeros(target_count, dtype=np.float32)
    indices = [target for target, _ in known]
    weights = np.asarray([weight for _, weight in known], dtype=np.float32)
    similarities = (sequence_matrix @ sequence_matrix[indices].T).toarray().astype(np.float32, copy=False)
    return np.max(similarities * weights[np.newaxis, :], axis=1)


def weighted_morgan_transfer_scores(
    query_fingerprint: object,
    historical_compounds: list[str],
    train_by_compound: dict[str, list[tuple[int, float]]],
    historical_fingerprints: dict[str, object],
    target_count: int,
) -> np.ndarray:
    scores = np.zeros(target_count, dtype=np.float32)
    if not historical_compounds:
        return scores
    similarities = DataStructs.BulkTanimotoSimilarity(
        query_fingerprint,
        [historical_fingerprints[compound] for compound in historical_compounds],
    )
    for compound, similarity in zip(historical_compounds, similarities):
        for target_idx, weight in train_by_compound[compound]:
            scores[target_idx] = max(scores[target_idx], float(similarity) * float(weight))
    return scores


def precompute_sequence_topk(
    sequence_matrix: object,
    historical_target_indices: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    require(historical_target_indices.size > 0, "At least one historical target is required")
    k = min(PAIR_NEIGHBOR_TOP_K, int(historical_target_indices.size))
    similarities = (sequence_matrix @ sequence_matrix[historical_target_indices].T).toarray().astype(
        np.float32, copy=False
    )
    columns = np.argpartition(similarities, -k, axis=1)[:, -k:]
    values = np.take_along_axis(similarities, columns, axis=1)
    return columns.astype(np.int32, copy=False), values.astype(np.float32, copy=False)


def pair_neighbor_scores(
    query_fingerprint: object,
    historical_compounds: list[str],
    train_by_compound: dict[str, list[tuple[int, float]]],
    historical_fingerprints: dict[str, object],
    target_to_historical_column: dict[int, int],
    sequence_topk_columns: np.ndarray,
    sequence_topk_similarities: np.ndarray,
) -> np.ndarray:
    if not historical_compounds:
        return np.zeros(sequence_topk_columns.shape[0], dtype=np.float32)
    chemical_similarities = DataStructs.BulkTanimotoSimilarity(
        query_fingerprint,
        [historical_fingerprints[compound] for compound in historical_compounds],
    )
    activation = np.zeros(len(target_to_historical_column), dtype=np.float32)
    for compound, similarity in zip(historical_compounds, chemical_similarities):
        for target_idx, evidence_weight in train_by_compound[compound]:
            column = target_to_historical_column[target_idx]
            activation[column] = max(
                activation[column], float(similarity) * float(evidence_weight)
            )
    return np.max(
        sequence_topk_similarities * activation[sequence_topk_columns], axis=1
    ).astype(np.float32, copy=False)


def _summarize_score_diagnostics(values: dict[str, list[int]]) -> dict[str, object]:
    positive = np.asarray(values["positive_target_count"], dtype=float)
    unique = np.asarray(values["unique_score_count"], dtype=float)
    return {
        "query_count": int(positive.size),
        "all_zero_query_count": int(np.sum(positive == 0)),
        "positive_score_target_count_min": int(positive.min()),
        "positive_score_target_count_median": float(np.median(positive)),
        "positive_score_target_count_max": int(positive.max()),
        "unique_score_count_min": int(unique.min()),
        "unique_score_count_median": float(np.median(unique)),
        "unique_score_count_max": int(unique.max()),
    }


def run(args: argparse.Namespace) -> int:
    total_started = time.perf_counter()
    tracemalloc.start()
    phase_seconds: dict[str, float] = {}

    validation_started = time.perf_counter()
    receipt_path = assert_isolated_input(args.receipt)
    receipt = load_receipt(receipt_path)
    inputs = {
        "historical_pairs": assert_isolated_input(args.historical_pairs),
        "scoring_queries": assert_isolated_input(args.scoring_queries),
        "historical_compounds": assert_isolated_input(args.historical_compounds),
        "query_compounds": assert_isolated_input(args.query_compounds),
        "candidate_targets": assert_isolated_input(args.candidate_targets),
        "candidate_sequences": assert_isolated_input(args.candidate_sequences),
    }
    require(
        inputs["historical_compounds"] != inputs["query_compounds"],
        "Historical and query structure roles must be supplied as distinct files",
    )
    manifest_path = assert_isolated_input(args.input_manifest)
    input_manifest = load_input_manifest(manifest_path, inputs, SCORING_INPUT_KIND)
    output_dir = assert_new_output_dir(args.output_dir)

    history_fields, history = read_tsv_gz(inputs["historical_pairs"])
    query_fields, queries = read_tsv_gz(inputs["scoring_queries"])
    historical_structure_fields, historical_structure_rows = read_tsv_gz(
        inputs["historical_compounds"]
    )
    query_structure_fields, query_structure_rows = read_tsv_gz(inputs["query_compounds"])
    target_fields, target_rows = read_tsv_gz(inputs["candidate_targets"])
    sequences = read_fasta(inputs["candidate_sequences"])

    validate_historical_pairs(history_fields, history)
    validate_queries(query_fields, queries)
    target_ids = validate_targets(target_fields, target_rows)
    require(set(target_ids) == set(sequences), "Candidate target TSV and FASTA keysets differ")
    historical_compound_keys = {row["inchikey_full"] for row in history}
    query_compound_keys = {row["inchikey_full"] for row in queries}
    historical_smiles = validate_role_structures(
        historical_structure_fields,
        historical_structure_rows,
        historical_compound_keys,
        "historical",
    )
    query_smiles = validate_role_structures(
        query_structure_fields,
        query_structure_rows,
        query_compound_keys,
        "query",
    )
    phase_seconds["input_load_and_validation"] = time.perf_counter() - validation_started

    feature_started = time.perf_counter()
    target_index = {target: index for index, target in enumerate(target_ids)}
    train_by_compound, weighted_popularity = build_train_maps(history, target_index)
    historical_fingerprints = morgan_fingerprints(historical_smiles, "historical")
    query_fingerprints = morgan_fingerprints(query_smiles, "query")
    sequence_matrix = build_sequence_matrix(target_ids, sequences)
    historical_target_indices = np.asarray(
        sorted({target_idx for pairs in train_by_compound.values() for target_idx, _ in pairs}),
        dtype=np.int32,
    )
    top_columns, top_similarities = precompute_sequence_topk(
        sequence_matrix, historical_target_indices
    )
    historical_column = {
        target_idx: column for column, target_idx in enumerate(historical_target_indices)
    }
    historical_compounds = list(train_by_compound)
    phase_seconds["feature_construction"] = time.perf_counter() - feature_started

    expected_rows_per_baseline = sum(
        EXPECTED_CANDIDATE_TARGETS - len(train_by_compound.get(row["inchikey_full"], []))
        for row in queries
    )
    require(
        expected_rows_per_baseline * len(BASELINES) == EXPECTED_COMPLETE_RANK_ROWS,
        "Corrective candidate masks do not reproduce the frozen complete-rank row count",
    )

    output_dir.mkdir(parents=True, exist_ok=False)
    rank_path = output_dir / "corrective_prediction_ranks.tsv.gz"
    fields = [
        "protocol_id",
        "baseline",
        "query_id",
        "query_compound_inchikey_full",
        "target_uniprot_accession",
        "rank",
        "score",
        "eligible_candidate_target_count",
    ]
    baseline_seconds = {baseline: 0.0 for baseline in BASELINES}
    diagnostics = {
        baseline: {"positive_target_count": [], "unique_score_count": []}
        for baseline in BASELINES
    }
    row_count = 0
    scoring_started = time.perf_counter()
    with rank_path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(text, fieldnames=fields, delimiter="\t", extrasaction="raise")
                writer.writeheader()
                for query_row in queries:
                    query_id = query_row["query_id"]
                    query_compound = query_row["inchikey_full"]
                    allowed = np.ones(len(target_ids), dtype=bool)
                    for target_idx, _ in train_by_compound.get(query_compound, []):
                        allowed[target_idx] = False
                    require(bool(np.any(allowed)), f"Query {query_id} has no eligible candidate target")
                    query_fingerprint = query_fingerprints[query_compound]
                    scorers: dict[str, Callable[[], np.ndarray]] = {
                        "weighted_target_popularity": lambda: weighted_popularity.copy(),
                        "sequence_3mer_transfer": lambda: sequence_transfer_scores(
                            query_compound,
                            train_by_compound,
                            sequence_matrix,
                            len(target_ids),
                        ),
                        "weighted_morgan_transfer": lambda: weighted_morgan_transfer_scores(
                            query_fingerprint,
                            historical_compounds,
                            train_by_compound,
                            historical_fingerprints,
                            len(target_ids),
                        ),
                        "structure_sequence_pair_neighbor": lambda: pair_neighbor_scores(
                            query_fingerprint,
                            historical_compounds,
                            train_by_compound,
                            historical_fingerprints,
                            historical_column,
                            top_columns,
                            top_similarities,
                        ),
                    }
                    for baseline in BASELINES:
                        baseline_started = time.perf_counter()
                        scores = scorers[baseline]()
                        finite_scores(scores, f"{baseline} scores for {query_id}")
                        diagnostics[baseline]["positive_target_count"].append(
                            int(np.count_nonzero(scores[allowed] > 0.0))
                        )
                        diagnostics[baseline]["unique_score_count"].append(
                            int(np.unique(scores[allowed]).size)
                        )
                        _, ranks = rank_scores(
                            scores, allowed, query_id, target_ids, LEGACY_TIE_SALT
                        )
                        require(
                            bool(np.all(ranks[allowed] >= 1) and np.all(ranks[~allowed] == -1)),
                            f"Rank-mask integrity failure for {baseline}, {query_id}",
                        )
                        candidate_count = int(allowed.sum())
                        for target_idx, target_id in enumerate(target_ids):
                            if not allowed[target_idx]:
                                continue
                            writer.writerow(
                                {
                                    "protocol_id": PROTOCOL_ID,
                                    "baseline": baseline,
                                    "query_id": query_id,
                                    "query_compound_inchikey_full": query_compound,
                                    "target_uniprot_accession": target_id,
                                    "rank": int(ranks[target_idx]),
                                    "score": f"{float(scores[target_idx]):.17g}",
                                    "eligible_candidate_target_count": candidate_count,
                                }
                            )
                            row_count += 1
                        baseline_seconds[baseline] += time.perf_counter() - baseline_started
    phase_seconds["scoring_ranking_and_rank_write"] = time.perf_counter() - scoring_started
    require(row_count == EXPECTED_COMPLETE_RANK_ROWS, "Corrective complete-rank row count changed")

    python_current, python_peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    peak_rss, peak_rss_method = peak_rss_bytes()
    phase_seconds["total_before_manifest_write"] = time.perf_counter() - total_started

    shared_rank_path = WORKSPACE / "scripts" / "pu_retrieval_metrics.py"
    executed_entrypoint = Path(sys.argv[0]).resolve()
    executed_code = [Path(__file__), ROOT / "scripts" / "asof_successor_common.py", shared_rank_path]
    if executed_entrypoint.is_file() and executed_entrypoint not in {path.resolve() for path in executed_code}:
        executed_code.append(executed_entrypoint)
    manifest_output = output_dir / "corrective_score_manifest.json"
    write_json(
        manifest_output,
        {
            "schema_version": "1.0",
            "protocol_id": PROTOCOL_ID,
            "run_id": input_manifest.get("run_id"),
            "stage": "corrective_score",
            "execution_mode": RUN_MODE,
            "claim_boundary": "Author-run non-independent post hoc temporal-purity correction; no blinded or external-validation claim.",
            "command_argv": list(sys.argv),
            "receipt": {"path": str(receipt_path), "sha256": sha256(receipt_path)},
            "input_manifest": {"path": str(manifest_path), "sha256": sha256(manifest_path)},
            "inputs": {
                name: {"path": str(path), "sha256": sha256(path)} for name, path in inputs.items()
            },
            "role_separation": {
                "historical_and_query_files_distinct": True,
                "historical_compound_count": len(historical_smiles),
                "query_compound_count": len(query_smiles),
                "shared_inchikey_count": len(set(historical_smiles).intersection(query_smiles)),
                "pooled_structure_map_used": False,
            },
            "rank_output": {"path": str(rank_path), "sha256": sha256(rank_path)},
            "target_count": len(target_ids),
            "historical_target_count": int(historical_target_indices.size),
            "query_count": len(queries),
            "row_count": row_count,
            "baselines": BASELINES,
            "parameters": {
                "strict_tier_weights": STRICT_TIERS,
                "morgan_radius": MORGAN_RADIUS,
                "morgan_bits": MORGAN_BITS,
                "sequence_feature": "character_3mer_tfidf_cosine",
                "pair_neighbor_sequence_top_k": PAIR_NEIGHBOR_TOP_K,
                "tie_salt": LEGACY_TIE_SALT,
                "mask_historical_targets_for_same_query": True,
            },
            "score_diagnostics": {
                baseline: _summarize_score_diagnostics(diagnostics[baseline])
                for baseline in BASELINES
            },
            "runtime": {
                "wall_seconds_by_phase": phase_seconds,
                "wall_seconds_by_baseline_including_rank_write": baseline_seconds,
                "python_tracemalloc_current_bytes": python_current,
                "python_tracemalloc_peak_bytes": python_peak,
                "process_peak_rss_bytes": peak_rss,
                "process_peak_rss_method": peak_rss_method,
            },
            "environment": environment_receipt(),
            "code": code_hashes(executed_code),
            "endpoint_file_supplied_to_score_command": False,
            "endpoint_read_by_score_engine": False,
            "legacy_outer_or_result_read": False,
            "figures_generated": False,
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
        },
    )
    print(
        json.dumps(
            {
                "output_dir": str(output_dir),
                "rank_rows": row_count,
                "endpoint_read": False,
                "figures_generated": False,
            }
        )
    )
    return 0


def main(argv: list[str] | None = None) -> int:
    return run(build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
