"""Generate strict A/B successor predictions without endpoint access.

This runner is deliberately separate from the legacy strict runner. It scores all
four locked baselines over every query/candidate target and never reads sealed
future targets, scaffold masks, homology masks, old outer ledgers, or results.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
sys.path.insert(0, str(WORKSPACE / "scripts"))

from pu_baseline_io import read_fasta, read_tsv_gz, sha256 as legacy_sha256
from pu_pair_neighbor_transfer_scores_v1_2 import pair_neighbor_scores, precompute_sequence_topk
from pu_retrieval_metrics import rank_scores
from pu_retrieval_scores import (
    build_sequence_kmer_matrix,
    build_train_maps,
    morgan_fingerprints,
    sequence_transfer_scores,
    software_versions,
    tanimoto_transfer_scores,
)
from successor_common import (
    PROTOCOL_ID,
    STRICT_TIERS,
    assert_isolated_input,
    finite_scores,
    input_manifest,
    load_authorized_manifest,
    load_receipt,
    validate_compounds,
    validate_historical_pairs,
    validate_queries,
    validate_targets,
)


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--receipt", required=True, type=Path)
    result.add_argument("--authorized-input-manifest", required=True, type=Path)
    result.add_argument("--historical-pairs", required=True, type=Path)
    result.add_argument("--scoring-queries", required=True, type=Path)
    result.add_argument("--compounds", required=True, type=Path)
    result.add_argument("--candidate-targets", required=True, type=Path)
    result.add_argument("--candidate-sequences", required=True, type=Path)
    result.add_argument("--output-dir", required=True, type=Path)
    return result


def main() -> int:
    args = parser().parse_args()
    receipt = load_receipt(args.receipt, "score")
    if receipt.get("actor_role") not in {"analysis_executor", "AI_assisted_analysis_executor"}:
        raise ValueError("Score receipt actor_role must be analysis_executor or AI_assisted_analysis_executor")

    inputs = {
        "historical_pairs": assert_isolated_input(args.historical_pairs),
        "scoring_queries": assert_isolated_input(args.scoring_queries),
        "compounds": assert_isolated_input(args.compounds),
        "candidate_targets": assert_isolated_input(args.candidate_targets),
        "candidate_sequences": assert_isolated_input(args.candidate_sequences),
    }
    input_authorization = load_authorized_manifest(
        assert_isolated_input(args.authorized_input_manifest),
        inputs,
        {
            "authorized_internal_successor_use": True,
            "sealed_endpoint_included": False,
            "legacy_outer_or_result_input": False,
        },
    )
    output_dir = assert_isolated_input(args.output_dir)
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite successor prediction directory: {output_dir}")

    history_fields, history = read_tsv_gz(inputs["historical_pairs"])
    query_fields, queries = read_tsv_gz(inputs["scoring_queries"])
    compound_fields, compounds = read_tsv_gz(inputs["compounds"])
    target_fields, target_rows = read_tsv_gz(inputs["candidate_targets"])
    sequences = read_fasta(inputs["candidate_sequences"])

    validate_historical_pairs(history_fields, history)
    validate_queries(query_fields, queries)
    validate_compounds(compound_fields, compounds)
    target_ids = validate_targets(target_fields, target_rows)
    if len(target_ids) != 4123:
        raise ValueError(f"Candidate target count must be 4123, found {len(target_ids)}")
    if set(target_ids) != set(sequences):
        raise ValueError("Candidate target TSV and FASTA accession sets differ")

    all_smiles = {row["inchikey_full"]: row["representative_smiles"] for row in compounds}
    required_compounds = {row["inchikey_full"] for row in history} | {row["inchikey_full"] for row in queries}
    missing_compounds = sorted(required_compounds.difference(all_smiles))
    if missing_compounds:
        raise ValueError(f"Compound structure table lacks required InChIKeys; first: {missing_compounds[:5]}")
    smiles = {compound: all_smiles[compound] for compound in required_compounds}

    target_index = {target: index for index, target in enumerate(target_ids)}
    helper_history = [
        {
            "inchikey_full": row["inchikey_full"],
            "uniprot_canonical_accession": row["uniprot_canonical_accession"],
            "best_strict_evidence_tier_v1_1": row["best_strict_evidence_tier"],
        }
        for row in history
    ]
    train_by_compound, weighted_popularity = build_train_maps(helper_history, target_index, STRICT_TIERS)
    fingerprints = morgan_fingerprints(smiles, 2, 2048)
    _, sequence_matrix = build_sequence_kmer_matrix(target_ids, sequences)
    historical_target_indices = np.asarray(
        sorted({target_index[row["uniprot_canonical_accession"]] for row in history}),
        dtype=np.int32,
    )
    top_columns, top_similarities = precompute_sequence_topk(sequence_matrix, historical_target_indices, 100)
    historical_column = {target_index: column for column, target_index in enumerate(historical_target_indices)}

    scorers: dict[str, Callable[[str], np.ndarray]] = {
        "weighted_target_popularity": lambda query: weighted_popularity.copy(),
        "sequence_3mer_transfer": lambda query: sequence_transfer_scores(
            query, train_by_compound, sequence_matrix, len(target_ids)
        ),
        "weighted_morgan_transfer": lambda query: tanimoto_transfer_scores(
            query, train_by_compound, fingerprints, len(target_ids)
        ),
        "structure_sequence_pair_neighbor": lambda query: pair_neighbor_scores(
            query,
            train_by_compound,
            fingerprints,
            historical_column,
            top_columns,
            top_similarities,
        ),
    }

    output_dir.mkdir(parents=True, exist_ok=False)
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
    rank_path = output_dir / "blind_prediction_ranks.tsv.gz"
    row_count = 0
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
                    if not np.any(allowed):
                        raise ValueError(f"Query {query_id} has no eligible candidate target")
                    for baseline, scorer in scorers.items():
                        scores = scorer(query_compound)
                        finite_scores(scores, f"{baseline} scores for {query_id}")
                        _, ranks = rank_scores(
                            scores,
                            allowed,
                            query_id,
                            target_ids,
                            "npass_strict_ab_doublecold_successor_v1_20260719",
                        )
                        if np.any(ranks[allowed] < 1) or np.any(ranks[~allowed] != -1):
                            raise ValueError(f"Rank-mask integrity failure for {baseline}, {query_id}")
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
                                    "eligible_candidate_target_count": int(allowed.sum()),
                                }
                            )
                            row_count += 1
    run_manifest = input_manifest(inputs, receipt)
    run_manifest.update(
        {
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
            "stage": "blind_prediction",
            "protocol_id": PROTOCOL_ID,
            "script": str(Path(__file__).resolve()),
            "script_sha256": legacy_sha256(Path(__file__)),
            "authorized_input_manifest": {
                "path": input_authorization["_path"],
                "sha256": legacy_sha256(Path(input_authorization["_path"])),
            },
            "rank_output": {"path": str(rank_path), "sha256": legacy_sha256(rank_path)},
            "target_count": len(target_ids),
            "query_count": len(queries),
            "row_count": row_count,
            "baselines": list(scorers),
            "software_versions": software_versions(),
            "environment": {"python": sys.version, "platform": platform.platform()},
            "sealed_endpoint_read": False,
            "legacy_outer_or_result_read": False,
        }
    )
    manifest_path = output_dir / "blind_prediction_manifest.json"
    manifest_path.write_text(json.dumps(run_manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"output_dir": str(output_dir), "rank_rows": row_count, "endpoint_read": False}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
