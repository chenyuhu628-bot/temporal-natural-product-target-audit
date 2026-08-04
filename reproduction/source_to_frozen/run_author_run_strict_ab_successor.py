"""Run a transparent, non-independent strict A/B successor calculation.

This runner is intentionally distinct from the role-separated successor route.
It validates the same strict A/B schemas, four fixed baselines, five scopes,
five metrics, and pre-fixed bootstrap, but makes no claim of a human gate,
independent evaluator, endpoint access control, or blinded external validation.
It produces no figures.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import platform
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent
sys.path.insert(0, str(WORKSPACE / "scripts"))

from pu_baseline_io import read_fasta, read_tsv_gz, sha256 as legacy_sha256
from pu_pair_neighbor_transfer_scores_v1_2 import pair_neighbor_scores, precompute_sequence_topk
from pu_retrieval_metrics import macro_average, query_metrics, rank_scores
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
    parse_bool,
    require_fields,
    require_unique,
    sha256,
    validate_compounds,
    validate_historical_pairs,
    validate_queries,
    validate_targets,
)
from evaluate_successor_sealed import (
    BASELINES,
    METRICS,
    SCOPES,
    bootstrap_baselines,
    bootstrap_pairs,
    build_scope_relevance,
    load_endpoint,
    load_prediction_ranks,
    read_bool_map,
    write_tsv_gz,
)


RUN_MODE = "author_run_non_independent"
SALT = "npass_strict_ab_doublecold_successor_v1_20260719"


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)

    score = commands.add_parser("score", help="Score the four fixed baselines without an endpoint argument.")
    score.add_argument("--receipt", required=True, type=Path)
    score.add_argument("--input-manifest", required=True, type=Path)
    score.add_argument("--historical-pairs", required=True, type=Path)
    score.add_argument("--scoring-queries", required=True, type=Path)
    score.add_argument("--compounds", required=True, type=Path)
    score.add_argument("--candidate-targets", required=True, type=Path)
    score.add_argument("--candidate-sequences", required=True, type=Path)
    score.add_argument("--output-dir", required=True, type=Path)

    evaluate = commands.add_parser("evaluate", help="Evaluate the fixed score artifact in author-run mode.")
    evaluate.add_argument("--receipt", required=True, type=Path)
    evaluate.add_argument("--input-manifest", required=True, type=Path)
    evaluate.add_argument("--score-dir", required=True, type=Path)
    evaluate.add_argument("--endpoint", required=True, type=Path)
    evaluate.add_argument("--scaffold-audit", required=True, type=Path)
    evaluate.add_argument("--homology-0-30", required=True, type=Path)
    evaluate.add_argument("--homology-0-50", required=True, type=Path)
    evaluate.add_argument("--homology-0-70", required=True, type=Path)
    evaluate.add_argument("--output-dir", required=True, type=Path)
    return result


def load_author_run_receipt(path: Path) -> dict[str, Any]:
    receipt = json.loads(path.read_text(encoding="utf-8"))
    if receipt.get("protocol_id") != PROTOCOL_ID:
        raise ValueError("Author-run receipt protocol ID mismatch")
    if receipt.get("execution_mode") != RUN_MODE:
        raise ValueError("Receipt does not declare the author-run non-independent mode")
    if receipt.get("project_lead_authorized_internal_use") is not True:
        raise ValueError("Receipt lacks the recorded project-lead internal-use authorization")
    if receipt.get("human_gate_status") != "not_claimed":
        raise ValueError("Author-run receipt may not claim completed human gates")
    if receipt.get("release_status") != "restricted_internal_no_public_release":
        raise ValueError("Author-run receipt lacks the restricted-release boundary")
    receipt["_path"] = str(path)
    return receipt


def load_author_run_manifest(
    path: Path,
    paths: dict[str, Path],
    expected_kind: str,
) -> dict[str, Any]:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("protocol_id") != PROTOCOL_ID:
        raise ValueError("Author-run input manifest protocol ID mismatch")
    if manifest.get("execution_mode") != RUN_MODE:
        raise ValueError("Author-run input manifest mode mismatch")
    if manifest.get("project_lead_authorized_internal_use") is not True:
        raise ValueError("Author-run input manifest lacks internal-use authorization")
    if manifest.get("input_kind") != expected_kind:
        raise ValueError(f"Author-run input manifest kind is not {expected_kind!r}")
    if manifest.get("legacy_outer_or_result_input") is not False:
        raise ValueError("Author-run input manifest does not exclude legacy result/outer inputs")
    declared = manifest.get("file_sha256", {})
    for name, source in paths.items():
        if declared.get(name) != sha256(source):
            raise ValueError(f"Author-run input manifest hash mismatch for {name}")
    manifest["_path"] = str(path)
    return manifest


def write_json(path: Path, payload: dict[str, Any]) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite output: {path}")
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def score(args: argparse.Namespace) -> int:
    receipt = load_author_run_receipt(assert_isolated_input(args.receipt))
    inputs = {
        "historical_pairs": assert_isolated_input(args.historical_pairs),
        "scoring_queries": assert_isolated_input(args.scoring_queries),
        "compounds": assert_isolated_input(args.compounds),
        "candidate_targets": assert_isolated_input(args.candidate_targets),
        "candidate_sequences": assert_isolated_input(args.candidate_sequences),
    }
    input_manifest = load_author_run_manifest(
        assert_isolated_input(args.input_manifest), inputs, "scoring_without_endpoint_file"
    )
    output_dir = assert_isolated_input(args.output_dir)
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite author-run score directory: {output_dir}")

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
        raise ValueError(f"Candidate target count must be 4,123, found {len(target_ids)}")
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
        sorted({target_index[row["uniprot_canonical_accession"]] for row in history}), dtype=np.int32
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
    rank_path = output_dir / "author_run_prediction_ranks.tsv.gz"
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
                        _, ranks = rank_scores(scores, allowed, query_id, target_ids, SALT)
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
    manifest_path = output_dir / "author_run_score_manifest.json"
    write_json(
        manifest_path,
        {
            "protocol_id": PROTOCOL_ID,
            "run_id": input_manifest["run_id"],
            "stage": "author_run_score",
            "execution_mode": RUN_MODE,
            "claim_boundary": "No independent evaluator, human-gate completion, endpoint access control, or blinded external-validation claim.",
            "receipt": {"path": receipt["_path"], "sha256": sha256(Path(receipt["_path"]))},
            "input_manifest": {"path": input_manifest["_path"], "sha256": sha256(Path(input_manifest["_path"]))},
            "inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in inputs.items()},
            "rank_output": {"path": str(rank_path), "sha256": sha256(rank_path)},
            "target_count": len(target_ids),
            "query_count": len(queries),
            "row_count": row_count,
            "baselines": list(scorers),
            "software_versions": software_versions(),
            "environment": {"python": sys.version, "platform": platform.platform()},
            "endpoint_file_supplied_to_score_command": False,
            "endpoint_read_by_score_engine": False,
            "endpoint_access_control_claim": "none",
            "legacy_outer_or_result_read": False,
            "figures_generated": False,
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
        },
    )
    print(json.dumps({"output_dir": str(output_dir), "rank_rows": row_count, "figures_generated": False}))
    return 0


def evaluate(args: argparse.Namespace) -> int:
    receipt = load_author_run_receipt(assert_isolated_input(args.receipt))
    input_paths = {
        "endpoint": assert_isolated_input(args.endpoint),
        "scaffold_audit": assert_isolated_input(args.scaffold_audit),
        "homology_0_30": assert_isolated_input(args.homology_0_30),
        "homology_0_50": assert_isolated_input(args.homology_0_50),
        "homology_0_70": assert_isolated_input(args.homology_0_70),
    }
    evaluation_input_manifest = load_author_run_manifest(
        assert_isolated_input(args.input_manifest), input_paths, "author_run_evaluation_endpoint"
    )
    score_dir = assert_isolated_input(args.score_dir)
    output_dir = assert_isolated_input(args.output_dir)
    if output_dir.exists():
        raise FileExistsError(f"Refusing to overwrite author-run evaluation directory: {output_dir}")
    score_manifest_path = score_dir / "author_run_score_manifest.json"
    rank_path = score_dir / "author_run_prediction_ranks.tsv.gz"
    if not score_manifest_path.is_file() or not rank_path.is_file():
        raise FileNotFoundError("Author-run score directory is incomplete")
    score_manifest = json.loads(score_manifest_path.read_text(encoding="utf-8"))
    if (
        score_manifest.get("protocol_id") != PROTOCOL_ID
        or score_manifest.get("stage") != "author_run_score"
        or score_manifest.get("execution_mode") != RUN_MODE
        or score_manifest.get("endpoint_file_supplied_to_score_command") is not False
        or score_manifest.get("endpoint_read_by_score_engine") is not False
        or score_manifest.get("endpoint_access_control_claim") != "none"
    ):
        raise ValueError("Author-run score manifest does not meet the transparent no-endpoint-argument contract")
    if score_manifest.get("rank_output", {}).get("sha256") != sha256(rank_path):
        raise ValueError("Author-run score-rank hash mismatch")

    endpoint = load_endpoint(input_paths["endpoint"])
    scaffold = read_bool_map(
        input_paths["scaffold_audit"],
        "canonical_pair_key",
        "audit_scaffold_cold_under_selected_policy",
        "audit_outcome",
        "author-run scaffold audit",
    )
    homology = {
        "0_30": read_bool_map(
            input_paths["homology_0_30"],
            "uniprot_canonical_accession",
            "is_future_target_homology_cold_candidate",
            "future_target_coldness_status",
            "author-run homology 0.30",
        ),
        "0_50": read_bool_map(
            input_paths["homology_0_50"],
            "uniprot_canonical_accession",
            "is_future_target_homology_cold_candidate",
            "future_target_coldness_status",
            "author-run homology 0.50",
        ),
        "0_70": read_bool_map(
            input_paths["homology_0_70"],
            "uniprot_canonical_accession",
            "is_future_target_homology_cold_candidate",
            "future_target_coldness_status",
            "author-run homology 0.70",
        ),
    }
    scoped = build_scope_relevance(endpoint, scaffold, homology)
    needed_targets: dict[str, set[str]] = defaultdict(set)
    for endpoint_row in endpoint:
        needed_targets[endpoint_row["query_id"]].add(endpoint_row["uniprot_canonical_accession"])
    ranks, rank_audit, query_compounds = load_prediction_ranks(
        rank_path, int(score_manifest["row_count"]), needed_targets
    )
    for endpoint_row in endpoint:
        query_id = endpoint_row["query_id"]
        if query_compounds.get(query_id) != endpoint_row["inchikey_full"]:
            raise ValueError(f"Endpoint/query compound mismatch for {query_id}")

    scope_audit_rows: list[dict[str, object]] = []
    for scope in SCOPES:
        scope_rows = [row for rows in scoped[scope].values() for row in rows]
        scope_audit_rows.append(
            {
                "scope": scope,
                "candidate_relation_count": len(scope_rows),
                "query_compound_count": len(scoped[scope]),
                "target_count": len({row["uniprot_canonical_accession"] for row in scope_rows}),
                "A_affinity_candidate_count": sum(
                    row["best_strict_evidence_tier"] == "A_affinity_candidate" for row in scope_rows
                ),
                "B_quantitative_functional_candidate_count": sum(
                    row["best_strict_evidence_tier"] == "B_quantitative_functional_candidate" for row in scope_rows
                ),
                "status": "estimable" if scoped[scope] else "not_estimable",
            }
        )

    metric_rows: list[dict[str, object]] = []
    baseline_bootstrap_rows: list[dict[str, object]] = []
    bootstrap_rows: list[dict[str, object]] = []
    generator = np.random.Generator(np.random.PCG64(20260719))
    for scope in SCOPES:
        relevant = scoped[scope]
        queries = sorted(relevant)
        arrays: dict[str, np.ndarray] = {}
        for baseline in BASELINES:
            per_query: list[dict[str, float]] = []
            candidate_counts: list[int] = []
            for query_id in queries:
                positive_ranks: list[int] = []
                for endpoint_row in relevant[query_id]:
                    target = endpoint_row["uniprot_canonical_accession"]
                    try:
                        rank = ranks[baseline][query_id][target]
                    except KeyError as exc:
                        raise ValueError(
                            f"Endpoint target lacks an author-run score for {baseline}, {query_id}, {target}"
                        ) from exc
                    positive_ranks.append(rank)
                per_query.append(query_metrics(positive_ranks, (10, 50)))
                candidate_counts.append(rank_audit[(baseline, query_id)]["candidate_count"])
            if per_query:
                summary = macro_average(per_query)
                arrays[baseline] = np.asarray(
                    [[row[metric] for metric in METRICS] for row in per_query], dtype=float
                )
                zero_at_10 = sum(row["Recall@10"] == 0.0 for row in per_query)
                zero_at_50 = sum(row["Recall@50"] == 0.0 for row in per_query)
                zero_mrr = sum(row["MRR"] == 0.0 for row in per_query)
                candidate_values = np.asarray(candidate_counts, dtype=float)
                metric_rows.append(
                    {
                        "scope": scope,
                        "baseline": baseline,
                        "evaluable_query_count": len(per_query),
                        "candidate_relation_count": sum(len(relevant[item]) for item in queries),
                        "candidate_target_median": f"{np.median(candidate_values):.17g}",
                        "candidate_target_iqr": f"{(np.percentile(candidate_values, 75) - np.percentile(candidate_values, 25)):.17g}",
                        "candidate_target_min": int(candidate_values.min()),
                        "candidate_target_max": int(candidate_values.max()),
                        **{metric: f"{summary[metric]:.17g}" for metric in METRICS},
                        "zero_recall_at_10_queries": zero_at_10,
                        "zero_recall_at_50_queries": zero_at_50,
                        "zero_mrr_queries": zero_mrr,
                        "status": "estimable" if len(per_query) > 1 else "n=1_descriptive_only",
                    }
                )
            else:
                arrays[baseline] = np.empty((0, len(METRICS)), dtype=float)
                metric_rows.append(
                    {
                        "scope": scope,
                        "baseline": baseline,
                        "evaluable_query_count": 0,
                        "candidate_relation_count": 0,
                        "candidate_target_median": "",
                        "candidate_target_iqr": "",
                        "candidate_target_min": "",
                        "candidate_target_max": "",
                        **{metric: "" for metric in METRICS},
                        "zero_recall_at_10_queries": 0,
                        "zero_recall_at_50_queries": 0,
                        "zero_mrr_queries": 0,
                        "status": "not_estimable",
                    }
                )
        n_scope_queries = len(queries)
        indices = (
            generator.integers(0, n_scope_queries, size=(10000, n_scope_queries), endpoint=False)
            if n_scope_queries > 1
            else None
        )
        baseline_bootstrap_rows.extend(bootstrap_baselines(arrays, scope, indices))
        bootstrap_rows.extend(bootstrap_pairs(arrays, scope, indices))

    output_dir.mkdir(parents=True, exist_ok=False)
    metric_path = output_dir / "author_run_aggregate_metrics.tsv.gz"
    scope_audit_path = output_dir / "author_run_scope_denominator_audit.tsv.gz"
    baseline_bootstrap_path = output_dir / "author_run_baseline_bootstrap_metrics.tsv.gz"
    bootstrap_path = output_dir / "author_run_paired_bootstrap_contrasts.tsv.gz"
    write_tsv_gz(metric_path, list(metric_rows[0]), metric_rows)
    write_tsv_gz(scope_audit_path, list(scope_audit_rows[0]), scope_audit_rows)
    write_tsv_gz(baseline_bootstrap_path, list(baseline_bootstrap_rows[0]), baseline_bootstrap_rows)
    contrast_fields = list(bootstrap_rows[0]) if bootstrap_rows else [
        "scope",
        "left_baseline",
        "right_baseline",
        "metric",
        "query_count",
        "mean_difference_left_minus_right",
        "ci95_low",
        "ci95_high",
        "status",
    ]
    write_tsv_gz(bootstrap_path, contrast_fields, bootstrap_rows)
    manifest_path = output_dir / "author_run_evaluation_manifest.json"
    write_json(
        manifest_path,
        {
            "protocol_id": PROTOCOL_ID,
            "run_id": evaluation_input_manifest["run_id"],
            "stage": "author_run_evaluation",
            "execution_mode": RUN_MODE,
            "claim_boundary": "Author-run non-independent calculation; no independent evaluator, human-gate completion, endpoint access control, or blinded external-validation claim.",
            "receipt": {"path": receipt["_path"], "sha256": sha256(Path(receipt["_path"]))},
            "evaluation_input_manifest": {
                "path": evaluation_input_manifest["_path"],
                "sha256": sha256(Path(evaluation_input_manifest["_path"])),
            },
            "score_manifest": {"path": str(score_manifest_path), "sha256": sha256(score_manifest_path)},
            "score_rank": {"path": str(rank_path), "sha256": sha256(rank_path)},
            "evaluation_inputs": {name: {"path": str(path), "sha256": sha256(path)} for name, path in input_paths.items()},
            "outputs": {
                "aggregate_metrics": {"path": str(metric_path), "sha256": sha256(metric_path)},
                "scope_denominator_audit": {"path": str(scope_audit_path), "sha256": sha256(scope_audit_path)},
                "baseline_bootstrap": {"path": str(baseline_bootstrap_path), "sha256": sha256(baseline_bootstrap_path)},
                "paired_bootstrap": {"path": str(bootstrap_path), "sha256": sha256(bootstrap_path)},
            },
            "all_scopes_reported": SCOPES,
            "all_baselines_reported": BASELINES,
            "all_metrics_reported": METRICS,
            "bootstrap": {"replicates": 10000, "prng": "PCG64", "seed": 20260719, "interval": "95% percentile"},
            "endpoint_access_control_claim": "none",
            "figures_generated": False,
            "created_at_utc": datetime.now(timezone.utc).isoformat(),
        },
    )
    print(json.dumps({"output_dir": str(output_dir), "figures_generated": False}))
    return 0


def main() -> int:
    args = parser().parse_args()
    return score(args) if args.command == "score" else evaluate(args)


if __name__ == "__main__":
    raise SystemExit(main())
