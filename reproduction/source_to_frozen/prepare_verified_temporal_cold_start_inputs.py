#!/usr/bin/env python3
"""Prepare reproducible cold-start inputs from the frozen strict temporal tables.

This is intentionally separate from the earlier broad-P1 cold-start input
preparation.  It starts only from the date-verified pair tables tagged
``v1_1_pmid_verified``: 4,990 pre-cutoff training pairs and 442 post-cutoff
future pairs.  It never creates a negative label and never changes raw data.

For historical pairs, structural provenance is restricted to NPASS v2.  For
future pairs it is restricted to NPASS v3.  The source table's archived
UniProt FASTA is filtered locally; this script makes no network request.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import shutil
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


TAG = "v1_1_pmid_verified"
PRIMARY_TIERS = {"A_affinity_candidate", "B_quantitative_functional_candidate"}
P1 = "P1_npass_raw_exact_candidate"
STRICT_MAPPING = "strict_one_to_one_reviewed_human"
EXPECTED_COUNTS = {"training": 4990, "future": 442}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def membership_digest(keys: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for key in sorted(keys):
        digest.update(key.encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest()


def read_delimited_gz(path: Path, delimiter: str) -> Iterable[dict[str, str]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        yield from csv.DictReader(handle, delimiter=delimiter)


def write_tsv_gz(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> None:
    """Write byte-stable gzip TSV output (mtime=0) for a deterministic run."""
    with path.open("wb") as binary:
        with gzip.GzipFile(filename="", mode="wb", fileobj=binary, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(text, fieldnames=fields, delimiter="\t", extrasaction="ignore")
                writer.writeheader()
                writer.writerows(rows)


def parse_bool(value: str) -> bool:
    return value.strip().casefold() in {"true", "1", "yes"}


def canonical_key(row: dict[str, str]) -> str:
    return f"{row.get('inchikey_full', '').strip()}|{row.get('uniprot_canonical_accession', '').strip()}"


def normalize_leakage_gated_future_row(row: dict[str, str], split: str) -> dict[str, str]:
    """Restore canonical temporal column names from the immutable C31 gate.

    The leakage gate prefixes frozen-future columns with ``future_`` so its
    provenance cannot be confused with its audit columns.  Downstream source
    selection needs the original names, so this adapter copies (never edits)
    those fields into a derived in-memory row.
    """
    if split != "future" or row.get("decision", "").strip() or not row.get("future_decision", "").strip():
        return row
    aliases = {
        "source_versions": "future_source_versions",
        "source_pair_keys": "future_source_pair_keys",
        "source_np_ids": "future_source_np_ids",
        "source_target_ids": "future_source_target_ids",
        "v2_all_record_count": "future_v2_all_record_count",
        "v3_all_record_count": "future_v3_all_record_count",
        "v2_primary_A_B_P1_record_count": "future_v2_primary_A_B_P1_record_count",
        "v3_primary_A_B_P1_record_count": "future_v3_primary_A_B_P1_record_count",
        "v2_strict_entity_primary_record_count": "future_v2_strict_entity_primary_record_count",
        "v3_strict_entity_primary_record_count": "future_v3_strict_entity_primary_record_count",
        "v2_day_precise_pre_cutoff_primary_record_count": "future_v2_day_precise_pre_cutoff_primary_record_count",
        "primary_evidence_tiers": "future_primary_evidence_tiers",
        "primary_activity_types": "future_primary_activity_types",
        "primary_references": "future_primary_references",
        "v3_cross_version_statuses": "future_v3_cross_version_statuses",
        "v3_temporal_screen_statuses": "future_v3_temporal_screen_statuses",
        "decision": "future_decision",
        "decision_rationale": "future_decision_rationale",
        "exclusion_or_holdout_reasons": "future_exclusion_or_holdout_reasons",
        "label_status": "future_label_status",
        "unrecorded_pair_policy": "future_unrecorded_pair_policy",
    }
    normalized = dict(row)
    for canonical, prefixed in aliases.items():
        normalized[canonical] = row.get(prefixed, "")
    return normalized


def load_strict_table(
    path: Path,
    split: str,
    expected_decision: str,
    expected_count: int | None,
) -> tuple[list[str], dict[str, dict[str, str]]]:
    rows = [normalize_leakage_gated_future_row(row, split) for row in read_delimited_gz(path, ",")]
    if expected_count is not None and len(rows) != expected_count:
        raise ValueError(f"{path.name}: expected {expected_count} {split} rows, found {len(rows)}")
    if not rows:
        raise ValueError(f"{path.name}: no rows")
    fields = list(rows[0])
    indexed: dict[str, dict[str, str]] = {}
    for row in rows:
        key = row.get("canonical_pair_key", "").strip()
        if not key or key != canonical_key(row):
            raise ValueError(f"{path.name}: malformed canonical_pair_key {key!r}")
        if row.get("decision", "") != expected_decision:
            raise ValueError(f"{path.name}: {key} has unexpected decision {row.get('decision')!r}")
        if key in indexed:
            raise ValueError(f"{path.name}: duplicate canonical pair {key}")
        indexed[key] = row
    return fields, indexed


def fasta_accession(header: str) -> str:
    text = header[1:].strip()
    parts = text.split("|", 2)
    return parts[1].strip() if len(parts) >= 2 else text.split(None, 1)[0]


def read_archived_fasta(path: Path) -> dict[str, str]:
    sequences: dict[str, str] = {}
    accession = ""
    chunks: list[str] = []
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if accession:
                    sequence = "".join(chunks)
                    if sequences.get(accession, sequence) != sequence:
                        raise ValueError(f"Conflicting archived sequences for {accession}")
                    sequences[accession] = sequence
                accession, chunks = fasta_accession(line), []
            else:
                chunks.append(line)
    if accession:
        sequence = "".join(chunks)
        if sequences.get(accession, sequence) != sequence:
            raise ValueError(f"Conflicting archived sequences for {accession}")
        sequences[accession] = sequence
    return sequences


def is_strict_primary(row: dict[str, str]) -> bool:
    return (
        row.get("automatic_verification_level", "") == P1
        and row.get("evidence_tier_v1_1", "") in PRIMARY_TIERS
        and row.get("activity_relation", "").strip() == "="
        and row.get("mapping_status", "") == STRICT_MAPPING
        and parse_bool(row.get("sequence_found", ""))
    )


def compact_source_row(row: dict[str, str], split: str) -> dict[str, str]:
    keep = [
        "source_version", "source_np_id", "source_target_id", "inchikey_full", "inchikey_connectivity", "smiles",
        "uniprot_raw", "uniprot_canonical_accession", "target_name", "target_type", "target_tax_id",
        "activity_relation", "activity_type", "activity_value", "activity_units", "activity_value_molar", "p_activity",
        "evidence_tier_v1_1", "evidence_weight_v1_1", "measurement_class", "ref_id_type", "ref_id",
        "pair_key", "mapping_status", "sequence_found", "sequence_md5",
    ]
    return {"split": split, "canonical_pair_key": canonical_key(row), **{field: row.get(field, "") for field in keep}}


def values_with_counts(values: Counter[str]) -> str:
    return json.dumps(dict(sorted(values.items())), ensure_ascii=False, separators=(",", ":"))


def safe_prepare_output(output_dir: Path, allowed_root: Path, overwrite: bool) -> None:
    allowed_root = allowed_root.resolve()
    output_dir = output_dir.resolve()
    try:
        output_dir.relative_to(allowed_root)
    except ValueError as exc:
        raise ValueError(f"Output must be inside {allowed_root}: {output_dir}") from exc
    if output_dir.exists():
        if not overwrite:
            raise FileExistsError(f"Refusing to overwrite existing output directory: {output_dir}")
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)


def run(
    root: Path,
    output_dir: Path,
    overwrite: bool,
    training_pairs_path: Path | None,
    future_pairs_path: Path | None,
    tag: str,
    archive_fasta_path: Path | None,
    uniprot_summary_path: Path | None,
) -> dict[str, object]:
    processed = root / "data" / "processed"
    source_paths = {
        "v2": processed / "npass_v2_evidence_records_v1_1_uniprot_mapped.tsv.gz",
        "v3": processed / "npass_v3_evidence_records_v1_1_uniprot_mapped.tsv.gz",
    }
    strict_paths = {
        "training": (training_pairs_path or processed / f"strict_temporal_training_candidates_{TAG}.csv.gz").resolve(),
        "future": (future_pairs_path or processed / f"strict_temporal_future_candidates_{TAG}.csv.gz").resolve(),
    }
    archive_fasta = (
        archive_fasta_path
        or root / "data/raw/uniprot/id_mapping_20260715T083054Z/strict_reviewed_human_sequences.fasta.gz"
    ).resolve()
    uniprot_summary = (
        uniprot_summary_path or root / "results/uniprot_target_mapping_summary.json"
    ).resolve()
    required = [*source_paths.values(), *strict_paths.values(), archive_fasta, uniprot_summary]
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)

    training_fields, training = load_strict_table(
        strict_paths["training"], "training", "strict_pre_cutoff_training_candidate", EXPECTED_COUNTS["training"]
    )
    future_fields, future = load_strict_table(
        strict_paths["future"],
        "future",
        "strict_post_cutoff_future_candidate",
        EXPECTED_COUNTS["future"] if future_pairs_path is None else None,
    )
    if set(training) & set(future):
        raise ValueError("Frozen training and future pair tables overlap")

    # ``allowed_pairs`` represents all source versions named in the frozen pair
    # table. ``active_pairs`` restricts model-input structural provenance to v2
    # for training and v3 for future.
    allowed_pairs: dict[str, set[str]] = {"v2": set(), "v3": set()}
    active_pairs: dict[str, dict[str, str]] = {"v2": {}, "v3": {}}
    split_pairs = {"training": training, "future": future}
    for split, table in split_pairs.items():
        active_version = "v2" if split == "training" else "v3"
        for key, row in table.items():
            versions = {value for value in row.get("source_versions", "").split(";") if value}
            if active_version not in versions:
                raise ValueError(f"{split}: {key} lacks expected active source version {active_version}")
            for version in versions:
                if version not in allowed_pairs:
                    raise ValueError(f"{split}: {key} lists unsupported source version {version!r}")
                allowed_pairs[version].add(key)
            active_pairs[active_version][key] = split

    # All permitted source rows are retained only for a SMILES-consistency audit.
    # Strict primary records from the temporally appropriate source version form
    # the actual model-input evidence membership and SMILES resolution set.
    all_source_smiles: dict[str, Counter[str]] = defaultdict(Counter)
    all_source_record_counts: Counter[str] = Counter()
    active_smiles: dict[str, Counter[str]] = defaultdict(Counter)
    active_source_records: dict[str, list[dict[str, str]]] = defaultdict(list)
    active_record_counts: Counter[str] = Counter()
    source_audit_counts: dict[str, Counter[str]] = {"v2": Counter(), "v3": Counter()}

    for version, source_path in source_paths.items():
        for row in read_delimited_gz(source_path, "\t"):
            key = canonical_key(row)
            if key not in allowed_pairs[version]:
                continue
            source_audit_counts[version]["all_allowed_pair_rows"] += 1
            smiles = row.get("smiles", "").strip()
            if smiles:
                all_source_smiles[row["inchikey_full"].strip()][smiles] += 1
            all_source_record_counts[key] += 1
            split = active_pairs[version].get(key)
            if split is None or not is_strict_primary(row):
                continue
            active_source_records[split].append(compact_source_row(row, split))
            active_record_counts[key] += 1
            if smiles:
                active_smiles[row["inchikey_full"].strip()][smiles] += 1
            source_audit_counts[version]["active_strict_primary_rows"] += 1

    missing_active = [key for version in active_pairs.values() for key in version if not active_record_counts[key]]
    if missing_active:
        raise ValueError(f"No strict A/B P1 source record for {len(missing_active)} frozen pairs; first: {missing_active[:5]}")
    strict_compounds = {row["inchikey_full"].strip() for table in split_pairs.values() for row in table.values()}
    missing_smiles = sorted(compound for compound in strict_compounds if not active_smiles[compound])
    if missing_smiles:
        raise ValueError(f"No active-source SMILES for {len(missing_smiles)} compounds; first: {missing_smiles[:5]}")

    # Choose the most frequent active-source SMILES; lexical ordering breaks
    # ties. This does not chemically canonicalize or silently erase conflicts.
    representative: dict[str, str] = {}
    for inchikey, counts in active_smiles.items():
        best_count = max(counts.values())
        representative[inchikey] = sorted(smiles for smiles, count in counts.items() if count == best_count)[0]

    train_compounds = {row["inchikey_full"].strip() for row in training.values()}
    future_compounds = {row["inchikey_full"].strip() for row in future.values()}
    train_targets = {row["uniprot_canonical_accession"].strip() for row in training.values()}
    future_targets = {row["uniprot_canonical_accession"].strip() for row in future.values()}
    sequences = read_archived_fasta(archive_fasta)
    target_union = train_targets | future_targets
    missing_sequences = sorted(target for target in target_union if not sequences.get(target))
    if missing_sequences:
        raise ValueError(f"Archived UniProt FASTA lacks {len(missing_sequences)} selected targets; first: {missing_sequences[:5]}")

    safe_prepare_output(output_dir, root / "data" / "interim", overwrite)

    pair_extra_fields = [
        "split", "representative_smiles", "representative_smiles_selection_rule",
        "active_source_primary_record_count", "all_allowed_source_record_count",
        "active_distinct_smiles_count", "all_allowed_distinct_smiles_count",
        "smiles_conflict_in_active_source", "smiles_conflict_in_any_allowed_source",
        "is_compound_unseen_vs_training", "is_target_unseen_vs_training", "is_date_entity_double_cold_candidate",
    ]
    pair_outputs: dict[str, list[dict[str, object]]] = {"training": [], "future": []}
    for split, table in split_pairs.items():
        for key, original in sorted(table.items()):
            inchikey = original["inchikey_full"].strip()
            target = original["uniprot_canonical_accession"].strip()
            compound_new = split == "future" and inchikey not in train_compounds
            target_new = split == "future" and target not in train_targets
            pair_outputs[split].append({
                **original,
                "split": split,
                "representative_smiles": representative[inchikey],
                "representative_smiles_selection_rule": "mode_of_active_source_strict_P1_rows__lexical_tie_break",
                "active_source_primary_record_count": active_record_counts[key],
                "all_allowed_source_record_count": all_source_record_counts[key],
                "active_distinct_smiles_count": len(active_smiles[inchikey]),
                "all_allowed_distinct_smiles_count": len(all_source_smiles[inchikey]),
                "smiles_conflict_in_active_source": len(active_smiles[inchikey]) > 1,
                "smiles_conflict_in_any_allowed_source": len(all_source_smiles[inchikey]) > 1,
                "is_compound_unseen_vs_training": compound_new,
                "is_target_unseen_vs_training": target_new,
                "is_date_entity_double_cold_candidate": compound_new and target_new,
            })

    compound_rows: list[dict[str, object]] = []
    conflict_rows: list[dict[str, object]] = []
    for inchikey in sorted(strict_compounds):
        active_counts = active_smiles[inchikey]
        all_counts = all_source_smiles[inchikey]
        row = {
            "inchikey_full": inchikey,
            "representative_smiles": representative[inchikey],
            "representative_smiles_selection_rule": "mode_of_active_source_strict_P1_rows__lexical_tie_break",
            "active_source_smiles_count": len(active_counts),
            "all_allowed_source_smiles_count": len(all_counts),
            "smiles_conflict_in_active_source": len(active_counts) > 1,
            "smiles_conflict_in_any_allowed_source": len(all_counts) > 1,
            "active_source_smiles_with_counts": values_with_counts(active_counts),
            "all_allowed_source_smiles_with_counts": values_with_counts(all_counts),
            "training_pair_count": sum(row["inchikey_full"].strip() == inchikey for row in training.values()),
            "future_pair_count": sum(row["inchikey_full"].strip() == inchikey for row in future.values()),
            "is_future_compound_unseen_vs_training": inchikey in future_compounds and inchikey not in train_compounds,
        }
        compound_rows.append(row)
        if row["smiles_conflict_in_active_source"] or row["smiles_conflict_in_any_allowed_source"]:
            conflict_rows.append(row)

    target_rows: list[dict[str, object]] = []
    for accession in sorted(target_union):
        sequence = sequences[accession]
        target_rows.append({
            "uniprot_canonical_accession": accession,
            "sequence_length": len(sequence),
            "sequence_md5": hashlib.md5(sequence.encode("ascii")).hexdigest(),
            "training_pair_count": sum(row["uniprot_canonical_accession"].strip() == accession for row in training.values()),
            "future_pair_count": sum(row["uniprot_canonical_accession"].strip() == accession for row in future.values()),
            "is_future_target_unseen_vs_training": accession in future_targets and accession not in train_targets,
            "sequence_origin": "archived_strict_reviewed_human_uniprot_fasta_20260715T083054Z",
        })

    files = {
        "training_pairs": output_dir / "strict_temporal_training_pairs.tsv.gz",
        "future_pairs": output_dir / "strict_temporal_future_pairs.tsv.gz",
        "training_evidence": output_dir / "training_active_source_strict_primary_evidence.tsv.gz",
        "future_evidence": output_dir / "future_active_source_strict_primary_evidence.tsv.gz",
        "compounds": output_dir / "compounds.tsv.gz",
        "smiles_conflicts": output_dir / "compound_smiles_conflicts.tsv.gz",
        "targets": output_dir / "targets.tsv.gz",
        "sequences": output_dir / "archived_strict_uniprot_sequences.fasta",
    }
    write_tsv_gz(files["training_pairs"], training_fields + pair_extra_fields, pair_outputs["training"])
    write_tsv_gz(files["future_pairs"], future_fields + pair_extra_fields, pair_outputs["future"])
    evidence_fields = list(active_source_records["training"][0])
    write_tsv_gz(files["training_evidence"], evidence_fields, sorted(active_source_records["training"], key=lambda x: (x["canonical_pair_key"], x["ref_id_type"], x["ref_id"], x["activity_type"], x["smiles"])))
    write_tsv_gz(files["future_evidence"], evidence_fields, sorted(active_source_records["future"], key=lambda x: (x["canonical_pair_key"], x["ref_id_type"], x["ref_id"], x["activity_type"], x["smiles"])))
    compound_fields = list(compound_rows[0])
    write_tsv_gz(files["compounds"], compound_fields, compound_rows)
    write_tsv_gz(files["smiles_conflicts"], compound_fields, conflict_rows)
    write_tsv_gz(files["targets"], list(target_rows[0]), target_rows)
    with files["sequences"].open("wt", encoding="utf-8", newline="\n") as handle:
        for accession in sorted(target_union):
            handle.write(f">{accession}\n{sequences[accession]}\n")

    future_new_compound_pairs = sum(row["is_compound_unseen_vs_training"] for row in pair_outputs["future"])
    future_new_target_pairs = sum(row["is_target_unseen_vs_training"] for row in pair_outputs["future"])
    future_double_cold_pairs = sum(row["is_date_entity_double_cold_candidate"] for row in pair_outputs["future"])
    uniprot_info = json.loads(uniprot_summary.read_text(encoding="utf-8"))
    manifest = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "tag": tag,
        "purpose": "Cold-start input preparation from frozen date-verified strict temporal A/B P1 pair membership only.",
        "label_caveat": "P1 candidates are not final positives; P2 assay/source-paper review remains required.",
        "no_negative_label_policy": "Unrecorded compound-target pairs are not represented or emitted as negatives.",
        "input_rules": {
            "training": "exactly the frozen strict pre-cutoff table; structural/evidence source restricted to NPASS v2",
            "future": "exactly the frozen strict post-cutoff table; structural/evidence source restricted to NPASS v3",
            "representative_smiles": "most frequent active-source strict A/B P1 SMILES; lexical tie-break; conflicts retained in a separate table",
            "sequences": "local archived strict reviewed-human UniProt FASTA only; no network retrieval",
        },
        "input_sha256": {str(path): sha256(path) for path in required},
        "frozen_pair_membership_sha256": {
            "training": membership_digest(training), "future": membership_digest(future),
        },
        "archived_uniprot": {
            "fasta": str(archive_fasta),
            "release": uniprot_info.get("fasta_headers", {}).get("X-UniProt-Release", "unknown"),
            "release_date": uniprot_info.get("fasta_headers", {}).get("X-UniProt-Release-Date", "unknown"),
            "retrieval_utc": uniprot_info.get("queried_at", "unknown"),
        },
        "counts": {
            "training_pairs": len(training), "future_pairs": len(future),
            "training_active_source_evidence_rows": len(active_source_records["training"]),
            "future_active_source_evidence_rows": len(active_source_records["future"]),
            "compound_union": len(strict_compounds), "target_union": len(target_union),
            "training_compounds": len(train_compounds), "future_compounds": len(future_compounds),
            "training_targets": len(train_targets), "future_targets": len(future_targets),
            "future_compounds_unseen_vs_training": len(future_compounds - train_compounds),
            "future_targets_unseen_vs_training": len(future_targets - train_targets),
            "future_pairs_with_unseen_compound": future_new_compound_pairs,
            "future_pairs_with_unseen_target": future_new_target_pairs,
            "future_date_entity_double_cold_pairs": future_double_cold_pairs,
            "smiles_conflicts_in_active_source": sum(len(active_smiles[key]) > 1 for key in strict_compounds),
            "smiles_conflicts_in_any_allowed_source": sum(len(all_source_smiles[key]) > 1 for key in strict_compounds),
            "archived_sequences_selected": len(target_union),
            "source_audit_counts": {version: dict(counter) for version, counter in source_audit_counts.items()},
        },
        "outputs": {name: str(path) for name, path in files.items()},
        "outputs_sha256": {name: sha256(path) for name, path in files.items()},
        "limitations": [
            "Date-cold flags mean unseen in the frozen training pair table; scaffold and homology isolation are not yet applied.",
            "The archived UniProt release was retrieved in 2026 and is entity-resolution provenance, not proof of a historical sequence annotation at the 2022 cutoff.",
            "SMILES are not chemically canonicalized by this script; all detected alternatives remain auditable.",
        ],
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    result_path = root / "results" / f"cold_start_inputs_{tag}_manifest.json"
    result_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return {"output_dir": str(output_dir), "result_path": str(result_path), "manifest": manifest}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output-dir", type=Path, default=None)
    parser.add_argument("--training-pairs", type=Path, default=None, help="override frozen strict historical pair table")
    parser.add_argument("--future-pairs", type=Path, default=None, help="override frozen strict future pair table after an explicit leakage gate")
    parser.add_argument("--uniprot-fasta", type=Path, default=None, help="archived reviewed-human UniProt FASTA.GZ")
    parser.add_argument("--uniprot-summary", type=Path, default=None, help="receipt JSON for the archived UniProt request")
    parser.add_argument("--tag", default=TAG, help="output tag used only for derived manifests")
    parser.add_argument("--overwrite", action="store_true", help="Replace only an existing directory under data/interim.")
    args = parser.parse_args()
    root = args.project_root.resolve()
    output_dir = (args.output_dir or root / "data/interim" / f"cold_start_inputs_{args.tag}").resolve()
    result = run(
        root,
        output_dir,
        args.overwrite,
        args.training_pairs,
        args.future_pairs,
        args.tag,
        args.uniprot_fasta,
        args.uniprot_summary,
    )
    print(result["result_path"])


if __name__ == "__main__":
    main()
