#!/usr/bin/env python3
"""Cluster union target sequences and audit future-target coldness with MMseqs2.

This is an audit utility, not an interaction-label generator.  It clusters a
single union FASTA at a prespecified identity grid and, at each threshold,
checks future targets directly against historical targets using exactly the
same sequence-identity and coverage rule.  Its only target-level outputs are
membership and coldness *candidates*; it never writes negative labels or
changes any raw input.

The default rule is deliberately fixed across all three pressure tests:

* identity: >= 0.30, 0.50, or 0.70 of aligned residues (MMseqs
  ``--seq-id-mode 0``);
* coverage: >= 0.80 of both query and target (MMseqs ``--cov-mode 0``);
* clustering: MMseqs set-cover clustering (``--cluster-mode 0``).

Inputs can point to a later frozen P2/P3 primary table.  The pair files only
need a target-accession column; no positive, negative, or PU label column is
read by this script.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import platform
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

from reproducible_io import deterministic_gzip_text


DEFAULT_THRESHOLDS = (0.30, 0.50, 0.70)
DEFAULT_MIN_COVERAGE = 0.80
DEFAULT_CLUSTER_MODE = 0
DEFAULT_COVERAGE_MODE = 0
DEFAULT_SEQUENCE_ID_MODE = 0
DEFAULT_SEARCH_SENSITIVITY = 7.5
DEFAULT_MAX_SEQS = 1_000_000


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    """Return a SHA-256 digest without loading the whole file into memory."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path, mode: str):
    if "b" in mode:
        raise ValueError("open_text only accepts text modes")
    if path.suffix.casefold() == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def accession_from_header(header: str) -> str:
    """Read a canonical accession from a simple or UniProt-style FASTA header."""
    token = header[1:].strip().split(None, 1)[0] if header[1:].strip() else ""
    parts = token.split("|")
    if len(parts) >= 2 and parts[0].casefold() in {"sp", "tr", "up"} and parts[1].strip():
        return parts[1].strip()
    return token.strip()


def read_fasta(path: Path) -> dict[str, str]:
    """Return accession -> sequence, rejecting ambiguous or malformed FASTA."""
    sequences: dict[str, str] = {}
    accession: str | None = None
    chunks: list[str] = []

    def save_current() -> None:
        nonlocal accession, chunks
        if accession is None:
            return
        # Preserve the submitted residue case in the derived subset FASTAs so
        # both MMseqs calls see the same biological sequence representation.
        sequence = "".join(chunks)
        if not sequence:
            raise ValueError(f"FASTA record {accession!r} has no sequence")
        if not all(character.isalpha() or character == "*" for character in sequence):
            raise ValueError(f"FASTA record {accession!r} contains unsupported non-letter residues")
        if accession in sequences:
            raise ValueError(f"Duplicate FASTA accession: {accession}")
        sequences[accession] = sequence

    with open_text(path, "rt") as handle:
        for number, raw in enumerate(handle, start=1):
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                save_current()
                accession = accession_from_header(line)
                if not accession:
                    raise ValueError(f"Empty FASTA accession at line {number}")
                chunks = []
            elif accession is None:
                raise ValueError(f"Sequence before the first FASTA header at line {number}")
            else:
                chunks.append(line)
    save_current()
    if not sequences:
        raise ValueError(f"No sequences found in {path}")
    return sequences


def read_target_membership(path: Path, target_column: str) -> tuple[set[str], int]:
    """Read unique target accessions from a TSV/TSV.GZ pair table."""
    targets: set[str] = set()
    row_count = 0
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or target_column not in reader.fieldnames:
            raise ValueError(
                f"{path} must be a tab-separated table with column {target_column!r}; "
                f"found {reader.fieldnames!r}"
            )
        for number, row in enumerate(reader, start=2):
            row_count += 1
            accession = (row.get(target_column) or "").strip()
            if not accession:
                raise ValueError(f"Empty {target_column!r} value in {path} at row {number}")
            targets.add(accession)
    if not targets:
        raise ValueError(f"No target accessions found in {path}")
    return targets, row_count


def write_fasta(path: Path, accessions: Iterable[str], sequences: dict[str, str]) -> None:
    with path.open("wt", encoding="utf-8", newline="\n") as handle:
        for accession in sorted(accessions):
            handle.write(f">{accession}\n{sequences[accession]}\n")


def write_tsv_gz(path: Path, fieldnames: list[str], rows: Iterable[dict[str, object]]) -> None:
    with deterministic_gzip_text(path) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def threshold_key(threshold: float) -> str:
    return f"identity_{threshold:.2f}".replace(".", "_")


def require_descendant(path: Path, ancestor: Path) -> None:
    """Fail closed before recursively removing a generated temporary directory."""
    resolved_path = path.resolve()
    resolved_ancestor = ancestor.resolve()
    try:
        resolved_path.relative_to(resolved_ancestor)
    except ValueError as exc:
        raise ValueError(f"Refusing to remove path outside output directory: {resolved_path}") from exc


def remove_generated_temp(path: Path, output_dir: Path) -> None:
    if path.exists():
        require_descendant(path, output_dir)
        shutil.rmtree(path)


def command_line(argv: list[str]) -> str:
    """Produce a Windows-compatible rendering while preserving the argv list too."""
    return subprocess.list2cmdline(argv)


def run_mmseqs(
    *,
    argv: list[str],
    label: str,
    log_directory: Path,
    commands: list[dict[str, object]],
) -> None:
    """Run one command, preserve its stdout/stderr, and fail without ambiguity."""
    stdout_path = log_directory / f"{label}.stdout.log"
    stderr_path = log_directory / f"{label}.stderr.log"
    started_at = utc_now()
    result = subprocess.run(
        argv,
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    stdout_path.write_text(result.stdout, encoding="utf-8")
    stderr_path.write_text(result.stderr, encoding="utf-8")
    record: dict[str, object] = {
        "label": label,
        "started_at_utc": started_at,
        "finished_at_utc": utc_now(),
        "argv": argv,
        "command_line_windows": command_line(argv),
        "returncode": result.returncode,
        "stdout_log": str(stdout_path),
        "stderr_log": str(stderr_path),
    }
    commands.append(record)
    if result.returncode != 0:
        raise RuntimeError(
            f"MMseqs2 command failed ({label}, exit {result.returncode}). "
            f"See {stderr_path} and {stdout_path}."
        )


def mmseqs_version(mmseqs: Path) -> str:
    result = subprocess.run(
        [str(mmseqs), "version"],
        check=False,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    if result.returncode != 0:
        raise RuntimeError(f"Could not run MMseqs2 version command at {mmseqs}: {result.stderr.strip()}")
    version = result.stdout.strip()
    if not version:
        raise RuntimeError(f"MMseqs2 version command returned no version at {mmseqs}")
    return version


def tool_hashes(mmseqs: Path) -> dict[str, dict[str, str]]:
    """Hash the launcher and, for the bundled Windows layout, its executable."""
    artifacts = {"launcher": mmseqs}
    bundled_executable = mmseqs.parent / "bin" / "mmseqs.exe"
    if bundled_executable.is_file():
        artifacts["mmseqs_executable"] = bundled_executable
    return {
        label: {"path": str(path), "sha256": sha256(path)}
        for label, path in artifacts.items()
    }


def locate_cluster_tsv(prefix: Path, temporary_dir: Path) -> Path:
    """Locate MMseqs2's membership table, including the Windows-bundle fallback.

    The pinned Windows bundle can complete ``easy-cluster`` while failing its
    optional final ``result2flat`` step. In that case the authoritative
    ``cluster.tsv`` still exists under the command's unique temporary-run
    directory. The caller copies it into the retained output before cleanup.
    """
    expected = Path(f"{prefix}_cluster.tsv")
    if expected.is_file():
        return expected
    candidates = sorted(prefix.parent.glob(f"{prefix.name}*cluster*.tsv"))
    if len(candidates) == 1:
        return candidates[0]
    temporary_candidates = sorted(
        candidate
        for candidate in temporary_dir.glob("*/cluster.tsv")
        if candidate.is_file()
    )
    if len(temporary_candidates) == 1:
        return temporary_candidates[0]
    raise FileNotFoundError(
        "MMseqs2 did not create a usable cluster membership table. "
        f"Expected {expected}; output candidates: {candidates}; temporary candidates: {temporary_candidates}"
    )


def mmseqs_tab_lines(path: Path) -> Iterable[tuple[int, str]]:
    """Yield tabular MMseqs output despite a Cygwin CR-before-tab quirk.

    The pinned Windows bundle can write ``representative\\r\\tmember\\r\\n``
    rather than conventional ``representative\\tmember\\r\\n``. Reading that
    as universal-newline text turns the two fields into separate lines. Input
    FASTA identifiers cannot contain CR, so normalizing only this byte pattern
    is lossless for the documented MMseqs outputs.
    """
    raw = path.read_bytes().replace(b"\r\t", b"\t").replace(b"\r\n", b"\n")
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise ValueError(f"MMseqs2 output is not UTF-8: {path}") from exc
    for line_number, line in enumerate(text.split("\n"), start=1):
        if line:
            yield line_number, line


def parse_cluster_table(cluster_path: Path, fasta_accessions: set[str]) -> dict[str, str]:
    """Return member -> representative, validating complete one-cluster membership."""
    assignments: dict[str, str] = {}
    representatives: set[str] = set()
    for line_number, line in mmseqs_tab_lines(cluster_path):
        parts = line.split("\t")
        if len(parts) < 2:
            raise ValueError(f"Malformed MMseqs2 cluster line {line_number} in {cluster_path}")
        representative, member = parts[0].strip(), parts[1].strip()
        if representative not in fasta_accessions or member not in fasta_accessions:
            raise ValueError(
                "MMseqs2 returned an identifier absent from the input FASTA at "
                f"{cluster_path}:{line_number}: {representative!r}, {member!r}"
            )
        representatives.add(representative)
        earlier = assignments.get(member)
        if earlier is not None and earlier != representative:
            raise ValueError(
                f"Target {member!r} belongs to two MMseqs2 clusters: {earlier!r} and {representative!r}"
            )
        assignments[member] = representative

    # MMseqs output normally includes a representative's self row.  Add one
    # only when it is absent, so singleton clusters remain explicit.
    for representative in representatives:
        earlier = assignments.get(representative)
        if earlier is None:
            assignments[representative] = representative
        elif earlier != representative:
            raise ValueError(
                f"MMseqs2 representative {representative!r} is assigned under {earlier!r}; cannot audit safely"
            )

    missing = fasta_accessions.difference(assignments)
    if missing:
        preview = ", ".join(sorted(missing)[:20])
        raise ValueError(f"MMseqs2 cluster table omitted {len(missing)} FASTA target(s): {preview}")
    extra = set(assignments).difference(fasta_accessions)
    if extra:
        raise ValueError(f"Unexpected MMseqs2 cluster members: {sorted(extra)[:20]}")
    return assignments


def parse_direct_hits(
    path: Path,
    future_targets: set[str],
    historical_targets: set[str],
) -> dict[str, list[dict[str, object]]]:
    """Parse the fixed MMseqs2 direct-search output and retain every detected hit."""
    hits: dict[str, list[dict[str, object]]] = defaultdict(list)
    for line_number, line in mmseqs_tab_lines(path):
        parts = line.split("\t")
        if len(parts) != 8:
            raise ValueError(
                f"Expected eight columns from fixed MMseqs2 direct search at {path}:{line_number}, got {len(parts)}"
            )
        query, target, pident, alignment_length, query_coverage, target_coverage, evalue, bits = parts
        if query not in future_targets or target not in historical_targets:
            raise ValueError(
                f"Unexpected direct-search identifier at {path}:{line_number}: {query!r} -> {target!r}"
            )
        try:
            hits[query].append(
                {
                    "historical_target_accession": target,
                    "pident": float(pident),
                    "alignment_length": int(float(alignment_length)),
                    "query_coverage": float(query_coverage),
                    "target_coverage": float(target_coverage),
                    "evalue": evalue,
                    "bits": float(bits),
                }
            )
        except ValueError as exc:
            raise ValueError(f"Could not parse direct-search values at {path}:{line_number}") from exc
    return hits


def output_hashes(output_dir: Path) -> dict[str, str]:
    """Hash all retained outputs except the manifest that contains this collection."""
    hashes: dict[str, str] = {}
    for path in sorted(item for item in output_dir.rglob("*") if item.is_file()):
        relative = path.relative_to(output_dir).as_posix()
        if relative == "run_manifest.json":
            continue
        hashes[relative] = sha256(path)
    return hashes


def target_coldness_status(
    accession: str,
    historical_targets: set[str],
    direct_hits: dict[str, list[dict[str, object]]],
    historical_members_in_cluster: list[str],
) -> str:
    """Return a descriptive audit status, never a biological interaction label."""
    if accession in historical_targets:
        return "future_target_accession_seen_in_historical_membership"
    if direct_hits.get(accession):
        return "future_target_direct_historical_alignment_detected"
    if historical_members_in_cluster:
        return "future_target_cluster_shares_historical_membership"
    return "future_target_homology_cold_candidate"


def make_assignments_and_audit(
    *,
    assignments: dict[str, str],
    historical_targets: set[str],
    future_targets: set[str],
    direct_hits: dict[str, list[dict[str, object]]],
) -> tuple[list[dict[str, object]], list[dict[str, object]], dict[str, object]]:
    members_by_representative: dict[str, list[str]] = defaultdict(list)
    for member, representative in assignments.items():
        members_by_representative[representative].append(member)
    for members in members_by_representative.values():
        members.sort()

    assignment_rows: list[dict[str, object]] = []
    future_audit_rows: list[dict[str, object]] = []
    for accession in sorted(assignments):
        representative = assignments[accession]
        members = members_by_representative[representative]
        historical_members = [member for member in members if member in historical_targets]
        future_members = [member for member in members if member in future_targets]
        direct_hit_rows = direct_hits.get(accession, [])
        assignment_rows.append(
            {
                "uniprot_canonical_accession": accession,
                "mmseqs_cluster_representative": representative,
                "cluster_size": len(members),
                "is_cluster_representative": accession == representative,
                "historical_target_present": accession in historical_targets,
                "future_target_present": accession in future_targets,
                "historical_target_count_in_cluster": len(historical_members),
                "future_target_count_in_cluster": len(future_members),
                "cluster_has_historical_target": bool(historical_members),
            }
        )
        if accession not in future_targets:
            continue
        status = target_coldness_status(accession, historical_targets, direct_hits, historical_members)
        best_hit = max(direct_hit_rows, key=lambda item: float(item["pident"]), default=None)
        future_audit_rows.append(
            {
                "uniprot_canonical_accession": accession,
                "historical_target_present": accession in historical_targets,
                "mmseqs_cluster_representative": representative,
                "cluster_size": len(members),
                "historical_target_count_in_cluster": len(historical_members),
                "historical_targets_in_cluster": "|".join(historical_members),
                "direct_historical_alignment_hit_count": len(direct_hit_rows),
                "direct_historical_targets_detected": "|".join(
                    sorted({str(item["historical_target_accession"]) for item in direct_hit_rows})
                ),
                "maximum_detected_direct_pident": "" if best_hit is None else best_hit["pident"],
                "maximum_detected_direct_query_coverage": "" if best_hit is None else best_hit["query_coverage"],
                "maximum_detected_direct_target_coverage": "" if best_hit is None else best_hit["target_coverage"],
                "future_target_coldness_status": status,
                "is_future_target_homology_cold_candidate": status == "future_target_homology_cold_candidate",
            }
        )

    status_counts = Counter(str(row["future_target_coldness_status"]) for row in future_audit_rows)
    cluster_sizes = [len(members) for members in members_by_representative.values()]
    summary: dict[str, object] = {
        "union_target_count": len(assignments),
        "historical_target_count": len(historical_targets),
        "future_target_count": len(future_targets),
        "cluster_count": len(members_by_representative),
        "cluster_size": {
            "minimum": min(cluster_sizes),
            "median": sorted(cluster_sizes)[len(cluster_sizes) // 2],
            "maximum": max(cluster_sizes),
        },
        "future_target_coldness_status_counts": dict(sorted(status_counts.items())),
        "future_target_homology_cold_candidate_count": sum(
            bool(row["is_future_target_homology_cold_candidate"]) for row in future_audit_rows
        ),
    }
    return assignment_rows, future_audit_rows, summary


def validate_arguments(args: argparse.Namespace) -> tuple[Path, Path, Path, Path, Path]:
    fasta = args.fasta.resolve()
    historical_pairs = args.historical_pairs.resolve()
    future_pairs = args.future_pairs.resolve()
    output_dir = args.output_dir.resolve()
    mmseqs = args.mmseqs.resolve()
    for path in (fasta, historical_pairs, future_pairs, mmseqs):
        if not path.is_file():
            raise FileNotFoundError(path)
    if output_dir.exists():
        raise FileExistsError(
            f"Output directory already exists: {output_dir}. Use a new run directory to preserve prior audit outputs."
        )
    if not 0 < args.min_coverage <= 1:
        raise ValueError("--min-coverage must be in (0, 1]")
    if args.threads < 1:
        raise ValueError("--threads must be >= 1")
    if args.search_sensitivity < 1:
        raise ValueError("--search-sensitivity must be >= 1")
    if args.max_seqs < 1:
        raise ValueError("--max-seqs must be >= 1")
    thresholds = tuple(float(value) for value in args.thresholds)
    if not thresholds or any(not 0 < threshold <= 1 for threshold in thresholds):
        raise ValueError("Every --thresholds value must be in (0, 1]")
    if len(set(thresholds)) != len(thresholds):
        raise ValueError("--thresholds must not contain duplicates")
    if tuple(sorted(thresholds)) != thresholds:
        raise ValueError("--thresholds must be supplied in ascending order")
    args.thresholds = thresholds
    return fasta, historical_pairs, future_pairs, output_dir, mmseqs


def run(args: argparse.Namespace) -> Path:
    fasta, historical_pairs, future_pairs, output_dir, mmseqs = validate_arguments(args)
    sequences = read_fasta(fasta)
    historical_targets, historical_row_count = read_target_membership(historical_pairs, args.target_column)
    future_targets, future_row_count = read_target_membership(future_pairs, args.target_column)
    fasta_targets = set(sequences)
    expected_union_targets = historical_targets | future_targets
    absent_from_fasta = expected_union_targets.difference(fasta_targets)
    if absent_from_fasta:
        raise ValueError(
            f"{len(absent_from_fasta)} target accession(s) in pair tables are absent from the union FASTA: "
            f"{', '.join(sorted(absent_from_fasta)[:20])}"
        )
    extra_fasta_targets = fasta_targets.difference(expected_union_targets)
    if extra_fasta_targets and not args.allow_extra_fasta_targets:
        raise ValueError(
            f"The FASTA contains {len(extra_fasta_targets)} target(s) outside historical/future membership. "
            "Pass --allow-extra-fasta-targets only when that inclusion is intentional: "
            f"{', '.join(sorted(extra_fasta_targets)[:20])}"
        )

    output_dir.mkdir(parents=True)
    input_dir = output_dir / "derived_target_subsets"
    input_dir.mkdir()
    logs_dir = output_dir / "logs"
    logs_dir.mkdir()
    write_fasta(input_dir / "historical_targets.fasta", historical_targets, sequences)
    write_fasta(input_dir / "future_targets.fasta", future_targets, sequences)

    manifest: dict[str, object] = {
        "schema_version": "1.0",
        "created_at_utc": utc_now(),
        "purpose": (
            "Temporal target-membership and homology audit only. No interaction positive, negative, "
            "or unlabeled labels are generated by this run."
        ),
        "status": "running",
        "tool": {
            "mmseqs_launcher": str(mmseqs),
            "mmseqs_version": mmseqs_version(mmseqs),
            "artifacts": tool_hashes(mmseqs),
        },
        "parameters": {
            "identity_thresholds": list(args.thresholds),
            "minimum_coverage": args.min_coverage,
            "coverage_mode": DEFAULT_COVERAGE_MODE,
            "coverage_definition": "at least the fixed coverage fraction of both query and target",
            "sequence_identity_mode": DEFAULT_SEQUENCE_ID_MODE,
            "sequence_identity_definition": "identical residues divided by alignment length",
            "cluster_mode": DEFAULT_CLUSTER_MODE,
            "cluster_definition": "MMseqs2 set-cover clustering",
            "search_sensitivity": args.search_sensitivity,
            "max_sequences_per_query": args.max_seqs,
            "threads": args.threads,
            "keep_mmseqs_temp": args.keep_mmseqs_temp,
        },
        "inputs": {
            "union_target_fasta": {
                "path": str(fasta),
                "sha256": sha256(fasta),
                "sequence_count": len(sequences),
                "expected_historical_future_union_count": len(expected_union_targets),
                "extra_fasta_target_count": len(extra_fasta_targets),
                "allow_extra_fasta_targets": args.allow_extra_fasta_targets,
            },
            "historical_pairs": {
                "path": str(historical_pairs),
                "sha256": sha256(historical_pairs),
                "row_count": historical_row_count,
                "unique_target_count": len(historical_targets),
            },
            "future_pairs": {
                "path": str(future_pairs),
                "sha256": sha256(future_pairs),
                "row_count": future_row_count,
                "unique_target_count": len(future_targets),
            },
            "target_column": args.target_column,
            "derived_subset_fastas": {},
        },
        "environment": {"python": sys.version, "platform": platform.platform()},
        "commands": [],
        "threshold_runs": [],
        "limitations": [
            "Cluster membership and detected alignments are target-level audit data, not interaction labels.",
            "A direct-search non-hit is reported as no detected hit under this frozen MMseqs2 run, not proof that no biological homology exists.",
            "A future target is called a homology-cold candidate only when it is accession-unseen, its union cluster has no historical target, and no direct historical alignment is detected under the recorded rule.",
            "A final double-cold evaluation also requires the separate chemical-scaffold isolation and evidence/temporal eligibility protocol.",
        ],
    }
    subset_fastas = [input_dir / "historical_targets.fasta", input_dir / "future_targets.fasta"]
    manifest["inputs"]["derived_subset_fastas"] = {
        path.name: {"path": str(path), "sha256": sha256(path)} for path in subset_fastas
    }

    try:
        for threshold in args.thresholds:
            label = threshold_key(threshold)
            threshold_dir = output_dir / label
            threshold_dir.mkdir()
            cluster_prefix = threshold_dir / "union_cluster"
            cluster_tmp = threshold_dir / "union_cluster_tmp"
            direct_hits_path = threshold_dir / "future_vs_historical_direct_hits.tsv"
            search_tmp = threshold_dir / "direct_search_tmp"
            commands: list[dict[str, object]] = manifest["commands"]  # type: ignore[assignment]

            cluster_command = [
                str(mmseqs),
                "easy-cluster",
                str(fasta),
                str(cluster_prefix),
                str(cluster_tmp),
                "--min-seq-id",
                f"{threshold:.6f}",
                "-c",
                f"{args.min_coverage:.6f}",
                "--cov-mode",
                str(DEFAULT_COVERAGE_MODE),
                "--seq-id-mode",
                str(DEFAULT_SEQUENCE_ID_MODE),
                "--cluster-mode",
                str(DEFAULT_CLUSTER_MODE),
                "--filter-hits",
                "1",
                "--threads",
                str(args.threads),
            ]
            run_mmseqs(
                argv=cluster_command,
                label=f"{label}_easy_cluster",
                log_directory=logs_dir,
                commands=commands,
            )
            located_cluster_table = locate_cluster_tsv(cluster_prefix, cluster_tmp)
            # Keep the exact source table even when the Windows bundle writes
            # it only beneath its temporary run directory.
            cluster_table = threshold_dir / "union_cluster_membership.tsv"
            shutil.copyfile(located_cluster_table, cluster_table)
            assignments = parse_cluster_table(cluster_table, fasta_targets)

            direct_search_command = [
                str(mmseqs),
                "easy-search",
                str(input_dir / "future_targets.fasta"),
                str(input_dir / "historical_targets.fasta"),
                str(direct_hits_path),
                str(search_tmp),
                "--min-seq-id",
                f"{threshold:.6f}",
                "-c",
                f"{args.min_coverage:.6f}",
                "--cov-mode",
                str(DEFAULT_COVERAGE_MODE),
                "--seq-id-mode",
                str(DEFAULT_SEQUENCE_ID_MODE),
                "--filter-hits",
                "1",
                "-s",
                f"{args.search_sensitivity:.6f}",
                "--max-seqs",
                str(args.max_seqs),
                "--format-output",
                "query,target,pident,alnlen,qcov,tcov,evalue,bits",
                "--threads",
                str(args.threads),
            ]
            run_mmseqs(
                argv=direct_search_command,
                label=f"{label}_easy_search",
                log_directory=logs_dir,
                commands=commands,
            )
            direct_hits = parse_direct_hits(direct_hits_path, future_targets, historical_targets)
            assignment_rows, future_audit_rows, summary = make_assignments_and_audit(
                assignments=assignments,
                historical_targets=historical_targets,
                future_targets=future_targets,
                direct_hits=direct_hits,
            )
            assignment_path = threshold_dir / "cluster_assignments.tsv.gz"
            future_audit_path = threshold_dir / "future_target_coldness_audit.tsv.gz"
            write_tsv_gz(
                assignment_path,
                [
                    "uniprot_canonical_accession",
                    "mmseqs_cluster_representative",
                    "cluster_size",
                    "is_cluster_representative",
                    "historical_target_present",
                    "future_target_present",
                    "historical_target_count_in_cluster",
                    "future_target_count_in_cluster",
                    "cluster_has_historical_target",
                ],
                assignment_rows,
            )
            write_tsv_gz(
                future_audit_path,
                [
                    "uniprot_canonical_accession",
                    "historical_target_present",
                    "mmseqs_cluster_representative",
                    "cluster_size",
                    "historical_target_count_in_cluster",
                    "historical_targets_in_cluster",
                    "direct_historical_alignment_hit_count",
                    "direct_historical_targets_detected",
                    "maximum_detected_direct_pident",
                    "maximum_detected_direct_query_coverage",
                    "maximum_detected_direct_target_coverage",
                    "future_target_coldness_status",
                    "is_future_target_homology_cold_candidate",
                ],
                future_audit_rows,
            )
            summary.update(
                {
                    "threshold": threshold,
                    "cluster_table_source_before_retention": str(located_cluster_table),
                    "cluster_table_retrieval": (
                        "standard_easy_cluster_output"
                        if located_cluster_table == Path(f"{cluster_prefix}_cluster.tsv")
                        else "mmseqs_temporary_run_fallback"
                    ),
                    "cluster_table": str(cluster_table),
                    "cluster_table_sha256": sha256(cluster_table),
                    "direct_hits": str(direct_hits_path),
                    "direct_hits_sha256": sha256(direct_hits_path),
                    "cluster_assignments": str(assignment_path),
                    "cluster_assignments_sha256": sha256(assignment_path),
                    "future_target_coldness_audit": str(future_audit_path),
                    "future_target_coldness_audit_sha256": sha256(future_audit_path),
                }
            )
            threshold_summary_path = threshold_dir / "summary.json"
            threshold_summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
            manifest["threshold_runs"].append(summary)  # type: ignore[union-attr]
            if not args.keep_mmseqs_temp:
                remove_generated_temp(cluster_tmp, output_dir)
                remove_generated_temp(search_tmp, output_dir)

        manifest["status"] = "complete"
        manifest["finished_at_utc"] = utc_now()
        manifest["output_sha256"] = output_hashes(output_dir)
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["failed_at_utc"] = utc_now()
        manifest["error"] = f"{type(exc).__name__}: {exc}"
        (output_dir / "run_manifest.json").write_text(
            json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
        )
        raise

    manifest_path = output_dir / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest_path


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", required=True, type=Path, help="Union target FASTA: exactly one sequence per target accession")
    parser.add_argument("--historical-pairs", required=True, type=Path, help="Historical pair table (TSV or TSV.GZ)")
    parser.add_argument("--future-pairs", required=True, type=Path, help="Future pair table (TSV or TSV.GZ)")
    parser.add_argument("--output-dir", required=True, type=Path, help="New, non-existing directory for this immutable audit run")
    parser.add_argument(
        "--mmseqs",
        type=Path,
        default=root / "tools/mmseqs2/18-8cc5c/mmseqs/mmseqs.bat",
        help="MMseqs2 launcher; default is the project's pinned Windows bundle",
    )
    parser.add_argument("--target-column", default="uniprot_canonical_accession")
    parser.add_argument(
        "--allow-extra-fasta-targets",
        action="store_true",
        help="Permit union FASTA accessions absent from both pair tables; disabled by default to prevent hidden clustering inputs",
    )
    parser.add_argument(
        "--thresholds",
        nargs="+",
        type=float,
        default=list(DEFAULT_THRESHOLDS),
        help="Ascending identity grid; default: 0.30 0.50 0.70",
    )
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=DEFAULT_MIN_COVERAGE,
        help="One fixed symmetric coverage cutoff shared by every threshold; default: 0.80",
    )
    parser.add_argument("--threads", type=int, default=8, help="MMseqs2 worker threads; default: 8")
    parser.add_argument(
        "--search-sensitivity",
        type=float,
        default=DEFAULT_SEARCH_SENSITIVITY,
        help="MMseqs2 direct-search sensitivity; default: 7.5",
    )
    parser.add_argument(
        "--max-seqs",
        type=int,
        default=DEFAULT_MAX_SEQS,
        help="Maximum historical matches retained per future query; default: 1000000",
    )
    parser.add_argument(
        "--keep-mmseqs-temp",
        action="store_true",
        help="Keep MMseqs2 temporary directories inside the new output directory for troubleshooting",
    )
    args = parser.parse_args()
    manifest = run(args)
    print(manifest)


if __name__ == "__main__":
    main()
