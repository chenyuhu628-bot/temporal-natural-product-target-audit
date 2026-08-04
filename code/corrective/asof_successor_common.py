"""Shared fail-closed helpers for the as-of-cutoff corrective successor.

This module contains only execution contracts and deterministic I/O helpers.
It does not know raw NPASS paths and never creates or interprets labels.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
import math
import os
import platform
import sys
from pathlib import Path
from typing import Any, Iterable


PROTOCOL_ID = "npass_strict_ab_asof_cutoff_corrective_successor_v1_20260728"
RUN_ID = "npass_strict_ab_asof_cutoff_author_run_v1_20260728"
RUN_MODE = "author_run_non_independent_corrective_successor"
LEGACY_TIE_SALT = "npass_strict_ab_doublecold_successor_v1_20260719"

STRICT_TIERS = {
    "A_affinity_candidate": 1.0,
    "B_quantitative_functional_candidate": 0.7,
}
HISTORY_DECISION = "strict_pre_cutoff_training_candidate"
ENDPOINT_DECISION = "strict_post_cutoff_future_candidate"
UNLABELED_POLICY = "unlabeled_not_negative"

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

EXPECTED_HISTORY_PAIRS = 4_990
EXPECTED_QUERIES = 222
EXPECTED_ENDPOINT_RELATIONS = 358
EXPECTED_ENDPOINT_TARGETS = 156
EXPECTED_CANDIDATE_TARGETS = 4_123
EXPECTED_COMPLETE_RANK_ROWS = 3_658_128
PAIR_NEIGHBOR_TOP_K = 100
MORGAN_RADIUS = 2
MORGAN_BITS = 2_048
BOOTSTRAP_REPLICATES = 10_000
BOOTSTRAP_SEED = 20_260_719

ROOT = Path(__file__).resolve().parents[1]
WORKSPACE = ROOT.parent


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def require_fields(fields: Iterable[str], required: set[str], label: str) -> None:
    missing = sorted(required.difference(fields))
    if missing:
        raise ValueError(f"{label} lacks required fields: {missing}")


def require_unique(rows: list[dict[str, str]], fields: tuple[str, ...], label: str) -> None:
    seen: set[tuple[str, ...]] = set()
    for row in rows:
        key = tuple(row.get(field, "") for field in fields)
        if not all(key):
            raise ValueError(f"{label} has an empty key for {fields}")
        if key in seen:
            raise ValueError(f"{label} has a duplicate key for {fields}: {key}")
        seen.add(key)


def parse_bool(value: str, label: str) -> bool:
    normalized = value.strip().casefold()
    if normalized in {"1", "true", "yes"}:
        return True
    if normalized in {"0", "false", "no"}:
        return False
    raise ValueError(f"{label} is not a boolean: {value!r}")


def parse_positive_int(value: str, label: str) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ValueError(f"{label} is not an integer: {value!r}") from exc
    if parsed < 1:
        raise ValueError(f"{label} must be positive: {parsed}")
    return parsed


def assert_isolated_input(path: Path) -> Path:
    resolved = path.resolve()
    if not resolved.is_file() and not resolved.is_dir():
        raise FileNotFoundError(f"Required successor input is absent: {resolved}")
    for blocked in (WORKSPACE / "data", WORKSPACE / "results", WORKSPACE / "manifests"):
        try:
            resolved.relative_to(blocked.resolve())
        except ValueError:
            continue
        raise ValueError(f"Corrective successor input must not be read directly from legacy tree: {blocked}")
    return resolved


def assert_new_output_dir(path: Path) -> Path:
    resolved = path.resolve()
    if resolved.exists():
        raise FileExistsError(f"Refusing to overwrite corrective output directory: {resolved}")
    for blocked in (WORKSPACE / "data", WORKSPACE / "results", WORKSPACE / "manifests"):
        try:
            resolved.relative_to(blocked.resolve())
        except ValueError:
            continue
        raise ValueError(f"Corrective output may not be written inside legacy tree: {blocked}")
    return resolved


def read_tsv_gz(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"Missing TSV header: {path}")
        return list(reader.fieldnames), list(reader)


def write_tsv_gz(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite corrective output: {path}")
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(text, fieldnames=fields, delimiter="\t", extrasaction="raise")
                writer.writeheader()
                writer.writerows(rows)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite corrective output: {path}")
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    accession = ""
    pieces: list[str] = []
    for raw in path.read_text(encoding="ascii").splitlines():
        if raw.startswith(">"):
            if accession:
                if accession in records:
                    raise ValueError(f"Duplicate FASTA accession: {accession}")
                records[accession] = "".join(pieces)
            accession = raw[1:].split("|")[0].strip()
            pieces = []
        else:
            pieces.append(raw.strip())
    if accession:
        if accession in records:
            raise ValueError(f"Duplicate FASTA accession: {accession}")
        records[accession] = "".join(pieces)
    if not records or any(not key or not value for key, value in records.items()):
        raise ValueError(f"FASTA contains an empty identifier or sequence: {path}")
    return records


def load_receipt(path: Path) -> dict[str, Any]:
    receipt = json.loads(path.read_text(encoding="utf-8"))
    require(receipt.get("protocol_id") == PROTOCOL_ID, "Corrective receipt protocol ID mismatch")
    require(receipt.get("run_id") == RUN_ID, "Corrective receipt run ID mismatch")
    require(receipt.get("execution_mode") == RUN_MODE, "Corrective receipt execution mode mismatch")
    authorized = (
        receipt.get("project_lead_authorized_internal_use") is True
        or receipt.get("execution_authorized") is True
    )
    require(authorized, "Corrective receipt lacks project-lead internal-use authorization")
    require(
        receipt.get("public_release_authorized") is not True,
        "Execution receipt must not assert public-release authorization",
    )
    for label in ("protocol_lock", "code_lock", "spec"):
        item = receipt.get(label)
        require(isinstance(item, dict), f"Corrective receipt lacks {label}")
        recorded_path = Path(str(item.get("path", "")))
        require(recorded_path.is_file(), f"Corrective receipt {label} file is absent")
        require(item.get("sha256") == sha256(recorded_path), f"Corrective receipt {label} hash mismatch")
    receipt["_path"] = str(path)
    return receipt


def load_input_manifest(
    path: Path,
    inputs: dict[str, Path],
    expected_kind: str,
) -> dict[str, Any]:
    manifest = json.loads(path.read_text(encoding="utf-8"))
    require(manifest.get("protocol_id") == PROTOCOL_ID, "Corrective input manifest protocol ID mismatch")
    require(manifest.get("run_id") == RUN_ID, "Corrective input manifest run ID mismatch")
    require(manifest.get("execution_mode") == RUN_MODE, "Corrective input manifest mode mismatch")
    require(manifest.get("input_kind") == expected_kind, f"Input manifest kind is not {expected_kind!r}")
    require(
        manifest.get("project_lead_authorized_internal_use") is True,
        "Corrective input manifest lacks internal-use authorization",
    )
    require(manifest.get("legacy_outer_or_result_input") is False, "Manifest permits a legacy result input")
    expected_endpoint_flag = expected_kind == "corrective_evaluation_endpoint"
    require(
        manifest.get("endpoint_file_included") is expected_endpoint_flag,
        "Corrective input manifest endpoint-inclusion flag is inconsistent with its stage",
    )
    declared = manifest.get("file_sha256")
    require(isinstance(declared, dict), "Corrective input manifest lacks file_sha256")
    for name, source in inputs.items():
        require(declared.get(name) == sha256(source), f"Input manifest hash mismatch for {name}")
    manifest["_path"] = str(path)
    return manifest


def validate_historical_pairs(fields: list[str], rows: list[dict[str, str]]) -> None:
    require_fields(
        fields,
        {
            "canonical_pair_key",
            "inchikey_full",
            "uniprot_canonical_accession",
            "best_strict_evidence_tier",
            "eligible_pre_cutoff_v2_row_count",
            "decision",
            "unrecorded_pair_policy",
        },
        "corrected historical pairs",
    )
    require(len(rows) == EXPECTED_HISTORY_PAIRS, f"Historical pair count is not {EXPECTED_HISTORY_PAIRS}")
    require_unique(rows, ("canonical_pair_key",), "corrected historical pairs")
    require_unique(rows, ("inchikey_full", "uniprot_canonical_accession"), "corrected historical pairs")
    for row in rows:
        require(row["best_strict_evidence_tier"] in STRICT_TIERS, "Historical pair has a non-strict tier")
        require(row["decision"] == HISTORY_DECISION, "Historical pair has an invalid decision")
        require(row["unrecorded_pair_policy"] == UNLABELED_POLICY, "Historical pair creates a negative label")
        parse_positive_int(
            row["eligible_pre_cutoff_v2_row_count"],
            f"eligible_pre_cutoff_v2_row_count for {row['canonical_pair_key']}",
        )


def validate_queries(fields: list[str], rows: list[dict[str, str]]) -> None:
    require_fields(fields, {"query_id", "inchikey_full"}, "corrective scoring queries")
    require(len(rows) == EXPECTED_QUERIES, f"Query count is not {EXPECTED_QUERIES}")
    require_unique(rows, ("query_id",), "corrective scoring queries")
    require_unique(rows, ("inchikey_full",), "corrective scoring queries")


def validate_role_structures(
    fields: list[str],
    rows: list[dict[str, str]],
    expected_compounds: set[str],
    role: str,
) -> dict[str, str]:
    require_fields(fields, {"inchikey_full", "representative_smiles"}, f"{role} structures")
    require_unique(rows, ("inchikey_full",), f"{role} structures")
    observed = {row["inchikey_full"] for row in rows}
    require(observed == expected_compounds, f"{role} structure keyset differs from its required compound role")
    require(all(row["representative_smiles"].strip() for row in rows), f"{role} structure contains empty SMILES")
    return {row["inchikey_full"]: row["representative_smiles"] for row in rows}


def validate_targets(fields: list[str], rows: list[dict[str, str]]) -> list[str]:
    require_fields(fields, {"uniprot_canonical_accession"}, "corrective candidate targets")
    require_unique(rows, ("uniprot_canonical_accession",), "corrective candidate targets")
    targets = [row["uniprot_canonical_accession"] for row in rows]
    require(len(targets) == EXPECTED_CANDIDATE_TARGETS, "Candidate target count is not 4,123")
    require(targets == sorted(targets), "Candidate targets are not in ascending deterministic order")
    return targets


def finite_scores(values: Any, label: str) -> None:
    if not all(math.isfinite(float(item)) for item in values):
        raise ValueError(f"{label} contains a non-finite score")


def peak_rss_bytes() -> tuple[int | None, str]:
    """Return process-lifetime peak RSS using only the standard library."""

    if os.name == "nt":
        try:
            import ctypes
            from ctypes import wintypes

            class ProcessMemoryCounters(ctypes.Structure):
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
            get_memory = psapi.GetProcessMemoryInfo
            get_memory.argtypes = [
                wintypes.HANDLE,
                ctypes.POINTER(ProcessMemoryCounters),
                wintypes.DWORD,
            ]
            get_memory.restype = wintypes.BOOL
            counters = ProcessMemoryCounters()
            counters.cb = ctypes.sizeof(counters)
            process = kernel32.GetCurrentProcess()
            ok = get_memory(process, ctypes.byref(counters), counters.cb)
            if ok:
                return int(counters.PeakWorkingSetSize), "windows_GetProcessMemoryInfo_PeakWorkingSetSize"
        except (AttributeError, OSError, TypeError):
            pass
    try:
        import resource

        value = int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        factor = 1 if sys.platform == "darwin" else 1024
        return value * factor, "resource_getrusage_ru_maxrss"
    except (ImportError, OSError, ValueError):
        return None, "not_available"


def environment_receipt() -> dict[str, Any]:
    versions: dict[str, str] = {}
    for module_name in ("numpy", "rdkit", "sklearn"):
        try:
            module = __import__(module_name)
            versions[module_name] = str(getattr(module, "__version__", "unknown"))
        except ImportError:
            versions[module_name] = "not_importable"
    return {
        "python": sys.version,
        "python_executable": sys.executable,
        "platform": platform.platform(),
        "packages": versions,
    }


def code_hashes(paths: Iterable[Path]) -> dict[str, dict[str, str]]:
    result: dict[str, dict[str, str]] = {}
    for path in paths:
        resolved = path.resolve()
        result[resolved.name] = {"path": str(resolved), "sha256": sha256(resolved)}
    return result

