"""Shared fail-closed utilities for the aggregate-only corrective audit suite."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import math
import re
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Iterator, Mapping, Sequence


PROTOCOL_ID = "npass_strict_ab_asof_cutoff_corrective_successor_v1_20260728"
SUITE_ID = "npass_strict_ab_asof_cutoff_aggregate_audit_suite_v1_20260728"
EXPECTED_PROTOCOL_LOCK_SHA256 = (
    "96befee13ae1d41ad433c8697fac92ccd30fb25e24c3cf1279c6b4b7e040abd9"
)
INCHIKEY_PATTERN = re.compile(r"\b[A-Z]{14}-[A-Z]{10}-[A-Z]\b")
FORBIDDEN_OUTPUT_KEYS = {
    "query_id",
    "pair_key",
    "canonical_pair_key",
    "inchikey_full",
    "query_compound_inchikey_full",
    "uniprot_canonical_accession",
    "target_uniprot_accession",
    "ref_id",
    "pmid",
    "source_document_id",
    "rank_vector",
    "score_vector",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def require_protocol_lock(lock_path: Path) -> dict[str, Any]:
    """Verify the exact pre-result protocol lock and every file it seals."""

    lock_path = lock_path.resolve()
    if not lock_path.is_file():
        raise FileNotFoundError(f"Protocol lock is absent: {lock_path}")
    observed_lock_hash = sha256(lock_path)
    if observed_lock_hash != EXPECTED_PROTOCOL_LOCK_SHA256:
        raise ValueError(
            "Protocol lock hash mismatch; refusing post-lock protocol drift: "
            f"{observed_lock_hash}"
        )
    payload = json.loads(lock_path.read_text(encoding="utf-8"))
    if payload.get("protocol_id") != PROTOCOL_ID:
        raise ValueError("Unexpected protocol_id in protocol lock")
    if payload.get("lock_timing") != "before corrective result generation":
        raise ValueError("Protocol was not certified as locked before result generation")
    if payload.get("legacy_artifacts_mutated") is not False:
        raise ValueError("Protocol lock does not certify preservation of legacy artifacts")
    if payload.get("public_release_authorized") is not False:
        raise ValueError("This suite expects restricted internal execution, not public release")

    successor_root = lock_path.parent.parent
    sealed = payload.get("files_sha256")
    if not isinstance(sealed, dict) or not sealed:
        raise ValueError("Protocol lock lacks its sealed file inventory")
    for relative_name, expected_hash in sorted(sealed.items()):
        candidate = (successor_root / relative_name).resolve()
        try:
            candidate.relative_to(successor_root.resolve())
        except ValueError as error:
            raise ValueError(f"Locked path escapes successor root: {relative_name}") from error
        if not candidate.is_file():
            raise FileNotFoundError(f"Locked protocol file is absent: {relative_name}")
        if sha256(candidate) != expected_hash:
            raise ValueError(f"Locked protocol file drifted: {relative_name}")
    return payload


def require_new_output_dir(path: Path) -> Path:
    path = path.resolve()
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite audit output: {path}")
    return path


@contextmanager
def open_text(path: Path) -> Iterator[Any]:
    if path.suffix.lower() == ".gz":
        with gzip.open(path, "rt", encoding="utf-8-sig", newline="") as handle:
            yield handle
    else:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            yield handle


@contextmanager
def open_dict_reader(path: Path) -> Iterator[csv.DictReader]:
    with open_text(path) as handle:
        first_line = handle.readline()
        if not first_line:
            raise ValueError(f"Input is empty: {path}")
        delimiter = "\t" if first_line.count("\t") >= first_line.count(",") else ","
        handle.seek(0)
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Input lacks a header: {path}")
        yield reader


def choose_field(fieldnames: Sequence[str], candidates: Sequence[str], role: str) -> str:
    for candidate in candidates:
        if candidate in fieldnames:
            return candidate
    raise ValueError(f"{role} lacks any accepted field: {', '.join(candidates)}")


def parse_bool(value: Any) -> bool:
    normalized = str(value).strip().casefold()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no", ""}:
        return False
    raise ValueError(f"Cannot parse boolean value: {value!r}")


def percentile(values: Sequence[float], probability: float) -> float | None:
    if not values:
        return None
    if not 0.0 <= probability <= 1.0:
        raise ValueError("Percentile probability must lie in [0, 1]")
    ordered = sorted(float(value) for value in values)
    position = (len(ordered) - 1) * probability
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def distribution(values: Sequence[float]) -> dict[str, float | int | None]:
    numeric = [float(value) for value in values]
    if not numeric:
        return {
            "n": 0,
            "min": None,
            "q25": None,
            "median": None,
            "q75": None,
            "max": None,
            "mean": None,
        }
    return {
        "n": len(numeric),
        "min": min(numeric),
        "q25": percentile(numeric, 0.25),
        "median": percentile(numeric, 0.50),
        "q75": percentile(numeric, 0.75),
        "max": max(numeric),
        "mean": sum(numeric) / len(numeric),
    }


def gini_nonnegative(values: Sequence[int]) -> float | None:
    ordered = sorted(float(value) for value in values)
    if not ordered:
        return None
    total = sum(ordered)
    if total == 0:
        return 0.0
    n = len(ordered)
    weighted = sum((index + 1) * value for index, value in enumerate(ordered))
    return (2.0 * weighted) / (n * total) - (n + 1.0) / n


def input_descriptor(role: str, path: Path) -> dict[str, Any]:
    path = path.resolve()
    if not path.is_file():
        raise FileNotFoundError(f"Required {role} input is absent: {path}")
    return {
        "role": role,
        "basename": path.name,
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
    }


def _walk_output(value: Any, path: tuple[str, ...] = ()) -> Iterable[tuple[tuple[str, ...], Any]]:
    if isinstance(value, Mapping):
        for key, child in value.items():
            key_string = str(key)
            if key_string in FORBIDDEN_OUTPUT_KEYS or key_string.endswith("_ids"):
                raise ValueError(f"Identifier-bearing output field is forbidden: {'.'.join(path + (key_string,))}")
            yield from _walk_output(child, path + (key_string,))
    elif isinstance(value, (list, tuple)):
        for index, child in enumerate(value):
            yield from _walk_output(child, path + (str(index),))
    else:
        yield path, value


def assert_aggregate_only(payload: Any) -> None:
    """Reject obvious row-level identifiers in a proposed output payload."""

    for field_path, value in _walk_output(payload):
        if isinstance(value, str) and INCHIKEY_PATTERN.search(value):
            raise ValueError(f"Full InChIKey leaked into aggregate output at {'.'.join(field_path)}")


def write_json_new(path: Path, payload: Any) -> None:
    assert_aggregate_only(payload)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="\n") as handle:
        json.dump(payload, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")


def write_tsv_new(path: Path, fieldnames: Sequence[str], rows: Iterable[Mapping[str, Any]]) -> None:
    if any(field in FORBIDDEN_OUTPUT_KEYS or field.endswith("_ids") for field in fieldnames):
        raise ValueError("Identifier-bearing TSV field requested")
    materialized = [dict(row) for row in rows]
    assert_aggregate_only(materialized)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(materialized)


def finalize_manifest(
    *,
    output_dir: Path,
    audit_id: str,
    script_path: Path,
    inputs: Sequence[dict[str, Any]],
    output_names: Sequence[str],
    extra: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    outputs = []
    for name in output_names:
        path = output_dir / name
        outputs.append({"basename": name, "bytes": path.stat().st_size, "sha256": sha256(path)})
    payload: dict[str, Any] = {
        "audit_id": audit_id,
        "audit_suite_id": SUITE_ID,
        "protocol_id": PROTOCOL_ID,
        "created_at_utc": utc_now(),
        "aggregate_only": True,
        "identifier_bearing_output": False,
        "alternative_rank_ledger_written": False,
        "script": {"basename": script_path.name, "sha256": sha256(script_path)},
        "protocol_lock_sha256": EXPECTED_PROTOCOL_LOCK_SHA256,
        "inputs": list(inputs),
        "outputs": outputs,
    }
    if extra:
        payload.update(dict(extra))
    return payload


