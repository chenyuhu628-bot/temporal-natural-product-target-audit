"""Shared deterministic I/O helpers for the as-of-cutoff corrective run."""

from __future__ import annotations

import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
from typing import Any, Iterable, Iterator


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def read_rows(path: Path, delimiter: str = "\t") -> Iterator[dict[str, str]]:
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Missing table header: {path}")
        yield from reader


def read_table(path: Path, delimiter: str = "\t") -> tuple[list[str], list[dict[str, str]]]:
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if not reader.fieldnames:
            raise ValueError(f"Missing table header: {path}")
        return list(reader.fieldnames), list(reader)


def write_tsv_gz(path: Path, fields: list[str], rows: Iterable[dict[str, object]]) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite output: {path}")
    with path.open("wb") as binary:
        with gzip.GzipFile(filename="", mode="wb", fileobj=binary, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                writer = csv.DictWriter(text, fieldnames=fields, delimiter="\t", extrasaction="raise")
                writer.writeheader()
                writer.writerows(rows)


def write_json(path: Path, payload: dict[str, Any]) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite output: {path}")
    path.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def parse_bool(value: str) -> bool:
    normalized = value.strip().casefold()
    if normalized in {"true", "1", "yes"}:
        return True
    if normalized in {"false", "0", "no"}:
        return False
    raise ValueError(f"Invalid boolean value: {value!r}")


def require_fields(fields: Iterable[str], required: set[str], label: str) -> None:
    missing = sorted(required.difference(fields))
    if missing:
        raise ValueError(f"{label} lacks required fields: {missing}")


def assert_unique(rows: Iterable[dict[str, str]], fields: tuple[str, ...], label: str) -> None:
    observed: set[tuple[str, ...]] = set()
    for row in rows:
        key = tuple(row.get(field, "").strip() for field in fields)
        if not all(key):
            raise ValueError(f"{label} contains an empty key for {fields}")
        if key in observed:
            raise ValueError(f"{label} contains a duplicate key for {fields}: {key}")
        observed.add(key)


def membership_sha256(values: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for value in sorted(set(values)):
        digest.update(value.encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()


