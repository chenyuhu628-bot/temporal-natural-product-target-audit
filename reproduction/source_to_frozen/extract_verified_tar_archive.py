#!/usr/bin/env python3
"""Verify a tar archive hash and extract it without overwriting any directory.

This is used for the frozen ChEMBL SQLite snapshot. It rejects path traversal,
links, and special files before extraction, leaves raw archives unchanged, and
records the derived extraction manifest separately. Every extracted SQLite file
is also hashed before the derived directory is committed, allowing later
read-only audits to detect post-extraction replacement or truncation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tarfile
from datetime import datetime, timezone
from pathlib import Path, PurePosixPath


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_relative_member(name: str) -> Path:
    normalised = name.replace("\\", "/")
    member = PurePosixPath(normalised)
    if member.is_absolute() or ".." in member.parts:
        raise ValueError(f"Unsafe archive member path: {name!r}")
    parts = [part for part in member.parts if part not in {"", "."}]
    if not parts:
        raise ValueError(f"Empty archive member path: {name!r}")
    return Path(*parts)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--destination", required=True, type=Path)
    parser.add_argument("--expected-sha256", required=True)
    parser.add_argument("--manifest", required=True, type=Path)
    args = parser.parse_args()

    archive = args.archive.resolve()
    destination = args.destination.resolve()
    manifest = args.manifest.resolve()
    expected = args.expected_sha256.strip().casefold()
    if not archive.is_file():
        raise FileNotFoundError(archive)
    if destination.exists():
        raise FileExistsError(f"Refusing to overwrite extraction destination: {destination}")
    if manifest.exists():
        raise FileExistsError(f"Refusing to overwrite extraction manifest: {manifest}")
    actual = sha256(archive)
    if actual.casefold() != expected:
        raise RuntimeError(f"Archive SHA-256 mismatch: got {actual}, expected {expected}")

    temporary = destination.with_name(destination.name + ".extracting")
    if temporary.exists():
        raise FileExistsError(f"Refusing to reuse existing temporary extraction directory: {temporary}")
    temporary.mkdir(parents=True)
    extracted_files: list[dict[str, int | str]] = []
    with tarfile.open(archive, "r:gz") as bundle:
        members = bundle.getmembers()
        for item in members:
            relative = safe_relative_member(item.name)
            if item.issym() or item.islnk() or item.isdev():
                raise ValueError(f"Refusing link or special archive member: {item.name!r}")
            target = (temporary / relative).resolve()
            try:
                target.relative_to(temporary.resolve())
            except ValueError as exc:
                raise ValueError(f"Archive member escapes extraction directory: {item.name!r}") from exc
            if item.isdir():
                target.mkdir(parents=True, exist_ok=True)
            elif item.isfile():
                target.parent.mkdir(parents=True, exist_ok=True)
                source = bundle.extractfile(item)
                if source is None:
                    raise RuntimeError(f"Could not read archive member: {item.name!r}")
                with source, target.open("wb") as handle:
                    while block := source.read(1024 * 1024):
                        handle.write(block)
                extracted_files.append({"path": str(relative), "bytes": target.stat().st_size})
            else:
                raise ValueError(f"Unsupported archive member type: {item.name!r}")

    sqlite_files = [
        {
            "path": str(path.relative_to(temporary)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
        for suffix in ("*.db", "*.sqlite", "*.sqlite3")
        for path in sorted(temporary.rglob(suffix))
    ]
    os.replace(temporary, destination)
    payload = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "archive": str(archive),
        "archive_sha256": actual,
        "destination": str(destination),
        "extracted_file_count": len(extracted_files),
        "extracted_total_bytes": sum(int(item["bytes"]) for item in extracted_files),
        "sqlite_candidates": [item["path"] for item in sqlite_files],
        "sqlite_files": sqlite_files,
        "files": extracted_files,
    }
    manifest.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(manifest)


if __name__ == "__main__":
    main()

