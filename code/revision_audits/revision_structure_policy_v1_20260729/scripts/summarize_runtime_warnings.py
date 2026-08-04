#!/usr/bin/env python
"""Create an aggregate addendum for transient RDKit runtime warnings."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


OUTPUT_DIR = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: object) -> None:
    if path.exists():
        raise FileExistsError(f"Create-once output already exists: {path.name}")
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )


def main() -> None:
    stderr_path = OUTPUT_DIR / "build.stderr.log"
    stdout_path = OUTPUT_DIR / "build.stdout.log"
    if not stderr_path.is_file() or not stdout_path.is_file():
        raise FileNotFoundError("Transient build logs are absent")
    lines = stderr_path.read_text(encoding="utf-8").splitlines()
    counts = {
        "max_tautomers_warning_event_count": sum(
            "max tautomers reached" in line for line in lines
        ),
        "max_transforms_warning_event_count": sum(
            "max transforms reached" in line for line in lines
        ),
        "kekulization_warning_event_count": sum("Can't kekulize" in line for line in lines),
    }
    categorized = sum(counts.values())
    summary = {
        "schema_version": "structure_policy_runtime_warning_summary_v1",
        "status": "AGGREGATED",
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "rdkit_version": "2026.03.4",
        "warning_event_count": len(lines),
        **counts,
        "uncategorized_warning_event_count": len(lines) - categorized,
        "source_stderr_sha256": sha256(stderr_path),
        "source_stdout_sha256": sha256(stdout_path),
        "source_stdout_byte_count": stdout_path.stat().st_size,
        "temporary_logs_deleted_after_aggregation": True,
        "interpretation": (
            "Counts are runtime warning events, not distinct molecules. Default "
            "RDKit tautomer limits bounded the canonicalization search; a returned "
            "molecule was treated as a deterministic policy result, not as proof "
            "that every possible tautomer was enumerated. Kekulization warnings "
            "were emitted internally; all role-policy transformations nevertheless "
            "returned sanitized nonempty molecules and no imputation was used."
        ),
    }
    summary_path = OUTPUT_DIR / "runtime_warning_summary.json"
    write_json(summary_path, summary)
    manifest_addendum = {
        "schema_version": "structure_policy_warning_manifest_addendum_v1",
        "status": "COMPLETE",
        "primary_manifest_sha256": sha256(OUTPUT_DIR / "manifest.json"),
        "aggregate_only": True,
        "identifier_bearing_outputs_retained": 0,
        "outputs": {
            "runtime_warning_summary.json": sha256(summary_path),
        },
        "scripts": {
            "scripts/summarize_runtime_warnings.py": sha256(Path(__file__).resolve()),
            "scripts/validate_warning_addendum.py": sha256(
                OUTPUT_DIR / "scripts/validate_warning_addendum.py"
            ),
        },
        "transient_log_hashes_before_deletion": {
            "stderr": summary["source_stderr_sha256"],
            "stdout": summary["source_stdout_sha256"],
        },
    }
    write_json(OUTPUT_DIR / "warning_manifest_addendum.json", manifest_addendum)
    stderr_path.unlink()
    stdout_path.unlink()


if __name__ == "__main__":
    main()
