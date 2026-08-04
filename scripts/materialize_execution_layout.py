#!/usr/bin/env python3
"""Materialize the historical execution layout from the portable source tree."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
STRICT = Path("strict_ab_asof_cutoff_successor_v1_20260728")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path, help="New empty work directory to create")
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"Refusing to overwrite existing path: {output}")

    chain = json.loads((ROOT / "manifests" / "CURRENT_EXECUTION_CHAIN.json").read_text(encoding="utf-8"))
    copies: list[tuple[Path, Path]] = []
    for step in chain["authoritative_steps"]:
        source = ROOT / step["package_path"]
        observed = sha256(source)
        if observed != step["sha256"]:
            raise ValueError(f"Authoritative source hash mismatch: {step['package_path']}")
        source_project_path = Path(step["source_project_path"])
        if source_project_path.parts[0] == STRICT.parts[0]:
            destination = output / source_project_path
        else:
            destination = output / source_project_path
        copies.append((source, destination))

    support = [
        (ROOT / "code" / "corrective" / "asof_common.py", output / STRICT / "scripts" / "asof_common.py"),
        (ROOT / "code" / "corrective" / "asof_successor_common.py", output / STRICT / "scripts" / "asof_successor_common.py"),
        (ROOT / "code" / "audits" / "audit_common.py", output / STRICT / "audit_suite_v1_20260728" / "scripts" / "audit_common.py"),
        (ROOT / "scripts" / "pu_baseline_io.py", output / "scripts" / "pu_baseline_io.py"),
        (ROOT / "scripts" / "pu_retrieval_metrics.py", output / "scripts" / "pu_retrieval_metrics.py"),
        (ROOT / "scripts" / "pu_retrieval_scores.py", output / "scripts" / "pu_retrieval_scores.py"),
    ]
    for path in sorted((ROOT / "manifests").glob("*.json")):
        support.append((path, output / STRICT / "manifests" / path.name))
    copies.extend(support)

    output.mkdir(parents=True)
    receipt_files: list[dict[str, object]] = []
    for source, destination in copies:
        if not source.is_file():
            raise FileNotFoundError(source)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source, destination)
        receipt_files.append(
            {
                "path": destination.relative_to(output).as_posix(),
                "sha256": sha256(destination),
                "bytes": destination.stat().st_size,
            }
        )

    receipt = {
        "schema_version": "portable_execution_layout_v1",
        "status": "MATERIALIZED",
        "authoritative_chain": chain["chain_id"],
        "restricted_inputs_included": False,
        "third_party_data_included": False,
        "files": sorted(receipt_files, key=lambda row: str(row["path"])),
        "next_step": "Supply locally acquired restricted inputs under the paths documented in the runbook; do not commit them.",
    }
    (output / "MATERIALIZATION_RECEIPT.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps({"status": "PASS", "files": len(receipt_files)}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
