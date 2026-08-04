"""Create the immutable pre-result protocol-lock manifest."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "manifests" / "protocol_lock_manifest_v1.json"
LOCKED = [
    "PROTOCOL_STATUS.md",
    "configs/asof_cutoff_rebuild_spec_v1.json",
    "plan/experiment-protocol.md",
    "plan/review/method-experiment-traceability.md",
    "plan/stage-gates.md",
    "tables/table-schema.md",
    "figures/data-manifest.md",
    "execution/execution-runbook.md",
    "governance/project_lead_authorization_20260728.md",
    "scripts/lock_protocol.py",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    if OUTPUT.exists():
        raise FileExistsError(f"Refusing to overwrite protocol lock: {OUTPUT}")
    missing = [item for item in LOCKED if not (ROOT / item).is_file()]
    if missing:
        raise FileNotFoundError(f"Protocol lock inputs missing: {missing}")
    payload = {
        "schema_version": "1.0",
        "protocol_id": "npass_strict_ab_asof_cutoff_corrective_successor_v1_20260728",
        "locked_at_utc": datetime.now(timezone.utc).isoformat(),
        "lock_timing": "before corrective result generation",
        "execution_mode": "author_run_non_independent_corrective_successor",
        "files_sha256": {item: sha256(ROOT / item) for item in LOCKED},
        "legacy_artifacts_mutated": False,
        "public_release_authorized": False,
    }
    OUTPUT.write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(OUTPUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

