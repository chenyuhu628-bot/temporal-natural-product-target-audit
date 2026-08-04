"""Run safe verification modes for the local release candidate."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_chain() -> int:
    chain = json.loads((ROOT / "manifests/CURRENT_EXECUTION_CHAIN.json").read_text(encoding="utf-8"))
    rows = []
    for step in chain["authoritative_steps"]:
        path = ROOT / step["package_path"]
        observed = digest(path) if path.is_file() else None
        rows.append({"step": step["step_id"], "path": step["package_path"], "expected": step["sha256"], "observed": observed, "match": observed == step["sha256"]})
    print(json.dumps(rows, indent=2))
    return 0 if all(row["match"] for row in rows) else 1


def smoke() -> int:
    environment = dict(**__import__("os").environ)
    environment["PYTHONPATH"] = __import__("os").pathsep.join(
        [str(ROOT / "scripts"), str(ROOT / "code/corrective"), str(ROOT / "code/audits"), environment.get("PYTHONPATH", "")]
    )
    return subprocess.call([sys.executable, "-m", "unittest", "tests.test_synthetic_smoke", "-v"], cwd=ROOT, env=environment)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", required=True, choices=("verify-chain", "smoke"))
    args = parser.parse_args()
    return verify_chain() if args.mode == "verify-chain" else smoke()


if __name__ == "__main__":
    raise SystemExit(main())

