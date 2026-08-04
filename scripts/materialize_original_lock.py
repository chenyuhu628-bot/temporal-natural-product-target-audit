"""Materialize the original CRLF protocol-lock bytes in a temporary work area."""

from __future__ import annotations

import argparse
import hashlib
from pathlib import Path

EXPECTED = "96befee13ae1d41ad433c8697fac92ccd30fb25e24c3cf1279c6b4b7e040abd9"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path("manifests/protocol_lock_manifest_v1.json"))
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"Refusing to overwrite {args.output}")
    text = args.source.read_text(encoding="utf-8").replace("\r\n", "\n")
    payload = text.replace("\n", "\r\n").encode("utf-8")
    observed = hashlib.sha256(payload).hexdigest()
    if observed != EXPECTED:
        raise ValueError(f"Materialized protocol lock hash mismatch: {observed}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_bytes(payload)
    print(f"materialized_sha256={observed}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

