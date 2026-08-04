"""Pin recorded gzip header metadata without changing compressed payloads."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import struct
from pathlib import Path


def digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def header_filename(data: bytes) -> str:
    if data[:3] != b"\x1f\x8b\x08" or not data[3] & 8:
        raise ValueError("Expected a gzip file with an original-filename header")
    start = 10
    end = data.index(0, start)
    return data[start:end].decode("latin-1")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve()
    config = json.loads(args.config.read_text(encoding="utf-8"))
    receipts = []
    for entry in config["files"]:
        path = (root / entry["path"]).resolve()
        path.relative_to(root)
        data = path.read_bytes()
        observed_name = header_filename(data)
        if observed_name != entry["header_filename"]:
            raise ValueError(f"Unexpected gzip header filename for {path}: {observed_name!r}")
        observed_content = digest(gzip.decompress(data))
        if observed_content != entry["decompressed_sha256"]:
            raise ValueError(f"Decompressed-content mismatch for {path}")
        patched = data[:4] + struct.pack("<I", int(entry["mtime"])) + data[8:]
        observed_gzip = digest(patched)
        if observed_gzip != entry["expected_gzip_sha256"]:
            raise ValueError(f"Pinned gzip hash mismatch for {path}: {observed_gzip}")
        temporary = path.with_suffix(path.suffix + ".header-pin.tmp")
        if temporary.exists():
            raise FileExistsError(temporary)
        temporary.write_bytes(patched)
        temporary.replace(path)
        receipts.append(
            {
                "path": entry["path"],
                "mtime": entry["mtime"],
                "header_filename": observed_name,
                "decompressed_sha256": observed_content,
                "gzip_sha256": observed_gzip,
            }
        )
    print(json.dumps(receipts, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
