"""Check or download link-only third-party sources with fail-closed hashes."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import urllib.error
import urllib.request
from pathlib import Path

LARGE_LIMIT = 500_000_000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def check(row: dict[str, str]) -> dict[str, object]:
    method = row.get("status_method", "HEAD") or "HEAD"
    request = urllib.request.Request(
        row["url"], method=method, headers={"User-Agent": "JCheminform-source-check/1.0"}
    )
    try:
        with urllib.request.urlopen(request, timeout=30) as response:
            if method == "GET":
                response.read(65536)
            return {
                "source": row["source"],
                "status": response.status,
                "auth_challenge": response.headers.get("WWW-Authenticate"),
                "content_length": response.headers.get("Content-Length"),
            }
    except urllib.error.HTTPError as error:
        return {"source": row["source"], "status": error.code, "error": str(error)}
    except Exception as error:
        return {"source": row["source"], "status": None, "error": f"{type(error).__name__}: {error}"}


def download(row: dict[str, str], output: Path, allow_large: bool) -> dict[str, object]:
    expected_hash = row.get("sha256", "").strip()
    expected_bytes = int(row["bytes"]) if row.get("bytes", "").strip() else None
    if row.get("direct_download") != "true" or not expected_hash:
        return {"source": row["source"], "status": "SKIPPED_NOT_STATIC_HASH_PINNED"}
    if expected_bytes and expected_bytes > LARGE_LIMIT and not allow_large:
        return {"source": row["source"], "status": "BLOCKED_REQUIRES_ALLOW_LARGE"}
    target = output / row["filename"]
    part = target.with_suffix(target.suffix + ".part")
    if target.exists() or part.exists():
        raise FileExistsError(f"Refusing to overwrite {target} or {part}")
    request = urllib.request.Request(row["url"], headers={"User-Agent": "JCheminform-source-download/1.0"})
    with urllib.request.urlopen(request, timeout=120) as response, part.open("xb") as handle:
        for block in iter(lambda: response.read(1024 * 1024), b""):
            handle.write(block)
    if expected_bytes is not None and part.stat().st_size != expected_bytes:
        raise ValueError(f"Byte-count mismatch for {row['source']}")
    observed = sha256(part)
    if observed != expected_hash:
        raise ValueError(f"SHA-256 mismatch for {row['source']}: {observed}")
    part.replace(target)
    return {"source": row["source"], "status": "VERIFIED", "path": str(target), "sha256": observed}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", required=True, type=Path)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--check-only", action="store_true")
    mode.add_argument("--download", action="store_true")
    parser.add_argument("--download-dir", type=Path)
    parser.add_argument("--allow-large", action="store_true")
    args = parser.parse_args()
    manifest_rows = rows(args.manifest)
    if args.check_only:
        result = [check(row) for row in manifest_rows]
    else:
        if args.download_dir is None:
            parser.error("--download-dir is required with --download")
        args.download_dir.mkdir(parents=True, exist_ok=False)
        result = [download(row, args.download_dir, args.allow_large) for row in manifest_rows]
    print(json.dumps(result, indent=2, ensure_ascii=False))
    bad = [item for item in result if item.get("status") in (None, "BLOCKED_REQUIRES_ALLOW_LARGE") or isinstance(item.get("status"), int) and int(item["status"]) >= 400]
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())

