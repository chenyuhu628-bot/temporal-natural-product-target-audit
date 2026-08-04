"""Deterministic compressed-text writers used by source-to-frozen steps."""

from __future__ import annotations

import gzip
import io
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, TextIO


PANDAS_GZIP = {"method": "gzip", "compresslevel": 9, "mtime": 0}


@contextmanager
def deterministic_gzip_text(path: Path, encoding: str = "utf-8") -> Iterator[TextIO]:
    """Write one gzip member with a blank filename and fixed timestamp."""

    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding=encoding, newline="") as text:
                yield text
