"""Identifier-free contract tests for every source-to-frozen stage."""

from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import io
import json
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CONTRACT = json.loads((ROOT / "tests/fixtures/source_to_frozen_contract.json").read_text(encoding="utf-8"))


class SourceToFrozenSyntheticTest(unittest.TestCase):
    def test_fixture_is_identifier_free(self) -> None:
        text = json.dumps(CONTRACT, sort_keys=True)
        self.assertNotIn("NPC", text)
        self.assertNotIn("PMID", text)
        self.assertNotIn("InChIKey", text)
        self.assertTrue(CONTRACT["synthetic_only"])

    def test_every_stage_has_a_compilable_script_and_contract(self) -> None:
        stages = CONTRACT["stages"]
        self.assertEqual([row["id"] for row in stages], [f"{value:02d}" for value in range(19)])
        for stage in stages:
            with self.subTest(stage=stage["id"]):
                path = ROOT / stage["script"]
                self.assertTrue(path.is_file(), path)
                compile(path.read_text(encoding="utf-8"), str(path), "exec")
                self.assertGreater(len(stage["contract"]), 20)

    def test_reproducible_gzip_writer(self) -> None:
        module_path = ROOT / "reproduction/source_to_frozen/reproducible_io.py"
        spec = importlib.util.spec_from_file_location("reproducible_io", module_path)
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        with tempfile.TemporaryDirectory() as temporary:
            first = Path(temporary) / "first.tsv.gz"
            second = Path(temporary) / "second.tsv.gz"
            for path in (first, second):
                with module.deterministic_gzip_text(path) as handle:
                    writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
                    writer.writerow(["synthetic_compound", "synthetic_target"])
                    writer.writerow(["SYNTH_COMPOUND_A", "SYNTH_TARGET_A"])
            self.assertEqual(first.read_bytes(), second.read_bytes())
            self.assertEqual(gzip.decompress(first.read_bytes()), gzip.decompress(second.read_bytes()))
            self.assertEqual(hashlib.sha256(first.read_bytes()).hexdigest(), hashlib.sha256(second.read_bytes()).hexdigest())

    def test_reconstruction_matrix_covers_sixteen_inputs(self) -> None:
        with (ROOT / "reproduction/INPUT_RECONSTRUCTION_MATRIX.tsv").open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(len(rows), 16)
        self.assertEqual(len({row["input_name"] for row in rows}), 16)
        self.assertTrue(all(len(row["reference_sha256"]) == 64 for row in rows))


if __name__ == "__main__":
    unittest.main()
