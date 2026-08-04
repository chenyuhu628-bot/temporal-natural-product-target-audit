#!/usr/bin/env python3
"""Derive Bemis-Murcko scaffold keys from prepared compound structures.

RDKit is intentionally an explicit runtime requirement.  The preparatory input
script can run without it; this script must not silently substitute a different
definition of chemical scaffold.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

from reproducible_io import deterministic_gzip_text


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return deterministic_gzip_text(path) if mode.startswith("w") else gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="compound_structures_for_rdkit.tsv.gz")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()
    try:
        from rdkit import Chem, rdBase
        from rdkit.Chem.Scaffolds import MurckoScaffold
    except ImportError as exc:
        raise SystemExit(
            "RDKit is required for true Bemis-Murcko scaffolds but is not available in this Python environment. "
            "Install or select an approved RDKit environment, then rerun this unchanged script."
        ) from exc
    if args.output.exists() or args.summary.exists():
        raise FileExistsError("Refusing to overwrite an existing scaffold output or summary")

    rows, invalid, acyclic = [], 0, 0
    with open_text(args.input, "rt") as source:
        reader = csv.DictReader(source, delimiter="\t")
        required = {"inchikey_full", "representative_smiles"}
        if not required.issubset(reader.fieldnames or set()):
            raise ValueError(f"Input missing required columns: {sorted(required)}")
        for item in reader:
            smiles = item["representative_smiles"].strip()
            molecule = Chem.MolFromSmiles(smiles)
            result = {
                "molecule_id": item.get("molecule_id", item["inchikey_full"]),
                "inchikey_full": item["inchikey_full"],
                "representative_smiles": smiles,
            }
            if molecule is None:
                invalid += 1
                result.update({"scaffold_status": "invalid_smiles", "bemis_murcko_smiles": "", "scaffold_key": ""})
            else:
                scaffold = MurckoScaffold.GetScaffoldForMol(molecule)
                scaffold_smiles = Chem.MolToSmiles(scaffold, canonical=True) if scaffold.GetNumAtoms() else ""
                if scaffold_smiles:
                    result.update({"scaffold_status": "ok", "bemis_murcko_smiles": scaffold_smiles, "scaffold_key": scaffold_smiles})
                else:
                    acyclic += 1
                    result.update({"scaffold_status": "acyclic_or_empty_bemis_murcko", "bemis_murcko_smiles": "", "scaffold_key": ""})
            rows.append(result)
    fieldnames = ["molecule_id", "inchikey_full", "representative_smiles", "scaffold_status", "bemis_murcko_smiles", "scaffold_key"]
    with open_text(args.output, "wt") as target:
        writer = csv.DictWriter(target, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "rdkit_version": rdBase.rdkitVersion,
        "input": str(args.input),
        "input_sha256": sha256(args.input),
        "output": str(args.output),
        "output_sha256": sha256(args.output),
        "counts": {"molecules": len(rows), "invalid_smiles": invalid, "acyclic_or_empty_bemis_murcko": acyclic, "nonempty_scaffold": len(rows) - invalid - acyclic},
        "warning": "Empty Bemis-Murcko scaffolds require an explicit acyclic-molecule split policy; do not treat empty strings as one shared scaffold group.",
    }
    args.summary.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(args.summary)


if __name__ == "__main__":
    main()
