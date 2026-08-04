"""Compatibility entrypoint for implementation amendment 2026-07-28-01."""

from __future__ import annotations

import rdkit
from rdkit.Chem import rdFingerprintGenerator

# The locked v1 builder imports this module from the rdkit package root. Expose
# the exact installed rdkit.Chem module at that location without changing any
# chemistry or builder logic.
rdkit.rdFingerprintGenerator = rdFingerprintGenerator

from build_asof_corrective_bundle import main


if __name__ == "__main__":
    raise SystemExit(main())

