# Source-aware date-precision policy audit of temporal natural-product–target retrieval

**Author:** Chenyu Hu  
**Affiliation and contact:** to be confirmed by the author before journal submission

## Purpose

This repository contains version 1.0.0 of the reproducibility software for
“Source-aware date-precision policy audit of temporal natural-product–target
retrieval”. The work is a post hoc audit of four fixed retrieval baselines under
positive–unlabeled semantics. It is not a new predictive model.

The study examined how source-aware publication-date precision, ties, source
dependence, scaffold coldness, and target homology affect evaluation against
later-recorded natural-product–target relations. The current analysis entry
points are locked by `manifests/CURRENT_EXECUTION_CHAIN.json`.

Author-side clean-environment reproduction is **COMPLETED**: the 19-stage
source-to-frozen chain matched 16/16 locked inputs by SHA-256, the 12-step
analysis chain passed, seven revision/v4 validators passed, and the manuscript
aggregate checker passed 17/17 checks. Independent third-party reproduction is
**NOT YET PERFORMED**. See `reproduction/clean_environment/`.

## Sources and the NPASS boundary

The analysis uses NPASS 2.0 and 3.0, ChEMBL 31, PubMed date metadata, UniProt
mappings and sequences, RDKit 2026.03.4, and MMseqs2 18-8cc5c. Source URLs,
recorded acquisition dates, sizes, and checksums are in
`manifests/source_download_manifest.tsv`.

NPASS is handled as `LINK_ONLY_NO_REDISTRIBUTION`. No NPASS file, identifier
subset, structure, activity record, target–compound–PMID mapping, or other
row-level derivative is included. Obtain NPASS files from the official download
page: https://bidd.group/NPASS/downloadnpass.html. The version-specific papers
are cited in `docs/NPASS_CITATION_AND_PROVENANCE.md`.

This repository is not affiliated with or endorsed by the NPASS database or
its maintainers.

## Repository map

- `code/`: the 12-step authoritative analysis chain and necessary import-only
  dependencies.
- `reproduction/source_to_frozen/`: evidence-backed reconstruction steps that
  build the 16 local frozen inputs from official sources.
- `code/revision_audits/`: transparent post hoc sensitivity-audit sources;
  these do not extend the authoritative 12-step chain.
- `scripts/`: release checks, verified downloading, orchestration, and
  manuscript aggregate validation.
- `configs/` and `manifests/`: source, policy, path-template, and code-lock
  metadata.
- `tests/fixtures/`: identifier-free synthetic fixtures.
- `results/aggregate/`, `figures/`, and `tables/`: only cleared,
  non-reconstructive project aggregates.

`code/corrective/build_asof_corrective_bundle.py` and
`code/audits/audit_document_component_bootstrap_v1.py` are custody-only import
dependencies. They are not current workflow entry points.

## Environment

The recorded Windows environment used Python 3.11.15, NumPy 2.4.6,
scikit-learn 1.9.0, RDKit 2026.03.4, pandas 3.0.3, Biopython 1.87, and
Matplotlib 3.11.0. For the exact Windows solve:

```text
conda create --name npass_temporal_release --file environment/conda-win-64-explicit.txt
conda activate npass_temporal_release
```

For a lighter cross-platform environment:

```text
conda env create -f environment/environment.yml
conda activate npass_temporal_release
```

MMseqs2 is not bundled. Install release 18-8cc5c from its official project.

## Data acquisition and reconstruction

Check the recorded official endpoints without downloading data:

```text
python scripts/download_sources.py --manifest manifests/source_download_manifest.tsv --check-only
```

Download hash-pinned direct sources into an ignored local directory:

```text
python scripts/download_sources.py --manifest manifests/source_download_manifest.tsv --download --download-dir work/sources --allow-large
```

Review `docs/SOURCE_AT_ACQUISITION.md`, `docs/NPASS_LINK_ONLY_POLICY.md`, and
`reproduction/RECONSTRUCTION_REPORT.md` before processing. Copy
`configs/reproduction_paths.example.json` to the ignored
`configs/reproduction_paths.local.json`, enter only local paths, then run:

```text
python scripts/rebuild_analysis.py --config configs/reproduction_paths.local.json
```

The input gate fails closed on any absent prerequisite or mismatch against the
16 locked hashes. To build the historical directory layout in a new work area:

```text
python scripts/materialize_execution_layout.py --output work/execution-layout
```

No restricted input is copied by the materializer. Follow
`reproduction/protocol/execution-runbook.md` and the stage commands recorded in
`reproduction/SOURCE_TO_FROZEN_PROVENANCE.json`. PubMed
requests record the requested-ID-list hash, UTC request times, batch sizes,
retry attempts, and response receipts. Gzip products are written with fixed
metadata and sorted rows.

## Safe verification without third-party records

```text
python scripts/run_reproduction.py --mode verify-chain
python scripts/run_reproduction.py --mode smoke
python -m unittest discover -s tests -v
```

These checks verify the 12 locked entry-point hashes and run identifier-free
synthetic tests. They do not download or expose third-party data.

## Expected aggregate results

The frozen evaluation contains 4,123 candidate human single-protein targets,
222 queries, 358 later-recorded strict A/B relations, and 4,990 historical
relations. Under the conservative day-only policy, 13,885 of 20,647 source
rows were admitted. Interval censoring classified all 20,455 date-resolved rows
as certainly before the cutoff, changing 141 evidence tiers without changing
historical membership or cold-scope denominators. Manuscript number checks and
figure/table regeneration are reported in the release audit; the README does
not substitute for those machine-readable receipts.

## Known limitations

- Later recorded does not mean first biological discovery.
- Unrecorded compound–target pairs are unlabeled, not confirmed negatives.
- The audit is outcome-visible, author-run, post hoc, and not independent.
- PubMed and UniProt services are dynamic; requests are receipted, and exact
  source versions or responses may need archival recovery for byte-identical
  reconstruction.
- ChEMBL 31 requires roughly 4.51 GB compressed and 23.75 GB extracted.
- Independent third-party reproduction is not claimed before it occurs.
- Five regenerated PNGs differ at the rendering layer under the fresh plotting
  stack; all six figure source-data TSVs are exact and Figure 5 is pixel-exact.

## Licence boundary

Project-authored software is MIT licensed, Copyright (c) 2026 Chenyu Hu.
Only aggregate files explicitly cleared in the dataset manifest are CC BY 4.0.
Neither licence covers NPASS, ChEMBL, PubMed, UniProt, identifiers or records,
third-party binaries, software packages, weights, or caches. See
`LICENSE_SCOPE.md` and `THIRD_PARTY_NOTICES.md`.

## Public release identifiers

- Repository: https://github.com/chenyuhu628-bot/temporal-natural-product-target-audit
- GitHub release: https://github.com/chenyuhu628-bot/temporal-natural-product-target-audit/releases/tag/v1.0.0
- Software archive DOI: https://doi.org/10.5281/zenodo.21788846
- Aggregate dataset DOI: https://doi.org/10.5281/zenodo.21788854

All four identifiers were verified after authenticated publication on 4 August
2026. The Zenodo file checksums match the local release packages.
