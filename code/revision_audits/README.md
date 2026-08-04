# Post hoc revision audit source

These project-authored scripts implement the date-precision, tie-aware,
structure-policy, weight-policy, unresolved-entity-bound, and reviewer-matrix
sensitivity analyses used for manuscript v4. They are retained as transparent
post hoc audit utilities and are **not** additional steps in the 12-step
authoritative execution chain in `manifests/CURRENT_EXECUTION_CHAIN.json`.

The clean reproduction used the locked scientific definitions and fresh local
manifests. One release-local portability adjustment in
`revision_unresolved_bounds_v1_20260729/scripts/run_unresolved_bounds.py`
updates expected gzip-container hashes only after decompressed content was
verified exact. The statistical calculations and aggregate TSV outputs are
unchanged. Details and comparisons are in
`reproduction/clean_environment/CLEAN_ENVIRONMENT_REPRODUCTION_REPORT.md`.

No revision-audit output directory, entity identifier, source record, sequence,
structure, rank, or per-query result is included here.
