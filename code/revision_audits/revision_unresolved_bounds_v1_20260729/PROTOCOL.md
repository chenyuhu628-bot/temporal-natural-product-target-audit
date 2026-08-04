# Entity-unresolved endpoint bounds protocol v1

Analysis ID: `revision_unresolved_bounds_v1_20260729`

Parent protocol: `npass_strict_ab_major_revision_v4_20260729`

Parent protocol SHA-256: `bd4902476160cc7c5cbacaf0cfd0f1a28c5300bde22232b3a2cc6c1f3c143dc2`

Status: locked before calculation of new bounds.

## Claim boundary

This is an author-run, outcome-visible, post hoc descriptive missing-endpoint
analysis. The 65 frozen entity-unresolved relations remain excluded from the
primary endpoint and are never interpreted as negatives, remapped, readmitted,
or assigned invented query/target ranks.

## Locked restricted inputs

| Role | Relative input | SHA-256 |
|---|---|---|
| Initial 442-relation C31 decision ledger | `results/strict_temporal_future_v1_1_pmid_verified_chembl31_leakage_decision_ledger.csv.gz` | `15a85f1994e78d07856a1d5ebed2ca790f3ba12592e51634e0f8a7c00893b33b` |
| Frozen 65-relation unresolved set | `data/processed/strict_temporal_future_v1_1_pmid_verified_chembl31_C31_entity_unresolved.csv.gz` | `bf6e922568954cef5a750df6bfddaafa0f556e794a556ad9e5dc8dbc5e93bec0` |
| Preliminary C31 mapping ledger | `data/interim/chembl_31_future_candidate_entity_mapping.csv.gz` | `246f04875bdeb3cd087419c92e3d660f2cbb9060272348ac9cd2152976f89b6d` |
| SQLite mapping validation ledger | `data/interim/chembl_31_future_candidate_sqlite_entity_validation.csv.gz` | `7d66328e4cca3fca567235efdcd448dbe4e5929c1934fc3fe73180f9b989bd2d` |
| Frozen 358-relation endpoint | `author_run_strict_ab_asof_cutoff_execution_v1_20260728/evaluation_inputs/evaluation_pairs.tsv.gz` | `09296b066a23197a7c178f00514f2b3d9ed7e6f3c459ea92a55e01a6010d1132` |
| Frozen complete ranks | `author_run_strict_ab_asof_cutoff_execution_v1_20260728/score/corrective_prediction_ranks.tsv.gz` | `87739aa818744c7084088d13c386444aa41bbef38c257083325298003181479e` |

## Analyses

1. Partition the initial candidate ledger into 358 resolved/no-historical-
   activity relations, 19 resolved/historical-activity exclusions, and 65
   entity-unresolved relations. Reproduce the preliminary failure strata.
2. Compare only observable source-record characteristics across the 65
   unresolved relations, the 358 retained resolved relations, and all 377
   entity-resolved relations. Report counts, proportions, and distribution
   summaries without p-values or causal missingness claims.
3. Report the identified endpoint cardinality interval `358–423`. The upper
   endpoint adds all 65 unresolved relations but retains the 19 definitive C31
   historical-activity exclusions.
4. For each frozen baseline, count observed endpoint relations ranked at most
   50. Over the 423 potentially eligible relations, report the sharp
   assumption bounds:
   - all unresolved fail: `observed_hits / 423`;
   - all unresolved succeed: `(observed_hits + 65) / 423`.

These are relation-weighted temporal top-50 bounds. They are not query-macro
Recall@50, NDCG, MRR, scaffold-cold, or joint scaffold-homology bounds.

## Non-identifiability rule

Query-macro and scope-specific retrieval estimands are not reported because
the unresolved mappings can change the query set, per-query relation counts,
target identity, candidate masks, ranks, scaffold membership, and homology
membership. Inventing any of these values would not be a mathematical bound on
the frozen estimand.

## Output contract

Only aggregate TSV/JSON/Markdown artifacts are emitted. No pair, compound,
target, document, query, SMILES, or absolute-path field is written. All source
and output hashes, tests, and execution/validation receipts are retained.

