# Experiment protocol: row-level as-of-cutoff corrective successor

## 1. Objective and defect being corrected

The run repairs historical feature construction for the strict A/B temporal
retrieval study. The prior implementation admitted a historical pair when at
least one v2 row had a day-precise PMID on or before 2022-08-31, but then allowed
all strict v2 rows for that pair, and cross-version tier aggregation, to affect
historical evidence and structures. The correction constructs a row-level
as-of ledger and permits only proven pre-cutoff v2 rows to contribute.

This is an author-run, non-independent, post hoc corrective successor. The
future endpoint was already visible in the project and is frozen rather than
redesigned. The run cannot restore blindness or establish external validation.

## 2. Dataset and split strategy

### Historical side

The historical keyset is exactly the existing 4,990 canonical pairs. Each pair
must contain at least one v2 strict A/B P1 row satisfying all of:

1. full InChIKey, canonical reviewed-human UniProt mapping, and sequence found;
2. `ref_id_type=PMID` with a numeric PMID found in the archived PubMed table;
3. publication date precise to day and no later than 2022-08-31.

Every v2 strict row is assigned exactly one status defined in the JSON spec.
Only `eligible_pre_cutoff` rows determine pair tier, A/B weight, representative
historical SMILES, Morgan fingerprint, and historical scaffold. Missing dates,
non-day precision, non-PMID references, PubMed misses, and after-cutoff dates
are excluded without imputation. No v3 row may contribute to history.

### Future side

The exact 358 relations / 222 queries / 156 targets from the prior successor
endpoint are retained byte-for-byte by hash. Their strict-v3-only tier has
already been checked against the current tier with zero differences. The
endpoint is not reconstructed, expanded, contracted, or filtered in this run.

### Candidate universe and masking

All methods rank the same fixed 4,123 C31 human SINGLE PROTEIN targets. For each
query, targets already associated historically with that compound are masked.
Unrecorded pairs are unlabeled and never treated as negatives.

## 3. Structure construction and role separation

Historical and query structures are separate role-specific maps. Historical
representative SMILES is the modal nonempty SMILES among eligible v2 rows, with
lexical tie breaking. Query structures use the already frozen endpoint-side
representation. A shared InChIKey must not cause the two maps to be pooled.

If a required historical structure cannot be parsed, repair is allowed only
from v2 source records whose `source_np_id` belongs to an eligible row. Raw
InChI must regenerate the requested full InChIKey, conversion back to canonical
isomeric SMILES must regenerate the same key, and all valid candidates must
agree uniquely. Otherwise the run fails closed.

## 4. Locked baselines and fairness

The four baselines remain unchanged:

| Baseline | Locked role |
|---|---|
| weighted target popularity | deterministic target-frequency reference with A=1.0/B=0.7 |
| sequence 3-mer transfer | maximum weighted character-3-mer TF-IDF cosine from historical targets |
| weighted Morgan transfer | maximum A/B-weighted Morgan Tanimoto, radius 2 and 2,048 bits |
| structure–sequence pair neighbour | maximum product over the candidate target's fixed top-100 historical sequence neighbours |

No new predictor, tuning, changed feature, changed weight, or result-dependent
threshold is allowed. The pair-neighbour primary implementation remains
top-100. An exhaustive 1,131-history-target calculation is prespecified only as
a fidelity sensitivity and cannot replace the primary result.

## 5. Metrics and imbalance handling

The unit is query compound and aggregation is macro average. Recall@50 is
primary; Recall@10, NDCG@10, NDCG@50, and MRR are secondary. All five frozen
scopes and all four baselines are reported. Query-cluster bootstrap uses PCG64,
10,000 replicates, seed 20260719, and percentile 95% intervals. The analysis
uses retrieval metrics because absent/unrecorded relations are unlabeled.
AUROC, ordinary AUPR, calibration, and binary accuracy remain prohibited.

## 6. Main comparison and correction estimand

The primary correction question is whether row-level temporal purification
changes historical evidence, rankings, or aggregate conclusions. The full
before/after comparison reports pair keysets, row eligibility, tier/weight,
role-specific SMILES, fingerprints, scaffolds, masks, ranks, metrics, and
bootstrap contrasts. All differences are reported; the corrected result is not
selected based on direction or apparent performance.

## 7. Robustness and audit families

Required audits are:

- per-baseline all-zero vectors, positive-score coverage, unique score levels,
  top-10/top-50 boundary tie blocks, and deterministic-tie dependence;
- top-100 versus exhaustive pair-neighbour fidelity without promoting the
  exhaustive variant to the primary analysis;
- PMID counts, date-source/precision attrition, document contribution
  concentration, and query–PMID connected-component summaries;
- aggregate-only decomposition of the frozen 65 entity-unresolved relations;
- measured runtime and peak memory for the corrective run where instrumentation
  is available; no retroactive runtime is fabricated.

This benchmark does not propose a new model, so no model-module ablation or XAI
claim is applicable. The top-100 fidelity calculation is a method audit, not an
ablation claiming improvement.

## 8. Reporting and claim boundaries

The correction may support claims only about construction and retrieval of
later-recorded strict A/B candidate relations under the stated frozen universe.
`Later recorded` is not `first biologically discovered`; `unrecorded` is not
`negative`; C31 and C37 are not external validation sources; and PMID overlap
is reported as provenance dependence rather than independent confirmation.

All identifier-rich ledgers, ranks, structures, sequences, and source-derived
records remain restricted pending field-level rights review. Aggregate outputs
may enter the manuscript only after integrity and release gates pass.


