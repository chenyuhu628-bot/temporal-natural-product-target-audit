# Date-precision policy sensitivity protocol v1

Analysis ID: `revision_date_policy_v1_20260729`

Parent frozen protocol: `npass_strict_ab_major_revision_v4_20260729`

Parent protocol SHA-256: `bd4902476160cc7c5cbacaf0cfd0f1a28c5300bde22232b3a2cc6c1f3c143dc2`

Status: implementation lock; written before computing scenario-level retrieval results.

## Claim boundary

This is an author-run, outcome-visible, post hoc descriptive sensitivity analysis.
It does not alter the frozen v3 manuscript, endpoint, target universe, baseline
definitions, homology thresholds, deterministic tie salt, or biological labels.
It is not independent validation and must not be used to select a favourable date
policy.

## Locked inputs

Inputs are read from the isolated corrective execution directory. Their exact
hashes are verified against the parent protocol before scoring. No row-level or
identifier-bearing output is written.

## Calendar-interval rules

The cutoff is 31 August 2022, inclusive. A day is represented as `[d, d]`, a
month as its first through last calendar day, and a year as 1 January through
31 December. Missing or non-PMID dates are not assigned an interval and remain
ineligible.

Three predeclared policies are evaluated:

1. `day_only_conservative`: include only PubMed-verified, day-precision rows on
   or before the cutoff; this exactly reproduces the frozen primary history.
2. `interval_certain_pre_cutoff`: include a PubMed-verified row only when the
   closed interval's upper bound is on or before the cutoff.
3. `interval_earliest_bound`: include a PubMed-verified row when the interval's
   lower bound is on or before the cutoff. This is explicitly liberal and could
   include cutoff-crossing rows if any existed.

All interval states (`definitely_before_or_on`, `crossing_cutoff`,
`definitely_after`, `unresolved_interval`) are counted before scenario scoring.
If predeclared scenarios collapse to identical selected row sets, each name is
retained in outputs and the equivalence is reported rather than hidden.

## Reconstruction and evaluation

For each policy, the analysis rebuilds:

- the best strict A/B evidence tier per historical pair (A dominates B) and the
  associated A=1.0/B=0.7 weight;
- the deterministic modal historical SMILES per full InChIKey, using
  lexicographic resolution of modal ties and the already validated v2 InChI
  repair only if the selected modal SMILES is unparsable;
- radius-2, 2,048-bit Morgan features and exact Bemis-Murcko scaffold keys;
- the four frozen baselines and salted complete candidate ranks;
- the frozen temporal and homology masks, with scaffold membership recomputed
  from the scenario-specific historical representatives.

Results are reported for five provenance scopes. For display, the identical
0.50 and 0.70 masks are combined as `joint_scaffold_homology_0.50/0.70`, while
both source hashes and separate internal denominator checks are retained.

Aggregate outputs include history/tier/structure changes, complete-score and
rank-change counts, endpoint-rank changes, scope denominators, and macro
Recall@50. No query-level ranks, pair identifiers, InChIKeys, accessions, SMILES,
or absolute filesystem paths are emitted.

## Failure rules

The run fails if any locked hash differs, any policy loses one of the frozen
4,990 historical pair keys, any required structure cannot be deterministically
resolved, candidate masks cease to be permutations, the primary policy fails to
reproduce frozen Recall@50 within 1e-12, or the 0.50 and 0.70 mask hashes differ.

