# Tier B weight-policy sensitivity

Analysis ID: `revision_weight_policy_v1_20260729`

This post hoc, author-run sensitivity implements the Tier B portion of
Analysis D in the frozen major-revision protocol. Tier A is fixed at 1.0 and
Tier B is evaluated at 0.5, the frozen primary value 0.7, and 1.0. The
1.0/1.0 condition is the all-equal policy; it is not an additional fourth
variant.

The corrected day-only historical relations, endpoint, role-separated
structures, target universe, query masks, scaffold and homology masks, four
baselines, pair-neighbour top-100 approximation, and SHA-256 tie salt remain
fixed. No weight is selected using the resulting metrics.

All three weights are evaluated in one process. Molecular fingerprints,
per-query chemical similarities, sequence representation, query-specific
sequence similarities, and pair-neighbour target similarities are reused
across weights. Every eligible candidate receives a complete score and rank
under every weight. The 0.7 scores and ranks are checked row by row against
the frozen 3,658,128-row primary ledger.

Only aggregate results are written: complete-rank and top-50 change counts,
scope cardinality invariance, macro retrieval metrics, and metric differences
from 0.7. Identifier-bearing score/rank rows are not released. The analysis is
descriptive and does not add p-values, winner selection, external validation,
or biological evidence.

The byte-identical 0.50 and 0.70 homology masks are displayed once as a joint
0.50/0.70 scope while their two input hashes remain separately receipted.
