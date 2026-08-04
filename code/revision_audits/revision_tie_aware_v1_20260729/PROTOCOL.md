# Tie-aware retrieval and dependence sensitivity

Analysis ID: `revision_tie_aware_v1_20260729`

Status: post hoc, reviewer-requested, author-run sensitivity.

This create-once analysis implements Analysis A and Analysis B of the frozen
major-revision protocol
`npass_strict_ab_major_revision_v4_20260729` (SHA-256
`bd4902476160cc7c5cbacaf0cfd0f1a28c5300bde22232b3a2cc6c1f3c143dc2`).
It reads the frozen complete score/rank ledger, strict A/B endpoint, scope
masks, and local source-evidence ledger. It does not alter the v3 manuscript,
endpoint, scores, ranks, or any existing output.

## Exact tie treatment

Distinct score levels retain their frozen order. Within each exact-score block,
all permutations are treated as equally likely. For Recall@10/50 and
NDCG@10/50 the analysis reports:

- the fixed SHA-256-salted realization;
- the uniform within-tie expectation (fractional metric);
- exact attainable worst and best bounds;
- counts of relevant relations whose top-k membership is
  `score_identifiable`, `boundary_tie_dependent`, or `not_retrieved`;
- counts of queries whose Recall or NDCG is score-identifiable versus
  tie-dependent.

The 3-mer baseline is additionally split into score-operational queries (at
least one positive eligible score) and structural all-zero queries. Metrics
for the latter describe tie allocation only and are not interpreted as
algorithmic ranking performance.

## Uncertainty estimands

The query bootstrap resamples observed queries and targets the equally
query-weighted macro mean under a query-independence approximation. The
PMID-component bootstrap constructs a scope/subset-specific undirected
query–PMID graph, resamples its non-overlapping connected components, retains
all queries in each selected component, and normalizes by the sampled query
count. It targets the same macro mean while diagnosing source-document
dependence. It is not document-disjoint validation.

Both bootstraps use 10,000 PCG64 replicates with independently derived,
cell-stable seeds. The base seeds are 2026072901 (query) and 2026072902
(PMID-component).

For each displayed project-defined joint scaffold–homology cold scope and
baseline, the analysis also reports the exact one-sided 95% Clopper–Pearson
upper bound for the probability that a query has at least one later-recorded
target in the fixed-salt top 50. A zero-width empirical resampling interval is
labelled `empirical_degenerate_zero_hits`.

Because the 0.50 and 0.70 homology inputs are byte-identical, they are displayed
once as `joint_scaffold_homology_cold_0_50_0_70_identical_mask`; both input
hashes remain separately receipted.

## Release boundary

All outputs are aggregate-only. No query, compound, target, pair, PMID, rank
ledger, or connected-component membership is written. Paths in outputs and
receipts are basenames or project-relative names, never absolute paths.
