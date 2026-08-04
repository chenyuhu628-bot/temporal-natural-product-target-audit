# Corrective successor table schema

These are locked data contracts. Values must be generated mechanically from
the versioned execution and may not be hand-entered.

| Table | Purpose | Rows | Required fields | Real data source | Replacement owner |
|---|---|---|---|---|---|
| T1 | temporal repair flow and frozen-set verification | eligibility status and keyset checks | source rows; eligible/excluded counts; affected pairs; 4,990/358/222/156/4,123 equality; hashes | corrective row ledger and input receipt | execution script |
| T2 | historical evidence before/after audit | change category | tier/weight changes; SMILES/fingerprint/scaffold changes; date precision; PMID concentration | difference and provenance audits | audit suite |
| T3 | corrected aggregate performance | 5 scopes × 4 baselines | denominators; Recall@10/50; NDCG@10/50; MRR | corrected evaluator output | evaluator script |
| T4 | corrected bootstrap summaries | baseline/scopes/comparisons | estimate; 95% CI; replicates; seed; estimability | corrected evaluator output | evaluator script |
| T5 | score degeneracy and tie audit | 5 scopes or query set × 4 baselines | all-zero queries; positive-score targets; unique scores; top-10/50 boundary tie blocks; tie-determined counts | aggregate-only score audit | audit suite |
| S1 | scope, mask and integrity verification | scope/check | pair/query/target counts; scaffold/homology mask hashes; rank completeness; failures | integrity audit | audit suite |
| S2 | top-100 versus exhaustive fidelity | scope/metric | maximizer outside top-100; score error; rank correlation; top-10/50 change; metric difference | prespecified sensitivity | audit suite |
| S3 | full zero/failure accounting | scope × baseline | zero Recall/MRR; incomplete ranks; nonfinite scores; not-estimable reasons | corrected evaluator | evaluator script |
| S4 | PMID/source-document dependence | aggregate statistic | unique PMIDs; relation/query contribution; largest component/document share; date source and precision | provenance audit | audit suite |
| S5 | frozen unresolved exclusions | aggregate reason | compound/target/both/validation reason counts; denominator | unresolved audit ledger | audit suite |
| S6 | reproducibility and release | artifact/check | input/code/output hashes; software; runtime; memory; path and content scan; shareability tier | manifests and release gate | package builder |

## Aggregation contract

Metrics are per query and macro averaged. All 5 scopes, 4 baselines and 5
metrics must appear. Bootstrap is query-cluster paired where applicable. No row
may be omitted because its result is zero, unfavorable, tied, or imprecise.
Identifier-rich inputs and ranks remain restricted; manuscript tables contain
only approved aggregate values.


