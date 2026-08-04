# Method–experiment traceability

| Contribution or issue | Method module | Experiment/audit | Table/Figure | Allowed claim | Evidence status before run |
|---|---|---|---|---|---|
| Historical evidence must be temporally pure | row-level v2 PMID/date eligibility ledger | full before/after row, tier and weight audit | T1, T2, Fig. 1B | corrected history uses only proven pre-cutoff v2 evidence | locked; pending execution |
| Endpoint must not be redesigned after opening | byte-hashed frozen endpoint and C31 universe | keyset/count/hash equality checks | T1, S1 | future endpoint and universe were unchanged during correction | locked; pending verification |
| Retrieval comparison must remain fair | same four baselines, masks, tie salt and candidate universe | complete rescoring of all methods/scopes | T3, T4 | corrected baselines were compared under the same frozen task | locked; pending execution |
| Approximation fidelity requires disclosure | top-100 primary plus exhaustive sensitivity | score/rank/top-K/metric difference audit | S2 | top-100 approximation fidelity was quantified descriptively | locked; pending execution |
| Sparse scores and ties require accounting | degeneracy and boundary-tie audit | all-zero, coverage, unique-score and tie-block counts | T5, S3 | deterministic tie dependence was explicitly measured | locked; pending execution |
| Source-document dependence limits independence | PMID provenance graph summaries | attrition, concentration and connected-component audit | T2, S4 | document dependence was measured as provenance, not validation | locked; pending execution |
| Entity exclusions may define selection boundaries | frozen-65 reason decomposition | aggregate reason counts only | S5 | unresolved entity exclusions were characterized without re-entry | locked; pending execution |
| Reproducibility must respect rights | manifests, environment lock and review-safe package | hash/path/sensitive-content audits | S6 | code and aggregate reconstruction materials are available within approved limits | pending release review |

## Claims that must not survive

- restored blindness, independent external validation, or prospective proof;
- biological discovery, direct binding, inactivity, or true-negative claims;
- a new algorithm or performance improvement selected after endpoint access;
- using exhaustive pair-neighbour or a different tie rule as the primary result;
- using C31, C37, shared PMID, or no-hit records as independent labels;
- public release of restricted fields without the separate rights gate.


