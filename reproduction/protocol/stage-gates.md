# Corrective successor stage gates

| Gate | Requirement | Evidence | Status at protocol lock |
|---|---|---|---|
| R0 | project-lead local repair authorization and claim boundary | governance authorization record | complete |
| R1 | protocol, traceability, table and figure contracts hashed before results | protocol lock manifest | pending manifest generation |
| R2 | implementation and environment locked | code lock and environment receipt | pending |
| R3 | row ledger and role-separated structures pass fail-closed validation | input manifest and validation report | pending |
| R4 | endpoint/universe/method hashes equal frozen values | input integrity report | pending |
| R5 | four baselines and all scopes fully rerun | score/evaluation manifests | pending |
| R6 | before/after, tie, PMID, approximation and unresolved audits complete | aggregate audit receipt | pending |
| R7 | manuscript/review package consistency and release review pass | integration and release reports | pending |

## Fail-closed conditions

The run stops if any of the following occurs:

1. an existing output or legacy artifact would be overwritten;
2. a v3 row contributes to historical tier, weight, SMILES, fingerprint,
   scaffold, or InChI repair;
3. a missing/non-day/non-PMID/PubMed-missing/after-cutoff row is treated as
   eligible history;
4. any retained historical pair lacks an eligible pre-cutoff v2 row;
5. the historical keyset is not exactly the frozen 4,990 keys;
6. the endpoint file/hash, 358/222/156 denominator, 4,123-target universe, or
   target/sequence hashes change;
7. historical and query structure maps are pooled;
8. a required historical structure fails parsing and lacks a unique, verified,
   eligible-v2-only InChI repair;
9. a scaffold or homology mask changes without a new audited mask derivation;
10. a baseline, A/B weight, top-100 rule, mask, tie salt, metric, bootstrap
    seed, threshold, or candidate universe changes;
11. a result is viewed and then a locked rule is modified;
12. an identifier-rich or source-derived artifact is copied into a release
    package without field-level approval.

An implementation defect found before result interpretation requires a dated,
append-only amendment and a new code-lock manifest. A scientific result that is
unfavorable is not a defect and must be retained.


