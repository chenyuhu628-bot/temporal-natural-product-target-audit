# Frozen structure-standardization sensitivity protocol

Protocol ID: `npass_structure_policy_sensitivity_v1_20260729`

Status: frozen before inspection of any standardized-structure result.

Parent protocol: `npass_strict_ab_major_revision_v4_20260729`

Parent protocol SHA-256: `bd4902476160cc7c5cbacaf0cfd0f1a28c5300bde22232b3a2cc6c1f3c143dc2`

## Claim boundary

This is an outcome-visible, author-run, descriptive Analysis D sensitivity.
It does not alter the frozen endpoint, candidate target universe, historical
relations, evidence weights, target sequences, candidate masks, tie salt, or
homology masks. It is not preregistered confirmation, independent validation,
chemical identity adjudication, or biological validation.

Historical and query structures remain in distinct role-specific maps for
every policy. A structure appearing in both roles is transformed independently
in each map. No pooled structure map is permitted.

## Locked inputs

| Input | SHA-256 |
|---|---|
| Parent v4 protocol | `bd4902476160cc7c5cbacaf0cfd0f1a28c5300bde22232b3a2cc6c1f3c143dc2` |
| Historical pairs | `cef748ae8ac277e49784d7e1fbf08e085beb14a84fc6c651a0fb8d99e88710d7` |
| Scoring queries | `0e6068d2e25cb3ea325656fb3517563788cd496e88cfaa3de761890fec9e9318` |
| Historical compounds | `f1f82793b5c652007a042699c19cd5640a8b68e8b1f0d4f94e4dc4f54045060c` |
| Query compounds | `f51670fffd21d2e9109b4376dd53aab55bbacb09d2e4795dfa515fcaca98b113` |
| Candidate targets | `0ee86746b306fb388a1f74a6b88ce4d1eba01b7a4eb473315f6b3def57145cdc` |
| Candidate sequences | `a83421dba2482f236fe18340dd592cc7d5ed22c98c4fc39435c40f04f289b442` |
| Evaluation pairs | `09296b066a23197a7c178f00514f2b3d9ed7e6f3c459ea92a55e01a6010d1132` |
| Frozen scaffold audit | `fa0029ef5b7822ad5ca93f7bd93ac808f85f1e0c02e827fa91be375031b2d7af` |
| Homology mask 0.30 | `3a8247ed8f683fe6fce5fb345f56e3ec73a872b065eca922e92e494f084a1793` |
| Homology mask 0.50 | `ec3bbd435f73bc1c724efdfd94ac10c32d6f9a55bd6c9a9349233a01e5dc7f5b` |
| Homology mask 0.70 | `ec3bbd435f73bc1c724efdfd94ac10c32d6f9a55bd6c9a9349233a01e5dc7f5b` |
| Frozen complete ranks | `87739aa818744c7084088d13c386444aa41bbef38c257083325298003181479e` |
| Frozen aggregate metrics | `fac75f28185cc4e2d320ab45ddba5d07d4857c85ce477e62ec4ae960ff656ea8` |
| Frozen scorer implementation | `7b8263828ab6eaf4756307f315df4502b1615a40bffa8219284612585ad9bdc8` |
| Frozen rank/metric implementation | `7e40e80c3203a4a6cf95ca675cc9c333168f8feb234449856f26b5e132ec3165` |

## Version lock and molecule conversion

The approved runtime must report RDKit `2026.03.4`; any other version fails
closed. Each input is parsed with `Chem.MolFromSmiles()` and sanitized by
RDKit. Each successful policy result is sanitized again and represented by
isomeric canonical SMILES for internal deterministic comparison only. No
identifier-bearing standardized structure is retained.

The following five policies are evaluated without result-based selection:

1. `raw_primary`: the sanitized molecule parsed from the locked role-specific
   primary SMILES; no MolStandardize transformation.
2. `cleanup_fragment_parent`: `rdMolStandardize.Cleanup`, followed by
   `rdMolStandardize.FragmentParent(..., skipStandardize=True)`.
3. `cleanup_charge_normalized`: `rdMolStandardize.Cleanup`, followed by the
   default canonical `rdMolStandardize.Uncharger().uncharge` operation.
4. `cleanup_canonical_tautomer`: `rdMolStandardize.Cleanup`, followed by
   `rdMolStandardize.TautomerEnumerator().Canonicalize`.
5. `cleanup_parent_charge_tautomer`: Cleanup, FragmentParent with
   `skipStandardize=True`, Uncharger, then canonical tautomer enumeration in
   that order.

No policy is allowed to fall back to a different molecule when parsing,
transformation, sanitization, or InChIKey generation fails. Such a policy is
reported as blocked for downstream scoring; no structure or score is imputed.

## Structure and scaffold summaries

For each role and policy, retain only aggregate counts of parse/transformation
success and failure; canonical structure, Morgan fingerprint, full InChIKey,
connectivity layer, and Bemis–Murcko scaffold changes; nonempty versus empty
scaffolds; and collision counts. Full identifiers and standardized SMILES are
never written.

Scaffold coldness is defined exactly as in the frozen audit: an endpoint
relation is scaffold-cold only when its query has a valid nonempty Bemis–Murcko
scaffold absent from all eligible historical nonempty scaffold groups. Empty or
acyclic scaffolds are excluded rather than pooled. Report relation and query
entries, exits, and symmetric membership changes relative to `raw_primary`.

## Raw fail-closed calibration

Before any perturbed result is accepted:

- all 358 raw endpoint scaffold flags must reproduce the frozen scaffold audit;
- all four raw baseline score strings and integer ranks must reproduce every
  eligible cell of the 3,658,128-row frozen complete-rank ledger;
- raw Recall@50 must reproduce every frozen baseline-by-scope cell;
- the 0.50 and 0.70 masks must be identical both by locked file hash and joined
  relation membership.

Any mismatch stops the run and prevents a PASS receipt.

## Scoring and reporting

The frozen target universe, target masking, Tier A/B weights 1.0/0.7, Morgan
radius 2 with 2,048 bits, native target 3-mer representation, pair-neighbour
sequence top-100 rule, and deterministic tie salt remain unchanged. Structure
policies can directly affect Morgan transfer and structure–sequence pair
neighbour scores. Popularity and sequence-transfer baselines are nevertheless
retained for every policy and must be verified invariant rather than omitted.

For each policy and all four baselines, report aggregate complete-ledger score,
rank, and top-50 membership changes relative to raw. For each of the following
nonduplicated scopes, report query-macro Recall@50 and changes relative to raw:

- temporal strict A/B;
- policy-specific scaffold-cold;
- policy-specific joint scaffold–homology cold at identity 0.30;
- policy-specific joint scaffold–homology cold for the identical 0.50/0.70 mask.

For a policy-specific scope, separate the score/rank effect on the same
membership from the total change relative to the frozen raw scope. Both
favourable and unfavourable changes are retained. No p-values or superiority
claims are permitted.

## Output contract

The create-once directory contains this protocol, aggregate-only TSV/JSON
outputs, scripts, input/output hashes, a manifest, tests, a validator, and an
execution receipt. It contains no compound, target, pair, query, structure,
scaffold, or source-document identifiers; no standardized SMILES; and no
absolute paths.
