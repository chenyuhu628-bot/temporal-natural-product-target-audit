# Corrective execution runbook

1. Generate and verify the protocol lock manifest. Do not run builders before
   it exists.
2. Implement the row ledger, role-separated structure builder, bundle builder,
   scorer/evaluator wrapper, and audit suite in this successor directory.
3. Create an immutable code lock containing every executable source hash and an
   exact software/environment receipt.
4. Build a new restricted execution directory. Existing directories are fatal.
5. Validate source hashes, classify v2 strict rows, and write the restricted
   row ledger plus aggregate eligibility summary.
6. Rebuild the 4,990 historical pair table and separate historical/query
   structure tables. Validate the endpoint and universe hashes before scoring.
7. Re-derive fingerprints and scaffold assignments. Reuse a cold mask only
   after keyset and hash/diff evidence proves it is unchanged; otherwise derive
   and audit a new mask.
8. Run all four baselines with measured wall time and peak memory where
   supported. Write ranks only inside the restricted execution directory.
9. Evaluate all five scopes and five metrics with the locked bootstrap.
10. Run the prespecified aggregate audits, then the integrity audit. Do not
    inspect or write manuscript conclusions until all required checks pass.
11. Build the review-safe package from approved code and aggregate outputs only;
    scan both filenames and file contents for absolute paths and restricted
    fields.
12. Create a new manuscript v3 directory and mechanically integrate verified
    outputs. Preserve manuscript v2 unchanged.


