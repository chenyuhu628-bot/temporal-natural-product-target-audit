# Source-at-acquisition reconstruction

This package uses a link-only reconstruction design for third-party data. It does not redistribute source files or source-derived identifier-bearing rows.

## Recorded sources

NPASS v2.0 and v3.0 static files, ChEMBL 31 structures/mapping/SQLite, the PubMed ESummary service, UniProt REST metadata, MMseqs2 18-8cc5c, and the RDKit 2026.03.4 conda artifact are listed in `manifests/source_download_manifest.tsv`. Static files retain the project-recorded acquisition timestamp, byte count, and SHA-256. Dynamic API services record the acquisition/query date and release evidence where available, but cannot be reproduced from a URL alone without the original identifier set.

## Fail-closed procedure

1. Check endpoint status with `python scripts/download_sources.py --manifest manifests/source_download_manifest.tsv --check-only`.
2. Download only to a new local directory. Do not place downloads under version control.
3. For every static hash-pinned file, verify byte count and SHA-256 before use. A mismatch is fatal.
4. Preserve download time, response headers, source release, and the exact request definition.
5. Reconstruct intermediate tables locally. Never copy NPASS, ChEMBL, PubMed, or UniProt rows into the release package unless a later field-level rights decision explicitly permits it.

## Source-specific limitations

- **NPASS:** Database-content redistribution terms remain unclear. The original NPASS files are link-only. The licences of articles describing NPASS do not grant a database-content licence. Written permission is required before releasing NPASS-derived row-level fields.
- **ChEMBL 31:** The official FAQ reports CC BY-SA 3.0, while EMBL-EBI terms also require attention to resource-specific and third-party rights. Exact snapshots, fields, attribution, and share-alike implications require review.
- **PubMed:** Only identifiers and date metadata were used. Abstracts and full text are excluded. Reproduction requires a locally reconstructed PMID request set derived from lawfully acquired NPASS records.
- **UniProt:** The project mapping receipt recorded UniProt release 2026_02, while the candidate-target sequence archive was identified as ChEMBL component-sequence database version 2022_02. The exact accession request set is not released. This difference must be documented and independently reconstructed.
- **RDKit:** Locked version 2026.03.4.
- **MMseqs2:** Locked release/commit identifier 18-8cc5c with 0.80 coverage and identity thresholds 0.30, 0.50, and 0.70. The executable is not shipped.

HTTP 200 and absence of a `WWW-Authenticate` header show only technical accessibility at the check time. They do not establish copyright, database rights, licence scope, or permission to redistribute.

