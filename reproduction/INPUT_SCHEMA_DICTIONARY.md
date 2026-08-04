# Input schema dictionary

The 16 inputs are local, row-level reconstruction products and are never
release data. Only schemas, rules, counts, missingness summaries, and hashes
are public. Full paths belong in ignored local configuration files.

## Evidence records

`npass_v2_evidence` and `npass_v3_evidence` contain source-version, source
entity/reference fields, compound identity/structure fields, target and assay
fields, quantitative activity fields, canonical pair key, evidence tier and
weight, temporal status, UniProt mapping status, and sequence length/MD5. The
complete ordered 43-field schema is enforced by `audit_npass_core.py`,
`build_evidence_tiers.py`, and `map_uniprot_target_sequences.py`.

`pubmed_v2_metadata` contains `pmid`, found flag, raw PubMed date strings,
normalized publication date, date precision and source, plus optional DOI,
journal, and title. Missing date fields remain missing; no day is invented.

`npass_v2_structures` is the unmodified official four-column NPASS 2.0 source
table (`np_id`, `InChI`, `InChIKey`, `SMILES`). It is used locally and is not
redistributed.

## Scoring inputs

- `frozen_historical_pairs`: canonical pair key, full compound key, canonical
  target accession, best strict evidence tier, decision, and PU policy.
- `frozen_scoring_queries`: deterministic query label and full compound key.
- `frozen_compounds`: full compound key and representative structure string.
- `candidate_targets`: one canonical target accession per row.
- `candidate_sequences`: one FASTA record per candidate target.

## Evaluation inputs

- `frozen_endpoint`: canonical pair/query/compound/target, best tier,
  temporal decision, and ChEMBL 31 gate status.
- `scaffold_audit`: pair key, selected-policy scaffold-cold flag, outcome, and
  explicit eligibility/exclusion reason.
- `homology_0_30`, `homology_0_50`, `homology_0_70`: canonical target,
  homology-cold candidate flag, and status.

## Candidate-universe inputs

`chembl31_target_catalogue` contains canonical accession, internal ChEMBL
target references, name/taxon, sequence length/MD5, source/version, and the
candidate-universe role. `chembl31_component_sequences` is the corresponding
one-record-per-target FASTA.

Exact reference rows, missingness, hashes, producer stages, sorting rules,
software, and seeds are listed in `INPUT_RECONSTRUCTION_MATRIX.tsv`.
