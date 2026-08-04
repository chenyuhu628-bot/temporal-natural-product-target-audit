# NPASS citation and provenance

The analysis used NPASS 2.0 and NPASS 3.0 as separately versioned source
snapshots. The files, acquisition timestamps, byte counts, hashes, purposes,
and transformations are recorded in `manifests/npass_source_at_acquisition.tsv`.

Version citations:

- NPASS 2.0: Zhao et al. (2023), “NPASS database update 2023: quantitative
  natural product activity and species source database for biomedical
  research”, *Nucleic Acids Research* 51(D1):D621–D628,
  https://doi.org/10.1093/nar/gkac1069.
- NPASS 3.0: Lin et al. (2026), “NPASS 3.0: an updated natural product activity
  and species source database with integrated AI-enabled insights”, *Nucleic
  Acids Research* 54(D1):D1519–D1528,
  https://doi.org/10.1093/nar/gkaf1196.

Official download page: https://bidd.group/NPASS/downloadnpass.html

The project records provenance facts but does not redistribute or relicense
NPASS database content. The public transformation chain operates on locally
downloaded, hash-verified files and emits row-level products only into ignored
local working directories.
