# V2 provenance

- `v1_snapshot.tsv` pins the immutable V1 repository state and core hashes.
- `planning_diagnostics.tsv` records lightweight read-only audit findings.
- `external_resource_audit.tsv` records every external dataset or paper
  considered for the plan, including access and statistics coverage.

- `planning_validation.tsv` is the manual/read-only QC register for the final
  plan, with each row naming its validation method.

Path convention for `planning_diagnostics.tsv`: semicolon-separated evidence
paths beginning with `models/` are repository-root-relative. Other paths are
relative to `diseases/COPD/`; a bare filename following a directory-qualified
entry inherits that preceding directory. Result IDs and decision IDs are
explicit non-path references.
