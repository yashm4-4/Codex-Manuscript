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

The completed RC audit adds `COPD-V2-RC_*` specification, preflight/scoring/
analysis/figure manifests, original-freeze checks, input and output hash ledgers,
232-file tracked-V1 baseline and 3,841-entry filesystem baseline, forward-only
p99 cutoffs, and independent final validation. The locked specification hash
must remain unchanged. Failed first numerical-verification artifacts are
preserved alongside the successful corrected run. `COPD-V2-RC_runtime/` is an
ignored cache directory, not an analysis result.

Shared V2 registers/documentation evolve by authorized module additions; prior
PHENO checksum entries for these shared files describe their earlier completed
state. PHENO-specific outputs and the prior checksum ledger remain preserved.

The external benchmark adds `COPD-V2-BENCH_*`: protocol lock, 316-entry immutable
input baseline, assembly manifest, independent source/identity reviews, a
200-file benchmark freeze, comparison run history, interpretation review,
independent QC and completed-artifact checksums. The master was frozen before
new benchmark-wide model extraction. No frozen benchmark membership or labels
were changed after comparison. Implementation-only validation notes are separate.

`COPD-V2-BENCH_external_evaluation_firewall.md` is mandatory for future use:
labels are not training data and outcomes cannot select models or settings.
The current V1 analysis is descriptive, not an untouched V1 selection test.
Prior PHENO/RC-specific files and ledgers remain byte-identical; their hashes of
shared V2 registers/documentation describe earlier historical completion states.
