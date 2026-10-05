# Section 5: Computational Validation

This section consumes the final Section 4 predicted-causal candidate table. It
keeps variant-level, regional experimental, and gene-level contextual evidence
separate so that database overlap is not mistaken for causal proof.

## Reproducible workflow

Run with the system Python environment, which provides `pyarrow` 6.0.1 and
SQLite support:

```bash
python3 diseases/COPD/05_computational_validation/scripts/run_section5.py
```

The input resolver prefers
`04_modeling/results/COPD-S4-R005_predicted_causal_population_genetics.tsv`
because it includes exact Ensembl/dbSNP identifiers. It can fall back to the
R004 causal or prioritized table. For isolated validation tests, set the
task-specific `COPD_S5_CANDIDATE_INPUT` environment variable to an alternate
TSV; production runs do not set it.

The access audit does not require Section 4 and can be run independently:

```bash
python3 diseases/COPD/05_computational_validation/scripts/run_section5.py --audit-only
```

## Result identifiers

- `COPD-S5-R001`: exact GRCh38 chromosome-position-REF-ALT matches to GTEx v10
  Lung significant cis-eQTL pairs.
- `COPD-S5-R002`: MPRAbase v4.9.3 matches. GRCh38 candidate bases are lifted to
  hg19 as one-base intervals and must return exactly on reverse liftover.
- `COPD-S5-R003`: HGNC resolution and direct Open Targets evidence at the exact
  COPD node `MONDO_0005002`.
- `COPD-S5-R004`: explicit local/public/controlled resource and biobank audit.
- `COPD-S5-R005`: integrated one-row-per-candidate matrix and Section 5 report;
  Section 4 rank order is preserved rather than recomputed.

`results/COPD-S5_validation_checks.tsv` contains cross-file consistency checks.
Each result has a JSON manifest with input/output hashes and interpretation
guardrails. No source or decision identifiers are assigned in this section;
the final traceability merge is left to the coordinating agent.

## Interpretation boundaries

GTEx is exact molecular association evidence from bulk non-diseased lung.
MPRAbase coordinate containment shows that an assayed sequence covered the
candidate base, but it is not assumed to be an allelic test. Open Targets is
gene-disease context and does not validate a candidate variant or establish an
enhancer-gene connection. Missing controlled-access or biobank data are logged
as availability gaps, never as failed replication.
