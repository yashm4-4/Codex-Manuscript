# Section 6: Experimental and Biological Validation

Status: computational experimental design and collaborator handoff complete.
Wet-lab execution is pending; no experimental result or causal proof is
claimed.

Primary narrative:
`results/section_6_experimental_validation.md`

## Reproducible workflow

From the repository root, run:

```bash
python3 diseases/COPD/06_experimental_validation/scripts/run_section6.py
```

Production execution requires the canonical Section 4 R010 list and its
R004-R008 supporting outputs, the final Section 5 integrated candidate table,
and both validated 2,001-bp candidate FASTAs. The workflow fails when the R006 target maps are missing. Upstream
paths and SHA-256 checksums are frozen in
`data/COPD-S6-R005_upstream_inputs.tsv`.

## Result identifiers

- `COPD-S6-R001`: 12-candidate experimental shortlist with exact model,
  population, target, motif, and public-validation evidence columns.
- `COPD-S6-R002`: all-337-candidate eligibility, linked-tag redundancy, and
  selection audit.
- `COPD-S6-R003`: 48 allele-specific MPRA inserts, covering REF and ALT in
  both orientations with 100 genomic bases on each side.
- `COPD-S6-R004`: candidate-specific MPRA/CRISPRi/editing plans, a TF-first
  perturbation plan, and COPD-relevant cell-context/control matrix.
- `COPD-S6-R005`: collaborator README, member manifest, checksum-verified
  deterministic archive, and upstream input manifest.
- `COPD-S6-R006`: integrated Section 6 report and analysis manifest.

`results/COPD-S6_validation_checks.tsv` contains executable reconciliation and
sequence checks. All experimental work-plan rows are explicitly labeled
`proposed_not_performed`.

## Interpretation boundary

The R003 records are genomic inserts, not complete ordering-ready
oligonucleotides. Vector, promoter, adapters, cloning sites, barcode design,
and sequence remediation remain platform decisions. Each insert isolates one
nominated allele in GRCh38 reference-sequence flanks; local phased donor
haplotypes and nearby linked alleles are not represented. Relevant phased
haplotypes should be tested where biological or prior functional evidence
warrants it. MPRA tests episomal sequence activity; CRISPRi tests an endogenous
region; motif matching is not TF binding evidence. Allele-level causal support
requires a genotype-verified endogenous allele perturbation with replicated
molecular consequences and appropriate controls.
