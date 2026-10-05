# COPD Section 6 collaborator handoff

Package status: computational design complete; all experiments are proposed and none has been performed.

## What is included

- `COPD-S6-R001_candidate_shortlist.tsv`: 12 nonredundant candidates with exact upstream evidence fields and explicit selection reasons.
- `COPD-S6-R002_selection_audit.tsv.gz`: every Section 4 predicted-causal record, LD-component assignment, eligibility, and selection outcome.
- `COPD-S6-R003_MPRA_constructs.tsv`: 48 validated insert designs, representing each shortlisted variant's REF and ALT alleles in both orientations.
- `COPD-S6-R004_candidate_validation_plan.tsv`: staged MPRA, endogenous CRISPRi, and allele-editing proposals.
- `COPD-S6-R004_TF_first_plan.tsv`: one leading sequence-motif hypothesis for each candidate/model combination and a binding/perturbation test plan.
- `COPD-S6-R004_cell_context_controls.tsv`: airway, alveolar, fibroblast, and technical-optimization contexts with controls and limitations.
- `COPD-S6-R005_upstream_inputs.tsv`: immutable paths, sizes, and SHA-256 checksums for every consumed upstream artifact.
- `section_6_experimental_validation.md`: integrated interpretation and design summary.
- `COPD-S6_validation_checks.tsv`: machine-readable reconciliation checks.

## Sequence-design boundary

The R003 sequences are genomic inserts, not complete synthesis oligonucleotides. Each forward insert contains 100 genomic bases, the normalized VCF allele, and 100 genomic bases. Reverse-complement designs are exact reverse complements. REF and ALT indel constructs therefore can differ in length. Each construct isolates one nominated allele in GRCh38 reference-sequence flanks; local phased donor haplotypes and nearby linked alleles are not represented. Effects can differ on native haplotypes, so relevant phased haplotypes should be tested when biological or published evidence warrants it. The receiving laboratory must choose the MPRA backbone, promoter, adapters, cloning sites, barcode architecture, minimum barcode count, and any restriction-site remediation before ordering.

## Evidence boundary

The shortlist preserves the canonical Section 4 R010 computational priority and adds auditable diversity constraints. A TREDNet call is a model prediction. TF-MoDISco/JASPAR and PWM results are sequence hypotheses, not TF binding evidence. GTEx is an exact variant-gene association in bulk lung. MPRAbase coordinate containment is regional assay coverage unless allelic testing is independently established. R006 target mappings are provisional. Neither selection nor database overlap constitutes causal proof.

## Recommended staged handoff

1. Platform review of insert architecture and sequence flags.
2. Allele-specific MPRA in pre-screened airway, alveolar, and fibroblast contexts.
3. Endogenous region-level CRISPRi only in contexts with accessible chromatin, target expression, and interpretable reporter or model evidence.
4. Genotype-verified allele editing for candidates with reproducible reporter and region-level effects.
5. TF binding and perturbation tests using the R004 TF-first hypotheses, including an allele-by-TF interaction test.

All statistical models, exclusion rules, biological replicates, and progression criteria should be preregistered with the receiving laboratory before unblinding primary readouts.
