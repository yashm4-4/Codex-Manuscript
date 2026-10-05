# COPD V2 gap-closure workspace

This directory is the isolated planning and execution namespace for the COPD
V2 gap-closure phase. The existing COPD Sections 1–6, manuscript, result IDs,
model artifacts, thresholds, and 337-candidate ordering are the frozen V1
record and must not be edited in place.

The investigator-authorized **GWAS phenotype-robustness module is complete**.
Read the [phenotype report](results/COPD-V2-PHENO_phenotype_robustness_report.md)
and [337-candidate support table](results/COPD-V2-PHENO-R005_frozen_337_phenotype_support.tsv).
Direct COPD susceptibility support retains 184/337 candidates, 72/153 frozen
components and 6/12 shortlist entries. This is phenotype dependence, not a
test of biological falsity. The original V1 ranking remains unchanged.

The investigator-authorized **frozen-model reverse-complement audit is also complete**.
Read the [RC audit report](results/COPD-V2-RC_reverse_complement_robustness_report.md)
and [rank-preserved 337 comparison](results/COPD-V2-RC-R005_frozen_337_comparison.tsv).
The prespecified severity is scientifically material: 153/337 retain a union
call, 118/337 retain exact model support, and 65/152 well-separated union
positives lose the call. All 15,303 frozen scorable pairs were audited.

The phenotype module was published at commit
`0ed50782139c9af113420e3e0a85518b620e9215`. The RC audit is local and stops for
investigator review. All other V2 modules await authorization; no retraining,
external functional benchmark, new cell context, fine-mapping, target analysis,
reranking or manuscript revision was performed by this audit.

## Identifier namespace

V2 records use `COPD-V2-*` identifiers. They must never reuse or reinterpret
`COPD-S1-*` through `COPD-S6-*` result identifiers.

## Contents

- `GAP_CLOSURE_PLAN.md`: investigator-facing scientific audit and execution plan.
- `gap_closure_decision_register.tsv`: V2 decisions and guardrails.
- `gap_closure_result_register.tsv`: V2 planning, phenotype and RC result registry.
- `activity_log.tsv`: V2 activity history.
- `data/`: reviewed accession-level phenotype inputs and explicit amendments.
- `results/`: phenotype and RC diagnostic tables, figures, QC and detailed reports.
- `scripts/`: phenotype reconstruction plus frozen-model RC scoring, analysis, figures and validation.
- `logs/`: phenotype and RC run/summary logs.
- `manuscript/`: future V2 supplements or amendment drafts, only after review.
- `provenance/`: frozen-V1 snapshot, planning audit, locked phenotype/RC specifications, manifests and checksums.

## Evidence policy

GWAS association, LD, statistical fine-mapping, model predictions, molecular
QTL association, GWAS–QTL colocalization, predictive enhancer–gene links,
physical contacts, motif similarity, occupancy, reporter assays, regional
perturbation, endogenous allele editing, and disease phenotypes remain separate
evidence classes. Missing coverage is not a negative result.
