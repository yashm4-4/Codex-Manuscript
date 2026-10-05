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

All other V2 modules await investigator review. No model training, new allele
scoring, fine-mapping, target mapping, large public-data download, candidate
reranking, manuscript rewrite or new GitHub push has been performed.

## Identifier namespace

V2 records use `COPD-V2-*` identifiers. They must never reuse or reinterpret
`COPD-S1-*` through `COPD-S6-*` result identifiers.

## Contents

- `GAP_CLOSURE_PLAN.md`: investigator-facing scientific audit and execution plan.
- `gap_closure_decision_register.tsv`: V2 decisions and guardrails.
- `gap_closure_result_register.tsv`: V2 planning and phenotype-module result registry.
- `activity_log.tsv`: V2 activity history.
- `data/`: reviewed accession-level phenotype inputs and explicit amendments.
- `results/`: phenotype-support tables, source ledger, checks and detailed report.
- `scripts/`: offline phenotype reconstruction and final validation; V1 code is unchanged.
- `logs/`: small phenotype-module runtime/summary log.
- `manuscript/`: future V2 supplements or amendment drafts, only after review.
- `provenance/`: frozen-V1 snapshot, planning audit, phenotype specification, manifests and checksums.

## Evidence policy

GWAS association, LD, statistical fine-mapping, model predictions, molecular
QTL association, GWAS–QTL colocalization, predictive enhancer–gene links,
physical contacts, motif similarity, occupancy, reporter assays, regional
perturbation, endogenous allele editing, and disease phenotypes remain separate
evidence classes. Missing coverage is not a negative result.
