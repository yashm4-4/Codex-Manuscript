# COPD V2 gap-closure workspace

This directory is the isolated planning and execution namespace for the COPD
V2 gap-closure phase. The existing COPD Sections 1–6, manuscript, result IDs,
model artifacts, thresholds, and 337-candidate ordering are the frozen V1
record and must not be edited in place.

The present state is **planning only**. No V2 model training, genome-wide
fine-mapping, large public-data download, candidate reranking, or manuscript
rewrite has been performed.

## Identifier namespace

V2 records use `COPD-V2-*` identifiers. They must never reuse or reinterpret
`COPD-S1-*` through `COPD-S6-*` result identifiers.

## Contents

- `GAP_CLOSURE_PLAN.md`: investigator-facing scientific audit and execution plan.
- `gap_closure_decision_register.tsv`: V2 decisions and guardrails.
- `gap_closure_result_register.tsv`: V2 result registry; currently planning artifacts only.
- `activity_log.tsv`: V2 activity history.
- `data/`: future V2 inputs and input manifests only.
- `results/`: future V2 analytical outputs only.
- `scripts/`: future V2 code only; V1 code is not to be edited.
- `logs/`: future V2 run logs.
- `manuscript/`: future V2 supplements or amendment drafts, only after review.
- `provenance/`: frozen-V1 snapshot and planning/source audit records.

## Evidence policy

GWAS association, LD, statistical fine-mapping, model predictions, molecular
QTL association, GWAS–QTL colocalization, predictive enhancer–gene links,
physical contacts, motif similarity, occupancy, reporter assays, regional
perturbation, endogenous allele editing, and disease phenotypes remain separate
evidence classes. Missing coverage is not a negative result.
