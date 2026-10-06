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
`0ed50782139c9af113420e3e0a85518b620e9215`; the RC audit was subsequently published
at `6343deb2b8dbacb31ea3e64dcd68e5e6e9b6d8e1`.

The separately authorized **external functional benchmark is complete and frozen**.
Read the [benchmark report](results/COPD-V2-BENCH_external_functional_benchmark_report.md)
and [frozen master](results/COPD-V2-BENCH-R003_frozen_benchmark_master.tsv).
The benchmark has 14,025 evidence rows and 1,731 reported rsIDs. Of 40 variants
with any in-scope positive assay, 23 have exact existing V1 sequence pairs:
one is recovered forward (rs2013701), none in RC. These are selected case-series
counts, not sensitivity. Castaldi GEO count/design data and Gong public null
results are preserved, but neither complete labeled/QC denominator passes.
The [future-evaluation firewall](provenance/COPD-V2-BENCH_external_evaluation_firewall.md)
prohibits training or model selection using benchmark labels/outcomes.

This benchmark remains local and stops for investigator review. No retraining,
new sequence/model inference, new cell-context model, matched controls,
fine-mapping, target analysis, reranking or manuscript revision was performed.
No next module, commit or push is authorized by this completion.

## Identifier namespace

V2 records use `COPD-V2-*` identifiers. They must never reuse or reinterpret
`COPD-S1-*` through `COPD-S6-*` result identifiers.

## Contents

- `GAP_CLOSURE_PLAN.md`: investigator-facing scientific audit and execution plan.
- `gap_closure_decision_register.tsv`: V2 decisions and guardrails.
- `gap_closure_result_register.tsv`: V2 planning, phenotype, RC and benchmark result registry.
- `activity_log.tsv`: V2 activity history.
- `data/`: reviewed phenotype inputs plus checksummed public benchmark sources and adjudications.
- `results/`: phenotype, RC and benchmark tables, figures where generated, QC and detailed reports.
- `scripts/`: module-specific reconstruction, source collection, frozen-score comparisons and validation.
- `logs/`: phenotype and RC run/summary logs; benchmark run provenance is under `provenance/`.
- `manuscript/`: future V2 supplements or amendment drafts, only after review.
- `provenance/`: frozen snapshots, locked specifications, benchmark data freeze/firewall, manifests and checksums.

## Evidence policy

GWAS association, LD, statistical fine-mapping, model predictions, molecular
QTL association, GWAS–QTL colocalization, predictive enhancer–gene links,
physical contacts, motif similarity, occupancy, reporter assays, regional
perturbation, endogenous allele editing, and disease phenotypes remain separate
evidence classes. Missing coverage is not a negative result.
