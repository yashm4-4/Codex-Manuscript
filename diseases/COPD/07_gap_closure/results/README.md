# V2 results

Start with `COPD-V2-PHENO_phenotype_robustness_report.md`.

R001 is the 104-accession adjudication register; R002/R002B cohort provenance;
R003 study/tag/panel support; R004 all 15,389 records; R005 the unchanged-order
337 candidates; R006 frozen components; R007 the unchanged shortlist; R008
stratum counts; R009 exact versus aggregate ML dependence; R010 all tags;
R011 frozen locus labels with 140/152 universe flags; R012 nonsignificant
aggregate links; R013 unresolved tag identities; R014 frozen biological
annotation summaries; R015 single-accession classification influence; and
R016 leave-one-publication influence. Source and validation ledgers accompany
them. All IDs start `COPD-V2-PHENO` and do not replace V1 result IDs.

Support columns are nonexclusive. `primary_retained` tests the direct COPD
definition; `secondary_ehr_supported` and `ml_supported` remain distinct.
`gws_supporting_studies` requires valid significant study/tag/panel evidence;
`v1_aggregate_linked_studies` preserves original broader provenance.
`GCST90244098_only_study` is sole-accession support, not a biological negative.
Unresolved reference matching and unresolved phenotype adjudication differ.
An empty set-valued support field means no supporting member of that set,
not unavailable evidence; absent design information is labeled explicitly in
the accession register. Every V1 S4-R010 candidate retains its original rank.

## Frozen-model reverse-complement audit

Start with [the RC report](COPD-V2-RC_reverse_complement_robustness_report.md).
R001 contains raw RC scores (the original R003-compatible score column names
refer to RC values in this file). R002/R003 record transformations and all-record
reconciliation. R004 compares all 15,303 scored pairs; R005/R006 preserve the
337 and 12 original ranks. R007/R008 retain original literature coverage and
exact-variant comparisons. R009–R012 contain score, delta, sign and forward-tail
metrics; R013/R014 call and context transitions; R015/R016 threshold proximity;
R017 subset summaries; R018 locked severity; R019 gate failures; R020 analysis QC.
F001–F004 are the corresponding PNG/PDF figures.

The overall result is scientifically material. Of the frozen 337, 153 retain
any call and 118 retain exact model support. RC-only positives are diagnostic
observations, not a replacement candidate list. H3K27me3-associated is the
interpretation of the internal `silencer` model name. PHENO strata do not change
model decisions. Undefined metrics are explicitly `not_evaluable`.

## Frozen external functional benchmark

Start with [the benchmark report](COPD-V2-BENCH_external_functional_benchmark_report.md).
`COPD-V2-BENCH-R001` through `R006B` are the source inventory, exact identities,
frozen master, assay classifications, separate contextual/excluded evidence,
denominator/search audit and source-locus grouping. R007/R008 preserve all
14,025 evidence rows in the forward and forward/RC comparisons. R009–R011 are
study/context, mechanism and locus summaries; R012 records why full-panel metrics
are not calculated; R013/R014 contain descriptive distributions and coverage;
R015 deduplicates variants without hiding context conflicts; R016 is analysis QC;
R017/R018 show original gate failures and orientation categories.

Empty model scores/calls mean unevaluable, not zero or negative. Explicit assay
`null` means reported nonsignificance under that study's rule, not biological
inactivity. Castaldi count-only records have unavailable binary labels even
though measurements are present. Original source rsIDs and exact biological
alleles remain distinct from model coverage. `silencer` means H3K27me3-associated,
not demonstrated functional silencing. No orientation is selected using labels.

Forty variants have any in-scope positive evidence; 23 are exactly scored, one
recovered forward and none RC. Do not interpret 1/23 as population sensitivity.
All frozen source/master files and future-use restrictions are pinned in
`../provenance/COPD-V2-BENCH_benchmark_freeze.json` and the firewall document.
