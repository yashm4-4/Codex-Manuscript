# V2 data

The three `adjudication_*.json` files are publication-batch reviews of all
104 frozen core accessions. They were classified before retention counts were
examined. `phenotype_adjudication_amendments.json` records explicit root-review
wording and composite-eligibility corrections, not silent source replacement.

The pipeline preserves reviewer inputs, normalizes documented cohort aliases,
clarifies Catalog availability wording, and outputs the complete accession
register under `results/`. That register includes ancestry/stage metadata,
publication/DOI, source URLs, full-statistics status, rationale and uncertainty.
Hashes of every input are in the phenotype run manifest. These are small
manual provenance records, not new downloaded GWAS data. No large analytical
dataset has been downloaded for this module.

## Frozen external functional benchmark

`COPD-V2-BENCH/` contains separately authorized public primary-study downloads,
Castaldi GEO GSE109452 counts/barcode design, Gong supplementary design/results,
mechanistic papers/figures, retrieval failures, authoritative allele mappings,
and three independently assembled source-evidence JSON files. These are source
evidence, not new GWAS or participant-level data. Binary labels are distinguished
from count measurements; absent labels are never inferred nulls.

Every source file is pinned by `../provenance/COPD-V2-BENCH_benchmark_freeze.json`.
Do not edit or regenerate them in place after viewing model performance.
This is external evaluation data: neither labels nor results may guide training,
architecture, thresholds, controls, seeds, orientation handling or run selection.
