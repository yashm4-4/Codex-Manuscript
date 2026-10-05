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
