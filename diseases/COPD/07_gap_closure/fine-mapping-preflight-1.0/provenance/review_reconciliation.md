# Pre-freeze review reconciliation

Independent read-only reviews identified integration issues before final validation. The saved source bytes were not changed to resolve them.

- Catalog-only risk, limitation, ancestry/sample and LD-requirement tables are retained under `catalog_*` names. Their authoritative counterparts now integrate verified original-provider evidence, including all three original BBJ accessions. A missing Catalog mirror no longer overrides successful original-provider access.
- The Catalog P-field recognizer explicitly retains `neg_log_10_p_value` and `neglog10_pval`. These encodings are not treated as ordinary P values.
- The Pan-UKB source crosswalk resolves GCST90691934 to J44 EUR and GCST90692407 to J44 CSA. GCST90692991 remains the EUR+CSA meta-analysis, with no single-European LD substitution.
- Harmonized allele-column names are reported from each actual header; the presence of a harmonized file is not interpreted as proof of any one fixed column name.
- The download-scope declaration distinguishes complete small public dbGaP display exports from complete genome-wide GWAS releases; no complete genome-wide GWAS or LD matrix values were acquired.
- BBJ's pooled Japanese population-control contribution is explicit in the report and overlap contract. BBJ ancestry alone does not prove exact cohort-matched LD.
- The first independent validator run retained 4 failures among 4,805 checks: an ambiguously named provider-override flag was incorrectly expected to be true for Kim's Catalog-hosted files as well as BBJ's provider-only files. The flag is now explicitly named `catalog_missing_source_access_resolved_by_provider` (true for BBJ). Revised validation directly checks public signed-source evidence and absence of stale missing-source blockers for all seven, and separately checks provider-override semantics. Original failed receipts and the validator version are retained under `validation_attempt_01/`; no saved source bytes or scientific conclusion changed.

All changes were source-metadata/schema integration or documentation corrections. No association scan, statistical fine-mapping, model inference, candidate lookup, external benchmark evaluation or earlier-stage scientific modification was involved. The final independent validation and freeze bind the corrected artifacts.
