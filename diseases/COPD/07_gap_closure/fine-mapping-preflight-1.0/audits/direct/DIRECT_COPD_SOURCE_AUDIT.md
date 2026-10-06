# Direct-COPD source and signed-statistics access audit

Audit date: 2026-10-06 UTC. This is a preflight, not fine-mapping or locus selection.

All 27 direct-COPD accessions were selected from the frozen 104-study inventory solely by the phenotype adjudication. The candidate-blind projection is `inputs/frozen_study_metadata.tsv`. Candidate membership, candidate-support counts, functional benchmark outcomes and model predictions did not determine source selection or eligibility. The machine-readable companion is `audits/direct/direct_source_audit.json`; every accession also has a separately acquired current Catalog API/FTP audit under `audits/catalog/`.

## Main findings

| Source | Actually inspected public material | Signed, complete regional statistics? | Preflight disposition |
|---|---|---|---|
| Sakornsakolpat 2019, GCST007692 | Primary paper, consortium announcement, dbGaP analysis page, and entire small `pha004766` public export | **No.** Export has 19,373 records versus 6,224,355 tested variants; beta is explicitly absolute; effect-allele display has an inconsistent first record | Request original complete signed data and ancestry-specific/matched-LD information; not executable from this export |
| Hobbs 2017, GCST004147 | Primary paper and entire `pha004496/497/498` public exports | **No.** 25,412 / 25,338 / 25,350 records, coefficients explicitly absolute; some coordinates absent | Request original signed Stage-1 statistics; selected Stage-2 replication is not a dense final genome-wide release |
| Kim 2021, GCST90016588/589/593/594 | Raw GRCh37, formatted and harmonized GRCh38 headers and representative records, source metadata, primary methods | Actual genome-wide files are anonymously retrievable; explicit effect/other alleles, OR, SE and P. Whole-file and locus coverage/QC are not yet verified | Promising direct-COPD smoking-stratum inputs, conditional on source-scale confirmation, matched signed LD and full-file QC |
| Ishigaki 2020, GCST90013709/746/781 | Primary article/methods and frozen phenotype metadata; source-specific NBDC/BBJ retrieval independently audited under `audits/alternatives/` | Do not infer absence from Catalog FTP404. Use the actual provider-specific file/allele audit | Japanese direct-COPD option; integrate alternatives and LD audits before execution recommendation |
| Joo 2022, GCST90103984/985 | Primary paper and current Catalog API/FTP attempts | No complete signed sex-specific file obtained; paper lead/gene tables do not suffice | Author/source request and matched UKB LD needed |
| Cho 2010/2011/2014, GCST000603/001321/002350/002351 | Primary papers plus current Catalog API/FTP attempts | No complete signed public source file verified | Request original discovery release; preserve hg18 versus later builds and selected-replication scope |
| Moll 2021, GCST011766 | Primary paper and figshare API file inventory | Exonic ascertainment prevents dense locus completeness even if all selected exonic results were available; figshare lists a DOCX appendix | Unsuitable as a dense regional primary fine-mapping input |
| Other 11 direct-COPD accessions | Frozen source adjudication and current accession-level API/FTP checks | No source-specific full signed file verified in this sub-audit | Explicitly unresolved; no claim of global unavailability or invented schema |

No dataset is execution-cleared by this direct-source sub-audit alone. A retrievable signed summary file is necessary but does not establish usable matched LD or locus-level adequacy. The separate alternatives audit may verify additional provider files; its actual evidence must be integrated, not overwritten by a Catalog-only availability flag.

## Sakornsakolpat: public access does not mean usable full signed statistics

The [primary study](https://www.nature.com/articles/s41588-018-0342-2) analyzed 35,735 cases and 222,076 controls across UK Biobank and ICGC, testing 6,224,355 variants. It reported 82 loci: 47 previously known and 35 novel, using ±1 Mb regions. These are published discovery counts, not prospective fine-mapping eligibility counts.

The paper names dbGaP `phs000179.v5.p2` and UK Biobank in data availability. The [consortium announcement](https://sites.google.com/a/channing.harvard.edu/icopd-genetics-consortium/announcements/icgc-biobank-paper-published) links the [specific dbGaP analysis](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/analysis.cgi?study_id=phs000179.v6.p2&pha=4766), under the later v6 study version. That page identifies the analysis as the ICGC/UKB COPD meta-analysis and explicitly lists European, African American, Hispanic and East Asian populations. Thus the main-result ancestry assessment is not inferred solely from aggregate Catalog labels.

The [public parent analyses directory](https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs000179/analyses/) is accessible even though the version-specific study directories do not contain the association exports and the current Catalog summary-stat directory returns HTTP404. The [public `pha004766` export](https://ftp.ncbi.nlm.nih.gov/dbgap/studies/phs000179/analyses/phs000179.pha004766.txt) was retrieved in full because it is only 1,749,842 bytes. Its complete SHA-256 is `9c94aa9086ce658c154f9cccfde43146734ddbc90391a7236d4d9509fe103574`.

This is **not** the required original complete signed GWAS release:

- It contains 19,373 data records, not the 6,224,355 tested variants.
- Its `&beta;` dictionary says “Absolute value of regression coefficient”; all effect strings in this export are nonnegative.
- Its first record, `rs12045923`, displays `Allele1=C`, `Allele2=G`, but `Coded Allele=T`. This cannot be silently harmonized into a valid tested-allele pair.
- The `Sample size` column is entirely blank; a `Minor allele` label is not an allele-frequency field.
- The display export labels its coordinates genome build 38, whereas the submitted association methods describe hg19 coordinates. Original association and display coordinate/allele fields cannot be mixed.

Neither a coded/minor allele label nor an absolute coefficient permits recovery of the lost sign. It is not legitimate to infer a disease-increasing allele or to reconstruct a regional denominator from these records or Catalog significant hits. An author/ICGC request or an appropriately authorized original-data route is required for the missing complete signed file; this audit does not establish that the required original file is already available to this account through controlled access.

The paper's 10,000-participant unrelated-UKB LD reference for GCTA-COJO is a methodological description, not a verified public signed LD download. Its adequacy for a mixed-ancestry combined meta-analysis is a separate question. Prefer ancestry-specific summaries with corresponding LD, or a justified meta-analysis-aware design; do not substitute a European reference merely because it is accessible.

## Hobbs and Cho: keep discovery versus selected follow-up distinct

The [Hobbs study](https://www.nature.com/articles/ng.3752) explicitly included non-European populations in its main Stage-1 fixed-effects analysis and separately examined European-only results. It combined 22 genome-wide and four custom-content cohorts. The paper's 15,256 cases / 47,936 controls at Stage 1 and 9,498 / 9,748 at Stage 2 must not be silently interchanged with differently staged Catalog sample descriptions.

The paper reported 13 Stage-1 genome-wide loci and 22 after selected Stage-2 follow-up. Stage 2 tested selected top Stage-1 results, not every genome-wide variant. Public exports are incomplete and explicitly unsigned:

| dbGaP analysis | Public export records | Bytes | Complete small-export SHA-256 |
|---|---:|---:|---|
| pha004496 | 25,412 | 2,528,808 | `443f6aac27a9799ba9a2e31f8375cfd1969eb78f340e540373b60018af3c6add` |
| pha004497 | 25,338 | 2,519,756 | `3705b3d4655d7308fe58bd99f7ee6362115ca184fc1dd537582d1f3ad9ec0131` |
| pha004498 | 25,350 | 2,522,624 | `8e8d01938011dcf78c39cfb99c021aaab605e25e23ef16412d1079e7afdcb8af` |

These are actual full small export hashes, **not** hashes of an original full GWAS. The export field `|&beta;|` is defined as an absolute coefficient. `dbgap_public_export_schema_audit.json` records exact headers, row counts, field missingness, first records, source methods and reasons for rejection.

Cho's [2010](https://pmc.ncbi.nlm.nih.gov/articles/PMC2828499/) and [2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC3298111/) papers use hg18/NCBI36 coordinates. The [2014 study](https://pmc.ncbi.nlm.nih.gov/articles/PMC4176924/) analyzes moderate-to-severe and severe COPD separately, combining European and African-American discovery cohorts and selected ICGN follow-up. Its six headline loci span related analyses and must not be assigned as the eligible-locus count of either accession. Shared COPDGene/GenKOLS/NETT-NAS/ECLIPSE/ICGN participants across publications are not independent replication datasets for a future meta-analysis.

## Kim: a usable signed schema, but missing fields and an unverified denominator remain

The [Kim methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC8096488/) distinguish the smoking-stratified marginal logistic GWAS from interaction and 2-df joint tests. The four direct-COPD accessions are:

| Accession | Stratum | Cases | Controls |
|---|---|---:|---:|
| GCST90016588 | Ever smokers | 12,446 | 59,145 |
| GCST90016589 | Never smokers | 8,631 | 120,544 |
| GCST90016593 | Current smokers | 4,589 | 10,001 |
| GCST90016594 | Noncurrent smokers | 16,488 | 169,688 |

The raw GRCh37 files explicitly contain `variant_id`, `p_value`, `chromosome`, `base_pair_location`, `effect_allele`, `other_allele`, `odds_ratio`, `standard_error`. Representative non-significant records demonstrate this is not simply a significant-hit table. Current Catalog metadata reports 9,886,854 tested variants per accession, but a bounded prefix does not verify all rows or adequate coverage at every locus. Original and harmonized prefix bytes, sizes and hashes are in each accession's `audits/catalog/` directory.

The paper reports PLINK 2 logistic regression. [Official PLINK output documentation](https://www.cog-genomics.org/plink/2.0/formats) defines logistic standard errors on the beta/log-odds scale, not the OR scale. This supports a future log(OR)/SE contract. However, these older Catalog files use a generic `standard_error` field and do not retain the original PLINK header or full transformation lineage; a future execution stage must verify the source scale and numerical consistency rather than treating the renamed column as conclusive. No such statistical checks were computed in this preflight.

The raw files omit AF, per-variant N and per-variant imputation quality. The formatted/harmonized records have placeholder AF/beta fields, with NA in inspected examples; merely having a column does not provide the values. Study-level filters were MAF ≥0.01 and imputation r² ≥0.5. These filters do not establish current locus-by-locus quality or justify permissive future gates.

Do not mix `hm_*` fields with original fields. The first ever-smoker raw record has effect allele C and OR 0.97009; its harmonized representation has effect allele G and OR 1.0308321908276552. Both representations can encode the same contrast only when the allele pair, effect orientation and build are kept together. The example is a schema observation, not an inferred risk-allele analysis.

The published **48 / 55** loci are explicitly from the **2-df joint interaction tests** and are not the four marginal-GWAS locus counts. No stratum-specific eligible-locus count is available from this preflight. Ever/never is one partition and current/noncurrent another; the four releases overlap and must not be pooled as independent cohorts. An ever-smoker primary direct-COPD analysis, if selected for its prespecified disease estimand and satisfactory LD/QC, should retain never-smoker and alternative-partition results as distinct analyses, not pretend it reconstructs the broader Sakornsakolpat meta-analysis.

## Other direct sources

The [Ishigaki paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC7968075/) uses physician-diagnosed Japanese BBJ COPD and reports unrestricted disease summary statistics through JENGER/NBDC hum0014. This direct clinical phenotype is distinct from the later expanded BBJ EHR-derived releases. The independent alternatives audit owns exact provider retrieval and allele2/ALT effect interpretation. Methods specify GRCh37, 1000 Genomes phase 3 v5 imputation, and exclusion of imputation Rsq below 0.7. Imputation-reference availability is not equivalent to cohort-matched signed LD availability. The study's all-disease locus totals must not be reported as COPD-specific counts.

The [Joo paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC8861659/) reports 17 male and 14 female genome-wide loci (nine shared) in white-British UKB direct spirometric COPD. Its female-control count differs between abstract/Catalog (123,714) and methods (123,741); that discrepancy is retained. Without a complete signed release, lead tables and gene-based output do not establish a regional denominator. Its UKB participants overlap other UKB COPD resources.

The [Moll exonic study](https://pmc.ncbi.nlm.nih.gov/articles/PMC8321852/) retained 20,536 functional exonic variants after selecting from 109,036 array variants; it was not designed as a dense noncoding-locus assay. The [figshare record](https://api.figshare.com/v2/articles/14538222) lists only a 6,177,157-byte DOCX supplementary appendix, not a complete genome-wide summary-statistics file. Its 80 exome-significant variants and 35 clumped lead variants are not a valid full regional fine-mapping denominator. A better LD reference cannot recover missing association statistics for the unassayed variant set.

## Provenance, failures and scope of checks

`source_access_ledger.json` consolidates 37 public acquisition records with exact URLs, effective URLs, UTC times, HTTP status/error, retained-body path, byte count, SHA-256 and truncation information. Original `access.json` and source bytes remain under `acquired/`. Current Catalog failures (including HTTP404) remain in the separate accession audit and are not converted into assertions of global unavailability.

Europe PMC XML attempts returned HTTP500 in several cases; NCBI BioC sometimes returned HTTP200 with a 6,148-byte HTML error/landing document instead of XML. These failed scientific-payload attempts are retained. Successful primary article XML or PMC/Nature HTML was obtained separately where reported above. HTTP200 alone was never treated as proof that a scientific file or usable summary table had been received.

For bounded requests, checksums cover only the acquired prefix. The helper's original `read_truncated` records whether the client byte cap truncated its HTTP response; an HTTP206 response can itself finish normally while still representing only a remote-file prefix. The consolidated ledger therefore separately marks `remote_object_partial`. Only the four small dbGaP exports were fully acquired and explicitly checked against HTTP Content-Length; completeness of a **display export** does not imply completeness of the underlying GWAS.

The audit computed no association statistics, PIPs, credible sets, LD matrices, model scores or candidate rankings. It counted rows and checked field presence/sign strings only to classify the small public exports and release schemas. Published locus counts are source descriptions, not outcome-based locus selection. Expected eligible loci remain **unknown**, not zero and not the number of Catalog hits. The stage must stop for investigator review before any full analysis.
