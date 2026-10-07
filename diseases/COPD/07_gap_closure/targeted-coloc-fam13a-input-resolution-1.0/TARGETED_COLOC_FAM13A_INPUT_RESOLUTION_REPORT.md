# FAM13A targeted coloc.abf input resolution — frozen review report

**Final readiness: NOT CLEARED. No colocalization was run.** This is a bounded input-resolution stage under the published target-gene execution freeze `0adeb7b90b762f15d2c6f940bf361c0300978ce3`. The 337 candidates, FAM13A Tier B, all prior freezes, Track-A GWAS statistics, and manuscript remain unchanged. The locus is the already frozen GWAS-defined **GRCh37 chr4:88,300,892–91,560,531**, inclusive. `4:88963935:G:T` is only a previously known biological reference; it was never used to define boundaries, variants, filters, priors, single-signal judgment, or clearance.

## Stop-condition answers

| Question | Finding |
|---|---|
| Compatible, sufficiently complete statistics? | **No.** The QTL has beta/SE/MAF and tested nulls, but GWAS ABF precision/N and significant-variant coverage are unresolved. The pinned v8 QTL has no robust FAM13A eQTL signal. |
| Source and intersection counts? | 12,094 native Track-A locus rows; 7,264 QTL source rows = 7,010 distinct tested variants after 254 consistent rsID alias rows; 5,915 coordinate/REF/ALT candidate GWAS–QTL row matches; **5,065 fully eligible exact variants** (one GWAS row per QTL variant). |
| GWAS beta, variance, N and case fraction? | ln(OR) and source SE are supported by the submission/PLINK convention but not certified for every row; 32 regional OR/SE/P precision exceptions, including 15 exact QTL overlaps. Nominal 12,446/59,145 gives N 71,591 and case fraction 0.1738486681; per-variant N is absent. **Not defensible as a complete ABF likelihood.** |
| GTEx beta, variance, N and MAF? | Quantitative ALT beta, supplied SE, nominal N 510, per-variant AN/2 = 485–510, and MAF are available. Variant-level r² is entirely NA. These fields support a QTL statistics contract in principle, with the AN variation requiring a rule before scalar-N use. |
| Single-signal assumption? | **INCONCLUSIVE.** No validated conditional Track-A decomposition; 50 of 51 significant GWAS rows cluster in one 100-kb bin and one is separated. The pinned v8 FAM13A slice has no official credible set, so its signal count is not inferable as one. |
| Could exclusions bias the comparison? | **Yes.** Only 5,065/12,094 GWAS rows and 5,065/7,010 QTL variants pass strict exact matching. Twelve of 51 significant GWAS rows, including the strongest source lead, are absent. This can materially distort regional ABFs or H3/H4 comparison. |
| Decision? | **NOT CLEARED**, freeze and stop for investigator review. |

## Gate audit

| Gate | State | Evidence |
|---|---|---|
| GWAS statistical contract | FAIL | 32 precision exceptions; 15 overlapping; genuine rowwise SE/N not independently certified. |
| QTL statistical contract | PASS WITH LIMITS | Full `.all` cis rows, beta/SE/nominal P/MAF present; 99 variants have AN below 1020 and r² is NA. |
| Cross-build identity | PARTIAL | 11,331 provider rows pass unique forward/reverse UCSC chain coordinates and independent hg38 REF checks; 676 GWAS rows lack a provider map, 87 provider rows are unoriented. |
| Allele harmonization | PARTIAL | 5,065 strict exact rows; 1,625 palindromic rows rejected without GWAS AF; 49 same-position/different-ALT situations are not coerced into matches. No rsID-only or proxy joins. |
| Regional overlap/coverage | FAIL | 5,065/12,094 native GWAS rows and 5,065/7,010 distinct QTL variants; 39/51 significant GWAS rows. |
| Required N/MAF information | FAIL FOR GWAS | GWAS has only nominal N and no AF or rowwise imputation QC. QTL has AN/2 and MAF. |
| Single-signal plausibility | INCONCLUSIVE | GWAS spatial profile is not conditional evidence; no FAM13A official credible set in the pinned v8 release. |
| Source provenance | PASS WITH VERSION LIMIT | All local inputs and targeted source artifacts have hashes. QTL association records give unversioned `ENSG00000138640`; exact Ensembl version is absent. Dataset ID was reused for v10 in pre-release 8, so file date/release are pinned. |

Rejection reasons in `regional_coverage_summary.tsv` are **nonexclusive**; a row may fail both identity and overlap. There are 5,416 GWAS rows with no exact QTL allele, 1,732 with unresolved native identity, 1,625 unresolved palindromes, 676 without unique provider mapping, 87 unoriented provider rows, and 32 precision failures. Direct REF/ALT conflicts and chain round-trip conflicts are both zero among interpretable records; 49 rows share a position with a different QTL ALT and are not forced into matches. The QTL has no REF conflicts or inconsistent rsID aliases. Both source denominators and every rowwise state are in the crosswalk and rejection table.

## Dataset identity and scientific interpretation

The [Catalogue release 7 metadata](https://github.com/eQTL-Catalogue/eQTL-Catalogue-resources/blob/master/data_tables/dataset_metadata_r7.tsv) identifies `QTS000015/QTD000271` as uniformly reprocessed **GTEx v8 Lung, n=510**, and the retrieved `.all` file is dated 2023-04-06. The [current release notes](https://www.ebi.ac.uk/eqtl/Release_notes/) and pre-release 8 metadata reuse that accession for **GTEx v10 Lung, n=601**. The pinned file date, metadata and hashes prevent accidental release substitution. The [Catalogue data-access page](https://www.ebi.ac.uk/eqtl/Data_access/) establishes that `.all` includes complete cis statistics, ALT is the effect allele, and rsID aliases require deduplication. The official permutation record reports 7,010 tests, smallest nominal P ≈ 5.59×10⁻⁴, empirical P ≈ 0.598 and adjusted P ≈ 0.605; the downloaded official credible-set file contains zero FAM13A entries. This v8 absence should not be read as a refutation of the separately frozen GTEx v10 exact significant pair or edited G→T allele, and those observations did not change any gate.

The future ABF specification is saved solely as a prespecified, dormant contract. Its minimum 90% reciprocal regional coverage and complete significant-GWAS capture are not met. No posterior probabilities, fine mapping, regulatory rescoring, candidate reranking, or manuscript edits occurred. The frozen stage ends here for investigator review.
