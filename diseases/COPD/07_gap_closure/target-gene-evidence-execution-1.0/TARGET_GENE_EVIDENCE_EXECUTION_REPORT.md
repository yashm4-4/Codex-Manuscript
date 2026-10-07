# COPD target-gene evidence execution 1.0

**Status: execution complete; frozen for investigator review.** The 337 V1
candidates retain their original order. This stage characterizes target-gene
evidence; it makes no causal-gene declaration, reranking, model change,
manuscript edit, or de novo colocalization.

## Answers to the stop questions

1. **Exact usable GWAS/eQTL allele join:** **36/337 candidates**. The frozen
   local A/B/C signed GWAS locus ledger contains an exact, reference-checked
   signed allele for 42 candidates (41 direct-support, one EHR-support), and
   36 of those have at least one exact GTEx v10 Lung significant pair. All
   36 usable joins are in track A and the direct clinical/spirometric
   phenotype-support stratum. Six exact signed candidates lack a significant
   GTEx pair in the frozen archive. Track B and C contribute zero exact
   candidates in the *available locus ledger*. This does not mean that the
   complete original GWAS files lack those alleles.
2. **Risk-to-expression direction:** **387 candidate–gene pairs** among the
   36 candidates, one verified A source row per pair. There are 178 increased
   and 209 decreased GTEx bulk-lung normalized-expression associations per
   A COPD disease-increasing allele. These are allele-arithmetic statements;
   none implies mediation, an exact causal variant, or colocalization.
3. **Independent convergence:** **FAM13A at `4:88963935:G:T`** has an exact
   GTEx v10 Lung eQTL and the already curated rs2013701 G-to-T endogenous
   edit in 16HBE bronchial epithelial cells, with a measured FAM13A
   expression effect. A's disease-increasing T allele aligns with increased
   GTEx FAM13A expression (ALT slope `+0.132724151015`; A's G-effect beta
   `-0.0773919739743477`, SE `0.014861`). The independent edited-cell assay
   supports the same G/T expression ordering. This pair reaches **Tier B**;
   the tissue contexts and experimental units differ. No pair reaches Tier A.
4. **Locus or predictive evidence:** Saferali Table 5 reports 38
   highest-posterior Moloc window/gene/QTL records, and Table S2 lists 57
   additional gene/QTL/PPA labels. Table S1 supplies 89 tested windows. We
   uniquely identify the reported window for 34 Table 5 records using the
   supplement; four overlapping MHC windows remain ambiguous and are not
   tiered by their uncertain window assignment. Candidate overlap with a
   uniquely identified window is only locus context. The ENCODE lung
   thresholded rE2G file yields 45 predictive element–gene overlaps at 33
   candidates and 35 genes. None is physical contact or perturbation.
5. **Frozen tier hierarchy:** At the distinct candidate–gene–phenotype
   stratum level, **Tier A 0, Tier B 1, Tier C 1,216, Unresolved 3,218**
   (4,435 pairs). At candidate–gene–stratum–context row level, the counts
   are **0, 1, 1,389, and 3,611** (5,001 rows). The B row explicitly joins
   GTEx Lung and 16HBE contexts; each source-specific row retains its own
   evidence tier. Tier C often denotes a significant eQTL, published *window*
   Moloc, or predictive rE2G link, not a demonstrated candidate target.
6. **Remaining contradictions and unresolved joins:** One candidate has a
   GRCh37/GRCh38 reference mismatch across the three track checks; six A
   candidate indel rows fail the frozen GWAS identity status; and the four
   MHC Moloc window identities remain ambiguous. There are 159 candidates
   with two or more genes having typed suggestive-or-stronger support. At
   the FAM13A candidate, the published direct-COPD Moloc window highlights
   **PKD2** in COPDGene whole blood sQTL data (PPA 0.832), whereas the exact
   edited allele and lung eQTL support **FAM13A**; the former is a window
   finding, not a candidate-level contradiction to the latter. No
   sign-arithmetic discrepancy was found among the 387 reported directions.
7. **Separate `coloc.abf` input-resolution stage:** **Yes, for targeted
   input resolution only.** The exact A–GTEx joins and the FAM13A/PKD2
   locus-level competition justify auditing full regional GTEx statistics
   against a prespecified GWAS-defined A locus. No current package is
   execution-cleared for `coloc.abf`: the GTEx v10 archive has only
   significant pairs, and regional coverage, effect/variance calibration,
   MAF/N and single-signal assumptions remain unresolved. This stage did
   not run formal colocalization.

## Frozen inputs and exact-join denominators

The R010, R006, GTEx, phenotype and signed-ledger SHA256 digests match the
frozen preflight receipts. The three contracts match the preflight checksum
record. `source_acquisition_receipts.tsv` lists
the hashes of every local reference and acquired small source. The V1
GRCh38 2,001-base reference windows validate candidate REF; the GRCh37
FASTA validates mapped source REF. All 337 candidate positions have the
earlier audited unique one-base liftOver and exact reverse round trip. This
execution additionally requires the signed GWAS REF/ALT and effect/other
alleles to agree with the mapped allele pair. Palindromic orientation is
withheld if encountered. Six local indel rows fail the frozen source
harmonization status, so no indel direction is rescued by coordinate or
rsID matching.

The primary support strata remain **184 direct**, **29 additional EHR**, and
**124 ML/surrogate-only** candidates. GTEx significant pairs are observed at
185/337 candidates (915 pairs), irrespective of GWAS signed coverage.
`crosswalk_denominators.tsv` reports every stratum × A/B/C count, including
separate absent-from-GWAS-ledger and absent-from-significant-GTEx flags.
Across the 1,011 candidate × track checks there are 42 exact eligible, 960
absent from the local GWAS ledger, six otherwise rejected, and three
build/reference mismatch statuses; ambiguity, unresolved palindromes, and
duplicate conflicts are zero in this ledger. The three mismatch statuses
refer to the same candidate checked against A, B, and C. Absence from a
significant-pair file is not evidence of no cis-eQTL.

A is direct ever-smoker spirometric COPD, B is direct Japanese clinical COPD,
and C is secondary EHR J44 COPD. Their source effects and phenotypes remain
separate. A and C share UKB participants and cannot be independent
replication. Track B's ancestry and SPA likelihood limits, and C's SPA
effect/SE limits, remain as frozen in the preflight.

## Published and predictive sources

The [Saferali paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC12481885/)
reports **33** COPD GWAS loci in its summary and **32** in its results text;
both claims are retained, with 38 reported windows. Table 5 is the maximum
PPA per window, so it is not a complete list of all tested gene/feature
combinations. Table S2 adds 57 published gene/QTL/PPA labels but does not
give full feature identifiers or regional summary data for those labels.
The structured tables retain the original QTL class (eQTL, sQTL, apaQTL),
lung LTRC or whole-blood COPDGene context, molecular feature where printed,
posterior and source COPD GWAS. Their Moloc result is a published
window/gene association using the Sakornsakolpat COPD GWAS; its molecular
cohorts differ from GTEx, but it is not independent disease-association
replication. Neither a best-colocalized SNP nor a window overlap assigns
an exact frozen candidate to the gene.

[ENCODE ENCSR528UQX](https://www.encodeproject.org/annotations/ENCSR528UQX/)
is a GRCh38 DNase-based lung annotation from donor ENCDO528BHB. The selected
thresholded rE2G file is `ENCFF324XYW`, generated by
`distal-regulation-encode_re2g 1.0.0`. Its 81,638 supplied links have an
observed score range of about 0.201–1.0; the lower observed score is not
treated as a derived cutoff. The [rE2G pipeline](https://github.com/EngreitzLab/ENCODE_rE2G/blob/main/README.md)
describes thresholding at a model operating point calibrated to 70% recall
of a CRISPR benchmark. The predictions use an ABC-derived feature lineage;
we do not count rE2G and ABC as independent evidence or relabel the
`3DContact.Feature` predictor as a measured lung contact in a candidate.

The one exact functional promotion is the already curated
[rs2013701 FAM13A editing study](https://pmc.ncbi.nlm.nih.gov/articles/PMC6353020/).
The allele pair is G/T at the exact V1 GRCh38 coordinate, confirmed against
the Ensembl variation mapping; the rsID alone was not used for identity.
The remaining seven curated records remain linked, nearby, haplotype, or
otherwise unresolved relative to a frozen candidate, and provide no exact
candidate–gene promotion. There was no V1 physical-contact map to import.

## Output and validation index

- `candidate_gene_evidence_matrix.tsv` — typed candidate × gene × stratum ×
  context rows, including proximity, Catalog, eQTL, direction, published
  Moloc, rE2G, coding and contradiction fields.
- `candidate_gene_tier_summary.tsv` and
  `target_evidence_tier_assignments.tsv` — pair and context tier views.
- `candidate_gwas_gtx_exact_crosswalk.tsv`, `crosswalk_denominators.tsv`,
  `risk_expression_direction_results.tsv`, and
  `risk_expression_direction_rejections.tsv` — identity and sign ledgers.
- `saferali_moloc_crosswalk.tsv`,
  `saferali_additional_moloc_signals.tsv`,
  `saferali_supplementary_table2_audit.tsv`,
  `encode_lung_re2g_crosswalk.tsv`, and
  `functional_target_evidence_audit.tsv` — external evidence with source
  scope preserved.
- `phenotype_stratum_target_summary.tsv` and
  `contradiction_unresolved_ledger.tsv` — stratum totals, competing genes,
  failed joins and unresolved assignments.
- `execution_contract.md`, `source_acquisition_receipts.tsv`,
  `validation_results.json`, `validation_audit.json`,
  `synthetic_validation_records.tsv`, `prior_freeze_integrity_audit.json`,
  scripts, and `provenance/freeze.json`.

The independent validator rechecked frozen order, raw signed GWAS beta/SE,
GTEx ALT slopes, all 387 allele signs, the reported intervals, source type
firewalls, tier rules and phenotype strata. Synthetic same-ALT, REF-reversal,
negative-slope, complement, palindrome and indel-mismatch cases passed.
The checksum manifest records the exact stage payload. Large previously
acquired resources remain outside Git; only the small Saferali and ENCODE
source files were downloaded to the ignored, reconstructable `sources/`
directory. No older frozen payload was changed. The historical preflight
freeze verifier passed before these authorized register additions. It binds
the *old register endings* and therefore fails its strict register-hash check
after the append-only rows. `prior_freeze_integrity_audit.json` separately
verifies all 15 preflight payload hashes and the exact old register byte
prefixes; only the six new register rows differ from the base commit.
