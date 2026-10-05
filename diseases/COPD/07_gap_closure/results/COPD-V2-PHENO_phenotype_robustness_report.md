# COPD V2: frozen GWAS phenotype-support robustness

Analysis/report ID: **COPD-V2-PHENO-REPORT-001**. Date: 2026-10-05. Scope: the investigator-authorized phenotype module only. V1, its manuscript, thresholds, models, candidate order, and experimental shortlist are unchanged.

## 1. Main conclusion

Phenotype heterogeneity is consequential for candidate-specific genetic interpretation, not merely a wording caveat. Direct clinically/spirometrically defined COPD susceptibility supports **184/337 frozen candidates (54.60%)**, **72/153 frozen redundancy components (47.06%)**, and **6/12 experimental-shortlist entries**. Another 29 candidates have EHR COPD support without direct-phenotype support; the remaining 124 are supported exclusively by the ML-liability accession GCST90244098. The original highest-ranked two candidates are in that ML-only group.

This does **not** establish that ML-supported candidates are biologically false. The analysis tests the phenotype provenance of their significant genetic support, not the validity of their regulatory sequence predictions or the existence of COPD biology. The noncoding predominance and a substantial directly supported candidate core survive. The claim that all 337 share equivalent direct COPD-susceptibility support does not.

Primary numerical sources: [stratum summary, R008](COPD-V2-PHENO-R008_stratum_summary.tsv), [all 337 records, R005](COPD-V2-PHENO-R005_frozen_337_phenotype_support.tsv), [GCST90244098 dependence, R009](COPD-V2-PHENO-R009_GCST90244098_dependence.tsv). Counts describe frozen V1 ascertainment, not all evidence that might exist in complete modern GWAS statistics.

## 2. Frozen inputs and prespecified method

### Input provenance

The core universe is the 104 accessions/42 publications in [V1 S2-R001](../../02_gwas/results/COPD-S2-R001_studies_core.tsv). The V1 Catalog freeze is 15 September 2026. Significant study–tag support comes from the 827 rows in [S2-R002 GWS associations](../../02_gwas/results/COPD-S2-R002_gws_associations.tsv), at the unchanged P ≤ 5×10⁻⁸ threshold, covering 660 tag identifiers. Only 35 of the 104 accessions contribute at least one such row. Zero curated GWS contribution is **not a negative full GWAS**.

Candidate identities and LD are reused from S2-R006A/B/C/E, including [focal-panel audit](../../02_gwas/results/COPD-S2-R006B_focal_panel_audit.tsv), [LD links](../../02_gwas/results/COPD-S2-R006C_ld_pairs.tsv.gz), and [candidate union](../../02_gwas/results/COPD-S2-R006E_candidate_variants_grch38.tsv.gz). The downstream records are the original [337 R010 candidates](../../04_modeling/results/COPD-S4-R010_THE_LIST.tsv), [S6 component/selection audit](../../06_experimental_validation/results/COPD-S6-R002_selection_audit.tsv.gz), and [12-candidate shortlist](../../06_experimental_validation/results/COPD-S6-R001_candidate_shortlist.tsv). Consumed files and SHA-256 values are in the [run manifest](../provenance/COPD-V2-PHENO_run_manifest.json).

### Adjudication and missingness

Each accession was classified from its actual tested outcome and ascertainment, not its publication title. The archived Catalog study/ancestry records were reconciled with primary publications and accessible supplementary evidence. Publication batches were adjudicated before downstream retention counts were inspected. The [specification](../provenance/COPD-V2-PHENO_analysis_specification.md), three reviewer-input JSON files, and explicit [amendments](../data/phenotype_adjudication_amendments.json) make judgments auditable. Source URLs, PMID/DOI, checked date, evidence paraphrases and uncertainties are provided per accession in [R001](COPD-V2-PHENO-R001_accession_adjudication.tsv) and the [source ledger](COPD-V2-PHENO-sources.tsv).

Primary eligibility requires COPD cases versus non-COPD controls with direct clinical or spirometric ascertainment. Severity-selected susceptibility and ordinary smoking-stratified susceptibility remain eligible. Different spirometric thresholds, pre/post-bronchodilator protocols, control screening, and cohort recruitment are recorded, not assumed identical. Secondary eligibility is ICD/PheCode/comparable EHR COPD susceptibility. Documented EHR-plus-self-report composites are separate, with an inclusive secondary sensitivity. Questionnaire-only or insufficiently resolved definitions do not become directly adjudicated COPD by assumption.

`not_reported` is not global unavailability; `unclear` is not a negative result; `not_applicable` denotes an inapplicable design field. Negative replication is retained where published. Some supplement endpoints were inaccessible: the precise uncertainty is in R001, not silently filled in. Catalog “full statistics” means a reported release for the study's tested units, not necessarily a complete common-SNP scan. Publication-wide availability does not prove a sex-specific accession file exists. Only one small 2,048-byte statistics prefix plus metadata was inspected for Kim's coefficient interpretation; no full GWAS file was downloaded for this analysis.

### Exact support reconstruction

For a directly represented tag, retain its significant study support even when reference matching or LD expansion failed. For an LD proxy, expand `supporting_studies` in R006B; select eligible study × focal-tag × ancestry-panel rows; then join the surviving tag–panel pairs to R006C. R006C has no study identity. This yielded 1,469 study–tag–panel rows over 1,245 original tag–panel combinations, covering all 660 tags; 33,378 frozen nonself LD links were reused. No LD was recalculated or borrowed from an excluded study's panel. [R003](COPD-V2-PHENO-R003_study_tag_ancestry_support.tsv) and [R004](COPD-V2-PHENO-R004_all_candidate_phenotype_support.tsv.gz) expose the routes.

A record is retained if any valid direct-tag or proxy route survives. Components use the **existing** S6 memberships and any-member support, with all/some/no-member states; components are not reconstructed. They remain shared-tag redundancy heuristics, not formal independent LD blocks. Original candidate ranks are copied unchanged; records outside R010 have rank `not_applicable`.

The phenotype sensitivity is an intersection with the frozen candidate list. It does not rescore variants, recompute the broad-universe delta percentiles, recalibrate thresholds, or estimate a new false-discovery rate. Study counts, correlated tags, and components are not treated as independent observations. All comparisons are descriptive; there are no new enrichment P values or causal probabilities.

## 3. Accession adjudication results

| Class | Accessions | Accessions with ≥1 V1 GWS association |
|---|---:|---:|
| Direct COPD susceptibility | 27 | 19 |
| EHR COPD susceptibility | 41 | 10 |
| ML spirogram liability | 1 | 1 |
| Chronic bronchitis | 6 | 2 |
| Conditional smoking interaction-model coefficients | 2 | 2 |
| Documented EHR/self-report composite susceptibility | 2 | 0 |
| Exacerbation/infection subtype endpoints | 2 | 0 |
| Gene/burden/CNV designs | 16 | 0 |
| Progression/decline | 1 | 0 |
| Within-COPD severity | 1 | 0 |
| Unresolved under the requested hierarchy | 5 | 1 |
| Total | 104 | 35 |

The title-level heterogeneity noted during planning overstated some exclusions. GCST002624 is COPD **without chronic bronchitis versus normal controls**, not the paper's within-case bronchitis contrast; GCST002795 similarly concerns COPD without pulmonary artery enlargement versus normal controls. GCST005417 tests binary post-bronchodilator COPD despite a lung-function-focused title. These remain primary with explicit attributes. The exact subtype contrasts were reconciled through sample counts and reported associations; inaccessible S8/E4 files are disclosed rather than claimed retrieved. See their R001/source-ledger rows and the primary [Lee chronic-bronchitis paper](https://doi.org/10.1186/s12931-014-0113-2) and [pulmonary-artery study](https://doi.org/10.1165/rcmb.2014-0210OC).

Joo (GCST90103984/985) and Moll (GCST011766) contribute individual-variant susceptibility tests even though their publications emphasize gene-based/exome approaches. Exome ascertainment alone is not a burden test. Conversely, actual burden/CNV accessions remain separate. These distinctions materially preserve direct support: excluding Moll alone would lose 20 of the 184 candidates; excluding Joo's two-accession publication would lose 10. These are classification-influence demonstrations, not recommendations to exclude either publication. [R015](COPD-V2-PHENO-R015_classification_influence.tsv), [R016](COPD-V2-PHENO-R016_publication_influence.tsv).

Four ordinary smoking-stratified Kim susceptibility accessions remain primary. Two interaction-model main-effect coefficients (GCST90016586/591) remain separate because they condition on reference exposure rather than represent the pooled marginal comparison. Their inclusion bound changes **none** of the five primary retained universes: their relevant routes already have primary support. Excluding the four ordinary smoking strata together would lose 12 candidates, so automatically excluding every smoking-stratified analysis would be consequential and contrary to the investigator's hierarchy. See [Kim et al.](https://doi.org/10.1093/aje/kwaa227), R001 and R016.

Ishigaki's original physician-diagnosed BBJ target recruitment is accepted as direct clinical susceptibility, with moderate confidence and unreported exact COPD spirometry; Sakaue's expanded medical-history/text-mined EHR definition stays secondary. Excluding the three Ishigaki accessions removes only 6 tags, 60 records, 2 candidates and 2 components; no shortlist entry changes. Those candidates are frozen ranks 130 (`19:40718299:G:C`) and 216 (`19:40908576:G:A`). The all-high-confidence-primary bound is therefore 219 tags, 8,152 records, 182 candidates, 70 components and 6 shortlist entries. Mathur's unscreened COGEND-control caveat and Do's broader replication ascertainment are recorded; neither uniquely retains an R010 candidate. [Ishigaki](https://doi.org/10.1038/s41588-020-0640-3), [Sakaue](https://doi.org/10.1038/s41588-021-00931-x), R001/R015.

Unresolved GCST007996 contributes one unmatched tag/record and no R010 candidate. Adding it to secondary EHR changes 96→97 tags and 1,977→1,978 records, but leaves 31 candidates, 15 components and 2 shortlist entries unchanged. The other four unresolved accessions contribute no V1 GWS row. Newbury ADO and Nam's documented EHR/self-report composites also contribute no V1 GWS row; including them with EHR changes no reported support count. Thus these uncertainties cannot explain the main 337-candidate result. The [Newbury algorithm paper](https://doi.org/10.1038/s41746-025-01815-8) and accession-specific R001 rows establish the distinction between known composite ascertainment and unresolved exact definitions.

## 4. Retention across frozen universes

| Frozen universe | Broad V1 | Direct primary | EHR secondary | Primary OR EHR | ML-liability support |
|---|---:|---:|---:|---:|---:|
| GWS tag identifiers | 660 | 225 (34.09%) | 96 (14.55%) | 315 (47.73%) | 356 (53.94%) |
| Tag/proxy records | 15,389 | 8,212 (53.36%) | 1,977 (12.85%) | 9,906 (64.37%) | 7,240 (47.05%) |
| Frozen regulatory candidates | 337 | 184 (54.60%) | 31 (9.20%) | 213 (63.20%) | 173 (51.34%) |
| Frozen components | 153 | 72 (47.06%) | 15 (9.80%) | 85 (55.56%) | 88 (57.52%) |
| Experimental shortlist | 12 | 6 (50.00%) | 2 (16.67%) | 8 (66.67%) | 6 (50.00%) |

Columns overlap and must not be added. These are exact GWS/tag/panel counts from R008, not aggregate linked-study annotations.

Among the 337, exclusive display states are 184 primary-retained; 28 EHR-only; 124 ML-only; and 1 EHR+ML without primary support (rank 244, `3:169003792:T:C`, GCST90399695 and GCST90244098). Two of the 31 EHR-supported candidates also have primary support. All 72 primary-supported components retain **all** their frozen members; 81 retain none and none is partially retained. This does not make the 72 independent signals.

All 85 originally unresolved reference-match tags are explicitly retained in [R013](COPD-V2-PHENO-R013_unresolved_tag_support.tsv): 12 primary, 31 EHR-only, 39 ML-only, 2 other-only, and 1 unresolved phenotype. Their match-failure statuses remain frozen. “Unresolved tag identity” and “unresolved study phenotype” are separate fields and separate limitations.

## 5. GCST90244098 dependence and planning-count reconciliation

The primary paper establishes that GCST90244098 tests a continuous model-derived liability from raw spirograms, not binary clinically adjudicated COPD. Record/self-report labels trained the model, but those labels are not the GWAS outcome. This remains COPD-relevant genetic evidence with a distinct estimand. [Cosentino et al., 2023](https://doi.org/10.1038/s41588-023-01372-4).

| Universe | ML-linked, original aggregate provenance | ML-linked, exact support | ML sole accession, exact support | Exact ML-linked also primary | Exact ML-linked not primary |
|---|---:|---:|---:|---:|---:|
| 660 tags | 356 | 356 | 339 | 16 | 340 |
| 15,389 records | 7,446 | 7,240 | 5,470 | 1,725 | 5,515 |
| 337 candidates | 177 | 173 | 124 | 48 | 125 |
| 153 components | 90 | 88 | 68 | 19 | 69 |
| 12 shortlist entries | 6 | 6 | 4 | 2 | 4 |

The planning observation is independently reproduced **under its original aggregate definition**: 356 tags, 337 exclusive tags, 177 linked R010 candidates, and 124 sole-accession candidates. Exact study-specific significance changes exclusive tags from 337 to 339 because the aggregate table includes nonsignificant reports for rs11594905 and rs12894780. In total, 37 tags have at least one such extra nonsignificant accession link; [R012](COPD-V2-PHENO-R012_aggregate_nonsignificant_provenance_audit.tsv) enumerates them. This is a provenance distinction, not a changed significance threshold or deletion from V1.

The 177→173 candidate difference is ancestry-panel provenance, not phenotype reclassification. Frozen ranks 62, 99, 138 and 193 have EAS-only proxy links to their relevant tags. GCST007692 supports those EAS tag–panel pairs, whereas GCST90244098 supports EUR, with no corresponding EUR proxy link to these four candidates. All four remain primary-retained. Original links remain in R005's `v1_aggregate_linked_studies`; exact routes are separate. The naive alternative of retaining every panel for a surviving primary tag would falsely retain four records outside R010; the implemented panel filter avoids that error.

None of the 124 sole-ML candidates retains primary support; 125 ML-linked candidates lack primary support because rank 244 has EHR as well as ML support. These results justify a phenotype-dependent genetic label, **not** a biological-negative label, removal from V1, or retroactive demotion of the frozen ranking.

## 6. Which candidates and loci are robust?

The full answer is the unchanged-order R005 table, not a replacement top list. Illustrative directly supported candidates include the three frozen FAM13A assignments at ranks 111 (`4:88964563:T:C`), 210 (`4:88963935:G:T`, the V1 rs2013701 record), and 310 (`4:88965146:C:T`); RIN3 ranks/records in R005; and inherited AGER, CHRNA3/CHRNA5 and IREB2 locus assignments. These are phenotype-support statements, not fresh functional validation or effector-gene inference.

Locus denominators must be kept distinct. V1's manuscript reports **85/140** candidate-bearing repeated-Catalog-mapping gene loci. R010 contains **88/152** all-selected locus labels, because VPS54, PELI1 and ALDH2 were additionally selected by reported effect size. There is no contradictory V1 count. Under primary support, **72/85** originally candidate-bearing repeated-mapping loci remain represented (72/140 in the fixed full set); the all-selected comparison is **74/88** (74/152 overall). These counts and selection flags are in [R011](COPD-V2-PHENO-R011_frozen_locus_support.tsv), using [S3-R002 definitions](../../03_regulatory_landscape/results/COPD-S3-R002_gene_loci_grch38.tsv) and the [V1 R004 locus summary](../../04_modeling/results/COPD-S4-R004_locus_summary.tsv). Locus assignment propagates existing tag-gene mappings or coordinate windows; it is not a new target-gene analysis, causal locus definition, or independent replication test.

### Unchanged experimental-shortlist membership

| Original shortlist position | Original R010 rank | Candidate GRCh38 ID | Exact phenotype support |
|---|---:|---|---|
| 1 | 1 | 11:62567436:G:C | ML only |
| 2 | 2 | 1:3528722:G:C | ML only |
| 3 | 3 | 14:92637384:G:A | Primary + ML |
| 4 | 4 | 6:31009903:G:A | Primary |
| 5 | 5 | 6:27556090:G:A | Primary |
| 6 | 6 | 16:75478398:G:GC | Primary |
| 7 | 7 | 15:67322629:G:A | EHR only |
| 8 | 8 | 16:28602644:A:G | Primary |
| 9 | 9 | 17:40058327:G:GCCCAGAC | EHR only |
| 10 | 38 | 6:4577675:T:A | ML only |
| 11 | 39 | 15:67150258:C:T | ML only |
| 12 | 77 | 11:13140768:T:C | Primary + ML + conditional smoking coefficient |

Source: [R007](COPD-V2-PHENO-R007_shortlist_phenotype_support.tsv). No candidate was replaced, reranked or experimentally tested. A future study may deliberately sample across these strata, but this module does not authorize redesign of the panel.

## 7. Cohort overlap and robustness limits

[R002](COPD-V2-PHENO-R002_cohort_overlap_provenance.tsv) records cohort provenance and [R002B](COPD-V2-PHENO-R002B_cohort_overlap_pairs.tsv) identifies shared-cohort accession pairs. `BEOCOPD` is normalized to `EOCOPD`; GBMI's `EB` to `EstBB`, with original labels retained. ICGC, CHARGE and UKECC are consortium umbrellas, not additional independent cohorts. ICGC is not ICGN. Discovery, replication and contextual-follow-up roles are described in the accession register; a shared label does not establish identical individuals or positive independent replication. Some consortium lists are documented overlap anchors rather than exhaustive rosters, so absence from R002B does not establish independence.

Both the direct and ML strata reuse UKB; COPDGene and other classical cohorts also recur. Thus primary+ML support is cross-phenotype support, **not necessarily independent participant replication**. Leave-one-publication sensitivity shows that removing GCST007692 would remove 69 of the 184 primary-supported candidates, 28 components and one shortlist entry. Moll uniquely supports another 20 candidates; Kim's ordinary smoking strata 12; Joo's publication 10. The overlapping cohorts behind these publications preclude interpreting these figures as counts of independent discoveries. [R016](COPD-V2-PHENO-R016_publication_influence.tsv).

This analysis inherits V1's Catalog significance censoring, tag normalization, ancestry descriptors, panel assignment and reference matching. It does not revisit whether every original ancestry assignment or rare/exonic LD estimate is optimal. Complete summary statistics might establish support absent from curated GWS hits; failure to retain is not a nonreplication test. Power, sample size, LD structure and phenotype ascertainment can all affect support breadth. The observed difference cannot be causally attributed to diagnostic heterogeneity alone, and the ML study's wider discovery is not itself evidence of bias.

## 8. Biological conclusions and future manuscript qualifications

The broad regulatory rationale survives, but its numerical and phenotype scope narrow. Among direct-primary tags, 187/220 consequence-annotated tags are noncoding (**85.00%**; 5/225 unknown), versus V1's 533/582 (**91.58%**; 78/660 unknown). Splice annotations remain a separate category in both calculations. Among 184 retained candidates, 98 retain the frozen enhancer call, 107 the frozen silencer call, and 21 both; 7 overlap coding CDS, 28 a donor-refined region and 130 a known regulatory catalog. These are descriptions of unchanged V1 annotations, not validation or new mechanism assignments. [R014](COPD-V2-PHENO-R014_frozen_biological_feature_summary.tsv).

The principal scientific defect exposed is **semantic overbreadth of the genetic phenotype claim**, accompanied by aggregate provenance that cannot be used directly as ancestry-specific significant support. It is substantial for prioritization: 45.40% of R010 lacks direct-primary support, half the proposed experimental panel does, and the top two original ranks are ML-only. It does not invalidate the broad-universe sequence computation or demonstrate that its predictions fail. V1's mechanistic literature, cell-context caveats, lack of colocalization/fine-mapping, and proposed-not-performed experimental status are not overturned or newly validated by this analysis.

The following future V2 manuscript qualifications are warranted; **the V1 manuscript has not been edited**:

| V1 location/statement | Qualification supported by this module |
|---|---|
| Title, abstract and opening Results framing as COPD susceptibility loci | Identify an ontology-defined, heterogeneous COPD-associated discovery universe; distinguish direct, EHR and ML-liability support. |
| GWAS Methods: ontology scope followed by manual phenotype review | Cite the new accession-level adjudication as V2 work; do not imply this register existed in V1. |
| Abstract/Results: 660 tags, 15,389 records, 337 candidates | Preserve those V1 counts and add the direct-primary intersections 225, 8,212 and 184, without redefining the frozen list. |
| Discussion: >91% noncoding significant tags | Restrict 91.58% to the broad universe; direct-primary annotated tags yield 85.00%. Neither percentage proves regulatory causation. |
| 85/140 candidate-bearing repeated-mapping gene loci | Add 72/140 under primary support; do not substitute the 152-label expanded universe or call repeated accessions independent replication. |
| Figure 3 and experimental-panel narrative | Preserve the original selection and identify six primary, two EHR-only and four ML-only entries. |
| Repeated gene mapping, conservation, chromosome clustering, target annotations and 185 GTEx matches | Keep their broad-universe scope. Those tests and mappings were not rerun or revalidated in a primary-only universe; exact eQTL overlap still is not colocalization. |

Source for the statements being qualified: [frozen V1 manuscript](../../manuscript/COPD_regulatory_genomics_manuscript.md), together with result IDs linked above. No disease risk direction is inferred from REF/ALT, allele frequency or unsigned r².

## 9. Reproducibility and validation

Run from the repository root:

```bash
python3 diseases/COPD/07_gap_closure/scripts/run_phenotype_robustness.py
python3 diseases/COPD/07_gap_closure/scripts/validate_phenotype_module.py
```

The pipeline is Python standard-library only, offline after adjudication. It uses one CPU, no GPU or scheduler job, approximately 3 seconds and under 0.5 GiB peak resident memory in the recorded execution; the publication audit, not computation, dominated effort. Exact software/platform/runtime details are in the run manifest and [run log](../logs/COPD-V2-PHENO_run_log.json).

The [53 analysis checks](COPD-V2-PHENO-validation_checks.tsv) cover 104-accession completeness, allowed classes, eligibility, GWS threshold, valid matched-polymorphic LD routes, original universe recovery, unchanged R010 order, exact component memberships, unchanged shortlist, all 85 unmatched tags, independent study-first reconstruction for every stratum, input hashes and TSV round-trip rectangularity. A separate reviewer reproduced the headline counts and diagnosed the four ML ancestry-panel differences. The [final validator](../scripts/validate_phenotype_module.py) checks the original freeze hashes and tracked V1 files against commit `875e995406ddf5d00ce77920d9e650ff128604c0`, analytical checksums and per-record summary agreement. Results are in [final validation](../provenance/COPD-V2-PHENO_final_validation.tsv); [artifact checksums](../provenance/COPD-V2-PHENO_artifact_checksums.tsv) cover the completed module. Structural validation is not proof that every literature judgment is biologically definitive.

## 10. Recommendation and stop

Phenotype heterogeneity has proved consequential. Retain the broad V1 record, carry explicit primary/EHR/ML support labels into any later V2 evidence comparison, and avoid presenting all candidates as equally direct-COPD-supported. The only moderate-confidence primary judgment affecting R010 is original BBJ clinical ascertainment (two candidates); the major loss is not driven by unresolved accessions or discretionary interaction exclusions.

This result does **not** justify immediate multi-context retraining, genome-wide fine-mapping, or a comprehensive smoking project. It strengthens the reason to conduct the already proposed **cheap frozen-model reverse-complement robustness check next**, followed by an external functional benchmark with phenotype/mechanism strata. Those checks address whether the sequence predictions themselves are stable and informative; they should include relevant direct and ML-only examples without tuning on positive controls. Additional context modeling remains conditional on those findings and coherent input data. This is a recommendation only.

**STOP: no reverse-complement scoring, functional benchmark, retraining, new target-gene analysis, fine-mapping, large download, manuscript revision, candidate reranking, or GitHub push has been performed by this module. Investigator review is required before another V2 module.**
