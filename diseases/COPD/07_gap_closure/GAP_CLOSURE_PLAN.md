# COPD V2 Gap-Closure Plan

**Plan ID:** COPD-V2-RES-PLAN-001
**Prepared:** 2026-10-05
**Scope:** scientific audit and planning only
**V1 snapshot:** git commit 875e995406ddf5d00ce77920d9e650ff128604c0
**Execution state:** no V2 model training, large download, genome-wide fine-mapping, candidate reranking, wet-lab experiment, or manuscript rewrite has been performed

This plan treats every COPD artifact outside 07_gap_closure as the immutable V1 record. The V1 final table remains 337 records in its original order. V2 may read V1 files and compare new evidence with them, but it must write only under 07_gap_closure, use COPD-V2 identifiers, and never present a V2 view as the original ranking or as a posterior probability unless a valid statistical model produced that probability.

Recommendation labels used below are:

- **GO:** scientifically useful with prerequisites already adequate for a planning-stage recommendation.
- **CONDITIONAL GO:** useful only after a named data or validity gate passes.
- **NO-GO:** do not execute in the proposed form.
- **Evidence type:** independent functional evidence, sensitivity evidence, association evidence, predictive annotation, physical-contact evidence, or synthesis. These types are not interchangeable.

The audit covered the COPD top-level registers; all Section 1–6 READMEs and integrated reports; important result tables; GWAS extraction, ancestry assignment, and LD code; regulatory construction; TREDNet training/evaluation/scoring/prioritization; DeepFootprinting and TFBS logic; target mapping; computational validation; experimental design; the final list; the manuscript, figure source data, rendering code, result-to-manuscript map, software records, section checks, and final 81-check audit. The exact frozen hashes are in [v1_snapshot.tsv](provenance/v1_snapshot.tsv); lightweight audit calculations are in [planning_diagnostics.tsv](provenance/planning_diagnostics.tsv); and checked external resources are in [external_resource_audit.tsv](provenance/external_resource_audit.tsv).

## 1. V1 scientific reconstruction

### 1.1 Central V1 question and claim

The central question was whether public COPD genetic associations could be converted into an auditable set of sequence-level regulatory hypotheses and prospective experiments while preserving phenotype, LD, regulatory, modeling, QTL, and experimental evidence boundaries.

The defensible central claim is narrower than causal discovery:

> V1 produced a transparent, reproducible framework that connects a broad, ontology-defined COPD GWAS universe to 337 enhancer/silencer-model hypotheses and an experimental design package; it did not establish 337 causal variants, causal genes, or COPD mechanisms.

This is the meaning recorded in COPD-DEC-029 and in the Discussion and limitations of the [V1 manuscript](../manuscript/COPD_regulatory_genomics_manuscript.md). It is also consistent with the rank boundary in [COPD-S4-R010_THE_LIST.tsv](../04_modeling/results/COPD-S4-R010_THE_LIST.tsv), which calls the order computational prioritization rather than causal proof.

### 1.2 Major analytical steps and exact chain of evidence

1. **Disease framing and source curation.** Section 1 distinguished COPD susceptibility from lung function, emphysema, chronic bronchitis, smoking, and related lung disease in prose and source records. Evidence: [Section 1 README](../01_background/README.md), [Section 1 report](../01_background/results/section_1_background.md), [sources.tsv](../sources.tsv), and [decisions.tsv](../decisions.tsv).

2. **Ontology-defined GWAS extraction.** The core set required primary mapped URI MONDO:0005002 or legacy EFO:0000341. This yielded 104 study accessions from 42 publications, 995 association rows, and 760 normalized tag variants. At P ≤ 5×10^-8, V1 retained 827 association rows and 660 unique tags. Evidence: [extraction script](../02_gwas/scripts/01_extract_copd_gwas.py), [core study table](../02_gwas/results/COPD-S2-R001_studies_core.tsv), [GWS associations](../02_gwas/results/COPD-S2-R002_gws_associations.tsv), [GWS tags](../02_gwas/results/COPD-S2-R002_gws_unique_tag_variants.tsv), and [Section 2 report](../02_gwas/results/section_2_gwas.md).

3. **Coding/noncoding annotation.** Of 582 consequence-resolved tags, 533 were noncoding, 46 coding, and three splice; 78 tags remained consequence-unknown. The strongest descriptive conclusion is therefore that reported COPD tags were predominantly noncoding, not that noncoding tags were causal. Evidence: [variant consequences](../02_gwas/results/COPD-S2-R004_variant_consequences.tsv) and [coding/noncoding summary](../02_gwas/results/COPD-S2-R004_coding_noncoding_summary.tsv).

4. **Ancestry-aware LD expansion.** Study-supported ancestry categories selected 1000 Genomes Phase 3 20190312 GRCh38 panels without an automatic European fallback. Tags were matched by coordinate and reported allele where available, with explicitly audited limited fallbacks. PLINK 1.9 expanded within ±500 kb at r² ≥ 0.8. Of 660 GWS tags, 575 matched at least one panel and 85 remained unresolved; 33,378 nonself tag-panel links produced 14,913 unique proxies. The tag-plus-proxy union was 15,389 records. This was a careful candidate-expansion procedure, not statistical fine-mapping. Evidence: [LD preparation](../02_gwas/scripts/03_prepare_ld.py), [LD execution](../02_gwas/scripts/04_run_ld.py), [LD summarization](../02_gwas/scripts/05_summarize_ld.py), [focal audit](../02_gwas/results/COPD-S2-R006A_focal_match_audit.tsv), [panel audit](../02_gwas/results/COPD-S2-R006B_focal_panel_audit.tsv), and [LD summary](../02_gwas/results/COPD-S2-R006F_summary.tsv).

5. **Regulatory landscape.** Disease-matched peaks came from three lobes of one ENCODE donor, ENCDO520EJG: a 60-year-old European male described as having severe emphysema, without repository evidence of clinically adjudicated COPD or spirometry. V1 used ATAC plus H3K27ac for enhancer-like regions and ATAC plus H3K27me3 for silencer-like/repressed regions, alongside public regulatory annotations, repeats, blacklist, and CDS. Section 3 correctly made refined annotation overlap same-lobe before union. Evidence: [Section 3 README](../03_regulatory_landscape/README.md), [ENCODE selection manifest](../03_regulatory_landscape/data/encode_peaks/manifest.tsv), [regulatory construction](../03_regulatory_landscape/scripts/03_gene_locus_regulatory_burden.py), [candidate classification](../03_regulatory_landscape/scripts/05_classify_candidate_variants.py), and [Section 3 report](../03_regulatory_landscape/results/section_3_regulatory_landscape.md).

6. **TREDNet training and internal evaluation.** Separate enhancer and silencer phase-II models used the frozen 565-MB phase-I representation. Phase-II chromosome partitions were train chr1–6,10–22,X,Y; validation chr7; and test chr8–9, with one deterministic seed. The test ROC AUC/PR AUC values were 0.938/0.885 for enhancer and 0.967/0.938 for silencer, with Brier scores 0.091 and 0.061. The ROC/PR discrimination is held out with respect to phase-II training, but the repository does not establish the upstream phase-I training chromosomes, and the nominal FPR thresholds were selected and evaluated on the same chr8–9 set. Brier/ECE also reflect an approximately 1:2 sampled positive:control mixture derived from the control-sampling design rather than genomic prevalence. These are strong internal phase-II split results, not an untouched threshold evaluation, an end-to-end chromosome holdout, or external donor/cell-context validation. Evidence: [training partitions](../04_modeling/results/COPD-S4-R001_training_partitions.tsv), [performance](../04_modeling/results/COPD-S4-R001_model_performance.tsv), [threshold performance](../04_modeling/results/COPD-S4-R001_test_threshold_performance.tsv), [training code](../04_modeling/trednet/TREDNet_v2_seeded.py), model-specific training manifests under 04_modeling/trednet/models_output, and [Section 4 report](../04_modeling/results/section_4_modeling.md).

7. **Allele scoring and candidate definition.** V1 scored 15,303 of 15,389 candidate allele pairs; 86 were unresolved because of allele or coordinate problems. A predicted candidate required both a region score above the internal 5% control-FPR operating point selected from the phase-II test ROC and an absolute REF–ALT delta at or above the 95th percentile within the SNV or indel/complex class of the eligible COPD candidate universe. This yielded 175 enhancer candidates, 199 silencer candidates, 37 in both, and 337 unique records; 297 were SNVs and 40 indels/complex. Twenty-one were tags and 324 were proxies; these statuses are not mutually exclusive because eight records were both. Evidence: [sequence audit](../04_modeling/results/COPD-S4-R002_candidate_sequence_audit.tsv.gz), [allele scores](../04_modeling/results/COPD-S4-R003_candidate_allele_scores.tsv.gz), [delta thresholds](../04_modeling/results/COPD-S4-R004_delta_thresholds.tsv), [variant summary](../04_modeling/results/COPD-S4-R004_variant_summary.tsv), and [final list](../04_modeling/results/COPD-S4-R010_THE_LIST.tsv).

8. **Prioritization.** The 337-record order used regulatory tier, number of positive models, within-class delta percentile, region-threshold ratio, maximum linked r², and stable candidate ID. It is neither a causal posterior nor an experimentally calibrated rank. Evidence: [R010 manifest](../04_modeling/results/COPD-S4-R010_analysis_manifest.json), [prioritization script](../04_modeling/scripts/09_finalize_section4.py), and the rank-basis field in R010.

9. **DeepFootprinting and motif analysis.** TF-MoDISco summarized sequence patterns from 1,000 positives and 100 controls per model. JASPAR similarity and a relative PWM threshold of 0.8 generated allele-compatible motif hypotheses. Motifs were also common outside predicted regions, and V1 correctly did not call them occupancy. Evidence: [R007 summary](../04_modeling/results/COPD-S4-R007_summary.tsv), [R008 variant fractions](../04_modeling/results/COPD-S4-R008_variant_fractions.tsv), [R008 candidate summary](../04_modeling/results/COPD-S4-R008_candidate_tfbs_summary.tsv.gz), and [Section 4 report](../04_modeling/results/section_4_modeling.md).

10. **Target-gene hypotheses.** V1 created 4,779 candidate–gene–method rows covering 1,635 gene labels from nearest TSS, all TSS within 100 kb, Catalog mapped genes propagated through tag–proxy links, selected locus membership, and CDS overlap. These are deliberately broad hypotheses, not effector-gene assignments. Evidence: [target evidence](../04_modeling/results/COPD-S4-R006_candidate_target_evidence.tsv), [target summaries](../04_modeling/results/COPD-S4-R006_candidate_target_summary.tsv), [target-ranking table](../04_modeling/results/COPD-S4-R006_target_gene_ranking.tsv), and [mapping code](../04_modeling/scripts/07_map_candidate_targets.py).

11. **Computational validation.** Exact GTEx v10 Lung significant-pair lookup covered 185 of 337 candidates; 72 had an eGene matching a GWAS/selected-locus gene and 144 matched the broad R006 target set. This was bulk, predominantly non-diseased lung and not colocalization. MPRAbase added regional HepG2 element evidence for one candidate (two nested assay elements), with no proven exact biallelic test. Open Targets provided broad gene context for 297 candidates but was excluded from the exact variant-level union. The integrated exact/regional public union covered 186 candidates and left 151 uncovered, not negative. Evidence: [Section 5 report](../05_computational_validation/results/section_5_computational_validation.md), [GTEx summary](../05_computational_validation/results/COPD-S5-R001_GTEx_v10_Lung_summary.tsv), [MPRAbase results](../05_computational_validation/results/COPD-S5-R002_MPRAbase_summary.tsv), [Open Targets results](../05_computational_validation/results/COPD-S5-R003_gene_catalog_summary.tsv), and [integrated summary](../05_computational_validation/results/COPD-S5-R005_integrated_summary.tsv).

12. **Prospective experimental package.** The 337 candidates collapsed into 153 shared-tag/source-tag components. A diversity-aware 12-candidate panel led to 48 REF/ALT × forward/reverse reporter constructs, candidate plans, and TF-first plans. All remained proposed_not_performed; the shortlist was not a top-causal-12 claim. Evidence: [Section 6 report](../06_experimental_validation/results/section_6_experimental_validation.md), [shortlist](../06_experimental_validation/results/COPD-S6-R001_candidate_shortlist.tsv), [selection audit](../06_experimental_validation/results/COPD-S6-R002_selection_audit.tsv.gz), [constructs](../06_experimental_validation/results/COPD-S6-R003_MPRA_constructs.tsv), and [collaborator README](../06_experimental_validation/results/COPD-S6-R005_collaborator_README.md).

13. **Traceability and validation.** The manuscript maps claims to result files, figures use registered source data, and the final audit passed 81 of 81 checks. This establishes artifact coherence and reproducibility, not biological truth. Evidence: [result-to-manuscript map](../manuscript/result_to_manuscript_map.tsv), [supplementary manifest](../manuscript/supplementary_file_manifest.tsv), [software environment](../manuscript/software_environment.tsv), [version capture](../manuscript/software_version_capture.txt), and [workflow validation](../manuscript/results/COPD-workflow-validation.tsv).

### 1.3 Strongest evidence and conclusions

The strongest V1 evidence is methodological and descriptive:

- The source-to-result chain is unusually auditable and the frozen artifacts are internally coherent.
- Reported, consequence-resolved COPD tags are predominantly noncoding.
- The LD expansion preserves tag, panel, ancestry, matching, and unresolved-record provenance and is appropriate for generating a broad test universe.
- Both sequence models discriminate their phase-II held-out chromosome examples strongly; upstream phase-I exposure and an untouched threshold-evaluation set are not established.
- The prespecified rules reproducibly nominate 337 testable sequence hypotheses.
- Evidence boundaries are mostly stated correctly: LD is not fine-mapping, exact eQTL overlap is not colocalization, motifs are not occupancy, and prospective experiments are not completed validation.

These conclusions are supported by Sections 2–6, R010, the manuscript claim map, and the final workflow audit cited above.

### 1.4 Conclusions that depend on assumptions

The following claims are materially assumption-dependent:

- **Which variants belong in the COPD discovery universe.** The result depends on a primary ontology URI being treated as sufficient to define core COPD despite heterogeneous study designs.
- **Whether model scores measure enhancer/silencer sequence rather than accessibility/source differences.** Most negative controls do not overlap the disease donor's ATAC peaks.
- **Whether one severe-emphysema bulk-lung donor transports to COPD-relevant airway, epithelial, mesenchymal, immune, and alveolar contexts.**
- **Whether forward-orientation scores are stable to reverse complementation, initialization, control sampling, and chromosome split.**
- **Whether H3K27me3-accessible sequence should be called silencer sequence.**
- **Which gene a candidate affects.** Proximity and propagated mapped genes are hypotheses, and the GTEx overlap comparison uses a structurally broad target set.
- **Which allele increases COPD risk.** V1 REF/ALT and unsigned LD do not generally provide disease-risk direction.
- **Whether the 95th-percentile delta rule identifies biologically unusual variants rather than the top tail of an ascertained input set.**

Exact audit values and source files are recorded in [planning_diagnostics.tsv](provenance/planning_diagnostics.tsv).

### 1.5 Major limitations already acknowledged by V1

V1 already acknowledges a single severe-emphysema donor; bulk-tissue and context limitations; lack of statistical fine-mapping; absence of GWAS–QTL colocalization and matched contact maps; sparse external functional recovery; motif-versus-occupancy limits; provisional target genes; no completed wet-lab validation; and the fact that the manuscript is not submission-ready. These appear in the [V1 manuscript](../manuscript/COPD_regulatory_genomics_manuscript.md), [Section 4 report](../04_modeling/results/section_4_modeling.md), [Section 5 report](../05_computational_validation/results/section_5_computational_validation.md), and [Section 6 report](../06_experimental_validation/results/section_6_experimental_validation.md).

The audit therefore does not treat every acknowledged caveat as a V2 work item. The question is which ones can change the central candidate and mechanistic conclusions.

## 2. Gap audit

### 2.1 Ranked gaps by potential impact

| Rank | Actual gap | Classification | Could it alter V1 conclusions? | Why it matters | Primary V1 evidence |
|---:|---|---|---|---|---|
| 1 | Core GWAS phenotype heterogeneity and missing auditable accession-level adjudication | **Central weakness; GO sensitivity analysis** | Yes, directly | GCST90244098 alone supplies 356/660 GWS tags, links to 177/337 candidates, and is the sole linked study accession for 124/337; smoking models, PheCodes, progression, and other designs also enter core | S2-R001, S2-R002, R010, extraction script, decisions COPD-DEC-004 |
| 2 | Negative-control accessibility/source confounding | **Central model-validity weakness; newly identified GO** | Yes, potentially | Only 12.31% of enhancer controls and 15.51% of silencer controls overlap donor ATAC while every positive does; high AUC may partly measure donor accessibility | TREDNet input BEDs, R001 partitions/performance, make_input_training_data.py |
| 3 | Single-donor, bulk severe-emphysema regulatory context | **Central transportability weakness; conditional sensitivity analysis** | Yes for candidate/context claims; no for workflow existence | A candidate may be specific to or absent from airway, AT2, fibroblast, normal lung, or other contexts; severe emphysema is not adjudicated COPD | S3 ENCODE manifest, Section 3/4 reports, R001 |
| 4 | No robust target/mechanism chain, including omitted splicing and polyadenylation evidence | **Central mechanistic weakness; staged GO** | Yes for gene/mechanism claims | V1 target mappings are broad; V1 used significant eQTL overlap only; a now-published COPD-specific lung eQTL/sQTL/apaQTL colocalization resource directly addresses this and shows non-expression mechanisms | R006, S5-R001/R005, Section 5, R010 |
| 5 | No statistical fine-mapping and no study-specific risk-allele orientation | **Important conditional extension** | Yes for selected loci and gain/loss interpretation; not required for V1's narrow framework claim | 316/337 candidates are proxy-only and V1 stores unsigned r²; credible sets/PIPs need compatible full statistics and LD | S2 LD artifacts, S2-R003B, R010 |
| 6 | Sparse and ascertained external functional benchmark | **Useful sensitivity analysis; GO** | Yes for model-sensitivity and scope claims | R009 has only eight literature rsIDs; one exact known functional variant is nominated, one is scorable but below thresholds, and six lack exact identity | S4-R009, R003/R004, R010 |
| 7 | Training-label lobe union and “silencer” semantics | **Major claim-strength limitation plus useful sensitivity; partly unavoidable** | Yes for interpretation of 199/337 silencer-positive and 162 silencer-only calls | All-lobe unions can permit cross-lobe ATAC/histone combinations; accessible H3K27me3 is not a functional silencer assay. Terminology can be corrected, but functional silencer identity remains unresolved without appropriate experiments | input-building code, R001, R004, S3 same-lobe code |
| 8 | Reverse-complement, seed, split, and control-sampling robustness unmeasured | **Useful but potentially consequential model audit; GO/conditional** | Possibly | The CNN is not explicitly strand-invariant; one seed and one split do not characterize uncertainty | TREDNet code, R001 partitions/manifests |
| 9 | Candidate thresholds and calibration are internal, not untouched-test or genomic error control | **Useful sensitivity/interpretation correction** | It can alter candidate-set size, not the existence of signal | The nominal FPR operating point was selected and evaluated on the same phase-II test set; the 95th-percentile delta is conditional on the COPD tag+LD universe; Brier/ECE use an approximately 1:2 sampled mixture from the control design. None is FDR, external calibration, or posterior causality | TREDNet training/evaluation code, R001 threshold performance, R004 thresholds |
| 10 | Upstream phase-I provenance and chromosome exposure are absent locally | **Reproducibility and leakage-interpretation gap** | Possibly for performance interpretation; unlikely to erase candidate existence | The frozen phase-I model has architecture/hash/features but no local training-data or chromosome manifest, so “chromosome holdout” is proven only for phase II | phase-I model, TREDNet code, manifests |
| 11 | Ancestry assignment and reference-LD resolution | **Important for selected loci; partly unavoidable for V1** | Possibly at multi-ancestry loci | V1 assigns all supported superpopulations at study level; 1000 Genomes panels are small for modern fine-mapping | S2 ancestry and LD audit files |
| 12 | TF motif hypotheses lack expression/occupancy confirmation | **Optional extension** | Unlikely to change the central framework claim | Motif similarity is deliberately weak evidence and is common outside model-positive regions | R007/R008 and Section 4 |
| 13 | Absolute paths and incomplete environment locking in attribution steps | **Reproducibility improvement, not scientific gap** | No, unless rerunning fails | The V1 audit is internally complete, but some attribution resources are external absolute paths and environments are not fully lockfile-pinned | software records and R008 manifest |

### 2.2 What is central versus unavoidable

The phenotype universe and model-control construction are the two gaps most capable of changing which candidates appear supported. They should be addressed before expensive context models or mechanistic storytelling.

The single-donor context and broad target mapping are central limits on biological generalization, but public data may not permit a clean complete solution. Their proper V2 outcome may be a bounded sensitivity result or an explicit unresolved status.

Sparse functional examples, limited ancestry references, unavailable in-sample LD, and lack of direct endogenous allele experiments are partly unavoidable. They warrant careful reporting, not manufactured completeness.

Cosmetic changes such as new figure styling, broad annotation accumulation, a larger undifferentiated source list, or another nearest-gene table have low information value and should not be executed.

## 3. Evaluation of proposed modules A–F

### A. External functional benchmark

**Recommendation: GO, as a prespecified frozen-model diagnostic; do not call it a population-wide validation rate.**

| Question | Assessment |
|---|---|
| Exact V1 limitation | S4-R009 contains only eight literature-tested rsIDs and does not use the full tested/null panels from COPD-specific functional studies. |
| Centrality | Important for sensitivity and context dependence, but not fatal to the narrow claim that V1 generated hypotheses. |
| V1 result tested | Frozen allele scores in S4-R003, frozen thresholds/calls in S4-R004, literature recovery in S4-R009, and the original R010 rank. |
| Information gained | Sensitivity to experimentally active alleles, false-negative patterns by mechanism/cell type, score distributions for assayed nulls, and whether failures cluster by accessibility/context or mechanism. This is independent functional evidence when the assay is independent; score concordance is a model diagnostic. |
| Local inputs | R003 scores, R004 thresholds/calls, R009, R010, source records COPD-SRC-044/045/047/048, candidate sequence audit, and frozen model files. |
| Public inputs | Castaldi reports a 606-variant design and 45 unique MPRA allele-effect positives; use the full assayed/null results only if a complete per-variant table is retrieved and verified. Gong reports a 1,120-variant five-locus/three-cell-type design and 25 allele-specific effects; include its assayed nulls only if the supplement exposes a complete results table. HHIP selected functional variants/haplotype, Stuart regional CRISPRi results, Saferali NPNT splice evidence, and Benway's 7,285 published locus candidates/PICS annotations provide typed case-series or association-derived context, not an unbiased benchmark denominator. Details and access checks are in external_resource_audit.tsv. |
| Access | Articles and available supplements are public; no controlled participant data should be needed. A complete per-variant Castaldi tested/null table and the Gong raw/full-null table have not yet been verified, so denominator-level performance claims are gated. |
| Statistical prerequisites | Freeze models, sequences, thresholds, and variant inclusion before looking at benchmark labels. Harmonize exact allele/build/strand. Define denominators from all assayed variants, include reported nulls, stratify by assay and mechanism, and cluster correlated variants/loci. Report confidence intervals and descriptive rank/score distributions; do not tune on these variants. |
| Compute | 1–4 CPU, less than 8 GB RAM, minutes to a few hours; no GPU/BioWulf. Rescoring genuinely new exact alleles with frozen models may use one GPU for less than a few hours but must not retrain. |
| Failure modes | Publication/selection bias; incomplete null tables; correlated variants; different sequence windows and cell types; reporter activity mistaken for endogenous regulation; CRISPRi region effects mistaken for nucleotide effects; build/allele errors; an assay panel drawn from the same GWAS evidence mistaken for independent genetic evidence. |
| Likely manuscript effect | Can support a bounded sensitivity statement, identify context/mechanism failure, and explain why literature-supported examples are or are not recovered. It cannot supply a general sensitivity/specificity claim unless full unbiased denominators are available. |

**FAM13A interpretation.** rs2013701 is not a V1 failure: exact record 4:88963935:G:T is enhancer-positive at frozen rank 210, with enhancer region score 0.9101 and SNV absolute-delta percentile 95.80. It ranks modestly because R010 does not encode external literature confidence and the variant is a proxy supported in one model/tier, not because thresholds were changed against it. The result is encouraging but one positive example is not calibration. Evidence: [R009](../04_modeling/results/COPD-S4-R009_literature_functional_variant_recovery.tsv) and [R010](../04_modeling/results/COPD-S4-R010_THE_LIST.tsv).

**Mechanism-stratified negative control.** rs34712979 is an exact, highly associated NPNT splice-region tag but falls below both model region/delta criteria. A 2025 functional study supports a splice mechanism. That is expected outside an enhancer/H3K27me3 sequence model and must be scored as “mechanism outside model scope,” not as a reason to relax thresholds. Evidence: [GWS tag table](../02_gwas/results/COPD-S2-R002_gws_unique_tag_variants.tsv), [R003](../04_modeling/results/COPD-S4-R003_candidate_allele_scores.tsv.gz), R004 thresholds, and external resource COPD-V2-EXT-015.

**Execution contract.**

1. Register every assayed allele before comparing labels to V1.
2. Use exact allele identity as the primary comparison; report locus/LD containment separately.
3. Retain active, null, discordant, unavailable, and not-assayed states.
4. Report assay class and cell context separately.
5. Never modify V1 thresholds, model weights, or R010 order.

### B. Regulatory-context robustness

**Recommendation: CONDITIONAL GO for new contexts; GO first for a same-donor matched-control validity model.**

| Question | Assessment |
|---|---|
| Exact V1 limitation | V1's positive labels come from one severe-emphysema bulk-lung donor. More importantly, its external allDHS controls usually do not overlap that donor's ATAC peaks, so performance may partly reflect accessibility/source. |
| Centrality | Control confounding is central to model validity. Single-donor context is central to transportability. |
| V1 result tested | R001 test metrics and thresholds; R003 allele scores; R004's 337 calls; model-specific R010 ranks and the 12 S6 selections. |
| Information gained | First, whether the model distinguishes H3K27ac/H3K27me3 state among regions accessible in the same donor rather than donor accessibility. Second, whether frozen candidates/scores transport across normal lung and relevant primary cell contexts. This is sensitivity evidence, not independent COPD validation. |
| Existing local resources | Selected ENCDO520EJG lobe-specific and union peaks; V1 positive/control BEDs; frozen partitions, weights, predictions, and input-building code; ENCODE lung inventory; Roadmap resources considered in S3. |
| Public resources | GSE152779 provides ATAC in normal primary bronchial epithelial, small-airway epithelial, AT2, fibroblast, and 16HBE samples, but no matched H3K27ac/H3K27me3 triplets. ENCODE has adult lung/bronchial/fibroblast records, but the local metadata audit found incomplete or donor-incompatible assay combinations. E096/E128 are older-build normal contexts. |
| Access | Processed resources are public. Controlled access is not necessary for the proposed peak-level sensitivity work. |
| Statistical/technical prerequisites | Same genome build and peak type; matched donor and biological sample where a model is trained; no mixing donors or cell contexts to manufacture a triplet; accessible negative controls matched for chromosome, GC, mappability, length, and ideally peak intensity; prespecified seed ensemble; same-lobe positive sensitivity; thresholds selected on validation/calibration data and evaluated once on untouched test chromosomes; fixed class sampling or prevalence-adjusted calibration when Brier/ECE are compared; fixed reporting metrics and candidate comparisons. Also reproduce V1's exact threshold convention separately for apples-to-apples candidate sensitivity. |
| Compute | Metadata/accessibility overlap: 1–4 CPU, under 16 GB, minutes-hours. One frozen scoring audit: CPU/GPU hours. A two-model, 3–5-seed retraining sensitivity: one GPU/run, roughly 16–64 GB host RAM, hours to days; BioWulf recommended. Context models scale similarly and can run one context/seed per GPU concurrently after inputs are frozen. |
| Failure modes | Incompatible assays/donors; small positive sets; donor/batch effects; normal tissue mislabeled disease validation; cell culture effects; cross-lobe label leakage; changing thresholds after seeing known variants; interpreting model disagreement as biological truth rather than combined biological/technical context. |
| Likely manuscript effect | If matched-control and context results are stable, the sequence-model conclusion strengthens materially. If unstable, the 337 remain V1 hypotheses but context-general claims must contract and V2 should identify context-sensitive subsets without reranking V1. |

**Why the control gate precedes context expansion.** V1 enhancer controls overlap donor ATAC at 12.31% and silencer controls at 15.51%, while positives are accessible by definition. Before spending GPU time on other donors, V2 should rebuild controls from the same donor-accessible universe, excluding positives/promoters/blacklist and matching basic sequence/peak properties. It should train new V2 sensitivity models in isolation, never replace V1 weights, and compare discrimination, scores, and candidate-call concordance with V1. Calibration metrics may be compared only under the same sampled class prevalence or after explicit prevalence adjustment; V2 operating points must be learned without reusing the untouched test chromosomes.

**Context feasibility judgment.**

- **Normal primary lung-cell accessibility:** feasible now as an overlap/transport annotation using GSE152779; not a TREDNet enhancer/silencer retraining dataset by itself.
- **Normal adult bulk lung:** potentially feasible only if a same-donor, technically compatible accessibility plus H3K27ac/H3K27me3 set passes metadata and QC.
- **Fibroblast:** local ENCODE inventory contains relevant assay types but donor matching/file compliance is unresolved; conditional.
- **Airway/bronchial epithelium:** accessibility data are feasible; a complete matched histone triplet was not established; conditional.
- **AT2/alveolar:** GSE152779 accessibility is feasible; matched histone model inputs were not established; conditional.
- **Other cell types:** no-go unless genetic/functional results first implicate the lineage and a coherent matched dataset exists.

No model should be created by combining incompatible donors, assays, or cell states.

### C. GWAS phenotype robustness

**Recommendation: GO; highest-value fast analysis and first scientific execution step.**

| Question | Assessment |
|---|---|
| Exact V1 limitation | Core classification uses primary ontology URIs/text but lacks a reproducible accession-level phenotype/design adjudication record. The core contains materially different COPD-related designs. |
| Centrality | Central: this defines the input tags and provenance of the final candidates. |
| V1 result tested | The 660 GWS tags, 15,389 tag/proxy union, all 337 R010 records, 153 redundancy components, and 12 S6 candidates. |
| Information gained | Which V1 findings are supported by strict susceptibility case-control studies versus ML liability, EHR/PheCode, smoking-stratified/joint, progression/severity, chronic bronchitis/emphysema, or other designs. This is phenotype sensitivity evidence, not a new comprehensive GWAS. |
| Existing local resources | Full core study register, association-level rows, tag provenance, ancestry assignments, LD links, candidate scores/calls, R010, and S6 selection audit. No new model scoring or LD calculation is necessary for the first pass. |
| Public inputs | Study methods/phenotype definitions for ambiguous accessions; representative full-stat metadata for GCST90244098, GCST90399694/5, smoking-stratified studies, MVP PheCode, and Pan-UKB J44. |
| Access | Catalog/publication metadata are public. Complete Sakornsakolpat meta-analysis statistics are not verified public and may require author/consortium contact; constituent cohort participant data are controlled. Adjudicating V1 provenance requires neither. |
| Statistical prerequisites | Prespecify a phenotype taxonomy and strict subset before candidate comparison. Keep spirometry/clinical susceptibility separate from ICD/PheCode EHR; keep progression/severity, lung function, emphysema, chronic bronchitis, smoking-stratified/joint, and cross-disease separate. Record cohort/sample overlap. Do not treat repeated UK Biobank-derived accessions as independent replication. |
| Compute | 1–4 CPU, under 16 GB, minutes to low hours; no GPU or BioWulf. |
| Failure modes | Post hoc strict-set definition; conflating smoking-stratified COPD with smoking behavior; treating EHR and spirometry definitions as identical; dropping unresolved studies; recomputing R010 order; counting overlapping cohorts as replication. |
| Likely manuscript effect | Could materially contract which candidates are described as COPD-susceptibility supported. The generic workflow claim would survive; claims about the breadth and independence of genetic support may change. |

**Planning diagnostic.** One accession-level result is reproducible directly from frozen V1 provenance: ML spirogram-liability study GCST90244098 contributes 356/660 tags, 337 exclusively; it links to 177/337 candidates and is the sole linked study accession for 124/337. Qualitative inspection also found broad case-control, EHR/PheCode, smoking-stratified/joint, within-COPD modifier, and cross/secondary designs among the 35 accessions with GWS association rows. The GWS-tag provenance represents 36 accessions because one additional accession contributes a non-GWS row for a tag made GWS by another study. Exploratory bucket totals used during planning are deliberately not treated as reproducible results because no prespecified accession-to-category register yet exists; Module C must create that register before reporting category aggregates. Evidence: S2-R001, S2-R002, R010, and [planning_diagnostics.tsv](provenance/planning_diagnostics.tsv).

**Proposed analysis.**

1. Build an accession register for all 104 core studies with phenotype class, case definition, ascertainment, analysis type, ancestry, cohort, overlap, full-stat availability, and adjudication rationale.
2. Primary strict stratum: clinically/spirometrically defined COPD susceptibility case-control.
3. Secondary stratum: ICD/PheCode/EHR COPD susceptibility, kept separate.
4. Keep ML/quantitative surrogate, severity/progression, lung function, emphysema, chronic bronchitis, smoking-stratified/joint, shared disease, and gene/burden/CNV tests separate.
5. Reuse frozen model scores and existing LD only after filtering at **study × focal tag × ancestry panel** level. Parse `supporting_studies` in S2-R006B for the retained study set, retain only those tag-panel pairs, and then join them to S2-R006C; R006C itself lacks study identity. This prevents retaining proxies from a panel supported only by an excluded heterogeneous study. No LD rerun is needed when the required tag-panel pair was already computed. Label each V1 candidate as retained, phenotype-dependent, multi-class, or unresolved; do not rerun models or R010 ranking.
6. Compare retention among 660 tags, 15,389 records, 337 candidates, 153 components, and 12 S6 selections.
7. Report cohort overlap and do not infer independent replication from accessions.

The investigator must approve the exact boundary between the primary strict and secondary EHR strata before execution.

### D. Statistical fine-mapping and risk-allele interpretation

**Recommendation: CONDITIONAL GO for selected phenotype-compatible loci; NO-GO for fine-mapping the heterogeneous V1 union. Study-specific effect-allele harmonization is GO once a study is chosen.**

| Question | Assessment |
|---|---|
| Exact V1 limitation | V1 used genome-wide-significant tags plus thresholded unsigned LD. It has no credible sets/PIPs and usually no disease-risk direction for a proxy. |
| Was V1 implementation appropriate? | Yes for broad, ancestry-aware candidate expansion. No implementation defect was found. It is explicitly not fine-mapping and should not be reinterpreted as such. |
| Centrality | Important for causal narrowing and directional mechanism claims, but V1 did not require it for the narrow framework claim. |
| V1 result tested | Whether V1 tags/proxies and the 337 candidates fall in study-specific 95% credible sets; whether the disease-increasing allele aligns to predicted enhancer/H3K27me3-model gain/loss. |
| Information gained | Association evidence narrowed within a specific phenotype/ancestry/study; PIPs/credible sets where valid; study-specific directional consistency. It is not functional validation. |
| Existing local resources | V1 tags, GRCh38 coordinates/alleles, ancestry panels, LD links, R003 scores, R004 calls, R010, and population annotations. The local 1000 Genomes data are not automatically suitable fine-mapping LD. |
| Public resources | Best current pilot: Pan-UKB European ICD-10 J44 GCST90691934 with complete summary data and matching public cohort LD. GBMI, MVP PheCode, ML liability, East Asian, and WGS datasets are useful comparisons but have phenotype/LD constraints. The strict Sakornsakolpat COPD anchor lacks a verified unrestricted complete summary-stat file. |
| Access | Pan-UKB aggregate statistics and LD are public; participant data are controlled. Fine-mapping GCST007692 as published would require complete meta-analysis summary statistics from the authors/consortium plus defensible matched LD. COPDGene or UK Biobank controlled access alone could support a new cohort-specific analysis, not reconstruct the multi-cohort Sakornsakolpat meta-analysis. |
| Statistical prerequisites | One phenotype and ancestry per analysis; complete locus statistics; effect allele/beta/SE or reliable z; per-variant N/AF where required; matching build and allele coding; sufficiently matched signed LD; locus QC; multi-signal method; window/prior sensitivity; MHC/complex-locus handling; convergence and summary-LD diagnostics. |
| Compute | Header/locus lookup: 1–4 CPU, under 16 GB, minutes. Pan-UKB LD extraction: Hail/Spark/cloud or BioWulf, about 4–16 CPU and 32–64 GB, hours. Each typical locus: 1–4 CPU and 8–32 GB, minutes-hours; dense loci more. No GPU. Loci can run concurrently after inputs/QC freeze. |
| Failure modes | Fine-mapping significant hits only; mismatching GRCh38 stats with GRCh37 LD; using small 1000 Genomes LD for a very large cohort without sensitivity checks; mixing ancestries; ignoring cohort overlap; ambiguous palindromes/indels; incomplete variants; MHC complexity; treating PIP as generic COPD causality; transferring tag risk direction through unsigned r². |
| Likely manuscript effect | Valid credible sets can strengthen selected loci and properly limit causal language. Failed prerequisites should produce an explicit unresolved result, not a forced analysis. |

**Feasible pilot.** Pan-UKB GCST90691934 is an EHR ICD-10 J44 phenotype, not spirometry-adjudicated COPD. It has public complete GRCh37 statistics and matching European in-sample dosage LD in a cloud resource. Use original-build statistics with original-build LD, fine-map prespecified autosomal GWS loci, exclude or separately flag the MHC, then allele-validated map credible variants to V1 GRCh38. Label every PIP as Pan-UKB ICD-J44-specific. Do not download the approximately 43.3-TB whole LD resource; extract locus slices.

**Risk-allele direction.** Of 337 final candidates, 316 are proxy-only and V1 uses unsigned r². For a chosen full-stat study:

1. align effect allele, other allele, build, strand, indels, and REF/ALT;
2. define which allele increases the study's COPD phenotype from beta/OR;
3. orient V1 ALT-minus-REF prediction to disease-increasing minus disease-decreasing allele;
4. preserve ambiguous, conflicting, and unavailable states;
5. describe “association-direction-consistent predicted gain/loss,” not causal gain/loss.

Do not derive risk direction from REF/ALT, allele frequency, ancestral allele, or unsigned LD. For a nonlead proxy, direction requires its own effect estimate or appropriately phased signed correlation/conditional result, not r².

### E. Variant-to-target-gene evidence

**Recommendation: staged GO for published COPD-specific molecular-QTL crosswalk and lung rE2G; CONDITIONAL GO for new colocalization; optional conditional GO for lung promoter contacts.**

| Question | Assessment |
|---|---|
| Exact V1 limitation | R006 is proximity/propagation/locus based; S5 uses significant eQTL overlap and broad gene concordance, not GWAS–QTL colocalization or matched contact maps. |
| Centrality | Central to a disease-associated variant → regulatory effect → target gene → COPD biology chain, but not to the narrow nomination-framework claim. |
| V1 result tested | Candidate/locus target hypotheses in R006, GTEx evidence in S5-R001, integrated evidence in S5-R005, R010 target columns, and the S6 candidate/TF plans. |
| Information gained | Published COPD-specific QTL-colocalized genes/mechanisms; predictive lung enhancer-gene links; orthogonal physical contacts; locus-level convergence/contradiction. Evidence remains typed and scoped. |
| Existing local resources | 4,779 target-evidence rows, 337 candidate summaries, GTEx v10 significant pairs, Open Targets context, selected loci, and candidate coordinates. |
| Public resources | Saferali et al. COPD/LTRC/COPDGene Moloc results/browser; ENCODE-rE2G lung ENCSR528UQX; Jung lung PCHi-C GSE86189; eQTL Catalogue release 7 GTEx v8 Lung full cis-eQTL statistics; GTEx v11 SuSiE sets; Natri GSE227136 cell-type lung eQTL for selected unresolved loci. |
| Access | Published Moloc tables/browser, rE2G, PCHi-C, GTEx aggregate, eQTL Catalogue, and GSE227136 processed results are public. LTRC/COPDGene participant data and Natri genotypes are controlled. |
| Statistical prerequisites | For a new coloc: complete locus-level GWAS and QTL statistics, common build and alleles, adequate ancestry/sample information, locus coverage, and multi-signal/conditional treatment. For contacts: audited liftOver and fragment-level scope. For rE2G: retain predictive status. |
| Compute | Published crosswalk/rE2G: 1 CPU, under 4–8 GB, minutes-hours. PCHi-C: 1–4 CPU, 4–8 GB, hours. Targeted coloc: 1–4 CPU, 4–16 GB per locus, minutes-hours. Cell-type QTL: 4–16 CPU, 32–64 GB, hours-days; BioWulf advisable. No GPU. |
| Failure modes | Relabeling exact eQTL overlap as colocalization; treating a colocalized locus/window as an exact causal candidate; nearest-gene circularity; contact as proof of regulation; rE2G as experimental validation; normal/ILD lung as COPD validation; mixing blood and lung; double-counting correlated ENCODE annotations; ignoring multiple QTL/GWAS signals. |
| Likely manuscript effect | Can materially improve target/mechanism claims for a subset and explicitly leave others unresolved. It should not create a universal target assignment. |

**Immediate high-value actions.**

1. **Published COPD Moloc crosswalk — GO.** Saferali et al. used the Sakornsakolpat COPD GWAS with LTRC bulk-lung and COPDGene blood eQTL/sQTL and additional transcript-processing QTL evidence. Crosswalk the published window, QTL class, gene, and posterior summaries to V1 loci and exact candidates, preserving window/locus scope. The published summary reports 38 windows corresponding to 33 loci, while a results-level count has been reported as 32 in one representation; V2 must retain and resolve that discrepancy rather than silently select one count. Underlying participant data are controlled, but the published results/browser are public. This is import/crosswalk of a published colocalization, not a new V2 colocalization claim.

2. **ENCODE lung rE2G — GO.** Intersect frozen candidates/elements with thresholded ENCSR528UQX links. Label these predictive enhancer–gene links from a non-COPD bulk-lung donor. Compare with R006 and QTL targets; do not count agreement as multiple independent validations if annotations share inputs.

3. **Lung promoter-capture Hi-C — CONDITIONAL GO.** GSE86189 provides public hg19 processed lung promoter-other/promoter-promoter contacts. Use only after base-level liftOver validation and restriction-fragment mapping. Label physical contact, not regulatory causation.

4. **New de novo colocalization — CONDITIONAL GO.** eQTL Catalogue QTD000271 provides complete GTEx v8 Lung cis-eQTL statistics on GRCh38 and SuSiE credible sets. Run only at phenotype-compatible loci after Module C and D data gates. Use a multi-signal-capable approach where necessary.

5. **Cell-type lung QTL — second-line CONDITIONAL GO.** GSE227136 is normal/ILD, not COPD, and full processed data are roughly 100 GB. Query only preselected unresolved loci/lineages after bulk-lung evidence; do not pull all files for completeness.

**No-go target work.** Do not add more nearest-gene/Open Targets/significant-pair overlap and call it stronger target validation. Do not hunt generic HiChIP datasets without a locus-specific question. Do not re-map controlled participant-level QTL data in immediate V2 unless public/published results prove insufficient and access/value are separately approved.

### F. Final robustness synthesis

**Recommendation: GO only after evidence-producing modules; NO-GO as a stand-alone analysis.**

| Question | Assessment |
|---|---|
| Exact V1 limitation | V1 has separate files but no V2 architecture that jointly reports phenotype, model, fine-map, functional, target, and experimental evidence while preserving their distinctions. |
| Centrality | Necessary for honest reporting, but it produces no new biological evidence. |
| V1 result tested | Every R010 record and locus/component, plus manuscript claims and S6 plans. |
| Information gained | A structured view of what is robust, sensitive, strengthened, contradicted, uncovered, or unresolved. This is synthesis only. |
| Inputs | Frozen R010, V2 outputs that pass gates, evidence provenance, and scope/independence metadata. |
| Compute | 1 CPU, under 8 GB, minutes-hours; no BioWulf. |
| Failure modes | Omnibus “validated” flag, majority vote, composite pseudo-posterior, double-counting correlated annotations, missing-as-negative, locus evidence assigned to a nucleotide, or a V2 order presented as R010. |
| Likely manuscript effect | Supports a future amendment with calibrated claims. It must follow, not substitute for, actual robustness analyses. |

The synthesis must distinguish:

- robust to phenotype definition;
- sensitive to phenotype definition;
- robust/sensitive to control construction;
- robust/sensitive to regulatory context;
- included in a valid study-specific credible set/PIP;
- directionally harmonized or unresolved;
- supported by external functional evidence;
- supported by exact QTL association;
- supported by formal GWAS–QTL colocalization;
- linked by rE2G prediction;
- linked by physical contact;
- motif-compatible versus occupied;
- reporter-active;
- region-perturbation supported;
- endogenous allele-perturbation supported;
- disease-phenotype supported;
- negative, contradictory, unavailable, not reported, not analyzed, or not covered.

## 4. Newly identified gaps

The items below were not fully captured by proposed modules A–F.

| New issue | Class | Why it could matter | Recommended response |
|---|---|---|---|
| V1 controls mostly lack donor ATAC overlap | **A: true scientific gap** | Can inflate chromosome-held-out discrimination and alter candidate calls | Definite V2 same-donor accessible, matched-control sensitivity models before context proliferation |
| No reverse-complement invariance audit | **A/B: cheap potentially consequential gap** | Forward genomic orientation should not arbitrarily change regulatory predictions | Score forward and reverse-complement with frozen weights; report score/call concordance before retraining |
| One training seed and one chromosome split | **B: useful extension** | Determinism is not prediction uncertainty | If control/RC gates pass, use 3–5 seeds; cross-chromosome folds only if candidate instability remains material |
| All-lobe unions before model-positive overlap | **B: useful sensitivity** | ATAC in one lobe can pair with histone signal in another | Rebuild same-lobe-union positives in V2 and compare; do not edit V1 |
| “Silencer” means accessible H3K27me3, not functional silencing | **A/C: major claim-strength gap, largely unavoidable experimentally** | 199/337 are silencer-positive and 162/337 are silencer-only, so the label materially affects interpretation | Preserve V1 IDs; correct V2 terminology to H3K27me3-associated accessible/repressed-sequence prediction; keep functional silencing unresolved unless an appropriate assay establishes it |
| Splice/apa mechanisms outside model and validation scope | **A/B: material mechanistic gap** | NPNT rs34712979 is a strong exact example; published COPD QTL work shows sQTL/apaQTL relevance | Mechanism-stratified benchmark and published Moloc crosswalk; no model threshold relaxation |
| Test-derived operating point and candidate delta cutoff lack untouched calibration/error control | **B: useful sensitivity** | 337 depends on an internally selected FPR operating point and an ascertained top tail, not controlled causal discoveries | In V2 select thresholds on validation/calibration data, evaluate once on untouched test chromosomes, hold class sampling fixed or prevalence-adjust Brier/ECE, and compare V1 convention separately; use external null/FDR calibration only under a defensible prespecified design |
| Phase-I training data/chromosome provenance unavailable locally | **B/C: performance-interpretation and reproducibility gap** | Phase-II chromosomes are held out, but end-to-end representation leakage cannot be excluded from the local record | Recover the upstream manifest/source if possible; otherwise constrain every holdout claim to phase II and retain unavailable status |
| No cohort-overlap map | **A/B: material evidence-independence gap** | Multiple accessions, especially UK Biobank-derived, can look like replication | Record cohorts and overlap in phenotype register; never count overlapping accessions as independent support |
| Broad R006 makes target “concordance” easy | **A/B: mechanistic calibration gap** | 144/185 overlap with a large hypothesis set is not a validation rate | Compare against target-set size/baseline and prioritize QTL/contact convergence rather than claiming rate |
| Evidence dependence is not explicit | **B: useful synthesis gap** | ENCODE-derived model, rE2G, ABC, and accessibility are correlated | Add source-family/dependence fields and do not double-count |
| Attribution environment has absolute external paths/incomplete lock | **B/C: reproducibility improvement** | Can impede exact rerun but does not change current result meaning | Use a V2 lock/container and relative, checksummed manifests; do not repair V1 in place |
| Genome-wide TF occupancy fishing | **D: low-value expansion** | Motif families and context make untargeted ChIP annotation noisy | No-go; use expression gating and candidate-specific occupancy only if experiments proceed |

The most important newly identified gap is control construction, not another annotation layer. It directly challenges what the high held-out AUCs measure. The most important newly available scientific evidence is the COPD-specific lung/blood molecular-QTL colocalization study, which exposes splicing and transcript-processing mechanisms absent from V1's GTEx significant-pair validation.



## 5. Data-feasibility matrix

The matrix records feasibility as of 2026-10-05. “Summary statistics” means complete locus/genome-level association statistics suitable for the named method, not Catalog significant hits. Exact versions, URLs, file sizes, access status, and limitations are recorded in [external_resource_audit.tsv](provenance/external_resource_audit.tsv).

| Analysis | Required dataset | Phenotype | Ancestry | Genome build | Local? | Public? | Controlled access? | Summary statistics available? | Approximate size | Proposed method | Feasibility | Major limitation |
|---|---|---|---|---|---:|---:|---:|---|---:|---|---|---|
| C1 accession adjudication | V1 core studies, publications, association rows | All V1 COPD-labeled designs | All represented groups | Mixed source; V1 normalized GRCh38 | Yes | Metadata yes | No | Not needed | MB-scale tables | Prespecified manual/design taxonomy plus cohort-overlap register | **High; GO** | Requires investigator agreement on strict phenotype boundary |
| C2 strict-support sensitivity | Frozen tags, S2-R006B tag/panel/study support, S2-R006C tag/panel LD links, R003/R004/R010/S6 | Strict susceptibility primary; EHR secondary | Preserve each study-supported panel | GRCh38 | Yes | Not needed | No | Not needed | Existing V1 tables | Filter study × focal tag × panel in R006B, then join retained tag-panel pairs to R006C; compare retention, never rerank | **High; GO** | R006C lacks study identity; failure to filter R006B first can retain proxies supported only by excluded studies |
| A1 FAM13A benchmark | Castaldi 606-variant design; complete per-variant null/results table if retrievable | COPD FAM13A locus | Association cohorts mixed; assay not ancestry-specific | Verify supplement; harmonize to GRCh38 | No complete panel locally | Article/supplement public | No | Functional table, not GWAS stats | Small | Case-series recovery from verified positives; denominator-level benchmark only if full table passes QC | **High as case series; full benchmark conditional** | Selected locus, reporter context, correlated variants; complete null table not yet verified |
| A2 multi-locus MPRA benchmark | Gong 1,120-variant design; complete per-variant null/results table if retrievable | Five COPD GWAS loci | Not assay-ancestry-specific | Must verify supplement | Selected examples only | Article/supplement link public | No | Functional table; raw/full-null table not verified | Small–moderate | Exact allele/cell-type case series; denominator benchmark only after table gate; cluster by locus | **Conditional-high; GO after table check** | Very recent paper; full null/raw availability not yet verified |
| A3 mechanism benchmark | HHIP, Stuart CRISPRi, NPNT splicing evidence | COPD/shared lung-disease loci | Varies | Verify each source | Selected V1 rows | Public papers | Some underlying cohorts controlled | Not required | Small | Typed case-series benchmark | **High; GO** | Heterogeneous assays cannot be collapsed to one label |
| New M1 reverse-complement audit | Frozen weights and all scorable candidate sequences | Not phenotype-specific | Not applicable | GRCh38 sequences | Yes | Not needed | No | Not needed | 15,303 sequence pairs | Score forward and reverse complement with identical frozen models; report full-universe score metrics, call-transition table, and Jaccard; summarize frozen 337 separately | **High; GO** | Reveals orientation sensitivity but not its biological source |
| New M2 same-donor accessible-control model | ENCDO520EJG ATAC/histone peaks and V1 partitions | Severe-emphysema donor | Donor reported European | GRCh38 | Yes | ENCODE public | No | Not applicable | Local peaks/model inputs plus 15,303 scorable pairs | Matched accessible controls, validation-calibrated threshold, untouched test, 3–5 seeds; score full frozen universe and report transitions/Jaccard | **High; GO after design freeze** | New V2 model; GPU; cannot prove generalization |
| New M3 same-lobe label sensitivity | Lobe-specific ATAC and histone peaks | Severe-emphysema donor lobes | Donor reported European | GRCh38 | Yes | ENCODE public | No | Not applicable | Local peaks plus 15,303 scorable pairs | Build same-lobe positives before union; matched controls; score full frozen universe | **High; GO with M2** | Fewer positives; lobe samples are same donor |
| B1 normal lung-cell accessibility | GSE152779 processed ATAC peaks | Normal primary bronchial, small-airway, AT2, fibroblast, 16HBE | Donor metadata in GEO | GRCh38 reported; verify files | No | Yes | No for processed | Not applicable | ~1.7 GB supplements | Frozen-candidate OCR overlap and cell-context coverage | **High; GO as annotation** | Accessibility only; cannot train enhancer/silencer model |
| B2 alternative TREDNet contexts | Same-donor accessibility + H3K27ac/H3K27me3 peaks | Normal adult lung or relevant cell type | Donor-specific | Prefer GRCh38 | Incomplete metadata inventory local | Potentially | No for processed ENCODE | Not applicable | Varies | Identical V2 training/evaluation with context-specific controls | **Conditional** | No complete compatible matched triplet yet established |
| D1 study-specific effect direction | Pan-UKB GCST90691934 full statistics | ICD-10 J44 COPD | European | GRCh37 original | No | Yes | Participant data only | Yes: beta/SE/EAF/P | ~337 MB raw | Allele/build harmonization; disease-increasing vs decreasing score | **High; GO if EHR phenotype accepted** | EHR phenotype, UKB overlap; direction is study-specific |
| D2 strict COPD fine-map | Sakornsakolpat GCST007692 complete published-meta-analysis stats plus defensible matched LD | Spirometry-defined COPD susceptibility | Predominantly European/multi-cohort | Original GRCh37 | No | Catalog hits only | Authors/consortium likely required; constituent cohorts controlled | No verified complete file | Large/unknown | Multi-signal fine-map only at published-study scope; otherwise label a new cohort-specific analysis separately | **Low now; conditional** | COPDGene/UKB access alone cannot reconstruct the meta-analysis or exact LD |
| D3 Pan-UKB fine-map pilot | GCST90691934 stats and Pan-UKB EUR LD slices | ICD-10 J44 | European | GRCh37 stats + GRCh37 LD | No | Aggregate yes | Participant data no need | Yes | GWAS ~337 MB; full LD ~43.3 TB, locus slices only | SuSiE-RSS/SuSiE-inf or FINEMAP with QC/sensitivity | **Moderate-high; conditional GO** | EHR phenotype; cloud/Hail LD extraction; not generic COPD |
| D4 GBMI EUR association/fine-map | GCST90399694 stats and exact meta-analysis LD | EHR COPD | European | GRCh38 | No | Stats yes | No for stats | Yes | ~862 MB | Association comparison; fine-map only if suitable LD found | **Comparison high; fine-map low** | Multi-biobank mixture and no exact combined in-sample LD |
| E1 published COPD Moloc crosswalk | Saferali published table/browser | COPD susceptibility GWAS; LTRC lung; COPDGene blood | GWAS largely European; QTL multi-ancestry details recorded | GRCh38 | No | Results yes | Participant data controlled | Published posterior results yes; public full QTL matrix not verified | Small browser/table | Crosswalk window/QTL class/gene/posterior to V1; no new coloc label | **High; GO** | Bulk lung/blood; window evidence not exact nucleotide; reported count discrepancy |
| E2 lung rE2G | ENCSR528UQX thresholded links | Non-COPD bulk lung | One reported-European donor | GRCh38 | No | Yes | No | Not applicable | 3.5 MB thresholded; 316.6 MB full | Candidate/element overlap to predictive gene links | **High; GO** | Predictive, non-COPD, shares annotation dependencies |
| E3 lung promoter contacts | GSE86189 lung PCHi-C | Non-COPD primary bulk lung | Not a GWAS cohort | hg19 | No | Yes | No | Not applicable | ~101 MB processed lung files | Validated liftOver, restriction-fragment overlap, promoter link | **Moderate; conditional GO** | Fragment resolution; physical contact is not regulation |
| E4 de novo bulk-lung coloc | Compatible COPD GWAS stats + eQTL Catalogue QTD000271 | Chosen COPD phenotype; GTEx v8 normal lung | GWAS-specific; GTEx predominantly European | GRCh38 QTL; harmonize GWAS | No | Yes | Raw GTEx controlled | Yes for gene eQTL: ~3.9 GB full cis stats | ~3.9 GB plus small credible sets | Multi-signal coloc after allele/coverage QC | **Moderate; conditional GO** | Normal bulk lung; version differs from V1; complete GWAS/matched LD gate |
| E5 cell-type QTL sensitivity | GSE227136 selected lineage/locus files | Normal/ILD lung, not COPD | ~67% self-reported European | Appears GRCh38; verify | No | Processed yes | Genotypes controlled | Processed eQTL statistics yes | ~100 GB all; 0.9–10.3 GB per large file | Query only selected unresolved loci/lineages | **Low-moderate; second-line** | ILD/normal context, modest n, large files |
| F synthesis | Frozen R010 plus completed V2 evidence ledgers | All, kept separate | All, kept separate | Record-specific | Future V2 | Not applicable | No | Not applicable | Small | Long-form evidence ledger and categorical robustness views | **High after upstream work** | Adds no new evidence; vulnerable to double-counting |

### Data availability conclusions

- **Already local and sufficient for immediate work:** phenotype provenance audit, strict-support comparison, reverse-complement audit, same-donor control/label redesign, V1 external score comparison, and final evidence architecture.
- **Small public retrieval with high value:** published COPD Moloc results, thresholded lung rE2G, FAM13A/Gong functional tables if complete, and later lung PCHi-C.
- **Public but large:** GTEx/eQTL Catalogue complete cis statistics, Pan-UKB summary statistics, selected cell-type lung QTL files.
- **Controlled/restricted:** participant-level COPDGene, UK Biobank, LTRC, GTEx, and Natri data. Such access is not required for immediate high-value work; it could enable separately scoped cohort-specific analyses but cannot by itself reconstruct the published Sakornsakolpat meta-analysis or its exact LD.
- **Not verified available:** an unrestricted, complete, strictly spirometry-defined COPD summary-statistic dataset with a sufficiently matched LD reference; a complete alternative lung/cell-type accessibility plus H3K27ac plus H3K27me3 triplet from the same donor/sample; full public LTRC/COPDGene QTL matrices.
- **Negative feasibility finding:** significant Catalog hits or local GTEx significant pairs are insufficient for fine-mapping/colocalization.

## 6. Computational plan

No major job should start until the investigator reviews this plan and the relevant stop/go gate is signed in the V2 decision register.

| Work package | CPU / GPU | Memory | Runtime class | BioWulf? | Safe concurrency | Dependencies |
|---|---|---:|---|---|---|---|
| Accession phenotype register and cohort-overlap map | 1–4 CPU; no GPU | <16 GB | Minutes–hours plus manual review | No | Can run with all read-only fast tasks | Investigator phenotype taxonomy |
| Frozen strict-support comparison | 1–4 CPU; no GPU | <16 GB | Minutes–hours | No | Yes | Completed phenotype register |
| Full functional benchmark crosswalk | 1–4 CPU; optional one GPU for new frozen scoring | <8–16 GB | Minutes–hours | No | Yes | Full tested/null tables; allele/build harmonization |
| Reverse-complement frozen scoring | 4–16 CPU or one GPU | 16–32 GB | Hours | Optional | One job/model; enhancer and silencer can run concurrently if GPU capacity permits | Frozen sequences/weights; no retraining |
| Same-donor control and same-lobe input construction | 4–16 CPU | 16–64 GB | Hours | Optional | Enhancer/silencer preprocessing can run concurrently after design freeze | Control-matching specification |
| Matched-control multi-seed TREDNet sensitivity | One GPU per run; 4–16 CPU support | 32–64 GB host RAM; GPU memory per V1 architecture | Hours–days for 6–10 runs | **Yes** | Context/model/seed runs can be array jobs, subject to storage/GPU quotas | RC audit; frozen inputs; matching QC |
| Alternative-context metadata/peak QC | 1–4 CPU | <16 GB | Minutes–hours | No | Yes | Resource retrieval and donor/assay review |
| Alternative-context models | One GPU/run | 32–64 GB | Hours–days | **Yes** | One context/seed per GPU | Compatible same-donor assay triplet; matched-control method validated |
| Published Moloc crosswalk | 1 CPU | <4 GB | Minutes–hours | No | Yes | Public table/browser export |
| Lung rE2G intersection | 1 CPU | <4 GB | Minutes | No | Yes | 3.5-MB thresholded file |
| Lung PCHi-C intersection | 1–4 CPU | 4–8 GB | Hours | No | Yes after liftOver QA | hg19→GRCh38 mapping and fragment QC |
| Pan-UKB effect harmonization | 1–4 CPU | <16 GB | Minutes–hours/locus | No | Yes | Phenotype decision; full-stat retrieval |
| Pan-UKB LD extraction | 4–16 CPU, Hail/Spark/cloud | 32–64 GB | Hours | **BioWulf or cloud** | Locus slices can run in parallel within quotas | Approved compute route; prespecified loci |
| Fine-mapping | 1–4 CPU/locus | 8–32 GB, more for dense loci | Minutes–hours/locus | Useful for arrays | Loci parallel after LD QC | Complete stats, matched LD, harmonization, complexity gate |
| Targeted de novo colocalization | 1–4 CPU/locus | 4–16 GB | Minutes–hours/locus | Optional | Loci/genes can run concurrently after input QC | Module C/D study choice; full QTL/GWAS stats |
| Cell-type lung QTL follow-up | 4–16 CPU | 32–64 GB | Hours–days | **Yes** | Only selected lineage files; avoid broad parallel download | Unresolved selected loci after bulk evidence |
| Final synthesis/validation | 1–4 CPU | <8 GB | Minutes–hours | No | After all accepted modules | Frozen schemas and completed evidence ledgers |

### Resource discipline

- GPU work starts only after the phenotype-independent RC audit and matched-control design are frozen.
- Do not download whole Pan-UKB LD; use locus slices.
- Do not download the roughly 100-GB cell-type QTL collection until a selected-locus/lineage query is justified.
- Do not download the roughly 260-GB GTEx v10 all-association bucket for immediate V2; use the targeted public eQTL Catalogue route or published COPD Moloc first.
- Hash every retrieved file, retain source metadata, record failed retrievals, and distinguish not available from access denied.
- Each parallel array must write to a unique COPD-V2 result/log namespace and must never target V1 directories.

## 7. Recommended execution order

### 1. Accession-level phenotype adjudication and frozen support sensitivity

**Why first:** It can change the interpretation of more than half of the tag universe and 177 of 337 final candidates, costs little, and defines which full-stat study is legitimate for later fine-mapping/colocalization. It also creates the cohort-overlap map needed to avoid false replication.

**Deliverables:** a 104-accession phenotype register; prespecified strict and EHR strata; retained/phenotype-dependent/unresolved flags for the 660 tags, 15,389 records, 337 candidates, 153 components, and 12 S6 candidates. No new rank.

### 2. Frozen-model diagnostic bundle: reverse complement plus external functional benchmark

**Why second:** Both are fast and can reveal whether model sensitivity is limited by orientation, context, or mechanism before expensive retraining. They are largely independent of the phenotype audit and may run concurrently with it.

**Deliverables:** forward-versus-reverse-complement score metrics and complete call-transition/Jaccard results over all 15,303 frozen scorable pairs, with the frozen 337 summarized separately; full tested/null functional performance only where complete tables are verified, otherwise typed case-series recovery; explicit FAM13A and NPNT case studies.

### 3. Fast disease-relevant target evidence: published COPD Moloc crosswalk plus lung rE2G

**Why third:** These provide the highest target-gene information per unit cost. The Moloc work is directly COPD-focused and includes splicing; rE2G gives a current GRCh38 lung predictive map. They can run concurrently with Steps 1–2 but must preserve locus/candidate scope and evidence dependence.

**Deliverables:** typed locus/candidate/gene evidence rows, discrepancies retained, no new colocalization claim for V2 unless exact published scope permits it.

### 4. Same-donor accessible-control and same-lobe V2 model sensitivity

**Why fourth:** This is the most important newly discovered technical weakness and the highest-value deeper analysis. It follows the RC audit so orientation behavior is understood and follows benchmark assembly so comparisons are prespecified. It precedes new contexts because otherwise context results could reproduce the same control confounding.

**Deliverables:** new COPD-V2 models only; validation-calibrated and untouched-test performance; prevalence-aware calibration; complete 15,303-pair score/call transition tables and Jaccard, with the frozen 337 summarized separately; seed variability; same-lobe-label sensitivity. Never replace V1 weights/calls, rerank R010, or rank V2-only calls.

### 5. Study-specific risk-allele harmonization and selected fine-mapping pilot

**Why fifth:** Module C must establish the phenotype scope first. If EHR J44 is acceptable as a secondary sensitivity phenotype, Pan-UKB offers the cleanest public stats/LD pilot. If the investigator requires spirometry-defined COPD only and no full statistics are authorized, stop rather than substitute a weaker phenotype.

**Deliverables:** study-specific harmonized disease direction; locus QC; credible sets/PIPs for passing loci; unresolved/failed loci retained; overlap with V1 tags/proxies/337 reported without reranking.

### 6. Alternative regulatory-context modeling

**Why sixth:** It is expensive and data compatibility is uncertain. Execute only if Step 4 establishes a defensible control design and if Step 2 or Step 4 indicates context sensitivity worth resolving. Begin with accessibility coverage in GSE152779; train only on coherent same-donor assay triplets.

**Deliverables:** context-specific scores/calls and concordance, not a replacement candidate list.

### 7. Deeper target closure at unresolved, high-value loci

**Why seventh:** Use lung PCHi-C, targeted de novo GTEx/eQTL Catalogue colocalization, or selected GSE227136 cell types only where Steps 1, 3, and 5 leave a specific target question. This prevents a large generic annotation exercise.

**Deliverables:** typed predictive/contact/colocalization evidence with candidate/element/fragment/window scope and contradictions retained.

### 8. Final evidence synthesis

**Why last:** It depends on upstream results and has no independent evidentiary value. Produce robustness strata and evidence ledgers, not a universal score or new causal ranking.

Several fast read-only/CPU steps can safely run concurrently after their schemas are frozen: phenotype metadata collection, external benchmark table ingestion, Moloc crosswalk, and rE2G. GPU retraining and formal fine-mapping should not start concurrently with their prerequisite audits because their designs may change.

## 8. Stop/go decision points

1. **Strict-phenotype definition gate.** If the investigator does not approve whether spirometric/clinical studies alone form the primary stratum and EHR/PheCode studies form a secondary stratum, do not compute a “strict COPD” retention rate.

2. **Phenotype sufficiency gate.** If manual adjudication leaves too few compatible strict studies, report study/phenotype-stratified support rather than pooling an unstable subset.

3. **Cohort independence gate.** If accessions reuse UK Biobank or another cohort, do not count them as independent replication; report shared-cohort support.

4. **Functional-table gate.** If a paper exposes only selected hits and not the full tested/null panel, perform a case-series recovery analysis and do not report sensitivity/specificity.

5. **Model-decision specification and orientation gate.** Before viewing any scoring/retraining result, sign a V2 specification defining full-universe score correlation/absolute-difference metrics, call-transition/Jaccard tolerances, AUC/PR changes with uncertainty, threshold-calibration data, seed aggregation, and the numerical results that can pause context modeling. Without it, the first run is diagnostic only and cannot trigger broader retraining. If reverse-complement results cross the registered tolerance over all 15,303 scorable pairs, elevate architecture/orientation investigation; never average scores into V1.

6. **Control-validity gate.** If same-donor accessible-control models cross preregistered discrimination, uncertainty, or full-universe call-transition thresholds, treat V1 model evidence as control-sensitive and contract model-based claims before adding contexts.

7. **Seed/split gate.** If matched-control calls are unstable across 3–5 seeds under the preregistered metric, do not spend resources on a large context grid; first quantify uncertainty and consider ensemble/reporting changes.

8. **Context-data compatibility gate.** If accessibility and histone files do not come from the same donor/biological context with compatible processing, do not train a complete enhancer/silencer context model. Report not available.

9. **Cell-context priority gate.** If the benchmark or target evidence identifies a strong cell-type-specific failure, promote that context; otherwise do not model every lung lineage for completeness.

10. **GWAS summary-statistics gate.** If complete phenotype-compatible locus statistics with effect alleles/effects/SE and adequate coverage are absent, do not fine-map.

11. **LD gate.** If suitably ancestry/cohort-matched signed LD is absent or summary-LD diagnostics fail, do not report credible sets/PIPs. A 1000 Genomes annotation may be descriptive but not a substitute.

12. **Complex-locus gate.** If the MHC or another locus violates model assumptions, fails convergence, or has incompatible variants, mark it unresolved rather than forcing a credible set.

13. **Risk-direction gate.** If effect allele orientation is ambiguous/conflicting or a proxy lacks its own signed relation/effect, report direction unavailable. Never infer it from REF/ALT, frequency, ancestry, or r².

14. **QTL-statistics gate.** If only significant eQTL pairs are available, do not perform or label colocalization. Use “exact significant-pair overlap” only.

15. **Multiple-signal gate.** If GWAS or QTL data show multiple signals, use conditioning or a SuSiE-capable colocalization approach; otherwise leave the locus unresolved.

16. **Published-coloc scope gate.** A published colocalized 2-Mb window supports a locus/gene/QTL mechanism, not an exact V1 nucleotide, unless the candidate is explicitly in the shared credible signal.

17. **Contact scope gate.** A PCHi-C restriction fragment or rE2G element link stays at fragment/element scope. Contact or prediction alone does not prove regulation.

18. **Evidence-disagreement gate.** If QTL, rE2G, contact, and functional results disagree, retain the contradiction; do not apply majority vote.

19. **Large-download gate.** Do not retrieve whole GTEx/Pan-UKB/cell-QTL resources if targeted endpoints or locus slices answer the prespecified question.

20. **Synthesis gate.** Do not create a V2 integrated view until evidence classes, scopes, provenance, and dependence fields validate. No missing value may silently become zero/negative.

## 9. Proposed V2 evidence architecture

### 9.1 Immutable core and linkage keys

The immutable core is [COPD-S4-R010_THE_LIST.tsv](../04_modeling/results/COPD-S4-R010_THE_LIST.tsv), identified by SHA-256 8241fcb2296ecd80c12fe97277daa2d8c8de87c86a3d8c895628be11cd3b49ce. V2 never adds columns to this file.

All V2 evidence should live in long-form tables linked by one or more of:

- candidate_record_id for exact GRCh38 REF/ALT identity;
- V1 focal/tag ID and tag→proxy link ID;
- prespecified locus/window ID with coordinates and source study;
- regulatory element/peak ID;
- contact restriction-fragment ID;
- gene ID with annotation version;
- study/accession and phenotype-stratum ID.

A locus-level observation must not be copied into an exact candidate row without an explicit relationship and evidence_scope field.

### 9.2 Required evidence-ledger fields

At minimum, every row should contain:

- v2_evidence_id;
- candidate_record_id, locus_id, element_id, and gene_id where applicable;
- evidence_class and evidence_subclass;
- evidence_scope: exact_allele, variant, element, fragment, locus/window, gene, pathway, or phenotype;
- phenotype class and exact case definition;
- cohort/study accession and cohort-overlap group;
- ancestry and LD reference;
- tissue, cell type, donor count, disease state, and assay;
- genome build, coordinate-normalization method, REF/ALT, effect allele, other allele, and orientation status;
- source/version/date, URL/accession/DOI, local checksum, and access status;
- method/software/version and prespecified threshold;
- numerical metric and uncertainty;
- result_state;
- directionality;
- independence_group/source_family;
- comparison_to_v1;
- interpretation boundary and notes.

### 9.3 Evidence classes that remain separate

| Evidence class | Permitted interpretation | Prohibited shortcut |
|---|---|---|
| GWAS association | Variant/region associated in named phenotype/study | Generic COPD causality |
| LD | Correlation in named reference/population | Causal proxy or risk direction |
| Statistical fine-mapping | Study/locus-specific PIP or credible-set membership | Universal posterior across heterogeneous COPD |
| Genomic/regulatory overlap | Spatial overlap with named annotation | Functional activity |
| Sequence-model prediction | Model/context-specific score/call | Endogenous regulation or causality |
| Allele-specific model prediction | Model predicts sequence difference | Disease-risk direction unless harmonized |
| Molecular QTL association | Variant associated with molecular trait | Shared causal signal with GWAS |
| GWAS–QTL colocalization | Evidence for shared signal under named method/data | Exact causal nucleotide/gene without scope |
| rE2G/ABC | Predictive enhancer–gene link | Contact or experimental validation |
| Chromatin contact | Physical proximity in named context/resolution | Regulatory effect |
| Motif similarity/PWM change | Compatible sequence motif | TF occupancy |
| Measured occupancy | TF binding in named assay/context | Allele causality without allele-specific test |
| Reporter assay | Episomal/integrated reporter activity | Endogenous regulation |
| CRISPRi/deletion | Region contributes to regulation | Nearby nucleotide is causal |
| Endogenous allele edit | Allele affects measured endogenous outcome | Organismal COPD phenotype by itself |
| Disease phenotype evidence | Effect on COPD-relevant phenotype | Molecular mechanism without chain |

### 9.4 Result states

Use explicit states:

- positive;
- negative;
- contradictory;
- non-replicating;
- not covered;
- not available;
- access restricted;
- not reported;
- not analyzed;
- failed QC;
- unresolved;
- not applicable.

“Uncovered” and “not analyzed” must never be collapsed into negative. Negative results and failures remain registered.

### 9.5 Robustness outputs

V2 should eventually produce separate, versioned views:

1. **Phenotype-support view:** exact supporting studies and strata for every V1 candidate.
2. **Model-robustness view:** V1, reverse-complement, matched-control, seed, same-lobe, and approved context comparisons.
3. **Fine-map/direction view:** one row per study/locus/candidate with study-specific PIP/credible-set and allele orientation.
4. **Functional-benchmark view:** one row per assayed allele/assay/context, including nulls.
5. **Target-evidence view:** QTL, formal colocalization, rE2G, contact, and perturbation kept separate.
6. **Claim-robustness view:** each manuscript-level conclusion classified as robust, sensitive, strengthened, contradicted, or unresolved with exact evidence references.

A separate V2 noncoding mechanistic-priority view may be created only after review. If created, it should use transparent categorical strata, retain all 337 records, and carry the original R010 rank as a frozen reference column. It must not be called “the updated ranking,” “validated candidates,” or a posterior causal probability. Missing data must not lower a candidate by default. Coding and splice mechanisms must be displayed separately rather than forced into enhancer/silencer evidence.

### 9.6 Independence and double-counting

Every row needs an independence_group. Examples:

- ENCODE peak overlap, TREDNet trained on those peaks, and ENCODE-rE2G are correlated annotation/model evidence.
- Multiple UK Biobank GWAS accessions can share participants.
- Exact GTEx eQTL overlap and a GTEx-derived colocalization are related, not two fully independent validations.
- MPRA, follow-up reporter, and CRISPRi in the same study are distinct assays but share variant selection and context.
- Proximity, Catalog mapped gene, and Open Targets proximity-derived links are not independent target evidence.

Synthesis should display convergent classes and dependencies, not sum them as equal votes.

## 10. Final recommendation

### Analyses to definitely perform after investigator approval

1. **Phenotype-defined GWAS robustness:** build the 104-study adjudication/cohort-overlap register and compare frozen V1 support under primary strict susceptibility and secondary EHR strata. This is the highest-value fast analysis.

2. **Frozen-model validity and external sensitivity:** run the reverse-complement audit and assemble/crosswalk the mechanism-stratified COPD functional benchmark without changing V1 thresholds. Use full-panel performance metrics only if complete tested/null tables pass the functional-table gate; otherwise report typed case-series recovery. Include FAM13A rs2013701 as a recovered model candidate with published endogenous allele editing and NPNT rs34712979 as an expected splice-mechanism non-recovery.

3. **Same-donor matched-control/same-lobe model sensitivity:** create isolated V2 models with donor-accessible matched controls, same-lobe positive construction, and 3–5 seeds; compare—not replace—V1 performance and scores/call transitions across all 15,303 frozen scorable pairs, then summarize the 337 subset. This is the highest-value deeper analysis.

4. **Fast target closure:** crosswalk the published COPD-specific Moloc results and add current lung rE2G predictive links, preserving QTL class, tissue, locus/candidate scope, and evidence dependence.

### Analyses to perform only if prerequisites pass

- Pan-UKB ICD-J44 risk-allele harmonization and locus fine-mapping, only if the investigator accepts EHR COPD as a secondary sensitivity phenotype and approves the LD compute route.
- Strict spirometry-defined COPD fine-mapping if compatible complete statistics and matched LD become accessible.
- Alternative donor/cell-context TREDNet models only after a coherent matched assay triplet and matched-control design pass QC.
- De novo GWAS–QTL colocalization only with complete compatible GWAS and full regional QTL statistics and a multi-signal-capable design.
- Lung PCHi-C after liftOver/fragment QC.
- GSE227136 cell-type QTL only for selected unresolved loci and lineages.
- Cross-chromosome folds or larger seed ensembles only if the first matched-control runs show material uncertainty.

### Analyses not to perform

- Do not fine-map the heterogeneous 660-tag V1 union or call its thresholded 1000 Genomes LD links credible sets.
- Do not launch a comprehensive smoking-genetics project; keep smoking-stratified/joint COPD evidence separate only where it affects interpretation.
- Do not combine donors, assays, or cell types to manufacture complete model inputs.
- Do not tune models or thresholds using known functional variants.
- Do not download whole Pan-UKB LD, all GTEx v10 associations, or all cell-type lung QTL files without a targeted need.
- Do not add more nearest-gene, broad Open Targets, motif, or generic HiChIP annotations merely for completeness.
- Do not call significant eQTL overlap colocalization, rE2G/contact a validated target, motif similarity occupancy, reporter activity endogenous regulation, regional CRISPRi an allele edit, or normal/ILD lung independent COPD validation.
- Do not create an omnibus validation score, causal rank, or replacement for R010.
- Do not treat Module F as an evidence-producing analysis; retain it only as final reporting architecture.

### Single most important weakness

The single most important scientific weakness is that the ontology-defined “core COPD” universe is not a homogeneous COPD-susceptibility case-control set and lacks an auditable manual adjudication record. Because one ML spirogram-liability study links to 177 of 337 candidates and exclusively supports 124, phenotype robustness can materially change the interpretation of the central candidate set.

### Highest-value fast analysis

The phenotype accession/provenance audit and frozen strict-support comparison. It is inexpensive, requires no new model or large download, and directly tests whether the downstream regulatory findings depend on heterogeneous phenotype definitions.

### Highest-value deeper analysis

The same-donor accessible-control and same-lobe multi-seed TREDNet sensitivity analysis. It directly tests whether V1's strong internal discrimination and candidate calls reflect regulatory state rather than donor accessibility/source and whether calls are stable under a defensible control design.

### Module disposition

- A: **GO**
- B: **CONDITIONAL GO for new contexts; the same-donor matched-control sensitivity model is a separate definite GO**
- C: **GO**
- D: **CONDITIONAL GO at selected compatible loci; NO-GO at whole-V1 scope**
- E: **GO in staged form; new de novo colocalization remains conditional**
- F: **GO only as downstream synthesis; drop it as a stand-alone analysis**

### Investigator decisions genuinely required

1. Approve the primary phenotype boundary: recommended clinically/spirometrically defined susceptibility case-control as primary, with ICD/PheCode/EHR as a separate secondary stratum.
2. Decide whether a Pan-UKB ICD-10 J44 fine-map is scientifically acceptable as secondary sensitivity evidence or whether fine-mapping should wait for phenotype-compatible complete statistics and defensible matched LD.
3. Decide whether to request complete GCST007692 meta-analysis statistics and suitable LD from the authors/consortium; separately confirm whether COPDGene/UK Biobank controlled access exists only if a distinctly labeled cohort-specific analysis is worth the delay.
4. Approve the compute route and budget for Pan-UKB locus-LD extraction if that pilot proceeds.
5. Approve GPU scope for matched-control modeling: recommended two model types × 3–5 seeds before any broad context grid.
6. Decide whether future prose should retain the V1 term “silencer” only as a frozen identifier while describing new evidence as H3K27me3-associated accessible/repressed-sequence prediction.

No other investigator choice is needed before the fast read-only/CPU analyses, but no analysis should begin until this plan is reviewed.

---

## Appendix A. V1 evidence index

This index ties planning claims to frozen V1 artifacts. It supplements the inline citations above.

| Domain | Frozen V1 files/results used |
|---|---|
| Scope, claim, freeze | ../README.md; ../decisions.tsv; ../activity_log.tsv; ../results_register.tsv; ../sources.tsv; ../manuscript/COPD_regulatory_genomics_manuscript.md |
| Disease evidence distinctions | ../01_background/README.md; ../01_background/results/section_1_background.md |
| GWAS universe/phenotypes | ../02_gwas/README.md; ../02_gwas/scripts/01_extract_copd_gwas.py; ../02_gwas/results/COPD-S2-R001_studies_core.tsv; ../02_gwas/results/COPD-S2-R002_gws_associations.tsv; ../02_gwas/results/COPD-S2-R002_gws_unique_tag_variants.tsv |
| Effect fields | ../02_gwas/results/COPD-S2-R003B_gws_effect_sizes_all.tsv |
| Coding/noncoding/splice | ../02_gwas/results/COPD-S2-R004_variant_consequences.tsv; ../02_gwas/results/COPD-S2-R004_coding_noncoding_summary.tsv; ../02_gwas/results/COPD-S2-R004_gene_coding_profile.tsv |
| Ancestry/LD | ../02_gwas/scripts/03_prepare_ld.py; ../02_gwas/scripts/04_run_ld.py; ../02_gwas/scripts/05_summarize_ld.py; ../02_gwas/results/COPD-S2-R006A_focal_match_audit.tsv; ../02_gwas/results/COPD-S2-R006B_focal_panel_audit.tsv; ../02_gwas/results/COPD-S2-R006C_ld_pairs.tsv.gz; ../02_gwas/results/COPD-S2-R006D_ld_expansion_by_focal.tsv; ../02_gwas/results/COPD-S2-R006F_summary.tsv |
| Regulatory donor/context | ../03_regulatory_landscape/README.md; ../03_regulatory_landscape/data/encode_peaks/manifest.tsv; ../03_regulatory_landscape/data/encode_lung_experiment_summary.tsv; ../03_regulatory_landscape/scripts/01_download_copd_lung_peaks.py; ../03_regulatory_landscape/scripts/03_gene_locus_regulatory_burden.py; ../03_regulatory_landscape/scripts/05_classify_candidate_variants.py; ../03_regulatory_landscape/results/section_3_regulatory_landscape.md |
| Training/evaluation | ../04_modeling/README.md; ../04_modeling/scripts/01_build_training_inputs.py; ../../../models/TREDNET_v2/make_input_training_data.py; ../04_modeling/trednet/TREDNet_v2_seeded.py; ../04_modeling/trednet/TREDNet_v2_seeded_v1.py; ../04_modeling/results/COPD-S4-R001_training_partitions.tsv; ../04_modeling/results/COPD-S4-R001_model_performance.tsv; ../04_modeling/results/COPD-S4-R001_test_threshold_performance.tsv; ../04_modeling/results/COPD-S4-R001_evaluation_manifest.json; ../04_modeling/trednet/models_output/COPD_SevereEmphysema_Lung_Enhancer_DHS_x2/training_manifest.tsv; ../04_modeling/trednet/models_output/COPD_SevereEmphysema_Lung_Silencer_DHS_x2/training_manifest.tsv |
| Allele scoring/prioritization | ../04_modeling/scripts/03_score_candidate_alleles.py; ../04_modeling/results/COPD-S4-R002_candidate_sequence_audit.tsv.gz; ../04_modeling/results/COPD-S4-R002_candidate_sequence_manifest.json; ../04_modeling/results/COPD-S4-R003_candidate_allele_scores.tsv.gz; ../04_modeling/results/COPD-S4-R003_scoring_manifest.json; ../04_modeling/results/COPD-S4-R004_delta_thresholds.tsv; ../04_modeling/results/COPD-S4-R004_predicted_causal_regulatory_variants.tsv; ../04_modeling/results/COPD-S4-R004_variant_summary.tsv; ../04_modeling/results/COPD-S4-R004_analysis_manifest.json; ../04_modeling/results/COPD-S4-R010_THE_LIST.tsv; ../04_modeling/results/COPD-S4-R010_analysis_manifest.json |
| Targets | ../04_modeling/scripts/07_map_candidate_targets.py; ../04_modeling/results/COPD-S4-R006_candidate_target_evidence.tsv; ../04_modeling/results/COPD-S4-R006_candidate_target_summary.tsv; ../04_modeling/results/COPD-S4-R006_target_gene_ranking.tsv; ../04_modeling/results/COPD-S4-R006_analysis_manifest.json |
| DeepFootprinting/TFBS | ../04_modeling/deepexplainer; ../04_modeling/results/COPD-S4-R007_summary.tsv; ../04_modeling/results/COPD-S4-R007_jaspar_motif_matches.tsv; ../04_modeling/results/COPD-S4-R007_analysis_manifest.json; ../04_modeling/results/COPD-S4-R008_allele_specific_motif_scores.tsv.gz; ../04_modeling/results/COPD-S4-R008_candidate_tfbs_summary.tsv.gz; ../04_modeling/results/COPD-S4-R008_variant_fractions.tsv; ../04_modeling/results/COPD-S4-R008_analysis_manifest.json |
| External functional recovery | ../04_modeling/results/COPD-S4-R009_literature_functional_variant_recovery.tsv; ../04_modeling/results/COPD-S4-R009_summary.tsv |
| Computational validation | ../05_computational_validation/README.md; ../05_computational_validation/results/section_5_computational_validation.md; ../05_computational_validation/results/COPD-S5-R001_GTEx_v10_Lung_summary.tsv; ../05_computational_validation/results/COPD-S5-R002_MPRAbase_summary.tsv; ../05_computational_validation/results/COPD-S5-R003_gene_catalog_summary.tsv; ../05_computational_validation/results/COPD-S5-R005_integrated_candidate_validation.tsv; ../05_computational_validation/results/COPD-S5-R005_integrated_summary.tsv |
| Experimental design | ../06_experimental_validation/README.md; ../06_experimental_validation/results/section_6_experimental_validation.md; ../06_experimental_validation/results/COPD-S6-R001_candidate_shortlist.tsv; ../06_experimental_validation/results/COPD-S6-R002_selection_audit.tsv.gz; ../06_experimental_validation/results/COPD-S6-R003_MPRA_constructs.tsv; ../06_experimental_validation/results/COPD-S6-R004_candidate_validation_plan.tsv; ../06_experimental_validation/results/COPD-S6-R004_TF_first_plan.tsv; ../06_experimental_validation/results/COPD-S6-R005_collaborator_README.md; ../06_experimental_validation/results/COPD-S6_validation_checks.tsv |
| Manuscript/figures | ../manuscript/COPD_regulatory_genomics_manuscript.md; ../manuscript/result_to_manuscript_map.tsv; ../manuscript/figure_data/figure_1_source_data.tsv; ../manuscript/figure_data/figure_2_source_data.tsv; ../manuscript/figure_data/figure_3_design_counts.tsv; ../manuscript/figure_data/figure_3_shortlist_source_data.tsv; ../manuscript/scripts/build_figures.py |
| Provenance/validation | ../manuscript/software_environment.tsv; ../manuscript/software_version_capture.txt; ../manuscript/supplementary_file_manifest.tsv; ../manuscript/results/COPD-workflow-validation.tsv; ../05_computational_validation/results/COPD-S5_validation_checks.tsv; ../06_experimental_validation/results/COPD-S6_validation_checks.tsv |

Paths in the table are relative to 07_gap_closure. Exact planning-derived counts and their source files are separately registered in [planning_diagnostics.tsv](provenance/planning_diagnostics.tsv).

## Appendix B. External-source provenance

Every external dataset or paper evaluated for feasibility is recorded row-by-row in [external_resource_audit.tsv](provenance/external_resource_audit.tsv), including:

- source and version/date;
- URL/accession/DOI;
- access status and controlled-access identifiers;
- phenotype/context and ancestry;
- genome build;
- data type;
- whether complete statistics, significant results, credible sets, processed peaks, or selected results are available;
- approximate size and local status;
- proposed use and major limitation.

The checks were lightweight metadata, header, directory, and primary-paper checks. No large external dataset was downloaded. If a future retrieval changes any status, V2 must append a dated provenance row rather than silently replacing this record.

## Appendix C. Terminology and reporting guardrails

- **Fine-mapping** requires a statistical model over complete locus evidence and suitable LD; V1 LD expansion is not fine-mapping.
- **Colocalization** requires locus-level GWAS and QTL evidence; exact significant-pair overlap is not colocalization.
- **Risk allele** comes from a named association analysis after harmonization; REF/ALT, allele frequency, ancestry, and r² do not define it.
- **Target gene** remains a hypothesis unless evidence supports the specific variant/element-to-gene relation.
- **Occupancy** requires a binding/occupancy assay; motif similarity is not occupancy.
- **Endogenous allele perturbation** is distinct from reporter assay and regional CRISPRi/deletion.
- **COPD validation** requires a COPD-relevant disease phenotype; normal lung, severe emphysema, and ILD are contexts, not synonyms for clinically adjudicated COPD.
- **Negative** means the assay/analysis was performed with adequate coverage and did not support the effect. Not available, not reported, not analyzed, failed QC, and not covered remain different states.
- **V2 robustness** means comparison against the frozen V1 artifact. It does not authorize retroactive threshold changes, removal of nulls, or replacement of V1.
