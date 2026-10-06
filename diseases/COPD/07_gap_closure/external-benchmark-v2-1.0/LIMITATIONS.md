# Limitations and interpretation boundaries

Stage `external-benchmark-v2-1.0`; generated 2026-10-06T20:14:33.345171+00:00.

The corrected V2-C ensemble shows a modest one-variant descriptive regional gain over each frozen V1 orientation on the shared positives, not broad external functional recovery.
Absolute recovery remains weak: 35/39 exact in-scope positive variants are not recovered by either frozen region threshold. These are region-threshold non-recoveries in a context-mismatched case series, not asserted biological false negatives.
Common-identity union comparison: higher than V1_forward (3/23 versus 2/23); higher than V1_RC (3/23 versus 2/23).

1. This is a one-time, frozen-model, retrospective/descriptive external evaluation. The benchmark had already been examined under V1; rs2013701 was known in advance, not blinded validation.
2. Experimental ascertainment, missing labels, unresolved identities and incomplete QC denominators prevent sensitivity, specificity, FPR, AUROC, AUPRC and enrichment claims. No such metrics or significance tests were calculated.
3. Complete Castaldi GEO counts/design are not complete published positive/null labels. Unidentified reported hits remain unidentified; low counts, missing barcodes, missing results and unreported significance are not biological negatives.
4. Gong nulls are assay-context outcomes, not universal variant negatives. Positive and null contexts can coexist for one variant; repeated contexts, shared variants and linked loci are not independent observations.
5. 21 frozen reporting groups lack exact usable identities, including the in-scope positive rs57658727. They are excluded as unavailable, not scored through proxies, approximate alleles or post-open identity rescue.
6. Regulatory reporter, TF-binding and endogenous-expression support are not interchangeable with bulk-lung enhancer/H3K27me3 labels. Cell lineage, stimulation, insert orientation, insert length and episomal versus endogenous conditions limit transferability.
7. H3K27me3-associated model scores indicate a learned regional association, not direct repression or causal gene regulation. No frozen assay supports a valid H3-specific allelic direction denominator; general expression/reporter directions were not automatically inverted.
8. Strict enhancer direction uses only 11 observations and 8 unique consensus variants. Its 7 concordant observations do not establish general allelic accuracy; experimental contexts are correlated. Exact zero predictions count as non-concordant ties, and no delta cutoff is invented.
9. Region-level recovery uses whether either allele passes an unchanged regional threshold. It is not evidence that ALT−REF deltas reproduce causal effects. No formal V1-to-V2 allelic-direction improvement analysis was prespecified, so improved regional recovery cannot be generalized to allelic improvement.
10. V1 and V2 have different sequence coverage. Paired gains/losses use only shared exact positive identities. Newly V2-scorable identities are coverage expansion, not V1 negatives. V1 region-only gates are distinct from historical combined region/delta candidate calls.
11. Splice-only NPNT and other out-of-model evidence are preserved, not counted as regulatory false negatives. Contextual contacts or regional perturbations are not promoted into exact variant labels.
12. No internal-reference percentile distribution was prespecified or constructed. Descriptive score quantiles are internal to this benchmark and do not measure calibrated variant pathogenicity or general-population rarity.
13. The hard post-opening firewall forbids retraining, retuning, recalibration, threshold changes, seed selection, model replacement, new controls or a return to internal-stage evaluation. Any future redesign needs a new scientific version and cannot reuse this benchmark as untouched validation.
14. No broader COPD GWAS/candidate scoring, 337-candidate reranking, fine-mapping, target-gene analysis, commit or push was part of this request.

Observed common-identity region-level conclusion: The corrected V2-C ensemble shows a modest one-variant descriptive regional gain over each frozen V1 orientation on the shared positives, not broad external functional recovery.
The sole newly recovered paired union case is rs35421223 at EEFSEC. rs141807665 at GSTCD loses the V1-forward enhancer region call, but its H3K27me3-associated region call retains union recovery; this is not a lost union case. rs7962469 at ACVR1B is recovered among newly V2-scorable sequences and is a coverage expansion, not a paired V1 gain. Despite that regional call, its enhancer ALT-minus-REF direction is discordant with both experimental reporter cell contexts, 16HBE and Jurkat. Those are two assay observations of one variant, not two model-context errors; H3K27me3 direction remains unevaluable.
Detailed quantitative results, missingness, per-locus directions, actual-network RC audit and validation receipt are in `EXTERNAL_BENCHMARK_V2_REPORT.md` and the machine-readable `results/` tables.
