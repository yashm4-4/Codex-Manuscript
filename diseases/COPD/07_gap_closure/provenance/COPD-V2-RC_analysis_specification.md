# COPD V2 frozen-model reverse-complement audit: locked specification

Specification ID: COPD-V2-RC-SPEC-001. Authorized by attachment `45a070c7-f717-4908-94c9-522798b1851e/pasted-text.txt` on 2026-10-05. Written before reverse-complement (RC) inference or inspection of RC scores. Preflight records this specification's SHA-256 and UTC lock time. Scoring/analysis must verify that hash. Do not alter this specification after RC inference starts; document implementation corrections separately.

## Scope and sources

Audit only the exact frozen V1 enhancer and H3K27me3-associated models including shared phase I. Internal `silencer` names mean H3K27me3-associated activity. New files belong exclusively under `diseases/COPD/07_gap_closure/`, with COPD-V2-RC identifiers. V1 remains immutable. No orientation averaging into V1, recalibration, reranking, replacement candidates, retraining, external functional benchmarking, other V2 modules, manuscript changes, or GitHub push is authorized.

Use all 15,303 exact paired inputs in `04_modeling/data/COPD_candidate_variants_{ref,alt}_2001bp.fa`; authoritative forward R003 scores; all 15,389 R002 sequence-audit rows; R004 prioritized candidates and delta thresholds; R010 THE_LIST; S6 R001 shortlist; PHENO R005/R007; and V1 R009 literature recovery. Retain the 86 unscorable records with original reasons, without generating inputs.

Use the exact V1 loader `models/TREDNET_v2/TREDNet_v2_inference.py`: reconstruct/load frozen phase I; load each phase-II architecture archive AND best `*_phase_two_weights.weights.h5`. The architecture archive alone has pre-training weights. Verify frozen R003 manifest hashes of all weights, sequences, scores and loader; additionally hash architecture files. The accessible Python 3.13.7 interpreter may use the original read-only TF2.20.0/Keras3.14.1/NumPy2.5.0 packages in place of inaccessible Python3.13.0 only after the numerical gate passes. Record exact versions/device; suppress bytecode; use inference mode and deterministic operations where supported; keep caches/logs/temp files in V2.

## Transformation

Uppercase each full frozen 2,001-base sequence, reverse it, and complement A<->T, C<->G, N<->N. Biological REF remains REF; ALT remains ALT. Verify for every allele `RC(RC(s))==s` and `encode(RC(s))==encode(s)[::-1,::-1]`, using frozen float32 [A,C,G,T] one-hot encoding (N=four zeros). Hash forward/RC strings.

The allele's original span [1000,1000+k) maps to [1001-k,1001). Verify allele letters and spans independently for REF/ALT. The original anchor maps to index1000, but multibase RC alleles start earlier. Never recenter, regenerate, normalize, pad, crop or swap alleles. Preserve V1's independently cropped right flanks. Indel results diagnose the exact frozen windows including their original differing edge context; this audit cannot isolate architecture sensitivity from that context.

Use frozen SNV versus indel/complex classes. Retain N-containing scorable sequences using V1 encoding, flag N counts and summarize separately when present. Record allele lengths/spans and distance from edges. Missing/duplicate IDs, pair-order disagreement, wrong lengths, unsupported symbols, allele-span mismatch, failed transformation checks or frozen-input hash mismatch stop inference; no silent exclusion. Preserve existing unscorable boundary/ambiguity reasons.

## Numerical gate and frozen decisions

Before any RC inference, independently forward-score a deterministic verification panel: first64 lexicographically sorted SNV IDs, first64 indel/complex IDs, all12 shortlist IDs and both exact literature IDs, deduplicated. Existing scores remain authoritative. Seed20261001; prediction batch64, outer batch256; smaller batches only for documented resource constraints. Per model require max REF/ALT error<=1e-5 against frozen R003 and recomputed delta error<=2e-5. Failure pauses inference for diagnosis; tolerances cannot be relaxed.

Region=max(REF,ALT); delta=ALT-REF from serialized authoritative scores as in R004. Preserve R003's separately serialized delta and report rounding discrepancy. Call=`scored AND NOT encode_blacklist AND region>=T AND abs(delta)>=D`, inclusive inequalities. All15,303 enter score metrics;20 blacklisted records remain excluded from calls, leaving15,283 eligible (13,747 SNVs;1,536 indel/complex). Reproduce every frozen forward call and R010 ID before RC:175 enhancer,199 H3K27me3,37 both,337 union.

| Model | T | SNV D | Indel/complex D |
|---|---:|---:|---:|
| Enhancer | 0.643623 | 0.05706318769999998 | 0.04936093850000001 |
| H3K27me3-associated | 0.58505 | 0.028802613899999996 | 0.026183359500000003 |

Read and assert original cutoffs; never estimate them from RC. Union=OR. RC-only observations are diagnostic, not newly nominated candidates.

## Metrics and margins

For each model and REF/ALT/region: n; Pearson r; Spearman rho with average ranks for ties; signed mean difference; mean/median/p90/p95/p99/max absolute difference; fractions above1e-5,0.02,0.10. Quantiles use linear interpolation. Constant vectors or n<2 yield explicit not-evaluable correlations.

RC delta=ALT score-REF score. Report correlations, signed mean difference, and mean/median/p90/p95/p99/max of abs(delta_RC-delta_F), abs(abs(delta_RC)-abs(delta_F)), and normalized discrepancy q=abs(delta_RC-delta_F)/D. Report raw sign agreement among pairs with both deltas nonzero, complete gain/loss/zero transitions and gain->loss/loss->gain counts. Meaningful sign uses epsilon=max(2e-5,0.1D): gain above epsilon, loss below minus epsilon, near-zero otherwise. Retention among forward-meaningful effects keeps RC near-zero collapse in the denominator. Neither delta is disease-risk direction.

Forward upper tails: eligible abs(delta_F)>=D, plus a separate forward model/class99th-percentile tail. Freeze those forward-only p99 cutoffs in preflight. Report tail correlations/errors/sign retention/near-zero collapse/call retention with the whole forward-selected denominator. Never select tails using RC.

Forward score-borderline: abs(region_F-T)<=0.02. Forward delta-borderline: abs(abs(delta_F)/D-1)<=0.10. Report four proximity combinations, signed continuous margins, crossing rates, and failed gates (region/delta/both). Well-separated positive=eligible, region_F>=T+0.05, abs(delta_F)>=1.5D. This is distance from computational thresholds, not calibrated biological confidence. Summarize score errors by forward-region bins bounded by0,T-0.10,T-0.02,T+0.02,T+0.10,1 and report discrepancies away from thresholds.

For enhancer/H3K27me3/union: forward-positive, RC-positive, both, forward-only, reverse-only, neither, Jaccard=both/(forward+RC-both), forward retention=both/forward, all2x2 transitions. Empty denominators are not evaluable. Report4x4 contexts: neither/enhancer only/H3K27me3 only/both.

## Prespecified descriptive severity

These are conservative transparent audit conventions, not universal published equivalence margins, clinical thresholds, hypothesis tests or data-tuned claims. The symmetry expectation is supported by [Zhou, Shrikumar and Kundaje,2022](https://proceedings.mlr.press/v165/zhou22a.html); that paper does not supply these numerical bands. On the0-1 score scale,0.005 is half a percentage point,0.05 five percentage points, and0.10 tail error is large. Delta error of one original cutoff can directly overturn effect selection. Positive-set metrics avoid dilution by abundant negatives. LD-linked variants are not independent replicates.

Per model over the complete universe: material if ANY material condition holds; negligible only if ALL negligible conditions hold; otherwise modest. Overall severity is the worst domain across models. Always report continuous effects and class-specific failures.

| Domain | Negligible | Scientifically material |
|---|---|---|
| REF/ALT/region scores | Each MAE<=0.005 and p95 error<=0.02 | Any MAE>=0.05 or p95 error>=0.10 |
| Allele deltas | Median q<=0.05, p95 q<=0.25, meaningful direction retention in forward95th-percentile tail>=99% | Median q>=0.50 or p95 q>=1.0 or tail direction retention<90% |
| Each model and union calls | Jaccard>=0.98, forward retention>=99%, no well-separated-positive loss | Jaccard<0.80 or retention<90% or well-separated retention<95% with at least3 such losses |
| Frozen337 support vectors | >=99% retain exact support vector and no well-separated model-call loss | >5% lose union call or >10% lose at least one original model call |

Union well-separated means at least one well-separated model call;337 well-separated loss means any such original supporting model is lost. A modest classification is predominantly borderline only if>=90% of changed model decisions satisfy the forward proximity band; otherwise state modest but not confined to borderlines. Give counts/denominators and the proportion even when overall severity is material. Zero changes need no borderline explanation.

## Subsets and outputs

Preserve all337 R010 ranks and12 S6 ranks. Per record include all forward/RC scores/deltas/regions/abs-deltas, fixed thresholds, calls, signed/absolute disagreements, margins, proximity flags, sign/context changes, loss of any original model call, gain of an additional RC-only model call, loss of union, exact-context stability and union retention. Support switching=loss of an original model plus gain of the other; also count pure enhancer-only<->H3K27me3-only switches.

Strata: all15,303; SNV/indel-complex; eligible/blacklisted; N-containing/free when present;337;184 PHENO R005 primary_retained;124 GCST90244098_only_study;12 shortlist; exact V1 literature. The124 are not all ML-supported variants. Phenotype classes never alter decisions.

Two exact scorable V1 R009 literature IDs: rs2013701=4:88963935:G:T (R010 rank210), rs7671167=4:88962828:C:T (outside337). Preserve all register rows and mark absent exact IDs unevaluable; no LD substitution or new functional benchmark.

Required outputs: sequence/transformation QC;15,389-record reconciliation; raw15,303 RC scores and full comparisons; rank-preserved337 and12 comparisons; phenotype-subset and literature coverage/exact tables; score/delta/tail/sign metrics;2x2/4x4 transitions; threshold-proximity and severity tables; standalone score/delta/discrepancy/transition plots; detailed report; forward-verification/run logs; model/input/software/script manifests/checksums; final QC; updated V2 activity/decision/result registers. All IDs COPD-V2-RC.

## Provenance and stop

Before inference record Git HEAD, tracked V1 content hashes, original freeze hashes, consumed-file hashes/sizes, and stat metadata for the entire local V1 disease tree including ignored outputs. Resolve model symlinks explicitly. Afterwards recheck consumed/tracked hashes and tree metadata; investigate changes. Record exact commands, UTC stage times, software/device, batches, scripts, spec lock and outputs. Keep caches/logs/temp files in V2 and suppress bytecode.

Test synthetic SNV/insertion/deletion/multibase/N cases then every pair. Stop on integrity/identity failure, failed forward gate, nonfinite/out-of-range scores, incomplete15,303 scoring,337/12 mismatch or V1 mutation. A documented environment/implementation correction consistent with this unchanged spec may resume work; do not relax tolerances or remove records.

Report score/delta/candidate stability; threshold proximity versus well-separated changes; model/class differences;184/124/12/literature behavior; whether sensitivity weakens V1 interpretation and warrants future RC augmentation/invariant-architecture investigation; and the resulting priority of planned matched-control/retraining experiments. RC disagreement is not biological strand-specific regulation. Stop for investigator review after this audit.
