# Frozen V2-C external functional benchmark evaluation

Stage: `external-benchmark-v2-1.0`. Report generated: 2026-10-06T20:12:17.548833+00:00.

## Principal findings

The unchanged V2-C models scored 1,710 exact benchmark identities, including 39 exact in-scope positive variants. Enhancer region recovery was 3/39; H3K27me3-associated recovery was 1/39; their union was 4/39. These are descriptive case-series recovery fractions, not sensitivity.

The corrected V2-C ensemble shows higher descriptive external region-level union recovery than both frozen V1 orientations on the same evaluable positives.
On the common-identity union denominator, V2 is higher than V1_forward (3/23 versus 2/23); higher than V1_RC (3/23 versus 2/23).

Strict enhancer-direction concordance was 7/11 observations (4 discordant, 0 exactly tied), and 5/8 unique nonconflicting variants. This is a small, correlated direction-resolved subset. No H3K27me3 direction denominator is defensible from the frozen assays. Any region-recovery improvement does not establish improved REF–ALT effect prediction; no formal V1-versus-V2 direction-improvement analysis was prespecified.

## Prospective design, one-time opening and frozen model contract

| Event | UTC / evidence |
| --- | --- |
| Prospective specification freeze | 2026-10-06T19:55:51.105921+00:00 |
| Single benchmark opening event | 2026-10-06T19:57:29.781057+00:00 |
| Snapshot ingestion completion | 2026-10-06T19:57:37.215877+00:00 |
| Annotation/identity freeze before V2 predictions | 2026-10-06T20:08:58.985271+00:00 |
| Inference started | 2026-10-06T20:09:05.680041+00:00 |
| Inference completed | 2026-10-06T20:09:46.860827+00:00 |
| Descriptive evaluation completed | 2026-10-06T20:10:19.889067+00:00 |

The benchmark was opened once as a scientific evaluation after the prospective specification and independent pre-open review. Subsequent annotation and analysis used the unchanged snapshot, not a new source selection or repeated tuning cycle. No fresh experimental curation, web retrieval, model selection, threshold fitting, training or candidate-universe scoring occurred. The benchmark had already been examined under V1 and is not newly blinded evidence.

| Frozen artifact | SHA-256 |
| --- | --- |
| PROSPECTIVE_SPECIFICATION.md | a8bbedf23a68c36b592473d542a656fbeb855b776b1a583696672a4a9c3940b1 |
| specification/evaluation_specification.json | 32432be4a7698256b8ae3cb6fe4139eac8a1ad6c83c5b8de41b880ce771e11bc |
| provenance/prospective_freeze.json | d4725233251785d90b3a768a115452e125fe5323a58768df5047043f6f020830 |
| provenance/preopen_source_manifest.json | f1186dfe3cb00eb7753a29321313851d8fd3c6a9926095ef9b283f8150a2ac29 |
| inputs/benchmark_snapshot.json.gz | 94b8cd3960458180e25357e3068ec229b1f260ce9b30d12dff704694f8520e90 |
| provenance/annotation_freeze.json | a434653d42bbe9620c288784f7bb8a694b409ebfd16c280e0de46e7f0edb47ee |
| predictions/inference_manifest.json | 0df7a4f0515cc6bd9f2fdb4e19c6ed2b72c21b7e07ab8e88b4572f01877c5af2 |
| results/evaluation_manifest.json | 05a60575ebbbe14981932ac747b40a9ef62d0afbfb2c303a0b4dd9c2e0220fd8 |

The pre-open source manifest verified 217 original COPD-V2-BENCH1.0 payloads against their original ledger. The master benchmark, all prior stages and original model archives were inputs only. Exact dependency hashes and managed paths are preserved in the source manifest.

| Context | Seed | Original selected checkpoint SHA-256 |
| --- | --- | --- |
| enhancer | 104729 | 1cb670f778a2d148d0ac2b71bea5c97f5a421bb64f744af1dc5b373a60ea3114 |
| enhancer | 130363 | 0e42ecc9364f967220ecc6c3bb0007dc25aced5d9c360355f6e48027c39f0795 |
| enhancer | 155921 | 733ca1845061508c820240d04d4b7c30d8ba0fa7c1086419eaec0b5c2bb041d5 |
| h3k27me3 | 104729 | 503683295ebeb48587f42748181bf6c5815f6a8c79b8dfdd4b669c931109be81 |
| h3k27me3 | 130363 | 56a2d331b6ed27e49f103dce228dda3c26a80ff296e91c9c740e52b9c36841c6 |
| h3k27me3 | 155921 | 995dc295d163a703b80bb3e451b8fc3bd7af01fb6ce041099fc5ce3d9cc82560 |

For each REF or ALT allele, each seed independently evaluates forward and nucleotide reverse-complement inputs. The two network float32 probabilities are promoted to float64 and averaged. The three seed means are averaged in order 104729, 130363, 155921. Delta is ensemble ALT minus ensemble REF. No seed, archive, architecture or input control was changed.

| Context | Frozen decimal threshold | Exact float64 encoding | Rule |
| --- | --- | --- | --- |
| enhancer | 0.74848511815071117 | 0x1.7f39710000001p-1 | max(REF, ALT) >= threshold |
| h3k27me3 | 0.76960810025533055 | 0x1.8a0a12aaaaaacp-1 | max(REF, ALT) >= threshold |

No allele-delta cutoff, internal-reference percentile, enrichment test, significance test, AUROC, AUPRC, FPR, specificity or sensitivity was calculated. Frozen region thresholds are not variant-level false-positive guarantees.

## Population and identity accounting

| Population | Count |
| --- | --- |
| original_master_observations | 14025 |
| contextual_excluded_source_rows | 267 |
| all_benchmark_identity_keys | 1731 |
| exact_benchmark_identity_keys | 1710 |
| exact_V2_scorable_identity_keys_all_states_mechanisms | 1710 |
| experimentally_positive_keys_all_mechanisms_including_unresolved | 41 |
| exact_experimentally_positive_keys_all_mechanisms | 40 |
| splice_only_positive_keys_out_of_model_not_false_negatives | 1 |
| positive_keys_with_out_of_model_observations | 1 |
| in_scope_positive_identity_keys_including_unresolved | 40 |
| exact_in_scope_positive_identity_keys | 39 |
| exact_V2_scorable_in_scope_positive_primary | 39 |
| exact_V2_scorable_reporter_positive | 39 |
| exact_V2_scorable_endogenous_positive | 2 |
| V1_forward_RC_V2_common_positive | 23 |
| strict_enhancer_direction_observations | 11 |
| strict_enhancer_direction_unique_variants_nonconflicting | 8 |
| strict_H3K27me3_direction_observations | 0 |

Variant-level fractions deduplicate GRCh38 chromosome:position:REF:ALT identities. Every original assay-context row is retained separately. Reporter and endogenous-editing subsets overlap and must not be summed. A positive context is not erased by a null context. The 21 unresolved reporting groups are not asserted exact sequence identities. No LD proxy, approximate position, allele/strand rescue or post-open re-normalization was used.

Sequence QC uses the frozen 2,001-bp construction: 1,000 upstream bases, the explicit allele, and shared downstream reference sequence cropped at the right edge. The allele starts at index 1000 for both REF and ALT; indels are not recentered. Genome REF, lengths, alphabet, reverse-complement involution, allele spans and independent sequence hashes were checked. Missing identities or incomplete predictions are unavailable, never negatives. All exact scoring inputs and source-row links are in the identity, sequence-QC and annotation tables.

### Context-specific source states

| Experimental state | Frozen scope | Unique source keys | Exact V2-scorable keys |
| --- | --- | --- | --- |
| ambiguous | yes | 11 | 11 |
| null | yes | 1125 | 1105 |
| positive | no | 1 | 1 |
| positive | partial | 3 | 3 |
| positive | yes | 40 | 39 |
| unavailable | yes | 620 | 619 |

These state populations overlap across assays; they are not a disjoint positive/negative partition.

## Region-level recovery and frozen V1 comparison

### Each method on its own available positive identities

| Population | Method | Enhancer | H3K27me3-associated | Union |
| --- | --- | --- | --- | --- |
| all_in_scope_positive | V1_forward | 2/23 | 1/23 | 2/23 |
| all_in_scope_positive | V1_RC | 1/23 | 1/23 | 2/23 |
| all_in_scope_positive | V2_C_symmetric | 3/39 | 1/39 | 4/39 |
| reporter_positive | V1_forward | 2/23 | 1/23 | 2/23 |
| reporter_positive | V1_RC | 1/23 | 1/23 | 2/23 |
| reporter_positive | V2_C_symmetric | 3/39 | 1/39 | 4/39 |
| endogenous_positive | V1_forward | 1/2 | 0/2 | 1/2 |
| endogenous_positive | V1_RC | 1/2 | 0/2 | 1/2 |
| endogenous_positive | V2_C_symmetric | 1/2 | 0/2 | 1/2 |

Different coverage denominators cannot by themselves establish a paired performance gain. V1 forward and RC values are read from frozen R015 region gates and allele scores only. The historical combined candidate calls also imposed allele-delta/eligibility gates and are deliberately not used as region-only recovery. No V1 inference or threshold change was performed.

### Same identities for all three methods

| Population | Method | Enhancer | H3K27me3-associated | Union |
| --- | --- | --- | --- | --- |
| three_method_exact_positive_intersection | V1_forward | 2/23 | 1/23 | 2/23 |
| three_method_exact_positive_intersection | V1_RC | 1/23 | 1/23 | 2/23 |
| three_method_exact_positive_intersection | V2_C_symmetric | 2/23 | 1/23 | 3/23 |

### Paired V1-to-V2 changes on shared positive identities

| Baseline | Context | Shared n | Retained | Newly recovered | Lost | Still missed |
| --- | --- | --- | --- | --- | --- | --- |
| V1_forward | enhancer | 23 | 1 | 1 | 1 | 20 |
| V1_forward | h3k27me3 | 23 | 1 | 0 | 0 | 22 |
| V1_forward | union | 23 | 2 | 1 | 0 | 20 |
| V1_RC | enhancer | 23 | 1 | 1 | 0 | 21 |
| V1_RC | h3k27me3 | 23 | 1 | 0 | 0 | 22 |
| V1_RC | union | 23 | 2 | 1 | 0 | 20 |

Exact variant-level transitions are retained in `results/paired_V1_V2_variant_transitions.tsv`. The reused variants and assays are correlated; no independence-assuming significance test was performed.

### Coverage expansion (16 positive identities absent from V1 forward coverage)

| Variant | rsID | Locus | Enhancer | H3K27me3 | Union |
| --- | --- | --- | --- | --- | --- |
| 12:51954475:A:G | rs7962469 | ACVR1B_12q13 | yes | no | yes |
| 3:128239571:A:G | rs2999077 | EEFSEC | no | no | no |
| 4:105554837:C:A | rs17035753 | GSTCD | no | no | no |
| 4:105662669:G:A | rs7696836 | GSTCD | no | no | no |
| 4:105743054:A:T | rs7684442 | GSTCD | no | no | no |
| 4:88805132:A:C | rs2276936 | FAM13A;FAM13A_4q22 | no | no | no |
| 4:88808923:C:T | rs2167750 | FAM13A;FAM13A_4q22 | no | no | no |
| 4:88816407:C:G | rs7695177 | FAM13A;FAM13A_4q22 | no | no | no |
| 4:88885633:G:A | rs7687539 | FAM13A | no | no | no |
| 4:88985447:C:T | rs78681184 | FAM13A | no | no | no |
| 5:148524678:T:C | rs10055487 | HTR4 | no | no | no |
| 5:148526442:T:C | rs1368390 | HTR4 | no | no | no |
| 5:148534098:C:T | rs9325103 | HTR4 | no | no | no |
| 5:148535101:C:T | rs74624276 | HTR4 | no | no | no |
| 5:157450895:T:C | rs6874324 | ADAM19 | no | no | no |
| 5:157597731:A:C | rs72811373 | ADAM19 | no | no | no |

These are newly scorable V2 sequences, not paired V1 false negatives or paired gains. The separately preserved RC-coverage expansion table makes the analogous comparison with V1 RC.

### All exact in-scope positives plus the unresolved positive

Scores below are display-rounded; machine-readable tables retain full precision. E and H denote enhancer and H3K27me3-associated contexts.

| Variant / source key | rsID | Locus | E max | H max | E ALT−REF | H ALT−REF | E call | H call | Union | V1 coverage |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 12:51954475:A:G | rs7962469 | ACVR1B_12q13 | 0.9719231 | 0.06018495 | 0.001505087 | -0.001851624 | yes | no | yes | absent_unavailable |
| 3:128111036:C:A | rs35421223 | EEFSEC | 0.7890836 | 0.1054182 | -0.002670954 | -0.0003892599 | yes | no | yes | exact_forward_and_RC |
| 3:128239571:A:G | rs2999077 | EEFSEC | 0.307436 | 0.1785504 | 0.01501791 | -0.009943523 | no | no | no | absent_unavailable |
| 3:128242335:T:A | rs2955083 | EEFSEC | 0.6926523 | 0.2691105 | 0.0129443 | -0.006913669 | no | no | no | exact_forward_and_RC |
| 3:128273095:T:C | rs2811416 | EEFSEC | 0.2308705 | 0.2793926 | 0.005975567 | -0.004530047 | no | no | no | exact_forward_and_RC |
| 3:128288623:A:G | rs6439124 | EEFSEC | 0.3870711 | 0.253149 | -0.001105875 | -0.01025858 | no | no | no | exact_forward_and_RC |
| 3:128308394:C:T | rs2811397 | EEFSEC | 0.4429078 | 0.2074025 | 0.01530092 | -0.003211655 | no | no | no | exact_forward_and_RC |
| 4:105554837:C:A | rs17035753 | GSTCD | 0.3484845 | 0.3962075 | 0.002199247 | -0.003733625 | no | no | no | absent_unavailable |
| 4:105610689:C:T | rs17035960 | GSTCD | 0.3362972 | 0.3677099 | 0.003816982 | -0.002726018 | no | no | no | exact_forward_and_RC |
| 4:105662669:G:A | rs7696836 | GSTCD | 0.6097033 | 0.3418588 | -0.01158606 | 0.003230626 | no | no | no | absent_unavailable |
| 4:105743054:A:T | rs7684442 | GSTCD | 0.4633778 | 0.3088018 | 0.009435306 | -0.006667033 | no | no | no | absent_unavailable |
| 4:105759206:G:A | rs72671892 | GSTCD | 0.340795 | 0.3805732 | -0.005316218 | 0.003948311 | no | no | no | exact_forward_and_RC |
| 4:105814505:C:T | rs142285263 | GSTCD | 0.3249486 | 0.3816584 | -0.0015068 | 0.001064718 | no | no | no | exact_forward_and_RC |
| 4:105829978:A:G | rs76817161 | GSTCD | 0.5720588 | 0.3306912 | 0.0006174942 | 0.004924849 | no | no | no | exact_forward_and_RC |
| 4:105894906:T:TA | rs141807665 | GSTCD | 0.5542237 | 0.9733737 | 0.004619131 | -0.0006551643 | no | yes | yes | exact_forward_and_RC |
| 4:144566782:A:G | rs6537296 | HHIP_4q31 | 0.3483778 | 0.3738731 | 0.0004494339 | -8.245309e-07 | no | no | no | exact_forward_and_RC |
| 4:144567182:C:T | rs1542725 | HHIP_4q31 | 0.380536 | 0.3620727 | -0.002355849 | 0.0004333407 | no | no | no | exact_forward_and_RC |
| 4:88805132:A:C | rs2276936 | FAM13A;FAM13A_4q22 | 0.6134071 | 0.1630293 | -0.008024007 | 0.0007054135 | no | no | no | absent_unavailable |
| 4:88808923:C:T | rs2167750 | FAM13A;FAM13A_4q22 | 0.5061566 | 0.3828911 | -0.02677834 | 0.003475954 | no | no | no | absent_unavailable |
| 4:88816407:C:G | rs7695177 | FAM13A;FAM13A_4q22 | 0.5264781 | 0.3372407 | -0.001786391 | -0.004724776 | no | no | no | absent_unavailable |
| 4:88885633:G:A | rs7687539 | FAM13A | 0.4508038 | 0.3621416 | 0.00439088 | -0.0001138747 | no | no | no | absent_unavailable |
| 4:88939696:C:T | rs2464523 | FAM13A | 0.3025733 | 0.3780399 | 0.001811708 | 0.0003545433 | no | no | no | exact_forward_and_RC |
| 4:88951025:G:A | rs7674369 | FAM13A | 0.2490933 | 0.4217928 | -0.003160097 | -0.008286943 | no | no | no | exact_forward_and_RC |
| 4:88954758:C:T | rs1964516 | FAM13A | 0.4671135 | 0.2578722 | -0.007551014 | 0.01078759 | no | no | no | exact_forward_and_RC |
| 4:88962828:C:T | rs7671167 | FAM13A | 0.5111824 | 0.3866351 | 0.001667768 | 0.0001608183 | no | no | no | exact_forward_and_RC |
| 4:88963935:G:T | rs2013701 | FAM13A | 0.7906222 | 0.134931 | 8.08239e-05 | -0.006963827 | yes | no | yes | exact_forward_and_RC |
| 4:88985447:C:T | rs78681184 | FAM13A | 0.4127075 | 0.3504664 | -0.002218346 | -0.000865072 | no | no | no | absent_unavailable |
| 4:89009162:A:G | rs1795739 | FAM13A | 0.5592687 | 0.2851102 | -0.0103275 | 0.007867739 | no | no | no | exact_forward_and_RC |
| 5:148463704:T:G | rs4705259 | HTR4 | 0.3758016 | 0.3760685 | -0.001511678 | 0.0001985133 | no | no | no | exact_forward_and_RC |
| 5:148524678:T:C | rs10055487 | HTR4 | 0.2156968 | 0.5230268 | 0.003828441 | -0.003546601 | no | no | no | absent_unavailable |
| 5:148526442:T:C | rs1368390 | HTR4 | 0.2703528 | 0.4040736 | 0.0003191705 | 0.0004808207 | no | no | no | absent_unavailable |
| 5:148534098:C:T | rs9325103 | HTR4 | 0.2845571 | 0.3957541 | 0.003690278 | -0.002990549 | no | no | no | absent_unavailable |
| 5:148535101:C:T | rs74624276 | HTR4 | 0.2969035 | 0.3861045 | 0.01290263 | -0.005317882 | no | no | no | absent_unavailable |
| 5:157450895:T:C | rs6874324 | ADAM19 | 0.4490611 | 0.3829594 | -0.01319868 | 0.002673879 | no | no | no | absent_unavailable |
| 5:157501000:C:T | rs56168343 | ADAM19 | 0.151168 | 0.4246058 | 0.0009633352 | -0.0001621842 | no | no | no | exact_forward_and_RC |
| 5:157501815:G:A | rs10476063 | ADAM19 | 0.6436365 | 0.3208456 | -0.009951721 | 0.00292624 | no | no | no | exact_forward_and_RC |
| 5:157597731:A:C | rs72811373 | ADAM19 | 0.5184152 | 0.3199345 | -0.01001331 | 0.005480349 | no | no | no | absent_unavailable |
| 6:142424746:T:G | rs7753012 | ADGRG6 | 0.3156146 | 0.3901089 | -0.0001681795 | 8.269151e-05 | no | no | no | exact_forward_and_RC |
| 6:7562999:T:G | rs2076295 | DSP_6p24 | 0.4437026 | 0.489586 | -0.05624785 | 0.01977012 | no | no | no | exact_forward_and_RC |
| unresolved:rs57658727:-/TTTTATCAACC | rs57658727 | GSTCD | NA | NA | NA | NA | NA | NA | NA | absent_unavailable |

rs2013701 is a pre-existing known case, not a blinded validation case. It remains visibly flagged in all relevant tables rather than being removed post hoc.

## Allelic-effect direction: distinct from region recognition

Only exact alleles, unambiguous frozen experimental ALT-minus-REF activity directions, and reasonably in-scope direct readouts enter this analysis. GWAS risk alleles, generic TF-binding changes, ambiguous MPRA ratio conventions and gene-expression anecdotes do not resolve direction. Positive endogenous allele editing provides partial regulatory support, not a direct bulk-lung chromatin label. The H3K27me3-associated score is a learned state association, not demonstrated repression causality; general expression effects are not automatically inverted.

| Unit | Subset | Context | Eligible n | Concordant | Discordant | Exact ties | Strict fraction | Secondary non-tied fraction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| observation | all_strict_observations | enhancer | 11 | 7 | 4 | 0 | 7/11 | 7/11 |
| observation | reporter | enhancer | 9 | 5 | 4 | 0 | 5/9 | 5/9 |
| observation | endogenous_editing | enhancer | 2 | 2 | 0 | 0 | 2/2 | 2/2 |
| unique_variant | nonconflicting_consensus | enhancer | 8 | 5 | 3 | 0 | 5/8 | 5/8 |
| observation | all_strict_observations | h3k27me3 | 0 | 0 | 0 | 0 | 0/0 (unevaluable) | 0/0 (unevaluable) |
| observation | reporter | h3k27me3 | 0 | 0 | 0 | 0 | 0/0 (unevaluable) | 0/0 (unevaluable) |
| observation | endogenous_editing | h3k27me3 | 0 | 0 | 0 | 0 | 0/0 (unevaluable) | 0/0 (unevaluable) |
| unique_variant | nonconflicting_consensus | h3k27me3 | 0 | 0 | 0 | 0 | 0/0 (unevaluable) | 0/0 (unevaluable) |

An exact unrounded model delta of zero is a tie: it remains in the strict denominator but is not concordant. No epsilon or allelic cutoff converts small deltas to ties. An observation-level percentage is not an independent-variant accuracy estimate.

| Assay | rsID | Cell context | Higher-activity allele | Expected ALT−REF sign | REF score | ALT score | Delta | Result |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| CASTALDI2019_REPORTER_rs2013701_16HBE | rs2013701 | 16HBE; pGL4.23 luciferase; no exposure reported | T | 1 | 0.7905414 | 0.7906222 | 8.08239e-05 | concordant |
| CASTALDI2019_REPORTER_rs7671167_16HBE | rs7671167 | 16HBE; pGL4.23 luciferase; no exposure reported | C | -1 | 0.5095146 | 0.5111824 | 0.001667768 | discordant |
| CASTALDI2019_REPORTER_rs1795739_16HBE | rs1795739 | 16HBE; pGL4.23 luciferase; no exposure reported | A | -1 | 0.5592687 | 0.5489412 | -0.0103275 | concordant |
| CASTALDI2019_EDIT_rs2013701_16HBE_expression | rs2013701 | 16HBE; CRISPR homology-directed GG-to-TT editing; submerged clonal culture | T | 1 | 0.7905414 | 0.7906222 | 8.08239e-05 | concordant |
| ZHOU2012_rs6537296_forward | rs6537296 | BEAS-2B | G | 1 | 0.3479284 | 0.3483778 | 0.0004494339 | concordant |
| ZHOU2012_rs1542725_reverse | rs1542725 | BEAS-2B | T | 1 | 0.380536 | 0.3781801 | -0.002355849 | discordant |
| BOUEIZ2019_rs7962469_16HBE | rs7962469 | 16HBE | A | -1 | 0.970418 | 0.9719231 | 0.001505087 | discordant |
| BOUEIZ2019_rs7962469_Jurkat | rs7962469 | Jurkat | A | -1 | 0.970418 | 0.9719231 | 0.001505087 | discordant |
| HAO2020_rs2076295_reporter | rs2076295 | 16HBE14o- | T | -1 | 0.4437026 | 0.3874547 | -0.05624785 | concordant |
| HAO2020_rs2076295_HDR | rs2076295 | 16HBE14o- isogenic HDR clones | T | -1 | 0.4437026 | 0.3874547 | -0.05624785 | concordant |
| LIN2020_rs2276936_forward | rs2276936 | HepG2 | A | -1 | 0.6134071 | 0.6053831 | -0.008024007 | concordant |

### Unique-variant consensus and locus accounting

| Variant | rsID | Locus | Context | Expected sign | Predicted sign | Delta | Result |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 12:51954475:A:G | rs7962469 | ACVR1B_12q13 | enhancer | -1 | 1 | 0.001505087 | discordant |
| 4:144566782:A:G | rs6537296 | HHIP_4q31 | enhancer | 1 | 1 | 0.0004494339 | concordant |
| 4:144567182:C:T | rs1542725 | HHIP_4q31 | enhancer | 1 | -1 | -0.002355849 | discordant |
| 4:88805132:A:C | rs2276936 | FAM13A;FAM13A_4q22 | enhancer | -1 | -1 | -0.008024007 | concordant |
| 4:88962828:C:T | rs7671167 | FAM13A | enhancer | -1 | 1 | 0.001667768 | discordant |
| 4:88963935:G:T | rs2013701 | FAM13A | enhancer | 1 | 1 | 8.08239e-05 | concordant |
| 4:89009162:A:G | rs1795739 | FAM13A | enhancer | -1 | -1 | -0.0103275 | concordant |
| 6:7562999:T:G | rs2076295 | DSP_6p24 | enhancer | -1 | -1 | -0.05624785 | concordant |

| Locus | Context | Observations | Unique variants | Concordant | Discordant | Tied |
| --- | --- | --- | --- | --- | --- | --- |
| ACVR1B_12q13 | enhancer | 2 | 1 | 0 | 2 | 0 |
| DSP_6p24 | enhancer | 2 | 1 | 2 | 0 | 0 |
| FAM13A | enhancer | 4 | 3 | 3 | 1 | 0 |
| FAM13A_4q22 | enhancer | 1 | 1 | 1 | 0 | 0 |
| HHIP_4q31 | enhancer | 2 | 2 | 1 | 1 | 0 |

Conflicting experimental directions, if present, are retained at observation level and excluded only from the unique-variant consensus summary. The frozen annotation map contained 0 such conflicting variants. Unevaluable directions and their reasons are preserved in `results/direction_unevaluable_observations.tsv`.

## Source, assay, mechanism and context stratification

Each cell below counts unique exact scorable positive variants within that stratum. Rows and strata overlap; small subsets are descriptive and do not establish tissue-specific validation.

### V1_presence

| Stratum | Enhancer | H3K27me3-associated | Union |
| --- | --- | --- | --- |
| absent_unavailable | 1/16 | 0/16 | 1/16 |
| exact_forward_and_RC | 2/23 | 1/23 | 3/23 |


### assay_class

| Stratum | Enhancer | H3K27me3-associated | Union |
| --- | --- | --- | --- |
| MPRA_allele_effect | 2/35 | 1/35 | 3/35 |
| TF_binding | 0/1 | 0/1 | 0/1 |
| conventional_reporter | 2/8 | 0/8 | 2/8 |
| endogenous_allele_editing | 1/2 | 0/2 | 1/2 |


### cell_context

| Stratum | Enhancer | H3K27me3-associated | Union |
| --- | --- | --- | --- |
| 16HBE | 1/15 | 1/15 | 2/15 |
| 16HBE14o- | 0/1 | 0/1 | 0/1 |
| 16HBE14o- isogenic HDR clones | 0/1 | 0/1 | 0/1 |
| 16HBE; CRISPR homology-directed GG-to-TT editing; submerged clonal culture | 1/1 | 0/1 | 1/1 |
| 16HBE; pGL4.23 luciferase; no exposure reported | 1/3 | 0/3 | 1/3 |
| BEAS-2B | 0/2 | 0/2 | 0/2 |
| BEAS-2B nuclear extract | 0/1 | 0/1 | 0/1 |
| Beas-2B; source-reported any oligo/promoter context | 1/11 | 0/11 | 1/11 |
| HUVEC | 0/5 | 0/5 | 0/5 |
| HepG2 | 0/1 | 0/1 | 0/1 |
| Jurkat | 1/1 | 0/1 | 1/1 |
| MRC5 | 1/6 | 0/6 | 1/6 |


### evidence_unit

| Stratum | Enhancer | H3K27me3-associated | Union |
| --- | --- | --- | --- |
| endogenous_editing | 1/2 | 0/2 | 1/2 |
| other_in_scope_evidence | 0/1 | 0/1 | 0/1 |
| reporter | 3/39 | 1/39 | 4/39 |


### locus

| Stratum | Enhancer | H3K27me3-associated | Union |
| --- | --- | --- | --- |
| ACVR1B_12q13 | 1/1 | 0/1 | 1/1 |
| ADAM19 | 0/4 | 0/4 | 0/4 |
| ADGRG6 | 0/1 | 0/1 | 0/1 |
| DSP_6p24 | 0/1 | 0/1 | 0/1 |
| EEFSEC | 1/6 | 0/6 | 1/6 |
| FAM13A | 1/11 | 0/11 | 1/11 |
| FAM13A_4q22 | 0/1 | 0/1 | 0/1 |
| GSTCD | 0/8 | 1/8 | 1/8 |
| HHIP_4q31 | 0/2 | 0/2 | 0/2 |
| HTR4 | 0/5 | 0/5 | 0/5 |


### mechanism_in_model_scope

| Stratum | Enhancer | H3K27me3-associated | Union |
| --- | --- | --- | --- |
| partial | 1/3 | 0/3 | 1/3 |
| yes | 3/39 | 1/39 | 4/39 |


### regulatory_mechanism

| Stratum | Enhancer | H3K27me3-associated | Union |
| --- | --- | --- | --- |
| TF_binding_partial_scope | 0/1 | 0/1 | 0/1 |
| endogenous_expression_partial_scope | 1/2 | 0/2 | 1/2 |
| transcriptional_reporter_not_H3_specific | 3/39 | 1/39 | 4/39 |


### study_id

| Stratum | Enhancer | H3K27me3-associated | Union |
| --- | --- | --- | --- |
| BOUEIZ2019_ACVR1B | 1/1 | 0/1 | 1/1 |
| CASTALDI2019 | 1/11 | 0/11 | 1/11 |
| GONG2026 | 1/24 | 1/24 | 2/24 |
| HAO2020_DSP | 0/1 | 0/1 | 0/1 |
| LIN2020_FAM13A | 0/1 | 0/1 | 0/1 |
| ZHOU2012_HHIP | 0/2 | 0/2 | 0/2 |


### variant_type

| Stratum | Enhancer | H3K27me3-associated | Union |
| --- | --- | --- | --- |
| SNV | 3/38 | 0/38 | 3/38 |
| indel_or_multibase | 0/1 | 1/1 | 1/1 |


## Continuous region and delta distributions

Quantiles use the prespecified NumPy float64 linear convention. They summarize this benchmark only; no internal reference distribution or outcome-conditioned percentile was constructed.

| Subset | Context | Score | n | Minimum | Q25 | Median | Q75 | Maximum |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| all_scorable_exact | enhancer | ref_score | 1710 | 0.04403425 | 0.3430664 | 0.4153059 | 0.5228953 | 0.9984607 |
| all_scorable_exact | enhancer | alt_score | 1710 | 0.04419324 | 0.3423962 | 0.4148712 | 0.5246781 | 0.998335 |
| all_scorable_exact | enhancer | max_allele_score | 1710 | 0.04419324 | 0.3444059 | 0.4179041 | 0.5274795 | 0.9984607 |
| all_scorable_exact | enhancer | delta | 1710 | -0.05624785 | -0.003887424 | 7.659445e-05 | 0.00358911 | 0.08425103 |
| all_scorable_exact | enhancer | abs_delta | 1710 | 1.400709e-06 | 0.001480192 | 0.003709572 | 0.007527489 | 0.08425103 |
| all_scorable_exact | h3k27me3 | ref_score | 1710 | 0.02774032 | 0.3132235 | 0.3622567 | 0.3861706 | 0.9820754 |
| all_scorable_exact | h3k27me3 | alt_score | 1710 | 0.02602257 | 0.3140988 | 0.3622059 | 0.3856671 | 0.9809916 |
| all_scorable_exact | h3k27me3 | max_allele_score | 1710 | 0.02774032 | 0.3162606 | 0.3633533 | 0.3866842 | 0.9820754 |
| all_scorable_exact | h3k27me3 | delta | 1710 | -0.03509869 | -0.001839149 | 6.29127e-05 | 0.001896148 | 0.02919235 |
| all_scorable_exact | h3k27me3 | abs_delta | 1710 | 8.245309e-07 | 0.0006409536 | 0.001873545 | 0.004525201 | 0.03509869 |
| in_scope_positive | enhancer | ref_score | 39 | 0.1502046 | 0.3202816 | 0.4276068 | 0.5380413 | 0.970418 |
| in_scope_positive | enhancer | alt_score | 39 | 0.151168 | 0.3194441 | 0.4104891 | 0.5368164 | 0.9719231 |
| in_scope_positive | enhancer | max_allele_score | 39 | 0.151168 | 0.3202816 | 0.4429078 | 0.5403509 | 0.9719231 |
| in_scope_positive | enhancer | delta | 39 | -0.05624785 | -0.004238157 | 8.08239e-05 | 0.00375363 | 0.01530092 |
| in_scope_positive | enhancer | abs_delta | 39 | 8.08239e-05 | 0.001589723 | 0.003816982 | 0.009982516 | 0.05624785 |
| in_scope_positive | h3k27me3 | ref_score | 39 | 0.06018495 | 0.2783175 | 0.3621416 | 0.3862894 | 0.9733737 |
| in_scope_positive | h3k27me3 | alt_score | 39 | 0.05833333 | 0.2799864 | 0.3620727 | 0.3847972 | 0.9727186 |
| in_scope_positive | h3k27me3 | max_allele_score | 39 | 0.06018495 | 0.2822514 | 0.3621416 | 0.3863698 | 0.9733737 |
| in_scope_positive | h3k27me3 | delta | 39 | -0.01025858 | -0.003640113 | -0.0001138747 | 0.001869299 | 0.01977012 |
| in_scope_positive | h3k27me3 | abs_delta | 39 | 8.245309e-07 | 0.0005679925 | 0.003211655 | 0.005399115 | 0.01977012 |
| reporter_positive | enhancer | ref_score | 39 | 0.1502046 | 0.3202816 | 0.4276068 | 0.5380413 | 0.970418 |
| reporter_positive | enhancer | alt_score | 39 | 0.151168 | 0.3194441 | 0.4104891 | 0.5368164 | 0.9719231 |
| reporter_positive | enhancer | max_allele_score | 39 | 0.151168 | 0.3202816 | 0.4429078 | 0.5403509 | 0.9719231 |
| reporter_positive | enhancer | delta | 39 | -0.05624785 | -0.004238157 | 8.08239e-05 | 0.00375363 | 0.01530092 |
| reporter_positive | enhancer | abs_delta | 39 | 8.08239e-05 | 0.001589723 | 0.003816982 | 0.009982516 | 0.05624785 |
| reporter_positive | h3k27me3 | ref_score | 39 | 0.06018495 | 0.2783175 | 0.3621416 | 0.3862894 | 0.9733737 |
| reporter_positive | h3k27me3 | alt_score | 39 | 0.05833333 | 0.2799864 | 0.3620727 | 0.3847972 | 0.9727186 |
| reporter_positive | h3k27me3 | max_allele_score | 39 | 0.06018495 | 0.2822514 | 0.3621416 | 0.3863698 | 0.9733737 |
| reporter_positive | h3k27me3 | delta | 39 | -0.01025858 | -0.003640113 | -0.0001138747 | 0.001869299 | 0.01977012 |
| reporter_positive | h3k27me3 | abs_delta | 39 | 8.245309e-07 | 0.0005679925 | 0.003211655 | 0.005399115 | 0.01977012 |
| endogenous_positive | enhancer | ref_score | 2 | 0.4437026 | 0.5304123 | 0.617122 | 0.7038317 | 0.7905414 |
| endogenous_positive | enhancer | alt_score | 2 | 0.3874547 | 0.4882466 | 0.5890385 | 0.6898304 | 0.7906222 |
| endogenous_positive | enhancer | max_allele_score | 2 | 0.4437026 | 0.5304325 | 0.6171624 | 0.7038923 | 0.7906222 |
| endogenous_positive | enhancer | delta | 2 | -0.05624785 | -0.04216568 | -0.02808351 | -0.01400134 | 8.08239e-05 |
| endogenous_positive | enhancer | abs_delta | 2 | 8.08239e-05 | 0.01412258 | 0.02816433 | 0.04220609 | 0.05624785 |
| endogenous_positive | h3k27me3 | ref_score | 2 | 0.134931 | 0.2186522 | 0.3023734 | 0.3860947 | 0.4698159 |
| endogenous_positive | h3k27me3 | alt_score | 2 | 0.1279672 | 0.2183719 | 0.3087766 | 0.3991813 | 0.489586 |
| endogenous_positive | h3k27me3 | max_allele_score | 2 | 0.134931 | 0.2235948 | 0.3122585 | 0.4009222 | 0.489586 |
| endogenous_positive | h3k27me3 | delta | 2 | -0.006963827 | -0.0002803411 | 0.006403144 | 0.01308663 | 0.01977012 |
| endogenous_positive | h3k27me3 | abs_delta | 2 | 0.006963827 | 0.0101654 | 0.01336697 | 0.01656854 | 0.01977012 |


## Missing, unresolved and mechanistically out-of-model evidence

| Unresolved source key | rsID | In-scope positive? | Frozen state | Reason |
| --- | --- | --- | --- | --- |
| unresolved:rs10634924:-/AA | rs10634924 | no | null | unresolved_exact_identity |
| unresolved:rs10637017:-/AAG | rs10637017 | no | null | unresolved_exact_identity |
| unresolved:rs10699204:-/GCATC | rs10699204 | no | null | unresolved_exact_identity |
| unresolved:rs113955166:-/GAAAA | rs113955166 | no | null | unresolved_exact_identity |
| unresolved:rs142307498:-/A | rs142307498 | no | null | unresolved_exact_identity |
| unresolved:rs199856654:-/ATGTAAAAT | rs199856654 | no | null | unresolved_exact_identity |
| unresolved:rs199910359:-/TGA | rs199910359 | no | null | unresolved_exact_identity |
| unresolved:rs200077507:-/GAT | rs200077507 | no | null | unresolved_exact_identity |
| unresolved:rs200407928:-/AC | rs200407928 | no | null | unresolved_exact_identity |
| unresolved:rs200608396:-/A | rs200608396 | no | null | unresolved_exact_identity |
| unresolved:rs201292172:-/ACAC | rs201292172 | no | null | unresolved_exact_identity |
| unresolved:rs201529721:-/TA | rs201529721 | no | null | unresolved_exact_identity |
| unresolved:rs202157195:-/TGTG | rs202157195 | no | null | unresolved_exact_identity |
| unresolved:rs34671427:-/GA | rs34671427 | no | null | unresolved_exact_identity |
| unresolved:rs35053582:-/AG | rs35053582 | no | null | unresolved_exact_identity |
| unresolved:rs3835235:-/T | rs3835235 | no | null | unresolved_exact_identity |
| unresolved:rs55840827:-/TCACGCCT | rs55840827 | no | null | unresolved_exact_identity |
| unresolved:rs57658727:-/TTTTATCAACC | rs57658727 | yes | conflicting | unresolved_exact_identity |
| unresolved:rs71153122:-/AAAC | rs71153122 | no | null | unresolved_exact_identity |
| unresolved:rs72589169:-/AA | rs72589169 | no | null | unresolved_exact_identity |
| unresolved:rs76419734:C/T | rs76419734 | no | unevaluable | unresolved_exact_identity |

The unresolved rs57658727 positive remains explicit and does not enter the exact positive denominator. Its invalid/ambiguous deposited contrast is not rescued by another orientation or an approximate genomic identity. All remaining unresolved keys are likewise unavailable, not negatives.

### Splice-only and other out-of-model positive variants

| Variant | rsID | Locus | Mechanism | V2 sequence scorable? | Interpretation |
| --- | --- | --- | --- | --- | --- |
| 4:105897896:G:A | rs34712979 | GSTCD;NPNT_4q24 | enhancer_like_reporter;splice_only_out_of_model | yes | Sequence scoring retained for complete identity accounting; outside regulatory-positive recovery and direction denominators, not a false negative |

NPNT splice-only evidence is preserved as out of model. The 267 contextual/excluded observations remain separate, including regional perturbations and contacts that do not establish an exact variant label.

### Frozen source denominator limitations

| Study | Complete denominator pass | Preserved limitation |
| --- | --- | --- |
| CASTALDI2019 | False | GEO counts and design enumerate all 606 variants and their allele pairs; they do not contain P/FDR/significance or analytical QC columns. The inaccessible supplement prevents reconstruction of published labels and 45-hit identities. Raw count reanalysis would require unverified filtering/normalization/testing choices; no such analysis was performed. |
| GONG2026 | NA | IncompleteQC/resultdenominator; preserve recoverable nulls and unavailable states; no fullpanel sensitivity/specificity/AUROC claim authorized. |
| BOUEIZ2019_ACVR1B | NA | No complete screening denominator beyond deliberately selected mechanistic examples. HHIP enumerates all four isolated Figure3B tests (two selected SNPs in one locus, each in two insert orientations); per-orientation descriptive counts are supported, but these correlated case-series observations do not establish population-like discrimination. Other exact sources select one variant. NPNT is splice-mediated and out of scope. |
| HAO2020_DSP | NA | No complete screening denominator beyond deliberately selected mechanistic examples. HHIP enumerates all four isolated Figure3B tests (two selected SNPs in one locus, each in two insert orientations); per-orientation descriptive counts are supported, but these correlated case-series observations do not establish population-like discrimination. Other exact sources select one variant. NPNT is splice-mediated and out of scope. |
| LIN2020_FAM13A | NA | No complete screening denominator beyond deliberately selected mechanistic examples. HHIP enumerates all four isolated Figure3B tests (two selected SNPs in one locus, each in two insert orientations); per-orientation descriptive counts are supported, but these correlated case-series observations do not establish population-like discrimination. Other exact sources select one variant. NPNT is splice-mediated and out of scope. |
| SAFERALI2025_NPNT | NA | No complete screening denominator beyond deliberately selected mechanistic examples. HHIP enumerates all four isolated Figure3B tests (two selected SNPs in one locus, each in two insert orientations); per-orientation descriptive counts are supported, but these correlated case-series observations do not establish population-like discrimination. Other exact sources select one variant. NPNT is splice-mediated and out of scope. |
| ZHOU2012_HHIP | NA | No complete screening denominator beyond deliberately selected mechanistic examples. HHIP enumerates all four isolated Figure3B tests (two selected SNPs in one locus, each in two insert orientations); per-orientation descriptive counts are supported, but these correlated case-series observations do not establish population-like discrimination. Other exact sources select one variant. NPNT is splice-mediated and out of scope. |

GEO design/count completeness does not recover missing Castaldi significance labels or QC exclusions; unidentified hits are not nulls. Gong null assay-context rows are not universal variant negatives, and missing/failed contexts do not become null. No classification or enrichment statistic is manufactured from these incomplete denominators.

## Real-network RC invariance and compute provenance

RC audit status: PASS; 6 seed audits and 2 ensemble audits. Maximum absolute residual: 0. Tolerance: abs(a−b) <= 1e-06 + 1e-06*abs(b).

The audit independently repeated nucleotide conversion and actual phase-I and phase-II network calls in RC-first/forward-second order. It was not constructed by swapping stored probabilities or reversing phase-I feature columns.

| Context | Unit | Seed | Probabilities | Failures | Maximum absolute residual |
| --- | --- | --- | --- | --- | --- |
| enhancer | seed | 104729 | 3420 | 0 | 0 |
| enhancer | seed | 130363 | 3420 | 0 | 0 |
| enhancer | seed | 155921 | 3420 | 0 | 0 |
| enhancer | ensemble | all3 | 3420 | 0 | 0 |
| h3k27me3 | seed | 104729 | 3420 | 0 | 0 |
| h3k27me3 | seed | 130363 | 3420 | 0 | 0 |
| h3k27me3 | seed | 155921 | 3420 | 0 | 0 |
| h3k27me3 | ensemble | all3 | 3420 | 0 | 0 |

| Compute item | Recorded value |
| --- | --- |
| Variants | 1710 |
| Allele sequences | 3420 |
| Original selected checkpoints read, not copied/re-saved | 6 |
| Variant/context scores | 3420 |
| Seed/context/variant scores | 10260 |
| Phase-I representation | final4560sigmoid |
| Phase-I feature shape per orientation | [3420, 4560] |
| Phase-I actual sequence-orientation evaluations | 13680 |
| Phase-I wall seconds | 7.562443354050629 |
| Full inference wall seconds | 41.18218468106352 |
| Inference peak RSS KiB | 2431012 |
| Descriptive evaluation wall seconds | 3.9242368310224265 |

Exact recorded inference invocation:

```json
[
  "diseases/COPD/07_gap_closure/external-benchmark-v2-1.0/scripts/score_frozen_v2.py",
  "--repo",
  "/vf/users/Dcode/yash/Codex-Manuscript/disease-regulatory-genomics-workflow",
  "--stage",
  "diseases/COPD/07_gap_closure/external-benchmark-v2-1.0"
]
```

Recorded software/environment (from this run, not an inferred environment):

```json
{
  "environment": {
    "CUDA_VISIBLE_DEVICES": "0",
    "PYTHONHASHSEED": "104729",
    "SLURM_ARRAY_TASK_ID": null,
    "SLURM_CPUS_PER_TASK": "16",
    "SLURM_JOB_GPUS": "3",
    "SLURM_JOB_ID": "32057893",
    "TF_CUDNN_DETERMINISTIC": "1",
    "TF_DETERMINISTIC_OPS": "1",
    "TF_ENABLE_ONEDNN_OPTS": null,
    "TF_XLA_FLAGS": "--tf_xla_enable_xla_devices=false"
  },
  "executable": "/home/maheshwarany2/.local/share/uv/python/cpython-3.13.0-linux-x86_64-gnu/bin/python3.13",
  "keras_epsilon": 1e-07,
  "keras_floatx": "float32",
  "logical_devices": [
    "LogicalDevice(name='/device:CPU:0', device_type='CPU')",
    "LogicalDevice(name='/device:GPU:0', device_type='GPU')"
  ],
  "nvidia_smi": {
    "returncode": 0,
    "stderr": "",
    "stdout": "NVIDIA A100-SXM4-80GB, GPU-1f1b65b3-e559-713a-e355-9dd901bfd8e6, 580.173.02, 81920 MiB\n"
  },
  "packages": {
    "h5py": "3.16.0",
    "keras": "3.14.1",
    "numpy": "2.5.0",
    "scikit-learn": "1.9.0",
    "scipy": "1.18.0",
    "tensorflow": "2.20.0"
  },
  "platform": "Linux-4.18.0-425.19.3b.el8.x86_64-x86_64-with-glibc2.28",
  "python": "3.13.0",
  "tensorflow_build": {
    "cpu_compiler": "clang 18",
    "cuda_compute_capabilities": [
      "sm_60",
      "sm_70",
      "sm_80",
      "sm_89",
      "compute_90"
    ],
    "cuda_version": "12.5.1",
    "cudnn_version": "9",
    "is_cuda_build": true,
    "is_rocm_build": false,
    "is_tensorrt_build": false
  }
}
```

### Command and implementation-attempt ledger


| Command record | Return code | Wall seconds | Recorded argv |
| --- | --- | --- | --- |
| evaluator_synthetic_selftest.completed.json | 0 | 0.3054830179316923 | ["bash", "diseases/COPD/07_gap_closure/external-benchmark-v2-1.0/scripts/runtime.sh", "diseases/COPD/07_gap_closure/external-benchmark-v2-1.0/scripts/evaluate_outputs.py", "--self-test"] |
| frozen_V2_C_external_inference.completed.json | 0 | 44.421370810014196 | ["bash", "diseases/COPD/07_gap_closure/external-benchmark-v2-1.0/scripts/runtime.sh", "diseases/COPD/07_gap_closure/external-benchmark-v2-1.0/scripts/score_frozen_v2.py", "--repo", "/vf/users/Dcode/yash/Codex-Manuscript/disease-regulatory-genomics-workflow", "--stage", "diseases/COPD/07_gap_closure/external-benchmark-v2-1.0"] |
| independent_validation_pre_final.completed.json | 0 | 14.711761220009066 | ["bash", "diseases/COPD/07_gap_closure/external-benchmark-v2-1.0/scripts/runtime.sh", "diseases/COPD/07_gap_closure/external-benchmark-v2-1.0/scripts/validate_external_stage.py", "--phase", "pre-final"] |
| independent_validator_synthetic_tests.completed.json | 0 | 0.3647142550908029 | ["bash", "diseases/COPD/07_gap_closure/external-benchmark-v2-1.0/scripts/runtime.sh", "-m", "unittest", "discover", "-s", "diseases/COPD/07_gap_closure/external-benchmark-v2-1.0/tests", "-p", "test_independent_validation.py", "-v"] |
| prespecified_external_evaluation.completed.json | 0 | 4.531564614037052 | ["bash", "diseases/COPD/07_gap_closure/external-benchmark-v2-1.0/scripts/runtime.sh", "diseases/COPD/07_gap_closure/external-benchmark-v2-1.0/scripts/evaluate_outputs.py"] |

Exclusive-output guards and pre-output assertion failures were implementation checks, not opportunities to alter scientific rules. The annotation and sequence-preparation first attempts exposed the same unresolved reporting-key allele-order convention; sorting the display-key pair matched the already frozen source summary without rescuing sequence identity or changing eligibility. Failed script versions and attempt ledgers are retained.

Annotation attempt ledger:

```json
[
  {
    "attempt": 1,
    "cause": "Frozen R015 unresolved reporting keys sort tested alleles; first implementation retained reported allele order.",
    "command": "python3 diseases/COPD/07_gap_closure/external-benchmark-v2-1.0/scripts/annotate_benchmark.py",
    "error": "AssertionError: unresolved:rs201292172:ACAC/-",
    "exit_code": 1,
    "failed_script": {
      "bytes": 18576,
      "path": "scripts/history/annotate_benchmark_attempt001.py",
      "sha256": "751a5081a3c77043cd962d210e71ef479f94f538bf3b36e67787889a4b704d23"
    },
    "repair": "Sort the two tested-allele labels only when recreating an unresolved R015 reporting key; no sequence identity rescued, no eligibility/direction rule changed.",
    "status": "FAILED_BEFORE_OUTPUTS"
  },
  {
    "attempt": 2,
    "command": "python3 diseases/COPD/07_gap_closure/external-benchmark-v2-1.0/scripts/annotate_benchmark.py",
    "exit_code": 0,
    "status": "PASS"
  }
]
```

Sequence-preparation attempt ledger:

```json
{
  "stage": "external-benchmark-v2-1.0",
  "attempts": [
    {
      "attempt": 1,
      "status": "FAILED_BEFORE_OUTPUT",
      "exit_code": 1,
      "original_script_sha256": "58f3feb3aa657047a58b036a6a47ae1c4d398fc6c3c0c055cfbea18f95c3eb19",
      "archived_source": "scripts/history/prepare_external_inputs_initial.py",
      "command": "env PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONPATH=/vf/users/Dcode/yash/Codex-Manuscript/disease-regulatory-genomics-workflow/models/TREDNET_v2/.venv/lib/python3.13/site-packages /home/maheshwarany2/.local/share/uv/python/cpython-3.13.0-linux-x86_64-gnu/bin/python3.13 diseases/COPD/07_gap_closure/external-benchmark-v2-1.0/scripts/prepare_external_inputs.py",
      "error": "AssertionError: unresolved:rs201292172:ACAC/-; at original source line204 assert key in unique_keys",
      "cause": "Unresolved display-key construction initially retained reported tested-allele order; frozen R015 sorts that pair.",
      "repair": "Retain the exact already-frozen R015 unresolved display key using sorted reported pair. This does not resolve identity, swap genomic REF/ALT, change eligibility or alter sequence/model rules.",
      "scientific_outputs_created": false,
      "models_or_predictions_accessed": false,
      "prior_artifacts_changed": false
    }
  ],
  "successful_attempt_receipt": "provenance/sequence_construction_validation.json",
  "prospective_scientific_specification_changed": false
}
```

## Independent validation, preservation and freeze

Independent pre-final scientific validation status: **PASS**. Receipt: `provenance/independent_validation_pre-final.json`; SHA-256 `9ce53af29adedd422130275fbabb4799d2cca36481fc494d87c3b48f4d8d0b1b`.

Authoritative final validation is recorded separately in `provenance/independent_final_validation.json` and its accompanying check table after this report is generated. That receipt binds the completed report, results and append-only register records before the final seal. Consult that receipt for the final status and exact verified scope. Its hash/status is not embedded in this report, avoiding a circular report-to-validator hash dependency; the finished report is not rewritten after final validation.

All original benchmark, pretraining, internal-training, internal-test, V1 and checkpoint artifacts are read-only inputs. Preservation and append-only register assertions are documented in the final validation/freeze records. No broader COPD candidate universe, 337-candidate reranking, fine-mapping or target-gene analysis was accessed or executed. The six retained V2-C checkpoint archives were reused in place, not retrained, re-saved, slimmed or duplicated.

The final checksum ledger and freeze record bind this report, inputs, predictions, results, scripts, environment and provenance after validation. The final freeze is intentionally created after this report and therefore its self-hash is not embedded here. Only the three V2 registers may receive append-only stage updates. This evaluation request does not authorize a commit or push.

## Interpretation and remaining failure modes

The corrected V2-C ensemble shows higher descriptive external region-level union recovery than both frozen V1 orientations on the same evaluable positives.
The common-identity union result is higher than V1_forward (3/23 versus 2/23); higher than V1_RC (3/23 versus 2/23).

Region misses remain explicitly listed in the per-variant and transition tables. Direction discordance (4/11 observations) must be reported separately from region recovery. Small allelic deltas have no validated decision threshold here. Improved regional recognition, if observed, cannot establish better causal allelic ranking, repression mechanism, clinical prediction or generalization to unassayed COPD variants.

The major remaining limitations are incomplete experimental denominators; unresolved exact identities; cell/assay and episomal-versus-endogenous context mismatch; correlated loci and repeated contexts; a small direction-resolved subset; no H3-specific experimental direction denominator; V1-dependent retrospective awareness including rs2013701; and out-of-model mechanisms such as splicing. See `LIMITATIONS.md` for the explicit interpretation boundaries. This stage stops here; any future redesign is a separately versioned scientific project and cannot call this benchmark untouched validation.
