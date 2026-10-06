# COPD V2 frozen chr8–9 internal holdout evaluation

Version: `internal-test-1.0`. Models: unchanged published `internal-training-1.0`; test construction: unchanged `pretraining-1.1`.

## Decision and stop boundary

READY FOR INVESTIGATOR REVIEW BEFORE ANY EXTERNAL FUNCTIONAL EVALUATION. Both C contexts pass the prospective holdout rule and all five recorded independent validation phases pass. This is a review disposition only: separate authorization is still required for external evaluation or candidate/allele scoring.

| C context | Prospective holdout decision | Failed or inconclusive conditions |
| --- | --- | --- |
| Enhancer-associated | PASS | none |
| H3K27me3-associated | PASS | none |

Chr8–9 was held out from V2 training, V2 checkpoint selection and chr7 calibration. It is **not historically pristine**: V1 used these chromosomes and V2 construction/QC inspected outcome-blind covariates. This is a **fixed internal V2 region-label holdout, not independent external validation**. No result changes the model configuration, matching, seed list, checkpoint, threshold or prospective contract.

The external functional benchmark remains unopened and unparsed. No COPD variant universe, 337 V1 candidates, 12-shortlist variants, rs2013701 or other special functional variants, or arbitrary REF/ALT sequences were scored. No retraining, fine-mapping, target-gene analysis, manuscript revision or other V2 module was performed. No Git commit or push is authorized by this stage request.

## Prospective contract and provenance

| Record | Value |
| --- | --- |
| Prospective freeze UTC | 2026-10-06T18:59:34.340609+00:00 |
| Before any test inference | yes |
| Baseline published commit | 07d72dad0f57d605f5f50bbc7fa6d2c7843a94ba |
| Machine specification SHA-256 | 55e5a7f340596244dd45150ce393b645bcc000d64e34c19c7c517f784e41dcca |
| Readable specification SHA-256 | 5c279df5c09699f3dbb70a953502e7cd994ed4fe093a7b4fdacad8a1bc8d1a19 |
| Prior training freeze SHA-256 | c151e61dba04a9e5500534e507a392d51e1af0601a6e62565288bc6b276e6de4 |

The [prospective specification](TEST_SPECIFICATION.md), [machine contract](specification/test_specification.json), [pre-inference freeze](provenance/prospective_specification_freeze.json), and [preserved authorization](provenance/user_authorization.txt) define this stage. Hash-bound gates precede input preparation, phase-I extraction, checkpoint inference and one-time evaluation. Model-performance outcomes were not used to revise the contract.

## Frozen population, membership and sequence eligibility

| Context | Positive | Control | Rows | Global components | chr8 | chr9 | Input audit |
| --- | --- | --- | --- | --- | --- | --- | --- |
| Enhancer-associated | 12695 | 12695 | 25390 | 8754 | 11963 | 13427 | PASS |
| H3K27me3-associated | 1754 | 1754 | 3508 | 1809 | 1724 | 1784 | PASS |

The union contains 26,225 unique genomic intervals and 26,225 canonical encoded identities in 8,926 existing global overlap/encoded-sequence components. The two contexts share 2,673 interval identities. All selected source rows, cell strings, order, labels, roles, components and train-derived covariates are preserved; no matching was rebuilt and no row was substituted.

| Frozen consumed input | Bytes | SHA-256 |
| --- | --- | --- |
| diseases/COPD/07_gap_closure/pretraining-1.1/data/evaluation/enhancer_common_challenge_panel.tsv.gz | 4962271 | 4f374ebf993b1e3b6a5fff5aba51ac6e0efaf6859300bf752f56245a7db65c20 |
| diseases/COPD/07_gap_closure/pretraining-1.1/data/evaluation/h3k27me3_common_challenge_panel.tsv.gz | 719591 | 895979a83b12d0bb372dd93a1baaa1a83d174e96a57e17acc4f5f3bea8775c26 |
| diseases/COPD/07_gap_closure/pretraining-1.1/data/configurations/V2-C_enhancer_interval_manifest.tsv.gz | 44914883 | 5165080dbeba1a43decb57286c9ed44e820eaa73617f252cecda6ebe4321e72f |
| diseases/COPD/07_gap_closure/pretraining-1.1/data/configurations/V2-C_h3k27me3_interval_manifest.tsv.gz | 6293534 | 1d325e61ae5c41bfabe955c12e0b84f683a80472fb855dbac9a19fabd1cbc7ec |
| diseases/COPD/04_modeling/trednet/fasta/hg38.fa | 3273481150 | 5be01555d98347fdb3714dc84c6f77c9d8bc774adcf32c6f7a8fa06f5baf5e51 |
| diseases/COPD/04_modeling/trednet/fasta/hg38.fa.fai | 19381 | 3b425de206296a5c8053023fa5ca61da43cfe78c1737c12e58c83367c7e83c21 |

[Input integrity details](inputs/test_panel_integrity_audit.json) and the [input manifest](inputs/input_manifest.json) record 1,000-bp cores, `[core_start−501, core_end+500)` 2,001-bp inputs, exact hg38/FAI hashes, raw sequence hashes, canonical nucleotide/RC hashes, GC and ambiguous-base checks. FASTA access is limited to the exact frozen test intervals on chr8 and chr9.

## Frozen phase-I extraction and unchanged checkpoint inference

| Record | Observed |
| --- | --- |
| Sequence geometry / encoding | [26225, 2001] / uint8_ASCII |
| Per-orientation feature geometry / dtype | [26225, 4560] / float32 |
| Representation | final4560sigmoid |
| Independent nucleotide-orientation phase-I evaluations | 52450 |
| Missing sequences / extraction failures | 0 / 0 |
| Phase-I weight SHA-256 | 483b6c0cafd750ea8e97932cbc8a949276eaf1e77e7d8766f9a22eeddc2916d6 |

Non-ACGT bases are normalized to N and all-zero one-hot encoding. Each canonical nucleotide sequence and its nucleotide reverse complement are independently evaluated through the original frozen phase-I network; representation columns are never reversed. The orientation map restores genomic-forward and nucleotide-RC inputs for each interval.

| Test cache | Bytes | SHA-256 |
| --- | --- | --- |
| cache/features_canonical.npy | 478344128 | 6681952719517712d123cd753c388b04049b88251ac5c434875ad4a98d32ed6e |
| cache/features_rc.npy | 478344128 | 2b99b0a2509c07af7ba5141ec7f67ca49b34c97f78285d3be5a190376985dcaa |

All 18 original selected `.keras` archives were loaded unchanged with `compile=False` and verified before and after inference. No retraining, checkpoint replacement, re-saving, optimizer-state editing or favorable-seed selection occurred. Each seed uses `q_s(x) = [float64(p_s(x)) + float64(p_s(RC(x)))] / 2`; the final score is the equal-three-seed float64 mean in order 104729, 130363, 155921. Actual network operations remain float32.

| Configuration | Context | Seed | Original checkpoint SHA-256 |
| --- | --- | --- | --- |
| V2-A | enhancer | 104729 | 8eb380c37e92f5cbe9cb71256421565f5924d1fe373ce7ee507dbb3545272a60 |
| V2-A | enhancer | 130363 | fb1f9103e7c824ffa280aa98fb2274a0abdb5eebbfcc381b2e8ed6df030faeed |
| V2-A | enhancer | 155921 | a74ee692f8b7386678a37de827e5ecb1bd597f665c6b475ba19fddadea961168 |
| V2-A | h3k27me3 | 104729 | 67b6c55c2dc7c79b8e70331ce4ee17c2f6b228ae95198312c8476084723b6249 |
| V2-A | h3k27me3 | 130363 | d20a5e35ac9e348ea09b1418abf18fd3a4c5c040811a8324d25ff39911b25ae8 |
| V2-A | h3k27me3 | 155921 | b350d45b7d126adf329c8379e5340a8cdd8b95c4028c00a99e726513c7a03820 |
| V2-B | enhancer | 104729 | 1dc13e371414fd81072b07689f3ce64446a32692e917d236983db024f28cabca |
| V2-B | enhancer | 130363 | 97bbc39bc5e93bf49b51169ad015fb4710708c7f1d3e7aeb232d28d3fa40d9ab |
| V2-B | enhancer | 155921 | c0e2a96ed19c3b9d5ddbb51a85b0d4ffce7701e41f12bd25ce69bcf763d6300e |
| V2-B | h3k27me3 | 104729 | c74842f069069ef544fbf8ddab04690188282a07c2aa7ebf2a2a67994d567e0c |
| V2-B | h3k27me3 | 130363 | 935e054078b14a491d69e89719a68a2fa359beb42afce814a5a496d74c95cb72 |
| V2-B | h3k27me3 | 155921 | 2e2fa2ec39cbaf8fea4080bd0576d6a1da6d968e2481c349d737f3db641553f6 |
| V2-C | enhancer | 104729 | 1cb670f778a2d148d0ac2b71bea5c97f5a421bb64f744af1dc5b373a60ea3114 |
| V2-C | enhancer | 130363 | 0e42ecc9364f967220ecc6c3bb0007dc25aced5d9c360355f6e48027c39f0795 |
| V2-C | enhancer | 155921 | 733ca1845061508c820240d04d4b7c30d8ba0fa7c1086419eaec0b5c2bb041d5 |
| V2-C | h3k27me3 | 104729 | 503683295ebeb48587f42748181bf6c5815f6a8c79b8dfdd4b669c931109be81 |
| V2-C | h3k27me3 | 130363 | 56a2d331b6ed27e49f103dce228dda3c26a80ff296e91c9c740e52b9c36841c6 |
| V2-C | h3k27me3 | 155921 | 995dc295d163a703b80bb3e451b8fc3bd7af01fb6ce041099fc5ce3d9cc82560 |

## Actual-network reverse-complement invariance

Saved invariance status: **PASS**; 18 seed and 6 ensemble audits. Every frozen test-panel sequence is covered. Additional real-network calls with swapped nucleotide-orientation inputs use the same row order and batching; the second wrapper is not obtained by algebraically recycling saved probabilities. The frozen tolerance is `abs(Q(x)−Q(RC(x))) <= 1e−6 + 1e−6*abs(Q(RC(x)))`. All raw, repeated, seed-symmetric and ensemble predictions must remain finite within [0,1].

| Configuration | Context | Unit | Seed | Sequences | Status | Failed | Maximum absolute residual |
| --- | --- | --- | --- | --- | --- | --- | --- |
| V2-A | enhancer | seed | 104729 | 25390 | PASS | 0 | 0 |
| V2-A | enhancer | seed | 130363 | 25390 | PASS | 0 | 0 |
| V2-A | enhancer | seed | 155921 | 25390 | PASS | 0 | 0 |
| V2-A | enhancer | ensemble | all3 | 25390 | PASS | 0 | 0 |
| V2-B | enhancer | seed | 104729 | 25390 | PASS | 0 | 0 |
| V2-B | enhancer | seed | 130363 | 25390 | PASS | 0 | 0 |
| V2-B | enhancer | seed | 155921 | 25390 | PASS | 0 | 0 |
| V2-B | enhancer | ensemble | all3 | 25390 | PASS | 0 | 0 |
| V2-C | enhancer | seed | 104729 | 25390 | PASS | 0 | 0 |
| V2-C | enhancer | seed | 130363 | 25390 | PASS | 0 | 0 |
| V2-C | enhancer | seed | 155921 | 25390 | PASS | 0 | 0 |
| V2-C | enhancer | ensemble | all3 | 25390 | PASS | 0 | 0 |
| V2-A | h3k27me3 | seed | 104729 | 3508 | PASS | 0 | 0 |
| V2-A | h3k27me3 | seed | 130363 | 3508 | PASS | 0 | 0 |
| V2-A | h3k27me3 | seed | 155921 | 3508 | PASS | 0 | 0 |
| V2-A | h3k27me3 | ensemble | all3 | 3508 | PASS | 0 | 0 |
| V2-B | h3k27me3 | seed | 104729 | 3508 | PASS | 0 | 0 |
| V2-B | h3k27me3 | seed | 130363 | 3508 | PASS | 0 | 0 |
| V2-B | h3k27me3 | seed | 155921 | 3508 | PASS | 0 | 0 |
| V2-B | h3k27me3 | ensemble | all3 | 3508 | PASS | 0 | 0 |
| V2-C | h3k27me3 | seed | 104729 | 3508 | PASS | 0 | 0 |
| V2-C | h3k27me3 | seed | 130363 | 3508 | PASS | 0 | 0 |
| V2-C | h3k27me3 | seed | 155921 | 3508 | PASS | 0 | 0 |
| V2-C | h3k27me3 | ensemble | all3 | 3508 | PASS | 0 | 0 |

## Primary common-panel A/B/C performance

All A/B/C models use the identical frozen C-task test panel within each context. AP is stepwise tie-grouped average precision, not trapezoidal PR-AUC; AUROC is tie-aware. Brier is mean squared probability error and Brier skill is `1 − Brier/[prevalence × (1−prevalence)]`. Display rounding does not replace full-precision TSV values. Entries below are point estimates [95% component-bootstrap confidence intervals].

| Context | Configuration | AP [95% CI] | AUROC [95% CI] | Brier [95% CI] | Brier skill [95% CI] |
| --- | --- | --- | --- | --- | --- |
| enhancer | V2-A | 0.7320958 [0.7135907, 0.7508474] | 0.7206394 [0.7083017, 0.7325479] | 0.2448588 [0.2378633, 0.2517026] | 0.02056486 [-0.006956272, 0.04812761] |
| enhancer | V2-B | 0.7994381 [0.7835919, 0.8151797] | 0.8004623 [0.7902778, 0.8104777] | 0.1854083 [0.1815867, 0.189224] | 0.2583667 [0.2429022, 0.2734466] |
| enhancer | V2-C | 0.7998119 [0.7839603, 0.8154022] | 0.8006103 [0.7902516, 0.8106114] | 0.1836227 [0.1795429, 0.1877439] | 0.2655091 [0.2490104, 0.2816364] |
| h3k27me3 | V2-A | 0.7499031 [0.6993771, 0.7925926] | 0.7192883 [0.6889317, 0.7488481] | 0.3148171 [0.2930947, 0.338965] | -0.2592684 [-0.3621814, -0.1765648] |
| h3k27me3 | V2-B | 0.81105 [0.7715478, 0.8448802] | 0.7994179 [0.7728915, 0.825429] | 0.1823339 [0.1709968, 0.193585] | 0.2706642 [0.2251271, 0.3153813] |
| h3k27me3 | V2-C | 0.8096768 [0.7698547, 0.8438564] | 0.7989697 [0.7722317, 0.8257346] | 0.1837016 [0.1711269, 0.1960853] | 0.2651937 [0.2150257, 0.3148315] |

### All frozen seed-specific symmetric predictors

| Context | Configuration | Predictor | AP [95% CI] | AUROC [95% CI] | Brier [95% CI] | Brier skill [95% CI] |
| --- | --- | --- | --- | --- | --- | --- |
| enhancer | V2-A | seed:104729 | 0.7310448 [0.7118493, 0.7500443] | 0.7177824 [0.7051502, 0.7299501] | 0.2474273 [0.2403456, 0.2543065] | 0.01029079 [-0.01766999, 0.03815818] |
| enhancer | V2-A | seed:130363 | 0.7303865 [0.7118026, 0.7495124] | 0.7211938 [0.7089662, 0.7331011] | 0.2440302 [0.2371039, 0.2508407] | 0.02387938 [-0.00382816, 0.05139072] |
| enhancer | V2-A | seed:155921 | 0.7324883 [0.7144097, 0.7508738] | 0.721664 [0.7093951, 0.7333641] | 0.2441907 [0.2372343, 0.2509505] | 0.02323716 [-0.003961557, 0.0507643] |
| enhancer | V2-B | seed:104729 | 0.7978786 [0.7819047, 0.8138024] | 0.7984715 [0.7882875, 0.8085731] | 0.1865776 [0.182721, 0.1903441] | 0.2536896 [0.238347, 0.2688099] |
| enhancer | V2-B | seed:130363 | 0.7997912 [0.78409, 0.8154532] | 0.8012538 [0.7910138, 0.8111712] | 0.184747 [0.1809342, 0.1885634] | 0.2610119 [0.245605, 0.2760974] |
| enhancer | V2-B | seed:155921 | 0.7993391 [0.7834246, 0.815187] | 0.8001669 [0.7899869, 0.8101723] | 0.1853743 [0.1814604, 0.189195] | 0.2585027 [0.2429522, 0.2739349] |
| enhancer | V2-C | seed:104729 | 0.7994352 [0.7836748, 0.8149559] | 0.8001262 [0.7899023, 0.8100133] | 0.1847329 [0.1807486, 0.1886649] | 0.2610686 [0.2451673, 0.2768507] |
| enhancer | V2-C | seed:130363 | 0.7996945 [0.7841444, 0.8153015] | 0.8001076 [0.7898055, 0.8102037] | 0.1834402 [0.1793092, 0.187627] | 0.2662393 [0.2494672, 0.2826692] |
| enhancer | V2-C | seed:155921 | 0.7989609 [0.7828587, 0.8148049] | 0.7998799 [0.7897061, 0.8098758] | 0.1834917 [0.179366, 0.1876147] | 0.2660333 [0.2495198, 0.2823029] |
| h3k27me3 | V2-A | seed:104729 | 0.7478234 [0.696693, 0.7914425] | 0.717946 [0.6875113, 0.7474299] | 0.3161176 [0.2942032, 0.3404536] | -0.2644704 [-0.3670419, -0.1819438] |
| h3k27me3 | V2-A | seed:130363 | 0.7447603 [0.6929482, 0.7890607] | 0.716006 [0.6848753, 0.7455497] | 0.3201798 [0.2979357, 0.3449291] | -0.2807191 [-0.3857233, -0.19556] |
| h3k27me3 | V2-A | seed:155921 | 0.7522092 [0.7024753, 0.7943072] | 0.7209363 [0.690544, 0.7502852] | 0.3089573 [0.2877867, 0.3324125] | -0.2358293 [-0.3354135, -0.1557013] |
| h3k27me3 | V2-B | seed:104729 | 0.8119052 [0.7730829, 0.8454855] | 0.7998258 [0.7736241, 0.8261643] | 0.1824219 [0.1709989, 0.1933524] | 0.2703124 [0.2252755, 0.3150155] |
| h3k27me3 | V2-B | seed:130363 | 0.8110815 [0.7715307, 0.8449755] | 0.7995593 [0.7729403, 0.8256587] | 0.1821486 [0.1706336, 0.1935318] | 0.2714056 [0.2250846, 0.315965] |
| h3k27me3 | V2-B | seed:155921 | 0.8078711 [0.7684646, 0.8417437] | 0.7966082 [0.769882, 0.8225495] | 0.1834314 [0.172185, 0.1947865] | 0.2662743 [0.2205321, 0.3105475] |
| h3k27me3 | V2-C | seed:104729 | 0.8089028 [0.7686711, 0.8429799] | 0.7981135 [0.7714981, 0.8249224] | 0.1840044 [0.1714125, 0.1963872] | 0.2639823 [0.213546, 0.3140195] |
| h3k27me3 | V2-C | seed:130363 | 0.8097726 [0.7698093, 0.8439955] | 0.798694 [0.7717801, 0.8254348] | 0.1837615 [0.1713356, 0.1959391] | 0.2649538 [0.215548, 0.3141748] |
| h3k27me3 | V2-C | seed:155921 | 0.8082443 [0.7680562, 0.8432533] | 0.7975921 [0.7707572, 0.8243459] | 0.1844729 [0.1719816, 0.1968346] | 0.2621085 [0.2120334, 0.3113565] |

### Component-bootstrap validity and scope

| Context | Components | Attempted | Valid | Invalid | Invalid fraction | Status |
| --- | --- | --- | --- | --- | --- | --- |
| enhancer | 8754 | 2000 | 2000 | 0 | 0 | PASS |
| h3k27me3 | 1809 | 2000 | 2000 | 0 | 0 | PASS |

The bootstrap uses the existing global full-input overlap/encoded-sequence component IDs, not matching pairs. Within each context all seed and ensemble A/B/C predictors and fixed-threshold C operating rates share the same draws: lexical component ordering, fresh NumPy default_rng(314159), 2,000 valid replicates, at most 20,000 attempts, and invalid/attempted strictly below 0.10. CIs use linear 2.5/97.5 percentiles. AP-minus-prevalence and Brier skill use each replicate's own weighted prevalence. Uncertainty is conditional on this internal genomic-component/fixed-matching design, not independent-donor or population uncertainty.

### Paired diagnostic ablations

| Context | Paired difference | Metric | Point [95% paired CI] |
| --- | --- | --- | --- |
| enhancer | V2-B minus V2-A | AP | 0.06734227 [0.05922835, 0.07553441] |
| enhancer | V2-B minus V2-A | AUROC | 0.07982288 [0.07171799, 0.08799002] |
| enhancer | V2-B minus V2-A | Brier | -0.05945047 [-0.06445239, -0.05461306] |
| enhancer | V2-C minus V2-B | AP | 0.0003738273 [-0.0003168871, 0.001082361] |
| enhancer | V2-C minus V2-B | AUROC | 0.0001480488 [-0.0005483908, 0.000820196] |
| enhancer | V2-C minus V2-B | Brier | -0.001785583 [-0.002199116, -0.001364004] |
| enhancer | V2-C minus V2-A | AP | 0.06771609 [0.05929479, 0.07628726] |
| enhancer | V2-C minus V2-A | AUROC | 0.07997092 [0.07174558, 0.08814995] |
| enhancer | V2-C minus V2-A | Brier | -0.06123605 [-0.06632305, -0.05617044] |
| h3k27me3 | V2-B minus V2-A | AP | 0.06114691 [0.0409381, 0.08408886] |
| h3k27me3 | V2-B minus V2-A | AUROC | 0.0801296 [0.05645418, 0.1045383] |
| h3k27me3 | V2-B minus V2-A | Brier | -0.1324831 [-0.1561556, -0.1112947] |
| h3k27me3 | V2-C minus V2-B | AP | -0.001373203 [-0.004090736, 0.001410924] |
| h3k27me3 | V2-C minus V2-B | AUROC | -0.0004482343 [-0.003281388, 0.002596417] |
| h3k27me3 | V2-C minus V2-B | Brier | 0.001367636 [-0.0007687577, 0.00344409] |
| h3k27me3 | V2-C minus V2-A | AP | 0.05977371 [0.04045699, 0.08108416] |
| h3k27me3 | V2-C minus V2-A | AUROC | 0.07968137 [0.05666227, 0.1048968] |
| h3k27me3 | V2-C minus V2-A | Brier | -0.1311155 [-0.1559795, -0.1086824] |

Enhancer-associated: V2-B minus V2-A: AP 0.06734227 (entire interval above zero); AUROC 0.07982288 (entire interval above zero); Brier -0.05945047 (entire interval below zero). V2-C minus V2-B: AP 0.0003738273 (interval includes zero); AUROC 0.0001480488 (interval includes zero); Brier -0.001785583 (entire interval below zero). V2-C minus V2-A: AP 0.06771609 (entire interval above zero); AUROC 0.07997092 (entire interval above zero); Brier -0.06123605 (entire interval below zero).

H3K27me3-associated: V2-B minus V2-A: AP 0.06114691 (entire interval above zero); AUROC 0.0801296 (entire interval above zero); Brier -0.1324831 (entire interval below zero). V2-C minus V2-B: AP -0.001373203 (interval includes zero); AUROC -0.0004482343 (interval includes zero); Brier 0.001367636 (interval includes zero). V2-C minus V2-A: AP 0.05977371 (entire interval above zero); AUROC 0.07968137 (entire interval above zero); Brier -0.1311155 (entire interval below zero).

For AP/AUROC, positive differences mean higher discrimination; for Brier, negative differences mean lower squared error. These comparisons are diagnostic ablations only. B→C changes both positive-label definition and independently rematched controls, so it is not a label-only contrast. C is not required to outperform B; no ordering changes the already-frozen model choice.

## C seed stability and prospective holdout-generalization safeguards

| Context | Configuration | AP seed 104729 | AP seed 130363 | AP seed 155921 | AP sample SD | AP range |
| --- | --- | --- | --- | --- | --- | --- |
| enhancer | V2-A | 0.7310448 | 0.7303865 | 0.7324883 | 0.001075042 | 0.002101744 |
| enhancer | V2-B | 0.7978786 | 0.7997912 | 0.7993391 | 0.0009996029 | 0.001912558 |
| enhancer | V2-C | 0.7994352 | 0.7996945 | 0.7989609 | 0.0003720143 | 0.000733605 |
| h3k27me3 | V2-A | 0.7478234 | 0.7447603 | 0.7522092 | 0.003744012 | 0.007448981 |
| h3k27me3 | V2-B | 0.8119052 | 0.8110815 | 0.8078711 | 0.002131466 | 0.00403409 |
| h3k27me3 | V2-C | 0.8089028 | 0.8097726 | 0.8082443 | 0.0007666086 | 0.001528361 |

For C, sample SD (`ddof=1`) ≤0.03 and range ≤0.10 are conjunctive release safeguards. These are internal generalization safeguards, not clinically meaningful effect-size thresholds. A/B stability is descriptive.

| Context | Prospective condition | Status | Observed | Frozen criterion |
| --- | --- | --- | --- | --- |
| enhancer | all_three_C_seeds_evaluable | PASS | 3 | 3 fixed seeds; no seed selection |
| enhancer | real_network_invariance | PASS | PASS | hash-bound actual-network all seeds/ensembles PASS |
| enhancer | finite_valid_probabilities | PASS | PASS | all orientations/seeds/ensembles finite in [0,1] |
| enhancer | exact_test_positive | PASS | 12695 | exactly 12695 |
| enhancer | exact_test_control | PASS | 12695 | exactly 12695 |
| enhancer | exact_test_components | PASS | 8754 | exactly 8754 |
| enhancer | bootstrap_valid_replicates | PASS | 2000 | exactly 2000 valid; at most 20000 attempts |
| enhancer | bootstrap_invalid_fraction | PASS | 0.0 | strictly <0.10 |
| enhancer | lower95_AUROC | PASS | 0.790251570476658 | strictly >0.5 |
| enhancer | lower95_AP_gain | PASS | 0.2872608458402373 | strictly >0.0 |
| enhancer | lower95_BrierSkill | PASS | 0.2490104267151466 | strictly >0.0 |
| enhancer | seed_AP_sample_SD | PASS | 0.0003720143229471607 | <=0.03; ddof=1 |
| enhancer | seed_AP_range | PASS | 0.0007336050444008979 | <=0.10 |
| h3k27me3 | all_three_C_seeds_evaluable | PASS | 3 | 3 fixed seeds; no seed selection |
| h3k27me3 | real_network_invariance | PASS | PASS | hash-bound actual-network all seeds/ensembles PASS |
| h3k27me3 | finite_valid_probabilities | PASS | PASS | all orientations/seeds/ensembles finite in [0,1] |
| h3k27me3 | exact_test_positive | PASS | 1754 | exactly 1754 |
| h3k27me3 | exact_test_control | PASS | 1754 | exactly 1754 |
| h3k27me3 | exact_test_components | PASS | 1809 | exactly 1809 |
| h3k27me3 | bootstrap_valid_replicates | PASS | 2000 | exactly 2000 valid; at most 20000 attempts |
| h3k27me3 | bootstrap_invalid_fraction | PASS | 0.0 | strictly <0.10 |
| h3k27me3 | lower95_AUROC | PASS | 0.7722317312841127 | strictly >0.5 |
| h3k27me3 | lower95_AP_gain | PASS | 0.27738755698445244 | strictly >0.0 |
| h3k27me3 | lower95_BrierSkill | PASS | 0.2150256954778419 | strictly >0.0 |
| h3k27me3 | seed_AP_sample_SD | PASS | 0.0007666085587933202 | <=0.03; ddof=1 |
| h3k27me3 | seed_AP_range | PASS | 0.0015283613320178357 | <=0.10 |

Any FAIL or INCONCLUSIVE context requires stopping for investigator review without redesign, retraining or external access. Even two PASS decisions do not themselves authorize the next module.

## Frozen chr7 region-score thresholds on the test panels

Apply the original C ensemble thresholds using `score >= threshold`; no chr8–9 score is used for calibration or adjustment.

| Context | Original decimal threshold | Exact float64 hex | TP | FP | TN | FN |
| --- | --- | --- | --- | --- | --- | --- |
| enhancer | 0.74848511815071117 | 0x1.7f39710000001p-1 | 4050 | 595 | 12100 | 8645 |
| h3k27me3 | 0.76960810025533055 | 0x1.8a0a12aaaaaacp-1 | 563 | 61 | 1693 | 1191 |

| Context | Operating metric | Point [95% component CI] | Defined draws | Undefined draws |
| --- | --- | --- | --- | --- |
| enhancer | sensitivity | 0.3190232 [0.299998, 0.3366452] | 2000 | 0 |
| enhancer | specificity | 0.9531312 [0.9464525, 0.9589969] | 2000 | 0 |
| enhancer | precision | 0.8719053 [0.8532232, 0.8893024] | 2000 | 0 |
| enhancer | NPV | 0.5832731 [0.5691379, 0.5969796] | 2000 | 0 |
| enhancer | accuracy | 0.6360772 [0.6242442, 0.6473099] | 2000 | 0 |
| enhancer | negative_row_FPR | 0.04686885 [0.04100311, 0.05354749] | 2000 | 0 |
| h3k27me3 | sensitivity | 0.3209806 [0.2745405, 0.3683645] | 2000 | 0 |
| h3k27me3 | specificity | 0.9652223 [0.952222, 0.9771075] | 2000 | 0 |
| h3k27me3 | precision | 0.9022436 [0.8615125, 0.936137] | 2000 | 0 |
| h3k27me3 | NPV | 0.5870319 [0.5536626, 0.6214815] | 2000 | 0 |
| h3k27me3 | accuracy | 0.6431015 [0.6146163, 0.6715902] | 2000 | 0 |
| h3k27me3 | negative_row_FPR | 0.03477765 [0.0228925, 0.04777802] | 2000 | 0 |

| Context | Observed chr7 negative-row FPR | Observed test negative-row FPR | Test minus chr7 | Test FPR >5% | Recalibrated |
| --- | --- | --- | --- | --- | --- |
| enhancer | 0.04949239 | 0.04686885 | -0.00262354 | False | False |
| h3k27me3 | 0.04942966 | 0.03477765 | -0.01465201 | False | False |

Any increase is reported as calibration drift. A test FPR above 5% is not automatically a protocol failure: the chr7 empirical 5% rule never guaranteed a test/population FPR. No recalibration is performed. These are region-label thresholds only—not variant-level FPR control, allele-effect thresholds or causal probabilities. Precision/NPV intervals omit undefined denominators with explicit counts; those omissions do not invalidate otherwise valid metric-bootstrap draws. Confusion counts are descriptive; rate CIs retain component dependence rather than assuming independent rows.

## Fixed reliability bins

All 24 predictors retain ten fixed bins `[0,.1), …, [.9,1]`, including empty bins. The ensemble bins are shown below; the full per-seed and ensemble table is linked in the artifact index. Observed fractions on this balanced, matched region panel are not population disease-risk calibration.

| Context | Configuration | Bin | Lower | Upper | Upper inclusive | Rows | Mean score | Observed positive fraction |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| enhancer | V2-A | 0 | 0.0 | 0.1 | False | 1430 | 0.05574171 | 0.2447552 |
| enhancer | V2-A | 1 | 0.1 | 0.2 | False | 1712 | 0.149864 | 0.2523364 |
| enhancer | V2-A | 2 | 0.2 | 0.3 | False | 1859 | 0.249838 | 0.2867133 |
| enhancer | V2-A | 3 | 0.3 | 0.4 | False | 1765 | 0.3510159 | 0.3348442 |
| enhancer | V2-A | 4 | 0.4 | 0.5 | False | 1817 | 0.4505291 | 0.3665382 |
| enhancer | V2-A | 5 | 0.5 | 0.6 | False | 1894 | 0.5496913 | 0.3917635 |
| enhancer | V2-A | 6 | 0.6 | 0.7 | False | 2219 | 0.6532712 | 0.4407391 |
| enhancer | V2-A | 7 | 0.7 | 0.8 | False | 2987 | 0.7528519 | 0.50385 |
| enhancer | V2-A | 8 | 0.8 | 0.9 | False | 4348 | 0.8535643 | 0.6078657 |
| enhancer | V2-A | 9 | 0.9 | 1.0 | True | 5359 | 0.9503088 | 0.7939914 |
| enhancer | V2-B | 0 | 0.0 | 0.1 | False | 419 | 0.07821916 | 0.04057279 |
| enhancer | V2-B | 1 | 0.1 | 0.2 | False | 1543 | 0.1534876 | 0.1017498 |
| enhancer | V2-B | 2 | 0.2 | 0.3 | False | 2425 | 0.2544073 | 0.1641237 |
| enhancer | V2-B | 3 | 0.3 | 0.4 | False | 3238 | 0.3520822 | 0.2745522 |
| enhancer | V2-B | 4 | 0.4 | 0.5 | False | 4097 | 0.4519533 | 0.3817427 |
| enhancer | V2-B | 5 | 0.5 | 0.6 | False | 4360 | 0.5495201 | 0.5428899 |
| enhancer | V2-B | 6 | 0.6 | 0.7 | False | 3643 | 0.648548 | 0.6733461 |
| enhancer | V2-B | 7 | 0.7 | 0.8 | False | 2786 | 0.747381 | 0.8018665 |
| enhancer | V2-B | 8 | 0.8 | 0.9 | False | 1749 | 0.8454291 | 0.8713551 |
| enhancer | V2-B | 9 | 0.9 | 1.0 | True | 1130 | 0.9506914 | 0.9663717 |
| enhancer | V2-C | 0 | 0.0 | 0.1 | False | 608 | 0.07460369 | 0.04440789 |
| enhancer | V2-C | 1 | 0.1 | 0.2 | False | 1867 | 0.1515708 | 0.115158 |
| enhancer | V2-C | 2 | 0.2 | 0.3 | False | 2701 | 0.2526746 | 0.1932618 |
| enhancer | V2-C | 3 | 0.3 | 0.4 | False | 3380 | 0.3513122 | 0.2926036 |
| enhancer | V2-C | 4 | 0.4 | 0.5 | False | 3723 | 0.4512698 | 0.4106903 |
| enhancer | V2-C | 5 | 0.5 | 0.6 | False | 3642 | 0.5489428 | 0.5494234 |
| enhancer | V2-C | 6 | 0.6 | 0.7 | False | 3381 | 0.6483627 | 0.6660751 |
| enhancer | V2-C | 7 | 0.7 | 0.8 | False | 2771 | 0.7474272 | 0.785637 |
| enhancer | V2-C | 8 | 0.8 | 0.9 | False | 1984 | 0.8460535 | 0.859879 |
| enhancer | V2-C | 9 | 0.9 | 1.0 | True | 1333 | 0.9521302 | 0.9579895 |
| h3k27me3 | V2-A | 0 | 0.0 | 0.1 | False | 171 | 0.04233741 | 0.3274854 |
| h3k27me3 | V2-A | 1 | 0.1 | 0.2 | False | 104 | 0.1459387 | 0.2403846 |
| h3k27me3 | V2-A | 2 | 0.2 | 0.3 | False | 94 | 0.2531863 | 0.3085106 |
| h3k27me3 | V2-A | 3 | 0.3 | 0.4 | False | 119 | 0.347818 | 0.4033613 |
| h3k27me3 | V2-A | 4 | 0.4 | 0.5 | False | 125 | 0.4502652 | 0.344 |
| h3k27me3 | V2-A | 5 | 0.5 | 0.6 | False | 211 | 0.5508105 | 0.2654028 |
| h3k27me3 | V2-A | 6 | 0.6 | 0.7 | False | 255 | 0.6534591 | 0.3529412 |
| h3k27me3 | V2-A | 7 | 0.7 | 0.8 | False | 313 | 0.7542372 | 0.3386581 |
| h3k27me3 | V2-A | 8 | 0.8 | 0.9 | False | 563 | 0.8569403 | 0.4156306 |
| h3k27me3 | V2-A | 9 | 0.9 | 1.0 | True | 1553 | 0.9688173 | 0.6870573 |
| h3k27me3 | V2-B | 0 | 0.0 | 0.1 | False | 125 | 0.07156304 | 0.104 |
| h3k27me3 | V2-B | 1 | 0.1 | 0.2 | False | 351 | 0.1503289 | 0.1623932 |
| h3k27me3 | V2-B | 2 | 0.2 | 0.3 | False | 439 | 0.2507093 | 0.23918 |
| h3k27me3 | V2-B | 3 | 0.3 | 0.4 | False | 546 | 0.3508419 | 0.3406593 |
| h3k27me3 | V2-B | 4 | 0.4 | 0.5 | False | 469 | 0.449632 | 0.4498934 |
| h3k27me3 | V2-B | 5 | 0.5 | 0.6 | False | 444 | 0.5485407 | 0.5788288 |
| h3k27me3 | V2-B | 6 | 0.6 | 0.7 | False | 322 | 0.6489115 | 0.6677019 |
| h3k27me3 | V2-B | 7 | 0.7 | 0.8 | False | 287 | 0.7494832 | 0.7804878 |
| h3k27me3 | V2-B | 8 | 0.8 | 0.9 | False | 206 | 0.8534242 | 0.8640777 |
| h3k27me3 | V2-B | 9 | 0.9 | 1.0 | True | 319 | 0.9414768 | 0.9655172 |
| h3k27me3 | V2-C | 0 | 0.0 | 0.1 | False | 167 | 0.07357995 | 0.1197605 |
| h3k27me3 | V2-C | 1 | 0.1 | 0.2 | False | 464 | 0.1486044 | 0.1724138 |
| h3k27me3 | V2-C | 2 | 0.2 | 0.3 | False | 520 | 0.2469972 | 0.2923077 |
| h3k27me3 | V2-C | 3 | 0.3 | 0.4 | False | 501 | 0.3508984 | 0.3592814 |
| h3k27me3 | V2-C | 4 | 0.4 | 0.5 | False | 462 | 0.4487459 | 0.5151515 |
| h3k27me3 | V2-C | 5 | 0.5 | 0.6 | False | 333 | 0.5440295 | 0.6366366 |
| h3k27me3 | V2-C | 6 | 0.6 | 0.7 | False | 244 | 0.6497862 | 0.6680328 |
| h3k27me3 | V2-C | 7 | 0.7 | 0.8 | False | 252 | 0.744904 | 0.7579365 |
| h3k27me3 | V2-C | 8 | 0.8 | 0.9 | False | 215 | 0.8519374 | 0.855814 |
| h3k27me3 | V2-C | 9 | 0.9 | 1.0 | True | 350 | 0.9443604 | 0.9542857 |

## Frozen stratified diagnostics

GC, train-derived ATAC rank and repeat quintiles reuse the exact frozen training-positive cutpoints with `searchsorted(side='right')`; no test-derived cutpoints or new subgroup definitions are introduced. Exact ATAC-lobe signatures, same-lobe histone-support signatures, chromosome and ambiguous-base status retain their frozen definitions. The table shows C ensembles; all seeds and A/B/C ensembles are retained in the full artifact. Every subgroup is descriptive. Cells below 100 rows are explicitly flagged, and one-class/undefined metrics remain empty (shown as —).

| Context | Stratifier | Stratum | Rows | Positive/control | Components | n<100 | AP | AUROC | Brier | Brier skill |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| enhancer | gc_fraction | quintile_1 | 5647 | 2769/2878 | 2924 | False | 0.7882223 | 0.7933227 | 0.1881329 | 0.2471877 |
| enhancer | gc_fraction | quintile_2 | 5656 | 2758/2898 | 2573 | False | 0.7614819 | 0.774447 | 0.1939365 | 0.2237786 |
| enhancer | gc_fraction | quintile_3 | 5092 | 2433/2659 | 2091 | False | 0.7557697 | 0.7837754 | 0.1903557 | 0.2370745 |
| enhancer | gc_fraction | quintile_4 | 4568 | 2324/2244 | 1597 | False | 0.8244181 | 0.8183327 | 0.1750643 | 0.299528 |
| enhancer | gc_fraction | quintile_5 | 4427 | 2411/2016 | 1173 | False | 0.8611107 | 0.8367932 | 0.1657793 | 0.3315612 |
| enhancer | atac_signal_percentile_max_train_only | quintile_1 | 4985 | 2513/2472 | 3120 | False | 0.8070814 | 0.8054384 | 0.1813997 | 0.2743523 |
| enhancer | atac_signal_percentile_max_train_only | quintile_2 | 5269 | 2456/2813 | 3533 | False | 0.7611632 | 0.7857365 | 0.1913977 | 0.2308784 |
| enhancer | atac_signal_percentile_max_train_only | quintile_3 | 5312 | 2585/2727 | 3399 | False | 0.7799706 | 0.7874256 | 0.1906482 | 0.2368618 |
| enhancer | atac_signal_percentile_max_train_only | quintile_4 | 4918 | 2470/2448 | 3111 | False | 0.8020976 | 0.8070023 | 0.1817666 | 0.2729189 |
| enhancer | atac_signal_percentile_max_train_only | quintile_5 | 4906 | 2671/2235 | 2562 | False | 0.8583465 | 0.8332855 | 0.1717851 | 0.3073893 |
| enhancer | repeat_fraction_2001 | quintile_1 | 4394 | 2396/1998 | 1581 | False | 0.8434695 | 0.8194639 | 0.1739386 | 0.2984903 |
| enhancer | repeat_fraction_2001 | quintile_2 | 5153 | 2663/2490 | 2242 | False | 0.8192863 | 0.8118576 | 0.1793756 | 0.2816881 |
| enhancer | repeat_fraction_2001 | quintile_3 | 5224 | 2560/2664 | 2461 | False | 0.7881778 | 0.7936328 | 0.1867139 | 0.2528481 |
| enhancer | repeat_fraction_2001 | quintile_4 | 5944 | 2899/3045 | 2727 | False | 0.7914479 | 0.8041108 | 0.1828985 | 0.2679645 |
| enhancer | repeat_fraction_2001 | quintile_5 | 4675 | 2177/2498 | 2356 | False | 0.7434515 | 0.7714118 | 0.1948729 | 0.216816 |
| enhancer | atac_lobe_signature | lower_left | 9825 | 4911/4914 | 6361 | False | 0.7772537 | 0.7762549 | 0.1938198 | 0.2247207 |
| enhancer | atac_lobe_signature | lower_left;lower_right | 254 | 130/124 | 247 | False | 0.8582917 | 0.8431762 | 0.1653174 | 0.3383611 |
| enhancer | atac_lobe_signature | lower_left;lower_right;upper_right | 10 | 5/5 | 10 | True | 0.9266667 | 0.92 | 0.1367054 | 0.4531785 |
| enhancer | atac_lobe_signature | lower_left;upper_right | 329 | 164/165 | 325 | False | 0.808987 | 0.8252772 | 0.1751046 | 0.2995749 |
| enhancer | atac_lobe_signature | lower_right | 7724 | 3861/3863 | 4857 | False | 0.8153064 | 0.8152216 | 0.1773269 | 0.2906924 |
| enhancer | atac_lobe_signature | lower_right;upper_right | 236 | 118/118 | 232 | False | 0.8155197 | 0.8024275 | 0.1808604 | 0.2765584 |
| enhancer | atac_lobe_signature | upper_right | 7012 | 3506/3506 | 4684 | False | 0.8109533 | 0.8148238 | 0.1774927 | 0.2900293 |
| enhancer | same_lobe_mark_support | lower_left | 4988 | 4988/0 | 2540 | False | — | — | — | — |
| enhancer | same_lobe_mark_support | lower_left;lower_right | 104 | 104/0 | 98 | False | — | — | — | — |
| enhancer | same_lobe_mark_support | lower_left;lower_right;upper_right | 3 | 3/0 | 3 | True | — | — | — | — |
| enhancer | same_lobe_mark_support | lower_left;upper_right | 95 | 95/0 | 95 | True | — | — | — | — |
| enhancer | same_lobe_mark_support | lower_right | 3894 | 3894/0 | 1937 | False | — | — | — | — |
| enhancer | same_lobe_mark_support | lower_right;upper_right | 89 | 89/0 | 88 | True | — | — | — | — |
| enhancer | same_lobe_mark_support | none | 12695 | 0/12695 | 5877 | False | — | — | — | — |
| enhancer | same_lobe_mark_support | upper_right | 3522 | 3522/0 | 1833 | False | — | — | — | — |
| enhancer | chromosome | chr8 | 11963 | 5979/5984 | 4310 | False | 0.8062642 | 0.8039349 | 0.1822399 | 0.2710401 |
| enhancer | chromosome | chr9 | 13427 | 6716/6711 | 4444 | False | 0.7940114 | 0.7974667 | 0.1848548 | 0.2605809 |
| enhancer | ambiguous_base_status | absent | 25390 | 12695/12695 | 8754 | False | 0.7998119 | 0.8006103 | 0.1836227 | 0.2655091 |
| h3k27me3 | gc_fraction | quintile_1 | 752 | 371/381 | 604 | False | 0.8164671 | 0.800914 | 0.1846038 | 0.2614542 |
| h3k27me3 | gc_fraction | quintile_2 | 648 | 327/321 | 434 | False | 0.819162 | 0.8017091 | 0.1828926 | 0.2683667 |
| h3k27me3 | gc_fraction | quintile_3 | 703 | 348/355 | 413 | False | 0.8126907 | 0.8119637 | 0.1770777 | 0.2916189 |
| h3k27me3 | gc_fraction | quintile_4 | 671 | 332/339 | 350 | False | 0.7788463 | 0.7653979 | 0.2024215 | 0.1902257 |
| h3k27me3 | gc_fraction | quintile_5 | 734 | 376/358 | 258 | False | 0.8346675 | 0.8240892 | 0.1727223 | 0.308695 |
| h3k27me3 | atac_signal_percentile_max_train_only | quintile_1 | 737 | 366/371 | 501 | False | 0.8047031 | 0.790133 | 0.1909224 | 0.2362752 |
| h3k27me3 | atac_signal_percentile_max_train_only | quintile_2 | 717 | 361/356 | 604 | False | 0.8167697 | 0.7921504 | 0.1877979 | 0.2487718 |
| h3k27me3 | atac_signal_percentile_max_train_only | quintile_3 | 629 | 314/315 | 484 | False | 0.8008656 | 0.7884036 | 0.187927 | 0.2482899 |
| h3k27me3 | atac_signal_percentile_max_train_only | quintile_4 | 729 | 365/364 | 519 | False | 0.8313346 | 0.8202092 | 0.1728573 | 0.3085694 |
| h3k27me3 | atac_signal_percentile_max_train_only | quintile_5 | 696 | 348/348 | 485 | False | 0.7992417 | 0.8080328 | 0.1793752 | 0.2824992 |
| h3k27me3 | repeat_fraction_2001 | quintile_1 | 633 | 326/307 | 315 | False | 0.8701015 | 0.8395516 | 0.1614761 | 0.3535133 |
| h3k27me3 | repeat_fraction_2001 | quintile_2 | 767 | 376/391 | 402 | False | 0.8183341 | 0.8364464 | 0.1646536 | 0.3411336 |
| h3k27me3 | repeat_fraction_2001 | quintile_3 | 766 | 376/390 | 472 | False | 0.7957115 | 0.7879774 | 0.1880856 | 0.2474061 |
| h3k27me3 | repeat_fraction_2001 | quintile_4 | 824 | 416/408 | 565 | False | 0.7937252 | 0.778528 | 0.1945783 | 0.2216135 |
| h3k27me3 | repeat_fraction_2001 | quintile_5 | 518 | 260/258 | 415 | False | 0.7515734 | 0.7330203 | 0.2152807 | 0.1388642 |
| h3k27me3 | atac_lobe_signature | lower_left | 1438 | 719/719 | 1100 | False | 0.8135808 | 0.8046526 | 0.1829045 | 0.268382 |
| h3k27me3 | atac_lobe_signature | lower_left;lower_right | 28 | 14/14 | 28 | True | 0.7254209 | 0.7040816 | 0.2302207 | 0.07911705 |
| h3k27me3 | atac_lobe_signature | lower_left;upper_right | 50 | 25/25 | 48 | True | 0.8387395 | 0.8112 | 0.1810629 | 0.2757486 |
| h3k27me3 | atac_lobe_signature | lower_right | 1128 | 564/564 | 767 | False | 0.7934878 | 0.7806606 | 0.191303 | 0.2347882 |
| h3k27me3 | atac_lobe_signature | lower_right;upper_right | 22 | 11/11 | 22 | True | 0.8143521 | 0.785124 | 0.1945631 | 0.2217477 |
| h3k27me3 | atac_lobe_signature | upper_right | 842 | 421/421 | 622 | False | 0.8313316 | 0.8207018 | 0.1732055 | 0.307178 |
| h3k27me3 | same_lobe_mark_support | lower_left | 735 | 735/0 | 485 | False | — | — | — | — |
| h3k27me3 | same_lobe_mark_support | lower_left;lower_right | 8 | 8/0 | 8 | True | — | — | — | — |
| h3k27me3 | same_lobe_mark_support | lower_left;upper_right | 12 | 12/0 | 10 | True | — | — | — | — |
| h3k27me3 | same_lobe_mark_support | lower_right | 571 | 571/0 | 309 | False | — | — | — | — |
| h3k27me3 | same_lobe_mark_support | lower_right;upper_right | 5 | 5/0 | 5 | True | — | — | — | — |
| h3k27me3 | same_lobe_mark_support | none | 1754 | 0/1754 | 1157 | False | — | — | — | — |
| h3k27me3 | same_lobe_mark_support | upper_right | 423 | 423/0 | 249 | False | — | — | — | — |
| h3k27me3 | chromosome | chr8 | 1724 | 862/862 | 869 | False | 0.8156773 | 0.8168709 | 0.1755482 | 0.2978074 |
| h3k27me3 | chromosome | chr9 | 1784 | 892/892 | 940 | False | 0.8043118 | 0.7807001 | 0.1915808 | 0.2336768 |
| h3k27me3 | ambiguous_base_status | absent | 3508 | 1754/1754 | 1809 | False | 0.8096768 | 0.7989697 | 0.1837016 | 0.2651937 |

## Independent validation and historical preservation

| Validation phase | Status | Checks | Failures | Evidence |
| --- | --- | --- | --- | --- |
| inputs | PASS | 244105 | 0 | provenance/inputs_independent_validation.json |
| cache | PASS | 151 | 0 | provenance/cache_independent_validation.json |
| predictions | PASS | 465 | 0 | provenance/predictions_independent_validation.json |
| evaluation | PASS | 383885 | 0 | provenance/evaluation_independent_validation.json |
| preservation | PASS | 1974 | 0 | provenance/preservation_independent_validation.json |

Independent preservation validation is PASS: historical tracked Git object identities and worktree scope are unchanged. Authorized pretraining-1.1 and internal-training-1.0 payloads and reference inputs are also directly rehashed. Benchmark payloads were not opened, parsed or newly hashed; their preservation evidence is Git object identity and worktree status only.

## Environment, compute and execution ledger

### Phase-I extraction

| Runtime field | Recorded value |
| --- | --- |
| Python | 3.13.0 |
| Executable | /home/maheshwarany2/.local/share/uv/python/cpython-3.13.0-linux-x86_64-gnu/bin/python3.13 |
| h5py | 3.16.0 |
| keras | 3.14.1 |
| numpy | 2.5.0 |
| scikit-learn | 1.9.0 |
| scipy | 1.18.0 |
| tensorflow | 2.20.0 |
| GPU/driver record | NVIDIA A100-SXM4-80GB, GPU-c0f9ea40-1229-3a05-cb51-36356549035c, 580.173.02, 81920 MiB  |
| Logical devices | LogicalDevice(name='/device:CPU:0', device_type='CPU'); LogicalDevice(name='/device:GPU:0', device_type='GPU') |
| CUDA_VISIBLE_DEVICES | 0 |
| PYTHONHASHSEED | 104729 |
| SLURM_ARRAY_TASK_ID | — |
| SLURM_CPUS_PER_TASK | 16 |
| SLURM_JOB_GPUS | 1 |
| SLURM_JOB_ID | 31881465 |
| TF_CUDNN_DETERMINISTIC | 1 |
| TF_DETERMINISTIC_OPS | 1 |
| TF_ENABLE_ONEDNN_OPTS | — |
| TF_XLA_FLAGS | --tf_xla_enable_xla_devices=false |

### Checkpoint inference

| Runtime field | Recorded value |
| --- | --- |
| Python | 3.13.0 |
| Executable | /home/maheshwarany2/.local/share/uv/python/cpython-3.13.0-linux-x86_64-gnu/bin/python3.13 |
| h5py | 3.16.0 |
| keras | 3.14.1 |
| numpy | 2.5.0 |
| scikit-learn | 1.9.0 |
| scipy | 1.18.0 |
| tensorflow | 2.20.0 |
| GPU/driver record | NVIDIA A100-SXM4-80GB, GPU-c0f9ea40-1229-3a05-cb51-36356549035c, 580.173.02, 81920 MiB  |
| Logical devices | LogicalDevice(name='/device:CPU:0', device_type='CPU'); LogicalDevice(name='/device:GPU:0', device_type='GPU') |
| CUDA_VISIBLE_DEVICES | 0 |
| PYTHONHASHSEED | 104729 |
| SLURM_ARRAY_TASK_ID | — |
| SLURM_CPUS_PER_TASK | 16 |
| SLURM_JOB_GPUS | 1 |
| SLURM_JOB_ID | 31881465 |
| TF_CUDNN_DETERMINISTIC | 1 |
| TF_DETERMINISTIC_OPS | 1 |
| TF_ENABLE_ONEDNN_OPTS | — |
| TF_XLA_FLAGS | --tf_xla_enable_xla_devices=false |

| Stage operation | Wall seconds | User CPU seconds | System CPU seconds | Peak RSS KiB |
| --- | --- | --- | --- | --- |
| Sequence preparation | 15.59217 | not separately recorded | not separately recorded | 358512 |
| Paired phase-I extraction | 26.03784 | 14.41359 | 2.591425 | 2395568 |
| All frozen checkpoint/RC inference | 71.58277 | 63.29055 | 6.714552 | 3087964 |

The frozen runtime is CPython 3.13.0 / TensorFlow 2.20.0 / Keras 3.14.1 / NumPy 2.5.0, with deterministic settings and suppressed bytecode writes to historical stages. Phase-I batch size is 32; phase-II inference batch size is 256. Compute records describe the actual allocated execution, not a projected retraining cost.

| Command | Return code | Wall seconds | Started UTC | Completed UTC |
| --- | --- | --- | --- | --- |
| collect_execution_provenance | 0 | 12.00202 | 2026-10-06T19:16:29.636481+00:00 | 2026-10-06T19:16:41.642320+00:00 |
| freeze_evaluation_execution_gate | 0 | 13.74709 | 2026-10-06T19:11:42.616617+00:00 | 2026-10-06T19:11:56.375442+00:00 |
| freeze_inference_execution_gate | 0 | 9.107804 | 2026-10-06T19:08:29.510518+00:00 | 2026-10-06T19:08:38.626379+00:00 |
| freeze_phase_I_execution_gate | 0 | 3.413867 | 2026-10-06T19:06:09.730071+00:00 | 2026-10-06T19:06:13.148415+00:00 |
| frozen18_test_inference_and_real_RC | 0 | 72.76768 | 2026-10-06T19:09:32.655828+00:00 | 2026-10-06T19:10:45.428672+00:00 |
| independent_cache_validation | 0 | 2.346657 | 2026-10-06T19:07:36.707279+00:00 | 2026-10-06T19:07:39.063658+00:00 |
| independent_evaluation_validation | 0 | 29.93257 | 2026-10-06T19:13:07.039249+00:00 | 2026-10-06T19:13:36.979253+00:00 |
| independent_historical_preservation | 0 | 33.85901 | 2026-10-06T19:15:57.644191+00:00 | 2026-10-06T19:16:31.507282+00:00 |
| independent_input_validation | 0 | 11.46415 | 2026-10-06T19:05:51.885301+00:00 | 2026-10-06T19:06:03.356100+00:00 |
| independent_prediction_validation | 0 | 6.599327 | 2026-10-06T19:11:04.121517+00:00 | 2026-10-06T19:11:10.725902+00:00 |
| one_time_common_test_evaluation | 0 | 27.69564 | 2026-10-06T19:12:01.276686+00:00 | 2026-10-06T19:12:28.979517+00:00 |
| phase_I_test_forward_and_RC | 0 | 28.86425 | 2026-10-06T19:06:39.571828+00:00 | 2026-10-06T19:07:08.444510+00:00 |
| prepare_frozen_test_panels | 0 | 16.02668 | 2026-10-06T19:02:21.007358+00:00 | 2026-10-06T19:02:37.038892+00:00 |
| synthetic_final_preinference_tests | 0 | 2.796495 | 2026-10-06T19:05:53.043508+00:00 | 2026-10-06T19:05:55.844836+00:00 |
| synthetic_independent_validator_tests | 0 | 2.711104 | 2026-10-06T19:07:37.880944+00:00 | 2026-10-06T19:07:40.598019+00:00 |
| synthetic_preexecution_tests | 0 | 2.848 | 2026-10-06T19:02:53.738690+00:00 | 2026-10-06T19:02:56.592146+00:00 |
| verify_immutable_test_sources | 0 | 32.06498 | 2026-10-06T19:00:21.905331+00:00 | 2026-10-06T19:00:53.975928+00:00 |

Exact commands, environment settings, stdout/stderr and return codes are preserved under [provenance/commands](provenance/commands/). Nonzero command returns are retained below rather than suppressed. Evaluator return code 3 denotes the saved FAIL/INCONCLUSIVE scientific stop decision, not permission for an outcome-driven rerun.

No records.

## Artifact index and final stop

| Artifact | Stage-relative path |
| --- | --- |
| metrics.tsv | [metrics.tsv](results/chr8_9_evaluation/metrics.tsv) |
| bootstrap_summary.tsv | [bootstrap_summary.tsv](results/chr8_9_evaluation/bootstrap_summary.tsv) |
| paired_ablation.tsv | [paired_ablation.tsv](results/chr8_9_evaluation/paired_ablation.tsv) |
| seed_stability.tsv | [seed_stability.tsv](results/chr8_9_evaluation/seed_stability.tsv) |
| C_holdout_generalization.tsv | [C_holdout_generalization.tsv](results/chr8_9_evaluation/C_holdout_generalization.tsv) |
| holdout_decision.json | [holdout_decision.json](results/chr8_9_evaluation/holdout_decision.json) |
| fixed_threshold_operating.tsv | [fixed_threshold_operating.tsv](results/chr8_9_evaluation/fixed_threshold_operating.tsv) |
| fixed_threshold_intervals.tsv | [fixed_threshold_intervals.tsv](results/chr8_9_evaluation/fixed_threshold_intervals.tsv) |
| reliability_bins.tsv | [reliability_bins.tsv](results/chr8_9_evaluation/reliability_bins.tsv) |
| stratified_metrics.tsv | [stratified_metrics.tsv](results/chr8_9_evaluation/stratified_metrics.tsv) |
| bootstrap_audit.json | [bootstrap_audit.json](results/chr8_9_evaluation/bootstrap_audit.json) |
| evaluation_manifest.json | [evaluation_manifest.json](results/chr8_9_evaluation/evaluation_manifest.json) |

Complete per-seed/ensemble scores, repeated independent RC outputs, 2,000-replicate bootstrap records and score-linked components are preserved in the new-stage prediction/result tables. This report reads summary artifacts only and does not recompute model scores, bootstrap statistics or threshold calibration. The final checksum/freeze manifest records this report and all completed payloads; shared V2 activity, decision and result registers receive append-only stage-completion entries.

READY FOR INVESTIGATOR REVIEW BEFORE ANY EXTERNAL FUNCTIONAL EVALUATION. Both C contexts pass the prospective holdout rule and all five recorded independent validation phases pass. This is a review disposition only: separate authorization is still required for external evaluation or candidate/allele scoring.

**Stop. Do not open the external functional benchmark or begin candidate/allele scoring without another explicit authorization.**
