# COPD V2 internally controlled training and chr7 evaluation

Version: `internal-training-1.0`; immutable design: `pretraining-1.1`.

## Decision and stop boundary

READY FOR INVESTIGATOR REVIEW BEFORE OPENING CHR8–9. Both C contexts, conditional calibration and the recorded independent QC passed. Separate test/external authorization remains required.

| V2-C context | Saved absolute-adequacy decision | Failed/inconclusive conditions |
| --- | --- | --- |
| Enhancer-associated | PASS | none |
| H3K27me3-associated | PASS | none |

Only region-label modeling was evaluated. A and B are diagnostic ablations, never eligible fallbacks or competitors for final-model selection. Neither chr8–9 performance nor the external COPD-V2-BENCH benchmark was opened in this stage. The COPD variant universe and 337 historical candidates were not scored. No new candidate list, fine-mapping, target-gene work or manuscript revision was performed.

## Frozen design and execution

A preserves historical labels and controls while correcting orientation handling. B replaces historical controls with frozen donor-accessible matched controls. C uses same-lobe positives and independently rematched controls. Consequently B→C changes both label membership and control matching, not labels alone. All comparisons below use the identical frozen C-task chr7 selection panel within each context; native-task metrics across different class definitions are not used to choose a model.

All 18 fits use seeds `104729`, `130363`, `155921`, the unchanged V1 phase-II architecture, unweighted binary cross-entropy, Adadelta (learning rate 0.001, rho 0.95, epsilon 1e-7), batch size 256, maximum 50 epochs and patience 15. One nucleotide orientation per biological interval per epoch is sampled with probability 0.5. Frozen seed substreams assign orientations in genomic order before the independent epoch permutation. The earliest epoch with minimum native chr7 checkpoint-role symmetric BCE is retained. Selection/calibration rows cannot choose a checkpoint. The selected checkpoint contains all selected-epoch trainable and BatchNorm weights and is inference-only, not an optimizer-resumption checkpoint.

Checkpoint-packaging qualification: Keras preserved `compile_config` during JSON reconstruction, so the `.keras` containers are not literally uncompiled or optimizer-free. They include fresh ancillary Adadelta state, not the selected epoch's training-optimizer state. Actual prediction loads use `compile=False`, which ignores that ancillary optimizer state. Selected network weights were restored and checked by bitwise serialization round-trip. These artifacts are frozen inference checkpoints only and must not be used to resume training. Any original implementation comment describing the saved container as uncompiled is superseded by the append-only packaging clarification; the checkpoints and scientific procedure were not rewritten or rerun.

The final symmetric score is the float64 equal-weight mean, in the frozen seed order, of each seed's `(p_forward+p_nucleotide_RC)/2`. No representation-axis reversal, seed selection, seed weighting, mixed precision or hyperparameter search is permitted.

## Phase-I cache and sequence integrity

| Record | Value |
| --- | --- |
| Unique genomic intervals used | 719156 |
| Encoded sequence/RC identities | 719147 |
| Independent nucleotide-orientation network evaluations | 1438294 |
| Per-orientation cache geometry | 719147 × 4560 |
| Cache dtype / representation | float32 / final4560sigmoid |
| Missing sequences / extraction failures | 0 / 0 |
| chr8–9 sequence extractions | 0 |
| Phase-I weight SHA-256 | 483b6c0cafd750ea8e97932cbc8a949276eaf1e77e7d8766f9a22eeddc2916d6 |

Both canonical nucleotide sequence and its nucleotide reverse complement independently passed through the same frozen, inference-mode phase-I network. Each interval's forward-orientation flag maps the two caches back to genomic forward/RC. Ambiguous bases retain the frozen all-zero encoding after normalization; no new interval exclusion was introduced.

| Cache artifact | Bytes | SHA-256 |
| --- | --- | --- |
| cache/features_canonical.npy | 13117241408 | 8d7bcb61888150e9d44ba9eb1204a2825e421827500628ed6a2c69a560f1cd0a |
| cache/features_rc.npy | 13117241408 | 4fab1cc6c84eab0f27d8ab2f70d23081617966afb8282c94cbb377dbb536fa33 |

Chromosome interval counts (training chromosomes plus chr7 only):

| Chromosome | Intervals |
| --- | --- |
| chr1 | 68625 |
| chr10 | 41500 |
| chr11 | 38007 |
| chr12 | 36953 |
| chr13 | 23638 |
| chr14 | 25134 |
| chr15 | 24713 |
| chr16 | 25302 |
| chr17 | 29914 |
| chr18 | 20775 |
| chr19 | 18511 |
| chr2 | 66380 |
| chr20 | 21234 |
| chr21 | 10127 |
| chr22 | 13774 |
| chr3 | 54207 |
| chr4 | 41111 |
| chr5 | 45464 |
| chr6 | 45597 |
| chr7 | 41485 |
| chrX | 25304 |
| chrY | 1401 |

## All 18 fit outcomes

| Configuration | Context | Seed | Outcome | Epochs | Selected epoch | Native checkpoint BCE | Attempt | Retries |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| V2-A | enhancer | 104729 | COMPLETED | 34 | 19 | 0.2901492 | 1 | 0 |
| V2-A | enhancer | 130363 | COMPLETED | 46 | 31 | 0.2861099 | 1 | 0 |
| V2-A | enhancer | 155921 | COMPLETED | 43 | 28 | 0.2862209 | 1 | 0 |
| V2-A | h3k27me3 | 104729 | COMPLETED | 42 | 27 | 0.1779631 | 1 | 0 |
| V2-A | h3k27me3 | 130363 | COMPLETED | 45 | 30 | 0.1764181 | 1 | 0 |
| V2-A | h3k27me3 | 155921 | COMPLETED | 43 | 28 | 0.1782582 | 1 | 0 |
| V2-B | enhancer | 104729 | COMPLETED | 21 | 6 | 0.5809753 | 1 | 0 |
| V2-B | enhancer | 130363 | COMPLETED | 23 | 8 | 0.5766301 | 1 | 0 |
| V2-B | enhancer | 155921 | COMPLETED | 22 | 7 | 0.5792029 | 1 | 0 |
| V2-B | h3k27me3 | 104729 | COMPLETED | 41 | 26 | 0.5222496 | 1 | 0 |
| V2-B | h3k27me3 | 130363 | COMPLETED | 34 | 19 | 0.5205834 | 1 | 0 |
| V2-B | h3k27me3 | 155921 | COMPLETED | 27 | 12 | 0.5220316 | 1 | 0 |
| V2-C | enhancer | 104729 | COMPLETED | 25 | 10 | 0.5572999 | 1 | 0 |
| V2-C | enhancer | 130363 | COMPLETED | 23 | 8 | 0.5540302 | 1 | 0 |
| V2-C | enhancer | 155921 | COMPLETED | 22 | 7 | 0.5552501 | 1 | 0 |
| V2-C | h3k27me3 | 104729 | COMPLETED | 29 | 14 | 0.4969136 | 1 | 0 |
| V2-C | h3k27me3 | 130363 | COMPLETED | 39 | 24 | 0.4916464 | 1 | 0 |
| V2-C | h3k27me3 | 155921 | COMPLETED | 23 | 8 | 0.4989785 | 1 | 0 |

Complete histories, epoch RNG hashes, checkpoint-improvement ledgers, selected-checkpoint records and serialized-weight round-trip checks remain under [runs](runs/). The all-18 checkpoint release and exact hashes are in [checkpoint_freeze.json](provenance/checkpoint_freeze.json).

### Selected-checkpoint container packaging

Independent final archive audit: **PASS**; 18 checkpoint containers recorded. [checkpoint_packaging_final_validation.json](provenance/checkpoint_packaging_final_validation.json) and its linked tensor-level ledger document the original, unchanged containers.

| Run | Packaging checks | Model tensors | Optimizer tensors | Optimizer slots | Iteration | All slots zero |
| --- | --- | --- | --- | --- | --- | --- |
| V2-A_enhancer_seed104729 | PASS | 14 | 26 | 24 | 0 | yes |
| V2-A_enhancer_seed130363 | PASS | 14 | 26 | 24 | 0 | yes |
| V2-A_enhancer_seed155921 | PASS | 14 | 26 | 24 | 0 | yes |
| V2-A_h3k27me3_seed104729 | PASS | 14 | 26 | 24 | 0 | yes |
| V2-A_h3k27me3_seed130363 | PASS | 14 | 26 | 24 | 0 | yes |
| V2-A_h3k27me3_seed155921 | PASS | 14 | 26 | 24 | 0 | yes |
| V2-B_enhancer_seed104729 | PASS | 14 | 26 | 24 | 0 | yes |
| V2-B_enhancer_seed130363 | PASS | 14 | 26 | 24 | 0 | yes |
| V2-B_enhancer_seed155921 | PASS | 14 | 26 | 24 | 0 | yes |
| V2-B_h3k27me3_seed104729 | PASS | 14 | 26 | 24 | 0 | yes |
| V2-B_h3k27me3_seed130363 | PASS | 14 | 26 | 24 | 0 | yes |
| V2-B_h3k27me3_seed155921 | PASS | 14 | 26 | 24 | 0 | yes |
| V2-C_enhancer_seed104729 | PASS | 14 | 26 | 24 | 0 | yes |
| V2-C_enhancer_seed130363 | PASS | 14 | 26 | 24 | 0 | yes |
| V2-C_enhancer_seed155921 | PASS | 14 | 26 | 24 | 0 | yes |
| V2-C_h3k27me3_seed104729 | PASS | 14 | 26 | 24 | 0 | yes |
| V2-C_h3k27me3_seed130363 | PASS | 14 | 26 | 24 | 0 | yes |
| V2-C_h3k27me3_seed155921 | PASS | 14 | 26 | 24 | 0 | yes |


## Actual-network reverse-complement invariance

Saved gate: **PASS**. Actual-network evidence: yes; 18 seed-level audits and 6 ensemble audits. Every frozen internal chr7 validation sequence in the context-specific native/common union is covered. The second wrapper pass calls the trained network again on swapped nucleotide-orientation inputs, using the same frozen row order and batch size 256; original probabilities are not reused as the purported second network pass. The frozen comparison is `abs(Q(x)-Q(RC(x))) <= 1e-6 + 1e-6*abs(Q(RC(x)))`, with finite [0,1] checks on all probabilities.

| Configuration | Context | Unit | Seed | Sequences | Status | Failed | Maximum absolute difference |
| --- | --- | --- | --- | --- | --- | --- | --- |
| V2-A | enhancer | seed | 104729 | 37494 | PASS | 0 | 0 |
| V2-A | enhancer | seed | 130363 | 37494 | PASS | 0 | 0 |
| V2-A | enhancer | seed | 155921 | 37494 | PASS | 0 | 0 |
| V2-A | enhancer | ensemble | all3 | 37494 | PASS | 0 | 0 |
| V2-B | enhancer | seed | 104729 | 37494 | PASS | 0 | 0 |
| V2-B | enhancer | seed | 130363 | 37494 | PASS | 0 | 0 |
| V2-B | enhancer | seed | 155921 | 37494 | PASS | 0 | 0 |
| V2-B | enhancer | ensemble | all3 | 37494 | PASS | 0 | 0 |
| V2-C | enhancer | seed | 104729 | 37494 | PASS | 0 | 0 |
| V2-C | enhancer | seed | 130363 | 37494 | PASS | 0 | 0 |
| V2-C | enhancer | seed | 155921 | 37494 | PASS | 0 | 0 |
| V2-C | enhancer | ensemble | all3 | 37494 | PASS | 0 | 0 |
| V2-A | h3k27me3 | seed | 104729 | 6692 | PASS | 0 | 0 |
| V2-A | h3k27me3 | seed | 130363 | 6692 | PASS | 0 | 0 |
| V2-A | h3k27me3 | seed | 155921 | 6692 | PASS | 0 | 0 |
| V2-A | h3k27me3 | ensemble | all3 | 6692 | PASS | 0 | 0 |
| V2-B | h3k27me3 | seed | 104729 | 6692 | PASS | 0 | 0 |
| V2-B | h3k27me3 | seed | 130363 | 6692 | PASS | 0 | 0 |
| V2-B | h3k27me3 | seed | 155921 | 6692 | PASS | 0 | 0 |
| V2-B | h3k27me3 | ensemble | all3 | 6692 | PASS | 0 | 0 |
| V2-C | h3k27me3 | seed | 104729 | 6692 | PASS | 0 | 0 |
| V2-C | h3k27me3 | seed | 130363 | 6692 | PASS | 0 | 0 |
| V2-C | h3k27me3 | seed | 155921 | 6692 | PASS | 0 | 0 |
| V2-C | h3k27me3 | ensemble | all3 | 6692 | PASS | 0 | 0 |

Detailed actual per-orientation and independent reverse-wrapper outputs are in [predictions](predictions/); [real_network_invariance.json](predictions/real_network_invariance.json) binds the checkpoints, inputs and implementation.

## Common C-task chr7 selection performance

The following are descriptive comparisons on the same frozen panel per context. AP is stepwise average precision; Brier skill uses the constant-prevalence baseline. Only C's predeclared absolute gates determine adequacy. Rounded display values below do not replace the full-precision TSV/JSON records.

| Context | Configuration | Predictor | Positive / control | Components | AP | AUROC | Brier | Brier skill |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| enhancer | V2-A | seed:104729 | 1755 / 1755 | 1295 | 0.7173839 | 0.6937791 | 0.2599604 | -0.03984156 |
| enhancer | V2-A | seed:130363 | 1755 / 1755 | 1295 | 0.714261 | 0.6963801 | 0.2568344 | -0.02733746 |
| enhancer | V2-A | seed:155921 | 1755 / 1755 | 1295 | 0.7235207 | 0.7008105 | 0.2557087 | -0.02283494 |
| enhancer | V2-A | ensemble | 1755 / 1755 | 1295 | 0.7195443 | 0.697293 | 0.257118 | -0.02847196 |
| enhancer | V2-B | seed:104729 | 1755 / 1755 | 1295 | 0.8130222 | 0.8042061 | 0.1846268 | 0.2614927 |
| enhancer | V2-B | seed:130363 | 1755 / 1755 | 1295 | 0.8146214 | 0.8056681 | 0.1828553 | 0.268579 |
| enhancer | V2-B | seed:155921 | 1755 / 1755 | 1295 | 0.8124065 | 0.8035876 | 0.1841034 | 0.2635862 |
| enhancer | V2-B | ensemble | 1755 / 1755 | 1295 | 0.8138468 | 0.8049892 | 0.1836989 | 0.2652046 |
| enhancer | V2-C | seed:104729 | 1755 / 1755 | 1295 | 0.8165024 | 0.8077662 | 0.1820957 | 0.2716171 |
| enhancer | V2-C | seed:130363 | 1755 / 1755 | 1295 | 0.8145257 | 0.8063292 | 0.180956 | 0.276176 |
| enhancer | V2-C | seed:155921 | 1755 / 1755 | 1295 | 0.8118307 | 0.8026529 | 0.1822093 | 0.271163 |
| enhancer | V2-C | ensemble | 1755 / 1755 | 1295 | 0.8147664 | 0.8061256 | 0.1814778 | 0.2740888 |
| h3k27me3 | V2-A | seed:104729 | 271 / 271 | 244 | 0.8519439 | 0.8207949 | 0.3076136 | -0.2304544 |
| h3k27me3 | V2-A | seed:130363 | 271 / 271 | 244 | 0.8451828 | 0.8153348 | 0.3116545 | -0.2466178 |
| h3k27me3 | V2-A | seed:155921 | 271 / 271 | 244 | 0.859516 | 0.8266363 | 0.297926 | -0.191704 |
| h3k27me3 | V2-A | ensemble | 271 / 271 | 244 | 0.8548062 | 0.8225378 | 0.3055108 | -0.2220434 |
| h3k27me3 | V2-B | seed:104729 | 271 / 271 | 244 | 0.8900858 | 0.8793589 | 0.145131 | 0.4194761 |
| h3k27me3 | V2-B | seed:130363 | 271 / 271 | 244 | 0.8865024 | 0.8790049 | 0.1455128 | 0.4179489 |
| h3k27me3 | V2-B | seed:155921 | 271 / 271 | 244 | 0.8892949 | 0.8797947 | 0.1453605 | 0.4185578 |
| h3k27me3 | V2-B | ensemble | 271 / 271 | 244 | 0.8893853 | 0.8801215 | 0.1450633 | 0.4197469 |
| h3k27me3 | V2-C | seed:104729 | 271 / 271 | 244 | 0.8906963 | 0.880598 | 0.1408204 | 0.4367183 |
| h3k27me3 | V2-C | seed:130363 | 271 / 271 | 244 | 0.8864167 | 0.8766357 | 0.1430001 | 0.4279996 |
| h3k27me3 | V2-C | seed:155921 | 271 / 271 | 244 | 0.8930045 | 0.8809112 | 0.1400247 | 0.4399012 |
| h3k27me3 | V2-C | ensemble | 271 / 271 | 244 | 0.8904478 | 0.8800261 | 0.1409546 | 0.4361816 |

### Fixed-seed stability

| Context | Configuration | AP 104729 | AP 130363 | AP 155921 | AP sample SD | AP range |
| --- | --- | --- | --- | --- | --- | --- |
| enhancer | V2-A | 0.7173839 | 0.714261 | 0.7235207 | 0.004710876 | 0.009259684 |
| enhancer | V2-B | 0.8130222 | 0.8146214 | 0.8124065 | 0.001143261 | 0.0022149 |
| enhancer | V2-C | 0.8165024 | 0.8145257 | 0.8118307 | 0.002345024 | 0.004671672 |
| h3k27me3 | V2-A | 0.8519439 | 0.8451828 | 0.859516 | 0.007170428 | 0.01433321 |
| h3k27me3 | V2-B | 0.8900858 | 0.8865024 | 0.8892949 | 0.001882597 | 0.003583466 |
| h3k27me3 | V2-C | 0.8906963 | 0.8864167 | 0.8930045 | 0.003342712 | 0.006587828 |

C requires sample SD ≤0.03 (`ddof=1`) and range ≤0.10 for each context. These limits were not revised after outcomes. A/B seed behavior is descriptive only.

### Paired component-bootstrap uncertainty

| Context | Components | Attempts | Valid | Invalid | Invalid fraction | Status |
| --- | --- | --- | --- | --- | --- | --- |
| enhancer | 1295 | 2000 | 2000 | 0 | 0 | PASS |
| h3k27me3 | 244 | 2000 | 2000 | 0 | 0 | PASS |

The frozen bootstrap uses seed 314159 and 2,000 valid genomic overlap/encoded-identity component draws, at most 20,000 attempts and invalid fraction <0.10. A/B/C share each draw and component multiplicities. The uncertainty is conditional on this donor and fixed matching design; it is not independent-donor or population generalization. AP minus prevalence and Brier skill are computed within each replicate before confidence limits, not by substituting the original prevalence afterward.

| Context | Configuration | Metric | Point [95% component-bootstrap CI] |
| --- | --- | --- | --- |
| enhancer | V2-A | AP | 0.7195443 [0.6686641, 0.7658479] |
| enhancer | V2-A | AUROC | 0.697293 [0.6647352, 0.7298917] |
| enhancer | V2-A | Brier | 0.257118 [0.2392257, 0.2757007] |
| enhancer | V2-A | BrierSkill | -0.02847196 [-0.1064087, 0.04128843] |
| enhancer | V2-A | AP_gain | 0.2195443 [0.1842824, 0.2553022] |
| enhancer | V2-B | AP | 0.8138468 [0.7735628, 0.8496027] |
| enhancer | V2-B | AUROC | 0.8049892 [0.7804997, 0.8289807] |
| enhancer | V2-B | Brier | 0.1836989 [0.173993, 0.1930351] |
| enhancer | V2-B | BrierSkill | 0.2652046 [0.2256113, 0.3026127] |
| enhancer | V2-B | AP_gain | 0.3138468 [0.2853406, 0.341062] |
| enhancer | V2-C | AP | 0.8147664 [0.7756742, 0.8500979] |
| enhancer | V2-C | AUROC | 0.8061256 [0.7819501, 0.8298224] |
| enhancer | V2-C | Brier | 0.1814778 [0.1713644, 0.1910657] |
| enhancer | V2-C | BrierSkill | 0.2740888 [0.2335666, 0.3127776] |
| enhancer | V2-C | AP_gain | 0.3147664 [0.2868651, 0.3412591] |
| h3k27me3 | V2-A | AP | 0.8548062 [0.7665115, 0.9173454] |
| h3k27me3 | V2-A | AUROC | 0.8225378 [0.7587103, 0.8775839] |
| h3k27me3 | V2-A | Brier | 0.3055108 [0.2469508, 0.3685997] |
| h3k27me3 | V2-A | BrierSkill | -0.2220434 [-0.5122831, -0.004616087] |
| h3k27me3 | V2-A | AP_gain | 0.3548062 [0.2936653, 0.4178093] |
| h3k27me3 | V2-B | AP | 0.8893853 [0.8217259, 0.9355439] |
| h3k27me3 | V2-B | AUROC | 0.8801215 [0.8326164, 0.9193523] |
| h3k27me3 | V2-B | Brier | 0.1450633 [0.1237513, 0.1676029] |
| h3k27me3 | V2-B | BrierSkill | 0.4197469 [0.3226465, 0.4995027] |
| h3k27me3 | V2-B | AP_gain | 0.3893853 [0.3207943, 0.4575918] |
| h3k27me3 | V2-C | AP | 0.8904478 [0.8214401, 0.9376147] |
| h3k27me3 | V2-C | AUROC | 0.8800261 [0.8314792, 0.921857] |
| h3k27me3 | V2-C | Brier | 0.1409546 [0.1181795, 0.1656096] |
| h3k27me3 | V2-C | BrierSkill | 0.4361816 [0.3324161, 0.5239712] |
| h3k27me3 | V2-C | AP_gain | 0.3904478 [0.3234105, 0.4580727] |

### Diagnostic A→B, B→C and A→C ablations

| Context | Difference | Metric | Point [95% paired CI] |
| --- | --- | --- | --- |
| enhancer | V2-B minus V2-A | AP | 0.09430253 [0.07330793, 0.1183187] |
| enhancer | V2-B minus V2-A | AUROC | 0.1076962 [0.0870592, 0.1297424] |
| enhancer | V2-B minus V2-A | Brier | -0.07341913 [-0.08664281, -0.06118751] |
| enhancer | V2-C minus V2-B | AP | 0.0009195484 [-0.001008576, 0.00291683] |
| enhancer | V2-C minus V2-B | AUROC | 0.001136354 [-0.0007484373, 0.003091818] |
| enhancer | V2-C minus V2-B | Brier | -0.002221062 [-0.003139031, -0.001282265] |
| enhancer | V2-C minus V2-A | AP | 0.09522208 [0.07344465, 0.1198952] |
| enhancer | V2-C minus V2-A | AUROC | 0.1088326 [0.08742173, 0.1307639] |
| enhancer | V2-C minus V2-A | Brier | -0.07564019 [-0.08907435, -0.06293641] |
| h3k27me3 | V2-B minus V2-A | AP | 0.03457912 [-0.008866655, 0.08480668] |
| h3k27me3 | V2-B minus V2-A | AUROC | 0.05758364 [0.00940841, 0.1039464] |
| h3k27me3 | V2-B minus V2-A | Brier | -0.1604476 [-0.2203053, -0.09988904] |
| h3k27me3 | V2-C minus V2-B | AP | 0.001062498 [-0.004448166, 0.006854882] |
| h3k27me3 | V2-C minus V2-B | AUROC | -9.531461e-05 [-0.006442071, 0.006118156] |
| h3k27me3 | V2-C minus V2-B | Brier | -0.004108668 [-0.008256547, 0.0002582189] |
| h3k27me3 | V2-C minus V2-A | AP | 0.03564162 [-0.004981273, 0.08262111] |
| h3k27me3 | V2-C minus V2-A | AUROC | 0.05748832 [0.01063244, 0.1045486] |
| h3k27me3 | V2-C minus V2-A | Brier | -0.1645562 [-0.2274295, -0.102435] |

Negative Brier differences mean lower error for the first configuration. Comparative superiority, equivalence or noninferiority is not a retention gate. A favorable A/B result cannot rescue a failed or inconclusive C context.

### Reliability and fixed strata

[reliability_bins.tsv](results/chr7_evaluation/reliability_bins.tsv) reports the ten prespecified bins `[0,.1), …, [.9,1]` for every seed and ensemble. [stratified_metrics.tsv](results/chr7_evaluation/stratified_metrics.tsv) reports GC, train-derived ATAC rank and repeat quintiles, ATAC lobe signature, same-lobe mark support, chromosome and ambiguous-base status. Quintile cutpoints come only from retained C training positives and are fixed across configurations. All strata are descriptive; cells with fewer than 100 rows are explicitly flagged. No fitted recalibration, subgroup selection or threshold tuning was performed from these displays.

## V2-C absolute-adequacy gates

All conditions are conjunctive within each C context, and both contexts must pass. Boundary equality cannot pass the strict lower-confidence-bound gates. An inconclusive support/bootstrap result blocks release just as a failed gate does.

| Context | Gate | Status | Saved observation | Frozen criterion |
| --- | --- | --- | --- | --- |
| enhancer | frozen_construction | PASS | PASS | B and C construction gates PASS |
| enhancer | all_three_C_seeds_completed | PASS | 3 | 3 fixed seeds; all checkpoints frozen |
| enhancer | real_network_invariance | PASS | PASS | actual-network all seeds and ensembles PASS |
| enhancer | finite_valid_probabilities | PASS | PASS | all orientation/seed/ensemble probabilities finite in [0,1] |
| enhancer | selection_positives | PASS | 1755 | >=100 |
| enhancer | selection_controls | PASS | 1755 | >=200 |
| enhancer | selection_components | PASS | 1295 | >=30 |
| enhancer | bootstrap_valid_replicates | PASS | 2000 | exactly2000valid; at most20000attempts |
| enhancer | bootstrap_invalid_fraction | PASS | 0.0 | strictly<0.10 |
| enhancer | lower95_AUROC | PASS | 0.781950099008121 | strictly>0.5 |
| enhancer | lower95_AP_gain | PASS | 0.28686507905101954 | strictly>0.0 |
| enhancer | lower95_BrierSkill | PASS | 0.23356659022794932 | strictly>0.0 |
| enhancer | seed_AP_sample_SD | PASS | 0.0023450235782206366 | <=0.03; ddof=1 |
| enhancer | seed_AP_range | PASS | 0.004671672459874321 | <=0.10 |
| h3k27me3 | frozen_construction | PASS | PASS | B and C construction gates PASS |
| h3k27me3 | all_three_C_seeds_completed | PASS | 3 | 3 fixed seeds; all checkpoints frozen |
| h3k27me3 | real_network_invariance | PASS | PASS | actual-network all seeds and ensembles PASS |
| h3k27me3 | finite_valid_probabilities | PASS | PASS | all orientation/seed/ensemble probabilities finite in [0,1] |
| h3k27me3 | selection_positives | PASS | 271 | >=100 |
| h3k27me3 | selection_controls | PASS | 271 | >=200 |
| h3k27me3 | selection_components | PASS | 244 | >=30 |
| h3k27me3 | bootstrap_valid_replicates | PASS | 2000 | exactly2000valid; at most20000attempts |
| h3k27me3 | bootstrap_invalid_fraction | PASS | 0.0 | strictly<0.10 |
| h3k27me3 | lower95_AUROC | PASS | 0.8314791998473509 | strictly>0.5 |
| h3k27me3 | lower95_AP_gain | PASS | 0.32341049138037337 | strictly>0.0 |
| h3k27me3 | lower95_BrierSkill | PASS | 0.33241614514708556 | strictly>0.0 |
| h3k27me3 | seed_AP_sample_SD | PASS | 0.0033427115344803146 | <=0.03; ddof=1 |
| h3k27me3 | seed_AP_range | PASS | 0.006587828199550949 | <=0.10 |

Saved next action: `FREEZE_C_CHECKPOINTS_ENSEMBLE_AND_ADEQUACY_THEN_CALIBRATE_CHR7_NEGATIVES_ONLY`.

## Conditional chr7 negative-control calibration

Saved calibration status: **PASS**. Calibration used only dedicated common chr7 calibration negative controls, after both C contexts and the six C checkpoints/ensemble were frozen. No positive score chose a threshold. With `k=floor(0.05*n)` and descending negative scores `d`, the rule is `t=nextafter(d[k], +infinity)` and calls use `score>=t`.

| Context | Controls | Control components | Boundary (17 digits) | Threshold (17 digits) | Threshold (hex) | Called / controls | Observed row FPR | Boundary ties |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| enhancer | 1576 | 809 | 0.74848511815071106 | 0.74848511815071117 | 0x1.7f39710000001p-1 | 78 / 1576 | 0.049492385786802033 | 1 |
| h3k27me3 | 263 | 176 | 0.76960810025533044 | 0.76960810025533055 | 0x1.8a0a12aaaaaacp-1 | 13 / 263 | 0.049429657794676805 | 1 |

All observations tied at the boundary are excluded; ties are never split. Exact decimal/hex serialization, empirical-row FPR checks and any threshold-above-one edge case are recorded in [C_region_thresholds.json](results/chr7_calibration/C_region_thresholds.json). This bounds the observed calibration-row fraction only; it does not establish a population-FPR guarantee, variant FPR, allele-effect FDR or causal probability. No allele-delta threshold was defined.

## Software, command provenance and compute

[provenance/commands](provenance/commands/) retains exact commands, launcher hashes, environment fields, streamed logs, return codes and elapsed times. Per-fit/cache/inference environment records preserve Python/package versions, deterministic settings, GPU identity and TensorFlow build metadata. The following are recorded process elapsed times, including I/O, hashing and serialization—not measured active GPU-kernel time or full scheduler billing.

| Process | Recorded wall seconds | Peak RSS (KiB) |
| --- | --- | --- |
| Phase-I cache extraction | 393.5700578 | 28448564 |
| V2-A_enhancer_seed104729 | 1943.969866 | 27342192 |
| V2-A_enhancer_seed130363 | 2609.422924 | 27422424 |
| V2-A_enhancer_seed155921 | 2629.701514 | 27345980 |
| V2-A_h3k27me3_seed104729 | 481.4463981 | 11427812 |
| V2-A_h3k27me3_seed130363 | 515.3140678 | 11496304 |
| V2-A_h3k27me3_seed155921 | 481.8152967 | 11621096 |
| V2-B_enhancer_seed104729 | 933.7876792 | 21607060 |
| V2-B_enhancer_seed130363 | 999.565644 | 21563272 |
| V2-B_enhancer_seed155921 | 1067.325572 | 21395392 |
| V2-B_h3k27me3_seed104729 | 368.8455849 | 7456832 |
| V2-B_h3k27me3_seed130363 | 293.0640669 | 7638236 |
| V2-B_h3k27me3_seed155921 | 242.3049678 | 7634052 |
| V2-C_enhancer_seed104729 | 809.1107897 | 20288360 |
| V2-C_enhancer_seed130363 | 809.3173237 | 20021268 |
| V2-C_enhancer_seed155921 | 688.8147052 | 20268948 |
| V2-C_h3k27me3_seed104729 | 194.4108558 | 6264048 |
| V2-C_h3k27me3_seed130363 | 192.2154948 | 6333432 |
| V2-C_h3k27me3_seed155921 | 165.1898356 | 6280372 |
| All-context chr7 inference and RC audit | 110.8607309 | 3656372 |

Sum of recorded successful one-GPU process walltimes: 4.425014826 GPU-process-hours. This accounting label does not imply 100% GPU utilization; failed-attempt allocations, queue wait and scheduler overhead are not included.

| Process | Python | TensorFlow | Keras | NumPy | Slurm job |
| --- | --- | --- | --- | --- | --- |
| phase_I | 3.13.0 | 2.20.0 | 3.14.1 | 2.5.0 | 31881465 |
| V2-A_enhancer_seed104729 | 3.13.0 | 2.20.0 | 3.14.1 | 2.5.0 | 32021390 |
| V2-A_enhancer_seed130363 | 3.13.0 | 2.20.0 | 3.14.1 | 2.5.0 | 32021391 |
| V2-A_enhancer_seed155921 | 3.13.0 | 2.20.0 | 3.14.1 | 2.5.0 | 32021392 |
| V2-A_h3k27me3_seed104729 | 3.13.0 | 2.20.0 | 3.14.1 | 2.5.0 | 32021393 |
| V2-A_h3k27me3_seed130363 | 3.13.0 | 2.20.0 | 3.14.1 | 2.5.0 | 32022534 |
| V2-A_h3k27me3_seed155921 | 3.13.0 | 2.20.0 | 3.14.1 | 2.5.0 | 32023281 |
| V2-B_enhancer_seed104729 | 3.13.0 | 2.20.0 | 3.14.1 | 2.5.0 | 32024685 |
| V2-B_enhancer_seed130363 | 3.13.0 | 2.20.0 | 3.14.1 | 2.5.0 | 32026106 |
| V2-B_enhancer_seed155921 | 3.13.0 | 2.20.0 | 3.14.1 | 2.5.0 | 32028256 |
| V2-B_h3k27me3_seed104729 | 3.13.0 | 2.20.0 | 3.14.1 | 2.5.0 | 32028257 |
| V2-B_h3k27me3_seed130363 | 3.13.0 | 2.20.0 | 3.14.1 | 2.5.0 | 32028258 |
| V2-B_h3k27me3_seed155921 | 3.13.0 | 2.20.0 | 3.14.1 | 2.5.0 | 32029656 |
| V2-C_enhancer_seed104729 | 3.13.0 | 2.20.0 | 3.14.1 | 2.5.0 | 32029844 |
| V2-C_enhancer_seed130363 | 3.13.0 | 2.20.0 | 3.14.1 | 2.5.0 | 32029845 |
| V2-C_enhancer_seed155921 | 3.13.0 | 2.20.0 | 3.14.1 | 2.5.0 | 32030334 |
| V2-C_h3k27me3_seed104729 | 3.13.0 | 2.20.0 | 3.14.1 | 2.5.0 | 32031469 |
| V2-C_h3k27me3_seed130363 | 3.13.0 | 2.20.0 | 3.14.1 | 2.5.0 | 32032032 |
| V2-C_h3k27me3_seed155921 | 3.13.0 | 2.20.0 | 3.14.1 | 2.5.0 | 32021359 |
| chr7_inference | 3.13.0 | 2.20.0 | 3.14.1 | 2.5.0 | 31881465 |

### Scheduler accounting

[slurm_accounting.tsv](provenance/slurm_accounting.tsv) preserves scheduler allocation and job-step accounting. Allocation and `.batch`/other step rows overlap and must not be summed together. Allocated GPU time is reserved resource time, not measured GPU-kernel utilization.

| JobID | State | ExitCode | ElapsedRaw | AllocTRES | TotalCPU | MaxRSS |
| --- | --- | --- | --- | --- | --- | --- |
| 32021359_0 | COMPLETED | 0:0 | 1951 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 26:06.251 | — |
| 32021359_0.batch | COMPLETED | 0:0 | 1951 | cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 26:06.251 | 26680639K |
| 32021359_0.extern | COMPLETED | 0:0 | 1951 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 00:00:00 | — |
| 32021359_1 | COMPLETED | 0:0 | 2617 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 35:11.611 | — |
| 32021359_1.batch | COMPLETED | 0:0 | 2617 | cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 35:11.611 | 26854034K |
| 32021359_1.extern | COMPLETED | 0:0 | 2617 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 00:00:00 | — |
| 32021359_2 | COMPLETED | 0:0 | 2637 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 36:49.411 | — |
| 32021359_2.batch | COMPLETED | 0:0 | 2637 | cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 36:49.411 | 26660596K |
| 32021359_2.extern | COMPLETED | 0:0 | 2637 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 00:00:00 | — |
| 32021359_3 | COMPLETED | 0:0 | 488 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 06:09.505 | — |
| 32021359_3.batch | COMPLETED | 0:0 | 488 | cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 06:09.505 | 10745391K |
| 32021359_3.extern | COMPLETED | 0:0 | 488 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 00:00:00 | — |
| 32021359_4 | COMPLETED | 0:0 | 520 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 06:37.788 | — |
| 32021359_4.batch | COMPLETED | 0:0 | 520 | cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 06:37.788 | 10811503K |
| 32021359_4.extern | COMPLETED | 0:0 | 520 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 00:00:00 | — |
| 32021359_5 | COMPLETED | 0:0 | 489 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 06:13.833 | — |
| 32021359_5.batch | COMPLETED | 0:0 | 489 | cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 06:13.833 | 10937381K |
| 32021359_5.extern | COMPLETED | 0:0 | 489 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 00:00:00 | — |
| 32021359_6 | COMPLETED | 0:0 | 939 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 12:30.208 | — |
| 32021359_6.batch | COMPLETED | 0:0 | 939 | cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 12:30.208 | 20917199K |
| 32021359_6.extern | COMPLETED | 0:0 | 939 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 00:00:00 | — |
| 32021359_7 | COMPLETED | 0:0 | 1006 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 13:14.085 | — |
| 32021359_7.batch | COMPLETED | 0:0 | 1006 | cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 13:14.085 | 21224270K |
| 32021359_7.extern | COMPLETED | 0:0 | 1006 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 00:00:00 | — |
| 32021359_8 | COMPLETED | 0:0 | 1075 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 11:57.934 | — |
| 32021359_8.batch | COMPLETED | 0:0 | 1075 | cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 11:57.934 | 20886710K |
| 32021359_8.extern | COMPLETED | 0:0 | 1075 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 00:00:00 | — |
| 32021359_9 | COMPLETED | 0:0 | 375 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 04:06.191 | — |
| 32021359_9.batch | COMPLETED | 0:0 | 375 | cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 04:06.191 | 6767534K |
| 32021359_9.extern | COMPLETED | 0:0 | 375 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 00:00:00 | — |
| 32021359_10 | COMPLETED | 0:0 | 299 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 03:52.391 | — |
| 32021359_10.batch | COMPLETED | 0:0 | 299 | cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 03:52.391 | 6953041K |
| 32021359_10.extern | COMPLETED | 0:0 | 299 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 00:00:00 | — |
| 32021359_11 | COMPLETED | 0:0 | 248 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 02:53.153 | — |
| 32021359_11.batch | COMPLETED | 0:0 | 248 | cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 02:53.153 | 6945216K |
| 32021359_11.extern | COMPLETED | 0:0 | 248 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 00:00:00 | — |
| 32021359_12 | COMPLETED | 0:0 | 814 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 11:29.946 | — |
| 32021359_12.batch | COMPLETED | 0:0 | 814 | cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 11:29.946 | 19607863K |
| 32021359_12.extern | COMPLETED | 0:0 | 814 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 00:00:00 | — |
| 32021359_13 | COMPLETED | 0:0 | 816 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 09:48.507 | — |
| 32021359_13.batch | COMPLETED | 0:0 | 816 | cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 09:48.507 | 19337703K |
| 32021359_13.extern | COMPLETED | 0:0 | 816 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 00:00:00 | — |
| 32021359_14 | COMPLETED | 0:0 | 694 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 09:33.125 | — |
| 32021359_14.batch | COMPLETED | 0:0 | 694 | cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 09:33.125 | 19607824K |
| 32021359_14.extern | COMPLETED | 0:0 | 694 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 00:00:00 | — |
| 32021359_15 | COMPLETED | 0:0 | 201 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 02:16.576 | — |
| 32021359_15.batch | COMPLETED | 0:0 | 201 | cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 02:16.576 | 5581058K |
| 32021359_15.extern | COMPLETED | 0:0 | 201 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 00:00:00 | — |
| 32021359_16 | COMPLETED | 0:0 | 196 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 02:49.728 | — |
| 32021359_16.batch | COMPLETED | 0:0 | 196 | cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 02:49.728 | 5651954K |
| 32021359_16.extern | COMPLETED | 0:0 | 196 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 00:00:00 | — |
| 32021359_17 | COMPLETED | 0:0 | 172 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 01:51.065 | — |
| 32021359_17.batch | COMPLETED | 0:0 | 172 | cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 01:51.065 | 5595001K |
| 32021359_17.extern | COMPLETED | 0:0 | 172 | billing=24,cpu=8,gres/gpu:a100=1,gres/gpu=1,mem=48G,node=1 | 00:00:00 | — |


### Infrastructure failures and retries

The checkpoint freeze records no phase-II infrastructure retries; each prescribed seed/configuration/context has exactly one successful fit.

Nonzero command exits (including administrative failures or scientific STOP codes, distinct from failed model fits):

| Command | Exit code | Preserved completion record |
| --- | --- | --- |
| freeze_C_release_for_calibration | 2 | provenance/commands/freeze_C_release_for_calibration.completed.json |

The original calibration-gate command stopped before creating a release or calibration output because the evaluation manifest used repository-relative artifact paths while the gate helper resolved them from the stage directory. An independently reviewed, append-only metadata-path adapter verified the original bytes and hashes and resolved those references without modifying the frozen helper, evaluation manifest, metrics, selected checkpoints or scientific design. Training, inference and evaluation were not rerun. The correction and exact original-to-resolved mappings are preserved in [calibration_gate_path_resolution.json](provenance/calibration_gate_path_resolution.json).

## Independent validation and preservation

| Audit record | Phase | Status | Checks | Failures |
| --- | --- | --- | --- | --- |
| provenance/calibration_final_validation.json | calibration | PASS | 19 | 0 |
| provenance/checkpoint_packaging_final_validation.json | packaging | PASS | 222 | 0 |
| provenance/evaluation_final_validation.json | evaluation | PASS | 14545 | 0 |
| provenance/input_integrity_verified_after.json | inputs | PASS | 518 | 0 |
| provenance/input_integrity_verified_before.json | inputs | PASS | 518 | 0 |
| provenance/phaseI_cache_independent_validation.json | cache | PASS | 21 | 0 |
| provenance/predictions_final_validation.json | predictions | PASS | 165 | 0 |
| provenance/prepared_input_independent_validation.json | prepared | PASS | 78 | 0 |
| provenance/preservation_final_validation.json | preservation | PASS | 837 | 0 |
| provenance/runs_final_validation.json | runs | PASS | 4684 | 0 |
| provenance/tracked_baseline_before.json | baseline | PASS | 835 | 0 |

Phase-level audit applicability and status used for the report's review-readiness wording. Calibration QC is not required when either C context fails or is inconclusive because calibration is scientifically prohibited:

| Phase | Latest applicable independent audit |
| --- | --- |
| inputs | PASS |
| prepared | PASS |
| cache | PASS |
| runs | PASS |
| packaging | PASS |
| predictions | PASS |
| evaluation | PASS |
| calibration | PASS |
| preservation | PASS |

Applicable execution QC complete/PASS: yes. A valid scientific STOP bundle can have complete execution QC without being ready for holdout evaluation.

Historical pretraining-1.0 and pretraining-1.1 remain immutable inputs. Shared-register activity, decision and result updates are recorded separately from the frozen scientific artifacts. The independent validator checks exact interval/role membership, source hashes, all prescribed seeds, RNG histories, earliest-minimum checkpoints, serialized outputs, probability reductions, component-bootstrap calculations, frozen adequacy rules and conditional calibration. A PASS artifact audit does not convert failed scientific adequacy into PASS.

## Interpretation limits and required stop

This is an internally controlled, single-donor region-label experiment, not independent population validation. Bulk same-lobe provenance does not establish matched cells or aliquots. H3K27me3-associated labels are not experimentally proven silencer labels. Matching reduces specified measured imbalances but cannot remove unmeasured confounding, and the repaired design does not imply local-caliper support for every matched pair. Improved region AP/AUROC/Brier does not establish biologically valid REF–ALT effects, causal variants, fine-mapping probabilities or variant-level false-positive control.

Stop after this internal stage. Do not score chr8–9, open the external functional benchmark, score the candidate universe, generate a candidate list or redesign/rerun this version in response to model outcomes. Any future test/external stage and any scientific redesign require separate investigator authorization. The frozen checkpoints, ensemble implementation and any permissible region thresholds must remain fixed.

## Summary-source integrity

This report formats existing stored summaries only. It does not rerun models, recompute metrics, bootstrap, recalibrate thresholds or read row-level scores. Display rounding is descriptive; all decisions and exact threshold values originate from the linked frozen records.

| Summary source | SHA-256 |
| --- | --- |
| provenance/checkpoint_freeze.json | d2330e4881c477cb7d8e5b2b621940674761133b779c89d0db1ad88635eee379 |
| inputs/input_manifest.json | 4c56f1aa34bc69da121e825d31c70c46d7b3de9be769d5b26a9c426d476bb538 |
| cache/cache_manifest.json | 671cf1161763a712a89c186200b449f0cd30892e2280018ba1427db5de55ce66 |
| predictions/real_network_invariance.json | 17cc6c3c247560554936555dbad5628c5aab70c27dbf36d410d664c2b8009940 |
| results/training_outcomes.tsv | 3ac836ac34f1a45816dbb5ea3cf8cffeb045ba304bb53d700a0e1e1dc84fe084 |
| results/chr7_evaluation/adequacy_decision.json | 267b7eeec486801fc5e0716089710d5c61dbb3ce25a64457dac39fed50adaf5e |
| results/chr7_evaluation/metrics.tsv | a38c1576534d6f09e61915d59287e0efa49a967361a4c7e41947aac24a030bf0 |
| results/chr7_evaluation/seed_stability.tsv | 7854bbab9c969118de3d9717f3fdd9be6c66d446a0489d148f85cdb1ac8ce2d1 |
| results/chr7_evaluation/bootstrap_summary.tsv | 1fa4756614c970950c8099b85a942d1472ade3a574fd7ec48349201751139c8c |
| results/chr7_evaluation/paired_ablation.tsv | 4659df59cea046b36b9d7383b4ace0e38f7dffea33b27e5d8613a93e33461a14 |
| results/chr7_evaluation/C_absolute_adequacy.tsv | 02a59ed7c198802fa1286fa3ee7bee62ef49daa2446124b023a9da8a6cc27b2a |
| results/chr7_evaluation/bootstrap_audit.json | 8f701f55ba189b0464aa2a2316f51950f7c9a09af649ce4baded341eb35829a8 |
| results/chr7_calibration/C_region_thresholds.json | 48f9a7a742a32ebb39d21ecfcde71c340cb2702ef124467d84ce40d0e33037a8 |
| provenance/checkpoint_packaging_final_validation.json | 50f731d56edbfdfca0c2418fd1dd4e33214ab6ed146a22c92c370707c871c8d3 |
