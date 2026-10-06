# COPD V2 model-redesign preflight: pretraining-1.1

Versioned, outcome-blind data-construction repair of the failed `pretraining-1.0`
snapshot. Baseline Git commit: `71f7b3105184d7e4425aa39a6d1d441b20bdc473`.

**Freeze disposition: READY FOR INVESTIGATOR TRAINING REVIEW. No training is
authorized.** The authoritative local freeze and payload hashes are in
[freeze.json](provenance/freeze.json) and
[artifact_checksums.tsv](provenance/artifact_checksums.tsv); the explicit
construction-versus-execution gates are in
[training_readiness.tsv](provenance/training_readiness.tsv).

All historical files, including the original preflight report, failed readiness
tables, matching specification, scripts, checksum ledgers and freeze, remain
unchanged. This version occupies only `pretraining-1.1/`. No GPU training,
phase-I feature extraction, model inference, candidate scoring, external
benchmark evaluation, fine-mapping, target-gene analysis, manuscript revision,
commit or push is part of this repair.

## 1. What changed, and what did not

The repair retains all V1 B positives and all same-lobe C positives, selects
unique 1:1 accessible controls, and rematches C independently. It changes matching
from mandatory local-caliper pairs to a prespecified **distributional** procedure
with unchanged SMD, chromosome and lobe gates. It also fits ATAC normalization
using training chromosomes only and removes bounded-neighbor tie ambiguity.

A and B are diagnostic ablations. C is the intended corrected model, contingent
on future absolute internal adequacy, seed stability, invariance and calibration
in both model contexts. C does not have to beat or be noninferior to B. Failure
of C means stop, not automatic fallback to A or B.

Unchanged: frozen TREDNet phase-I weights/representation; phase-II architecture,
optimizer, batch/epoch/checkpoint framework; RC augmentation and exact symmetric
forward/RC probability averaging; seeds **104729, 130363, 155921**; donor
ENCDO520EJG severe-emphysema context; chromosome and chr7 component-role assignments;
validation-only calibration; benchmark firewall; and absence of a new
allele-delta cutoff. H3K27me3 association is not relabeled experimentally proven
silencer function. The source/architecture records in 1.0 are immutable inputs.

## 2. Matching development and deterministic semantics

The specification was written and hashed before attempt 001 executed. Its inputs
were construction annotations only: GC, training-normalized ATAC peak strength,
repeat coverage, chromosome/role, lobe provenance, sequence/QC and eligibility.
No model score, GWAS priority, known functional-variant list, benchmark label,
AP/AUROC/Brier result or future model outcome entered matching.

Attempt 001 has two recorded stages, not hidden alternative implementations:

1. Match every target to a unique eligible control. Process targets by frozen
   `SHA256(271828|model|configuration|interval_id)`. Prefer the same
   role/lobe-signature/chromosome pool, then same role/signature, then same role.
   Roles never mix train, chr7 checkpoint/selection/calibration or chr8–9 test.
   Lobe signatures and chromosomes may differ pairwise only in the explicit
   fallback tiers; their distributional gaps remain hard-gated.
2. If needed, exchange selected and unused controls within the same
   role/chromosome/signature, preserving those margins, to improve the five
   covariate mean residuals. The complete deterministic proposal ordering and
   actual-loss acceptance rule are specified in
   [matching_attempt_001.json](specification/matching_attempt_001.json).
   The method is a reproducible heuristic, not a claimed globally optimal
   assignment or globally optimal Cartesian exchange solution.

Distances use float64 Chebyshev distance on GC/0.05, ATAC/0.20 and repeat/0.20.
These denominators are **distance scales, not retained hard pair calipers**.
Initial nearest matching is global within the first available priority tier:
after locating a nearest available distance, enumerate the complete expanded
radius, recompute exact distances and break exact ties by natural chromosome,
numeric start/end, then interval ID. A bounded-k query locates the radius but
never resolves ties. All exchange projection/proposal ties likewise use complete
global coordinate ordering.

The initial full-population selection failed four enhancer-B partition-level
SMD checks: train GC 0.112010, train repeat −0.125824, test GC 0.125912 and test
repeat −0.119250. The prespecified exchange stage made **1,371** accepted
improvements (1,175 train; 196 test). C and both H3K27me3 constructions required
no exchanges. Initial pairs, initial failures, final pairs, every accepted
exchange, round summaries and both diagnostic stages remain available under
[attempt 001](attempts/001_full_population/).

No alternative matching implementation or common-support restriction was needed.
The positive-exclusion table is intentionally header-only. B still represents
the complete V1 positive population. The initial/final stages and historical
1.0 failure are distinguished in the attempt ledger; 1.0 was not rerun.

## 3. Retention and class counts

Every B/C row below has an equal number of **unique** controls; none is duplicated
to manufacture coverage. Retention is 100% relative to each configuration's
positive universe, including separately in every chromosome partition.

| Model | Configuration | Train positives / controls | chr7 positives / controls | chr8–9 positives / controls | Total positives / controls |
|---|---|---:|---:|---:|---:|
| Enhancer | B | 156,789 / 156,789 | 8,460 / 8,460 | 16,103 / 16,103 | 181,352 / 181,352 |
| Enhancer | C | 125,874 / 125,874 | 6,631 / 6,631 | 12,695 / 12,695 | 145,200 / 145,200 |
| H3K27me3-associated | B | 26,326 / 26,326 | 1,612 / 1,612 | 2,646 / 2,646 | 30,584 / 30,584 |
| H3K27me3-associated | C | 17,132 / 17,132 | 1,096 / 1,096 | 1,754 / 1,754 | 19,982 / 19,982 |

A keeps exact V1 positive/control membership: enhancer 181,352 / 362,704 and
H3K27me3-associated 30,584 / 61,168. Its known historical task defects are not
silently repaired. Six rich manifests and 36 class/split BED files are provided.

The direct same-lobe label correction removes 36,152 enhancer positives
(19.935%) and 10,602 H3K27me3-associated positives (34.665%) from B's universe;
this is unchanged biological membership evidence, not a new matching exclusion.
The B→C model contrast is explicitly **same-lobe label correction plus induced
control rematching**, not a pure label-only contrast. C has 3,282 enhancer
controls and 606 H3K27me3 controls not selected in B; its control sets were not
copied unchanged. Shared controls across configurations are permissible staged
ablation reuse, not duplication within a configuration.

## 4. Complete balance metrics and gate interpretation

The hard limits remain `|SMD| <= 0.10`, chromosome proportion gap `<= 0.02`,
and ATAC-lobe-signature proportion gap `<= 0.05`. Individual-lobe support gaps
are additionally checked against 0.05. All five SMD variables and all represented
chromosome/signature/lobe categories pass in train, full chr7 and chr8–9.
The chr7 checkpoint/selection/calibration roles also pass separately.

Signed SMD is `(positive mean - control mean) / pooled sample SD`, `ddof=1`.
Values below are rounded for display; gates use unrounded values with no epsilon
relaxation.

| Model/configuration | Partition | GC SMD | Training-only ATAC SMD | Repeat SMD | Full-input blacklist-any SMD | Non-ACGT-any SMD |
|---|---|---:|---:|---:|---:|---:|
| Enhancer B | train | 0.095965 | 0.076598 | −0.099993 | −0.002789 | −0.001458 |
| Enhancer B | chr7 | 0.035039 | 0.025056 | −0.029606 | −0.015376 | −0.026634 |
| Enhancer B | chr8–9 | 0.099992 | 0.075265 | −0.074358 | −0.011145 | 0 |
| Enhancer C | train | 0.089785 | 0.051493 | −0.092628 | −0.008722 | −0.015377 |
| Enhancer C | chr7 | 0.032731 | 0.016379 | −0.020353 | −0.017367 | −0.030085 |
| Enhancer C | chr8–9 | 0.099411 | 0.051747 | −0.088504 | −0.025106 | 0 |
| H3K27me3 B | train | 0.026034 | 0.002921 | 0.001701 | −0.003899 | 0.014800 |
| H3K27me3 B | chr7 | 0.038280 | −0.001440 | −0.004148 | 0.035223 | −0.049829 |
| H3K27me3 B | chr8–9 | 0.016172 | 0.002073 | 0.005922 | 0 | 0 |
| H3K27me3 C | train | 0.030121 | 0.002605 | 0.002396 | 0.003258 | 0.013669 |
| H3K27me3 C | chr7 | 0.042454 | 0.003607 | −0.005524 | 0.042718 | −0.042718 |
| H3K27me3 C | chr8–9 | 0.013267 | 0.003166 | 0.007863 | 0 | 0 |

Maximum absolute category gap within each partition:

| Model/configuration | Partition | Chromosome | Lobe signature | Individual lobe support |
|---|---|---:|---:|---:|
| Enhancer B | train | 0.003323 | 0 | 0 |
| Enhancer B | chr7 | 0 | 0.000355 | 0.000236 |
| Enhancer B | chr8–9 | 0.000559 | 0.000621 | 0.000621 |
| Enhancer C | train | 0.000802 | 0 | 0 |
| Enhancer C | chr7 | 0 | 0.000302 | 0.000151 |
| Enhancer C | chr8–9 | 0.000394 | 0.000473 | 0.000315 |
| H3K27me3 B and C | all three partitions | 0 | 0 | 0 |

C rematching resolves the former GC/repeat/lobe failures under the same gates.
However, enhancer B test GC, B train repeat and C test GC are close to 0.10.
Passing is not a large safety margin, proof of distributional identity or proof
that unmeasured confounding is absent. No gate was loosened or redefined after
the result. Quantiles, dispersions, missingness and source counts are reported
alongside mean-balance measures.

The full unrounded values, every category and every chr7 role are in
[covariate_balance.tsv](attempts/001_full_population/covariate_balance.tsv),
[covariate_distributions.tsv](results/covariate_distributions.tsv) and
[lobe_provenance_balance.tsv](results/lobe_provenance_balance.tsv).

## 5. Local-distance qualification

Full positive retention is **not** a claim of 100% local-caliper common support
or pairwise biological exchangeability. The user-prioritized distributional
construction deliberately removes the old pair-distance limits while preserving
all scientific balance/contamination gates.

| Model/configuration | Final pairs outside at least one old distance scale | Fraction | Maximum normalized distance |
|---|---:|---:|---:|
| Enhancer B | 17,094 / 181,352 | 9.4259% | 7.8761 |
| Enhancer C | 9,478 / 145,200 | 6.5275% | 6.5967 |
| H3K27me3 B | 180 / 30,584 | 0.5885% | 5.1374 |
| H3K27me3 C | 143 / 19,982 | 0.7156% | 3.0585 |

Enhancer B/C have respectively 1,136/258 cross-chromosome pairs and 13/8
cross-signature pairs. H3K27me3 has none. All remain within the identical fixed
partition/chr7 role. These are bookkeeping pair links, not biological replicate
or bootstrap units. The final distributional exchange can increase an individual
pair's distance while improving the relevant marginal balance; this tradeoff is
retained, not hidden. Initial/final and per-feature/partition/role distance tails
are in [matched_pair_distance_summary.tsv](results/matched_pair_distance_summary.tsv).

## 6. Training-derived ATAC normalization

For each source, fit the midpoint empirical CDF on training chromosomes only:
`F(x) = (count(training signal < x) + 0.5*count(training signal = x))/n_train`.
Ties share the midpoint rank; interior unseen values use the strict-less count;
values below/above the training range map to 0/1. No interpolation, chr7 fit,
test fit or held-out refit occurs. Identical interval anchors take the maximum
of their independently mapped source-peak values. This is matching metadata,
not a neural-network input.

Training denominators: ENCFF189LZA 189,831; ENCFF899TQV 239,613;
ENCFF906HOT 212,067. All 745,045 source peaks and 1,143,370 master interval IDs
are reconciled, including 725,105 anchored intervals. Four test peaks exceed
their source's training maximum and correctly map to 1. Missing rank values in
non-anchored historical A controls are preserved as missing, not imputed.

The full frozen CDF knots, sidecar, specification, denominator table and raw
source reconciliation are supplied. Original 1.0 rank strings remain in the
rich manifests for provenance; only the explicitly named
`atac_signal_percentile_max_train_only` drives repaired matching and balance.

## 7. Common C-task panels and unchanged role separation

The common panels use final same-lobe C positives and independently C-matched
controls. They contain selection, calibration and test roles only; checkpoint
rows stay in native training/checkpoint manifests. All A/B/C models may later
be evaluated on these same panels for paired comparability. Native-task results
remain descriptive.

| Model | Role | Positives | Controls | All components | Control components |
|---|---|---:|---:|---:|---:|
| Enhancer | chr7 selection | 1,755 | 1,755 | 1,295 | 885 |
| Enhancer | chr7 calibration | 1,576 | 1,576 | 1,157 | 809 |
| Enhancer | chr8–9 test | 12,695 | 12,695 | 8,754 | 5,877 |
| H3K27me3 | chr7 selection | 271 | 271 | 244 | 160 |
| H3K27me3 | chr7 calibration | 263 | 263 | 259 | 176 |
| H3K27me3 | chr8–9 test | 1,754 | 1,754 | 1,809 | 1,157 |

Each selection panel exceeds the preserved 100-positive/200-control/30-component
minimum. Calibration exceeds 200 controls and an explicitly clarified minimum
of 30 independent control components. Components retain the original global
genomic-overlap/encoded-RC grouping and hashed chr7 assignment; matching pairs
are not used to invent independent biological units or reassign roles.

Chr8–9 remains a fixed internal holdout from V2 training/model selection, but is
**not historically untouched**. This repair uses its permitted construction
covariates and QC, including matching exchanges. No test model score, performance
metric or threshold is computed or used. Final construction/specifications are
frozen before any V2 inference.

## 8. Prospective model adequacy and calibration

The complete rule is [selection_adequacy.json](specification/selection_adequacy.json).
For both C model contexts, require all three seeds, valid finite probabilities,
real-network symmetric-inference invariance, lower 95% component-bootstrap CI
for AUROC >0.5, AP minus replicate prevalence >0 and Brier skill >0. Seed AP
sample SD >0.03 or range >0.10 stops retention. The latter range stop is redundant
for three seeds but retained, not loosened. Require 2,000 valid bootstrap draws,
the frozen invalid-draw safeguards and sample/component minima.

A/B/C differences retain paired component-bootstrap uncertainty, but comparative
AP/AUROC/Brier margins no longer choose the final configuration. C must satisfy
absolute criteria; if it fails, stop. No A/B or mixed-head fallback, favorable
seed subset, test tie breaker or external tie breaker is allowed.

Threshold calibration uses only the dedicated common chr7 calibration controls
after the retained C ensemble is fixed: with `k=floor(0.05*n)`, take
`nextafter(descending_negative_scores[k], +infinity)` and call `score >= threshold`.
Ties are conservative. No calibration positive, test or external score chooses
the threshold. This is a region-label operating threshold, not a calibrated
allele-effect FDR, causal-variant list or new allele-delta cutoff.

The CPU specification review passes 95 invariant/synthetic checks. No model
adequacy, real-network invariance or calibration result has been generated;
those remain future gates after separate authorization.

## 9. Sequence, contamination, provenance and independent validation

Controls retain donor ATAC anchoring, absence of relevant mark in every original
supporting peak and the full 2,001-bp input across lobes, exclusion of same-model
known positive full inputs, historical core blacklist/TSS filters and available
sequence. Full-input flank blacklist/promoter annotations remain audited rather
than retroactively changing A or positive definitions. Bulk-lobe support does
not establish matched-cell/aliquot biology or healthy-donor generalization.

The selected union has 792,681 unique intervals and 792,672 encoded sequence
identities. Nine encoded forward/RC duplicate groups remain within training;
none crosses chromosome partitions or chr7 roles. No selected raw exact-duplicate
group crosses a boundary. The independent audit checks original peak/mark/source
annotations, full-input contamination, all selected FASTA sequences, hashes,
roles/components, complete initial matching ties and the full exchange ledger.

Synthetic nearest-tie conformance covers equal-distance sets larger than the
bounded query, exhausted candidates, shuffled coordinate indices and nearby
unequal distances. A bounded constructor robustness caveat outside valid inputs
is documented: its SMD helper can mask undefined singleton/nonfinite variance.
Readiness therefore additionally requires the independent audit to establish
finite features/pooled variances and at least two rows per class in every B/C
scope; the actual smallest constructed class/role has 263 observations.

The corrected independent final audit passes **1,154 checks with zero artifact
failures or readiness blockers**. It replays all **377,118** initial assignments,
including **6,868** exact global ties, and all **1,371** accepted exchanges;
checks every selected source/FASTA sequence; and independently verifies all
balance, contamination, partition/role, sample/component and historical-input
conditions. See
[validation_corrected_final_manifest.json](provenance/validation_corrected_final_manifest.json)
and its linked check, balance, readiness and tie tables.

A separate packaging audit passes **625 checks**, and a supplementary pair-link
audit passes **57 checks** covering all 754,236 B/C manifest endpoints; all
635,808 A rows correctly have empty pair fields. These independently reconcile
the six rich manifests, 36 BEDs, 18-run matrix, distributions/distances, compute
projection, duplicate accounting and the finite-SMD prerequisites. All **708**
historical tracked files remain byte-identical, including 1.0 and the benchmark.

The first independent validation attempt is preserved, not silently replaced:
it had 1,071 checks and two validator-only failures. One incorrectly compared an
absent old `external_recalibration` field with the new explicit `false`
prohibition. The other used Python decimal parsing rather than the frozen
constructor's pandas 2.0.3 input decoding, turning a one-ulp distance difference
into an artificial tie for one enhancer-B target. The corrected validator
independently reconstructs the source mapping, then replays the exact constructor
numeric inputs, with strict tie identities and **no tolerance relaxation**.
Its complete repeat audit passes. The original failed outputs, original
validator source, corrected source, exact float/operand diagnostic and resolution
record are all retained. No matching, model specification, scientific data or
gate was changed to make the validator pass.

Intermediate assembly/attempt manifests truthfully retain their historical
"pending independent validation" state. The final corrected audit and local
freeze are authoritative for this completed version. This is not a revision or
reinterpretation of the permanently preserved failed 1.0 freeze.

## 10. Compute projection, readiness and stop

The design still contains 18 planned fits: three configurations × two contexts
× three seeds. None is launched or authorized. Count-scaled fit-only estimate:
**5.71305 A100 GPU-hours**. Retain the conservative **12–24 GPU-hour** total
planning range for two-orientation representation extraction, symmetric
validation, I/O and overhead. This is an estimate, not a measured V2 runtime.

Planned per fit: one A100, 8 CPUs, 48 GiB RAM and 4-hour walltime cap, at most
four concurrent fits. The aggregate reservation ceiling is 72 GPU-hours, not
expected use. Both-orientation float32 caches for 792,672 encoded identities
would occupy 28,916,674,560 bytes (**26.9308 GiB**); scratch budget remains
120 GiB. Cache train/chr7 first only after authorization; defer test model
execution until the retained-model freeze.

Observed CPU-only stages: ATAC normalization 41.67 s, matching 100.76 s,
assembly 98.98 s. Independent validation costs are separately recorded in its
manifest. No phase-I features or model predictions were produced.

Construction gates, independent source/sequence/tie/hash validation and packaging
all pass with full retention and independent C rematching. The local checksum
freeze preserves the complete repaired version. The remaining investigator decision is review of the distributional
matching tradeoff, narrow enhancer balance margins and C-only absolute-adequacy
rule, followed—if accepted—by **separate explicit GPU authorization**. Stop after
this local preflight; do not launch the 18 fits.
