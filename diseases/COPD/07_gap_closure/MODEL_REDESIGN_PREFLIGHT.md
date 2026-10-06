# COPD V2 model-robustness pretraining preflight

Module: `COPD-V2-PREFLIGHT`, experimental version `pretraining-1.0`.
Scope: design, CPU-only data construction, provenance/sequence audits, and an
immutable pretraining freeze. **GPU training and phase-I feature extraction are
not authorized.** No new candidate ranking, external evaluation, retraining,
architecture search, or Git publication is part of this module.

This document and the machine-readable construction/selection specifications
must be frozen together before training can be approved. Numerical construction
results and readiness are recorded in the final construction/QC manifests.

**Preflight disposition: constructed and frozen for review; NOT READY FOR
TRAINING.** Matching fails its prespecified covariate-balance gates. This is a
data-construction finding, not a model-performance result. Freezing preserves
the complete candidate experiment and its failures; it does not certify a
successful matched experiment or authorize GPU use. No matching rule or model
configuration has been retuned to make these diagnostics pass.

## 1. Scientific hypotheses and limits

H1: RC augmentation plus mathematically symmetric inference removes orientation
dependence without replacing the frozen TREDNet representation or phase-II
model family. Augmentation may improve what is learned; symmetric inference
enforces output invariance independently of whether augmentation succeeds.

H2: changing the historical multi-source DHS controls to donor-accessible,
relevant-mark-absent controls changes the learned contrast from a partly
accessibility/source-discriminating task to regulatory-mark classification
conditional on donor accessibility. An AUROC decrease on this harder task is
not, by itself, evidence that the correction failed.

H3: requiring same-lobe ATAC/mark co-occurrence materially changes positive
labels relative to the V1 all-lobe union. The experiment tests this label
construction effect; it cannot identify a cell-intrinsic causal histone effect.

H4: three independently initialized phase-II fits expose optimization/seed
instability hidden by the one-seed historical reference. Three seeds do not
estimate the full distribution of possible training runs precisely.

The donor remains ENCDO520EJG (severe-emphysema bulk lung). Generalization to
other donors, healthy controls, or cell types is not established by this design.
H3K27me3 association is not interchangeable with experimentally proven silencer
function. Regulatory-region labels are not allele-effect labels.

## 2. Staged configuration matrix

| Configuration | Positive labels | Negative controls | Orientation | Fits |
|---|---|---|---|---:|
| Frozen V1 | Exact frozen V1 | Exact frozen V1 | Original frozen forward scores | 0 |
| V1 symmetric diagnostic | Unchanged reference | Unchanged reference | Equal forward/RC probability mean; separate diagnostic output only | 0 |
| V2-A | Exact V1 | Exact V1, including known construction defects | RC augmentation + symmetric inference | 2 models × 3 seeds = 6 |
| V2-B | Exact V1 | Same-donor accessible matched, relevant-mark-absent controls | Same as A | 6 |
| V2-C | Same-lobe-supported subset of V1 positives | Prespecified accessible-control procedure; B/C identity preferred to isolate labels | Same as A | 6 |

The primary matrix is 18 fits, not an architecture/hyperparameter search.
All rows preserve the frozen phase-I weights, phase-II layer structure,
chromosome framework, donor, sequence geometry, optimizer, and basic batch
semantics. Test/external outcomes cannot add configurations or choose seeds.
No V1 symmetric scores need to be computed during this pretraining phase.

V1→A is an orientation/validation-robustness package, not a pure augmentation
effect: seed aggregation, checkpoint-role separation and threshold-leakage
repair also change. Frozen V1 forward versus symmetric inference is the isolated
inference diagnostic. B→C keeps selected controls identical in this construction,
making the label change interpretable but subject to a matching-readiness gate.

## 3. Frozen data definitions and V1 preservation

BED coordinates are zero-based half-open. The V1 positive definition is at
least 1 bp overlap of an **original ATAC peak** with a relevant histone peak,
not overlap of an arbitrarily recentered sequence. Its core window is
`[summit-500, summit+500)`. The exact historical model input is
`[core_start-501, core_end+500)`, length 2,001 bp. This asymmetric expansion
places input index 1,000 one base before the original summit; it is retained
across A/B/C, not silently corrected as an extra factor.

V1 exact duplicate core coordinates were removed before training. A retains
the exact four frozen V1 BED sets, their labels, and their original controls.
Ambiguous reference bases remain zero-valued one-hot rows as in V1; availability
and ambiguity fractions are audited rather than silently filtering A.

All nine ENCODE peak files, the historical BEDs, hg38 sequence and index,
GENCODE v50 promoter definition, hg38 blacklist v2, RepeatMasker track, scripts,
and model weights are identified by path and SHA-256. Source files and V1
outputs are read-only. The completed benchmark is hash-checked as opaque bytes
only; no benchmark outcomes are parsed for design.

## 4. Accessible-control construction and confounding audit

For each model, candidate controls must be summit windows anchored in this
donor's ATAC peaks. They must pass the historical core-window blacklist and
GENCODE gene-TSS ±2-kb exclusion. A candidate is rejected if the relevant mark
overlaps any original supporting ATAC peak or any base of its full 2,001-bp
model input in **any** sampled lobe. It must also not overlap a V1-positive
model-input window for that model. A cross-lobe-only positive is excluded from
C positives, never knowingly reassigned to C negatives.

The full-input negative exclusion is a documented part of the B control
correction. It prevents treating a known marked flank as a negative even when
the historical 1-kb core is mark-free. The corresponding defects in A remain
audited historical defects, not repaired in place.

The matching specification is fixed before any model outcomes. It must use
chromosome, GC, source-normalized ATAC peak strength, repeat coverage, identical
window length, and recorded genomic/lobe provenance. ATAC signalValue is
available from narrowPeak column 7; ranks are calculated within each source
file before combining identical anchors. It is a peak-strength proxy, not a
depth-normalized cross-experiment bigWig measurement. Raw signal and all lobe
support remain available for audit. Do not match on motifs, histone intensity,
model predictions, GWAS associations, or benchmark outcomes: these could
remove the biological signal or breach the evaluation firewall.

The implemented target is one control per B positive for both models. The
enhancer pool has 286,943 eligible windows for 181,352 positives, insufficient
for the historical 2:1 ratio even before matching. A common 1:1 target avoids
duplicate reuse and unequal model-specific policies. Exact strata are chromosome,
global chr7 role, ATAC lobe-support signature, any full-input blacklist overlap,
and any non-ACGT base. Process positive targets by ascending
SHA256(`271828|model|interval_id`); select the nearest unused control by Chebyshev
distance after dividing GC difference by 0.05, ATAC signal-percentile difference
by 0.20, and repeat-fraction difference by 0.20. All three calipers must hold.
Within the returned nearest-neighbor candidate batch, coordinate-ID order
resolves ties. The exact bounded-query and fallback semantics, software versions,
and role-hash truncation are recorded in
`provenance/COPD-V2-PREFLIGHT_matching_implementation_notes.json`; this clarifies
the initial abbreviated specification rather than changing selected intervals.
These are transparent design
tolerances, not learned or benchmark-optimized values.

Keep every unmatched positive and record its failure reason. Clean matching
requires at least 90% B-target coverage, absolute SMD ≤0.10, chromosome-proportion
gap ≤0.02, and lobe-signature proportion gap ≤0.05 in each chromosome partition.
Actual ratios may be below target. C retains B controls and therefore has a
different ratio after positive restriction. Unweighted BCE is preserved; changed
class prevalence is reported, not hidden by unregistered loss weighting.

No duplicate control rows, fabricated 2:1 ratio, or silent positive subsampling
is allowed to conceal lack of common support. Matching diagnostics must show
requested versus realized ratios, unmatched strata/targets, distributions,
standardized mean differences, and any balance failures. Absolute SMD >0.10
is a predeclared imbalance flag, not a universal proof of confounding below
that value. Source overlap is established by construction, while technical
cleanliness additionally requires sequence, exclusion, independence, and
matching QC. Blacklist/promoter overlap is audited on both core and full input;
flank findings cannot be silently erased from A or B positives.

## 5. Same-lobe labels and technical compatibility

For each original ATAC peak, histone support must come from the same donor,
anatomical lobe, GRCh38 assembly, and compatible released peak-call context.
The available pairs are upper-right, lower-right, and lower-left lung.
ATAC and ChIP biosample accessions differ: these are donor/lobe-matched bulk
observations, not matched aliquots, cells, or a joint assay.

A genomic summit window is positive in C if at least one of its original ATAC
observations has same-lobe support. Combine the evidence as a union and retain
one row per exact genomic window with all contributing source IDs. Do not give
a window threefold weight because three lobes support it. Nearby overlapping
windows are dependent observations, not duplicate independent validations;
their overlap components are recorded and kept together for validation roles
and uncertainty estimation.

ENCODE warnings are retained, including lower-right H3K27me3 low read depth and
missing-control warnings. No lobe is selected after seeing model performance.
A lobe-quality exclusion would be a separately documented conditional version,
not an automatic fourth primary configuration.

## 6. Orientation handling and mathematical contract

Let `g` be frozen phase I, `h[c,m,s]` the phase-II network, and `R` nucleotide
reverse complement. For every sequence, model, configuration, and seed:

`f_s(x) = h_s(g(x))`

`q_s(x) = (f_s(x) + f_s(R(x))) / 2`

`Q(x) = (q_104729(x) + q_130363(x) + q_155921(x)) / 3`.

Because `R(R(x)) = x`, `q_s(R(x)) = q_s(x)` and therefore `Q(R(x)) = Q(x)`.
RC is applied to the nucleotide sequence **before** phase I. Reversing the
4,560-dimensional epigenomic embedding is not reverse complementation.
Average final probabilities, not embeddings before a nonlinear phase-II layer.

Training augmentation draws one Bernoulli(0.5) orientation for each original
interval in each epoch from a recorded seed-derived RNG stream. It preserves
one example weight per interval and the original number of steps per epoch;
do not double count paired views. Forward/RC frozen phase-I caches are reused
across configurations and seeds, keyed by sequence and model-weight hashes.

Inference disables dropout and BatchNorm updates, uses fixed ordered float64
means, and retains raw orientation scores. REF and ALT are transformed
independently without swapping allele identities; index i maps to 2,000−i.
Required numerical tolerance is `abs(Q(x)-Q(R(x))) <= 1e-6 + 1e-6*abs(Q(x))`,
with finite probabilities. The CPU contract tests use synthetic sequences and
orientation-sensitive mock callbacks only. They verify the wrapper, not any
untrained V2 network. Real network invariance remains a mandatory later gate.

## 7. Seed, architecture, optimization, and checkpoint policy

Training seeds are exactly **104729, 130363, 155921** for every new configuration
and both models. Matching seed is 271828; bootstrap seed is 314159. These values
were fixed without model outcomes. Set Python/NumPy/TensorFlow RNGs, shuffle and
augmentation substreams, and deterministic operations explicitly. All three
seed-level outputs are retained. Infrastructure retries use the same seed and
are logged; failed seeds cannot be replaced by more favorable seeds.

Phase I remains the exact frozen 4,560-output representation. Phase II remains
Conv1D(64,k4)→BatchNorm→pool2→dropout0.4→Conv1D(128,k2)→dropout0.4→flatten→
Dense100→Dense50→sigmoid, with the historical activations. Preserve unweighted
binary cross-entropy, Adadelta learning rate0.001/rho0.95/epsilon1e-7, batch256,
maximum50epochs, and patience15. No adaptive learning-rate, class weighting,
architecture, or optimizer search is introduced silently.

Per seed choose the earliest epoch attaining the minimum **symmetric** BCE on
the checkpoint-only chr7 subset; ties retain the earlier epoch, min_delta=0.
The checkpoint callback evaluates final forward/RC probabilities, not a
forward-only loss. Restore that one checkpoint. All three successful seed
checkpoints are equally averaged; no seed/cherry-picked checkpoint ensemble.

Checkpoint loss uses each configuration's native labels/controls restricted to
the global chr7 checkpoint role. Neither common selection labels nor
calibration-role scores enter that callback. Derive shuffle and augmentation
RNGs with NumPy SeedSequence([training_seed,1,epoch]) and
SeedSequence([training_seed,2,epoch]), respectively. Epoch is zero-based; start
from the frozen coordinate-sorted interval ordering. Phase-II initialization
uses the training seed. Frozen phase-I weights are never randomized or optimized.

## 8. Chromosome independence and validation roles

Train: chr1–6, chr10–22, chrX, chrY. Validation/calibration: chr7 only.
Untouched final internal test: chr8–9. Construction counts/availability and
prespecified QC may be audited on test intervals before training; no model
scores, test metrics, threshold fitting, or performance-guided redesign is
permitted. Fixed prespecified matching calipers and strata are applied mechanically
to every partition. ATAC strength percentiles in this version are assay-wide
annotations ranked within each source file over all canonical chromosomes,
including chr7–9; they are not train-fitted preprocessing. Held-out covariate
distributions therefore contribute to interval construction, while held-out
predictive scores and performance remain unopened. These ranks are control-
selection metadata, not neural-network inputs or normalization of model features.
For an approved replacement design, a per-source training-chromosome-only
empirical CDF, frozen and applied to chr7–9 with fixed tie/extreme handling,
would strengthen strictly inductive preprocessing. That is a prospective
versioned improvement, not a claim that the current annotation was train-only.

Freeze one global assignment of chr7 full-input overlap components to three
roles: approximately 50% checkpoint, 25% configuration selection, and 25%
threshold calibration, using a fixed hash rule. Exact or RC-equivalent encoded
sequence duplicates must share a role, even when genomically disjoint.
Cross-chromosome derived duplicates are a hard audit condition: do not quietly
alter A's exact frozen inputs to make its independence checks pass. Such a
conflict requires an explicit conditional safety version and investigator
approval before training, not an inaccurate claim of clean independence.

Freeze a common same-lobe-positive/accessibility-control chr7 challenge panel
for cross-configuration comparison, and an analogous chr8–9 panel for final
testing. Native A/B/C panels remain descriptive staged-task evaluations.
Native AUPRC values cannot rank configurations with different labels, controls,
or prevalence. All model-independent labels, panels, roles, and hashes are fixed
before training. Matching pairs remain in one validation role by exact matching.
A matching link between distant loci is a fixed sampling-design link, not a
shared biological observation. Uncertainty estimation conditions on the frozen
matching design and resamples genomic/encoded-sequence components rather than
merging distant components by matching links.

## 9. Validation-only threshold calibration

After internal configuration selection and checkpoint/ensemble freeze, use
only negatives in the common chr7 calibration panel. With n negative scores
sorted descending as d1…dn and k=floor(0.05*n), set
`t = nextafter(d[k+1], +infinity)` (one-based indexing). Call `score >= t`.
This guarantees empirical calibration FPR ≤5%, including ties; a percentile
interpolation or nearest-ROC-point procedure is not equivalent. Require n≥200
negative controls; report uncertainty and do not promise 5% population/test FPR.
Store the binary float and at least 17 significant decimal digits.

Apply the frozen threshold unchanged to chr8–9. Never optimize thresholds on
test negatives as V1 did. Any historical V1 cutoff is identified as historical
and test-derived, not leak-free V2 calibration. Per-seed thresholds may be
reported descriptively, but the operative threshold belongs to the fixed
three-seed symmetric ensemble. No post-hoc external threshold tuning is allowed.

## 10. Internal configuration selection before test/external opening

Primary comparison: average precision (AP, the stepwise definition, not a
trapezoidal PR area) on the unchanged common chr7 selection panel. AUROC and
Brier are separate guardrails; no weighted omnibus score. Report native class
ratios, common-panel prevalence, reliability by ten fixed probability bins,
Brier skill relative to a constant-prevalence predictor, and seed-level metrics.

Use 2,000 paired overlap-component bootstrap replicates, seed314159. Proposed
absolute margins are 0.02 for AP, AUROC, and Brier: pragmatic pretraining
tolerances, not benchmark-derived or empirically established clinical margins.
Report all point estimates and uncertainty, including insufficient independent
components or undefined replicates. At least 30 independent components and
100 positives/200 controls are required for a definitive selection comparison.

Superior: lower95%CI for ΔAP >0, point ΔAP≥0.02, lower95%CI ΔAUROC≥−0.02,
and upper95%CI ΔBrier≤+0.02. Equivalent: complete paired90%CIs for all three
differences lie within ±0.02. Noninferiority alone is not equivalence; failure
to reject a difference is not equivalence. Noninferiority is separately defined:
lower 95% CI ΔAP ≥−0.02, lower 95% CI ΔAUROC ≥−0.02, and upper 95% CI ΔBrier
≤+0.02. A result can be noninferior without being superior or equivalent.
Remaining results are inconclusive. Apply the rules separately to enhancer and
H3K27me3; full C retention requires both models, not an unregistered mixed choice.

Bootstrap draws use the same component IDs and multiplicities for every
configuration and carry every row of each sampled component. Compute paired
differences per draw; use percentile intervals at 2.5/97.5% or 5/95% for 95% or
90% intervals. Reject one-class draws, allow at most 20,000 attempts, and require
2,000 valid draws with fewer than 10% rejected draws. Otherwise comparison is
inconclusive. Require 30 represented components separately for each model's
common selection panel. These conditional genomic-component intervals do not
establish independent-donor generalization.

Adequacy additionally requires valid construction/independence QC, finite
predictions and network invariance, all three seeds completed, ensemble lower
95% AUROC bound >0.5, lower 95% bound of AP minus replicate-specific prevalence
above zero, and lower 95% bound of Brier skill above zero. Bootstrap these
quantities directly. Seed instability is a hard review stop if AP sample
standard deviation >0.03 or range >0.10;
do not fix it by selecting a favorable seed. Report stratified GC, accessibility,
repeat, chromosome/lobe-support and ambiguous-base results; strata too small
for stable inference remain descriptive and cannot drive a hidden search.

The scientifically preferred retained configuration is C only if adequate and
noninferior to B on the common selection task with the above guardrails. If C
is inconclusive/inadequate, stop for review. A and B are diagnostics retaining
known construction defects, not automatic substitutes for a fully corrected
final V2 model. A later architecture redesign or alternative data definition
requires a new registered version, never the test/external set as tie-breaker.

## 11. Candidate-scoring plan, not a new causal list

After approved training and model freeze, score the unchanged candidate
universe with per-model/per-seed forward and RC REF and ALT probabilities;
report Q_REF, Q_ALT, symmetric region score=max(Q_REF,Q_ALT), delta=Q_ALT−Q_REF,
seed SD/range, raw orientation differences, numerical invariance status, and
validation-calibrated region-threshold status. Unscorable alleles retain an
explicit reason. Averaging never changes REF/ALT identity. No candidate scores,
new ranking, or causal list are produced by this pretraining module.

The 5% operating point is calibrated for a single genomic-sequence region
score. Taking max(REF,ALT) does **not** inherit a guaranteed 5% variant-level FPR.
Region-threshold status is a comparator, not an allele-effect or variant-level
error guarantee.

## 12. Allele-effect calibration is a separate unresolved scientific task

Improved region AUROC does not validate REF–ALT deltas. The V1 candidate-universe
95th-percentile |delta| rule is a historical descriptive comparator only: a
GWAS/LD-enriched universe is not an experimental null, and its percentile is
not a false-discovery or causal-probability calibration. Do not make it the
definitive V2 rule and do not replace it with a benchmark-optimized cutoff.

An eventual allele-effect calibration requires independently justified labels
or nulls, strand-consistent allele representation, assay/mechanism awareness,
and a separately frozen development/calibration split. Internally generated
mutations can assess numerical stability but are not biological negatives.
Without a validated entirely internal method, retain continuous signed deltas
and uncertainty; authorize a separate calibration protocol before binary
allele-effect claims. The frozen external benchmark cannot serve as calibration.

## 13. External benchmark firewall

`COPD-V2-BENCH` version1.0 remains unchanged frozen external evaluation data.
No outcome table, label, recovery rate, or benchmark-derived preference is read
or used for architecture, controls, lobes, augmentation, inference, seeds,
hyperparameters, checkpoint, thresholds, or configuration decisions. Opaque
file hashes are allowed solely to verify preservation. Scripts whitelist V1
training/source inputs and reject external benchmark data as design inputs.

Freeze retained model, all three checkpoints, inference contract, thresholds,
and the internal-selection record before opening the benchmark once for final
external comparison. A poor external result is reported, not used to revise
this version. Freeze before final internal-test opening as well. External
functional labels are distinct from region labels; external genomic loci may
overlap training chromosomes, so genomic independence cannot be assumed from
the word external. Stratify overlap only after model freeze without retuning.

## 14. Compute/runtime plan and parallelization

Historical fits used one A10080GB, not four GPUs: the V1 constant GPUS=4 only
multiplied batch64 to batch256. Both fits ran50epochs; epoch sums were44.3min
enhancer and7.57min H3K27me3. Eighteen full-V1-sized fits imply about7.78GPU-hours
before extra symmetric-validation/cache/I/O work. Budget12–24GPU-hours for the
full approved matrix, to be refined from frozen counts, not external results.
Original MaxRSS was not recovered; runtime estimates are not guarantees.

Count scaling for the constructed datasets gives **5.80 GPU-hours of fitting**
before overhead; the conservative total budget remains **12–24 A100 GPU-hours**.
There are 780,226 selected genomic intervals and 780,217 unique encoded
sequence/RC identities. Two-orientation float32 features occupy approximately
28.46 GB (26.51 GiB), excluding metadata. Reserve about 120 GiB scratch for
cache, immutable checkpoints, predictions, temporary files and manifests.
The 18 separate four-hour caps sum to 72 GPU-hours of maximum reservations,
not the projected utilization. Four-way concurrency suggests roughly 4–8 hours
of execution including cache/overheads, excluding queue waits, if approved.

Observed CPU construction: 160.17 seconds, peak RSS 1.85 GiB, for the complete
source/sequence feature audit; 142.72 seconds for configuration construction.
Independent validation is additional. An 8-CPU/16-GiB/one-hour CPU request is
a conservative reproduction allocation; these stages used no GPU. The current
feature table is 125.6 MB compressed; this local preflight is not being pushed,
and publication packaging must respect GitHub's per-file limit without dropping
rows or changing the frozen data identity.

Future request per fit: one A100,8CPUs,48GBRAM,4h wall-time cap; initially at
most4 concurrent fits. Cache forward/RC phase-I features once per unique encoded
sequence, frozen model hash, and preprocessing version; each pair of float32
4,560-vectors costs36,480bytes before metadata. Cache extraction and training
are separate stages. Do not materialize all one-hot sequences in RAM.

CPU interval construction uses the existing compute-node allocation only,
streaming reference sequence and interval indexes. Parallelize chromosome
audits and later independent model/seed fits; do not permit concurrent writes
to a cache. Atomically finalize immutable caches before dependent fits. Retain
CPU/RAM/time measurements for this preflight and observed accounting for later
jobs. A scheduler template is a plan, not permission to submit it.

Current request syntax and allocation constraints are described in the
[official BioWulf user guide](https://hpc.nih.gov/docs/userguide.html).
Exact future resource availability must be rechecked when execution is approved.

## 15. Expected outputs and freeze

Deliverables include full source/peak provenance, master interval/sequence
features, exact per-configuration/model/partition BEDs and manifests, chr7 role
assignments, common challenge panels, duplicates/overlap audits, matching pairs
and diagnostics, class/ratio tables, configuration overlaps,18-row run matrix,
source/code/data SHA-256 ledgers, CPU orientation/threshold contract tests,
independent validation, the pretraining specification, and a freeze manifest.

The freeze records exact seeds, augmentation and inference, checkpoint rule,
threshold rule, selection rule, source/input identities and every interval set.
No post-training edits within this version are permitted. The freeze describes
the constructed experimental version; it does not by itself authorize training
or waive any failed readiness gate.

## 16. Failure and stop conditions

Stop before training for input-hash drift, benchmark outcome access, missing
phase-I weights, unreconciled V1 membership, unavailable/incorrect-length
sequence, B/C positive-mark negative contamination, cross-partition genomic or
encoded-sequence leakage, inconsistent same-lobe provenance, inadequate
matching/common support, insufficient independent validation/calibration data,
or undefined construction semantics. Do not fall back to random phase-I
weights if a weight file is missing. Ambiguous bases and known frozen A defects
must be explicit rather than silently presented as clean.

Later stop conditions include any real-network invariance failure, incomplete
three-seed ensemble, major seed instability, inadequate common-panel performance,
ambiguous internal selection, or unexpected compute/resource failure. Fixes that
change frozen data/rules require a new experimental version. Test/external
performance never grants permission to revise this version.

## 17. Decisions requiring investigator approval

The mandatory decision is whether to approve this frozen preflight and authorize
the proposed experiment **after its matching blockers are resolved**. The current
freeze must not be approved directly for GPU training. No job starts automatically.
The final readiness/QC
record identifies any additional actual data limitations that require a choice;
an unresolved hard gate cannot be waived implicitly. Proposed0.02 performance
margins and seed-stability tolerances are explicitly part of this pretraining
review, not claimed as uniquely correct scientific constants. A new architecture,
alternative lobe exclusion, safety-decontaminated A, or allele-effect calibration
is conditional and separately authorized, not part of the current run matrix.

Two actual design choices now need review:

1. **B common support:** keep the full V1 positive population with explicitly
   residual GC imbalance, or authorize a revised ratio/matching/common-support
   construction. I recommend further data-design work, not training the current
   imbalanced set. A common-support positive restriction changes the estimand
   and breaks the simple exact-V1-positive B contrast; if selected, apply and
   audit the corresponding restriction in a conditional matched A reference
   rather than quietly attributing every change to controls.
2. **C controls:** identical B/C controls isolate the label change but are not
   balanced to the same-lobe positive population. I recommend authorizing a
   C-specific rematching design in a new pretraining version, explicitly calling
   its contrast the label change plus induced control-rematching effect. The
   present fixed-control contrast is a frozen diagnostic, not proof of adequate
   accessible matching. No additional configuration is launched automatically.

These are genuine population/contrast choices, not a request to choose a seed,
architecture, or threshold after outcomes. Whichever correction is approved
must produce a replacement audited data/rule freeze **before** a separate GPU
authorization. Merely loosening the SMD gates until the current data pass is
not recommended.

## Constructed datasets and quantitative readiness findings

Counts below are positive/control, before any training. Validation counts include
all three disjoint chr7 roles; their exact IDs and per-role counts are frozen.

| Model | Configuration | Train | chr7 validation | chr8–9 test | All partitions |
|---|---|---:|---:|---:|---:|
| Enhancer | A | 156,789 / 307,473 | 8,460 / 20,328 | 16,103 / 34,903 | 181,352 / 362,704 |
| Enhancer | B | 156,789 / 146,031 | 8,460 / 8,259 | 16,103 / 15,093 | 181,352 / 169,383 |
| Enhancer | C | 125,874 / 146,031 | 6,631 / 8,259 | 12,695 / 15,093 | 145,200 / 169,383 |
| H3K27me3 | A | 26,326 / 51,839 | 1,612 / 3,427 | 2,646 / 5,902 | 30,584 / 61,168 |
| H3K27me3 | B | 26,326 / 26,162 | 1,612 / 1,599 | 2,646 / 2,644 | 30,584 / 30,405 |
| H3K27me3 | C | 17,132 / 26,162 | 1,096 / 1,599 | 1,754 / 2,644 | 19,982 / 30,405 |

Same-lobe construction removes 36,152 enhancer positives (19.935%) and 10,602
H3K27me3 positives (34.665%). This is materially different label construction,
not a negligible deduplication. The removed windows remain in the audit with
their original V1 labels; they are not recast as negatives.

V1 donor-accessible controls at the 1-kb core: enhancer 44,639/362,704 (12.307%);
H3K27me3 9,490/61,168 (15.515%). At the 2,001-bp input these become 67,244
(18.540%) and 13,893 (22.713%). At the original summit they are 17,049 (4.701%)
and 3,787 (6.191%). These are explicitly different overlap definitions.
V1 mark-contaminated expanded control windows number 5,966 enhancer and 1,041
H3K27me3, despite zero corresponding-mark overlap at the original core.

The new eligible pools contain 286,943 enhancer and 446,299 H3K27me3 controls,
**all donor-ATAC anchored by construction**. Matching selects 169,383 and 30,405
unique controls, respectively: 93.400% and 99.415% of the 1:1 targets.
Unmatched positive targets are 11,969 and 179; none was silently removed.
Every selected B/C negative passes the full-input same-mark exclusion.
All 1,143,370 master intervals have available 2,001-bp sequences. No exact or
model-encoded RC-equivalent sequence group crosses chromosome partitions.
Nine selected encoded-identity groups (18 genomic rows) occur within the
enhancer B/C training partition only; they remain explicit dependent genomic
observations, not independent validation replicates. Exact repeated genomic
coordinates from multiple lobes do not receive duplicate rows.

**Accessible provenance is clean; covariate matching is not fully clean.** There
are 18 prespecified B/C balance failures. B enhancer training GC SMD is 0.116
(limit 0.10). C enhancer training GC/repeat SMDs are 0.180/−0.123; C H3K27me3
training GC/repeat SMDs are 0.178/−0.111. C H3K27me3 lower-left and upper-right
lobe-support proportion gaps are +0.095 and −0.105 (limit 0.05). Additional
held-out construction imbalance is tabulated, not used as model-performance
feedback. No model scores have been generated. A diagnostic ATAC-strength SMD
uses only exact source-anchor matches and must not be interpreted as an
all-control estimate; the explicit observed/missing denominator table resolves
that limitation. Interval-overlap fractions above use every V1 control.

The common selection panels contain 1,755 enhancer positives/2,171 controls in
1,405 components and 271 H3K27me3 positives/379 controls in 303 components.
Common calibration has 1,576/1,944 enhancer and 263/391 H3K27me3 examples. These
meet the prespecified count minima but do not override the matching failures.

Machine-readable outcomes are in `results/COPD-V2-PREFLIGHT_*` and
`provenance/COPD-V2-PREFLIGHT_final_validation_manifest.json`; all constructed
intervals, BEDs, roles, matching pairs, unmatched targets and the 18-run matrix
are under `data/COPD-V2-PREFLIGHT/`. A complete artifact inventory and SHA-256
ledger accompany the final freeze. Independent artifact-validation success
means these records faithfully describe the constructed data, **not** that
the failed matching-readiness gates passed.

## Methodological references

[Zhou, Shrikumar and Kundaje (2022)](https://proceedings.mlr.press/v165/zhou22a.html)
distinguish RC augmentation from post-hoc conjoined inference; the invariant
averaging rule here is also established directly by the algebra in section6.
[Average precision definition](https://scikit-learn.org/stable/modules/generated/sklearn.metrics.average_precision_score.html)
and [validation-monitored early stopping](https://www.tensorflow.org/api_docs/python/tf/keras/callbacks/EarlyStopping)
specify metric/callback semantics. Frozen V1 code and independently audited
source provenance are the primary evidence for this particular experiment.
