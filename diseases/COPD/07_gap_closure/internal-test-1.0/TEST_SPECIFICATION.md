# COPD V2 internal-test-1.0: prospective chr8–9 holdout contract

This specification is recorded before any chr8–9 model inference. The complete
machine-readable contract is `specification/test_specification.json`, with
its SHA-256 and this document's SHA-256 in
`provenance/prospective_specification_freeze.json`. The investigator request is
preserved in `provenance/user_authorization.txt`. No test-performance outcome
may alter this specification.

## Interpretation and firewall

Chr8–9 was held out from V2 training, checkpoint selection and chr7 calibration.
It is **not historically pristine**: V1 used these chromosomes and V2
construction/QC inspected outcome-blind covariates. This is a fixed internal
V2 region-label holdout, **not independent external validation**.

The published `internal-training-1.0`, both pretraining stages and all other
historical artifacts remain immutable. The external functional benchmark must
remain unopened and unparsed, including labels and outcomes. No COPD variant
universe, 337 candidates, 12 shortlist variants, rs2013701 or other special-case
variants, or arbitrary REF/ALT sequences may be scored. There is no training,
checkpoint replacement, seed selection, threshold calibration, matching change,
fine-mapping, target-gene analysis, manuscript revision, or automatic publication.

## Fixed inputs and inference

Use only the frozen pretraining-1.1 common C-task test rows, preserving all
membership, labels, roles and existing global component IDs. Enhancer has
12,695 positives and 12,695 controls in 8,754 components; H3K27me3-associated
has 1,754 positives and 1,754 controls in 1,809 components. Their union is
26,225 unique intervals/encoded identities. Exact source hashes and counts are
locked in the machine contract. Any mismatch is a hard stop, never permission
to substitute rows or reconstruct matching.

Verify the original hg38 FASTA/FAI, frozen raw and canonical sequence hashes,
1,000-bp core and 2,001-bp input geometry. Normalize non-ACGT to zero-encoded N;
independently evaluate canonical nucleotide and nucleotide-RC one-hot arrays
through the unchanged phase-I model to obtain float32 [26225,4560] caches.
Never reverse representation columns. Phase-I batch size is 32.

Load all 18 exact original selected `.keras` archives with `compile=False`
after frozen SHA-256 verification. Use configurations A/B/C, both contexts,
and fixed seeds 104729, 130363, 155921. Phase-II batch size is 256;
network operations remain float32. Each seed score is the float64 mean of
forward and nucleotide-RC probabilities; the ensemble is the fixed equal-three-
seed float64 mean in that seed order. No archive is altered or re-saved.

For every test sequence and every seed and ensemble, make independent repeated
actual-network calls with swapped orientation inputs. The frozen RC gate is
`abs(q(x)-q(RC(x))) <= 1e-6 + 1e-6*abs(q(RC(x)))`. Check finite [0,1]
probabilities at every raw/repeated/symmetric/ensemble step. Any failure stops.

## One-time evaluation and release gate

Evaluate all seeds and fixed ensembles on identical C-task panels: stepwise
tie-aware AP, tie-aware AUROC, Brier and Brier skill; retain ten fixed reliability
bins [0,.1), …, [.9,1]. Report 95% confidence intervals using the inherited
component bootstrap: lexical component order, fresh NumPy default_rng(314159)
per context, sample the observed number of components with replacement, carry
all rows with component multiplicities, exactly 2,000 valid draws, at most
20,000 attempts, invalid/attempted strictly below 0.10. All A/B/C predictors
share draws within each context. Percentiles use linear interpolation.
AP-minus-prevalence and Brier skill use each replicate's own prevalence.
Report paired B−A, C−B and C−A AP/AUROC/Brier differences. These are diagnostic
ablations; C need not outperform B and model choice cannot change.

Separately for each C context, PASS requires all three seeds evaluable,
real-network RC PASS, finite valid probabilities, exact sample/component counts,
2,000 valid bootstrap draws with invalid fraction below 0.10, lower 95% AUROC
bound above 0.5, lower 95% AP-minus-replicate-prevalence above zero, and lower
95% Brier-skill bound above zero. Three-seed AP sample SD (`ddof=1`) must be
at most 0.03 and range at most 0.10. All conditions are conjunctive. Missing
evidence is INCONCLUSIVE; observed failed conditions are FAIL. Either non-PASS
stops downstream work for investigator review, without redesign or retraining.
These are internal generalization safeguards, not clinical effect-size criteria.

## Frozen thresholds and descriptive diagnostics

Apply `score >= threshold` to C ensembles only, without recalibration:

- Enhancer: `0.74848511815071117` (`0x1.7f39710000001p-1`).
- H3K27me3-associated: `0.76960810025533055` (`0x1.8a0a12aaaaaacp-1`).

Report TP/FP/TN/FN, sensitivity, specificity, precision, NPV, accuracy and
observed negative-row FPR. Rate CIs use the same component draws at fixed
thresholds; undefined ratio denominators have explicit omitted/valid counts.
A test FPR above 5% is calibration drift, not automatic protocol failure.
These are region-label thresholds, not variant FPR control, allele-effect
thresholds or causal probabilities.

Reuse the frozen training-derived GC, ATAC-rank and repeat quintile cutpoints
verbatim, with searchsorted(side=right). Use only already-specified exact ATAC
lobe signatures, context-specific same-lobe support signatures, chromosome and
non-ACGT status. Flag n<100 cells; undefined one-class metrics remain empty.
All subgroup results are descriptive and cannot define new subgroup searches.

## Execution, validation and stop

Use the frozen CPython 3.13.0 / TensorFlow 2.20.0 / Keras 3.14.1 / NumPy 2.5.0
environment and deterministic settings. Read-only helper imports must suppress
bytecode writes. Hash-bound gates precede sequence preparation, phase-I
extraction, checkpoint inference and evaluation. Record commands, resources,
software, hashes and any failure; no outcome-responsive reruns or changes.

Produce all requested results, an independent validation, append-only V2
activity/decision/result entries, the detailed holdout report, and final
checksums/freeze. No Git commit/push is authorized in this request. Report both
C decisions and readiness for investigator review, then stop. A PASS never
authorizes external benchmark access or candidate/allele scoring.
