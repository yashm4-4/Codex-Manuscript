# Frozen external evaluation data: mandatory future-use boundary

Benchmark version: 1.0. Frozen 2026-10-05 23:31:15 UTC.

Master: `results/COPD-V2-BENCH-R003_frozen_benchmark_master.tsv`, relative to
the V2 workspace. SHA-256:
`a69a8b1d5c82299075b2d5203fee3d4a23105572c6cf456bab834b82ef13bbdd`.
The authoritative file/label/source boundary is
`COPD-V2-BENCH_benchmark_freeze.json` in this directory.

This is **FROZEN EXTERNAL EVALUATION DATA**, not a training, calibration or
model-selection resource. The current comparison describes already frozen V1
models; it does not make V1 an untouched test of model selection.

Future V2 development may know that the benchmark exists, but:

- Do not use its labels for training.
- Do not use its outcomes to select architecture, hyperparameters or thresholds.
- Do not use its outcomes to construct controls or training labels.
- Do not use its outcomes to select seeds, orientation handling, cell models,
  checkpoints, or which run is retained.
- Select models exclusively with prespecified internal validation/calibration
  criteria, independent of this benchmark's results.
- Freeze the future model design and selection procedure before final external
  evaluation. Record that freeze and evaluate the prespecified retained model.
- Preserve all source rows, including nulls, missing labels, ambiguous results,
  context conflicts and identity exclusions. Never add/remove cases to improve
  external performance.

Experimental assay orientation and computational sequence orientation are
distinct. Neither benchmark labels nor outcomes choose a correct orientation.
Do not average forward/RC scores or invent a benchmark-specific call rule.

The current set has incomplete labeled/QC screening denominators and incomplete
existing-V1 sequence coverage. It supports typed descriptive evaluation, not
unqualified sensitivity, specificity, AUROC or population error estimates.
Future improved access must be documented as a new evidence version before
evaluation, never a silent mutation of this frozen version.

No model development, next V2 module, manuscript revision or publication action
is authorized by completion of this benchmark.
