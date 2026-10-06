# COPD V2 internal-training-1.0: prospective execution contract

This stage implements the explicit investigator authorization attached as
`933794dc-e228-4aa9-a4d9-db8cf11e31e7/pasted-text.txt`. The immutable scientific
contract is pretraining-1.1, published at repository commit
`7c0caf4af5f8f346089c92f43b88161331ca420b`. Its historical prohibition on training
is superseded only by this new authorization; no frozen artifact is edited.

## Scope and gates

Only training chromosomes (1–6, 10–22, X, Y) and chr7 are eligible for sequence
fetch, feature extraction, or inference. Reading frozen interval metadata to
filter its train/validation rows does not authorize extraction of chr8–9.
Benchmark files remain opaque; no benchmark, candidate, test, or manuscript
analysis is performed. Frozen 1.0 and 1.1 payloads and all consumed-source
hashes must pass the independent integrity gate before execution.

1. Prepare exact frozen-order training and native chr7 checkpoint tables;
   verify raw/canonical sequence hashes, 1000-bp core and 2001-bp geometry.
2. Independently evaluate both nucleotide orientations through the strict,
   loaded frozen phase-I network. Freeze both float32 cache matrices and hashes.
3. Execute precisely the 18 prescribed phase-II fits (A/B/C × two contexts ×
   seeds 104729, 130363, 155921). Training workers open only train and native
   checkpoint labels; no selection/calibration performance is inspected.
4. Freeze all 18 selected checkpoints before any common-panel inference.
   Run actual-network orientation-swapped passes and seed/ensemble invariance.
   Any numeric failure stops this stage.
5. Evaluate only the common C-task chr7 selection panels. Apply all absolute
   C adequacy and seed-stability gates independently to both contexts. A/B
   comparisons are diagnostic and cannot select the final model.
6. Only if BOTH C contexts PASS, freeze the six C checkpoints, inference code,
   equal-seed ensemble and adequacy decision; then calibrate C region-score
   thresholds using dedicated chr7 calibration negatives alone. Otherwise STOP.
7. Complete independent QC, report and freeze this internal stage. No chr8–9
   or external benchmark access follows automatically, even after a PASS.

## Prospective numerical implementation details

Every non-ACGT nucleotide maps to the zero-encoded N symbol. Cache rows dedupe
the frozen canonical encoded forward/RC identity. The canonical nucleotide
sequence and its nucleotide RC each independently enter phase I; a mapping
restores each interval's forward identity. No 4560-feature axis is reversed.
Inference uses batch size 32 for phase I and 256 for phase II, float32 network
operations, and float64 orientation/seed aggregation in frozen seed order.

The exact frozen V1 phase-II architecture is independently initialized for
each seed. Unweighted BCE, Adadelta(0.001, rho=0.95, epsilon=1e-7), batch256,
max50 epochs, patience15 and min_delta0 remain unchanged. Epoch indices are
zero-based. Frozen substream1 defines the permutation; substream2 uniform
draws <0.5 assign one RC view per interval in frozen genomic order BEFORE
applying the permutation. No replacement seeds or seed weighting are allowed.

Native checkpoint BCE is the float64 mean Bernoulli log loss of symmetric
probabilities, clipped to [1e-7,1-1e-7] using the original Keras epsilon.
Strict numerical decrease alone replaces the best checkpoint, so ties keep
the earliest epoch. Selected full trainable and BatchNorm weights are restored,
and restored BCE must match exactly. Checkpoints are immutable inference-only
Keras files, not resumable optimizer snapshots.

The real-network gate covers every sequence in the union of native chr7 and
common chr7 panels, for each model/configuration/seed and equal-seed ensemble.
Independent second network calls swap nucleotide-orientation cache inputs;
they do not reuse an algebraically commuted saved probability pair. Calibration
rows may participate in this label-free invariance gate, but calibration
performance/thresholds cannot be examined before both C adequacy decisions PASS.

Evaluation conventions, exact bootstrap implementation, fixed strata and
conditional calibration schema are locked in
`specification/internal_evaluation_implementation.json` before outcomes.
Bootstrap seed314159 is reset per context, with the same sampled components
across A/B/C in each replicate, exactly2000 valid replicates, max20000 attempts,
invalid fraction<0.10. No result-responsive design changes or margin searches.

## Runtime and infrastructure

The original venv's interpreter points inside another user's private home.
Its accessible pyvenv.cfg currently records Python3.13.13, whereas the frozen
training manifest and architecture audit prescribe Python3.13.0. Execution
therefore uses a newly accessible uv-installed CPython3.13.0 interpreter with
the original, read-only site-packages via explicit PYTHONPATH. Verified versions:
TensorFlow2.20.0, Keras3.14.1, NumPy2.5.0. No original environment file is changed.
`runtime.sh` records the deterministic settings and shared NVIDIA library paths;
environment manifests capture executed interpreter/package metadata hashes.

Phase-I extraction uses the already allocated idle A100 on the interactive
allocation (job31881465), bounded to 8 compute threads. Training is a Slurm
array of exactly18 elements with at most4 concurrent fits, each one A100,
8 CPUs, 48 GiB RAM and 4-hour walltime, matching the frozen compute plan.
No ad-hoc pilot fit is added. A failed run is never replaced for performance;
only an explicitly documented identical same-seed infrastructure retry is legal.
Scientific/sequence/hash/numeric failures stop and require reporting.

All substantive commands have exclusive started/completed records in
`provenance/commands/`; the earlier read-only diagnostics and CPython access
restoration are described in runtime provenance. Outputs, input lineage, code
and checkpoints are hashed at the relevant pre-execution and final gates.

This is region-label modeling only. Neither a region-classification improvement
nor a calibrated region threshold establishes REF–ALT effect validity or a
causal variant. Investigator review is required before opening chr8–9.
