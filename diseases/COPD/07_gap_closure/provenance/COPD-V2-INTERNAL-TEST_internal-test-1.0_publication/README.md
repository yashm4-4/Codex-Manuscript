# Frozen internal-test-1.0 publication storage addendum

This is an administrative storage/publication addendum outside the immutable scientific stage. It does not change the original report, specifications, results, checksum ledger, freeze record, prior stages or model archives, and does not authorize any new scientific execution.

The target is `yashm4-4/Codex-Manuscript`, branch `main`, from historical baseline `07d72dad0f57d605f5f50bbc7fa6d2c7843a94ba`. The original [scientific freeze](../../internal-test-1.0/provenance/freeze.json) and [checksum ledger](../../internal-test-1.0/provenance/artifact_checksums.tsv) remain byte-for-byte unchanged. Their historical statements that no commit/push had yet occurred describe the time of scientific freezing; publication does not rewrite them.

## Complete storage accounting

| Storage class | Count | Exact logical bytes | Treatment |
|---|---:|---:|---|
| Git ordinary: frozen scientific payloads | 170 | 136,520,145 | Publish unchanged |
| Git ordinary: original checksum ledger and freeze record | 2 | 22,122 | Publish unchanged |
| External reconstructable phase-I test intermediates | 2 | 956,688,256 | Preserve at existing managed paths; no Git or LFS |
| New Git LFS checkpoint objects | 0 | 0 | No checkpoint duplication |

The scientific freeze contains **172 payloads totaling 1,093,208,401 bytes**: 170 Git-hosted payloads plus two external caches. Including the original ledger and freeze record, the Git-hosted stage contains **172 files totaling 136,542,267 bytes**. These are different counts with intentionally different scopes.

The three completed append-only V2 register files are also published unchanged from their frozen update records (70,126 full-file bytes; 8,665 newly appended bytes). Their existing historical prefixes are preserved. Publication policy additions and this addendum are outside the scientific freeze and outside the stage-byte totals.

[storage_inventory.tsv](storage_inventory.tsv) lists all 172 frozen scientific payloads, the original ledger/freeze, and the three register files (177 records), with exact repository paths, byte sizes, SHA-256 values, storage class and managed-storage path. [storage_manifest.json](storage_manifest.json) contains the full machine-readable policy, dependencies, environments, command provenance and reconstruction contract. Publication-only metadata and root policy changes are bound by the Git commit, rather than introducing recursive self-checksums.

All report/specification files, test panels, sequence input, provenance, predictions, RC audits, A/B/C and bootstrap results, comparisons, holdout-gate results, frozen-threshold diagnostics, environment/compute records, independent validation and scripts remain ordinary Git artifacts. No scientific result is removed from the original checksum ledger.

## Oversized-artifact classification

1. **Scientifically essential non-regenerable artifacts:** no new oversized artifact in this stage. The 18 original selected `.keras` archives are already preserved by Git LFS in `internal-training-1.0`; their existing paths/objects are reused, unchanged and without copies in this stage.
2. **Reconstructable caches/intermediates:** exactly the two arrays below. They are explicitly **external reconstructable intermediates, not part of the Git-hosted frozen payload**.
3. **Other artifacts requiring a storage decision:** none.

Both external arrays are NumPy NPY v1.0 files, little-endian float32 (`<f4`), C order, shape **[26,225, 4,560]**, with 128-byte headers. Each is **478,344,128 bytes**.

| Stage-relative path | SHA-256 |
|---|---|
| `cache/features_canonical.npy` | `6681952719517712d123cd753c388b04049b88251ac5c434875ad4a98d32ed6e` |
| `cache/features_rc.npy` | `2b99b0a2509c07af7ba5141ec7f67ca49b34c97f78285d3be5a190376985dcaa` |

They remain at these exact managed shared-filesystem locations:

- `/vf/users/Dcode/yash/Codex-Manuscript/disease-regulatory-genomics-workflow/diseases/COPD/07_gap_closure/internal-test-1.0/cache/features_canonical.npy`
- `/vf/users/Dcode/yash/Codex-Manuscript/disease-regulatory-genomics-workflow/diseases/COPD/07_gap_closure/internal-test-1.0/cache/features_rc.npy`

Neither file is added to Git LFS, moved, transformed or regenerated. Documentation of current managed storage is not a new backup, replication or retention guarantee.

The canonical array contains the original frozen phase-I network's final 4,560 sigmoid activations. The RC array contains a separate network evaluation of nucleotide reverse-complement inputs, not a reversal of the feature axis. Rows are bound to the unchanged `inputs/cache_index.tsv.gz` identity order and the frozen union of the common chr8–9 test panels.

The **52,476,353-byte** `inputs/canonical_sequences.npy` is retained in ordinary Git, unchanged: uint8 ASCII A/C/G/T/N, shape [26,225, 2,001], SHA-256 `589ce6195f339aec072a7f20453cb13071d9d37a5bbf3f3705aa763cac616975`. It is not an excluded intermediate.

## Reconstruction contract: documentation only

**No reconstruction is performed or authorized by publication.** Future cache recovery must occur only under separate authorization in a new scratch output stage, never in the frozen stage.

The full procedural and dependency records are in `shared_cache_reconstruction` in [storage_manifest.json](storage_manifest.json), backed by the unchanged [input manifest](../../internal-test-1.0/inputs/input_manifest.json), [phase-I cache manifest](../../internal-test-1.0/cache/cache_manifest.json), [execution source manifest](../../internal-test-1.0/provenance/execution_source_manifest.json) and [phase-I execution gate](../../internal-test-1.0/provenance/phase_I_execution_gate.json).

Direct recovery requires the published canonical-sequence array and row index, frozen stage-side input/prospective gate records and specifications, original `scripts/extract_test_phase_one.py`, `scripts/test_contract.py`, `scripts/runtime.sh`, and the already-published immutable `internal-training-1.0/scripts/phase_two_contract.py`. The stage-side gate enumerates all checked records. The reconstruction does **not** require any of the 18 selected phase-II checkpoints.

The original phase-I weights must remain available unchanged:

- Repository dependency: `diseases/COPD/04_modeling/trednet/model_phase_I/phase_one_weights.h5`.
- Resolved managed-storage path: `/vf/users/Dcode/gaetano/projects/Writing_Scientific_Manuscript_Codex/models/TREDNET_v2/model_phase_I/phase_one_weights.h5`.
- Size: **565,399,080 bytes**.
- SHA-256: `483b6c0cafd750ea8e97932cbc8a949276eaf1e77e7d8766f9a22eeddc2916d6`.

The JSON also records the exact hg38 FASTA/FAI and pretraining-1.1 common C panel/interval manifests used for the original sequence construction. They are not needed for direct cache recovery while the Git-hosted sequence input remains available.

The preserved extractor requires a previously nonexistent `cache` output directory, verifies frozen inputs and weights, uses batch size 32 and two independent nucleotide-orientation calls with `training=False`, and writes both arrays in unchanged row order. For separately authorized recovery, prepare an isolated stage holding the exact gate-referenced inputs and provenance but no cache directory; use the preserved extractor/runtime with `--repo` pointing to the frozen repository and `--stage` pointing only to that new stage. The original command's exact argv is preserved in `provenance/commands/phase_I_test_forward_and_RC.started.json` and the publication manifest. Do not run that historical argv unchanged against the original stage.

The strict original runtime is CPython **3.13.0**, TensorFlow **2.20.0**, Keras **3.14.1**, NumPy **2.5.0**, float32 policy, and exactly one allocated GPU. Recorded hardware was an NVIDIA **A100-SXM4-80GB**, driver **580.173.02**, TensorFlow CUDA build **12.5.1** / cuDNN **9**. Set `PYTHONHASHSEED=104729` before interpreter startup; retain deterministic TensorFlow/cuDNN operations, disabled XLA JIT, and the original thread/runtime/library settings. The JSON retains the entire recorded environment and package versions, as well as the managed interpreter and site-packages locations.

Recorded cost for **both caches jointly** was **26.03783931955695 seconds** of extractor wall time (logged wrapper **28.86424876190722 seconds**) on one A100; peak host RSS was **2,395,568 KiB**. These are original observations, not a new regeneration measurement or cost guarantee. Optional original sequence preparation took **15.592166801914573 seconds**, peak RSS **358,512 KiB**. Queue time, transfers, setup and additional source verification are not included.

Reconstructability is conditional on continued availability of exact inputs, source weights and a compatible exact software/hardware stack. Deterministic settings were used, but independent same-byte regeneration has **not** been demonstrated. Compare any future recovered files against the frozen sizes, headers and SHA-256 values; do not accept a mismatch as equivalent or rewrite the freeze.

## Publication verification and boundaries

The publication workflow independently verifies the remote commit, all Git-hosted frozen payload hashes, exact external-cache exclusions and accounting, original local frozen hashes, append-only register prefixes, and unchanged historical stage/checkpoint references. A separate `remote_verification.json` may record completed publication checks in this addendum directory; this document alone does not assert that publication has succeeded.

No inference, feature extraction, metric recomputation, recalibration, retraining, external benchmark outcome access, COPD variant/candidate scoring or other V2 analysis is performed as part of publication. Earlier V1, benchmark, pretraining and internal-training stages remain unchanged. The scientific stage remains frozen.
