# internal-training-1.0 publication storage addendum

This directory records publication of the **unchanged** frozen internal training result to `yashm4-4/Codex-Manuscript`, branch `main`. It is outside the frozen stage and authorizes no new scientific computation.

## Storage classes

| Class | Stage files | Original content bytes | Publication policy |
|---|---:|---:|---|
| Ordinary Git artifacts | 451 | 295,006,279 | All 449 ordinary frozen payloads plus the original checksum ledger and freeze record |
| Git LFS original checkpoints | 18 | 6,301,157,652 | Every selected original `.keras` archive, unchanged |
| External reconstructable intermediates | 3 | 27,673,496,091 | Retained at recorded managed-storage paths; neither Git nor LFS |

The Git-hosted stage comprises **469 files and 6,596,163,931 original content bytes**. This excludes publication metadata, the three shared registers, root storage-policy files, Git history/pack overhead and LFS pointer overhead. These are original-content accounting totals, not network-transfer estimates.

[storage_inventory.tsv](storage_inventory.tsv) inventories all 472 stage files, the three append-only V2 registers and the two root storage-policy files: **477 rows**, with exact repository-relative path, bytes, SHA-256 and storage class. [storage_manifest.json](storage_manifest.json) records reconstruction and preservation requirements. Publication-only metadata is bound by Git commits rather than creating a circular self-hash manifest.

## Preserve all 18 original checkpoints

The six V2-C archives are retained-model ensemble members; the twelve V2-A/V2-B archives preserve the actual ablation models. All are scientific-result artifacts and all are assigned to LFS. Their SHA-256 values are both the original frozen content hashes and their LFS object IDs.

Do **not** retrain, re-save, slim, strip optimizer state or repackage these archives. Existing packaging validation shows valid selected inference weights plus fresh ancillary Adadelta state; inference loads with `compile=False`. They are **not resumable training snapshots**. Their archive hashes remain the identity of the preserved result.

The LFS attributes enumerate these 18 paths specifically. The existing unrelated LFS rule is preserved. [preflight.json](preflight.json) records installation/configuration, repository access and known storage-blocker checks; an authorized upload response is not a guarantee of unused account quota. A quota/storage rejection is a hard publication stop, never permission to exclude or modify checkpoints.

Post-publication verification **passed** for artifact commit `276a7881ae970fbdb19c78fd17dd1453c8b6f58f`. [remote_verification.json](remote_verification.json) records fresh downloads of all 18 remote LFS objects into a previously empty isolated object store: every downloaded size and SHA-256 matches its frozen local original. A second read-only audit independently confirmed all 18 downloaded hashes. All 470 original scientific payloads were rehashed unchanged; the 451 ordinary stage files, three excluded intermediates and historical artifact preservation checks also passed.

The checkpoint upload and fresh download both succeeded without quota/storage rejection. Billing usage visibility was unavailable (HTTP 404), so no claim is made about remaining account quota. GitHub accepted the ordinary Git payload with non-blocking recommended-size warnings for `inputs/cache_index.tsv.gz` and `inputs/configurations/V2-A_enhancer_train.tsv.gz`; neither file was changed. The receipt records exact original-content payload accounting at the artifact commit. Its later containing commit adds only publication verification metadata/documentation. [verify_publication.py](verify_publication.py) is the read-only byte/hash/tree verification helper; it does not load models or perform scientific analysis.

## Three external reconstructable intermediates

These are **not part of the Git-hosted frozen payload** and are explicitly ignored, not LFS-tracked:

| Stage-relative path | Bytes | Schema |
|---|---:|---|
| `cache/features_canonical.npy` | 13,117,241,408 | float32 [719147, 4560], canonical nucleotide phase-I representations |
| `cache/features_rc.npy` | 13,117,241,408 | float32 [719147, 4560], independently evaluated nucleotide-RC representations |
| `inputs/canonical_sequences.npy` | 1,439,013,275 | uint8 ASCII [719147, 2001], canonicalized normalized nucleotide strings |

All remain in place beneath:

```text
/vf/users/Dcode/yash/Codex-Manuscript/disease-regulatory-genomics-workflow/diseases/COPD/07_gap_closure/internal-training-1.0/
```

Their exact individual absolute/resolved paths, sizes, SHA-256 values, ordering/index mappings, reconstruction inputs and script/specification hashes are in the storage manifest. Merely recording a filesystem path is not a new backup or retention guarantee.

Sequence reconstruction requires the exact frozen hg38 FASTA/FAI and pretraining-1.1 interval/common-panel manifests. Phase-I feature reconstruction additionally requires the original phase-I weights, canonical sequences, architecture, runtime and gates. Those upstream inputs must remain available as actual files: **hashes are not backups**. The manifest identifies the existing managed reference/model realpaths and frozen hashes without altering or republishing them.

The original complete input-preparation step took **174.34796173404902 CPU wall seconds**. Both orientation caches together took **393.5700578056276 A100 wall seconds**, with 28,448,564 KiB peak host RSS. The paired extraction cost is shared and must not be counted once per cache. Costs are historical observations, not guarantees.

The recorded scientific runtime is CPython 3.13.0, TensorFlow 2.20.0, Keras 3.14.1 and NumPy 2.5.0; exact package/build/environment and source hashes are linked in the manifest. The procedure is preserved, but independent byte-identical numerical regeneration was not demonstrated. Canonical sequence extraction has deterministic rules; exact NPY bytes additionally depend on the writer format. **Do not derive RC features by reversing feature columns.** No reconstruction was performed for publication.

Frozen scripts reject replacement of existing outputs. Any future authorized reconstruction must use a separately prepared output location and the required gates/dependencies; do not run reconstruction against this immutable stage.

## Obtaining the published model files

For a separate checkout, ordinary Git retrieval followed by an LFS pull is sufficient for the 18 preserved archives:

```bash
GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/yashm4-4/Codex-Manuscript.git
cd Codex-Manuscript
git lfs install --local
git lfs pull --include="diseases/COPD/07_gap_closure/internal-training-1.0/runs/*/attempt-001/selected_checkpoint.keras"
```

These are documentation-only commands, not analysis instructions. Obtain external intermediates separately only if full train/chr7 replay is required.

## Frozen versus published coverage

The original checksum ledger covers **470 scientific payload files**, including the three external intermediates. It and the original freeze record remain unchanged. A Git clone plus LFS does **not** claim to contain the entire full-bundle ledger: it contains 467 of those payloads plus the two freeze metadata files, while the three excluded files are explicitly accounted for here.

The original freeze's `commit_or_push_performed:false` is a historical statement about the earlier scientific freeze. It is intentionally unchanged; this separate addendum records the later publication authorization.

V1, the external functional benchmark, pretraining-1.0 and pretraining-1.1 must remain unchanged. Shared V2 registers retain their already-completed append-only internal-training entries. This publication does not authorize chr8–9 execution, benchmark evaluation, candidate scoring, retraining or another V2 module.
