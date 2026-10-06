# V2 scripts

`run_phenotype_robustness.py` reads frozen V1 tables and reviewed V2 JSON inputs,
then reconstructs exact significant study/tag/panel support. It is offline,
standard-library-only Python and writes exclusively within V2.

`validate_phenotype_module.py` verifies original V1 freeze hashes, tracked-file
immutability, input/output checksums and per-record/summary agreement; it then
writes final validation and artifact checksums. Run it after any authorized
V2 report/register edits to refresh the final artifact manifest.

Run both from the repository root as shown in the report. Neither script
trains models, scores alleles, recalculates LD, maps targets or reranks V1.

## Frozen-model reverse-complement audit

- `run_reverse_complement_scoring.py`: locked preflight, all-sequence checks,
  deterministic forward verification, then exact frozen-model RC inference.
  It refuses to overwrite an existing RC score file.
- `analyze_reverse_complement_audit.py`: original-threshold comparisons,
  rank-preserved subsets and prespecified quantitative/severity tables.
- `plot_reverse_complement_audit.py`: four standalone PNG/PDF figure families.
- `validate_reverse_complement_audit.py`: independent source-based result
  reconstruction, V1/spec integrity checks and final artifact checksums.

The RC report records the compatible Python/library invocation. The first
verification attempt exposed partial-batch numerical dependence; the corrected
verification pads execution to full minibatches while retaining exactly the
same 142 assessed IDs and fixed tolerances. RC scoring preserves original V1
batch shapes. Every new output is confined to V2; no script retrains or modifies V1.

## External functional benchmark

- `initialize_external_benchmark.py`: one-time protocol/baseline lock; refuses overwrite.
- `collect_benchmark_castaldi.py`, `collect_benchmark_gong.py`,
  `collect_benchmark_mechanisms.py`: public retrieval and source-only evidence extraction.
- `assemble_external_benchmark.py`: authoritative allele mapping, construct/reference
  checks, source-only tables and benchmark freeze; refuses any operation after freezing.
- `compare_external_benchmark.py`: verifies frozen source hashes before exact joins to
  existing V1/RC scores; preserves original calls/ranks and descriptive denominators.
- `validate_external_benchmark.py`: independent source, score, call, summary, freeze,
  boundary and final-checksum validation; requires `--comparison-complete`.

Do not rerun source collectors into this frozen namespace. Post-freeze comparison
can be reproduced as documented in the report; it executes no model and creates
no new model sequences. Its history records implementation-only refinements.
Source/runtime package details are in collection and run manifests. Benchmark
labels and results must never be used for training or model selection.
