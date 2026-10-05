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
