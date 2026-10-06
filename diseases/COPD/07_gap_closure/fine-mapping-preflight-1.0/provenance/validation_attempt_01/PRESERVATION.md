# Preserved first validation attempt

The 2026-10-06T20:57:21Z run completed 4,805 checks and reported four failures, one for each Kim smoking stratum. Every failure concerned an incorrectly generalized bookkeeping flag: the validator required a provider-resolution flag to be true for all seven accessible direct sources, although that flag applied only to the three BBJ sources resolved outside Catalog.

The actual Kim public-access fields and remaining LD/scale/missingness blockers were correct. No source hash, header, historical-content, register or reported-headline check failed. This was a prefreeze schema/test semantics mismatch, not an observed scientific incompatibility or evidence of completed locus QC.

The root agent clarified the column name to `catalog_missing_source_access_resolved_by_provider`. The revised validator checks public access and absence of stale source/schema blockers for all seven sources, and separately checks that the provider-override flag is true only for BBJ. No check was silently suppressed, and no source bytes or scientific results were altered.

These three files are preserved byte-for-byte from the failed run:

| File | SHA-256 |
|---|---|
| `validate_preflight.py` | `7ec74a00f3aa72e7e662310c892934953169e650e462ac33b1b1213e1df4a8ef` |
| `independent_validation.json` | `7067031852d2ad52c470beb81938c3bdac0e71c251fb18271cb29dcc4cb5c868` |
| `independent_validation.tsv` | `4582dd3ef412024122261268609e265fcf53088fdd0bb55e9d824fb71d59f7e8` |

The later canonical receipt is written separately to `provenance/independent_validation.json`; the final freeze must retain this first attempt as well.
