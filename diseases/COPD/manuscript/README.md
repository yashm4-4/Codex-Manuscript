# COPD manuscript bundle

This directory contains the audited manuscript, figures, plotted source data,
claim-level result map, software-version capture, and final workflow validation.

## Primary files

- `COPD_regulatory_genomics_manuscript.md`: submission-format scientific draft.
- `result_to_manuscript_map.tsv`: claim-to-primary-result traceability.
- `supplementary_file_manifest.tsv`: size and SHA-256 inventory for the static data and reproducibility bundle.
- `software_environment.tsv` and `software_version_capture.txt`: recorded and recovered software versions; these are not a complete environment lockfile.
- `figures/` and `figure_data/`: reproducible figures and their plotted values.
- `results/COPD-workflow-validation.tsv`: whole-workflow audit table.

The validation table and log are dynamic audit artifacts and are intentionally
excluded from the supplementary checksum manifest; including files rewritten by
the audit that validates the manifest would create a checksum cycle.

## Rebuild and audit

Run from the repository root:

```bash
python3 diseases/COPD/manuscript/scripts/build_figures.py
python3 diseases/COPD/manuscript/scripts/renumber_citations.py
python3 diseases/COPD/manuscript/scripts/build_supplementary_manifest.py
python3 diseases/COPD/manuscript/scripts/validate_complete_workflow.py --strict-final
```

The citation script should be run only after substantive prose edits. Rebuild the
supplementary manifest after changing any file listed in that manifest.

## Submission status

The computational workflow is complete. The manuscript is not ready for journal
submission until named authors approve the science and supply authorship,
affiliations, correspondence, CRediT roles, funding, competing interests,
acknowledgments, institutional ethics wording, journal-specific AI disclosure,
and a persistent public data-and-code archive.
