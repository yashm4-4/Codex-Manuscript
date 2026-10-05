# COPD V2 phenotype robustness: analysis specification

Specification ID: COPD-V2-PHENO-SPEC-001. Authorized 2026-10-05 by the investigator's attached phenotype-module instructions. This specification is written before phenotype-retention results are calculated.

## Scope and phenotype rules

Classify all 104 frozen core study accessions by their actual association question and diagnostic method using the Catalog record and primary publication. Primary eligibility means direct clinically or spirometrically defined COPD cases versus non-COPD controls. Moderate/severe case selection and ordinary smoking restriction/stratification remain eligible, with those attributes explicit. Secondary eligibility means ICD/PheCode or comparable EHR-defined COPD case-control. Keep mixed direct/EHR meta-analyses, quantitative/ML liability, within-case outcomes, subtype contrasts, interaction coefficients, and nonstandard gene/burden/CNV tests separately identified. A title alone cannot establish the tested phenotype or genetic-test unit. Unresolved definitions remain unresolved.

The publication reviewers will record classifications before consulting downstream retention counts. Multiple association tests within one accession must be disclosed. If no single eligible phenotype can be justified, exclude that accession from the primary estimate and show an explicit bounded sensitivity for an arguable alternative. Do not retrospectively loosen the primary rule to retain a desired candidate.

Missing states are `not_reported`, `not_applicable`, `unclear`, and `negative` (an observed negative only); none is silently converted into another. Cohort overlap is documented at cohort level without participant deduplication.

Before the first retention calculation, publication review identified two documented EHR-plus-self-report composites (GCST90668057 and GCST90837197). These are classified `mixed_ehr_self_report_susceptibility`, not mislabeled direct clinical/EHR mixtures. They are excluded from the strict EHR comparison and included in an explicitly labeled EHR-plus-composite bound. GCST007996 remains unresolved because its exact COPD supplement rule was not retrieved; a separate secondary-inclusion bound is reported. Also report high-confidence-primary-only, primary plus conditional smoking coefficients, and leave-one-publication support loss. None changes the primary hierarchy. Reviewer-input wording corrections and the two composite eligibility judgments are recorded in `data/phenotype_adjudication_amendments.json` and transparent normalization rules in the script, before retention counts were examined.

## Frozen inputs and statistical support

Use S2-R002 GWS association rows at the unchanged P <= 5e-8 cutoff to define study-specific support for each of the 660 tag identifiers. The aggregated tag table's `study_accessions` also includes some nonsignificant reports; preserve those as original provenance but do not promote them to significant support. Explicitly audit the difference.

For each LD relationship, expand `supporting_studies` in S2-R006B into study × tag × ancestry-panel rows, confirm the study/tag has a GWS association, retain the rows eligible for the stratum, and join the surviving tag-panel pairs to S2-R006C. The LD pair table has no study field. Preserve exact tag support regardless of whether the tag could be reference matched or LD expanded. Retain all 85 unmatched tags with their frozen reasons. Do not create missing LD links or silently borrow a panel supported only by an excluded study.

Compare all 660 tags; the frozen 15,389 records in S2-R006E; 337 R010 candidates in original rank order; 153 frozen S6 components; and 12 frozen S6 shortlist entries. A component is retained if any existing member is retained; also report all/partial/no-member retention. Do not rebuild components. Include the original rank for every R010 record and `not_applicable` for records outside R010.

## Support states and outputs

Report nonexclusive primary/EHR/ML/other/unresolved indicators and exact supporting accessions, classes, tags, panels, and cohorts. The exclusive display state follows: primary-supported; otherwise multiple classes; otherwise EHR-only; ML-only; other-only; unresolved. Explicit mixed-class indicators accompany the primary-supported state. Class-specific counts need not sum to the broad count.

GCST90244098 dependence is measured twice: the frozen, aggregated V1 linked-study provenance (to reproduce planning counts) and the exact study/tag/panel reconstruction used for the main sensitivity. Explain any differences. `ML-only` at phenotype-class level and `GCST90244098 sole supporting accession` are separately labeled. Loss of primary support is a limitation of genetic phenotype provenance, not biological falsity.

## Interpretation and validation

Use descriptive counts/proportions rather than treating correlated variants or overlapping accessions as independent samples. No enrichment significance test is planned. Summarize frozen consequence, model-context, annotation, and selected-locus labels solely to assess the scope of V1 claims. Do not derive new target genes, model scores, ranks, or causal probabilities.

Validate complete accession coverage, legal/mutually compatible classes and eligibility, panel join integrity, broad-universe recovery, stable R010 ordering, frozen component memberships, shortlist identity, explicit unresolved tags, independent set-based reconstruction, and all input hashes. Write scripts, evidence registers, result tables, the report, manifests, and checks exclusively under `07_gap_closure/`.
