# COPD regulatory genomics project

## Target disease

Chronic obstructive pulmonary disease (COPD) is the primary phenotype. GWAS
Catalog records mapped directly to COPD ontology terms are included in the core
analysis. Quantitative lung-function traits, emphysema-only studies, chronic
bronchitis-only studies, asthma-COPD overlap, exacerbation traits, and smoking
behavior are retained as contextual or related evidence and are not merged into
the core COPD estimates unless the source explicitly defines a COPD cohort.

## Framework

This workspace follows `frameworks/DiseaseFramework_SL091626.md` in order.
Each numbered directory corresponds to one framework section. Scientific claims
must have a stable result identifier and a verified entry in `sources.tsv`.

## Analysis conventions

- Human genome build: GRCh38, with source builds and coordinate conversion recorded.
- Population scope: all ancestries, ages, sexes, and geographies, preserving source definitions and missing strata.
- Evidence cutoff: 2026-10-01 for the initial run; later updates must record their own access date.
- Disease-specific inputs and all generated outputs stay inside `diseases/COPD/`.
- Shared `data/`, `models/`, and `reference_papers/` assets are read-only.

## Section 1 plan

1. Define COPD phenotype boundaries and diagnostic terminology.
2. Compile epidemiology by age, sex, geography, and population where available.
3. Summarize family, twin, and SNP-heritability evidence by ancestry.
4. Identify replicated genes, pathways, affected lung tissues and cell types, and relevant ENCODE/Roadmap biosamples.
5. Review comorbidities, childhood origins, progression, diagnosis, and treatments.
6. Record unresolved questions and evidence gaps with verified citations.

## Current status

The 2026-10-01 public-data workflow is complete through Sections 1 to 6 and the
manuscript-level internal audit.
It includes a 104-study core COPD GWAS inventory, ancestry-aware LD expansion,
single-donor severe-emphysema lung regulatory maps, held-out enhancer and
silencer models, a 337-candidate ranked list, public computational validation,
and a checksum-verified 12-candidate experimental handoff. Wet-lab work was
designed but not performed. The manuscript, claim map, figure source data,
supplementary checksum manifest, and strict final audit are under `manuscript/`.
The draft still requires named-author review, declarations, institutional ethics
wording, journal-specific disclosure, and a persistent public archive before
journal submission.
