# COPD fine-mapping and risk-direction preflight 1.0

This stage audits public data access, phenotype fidelity, signed summary-statistic schemas, ancestry/LD compatibility and prospective statistical prerequisites. It does not execute statistical fine-mapping, model inference, candidate scoring, candidate reranking, colocalization, gene assignment or manuscript changes. No commit or push is authorized.

The starting population is all 104 accessions in the frozen COPD core study inventory. The frozen phenotype adjudication supplies the biological hierarchy; candidate support counts are not selection inputs. Source-specific alternatives may be added only because they measure relevant COPD phenotypes and offer potentially usable statistics/LD, not because of candidate or benchmark overlap.

Direct clinical/spirometric case-control COPD is preferred for a primary analysis if statistical and LD prerequisites can be demonstrated. EHR/ICD/PheCode, mixed EHR/self-report, ML liability, within-case severity/progression, chronic bronchitis/subphenotypes, smoking-interaction coefficients and non-COPD lung-function traits remain distinct estimands. No pooling across these classes is proposed.

Access is tested using public metadata and bounded file prefixes/headers. A successful HTTP response or Catalog availability flag is not evidence of whole-file completeness, valid signed effects, matched LD, or locus eligibility. Partial samples are explicitly labelled; no full-file checksum is invented from a prefix. Failed or restricted endpoints remain documented, not silently replaced.

Readiness has separate levels: source advertised; metadata inspected; header/representative records verified; whole-file integrity/coverage verified; matched-LD inputs available; locus-level QC passed. This preflight may establish the earlier levels but must not claim later unperformed checks passed. Unknown prerequisites block execution rather than being filled by assumption.

The source selection firewall prohibits accessing candidate membership or external functional benchmark outcomes to define datasets, loci or eligibility. Reading the historical phenotype report establishes its hierarchy only; its inherited candidate examples/counts are not consumed by the preflight scripts. Loci must later arise from the chosen GWAS's own association evidence under prospective rules.

Earlier stages, all model checkpoints and original registers are immutable except new append-only rows in the three V2 registers. Final validation will compare historical Git tree entries and original register prefixes with baseline commit `1f41df67dcf18d626a3862a38877919d2d3114a4`.
