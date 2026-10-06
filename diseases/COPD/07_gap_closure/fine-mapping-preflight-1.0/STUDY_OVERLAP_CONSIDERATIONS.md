# Cohort overlap and estimand separation

The 104-row ancestry/cohort table retains frozen ascertainment and current Catalog discovery/replication sample descriptions. The original phenotype-stage cohort and pairwise-overlap tables are hash-bound inputs, not recalculated independent-replication counts. Shared cohort names indicate possible participant overlap; absence of a shared name does not establish independence. No participant-level overlap estimate is available from this audit.

## UK Biobank analyses

Kim's ordinary ever/never-smoker and current/noncurrent-smoker COPD analyses are four strata, not four independent discovery cohorts. Ever and never partition one exposure definition; current and noncurrent partition another, with substantial overlap across the two partitions. Do not combine all four by ordinary independent-effects meta-analysis. A primary ever-smoker analysis and a separately labelled never-smoker analysis preserve distinct conditional population estimands; current/noncurrent results can only be additional overlapping sensitivity analyses, not independent replication.

The two Kim SNP main-effect coefficients from models containing smoking interaction terms remain separate from the four ordinary stratum GWAS. They are conditional on the encoded reference exposure. The paper's 48 and 55 joint-test loci do not enumerate loci for the four marginal smoking-stratum GWAS.

Pan-UKB J44, Pan-UKB PheCode 496/496.21, earlier UKB SAIGE/PheWeb, UKB exome and whole-genome sequencing, later UKB phenotyping algorithms, ML spirogram liability, and the UKB contribution to Sakornsakolpat/GBMI reuse participants to varying degrees. Comparisons across their phenotypes are not independent cohort replication. Similar UKB origin alone does not prove identical LD samples or covariate projection. Sequencing-derived rare variants may be absent from array/imputation LD resources.

## BioBank Japan and other meta-analyses

Ishigaki original physician-diagnosed COPD and later Sakaue EHR/text-mined definitions are not interchangeable. Their BBJ cohorts overlap. Male, female and combined-sex releases also overlap by design; the combined release must not be meta-analyzed with its sex-specific constituents. Japanese ancestry does not make every East Asian reference equally appropriate. Original study controls included ToMMo/IMM/JPHC/J-MICC population cohorts as well as BBJ participants without relevant diagnoses; an LD matching contract must resolve that pooled control/sample composition rather than treating all controls as BBJ.

GBMI includes multiple biobanks and ancestry-specific as well as pooled releases. UKB/BBJ/FinnGen and other constituent overlap must be mapped to each exact GBMI release before any cross-study synthesis. A pooled all-ancestry result cannot use a single-European LD matrix. Even an ancestry-specific meta-analysis requires consideration of cohort ancestry structure, differing sample size across variants, and differing covariate/phenotype models.

MVP PheCode COPD/airway-obstruction files provide useful secondary ancestry-specific resources but do not inherit UKB cohort-matched LD. Older COPD consortia reuse COPDGene, ECLIPSE, GenKOLS, NETT/NAS and other cohorts. Replication-stage sample counts may relate to selected follow-up variants rather than the full discovery GWAS and must not be added to every variant's N.

## Sample-size and pooling rules

- Preserve reported cases and controls separately from total N, effective N and per-variant N.
- Do not sum duplicate sample entries in Catalog YAML metadata. The audited Sakaue metadata include repeated ancestry/sample-size entries; these are not demonstrated additional participants.
- Before an RSS analysis of a binary outcome, document the method's sample-size approximation and any unbalanced-case/rare-variant limitations. A case/control effective-size formula is a modelling choice, not permission to replace unexplained per-variant sample sizes.
- Do not assume a heterogeneous meta-analysis has the same genotype covariance as its largest constituent cohort. Obtain compatible ancestry-specific statistics/LD or a justified multi-population method and its actual required inputs.
- Do not pool across direct COPD, EHR/ICD, ML liability, within-case decline/severity, bronchitis, smoking interactions or non-COPD lung-function traits to increase apparent power.

No meta-analysis, overlap adjustment, locus selection, fine-mapping, colocalization or target-gene analysis was performed here.
