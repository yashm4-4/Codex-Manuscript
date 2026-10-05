# Section 1: COPD background

Evidence cutoff: 2026-10-01. Stable outputs are registered as COPD-S1-R001
through COPD-S1-R008. Detailed values and provenance are in the accompanying
TSV files and `../sources.tsv`.

## Phenotype and diagnostic boundary

Chronic obstructive pulmonary disease is a heterogeneous lung condition with
chronic respiratory symptoms caused by airway and/or alveolar abnormalities
that produce persistent, often progressive airflow obstruction. The operational
GOLD 2026 diagnosis requires an appropriate clinical context and a
post-bronchodilator FEV1/FVC ratio below 0.70 [COPD-SRC-001]. Pre-COPD and
preserved-ratio impaired spirometry are risk states, not established COPD. This
study therefore treats GWAS Catalog records mapped only to COPD as the core
phenotype. Emphysema, chronic bronchitis, quantitative lung-function traits,
asthma-COPD overlap, exacerbation, biomarker, treatment-response, and smoking
interaction records are retained separately.

Routine molecular evaluation does not replace spirometry. Alpha-1 antitrypsin
deficiency is the principal established etiologic molecular subtype for which
testing is recommended broadly in COPD. Testing begins with serum alpha-1
antitrypsin and at least the SERPINA1 S and Z alleles, with phenotype or expanded
sequencing when clinical suspicion persists [COPD-SRC-023]. Blood eosinophils
guide anti-inflammatory treatment response rather than establish the diagnosis
[COPD-SRC-024]. No polygenic score or other omics assay is validated to replace
clinical assessment and post-bronchodilator spirometry.

## Epidemiology and population boundaries

GBD 2021 estimated 213.39 million prevalent cases, 16.9 million incident cases,
3.72 million deaths, and 79.78 million disability-adjusted life years worldwide
[COPD-SRC-002]. WHO reported a later estimate of 3.4 million deaths in 2023 and
that nearly 90% of deaths below age 70 occurred in low- and middle-income
countries [COPD-SRC-003]. These estimates use different reference years and
modeling systems and are not pooled.

Case definition materially changes prevalence. Among adults age 30 to 79 in
2019, spirometry-based modeling estimated 391.9 million cases using a fixed
FEV1/FVC ratio, compared with 292.0 million using the lower limit of normal
[COPD-SRC-004]. In the 2023 US National Health Interview Survey, age-adjusted
self-reported diagnosed prevalence was 3.8%, with estimates of 4.1% in women and
3.4% in men. Prevalence increased from 0.4% at age 18 to 24 to 10.5% at age 75
or older [COPD-SRC-005]. Survey race and Hispanic-origin categories are social
and administrative categories, not genetic ancestry. Detailed denominators,
uncertainty, geography, and limitations are in `epidemiology.tsv`.

## Heritability and genetic architecture

COPD is heritable, but estimates depend on phenotype and scale. A 2025 Danish
twin study estimated clinical COPD heritability at 47% with a wide 95% confidence
interval of 16% to 78%; the fixed-ratio sensitivity estimate was approximately
55% [COPD-SRC-006]. Nordic registry studies of hospitalized COPD estimated
additive genetic contributions near 61% to 63%, but captured severe disease
[COPD-SRC-007]. Family studies of severe probands also found excess airflow
obstruction among relatives [COPD-SRC-008, COPD-SRC-009].

Common-variant estimates are lower and cannot be compared directly with twin
estimates. COPDGene GREML liability estimates were 37.7% (SE 7.4%) in
non-Hispanic White smokers and 37.9% (SE 20.4%) in African-American smokers
[COPD-SRC-010]. The second estimate was imprecise, and a null between-group test
does not establish equivalent architecture. UK Biobank LDSC yielded 20.11%
(SE 2.61%) on the observed scale [COPD-SRC-011]. Evidence remains sparse for
continental African, Indigenous, South Asian, Latin American, non-smoking, and
biomass-exposed populations.

Replicated regions include HHIP, FAM13A, and CHRNA3-CHRNA5-IREB2, with additional
support for RIN3, MMP12, and TGFB2 [COPD-SRC-012]. A large COPD GWAS identified
82 loci in 35,735 cases and 222,076 controls; the lead variants explained at
most about 7% of liability [COPD-SRC-013]. Lung development, WNT, MAPK/ERK,
nerve-growth-factor, extracellular-matrix, protease-antiprotease, epithelial,
and inflammatory programs recur. Lung-function GWAS provide useful context but
are not relabeled as COPD case-control evidence [COPD-SRC-014].

## Tissues, cells, and regulatory resources

The principal compartments are terminal and small conducting airways, alveolar
parenchyma, and pulmonary vasculature. Relevant epithelial cells include basal,
goblet, club/secretory, ciliated, AT1, and AT2 cells. Relevant stromal and
vascular states include alveolar and CTHRC1-positive fibroblasts, airway and
vascular smooth muscle, and capillary, aerocyte, and arterial endothelium.
Alveolar/interstitial macrophages, T cells, B cells, NK cells, dendritic cells,
and neutrophils contribute heterogeneous inflammatory programs
[COPD-SRC-015, COPD-SRC-016].

Integration of 7,285 fine-mapped variants across 82 loci found 250 overlaps with
open chromatin and enrichment in primary lung fibroblasts, AT2 cells,
small-airway epithelial cells, and 16HBE cells [COPD-SRC-017]. ENCODE metadata
inspection identified a particularly relevant model source: donor ENCDO520EJG,
whose record states severe emphysema, has matched ATAC-seq, H3K27ac, and
H3K27me3 experiments in three lung lobes [COPD-SRC-018]. This donor-matched
disease tissue is the primary Section 3 and 4 context. Its single-donor,
bulk-tissue, age, sex, ancestry, and comorbidity limitations require normal-lung
and fibroblast sensitivity analyses. Exact accessions are in
`cell_tissue_resources.tsv`.

## Comorbidity, development, and progression

COPD progression is heterogeneous. In a 22-year study, COPD developed in 26%
of participants who began young adulthood with FEV1 below 80% predicted and 7%
of those with normal FEV1. Approximately half of incident cases followed a
rapid-decline trajectory, while the remainder began with low maximally attained
lung function [COPD-SRC-020]. COPD is rarely an established pediatric diagnosis,
but prematurity, childhood infection, persistent asthma, smoke and pollution,
and impaired lung growth can establish adult susceptibility. Improved air
quality was associated with fewer children having low FEV1 at age 15
[COPD-SRC-022].

Common comorbidities include cardiovascular disease, lung cancer,
bronchiectasis, anxiety and depression, osteoporosis, sarcopenia, diabetes,
sleep apnea, gastroesophageal reflux, anemia, and chronic kidney disease.
Published prevalence ranges are broad and source-dependent, so they are not
treated as universal estimates [COPD-SRC-021].

## Treatment and pharmacogenomics

Exposure reduction, vaccination, rehabilitation, exercise, selected oxygen or
noninvasive ventilation, and appropriate procedures remain foundational.
Bronchodilators target ADRB2 or muscarinic receptors. Inhaled corticosteroids
bind the nuclear receptor NR3C1 and alter transcription, but do not edit DNA.
Roflumilast inhibits PDE4. Ensifentrine, approved in 2024, inhibits PDE3 and
PDE4 [COPD-SRC-025]. Dupilumab, approved for eosinophilic COPD in 2024, blocks
IL-4 receptor alpha [COPD-SRC-026]. Mepolizumab, approved for inadequately
controlled eosinophilic COPD in May 2025, neutralizes IL-5
[COPD-SRC-027, COPD-SRC-028]. Alpha-1 proteinase inhibitor augmentation is
restricted to severe SERPINA1-deficiency emphysema [COPD-SRC-023].

Tozorakimab, itepekimab, tezepelumab, and astegolimab remain investigational for
COPD at the evidence cutoff, with positive, discordant, or negative programs as
specified in `treatments.tsv` [COPD-SRC-029 through COPD-SRC-032]. No approved
therapy directly targets a COPD enhancer or silencer. ADRB2 candidate-variant
results are inconsistent, and bronchodilator-response GWAS findings have not
produced routine genotype-guided prescribing [COPD-SRC-033, COPD-SRC-034].

## Open questions

The primary unresolved questions are which pre-COPD states progress, which
noncoding GWAS variants are causal, which cell states mediate each locus, and
how regulation changes during injury, repair, and exacerbation. Donor-matched
chromatin data are needed across small-airway epithelium, alveolar cells,
endothelium, macrophages, fibroblasts, sex, ancestry, exposure, and disease
stage. Polygenic and molecular endotypes require validation across populations.
Therapeutic priorities include predicting decline and biologic response,
identifying resilient exposed individuals, restoring terminal airways and
alveoli, and improving diagnosis and exposure reduction in high-burden settings.

## Evidence gaps

No reliable global incidence or prevalence table simultaneously stratified by
age, sex, geography, and genetic ancestry was identified. GWAS ancestry is
incompletely represented, and epidemiologic race categories cannot substitute
for genetic ancestry. Comprehensive donor-matched disease-state ATAC/DNase,
H3K27ac, and H3K27me3 maps are unavailable for every prioritized lung cell
type. Controlled biobank data, participant-level clinical phenotypes, and
wet-lab validation are outside the public-data scope of this section.
