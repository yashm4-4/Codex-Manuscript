# COPD enhancer and silencer modeling

Evidence cutoff: 2026-10-01. Genome build: GRCh38. This report integrates
`COPD-S4-R001` through `COPD-S4-R010`. It describes predictions from one
severe-emphysema bulk-lung donor. It does not establish variant causality,
transcription-factor occupancy, or clinical utility.

## Disease-relevant models and held-out performance

Two TREDNet v2 phase-II models were fine-tuned from the same frozen phase-I
representation (COPD-SRC-060). Positives were 1-kb windows centered on
three-lobe donor-union ATAC peaks that overlapped H3K27ac for the enhancer
model or H3K27me3 for the silencer model. Controls came from the reusable
non-promoter, non-exon DHS pool after removal of the ENCODE blacklist,
model-matched histone peaks, and positive windows. Canonical chromosomes only
were retained. Sampling used seed 20261001 and an exact 1:2
positive-to-control ratio overall.

| Model | Training positives | Training controls | Validation positives | Validation controls | Test positives | Test controls |
|---|---:|---:|---:|---:|---:|---:|
| Enhancer | 156,789 | 307,473 | 8,460 | 20,328 | 16,103 | 34,903 |
| Silencer | 26,326 | 51,839 | 1,612 | 3,427 | 2,646 | 5,902 |

Training used chr1-6, chr10-22, chrX, and chrY; chr7 was validation; chr8-9
were held out for testing. Both models ran 50 epochs with deterministic
TensorFlow operations, a nominal batch size of 256, binary cross-entropy,
Adadelta, and best-validation-loss checkpoint restoration. The selected
enhancer checkpoint was epoch 42 and the selected silencer checkpoint was
epoch 37.

| Model | Held-out n, positive/control | ROC AUC (95% bootstrap CI) | PR AUC (95% bootstrap CI) | Brier score | Log loss | ECE |
|---|---|---:|---:|---:|---:|---:|
| Enhancer | 51,006, 16,103/34,903 | 0.938272 (0.936441-0.940014) | 0.885103 (0.881510-0.888851) | 0.090628 | 0.292544 | 0.018553 |
| Silencer | 8,548, 2,646/5,902 | 0.967243 (0.963408-0.970706) | 0.937608 (0.931289-0.943635) | 0.061341 | 0.208854 | 0.011078 |

Confidence intervals used 500 stratified bootstrap replicates, separately
resampling positives and controls. The score thresholds at nominal held-out
false-positive rates of 10%, 5%, 3%, and 1% were, respectively, 0.464482,
0.643623, 0.738453, and 0.861629 for the enhancer model, and 0.294698,
0.585050, 0.720549, and 0.864801 for the silencer model. Observed
false-positive rates at the selected 5% thresholds were 0.049938 and 0.050152.

The biosample is ENCODE donor `ENCDO520EJG`, a 60-year-old man of reported
European ancestry with severe emphysema (COPD-SRC-018). ENCODE does not
provide post-bronchodilator spirometry or an adjudicated COPD diagnosis. Both
models therefore capture one bulk-lung epigenomic context, not COPD
case-control differences.

## Allele-specific scoring identifies 337 candidate regulatory variants

Variant-centered 2,001-bp REF and ALT sequences were constructed for the
15,389 frozen GWAS-tag and ancestry-matched LD-proxy records. Both alleles were
validated and scored for 15,303 records (99.44%). Eighty-two unresolved
alleles, three missing coordinates, and one ambiguous ALT were not scored.

A predicted causal call required all of the following: a successfully scored
sequence, no ENCODE blacklist overlap, a maximum REF/ALT region score at or
above the model's held-out 5% FPR threshold, and an absolute REF-ALT score
difference at or above the 95th percentile among eligible candidates in the
same broad allele class. SNVs and indel/complex variants were calibrated
separately. The absolute-delta cutoffs were 0.0570632 and 0.0493609 for
enhancer SNVs and indels, and 0.0288026 and 0.0261834 for silencer SNVs and
indels.

The 5% FPR region threshold identified 1,392 enhancer and 1,140 silencer
records, 1,830 in union. The joint region-and-delta definition identified 175
enhancer and 199 silencer candidates, with 37 shared between models and 337 in
union. The union comprised 297 SNVs and 40 indel/complex records. Twenty-one
were GWS tag records and 324 were LD proxies; these roles overlap. No
predicted causal record overlapped the ENCODE blacklist.

Among the 337 candidates, 138 were enhancer-only, 162 were silencer-only, and
37 were called by both models. Enhancer calls included 91 predicted losses and
84 gains; silencer calls included 107 losses and 92 gains. Fifteen candidates
overlapped a GENCODE CDS, 50 overlapped a refined donor enhancer or silencer,
216 overlapped SCREEN, Ensembl, or FANTOM5 regulatory annotations, and 220
overlapped either donor-refined or general known-regulatory annotations. Those
overlaps are contextual and are not independent functional validation.

## Locus classification suggests context-specific missing regulation

Eighty-five of 140 replicated COPD GWAS-gene loci (60.71%) contained at least
one predicted causal regulatory candidate. The mutually exclusive locus
classes were 45 coding plus predicted-causal regulatory, 15 coding only, 40
predicted-causal regulatory only, and 40 other. Six loci (4.29%) had neither
an observed regulatory annotation nor a 5% FPR prediction. The 40-locus other
class includes loci with noncausal model calls or observed annotations and is
not equivalent to the six no-regulatory loci. Additional donors, airway and
alveolar cell types, disease stages, exposure states, and single-cell models
are warranted for this class.

## DeepFootprinting nominates distinct enhancer and silencer motif programs

DeepExplainer attribution and TF-MoDISco were run on 1,000 positive sequences
per model (COPD-SRC-061; COPD-SRC-063; COPD-SRC-071). Enhancer analysis yielded 463 input
positive seqlets, five retained patterns, and 249 pattern-assigned seqlets.
Six of 15 reported JASPAR 2024 matches passed q <= 0.05. Silencer analysis
yielded 343 input positive seqlets, four patterns, and 256 pattern-assigned
seqlets; all 12 reported JASPAR matches passed q <= 0.05 (COPD-SRC-062).

Significant enhancer matches nominated EHF, ELF3, IKZF3, BATF::JUN, BATF3,
and BATF. Significant silencer matches nominated CTCF-family, PATZ1, KLF12,
KLF16, SPIB, SPI1, ERG, FOS::JUN, FOSL2, and FOSL2::JUN motifs. These are
sequence-similarity hypotheses, not direct binding measurements.

A motif-compatible site required a best two-strand PWM score at least 0.8 of
the theoretical motif range, and creation or disruption required an allele to
cross that threshold. Among 175 enhancer candidates, 49 (28.00%) had a
compatible site, 21 (12.00%) disrupted one, 17 (9.71%) created one, and 38
(21.71%) disrupted or created at least one. The leading altered-site
hypotheses were IKZF3 (20 candidates), EHF (18), and ELF3 (16). Among 199
silencer candidates, 144 (72.36%) had a compatible site, 62 (31.16%)
disrupted one, 55 (27.64%) created one, and 103 (51.76%) disrupted or created
at least one. The leading hypotheses were KLF16 (37), CTCF (31), and KLF12
(25).

Motif compatibility was common outside model-positive regions: 1,726
candidates outside enhancer predictions contained a significant
enhancer-model motif and 7,346 outside silencer predictions contained a
significant silencer-model motif. Across both models, 7,151 of 13,559
candidates outside either 5% FPR region call contained at least one such site.
This shows why motif occurrence alone cannot substitute for a regulatory
element prediction or occupancy assay.

## Population frequencies and ancestral states

Exact alternate-allele frequencies were recalculated from 1000 Genomes Phase
3 GRCh38 genotypes for all 337 candidates across AFR, AMR, EAS, EUR, and SAS
superpopulations (COPD-SRC-037). A variant was classified common if its
maximum superpopulation ALT frequency was at least 0.05, low frequency if at
least 0.01 but below 0.05, and rare if below 0.01. There were 323 common, 11
low-frequency, and three rare candidates. None met the strict
population-specific rule of ALT AF at least 0.01 in exactly one
superpopulation and below 0.001 in each other superpopulation.

Ensembl ancestral alleles could be mapped to REF or ALT for 310 candidates.
Seventy-one of 310 (22.90%) had global derived-allele frequency above 0.5, and
131 of 310 (42.26%) exceeded 0.5 in at least one superpopulation. The
remaining 27 candidates have unresolved ancestral state and were excluded
from those denominators.

## Target-gene hypotheses remain provisional

GENCODE v50 proximity, GWAS Catalog mapped-gene links propagated through
r2 >= 0.8 LD, selected-locus membership, and direct CDS overlap produced
4,779 candidate-gene-method rows for 337 candidates and 1,635 distinct gene
labels (COPD-SRC-035; COPD-SRC-043). Each candidate has nearest-any-gene and
nearest-protein-coding TSS hypotheses, plus all TSSs within 100 kb where
present. No matched severe-emphysema donor rE2G, HiChIP, or equivalent
chromatin-contact map was available. The large target set, including many
unnamed GENCODE entries and proximity-only links, must not be interpreted as
a set of proven targets. Section 5 evaluates these hypotheses against GTEx
v10 Lung eQTLs.

## Comparison with experimentally studied COPD variants

Section 3 recorded eight functional evidence rows spanning eight distinct
rsIDs. Exact identity matching used only GWAS tag IDs and Ensembl variation
aliases; LD source-focal labels were not treated as identity. Two of eight
tested rsIDs were exact frozen candidates. The endogenous-edited FAM13A
variant rs2013701 was recovered as predicted-causal enhancer candidate
`4:88963935:G:T`, rank 210. The reporter-tested GWAS tag rs7671167 was
scorable but its enhancer and silencer region scores were both below their
held-out 5% FPR thresholds, so it was not called causal. The other six tested
rsIDs were not exact identities in the frozen resolved set, making model
recovery non-estimable rather than negative. Thus, one of 337 predicted
causal candidates (0.30%) had exact prior literature evidence at an
experimentally tested rsID. Section 5 separately identified one regional
MPRAbase overlap, which did not test both candidate alleles.

Reasons for non-recovery include absence from the frozen tag-plus-LD candidate
identity set, cell-context differences between bronchial epithelial assays
and bulk emphysematous lung, region scores below calibrated thresholds, and
the requirement for a top-5% allele delta. Recovery requires deliberate
inclusion of literature variants, additional cell-type and donor models,
region and allelic assays, and contact-informed target mapping. Thresholds
should not be relaxed after inspecting known variants without a prespecified
sensitivity analysis.

## Ranked experimental list

`COPD-S4-R010_THE_LIST.tsv` is the complete 337-row ranked list. The order
was frozen in R004: regulatory evidence tier, number of causal model calls,
matched-class absolute-delta percentile, region-threshold ratio, maximum
tag-proxy r2, and stable candidate identifier. Population, motif, and target
annotations were added without reranking. The table is a prioritization
instrument, not a posterior causal probability.

## Reproducibility

Production scripts, logs, model checkpoints, prediction tables, input and
output SHA-256 values, and software versions are retained under
`04_modeling/`. R001 contains held-out predictions and thresholds; R002
audits sequence construction; R003 contains paired allele scores; R004
contains causal and locus classifications; R005 contains population genetics;
R006 contains target hypotheses; R007 and R008 contain attribution, motif,
and allele-specific TFBS results; R009 audits literature recovery; and R010
is THE LIST. The streaming phase-I embedding implementation reproduced the
legacy implementation exactly for the prespecified test batch (maximum
absolute difference 0.0).

## Limitations

The strongest limitation is biological replication: all cell-matched training
peaks came from one bulk-lung severe-emphysema donor. Anatomy, smoking
exposure, cell mixture, treatment, comorbidity, and inter-individual
differences cannot be separated. H3K27me3 is a broad repressive mark and does
not itself prove silencer activity. Chromosome holdouts reduce local sequence
leakage but do not provide an external donor, disease-control, or assay
benchmark. The 337 candidates were generated by GWAS tag and r2-based LD
expansion, not statistical fine-mapping. Thresholds define high-scoring model
perturbations, not causal effects. Motif names can represent TF families and
do not demonstrate occupancy. Proximity and GWAS-mapped genes are hypotheses,
and population frequencies do not establish selection or disease effect.
Direct allele, chromatin, expression, and phenotype perturbations are needed.
