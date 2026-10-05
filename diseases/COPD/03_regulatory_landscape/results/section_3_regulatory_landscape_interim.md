# COPD regulatory landscape: interim results

Superseded for the complete computational analysis by
`section_3_regulatory_landscape.md`; retained as the pre-R003 frozen snapshot.

Evidence cutoff: 2026-10-01. Genome build: GRCh38. This report covers results
`COPD-S3-R001` and `COPD-S3-R002`. Candidate-variant intersections are excluded
until the Section 2 tag plus LD-proxy set is frozen.

## Regulatory reference preparation (COPD-S3-R001)

Five local shared regulatory resources were normalized to a common six-column,
0-based, half-open BED representation restricted to `chr1-22`, `chrX`, and
`chrY`. The source files were not modified or downloaded again.

| Reference | Release | Canonical intervals retained | Noncanonical excluded | Invalid excluded |
|---|---:|---:|---:|---:|
| ENCODE SCREEN cCREs (COPD-SRC-038) | Registry V3 | 1,063,878 | 0 | 0 |
| Ensembl Regulatory Build (COPD-SRC-039) | 116 | 380,818 | 0 | 0 |
| FANTOM5 CAGE enhancers (COPD-SRC-040) | hg38_latest | 63,285 | 0 | 0 |
| UCSC RepeatMasker (COPD-SRC-041) | hg38 | 5,317,286 | 366,404 | 0 |
| ENCODE blacklist (COPD-SRC-042) | hg38 v2 | 636 | 0 | 0 |

The complete feature-class breakdown and checksums are in
`COPD-S3-R001_regulatory_reference_feature_counts.tsv` and
`COPD-S3-R001_regulatory_reference_inventory.tsv`. These references are ready
for candidate-variant annotation but no overlap fraction is reported before the
candidate set is final.

## Disease-relevant epigenomic context

The primary context uses ENCODE donor `ENCDO520EJG`, whose metadata describe a
60-year-old male donor of reported European ancestry with "lungs with severe
emphysema" (COPD-SRC-018). Released GRCh38 pseudoreplicated ATAC-seq, H3K27ac,
and H3K27me3 narrowPeak files were available for upper-right, lower-right, and
lower-left lung lobes. Assays are matched at the donor and anatomical-lobe
levels, not at the biosample-accession level. ENCODE provides neither
post-bronchodilator spirometry nor clinical COPD adjudication for this donor.

Preliminary enhancers were H3K27ac peaks and preliminary silencers were
H3K27me3 peaks. Refined elements were ATAC peaks overlapping the corresponding
same-lobe histone mark by at least 1 bp. Across lobes, preliminary peak counts
ranged from 60,451 to 67,566 for enhancers and 33,416 to 68,141 for silencers.
Refined counts ranged from 91,110 to 113,178 accessibility peaks for enhancers
and 13,981 to 22,432 accessibility peaks for silencers. Merging overlapping or
book-ended regions across lobes produced 77,164 preliminary enhancer regions,
78,899 preliminary silencer regions, 71,746 refined enhancer regions, and
15,699 refined silencer regions.

## GWAS-gene locus construction

Section 2 selected 155 gene symbols: 143 observed in at least two core COPD
study accessions and 12 additional effect-size-prioritized genes. Coordinate
analysis required one exact GENCODE v50 gene record. `Y_RNA` mapped to 756
records and `GUSBP5` and `CYP2B7P` each mapped to two, so all three were
excluded rather than assigning arbitrary coordinates. The analyzed sets
therefore contain 152 selected genes and 140 replicated genes. The replicated
gene-body plus 100 kb windows merged to 102 intervals spanning 42,874,207 bp;
the selected windows merged to 109 intervals spanning 47,285,591 bp.

## Regulatory burden in replicated COPD GWAS-gene loci

Counts below are distinct peak calls for individual lobes and distinct merged
regions for the donor union. The percentages in parentheses are the fractions
of the 140 unambiguous replicated gene loci with at least one refined element.

| Context | Preliminary enhancers in locus union | Preliminary silencers in locus union | Refined enhancers in locus union | Replicated loci with refined enhancer | Refined silencers in locus union | Replicated loci with refined silencer |
|---|---:|---:|---:|---:|---:|---:|
| Upper-right lobe | 1,549 | 657 | 2,230 | 128 (91.43%) | 286 | 69 (49.29%) |
| Lower-right lobe | 1,497 | 762 | 2,455 | 126 (90.00%) | 352 | 74 (52.86%) |
| Lower-left lobe | 1,652 | 1,087 | 2,684 | 127 (90.71%) | 435 | 93 (66.43%) |
| Three-lobe donor union | 1,892 | 1,201 | 1,815 | 132 (94.29%) | 282 | 94 (67.14%) |

In the donor union, refined enhancer and silencer densities were 42.333 and
6.577 regions per megabase of merged replicated-locus sequence, respectively.
The median per replicated gene locus was 11 refined enhancers and 2 refined
silencers. Preliminary and refined counts are not nested because they use
different units: histone peaks for the preliminary definition and
accessibility peaks for the refined definition.

## Regulatory burden in the full selected gene set

Among the 152 unambiguous selected loci, the donor union contained 2,027
preliminary enhancer regions, 1,241 preliminary silencer regions, 1,953 refined
enhancer regions, and 298 refined silencer regions. At least one refined
enhancer occurred in 143 loci (94.08%), and at least one refined silencer
occurred in 98 loci (64.47%). The corresponding densities were 41.302 refined
enhancers and 6.302 refined silencers per megabase of merged selected-locus
sequence. Full lobe-specific results and per-gene values are in
`COPD-S3-R002_lobe_union_summary.tsv` and
`COPD-S3-R002_per_gene_burden.tsv`.

## Interpretation and limitations

Most replicated COPD GWAS-mapped gene windows contain accessible,
H3K27ac-marked regions in this severe-emphysema lung donor, while accessible,
H3K27me3-marked regions are less ubiquitous. This is a descriptive overlap and
does not establish variant function, regulatory target genes, or disease
causality. The one-donor bulk-lung design cannot separate disease effects from
individual, anatomic, exposure, treatment, or cell-composition effects. The
donor metadata also report hypertension. H3K27ac includes promoter-associated
activity, and H3K27me3 is not direct evidence of a silencer. No promoter or
blacklist exclusion was applied. Finally, a 100 kb flank is an operational
window; chromatin-contact and expression-QTL evidence are required to assign
targets.

## Reproducibility checks

All nine ENCODE peak files matched the SHA-256 values in the download manifest.
All derived intervals use canonical GRCh38 chromosomes and valid half-open
coordinates. An independent bedtools 2.31.1 intersection reproduced the six
same-lobe refined counts exactly. Reruns from `/tmp` reproduced every result TSV
and compressed BED file byte-for-byte. Replay commands, versions, definitions,
input checksums, and output checksums are recorded in the R001 and R002
manifests.
