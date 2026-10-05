# COPD regulatory landscape: computational results

Evidence cutoff: 2026-10-01. Genome build: GRCh38. This report covers
`COPD-S3-R001` through `COPD-S3-R003`: shared regulatory-reference
normalization, severe-emphysema lung regulatory-element burden in COPD
GWAS-gene loci, and regulatory classification of the frozen Section 2 tag plus
LD-proxy candidate set. Literature synthesis and normal-lung or cell-type
sensitivity analyses are outside this computational report.

## Summary

The frozen candidate set contains 15,389 unique reference records: 476
genome-wide-significant (GWS) tag-only records, 176 records that are both a GWS
tag and an LD proxy, and 14,737 proxy-only records. Of these, 15,386 have
GRCh38 BED coordinates and three unlocalized tags do not. Among BED-eligible
records, 257 (1.67%) overlap a GENCODE v50 coding sequence, 1,956 (12.71%)
overlap either preliminary donor regulatory-element class, 639 (4.15%) overlap
either refined donor class, 2,638 (17.15%) overlap SCREEN, Ensembl, or FANTOM5,
8,131 (52.85%) overlap RepeatMasker, and 21 (0.14%) overlap the ENCODE hg38
blacklist. These categories are mutually nonexclusive.

The donor-specific signal is descriptive rather than causal: the ENCODE
context is one bulk-lung donor with reported severe emphysema, not a COPD
case-control cohort. Similarly, source-wise and ancestry-panel percentages are
annotations of a deliberately ascertained tag-plus-LD set, not enrichment
tests.

## Regulatory-reference preparation (COPD-S3-R001)

Five immutable shared resources were normalized to a six-column, 0-based,
half-open BED representation restricted to `chr1-22`, `chrX`, and `chrY`.

| Reference | Source ID | Release | Canonical intervals retained | Noncanonical excluded | Invalid excluded |
|---|---|---:|---:|---:|---:|
| ENCODE SCREEN cCREs | COPD-SRC-038 | Registry V3 | 1,063,878 | 0 | 0 |
| Ensembl Regulatory Build | COPD-SRC-039 | 116 | 380,818 | 0 | 0 |
| FANTOM5 CAGE enhancers | COPD-SRC-040 | hg38_latest | 63,285 | 0 | 0 |
| UCSC RepeatMasker | COPD-SRC-041 | hg38 | 5,317,286 | 366,404 | 0 |
| ENCODE blacklist | COPD-SRC-042 | hg38 v2 | 636 | 0 | 0 |

Source paths, access dates, source and normalized-file SHA-256 values, and
feature-class counts are in
`COPD-S3-R001_regulatory_reference_inventory.tsv` and
`COPD-S3-R001_regulatory_reference_feature_counts.tsv`.

## Severe-emphysema donor regulatory context (COPD-S3-R002)

The primary epigenomic context is ENCODE donor `ENCDO520EJG`, a 60-year-old
male donor of reported European ancestry whose metadata state "lungs with
severe emphysema" (COPD-SRC-018). Released GRCh38 pseudoreplicated ATAC-seq,
H3K27ac, and H3K27me3 narrowPeak files were available for upper-right,
lower-right, and lower-left lung lobes. ENCODE supplies neither
post-bronchodilator spirometry nor a clinically adjudicated COPD diagnosis for
this donor.

The operational element definitions were:

- preliminary enhancer: an H3K27ac peak;
- preliminary silencer: an H3K27me3 peak;
- refined enhancer: an ATAC peak overlapping an H3K27ac peak from the same
  donor and lobe by at least 1 bp;
- refined silencer: an ATAC peak overlapping an H3K27me3 peak from the same
  donor and lobe by at least 1 bp; and
- donor union: overlapping or book-ended intervals of one definition and
  element type merged across the three lobes.

Preliminary and refined counts use different units: histone peaks and
accessibility peaks, respectively. A broad histone interval can overlap
several ATAC peaks. The donor union contained 77,164 preliminary enhancer,
78,899 preliminary silencer, 71,746 refined enhancer, and 15,699 refined
silencer regions.

Section 2 selected 155 gene symbols. Coordinate analysis required one exact
GENCODE v50 gene record; `Y_RNA`, `GUSBP5`, and `CYP2B7P` were excluded because
their exact symbols were ambiguous. The analyzed sets therefore contain 152
selected and 140 replicated genes, represented as gene bodies plus 100,000 bp
on each side.

| Context | Preliminary enhancers in replicated-locus union | Preliminary silencers | Refined enhancers | Replicated loci with refined enhancer | Refined silencers | Replicated loci with refined silencer |
|---|---:|---:|---:|---:|---:|---:|
| Upper-right lobe | 1,549 | 657 | 2,230 | 128 (91.43%) | 286 | 69 (49.29%) |
| Lower-right lobe | 1,497 | 762 | 2,455 | 126 (90.00%) | 352 | 74 (52.86%) |
| Lower-left lobe | 1,652 | 1,087 | 2,684 | 127 (90.71%) | 435 | 93 (66.43%) |
| Three-lobe donor union | 1,892 | 1,201 | 1,815 | 132 (94.29%) | 282 | 94 (67.14%) |

In the full 152-gene selected set, the donor union contained 2,027 preliminary
enhancer, 1,241 preliminary silencer, 1,953 refined enhancer, and 298 refined
silencer regions. At least one refined enhancer occurred in 143 loci (94.08%)
and at least one refined silencer in 98 (64.47%). Complete lobe, union, and
per-gene results are in the R002 tables.

## Frozen candidate-set provenance (COPD-S3-R003)

`COPD-S2-R006E_candidate_variants_grch38.tsv.gz` contains 15,389 unique
candidate reference records derived from the GWAS Catalog (COPD-SRC-035) and
the ancestry-aware 1000 Genomes Phase 3 GRCh38 LD expansion (COPD-SRC-037).
The record-level roles deliberately overlap:
652 records represent a GWS tag and 14,913 represent an ancestry-matched LD
proxy, with 176 records in both roles. Expanding the tag records by
`gws_tag_ids` recovers all 660 unique Catalog tag identifiers; eight pairs of
tag labels collapse to the same reference record. Candidate origin, which is
exclusive, is distributed as follows.

| Candidate origin | Records | BED-eligible |
|---|---:|---:|
| GWS tag only | 476 | 473 |
| GWS tag and LD proxy | 176 | 176 |
| LD proxy only | 14,737 | 14,737 |
| **Total** | **15,389** | **15,386** |

The set comprises 13,764 SNVs, 1,540 indels or complex alleles, and 85
unresolved tags. Eighty-two unresolved tags have usable reported positions and
were annotated positionally despite failing local reference matching. The
three tags without BED-eligible coordinates are retained in every primary
table as `bed_ineligible`, with feature flags set to
`NA_not_bed_eligible`:

- `unmatched_tag:15q25.1`
- `unmatched_tag:rs139284640`
- `unmatched_tag:rs751872749`

LD proxy records retain their AFR, AMR, EAS, EUR, and SAS panel provenance from
Section 2. Panel membership is nonexclusive because one variant may be a proxy
in more than one ancestry panel. The tag-supported and tag-polymorphic panel
fields are retained separately and are not treated as proxy provenance.

## Candidate-classification method

GENCODE v50 CDS records were converted from one-based inclusive GTF
coordinates to GRCh38 BED and deduplicated at the gene-aware interval level,
yielding 387,360 canonical-chromosome intervals (COPD-SRC-043). Each eligible
candidate BED interval starts at `position_grch38 - 1` and spans
`max(1, length(ref))` bases. bedtools 2.31.1 required at least 1 bp of overlap
with each of the following:

1. GENCODE v50 CDS;
2. donor-union preliminary enhancer and silencer regions;
3. donor-union refined enhancer and silencer regions;
4. SCREEN V3 cCREs;
5. Ensembl Regulatory Build 116 features;
6. FANTOM5 CAGE enhancers;
7. UCSC hg38 RepeatMasker intervals; and
8. the ENCODE hg38 blacklist v2.

All original Section 2 columns are preserved. Per-source flags, hit counts,
feature identifiers, feature types, and source details are followed by derived
union flags and three explicitly defined exclusive classifications. Normalized
raw intersection tables retain every source hit rather than only boolean
flags.

## Mutually nonexclusive candidate annotations

Percentages below use the 15,386 BED-eligible records as the denominator. An
eligible record may occur in any number of rows.

| Annotation | Overlapping records | Percent eligible |
|---|---:|---:|
| GENCODE v50 CDS | 257 | 1.670% |
| Preliminary donor enhancer | 1,377 | 8.950% |
| Preliminary donor silencer | 612 | 3.978% |
| Either preliminary donor class | 1,956 | 12.713% |
| Refined donor enhancer | 588 | 3.822% |
| Refined donor silencer | 86 | 0.559% |
| Either refined donor class | 639 | 4.153% |
| SCREEN V3 cCRE | 2,247 | 14.604% |
| Ensembl Regulatory Build 116 | 1,276 | 8.293% |
| FANTOM5 CAGE enhancer | 162 | 1.053% |
| SCREEN, Ensembl, or FANTOM5 | 2,638 | 17.145% |
| RepeatMasker | 8,131 | 52.847% |
| ENCODE blacklist | 21 | 0.136% |
| Any tested biological annotation, including repeat | 10,649 | 69.212% |
| No tested biological annotation | 4,737 | 30.788% |

There are 33 eligible records called both a preliminary enhancer and
preliminary silencer, and 35 called both a refined enhancer and refined
silencer. The separate flags preserve those overlaps. Blacklist status is a
QC flag and is not included in the `any_tested_biological_annotation` union.

### Tag and proxy reporting units

The two record-level role strata overlap by the 176 records that are both a
tag and a proxy, so their counts must not be added. The tag-identifier stratum
expands the 652 tag records to the 660 original GWS identifiers and is the
appropriate denominator for statements specifically about Catalog tag labels.

| Unit | Total | BED-eligible | Coding CDS | Refined enhancer | Refined silencer | SCREEN/Ensembl/FANTOM5 | RepeatMasker | Blacklist |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| GWS tag records | 652 | 649 | 49 (7.55%) | 48 (7.40%) | 8 (1.23%) | 168 (25.89%) | 258 (39.75%) | 1 (0.15%) |
| GWS tag identifiers | 660 | 657 | 53 (8.07%) | 48 (7.31%) | 9 (1.37%) | 173 (26.33%) | 259 (39.42%) | 1 (0.15%) |
| LD proxy records | 14,913 | 14,913 | 228 (1.53%) | 548 (3.67%) | 80 (0.54%) | 2,513 (16.85%) | 7,944 (53.27%) | 20 (0.13%) |

### Ancestry-panel provenance of proxy records

Every panel stratum below contains LD proxy records and is BED-eligible.
Percentages are within panel. A proxy can occur in multiple rows, so neither
record counts nor percentages should be summed across panels or interpreted
as independent ancestry-specific enrichment.

| LD panel | Proxy records | Coding CDS | Refined enhancer | Refined silencer | SCREEN/Ensembl/FANTOM5 | RepeatMasker |
|---|---:|---:|---:|---:|---:|---:|
| AFR | 1,819 | 33 (1.81%) | 61 (3.35%) | 9 (0.49%) | 271 (14.90%) | 992 (54.54%) |
| AMR | 4,515 | 72 (1.59%) | 171 (3.79%) | 16 (0.35%) | 758 (16.79%) | 2,467 (54.64%) |
| EAS | 5,315 | 83 (1.56%) | 177 (3.33%) | 22 (0.41%) | 868 (16.33%) | 2,890 (54.37%) |
| EUR | 12,697 | 175 (1.38%) | 464 (3.65%) | 65 (0.51%) | 2,138 (16.84%) | 6,790 (53.48%) |
| SAS | 722 | 13 (1.80%) | 37 (5.12%) | 1 (0.14%) | 156 (21.61%) | 424 (58.73%) |

## Mutually exclusive reporting hierarchies

Exclusive classes are deterministic reporting summaries, not assertions that
one annotation is biologically more important than another. Each hierarchy
first separates `bed_ineligible`, then assigns every eligible record to the
first matching class:

- `framework_preliminary`: coding CDS > preliminary enhancer > preliminary
  silencer > RepeatMasker > other;
- `framework_refined`: coding CDS > refined enhancer > refined silencer >
  RepeatMasker > other; and
- `comprehensive`: coding CDS > refined enhancer > refined silencer >
  preliminary enhancer > preliminary silencer > other known regulatory
  (SCREEN, Ensembl, or FANTOM5) > RepeatMasker > other.

The overall comprehensive result is:

| Exclusive class | Records | Percent of BED-eligible records |
|---|---:|---:|
| BED-ineligible | 3 | not applicable |
| Coding CDS | 257 | 1.670% |
| Refined enhancer | 561 | 3.646% |
| Refined silencer | 46 | 0.299% |
| Preliminary enhancer | 835 | 5.427% |
| Preliminary silencer | 518 | 3.367% |
| Other known regulatory | 1,515 | 9.847% |
| Repetitive only in hierarchy | 6,917 | 44.956% |
| Other | 4,737 | 30.788% |

The category counts sum to all 15,389 records. The complete exclusive table
also provides origin, role, variant-class, reference-status, proxy-panel,
tag-supported-panel, tag-polymorphic-panel, and 660-tag-identifier strata.

## Reproducibility and validation

From the repository root:

```bash
python diseases/COPD/03_regulatory_landscape/scripts/02_prepare_regulatory_references.py
python diseases/COPD/03_regulatory_landscape/scripts/03_gene_locus_regulatory_burden.py
bash diseases/COPD/03_regulatory_landscape/scripts/04_validate_regulatory_burden.sh
python diseases/COPD/03_regulatory_landscape/scripts/05_classify_candidate_variants.py
bash diseases/COPD/03_regulatory_landscape/scripts/06_validate_candidate_classification.sh
```

The production scripts locate repository inputs relative to their own paths.
The R003 manifest records input, script, feature, raw-intersection, and output
hashes; software versions; exact replay command; denominators; and all three
unlocalized identifiers. The independent validator recomputed all ten overlap
counts directly with bedtools, verified exact preservation of the 15,389 input
records and their original columns, checked the 660-identifier tag expansion,
confirmed the explicit ineligible states, validated all 81 exclusive
stratum-by-hierarchy partitions, and matched manifest output hashes. All checks
passed.

## Interpretation and limitations

These are physical interval overlaps, not functional validation, target-gene
assignment, fine-mapping, or causal classification. The 82 positionally
annotated reference-unmatched tags have usable reported coordinates but are
not sequence-normalized against the local reference. The three remaining tags
cannot be annotated without resolvable GRCh38 coordinates.

The donor-specific annotations come from one severe-emphysema bulk-lung donor
and cannot separate COPD effects from individual, anatomical, exposure,
treatment, comorbidity, or cell-composition effects. H3K27ac is not
enhancer-specific, and H3K27me3 is a broad repressive mark rather than direct
evidence of silencer function. The generic SCREEN, Ensembl, FANTOM5,
RepeatMasker, and blacklist resources are not COPD-specific. The high repeat
fraction is reported descriptively and must not be interpreted as disease
enrichment. The 21 blacklist-overlapping candidates should be flagged in
downstream analyses. Finally, ancestry-panel percentages reflect different,
overlapping LD-derived candidate sets and were not tested for comparative
enrichment.
