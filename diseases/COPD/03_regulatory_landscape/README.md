# Section 3: Regulatory Landscape

Status: the computational analyses are complete through frozen candidate-set
classification (`COPD-S3-R001` through `COPD-S3-R003`), and the independently
traceable literature evidence synthesis is complete as `COPD-S3-R004`. The
computational report is `results/section_3_regulatory_landscape.md`, the
evidence report is `results/section_3_evidence_synthesis.md`, and the earlier
pre-R003 snapshot is retained as
`results/section_3_regulatory_landscape_interim.md`.

## Primary context

The primary regulatory context is ENCODE donor `ENCDO520EJG`, a 60-year-old
male donor of reported European ancestry whose health-status metadata state
"lungs with severe emphysema". Upper-right, lower-right, and lower-left lung
lobes each have released GRCh38 pseudoreplicated ATAC-seq, H3K27ac, and
H3K27me3 narrowPeak files. The experiments are matched by donor and anatomic
site, but the assay-specific biosample accessions differ. ENCODE does not
provide post-bronchodilator spirometry or a clinically adjudicated COPD
diagnosis for this donor, so this is a severe-emphysema regulatory context and
not a COPD case cohort.

## Completed outputs

- `COPD-S3-R001_regulatory_reference_inventory.tsv`: normalized GRCh38 SCREEN
  V3, Ensembl Regulatory Build 116, FANTOM5 enhancer, UCSC RepeatMasker, and
  ENCODE blacklist resources with source and output SHA-256 values.
- `COPD-S3-R001_regulatory_reference_feature_counts.tsv`: interval counts by
  source and feature class.
- `COPD-S3-R002_gene_loci_grch38.tsv`: 152 unambiguous selected COPD GWAS genes
  represented as GENCODE v50 gene bodies plus 100 kb on each side.
- `COPD-S3-R002_locus_exclusions.tsv`: three ambiguous gene symbols excluded
  from coordinate-level analysis.
- `COPD-S3-R002_donor_element_sets.tsv`: lobe-specific and donor-union element
  counts, sources, paths, and checksums.
- `COPD-S3-R002_lobe_union_summary.tsv`: distinct element counts and locus
  coverage for the 152 selected and 140 replicated, unambiguous gene sets.
- `COPD-S3-R002_per_gene_burden.tsv`: element count and count per megabase for
  every analyzed gene, lobe or donor union, definition, and element type.
- `COPD-S3-R002_analysis_manifest.json`: exact definitions, software versions,
  input hashes, replay command, and derived-file metadata.
- `COPD-S3-R003_candidate_classification.tsv.gz`: all 15,389 frozen tag and
  LD-proxy records, including three explicit BED-ineligible tags, with source
  overlaps, details, derived flags, and three exclusive classes.
- `COPD-S3-R003_gws_tag_identifier_classification.tsv.gz`: the corresponding
  660-identifier GWS tag view, preserving tag labels that share a reference
  record.
- `COPD-S3-R003_nonexclusive_summary.tsv`: overlap counts and denominators by
  candidate origin, tag/proxy role, variant class, reference status, ancestry
  panel provenance, and tag identifier.
- `COPD-S3-R003_exclusive_summary.tsv`: complete partitions for preliminary,
  refined, and comprehensive reporting hierarchies.
- `COPD-S3-R003_feature_definitions.tsv` and
  `COPD-S3-R003_hierarchy_definitions.tsv`: exact overlap and precedence rules.
- `COPD-S3-R003_analysis_manifest.json`: input, script, feature,
  raw-intersection, and output hashes plus exact software and replay metadata.
- `COPD-S3-R004_expression_evidence.tsv`: 12 traceable human expression rows
  with explicit phenotype, tissue, composition, and perturbation qualifiers.
- `COPD-S3-R004_regulatory_element_evidence.tsv`: 12 element-to-target rows
  graded from association or prediction through endogenous perturbation.
- `COPD-S3-R004_functional_variant_evidence.tsv`: eight tested-variant rows
  distinguishing endogenous allele editing, allele-specific reporter assays,
  and variant-containing-region CRISPRi.
- `section_3_evidence_synthesis.md`: integrated interpretation of COPD versus
  lung-function or emphysema evidence, association versus perturbation, and
  bulk composition versus cell-intrinsic expression.

The normalized shared references are in `data/regulatory_references/`. The 16
derived severe-emphysema element sets are in
`data/copd_donor_regulatory_elements/`. The normalized GENCODE v50 CDS BED and
the ten source-wise raw candidate-intersection tables are in
`data/candidate_classification/` and `data/candidate_intersections/`. The R004
search audit and exact proposed `COPD-SRC-044` through `COPD-SRC-059` source
rows are in `data/COPD-S3-R004_literature_search_log.tsv` and
`data/COPD-S3-R004_source_rows_proposed.tsv`; the global source register was
not edited during concurrent work.

## Element and locus definitions

- Preliminary enhancer: one H3K27ac narrowPeak.
- Preliminary silencer: one H3K27me3 narrowPeak.
- Refined enhancer: one ATAC narrowPeak with at least 1 bp overlap of an
  H3K27ac narrowPeak from the same donor and lobe.
- Refined silencer: one ATAC narrowPeak with at least 1 bp overlap of an
  H3K27me3 narrowPeak from the same donor and lobe.
- Donor union: overlapping or book-ended intervals of one definition and
  element type merged across all three lobes.
- Gene locus: one unambiguous GENCODE v50 gene body plus 100,000 bp on both
  sides, clipped at coordinate zero, on GRCh38 canonical chromosomes.

Preliminary and refined counts do not use the same counting unit. Preliminary
counts enumerate histone peaks, whereas refined counts enumerate accessibility
peaks. A broad histone interval can overlap several ATAC peaks, so a refined
count can exceed its preliminary count. Donor-union counts enumerate merged
genomic regions rather than the sum of lobe-specific peak calls.

## Reproduction

From the repository root:

```bash
python diseases/COPD/03_regulatory_landscape/scripts/02_prepare_regulatory_references.py
python diseases/COPD/03_regulatory_landscape/scripts/03_gene_locus_regulatory_burden.py
bash diseases/COPD/03_regulatory_landscape/scripts/04_validate_regulatory_burden.sh
python diseases/COPD/03_regulatory_landscape/scripts/05_classify_candidate_variants.py
bash diseases/COPD/03_regulatory_landscape/scripts/06_validate_candidate_classification.sh
```

The production scripts locate inputs relative to their own files and run from
outside the repository. R001 and R002 output tables and compressed BED files
reproduced byte-for-byte in an independent rerun from `/tmp`. The R002
same-lobe refined counts and all ten R003 candidate-overlap counts were
recomputed independently with bedtools 2.31.1. R003 validation also checks
source-row preservation, explicit handling of unlocalized tags, tag-identifier
expansion, hierarchy partition sums, source identifiers, and manifest hashes.

## Deferred analyses

- Repeat the peak-burden analysis in normal adult lung and relevant lung cell
  types as sensitivity analyses.
- Functionally validate or fine-map candidate variants; interval overlap alone
  is not evidence of causal regulatory activity.

## Limitations

The primary data come from one severe-emphysema donor and three anatomical
samples, not a case-control series. Bulk lung cannot assign elements to airway,
alveolar, immune, endothelial, fibroblast, or smooth-muscle cell types.
H3K27ac is not enhancer-specific, and H3K27me3 is a broad repressive mark rather
than direct evidence of silencer activity. No promoter or blacklist exclusion
was applied to these descriptive peak counts. The 100 kb flank is an
operational window and does not establish a target relationship. GWAS Catalog
mapped genes are annotation candidates rather than proven causal genes. The
candidate set is LD-ascertained and ancestry-panel memberships overlap, so raw
panel percentages are descriptive rather than comparative enrichment.
Eighty-two reference-unmatched tags were annotated positionally, while three
tags had no BED-eligible GRCh38 coordinate. Generic reference overlaps are not
lung-specific, and the 21 blacklist-overlapping candidates require downstream
QC.
