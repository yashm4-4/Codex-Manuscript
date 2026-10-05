# COPD Section 3 evidence synthesis

Evidence cutoff: 2026-10-01. This synthesis is deliberately separate from the
candidate-variant overlap computation. It evaluates published human COPD
expression and regulatory evidence and does not treat physical overlap with an
annotation as functional validation.

## Traceable outputs

- `COPD-S3-R004_expression_evidence.tsv` distinguishes cell-state expression,
  cell abundance, bulk composition, phenotype scope, and perturbational
  follow-up.
- `COPD-S3-R004_regulatory_element_evidence.tsv` records element type, lung
  context, target linkage, disease linkage, functional test, and causal
  strength.
- `COPD-S3-R004_functional_variant_evidence.tsv` separates endogenous
  allele tests from episomal allele tests and variant-containing-region
  perturbations.
- `../data/COPD-S3-R004_literature_search_log.tsv` records the searches,
  screening decisions, and important downgrades.
- `../data/COPD-S3-R004_source_rows_proposed.tsv` contains source rows that can
  be appended to the global register after concurrent work is reconciled.

The tables also reuse already registered `COPD-SRC-015` through
`COPD-SRC-017`; they do not duplicate those source records. New literature
source IDs are `COPD-SRC-044` through `COPD-SRC-059`. No global register was
edited.

## Evidence scale used here

The highest tier is an alternate-allele perturbation at the endogenous locus.
The next tier is endogenous deletion or CRISPRi of a regulatory element, which
can establish region-to-gene function but not necessarily identify a causal
nucleotide. Allele-specific reporter and binding assays show sequence-dependent
activity outside the native chromosome. Disease-associated accessibility,
methylation, expression, eQTLs, chromatin contact, and sequence-effect models
are supportive associations unless paired with a perturbation. These categories
are recorded explicitly rather than collapsed into a single "functional"
label.

## Main findings

### Expression changes are strongly compartment and phenotype dependent

The most robust conclusion is not a universal COPD expression signature, but a
set of context-specific programs.

- Large contemporary single-cell studies resolve small-airway and lung-wide
  cellular states and also quantify changes in cell abundance
  (`COPD-SRC-015`, `COPD-SRC-016`). These two signal types are not
  interchangeable.
- In advanced COPD with radiographic emphysema, AT2 cells showed lower
  `NUPR1`, capillary endothelium showed increased `CXCL12`-centered signaling,
  and high-metallothionein macrophages were enriched (`COPD-SRC-057`). The
  first two are within-state expression findings; the macrophage result is in
  part a state-abundance result. It should be labeled advanced-emphysema
  evidence, not all-COPD evidence.
- A smaller severe-COPD single-cell study localized many changes to monocytes,
  macrophages, and ciliated epithelium and supported altered `IGFBP5` and `QKI`
  protein (`COPD-SRC-056`). Its three-versus-three, age- and smoking-mismatched
  design warrants caution. The paper contains conflicting narrative statements
  about the direction of QKI protein, so the evidence table intentionally does
  not assign one.
- Bronchial and small-airway brushing studies provide epithelial-enriched, not
  pure-cell, profiles. The 98-gene Steiling signature was selected for joint
  association with COPD status and spirometry, and also tracked emphysema
  (`COPD-SRC-054`). The distal-to-proximal airway program was explicitly
  smoking-dependent and was reproduced with EGF (`COPD-SRC-055`). Neither
  should be reported as COPD-exclusive.
- Purified, cultured COPD fibroblasts retained widespread methylation and
  expression differences, supporting a cell-intrinsic component despite the
  small discovery series (`COPD-SRC-050`). COPD-derived differentiated airway
  epithelium likewise retained lower `TET1` and `CDH1`, coupled to enhancer-D
  hypermethylation (`COPD-SRC-049`).
- Bulk lung and BAL studies remain informative but composition-sensitive.
  `EPAS1` was supported by bulk multi-omics, lower lung protein, and endothelial
  knockdown (`COPD-SRC-052`), but the original lung association can reflect
  endothelial abundance as well as regulation. BAL methylation
  (`COPD-SRC-053`) is a mixed immune-cell signal and its expression correlation
  used a separate cohort.

WNT findings illustrate why tissue labels matter. `FZD4` and canonical WNT
repair signaling were reduced in COPD alveolar epithelium, principally an
emphysema/repair result (`COPD-SRC-058`). Canonical WNT activity was increased
in COPD bronchial epithelium and experimentally promoted barrier and
differentiation defects (`COPD-SRC-059`). These are different compartments and
phenotypes, not a single contradictory whole-lung direction.

### Regulatory elements with the strongest target evidence

1. **FAM13A intronic regulatory sequence.** A 606-variant MPRA at the COPD
   locus identified 45 allele-sensitive sequences. Three variants validated in
   conventional reporter assays, and `rs2013701` contacted the `FAM13A`
   promoter and was edited at the endogenous locus. The edit changed `FAM13A`
   expression and cellular proliferation (`COPD-SRC-044`). This is the clearest
   endogenous allele-specific result in the reviewed set.
2. **HHIP upstream enhancer carrying a COPD-risk haplotype.** An approximately
   2.4-kb region about 85 kb upstream of `HHIP`, narrowed to an approximately
   500-bp element, contacted the promoter. The `rs6537296-A` plus
   `rs1542725-C` COPD-risk haplotype lowered reporter activity, and
   `rs1542725-C` bound the Sp3 repressor more strongly (`COPD-SRC-045`). This is
   strong allele-specific reporter and binding evidence, but not an endogenous
   allele edit.
3. **HHIP T2 distal enhancer.** A human epithelial H3K27ac-marked element about
   60 kb from `HHIP` contacted its promoter. Deleting less than 500 bp reduced
   `HHIP`, disrupted local topology, and weakened the TGF-beta/SMAD3 response
   in BEAS-2B, primary bronchial, and iPSC-derived AT2 contexts
   (`COPD-SRC-046`). This is strong element-to-target perturbation. The paper
   explicitly reports that no COPD GWAS variant lies within T2, so it is not a
   functionally tested disease variant.
4. **`rs35421223` region.** A five-locus COPD MPRA identified allele-sensitive
   activity. CRISPRi of the `rs35421223`-containing open, H3K27ac-enriched
   region changed `RUVBL1` and `RAB7A` in 16HBE and primary NHBE cells
   (`COPD-SRC-047`). The target link is endogenous and reproducible, but the
   CRISPRi experiment represses the region rather than comparing endogenous
   alleles.
5. **CDH1 enhancer D.** COPD-derived airway epithelium showed enhancer-D
   hypermethylation, lower RNA-polymerase-II occupancy, and reduced `CDH1`.
   Demethylating treatment restored barrier-related phenotypes and improved
   damage measures in COPD lung slices (`COPD-SRC-049`). The intervention is
   genome-wide, so it does not isolate enhancer D as the sole mediator.

The regions containing `rs1800469` and `rs2241712` regulate `TGFB1` and, for
`rs2241712`, `B9D2` and `TMEM91` by CRISPRi (`COPD-SRC-048`). They are useful
region-to-gene evidence, but the loci are shared asthma/CF/COPD or asthma/COPD
modifiers and the alternate alleles were not tested. They must not be promoted
to COPD-specific causal variants.

### Association-only regulatory evidence remains useful for prioritization

ATAC-seq in nondiseased primary lung cell types placed 250 fine-mapped COPD
variants in open chromatin and prioritized 22 variants by OCR overlap plus
deltaSVM (`COPD-SRC-017`). This is a cell-context prediction resource, not a
COPD case-control chromatin study. Similarly, airway, fibroblast, whole-lung,
and BAL methylation studies identify plausible disease-associated regulatory
programs (`COPD-SRC-050` through `COPD-SRC-053`) but do not establish that an
individual CpG or enhancer causes the paired expression change.

## Integration rules for the Section 3 computational results

- Keep the ENCODE severe-emphysema donor interpretation unchanged: it is one
  donor with reported severe emphysema, not a clinically adjudicated COPD case
  cohort.
- A candidate variant overlapping an enhancer, silencer, OCR, methylated site,
  or predicted motif receives an annotation, not a causal label.
- Promote a variant to "allele-function tested" only when alternate alleles
  were compared. Promote it to "endogenous allele-function tested" only when
  the endogenous nucleotide was edited; in this reviewed set, that highest
  label applies to `rs2013701`.
- CRISPR deletion or CRISPRi supports the perturbed region and measured target.
  It does not prove that a nearby SNP is the functional nucleotide.
- Preserve phenotype qualifiers in downstream tables: `COPD`,
  `COPD-plus-lung-function`, `smoking-dependent`, `advanced emphysema`, and
  `shared lung-disease modifier` are not interchangeable.
- Preserve tissue qualifiers. Airway epithelium, alveolar epithelium,
  fibroblasts, endothelium, macrophages, BAL, and homogenized lung can show
  different or opposing programs.
- For bulk or mixed-cell assays, record cellular composition as an alternative
  explanation. For single-cell data, keep within-state differential expression
  separate from shifts in state abundance.

## Bottom line

The most defensible COPD regulatory chains are
`rs2013701 -> FAM13A`, the `rs6537296/rs1542725 -> HHIP` enhancer haplotype,
the non-variant T2 enhancer `-> HHIP`, and the
`rs35421223-containing region -> RUVBL1/RAB7A`. Only the first includes an
endogenous alternate-allele edit. `CDH1` enhancer-D methylation is strong
disease-altered element evidence with a global epigenetic rescue. Other
epigenomic and transcriptomic results should guide cell-type and target
prioritization while retaining their phenotype, composition, and causal limits.
