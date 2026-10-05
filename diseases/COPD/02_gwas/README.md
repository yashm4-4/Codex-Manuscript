# Section 2: GWAS

Status: Section 2 complete (2026-10-01): core COPD GWAS analysis,
ancestry-aware 1000 Genomes LD expansion, replicated-gene locus properties,
conservation, and non-circular Open Targets disease links. No live external API
call is required; the disease-link analysis uses a dated local snapshot.

The core set is restricted to GWAS Catalog studies and associations mapped
directly to COPD ontology terms. Related respiratory phenotypes are exported
separately for audit and are not included in core counts.

Scientific outputs use identifiers beginning `COPD-S2-`.

## Reproduce the core analysis

```bash
python scripts/02_core_gwas_analysis.py
```

The script derives the Section 2 and repository roots from its own location,
so it can be launched from any working directory. It is offline: inputs are
the step-01 COPD extracts plus the repository-local GENCODE v50 and 1000
Genomes GRCh38 references. It makes no web/API calls and does no LD expansion.

The exact replay command, working directory, software/reference versions,
input checksums, script checksum, and output row counts are recorded in
`logs/COPD-S2-core_analysis_manifest.tsv`; the complete run log is
`logs/COPD-S2-core_analysis.log`.

## Reproduce locus-property and other-disease analyses

The replicated-gene analysis uses local GENCODE v50 and UCSC hg38 100-way
conservation tracks. UCSC command-line tools are needed to score the bigWigs.
The disease-link analysis uses the local Open Targets Platform snapshot
downloaded from the official `latest` endpoint on 2026-09-21 and subsequently
resolved byte-for-byte to the stable Open Targets Platform 26.06 release.

```bash
source /etc/profile.d/modules.sh
module load ucsc/503
python scripts/06_locus_properties.py --permutations 100000 --seed 20261001
python scripts/07_other_disease_links.py
```

Both scripts derive repository-local defaults from their own paths and accept
command-line overrides for inputs and output directories. Exact software,
reference, seed, checksum, and row-count provenance is in
`logs/COPD-S2-locus_properties_manifest.tsv` and
`logs/COPD-S2-other_disease_links_manifest.tsv`.

## Reproduce the LD expansion

The shared reference is the 1000 Genomes Phase 3 phased biallelic SNV/indel
release remapped to GRCh38 (20190312). The GWAS Catalog coordinates are also
GRCh38, so no liftOver is performed and no build-conversion loss is introduced.

```bash
python scripts/03_prepare_ld.py
source /etc/profile.d/modules.sh
module load bcftools/1.23 plink/1.9.0-beta4.4
python scripts/04_run_ld.py --jobs 5 --window-bp 500000 --r2 0.8
python scripts/05_summarize_ld.py --window-bp 500000 --r2 0.8
```

Each tag is assigned only to 1000 Genomes super-populations represented in the
initial-stage ancestry rows for the study or studies in which that tag is
genome-wide significant. The explicit Catalog-to-panel mapping is European to
EUR, East Asian to EAS, African American/Afro-Caribbean or African unspecified
to AFR, Hispanic/Latin American to AMR, and South Asian to SAS. Composite
Catalog labels contribute each named panel. There is no default-EUR fallback.

Focal matching is deliberately strict. A tag must have a valid Catalog GRCh38
coordinate and either (i) a unique reference record with a reported REF/ALT
allele match, (ii) a unique record at that coordinate when no sequence allele
was reported, or (iii) a unique indel one base to either side whose VCF REF or
ALT matches the reported allele, reflecting anchor normalization. Neighboring
SNVs and allele-discordant or ambiguous records are not accepted. All failures
remain in the match audit and candidate table.

The VCF header was checked on every used chromosome. Of the 2,504 samples in
the standard Phase 3 population panel, one AFR metadata sample (`NA18498`) is
absent from this GRCh38 VCF release. The effective panels are therefore AFR
660, AMR 347, EAS 504, EUR 503, and SAS 489 samples. The exact exclusion is in
`data/ld_work/reference_sample_exclusions.tsv`.

## Stable result sets

- `COPD-S2-R001_*`: core study, ancestry, and sample summaries. Parsed
  case/control and ancestry-individual totals are descriptive sums across
  overlapping GWAS cohorts, not counts of unique people.
- `COPD-S2-R002_*`: normalized unique tag variants, local sequence-variant
  types, GWAS Catalog consequences, and the genome-wide-significant subset
  (`P <= 5e-8`). Sequence type is assigned only when a variant is found at its
  exact Catalog GRCh38 coordinate in the local 1000 Genomes panel; unresolved
  variants are retained explicitly.
- `COPD-S2-R003A_*`: genes mapped to genome-wide-significant tag variants.
  The multi-study definition is at least two distinct GWAS Catalog study
  accessions; the tables also report distinct-publication counts.
- `COPD-S2-R003B_*`: all reported numeric effects and the top 15 per inferred
  measure class. The Catalog combines OR and beta in one field, so unit text
  identifies beta/other effects; remaining effects are labelled inferred ORs
  and ranked symmetrically as `max(OR, 1/OR)`.
- `COPD-S2-R004_*`: coding, splice, noncoding, and unknown consequence
  summaries for unique tag variants and per-gene profiles.
- `COPD-S2-R005_locus_table_gencode_v50.tsv`: GENCODE v50 GRCh38 gene bodies
  and nearest non-overlapping protein-coding/lncRNA flanking genes for the
  replicated genes plus genes in the top-effect table.
- `COPD-S2-R006A_focal_match_audit.tsv`: all 660 significant tag variants,
  including strict match class, coordinate shift, alleles, and every explicit
  exclusion.
- `COPD-S2-R006B_focal_panel_audit.tsv`: all 1,245 ancestry-supported
  tag-panel assignments, panel MAF, and polymorphic/monomorphic/not-run status.
- `COPD-S2-R006C_ld_pairs.tsv.gz`: validated non-self tag-proxy links at
  `r2 >= 0.8` and absolute distance no greater than 500 kb, with panel and
  GRCh38 allele provenance.
- `COPD-S2-R006D_ld_expansion_by_focal.tsv`: proxy count and inclusive block
  extent for every tag-panel assignment, including zero-proxy and unresolved
  rows.
- `COPD-S2-R006E_candidate_variants_grch38.tsv.gz` and `.bed`: union of the
  original GWS tags and LD proxies for downstream regulatory analysis. The TSV
  retains the three coordinate-unresolved tags that cannot appear in BED.
- `COPD-S2-R006F_summary.tsv` and `_panel_summary.tsv`: overall and
  population-specific counts. Intermediate chromosome-panel VCF, PLINK, MAF,
  and LD files are retained under `data/ld_work/jobs/` for auditability.
- `COPD-S2-R007A_*`: the 140 unambiguous replicated genes, gene-biotype-matched
  GENCODE v50 background, gene-body and gene-plus-100-kb length tests, and
  per-gene conservation metrics. The three ambiguous replicated names are
  retained in `COPD-S2-R007_excluded_ambiguous_replicated_genes.tsv`.
- `COPD-S2-R007B_*`: chromosome-level observed and matched-null counts, the
  global Monte Carlo goodness-of-fit result, and BH correction across all 24
  canonical chromosomes.
- `COPD-S2-R007C_conservation_tests.tsv`: phyloP100way and phastCons100way
  gene-body mean0 comparisons, with four-test BH correction.
- `COPD-S2-R008A_*` through `R008E_*`: per-gene, compact top, complete, COPD
  ontology/lexical exclusion-audit, and summary tables for direct Open Targets
  links to other MONDO diseases.
- `section_2_gwas.md`: integrated Section 2 interpretation, exact tests,
  effect sizes, limitations, and downstream handoff.

## Core counts

- 104 core COPD study accessions from 42 publications; 45 accessions have
  curated association rows and 69 report full summary statistics.
- 995 association rows and 760 normalized unique tag variants.
- 827 associations at `P <= 5e-8`, representing 660 unique tag variants.
- Local 1000 Genomes lookup resolved 644/760 variants: 637 SNVs and 7 indels;
  115 were unresolved and one record was a cytogenetic region rather than a
  single sequence variant.
- 143 mapped genes were genome-wide significant in at least two study
  accessions; 117 of these also occurred across at least two publications.
- Of 660 significant tag variants, 582 had a Catalog consequence: 46 coding,
  3 splice, and 533 noncoding (91.6% noncoding among annotated variants); 78
  lacked consequence annotation.
- 588 genes were mapped to significant variants: 535 had only noncoding tag
  variants, 51 had a coding tag variant, and 2 had splice but no coding tag
  variant.
- All 155 selected locus-table genes matched GENCODE v50 by exact gene name.

## LD expansion counts

- All 660 significant tags received at least one study-supported population
  assignment: EUR 650, AFR 193, EAS 190, AMR 173, and SAS 39 (1,245 total).
- 575 tags matched the local reference: 539 by exact coordinate plus reported
  allele, 7 by a unique exact-position record when no allele was reported, and
  29 by an explicitly flagged adjacent indel-anchor match. The remaining 85
  tags are retained but were not LD-expanded.
- Eight pairs of Catalog identifiers resolve to the same VCF record, leaving
  567 unique matched reference records.
- Of 1,132 matched tag-panel assignments, 1,108 were polymorphic and 24 were
  monomorphic in the assigned panel. The 113 assignments belonging to the 85
  unmatched tags were not run.
- PLINK emitted 34,208 rows, including 1,088 focal self-pairs. After removing
  self-pairs and expanding shared VCF records back to their Catalog tag labels,
  33,378 tag-specific links remain.
- The union contains 14,913 unique LD proxy records: 13,403 SNVs and 1,510
  indel/complex records. Combining these with all tag records yields 15,389
  candidate records, of which 15,386 have a usable GRCh38 coordinate in BED.

## LD limitations

The panel mapping reflects the ancestry metadata attached to each discovery
study, not participant-level ancestry or association-stratum labels, which are
usually unavailable in the Catalog association rows. A composite or
multi-ancestry study therefore contributes all locally represented panels.
The 1000 Genomes super-populations are proxies for LD in the reported cohorts,
not exact cohort matches. The seven positional-only focal matches have lower
allelic confidence than coordinate-plus-allele matches. LD proxies identify
correlated variants, not causal variants, and variants absent from this Phase 3
release cannot be expanded even when the Catalog supplies a coordinate.

## Replicated-gene locus-property results

- Of 143 replicated Catalog-mapped gene names, 140 map unambiguously to one
  GENCODE v50 record. Y_RNA, GUSBP5, and CYP2B7P are retained as explicit
  ambiguous exclusions.
- The background contains 70,016 unique-name GENCODE genes on canonical
  chromosomes with the same eight biotypes. One hundred thousand random draws
  preserve the exact observed biotype counts.
- Mean gene-body length is 149,940.3 bp versus a matched-null expectation of
  54,756.506 bp (ratio 2.73831; two-sided empirical `P = q = 1/100001`). Mean
  gene-plus-100-kb locus length is 349,940.3 versus 254,666.104 bp (ratio
  1.37411; `P = q = 1/100001`). Median tests agree.
- The chromosome distribution is non-random under the matched null (Monte
  Carlo Pearson statistic 93.249672; `P = 1/100001`). Only chromosome 6
  survives BH correction across 24 chromosomes: 27 observed versus 7.49505
  expected genes, ratio 3.60238, `q = 0.000240`.
- Conservation point estimates are lower than the matched expectation, but no
  phyloP100way or phastCons100way mean/median test survives four-test BH
  correction (smallest `q = 0.111506`). Mean track coverage across the 140 gene
  bodies is 0.996372; the minimum is 0.884328.

## Other-disease links

After restricting to MONDO terms and excluding the COPD ontology closure plus
explicit COPD/emphysema/chronic-bronchitis/bronchiectasis label safeguards,
125/140 genes have at least one direct Open Targets link to another disease.
The complete table has 43,636 unique gene-disease pairs across 8,549 diseases;
15,767 pairs have at least two evidence records, covering 111 genes. Eighty-five
genes have at least one link from a predefined human genetic or curated source.
The 32 excluded terms and 271 removed target-disease rows are fully audited.
Open Targets scores rank heterogeneous evidence; they are not P values and do
not establish pleiotropic causality or COPD comorbidity.

## Additional locus-property limitations

The analysis rows are mapped genes rather than independent fine-mapped loci;
neighboring genes from one association region can therefore amplify chromosome
counts. The matched null controls gene biotype, but not expression,
recombination, GC content, or study intensity. The chromosome test addresses
chromosome-level representation, not physical spacing within a chromosome.
Fixed 100-kb flanks are operational windows rather than causal regulatory
domains, and whole-gene conservation averages can mask constrained subregions.
