# Section 2: COPD GWAS and locus characterization

Evidence cutoff: 2026-10-01. Stable outputs are registered as COPD-S2-R001
through COPD-S2-R008. The core phenotype is restricted to GWAS Catalog studies
mapped directly to chronic obstructive pulmonary disease; related respiratory,
smoking, biomarker, and interaction traits remain outside the core set
[COPD-SRC-035].

## GWAS inventory and ancestry

The local GWAS Catalog release contains 104 core COPD study accessions from 42
publications. Forty-five accessions have curated association rows and 69 report
full summary statistics. The analysis includes 995 association rows and 760
normalized unique tag variants. At the prespecified genome-wide threshold of
`P <= 5e-8`, 827 association rows represent 660 unique tag variants.

Initial-stage ancestry metadata are dominated by European cohorts, but also
contain East Asian, African/African-American or Afro-Caribbean, Hispanic/Latin
American, and South Asian entries. The reported participant totals in R001 are
sums across Catalog cohort entries and are not unique-person counts because
cohorts recur across studies. Ancestry is study-level metadata rather than a
variant-stratum assignment, which limits population-specific inference.

## Genes, effects, and coding status

Genome-wide-significant variants map to 588 reported genes. Of these, 143 occur
in at least two distinct GWAS Catalog study accessions and 117 also occur in at
least two publications. Replication here means repeated Catalog mapping, not
independent fine-mapping or functional validation.

Reported effects are separated by inferred measure class. The strongest
inferred odds-ratio entries are rs2286351 (reported effect 6.367955) and
rs7709630 (4.045759), both from the same 2018 study; more widely replicated
signals include rs12914385 near CHRNA3 (1.39) and rs13141641 near HHIP-AS1
(1.39). These rankings reproduce Catalog values and units rather than
harmonized alleles or phenotypes. Beta/other-unit effects are ranked separately,
because their numerical magnitudes cannot be compared with odds ratios.

Among 660 significant tags, 582 have a Catalog consequence: 46 are coding, 3
splice-related, and 533 noncoding. Thus, 91.58% of annotated significant tags
are noncoding; 78 tags lack a consequence annotation. At gene level, 535 of 588
mapped genes have only noncoding significant tags, 51 have at least one coding
tag, and 2 have splice but no coding tag. These are tag-variant annotations and
do not yet include the LD proxies.

## Ancestry-aware LD expansion

The LD expansion uses the local phased, biallelic 1000 Genomes Phase 3 release
remapped to GRCh38 on 2019-03-12 [COPD-SRC-037]. GWAS Catalog coordinates are
also GRCh38, so no liftOver is used. Each significant tag is assigned only to
super-populations represented in its discovery-stage study metadata. Across
1,245 tag-panel assignments, the counts are EUR 650, AFR 193, EAS 190, AMR 173,
and SAS 39.

Of 660 tags, 575 match the reference under the documented strict rules: 539 by
exact coordinate and reported allele, 7 by a unique exact-position record when
no allele was reported, and 29 by an explicitly flagged adjacent indel-anchor
match. All 85 unresolved tags remain in the audit and candidate table. Among
1,132 matched tag-panel assignments, 1,108 are polymorphic and 24 monomorphic.

At `r2 >= 0.8` within an absolute distance of 500 kb, the validated expansion
contains 33,378 tag-specific non-self links and 14,913 unique proxy records
(13,403 SNVs and 1,510 indel/complex records). The frozen union contains 15,389
candidate records; 15,386 have valid GRCh38 coordinates and are available in
BED for downstream regulatory classification.

## Replicated-gene analysis set and comparator

The locus-property analysis starts from the 143 replicated mapped gene names.
Three names are ambiguous in GENCODE v50 and are retained in an exclusion
table: Y_RNA maps to 756 records, while GUSBP5 and CYP2B7P each map to two.
The remaining 140 genes resolve to one GENCODE v50 GRCh38 record
[COPD-SRC-068].

The comparison pool comprises 70,016 non-COPD GENCODE v50 genes on chromosomes
1-22, X, or Y with a unique gene name and one of the eight gene biotypes seen
in the COPD set. Every null draw samples without replacement within exact gene
biotype, preserving the observed composition (98 protein-coding, 26 lncRNA,
and 16 genes across six other biotypes). The operational locus is the gene body
plus 100 kb on each side, clipped to GRCh38 chromosome bounds. Tests use 100,000
Monte Carlo draws and seed 20261001. Two-sided empirical P values use an add-one
correction; BH correction is applied within each four-test metric family.

## Are COPD genes or operational loci longer?

Yes. Both gene bodies and gene-plus-100-kb operational loci are longer than the
gene-biotype-matched genome-wide expectation.

| Metric | Statistic | COPD observed (bp) | Null expectation (bp) | Ratio | Null-standardized difference | Empirical P | BH q |
|---|---:|---:|---:|---:|---:|---:|---:|
| Gene body | Mean | 149,940.300 | 54,756.506 | 2.73831 | 9.75748 | 1/100001 | 1/100001 |
| Gene body | Median | 67,274.000 | 18,230.427 | 3.69020 | 14.77128 | 1/100001 | 1/100001 |
| Gene + 100 kb per side | Mean | 349,940.300 | 254,666.104 | 1.37411 | 9.76238 | 1/100001 | 1/100001 |
| Gene + 100 kb per side | Median | 267,274.000 | 218,185.937 | 1.22498 | 14.78395 | 1/100001 | 1/100001 |

No one of the 100,000 matched draws was at least as extreme, so `1/100001` is
the resolution floor rather than an exact tail probability below that value.
The 100-kb flank attenuates, but does not remove, the length difference.

## Do replicated genes cluster by chromosome?

A gene-biotype-stratified Monte Carlo Pearson goodness-of-fit test rejects the
expected chromosome distribution (statistic 93.249672; 24 chromosomes;
sqrt(statistic/140) = 0.816131; empirical `P = 1/100001`). The expected counts
are the chromosome frequencies produced by the same exact-biotype matched null
draws used for the length analysis, rather than uniform counts or chromosome
length alone.

After BH correction across 24 chromosome tests, only chromosome 6 remains
significant: 27 observed genes versus 7.49505 expected, enrichment ratio
3.60238, two-sided empirical `P = 1/100001`, and `q = 0.000240`. Nominal signals
on chromosome 4 (12 versus 5.80732; `P = 0.013850`), chromosome X (0 versus
5.75646; `P = 0.015660`), and chromosome 15 (9 versus 4.33914;
`P = 0.029660`) do not survive correction. This is chromosome-level
overrepresentation, not a test of within-chromosome physical spacing.

## Are replicated gene bodies unusually conserved?

The local UCSC hg38 100-vertebrate phyloP and phastCons bigWigs are available
and were scored over the 140 gene bodies [COPD-SRC-069]. `mean0` averages over
the full gene body and assigns zero to positions without a bigWig value; track
coverage is reported separately. Mean coverage is 0.996372 and the minimum is
0.884328 for both tracks.

| Metric | Statistic | COPD observed | Null expectation | Ratio | Null-standardized difference | Empirical P | BH q |
|---|---:|---:|---:|---:|---:|---:|---:|
| phyloP100way mean0 | Mean | 0.254147 | 0.313007 | 0.81195 | -1.36289 | 0.169008 | 0.169008 |
| phyloP100way mean0 | Median | 0.111256 | 0.148835 | 0.74751 | -1.73264 | 0.083629 | 0.111506 |
| phastCons100way mean0 | Mean | 0.148951 | 0.169736 | 0.87755 | -1.86604 | 0.061999 | 0.111506 |
| phastCons100way mean0 | Median | 0.098077 | 0.110194 | 0.89004 | -1.93965 | 0.052549 | 0.111506 |

All point estimates are lower than the matched expectation, but none survives
BH correction. The supported conclusion is therefore no detectable difference
in gene-body conservation under this analysis, not that COPD genes are more or
less evolutionarily constrained.

## Links to other diseases

Other-disease links use the local Open Targets Platform direct target-disease
association snapshot downloaded from the official `latest` endpoint on
2026-09-21 [COPD-SRC-070]. Only positive MONDO disease associations are
retained. COPD itself, all Open Targets COPD ancestors and descendants, and an
additional case-insensitive label safeguard for COPD, chronic obstructive,
emphysema, chronic bronchitis, and bronchiectasis are excluded. The safeguard
removes compound pulmonary terms that fall outside the formal COPD ontology
closure; the full 32-term exclusion audit is R008D.

Of the 140 replicated unambiguous genes, 125 have at least one non-circular
direct MONDO association, 111 have at least one link supported by two or more
evidence records, and 85 have a link from a predefined human genetic or curated
source. The full table contains 43,636 unique gene-disease pairs spanning 8,549
MONDO diseases; 15,767 pairs have at least two evidence records. Examples among
the highest-scoring retained links include GM2A with Tay-Sachs disease AB
variant (score 0.811065; 211 evidence records), TET2 with myelodysplastic
syndrome (0.800389; 942), GLIS3 with neonatal diabetes and congenital
hypothyroidism (0.788320; 446), and SMAD3 with aneurysm-osteoarthritis syndrome
(0.786349; 393).

Open Targets association scores combine heterogeneous evidence and are ranking
scores, not P values. These links establish that mapped COPD genes have evidence
in other disease contexts; they do not demonstrate COPD comorbidity, shared
causality, or that the GWAS-mapped gene mediates its COPD locus. R008B supplies
up to five compact links per gene, requiring at least two evidence records;
R008C retains every eligible pair and its datasource-level provenance.

## Reproducibility and limitations

The complete replay commands, local resource paths and sizes, software
versions, random seed, script checksums, input checksums, output checksums, and
row counts are stored in the two R007/R008 manifests under `../logs/`. From the
repository root, the analyses can be rerun with:

```bash
source /etc/profile.d/modules.sh
module load ucsc/503
python diseases/COPD/02_gwas/scripts/06_locus_properties.py --permutations 100000 --seed 20261001
python diseases/COPD/02_gwas/scripts/07_other_disease_links.py
```

The 140 rows are mapped genes, not statistically independent fine-mapped loci;
multiple neighboring genes can represent one association region, particularly
on chromosome 6. Gene-name ambiguity excluded three replicated names, and the
matched background controls biotype but not expression, recombination, GC
content, ascertainment, or publication intensity. Operational loci can overlap
and their fixed flanks do not define causal regulatory domains. Conservation is
an average over whole gene bodies and can obscure exon- or element-level
constraint. Open Targets coverage varies by gene and datasource and includes
literature, model, somatic, and genetic evidence. These limitations require
candidate-level functional follow-up rather than causal interpretation.
