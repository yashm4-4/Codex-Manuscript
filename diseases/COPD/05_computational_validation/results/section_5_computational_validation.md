# Section 5: Computational Validation

## Scope and inputs

The analysis evaluated all 337 Section 4 predicted-causal regulatory candidate records without changing their pre-existing priority order. Three evidence types were kept separate: exact variant-level GTEx v10 Lung significant cis-eQTLs (COPD-S5-R001), MPRAbase v4.9.3 assayed-element evidence after audited genome-build conversion (COPD-S5-R002), and HGNC-resolved, exact-node Open Targets COPD gene associations (COPD-S5-R003). Resource and biobank availability was audited in COPD-S5-R004.

## Exact GTEx v10 Lung cis-eQTL evidence

Exact matching required the same GRCh38 chromosome, 1-based position, REF, and ALT. 185 of 337 predicted-causal candidates (54.90%) matched at least one row in the GTEx v10 Lung significant-pairs file. For 72 candidates, at least one significant eGene also matched an HGNC-resolved GWAS mapped gene or selected gene locus. A significant eGene matched at least one broad R006 provisional target for 144 of 185 GTEx-supported candidates with target hypotheses (77.84%). R006 includes nearest-TSS and all-TSS-within-100-kb hypotheses, so this concordance is not a chromatin-contact assignment. These are molecular association results from bulk, non-diseased lung and do not by themselves establish COPD causality or a specific causal lung cell type.

## MPRAbase evidence

Candidate sites were represented as one-base GRCh38 BED intervals, lifted to hg19 with the UCSC chain, and required to map uniquely and return to the original GRCh38 base on reverse liftover. 1 candidate (0.30%) had at least one MPRAbase element match. Coordinate containment was observed for 1 candidate and an exact candidate rsID in element metadata was observed for 0. The sole match was 11:13140768:T:C (rs11022677), covered by 2 elements in Klein_MPRA_HepG2 (HepG2; PMID 33046894). Coordinate containment shows that an assayed sequence covered the candidate base. It does not demonstrate that the REF and ALT alleles were both tested, and the observed HepG2 context is not a COPD lung-cell assay. Raw score scales cannot be compared across heterogeneous MPRA studies without study-specific calibration.

## Gene catalog context

At least one HGNC-resolved candidate-linked gene had an exact `MONDO_0005002` Open Targets association for 297 candidates (88.13%). The tested gene universe includes broad R006 proximity hypotheses, including all TSSs within 100 kb, so this high coverage is descriptive context rather than evidence of enrichment. It is not independent evidence for the candidate variant or proof that a candidate element regulates that gene. All Section 4 target assignments remain explicitly labeled provisional.

## Integrated result

186 candidates (55.19%) had exact GTEx or MPRAbase regional/identifier evidence; 0 had both evidence types. The remaining 151 candidates had no support in these two primary public queries. That category is not a negative causal call because available resources do not cover all lung cell types, disease states, perturbations, alleles, or regulatory assay designs.

QTLbase was not used as primary evidence. Its public endpoint is reserved for exact candidate-rsID requests with returned-count auditing; GWAS source focal-tag rsIDs are not substituted for LD-proxy identifiers. No framework-named biobank had local exact-candidate COPD statistics. Controlled-access and missing-summary resources are therefore recorded as validation gaps rather than null replications.

## Reproducible outputs

- `COPD-S5-R001_GTEx_v10_Lung_exact_significant_pairs.tsv.gz`: exact significant variant-gene pairs.
- `COPD-S5-R002_MPRAbase_v4_9_3_candidate_elements.tsv.gz`: candidate-element matches with assay metadata and evidence scope.
- `COPD-S5-R003_candidate_gene_context.tsv.gz`: HGNC resolution and exact COPD Open Targets context.
- `COPD-S5-R004_resource_access_audit.tsv` and `COPD-S5-R004_biobank_access_audit.tsv`: public and controlled-access audit.
- `COPD-S5-R005_integrated_candidate_validation.tsv`: one row per predicted-causal candidate, in Section 4 priority order.
- `COPD-S5-R005_integrated_summary.tsv`: denominators and evidence-overlap counts used above.
