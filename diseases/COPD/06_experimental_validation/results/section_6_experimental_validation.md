# Section 6: Experimental and Biological Validation

Result: COPD-S6-R006. Status: design and collaborator handoff complete; wet-lab execution is pending. No experimental result or causal proof is claimed. The computational shortlist is not causal proof.

## Scope

Section 6 consumed the canonical 337-row Section 4 R010 list and preserved its `predicted_causal_priority_rank`. R004-R008 supplied sequence eligibility, interval, population, target, and allele-specific motif detail; the analysis also consumed all 4,779 R006 target-evidence rows and the final, validated Section 5 integrated matrix. All Section 6 tables distinguish observed or computed upstream evidence from proposed future experiments.

## Transparent candidate selection

The starting set contained 337 predicted-causal records. Exact shared linked GWAS tags or source focal tags were joined transitively into 153 LD-redundancy components; 184 records were lower-priority members of a component. At most one candidate was selected per component. The policy first ensured coverage of model context, variant class, severe-emphysema-donor overlap, exact GTEx Lung eQTL evidence, public MPRAbase element evidence, and predicted motif creation/disruption when available. A unique diversity anchor could replace the default highest-priority member of its component; remaining slots were filled with default component representatives by the unchanged Section 4 rank. It did not create or optimize a new composite score.

The final panel has 9 both-model, 2 enhancer-only, and 1 silencer-only candidates. It includes 2 indel/complex candidate(s), 3 candidate(s) overlapping a refined regulatory interval from the severe-emphysema donor, 8 exact GTEx v10 Lung significant eQTL candidate(s), 1 candidate(s) with MPRAbase element evidence, and 7 candidate(s) for which an exact GTEx Lung eGene also occurs in the R006 provisional target set.

| S6 rank | Candidate | Model context | Class | Exact-concordant or nearest protein-coding target | Exact lung eQTL | Leading motif TF names | Selection rule |
|---:|---|---|---|---|---|---|---|
| 1 | `11:62567436:G:C` | both | SNV | EEF1G;EML3;INTS5;ROM1 | yes | EHF;Erg | both_model_prediction |
| 2 | `1:3528722:G:C` | both | SNV | ARHGEF16 | no | BATF;FOS::JUN | priority_fill_after_diversity_anchors |
| 3 | `14:92637384:G:A` | both | SNV | RIN3 | yes | BATF3;FOSL2::JUN | priority_fill_after_diversity_anchors |
| 4 | `6:31009903:G:A` | both | SNV | HCG22;SFTA2;VARS2 | yes | EHF;KLF12 | priority_fill_after_diversity_anchors |
| 5 | `6:27556090:G:A` | both | SNV | ZNF184 | yes | ELF3;PATZ1 | priority_fill_after_diversity_anchors |
| 6 | `16:75478398:G:GC` | both | indel_or_complex | ENSG00000261783 | yes | BATF;FOS::JUN | indel_or_complex |
| 7 | `15:67322629:G:A` | both | SNV | AAGAB;IQCH | yes | Erg;Ikzf3 | priority_fill_after_diversity_anchors |
| 8 | `16:28602644:A:G` | both | SNV | CDC37P1;EIF3C;SULT1A1 | yes | BATF;KLF12 | priority_fill_after_diversity_anchors |
| 9 | `17:40058327:G:GCCCAGAC` | both | indel_or_complex | ORMDL3 | yes | ELF3;PATZ1 | priority_fill_after_diversity_anchors |
| 10 | `6:4577675:T:A` | enhancer_only | SNV | CDYL | no | BATF;FOS::JUN | enhancer_only_prediction |
| 11 | `15:67150258:C:T` | silencer_only | SNV | SMAD3 | no | BATF::JUN;FOSL2 | silencer_only_prediction |
| 12 | `11:13140768:T:C` | enhancer_only | SNV | RASSF10 | no | BATF3;KLF12 | public_MPRAbase_element_overlap |

`COPD-S6-R002_selection_audit.tsv.gz` retains every candidate, component membership, eligibility, coverage strata, and exclusion or inclusion reason. Shared-tag components are an LD-redundancy heuristic tied to the available workflow links; they are not claims that every member is mutually correlated in every ancestry.

## Allele-specific MPRA design

R003 contains 48 insert designs: REF and ALT alleles in forward and reverse-complement orientations for each of 12 variants. Each forward insert has exactly 100 genomic bases on each side of the normalized allele. Insert lengths are 201 bp (n=44), 202 bp (n=2), 208 bp (n=2). Paired FASTA allele identity, shared genomic flanks, insert length, centered allele, and reverse-complement identity are validated computationally. 4 constructs carry a sequence-review flag for GC, homopolymers, or common Type IIS recognition sites. Each design isolates the nominated allele in GRCh38 reference-sequence flanks and does not encode local phased donor haplotypes or nearby linked alleles. Because reporter effects can depend on haplotype background, relevant phased haplotypes should be tested where biological or prior functional evidence warrants it.

These are inserts, not ordering-ready oligonucleotides. Vector, promoter, adapters, barcodes, restriction-site remediation, randomization, and assay-specific positive/negative controls must be chosen with the receiving laboratory. The proposed minimum is 10 independent barcodes per construct and 3 independent biological replicates per selected context, subject to prospective power analysis.

## Endogenous validation plan

The R004 candidate plan is staged so that evidence classes remain separate.

1. Allelic MPRA tests sequence-dependent reporter activity in pre-screened airway epithelium, alveolar AT2, and parenchymal lung fibroblast contexts.
2. CRISPRi tests the endogenous region with at least three non-overlapping guides. A region effect does not identify the causal nucleotide.
3. Base editing is considered only for compatible transitions after PAM, activity-window, and bystander review. Prime editing is the primary proposal for transversions and indels. Genotype-verified edited pools and independent clones, off-target review, and clone-aware controls are required.
4. Target-gene readouts combine exact GTEx Lung eGenes and R006 hypotheses with targeted and transcriptome-wide measurements. Target assignments remain provisional until endogenous perturbation supports them.

Primary COPD-relevant models are donor-derived airway epithelial cultures, primary or iPSC-derived AT2 cells, and low-passage parenchymal lung fibroblasts. 16HBE, BEAS-2B, or MRC5 may support technical optimization but cannot substitute for donor replication. COPD and smoking-history-matched control donors, cell identity, viability, passage or differentiation state, and exposure conditions are explicit covariates and controls.

## TF-first perturbation plan

R004 provides 21 candidate-by-model TF hypotheses. For each model-predicted candidate, it selects the strongest compatible-site creation/disruption first, then a retained site or score-only shift. The proposed sequence is expression gating, allele-paired biochemical binding, chromatin binding in a relevant heterozygous or edited context, TF perturbation with independent reagents and rescue, and a formal allele-by-TF interaction test. Motif similarity and PWM changes are not TF binding evidence or TF occupancy evidence.

## Interpretation and status

All 12 candidate plans and all 21 TF plans are labeled `proposed_not_performed`. MPRA can support allele-sensitive reporter activity; CRISPRi can support endogenous region function; target expression after region repression can support a region-to-gene link. Strong allele-level causal support requires a genotype-verified endogenous allele change with replicated molecular consequences and appropriate controls. No such Section 6 experiment has been conducted in this workflow.

## Reproducible outputs

- `COPD-S6-R001_candidate_shortlist.tsv`: exact candidate evidence and selection fields.
- `COPD-S6-R002_selection_audit.tsv.gz`: all-candidate LD redundancy and selection audit.
- `COPD-S6-R003_MPRA_constructs.tsv`: paired-allele, paired-orientation inserts.
- `COPD-S6-R004_candidate_validation_plan.tsv`: staged endogenous-validation plan.
- `COPD-S6-R004_TF_first_plan.tsv`: candidate-specific TF hypotheses and tests.
- `COPD-S6-R004_cell_context_controls.tsv`: cell models, controls, and limitations.
- `COPD-S6-R005_collaborator_README.md` and package manifest/archive: collaborator handoff.
- `COPD-S6-R005_upstream_inputs.tsv`: upstream checksums.
- `COPD-S6_validation_checks.tsv`: executable consistency checks.
