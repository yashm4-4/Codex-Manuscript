# COPD V2 external functional benchmark and frozen V1 comparison

Date: 2026-10-05. Module: `COPD-V2-BENCH`. Status: completed external-evidence
audit; final independent checks and artifact hashes are recorded in the linked
[validation ledger](../provenance/COPD-V2-BENCH_final_validation.tsv) and
[manifest](../provenance/COPD-V2-BENCH_final_validation_manifest.json).

## Main finding

Among **40 distinct variants with at least one reported in-scope positive
assay**, **23 have an exact existing frozen V1 sequence pair**. V1's original
forward rule recovers **1/23 (4.35%)**, rs2013701; the already completed RC audit
recovers **0/23**. The other 22 are negative in both orientations. This is
**descriptive, coverage-conditioned case-series recovery, not sensitivity**.
Seventeen positive cases cannot enter this comparison: 16 lack the exact
existing V1 record and one has unresolved deposited-construct identity.

The complete source benchmark contains 14,025 rows, 1,731 reported rsIDs,
13,899 exact-identity rows and 1,710 distinct exact GRCh38 variants. Existing
V1 scores cover 726 distinct exact variants. Measurements, selected positive
summaries, repeated constructs and distinct cell contexts are explicitly typed;
14,025 rows are not 14,025 independent variants or biological replicates.

The data support regulatory mechanisms at COPD-associated loci. They do not
establish general accuracy of the frozen binary nomination rule, repair its RC
instability, or justify changing a threshold. No retraining, new sequences,
new inference, orientation averaging, cell-context model, matched-control model,
fine-mapping, target-gene analysis, reranking or manuscript revision occurred.

## 1. Prespecification, independence and freeze

The [protocol](../provenance/COPD-V2-BENCH_analysis_specification.md) was locked
at **23:02:10 UTC**, before source adjudication or new benchmark-wide model
comparison. Its SHA-256 is
`46d3aac26b1f980d7c2dcd5034d6cc64566524ca7d99cc635f7d979e04f2de65`.
Source adjudicators did not inspect model outputs. Previously disclosed selected
examples, particularly rs2013701, were already known to the investigator and
lead analyst; this is not a claim of retrospective blinding to those examples.

After source and identity reviews, **200 source/protocol/assembly artifacts**
were frozen at **23:31:15 UTC**, before the new comparison started. The
[freeze manifest](../provenance/COPD-V2-BENCH_benchmark_freeze.json) pins every
constituent file. The [master table](COPD-V2-BENCH-R003_frozen_benchmark_master.tsv)
SHA-256 is
`a69a8b1d5c82299075b2d5203fee3d4a23105572c6cf456bab834b82ef13bbdd`.
No variants or labels were added, removed or changed after viewing performance.

The search was a documented targeted audit of named COPD functional studies,
the eight V1 literature examples, and directly encountered primary follow-ups;
it was not an exhaustive systematic review. Source searches, inaccessible
resources and exclusion reasons are retained. "External" denotes independent
experimental evidence, not a verified sequence-disjoint V1 training holdout.
No new phase-I training-corpus overlap analysis was performed.

## 2. Source and denominator audit

### Castaldi FAM13A: complete count/design data, incomplete published labels

The investigator identified [GEO GSE109452](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE109452).
Both public files were retrieved through standard NCBI HTTPS/FTP paths:
`GSE109452_FAM13A_MPRA.txt.gz` and
`GSE109452_secondary_design_barcodes.dat.gz`; raw bytes, URLs and checksums are
preserved under `data/COPD-V2-BENCH/castaldi/`.

Independent count reconciliation establishes 606 SNPs, 7,270 allele constructs,
3,635 paired SNP/window/orientation contexts, 239,910 designed barcodes and
218,526 observed barcode rows. Every designed construct is represented by at
least one observed barcode. The 21,384 absent designed tags are not null variants.
The file contains input/output counts for two promoter experiments, not
per-contrast P values, FDR decisions or the authors' complete QC/exclusion ledger.
It therefore recovers the tested identities and measurements but does **not**
establish a complete labeled positive/null denominator. No significance analysis
was invented from counts.

The [primary article](https://pmc.ncbi.nlm.nih.gov/articles/PMC6353020/) reports
45 MPRA-positive SNPs and 85 significant oligonucleotide contexts. Its selected
Table 1 exposes eight SNPs and 32 promoter-specific statistics. The numerical
decision rule could not be verified from the retrieved primary methods, so those
individual contexts remain ambiguous, not uniformly positive. Three further
hit identities are explicitly reported as reuse of the same experiment by
[Lin et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC7063182/); they are not
independent replications. Thus 11 of 45 reported hit identities are available.
Thirty-four remain unidentified at the hit-list level, not presumed null.

Castaldi contributes 7,289 rows: 7,270 individual MPRA contexts, 12 explicitly
typed variant-level summaries, six follow-up reporter contrasts and one
endogenous expression-editing contrast. The 12 summaries comprise 11 reported
positive identities and the explicitly negative rs4416442 observation; they are
not 12 additional independent experiments. The six reporters have three
reported positive and three explicit null contrasts. Reporter-only rs147089648
adds a 607th study-associated identity beyond the 606-SNP MPRA design.

The supplementary PDF remained inaccessible after documented official-endpoint
attempts; the successful GEO retrieval supersedes any suggestion that count data
were unavailable. Counts-only measurements remain `unavailable` **for a
published binary label**, not unavailable as raw measurements. MPRA effect signs
and effect alleles are retained, but their input/output-ratio convention was
not silently converted to higher-activity alleles.

### Gong multi-locus MPRA: explicit nulls, incomplete QC denominator

The [public supplement](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/ajrcmb/PAP/10.1093_ajrcmb_aanag179/1/aanag179_supplementary_data.zip)
contains design oligos, three-cell MPRA results, methods, figures and selected
follow-up tables. The [published abstract](https://pubmed.ncbi.nlm.nih.gov/42745457/)
provides study-level context; the subscriber-only main reporter results were
not invented from methods.

The deposited design has 1,120 experimental rsIDs plus 34 controls. The methods
require all three BH-adjusted P values—QuASAR, mpralm_sum and DESeq2—to be below
0.05. Applying that published conjunction yields 29 positive cell/orientation
observations representing 25 variants. A complete reported contrast failing
that rule is a study-specific null, not proof of biological inactivity.

| Cell | Designed cell/orientation contrasts | Positive | Explicit consensus-null | Unavailable |
|---|---:|---:|---:|---:|
| 16HBE | 2,240 | 16 | 2,222 | 2 |
| HUVEC | 2,240 | 6 | 2,220 | 14 |
| MRC5 | 2,240 | 7 | 2,229 | 4 |
| Total | 6,720 | 29 | 6,671 | 20 |

All six cell/orientation results for rs76419734 are absent, and its duplicated
design has conflicting allele-role/flank representations. Fourteen additional
reported observations lack a required method statistic. The missing-result/QC
reason cannot be reconstructed completely. Hence the locked full-denominator
gate fails despite extensive recoverable null data. The 34 controls remain in
the separate context/exclusion table. No GWAS effect allele is treated as a
measured reporter-effect direction.

### Other mechanism-stratified evidence

| Resource | Exact assay rows | Treatment |
|---|---:|---|
| Zhou HHIP | 5 | Four isolated reporter tests: two positive, two explicitly NS; one TF-binding assay is partial model scope |
| Boueiz ACVR1B | 2 | Positive reporters in 16HBE and Jurkat; no existing exact V1 sequence pair |
| Hao DSP | 2 | Reporter and endogenous HDR expression evidence kept separate |
| Saferali NPNT | 1 | Positive splice-mechanism evidence; outside enhancer/H3K27me3 scope |
| Lin FAM13A metabolic follow-up | 6 | One explicit HepG2 reporter positive, five ambiguous unmarked contrasts; not independent COPD disease validation |

[HHIP](https://pmc.ncbi.nlm.nih.gov/articles/PMC3284120/),
[ACVR1B](https://pmc.ncbi.nlm.nih.gov/articles/PMC6444627/),
[DSP](https://pmc.ncbi.nlm.nih.gov/articles/PMC7605184/) and
[NPNT](https://pmc.ncbi.nlm.nih.gov/articles/PMC11968218/) primary sources support
these distinctions. HHIP double-mutant haplotypes, chromatin contacts,
[Stuart regional perturbations](https://pmc.ncbi.nlm.nih.gov/articles/PMC7252574/),
[TGFB2 regional deletion](https://pmc.ncbi.nlm.nih.gov/articles/PMC6693893/),
HHIP follow-ups, QTL evidence and predicted TF effects remain contextual.
They are not promoted to isolated nucleotide validation. No exact alternative-
polyadenylation perturbation satisfying the criteria was recovered; that class
remains separate rather than being relabeled as enhancer evidence.

## 3. Exact alleles, source conflicts and exclusions

The [harmonization table](COPD-V2-BENCH-R002_exact_allele_harmonization.tsv)
preserves source alleles, coordinates, build statements, strand uncertainty,
palindromic flags, current primary-assembly mappings and reference checks.
Ensembl/dbSNP mappings and the existing GRCh38 FASTA were used. Database assembly
mapping is not described as liftOver; **no liftOver was used**, so round trips
are explicitly not applicable. Castaldi coordinates match GRCh37 database
positions, but the retrieved source did not explicitly state its build; that
inference is distinguished from reported metadata.

Deposited Gong allele inserts were checked against the reference and the
authoritative allele contrast, with reference-backed left normalization for
indels. Exact identities were resolved for 1,099 of 1,120 reported Gong IDs.
Twenty-one IDs, comprising 126 assay rows, remain identity-unevaluable. Published
experimental labels are retained even when identity fails.

The source audit identified 14 internally discordant reverse ALT insertion
constructs, duplicated rs76419734 design ambiguity, and other deposited-versus-
database contrast discrepancies. These are not silently repaired. For example,
rs57658727's deposited insertion does not match its current database contrast;
its positive experiment cannot validate a V1 nucleotide by rsID alone. A
single-SNP fallback resolves rs11732650 using an isolated C/G contrast and a
unique exact terminal anchor; its shared nonreference background is documented.
No approximate indel rescue was allowed. Explicit dbSNP merges are recorded,
not interpreted as LD proxies.

Published Gong Figure S1 locus ranges define source-region grouping, not new
target-gene assignments. Reconstructed counts are GSTCD 493, EEFSEC 224,
ADAM19 203, ADGRG6 109 and HTR4 91. Figure S1 states 492, 224, 204, 108 and 92;
both total 1,120. Both sets are retained without forcing agreement.

Two other source conflicts are explicit: the frozen V1 literature annotation
overstates rs1795739 promoter contact, whereas the primary Castaldi result is
negative for that region; and Gong's regional supplementary results contain
numeric/significance and CRISPRi-versus-active-Cas9 terminology inconsistencies.
V1 was not edited. These limitations cannot change nucleotide labels or become
model false negatives.

## 4. Model coverage and recovery

Only exact canonical alleles were joined to existing V1 sequences/scores and
the completed RC outputs. Missing sequences were not reconstructed and missing
scores were not replaced with zero. Original blacklist/eligibility gates,
class-specific cutoffs, 337-candidate membership and ranks remain unchanged.

| Source | Distinct reported assay variants | Exact identities | Existing scored variants | In-scope positive cases scored | Forward union recovery | RC union recovery |
|---|---:|---:|---:|---:|---:|---:|
| Castaldi | 607 | 607 | 75 | 6 | 1/6 | 0/6 |
| Gong | 1,120 | 1,099 | 648 | 14 | 0/14 | 0/14 |
| HHIP | 2 | 2 | 2 | 2 | 0/2 | 0/2 |
| DSP | 1 | 1 | 1 | 1 | 0/1 | 0/1 |
| ACVR1B | 1 | 1 | 0 | 0 | Unevaluable | Unevaluable |
| Lin metabolic follow-up | 3 | 3 | 0 | 0 | Unevaluable | Unevaluable |
| NPNT splice study | 1 | 1 | 1 | Not in scope | Not a false negative | Not a false negative |

Study counts overlap: the Lin variants already occur in Castaldi, and NPNT
also occurs in Gong. Deduplication yields 726 scored exact variants overall,
not the sum of study counts. The 23 in-scope positive cases span eight source
loci; 17 have positive/null coexistence across contexts. The summary state
`conflicting` denotes this heterogeneity, not necessarily irreproducibility.

Gong's 14 scored positive variants comprise EEFSEC 5, GSTCD 5, ADAM19 2,
ADGRG6 1 and HTR4 1; none meets the frozen union rule. Castaldi's six scored
positive identities are rs2464523, rs7674369, rs1964516, rs7671167, rs2013701
and rs1795739. Only rs2013701 passes. These linked observations do not provide
23 independent biological replications.

The [full forward table](COPD-V2-BENCH-R007_V1_forward_comparison.tsv) and
[orientation table](COPD-V2-BENCH-R008_forward_RC_comparison.tsv) retain every
row, including nulls and unevaluable states. Among Gong variants with at least
one reported null context, 648 are scored; seven have a forward union call and
eight an RC union call. Those are descriptive overlaps, **not false-positive
rates**: contexts repeat, some variants are positive in another context, and
the complete QC denominator is not established.

### Fixed-gate failure patterns

The enhancer region threshold is 0.643623; the H3K27me3-associated threshold is
0.58505. The SNV absolute-delta cutoffs are 0.0570631877 and 0.0288026139;
the indel/complex cutoffs are 0.0493609385 and 0.0261833595, respectively.
The original rule requires both region and delta gates plus original eligibility.

| Forward model, 23 positive cases | Recovered | Region gate only fails | Delta gate only fails | Both fail | Eligibility excluded |
|---|---:|---:|---:|---:|---:|
| Enhancer | 1 | 2 | 1 | 19 | 0 |
| H3K27me3-associated | 0 | 2 | 1 | 20 | 0 |

The two forward region-only failures in both models are rs72671892 and DSP
rs2076295. DSP has enhancer delta −0.199954 but maximum score 0.568677, below
0.643623. Its RC delta is −0.223099, with maximum score 0.460304; it remains
negative. Conversely, rs141807665 has high region scores but insufficient allele
deltas. These examples distinguish region classification from allelic sensitivity;
they are not an instruction to relax either gate.

For six scored, positive reporter contrasts with resolved activity direction,
the enhancer delta has the reported sign in all six: three FAM13A variants,
two HHIP variants and DSP. They span only three loci, are selected examples,
and include very small deltas. Sign agreement is not candidate recovery or
effect-size calibration. Gong has no deposited signed effect estimates, and
Castaldi MPRA ratio direction is unresolved, so neither is added to this direction
denominator. H3K27me3 changes are not called validated repression/activation.

## 5. Orientation sensitivity and rs2013701

| Union-call category | All 726 scored benchmark variants | 23 scored in-scope positive cases |
|---|---:|---:|
| Stable positive in both | 3 | 0 |
| Forward-only | 7 | 1 |
| RC-only | 6 | 0 |
| Negative in both | 710 | 22 |

Overall benchmark union calls change from 10 forward to nine RC, while the
single known-positive recovery disappears. Global negative agreement does not
establish functional-variant recovery. Experimental insert orientation is a
different variable from computational reversal of the exact frozen 2,001-bp
sequence; these were never equated or averaged.

### rs2013701: biological evidence persists when the model call changes

Exact identity: **4:88963935:G:T**. Original V1 enhancer-positive candidate
rank: **210**. The publication supplies MPRA evidence, a 16HBE reporter result,
regional promoter-contact evidence and endogenous allele-editing effects on
FAM13A expression. The T allele's higher expression in the reporter/editing
experiments is biological evidence independent of this computational score.
The contact assay localizes a region, not an isolated base.

| Frozen score | Forward | RC |
|---|---:|---:|
| Enhancer REF | 0.846710801 | 0.884276330 |
| Enhancer ALT | 0.910054922 | 0.937941968 |
| Enhancer ALT−REF | 0.063344121 | 0.053665638 |
| Enhancer maximum | 0.910054922 | 0.937941968 |
| Enhancer call | Positive | Negative |
| H3K27me3 REF | 0.036009070 | 0.047603641 |
| H3K27me3 ALT | 0.063762069 | 0.085098609 |
| H3K27me3 ALT−REF | 0.027752999 | 0.037494969 |
| H3K27me3 call | Negative | Negative |

The enhancer region score increases and the gain sign persists, but the RC
delta falls below the unchanged 0.0570631877 cutoff. The H3K27me3 RC delta
passes its cutoff but its region score remains below 0.58505. Neither rule is
modified. This is an orientation-sensitive model recovery, not contradictory
biological evidence or proof that the experiment is wrong.

rs2013701 is **exceptional within this evaluable case series**: it is the only
forward-positive functional case. It cannot represent broad recovery, and its
1/1 loss among recovered cases is too small and ascertained a denominator for a
general orientation-failure probability. The other 22 cases are not rescued
by switching orientation.

## 6. Mechanism and context interpretation

NPNT rs34712979 maps exactly to **4:105897896:G:A** and is negative under both
models in both orientations. The Saferali positive evidence concerns altered
splicing and transcript isoforms, not a matched enhancer/H3K27me3 experiment.
It is explicitly **outside model scope, not a sequence-model false negative**.
The same nucleotide has six in-scope reporter-null observations in Gong; these
are distinct assays and do not negate its splice mechanism. QTL associations,
contacts, regional deletions and haplotypes are likewise not interchangeable
with isolated allele-dependent reporter or endogenous-editing evidence.

The seven requested interpretive answers are:

1. **How often does V1 recover experimentally active regulatory variants?**
   One of 23 exact scored cases with any reported in-scope positive assay
   (4.35%), conditional on coverage and selection; not population sensitivity.
   Seventeen additional positive identities are not model-evaluable.
2. **Which failures might reflect context?** Low bulk-lung region scores at
   airway reporter loci, and the DSP region-gate failure despite substantial
   allele deltas, are compatible with context mismatch. Short episomal constructs,
   cell lines, promoter/window choice and shared nonreference backgrounds also
   differ from frozen genomic sequences. This audit does not establish which
   factor caused any individual failure; genuine model limitations remain possible.
3. **Which mechanisms are outside scope?** NPNT splicing is the clear positive
   example. Alternative polyadenylation would also require separate treatment;
   no eligible exact perturbation of that class was recovered. Regional/contact/
   QTL evidence is not an exact binary validation target in the first place.
4. **How much does orientation change recovery?** From 1/23 forward to 0/23 RC;
   no positive case is RC-only or stably recovered in both. This does not imply
   that every negative prediction is caused by orientation instability.
5. **Is rs2013701 representative?** No: it is the sole recovered forward case
   here, with unusually rich orthogonal evidence. Its score instability is
   illustrative, but cannot estimate a general error probability by itself.
6. **Does the evidence support a regulatory-model strategy?** The functional
   studies justify investigating regulatory mechanisms. Selected sign agreement
   suggests that some directional information is present. Sparse binary recovery
   and the loss of the sole recovered case do not validate this frozen predictor
   or justify confidence in its current candidate calls.
7. **Which assay/cell contexts are informative?** Allele-specific airway reporters
   paired with endogenous nucleotide editing provide direct mechanistic evidence;
   multi-cell MPRA exposes context dependence. Primary-airway regional perturbation
   adds context but not nucleotide specificity; HepG2/Jurkat/HUVEC results need
   their tissue caveats. These are descriptions of evidence quality, **not choices
   of future V2 architecture, cell model or retained run based on benchmark outcomes**.

The internal `silencer` field name is reported as **H3K27me3-associated**.
Active reporters are not a direct functional-silencing benchmark; failure of
that model to call an activating reporter is not automatically a mechanistically
matched false negative.

## 7. Quantitative limits and QC

No complete, exact-scored, prespecified positive/null screening panel satisfies
all locked gates. Consequently no sensitivity, specificity, AUROC, average
precision, false-positive rate, rank-enrichment test or inferential confidence
interval is reported. [Metric eligibility](COPD-V2-BENCH-R012_quantitative_metric_eligibility.tsv)
gives explicit reasons per study/assay/cell/orientation. Correlated loci and
repeated assays are not treated as independent binomial replicates. The
prespecified cluster bootstrap is not run when its underlying metric is invalid.
Score distributions and original candidate ranks remain available descriptively.

Independent source review reproduced all Gong statistics/labels and checked
HHIP/Lin figures, NPNT scope and the eight historical literature examples.
Independent identity review reproduced all Gong mappings and exclusions.
The lead analyst separately reconciled the original Castaldi count/design files.
The comparison pipeline reports **5,290 PASS checks**; a separate validator
reconstructs scores, calls, ranks, summaries, missingness and immutable hashes
without importing the comparison implementation. Final validation also verifies
the 316-entry frozen-input ledger and repository boundary.

Implementation-only refinements after comparison added explicit gate-failure,
orientation and null-overlap summaries and replaced an exemplar-assay direction
field in the unique-variant table with an assay-count summary. A validator-only
NPNT check was narrowed to the splice assay rather than every assay at that
nucleotide. Run history preserves these corrections. No benchmark labels,
membership, harmonization rules, frozen scores or thresholds changed.

## 8. Deliverables and reproducibility

| Artifact | Contents |
|---|---|
| R001 | Source/study inventory |
| R002 | Complete exact-allele harmonization, failures and strand/build flags |
| R003 | Frozen master benchmark, including unavailable/ambiguous/null states |
| R004 / R004B | Assay/mechanism classification / contextual and excluded evidence |
| R005 / R006 / R006B | Denominator audit / search audit / source-locus grouping |
| R007 / R008 | Every-row V1 forward / forward-versus-RC comparison |
| R009 / R010 / R011 | Study-context, mechanism and locus summaries |
| R012 / R013 / R014 | Metric gates, score distributions and coverage flow |
| R015 / R016 | Deduplicated variant summaries / analysis validation |
| R017 / R018 | Fixed-gate failures / orientation recovery counts |

All result names begin `COPD-V2-BENCH-`; V1 result IDs are unchanged. Raw primary
downloads, failed-response records, extraction ledgers and curated JSON are in
`../data/COPD-V2-BENCH/`. Collection scripts reproduce source extraction from
cached files. `initialize_external_benchmark.py` locks the protocol;
`assemble_external_benchmark.py` performs source-only assembly and refuses to
overwrite a frozen benchmark; `compare_external_benchmark.py` requires every
freeze hash to match before reading model values. `validate_external_benchmark.py`
performs independent final QC. Python/package versions and commands are recorded
in manifests; no GPU or model execution was needed.

From the repository root, post-freeze comparison can be reproduced with
`python3 -B diseases/COPD/07_gap_closure/scripts/compare_external_benchmark.py`.
This rewrites only derived comparison outputs and appends run history; it does
not amend the frozen source set. Final QC is then run with
`python3 -B diseases/COPD/07_gap_closure/scripts/validate_external_benchmark.py --comparison-complete`.
The [checksum ledger](../provenance/COPD-V2-BENCH_artifact_checksums.tsv) pins
the completed artifact bundle. Shared V2 documentation/registers evolve; earlier
PHENO/RC checksum ledgers retain their historical shared-file states unchanged.

## 9. Frozen external-evaluation firewall and stop

This benchmark is **FROZEN EXTERNAL EVALUATION DATA**. Its labels must not be
used for training. Its outcomes must not determine architecture, hyperparameters,
thresholds, controls, seeds, orientation rules, cell-model selection or which
V2 run is retained. Future model selection must use only prespecified internal
validation/calibration criteria. Freeze the future V2 model design before using
this set for its final external comparison; keep all frozen variants and labels,
including failures and unavailable records.

The current result is descriptive V1 characterization, not an untouched V1
model-selection test. This completed module stops for investigator review.
No subsequent V2 module, manuscript revision, commit or push was performed.
