# COPD V2 external functional benchmark: locked protocol

Protocol version: 1.0. Investigator authorization: attachment
`31c0f2f1-803e-4c5e-81ad-bd39763c2f0a/pasted-text.txt`, 2026-10-05.
Namespace: `COPD-V2-BENCH-*`. The machine-readable protocol-lock manifest
records the UTC lock time, SHA-256, repository baseline and frozen input ledger.
This document must not be edited after locking. Implementation corrections must
be documented separately; no outcome-dependent eligibility or rule changes.

## Scope and temporal separation

Assemble independently sourced COPD functional evidence, freeze the complete
benchmark with labels and exclusions, then join authoritative frozen V1 forward
scores and already completed RC scores. Existing V1, phenotype and RC outputs
are immutable. No retraining, augmentation, new sequences, new RC inference,
matched-control construction, alternative-context modeling, fine-mapping,
target-gene analysis, reranking, manuscript revision or publishing is authorized.
Use only the exact existing V1 allele pairs for model comparison. Assayed
variants without these pairs remain explicit model-unevaluable observations.

The investigator has already disclosed rs2013701's rank and orientation
behavior; prior reports also disclose other selected examples. This is not a
claim of retrospective blinding to those known observations. No new benchmark
aggregate or panel-level model recovery will be extracted until the complete
master benchmark and its constituent source/label tables have been frozen.
Source adjudicators must not inspect V1/RC score or candidate tables during
assembly. The benchmark is not selected from V1 candidates.

## Source search and eligibility

Review the named primary studies: Castaldi FAM13A, Gong multi-locus COPD MPRA,
Zhou HHIP and relevant HHIP functional follow-ups, rs2013701 endogenous editing,
Saferali rs34712979 NPNT splicing, and every source/variant previously cataloged
in V1's functional literature register. Also screen directly cited or retrieved
primary COPD functional studies for additional allele-resolved examples.
Record search queries, retrieval time, source URL/DOI/accession, supplement
availability, exclusions and unresolved access. Search closes at benchmark
freeze; no later addition may depend on model recovery. This is a documented,
targeted evidence audit, not an exhaustive systematic review of all COPD biology.

Eligible studies connect a tested human variant or explicitly defined region
to a COPD/lung-function GWAS locus and report molecular experiments or a
mechanistically resolved allele-dependent transcript outcome. Record COPD
versus related lung-function/shared-disease ascertainment separately. Pure
association, prediction, fine-mapping, motif or QTL significance alone is
supporting evidence, not an experimentally positive sequence-model label.
Clinical cohorts and expression/QTL results may support a tested mechanism.
Nonhuman mechanistic experiments are context, not human allele benchmarks.

Retrieve public article and supplements where possible. Preserve raw downloads
with URLs, hashes, file type and access failures. Do not bypass access controls
or seek controlled participant data. Public text and accessible metadata can
support a case-series extraction but cannot establish an unseen null panel.

## Evidence units and mechanism classes

The master unit is study x exact assayed variant/allele contrast x assay x
biological context. Multiple rows for the same nucleotide are retained and are
not independent variants. Preserve source-table row identifiers. Repeated
technical/biological measurements of one contrast are not separate variants.
Maintain a separate evidence-context table for regions, contacts, unresolved
identities and excluded designs; never promote a region or LD proxy to an
experimentally validated nucleotide.

Use distinct mechanism labels: `MPRA_allele_effect`,
`conventional_reporter`, `endogenous_allele_editing`,
`CRISPRi_region`, `CRISPR_deletion_region`, `chromatin_contact`,
`TF_binding`, `splicing`, `alternative_polyadenylation`,
`expression_QTL`, and `other_or_unresolved`. A study may contribute several
classes, but regional perturbation/contact/QTL does not inherit an allele label.
Haplotype assays are labeled `haplotype_only` unless single-nucleotide contrasts
or an explicitly isolating experiment are reported. Correlated variants are
retained with study/locus and haplotype identifiers; no new LD calculation.

`mechanism_in_model_scope` is `yes`, `partial`, `no`, or `unevaluable`, with a
row-specific rationale. Allele-specific transcriptional reporter/MPRA activity
is in broad regulatory-sequence scope (`yes`) but not necessarily the frozen
bulk-lung cell state. Endogenous editing affecting expression and isolated
TF-binding effects are `partial` unless a matching regulatory mechanism is
established. Regional perturbation/contact/expression-QTL alone is contextual,
not an exact positive/negative assay. Splicing and alternative polyadenylation
are `no` absent an independently demonstrated enhancer/H3K27me3 mechanism.
The H3K27me3 model is not a demonstrated functional silencer assay.

## Experimental states and direction

Labels are assigned from published assay-specific definitions before scores:
`positive`, `null`, `ambiguous`, `conflicting`, `unavailable`, `unevaluable`.
Positive requires a reported allele-dependent experimental effect, or a clearly
typed regional effect in the separate context table. Null requires an explicit
tested contrast with reported nonsignificance/no detected effect under the
study's QC and test rule. Nonsignificance is not proof of biological inactivity.
Untested, failed-QC, missing-result and inaccessible observations are never null.
Reporter activity versus empty vector is distinct from an allele-specific
effect; an active fragment without an allelic contrast is not a positive allele.
Retain source P/FDR/effect values and each study's threshold as reported; do not
derive unreported labels from hit lists or introduce a universal significance
cutoff. When statistical fields or testing rules are ambiguous, label ambiguous.

For multiple cells/assays retain every label. The unique-variant case-series
summary is positive if at least one eligible in-scope assay is positive, null
only if all reported eligible assays are null, and conflicting if explicit
positive/null or opposite-direction results coexist; expose the underlying
counts. Such an any-assay summary is descriptive, not an independent sensitivity
denominator. Do not collapse mechanistically different assays for AUROC.
Record reported tested alleles, higher-activity allele and effect orientation.
Direction comparison requires an unambiguous isolated contrast mapped to
GRCh38 REF/ALT; otherwise not evaluable. Risk-allele and activity directions
are not interchangeable. TF binding, transcript abundance, splicing and contact
do not become reporter-direction labels.

## Exact genomic identity and harmonization

Preserve source rsID, chromosome/coordinate/build, reported allele pair, strand,
sequence and source location. Normalize to GRCh38 forward-strand, 1-based
`chrom:pos:REF:ALT`; retain original representations. Resolve rsIDs with public
authoritative dbSNP/Ensembl records and verify primary-assembly mapping and
reference base with GRCh38 reference sequence or authoritative genomic SPDI.
Reference/alternate in a paper may not mean GRCh38 REF/ALT. Complement and/or
swap only with documented evidence; transform direction consistently.

For rsID-only studies, record source build as not reported if not established,
and distinguish current database mapping from validated source-build liftOver.
Accept an rsID-based exact contrast only if its tested allele pair can be
uniquely reconciled to the authoritative GRCh38 biallelic identity. Missing
tested alleles remain unresolved, even if a coordinate can be looked up.
For palindromic A/T or C/G contrasts, flag strand ambiguity explicitly; genomic
identity may be resolvable but allele-direction evaluation requires reported
forward strand, allele-resolved construct sequence or equivalent evidence.
Reject nonunique/multiallelic mappings unless the assayed alternate is resolved.
Do not replace assayed nucleotides with linked variants.

Coordinate-only observations require an explicit source build. If liftOver is
used, preserve chain/version/hash, require a unique interval mapping, a round
trip to the identical source interval and allele/reference agreement; failed
or ambiguous lifts remain unevaluable. Database accession mappings are not
described as liftOver. For indels require normalized reference-backed exact
representation and tested contrast; no reconstruction of new model sequences.

## Denominator gates and quantitative analysis

For Castaldi and Gong audit separately: designed variants, assayed variants,
QC-passing variants, all per-cell allele results, significant-hit subset, missing
rows and null availability. A complete denominator requires a recoverable
enumerated assay panel, explicit tested allele contrasts, known QC exclusions,
and reported result state for every retained contrast in the specified context.
A design/oligo table plus selected hits does not establish tested nulls unless
the authors explicitly confirm complete testing/QC and exhaustive results.

Study/context sensitivity/specificity-style metrics are allowed only when that
complete assay denominator, positive and null labels, exact identities and
frozen model coverage all pass. Otherwise report selected-positive/case-series
recovery or explicitly `model-covered subset of complete assay panel`; never
claim full-panel performance from an incomplete scored intersection. Show all
unavailable/ambiguous/out-of-scope/missing-model counts in the flow denominator.
False-positive rate against reporter nulls is assay-disagreement, not genomic
FPR or causal error control. No population prevalence or clinical claims.

For valid complete scored study/context panels report fixed-call 2x2 tables,
positive recovery, null call fraction, and distributions of frozen region score
and absolute allele delta by model. If at least five positives and five nulls
exist, compute descriptive AUROC and average precision (AUPRC as stepwise AP)
separately for absolute enhancer delta (primary MPRA discriminator), absolute
H3K27me3 delta and each model's max REF/ALT region score (secondary). No composite
score, threshold optimization or orientation averaging. AP baseline is the
panel's observed positive fraction. Rank enrichment is a descriptive positive
rank distribution within that same panel, not a significance search.

Variants from one locus/haplotype are not independent biological replicates.
Always show locus-level counts and per-locus recovery; do not pool studies as
independent if they repeat variants. For panels spanning at least five distinct
source-defined loci, use 2,000 locus-cluster bootstrap replicates, seed 20261005,
and percentile 95% intervals for eligible metrics, discarding replicates without
both labels for discrimination metrics and reporting valid replicate counts.
With fewer loci, provide descriptive fractions without inferential CI; do not
use a binomial variant-level CI to manufacture precision. Heterogeneous selected
case series receive counts/fractions, not sensitivity, specificity, AUROC or FPR.

## Frozen model and RC comparison after benchmark freeze

Hash the final master table, source inventory, allele harmonization and
assay/context evidence tables before any new aggregate model extraction.
Record freeze UTC and protocol hash. Comparison script must refuse to run
without a matching freeze manifest. No inclusion/label changes after comparison;
a genuine source correction must be appended transparently as a future version,
not silently overwrite this evaluation set.

Join exact canonical allele IDs to authoritative V1 forward scores and original
eligibility/calls/ranks; no rsID-only model join and no proxy. Report enhancer and
H3K27me3 REF, ALT, ALT-minus-REF, region, fixed class-specific thresholds and
calls. Calls retain the original scoreability and blacklist gates. Original
337 candidate rank is present only for actual members; absence is not a new rank.
Recheck calls from frozen scores and original rules without changing results.

For exact variants already in RC output, add frozen RC scores/calls and classify
each model and union as `stable_both`, `forward_only`, `RC_only`,
`negative_both`, or `unevaluable`. Retain distinct model context transitions.
Do not pick a correct orientation using assay labels. Model recovery and
experimental biological evidence remain separate columns and interpretations.
Compare reporter direction only with enhancer ALT-minus-REF where allele
direction is resolved; H3K27me3 changes are reported without calling them proof
of activating/repressing direction. No model-negative splice variant is counted
as an enhancer/H3K27me3 false negative.

Preserve the known rs2013701 `4:88963935:G:T` forward enhancer-positive rank 210
and RC delta-cutoff loss; show the separate MPRA, reporter, contact and endogenous
editing evidence without retuning. Assess whether its orientation category is
common within the eligible exact benchmark, with the sample-size/ascertainment
limitations. Context mismatch is a plausible explanation, not established
causation of a model failure.

## QC, reporting and future-data firewall

Deliver specification, source/search/retrieval inventory, exact harmonization,
frozen master, assay/mechanism context table, forward/RC comparisons,
study/mechanism/locus summaries, valid quantitative metrics or explicit reasons
for omission, detailed report, scripts, manifests/checksums and independent QC.
Validate source row counts, label provenance, identity/alleles, duplicates,
denominators, missingness, exact joins, score/rank fidelity, specification/freeze
hashes, frozen V1/phenotype/RC artifacts and repository boundary. Append V2
activity, decision and result registers only after verified progress.

Report all seven requested interpretive questions, including context and
mechanism limitations and orientation-sensitive recovery. Assay/cell contexts
may be characterized for biological informativeness, but these benchmark
outcomes must not select future model designs. Apparent biological support for
regulatory mechanisms does not establish accuracy or repair RC instability.

Mark the benchmark **frozen external evaluation data**. Its labels must not be
used for training; its outcomes must not determine architecture, hyperparameters,
thresholds, control construction, seed selection, orientation rules or retained
V2 runs. Future model selection must rely exclusively on prespecified internal
validation/calibration criteria. Freeze V2 design before final external use.
The present audit is descriptive V1 characterization, not an untouched V1 test.
At completion stop for investigator review. No next module or commit/push.
