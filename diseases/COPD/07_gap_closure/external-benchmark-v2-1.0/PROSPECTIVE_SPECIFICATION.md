# Prospective frozen external evaluation: external-benchmark-v2-1.0

This is a one-time, evaluation-only application of already selected V2-C models to the unchanged COPD-V2-BENCH1.0 benchmark. This specification must be hashed and frozen before any benchmark outcome or label is interpreted. Raw-byte hashing and checksum-ledger metadata inspection are not outcome opening. Exact source and model hashes are enumerated in `specification/evaluation_specification.json` and `provenance/preopen_source_manifest.json`. The first outcome-ingestion event must verify that freeze, be logged once, and refuse a second open event. Subsequent analysis of that immutable snapshot is the same evaluation, not a new benchmark selection or tuning cycle.

## Immutable inputs and models

Baseline repository commit: `c9515770db4dc38064d4c1fd7fa323788040ec37`. Preserve every earlier stage, source benchmark, original model archive, label, control, calibration and seed unchanged. The prior benchmark has already been examined under V1; it is not newly blinded experimental evidence. In particular, rs2013701 is a pre-existing known case.

Benchmark version: COPD-V2-BENCH1.0. Freeze the exact original master, allele-harmonization, assay/mechanism, study, denominator, V1 forward/RC comparison and related provenance files by their existing checksum ledger. This stage will not query new external sources, rescue unresolved identities, or modify the benchmark. All original benchmark observations, including contextual, unresolved and excluded observations, remain explicitly accounted for.

Use exactly V2-C seeds **104729, 130363, 155921**, in that order, for each of enhancer and H3K27me3-associated contexts. Six original selected `.keras` archives are hash-bound in the JSON specification. No seed may be excluded and no alternative checkpoint substituted. Load with `compile=False`; no optimizer use or model re-saving.

For each allele and seed, obtain real-network forward and nucleotide-RC probabilities. Promote the two float32 network outputs to float64 and calculate `q_seed = (p_forward + p_RC)/2`. The context score is the NumPy float64 mean of the three seed q values in frozen seed order. Compute ALT-minus-REF only after ensembling. Use the original phase-I weights and final 4,560 sigmoid activations, the hash-bound frozen phase-II helper arithmetic, phase-I batch size 32, phase-II batch size 256, `training=False`, no XLA JIT, and the exact deterministic original runtime. Repeat actual-network orientation calls for a separate invariance audit; do not manufacture the audit by swapping stored probabilities. Numerical invariance tolerance is `abs(a-b) <= 1e-6 + 1e-6*abs(b)`; also report exact residuals.

Frozen region thresholds and inclusive operators:

| Context | Decimal17g | Exact float64 encoding | Rule |
|---|---|---|---|
| enhancer | 0.74848511815071117 | `0x1.7f39710000001p-1` | score >= threshold |
| h3k27me3 | 0.76960810025533055 | `0x1.8a0a12aaaaaacp-1` | score >= threshold |

Use `float.fromhex` and verify the decimal encoding. These are region-label thresholds, not allele-effect thresholds or variant false-positive guarantees. No delta threshold will be introduced.

## Exact identity and sequence construction

The immutable identity is GRCh38 chromosome, 1-based position, REF, ALT from the frozen benchmark harmonization. Preserve reported IDs/aliases as provenance; collapse identical exact keys for sequence scoring and variant recovery while retaining every observation. Only explicit frozen exact aliases may map to the same key. No LD proxies, approximate positions, inferred alleles, allele swaps, strand-flip rescue, left-alignment or re-normalization may be introduced. An rsID without an exact sequence-resolved identity remains unresolved.

Use positive hg38 genomic orientation, never gene strand or risk-allele orientation. For `p0=pos1-1`, construct the same previous V1/RC geometry:

```
left  = hg38[chrom][p0-1000:p0]
right = hg38[chrom][p0+len(REF):p0+len(REF)+2001]
REF_sequence = (left + REF + right)[:2001]
ALT_sequence = (left + ALT + right)[:2001]
```

Check the genome's REF bases exactly. Each allele begins at index 1000. Indels share the same left and downstream reference sequence and crop the right end independently, without recentering. Require nonempty ACGT alleles of length at most 1001, sequence length exactly 2001, and ACGTN context. Normalize reference context case only; unsupported context symbols, boundaries, invalid coordinates, absent chromosome, ambiguous alleles or REF mismatch produce explicit unevaluable reasons. N encodes four zeros in A/C/G/T one-hot float32. Nucleotide RC complements ACGTN to TGCAN and reverses positions; an allele of length L occupies RC span `[1001-L,1001)`. Validate allele spans, RC involution, sequence hashes and one-hot orientation equivalence.

Score every distinct exact sequence-eligible benchmark identity regardless of label or mechanism. Sequence eligibility is not biological scope. Both alleles, both orientations, all three seeds and both contexts must be complete for the shared primary recovery denominator. Failed/missing predictions are unavailable, never zeros and never an invitation to average fewer seeds.

## Frozen populations and denominator rules

Keep the existing benchmark taxonomy. Before revealing any V2 predictions, derive and freeze an outcome/identity annotation map from the original benchmark fields and source evidence, documenting any field interpretation. This is a deterministic mapping of frozen evidence, not new experimental curation. The map must separately retain exact identities, positive in-scope observations, exact scorable positives, reporter positives, endogenous-editing positives, splice-only/other out-of-model observations, context-specific null/positive observations, and unresolved identities. NPNT splice-only evidence remains mechanistically out-of-model, not a false negative. Gong null-context observations are not universal negatives; unidentified Castaldi hits are not nulls.

A variant is an in-scope positive when at least one frozen positive observation supports an in-scope regulatory mechanism. Negative/null contexts do not erase its positive context; mixed mechanisms and conflicting contexts remain visible. Reporter and endogenous subsets may overlap and must not be summed as a partition. Primary case-series recovery counts unique exact, in-scope, experimentally positive variants with complete V2 scores; list all exclusions between the original positive population and this denominator. Include an observation-level table alongside deduplicated variant summaries.

No sensitivity, specificity, FPR, AUROC, AUPRC, enrichment p-value, or significance test is planned for this stage. Retain and report the original source denominator audit; do not fabricate completeness. All recovery fractions are descriptive case-series statistics. No internal-reference percentile comparison is planned, and no outcome-conditioned internal reference distribution will be constructed.

## Prespecified outputs and metrics

For each exact scorable identity and context: REF score, ALT score, delta=ALT-REF, absolute delta, unrounded sign (-1/0/+1), maximum allele region score, each allele's region call and either-allele call. Preserve seed/orientation outputs and inference provenance.

Region recovery uses either allele >= its frozen context threshold. Report enhancer, H3K27me3-associated and their OR union on the same primary positive denominator, plus overlapping reporter/endogenous subsets. Report continuous region-score and delta distributions with n, minimum, 25th percentile, median, 75th percentile and maximum, using NumPy float64 linear quantiles. Show individual loci rather than allowing small strata to imply generalization.

Allelic direction requires exact identity, explicit experimentally compared nucleotide alleles, unambiguous REF/ALT mapping and a mechanism supporting the interpreted readout. Experimental sign is ALT-minus-REF activity, never disease-risk direction, a proxy direction or a gene-expression anecdote. For enhancer effects require a direct activating/regulatory readout with a defensible sign. For H3K27me3 effects require explicit H3/repression-related experimental support and a defensible mapping to the model's repressive-state propensity: a greater measured repressive/H3 state maps positive, or a direct transcriptional effect in an explicitly established repressive mechanism maps oppositely. Do not automatically invert general reporter/expression effects to create an H3 denominator. Unsupported H3 direction is unevaluable, potentially leaving n=0.

Report direction observation-by-observation and by locus, with exact concordant/discordant/tied numerators and denominators. Exact model delta=0 is a tie, counted in the strict eligible denominator but not concordant; no epsilon or delta cutoff determines sign. Missing predictions are separate. For a unique-variant consensus summary exclude variants whose eligible experimental directions conflict, rather than majority-voting or selecting a favorable context. Report conflicts and strict versus non-tied fractions distinctly. State explicitly when n is small. Region recovery alone is not evidence of accurate allelic effects.

V1 comparisons use only unchanged frozen forward and RC benchmark region calls/scores and their frozen threshold rules; do not rescore V1. Report own coverage, three-method common-identity intersection, and each pair's shared evaluable-identity intersection. On each shared positive set report retained recovery, newly recovered, lost and still missed. A newly V2-scorable identity missing in V1 is coverage expansion, not a paired gain or V1 false negative. Preserve allele orientation and ID mapping. No independence-assuming significance tests. This comparison is retrospective/descriptive; explicitly flag rs2013701 as previously known.

Stratify, where present, by source study, cell type, assay type, reporter versus endogenous perturbation, enhancer-like versus repressor/H3-relevant mechanism, SNV versus indel, and exact V1 presence versus newly V2-scorable sequence. Missing strata remain explicit and small strata descriptive. Do not infer missing cell types or mechanism direction.

## Hard post-opening firewall and verification

After the single outcome-opening event, prohibit training, architecture/model changes, seed choice, checkpoint replacement, altered controls/labels, new calibration or threshold changes, and return to internal training/test evaluation. Read-only hash-bound reuse of frozen model helpers, thresholds and archives is permitted; execution of earlier stage entrypoints is not. No broader COPD GWAS/candidate universe, 337-candidate reranking, fine mapping or gene targeting is authorized.

Implementation errors must be logged and repaired only to implement this unchanged specification; preserve failed attempts and do not use outcomes to choose methods. Any invalid model/source contract is a stop condition. Independently validate input preservation, prospective chronology, identity and sequence geometry, complete six-model inference, ensemble/delta/threshold arithmetic, actual-network RC invariance, population accounting, direction mapping, paired comparisons, and exact denominators. Preserve all earlier artifacts and append only the three V2 registers. Freeze all final new-stage artifacts by SHA-256 with explicit self-reference exclusions.

Required final report: exact scorable positive count; enhancer/H3/union recovery; frozen V1 comparison; strict direction concordance; whether improvement is supported and whether it concerns region recovery only or also REF-ALT effects; unresolved failures and limitations. Report poor performance honestly. Stop after this stage. No commit/push is authorized by this evaluation request.
