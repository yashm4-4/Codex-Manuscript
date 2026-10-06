#!/usr/bin/env python3
"""Render the completed frozen-method external evaluation; never rerun science.

Reads only this new stage's results, prediction provenance and annotations.
Writes only EXTERNAL_BENCHMARK_V2_REPORT.md and LIMITATIONS.md. No inference,
metric generation, benchmark reopening, prior-stage entrypoint or model mutation.
"""
import argparse
from collections import defaultdict
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path

STAGE = Path(__file__).resolve().parents[1]
CONTEXTS = ('enhancer', 'h3k27me3', 'union')


def load(relative):
    return json.loads((STAGE / relative).read_text())


def rows(relative):
    with (STAGE / relative).open(newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def true(value):
    return value is True or str(value).lower() in {'true', '1'}


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def scalar(value):
    if value is None or value == '':
        return 'NA'
    if isinstance(value, bool):
        return 'yes' if value else 'no'
    return str(value).replace('|', '\\|').replace('\n', '<br>')


def num(value, precision=7):
    return 'NA' if value is None or value == '' else format(float(value), f'.{precision}g')


def call(value):
    return 'NA' if value is None or value == '' else 'yes' if true(value) else 'no'


def table(headers, data):
    if not data:
        return '_No eligible rows; denominator is zero._\n'
    return '\n'.join(['| ' + ' | '.join(headers) + ' |',
                      '| ' + ' | '.join('---' for _ in headers) + ' |'] +
                     ['| ' + ' | '.join(scalar(value) for value in row) + ' |' for row in data]) + '\n'


def fraction(numerator, denominator):
    return f'{numerator}/{denominator}' if int(denominator) else '0/0 (unevaluable)'


def recovery_table(data, denominator_name):
    groups = defaultdict(dict)
    for row in data:
        groups[(row.get('subset', row.get('population', '')), row['method'])][row['context']] = row
    output = []
    for (subset, method), group in groups.items():
        output.append([subset, method] + [fraction(group[c]['recovered'], group[c][denominator_name]) for c in CONTEXTS])
    return table(['Population', 'Method', 'Enhancer', 'H3K27me3-associated', 'Union'], output)


def relation(new, old):
    return 'higher' if new > old else 'lower' if new < old else 'unchanged'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--refresh', action='store_true',
                        help='Refresh only the two generated administrative Markdown reports before final validation; never modify scientific results')
    args = parser.parse_args()
    targets = [STAGE / 'EXTERNAL_BENCHMARK_V2_REPORT.md', STAGE / 'LIMITATIONS.md']
    if not args.refresh and any(path.exists() for path in targets):
        raise RuntimeError('Reports already exist; --refresh is explicit and report-only')
    if (STAGE / 'provenance/freeze.json').exists():
        raise RuntimeError('Final stage already frozen; do not rewrite reports')
    if (STAGE / 'provenance/independent_final_validation.json').exists():
        raise RuntimeError('Final validation already binds the report; do not rewrite it')
    prefinal_path = STAGE / 'provenance/independent_validation_pre-final.json'
    prefinal = json.loads(prefinal_path.read_text())
    assert prefinal['status'] == 'PASS', 'Independent pre-final validation must pass before report rendering'
    manifest = load('results/evaluation_manifest.json')
    inference = load('predictions/inference_manifest.json')
    invariance = load('predictions/real_network_invariance.json')
    prospective = load('provenance/prospective_freeze.json')
    specification = load('specification/evaluation_specification.json')
    opening = load('provenance/benchmark_opening_event.json')
    ingestion = load('provenance/benchmark_ingestion_complete.json')
    annotation = load('provenance/annotation_freeze.json')
    annotation_rules = load('provenance/annotation_rules.json')
    source_manifest = load('provenance/preopen_source_manifest.json')
    input_manifest = load('inputs/input_manifest.json')
    phase_one = load('predictions/phase_I_provenance.json')
    start = load('predictions/started.json')
    environment = load('predictions/environment.json')
    for value in (manifest, inference, invariance, prospective, annotation, source_manifest, input_manifest):
        assert value['status'] == 'PASS'
    assert inference['completed'] is True
    # Verify already-produced result bytes instead of regenerating any metric.
    for item in manifest['files']:
        path = STAGE / item['path']
        assert path.stat().st_size == item['bytes'] and sha(path) == item['sha256'], item['path']
    flow = {row['population']: int(row['count']) for row in manifest['population_accounting']}
    recovery = manifest['region_recovery']
    common = manifest['common_identity_recovery']
    paired = manifest['paired_comparisons']
    direction_summaries = manifest['direction_summaries']
    positive_rows = rows('results/region_recovery_per_variant.tsv')
    all_variants = rows('results/all_variant_evaluation.tsv')
    direction_rows = rows('results/direction_resolved_allelic_effects.tsv')
    consensus_rows = rows('results/direction_unique_variant_consensus.tsv')
    missing = rows('results/unevaluable_variants.tsv')
    expansion = rows('results/newly_V2_scorable_positive_variants.tsv')
    stratified = rows('results/source_mechanism_stratified_recovery.tsv')
    distributions = rows('results/continuous_score_distributions.tsv')
    denominator_rules = rows('results/frozen_denominator_rules.tsv')
    direction_loci = rows('results/direction_locus_summary.tsv')
    strict = next(row for row in direction_summaries if row['context'] == 'enhancer' and row['subset'] == 'all_strict_observations')
    unique_direction = next(row for row in direction_summaries if row['context'] == 'enhancer' and row['unit'] == 'unique_variant')
    primary = {row['context']: row for row in recovery if row['method'] == 'V2_C_symmetric' and row['subset'] == 'all_in_scope_positive'}
    common_union = {row['method']: row for row in common if row['context'] == 'union'}
    cnew = int(common_union['V2_C_symmetric']['recovered'])
    union_conclusion = '; '.join(
        f"{relation(cnew, int(common_union[method]['recovered']))} than {method} "
        f"({fraction(cnew, common_union['V2_C_symmetric']['denominator'])} versus "
        f"{fraction(common_union[method]['recovered'], common_union[method]['denominator'])})"
        for method in ['V1_forward', 'V1_RC'])
    if all(cnew > int(common_union[method]['recovered']) for method in ['V1_forward', 'V1_RC']):
        paired_union_gains = {cnew - int(common_union[method]['recovered']) for method in ['V1_forward', 'V1_RC']}
        improvement = ('The corrected V2-C ensemble shows a modest one-variant descriptive regional gain over each frozen V1 orientation on the shared positives, not broad external functional recovery.'
                       if paired_union_gains == {1} else
                       'The corrected V2-C ensemble has higher descriptive region-level union recovery on shared positives; the absolute recovery and exact paired gain sizes below limit the interpretation.')
    elif all(cnew == int(common_union[method]['recovered']) for method in ['V1_forward', 'V1_RC']):
        improvement = 'The corrected V2-C ensemble does not improve descriptive region-level union recovery on the shared positives; the union count is unchanged from both frozen V1 orientations.'
    elif all(cnew < int(common_union[method]['recovered']) for method in ['V1_forward', 'V1_RC']):
        improvement = 'The corrected V2-C ensemble has lower descriptive region-level union recovery than both frozen V1 orientations on the same evaluable positives.'
    else:
        improvement = 'External region-level recovery is mixed relative to the two frozen V1 orientations; a blanket improvement claim is not supported.'
    not_recovered = int(primary['union']['evaluable_positive_denominator']) - int(primary['union']['recovered'])
    absolute_recovery = (
        f"Absolute recovery remains weak: {fraction(not_recovered, primary['union']['evaluable_positive_denominator'])} "
        'exact in-scope positive variants are not recovered by either frozen region threshold. '
        'These are region-threshold non-recoveries in a context-mismatched case series, not asserted biological false negatives.')
    by_exact_id = {row['variant_id']: row for row in positive_rows if row['variant_id']}
    eefsec = by_exact_id['3:128111036:C:A']
    gstcd = by_exact_id['4:105894906:T:TA']
    acvr = by_exact_id['12:51954475:A:G']
    assert true(eefsec['V2_union_region_call']) and not true(eefsec['V1_forward_union_region_call']) and not true(eefsec['V1_RC_union_region_call'])
    assert true(gstcd['V1_forward_enhancer_region_call']) and not true(gstcd['V2_enhancer_either_allele_region_call'])
    assert true(gstcd['V2_h3k27me3_either_allele_region_call']) and true(gstcd['V2_union_region_call'])
    assert acvr['V1_presence'] == 'absent_unavailable' and true(acvr['V2_enhancer_either_allele_region_call'])
    acvr_direction = [row for row in direction_rows if row['rsid'] == 'rs7962469' and row['context'] == 'enhancer']
    assert len(acvr_direction) == 2 and all(row['concordance'] == 'discordant' for row in acvr_direction)
    case_interpretation = (
        'The sole newly recovered paired union case is rs35421223 at EEFSEC. '
        'rs141807665 at GSTCD loses the V1-forward enhancer region call, but its H3K27me3-associated region call retains union recovery; this is not a lost union case. '
        'rs7962469 at ACVR1B is recovered among newly V2-scorable sequences and is a coverage expansion, not a paired V1 gain. '
        'Despite that regional call, its enhancer ALT-minus-REF direction is discordant with both experimental reporter cell contexts, 16HBE and Jurkat. '
        'Those are two assay observations of one variant, not two model-context errors; H3K27me3 direction remains unevaluable.')
    generated = datetime.now(timezone.utc).isoformat()
    report = ['# Frozen V2-C external functional benchmark evaluation',
              '', f'Stage: `external-benchmark-v2-1.0`. Report generated: {generated}.',
              '', '## Principal findings', '',
              f"The unchanged V2-C models scored {flow['exact_V2_scorable_identity_keys_all_states_mechanisms']:,} exact benchmark identities, including "
              f"{flow['exact_V2_scorable_in_scope_positive_primary']} exact in-scope positive variants. "
              f"Enhancer region recovery was {fraction(primary['enhancer']['recovered'], primary['enhancer']['evaluable_positive_denominator'])}; "
              f"H3K27me3-associated recovery was {fraction(primary['h3k27me3']['recovered'], primary['h3k27me3']['evaluable_positive_denominator'])}; "
              f"their union was {fraction(primary['union']['recovered'], primary['union']['evaluable_positive_denominator'])}. "
              'These are descriptive case-series recovery fractions, not sensitivity.', '',
              improvement, f'On the common-identity union denominator, V2 is {union_conclusion}.',
              absolute_recovery, '',
              f"Strict enhancer-direction concordance was {fraction(strict['concordant'], strict['denominator'])} observations "
              f"({strict['discordant']} discordant, {strict['tied']} exactly tied), and "
              f"{fraction(unique_direction['concordant'], unique_direction['denominator'])} unique nonconflicting variants. "
              'This is a small, correlated direction-resolved subset. No H3K27me3 direction denominator is defensible from the frozen assays. '
              'Any region-recovery improvement does not establish improved REF–ALT effect prediction; no formal V1-versus-V2 direction-improvement analysis was prespecified.', '',
              '## Prospective design, one-time opening and frozen model contract', '',
              table(['Event', 'UTC / evidence'], [
                  ['Prospective specification freeze', prospective['created_utc']],
                  ['Single benchmark opening event', opening['opened_utc']],
                  ['Snapshot ingestion completion', ingestion.get('completed_utc', ingestion.get('created_utc', 'see provenance/benchmark_ingestion_complete.json'))],
                  ['Annotation/identity freeze before V2 predictions', annotation['created_utc']],
                  ['Inference started', start['started_utc']],
                  ['Inference completed', inference['resources']['finished_utc']],
                  ['Descriptive evaluation completed', manifest['completed_utc']],
              ]),
              'The benchmark was opened once as a scientific evaluation after the prospective specification and independent pre-open review. '
              'Subsequent annotation and analysis used the unchanged snapshot, not a new source selection or repeated tuning cycle. '
              'No fresh experimental curation, web retrieval, model selection, threshold fitting, training or candidate-universe scoring occurred. '
              'The benchmark had already been examined under V1 and is not newly blinded evidence.', '',
              table(['Frozen artifact', 'SHA-256'], [
                  [relative, sha(STAGE / relative)] for relative in [
                      'PROSPECTIVE_SPECIFICATION.md', 'specification/evaluation_specification.json',
                      'provenance/prospective_freeze.json', 'provenance/preopen_source_manifest.json',
                      'inputs/benchmark_snapshot.json.gz', 'provenance/annotation_freeze.json',
                      'predictions/inference_manifest.json', 'results/evaluation_manifest.json']]),
              f"The pre-open source manifest verified {source_manifest['benchmark_payloads_verified_against_original_ledger']} original COPD-V2-BENCH1.0 payloads against their original ledger. "
              'The master benchmark, all prior stages and original model archives were inputs only. Exact dependency hashes and managed paths are preserved in the source manifest.', '',
              table(['Context', 'Seed', 'Original selected checkpoint SHA-256'], [
                  [row['context'], row['seed'], row['sha256']] for row in specification['checkpoints']]),
              'For each REF or ALT allele, each seed independently evaluates forward and nucleotide reverse-complement inputs. '
              'The two network float32 probabilities are promoted to float64 and averaged. The three seed means are averaged in order 104729, 130363, 155921. '
              'Delta is ensemble ALT minus ensemble REF. No seed, archive, architecture or input control was changed.', '',
              table(['Context', 'Frozen decimal threshold', 'Exact float64 encoding', 'Rule'], [
                  [context, values['decimal17g'], values['float64_hex'], 'max(REF, ALT) >= threshold']
                  for context, values in specification['thresholds'].items()]),
              'No allele-delta cutoff, internal-reference percentile, enrichment test, significance test, AUROC, AUPRC, FPR, specificity or sensitivity was calculated. '
              'Frozen region thresholds are not variant-level false-positive guarantees.', '',
              '## Population and identity accounting', '',
              table(['Population', 'Count'], [[row['population'], row['count']] for row in manifest['population_accounting']]),
              'Variant-level fractions deduplicate GRCh38 chromosome:position:REF:ALT identities. Every original assay-context row is retained separately. '
              'Reporter and endogenous-editing subsets overlap and must not be summed. A positive context is not erased by a null context. '
              'The 21 unresolved reporting groups are not asserted exact sequence identities. No LD proxy, approximate position, allele/strand rescue or post-open re-normalization was used.', '',
              'Sequence QC uses the frozen 2,001-bp construction: 1,000 upstream bases, the explicit allele, and shared downstream reference sequence cropped at the right edge. '
              'The allele starts at index 1000 for both REF and ALT; indels are not recentered. Genome REF, lengths, alphabet, reverse-complement involution, allele spans and independent sequence hashes were checked. '
              'Missing identities or incomplete predictions are unavailable, never negatives. All exact scoring inputs and source-row links are in the identity, sequence-QC and annotation tables.', '',
              '### Context-specific source states', '',
              table(['Experimental state', 'Frozen scope', 'Unique source keys', 'Exact V2-scorable keys'], [
                  [row['experimental_state'], row['mechanism_in_model_scope'], row['unique_benchmark_keys'], row['exact_V2_scorable_keys']]
                  for row in rows('results/context_state_populations.tsv')]),
              'These state populations overlap across assays; they are not a disjoint positive/negative partition.', '',
              '## Region-level recovery and frozen V1 comparison', '',
              '### Each method on its own available positive identities', '',
              recovery_table(recovery, 'evaluable_positive_denominator'),
              'Different coverage denominators cannot by themselves establish a paired performance gain. '
              'V1 forward and RC values are read from frozen R015 region gates and allele scores only. '
              'The historical combined candidate calls also imposed allele-delta/eligibility gates and are deliberately not used as region-only recovery. '
              'No V1 inference or threshold change was performed.', '',
              '### Same identities for all three methods', '',
              recovery_table(common, 'denominator'),
              '### Paired V1-to-V2 changes on shared positive identities', '',
              table(['Baseline', 'Context', 'Shared n', 'Retained', 'Newly recovered', 'Lost', 'Still missed'], [
                  [row['baseline'], row['context'], row['shared_positive_denominator'], row['retained_recovery'],
                   row['newly_recovered'], row['lost'], row['still_missed']] for row in paired]),
              '### Case-level interpretation of the modest paired change', '',
              case_interpretation, '',
              'Exact variant-level transitions are retained in `results/paired_V1_V2_variant_transitions.tsv`. '
              'The reused variants and assays are correlated; no independence-assuming significance test was performed.', '',
              f"### Coverage expansion ({manifest['new_V2_coverage_positive_count']} positive identities absent from V1 forward coverage)", '',
              table(['Variant', 'rsID', 'Locus', 'Enhancer', 'H3K27me3', 'Union'], [
                  [row['benchmark_variant_key'], row['rsids'], row['loci'], call(row['enhancer_recovered']),
                   call(row['h3k27me3_recovered']), call(row['union_recovered'])] for row in expansion]),
              'These are newly scorable V2 sequences, not paired V1 false negatives or paired gains. '
              'The separately preserved RC-coverage expansion table makes the analogous comparison with V1 RC.', '',
              '### All exact in-scope positives plus the unresolved positive', '',
              'Scores below are display-rounded; machine-readable tables retain full precision. E and H denote enhancer and H3K27me3-associated contexts.', '',
              table(['Variant / source key', 'rsID', 'Locus', 'E max', 'H max', 'E ALT−REF', 'H ALT−REF', 'E call', 'H call', 'Union', 'V1 coverage'], [
                  [row['benchmark_variant_key'], row['rsids'], row['loci'],
                   num(row.get('V2_enhancer_max_allele_score')), num(row.get('V2_h3k27me3_max_allele_score')),
                   num(row.get('V2_enhancer_delta')), num(row.get('V2_h3k27me3_delta')),
                   call(row.get('V2_enhancer_either_allele_region_call')), call(row.get('V2_h3k27me3_either_allele_region_call')),
                   call(row.get('V2_union_region_call')), row['V1_presence']] for row in positive_rows]),
              'rs2013701 is a pre-existing known case, not a blinded validation case. It remains visibly flagged in all relevant tables rather than being removed post hoc.', '',
              '## Allelic-effect direction: distinct from region recognition', '',
              'Only exact alleles, unambiguous frozen experimental ALT-minus-REF activity directions, and reasonably in-scope direct readouts enter this analysis. '
              'GWAS risk alleles, generic TF-binding changes, ambiguous MPRA ratio conventions and gene-expression anecdotes do not resolve direction. '
              'Positive endogenous allele editing provides partial regulatory support, not a direct bulk-lung chromatin label. '
              'The H3K27me3-associated score is a learned state association, not demonstrated repression causality; general expression effects are not automatically inverted.', '',
              table(['Unit', 'Subset', 'Context', 'Eligible n', 'Concordant', 'Discordant', 'Exact ties', 'Strict fraction', 'Secondary non-tied fraction'], [
                  [row['unit'], row['subset'], row['context'], row['denominator'], row['concordant'], row['discordant'], row['tied'],
                   fraction(row['concordant'], row['denominator']), fraction(row['concordant'], row['non_tied_denominator'])]
                  for row in direction_summaries]),
              'An exact unrounded model delta of zero is a tie: it remains in the strict denominator but is not concordant. '
              'No epsilon or allelic cutoff converts small deltas to ties. An observation-level percentage is not an independent-variant accuracy estimate.', '',
              table(['Assay', 'rsID', 'Cell context', 'Higher-activity allele', 'Expected ALT−REF sign', 'REF score', 'ALT score', 'Delta', 'Result'], [
                  [row['assay_id'], row['rsid'], row['cell_context'], row['higher_activity_grch38_allele'],
                   row['expected_ALT_minus_REF_sign'], num(row['ref_score']), num(row['alt_score']), num(row['delta']), row['concordance']]
                  for row in direction_rows]),
              '### Unique-variant consensus and locus accounting', '',
              table(['Variant', 'rsID', 'Locus', 'Context', 'Expected sign', 'Predicted sign', 'Delta', 'Result'], [
                  [row['variant_id'], row['rsids'], row['loci'], row['context'], row['expected_ALT_minus_REF_sign'],
                   row['predicted_sign'], num(row['delta']), row['concordance']] for row in consensus_rows]),
              table(['Locus', 'Context', 'Observations', 'Unique variants', 'Concordant', 'Discordant', 'Tied'], [
                  [row['locus'], row['context'], row['denominator'], row['unique_variants'], row['concordant'], row['discordant'], row['tied']]
                  for row in direction_loci]),
              'Conflicting experimental directions, if present, are retained at observation level and excluded only from the unique-variant consensus summary. '
              f"The frozen annotation map contained {annotation_rules['counts']['enhancer_direction_conflicting_variants']} such conflicting variants. "
              'Unevaluable directions and their reasons are preserved in `results/direction_unevaluable_observations.tsv`.', '',
              '## Source, assay, mechanism and context stratification', '',
              'Each cell below counts unique exact scorable positive variants within that stratum. Rows and strata overlap; small subsets are descriptive and do not establish tissue-specific validation.', '']
    strata_by_axis = defaultdict(lambda: defaultdict(dict))
    for row in stratified:
        strata_by_axis[row['axis']][row['stratum']][row['context']] = row
    for axis, groups in strata_by_axis.items():
        report.extend([f'### {axis}', '', table(['Stratum', 'Enhancer', 'H3K27me3-associated', 'Union'], [
            [group] + [fraction(cells[c]['recovered'], cells[c]['unique_positive_scorable_variants']) for c in CONTEXTS]
            for group, cells in groups.items()]), ''])
    report.extend(['## Continuous region and delta distributions', '',
                   'Quantiles use the prespecified NumPy float64 linear convention. They summarize this benchmark only; no internal reference distribution or outcome-conditioned percentile was constructed.', '',
                   table(['Subset', 'Context', 'Score', 'n', 'Minimum', 'Q25', 'Median', 'Q75', 'Maximum'], [
                       [row['subset'], row['context'], row['score'], row['n']] +
                       [num(row[name]) for name in ['minimum', 'q25', 'median', 'q75', 'maximum']] for row in distributions]), '',
                   '## Missing, unresolved and mechanistically out-of-model evidence', '',
                   table(['Unresolved source key', 'rsID', 'In-scope positive?', 'Frozen state', 'Reason'], [
                       [row['benchmark_variant_key'], row['rsids'], call(row['positive_in_scope']),
                        row['frozen_evidence_summary_state'], row['V2_missing_reason']] for row in missing]),
                   'The unresolved rs57658727 positive remains explicit and does not enter the exact positive denominator. '
                   'Its invalid/ambiguous deposited contrast is not rescued by another orientation or an approximate genomic identity. '
                   'All remaining unresolved keys are likewise unavailable, not negatives.', '',
                   '### Splice-only and other out-of-model positive variants', '',
                   table(['Variant', 'rsID', 'Locus', 'Mechanism', 'V2 sequence scorable?', 'Interpretation'], [
                       [row['benchmark_variant_key'], row['rsids'], row['loci'], row['mechanism_groups'], call(row['v2_scorable']),
                        'Sequence scoring retained for complete identity accounting; outside regulatory-positive recovery and direction denominators, not a false negative']
                       for row in all_variants if true(row['out_of_model_positive']) or true(row['splice_only'])]),
                   'NPNT splice-only evidence is preserved as out of model. The 267 contextual/excluded observations remain separate, including regional perturbations and contacts that do not establish an exact variant label.', '',
                   '### Frozen source denominator limitations', '',
                   table(['Study', 'Complete denominator pass', 'Preserved limitation'], [
                       [row['study_id'], row.get('complete_denominator_pass', ''),
                        row.get('denominator_failure_reason') or row.get('reason') or row.get('omitted_from_metrics_reason') or row.get('warning') or row.get('permitted_analysis')]
                       for row in denominator_rules]),
                   'GEO design/count completeness does not recover missing Castaldi significance labels or QC exclusions; unidentified hits are not nulls. '
                   'Gong null assay-context rows are not universal variant negatives, and missing/failed contexts do not become null. '
                   'No classification or enrichment statistic is manufactured from these incomplete denominators.', '',
                   '## Real-network RC invariance and compute provenance', '',
                   f"RC audit status: {invariance['status']}; {invariance['n_seed_audits']} seed audits and "
                   f"{invariance['n_ensemble_audits']} ensemble audits. Maximum absolute residual: "
                   f"{num(invariance['maximum_absolute_residual'], 17)}. Tolerance: "
                   f"abs(a−b) <= {invariance['atol']} + {invariance['rtol']}*abs(b).", '',
                   'The audit independently repeated nucleotide conversion and actual phase-I and phase-II network calls in RC-first/forward-second order. '
                   'It was not constructed by swapping stored probabilities or reversing phase-I feature columns.', '',
                   table(['Context', 'Unit', 'Seed', 'Probabilities', 'Failures', 'Maximum absolute residual'], [
                       [row['context'], row['unit'], row['seed'], row['n_probabilities'], row['n_failed'], num(row['max_absolute_difference'], 17)]
                       for row in invariance['audit_rows']]),
                   table(['Compute item', 'Recorded value'], [
                       ['Variants', inference['n_variants']], ['Allele sequences', inference['n_allele_sequences']],
                       ['Original selected checkpoints read, not copied/re-saved', inference['n_original_checkpoints']],
                       ['Variant/context scores', inference['n_context_variant_scores']],
                       ['Seed/context/variant scores', inference['n_seed_variant_scores']],
                       ['Phase-I representation', phase_one['representation']],
                       ['Phase-I feature shape per orientation', phase_one['shape_per_orientation']],
                       ['Phase-I actual sequence-orientation evaluations', phase_one['orientation_sequence_evaluations']],
                       ['Phase-I wall seconds', phase_one['wall_seconds']],
                       ['Full inference wall seconds', inference['resources']['wall_seconds']],
                       ['Inference peak RSS KiB', inference['resources']['peak_RSS_KiB']],
                       ['Descriptive evaluation wall seconds', manifest['resources']['wall_seconds']],
                   ]),
                   'Exact recorded inference invocation:', '', '```json', json.dumps(start['argv'], indent=2), '```', '',
                   'Recorded software/environment (from this run, not an inferred environment):', '',
                   '```json', json.dumps(environment, indent=2, sort_keys=True), '```', '',
                   '### Command and implementation-attempt ledger', ''])
    command_rows = []
    for path in sorted((STAGE / 'provenance/commands').glob('*.completed.json')):
        item = json.loads(path.read_text())
        command_rows.append([path.name, item.get('returncode'), item.get('wall_seconds'), json.dumps(item.get('argv', []))])
    report.extend(['', table(['Command record', 'Return code', 'Wall seconds', 'Recorded argv'], command_rows),
                   'Exclusive-output guards and pre-output assertion failures were implementation checks, not opportunities to alter scientific rules. '
                   'The annotation and sequence-preparation first attempts exposed the same unresolved reporting-key allele-order convention; sorting the display-key pair matched the already frozen source summary without rescuing sequence identity or changing eligibility. '
                   'Failed script versions and attempt ledgers are retained.', '',
                   'Annotation attempt ledger:', '', '```json', json.dumps(annotation_rules.get('attempt_ledger', []), indent=2), '```', ''])
    sequence_attempt_path = STAGE / 'provenance/sequence_preparation_attempts.json'
    if sequence_attempt_path.exists():
        report.extend(['Sequence-preparation attempt ledger:', '', '```json', sequence_attempt_path.read_text().strip(), '```', ''])
    report.extend(['## Independent validation, preservation and freeze', ''])
    report.extend([f"Independent pre-final scientific validation status: **{prefinal['status']}**. "
                   f'Receipt: `provenance/independent_validation_pre-final.json`; SHA-256 `{sha(prefinal_path)}`.', '',
                   'Authoritative final validation is recorded separately in `provenance/independent_final_validation.json` and its accompanying check table after this report is generated. '
                   'That receipt binds the completed report, results and append-only register records before the final seal. '
                   'Consult that receipt for the final status and exact verified scope. Its hash/status is not embedded in this report, avoiding a circular report-to-validator hash dependency; the finished report is not rewritten after final validation.'])
    report.extend(['', 'All original benchmark, pretraining, internal-training, internal-test, V1 and checkpoint artifacts are read-only inputs. '
                   'Preservation and append-only register assertions are documented in the final validation/freeze records. '
                   'No broader COPD candidate universe, 337-candidate reranking, fine-mapping or target-gene analysis was accessed or executed. '
                   'The six retained V2-C checkpoint archives were reused in place, not retrained, re-saved, slimmed or duplicated.', '',
                   'The final checksum ledger and freeze record bind this report, inputs, predictions, results, scripts, environment and provenance after validation. '
                   'The final freeze is intentionally created after this report and therefore its self-hash is not embedded here. '
                   'Only the three V2 registers may receive append-only stage updates. This evaluation request does not authorize a commit or push.', '',
                   '## Interpretation and remaining failure modes', '', improvement,
                   f'The common-identity union result is {union_conclusion}.',
                   absolute_recovery, '', case_interpretation, '',
                   f"Region misses remain explicitly listed in the per-variant and transition tables. Direction discordance ({strict['discordant']}/{strict['denominator']} observations) "
                   'must be reported separately from region recovery. Small allelic deltas have no validated decision threshold here. '
                   'Improved regional recognition, if observed, cannot establish better causal allelic ranking, repression mechanism, clinical prediction or generalization to unassayed COPD variants.', '',
                   'The major remaining limitations are incomplete experimental denominators; unresolved exact identities; cell/assay and episomal-versus-endogenous context mismatch; correlated loci and repeated contexts; '
                   'a small direction-resolved subset; no H3-specific experimental direction denominator; V1-dependent retrospective awareness including rs2013701; and out-of-model mechanisms such as splicing. '
                   'See `LIMITATIONS.md` for the explicit interpretation boundaries. This stage stops here; any future redesign is a separately versioned scientific project and cannot call this benchmark untouched validation.', ''])
    limitations = [
        '# Limitations and interpretation boundaries', '',
        f'Stage `external-benchmark-v2-1.0`; generated {generated}.', '',
        improvement, absolute_recovery,
        f'Common-identity union comparison: {union_conclusion}.', '',
        '1. This is a one-time, frozen-model, retrospective/descriptive external evaluation. The benchmark had already been examined under V1; rs2013701 was known in advance, not blinded validation.',
        '2. Experimental ascertainment, missing labels, unresolved identities and incomplete QC denominators prevent sensitivity, specificity, FPR, AUROC, AUPRC and enrichment claims. No such metrics or significance tests were calculated.',
        '3. Complete Castaldi GEO counts/design are not complete published positive/null labels. Unidentified reported hits remain unidentified; low counts, missing barcodes, missing results and unreported significance are not biological negatives.',
        '4. Gong nulls are assay-context outcomes, not universal variant negatives. Positive and null contexts can coexist for one variant; repeated contexts, shared variants and linked loci are not independent observations.',
        f"5. {len(missing)} frozen reporting groups lack exact usable identities, including the in-scope positive rs57658727. They are excluded as unavailable, not scored through proxies, approximate alleles or post-open identity rescue.",
        '6. Regulatory reporter, TF-binding and endogenous-expression support are not interchangeable with bulk-lung enhancer/H3K27me3 labels. Cell lineage, stimulation, insert orientation, insert length and episomal versus endogenous conditions limit transferability.',
        '7. H3K27me3-associated model scores indicate a learned regional association, not direct repression or causal gene regulation. No frozen assay supports a valid H3-specific allelic direction denominator; general expression/reporter directions were not automatically inverted.',
        f"8. Strict enhancer direction uses only {strict['denominator']} observations and {unique_direction['denominator']} unique consensus variants. Its {strict['concordant']} concordant observations do not establish general allelic accuracy; experimental contexts are correlated. Exact zero predictions count as non-concordant ties, and no delta cutoff is invented.",
        '9. Region-level recovery uses whether either allele passes an unchanged regional threshold. It is not evidence that ALT−REF deltas reproduce causal effects. No formal V1-to-V2 allelic-direction improvement analysis was prespecified, so improved regional recovery cannot be generalized to allelic improvement.',
        '10. V1 and V2 have different sequence coverage. Paired gains/losses use only shared exact positive identities. Newly V2-scorable identities are coverage expansion, not V1 negatives. V1 region-only gates are distinct from historical combined region/delta candidate calls.',
        '11. Splice-only NPNT and other out-of-model evidence are preserved, not counted as regulatory false negatives. Contextual contacts or regional perturbations are not promoted into exact variant labels.',
        '12. No internal-reference percentile distribution was prespecified or constructed. Descriptive score quantiles are internal to this benchmark and do not measure calibrated variant pathogenicity or general-population rarity.',
        '13. The hard post-opening firewall forbids retraining, retuning, recalibration, threshold changes, seed selection, model replacement, new controls or a return to internal-stage evaluation. Any future redesign needs a new scientific version and cannot reuse this benchmark as untouched validation.',
        '14. No broader COPD GWAS/candidate scoring, 337-candidate reranking, fine-mapping, target-gene analysis, commit or push was part of this request.', '',
        'Observed common-identity region-level conclusion: ' + improvement,
        case_interpretation,
        'Detailed quantitative results, missingness, per-locus directions, actual-network RC audit and validation receipt are in `EXTERNAL_BENCHMARK_V2_REPORT.md` and the machine-readable `results/` tables.', '']
    for path, lines in zip(targets, [report, limitations]):
        with path.open('w' if args.refresh else 'x') as handle:
            handle.write('\n'.join(lines))
    print(json.dumps({'status': 'PASS', 'generated_utc': generated,
                      'independent_pre_final_validation_status': prefinal['status'],
                      'independent_final_validation': 'Recorded separately after finished report; report remains unchanged after that validation',
                      'report_headline_metrics': {
                          'primary_region_recovery': primary,
                          'common_identity_union': common_union,
                          'strict_enhancer_direction': strict,
                          'unique_variant_enhancer_direction': unique_direction,
                          'exact_scorable_variants': flow['exact_V2_scorable_identity_keys_all_states_mechanisms'],
                          'exact_scorable_in_scope_positive_variants': flow['exact_V2_scorable_in_scope_positive_primary'],
                          'unresolved_reporting_groups': len(missing),
                          'coverage_expansion_positive_count': manifest['new_V2_coverage_positive_count'],
                      },
                      'reports': [{'path': str(path.relative_to(STAGE)), 'bytes': path.stat().st_size, 'sha256': sha(path)} for path in targets],
                      'scientific_results_recomputed': False, 'original_benchmark_reopened': False}, indent=2))


if __name__ == '__main__':
    main()
