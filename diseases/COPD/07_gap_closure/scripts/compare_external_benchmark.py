#!/usr/bin/env python3
"""Compare a hash-locked source benchmark with existing frozen V1 and RC scores.

No inference, training, new sequence construction or threshold selection occurs.
"""
import csv
import gzip
import hashlib
import json
import math
import platform
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
BASE = ROOT / 'diseases/COPD/07_gap_closure'
RES, PROV = BASE / 'results', BASE / 'provenance'
PREFIX = 'COPD-V2-BENCH'
SPEC_SHA = '46d3aac26b1f980d7c2dcd5034d6cc64566524ca7d99cc635f7d979e04f2de65'
MODEL_DIR = ROOT / 'diseases/COPD/04_modeling/results'
REGION = {'enhancer': .643623, 'silencer': .58505}
DELTA = {'enhancer': {'SNV': .05706318769999998, 'indel_or_complex': .04936093850000001},
         'silencer': {'SNV': .028802613899999996, 'indel_or_complex': .026183359500000003}}
CHECKS = []


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write(path, rows, first=()):
    fields = list(first) + sorted(set().union(*(set(r) for r in rows)) - set(first))
    with Path(path).open('w', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        for r in rows:
            w.writerow({k: '' if v is None else json.dumps(v, sort_keys=True) if isinstance(v, (list, dict)) else v
                        for k, v in r.items()})


def boolean(value):
    if str(value).lower() in ('true', '1'):
        return True
    if str(value).lower() in ('false', '0'):
        return False
    raise ValueError('Unexpected boolean: ' + repr(value))


def check(name, passed, detail=''):
    CHECKS.append({'check': name, 'status': 'PASS' if passed else 'FAIL', 'detail': str(detail)})
    if not passed:
        raise AssertionError(name + ': ' + str(detail))


def category(a, b):
    if a is None or b is None:
        return 'unevaluable'
    return 'stable_both' if a and b else 'forward_only' if a else 'RC_only' if b else 'negative_both'


def context(e, s):
    return 'both' if e and s else 'enhancer_only' if e else 'H3K27me3_only' if s else 'neither'


def rate(n, d):
    return n / d if d else None


def key_variant(r):
    return r.get('canonical_variant_id') or 'unresolved:' + r.get('rsid', '') + ':' + '/'.join(sorted(
        [r.get('tested_allele1', ''), r.get('tested_allele2', '')]))


def quantiles(values):
    if not values:
        return {'n': 0}
    values = sorted(values)
    def q(p):
        at = p * (len(values) - 1)
        a = int(at)
        b = min(a + 1, len(values) - 1)
        return values[a] + (values[b] - values[a]) * (at - a)
    return {'n': len(values), 'min': values[0], 'p25': q(.25), 'median': q(.5),
            'p75': q(.75), 'max': values[-1], 'mean': sum(values) / len(values)}


def main():
    freeze_path = PROV / f'{PREFIX}_benchmark_freeze.json'
    frozen = json.loads(freeze_path.read_text())
    check('Protocol lock unchanged', sha(PROV / f'{PREFIX}_analysis_specification.md') == SPEC_SHA)
    check('Benchmark marked external evaluation only', frozen['status'] == 'FROZEN_EXTERNAL_EVALUATION_DATA')
    for path, expected in frozen['frozen_files'].items():
        check('Frozen source/table unchanged: ' + path, sha(ROOT / path) == expected)
    started = now()
    run_path = PROV / f'{PREFIX}_comparison_manifest.json'
    previous = json.loads(run_path.read_text()) if run_path.exists() else None
    manifest = {'comparison_started_utc': started, 'benchmark_frozen_utc': frozen['frozen_utc'],
                'benchmark_freeze_sha256': sha(freeze_path), 'protocol_sha256': SPEC_SHA,
                'script_sha256': sha(__file__), 'python': platform.python_version(),
                'previous_runs': (previous.get('previous_runs', []) + [
                    {k: v for k, v in previous.items() if k != 'previous_runs'}]) if previous else [],
                'new_model_inference': False, 'new_sequences': False, 'thresholds_changed': False}
    check('Benchmark frozen before new model comparison', frozen['frozen_utc'] < started)
    run_path.write_text(json.dumps(manifest, indent=2) + '\n')
    master = read(ROOT / frozen['master_path'])
    score_path = MODEL_DIR / 'COPD-S4-R003_candidate_allele_scores.tsv.gz'
    meta_path = MODEL_DIR / 'COPD-S4-R004_prioritized_candidates.tsv.gz'
    rank_path = MODEL_DIR / 'COPD-S4-R010_THE_LIST.tsv'
    rc_path = RES / 'COPD-V2-RC-R001_reverse_complement_scores.tsv.gz'
    rc_check_path = RES / 'COPD-V2-RC-R004_all_scorable_comparison.tsv.gz'
    scores = {r['candidate_record_id']: r for r in read(score_path)}
    metadata = {r['candidate_record_id']: r for r in read(meta_path)}
    ranks = {r['candidate_record_id']: int(r['predicted_causal_priority_rank']) for r in read(rank_path)}
    rc_scores = {r['candidate_record_id']: r for r in read(rc_path)}
    rc_authority = {r['candidate_record_id']: r for r in read(rc_check_path)}
    check('Frozen universe sizes', (len(scores), len(metadata), len(ranks), len(rc_scores)) == (15303, 15389, 337, 15303))
    check('Forward and completed RC exact identity universes', set(scores) == set(rc_scores) == set(rc_authority))
    computed = {}
    for variant in sorted({r['canonical_variant_id'] for r in master if r['canonical_variant_id']}):
        record = {'canonical_variant_id': variant, 'v1_model_evaluable': variant in scores,
                  'rc_model_evaluable': variant in rc_scores,
                  'v1_candidate_status': variant in ranks if variant in metadata else None,
                  'v1_candidate_rank': ranks.get(variant),
                  'v1_match_status': 'exact_existing_scored_pair' if variant in scores else
                  'exact_existing_unscorable_record' if variant in metadata else 'no_frozen_V1_record'}
        if variant in scores:
            meta = metadata[variant]
            eligible = boolean(meta['causal_call_eligible'])
            record.update({'v1_class_group': meta['variant_class_group'], 'v1_call_eligible': eligible,
                           'v1_blacklisted': boolean(meta['blacklisted'])})
            for model in REGION:
                record[model + '_region_threshold'] = REGION[model]
                record[model + '_abs_delta_threshold'] = DELTA[model][meta['variant_class_group']]
                for orientation, source in [('forward', scores), ('rc', rc_scores)]:
                    raw = source[variant]
                    ref, alt = float(raw[model + '_ref_score']), float(raw[model + '_alt_score'])
                    delta = alt - ref
                    region = max(ref, alt)
                    call = eligible and region >= REGION[model] and abs(delta) >= record[model + '_abs_delta_threshold']
                    prefix = model + '_' + orientation
                    record.update({prefix + '_ref_score': ref, prefix + '_alt_score': alt,
                                   prefix + '_delta': delta,
                                   prefix + '_serialized_delta': float(raw[model + '_delta_alt_minus_ref']),
                                   prefix + '_region_score': region, prefix + '_abs_delta': abs(delta),
                                   prefix + '_call': call,
                                   prefix + '_region_gate': region >= REGION[model],
                                   prefix + '_delta_gate': abs(delta) >= record[model + '_abs_delta_threshold']})
                    check(f'{variant} {model} {orientation} original audit call fidelity',
                          call == boolean(rc_authority[variant][prefix + '_call']))
                    if orientation == 'forward':
                        check(f'{variant} {model} authoritative V1 call fidelity', call == boolean(meta['predicted_causal_' + model]))
                record[model + '_orientation_category'] = category(record[model + '_forward_call'], record[model + '_rc_call'])
            for orientation in ['forward', 'rc']:
                e, s = record['enhancer_' + orientation + '_call'], record['silencer_' + orientation + '_call']
                record[orientation + '_union_call'] = e or s
                record[orientation + '_model_context'] = context(e, s)
            record['union_orientation_category'] = category(record['forward_union_call'], record['rc_union_call'])
            check(f'{variant} frozen 337 membership fidelity', record['forward_union_call'] == (variant in ranks))
        else:
            record['union_orientation_category'] = 'unevaluable'
            record['enhancer_orientation_category'] = 'unevaluable'
            record['silencer_orientation_category'] = 'unevaluable'
        computed[variant] = record
    compared = []
    for assay in master:
        row = dict(assay)
        row.update(computed.get(assay['canonical_variant_id'], {
            'v1_model_evaluable': False, 'rc_model_evaluable': False,
            'v1_match_status': 'unresolved_assayed_identity', 'v1_candidate_status': None,
            'v1_candidate_rank': None, 'union_orientation_category': 'unevaluable',
            'enhancer_orientation_category': 'unevaluable', 'silencer_orientation_category': 'unevaluable'}))
        state, scope = assay['experimental_state'], assay['mechanism_in_model_scope']
        if scope == 'no':
            concordance = 'outside_model_scope_not_a_sequence_model_false_negative'
        elif not row['v1_model_evaluable']:
            concordance = 'model_unevaluable_not_a_false_negative'
        elif scope != 'yes':
            concordance = 'partial_or_unresolved_mechanism_descriptive_only'
        elif state == 'positive':
            concordance = 'reported_active_forward_recovered' if row['forward_union_call'] else 'reported_active_forward_not_recovered'
        elif state == 'null':
            concordance = 'assay_null_model_positive' if row['forward_union_call'] else 'assay_null_model_negative'
        else:
            concordance = 'experimental_label_not_evaluable_for_binary_performance'
        row['concordance_category'] = concordance
        row['enhancer_reporter_direction_agreement'] = 'unevaluable'
        if (row['v1_model_evaluable'] and state == 'positive' and
            row['assay_class'] in ('MPRA_allele_effect', 'conventional_reporter') and
            row.get('reported_direction_alt_minus_ref') in ('1', '-1')):
            direction = int(row['reported_direction_alt_minus_ref'])
            model_sign = (row['enhancer_forward_delta'] > 0) - (row['enhancer_forward_delta'] < 0)
            row['enhancer_reporter_direction_agreement'] = 'same' if model_sign == direction else 'opposite' if model_sign else 'zero'
        row['H3K27me3_activity_direction_interpretation'] = 'not_a_validated_reporter_activation_or_repression_direction'
        compared.append(row)
    forward = [{k: v for k, v in r.items() if '_rc_' not in k and not k.startswith('rc_') and 'orientation_category' not in k}
               for r in compared]
    write(RES / f'{PREFIX}-R007_V1_forward_comparison.tsv', forward,
          ('assay_id', 'study_id', 'rsid', 'canonical_variant_id', 'experimental_state', 'v1_match_status'))
    write(RES / f'{PREFIX}-R008_forward_RC_comparison.tsv', compared,
          ('assay_id', 'study_id', 'rsid', 'canonical_variant_id', 'experimental_state', 'union_orientation_category'))
    # Each source assay remains available; summary observations below are not independent replicates.
    summaries = []
    group_specs = [('study_assay_context', ('study_id', 'assay_class', 'cell_context', 'assay_orientation')),
                   ('mechanism', ('assay_class', 'mechanism_in_model_scope')),
                   ('study_locus', ('study_id', 'locus'))]
    for kind, fields in group_specs:
        grouped = defaultdict(list)
        for row in compared:
            grouped[tuple(row.get(k, '') for k in fields)].append(row)
        for values, group in sorted(grouped.items()):
            out = dict(zip(fields, values))
            out.update(summarize(group))
            out['summary_unit'] = kind
            summaries.append(out)
    for kind, name in [('study_assay_context', 'R009_study_specific_summary.tsv'),
                       ('mechanism', 'R010_mechanism_specific_summary.tsv'),
                       ('study_locus', 'R011_locus_specific_summary.tsv')]:
        write(RES / f'{PREFIX}-{name}', [r for r in summaries if r['summary_unit'] == kind])
    unique = variant_summary(compared)
    write(RES / f'{PREFIX}-R015_unique_variant_summary.tsv', unique,
          ('benchmark_variant_key', 'canonical_variant_id', 'rsids', 'evidence_summary_state'))
    failures = []
    for study in ['ALL_DEDUPLICATED'] + sorted({r['study_id'] for r in compared}):
        selected = {r['canonical_variant_id']: r for r in compared
                    if (study == 'ALL_DEDUPLICATED' or r['study_id'] == study) and
                    r['experimental_state'] == 'positive' and r['mechanism_in_model_scope'] == 'yes' and r['v1_model_evaluable']}
        for model in REGION:
            for orientation in ['forward', 'rc']:
                counts = Counter()
                for r in selected.values():
                    prefix = model + '_' + orientation
                    state = 'recovered' if r[prefix + '_call'] else 'original_eligibility_exclusion' if not r['v1_call_eligible'] else \
                            'both_gates_fail' if not r[prefix + '_region_gate'] and not r[prefix + '_delta_gate'] else \
                            'region_only_fail' if not r[prefix + '_region_gate'] else 'delta_only_fail'
                    counts[state] += 1
                failures.append({'study_id': study, 'model': model, 'orientation': orientation,
                                 'in_scope_positive_scored_unique_variants': len(selected),
                                 **{k: counts[k] for k in ['recovered', 'original_eligibility_exclusion',
                                                         'both_gates_fail', 'region_only_fail', 'delta_only_fail']}})
    write(RES / f'{PREFIX}-R017_positive_case_gate_failures.tsv', failures)
    orientation_rows = []
    for subset, selected in [('all_exact_scored_variants', [r for r in unique if r['v1_model_evaluable']]),
                             ('any_in_scope_positive_scored_variants', [r for r in unique if r['v1_model_evaluable'] and r['any_in_scope_positive_assay']])]:
        for model in ['enhancer', 'silencer', 'union']:
            orientation_rows.append({'subset': subset, 'model': model, 'unique_variants': len(selected),
                                     **{cat: sum(r[model + '_orientation_category'] == cat for r in selected)
                                        for cat in ['stable_both', 'forward_only', 'RC_only', 'negative_both']}})
    write(RES / f'{PREFIX}-R018_orientation_recovery_summary.tsv', orientation_rows)
    distributions = []
    grouped = defaultdict(list)
    for r in compared:
        if r['v1_model_evaluable']:
            grouped[(r['study_id'], r['assay_class'], r['cell_context'], r.get('assay_orientation', ''), r['experimental_state'])].append(r)
    for values, group in sorted(grouped.items()):
        for model in REGION:
            for metric in ['region_score', 'abs_delta']:
                vals = list({r['canonical_variant_id']: r[f'{model}_forward_{metric}'] for r in group}.values())
                distributions.append(dict(zip(['study_id', 'assay_class', 'cell_context', 'assay_orientation', 'experimental_state'], values),
                                          model=model, score=metric, **quantiles(vals),
                                          interpretation='descriptive_model_covered_source_stratum_not_independent_variants'))
    write(RES / f'{PREFIX}-R013_score_distributions.tsv', distributions)
    flow = []
    for study in sorted({r['study_id'] for r in compared}):
        rows = [r for r in compared if r['study_id'] == study]
        flow.append(dict(study_id=study, **summarize(rows),
                         match_status_counts=dict(Counter(r['v1_match_status'] for r in rows))))
    write(RES / f'{PREFIX}-R014_benchmark_coverage_flow.tsv', flow, ('study_id',))
    gates = performance_gates(compared, read(RES / f'{PREFIX}-R005_assayed_denominator_audit.tsv'))
    write(RES / f'{PREFIX}-R012_quantitative_metric_eligibility.tsv', gates)
    check('Every frozen assay preserved in forward/RC comparison', len(compared) == len(master) and
          [r['assay_id'] for r in compared] == [r['assay_id'] for r in master])
    check('No nonexact model comparison', all(not r['v1_model_evaluable'] or r['identity_status'] == 'exact' for r in compared))
    check('Out-of-scope rows never called sequence model false negatives', all(
        r['concordance_category'] == 'outside_model_scope_not_a_sequence_model_false_negative'
        for r in compared if r['mechanism_in_model_scope'] == 'no'))
    write(RES / f'{PREFIX}-R016_analysis_validation.tsv', CHECKS, ('check', 'status', 'detail'))
    manifest.update({'completed_utc': now(), 'assay_rows': len(compared), 'unique_variant_records': len(unique),
                     'model_covered_unique_exact_variants': sum(r['v1_model_evaluable'] for r in computed.values()),
                     'analysis_checks': len(CHECKS),
                     'inputs': {str(p.relative_to(ROOT)): sha(p) for p in [score_path, meta_path, rank_path, rc_path, rc_check_path]},
                     'outputs': {str(p.relative_to(ROOT)): sha(p) for p in RES.glob(f'{PREFIX}-R*.tsv')},
                     'all_variant_summary': summarize(compared),
                     'unique_any_positive_in_scope_summary': summarize_unique_positive(unique)})
    run_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + '\n')
    print(json.dumps({k: manifest[k] for k in ['assay_rows', 'unique_variant_records',
          'model_covered_unique_exact_variants', 'analysis_checks', 'unique_any_positive_in_scope_summary']}, indent=2))


def summarize(group):
    out = {'assay_rows': len(group), 'unique_variant_keys': len({key_variant(r) for r in group}),
           'unique_reported_rsids': len({r['rsid'] for r in group if r.get('rsid')}),
           'unique_exact_variants': len({r['canonical_variant_id'] for r in group if r['canonical_variant_id']}),
           'model_evaluable_assay_rows': sum(r['v1_model_evaluable'] for r in group),
           'model_evaluable_unique_variants': len({r['canonical_variant_id'] for r in group if r['v1_model_evaluable']}),
           'in_scope_positive_case_series_only': True}
    for state in ['positive', 'null', 'ambiguous', 'conflicting', 'unavailable', 'unevaluable']:
        rows = [r for r in group if r['experimental_state'] == state]
        out[state + '_assay_rows'] = len(rows)
        out[state + '_unique_variant_keys'] = len({key_variant(r) for r in rows})
        covered = {r['canonical_variant_id']: r for r in rows if r['v1_model_evaluable']}
        out[state + '_model_covered_unique_variants'] = len(covered)
        out[state + '_forward_union_call_unique_variants'] = sum(r['forward_union_call'] for r in covered.values())
        out[state + '_rc_union_call_unique_variants'] = sum(r['rc_union_call'] for r in covered.values())
    active = {r['canonical_variant_id']: r for r in group if r['experimental_state'] == 'positive' and
              r['mechanism_in_model_scope'] == 'yes' and r['v1_model_evaluable']}
    out['in_scope_positive_model_evaluable_unique'] = len(active)
    for model, field in [('enhancer', 'enhancer_forward_call'), ('H3K27me3', 'silencer_forward_call'), ('union', 'forward_union_call')]:
        count = sum(r[field] for r in active.values())
        out['in_scope_positive_forward_' + model + '_count'] = count
        out['in_scope_positive_forward_' + model + '_case_recovery_fraction'] = rate(count, len(active))
    for cat in ['stable_both', 'forward_only', 'RC_only', 'negative_both']:
        out['in_scope_positive_union_' + cat] = sum(r['union_orientation_category'] == cat for r in active.values())
    out['direction_agreement_positive_assay_counts'] = dict(Counter(r['enhancer_reporter_direction_agreement'] for r in group
                                                                  if r['experimental_state'] == 'positive'))
    out['independence_limit'] = 'Repeated assay contexts and linked variants are not independent biological replicates; no binomial CI'
    return out


def variant_summary(rows):
    groups = defaultdict(list)
    for r in rows:
        groups[key_variant(r)].append(r)
    result = []
    for key, group in sorted(groups.items()):
        in_scope = [r for r in group if r['mechanism_in_model_scope'] == 'yes']
        labels = {r['experimental_state'] for r in in_scope}
        directions = {r['reported_direction_alt_minus_ref'] for r in in_scope if r['experimental_state'] == 'positive' and
                      r.get('reported_direction_alt_minus_ref') in ('1', '-1')}
        state = 'conflicting' if 'conflicting' in labels or {'positive', 'null'} <= labels or len(directions) > 1 else \
                'positive' if 'positive' in labels else 'null' if labels == {'null'} else 'unevaluable'
        exemplar = group[0]
        out = {k: v for k, v in exemplar.items() if k.startswith(('v1_', 'enhancer_', 'silencer_', 'forward_', 'rc_')) or
               k in ('canonical_variant_id', 'union_orientation_category')}
        out.pop('enhancer_reporter_direction_agreement', None)
        out.update({'benchmark_variant_key': key, 'rsids': ';'.join(sorted({r['rsid'] for r in group})),
                    'study_ids': ';'.join(sorted({r['study_id'] for r in group})),
                    'loci': ';'.join(sorted({r.get('locus', '') for r in group})),
                    'assay_classes': ';'.join(sorted({r['assay_class'] for r in group})),
                    'assay_rows': len(group), 'evidence_summary_state': state,
                    'any_in_scope_positive_assay': 'positive' in labels,
                    'in_scope_positive_assay_rows': sum(r['experimental_state'] == 'positive' for r in in_scope),
                    'in_scope_null_assay_rows': sum(r['experimental_state'] == 'null' for r in in_scope),
                    'opposite_activity_directions_across_contexts': len(directions) > 1,
                    'positive_reporter_direction_assay_counts': dict(Counter(r['enhancer_reporter_direction_agreement']
                        for r in group if r['experimental_state'] == 'positive' and r['assay_class'] in ('MPRA_allele_effect', 'conventional_reporter'))),
                    'evidence_summary_interpretation': 'Positive/null coexistence denotes context heterogeneity, not necessarily failed replication'})
        result.append(out)
    return result


def summarize_unique_positive(rows):
    selected = [r for r in rows if r['any_in_scope_positive_assay']]
    covered = [r for r in selected if r['v1_model_evaluable']]
    return {'any_positive_unique_variant_keys': len(selected), 'model_covered': len(covered),
            'model_unevaluable': len(selected) - len(covered),
            'forward_union_recovered': sum(r['forward_union_call'] for r in covered),
            'rc_union_recovered': sum(r['rc_union_call'] for r in covered),
            'orientation_categories': dict(Counter(r['union_orientation_category'] for r in covered)),
            'context_heterogeneous_among_covered': sum(r['evidence_summary_state'] == 'conflicting' for r in covered),
            'interpretation': 'Any-reported-positive descriptive case series, not sensitivity or independent validation replicates'}


def performance_gates(rows, denominators):
    # A complete observed table and a complete prespecified tested/QC denominator are different.
    approved = {r['study_id'] for r in denominators if str(r.get('complete_assayed_denominator_valid', '')).lower() == 'true'}
    grouped = defaultdict(list)
    for r in rows:
        grouped[(r['study_id'], r['assay_class'], r['cell_context'], r.get('assay_orientation', ''))].append(r)
    gates = []
    for (study, assay, cell, orientation), group in sorted(grouped.items()):
        reasons = []
        if study not in approved:
            reasons.append('complete_prespecified_tested_and_QC_denominator_not_verified_or_selected_case_series')
        if any(r['experimental_state'] not in ('positive', 'null') for r in group):
            reasons.append('missing_or_ambiguous_experimental_labels')
        if any(r['mechanism_in_model_scope'] != 'yes' for r in group):
            reasons.append('not_entirely_in_scope_allele_activity_assays')
        if any(not r['v1_model_evaluable'] for r in group):
            reasons.append('incomplete_exact_existing_V1_sequence_coverage')
        positives = {key_variant(r) for r in group if r['experimental_state'] == 'positive'}
        nulls = {key_variant(r) for r in group if r['experimental_state'] == 'null'}
        if not positives or not nulls:
            reasons.append('both_positive_and_null_classes_not_available')
        if positives & nulls:
            reasons.append('within_context_conflicting_variant_labels')
        if not reasons:
            raise RuntimeError('A complete valid panel was found; implement prespecified discrimination and cluster inference before completion')
        gates.append({'study_id': study, 'assay_class': assay, 'cell_context': cell, 'assay_orientation': orientation,
                      'positive_variant_keys': len(positives), 'null_variant_keys': len(nulls),
                      'full_panel_metrics_valid': False, 'AUROC': 'not_calculated', 'average_precision': 'not_calculated',
                      'sensitivity_specificity_FPR': 'not_calculated', 'confidence_interval': 'not_calculated',
                      'reason': ';'.join(reasons), 'permitted_output': 'descriptive_evaluable_case_recovery_counts_and_score_distributions'})
    return gates


if __name__ == '__main__':
    main()
