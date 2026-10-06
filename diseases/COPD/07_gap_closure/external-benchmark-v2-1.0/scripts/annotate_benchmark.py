#!/usr/bin/env python3
"""Freeze descriptive annotations from the single immutable benchmark snapshot.

No model, prediction, candidate universe, original benchmark, or scientific
runtime is opened. Full original rows remain in the hash-bound snapshot; row
locators and canonical-JSON hashes bind each compact annotation to its evidence.
"""
from collections import Counter, defaultdict
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path

STAGE = Path(__file__).resolve().parents[1]
MASTER = 'COPD-V2-BENCH-R003_frozen_benchmark_master.tsv'
CONTEXT = 'COPD-V2-BENCH-R004B_contextual_and_excluded_evidence.tsv'
UNIQUE = 'COPD-V2-BENCH-R015_unique_variant_summary.tsv'
REPORTER = {'MPRA_allele_effect', 'conventional_reporter'}
IN_SCOPE = {'yes', 'partial'}


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def true(value):
    return str(value).strip().lower() == 'true'


def false(value):
    return str(value).strip().lower() == 'false'


def row_hash(row):
    return hashlib.sha256(json.dumps(row, sort_keys=True, separators=(',', ':'),
                                     ensure_ascii=True).encode()).hexdigest()


def identity(row):
    reasons = []
    if row.get('identity_status') != 'exact':
        reasons.append('frozen_identity_status_not_exact')
    if not true(row.get('reference_verified')):
        reasons.append('frozen_reference_not_verified')
    if false(row.get('exact_identity_eligible')):
        reasons.append('frozen_exact_identity_eligible_false')
    if false(row.get('construct_exact_allele_identity_valid')):
        reasons.append('frozen_construct_exact_identity_false')
    parts = [row.get(k, '') for k in ['grch38_chrom', 'grch38_pos', 'grch38_ref', 'grch38_alt']]
    if not all(parts) or ':'.join(parts) != row.get('canonical_variant_id'):
        reasons.append('incomplete_or_inconsistent_frozen_canonical_key')
    return not reasons, ';'.join(reasons)


def benchmark_key(row):
    if row.get('canonical_variant_id'):
        return row['canonical_variant_id']
    allele1, allele2 = sorted([row['tested_allele1'], row['tested_allele2']])
    return f"unresolved:{row['rsid']}:{allele1}/{allele2}"


def mechanism_group(row):
    assay = row.get('assay_class', '')
    if assay == 'splicing':
        return 'splice_only_out_of_model'
    if assay in REPORTER:
        return 'enhancer_like_reporter'
    if assay == 'endogenous_allele_editing':
        return 'endogenous_expression_support_partial'
    if assay == 'TF_binding':
        return 'TF_binding_support_partial'
    return 'contextual_or_other_unresolved'


def direction(row, exact):
    reasons = []
    if not exact:
        reasons.append('identity_not_exact')
    if row['experimental_state'] != 'positive':
        reasons.append('not_frozen_positive_observation')
    if row['mechanism_in_model_scope'] not in IN_SCOPE:
        reasons.append('mechanism_out_of_scope_or_unresolved')
    if row['assay_class'] not in REPORTER | {'endogenous_allele_editing'}:
        reasons.append('not_direct_regulatory_activity_readout')
    if not true(row.get('direction_identity_resolved')):
        reasons.append('frozen_direction_identity_unresolved')
    higher = row.get('higher_activity_grch38_allele', '')
    sign = row.get('reported_direction_alt_minus_ref', '')
    ref, alt = row.get('grch38_ref', ''), row.get('grch38_alt', '')
    if not higher or higher not in {ref, alt} or sign not in {'-1', '1'}:
        reasons.append('explicit_frozen_activity_allele_or_signed_direction_unavailable')
    if not reasons:
        assert int(sign) == (1 if higher == alt else -1), row['assay_id']
        # Every currently direction-resolved observation has an identity transform.
        # Do not synthesize new complements or reconstruct an absent direction.
        assert row['allele_transform'] == 'identity', row['assay_id']
        assert {row['tested_allele1'], row['tested_allele2']} == {ref, alt}, row['assay_id']
    return not reasons, int(sign) if not reasons else '', ';'.join(reasons)


def write_tsv(path, rows):
    columns = list(rows[0])
    with path.open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def joined(rows, key):
    return ';'.join(sorted({str(row.get(key, '')) for row in rows if str(row.get(key, ''))}))


def file_record(path):
    return {'path': str(path.relative_to(STAGE)), 'bytes': path.stat().st_size, 'sha256': sha(path)}


def main():
    outputs = [STAGE / p for p in ['inputs/observation_annotations.tsv',
                                   'inputs/variant_annotations.tsv',
                                   'inputs/contextual_annotations.tsv',
                                   'provenance/annotation_rules.json']]
    assert not any(path.exists() for path in outputs), 'Refuse to overwrite annotations'
    assert not any(path.is_file() for path in (STAGE / 'predictions').rglob('*')), 'Annotations must precede predictions'
    freeze_path = STAGE / 'provenance/prospective_freeze.json'
    freeze = json.loads(freeze_path.read_text())
    for item in freeze['files']:
        path = STAGE / item['path']
        assert sha(path) == item['sha256'] and path.stat().st_size == item['bytes']
    snapshot_path = STAGE / 'inputs/benchmark_snapshot.json.gz'
    with gzip.open(snapshot_path, 'rt') as handle:
        snapshot = json.load(handle)
    master = snapshot['tables'][MASTER]['rows']
    contextual = snapshot['tables'][CONTEXT]['rows']
    unique = snapshot['tables'][UNIQUE]['rows']
    assert (len(master), len(contextual), len(unique)) == (14025, 267, 1731)
    assert len({row['assay_id'] for row in master}) == len(master)
    original_unique = {row['benchmark_variant_key']: row for row in unique}
    assert len(original_unique) == len(unique)
    observations = []
    groups = defaultdict(list)
    evidence_keys = ['study_id', 'rsid', 'resolved_rsid', 'locus', 'cell_context',
                     'assay_class', 'assay_orientation', 'experimental_state',
                     'mechanism_in_model_scope', 'mechanism_scope_rationale',
                     'evidence_role', 'identity_status', 'canonical_variant_id',
                     'grch38_chrom', 'grch38_pos', 'grch38_ref', 'grch38_alt',
                     'reference_verified', 'exact_identity_eligible',
                     'construct_exact_allele_identity_valid', 'tested_allele1',
                     'tested_allele2', 'allele_transform',
                     'higher_activity_grch38_allele', 'reported_direction_alt_minus_ref',
                     'direction_identity_resolved', 'reported_effect_direction',
                     'source_locator', 'source_url', 'notes']
    for index, row in enumerate(master):
        exact, reason = identity(row)
        key = benchmark_key(row)
        assert key in original_unique, key
        positive = row['experimental_state'] == 'positive'
        scoped = row['mechanism_in_model_scope'] in IN_SCOPE
        eligible, sign, direction_reason = direction(row, exact)
        out = {'observation_id': f'R003:{row["assay_id"]}', 'source_table': MASTER,
               'snapshot_row_index': index, 'assay_id': row['assay_id'],
               'variant_id': row['canonical_variant_id'] if exact else '',
               'benchmark_variant_key': key}
        out.update({key: row.get(key, '') for key in evidence_keys})
        out.update(identity_exact=exact, identity_exclusion_reason=reason,
                   experimental_positive=positive, positive_in_scope=positive and scoped,
                   reporter_assay=row['assay_class'] in REPORTER,
                   endogenous_assay=row['assay_class'] == 'endogenous_allele_editing',
                   reporter_positive=positive and scoped and row['assay_class'] in REPORTER,
                   endogenous_positive=positive and scoped and row['assay_class'] == 'endogenous_allele_editing',
                   splice_assay=row['assay_class'] == 'splicing',
                   out_of_model_positive=positive and not scoped,
                   mechanism_group=mechanism_group(row),
                   enhancer_direction_eligible=eligible, enhancer_expected_sign=sign,
                   h3k27me3_direction_eligible=False, h3k27me3_expected_sign='',
                   direction_exclusion_reason=direction_reason,
                   h3k27me3_direction_exclusion_reason='No frozen direct H3K27me3 or established repressive-state allelic readout; general activity direction is not automatically inverted',
                   preexisting_known_case=row['rsid'] == 'rs2013701',
                   frozen_row_json_sha256=row_hash(row))
        observations.append(out)
        groups[key].append(out)
    variants = []
    for key in sorted(groups):
        rows = groups[key]
        original = original_unique[key]
        ids = {row['variant_id'] for row in rows if row['identity_exact']}
        assert len(ids) <= 1
        variant_id = next(iter(ids)) if ids else ''
        positive = any(row['positive_in_scope'] for row in rows)
        assert positive == true(original['any_in_scope_positive_assay']), key
        assert len(rows) == int(original['assay_rows']), key
        assert variant_id == original['canonical_variant_id'], key
        signs = {row['enhancer_expected_sign'] for row in rows if row['enhancer_direction_eligible']}
        conflict = len(signs) > 1
        consensus = next(iter(signs)) if len(signs) == 1 else ''
        chrom, pos, ref, alt = variant_id.split(':') if variant_id else ('', '', '', '')
        variant_type = 'SNV' if len(ref) == len(alt) == 1 else 'indel' if ref and alt and len(ref) != len(alt) else 'MNV' if ref and alt else 'unresolved'
        variants.append({
            'variant_id': variant_id, 'benchmark_variant_key': key,
            'canonical_variant_id': variant_id, 'identity_exact': bool(variant_id),
            'chrom': chrom, 'pos1': pos, 'ref': ref, 'alt': alt, 'variant_type': variant_type,
            'rsids': joined(rows, 'rsid'), 'loci': joined(rows, 'locus'),
            'study_ids': joined(rows, 'study_id'), 'cell_contexts': joined(rows, 'cell_context'),
            'assay_classes': joined(rows, 'assay_class'), 'mechanism_groups': joined(rows, 'mechanism_group'),
            'mechanism_in_model_scope_values': joined(rows, 'mechanism_in_model_scope'),
            'experimental_states': joined(rows, 'experimental_state'),
            'frozen_evidence_summary_state': original['evidence_summary_state'],
            'experimental_positive': any(row['experimental_positive'] for row in rows),
            'positive_in_scope': positive,
            'reporter_positive': any(row['reporter_positive'] for row in rows),
            'endogenous_positive': any(row['endogenous_positive'] for row in rows),
            'splice_positive': any(row['experimental_positive'] and row['splice_assay'] for row in rows),
            'splice_only': any(row['splice_assay'] for row in rows) and not positive,
            'out_of_model_positive': any(row['out_of_model_positive'] for row in rows),
            'contextual_positive_null_conflict': any(row['experimental_state'] == 'positive' for row in rows) and any(row['experimental_state'] == 'null' for row in rows),
            'enhancer_direction_eligible': bool(signs) and not conflict,
            'enhancer_expected_sign': consensus,
            'enhancer_direction_consensus_sign': consensus,
            'enhancer_direction_conflict': conflict,
            'enhancer_direction_observation_count': sum(row['enhancer_direction_eligible'] for row in rows),
            'h3k27me3_direction_eligible': False, 'h3k27me3_expected_sign': '',
            'h3k27me3_direction_consensus_sign': '', 'h3k27me3_direction_conflict': False,
            'h3k27me3_direction_observation_count': 0,
            'observation_count': len(rows),
            'positive_observation_count': sum(row['experimental_positive'] for row in rows),
            'in_scope_positive_observation_count': sum(row['positive_in_scope'] for row in rows),
            'null_observation_count': sum(row['experimental_state'] == 'null' for row in rows),
            'preexisting_known_case': any(row['preexisting_known_case'] for row in rows),
            'frozen_R015_row_json_sha256': row_hash(original),
        })
    contexts = []
    for index, row in enumerate(contextual):
        out = {'observation_id': f'R004B:{index:04d}', 'source_table': CONTEXT,
               'snapshot_row_index': index, 'variant_id': '', 'identity_exact': False,
               'sequence_scoring_population': False,
               'population_exclusion_reason': 'Frozen contextual/excluded evidence; not an exact master variant-identity observation',
               'experimental_positive_context': row.get('experimental_state') == 'positive',
               'frozen_row_json_sha256': row_hash(row)}
        # Preserve every original contextual field, including exclusions/construct QC.
        out.update(row)
        contexts.append(out)
    positives = {row['variant_id'] for row in variants if row['positive_in_scope'] and row['identity_exact']}
    original_positives = {row['canonical_variant_id'] for row in unique if true(row['any_in_scope_positive_assay']) and row['canonical_variant_id']}
    assert positives == original_positives and len(positives) == 39
    assert len(variants) == 1731 and sum(row['identity_exact'] for row in variants) == 1710
    assert sum(row['positive_in_scope'] for row in variants) == 40
    assert sum(row['enhancer_direction_eligible'] for row in observations) == 11
    assert sum(row['enhancer_direction_eligible'] for row in variants) == 8
    assert not any(row['enhancer_direction_conflict'] for row in variants)
    assert sum(row['identity_exact'] for row in observations) == 13899
    for path, rows in zip(outputs[:3], [observations, variants, contexts]):
        write_tsv(path, rows)
    rules = {
        'status': 'PASS', 'created_utc': datetime.now(timezone.utc).isoformat(),
        'stage': 'external-benchmark-v2-1.0',
        'source': file_record(snapshot_path), 'script': file_record(Path(__file__).resolve()),
        'prospective_freeze_sha256': sha(freeze_path),
        'before_any_V2_predictions': True, 'V2_prediction_files_opened': False,
        'original_benchmark_files_reopened': False, 'candidate_universe_opened': False,
        'identity_rule': 'Frozen exact identity, reference_verified True, canonical coordinates/alleles consistent; explicit False identity/construct eligibility vetoes. Sequence validation is separate and cannot rescue identities.',
        'scope_values': sorted(IN_SCOPE),
        'positive_rule': 'experimental_state == positive and mechanism_in_model_scope in yes/partial, preserving frozen R015 any_in_scope_positive_assay for every unique benchmark key',
        'reporter_assay_classes': sorted(REPORTER),
        'endogenous_assay_class': 'endogenous_allele_editing',
        'direction_rule': 'Positive exact in-scope direct reporter or isolated endogenous expression observation, frozen direction_identity_resolved True, explicit higher_activity_grch38_allele in REF/ALT and matching signed ALT-minus-REF direction. No narrative rescue, GWAS/risk effect, missing MPRA ratio convention, or TF-binding-to-activity inference.',
        'h3k27me3_direction_rule': 'No frozen assay establishes a directly measured H3K27me3 or explicitly demonstrated repressive-state allele direction; denominator is empty rather than inverting general expression/reporter readouts.',
        'h3k27me3_eligible_assay_ids': [],
        'unique_variant_direction_rule': 'All eligible experimental signs must agree; conflicting variants excluded only from the consensus summary and retained in observation-level tables.',
        'contextual_rule': 'All 267 original R004B rows retained without promotion to exact-variant labels or negatives.',
        'unresolved_rule': '21 unresolved frozen benchmark keys retained as source reporting groups, not asserted sequence identities. variant_id blank; benchmark_variant_key and unique observation IDs preserve complete accounting.',
        'full_original_row_preservation': 'All exact original rows remain in the immutable benchmark_snapshot.json.gz. Annotation source_table and zero-based snapshot_row_index locate each row; frozen_row_json_sha256 uses json.dumps(sort_keys=True,separators=(comma,colon),ensure_ascii=True). Compact annotation tables do not duplicate the 62 MB master JSON.',
        'counts': {
            'master_observations': len(observations), 'contextual_observations': len(contexts),
            'unique_frozen_benchmark_keys': len(variants), 'exact_variant_identities': 1710,
            'unresolved_frozen_benchmark_keys': 21, 'exact_master_observations': 13899,
            'in_scope_positive_benchmark_keys_including_unresolved': 40,
            'exact_in_scope_positive_variants_before_sequence_QC': len(positives),
            'enhancer_direction_observations': 11, 'enhancer_direction_unique_consensus_variants': 8,
            'enhancer_direction_conflicting_variants': 0,
            'h3k27me3_direction_observations': 0,
            'experimental_state_observations': dict(Counter(row['experimental_state'] for row in observations)),
        },
        'direction_resolved_observations': [
            {key: row[key] for key in ['observation_id', 'assay_id', 'variant_id', 'study_id',
                                      'cell_context', 'assay_class', 'higher_activity_grch38_allele',
                                      'enhancer_expected_sign', 'source_locator']}
            for row in observations if row['enhancer_direction_eligible']
        ],
        'crosschecks': {
            'R015_unique_keys_and_observation_counts_match': True,
            'R015_all_positive_flags_match': True,
            'R015_exact_positive_key_set_match': True,
            'all_14025_master_and_267_contextual_rows_accounted': True,
        },
        'attempt_ledger': [
            {'attempt': 1, 'status': 'FAILED_BEFORE_OUTPUTS',
             'command': 'python3 diseases/COPD/07_gap_closure/external-benchmark-v2-1.0/scripts/annotate_benchmark.py',
             'exit_code': 1,
             'error': 'AssertionError: unresolved:rs201292172:ACAC/-',
             'cause': 'Frozen R015 unresolved reporting keys sort tested alleles; first implementation retained reported allele order.',
             'repair': 'Sort the two tested-allele labels only when recreating an unresolved R015 reporting key; no sequence identity rescued, no eligibility/direction rule changed.',
             'failed_script': file_record(STAGE / 'scripts/history/annotate_benchmark_attempt001.py')},
            {'attempt': 2, 'status': 'PASS', 'exit_code': 0,
             'command': 'python3 diseases/COPD/07_gap_closure/external-benchmark-v2-1.0/scripts/annotate_benchmark.py'},
        ],
        'outputs': [file_record(path) for path in outputs[:3]],
    }
    with outputs[3].open('x') as handle:
        json.dump(rules, handle, indent=2, sort_keys=True)
        handle.write('\n')
    print(json.dumps({'status': 'PASS', 'counts': rules['counts'],
                      'annotation_rules': file_record(outputs[3])}, indent=2))


if __name__ == '__main__':
    main()
