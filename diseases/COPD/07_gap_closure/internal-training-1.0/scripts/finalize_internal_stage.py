#!/usr/bin/env python3
"""Freeze completed internal-stage artifacts after independent output validation."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
from datetime import datetime, timezone


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


def load(path):
    return json.loads(path.read_text())


def write(path, value):
    with path.open('x') as f:
        json.dump(value, f, indent=2, sort_keys=True, allow_nan=False)
        f.write('\n')


def verify_refs(records, base):
    for rec in records:
        path = Path(rec['path'])
        if not path.is_absolute():
            path = base / path
        if sha(path) != rec['sha256'] or ('bytes' in rec and path.stat().st_size != rec['bytes']):
            raise RuntimeError('Producer-manifest hash/size mismatch: ' + str(path))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', required=True, type=Path)
    args = parser.parse_args()
    stage = args.stage.resolve()
    repo = stage.parents[3]
    provenance = stage / 'provenance'
    if (provenance / 'freeze.json').exists() or (provenance / 'artifact_checksums.tsv').exists():
        raise RuntimeError('Stage already frozen; overwriting is prohibited')
    decision = load(stage / 'results/chr7_evaluation/adequacy_decision.json')
    contexts = decision['model_contexts']
    if set(contexts) != {'enhancer', 'h3k27me3'} or any(x['status'] not in ('PASS', 'FAIL', 'INCONCLUSIVE') for x in contexts.values()):
        raise RuntimeError('Incomplete or invalid model-context adequacy decision')
    passed = decision['both_C_contexts_pass'] is True
    if not isinstance(decision['both_C_contexts_pass'], bool) or passed != all(x['status'] == 'PASS' for x in contexts.values()):
        raise RuntimeError('Inconsistent conjunctive C adequacy decision')
    if not passed and ((stage / 'results/chr7_calibration').exists() or (provenance / 'C_release_precalibration_freeze.json').exists()):
        raise RuntimeError('Calibration/release exists after failed or inconclusive adequacy')
    required_qc = ['input_integrity_verified_before', 'input_integrity_verified_after',
                   'prepared_input_independent_validation', 'phaseI_cache_independent_validation',
                   'runs_final_validation', 'predictions_final_validation',
                   'evaluation_final_validation', 'preservation_final_validation',
                   'checkpoint_packaging_final_validation']
    if passed:
        required_qc.append('calibration_final_validation')
    qc = []
    phases = dict(zip(required_qc, ['inputs', 'inputs', 'prepared', 'cache', 'runs',
                                   'predictions', 'evaluation', 'preservation', 'packaging'] + (['calibration'] if passed else [])))
    for name in required_qc:
        path = provenance / (name + '.json')
        item = load(path)
        if (item.get('status') != 'PASS' or item.get('failures') != 0 or
                item.get('phase') != phases[name] or item.get('validator') != 'independent-execution-audit' or item.get('checks', 0) <= 0):
            raise RuntimeError('Independent final QC did not pass: ' + name)
        table = repo / item['checks_path']
        if sha(table) != item['checks_sha256']:
            raise RuntimeError('Independent QC table changed: ' + name)
        with table.open(newline='') as f:
            check_rows = list(csv.DictReader(f, delimiter='\t'))
        if len(check_rows) != item['checks'] or any(row['status'] != 'PASS' for row in check_rows):
            raise RuntimeError('Independent QC count/status mismatch: ' + name)
        qc.append({'path': str(path.relative_to(stage)), 'sha256': sha(path), 'checks': item['checks']})
    prepared = load(stage / 'inputs/input_manifest.json')
    verify_refs(prepared['files'], stage)
    verify_refs(prepared['frozen_sources'], repo)
    verify_refs(load(stage / 'cache/cache_manifest.json')['files'], stage)
    for rec in load(provenance / 'execution_input_gate.json')['files']:
        path = stage / rec['path']
        if path.stat().st_size != rec['bytes'] or sha(path) != rec['sha256']:
            raise RuntimeError('Prospective execution input/code changed: ' + rec['path'])
    checkpoints = load(provenance / 'checkpoint_freeze.json')
    if checkpoints.get('status') != 'PASS' or len(checkpoints['checkpoints']) != 18:
        raise RuntimeError('Incomplete eighteen-fit checkpoint gate')
    for rec in checkpoints['checkpoints']:
        path = stage / rec['path']
        if sha(path) != rec['sha256']:
            raise RuntimeError('Selected checkpoint changed')
    invariance = load(stage / 'predictions/real_network_invariance.json')
    if invariance.get('status') != 'PASS' or invariance.get('n_seed_audits') != 18 or invariance.get('n_ensemble_audits') != 6:
        raise RuntimeError('Incomplete actual-network invariance gate')
    verify_refs(invariance['files'], stage)
    evaluation = load(stage / 'results/chr7_evaluation/evaluation_manifest.json')
    verify_refs(evaluation['artifacts'], repo)
    verify_refs([{'path': evaluation['input_manifest'], 'sha256': evaluation['input_manifest_sha256']}], repo)
    if passed:
        calibration = load(stage / 'results/chr7_calibration/calibration_manifest.json')
        verify_refs(calibration['artifacts'], repo)
        verify_refs([{'path': calibration['input_manifest'], 'sha256': calibration['input_manifest_sha256']}], repo)
        thresholds = load(stage / 'results/chr7_calibration/C_region_thresholds.json')
        if thresholds.get('status') != 'PASS' or set(thresholds['models']) != {'enhancer', 'h3k27me3'}:
            raise RuntimeError('Incomplete conditional threshold calibration')
    required = ['INTERNAL_TRAINING_REPORT.md', 'results/training_outcomes.tsv',
                'provenance/slurm_accounting.tsv', 'provenance/failure_retry_ledger.tsv',
                'provenance/shared_register_updates.json']
    if passed:
        required += ['provenance/C_release_precalibration_freeze.json',
                     'results/chr7_calibration/calibration_manifest.json']
    for name in required:
        if not (stage / name).is_file():
            raise RuntimeError('Required final artifact missing: ' + name)
    registers = load(provenance / 'shared_register_updates.json')
    for rec in registers['files']:
        path = repo / rec['path']
        if path.stat().st_size != rec['bytes'] or sha(path) != rec['sha256']:
            raise RuntimeError('Final shared-register hash mismatch')
    final_qc = {'status': 'PASS', 'created_utc': datetime.now(timezone.utc).isoformat(),
                'independent_validation_reports': qc, 'total_checks_across_reports': sum(r['checks'] for r in qc),
                'all_18_fits_completed': True, 'real_network_RC_invariance': 'PASS',
                'both_C_contexts_pass': passed, 'calibration_completed': passed,
                'scientific_decision_is_separate_from_execution_QC': True,
                'chr8_chr9_model_execution': False, 'external_benchmark_evaluation': False,
                'candidate_scoring': False, 'pretraining_1_0_and_1_1_unchanged': True,
                'tracked_changes_limited_to_append_only_V2_registers': True}
    write(provenance / 'final_validation.json', final_qc)
    write(provenance / 'finalization_command.json', {'argv': sys.argv, 'cwd': str(Path.cwd()),
          'invoked_utc': datetime.now(timezone.utc).isoformat(), 'script_sha256': sha(Path(__file__)),
          'note': 'Finalizer records invocation before hashing; final freeze is its completion evidence. No model execution.'})
    # No process may still be writing payloads. Exclude only Python bytecode and
    # the checksum/freeze records themselves to avoid a self-referential digest.
    rows = []
    for path in sorted(stage.rglob('*')):
        if not path.is_file() or '__pycache__' in path.parts or path.suffix == '.pyc':
            continue
        rows.append({'path': str(path.relative_to(stage)), 'bytes': path.stat().st_size, 'sha256': sha(path)})
    ledger = provenance / 'artifact_checksums.tsv'
    with ledger.open('x', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=['path', 'bytes', 'sha256'], delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    freeze = {'module': 'COPD-V2-INTERNAL-TRAINING', 'version': 'internal-training-1.0',
              'status': 'COMPLETE_INTERNAL_STAGE_STOP_FOR_INVESTIGATOR_REVIEW',
              'frozen_utc': datetime.now(timezone.utc).isoformat(),
              'payload_files': len(rows), 'payload_bytes': sum(r['bytes'] for r in rows),
              'artifact_checksums': {'path': str(ledger.relative_to(stage)), 'sha256': sha(ledger)},
              'pretraining_version': 'pretraining-1.1', 'publication_baseline_commit': '7c0caf4af5f8f346089c92f43b88161331ca420b',
              'adequacy_decision': {'path': 'results/chr7_evaluation/adequacy_decision.json', 'sha256': sha(stage / 'results/chr7_evaluation/adequacy_decision.json')},
              'both_C_contexts_pass': passed, 'threshold_calibration_completed': passed,
              'ready_for_investigator_review_before_test_opening': passed,
              'chr8_chr9_authorized_or_evaluated': False, 'external_benchmark_evaluated': False,
              'candidate_or_allele_scoring_performed': False, 'automatic_A_B_fallback': False,
              'commit_or_push_performed': False,
              'scope_limit': 'Region-label modeling only; no causal variant or REF-ALT validity claim.'}
    write(provenance / 'freeze.json', freeze)
    print(json.dumps(freeze, indent=2))


if __name__ == '__main__':
    main()
