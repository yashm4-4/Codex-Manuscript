#!/usr/bin/env python3
"""Exclusive hash-bound transitions for the frozen test stage; no model calls."""
import argparse
from pathlib import Path
from test_contract import (file_record, read_json, sha256, utcnow, verify_record,
                           verify_prospective, verify_test_checkpoints, write_json)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', type=Path, required=True)
    parser.add_argument('--repo', type=Path, required=True)
    parser.add_argument('--phase', choices=['phase-I', 'inference', 'evaluation'], required=True)
    args = parser.parse_args()
    stage, repo = args.stage.resolve(), args.repo.resolve()
    frozen = verify_prospective(stage)
    for record in frozen['sources']:
        path = repo / record['path']
        if sha256(path) != record['sha256'] or ('bytes' in record and path.stat().st_size != record['bytes']):
            raise RuntimeError('Immutable consumed source changed: ' + record['path'])
    paths = set()
    def bind(path):
        paths.add(path)
        return read_json(stage / path)
    def passed(path):
        value = bind(path)
        if value.get('status') != 'PASS':
            raise RuntimeError('Required independent gate is not PASS: ' + path)
        return value
    for path in ['provenance/prospective_specification_freeze.json',
                 'specification/test_specification.json',
                 'provenance/input_integrity_verified_before.json',
                 'provenance/checkpoint_verification.json']:
        bind(path)
    paths.update(['TEST_SPECIFICATION.md', 'provenance/user_authorization.txt'])
    passed('provenance/inputs_independent_validation.json')
    manifest = passed('inputs/input_manifest.json')
    if manifest.get('completed') is not True or manifest.get('encoded_identity_count') != 26225:
        raise RuntimeError('Incomplete frozen test sequence input manifest')
    for record in manifest['files']:
        verify_record(stage, record)
        paths.add(record['path'])
    for name in ['test_contract.py', 'extract_test_phase_one.py', 'predict_test.py',
                 'evaluate_test.py', 'prepare_test_inputs.py', 'prepare_test_gates.py',
                 'runtime.sh', 'run_logged.py']:
        paths.add('scripts/' + name)
    paths.update(str(p.relative_to(stage)) for p in (stage / 'tests').glob('test_*.py'))
    if args.phase == 'phase-I':
        if (stage / 'cache').exists() or (stage / 'predictions').exists():
            raise RuntimeError('Phase-I gate must precede first test model invocation')
        output = 'provenance/phase_I_execution_gate.json'
        authorization = {'phase_I_authorized': True, 'phase_II_inference_authorized': False}
    else:
        passed('provenance/cache_independent_validation.json')
        cache = passed('cache/cache_manifest.json')
        for record in cache['files']:
            verify_record(stage, record)
            paths.add(record['path'])
        phase_one_gate = passed('provenance/phase_I_execution_gate.json')
        for record in phase_one_gate['files']:
            verify_record(stage, record)
        verify_test_checkpoints(stage, repo)
        if args.phase == 'inference':
            if (stage / 'predictions').exists():
                raise RuntimeError('Inference gate must precede first phase-II test score')
            output = 'provenance/inference_execution_gate.json'
            authorization = {'phase_II_inference_authorized': True}
        else:
            passed('provenance/predictions_independent_validation.json')
            inference_gate = passed('provenance/inference_execution_gate.json')
            for record in inference_gate['files']:
                verify_record(stage, record)
            invariance = passed('predictions/real_network_invariance.json')
            if (not invariance['real_network'] or not invariance['all_configurations_all_seeds_and_ensembles']
                    or invariance['n_seed_audits'] != 18 or invariance['n_ensemble_audits'] != 6):
                raise RuntimeError('Incomplete real-network test RC gate')
            for record in invariance['files']:
                verify_record(stage, record)
                paths.add(record['path'])
            output = 'provenance/evaluation_execution_gate.json'
            authorization = {'one_time_test_evaluation_authorized': True}
    gate = {'status': 'PASS', 'created_utc': utcnow(), 'phase': args.phase,
            'scope': 'Frozen common C-task chr8-9 internal holdout only',
            'files': [file_record(stage, stage / path) for path in sorted(paths)],
            **authorization}
    write_json(stage / output, gate)
    if args.phase == 'evaluation':
        def absolute(path):
            path = Path(path)
            return {'path': str(path), 'bytes': path.stat().st_size, 'sha256': sha256(path)}
        training = stage.parent / 'internal-training-1.0'
        pretraining = stage.parent / 'pretraining-1.1'
        checkpoints = read_json(stage / 'provenance/checkpoint_verification.json')['checkpoints']
        ready = {
            'version': 'internal-test-1.0',
            'specification': absolute(stage / 'specification/test_specification.json'),
            'prospective_freeze': absolute(stage / 'provenance/prospective_specification_freeze.json'),
            'frozen_design': absolute(pretraining / 'specification/selection_adequacy.json'),
            'frozen_cutpoints': absolute(training / 'results/chr7_evaluation/stratification_cutpoints.json'),
            'frozen_thresholds': absolute(training / 'results/chr7_calibration/C_region_thresholds.json'),
            'evaluation_execution_gate': absolute(stage / output),
            'prerequisites': {
                'input_integrity': {'status': 'PASS', 'evidence': absolute(stage / 'provenance/inputs_independent_validation.json')},
                'checkpoint_verification': {'status': 'PASS', 'evidence': absolute(stage / 'provenance/checkpoint_verification.json')},
                'invariance': {'status': 'PASS', 'real_network': True,
                    'all_configurations_all_seeds_and_ensembles': True,
                    'evidence': absolute(stage / 'predictions/real_network_invariance.json')},
                'checkpoints': [dict(absolute(repo / r['path']), configuration=r['configuration'],
                    model=r['model'], seed=r['seed'], completed=True) for r in checkpoints]},
            'models': {context: {'test': absolute(stage / f'inputs/common/{context}_test.tsv.gz'),
                      'predictions': absolute(stage / f'predictions/{context}_chr8_9.tsv.gz')}
                       for context in ('enhancer', 'h3k27me3')}}
        write_json(stage / 'provenance/evaluation_input.json', ready)
    print('PASS: ' + output, flush=True)


if __name__ == '__main__':
    main()
