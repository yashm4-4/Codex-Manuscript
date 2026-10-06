#!/usr/bin/env python3
"""Hash-bind verified inputs and prospective implementation before GPU stages."""
import argparse
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys

from phase_two_contract import (CONFIGURATIONS, MODELS, file_record, read_json,
                                runtime_environment, sha256, verify_record,
                                write_json, utcnow)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', type=Path, required=True)
    parser.add_argument('--phase', choices=['phase-I', 'training'], required=True)
    args = parser.parse_args()
    stage = args.stage.resolve()
    repo = stage.parents[3]
    integrity = read_json(stage / 'provenance/input_integrity_verified_before.json')
    if integrity.get('status') != 'PASS' or integrity.get('failures') != 0 or integrity.get('checks', 0) < 518:
        raise RuntimeError('Independent input integrity gate failed')
    im = read_json(stage / 'inputs/input_manifest.json')
    if im.get('status') != 'PASS' or im.get('completed') is not True:
        raise RuntimeError('Input preparation gate failed')
    for rec in im['files']:
        verify_record(stage, rec)
    common = ['EXECUTION_SPECIFICATION.md', 'specification/internal_evaluation_implementation.json',
              'scripts/runtime.sh', 'scripts/run_logged.py', 'scripts/phase_two_contract.py',
              'scripts/prepare_internal_inputs.py', 'scripts/prepare_execution_gate.py',
              'inputs/input_manifest.json', 'provenance/input_integrity_verified_before.json',
              'provenance/input_integrity_verified_before.tsv']
    if args.phase == 'phase-I':
        packages = []
        for dist in importlib.metadata.distributions():
            metadata = []
            for name in ('METADATA', 'RECORD', 'WHEEL', 'INSTALLER'):
                path = Path(dist._path) / name
                if path.is_file():
                    metadata.append({'path': str(path.resolve()), 'bytes': path.stat().st_size, 'sha256': sha256(path)})
            packages.append({'name': dist.metadata['Name'], 'version': dist.version, 'metadata': metadata})
        env = runtime_environment()
        env.update({'created_utc': utcnow(), 'interpreter_sha256': sha256(sys.executable),
                    'packages_all': sorted(packages, key=lambda row: (row['name'] or '').lower()),
                    'runtime_script': file_record(stage, stage / 'scripts/runtime.sh'),
                    'python_install_command': ['uv', 'python', 'install', '3.13.0'],
                    'original_venv_python_symlink': os.readlink(repo / 'models/TREDNET_v2/.venv/bin/python'),
                    'original_pyvenv_cfg': (repo / 'models/TREDNET_v2/.venv/pyvenv.cfg').read_text(),
                    'interpretation': 'Match frozen training manifest Python3.13.0; preserve original read-only site-packages; original venv interpreter inaccessible and pyvenv.cfg patch differs.'})
        if (platform.python_version(), *(importlib.metadata.version(k) for k in ('tensorflow', 'keras', 'numpy'))) != ('3.13.0', '2.20.0', '3.14.1', '2.5.0'):
            raise RuntimeError('Runtime version mismatch')
        write_json(stage / 'provenance/executed_environment.json', env)
        common += ['scripts/extract_phase_one.py', 'provenance/executed_environment.json']
        target = 'provenance/phase_I_execution_gate.json'
        extra = {'phase_I_authorized': True, 'sequence_intervals': im['sequence_intervals'],
                 'encoded_identity_count': im['encoded_identity_count'], 'test_extraction_authorized': False}
    else:
        cache = read_json(stage / 'cache/cache_manifest.json')
        if cache.get('status') != 'PASS' or cache.get('completed') is not True or cache.get('failures') != 0:
            raise RuntimeError('Feature extraction gate failed')
        for rec in cache['files']:
            verify_record(stage, rec)
        for config in CONFIGURATIONS:
            for model in MODELS:
                for role in ('train', 'checkpoint'):
                    common.append(f'inputs/configurations/{config}_{model}_{role}.tsv.gz')
        common += [f'inputs/validation/{model}_chr7.tsv.gz' for model in MODELS]
        common += ['cache/cache_manifest.json', 'scripts/train_phase_two.py', 'scripts/predict_chr7.py',
                   'scripts/train_array.sh', 'scripts/freeze_checkpoints.py', 'scripts/evaluate_chr7.py',
                   'scripts/calibrate_chr7.py', 'scripts/prepare_evaluation_gates.py',
                   'provenance/phase_I_execution_gate.json', 'provenance/executed_environment.json']
        target = 'provenance/execution_input_gate.json'
        extra = {'training_authorized': True, 'input_hashes_verified': True,
                 'fits_authorized': 18, 'selection_performance_authorized_during_training': False,
                 'test_inference_authorized': False, 'external_evaluation_authorized': False}
    payload = {'status': 'PASS', 'created_utc': utcnow(), 'phase': args.phase, **extra,
               'authorization_attachment_sha256': sha256('/home/maheshwarany2/.codex/attachments/933794dc-e228-4aa9-a4d9-db8cf11e31e7/pasted-text.txt'),
               'files': [file_record(stage, stage / name) for name in common]}
    write_json(stage / target, payload)
    print(json.dumps({'status': 'PASS', 'gate': target, 'files': len(common)}))


if __name__ == '__main__':
    main()
