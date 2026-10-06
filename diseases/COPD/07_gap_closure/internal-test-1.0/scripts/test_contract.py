"""Test-only adapters around immutable numerical helpers; no training entrypoint."""
from pathlib import Path
import hashlib
import sys
sys.dont_write_bytecode = True

import numpy as np

TRAINING_STAGE = Path(__file__).resolve().parents[2] / 'internal-training-1.0'
if hashlib.sha256((TRAINING_STAGE / 'scripts/phase_two_contract.py').read_bytes()).hexdigest() != 'f50f226b901df6bbd66455462e4da922d8229a4d030f09613a8483a74f636e2f':
    raise RuntimeError('Immutable numerical helper SHA-256 mismatch before import')
sys.path.insert(0, str(TRAINING_STAGE / 'scripts'))
from phase_two_contract import (
    ATOL, RTOL, PHASE_I_SHA256, SEEDS, MODELS, CONFIGURATIONS, FeatureCache,
    array_digest, bool_value, file_record, inference_function, initialize_runtime,
    predict_orientation, read_json, read_rows, resource_record, runtime_environment,
    seed_ensemble, sha256, symmetric_probability, utcnow, verify_record, write_json,
)

EXPECTED_COUNTS = {'enhancer': (25390, 12695, 8754), 'h3k27me3': (3508, 1754, 1809)}


def verify_prospective(stage):
    stage = Path(stage).resolve()
    frozen = read_json(stage / 'provenance/prospective_specification_freeze.json')
    if frozen.get('status') != 'PASS' or frozen.get('before_any_test_model_inference') is not True:
        raise RuntimeError('Prospective specification was not frozen before inference')
    for key in ('specification', 'markdown', 'authorization'):
        verify_record(stage, frozen[key])
    integrity = read_json(stage / 'provenance/input_integrity_verified_before.json')
    if integrity.get('status') != 'PASS' or integrity.get('input_hashes_verified') is not True:
        raise RuntimeError('Source integrity gate is not PASS')
    verify_record(stage, integrity['prospective_specification_freeze'])
    return frozen


def load_test_table(path, context):
    rows = read_rows(path)
    n, positive, components = EXPECTED_COUNTS[context]
    if len(rows) != n or len({r['interval_id'] for r in rows}) != n:
        raise RuntimeError('Frozen test population count/identity mismatch')
    if any(r['chrom'] not in {'chr8', 'chr9'} or r['partition'] != 'test'
           or r['validation_role'] != 'test' or r['configuration'] != 'V2-C'
           or r['model'] != context for r in rows):
        raise RuntimeError('Test chromosome/role firewall breach')
    labels = np.array([int(r['label']) for r in rows], dtype=np.int8)
    if not np.isin(labels, [0, 1]).all() or labels.sum() != positive:
        raise RuntimeError('Frozen test label count mismatch')
    if len({r['component_id'] for r in rows}) != components:
        raise RuntimeError('Frozen test component count mismatch')
    return {'rows': rows, 'interval_id': np.array([r['interval_id'] for r in rows]),
            'cache_row': np.array([int(r['cache_row']) for r in rows], dtype=np.int64),
            'forward_is_canonical': np.array([bool_value(r['forward_is_canonical']) for r in rows])}


def verify_test_checkpoints(stage, repo):
    verify_prospective(stage)
    gate_path = stage / 'provenance/checkpoint_verification.json'
    gate = read_json(gate_path)
    if gate.get('status') != 'PASS' or gate.get('n_selected_original_checkpoints') != 18:
        raise RuntimeError('Original checkpoint verification is incomplete')
    verify_record(stage, gate['prospective_specification_freeze'])
    expected = {(c, m, s) for c in CONFIGURATIONS for m in MODELS for s in SEEDS}
    result = {}
    for rec in gate['checkpoints']:
        key = rec['configuration'], rec['model'], int(rec['seed'])
        path = repo / rec['path']
        if key not in expected or key in result or not path.resolve().is_relative_to(TRAINING_STAGE):
            raise RuntimeError('Unexpected checkpoint identity/source')
        exact = TRAINING_STAGE / f'runs/{key[0]}_{key[1]}_seed{key[2]}/attempt-001/selected_checkpoint.keras'
        if path != exact or rec.get('status') != 'PASS' or rec.get('preserved_original_archive') is not True:
            raise RuntimeError('Checkpoint configuration/model/seed path mapping changed')
        if path.stat().st_size != rec['bytes'] or sha256(path) != rec['sha256']:
            raise RuntimeError('Original checkpoint bytes changed')
        result[key] = path
    if set(result) != expected:
        raise RuntimeError('Missing frozen selected checkpoint')
    return result, gate_path


def external_checkpoint_record(repo, path):
    return {'path': str(path.relative_to(repo)), 'bytes': path.stat().st_size,
            'sha256': sha256(path), 'path_base': 'repository', 'original_archive_unchanged': True}
