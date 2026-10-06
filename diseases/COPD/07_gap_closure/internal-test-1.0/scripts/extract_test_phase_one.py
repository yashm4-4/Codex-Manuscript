#!/usr/bin/env python3
"""Frozen phase-I extraction for the prospective chr8–9 common test panels only."""
import argparse
import json
import os
from pathlib import Path
import time
import traceback

import numpy as np

from test_contract import (PHASE_I_SHA256, SEEDS, verify_prospective, file_record, initialize_runtime,
                                read_json, resource_record, runtime_environment,
                                sha256, verify_record, write_json, utcnow)


def build_phase_one(keras):
    max_norm = keras.constraints.MaxNorm(max_value=0.9, axis=0)
    glorot_uniform = keras.initializers.GlorotUniform()
    return keras.Sequential([
        keras.layers.Input(shape=(2001, 4)),
        keras.layers.Conv1D(320, 8, activation='relu', padding='valid', kernel_constraint=max_norm, kernel_initializer=glorot_uniform, name='conv1d_1'),
        keras.layers.Conv1D(320, 8, activation='relu', padding='valid', kernel_constraint=max_norm, kernel_initializer=glorot_uniform, name='conv1d_2'),
        keras.layers.Dropout(0.2, name='dropout_1'),
        keras.layers.MaxPooling1D(pool_size=6, name='max_pooling1d_1'),
        keras.layers.Conv1D(480, 8, activation='relu', padding='valid', kernel_constraint=max_norm, kernel_initializer=glorot_uniform, name='conv1d_3'),
        keras.layers.Conv1D(480, 8, activation='relu', padding='valid', kernel_constraint=max_norm, kernel_initializer=glorot_uniform, name='conv1d_4'),
        keras.layers.Dropout(0.2, name='dropout_2'),
        keras.layers.MaxPooling1D(pool_size=6, name='max_pooling1d_2'),
        keras.layers.Conv1D(640, 8, activation='relu', padding='valid', kernel_constraint=max_norm, kernel_initializer=glorot_uniform, name='conv1d_5'),
        keras.layers.Conv1D(640, 8, activation='relu', padding='valid', kernel_constraint=max_norm, kernel_initializer=glorot_uniform, name='conv1d_6'),
        keras.layers.Dropout(0.5, name='dropout_3'),
        keras.layers.Flatten(name='flatten_1'),
        keras.layers.Dense(4560, activation='relu', name='dense_1'),
        keras.layers.Dense(4560, activation='linear', name='dense_2'),
        keras.layers.Activation('sigmoid', name='activation_1'),
    ])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', type=Path, required=True)
    parser.add_argument('--repo', type=Path, required=True)
    args = parser.parse_args()
    stage, repo = args.stage.resolve(), args.repo.resolve()
    verify_prospective(stage)
    cache = stage / 'cache'
    cache.mkdir(exist_ok=False)
    started = time.monotonic()
    try:
        gate = read_json(stage / 'provenance/phase_I_execution_gate.json')
        if gate.get('status') != 'PASS' or gate.get('phase_I_authorized') is not True:
            raise RuntimeError('No approved phase-I input/code gate')
        for rec in gate['files']:
            verify_record(stage, rec)
        manifest_path = stage / 'inputs/input_manifest.json'
        manifest = read_json(manifest_path)
        if manifest.get('status') != 'PASS' or manifest.get('completed') is not True:
            raise RuntimeError('Sequence input preparation failed')
        if not set(manifest['chromosome_counts']) <= {'chr8', 'chr9'}:
            raise RuntimeError('Forbidden chromosome present in cache input')
        if manifest['sequence_intervals'] != 26225 or manifest['encoded_identity_count'] != 26225:
            raise RuntimeError('Exact frozen test union geometry mismatch')
        for rec in manifest['files']:
            verify_record(stage, rec)
        weights = repo / 'diseases/COPD/04_modeling/trednet/model_phase_I/phase_one_weights.h5'
        if sha256(weights) != PHASE_I_SHA256:
            raise RuntimeError('Phase-I weight hash mismatch; random fallback prohibited')
        sequences = np.load(stage / 'inputs/canonical_sequences.npy', mmap_mode='r', allow_pickle=False)
        if sequences.dtype != np.uint8 or sequences.shape != (manifest['encoded_identity_count'], 2001):
            raise RuntimeError('Sequence representation geometry failure')
        tf, keras = initialize_runtime(SEEDS[0])
        write_json(cache / 'environment.json', runtime_environment(tf, keras))
        model = build_phase_one(keras)
        model.load_weights(weights)
        model.trainable = False
        if model.output_shape != (None, 4560):
            raise RuntimeError('Phase-I representation geometry failure')
        write_json(cache / 'architecture.json', json.loads(model.to_json()))
        weight_digest_before = [sha256(weights)]
        # All weights loaded through strict native Keras, and no optimizer exists.
        for variable in model.weights:
            if not np.isfinite(variable.numpy()).all():
                raise RuntimeError('Nonfinite frozen model weight')
        infer = tf.function(lambda batch: model(batch, training=False),
                            input_signature=[tf.TensorSpec((None, 2001, 4), tf.float32)], jit_compile=False)
        canonical = np.lib.format.open_memmap(cache / 'features_canonical.npy', mode='w+', dtype=np.float32, shape=(len(sequences), 4560))
        reverse = np.lib.format.open_memmap(cache / 'features_rc.npy', mode='w+', dtype=np.float32, shape=(len(sequences), 4560))
        alphabet = np.frombuffer(b'ACGT', dtype=np.uint8)
        allowed = np.frombuffer(b'ACGTN', dtype=np.uint8)
        batch_size = 32
        write_json(cache / 'started.json', {'started_utc': utcnow(), 'batch_size': batch_size,
                   'shape': [len(sequences), 4560], 'phase_I_weights_sha256': PHASE_I_SHA256,
                   'input_manifest': file_record(stage, manifest_path),
                   'phase_I_gate': file_record(stage, stage / 'provenance/phase_I_execution_gate.json'),
                   'two_independent_nucleotide_orientation_passes': True,
                   'representation_axis_reversal': False, 'optimizer': None})
        for start in range(0, len(sequences), batch_size):
            chars = sequences[start:start + batch_size]
            if not np.isin(chars, allowed).all():
                raise RuntimeError('Non-normalized nucleotide cache input')
            onehot = (chars[:, :, None] == alphabet[None, None, :]).astype(np.float32)
            # This reverses nucleotide positions and A/C/G/T channels before phase I,
            # never the 4560 features. The second network call is independent.
            rc_onehot = np.ascontiguousarray(onehot[:, ::-1, ::-1])
            forward_features = np.asarray(infer(onehot), dtype=np.float32)
            reverse_features = np.asarray(infer(rc_onehot), dtype=np.float32)
            expected = (len(chars), 4560)
            for features in (forward_features, reverse_features):
                if features.shape != expected or not np.isfinite(features).all() or np.any(features < 0) or np.any(features > 1):
                    raise RuntimeError('Frozen representation failure')
            canonical[start:start + len(chars)] = forward_features
            reverse[start:start + len(chars)] = reverse_features
            if (start // batch_size + 1) % 200 == 0:
                canonical.flush()
                reverse.flush()
                print(json.dumps({'encoded_identities_completed': start + len(chars), 'total': len(sequences),
                                  'elapsed_seconds': time.monotonic() - started}), flush=True)
        canonical.flush()
        reverse.flush()
        if sha256(weights) != weight_digest_before[0]:
            raise RuntimeError('Phase-I source weights changed during extraction')
        outputs = [file_record(stage, cache / name) for name in ('features_canonical.npy', 'features_rc.npy', 'environment.json', 'architecture.json', 'started.json')]
        write_json(cache / 'cache_manifest.json', {'completed': True, 'status': 'PASS', 'completed_utc': utcnow(),
                   'sequence_intervals': manifest['sequence_intervals'], 'encoded_identity_count': len(sequences),
                   'shape': [len(sequences), 4560], 'dtype': 'float32', 'orientation_network_evaluations': len(sequences) * 2,
                   'phase_I_weights_sha256': PHASE_I_SHA256, 'files': outputs,
                   'representation': 'final4560sigmoid', 'nucleotide_rc_independently_evaluated': True,
                   'test_sequence_extractions': len(sequences), 'failures': 0, 'missing_sequences': 0, 'scope': 'frozen common C-task chr8-9 panels only',
                   'prospective_freeze': file_record(stage, stage / 'provenance/prospective_specification_freeze.json'),
                   'resources': resource_record(started)})
        print(json.dumps({'status': 'PASS', 'cache_manifest': str(cache / 'cache_manifest.json')}), flush=True)
    except BaseException as exc:
        write_json(cache / 'failure.json', {'status': 'FAIL', 'time_utc': utcnow(), 'error': repr(exc),
                   'traceback': traceback.format_exc(), 'resources': resource_record(started)})
        raise


if __name__ == '__main__':
    main()
