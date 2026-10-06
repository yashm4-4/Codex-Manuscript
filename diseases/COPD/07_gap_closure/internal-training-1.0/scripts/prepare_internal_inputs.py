#!/usr/bin/env python3
"""Derive allowed training/chr7 inputs without changing frozen labels or roles."""
import argparse
import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

CONFIGS = ('V2-A', 'V2-B', 'V2-C')
MODELS = ('enhancer', 'h3k27me3')
TRAIN_CHROMS = ['chr' + str(i) for i in range(1, 23) if i not in (7, 8, 9)] + ['chrX', 'chrY']
ALLOWED = set(TRAIN_CHROMS + ['chr7'])
CHROM_ORDER = {c: i for i, c in enumerate(['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY'])}
COLS = ['interval_id', 'chrom', 'core_start', 'core_end', 'input_start', 'input_end',
        'partition', 'validation_role', 'component_id', 'sequence_available', 'sequence_length',
        'sequence_sha256', 'canonical_rc_sequence_sha256', 'gc_fraction', 'non_acgt_fraction']
NORMALIZE = bytes(x if x in b'ACGT' else ord('N') for x in range(256))
COMPLEMENT = bytes.maketrans(b'ACGTN', b'TGCAN')


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(8 << 20), b''):
            h.update(block)
    return h.hexdigest()


def record(path, base):
    return {'path': str(path.relative_to(base)), 'bytes': path.stat().st_size, 'sha256': sha(path)}


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as f:
        json.dump(value, f, indent=2, sort_keys=True)
        f.write('\n')


def read(path, **kwargs):
    return pd.read_csv(path, sep='\t', dtype=str, keep_default_na=False, **kwargs)


def write_table(path, frame):
    if path.exists():
        raise RuntimeError(f'Refusing to overwrite {path}')
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep='\t', index=False, compression={'method': 'gzip', 'mtime': 0})


class IndexedFasta:
    """Zero-based half-open fetch using the already frozen samtools FAI."""
    def __init__(self, path):
        self.handle = path.open('rb')
        self.index = {}
        for line in Path(str(path) + '.fai').read_text().splitlines():
            name, length, offset, bases, width = line.split('\t')[:5]
            self.index[name] = tuple(map(int, (length, offset, bases, width)))

    def fetch(self, chrom, start, end):
        if chrom not in ALLOWED:
            raise ValueError('Forbidden chromosome sequence access: ' + chrom)
        key = chrom if chrom in self.index else chrom.removeprefix('chr')
        length, offset, bases, width = self.index[key]
        if not 0 <= start < end <= length:
            raise ValueError('Sequence bounds failure')
        first = offset + (start // bases) * width + start % bases
        last = offset + ((end - 1) // bases) * width + (end - 1) % bases
        self.handle.seek(first)
        return self.handle.read(last - first + 1).replace(b'\n', b'').replace(b'\r', b'').upper()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', type=Path, required=True)
    parser.add_argument('--repo', type=Path, required=True)
    args = parser.parse_args()
    stage, repo = args.stage.resolve(), args.repo.resolve()
    frozen = stage.parent / 'pretraining-1.1'
    integrity = json.loads((stage / 'provenance/input_integrity_verified_before.json').read_text())
    # The independent report must contain an explicit PASS; no model call occurs here.
    if integrity.get('status') != 'PASS':
        raise RuntimeError('Independent immutable-input gate is not PASS')
    if (stage / 'inputs/input_manifest.json').exists():
        raise RuntimeError('Inputs already prepared; immutable output cannot be replaced')
    started = time.monotonic()
    paths, frames = [], []
    for model in MODELS:
        for config in CONFIGS:
            path = frozen / f'data/configurations/{config}_{model}_interval_manifest.tsv.gz'
            frame = read(path, usecols=COLS)
            allowed = frame[frame.partition.isin(['train', 'validation'])].copy()
            if not allowed.chrom.isin(ALLOWED).all():
                raise RuntimeError('Forbidden chromosome in permitted partition')
            if not ((allowed.partition.eq('train') & allowed.chrom.isin(TRAIN_CHROMS) & allowed.validation_role.eq('train')) |
                    (allowed.partition.eq('validation') & allowed.chrom.eq('chr7') & allowed.validation_role.isin(['checkpoint', 'selection', 'calibration']))).all():
                raise RuntimeError('Chromosome/role mismatch')
            frames.append(allowed)
            paths.append(path)
    joined = pd.concat(frames, ignore_index=True)
    # All repeated interval IDs must have identical frozen geometry/sequence/role metadata.
    unique = joined.drop_duplicates()
    if unique.interval_id.duplicated().any():
        raise RuntimeError('Inconsistent repeated frozen interval')
    unique = unique.assign(_chrom=unique.chrom.map(CHROM_ORDER), _start=unique.core_start.astype(int))
    unique = unique.sort_values(['_chrom', '_start', 'interval_id'], kind='stable').drop(columns=['_chrom', '_start']).reset_index(drop=True)
    del joined, frames
    n_encoded = unique.canonical_rc_sequence_sha256.nunique()
    seq_path = stage / 'inputs/canonical_sequences.npy'
    seq_path.parent.mkdir(parents=True, exist_ok=True)
    if seq_path.exists():
        raise RuntimeError('Refusing to overwrite sequence array')
    seqs = np.lib.format.open_memmap(seq_path, mode='w+', dtype=np.uint8, shape=(n_encoded, 2001))
    fasta = IndexedFasta(repo / 'diseases/COPD/04_modeling/trednet/fasta/hg38.fa')
    encoded_map, cache_rows, forward_flags = {}, [], []
    for i, row in enumerate(unique.itertuples(index=False)):
        start, end, core_start, core_end = map(int, (row.input_start, row.input_end, row.core_start, row.core_end))
        if not (end - start == 2001 and core_end - core_start == 1000 and start == core_start - 501 and end == core_end + 500 and row.sequence_available == '1' and row.sequence_length == '2001'):
            raise RuntimeError('Frozen geometry failure: ' + row.interval_id)
        raw = fasta.fetch(row.chrom, start, end)
        normalized = raw.translate(NORMALIZE)
        reverse = normalized.translate(COMPLEMENT)[::-1]
        canonical = min(normalized, reverse)
        if (len(raw) != 2001 or hashlib.sha256(raw).hexdigest() != row.sequence_sha256 or
                hashlib.sha256(canonical).hexdigest() != row.canonical_rc_sequence_sha256 or
                abs((raw.count(b'G') + raw.count(b'C')) / 2001 - float(row.gc_fraction)) > 1e-9 or
                abs(normalized.count(b'N') / 2001 - float(row.non_acgt_fraction)) > 1e-9):
            raise RuntimeError('Frozen sequence/hash/covariate failure: ' + row.interval_id)
        key = row.canonical_rc_sequence_sha256
        if key not in encoded_map:
            encoded_map[key] = len(encoded_map)
            seqs[encoded_map[key]] = np.frombuffer(canonical, dtype=np.uint8)
        elif not np.array_equal(seqs[encoded_map[key]], np.frombuffer(canonical, dtype=np.uint8)):
            raise RuntimeError('Encoded hash collision')
        cache_rows.append(encoded_map[key])
        forward_flags.append(int(normalized == canonical))
        if (i + 1) % 100000 == 0:
            print(json.dumps({'sequence_intervals_checked': i + 1, 'encoded_identities': len(encoded_map)}), flush=True)
    seqs.flush()
    fasta.handle.close()
    if len(encoded_map) != n_encoded:
        raise RuntimeError('Encoded count mismatch')
    unique['cache_row'] = cache_rows
    unique['forward_is_canonical'] = forward_flags
    write_table(stage / 'inputs/cache_index.tsv.gz', unique)
    mapping = unique.set_index('interval_id')[['cache_row', 'forward_is_canonical']]
    counts, validation = [], {model: [] for model in MODELS}
    for model in MODELS:
        for config in CONFIGS:
            path = frozen / f'data/configurations/{config}_{model}_interval_manifest.tsv.gz'
            frame = read(path)
            frame = frame[frame.partition.isin(['train', 'validation'])].copy()
            frame = frame.join(mapping, on='interval_id', validate='many_to_one')
            if frame.cache_row.isna().any():
                raise RuntimeError('Missing feature mapping')
            frame['cache_row'] = frame.cache_row.astype(int)
            frame['forward_is_canonical'] = frame.forward_is_canonical.astype(int)
            for role in ('train', 'checkpoint'):
                sub = frame[frame.validation_role.eq(role)]
                out = stage / f'inputs/configurations/{config}_{model}_{role}.tsv.gz'
                write_table(out, sub)
                counts.append({'configuration': config, 'model': model, 'role': role, 'rows': len(sub),
                               'positive': int(sub.label.eq('1').sum()), 'negative': int(sub.label.eq('0').sum())})
            validation[model].append(frame[frame.chrom.eq('chr7')].drop(columns=['configuration', 'model', 'label', 'matched_positive_id', 'matched_control_id', 'matched_final_normalized_chebyshev_distance']))
            if config == 'V2-C':
                write_table(stage / f'inputs/common/{model}_train_C_positive_covariates.tsv.gz', frame[frame.validation_role.eq('train') & frame.label.eq('1')])
        native_union = pd.concat(validation[model], ignore_index=True).drop_duplicates()
        if native_union.interval_id.duplicated().any():
            raise RuntimeError('Conflicting validation metadata')
        native_union = native_union.assign(_start=native_union.core_start.astype(int)).sort_values(['_start', 'interval_id'], kind='stable').drop(columns='_start')
        write_table(stage / f'inputs/validation/{model}_chr7.tsv.gz', native_union)
        path = frozen / f'data/evaluation/{model}_common_challenge_panel.tsv.gz'
        paths.append(path)
        panel = read(path)
        for role in ('selection', 'calibration'):
            sub = panel[panel.validation_role.eq(role)].copy()
            if not sub.chrom.eq('chr7').all() or not set(sub.interval_id) <= set(native_union.interval_id):
                raise RuntimeError('Common-panel role/coverage firewall failure')
            write_table(stage / f'inputs/common/{model}_{role}.tsv.gz', sub)
    outputs = [record(p, stage) for p in sorted((stage / 'inputs').rglob('*')) if p.is_file()]
    manifest = {'status': 'PASS', 'completed': True, 'created_utc': datetime.now(timezone.utc).isoformat(),
                'sequence_intervals': len(unique), 'encoded_identity_count': n_encoded,
                'sequence_shape': [n_encoded, 2001], 'sequence_dtype': 'uint8_ASCII',
                'normalization': 'nonACGT to N; canonical=min(normalized_forward,nucleotide_RC)',
                'chromosome_counts': unique.chrom.value_counts().sort_index().to_dict(),
                'forbidden_sequence_fetches': 0, 'missing_sequences': 0, 'hash_failures': 0,
                'class_role_counts': counts, 'files': outputs,
                'frozen_sources': [record(p, repo) for p in paths], 'elapsed_seconds': time.monotonic() - started}
    write_json(stage / 'inputs/input_manifest.json', manifest)
    print(json.dumps({k: manifest[k] for k in ('status', 'sequence_intervals', 'encoded_identity_count', 'elapsed_seconds')}), flush=True)


if __name__ == '__main__':
    main()
