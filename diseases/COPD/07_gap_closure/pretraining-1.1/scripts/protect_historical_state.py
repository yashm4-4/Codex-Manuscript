#!/usr/bin/env python3
"""Opaque integrity snapshots only; never parse historical benchmark outcomes."""
import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
VERSION = ROOT/'diseases/COPD/07_gap_closure/pretraining-1.1'


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    out = VERSION/'provenance/historical_tracked_state.tsv'
    meta = VERSION/'provenance/historical_state_manifest.json'
    paths = subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT).decode().strip('\0').split('\0')
    if out.exists():
        with out.open() as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        assert set(paths) == {row['path'] for row in rows}, 'Tracked file scope changed'
        for row in rows:
            assert digest(ROOT/row['path']) == row['sha256'], row['path']
        assert not subprocess.check_output(['git', 'diff', '--name-only', 'HEAD'], cwd=ROOT).strip()
        print(f'PASS: all {len(rows)} historical tracked working-tree byte hashes unchanged')
        return
    assert not subprocess.check_output(['git', 'diff', '--name-only', 'HEAD'], cwd=ROOT).strip()
    rows = [{'path': path, 'bytes': (ROOT/path).stat().st_size, 'sha256': digest(ROOT/path)} for path in paths]
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['path','bytes','sha256'], delimiter='\t')
        writer.writeheader(); writer.writerows(rows)
    with meta.open('x') as handle:
        json.dump({'recorded_utc': datetime.now(timezone.utc).isoformat(), 'baseline_commit': subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT).decode().strip(), 'tracked_paths': len(rows), 'opaque_byte_hashes_only': True, 'benchmark_outcomes_parsed': False, 'ledger_sha256': digest(out)}, handle, indent=2)
        handle.write('\n')
    print(f'Protected {len(rows)} tracked historical files, including all pretraining-1.0 and benchmark artifacts')


if __name__ == '__main__':
    main()
