#!/usr/bin/env python3
"""Lock the benchmark protocol and inventory immutable inputs before assembly."""
import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
BASE = ROOT / 'diseases/COPD/07_gap_closure'
PROV = BASE / 'provenance'
PREFIX = 'COPD-V2-BENCH'


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    lock = PROV / f'{PREFIX}_protocol_lock.json'
    if lock.exists():
        raise SystemExit('Existing protocol lock: refusing to overwrite.')
    started = datetime.now(timezone.utc).isoformat()
    tracked = subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT).decode().split('\0')
    paths = {p for p in tracked if p.startswith('diseases/COPD/') and
             ('/07_gap_closure/' not in p or 'COPD-V2-PHENO' in p or 'COPD-V2-RC' in p)}
    with (PROV / 'COPD-V2-RC_consumed_file_hashes.tsv').open() as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            paths.add(row['path'] if 'path' in row else row['repo_path'])
    rows = []
    for path in sorted(paths):
        target = ROOT / path
        rows.append({'path': path, 'sha256': digest(target), 'bytes': target.stat().st_size,
                     'resolved_path': str(target.resolve())})
    ledger = PROV / f'{PREFIX}_frozen_input_hashes_before.tsv'
    with ledger.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['path', 'sha256', 'bytes', 'resolved_path'], delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    specification = PROV / f'{PREFIX}_analysis_specification.md'
    status = subprocess.check_output(['git', 'status', '--porcelain=v1', '-uall'], cwd=ROOT).decode()
    record = {'started_utc': started, 'protocol_locked_utc': datetime.now(timezone.utc).isoformat(),
              'protocol_path': str(specification.relative_to(ROOT)), 'protocol_sha256': digest(specification),
              'baseline_git_head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT).decode().strip(),
              'baseline_git_status': status, 'frozen_input_count': len(rows),
              'frozen_input_ledger': str(ledger.relative_to(ROOT)), 'frozen_input_ledger_sha256': digest(ledger),
              'new_benchmark_model_comparison_started': False,
              'known_prior_observations_disclosed': ['rs2013701 forward recovery and RC loss',
                                                    'prior V1/RC selected literature examples'],
              'future_external_evaluation_firewall': True}
    lock.write_text(json.dumps(record, indent=2) + '\n')
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
