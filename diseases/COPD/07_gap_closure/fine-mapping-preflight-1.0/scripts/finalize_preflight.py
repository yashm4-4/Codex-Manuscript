#!/usr/bin/env python3
"""Seal or read-only verify the completed preflight; no scientific operations.

The checksum ledger binds every stage file except itself and freeze.json.
The freeze binds that ledger, independent validation and shared register bytes.
Interpreter __pycache__/*.pyc products are explicitly not scientific payloads.
Neither sealing nor verification changes an earlier-stage artifact or Git state.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess

STAGE = Path(__file__).resolve().parents[1]
ROOT = STAGE.parents[3]
LEDGER = STAGE / 'provenance/artifact_checksums.tsv'
FREEZE = STAGE / 'provenance/freeze.json'


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def readj(path):
    return json.loads(path.read_text())


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT,
                                   env=dict(os.environ, GIT_OPTIONAL_LOCKS='0'))


def payload():
    result = []
    for path in sorted(STAGE.rglob('*')):
        if path.is_symlink():
            raise RuntimeError('Symlink not allowed in frozen payload: '+str(path))
        if not path.is_file() or path in {LEDGER, FREEZE}:
            continue
        if '__pycache__' in path.parts or path.suffix == '.pyc':
            continue
        result.append({'path': str(path.relative_to(STAGE)),
                       'bytes': path.stat().st_size, 'sha256': sha(path)})
    return result


def prerequisites():
    init = readj(STAGE/'provenance/initialization.json')
    baseline = init['baseline_commit']
    assert git('rev-parse', 'HEAD').decode().strip() == baseline, 'HEAD changed'
    assert git('ls-tree', '-r', baseline) == (STAGE/'provenance/baseline_git_tree.tsv').read_bytes(), 'Baseline tree mismatch'
    register_receipt = readj(STAGE/'provenance/shared_register_updates.json')
    registers = register_receipt['registers']
    allowed = {r['path'] for r in registers}
    changed = git('diff', '--name-only', '--no-ext-diff', baseline, '--').decode().splitlines()
    assert set(changed) == allowed, ('Unexpected historical changes', changed)
    assert not git('diff', '--cached', '--name-only').strip(), 'Unexpected staged files'
    for row in registers:
        path = ROOT/row['path']
        current = path.read_bytes()
        old = git('show', baseline+':'+row['path'])
        assert current.startswith(old), 'Register prefix changed'
        assert len(old) == row['before_bytes'] and hashlib.sha256(old).hexdigest() == row['before_sha256']
        assert len(current) == row['bytes'] and sha(path) == row['sha256'], 'Register changed since validation'
    validation = readj(STAGE/'provenance/independent_validation.json')
    assert validation['status'] == 'PASS' and validation['failed_checks'] == 0, 'Independent validation not PASS'
    for row in validation['bound_artifacts']:
        path = STAGE/row['path']
        assert path.stat().st_size == row['bytes'] and sha(path) == row['sha256'], ('Validated artifact changed', row['path'])
    for row in init['frozen_inputs']:
        path = ROOT/row['path']
        assert path.stat().st_size == row['bytes'] and sha(path) == row['sha256'], 'Historical input changed'
    return init, registers, validation


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--verify', action='store_true', help='Verify only; create or modify nothing')
    args = parser.parse_args()
    init, registers, validation = prerequisites()
    records = payload()
    if args.verify:
        freeze = readj(FREEZE)
        with LEDGER.open() as handle:
            saved = list(csv.DictReader(handle, delimiter='\t'))
        saved = [dict(r, bytes=int(r['bytes'])) for r in saved]
        assert records == saved, 'Payload membership/size/hash mismatch'
        assert sha(LEDGER) == freeze['artifact_checksums_sha256'], 'Ledger mismatch'
        assert sha(STAGE/'provenance/independent_validation.json') == freeze['independent_validation_sha256']
        assert freeze['shared_registers'] == registers
        assert freeze['payload_files_bound_by_ledger'] == len(records)
        assert freeze['payload_bytes_bound_by_ledger'] == sum(r['bytes'] for r in records)
        assert freeze['baseline_commit'] == init['baseline_commit']
        print(json.dumps({'status':'PASS', 'read_only_verification':True,
                          'payload_files':len(records),
                          'payload_bytes':sum(r['bytes'] for r in records),
                          'artifact_checksums_sha256':sha(LEDGER),
                          'freeze_sha256':sha(FREEZE),
                          'historical_git_content_unchanged_except_append_only_registers':True}, indent=2))
        return
    if FREEZE.exists() or LEDGER.exists():
        raise RuntimeError('Seal already exists or is partial; no overwriting permitted')
    with LEDGER.open('x') as handle:
        writer = csv.DictWriter(handle, ['path','bytes','sha256'], delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(records)
    freeze = {
        'stage':'fine-mapping-preflight-1.0',
        'frozen_utc':datetime.now(timezone.utc).isoformat(),
        'status':'FROZEN_PREFLIGHT_STOP_FOR_INVESTIGATOR_REVIEW',
        'baseline_commit':init['baseline_commit'],
        'payload_files_bound_by_ledger':len(records),
        'payload_bytes_bound_by_ledger':sum(r['bytes'] for r in records),
        'total_stage_payload_files_including_ledger_and_freeze':len(records)+2,
        'artifact_checksums_path':'provenance/artifact_checksums.tsv',
        'artifact_checksums_sha256':sha(LEDGER),
        'checksum_closure':'All stage payload files bound by artifact_checksums.tsv except ledger and freeze itself; freeze binds ledger. External shared registers bound separately.',
        'excluded_ephemeral_files':'Interpreter __pycache__ directories and .pyc files only; no acquired source or scientific artifact excluded.',
        'independent_validation_path':'provenance/independent_validation.json',
        'independent_validation_sha256':sha(STAGE/'provenance/independent_validation.json'),
        'independent_validation_checks':validation['checks'],
        'independent_validation_status':validation['status'],
        'shared_registers':registers,
        'historical_preservation_scope':'Baseline Git tree and current tracked content comparison, explicit frozen-input SHA256 and original-register prefix verification; no historical scientific data parsed.',
        'source_receipts':validation['source_receipts'],
        'saved_source_bytes':validation['saved_source_bytes'],
        'fine_mapping_executed':False,
        'model_or_candidate_scoring_executed':False,
        'external_benchmark_evaluation_executed':False,
        'manuscript_revised':False,
        'commit_or_push_performed':False,
        'future_execution_authorized':False,
    }
    with FREEZE.open('x') as handle:
        json.dump(freeze, handle, indent=2)
        handle.write('\n')
    print(json.dumps(freeze, indent=2))


if __name__ == '__main__':
    main()
