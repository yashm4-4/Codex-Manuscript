#!/usr/bin/env python3
"""Exclusive exact-command ledger; streams subprocess output without alteration."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time
from datetime import datetime, timezone


def now():
    return datetime.now(timezone.utc).isoformat()


def write(path, value):
    with path.open('x') as f:
        json.dump(value, f, indent=2, sort_keys=True)
        f.write('\n')


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--stage', type=Path, required=True)
    parser.add_argument('--label', required=True)
    parser.add_argument('command', nargs=argparse.REMAINDER)
    args = parser.parse_args()
    command = args.command[1:] if args.command[:1] == ['--'] else args.command
    if not command or '/' in args.label or '..' in args.label:
        raise ValueError('Invalid command/ledger label')
    ledger = args.stage.resolve() / 'provenance/commands'
    ledger.mkdir(parents=True, exist_ok=True)
    record = {'label': args.label, 'argv': command, 'cwd': os.getcwd(), 'started_utc': now(),
              'launcher': str(Path(__file__).resolve()),
              'launcher_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              'environment': {k: os.environ.get(k) for k in ('PYTHONHASHSEED', 'CUDA_VISIBLE_DEVICES', 'SLURM_JOB_ID', 'SLURM_ARRAY_TASK_ID')}}
    write(ledger / (args.label + '.started.json'), record)
    start = time.monotonic()
    with (ledger / (args.label + '.log')).open('x') as log:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in process.stdout:
            log.write(line)
            log.flush()
            print(line, end='', flush=True)
        code = process.wait()
    write(ledger / (args.label + '.completed.json'), {**record, 'returncode': code, 'completed_utc': now(), 'elapsed_seconds': time.monotonic() - start})
    sys.exit(code)


if __name__ == '__main__':
    main()
