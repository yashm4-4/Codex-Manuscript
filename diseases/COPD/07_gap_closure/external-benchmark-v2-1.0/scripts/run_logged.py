#!/usr/bin/env python3
"""Exclusive-create command provenance for this evaluation stage only."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

STAGE=Path(__file__).resolve().parents[1]

def now():return datetime.now(timezone.utc).isoformat()
def rec(path):
    data=path.read_bytes()
    return {'path':str(path.relative_to(STAGE)),'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
def put(path,data):
    with path.open('x') as h:json.dump(data,h,indent=2,sort_keys=True);h.write('\n')
def main():
    p=argparse.ArgumentParser();p.add_argument('--label',required=True);p.add_argument('command',nargs=argparse.REMAINDER);a=p.parse_args()
    assert a.label.replace('_','').replace('-','').isalnum()
    command=a.command[1:] if a.command and a.command[0]=='--' else a.command
    assert command
    directory=STAGE/'provenance/commands';directory.mkdir(parents=True,exist_ok=True)
    base=directory/a.label
    put(base.with_suffix('.started.json'),{'started_utc':now(),'argv':command,'cwd':os.getcwd(),
        'environment':{k:os.environ.get(k) for k in ['SLURM_JOB_ID','CUDA_VISIBLE_DEVICES','SLURM_JOB_GPUS','PYTHONHASHSEED']},
        'launcher':rec(Path(__file__).resolve())})
    begin=time.monotonic()
    stdout=base.with_suffix('.stdout.log');stderr=base.with_suffix('.stderr.log')
    with stdout.open('xb') as out,stderr.open('xb') as err:
        result=subprocess.run(command,stdout=out,stderr=err,check=False)
    put(base.with_suffix('.completed.json'),{'completed_utc':now(),'argv':command,'returncode':result.returncode,
        'wall_seconds':time.monotonic()-begin,'stdout':rec(stdout),'stderr':rec(stderr)})
    print(json.dumps({'label':a.label,'returncode':result.returncode,'wall_seconds':time.monotonic()-begin}))
    raise SystemExit(result.returncode)
if __name__=='__main__':main()
