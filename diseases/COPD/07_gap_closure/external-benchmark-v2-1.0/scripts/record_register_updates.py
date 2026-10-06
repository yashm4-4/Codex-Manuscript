#!/usr/bin/env python3
"""Record and verify already appended stage register updates; never edit registers."""
import csv
from datetime import datetime,timezone
import hashlib
import io
import json
from pathlib import Path
import subprocess

STAGE=Path(__file__).resolve().parents[1]
REPO=STAGE.parents[3]
BASELINE='c9515770db4dc38064d4c1fd7fa323788040ec37'
EXPECTED={'activity_log.tsv':7,'gap_closure_decision_register.tsv':5,'gap_closure_result_register.tsv':8}
def sha(data):return hashlib.sha256(data).hexdigest()
def main():
    original=json.loads((STAGE/'provenance/preopen_source_manifest.json').read_text())
    records=[]
    for old in original['register_baselines']:
        path=REPO/old['path']
        before=subprocess.check_output(['git','-C',str(REPO),'show',BASELINE+':'+old['path']])
        after=path.read_bytes()
        assert len(before)==old['bytes'] and sha(before)==old['sha256']
        assert after.startswith(before)
        suffix=after[len(before):]
        assert suffix.count(b'\n')==EXPECTED[path.name]
        parsed=list(csv.reader(io.StringIO(after.decode()),delimiter='\t'))
        assert all(len(row)==len(parsed[0]) for row in parsed)
        assert len({row[0] for row in parsed})==len(parsed)
        records.append({'path':old['path'],'bytes':len(after),'sha256':sha(after),'append_only':True,
            'before_bytes':len(before),'before_sha256':sha(before),'appended_bytes':len(suffix),'appended_rows':suffix.count(b'\n')})
    data={'status':'PASS','created_utc':datetime.now(timezone.utc).isoformat(),'stage':'external-benchmark-v2-1.0',
        'baseline_commit':BASELINE,'registers':records,'prior_bytes_changed':False,'commit_or_push_performed':False,
        'registers_updated_by':'Root append-only apply_patch; this script performs read-only verification'}
    target=STAGE/'provenance/shared_register_updates.json'
    with target.open('x') as h:json.dump(data,h,indent=2,sort_keys=True);h.write('\n')
    print(json.dumps(data,indent=2))
if __name__=='__main__':main()
