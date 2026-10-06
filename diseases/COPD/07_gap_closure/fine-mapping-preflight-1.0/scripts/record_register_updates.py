#!/usr/bin/env python3
"""Validate already appended register rows and record exact before/after bytes."""
import csv,datetime,hashlib,io,json,pathlib,subprocess
STAGE=pathlib.Path(__file__).resolve().parents[1];ROOT=STAGE.parents[3]
def sha(b):return hashlib.sha256(b).hexdigest()
def main():
    baseline=json.loads((STAGE/'provenance/initialization.json').read_text());records=[]
    expected={'activity_log.tsv':6,'gap_closure_decision_register.tsv':6,'gap_closure_result_register.tsv':8}
    for old in baseline['baseline_registers']:
        name=pathlib.Path(old['path']).name;before=(STAGE/'provenance'/('baseline_'+name)).read_bytes();now=(ROOT/old['path']).read_bytes()
        assert len(before)==old['bytes'] and sha(before)==old['sha256'] and now.startswith(before)
        assert before==subprocess.check_output(['git','show',baseline['baseline_commit']+':'+old['path']],cwd=ROOT)
        rows=list(csv.reader(io.StringIO(now.decode()),delimiter='\t'));assert all(len(r)==len(rows[0]) for r in rows)
        assert len({r[0] for r in rows[1:]})==len(rows)-1
        added=now[len(before):];assert len(added.splitlines())==expected[name]
        records.append({'path':old['path'],'before_bytes':len(before),'before_sha256':sha(before),'bytes':len(now),'sha256':sha(now),
                        'appended_bytes':len(added),'appended_rows':len(added.splitlines()),'append_only':True})
    result={'stage':'fine-mapping-preflight-1.0','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'baseline_commit':baseline['baseline_commit'],'status':'PASS','prior_bytes_changed':False,'registers':records,
            'updates_applied_with':'apply_patch; this script verifies only'}
    (STAGE/'provenance/shared_register_updates.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
