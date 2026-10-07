#!/usr/bin/env python3
"""Verify historical payloads while allowing documented register appends."""
import csv
import hashlib
import json
import subprocess
from pathlib import Path
S=Path(__file__).resolve().parents[1]
R=S.parents[3]
G=R/'diseases/COPD/07_gap_closure'
P=G/'target-gene-evidence-preflight-1.0'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def main():
 rows=list(csv.DictReader(open(P/'provenance/artifact_checksums.tsv'),delimiter='\t'))
 for r in rows:
  p=P/r['relative_path'];assert sha(p)==r['sha256'],r['relative_path']
 changed=subprocess.check_output(['git','diff','--name-only','HEAD'],cwd=R,text=True).splitlines()
 allowed={f'diseases/COPD/07_gap_closure/{n}' for n in ['activity_log.tsv','gap_closure_decision_register.tsv','gap_closure_result_register.tsv']}
 assert set(changed)<=allowed,changed
 for name in [p.split('/')[-1] for p in allowed]:
  current=(G/name).read_bytes()
  old=subprocess.check_output(['git','show','HEAD:diseases/COPD/07_gap_closure/'+name],cwd=R)
  assert current.startswith(old),name
 result=dict(status='PASS',preflight_payload_checksums=len(rows),tracked_changes=changed,
  prior_register_bytes_are_exact_prefix=True,
  historical_preflight_verifier_note='Its strict historical register digests intentionally fail after authorized append-only entries; prior payload hashes pass.')
 (S/'prior_freeze_integrity_audit.json').write_text(json.dumps(result,indent=2)+'\n')
 print('PASS:',len(rows),'preflight payloads; only append-only register changes')
if __name__=='__main__':main()
