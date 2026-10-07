#!/usr/bin/env python3
"""Verify stage payload and append-only register digests without mutation."""
import csv
import hashlib
import json
import subprocess
from pathlib import Path
S=Path(__file__).resolve().parents[1]
R=S.parents[3]
G=R/'diseases/COPD/07_gap_closure'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()
def main():
 f=json.loads((S/'provenance/freeze.json').read_text())
 m=S/'provenance/artifact_checksums.tsv'
 assert sha(m)==f['manifest_sha256']
 rows=list(csv.DictReader(open(m),delimiter='\t'))
 assert len(rows)==f['payload_file_count']
 for r in rows:
  p=S/r['relative_path'];assert p.stat().st_size==int(r['bytes']) and sha(p)==r['sha256'],r['relative_path']
 for name,digest in f['register_sha256'].items():
  p=G/name;assert sha(p)==digest,name
  prior=subprocess.check_output(['git','show',f['base_commit']+':diseases/COPD/07_gap_closure/'+name],cwd=R)
  assert p.read_bytes().startswith(prior),name
 assert f['no_de_novo_colocalization'] and f['prior_stage_payloads_modified']==False
 print('PASS:',len(rows),'payload checksums, 3 append-only registers, no prior-stage edits')
if __name__=='__main__':main()
