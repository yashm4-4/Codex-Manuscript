#!/usr/bin/env python3
"""Freeze a validated review-only proposal; never authorize/launch training."""
from __future__ import annotations
import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
BASE=ROOT/'diseases/COPD/07_gap_closure'
DATA=BASE/'data/COPD-V2-PREFLIGHT'
PROV=BASE/'provenance'
OUT=BASE/'results'
P='COPD-V2-PREFLIGHT'
SCRIPT_NAMES=['build_model_redesign_preflight.py','construct_preflight_configurations.py','preflight_orientation_contract.py','summarize_preflight_readiness.py','validate_model_redesign_preflight.py','validate_model_redesign_preflight_strict_tie_attempt.py','validate_preflight_ancillary.py','freeze_model_redesign_preflight.py']

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()

def read_json(path):return json.loads(path.read_text())

def table(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))

def write_table(path,rows,fields):
 with path.open('x',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)

def main():
 frozen=PROV/f'{P}_freeze.json'
 ledger=PROV/f'{P}_artifact_checksums.tsv'
 inventory=PROV/f'{P}_artifact_inventory.tsv'
 if any(x.exists() for x in (frozen,ledger,inventory)):
  raise RuntimeError('Preflight already frozen or partial freeze exists; do not overwrite')
 validation=read_json(PROV/f'{P}_final_validation_manifest.json')
 if validation['artifact_validation_status']!='PASS' or validation['artifact_failures']!=0:
  raise RuntimeError('Artifact validation has not passed')
 ancillary=read_json(PROV/f'{P}_ancillary_validation.json')
 if ancillary['status']!='PASS':
  raise RuntimeError('Ancillary independent validation has not passed')
 selection=read_json(PROV/f'{P}_selection_specification.json')
 assert selection['training_authorized'] is False
 assert selection['training_readiness']=='BLOCKED_BY_PRESPECIFIED_MATCHING_GATES'
 assert selection['seeds']==[104729,130363,155921]
 assert validation['training_readiness']=='BLOCKED_NOT_AUTHORIZED'
 review=read_json(PROV/f'{P}_design_review.json')
 assert review['methodology_review_status']=='PASS'
 assert review['dataset_training_readiness']=='NOT_READY'
 for item in review['reviewed_artifacts']:
  if sha(BASE/item['path'])!=item['sha256']:
   raise RuntimeError(f'Reviewed artifact changed after review: {item["path"]}')
 # Review JSON is preserved in full; the independently validated artifact status
 # is separate from scientific readiness, which remains blocked by this snapshot.
 preserved=[]
 for filename in [f'{P}_consumed_input_hashes.tsv',f'{P}_protected_tracked_hashes_before.tsv']:
  source=PROV/filename
  rows=table(source)
  for row in rows:
   path=ROOT/row['path']
   if not path.is_file() or path.stat().st_size!=int(row['bytes']) or sha(path)!=row['sha256']:
    raise RuntimeError(f'Protected/input file changed: {row["path"]}')
  preserved.append({'ledger':str(source.relative_to(ROOT)),'rows_verified':len(rows),'sha256':sha(source)})
 # No pre-existing tracked file may have changed, including every BENCH artifact
 # and the shared historical registers. BENCH bytes are hashed, never parsed.
 subprocess.run(['git','diff','--quiet','HEAD','--'],cwd=ROOT,check=True)
 paths={BASE/'MODEL_REDESIGN_PREFLIGHT.md'}
 paths.update(x for x in DATA.rglob('*') if x.is_file())
 paths.update(OUT.glob(f'{P}_*'))
 paths.update(PROV.glob(f'{P}_*'))
 paths.update(BASE/'scripts'/x for x in SCRIPT_NAMES)
 if any(not p.is_file() or p.is_symlink() for p in paths):raise RuntimeError('Missing or symlinked artifact')
 if any('__pycache__' in str(p) or p.suffix=='.pyc' for p in paths):raise RuntimeError('Runtime cache in payload')
 rows=[]
 for path in sorted(paths):
  rel=str(path.relative_to(ROOT))
  category='data' if '/data/' in rel else 'result' if '/results/' in rel else 'script' if '/scripts/' in rel else 'specification' if path.name=='MODEL_REDESIGN_PREFLIGHT.md' else 'provenance_or_QC'
  rows.append({'path':rel,'category':category,'bytes':path.stat().st_size,'sha256':sha(path)})
 write_table(inventory,rows,['path','category','bytes','sha256'])
 rows.append({'path':str(inventory.relative_to(ROOT)),'category':'artifact_inventory','bytes':inventory.stat().st_size,'sha256':sha(inventory)})
 write_table(ledger,[{k:r[k] for k in ['path','bytes','sha256']} for r in rows],['path','bytes','sha256'])
 request=Path('/home/maheshwarany2/.codex/attachments/d6278c56-b146-4a32-b6db-8cc8f006851d/pasted-text.txt')
 value={'module':P,'version':'pretraining-1.0','frozen_utc':datetime.now(timezone.utc).isoformat(),'git_baseline':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'status':'FROZEN_REVIEW_ONLY_TRAINING_BLOCKED','artifact_validation_status':validation['artifact_validation_status'],'artifact_checks':validation['artifact_checks'],'artifact_failures':0,'matching_balance_failures':len(validation['B_C_balance_failures']),'training_authorized':False,'training_started':False,'GPU_training_or_inference_performed':False,'benchmark_outcomes_parsed':False,'benchmark_version_1_unchanged':True,'protected_ledgers_verified':preserved,'training_seeds':selection['seeds'],'n_planned_fits':18,'payload_count':len(rows),'payload_bytes':sum(x['bytes'] for x in rows),'artifact_checksum_ledger':str(ledger.relative_to(ROOT)),'artifact_checksum_ledger_sha256':sha(ledger),'checksum_self_exclusions':['this freeze JSON','artifact checksum ledger itself'],'request_attachment_sha256':sha(request),'design_report':str((BASE/'MODEL_REDESIGN_PREFLIGHT.md').relative_to(ROOT)),'design_report_sha256':sha(BASE/'MODEL_REDESIGN_PREFLIGHT.md'),'selection_specification_sha256':sha(PROV/f'{P}_selection_specification.json'),'design_review_sha256':sha(PROV/f'{P}_design_review.json'),'interpretation':'Artifacts faithfully preserve an audited candidate experiment whose matching readiness FAILED. A new reviewed construction version is required to resolve matching; this freeze is not authority to train. No experiment has begun.','next_authority_required':'Investigator decision on control/common-support revision, followed by replacement data/rule freeze and separate explicit GPU authorization.','publication':'No commit or push in this module; complete artifacts exist locally including data paths possibly covered by repository ignore rules.'}
 value['ancillary_validation_status']=ancillary['status']
 value['ancillary_checks']=ancillary['n_checks']
 with frozen.open('x') as f:json.dump(value,f,indent=2,sort_keys=True);f.write('\n')
 for row in rows:
  if sha(ROOT/row['path'])!=row['sha256']:raise RuntimeError('Post-freeze payload drift')
 print(json.dumps(value,indent=2))
 print('Freeze SHA256:',sha(frozen))

if __name__=='__main__':main()
