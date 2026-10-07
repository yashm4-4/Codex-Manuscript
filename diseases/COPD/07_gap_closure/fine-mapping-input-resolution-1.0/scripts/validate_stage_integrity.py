"""Read-only scientific-output, historical-content and append-only verification."""
from pathlib import Path
import csv
import hashlib
import json
import subprocess

s=Path(__file__).resolve().parents[1]
repo=s.parents[3]
checks=[]


def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''): h.update(b)
    return h.hexdigest()


def check(name, ok, detail=None):
    checks.append({'check':name,'pass':bool(ok),'detail':detail})


init=json.loads((s/'provenance/initialization.json').read_text())
baseline=init['baseline_commit']
head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=repo,text=True).strip()
check('HEAD_unchanged_no_commit',head==baseline,head)
check('baseline_git_tree_unchanged',subprocess.check_output(['git','ls-tree','-r',head],cwd=repo)==(s/'provenance/baseline_git_tree.tsv').read_bytes())
changed=subprocess.check_output(['git','diff','--name-only',baseline,'--'],cwd=repo,text=True).splitlines()
allowed={x['path'] for x in init['baseline_registers']}
check('tracked_changes_only_authorized_registers',set(changed).issubset(allowed),changed)
registers=[]
for r in init['baseline_registers']:
    p=repo/r['path']; data=p.read_bytes(); prefix=data[:r['bytes']]
    valid=hashlib.sha256(prefix).hexdigest()==r['sha256']
    check('append_only_prefix_'+p.name,valid)
    check('register_has_new_stage_append_'+p.name,len(data)>r['bytes'] and b'fine-mapping-input-resolution-1.0' in data[r['bytes']:])
    rows=list(csv.reader(data.decode().splitlines(),delimiter='\t'))
    check('register_unique_ids_'+p.name,len({x[0] for x in rows[1:]})==len(rows)-1)
    check('register_column_counts_'+p.name,all(len(x)==len(rows[0]) for x in rows))
    registers.append(dict(path=r['path'],before_bytes=r['bytes'],before_sha256=r['sha256'],
        bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),appended_bytes=len(data)-r['bytes'],
        appended_rows=data[r['bytes']:].count(b'\n'),append_only=valid))
pre=repo/init['governing_preflight']
check('governing_preflight_freeze_unchanged',sha(pre/'provenance/freeze.json')==init['preflight_freeze_sha256'])
check('governing_preflight_ledger_unchanged',sha(pre/'provenance/artifact_checksums.tsv')==init['preflight_ledger_sha256'])
with (pre/'provenance/artifact_checksums.tsv').open() as f:
    prior=list(csv.DictReader(f,delimiter='\t'))
for r in prior:
    p=pre/r['path'];check('prior_preflight_payload_'+r['path'],p.is_file() and sha(p)==r['sha256'])
integration=json.loads((s/'tables/readiness_integration.json').read_text())
check('integration_final',integration['final'])
check('all16_obtainable_nonMHC_matrices_diagnosed',integration['diagnostic_matrices_completed']==16)
check('23_loci_21_nonMHC',integration['prospective_loci_total']==23 and integration['non_MHC_loci_total']==21)
check('zero_execution_cleared',integration['execution_cleared_loci']==0)
manifest=json.loads((s/'execution_inputs/manifest.json').read_text())
check('no_execution_package_mislabel',manifest['packages']==[] and not manifest['future_statistical_execution_authorized'])
extraction=json.loads((s/'ld/extraction_validation.json').read_text())
check('final_extraction_validation_PASS',extraction['mode']=='final' and extraction['status']=='PASS_EXTRACTION_AND_OUTPUT_VALIDATION_ONLY'
      and extraction['planned_loci']==16 and extraction['check_counts'].get('FAIL',0)==0
      and extraction['check_counts'].get('PENDING',0)==0)
check('final_report_exists',(s/'FINE_MAPPING_INPUT_RESOLUTION_REPORT.md').is_file())
untracked=subprocess.check_output(['git','ls-files','--others','--exclude-standard'],cwd=repo,text=True).splitlines()
rel=str(s.relative_to(repo))+'/'
check('new_files_within_authorized_stage',all(x.startswith(rel) for x in untracked),[x for x in untracked if not x.startswith(rel)])
result={'status':'PASS' if all(x['pass'] for x in checks) else 'FAIL','checks_passed':sum(x['pass'] for x in checks),
    'checks_failed':sum(not x['pass'] for x in checks),'baseline_commit':baseline,'current_HEAD':head,
    'shared_registers':registers,'prior_preflight_payloads_checked':len(prior),'checks':checks,
    'scope_attestation':{'fine_mapping_executed':False,'PIPs_or_credible_sets_computed':False,
        'model_or_candidate_scoring_executed':False,'target_genes_or_colocalization_executed':False,
        'manuscript_modified':False,'external_request_sent':False,'commit_or_push_performed_this_stage':False},
    'scope_evidence':'Authorized scripts, source/diagnostic logs, current outputs, independent review and unchanged historical Git tree. Future software archives/configuration are design artifacts, not inference runs.'}
(s/'provenance/stage_integrity_validation.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({k:result[k] for k in ['status','checks_passed','checks_failed','prior_preflight_payloads_checked']}))
if result['status']!='PASS':raise SystemExit(1)
