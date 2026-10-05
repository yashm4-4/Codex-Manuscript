#!/usr/bin/env python3
"""Final read-only V1 boundary audit and V2 artifact checksum manifest."""
import csv
import gzip
import hashlib
import json
import re
import subprocess
from pathlib import Path

V2=Path(__file__).resolve().parents[1]
ROOT=V2.parents[2]
PREFIX='COPD-V2-PHENO'


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''): h.update(b)
    return h.hexdigest()


def table(path):
    opener=gzip.open if path.suffix=='.gz' else open
    with opener(path,'rt',newline='') as f: return list(csv.DictReader(f,delimiter='\t'))


def main():
    checks=[]
    def check(name,value,detail):
        checks.append({'check':name,'status':'PASS' if value else 'FAIL','detail':detail})
        if not value: raise AssertionError((name,detail))
    snapshot=table(V2/'provenance/v1_snapshot.tsv')
    check('Original freeze snapshot SHA256',all(sha(V2.parent/r['artifact'])==r['sha256'] for r in snapshot),'Original R010, manuscript, final V1 validation unchanged')
    paths=subprocess.check_output(['git','diff','--name-only','875e995406ddf5d00ce77920d9e650ff128604c0','--','diseases/COPD'],cwd=ROOT,text=True).splitlines()
    check('Tracked V1 files unchanged from freeze commit',all(p.startswith('diseases/COPD/07_gap_closure/') for p in paths),f'{len(paths)} changed tracked paths, all inside V2')
    manifest=json.loads((V2/'provenance'/f'{PREFIX}_run_manifest.json').read_text())
    check('All consumed input hashes remain identical',all(sha(ROOT/r['path'])==r['sha256'] for r in manifest['inputs']),str(len(manifest['inputs'])))
    check('All analytical output hashes match run manifest',all(sha(ROOT/r['path'])==r['sha256'] for r in manifest['outputs']),str(len(manifest['outputs'])))
    qc=table(V2/'results'/f'{PREFIX}-validation_checks.tsv')
    check('All analysis checks PASS',all(r['status']=='PASS' for r in qc),str(len(qc)))
    summary=table(V2/'results'/f'{PREFIX}-R008_stratum_summary.tsv')
    lookup={(r['universe'],r['stratum']):int(r['supported']) for r in summary}
    files={'gws_tags':'R010_gws_tag_phenotype_support.tsv','all_15389':'R004_all_candidate_phenotype_support.tsv.gz','frozen_337':'R005_frozen_337_phenotype_support.tsv','components_153':'R006_component_phenotype_support.tsv','shortlist_12':'R007_shortlist_phenotype_support.tsv'}
    for universe,filename in files.items():
        rows=table(V2/'results'/f'{PREFIX}-{filename}')
        check(f'{universe}: summary agrees with per-record outputs',all(sum(r[col]=='True' for r in rows)==lookup[(universe,group)] for col,group in [('primary_retained','primary_direct'),('secondary_ehr_supported','secondary_ehr'),('ml_supported','ml_surrogate')]),f'{len(rows)} rows')
    check('Required detailed report exists',(V2/'results'/f'{PREFIX}_phenotype_robustness_report.md').is_file(),'Markdown report; V1 manuscript not rewritten')
    broken=[]
    report=V2/'results'/f'{PREFIX}_phenotype_robustness_report.md'
    for target in re.findall(r'\]\(([^)]+)\)',report.read_text()):
        if '://' not in target and not (report.parent/target.split('#')[0]).exists(): broken.append(target)
    check('Report local artifact links resolve',not broken,str(broken))
    registers=['gap_closure_result_register.tsv','gap_closure_decision_register.tsv','activity_log.tsv']
    for name in registers:
        with (V2/name).open(newline='') as f: rows=list(csv.reader(f,delimiter='\t'))
        check(f'Register integrity: {name}',len(set(map(len,rows)))==1 and len({r[0] for r in rows[1:]})==len(rows)-1,'Rectangular; identifiers unique')
    entries=table(V2/'gap_closure_result_register.tsv')
    check('Registered result paths exist',all((V2/r['primary_file']).is_file() for r in entries),'Planning and phenotype entries preserved')
    output=V2/'provenance'/f'{PREFIX}_final_validation.tsv'
    with output.open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=['check','status','detail'],delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(checks)
    artifact_paths=[]
    for path in V2.rglob('*'):
        if not path.is_file() or '__pycache__' in path.parts: continue
        if path.name.startswith(PREFIX) or path.name.startswith('adjudication_') or path.name in {'.gitignore','phenotype_adjudication_amendments.json','run_phenotype_robustness.py','validate_phenotype_module.py','gap_closure_result_register.tsv','gap_closure_decision_register.tsv','activity_log.tsv','README.md'}:
            if 'artifact_checksums' not in path.name:artifact_paths.append(path)
    checksum=V2/'provenance'/f'{PREFIX}_artifact_checksums.tsv'
    with checksum.open('w',newline='') as f:
        writer=csv.writer(f,delimiter='\t',lineterminator='\n');writer.writerow(['path','bytes','sha256'])
        for path in sorted(artifact_paths): writer.writerow([path.relative_to(ROOT),path.stat().st_size,sha(path)])
    print(json.dumps({'final_checks':len(checks),'all_pass':True,'checksummed_artifacts':len(artifact_paths),'v1_unchanged':True},indent=2))


if __name__=='__main__': main()
