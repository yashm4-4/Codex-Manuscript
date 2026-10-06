#!/usr/bin/env python3
"""Create provenance and a candidate-blind study metadata projection, not analysis."""
import csv, datetime, hashlib, io, json, pathlib, platform, subprocess

ROOT = pathlib.Path(__file__).resolve().parents[5]
STAGE = pathlib.Path(__file__).resolve().parents[1]
BASE = '1f41df67dcf18d626a3862a38877919d2d3114a4'
def sha(b): return hashlib.sha256(b).hexdigest()
def main():
    assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT).decode().strip()==BASE
    prov=STAGE/'provenance'; prov.mkdir(exist_ok=True)
    inputs=STAGE/'inputs'; inputs.mkdir(exist_ok=True)
    paths=['diseases/COPD/02_gwas/results/COPD-S2-R001_studies_core.tsv',
           'diseases/COPD/07_gap_closure/results/COPD-V2-PHENO-R001_accession_adjudication.tsv',
           'diseases/COPD/07_gap_closure/results/COPD-V2-PHENO-R002_cohort_overlap_provenance.tsv',
           'diseases/COPD/07_gap_closure/results/COPD-V2-PHENO-R002B_cohort_overlap_pairs.tsv',
           'diseases/COPD/02_gwas/data/copd_ancestries.tsv']
    records=[]
    for path in paths:
        b=(ROOT/path).read_bytes(); records.append({'path':path,'bytes':len(b),'sha256':sha(b)})
    raw=list(csv.DictReader((ROOT/paths[1]).open(),delimiter='\t'))
    excluded={'v1_gws_association_rows','v1_gws_tags','v1_candidate_records_supported','v1_frozen_337_supported','association_coverage_caveat'}
    fields=[k for k in raw[0] if k not in excluded]
    assert len(raw)==104 and len({r['study_accession'] for r in raw})==104
    with (inputs/'frozen_study_metadata.tsv').open('w') as out:
        w=csv.DictWriter(out,fields,delimiter='\t',lineterminator='\n');w.writeheader()
        w.writerows({k:r[k] for k in fields} for r in raw)
    registers=[]
    for name in ['activity_log.tsv','gap_closure_decision_register.tsv','gap_closure_result_register.tsv']:
        p='diseases/COPD/07_gap_closure/'+name;b=(ROOT/p).read_bytes()
        assert b==subprocess.check_output(['git','show',BASE+':'+p],cwd=ROOT)
        (prov/('baseline_'+name)).write_bytes(b)
        registers.append({'path':p,'bytes':len(b),'sha256':sha(b)})
    tree=subprocess.check_output(['git','ls-tree','-r',BASE],cwd=ROOT)
    (prov/'baseline_git_tree.tsv').write_bytes(tree)
    attachment=pathlib.Path('/home/maheshwarany2/.codex/attachments/34ad4eec-78fe-4009-a3b7-c6ed79b98d4b/pasted-text.txt')
    (prov/'user_authorization.txt').write_bytes(attachment.read_bytes())
    manifest={'stage':'fine-mapping-preflight-1.0','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'baseline_commit':BASE,'frozen_inputs':records,'baseline_registers':registers,'core_accessions':104,
        'candidate_support_columns_excluded':sorted(excluded),'candidate_membership_used':False,
        'python':platform.python_version(),'platform':platform.platform(),'scope_sha256':sha((STAGE/'PREFLIGHT_SCOPE.md').read_bytes())}
    (prov/'initialization.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'status':'INITIALIZED','core_accessions':104,'stage':str(STAGE)}))
if __name__=='__main__':main()
