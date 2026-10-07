#!/usr/bin/env python3
"""Validate, append register rows, and checksum-freeze this input-only stage."""
import csv
import gzip
import hashlib
import io
import json
import subprocess
from pathlib import Path

S=Path(__file__).resolve().parents[1]
G=S.parent
R=S.parents[3]
BASE='0adeb7b90b762f15d2c6f940bf361c0300978ce3'
DATE='2026-10-07'
A=G/'fine-mapping-input-resolution-1.0/tracks/A'

def sha(p):
    h=hashlib.sha256()
    with open(p,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def rows(p):
    op=gzip.open if str(p).endswith('.gz') else open
    with op(p,'rt',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def table(p,values,fields):
    with open(p,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(values)
def git(*args):return subprocess.check_output(['git',*args],cwd=R).decode().strip()

assert git('rev-parse','HEAD')==BASE
changed=set(git('diff','--name-only').splitlines()) | set(git('diff','--cached','--name-only').splitlines())
allowed_registers={'diseases/COPD/07_gap_closure/'+x for x in ['activity_log.tsv','gap_closure_decision_register.tsv','gap_closure_result_register.tsv']}
assert all(x in allowed_registers or x.startswith('diseases/COPD/07_gap_closure/targeted-coloc-fam13a-input-resolution-1.0/') for x in changed), 'Tracked prior payload changed'
v=json.loads((S/'validation_records.json').read_text())
g=rows(S/'fam13a_gwas_region.tsv');q=rows(S/'fam13a_gtex_v8_lung_full_cis.tsv');c=rows(S/'fam13a_exact_gwas_eqtl_crosswalk.tsv')
with gzip.open(S/'sources/QTD000271.permuted.tsv.gz','rt') as f:
    p=[x for x in csv.DictReader(f,delimiter='\t') if x['molecular_trait_id']=='ENSG00000138640']
with gzip.open(S/'sources/QTD000271.credible_sets.tsv.gz','rt') as f:
    cs=[x for x in csv.DictReader(f,delimiter='\t') if x['gene_id']=='ENSG00000138640']
assert len(g)==len(c)==12094 and len(q)==7010 and v['qtl_rows']==7264 and len(p)==1 and p[0]['n_variants']=='7010' and not cs
exact=[x for x in c if x['status']=='EXACT_PASS']
assert len(exact)==5065 and len(set(x['qtl_variant'] for x in exact))==5065
assert sum(float(x['native_p'])<5e-8 for x in c)==51
assert sum(float(x['native_p'])<5e-8 for x in exact)==39
assert sum(x['b38_reference_pass']=='False' for x in q)==0
assert sum(int(x['an'])<1020 for x in q)==99
assert all(x['chain_status']=='UNIQUE_ROUNDTRIP_EXACT' and x['reason']=='.' for x in exact)
assert not any('coloc.abf(' in x.read_text() for x in [S/'scripts/acquire_qtl.py',S/'scripts/resolve_inputs.py'])
checks=[('native_gwas_rows',12094),('qtl_source_rows',7264),('qtl_distinct_tested',7010),('official_permutation_n_variants',p[0]['n_variants']),('official_credible_set_rows',0),('permutation_p_beta',p[0]['p_beta']),('eligible_exact_intersection',5065),('unique_intersection_variants',5065),('gwas_significant_total',51),('gwas_significant_eligible',39),('qtl_ref_mismatches',0),('qtl_an_below_1020',99),('colocalization_executed','NO')]
table(S/'validation_checks.tsv',[{'check':k,'value':str(x),'status':'PASS'} for k,x in checks],['check','value','status'])

inputs=[
 ('frozen_track_a_native',A/'sources/GCST90016588_buildGRCh37.tsv'),
 ('frozen_track_a_harmonized_file',A/'sources/harmonised/33106845-GCST90016588-EFO_0006527.h.tsv.gz'),
 ('frozen_harmonized_locus_ledger',A/'results/harmonized_locus_ledger.tsv.gz'),
 ('frozen_locus_definition',A/'results/loci.tsv'),
 ('frozen_native_audit',A/'results/native_audit.json'),
 ('frozen_reference_checked_gwas',A/'results/locus_variants_reference_checked.tsv.gz'),
 ('grch37_p13_reference',G/'fine-mapping-input-resolution-1.0/references/GRCh37.p13.genome.fa'),
 ('prior_target_evidence_freeze',G/'target-gene-evidence-execution-1.0/provenance/freeze.json'),
]
inputs += [('targeted_'+p.name,p) for p in sorted((S/'sources').iterdir()) if p.is_file()]
table(S/'input_checksums.tsv',[{'source_id':n,'path':str(p.relative_to(R)),'bytes':p.stat().st_size,'sha256':sha(p)} for n,p in inputs],['source_id','path','bytes','sha256'])

updates={
 'activity_log.tsv':[
  ['COPD-V2-ACT-063',DATE,'complete_input_resolution','Acquire full FAM13A v8 Lung cis rows and audit frozen Track-A locus with exact reference-checked chain joins','frozen Track-A native/harmonized locus; Catalogue QTD000271 release-7 archive; UCSC chains/reference','targeted-coloc-fam13a-input-resolution-1.0/TARGETED_COLOC_FAM13A_INPUT_RESOLUTION_REPORT.md;targeted-coloc-fam13a-input-resolution-1.0/fam13a_exact_gwas_eqtl_crosswalk.tsv','No coloc or source-statistic changes.'],
  ['COPD-V2-ACT-064',DATE,'frozen_stop_for_investigator_review','Validate source denominators, statistical contracts, single-signal gate and checksum freeze','targeted-coloc-fam13a-input-resolution-1.0/validation_checks.tsv;targeted-coloc-fam13a-input-resolution-1.0/input_checksums.tsv','targeted-coloc-fam13a-input-resolution-1.0/provenance/freeze.json','NOT CLEARED; all prior frozen payloads unchanged.']],
 'gap_closure_decision_register.tsv':[
  ['COPD-V2-DEC-067',DATE,'targeted_fam13a_abf_not_cleared','Do not execute Track-A × GTEx v8 Lung FAM13A coloc.abf','Incomplete significant-GWAS capture, regional precision/N contract gaps, and unproven single-signal assumption; no gate is lowered for favorable FAM13A biology.','targeted-coloc-fam13a-input-resolution-1.0/TARGETED_COLOC_FAM13A_INPUT_RESOLUTION_REPORT.md;targeted-coloc-fam13a-input-resolution-1.0/future_coloc_abf_execution_specification.md']],
 'gap_closure_result_register.tsv':[
  ['COPD-V2-FAM13A-ABF-INPUT-001',DATE,'not_cleared_input_resolution','targeted_coloc_readiness','12,094 GWAS rows; 7,010 distinct QTL tests; 5,065 eligible exact overlaps; no FAM13A v8 credible set; no coloc.','targeted-coloc-fam13a-input-resolution-1.0/TARGETED_COLOC_FAM13A_INPUT_RESOLUTION_REPORT.md','no'],
  ['COPD-V2-FAM13A-ABF-INPUT-FREEZE-001',DATE,'frozen_stop_for_review','input_resolution_freeze','Checksum-bound targeted source and rowwise harmonization audit; prior freezes preserved.','targeted-coloc-fam13a-input-resolution-1.0/provenance/freeze.json','no']]
}
for name,add in updates.items():
    path=G/name
    baseline=subprocess.check_output(['git','show',BASE+':diseases/COPD/07_gap_closure/'+name],cwd=R)
    buf=io.StringIO(newline='');csv.writer(buf,delimiter='\t',lineterminator='\n').writerows(add)
    expected=baseline+buf.getvalue().encode()
    if path.read_bytes()==baseline:path.write_bytes(expected)
    else:assert path.read_bytes()==expected,name
    assert path.read_bytes().startswith(baseline)

prov=S/'provenance';prov.mkdir(exist_ok=True)
payload=[]
for p in sorted(S.rglob('*')):
    if not p.is_file() or 'sources' in p.parts or p.parent==prov:continue
    payload.append({'relative_path':str(p.relative_to(S)),'bytes':p.stat().st_size,'sha256':sha(p)})
table(prov/'artifact_checksums.tsv',payload,['relative_path','bytes','sha256'])
freeze={'base_commit':BASE,'client_date':DATE,'decision':'NOT CLEARED','no_colocalization_executed':True,'prior_stage_payloads_modified':False,'payload_file_count':len(payload),'manifest_sha256':sha(prov/'artifact_checksums.tsv'),'input_manifest_sha256':sha(S/'input_checksums.tsv'),'register_sha256':{name:sha(G/name) for name in updates},'raw_sources_excluded_from_git':True}
(prov/'freeze.json').write_text(json.dumps(freeze,indent=2)+'\n')
print(json.dumps(freeze,indent=2))
