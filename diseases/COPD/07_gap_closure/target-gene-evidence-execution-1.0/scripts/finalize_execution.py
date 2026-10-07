#!/usr/bin/env python3
"""Append stage records, then checksum-freeze the execution payload."""
import csv
import hashlib
import json
import subprocess
from pathlib import Path

STAGE=Path(__file__).resolve().parents[1]
ROOT=STAGE.parents[3]
GAP=ROOT/'diseases/COPD/07_gap_closure'
REL='target-gene-evidence-execution-1.0/'

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def append(name,row):
 p=GAP/name
 with open(p,newline='') as f:
  r=csv.reader(f,delimiter='\t');head=next(r);records=list(r)
 if any(x[0]==row[0] for x in records):return
 assert len(row)==len(head),(name,head,row)
 with open(p,'a',newline='') as f:
  csv.writer(f,delimiter='\t',lineterminator='\n').writerow(row)

def main():
 assert json.loads((STAGE/'validation_audit.json').read_text())['status']=='PASS'
 v=json.loads((STAGE/'validation_results.json').read_text())
 assert v['candidate_count']==337 and v['directed_candidate_gene_pairs']==387
 assert v['tiers_candidate_gene_pairs']=={'Unresolved':3218,'C':1216,'B':1}
 head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
 assert head=='dfb6390d24e26abf5b2755ee997d9b29776125f6',head
 # Verify every pre-existing register byte is an unchanged prefix of the
 # current append-only file before adding a new row.
 for name in ['activity_log.tsv','gap_closure_decision_register.tsv','gap_closure_result_register.tsv']:
  current=(GAP/name).read_bytes();old=subprocess.check_output(['git','show','HEAD:diseases/COPD/07_gap_closure/'+name],cwd=ROOT)
  assert current.startswith(old),f'PRIOR REGISTER CHANGED: {name}'
 append('activity_log.tsv',['COPD-V2-ACT-061','2026-10-07','complete_target_gene_execution',
  'Execute frozen exact candidate/GWAS/GTEx joins, published Moloc import, ENCODE lung rE2G and V1 functional re-audit',
  'frozen R010/R006/S5/PHENO/GWAS;target-gene preflight contracts;Saferali/ENCODE public tables',
  REL+'TARGET_GENE_EVIDENCE_EXECUTION_REPORT.md;'+REL+'candidate_gene_evidence_matrix.tsv',
  '337 rank unchanged; 36 exact usable GWAS/eQTL candidates, 387 directed pairs; no de novo coloc.'])
 append('activity_log.tsv',['COPD-V2-ACT-062','2026-10-07','frozen_stop_for_investigator_review',
  'Validate and checksum-freeze target-gene evidence execution; preserve failed joins and competing genes',
  REL+'validation_audit.json;'+REL+'source_acquisition_receipts.tsv',
  REL+'provenance/freeze.json;'+REL+'provenance/artifact_checksums.tsv',
  'No prior stage or manuscript edited; targeted coloc.abf input resolution remains a separate future stage.'])
 append('gap_closure_decision_register.tsv',['COPD-V2-DEC-065','2026-10-07','frozen_target_evidence_hierarchy',
  'Assign one FAM13A pair Tier B; keep published Moloc window and rE2G prediction at their original scopes',
  'Exact GTEx eQTL plus independently edited rs2013701 G/T allele supports FAM13A; no other pair meets Tier A/B gates.',
  REL+'target_evidence_tier_assignments.tsv;'+REL+'functional_target_evidence_audit.tsv'])
 append('gap_closure_decision_register.tsv',['COPD-V2-DEC-066','2026-10-07','targeted_abf_input_resolution_only',
  'A separate A-locus coloc.abf input-resolution stage is scientifically justified, but no coloc execution is cleared',
  'Exact signed GWAS/GTEx overlaps and FAM13A/PKD2 locus competition warrant full regional QTL and method-gate audit.',
  REL+'TARGET_GENE_EVIDENCE_EXECUTION_REPORT.md;'+REL+'risk_expression_direction_results.tsv'])
 append('gap_closure_result_register.tsv',['COPD-V2-TGE-EXEC-001','2026-10-07','complete_frozen_execution','target_gene_evidence_matrix',
  '36 exact usable GWAS/eQTL candidates, 387 risk-expression pairs; pair tiers A=0 B=1 C=1216 Unresolved=3218; no de novo coloc.',
  REL+'TARGET_GENE_EVIDENCE_EXECUTION_REPORT.md','no'])
 append('gap_closure_result_register.tsv',['COPD-V2-TGE-EXEC-FREEZE-001','2026-10-07','frozen_stop_for_review','target_gene_execution_freeze',
  'Checksum-bound typed evidence outputs and append-only registers; prior frozen stages unchanged.',
  REL+'provenance/freeze.json','no'])
 prov=STAGE/'provenance';prov.mkdir(exist_ok=True)
 paths=sorted(p for p in STAGE.rglob('*') if p.is_file() and 'sources' not in p.relative_to(STAGE).parts and 'provenance' not in p.relative_to(STAGE).parts and '__pycache__' not in p.parts)
 with open(prov/'artifact_checksums.tsv','w',newline='') as f:
  w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(['relative_path','bytes','sha256'])
  for p in paths:w.writerow([p.relative_to(STAGE),p.stat().st_size,sha(p)])
 freeze=dict(status='FROZEN_FOR_INVESTIGATOR_REVIEW',client_date='2026-10-07',base_commit=head,
  stage='target-gene-evidence-execution-1.0',manifest_sha256=sha(prov/'artifact_checksums.tsv'),
  payload_file_count=len(paths),source_receipts_sha256=sha(STAGE/'source_acquisition_receipts.tsv'),
  validation_sha256=sha(STAGE/'validation_audit.json'),results=v,
  register_sha256={name:sha(GAP/name) for name in ['activity_log.tsv','gap_closure_decision_register.tsv','gap_closure_result_register.tsv']},
  no_de_novo_colocalization=True,prior_stage_payloads_modified=False)
 (prov/'freeze.json').write_text(json.dumps(freeze,indent=2,sort_keys=True)+'\n')
 print(f'FROZEN: {len(paths)} payload files; manifest {freeze["manifest_sha256"]}')

if __name__=='__main__':main()
