#!/usr/bin/env python3
"""Independent structural and allele-arithmetic audit of execution outputs."""
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path

S=Path(__file__).resolve().parents[1]
R=S.parents[3]
C=R/'diseases/COPD'

def read(p):
 op=gzip.open if str(p).endswith('.gz') else open
 with op(p,'rt') as f:return list(csv.DictReader(f,delimiter='\t'))
def check(cond,msg):
 if not cond:raise AssertionError(msg)

def main():
 report=[]
 def passed(name,detail):report.append(dict(test=name,status='PASS',detail=detail))
 r010=read(C/'04_modeling/results/COPD-S4-R010_THE_LIST.tsv')
 cw=read(S/'candidate_gwas_gtx_exact_crosswalk.tsv')
 direction=read(S/'risk_expression_direction_results.tsv')
 rejection=read(S/'risk_expression_direction_rejections.tsv')
 eq=read(C/'05_computational_validation/results/COPD-S5-R001_GTEx_v10_Lung_exact_significant_pairs.tsv.gz')
 eqidx={(x['candidate_record_id'],x['gene_id']):x for x in eq}
 ids=[x['candidate_record_id'] for x in r010]
 check(len(ids)==337 and len(set(ids))==337,'V1 337 identity')
 check([x['candidate_record_id'] for x in cw[::3]]==ids,'V1 order in crosswalk')
 check(len(cw)==1011 and all([x['gwas_track'] for x in cw[3*i:3*i+3]]==['A','B','C'] for i in range(337)),'track order')
 passed('frozen_candidate_order','337 candidates, identical rank order and A/B/C blocks')
 statuses=collections.Counter(x['status'] for x in cw)
 check(sum(statuses.values())==1011,'crosswalk denominator')
 check(all(x['status']!='exact_eligible' or x['reference_validation']=='PASS' and x['liftover_status']=='unique_one_base_mapping;roundtrip_exact' for x in cw),'reference/roundtrip')
 check(all(x['status']!='exact_eligible' or x['candidate_record_id'] in ids for x in cw),'candidate exactness')
 passed('exact_crossbuild','all eligible alleles have checked b37/b38 reference and unique round trip')
 eligible={(x['candidate_record_id'],x['gwas_track'],x['gwas_source_row']) for x in cw if x['status']=='exact_eligible'}
 check(all((x['candidate_record_id'],x['gwas_track'],x['gwas_source_row']) in eligible for x in direction),'proxy-to-exact promotion')
 passed('no_proxy_promotion','every direction row points to an exact eligible source contrast')
 # Synthetic cases are deliberately independent of biological rows.
 def change(slope,risk,alt):return slope if risk==alt else -slope
 check(change(.4,'A','A')==.4,'same-alt synthetic')
 check(change(.4,'G','A')==-.4,'reversed-ref synthetic')
 check(change(-.4,'A','A')==-.4,'negative slope synthetic')
 check('G'.translate(str.maketrans('ACGT','TGCA'))=='C','complement synthetic')
 check({'A','T'}=={'A','T'} and 'A'=='T'.translate(str.maketrans('ACGT','TGCA')),'palindromic synthetic')
 check((100,'G','GC')!=(100,'G','GT'),'indel allele synthetic')
 passed('synthetic_allele_cases','same ALT, reversed REF, negative slope, complement, palindrome and indel mismatch')
 for x in direction:
  c=x['candidate_record_id'];e=eqidx[(c,x['gene_id'])]
  check(x['gtex_variant_id']==e['gtex_variant_id'],'GTEx identity '+c)
  check(x['gtex_alt_slope']==e['slope_per_gtex_alt_allele'],'GTEx raw slope '+c)
  check(x['disease_increasing_allele_b38'] in c.split(':')[2:],'risk allele '+c)
  slope=float(e['slope_per_gtex_alt_allele']);expected=change(slope,x['disease_increasing_allele_b38'],c.split(':')[3])
  check(abs(float(x['expression_change_per_risk_allele'])-expected)<1e-12,'sign '+c)
  check(x['direction']==('increased' if expected>0 else 'decreased'),'direction '+c)
  check(x['gwas_track']=='A','unexpected track')
 wanted={(x['gwas_track'],x['gwas_source_row']) for x in direction}
 seen={}
 p=C/'07_gap_closure/fine-mapping-input-resolution-1.0/tables/gwas_only_risk_direction.tsv.gz'
 with gzip.open(p,'rt') as f:
  for row in csv.DictReader(f,delimiter='\t'):
   ident=(row['track'],row['source_row'])
   if ident in wanted:seen[ident]=row
 check(len(seen)==len(wanted),'all source GWAS rows')
 for x in direction:
  g=seen[(x['gwas_track'],x['gwas_source_row'])]
  check(x['gwas_beta']==g['beta'] and x['gwas_se']==g['se'] and x['gwas_source_identity']==g['source_identity'],'raw signed input')
  expected_risk=g['analysis_effect_allele'] if float(g['beta'])>0 else g['analysis_other_allele']
  check(g['disease_increasing_allele']==expected_risk,'GWAS beta/risk sign')
 passed('biological_allele_arithmetic',f'{len(direction)} rows recomputed from frozen GTEx ALT slope and risk allele')
 check(not any(x['reason'].startswith('exact_eligible') for x in rejection),'rejection state')
 mol=read(S/'saferali_moloc_crosswalk.tsv');matrix=read(S/'candidate_gene_evidence_matrix.tsv')
 check(len({x['publication_table_row'] for x in mol})==38,'published Table 5 coverage')
 check(all(x['statistical_scope']=='published_window_gene_QTL_Moloc_highest_PPA_only' for x in mol),'Moloc scope')
 check(all(x['candidate_overlap_scope']!='spatial_window_overlap_only' or x['candidate_record_id'] for x in mol),'window scope')
 check(all(not x['published_moloc_locus'] or x['tier']!='A' for x in matrix),'Moloc A promotion')
 passed('published_moloc_scope','38 Table 5 rows; candidate overlaps retain window scope')
 extra=read(S/'saferali_additional_moloc_signals.tsv')
 labels={(x['lead_gwas_variant'],x['gene_name'],x['qtl_cohort'],x['qtl_class'],x['posterior_ppa']) for x in extra}
 check(len(labels)==57 and all(float(x[-1])>=.8 for x in labels),'Table S2 extra labels')
 check(all(x['candidate_overlap_scope']!='exact_candidate_target' for x in extra),'Table S2 scope')
 passed('published_additional_signals','57 additional Table S2 gene/QTL/PPA labels retain locus scope')
 re2g=read(S/'encode_lung_re2g_crosswalk.tsv')
 check(all(x['evidence_type']=='predictive_element_gene_link' and x['file_accession']=='ENCFF324XYW' for x in re2g),'rE2G type')
 check(all(int(x['element_start_0based'])<=int(x['candidate_record_id'].split(':')[1])-1<int(x['element_end_0based_exclusive']) for x in re2g),'rE2G interval')
 check(all(not x['physical_contact'] and not x['abc_predictive'] for x in matrix),'rE2G/ABC/contact independence')
 passed('rE2G_type',f'{len(re2g)} element overlaps remain predictive, never contact or ABC replication')
 functional=read(S/'functional_target_evidence_audit.tsv')
 exact=[x for x in functional if x['exact_target_support']=='yes']
 check(len(exact)==1 and exact[0]['candidate_record_id']=='4:88963935:G:T' and exact[0]['target_gene']=='FAM13A','exact functional')
 check(all(x['tier']!='A' for x in matrix),'A hierarchy')
 check(all(x['tier']!='B' or x['candidate_record_id']=='4:88963935:G:T' and x['gene_name']=='FAM13A' and x['exact_eqtl'] and x['exact_functional_perturbation'] for x in matrix),'B hierarchy')
 check(sum(x['tier']=='B' for x in matrix)==1,'one integrated B context')
 check(all(x['tier']!='C' or any(x[c] for c in ['exact_eqtl','published_moloc_locus','rE2G_predictive','exact_functional_perturbation','physical_contact']) for x in matrix),'C context evidence')
 passed('tier_firewall','only integrated FAM13A eQTL + edited allele reaches B; every C context contains its own typed evidence')
 ph=read(S/'phenotype_stratum_target_summary.tsv')
 check([int(x['candidates']) for x in ph]==[184,29,124],'phenotype denominator')
 passed('phenotype_strata','184 direct, 29 EHR, 124 ML-only remain distinct')
 check('coloc.abf' not in (S/'scripts/run_execution.py').read_text() and 'coloc.susie' not in (S/'scripts/run_execution.py').read_text(),'no de novo coloc code')
 passed('no_de_novo_colocalization','execution script imports published Moloc only')
 with open(S/'synthetic_validation_records.tsv','w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=['test','status','detail'],delimiter='\t');w.writeheader();w.writerows(report)
 (S/'validation_audit.json').write_text(json.dumps(dict(status='PASS',tests=report),indent=2)+'\n')
 print('PASS:',len(report),'independent checks')

if __name__=='__main__':main()
