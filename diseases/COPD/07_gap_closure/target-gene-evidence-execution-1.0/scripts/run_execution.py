#!/usr/bin/env python3
"""Deterministic, fail-closed target-gene evidence execution; no coloc."""
import collections
import csv
import gzip
import hashlib
import json
import math
import re
from pathlib import Path

import lxml.etree as ET
import pdfplumber
import pysam
from Bio import SeqIO

STAGE = Path(__file__).resolve().parents[1]
ROOT = STAGE.parents[3]
COPD = ROOT / 'diseases/COPD'
GAP = COPD / '07_gap_closure'
PRE = GAP / 'target-gene-evidence-preflight-1.0'
SRC = STAGE / 'sources'

INPUTS = {
 'LOCAL_R010': COPD/'04_modeling/results/COPD-S4-R010_THE_LIST.tsv',
 'LOCAL_R006': COPD/'04_modeling/results/COPD-S4-R006_candidate_target_evidence.tsv',
 'LOCAL_S5_GTEX': COPD/'05_computational_validation/results/COPD-S5-R001_GTEx_v10_Lung_exact_significant_pairs.tsv.gz',
 'LOCAL_PHENO': GAP/'results/COPD-V2-PHENO-R005_frozen_337_phenotype_support.tsv',
 'LOCAL_GWAS_RISK': GAP/'fine-mapping-input-resolution-1.0/tables/gwas_only_risk_direction.tsv.gz',
}
CONTRACTS = ['eqtl_harmonization_contract.md','risk_expression_direction_contract.md','proposed_target_evidence_hierarchy.md']
LIFT = COPD/'05_computational_validation/results/COPD-S5-R002_GRCh38_to_hg19_liftover_audit.tsv'
FA37 = GAP/'fine-mapping-input-resolution-1.0/references/GRCh37.p13.genome.fa'
FA38_WINDOWS = COPD/'04_modeling/data/COPD_candidate_variants_ref_2001bp.fa'
FUNC = COPD/'03_regulatory_landscape/results/COPD-S3-R004_functional_variant_evidence.tsv'

def sha(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
 return h.hexdigest()

def read(path):
 op=gzip.open if str(path).endswith('.gz') else open
 with op(path,'rt',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def write(name,rows,fields):
 with open(STAGE/name,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore',lineterminator='\n')
  w.writeheader();w.writerows(rows)

def comp(s):return s.translate(str.maketrans('ACGT','TGCA'))
def revcomp(s):return comp(s)[::-1]
def pal(a,b):return len(a)==len(b)==1 and a==comp(b)
def key(r):return r['candidate_record_id']
def stratum(r):
 if r['primary_retained']=='True':return 'direct_clinical_spirometric'
 if r['secondary_ehr_supported']=='True':return 'ehr_copd'
 return 'ml_surrogate_only'

def norm(pos,ref,alt,fetch):
 pos=int(pos);ref=ref.upper();alt=alt.upper()
 while len(ref)>1 and len(alt)>1 and ref[-1]==alt[-1]:ref=ref[:-1];alt=alt[:-1]
 while len(ref)>1 and len(alt)>1 and ref[0]==alt[0]:pos+=1;ref=ref[1:];alt=alt[1:]
 if len(ref)!=len(alt):
  for _ in range(1000):
   if pos<=1:break
   prior=fetch(pos-2,pos-1).upper()
   if len(prior)!=1 or ref[-1]!=alt[-1]:break
   pos-=1;ref=prior+ref[:-1];alt=prior+alt[:-1]
 return pos,ref,alt

def validate_inputs():
 receipts={x['receipt_id']:x for x in read(PRE/'provenance/acquisition_inspection_receipts.tsv')}
 rows=[]
 for code,path in INPUTS.items():
  digest=sha(path);expected=receipts[code]['sha256']
  if digest!=expected:raise RuntimeError(f'FROZEN INPUT DISCREPANCY: {code}: {digest} != {expected}')
  rows.append(dict(source_id=code,locator=str(path.relative_to(ROOT)),bytes=path.stat().st_size,sha256=digest,expected_sha256=expected,status='MATCH'))
 for name in CONTRACTS:
  p=PRE/name
  frozen={x['relative_path']:x['sha256'] for x in read(PRE/'provenance/artifact_checksums.tsv')}
  if sha(p)!=frozen[name]:raise RuntimeError('FROZEN CONTRACT DISCREPANCY: '+name)
  rows.append(dict(source_id=name,locator=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=sha(p),expected_sha256=frozen[name],status='MATCH'))
 for code,path in [('LIFTOVER_337',LIFT),('GRCH37_FASTA',FA37),('GRCH38_CANDIDATE_WINDOWS',FA38_WINDOWS),('V1_FUNCTIONAL',FUNC)]:
  rows.append(dict(source_id=code,locator=str(path.relative_to(ROOT)),bytes=path.stat().st_size,sha256=sha(path),expected_sha256='',status='READ_ONLY_REFERENCE'))
 for name,url in [('saferali_fulltext.xml','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12481885/fullTextXML'),('saferali_supplementary.zip','https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12481885/supplementaryFiles'),('saferali_mmc1.pdf','supplementaryFiles: mmc1.pdf'),('saferali_mmc1.txt','derived locally from supplementaryFiles: mmc1.pdf with pypdf'),('ENCFF324XYW.bed.gz','https://www.encodeproject.org/files/ENCFF324XYW/@@download/ENCFF324XYW.bed.gz'),('encode_ENCSR528UQX.json','https://www.encodeproject.org/annotations/ENCSR528UQX/?format=json&frame=embedded'),('ensembl_rs2013701.json','https://rest.ensembl.org/variation/human/rs2013701?content-type=application/json')]:
  p=SRC/name
  rows.append(dict(source_id=name,locator=url,bytes=p.stat().st_size,sha256=sha(p),expected_sha256='',status='PUBLIC_SOURCE_LOCAL_IGNORED'))
 for row in rows:row['client_date']='2026-10-07'
 write('source_acquisition_receipts.tsv',rows,['source_id','locator','client_date','bytes','sha256','expected_sha256','status'])

def parse_windows():
 text=(SRC/'saferali_mmc1.txt').read_text()
 table=text.split('Supplementary Table 1:',1)[1].split('Supplementary Table 2:',1)[0]
 # The PDF repeats rows on separate lines in the MHC. Parse each window
 # variant followed by its printed start/end, retaining the literal bounds.
 pat=re.compile(r'(\d+):(\d+):([ACGT]+):([ACGT]+)\s+(\d+)\s+(\d+)')
 out=[];seen=set()
 for m in pat.finditer(table):
  chrom,pos,ref,alt,start,end=m.groups();ident=f'{chrom}:{pos}:{ref}:{alt}'
  if (ident,start,end) in seen:continue
  seen.add((ident,start,end));out.append(dict(window_variant=ident,chrom=chrom,start=int(start),end=int(end)))
 # Some rows have an rsID between window variant and bounds. Capture those.
 pat2=re.compile(r'(\d+):(\d+):([ACGT]+):([ACGT]+)\s+rs\d+\s+\1:\2:\3:\4\s+(\d+)\s+(\d+)')
 for m in pat2.finditer(table):
  chrom,pos,ref,alt,start,end=m.groups();ident=f'{chrom}:{pos}:{ref}:{alt}'
  if (ident,start,end) not in seen:
   seen.add((ident,start,end));out.append(dict(window_variant=ident,chrom=chrom,start=int(start),end=int(end)))
 return out

def parse_additional_moloc(windows,candidates):
 # Supplementary Table 2 lists extra colocalized gene/QTL/PPA labels,
 # without a feature identifier or full regional statistical output.
 result=[];raw=[]
 with pdfplumber.open(SRC/'saferali_mmc1.pdf') as pdf:
  rows=[]
  for page in [8,9,10]:
   rows += pdf.pages[page].extract_tables()[0]
  for row in rows:
   if not row[0] or not re.search(r'\d+:\d+',row[0]):continue
   lead=(row[0] or '').replace('\n','').strip();best=(row[1] or '').replace('\n','').strip()
   if not re.fullmatch(r'\d+:\d+:[ACGT]+:[ACGT]+',lead):continue
   raw_cell=(row[5] or '').strip()
   raw.append(dict(lead_gwas_variant=lead,best_colocalized_snp=best,locus=row[2] or '',best_gene=row[3] or '',best_qtl_class=row[4] or '',
    other_colocalized_gene_labels_raw=raw_cell,gwas_identified_target_gene=row[6] or ''))
   fixed=raw_cell.replace('-\n','-').replace('_\n','_')
   fixed=re.sub(r'(0\.)\s*\n\s*(\d+)',r'\1\2',fixed)
   fixed=re.sub(r'(?<=_)\s*\n\s*','',fixed)
   tokens=re.findall(r'(LTRC|COPDGene)_(eQTL|sQTL|APAQTL)_([A-Za-z0-9.\-]+)_(0\.\d+)',fixed)
   # The PDF clips a few cells at column boundaries. Each correction is
   # transcribed from the same supplement's extractable page text (pages 9–11).
   corrections={
    ('5:157510035:A:G','COPDGene','sQTL','ADAM19'):'0.89',
    ('7:100032719:C:T','COPDGene','APAQTL','ZKSCAN1'):'0.97',
    ('7:100032719:C:T','COPDGene','APAQTL','STAG3L5P'):'0.96',
    ('17:38678826:C:T','COPDGene','eQTL','AC006449.6'):'0.88',
    ('1:45480964:G:T','COPDGene','APAQTL','MUTYH'):'0.93',
    ('15:49692513:G:C','COPDGene','sQTL','FAM227B'):'0.92',
    ('10:79946568:A:G','COPDGene','sQTL','FAM213A'):'0.95',
   }
   tokmap={(co,qtl,gene):ppa for co,qtl,gene,ppa in tokens}
   for (l,co,qtl,gene),ppa in corrections.items():
    if l==lead:tokmap[(co,qtl,gene)]=ppa
   tokens=[(*k,v) for k,v in tokmap.items()]
   lead_windows=[w for w in windows if w['window_variant']==lead]
   for cohort,qtl,gene,ppa in tokens:
    for w in lead_windows or [dict(window_variant=lead,chrom=lead.split(':')[0],start='',end='')]:
     overlaps=[c for c in candidates if w['start']!='' and c['chromosome_grch38']==w['chrom'] and w['start']<=int(c['position_grch38'])<=w['end']]
     for c in overlaps or [None]:
      result.append(dict(lead_gwas_variant=lead,published_best_snp=best,window_start_grch38=w['start'],window_end_grch38=w['end'],
       window_resolution='unique_Table_S1_lead' if len(lead_windows)==1 else 'ambiguous_or_unparsed_Table_S1_lead',
       gene_name=gene,posterior_ppa=ppa,qtl_cohort=cohort,qtl_class=qtl,tissue='lung' if cohort=='LTRC' else 'whole_blood',
       molecular_feature='not_reported_in_Table_S2',candidate_record_id=key(c) if c else '',
       candidate_overlap_scope='spatial_window_overlap_only' if c else 'no_frozen_candidate_in_window',
       statistical_scope='published_Table_S2_additional_window_gene_QTL_label',independent_disease_replication='no_reuses_Sakornsakolpat_GWAS'))
 write('saferali_supplementary_table2_audit.tsv',raw,list(raw[0]))
 write('saferali_additional_moloc_signals.tsv',result,list(result[0]) if result else ['lead_gwas_variant','gene_name','posterior_ppa'])
 return result

def table2_best_to_lead():
 out={}
 with pdfplumber.open(SRC/'saferali_mmc1.pdf') as pdf:
  for page in [8,9,10]:
   for row in pdf.pages[page].extract_tables()[0]:
    lead=(row[0] or '').replace('\n','').strip();best=(row[1] or '').replace('\n','').strip()
    if re.fullmatch(r'\d+:\d+:[ACGT]+:[ACGT]+',lead) and re.fullmatch(r'\d+:\d+:[ACGT]+:[ACGT]+',best):
     if best in out and out[best]!=lead:raise RuntimeError('Table S2 conflicting best-SNP lead mapping: '+best)
     out[best]=lead
 return out

def main():
 validate_inputs()
 candidates=read(INPUTS['LOCAL_R010']);pheno={key(x):x for x in read(INPUTS['LOCAL_PHENO'])}
 assert len(candidates)==337 and len(pheno)==337
 assert [key(x) for x in candidates]==[key(x) for x in read(INPUTS['LOCAL_PHENO'])]
 lift={key(x):x for x in read(LIFT)};assert len(lift)==337
 gtex=collections.defaultdict(list)
 for x in read(INPUTS['LOCAL_S5_GTEX']):gtex[key(x)].append(x)
 r006=collections.defaultdict(list)
 for x in read(INPUTS['LOCAL_R006']):r006[key(x)].append(x)
 refwin={x.id:str(x.seq).upper() for x in SeqIO.parse(FA38_WINDOWS,'fasta')}
 ref37=pysam.FastaFile(str(FA37))
 mapped={};at37=collections.defaultdict(list)
 for c in candidates:
  k=key(c);li=lift[k];ch=c['chromosome_grch38'];pos=int(c['position_grch38']);ref=c['ref'].upper();alt=c['alt'].upper()
  assert li['liftover_status']=='unique_one_base_mapping' and li['roundtrip_status']=='roundtrip_exact' and li['coordinate_match_eligible']=='true'
  assert li['chromosome_grch38']==ch and int(li['position_grch38'])==pos
  w=refwin[k];assert w[1000:1000+len(ref)]==ref
  p37=int(li['position_hg19']);ch37=li['chromosome_hg19']
  r37=ref37.fetch('chr'+ch37,p37-1,p37-1+len(ref)).upper()
  if r37==ref:expected=(ref,alt,'forward')
  elif r37==revcomp(ref):expected=(revcomp(ref),revcomp(alt),'reverse_complement')
  else:expected=None
  mapped[k]=dict(ch37=ch37,p37=p37,reference37=r37,expected=expected,reference38=ref)
  at37[(ch37,p37)].append(k)
 # Only the signed, frozen GWAS locus ledger is searched. Rows at the
 # mapped coordinate are retained even when their allele/QC state fails.
 atpos=collections.defaultdict(list)
 with gzip.open(INPUTS['LOCAL_GWAS_RISK'],'rt',newline='') as f:
  for row in csv.DictReader(f,delimiter='\t'):
   match=re.match(r'^(?:chr)?([^:]+):(\d+):([ACGT]+):([ACGT]+)$',row['analysis_identity'])
   if match and (match[1],int(match[2])) in at37:atpos[(match[1],int(match[2]),row['track'])].append(row)
 cw=[];directions=[];reject=[];passby=collections.defaultdict(list)
 for c in candidates:
  k=key(c);m=mapped[k];st=stratum(pheno[k]);glist=gtex[k]
  for track in 'ABC':
   raw=atpos.get((m['ch37'],m['p37'],track),[]);accepted=[];fail=[]
   for g in raw:
    parts=g['analysis_identity'].split(':');gr,ga=parts[2].upper(),parts[3].upper()
    reason=''
    if m['expected'] is None:reason='BUILD_REF_MISMATCH'
    elif ref37.fetch('chr'+m['ch37'],m['p37']-1,m['p37']-1+len(gr)).upper()!=gr:reason='GWAS_REF_MISMATCH'
    elif not (gr==m['expected'][0] and ga==m['expected'][1]):reason='ALLELE_CONFLICT_AT_POSITION'
    elif pal(c['ref'],c['alt']):reason='PALINDROMIC_UNRESOLVED'
    elif not g['beta'] or not math.isfinite(float(g['beta'])) or float(g['beta'])==0:reason='ZERO_OR_NONFINITE_GWAS_BETA'
    elif {g['analysis_effect_allele'].upper(),g['analysis_other_allele'].upper()}!={gr,ga}:reason='EFFECT_OTHER_ALLELE_CONFLICT'
    elif g['disease_increasing_allele'].upper()!=(g['analysis_effect_allele'].upper() if float(g['beta'])>0 else g['analysis_other_allele'].upper()):reason='RISK_SIGN_CONFLICT'
    elif g['qc_status']!='SIGNED_STATISTIC_PASS':reason='GWAS_SIGNED_QC_FAILURE'
    elif g['direction_status'] not in ('VERIFIED_GWAS_SIGNED_CONTRAST_NOT_CAUSALITY','VERIFIED_SIGNED_GWAS_CONTRAST','VERIFIED_SIGNED_GWAS_CONTRAST_ESTIMATED_DIRECTION'):reason='UNVERIFIED_SIGNED_CONTRAST'
    elif g['harmonization_status'] not in ('REFERENCE_VERIFIED_SNP','REFERENCE_VERIFIED_INDEL','REFERENCE_VERIFIED','PASS') and not g['harmonization_status'].startswith('REFERENCE_VERIFIED'):reason='UNVERIFIED_GWAS_IDENTITY'
    elif len(c['ref'])!=len(c['alt']):
     f37=lambda a,b:ref37.fetch('chr'+m['ch37'],a,b)
     f38=lambda a,b:refwin[k][1000+(a-int(c['position_grch38'])):1000+(b-int(c['position_grch38']))]
     n37=norm(m['p37'],gr,ga,f37);n38=norm(c['position_grch38'],c['ref'],c['alt'],f38)
     if (n37[1],n37[2])!=(n38[1],n38[2]) and m['expected'][2]=='forward':reason='INDEL_NORMALIZATION_CONFLICT'
    if reason:fail.append((g,reason))
    else:accepted.append(g)
   # A repeated source contrast is admissible only when all copies agree.
   sig={(x['analysis_identity'],x['effect_allele'],x['other_allele'],x['beta'],x['disease_increasing_allele']) for x in accepted}
   if len(sig)>1:
    fail += [(x,'DUPLICATE_CONFLICTING_SIGNED_ROWS') for x in accepted];accepted=[]
   elif len(accepted)>1:accepted=accepted[:1]
   if accepted:status='exact_eligible'
   elif m['expected'] is None:status='build_ref_mismatch'
   elif any(z=='PALINDROMIC_UNRESOLVED' for _,z in fail):status='palindromic_unresolved'
   elif any(z=='DUPLICATE_CONFLICTING_SIGNED_ROWS' for _,z in fail):status='duplicate_conflicting'
   elif raw:status='otherwise_rejected'
   else:status='absent_from_gwas_ledger'
   for g in accepted or [None]:
    cr=dict(candidate_record_id=k,v1_rank=c['predicted_causal_priority_rank'],phenotype_stratum=st,gwas_track=track,
      gwas_phenotype={'A':'direct_ever_smoker_spirometric_COPD','B':'direct_Japanese_clinical_COPD','C':'secondary_EHR_J44_COPD'}[track],
      grch38_identity=k,grch37_coordinate=f"{m['ch37']}:{m['p37']}",grch37_reference=m['reference37'],
      liftover_status=lift[k]['liftover_status']+';'+lift[k]['roundtrip_status'],reference_validation='PASS' if m['expected'] else 'FAIL',
      candidate_strand=m['expected'][2] if m['expected'] else '',gwas_source_row=g['source_row'] if g else '',gwas_source_identity=g['source_identity'] if g else '',
      gwas_analysis_identity=g['analysis_identity'] if g else '',gwas_effect_allele=g['effect_allele'] if g else '',gwas_other_allele=g['other_allele'] if g else '',
      disease_increasing_allele=g['disease_increasing_allele'] if g else '',gwas_beta=g['beta'] if g else '',gwas_se=g['se'] if g else '',
      exact_gtex_pair_count=len(glist),gtex_status='exact_significant_pair_observed' if glist else 'absent_from_significant_pairs',
      status=status,rejection_reasons=';'.join(sorted(set(z for _,z in fail))))
    cw.append(cr)
    if g:
     passby[(k,track)].append(g)
     for eq in glist:
      slope=float(eq['slope_per_gtex_alt_allele']);risk=g['disease_increasing_allele'].upper();orient=m['expected'][2]
      if orient=='reverse_complement':risk=revcomp(risk)
      if risk not in (c['ref'],c['alt']):raise RuntimeError('Risk allele conversion failed: '+k)
      if not math.isfinite(slope) or slope==0:
       reject.append(dict(candidate_record_id=k,gwas_track=track,gene_id=eq['gene_id'],reason='ZERO_OR_NONFINITE_GTEX_SLOPE'));continue
      multiplier=1 if risk==c['alt'] else -1
      d=dict(candidate_record_id=k,phenotype_stratum=st,gwas_track=track,gwas_phenotype=cr['gwas_phenotype'],gwas_source_row=g['source_row'],
       gwas_source_identity=g['source_identity'],gwas_beta=g['beta'],gwas_se=g['se'],gwas_p=g['p'],source_effect_allele=g['effect_allele'],
       source_other_allele=g['other_allele'],disease_increasing_allele_b38=risk,gtex_variant_id=eq['gtex_variant_id'],gene_id=eq['gene_id'],
       gene_id_versionless=eq['gene_id_versionless'],gene_name=eq['gene_name'],tissue='GTEx_v10_bulk_Lung',gtex_effect_allele=eq['alt'],
       gtex_alt_slope=eq['slope_per_gtex_alt_allele'],gtex_slope_se=eq['slope_se'],harmonization_operation=orient+';'+('same_alt' if multiplier==1 else 'reversed_ref'),
       multiplier=multiplier,expression_change_per_risk_allele=multiplier*slope,direction='increased' if multiplier*slope>0 else 'decreased',
       interpretation='allele_arithmetic_association_not_colocalization_or_mediation',ukb_overlap='A_C_overlap' if track in 'AC' else 'not_A_C_UKB')
      directions.append(d)
   if not accepted:
    for eq in glist or [None]:reject.append(dict(candidate_record_id=k,gwas_track=track,gene_id=eq['gene_id'] if eq else '',reason=status+';'+cr['rejection_reasons'] if glist else status+';NO_EXACT_SIGNIFICANT_GTEX_PAIR'))
 write('candidate_gwas_gtx_exact_crosswalk.tsv',cw,list(cw[0]))
 write('risk_expression_direction_results.tsv',directions,list(directions[0]) if directions else ['candidate_record_id','gene_id','direction'])
 write('risk_expression_direction_rejections.tsv',reject,['candidate_record_id','gwas_track','gene_id','reason'])

 # Published Moloc Table 5: one highest-PPA QTL record per reported window.
 xml=ET.parse(str(SRC/'saferali_fulltext.xml'));table=xml.xpath('.//table-wrap[@id="tbl5"]')[0]
 windows=parse_windows();best_to_lead=table2_best_to_lead();mol=[]
 for ix,tr in enumerate(table.xpath('.//tr')[2:],1):
  v=[''.join(c.itertext()).strip().replace('−','-') for c in tr.xpath('./th|./td')]
  mm=re.match(r'chr(\d+):(\d+)([ACGT]+)>([ACGT]+)',v[0]);assert mm and len(v)==9
  ch,pos,ref,alt=mm.groups();pos=int(pos);best=f'{ch}:{pos}:{ref}:{alt}'
  possible=[w for w in windows if w['chrom']==ch and w['start']<=pos<=w['end']]
  lead=best_to_lead.get(best)
  matching=[w for w in possible if w['window_variant']==lead] if lead else possible
  ws=matching or [dict(window_variant='',chrom=ch,start='',end='')]
  for win in ws:
   overlapping=[c for c in candidates if c['chromosome_grch38']==ch and win['start']!='' and win['start']<=int(c['position_grch38'])<=win['end']]
   if not overlapping:overlapping=[None]
   for c in overlapping:
    exact=bool(c and key(c)==f'{ch}:{pos}:{ref}:{alt}')
    mol.append(dict(publication_table_row=ix,published_best_snp=f'{ch}:{pos}:{ref}:{alt}',posterior_ppa=v[1],gene_name=v[2],
      qtl_class=v[3].split('_',1)[1],qtl_cohort=v[3].split('_',1)[0],tissue='lung' if v[3].startswith('LTRC') else 'whole_blood',
      molecular_feature=v[4],qtl_p=v[5],qtl_effect=v[6],gwas_p=v[7],gwas_effect=v[8],
      source_gwas='Sakornsakolpat_2019_COPD_case_control',window_variant=win['window_variant'],window_start_grch38=win['start'],window_end_grch38=win['end'],
      window_resolution='unique_Table_S1_window' if len(matching)==1 else 'possible_Table_S1_window_not_uniquely_identified',
      window_link_method='Table_S2_lead_plus_Table_S1_bounds' if lead and len(matching)==1 else 'Table_S1_coordinate_containment',
      candidate_record_id=key(c) if c else '',candidate_overlap_scope='best_snp_allele_identity_metadata_only' if exact else ('spatial_window_overlap_only' if c else 'no_frozen_candidate_in_window'),
      statistical_scope='published_window_gene_QTL_Moloc_highest_PPA_only',independent_disease_replication='no_reuses_Sakornsakolpat_GWAS'))
 write('saferali_moloc_crosswalk.tsv',mol,list(mol[0]))
 additional_moloc=parse_additional_moloc(windows,candidates)

 enc=json.loads((SRC/'encode_ENCSR528UQX.json').read_text())
 f=[x for x in enc['files'] if x['accession']=='ENCFF324XYW'][0]
 assert f['assembly']=='GRCh38' and f['output_type']=='thresholded element gene links' and 'encode_re2g' in ' '.join(f['aliases'])
 bychr=collections.defaultdict(list)
 for c in candidates:bychr['chr'+c['chromosome_grch38']].append(c)
 re2g=[]
 with gzip.open(SRC/'ENCFF324XYW.bed.gz','rt') as h:
  reader=csv.DictReader(h,delimiter='\t')
  for e in reader:
   for c in bychr.get(e['#chr'],[]):
    p=int(c['position_grch38'])-1
    if int(e['start'])<=p<int(e['end']):
     re2g.append(dict(candidate_record_id=key(c),phenotype_stratum=stratum(pheno[key(c)]),accession='ENCSR528UQX',file_accession='ENCFF324XYW',
      donor='ENCDO528BHB',context='non_COPD_bulk_lung_DNase',assembly='GRCh38',pipeline='distal-regulation-encode_re2g_1.0.0',element_chrom=e['#chr'],element_start_0based=e['start'],element_end_0based_exclusive=e['end'],
      element_class=e['class'],gene_id=e['TargetGeneEnsemblID'],gene_name=e['TargetGene'],score=e['Score'],threshold_semantics='ENCODE_supplied_thresholded_predictions;model_operating_point_70pct_CRISPR_link_recall;numeric_cutoff_not_rederived',
      evidence_type='predictive_element_gene_link',independence_family='ENCODE_rE2G_ABC_shared_lineage'))
 write('encode_lung_re2g_crosswalk.tsv',re2g,list(re2g[0]) if re2g else ['candidate_record_id','gene_id','score'])

 # Existing curated literature only. The single exact rs2013701 assay is
 # checked against an independent GRCh38 variation mapping and tested G/T.
 var=json.loads((SRC/'ensembl_rs2013701.json').read_text())
 assert any(x['assembly_name']=='GRCh38' and x['location']=='4:88963935-88963935' and 'G' in x['allele_string'] and 'T' in x['allele_string'] for x in var['mappings'])
 functional=[]
 for x in read(FUNC):
  ids=set(re.findall(r'rs\d+',x['variant_or_haplotype']))
  hits=[c for c in candidates if c['selected_ensembl_variation_id'] in ids]
  if not hits:hits=[None]
  for c in hits:
   exact=bool(c and x['variant_evidence_id']=='COPD-VAR-001' and key(c)=='4:88963935:G:T')
   functional.append(dict(evidence_id=x['variant_evidence_id'],source_id=x['source_id'],published_variant=x['variant_or_haplotype'],candidate_record_id=key(c) if c else '',
     target_gene=x['target_gene'],cell_context=x['cell_context'],tested_unit=x['tested_unit'],endogenous_allele_tested=x['endogenous_allele_tested'],
     scope='exact_candidate_exact_G_T_endogenous_allele' if exact else ('rsid_overlap_allele_or_build_unresolved' if c else 'linked_or_nearby_locus_only'),
     exact_target_support='yes' if exact else 'no',readout=x['functional_readouts'],caveat=x['caveat']))
 write('functional_target_evidence_audit.tsv',functional,list(functional[0]))

 # Union of GENCODE/Catalog candidates, exact eGenes, rE2G and published
 # Moloc genes. Symbol-only Moloc links remain separate from ENSG IDs.
 gene_rows={};name_to_ids=collections.defaultdict(set)
 for rows in r006.values():
  for x in rows:
   if x['gencode_gene_id']:name_to_ids[x['gene_name']].add(x['gencode_gene_id'])
 for rows in gtex.values():
  for x in rows:name_to_ids[x['gene_name']].add(x['gene_id_versionless'])
 for x in re2g:name_to_ids[x['gene_name']].add(x['gene_id'])
 def ensure(k,gid,name,context):
  z=(k,gid,name,context)
  if z not in gene_rows:gene_rows[z]=dict(candidate_record_id=k,phenotype_stratum=stratum(pheno[k]),gene_id=gid,gene_name=name,tissue_context=context,
   proximity_tss='',gwas_catalog_mapped_gene='',selected_locus_membership='',exact_eqtl='',risk_expression_direction='',published_moloc_locus='',rE2G_predictive='',abc_predictive='',
   physical_contact='',exact_functional_perturbation='',coding_overlap='',contradictions='',tier='Unresolved',tier_basis='proximity_or_no_usable_target_evidence')
  return gene_rows[z]
 for k,rows in r006.items():
  for x in rows:
   z=ensure(k,x['gencode_gene_id'],x['gene_name'],'unspecified')
   method=x['mapping_method']
   if 'TSS' in method or 'tss' in method:z['proximity_tss']+=';'+method+':'+x['absolute_tss_distance_bp']
   if method=='GWAS_Catalog_mapped_gene_via_linked_tag':z['gwas_catalog_mapped_gene']+=';'+method
   if method=='selected_COPD_GWAS_locus_membership':z['selected_locus_membership']='selected_COPD_GWAS_locus_membership'
 for k,rows in gtex.items():
  for x in rows:
   z=ensure(k,x['gene_id_versionless'],x['gene_name'],'GTEx_v10_bulk_Lung')
   z['exact_eqtl']=x['gene_id']+';slope='+x['slope_per_gtex_alt_allele']+';p='+x['pval_nominal']
 for x in directions:
  z=ensure(x['candidate_record_id'],x['gene_id_versionless'],x['gene_name'],'GTEx_v10_bulk_Lung')
  z['risk_expression_direction']+=';'+x['gwas_track']+':'+x['direction']
 for x in re2g:
  z=ensure(x['candidate_record_id'],x['gene_id'],x['gene_name'],'non_COPD_bulk_lung_DNase')
  z['rE2G_predictive']+=';ENCFF324XYW:'+x['element_chrom']+':'+x['element_start_0based']+'-'+x['element_end_0based_exclusive']+':'+x['score']
 for x in mol:
  k=x['candidate_record_id']
  if not k or x['window_resolution']!='unique_Table_S1_window':continue
  ids=name_to_ids[x['gene_name']]
  gid=next(iter(ids)) if len(ids)==1 else ''
  z=ensure(k,gid,x['gene_name'],x['tissue'])
  z['published_moloc_locus']+=';Table5:'+str(x['publication_table_row'])+':'+x['qtl_cohort']+'_'+x['qtl_class']+':PPA='+x['posterior_ppa']
  if len(ids)>1:z['contradictions']='symbol_to_ENSG_ambiguous'
 for x in additional_moloc:
  k=x['candidate_record_id']
  if not k or x['window_resolution']!='unique_Table_S1_lead':continue
  ids=name_to_ids[x['gene_name']];gid=next(iter(ids)) if len(ids)==1 else ''
  z=ensure(k,gid,x['gene_name'],x['tissue'])
  z['published_moloc_locus']+=';TableS2:'+x['qtl_cohort']+'_'+x['qtl_class']+':PPA='+x['posterior_ppa']
  if len(ids)>1:z['contradictions']='symbol_to_ENSG_ambiguous'
 for x in functional:
  if x['exact_target_support']!='yes':continue
  k=x['candidate_record_id'];ids=name_to_ids[x['target_gene']];gid=next(iter(ids)) if len(ids)==1 else ''
  z=ensure(k,gid,x['target_gene'],'16HBE_bronchial_epithelial')
  z['exact_functional_perturbation']=x['evidence_id']+':endogenous_G_to_T_FAM13A_expression'
 for c in candidates:
  k=key(c)
  for gene in filter(None,c['coding_CDS_gene_names'].split(';')):
   ids=name_to_ids[gene];gid=next(iter(ids)) if len(ids)==1 else ''
   z=ensure(k,gid,gene,'coding_annotation');z['coding_overlap']='exact_CDS_overlap_annotation_only'
 # A combined-context record explicitly carries the two independent lines.
 # Source-specific and proximity-only rows keep their own interpretation.
 k='4:88963935:G:T';fam=[z for z in gene_rows.values() if z['candidate_record_id']==k and z['gene_name']=='FAM13A']
 qtl=next((z for z in fam if z['exact_eqtl']),None);assay=next((z for z in fam if z['exact_functional_perturbation']),None)
 if qtl and assay:
  z=ensure(k,qtl['gene_id'],'FAM13A','GTEx_v10_bulk_Lung+16HBE_bronchial_epithelial')
  z['exact_eqtl']=qtl['exact_eqtl'];z['risk_expression_direction']=qtl['risk_expression_direction']
  z['exact_functional_perturbation']=assay['exact_functional_perturbation']
 for z in gene_rows.values():
  eq=bool(z['exact_eqtl']);fun=bool(z['exact_functional_perturbation'])
  if fun and eq and not z['contradictions']:
   z['tier']='B';z['tier_basis']='exact_candidate_eQTL_plus_independent_validated_allele_specific_target_assay;combined_context'
  elif eq or z['published_moloc_locus'] or z['rE2G_predictive'] or fun or z['physical_contact']:
   z['tier']='C';z['tier_basis']='suggestive_typed_molecular_or_predictive_evidence_in_this_context'
  if z['coding_overlap'] and not (eq or z['published_moloc_locus'] or z['rE2G_predictive'] or fun):
   z['tier']='Unresolved';z['tier_basis']='coding_overlap_not_regulatory_target'
 matrix=list(gene_rows.values());matrix.sort(key=lambda z:(next(i for i,c in enumerate(candidates) if key(c)==z['candidate_record_id']),z['gene_name'],z['tissue_context']))
 write('candidate_gene_evidence_matrix.tsv',matrix,list(matrix[0]))
 assignments=[{a:z[a] for a in ['candidate_record_id','phenotype_stratum','gene_id','gene_name','tissue_context','tier','tier_basis','contradictions']} for z in matrix]
 write('target_evidence_tier_assignments.tsv',assignments,list(assignments[0]))
 summary=[]
 for st in ['direct_clinical_spirometric','ehr_copd','ml_surrogate_only']:
  cs=[c for c in candidates if stratum(pheno[key(c)])==st];ks={key(c) for c in cs};cws=[x for x in cw if x['candidate_record_id'] in ks]
  summary.append(dict(phenotype_stratum=st,candidates=len(cs),candidates_with_exact_gwas= len({x['candidate_record_id'] for x in cws if x['status']=='exact_eligible'}),
    candidates_with_significant_gtex=len({k for k in ks if gtex[k]}),candidates_with_usable_gwas_gtex=len({x['candidate_record_id'] for x in directions if x['candidate_record_id'] in ks}),
    directed_candidate_gene_track_rows=sum(x['candidate_record_id'] in ks for x in directions),
    moloc_window_overlap_candidates=len({x['candidate_record_id'] for x in mol if x['candidate_record_id'] in ks and x['window_resolution']=='unique_Table_S1_window'}),
    re2g_overlap_candidates=len({x['candidate_record_id'] for x in re2g if x['candidate_record_id'] in ks}),
    tier_A_rows=sum(z['candidate_record_id'] in ks and z['tier']=='A' for z in matrix),tier_B_rows=sum(z['candidate_record_id'] in ks and z['tier']=='B' for z in matrix),
    tier_C_rows=sum(z['candidate_record_id'] in ks and z['tier']=='C' for z in matrix),tier_Unresolved_rows=sum(z['candidate_record_id'] in ks and z['tier']=='Unresolved' for z in matrix)))
 write('phenotype_stratum_target_summary.tsv',summary,list(summary[0]))
 denominators=[]
 all_statuses=['exact_eligible','absent_from_gwas_ledger','absent_from_significant_gtex_pairs','ambiguous','palindromic_unresolved','build_ref_mismatch','duplicate_conflicting','otherwise_rejected']
 for st in ['direct_clinical_spirometric','ehr_copd','ml_surrogate_only']:
  for track in 'ABC':
   subset=[x for x in cw if x['phenotype_stratum']==st and x['gwas_track']==track]
   n=collections.Counter(x['status'] for x in subset)
   record=dict(phenotype_stratum=st,gwas_track=track,candidates=len(subset),**{q:n[q] for q in all_statuses if q!='absent_from_significant_gtex_pairs'})
   record['absent_from_significant_gtex_pairs']=sum(int(x['exact_gtex_pair_count'])==0 for x in subset)
   record['exact_gwas_and_gtex']=sum(x['status']=='exact_eligible' and int(x['exact_gtex_pair_count'])>0 for x in subset)
   denominators.append(record)
 write('crosswalk_denominators.tsv',denominators,['phenotype_stratum','gwas_track','candidates']+all_statuses+['exact_gwas_and_gtex'])
 pair_tiers={};tier_order={'Unresolved':0,'C':1,'B':2,'A':3}
 for z in matrix:
  pair=(z['candidate_record_id'],z['gene_id'],z['gene_name'],z['phenotype_stratum'])
  if pair not in pair_tiers or tier_order[z['tier']]>tier_order[pair_tiers[pair]['tier']]:
   pair_tiers[pair]=dict(candidate_record_id=z['candidate_record_id'],gene_id=z['gene_id'],gene_name=z['gene_name'],phenotype_stratum=z['phenotype_stratum'],tier=z['tier'],tier_basis=z['tier_basis'],contexts='')
 for pair,row in pair_tiers.items():
  row['contexts']=';'.join(sorted({z['tissue_context'] for z in matrix if (z['candidate_record_id'],z['gene_id'],z['gene_name'],z['phenotype_stratum'])==pair}))
 write('candidate_gene_tier_summary.tsv',list(pair_tiers.values()),['candidate_record_id','gene_id','gene_name','phenotype_stratum','tier','tier_basis','contexts'])
 issues=[]
 for x in cw:
  if x['status']!='exact_eligible':issues.append(dict(candidate_record_id=x['candidate_record_id'],gene_name='',issue_type='crosswalk_'+x['status'],scope=x['gwas_track'],detail=x['rejection_reasons']))
 for z in matrix:
  if z['contradictions']:issues.append(dict(candidate_record_id=z['candidate_record_id'],gene_name=z['gene_name'],issue_type='gene_identity_conflict',scope=z['tissue_context'],detail=z['contradictions']))
  if z['tier']=='Unresolved':issues.append(dict(candidate_record_id=z['candidate_record_id'],gene_name=z['gene_name'],issue_type='unresolved_target',scope=z['tissue_context'],detail=z['tier_basis']))
 for k in [key(c) for c in candidates]:
  genes=sorted({z['gene_name'] for z in matrix if z['candidate_record_id']==k and z['tier'] in ('A','B','C')})
  if len(genes)>1:issues.append(dict(candidate_record_id=k,gene_name=';'.join(genes),issue_type='competing_supported_genes',scope='candidate',detail='Multiple typed suggestive or stronger links; no single target selected'))
 write('contradiction_unresolved_ledger.tsv',issues,['candidate_record_id','gene_name','issue_type','scope','detail'])
 counts=collections.Counter(x['status'] for x in cw)
 validation=dict(candidate_count=len(candidates),candidate_order_sha256=hashlib.sha256('\n'.join(key(c) for c in candidates).encode()).hexdigest(),
  exact_gwas_candidate_count=len({x['candidate_record_id'] for x in cw if x['status']=='exact_eligible'}),directed_candidate_count=len({x['candidate_record_id'] for x in directions}),
  directed_candidate_gene_pairs=len({(x['candidate_record_id'],x['gene_id_versionless']) for x in directions}),directed_rows=len(directions),
  gwas_crosswalk_status_counts=dict(counts),moloc_table5_rows=len(set(x['publication_table_row'] for x in mol)),moloc_window_parsed=len(windows),
  re2g_rows=len(re2g),saferali_additional_signal_labels=len({(x['lead_gwas_variant'],x['gene_name'],x['qtl_cohort'],x['qtl_class'],x['posterior_ppa']) for x in additional_moloc}),
  functional_exact_count=sum(x['exact_target_support']=='yes' for x in functional),tiers_context_rows=dict(collections.Counter(x['tier'] for x in matrix)),
  tiers_candidate_gene_pairs=dict(collections.Counter(x['tier'] for x in pair_tiers.values())),
  no_de_novo_colocalization=True,source_gwas_scope='prior_frozen_locus_ledger_only')
 (STAGE/'validation_results.json').write_text(json.dumps(validation,indent=2,sort_keys=True)+'\n')
 print(json.dumps(validation,indent=2))

if __name__=='__main__':main()
