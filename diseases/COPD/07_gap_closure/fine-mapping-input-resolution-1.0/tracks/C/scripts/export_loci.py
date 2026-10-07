from pathlib import Path
import pandas as pd,numpy as np,gzip,json,time,collections,math
B=Path(__file__).resolve().parents[1];R=B/'results'
while not (R/'loci.tsv').exists():time.sleep(5)
lo=pd.read_csv(R/'loci.tsv',sep='\t',dtype={'chrom':str});audit=json.loads((R/'whole_file_audit.json').read_text());assert audit['gzip_complete_decode_crc_pass'] and audit['row_count_matches_provider']
F=['chr','pos','ref','alt','af_cases_EUR','af_controls_EUR','beta_EUR','se_EUR','neglog10_pval_EUR','low_confidence_EUR'];rows=0;out=gzip.open(R/'locus_variants.source.tsv.gz','wt');first=True;counts=collections.Counter()
for d in pd.read_csv(B/'data/icd10-J44-both_sexes.tsv.bgz',compression='gzip',sep='\t',usecols=F,dtype=str,keep_default_na=False,chunksize=250000):
 p=pd.to_numeric(d['pos'],errors='coerce');sr=np.arange(rows+1,rows+len(d)+1);d['source_row']=sr
 for l in lo.to_dict('records'):
  sel=(d['chr']==l['chrom'])&(p>=l['start'])&(p<=l['end']);s=d.loc[sel].copy()
  if s.empty:continue
  o=pd.DataFrame(index=s.index);o['track']='C';o['locus_id']=l['locus_id'];o['chrom']=s['chr'];o['pos']=s['pos'];o['ref']=s['ref'];o['alt']=s['alt'];o['effect_allele']=s['alt'];o['other_allele']=s['ref']
  b=pd.to_numeric(s['beta_EUR'],errors='coerce');se=pd.to_numeric(s['se_EUR'],errors='coerce');nlp=pd.to_numeric(s['neglog10_pval_EUR'],errors='coerce');afc=pd.to_numeric(s['af_cases_EUR'],errors='coerce');afn=pd.to_numeric(s['af_controls_EUR'],errors='coerce')
  o['beta']=b;o['se']=se;o['z']=b/se;o['p']=10.**(-nlp);o['neglog10_p']=nlp;o['source_row']=s['source_row'];o['rsid']='';o['af']=(afc*11536+afn*408995)/420531;o['n']=420531
  valid=np.isfinite(b)&np.isfinite(se)&(se>0)&np.isfinite(nlp)&(nlp>=0)&np.isfinite(afc)&(afc>=0)&(afc<=1)&np.isfinite(afn)&(afn>=0)&(afn<=1)&s['ref'].str.fullmatch('[ACGT]+')&s['alt'].str.fullmatch('[ACGT]+')&(s['ref']!=s['alt'])
  o['qc_status']=np.select([~valid,s['low_confidence_EUR']=='true',s['low_confidence_EUR']!='false'],['UNRESOLVED_MISSING_OR_INVALID_GWAS','SOURCE_LOW_CONFIDENCE','UNRESOLVED_SOURCE_QC'],default='SOURCE_STATISTIC_QC_PASS');o['harmonization_status']='NATIVE_GRCH37_FORWARD_ALT_CONTRAST_REFERENCE_CHECK_PENDING';o['variant_id']='GRCh37:'+s['chr']+':'+s['pos']+':'+s['ref']+':'+s['alt'];o['locus_status']=l['locus_status'];o['source_identity']='GRCh37:'+s['chr']+':'+s['pos']+':'+s['ref']+':'+s['alt'];o['n_status']='PHENOTYPE_LEVEL_CONSTANT_PER_VARIANT_N_ABSENT';o['af_status']='CASE_CONTROL_WEIGHTED_AF_APPROXIMATION_FROM_ROUNDED_FIELDS'
  for c in F:o['original_'+c]=s[c]
  o.to_csv(out,sep='\t',index=False,header=first,na_rep='NA');first=False;counts[l['locus_id']]+=len(o)
 rows+=len(d)
 if rows%2000000==0:print('ROWS',rows,'EXPORTED',sum(counts.values()),flush=True)
out.close();assert rows==audit['rows'];(R/'locus_export.json').write_text(json.dumps({'complete_source_rows_revisited':rows,'locus_rows':dict(counts),'total_locus_rows':sum(counts.values()),'includes_nonsignificant_missing_invalid_and_deferred_rows':True},indent=2)+'\n');print('COMPLETE',dict(counts),flush=True)
