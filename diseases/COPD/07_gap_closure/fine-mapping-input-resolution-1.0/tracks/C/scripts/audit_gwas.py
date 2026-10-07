from pathlib import Path
import pandas as pd,numpy as np,gzip,csv,json,hashlib,datetime,math,collections,time
from scipy.special import log_ndtr
B=Path(__file__).resolve().parents[1];R=B/'results';R.mkdir(exist_ok=True)
config=B/'config/gwas_contract.json';freeze=json.loads((B/'config/gwas_contract.freeze.json').read_text());assert hashlib.sha256(config.read_bytes()).hexdigest()==freeze['sha256']
acq=json.loads((B/'sources/acquisition.json').read_text());assert acq['integrity_pass']
LENGTHS=dict(zip(map(str,range(1,23)),[249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566]));LENGTHS['X']=155270560
FIELDS=['chr','pos','ref','alt','af_cases_EUR','af_controls_EUR','beta_EUR','se_EUR','neglog10_pval_EUR','low_confidence_EUR']
numfields=['af_cases_EUR','af_controls_EUR','beta_EUR','se_EUR','neglog10_pval_EUR']
A={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_sha256':acq['sha256'],'source_md5':acq['md5'],'source_bytes':acq['actual_bytes'],'track':'C','expected_rows_from_provider_docs':28987534,'config_sha256':freeze['sha256'],'header':None,'rows':0,'missingness':collections.Counter(),'chromosome_rows':collections.Counter(),'low_confidence_EUR':collections.Counter(),'counts':collections.Counter(),'wald_by_chrom_pband':{},'malformed_rows':0,'duplicate_exact_identity_rows':0,'multiple_allelic_rows_same_position':0,'source_position_order_violations':0,'inference_performed':False}
rowbase=0;seeds=[];lastkey=None;lastpos=None;prevord=None;sample_counts=collections.Counter();t=time.time();out=gzip.open(R/'effect_test_sample.tsv.gz','wt');sw=None
# All source columns are parsed only for syntax/missingness; numerical statistics use EUR columns exclusively.
for d in pd.read_csv(B/'data/icd10-J44-both_sexes.tsv.bgz',sep='\t',compression='gzip',dtype=str,keep_default_na=False,chunksize=100000,on_bad_lines='error'):
 if A['header'] is None:A['header']=list(d.columns)
 n=len(d);rows=np.arange(rowbase+1,rowbase+n+1,dtype=np.int64);A['rows']+=n
 for col in d.columns:A['missingness'][col]+=int(d[col].isin(['NA','','nan','NaN']).sum())
 A['chromosome_rows'].update(d['chr'].value_counts().to_dict());A['low_confidence_EUR'].update(d['low_confidence_EUR'].value_counts().to_dict())
 vals={k:pd.to_numeric(d[k],errors='coerce').to_numpy(dtype=float) for k in numfields};pos=pd.to_numeric(d['pos'],errors='coerce').to_numpy(dtype=float);chrom=d['chr'];cs=chrom.to_numpy();ref=d['ref'];alt=d['alt'];b=vals['beta_EUR'];se=vals['se_EUR'];nlp=vals['neglog10_pval_EUR'];afc=vals['af_cases_EUR'];afn=vals['af_controls_EUR']
 lim=chrom.map(LENGTHS).fillna(0).to_numpy();identity=np.isfinite(pos)&(pos>=1)&(pos<=lim)&(pos==np.floor(pos))&ref.str.fullmatch('[ACGT]+').to_numpy()&alt.str.fullmatch('[ACGT]+').to_numpy()&(ref.to_numpy()!=alt.to_numpy())
 valid=np.isfinite(b)&np.isfinite(se)&(se>0)&np.isfinite(nlp)&(nlp>=0)&np.isfinite(afc)&(afc>=0)&(afc<=1)&np.isfinite(afn)&(afn>=0)&(afn<=1)&identity
 hq=valid&(d['low_confidence_EUR'].to_numpy()=='false');auto=chrom.isin(list(LENGTHS)[:22]).to_numpy()
 A['counts'].update({'valid_signed_rows':int(valid.sum()),'valid_signed_low_confidence_false_rows':int(hq.sum()),'valid_signed_autosomal_low_confidence_false_rows':int((hq&auto).sum()),'invalid_identity_rows':int((~identity).sum()),'indel_or_multibase_rows':int(((ref.str.len()!=1)|(alt.str.len()!=1)).sum()),'palindromic_snp_rows':int(((ref+alt).isin(['AT','TA','CG','GC'])).sum())})
 key=d['chr']+':'+d['pos']+':'+ref+':'+alt;dp=d.duplicated(['chr','pos','ref','alt'],keep='first').to_numpy();dp[0]|=key.iloc[0]==lastkey;A['duplicate_exact_identity_rows']+=int(dp.sum());lastkey=key.iloc[-1]
 pp=d['chr']+':'+d['pos'];dup_pos=pp.duplicated().to_numpy();dup_pos[0]|=pp.iloc[0]==lastpos;A['multiple_allelic_rows_same_position']+=int(dup_pos.sum());lastpos=pp.iloc[-1]
 ordv=chrom.map({**{str(i):i for i in range(1,23)},'X':23}).to_numpy()*1000000000+pos;A['source_position_order_violations']+=int((np.diff(ordv)<0).sum())+int(prevord is not None and ordv[0]<prevord);prevord=ordv[-1]
 with np.errstate(divide='ignore',invalid='ignore',over='ignore'):
  z=b/se;wald=-(math.log(2)+log_ndtr(-np.abs(z)))/math.log(10)
  # Numeric exponent equals printed exponent only for normalized d.ddd e notation, which is checked on every defined value below.
  widths={}
  for k,v in [('beta_EUR',b),('se_EUR',se),('neglog10_pval_EUR',nlp)]:
   fmt=d[k].str.fullmatch(r'-?[0-9]\.[0-9]{3}e[+-][0-9]+').to_numpy();defined=np.isfinite(v);bad=defined&~fmt;A['counts'][k+'_unexpected_numeric_print_format']+=int(bad.sum())
   ex=np.floor(np.log10(np.abs(v)));q=0.0005*np.power(10.,ex);zero=(v==0)
   if zero.any():q[zero]=0.0005*np.power(10.,d.loc[zero,k].str.extract(r'e([+-][0-9]+)$')[0].astype(float).to_numpy())
   q[bad]=np.nan;widths[k]=q
  blo=b-widths['beta_EUR'];bhi=b+widths['beta_EUR'];slo=se-widths['se_EUR'];shi=se+widths['se_EUR'];azlo=np.maximum(0,np.abs(b)-widths['beta_EUR'])/shi;azhi=(np.abs(b)+widths['beta_EUR'])/slo
  wlo=-(math.log(2)+log_ndtr(-azlo))/math.log(10);whi=-(math.log(2)+log_ndtr(-azhi))/math.log(10);plo=np.maximum(0,nlp-widths['neglog10_pval_EUR']);phi=nlp+widths['neglog10_pval_EUR'];comp=valid&(wlo<=phi+1e-10)&(whi>=plo-1e-10);delta=np.abs(wald-nlp)
 A['counts'].update({'wald_rounding_compatible_valid_rows':int(comp.sum()),'wald_rounding_discrepant_valid_rows':int((valid&~comp).sum()),'wald_rounding_discrepant_hq_rows':int((hq&~comp).sum()),'abs_z_gt_40_valid_rows':int((valid&(np.abs(z)>40)).sum())})
 bands=np.select([nlp>12,nlp>-math.log10(5e-8),nlp>-math.log10(.05)],['lt_1e-12','1e-12_to_5e-8','5e-8_to_0.05'],default='ge_0.05')
 samp=(rows%1000==0)&valid;sig=hq&auto&(nlp>-math.log10(5e-8));A['counts']['genomewide_significant_seed_rows']+=int(sig.sum());samp|=valid&(nlp>-math.log10(5e-8))
 for cc in np.unique(cs):
  for band in np.unique(bands):
   sel=valid&(cs==cc)&(bands==band);k=cc+'|'+band;rec=A['wald_by_chrom_pband'].setdefault(k,{'valid':0,'hq':0,'rounding_compatible':0,'max_abs_log10_p_difference':0});rec['valid']+=int(sel.sum());rec['hq']+=int((sel&hq).sum());rec['rounding_compatible']+=int((sel&comp).sum());rec['max_abs_log10_p_difference']=max(rec['max_abs_log10_p_difference'],float(np.max(delta[sel])) if sel.any() else 0)
   need=max(0,100-sample_counts[k]);ids=np.where(sel)[0][:need];samp[ids]=True;sample_counts[k]+=len(ids)
 s=d.loc[samp,FIELDS].copy();s['source_row']=rows[samp];s['z_beta_over_se']=z[samp];s['wald_neglog10_p']=wald[samp];s['rounding_compatible']=comp[samp];s['diagnostic_status']=np.where(comp[samp],'WALD_ROUNDING_COMPATIBLE','DIAGNOSTIC_SAIGE_VS_WALD');s.to_csv(out,sep='\t',index=False,header=(rowbase==0))
 for j in np.where(sig)[0]:seeds.append({'track':'C','chrom':cs[j],'pos':int(pos[j]),'ref':ref.iloc[j],'alt':alt.iloc[j],'source_row':int(rows[j]),'neglog10_p':float(nlp[j]),'beta':float(b[j]),'se':float(se[j]),'source_neglog10_p':d['neglog10_pval_EUR'].iloc[j]})
 rowbase+=n
 if rowbase%1000000==0:print('ROWS',rowbase,'SECONDS',round(time.time()-t,1),'SEEDS',len(seeds),flush=True)
out.close();A['gzip_complete_decode_crc_pass']=True;A['row_count_matches_provider']=A['rows']==28987534;A['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
for k in ['missingness','chromosome_rows','low_confidence_EUR','counts']:A[k]=dict(A[k])
(R/'whole_file_audit.json').write_text(json.dumps(A,indent=2,allow_nan=False)+'\n')
seeddf=pd.DataFrame(seeds);seeddf.to_csv(R/'gwas_significant_seeds.tsv',sep='\t',index=False)
merged=[]
for s in sorted(seeds,key=lambda s:(int(s['chrom']),s['pos'],s['ref'],s['alt'],s['source_row'])):
 start=max(1,s['pos']-1500000);end=min(LENGTHS[s['chrom']],s['pos']+1500000)
 if merged and merged[-1]['chrom']==s['chrom'] and start<=merged[-1]['end']:
  m=merged[-1];m['end']=max(m['end'],end);m['seed_rows'].append(s)
 else:merged.append({'chrom':s['chrom'],'start':start,'end':end,'seed_rows':[s]})
records=[]
for i,m in enumerate(merged,1):
 lead=sorted(m['seed_rows'],key=lambda s:(-s['neglog10_p'],int(s['chrom']),s['pos'],s['ref'],s['alt'],s['source_row']))[0]
 mhc=m['chrom']=='6' and m['start']<=36000000 and m['end']>=25000000
 records.append({'track':'C','locus_id':f'C_L{i:03d}','chrom':m['chrom'],'start':m['start'],'end':m['end'],'build':'GRCh37','width_bp':m['end']-m['start']+1,'n_significant_seeds':len(m['seed_rows']),'lead_variant':f"{lead['chrom']}:{lead['pos']}:{lead['ref']}:{lead['alt']}",'lead_pos':lead['pos'],'lead_neglog10_p':lead['neglog10_p'],'lead_source_row':lead['source_row'],'locus_status':'DEFERRED_MHC' if mhc else 'PROSPECTIVE_REQUIRES_LD_QC','chromosome_edge_clipped':m['start']==1 or m['end']==LENGTHS[m['chrom']]})
pd.DataFrame(records).to_csv(R/'loci.tsv',sep='\t',index=False)
print('COMPLETE',A['rows'],'LOCI',len(records),'NON_MHC',sum(r['locus_status']!='DEFERRED_MHC' for r in records),flush=True)
