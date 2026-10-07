"""Whole combined-sex autosomal source audit; no inference or LD construction."""
from pathlib import Path
from collections import Counter,defaultdict
import pandas as pd,numpy as np,json,gzip,hashlib,datetime,csv,subprocess,os
from scipy.special import log_ndtr
r=Path(__file__).resolve().parents[1];out=r/'results';cfg=json.loads((r/'config/tolerances.json').read_text());source=r/'acquisition/COPD.auto.rsq07.mac10.txt.gz'
a={'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'track':'B','accession':'GCST90013709','source':'acquisition/COPD.auto.rsq07.mac10.txt.gz','tolerance_sha256':hashlib.sha256((r/'config/tolerances.json').read_bytes()).hexdigest(),'rows':0,'chromosomes':Counter(),'missingness':Counter(),'qc_failures':Counter(),'ranges':{},'p_consistency':{},'counts':Counter(),'chunk_log':[]}
counts=Counter();sig=[];examples=[];prev=None;prevchr=None; seenchr=set(); perchr_ids=set();perchr_pos=defaultdict(set);duplicate_rows=[];multipos=[]; ranges={}; metrics=defaultdict(lambda:{'n':0,'fail':0,'max_error':0.0,'max_error_row':None})
num=['CHR','POS','AC_Allele2','AF_Allele2','N','BETA','SE','Tstat','p.value','p.value.NA','Is.SPA.converge','varT','varTstar','AF.Cases','AF.Controls','Rsq','MAC']
# Exact identity strings are stored per chromosome; release order is independently checked.
def vals(d):
 v={k:pd.to_numeric(d[k],errors='coerce').to_numpy() for k in num}
 return v

def checks(d,v):
 f={}
 for k in num:f['invalid_'+k]=~np.isfinite(v[k])
 f['nonpositive_SE']=v['SE']<=0;f['bad_p']=(v['p.value']<0)|(v['p.value']>1);f['bad_p_NA']=(v['p.value.NA']<0)|(v['p.value.NA']>1)
 f['N_unexpected']=v['N']!=cfg['n_expected'];f['invalid_AF']=(v['AF_Allele2']<0)|(v['AF_Allele2']>1);f['Rsq_below_07']=v['Rsq']<.7;f['MAC_below_10']=v['MAC']<10
 f['SPA_not_converged']=v['Is.SPA.converge']!=1
 f['allele_sequence_invalid']=~(d.Allele1.str.fullmatch('[ACGT]+') & d.Allele2.str.fullmatch('[ACGT]+')).to_numpy();f['same_alleles']=(d.Allele1==d.Allele2).to_numpy()
 expected=d.CHR+'_'+d.POS+'_'+d.Allele1+'_'+d.Allele2;f['SNPID_mismatch']=(expected!=d.SNPID).to_numpy()
 f['position_outside_chr']=np.array([not (c in cfg['chromosome_lengths_GRCh37'] and 1<=p<=cfg['chromosome_lengths_GRCh37'][c]) for c,p in zip(d.CHR,v['POS'])])
 return f

for idx,d in enumerate(pd.read_csv(source,sep=r'\s+',dtype=str,keep_default_na=False,chunksize=200000)):
 if idx==0:a['schema']=list(d.columns); assert len(d.columns)==20
 rownums=np.arange(a['rows']+1,a['rows']+len(d)+1);v=vals(d);f=checks(d,v)
 for k in d.columns:a['missingness'][k]+=int(d[k].isin(['','NA','NaN','nan','.']).sum())
 for k,x in f.items():a['qc_failures'][k]+=int(x.sum())
 for k,x in v.items():
  finite=x[np.isfinite(x)]
  if len(finite):
   old=ranges.get(k,[float('inf'),float('-inf')]);ranges[k]=[min(old[0],float(finite.min())),max(old[1],float(finite.max()))]
 a['chromosomes'].update(d.CHR)
 allgood=np.ones(len(d),dtype=bool)
 for x in f.values():allgood&=~x
 counts['rows_passing_source_stat_qc']+=int(allgood.sum());counts['p_zero']+=int((v['p.value']==0).sum());counts['p_NA_zero']+=int((v['p.value.NA']==0).sum())
 lengths=(d.Allele1.str.len()==1)&(d.Allele2.str.len()==1);counts['SNP']+=int(lengths.sum());counts['indel_or_MNV']+=int((~lengths).sum());counts['palindromic_SNP']+=int((lengths & ((d.Allele1+d.Allele2).isin(['AT','TA','CG','GC']))).sum())
 z=v['BETA']/v['SE'];logp=(np.log(2)+log_ndtr(-np.abs(z)))/np.log(10)
 for pk in ['p.value','p.value.NA']:
  good=np.isfinite(z)&(v[pk]>0)&(v[pk]<=1);errs=np.abs(logp[good]-np.log10(v[pk][good]));limit=cfg['wald_log10_p_absolute_tolerance']+cfg['wald_log10_p_relative_tolerance']*np.abs(np.log10(v[pk][good]));m=metrics[pk];m['n']+=len(errs);m['fail']+=int((errs>limit).sum())
  if len(errs) and float(errs.max())>m['max_error']: ix=np.flatnonzero(good)[errs.argmax()];m['max_error']=float(errs.max());m['max_error_row']=int(rownums[ix]);m['max_error_record']=d.iloc[ix].to_dict()
 for name,lhs,rhs,atol,rtol in [('AC_equals_2N_AF',v['AC_Allele2'],2*v['N']*v['AF_Allele2'],.001,1e-6),('BETA_equals_Tstat_over_varT',v['BETA'],v['Tstat']/v['varT'],1e-9,1e-6)]:
  good=np.isfinite(lhs)&np.isfinite(rhs);errs=np.abs(lhs[good]-rhs[good]);m=metrics[name];m['n']+=len(errs);m['fail']+=int((errs>(atol+rtol*np.abs(rhs[good]))).sum());m['max_error']=max(m['max_error'],float(errs.max(initial=0)))
 # Deterministic source-row sample plus all significant signals spans all chromosomes/P ranges.
 mask=(rownums%1000==1)|(v['p.value']<5e-8)
 for i in np.flatnonzero(mask):
  ex=d.iloc[i].to_dict();ex.update(source_row=int(rownums[i]),z=float(z[i]),wald_neglog10p=float(-logp[i]),source_stat_qc='PASS' if allgood[i] else ';'.join(k for k,x in f.items() if x[i]));examples.append(ex)
 for i in np.flatnonzero((v['p.value']<5e-8)&allgood):
  s=d.iloc[i].to_dict();s['source_row']=int(rownums[i]);sig.append(s)
 # Exact duplicates and multiallelic-coordinate representation; no hash-only identity joins.
 for c,p,k,a1,a2,sr in zip(d.CHR,d.POS,d.SNPID,d.Allele1,d.Allele2,rownums):
  coord=(int(c),int(p))
  if prev is not None and coord<prev:counts['coordinate_order_decreases']+=1
  if c!=prevchr:
   if c in seenchr: raise AssertionError('Chromosome repeated non-contiguously: exact duplicate audit requires disk-sort fallback')
   if prevchr is not None:multipos.extend({'chrom':prevchr,'pos':pp,'allelic_identities':';'.join(sorted(ks))} for pp,ks in perchr_pos.items() if len(ks)>1)
   seenchr.add(c);perchr_ids=set();perchr_pos=defaultdict(set);prevchr=c
  identity=(p,a1,a2)
  if identity in perchr_ids:duplicate_rows.append({'source_row':int(sr),'chrom':c,'pos':p,'ref':a1,'alt':a2,'SNPID':k})
  perchr_ids.add(identity);perchr_pos[p].add(k);prev=coord
 a['rows']+=len(d);a['chunk_log'].append({'rows':a['rows'],'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()});print('audited',a['rows'],flush=True)
multipos.extend({'chrom':prevchr,'pos':pp,'allelic_identities':';'.join(sorted(ks))} for pp,ks in perchr_pos.items() if len(ks)>1)
a['counts']=dict(counts);a['ranges']=ranges;a['p_consistency']=dict(metrics);a['duplicates_extra_rows']=len(duplicate_rows);a['multiallelic_positions']=len(multipos);a['significant_valid_rows']=len(sig);a['gzip_end_to_end_crc_pass']=True
pd.DataFrame(examples).to_csv(out/'effect_consistency_sample.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
pd.DataFrame(duplicate_rows,columns=['source_row','chrom','pos','ref','alt','SNPID']).to_csv(out/'duplicate_identities.tsv',sep='\t',index=False)
pd.DataFrame(multipos,columns=['chrom','pos','allelic_identities']).to_csv(out/'multiallelic_positions.tsv',sep='\t',index=False)
pd.DataFrame(sig,columns=a['schema']+['source_row']).to_csv(out/'significant_variants.tsv',sep='\t',index=False)
# Deterministic transitive window merge uses source P, chromosome, position, exact allele key, row.
sig=sorted(sig,key=lambda s:(int(s['CHR']),int(s['POS']),s['SNPID'],s['source_row']))
loci=[]
for s in sig:
 c=int(s['CHR']);p=int(s['POS']);start=max(1,p-1500000);end=min(cfg['chromosome_lengths_GRCh37'][str(c)],p+1500000)
 if loci and loci[-1]['chrom']==c and start<=loci[-1]['end']:
  loci[-1]['end']=max(loci[-1]['end'],end);loci[-1]['signals'].append(s)
 else:loci.append({'track':'B','chrom':c,'start':start,'end':end,'signals':[s]})
for l in loci:
 lead=min(l['signals'],key=lambda s:(float(s['p.value']),int(s['CHR']),int(s['POS']),s['SNPID'],s['source_row']));l.update(locus_id=f"B_GRCh37_chr{l['chrom']}_{l['start']}_{l['end']}",lead_variant=lead['SNPID'],lead_p=lead['p.value'],lead_source_row=lead['source_row'],significant_rows=len(l['signals']),mhc_deferred=l['chrom']==6 and l['start']<=36000000 and l['end']>=25000000,chromosome_edge_clipped=l['start']==1 or l['end']==cfg['chromosome_lengths_GRCh37'][str(l['chrom'])],summary_stat_gate='PENDING_REFERENCE_IDENTITY_AND_TEST_LIKELIHOOD',ld_gate='NOT_RUN_NO_DEFENSIBLE_DENSE_SIGNED_JAPANESE_LD',readiness='NOT CLEARED');del l['signals']
pd.DataFrame(loci,columns=['track','locus_id','chrom','start','end','lead_variant','lead_p','lead_source_row','significant_rows','mhc_deferred','chromosome_edge_clipped','summary_stat_gate','ld_gate','readiness']).to_csv(out/'loci.tsv',sep='\t',index=False)
a['prospective_merged_loci_before_MHC']=len(loci);a['prospective_loci']=sum(not l['mhc_deferred'] for l in loci);a['mhc_deferred_loci']=sum(l['mhc_deferred'] for l in loci)
# Second full streaming pass retains every row within merged intervals incl nonsignificant, without LD filtering.
col=['track','locus_id','chrom','pos','ref','alt','effect_allele','other_allele','beta','se','z','p','neglog10_p','source_row','rsid','af','n','qc_status','harmonization_status','variant_identity']+['original_'+k for k in a['schema']]
rows_by_locus=Counter();offset=0
with (out/'locus_variants.tsv.gz').open('wb') as raw:
 with gzip.GzipFile(fileobj=raw,mode='wb',mtime=0) as gz:
  import io
  fo=io.TextIOWrapper(gz);w=csv.DictWriter(fo,fieldnames=col,delimiter='\t',lineterminator='\n');w.writeheader()
  for d in pd.read_csv(source,sep=r'\s+',dtype=str,keep_default_na=False,chunksize=200000):
   v=vals(d);f=checks(d,v)
   for l in loci:
    hit=(v['CHR']==l['chrom'])&(v['POS']>=l['start'])&(v['POS']<=l['end'])
    for i in np.flatnonzero(hit):
     s=d.iloc[i];fails=[k for k,x in f.items() if x[i]];x=dict(track='B',locus_id=l['locus_id'],chrom=int(v['CHR'][i]),pos=int(v['POS'][i]),ref=s.Allele1,alt=s.Allele2,effect_allele=s.Allele2,other_allele=s.Allele1,beta=s.BETA,se=s.SE,z=v['BETA'][i]/v['SE'][i],p=s['p.value'],neglog10_p=-np.log10(v['p.value'][i]) if v['p.value'][i]>0 else 'inf',source_row=offset+i+1,rsid='',af=s.AF_Allele2,n=s.N,qc_status='PASS' if not fails else ';'.join(fails),harmonization_status='SOURCE_GRCh37_REF_ALT_VERIFIED_REFERENCE_SEQUENCE_PENDING',variant_identity=f'GRCh37:{s.CHR}:{s.POS}:{s.Allele1}:{s.Allele2}')
     x.update({'original_'+k:s[k] for k in a['schema']});w.writerow(x);rows_by_locus[l['locus_id']]+=1
   offset+=len(d)
  fo.flush()
a['locus_variant_rows']=dict(rows_by_locus);a['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();a['readiness']='NOT CLEARED';a['execution_cleared_loci']=0
(out/'whole_file_audit.json').write_text(json.dumps(a,indent=2)+'\n');print(json.dumps(a,indent=2),flush=True)
