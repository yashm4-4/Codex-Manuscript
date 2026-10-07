from pathlib import Path
import pandas as pd,numpy as np,json,math,hashlib,datetime,collections,gzip
from scipy.special import log_ndtr
R=Path(__file__).resolve().parents[1];O=R/'results';S=R/'sources';C=json.loads((R/'config/effect_contract.json').read_text())
def write(df,name):df.to_csv(O/name,sep='\t',index=False,compression={'method':'gzip','mtime':0} if name.endswith('.gz') else None)
def halfunit(x):
 x=x.lower();a,*b=x.split('e');ex=int(b[0]) if b else 0;places=len(a.split('.')[1]) if '.' in a else 0
 return 0.5*10.**(ex-places)
def logp(z):return math.log(2)+log_ndtr(-np.abs(z))
meta=json.loads((S/'GCST90016588_buildGRCh37.tsv.access.json').read_text());assert meta['complete_http'] and meta['status']==200
provider=[l.split()[0] for l in (S/'md5sum.txt').read_text().splitlines() if l.rstrip().endswith('GCST90016588_buildGRCh37.tsv')][0]
assert meta['md5']==provider
print('load',datetime.datetime.now().isoformat(),flush=True)
d=pd.read_csv(S/'GCST90016588_buildGRCh37.tsv',sep='\t',dtype=str,keep_default_na=False);N=len(d);d['source_row']=np.arange(1,N+1)
a={'track':'A','source':'GCST90016588','rows':N,'native_build':'GRCh37','source_md5_pass':True,'native_sha256':meta['sha256'],'compression':'none','source_columns':list(d.columns[:-1]),'missingness':{k:int(d[k].str.strip('"').isin(['','NA','NaN','.','nan']).sum()) for k in d.columns[:-1]},'chromosomes':dict(d.chromosome.value_counts().sort_index()),'af_fields':[],'n_fields':[],'imputation_qc_fields':[],'nominal_n':71591,'n_cases':12446,'n_controls':59145,'per_variant_N_status':'absent; no assumption that nominal N is per-variant observed N','provider_qc':'Publication states MAF >= 0.01 and imputation r2 >= 0.5; fields absent from file; row-level recomputation unavailable'}
a['chromosomes']={str(k):int(v) for k,v in a['chromosomes'].items()}
for col in ['p_value','odds_ratio','standard_error','base_pair_location']:
 d[col+'_numeric']=pd.to_numeric(d[col],errors='coerce')
p=d.p_value_numeric.to_numpy();o=d.odds_ratio_numeric.to_numpy();s=d.standard_error_numeric.to_numpy();pos=d.base_pair_location_numeric.to_numpy()
valid=np.isfinite(p)&(p>0)&(p<=1)&np.isfinite(o)&(o>0)&np.isfinite(s)&(s>0)&np.isfinite(pos)&(pos>=1)&(pos==np.floor(pos))
beta=np.log(o);z=beta/s
print('rounding intervals',flush=True)
oh=np.fromiter((halfunit(x) for x in d.odds_ratio),float,count=N);sh=np.fromiter((halfunit(x) for x in d.standard_error),float,count=N);ph=np.fromiter((halfunit(x) for x in d.p_value),float,count=N)
bl=np.log(np.maximum(o-oh,np.finfo(float).tiny));bh=np.log(o+oh);small=np.where((bl<=0)&(bh>=0),0,np.minimum(abs(bl),abs(bh)));large=np.maximum(abs(bl),abs(bh))
zlo=small/(s+sh);zhi=large/np.maximum(s-sh,np.finfo(float).tiny)
log_model_lo=logp(zhi);log_model_hi=logp(zlo);log_report_lo=np.log(np.maximum(p-ph,np.finfo(float).tiny));log_report_hi=np.log(np.minimum(p+ph,1))
passed=valid&(log_model_lo<=log_report_hi+C['log_p_numerical_slack'])&(log_report_lo<=log_model_hi+C['log_p_numerical_slack'])
point_delta=(logp(z)-np.log(p))/math.log(10)
d['beta_logOR']=beta;d['z']=z;d['wald_p_rounding_pass']=passed;d['rounding_direction_unresolved']=(o-oh<=1)&(o+oh>=1)
d['p_band']=np.select([p<1e-12,p<5e-8,p<.05],['lt_1e-12','1e-12_to_5e-8','5e-8_to_0.05'],default='ge_0.05')
a['statistic_valid_rows']=int(valid.sum());a['rounding_consistency_pass']=int(passed.sum());a['rounding_consistency_fail']=int((~passed).sum());a['p_zero_rows']=int((p==0).sum());a['point_log10p_abs_difference_max']=float(np.nanmax(abs(point_delta)));a['rounding_direction_unresolved_rows']=int(d.rounding_direction_unresolved.sum())
a['beta_range']=[float(np.nanmin(beta)),float(np.nanmax(beta))];a['z_range']=[float(np.nanmin(z)),float(np.nanmax(z))];a['p_min']=float(np.nanmin(p));a['p_max']=float(np.nanmax(p))
sel=(d.groupby(['chromosome','p_band'],sort=False).cumcount()<250)|(d.source_row%1000==0)|(p<5e-8)
sample=d.loc[sel].copy();sample['point_log10p_delta']=point_delta[sel];sample['model_log10p_low']=log_model_lo[sel]/math.log(10);sample['model_log10p_high']=log_model_hi[sel]/math.log(10);sample['reported_log10p_low']=log_report_lo[sel]/math.log(10);sample['reported_log10p_high']=log_report_hi[sel]/math.log(10)
write(sample,'effect_consistency_deterministic_sample.tsv.gz');write(d.loc[~passed],'effect_consistency_failures.tsv.gz');a['deterministic_sample_rows']=len(sample)
print('duplicates',flush=True)
identity=['chromosome','base_pair_location','effect_allele','other_allele'];dups=d.duplicated(identity,keep=False);a['directional_identity_duplicate_rows']=int(dups.sum());write(d.loc[dups],'duplicate_identity_records.tsv.gz')
a['exact_duplicate_extra_rows']=int(d.duplicated(list(d.columns[:8])).sum())
# Unordered allele pairs distinguish true multiallelic sites from duplicate contrasts.
x=np.minimum(d.effect_allele.to_numpy(dtype=str),d.other_allele.to_numpy(dtype=str)) if False else None
pairs=pd.DataFrame({'chrom':d.chromosome,'pos':d.base_pair_location,'ea':d.effect_allele,'oa':d.other_allele})
pairs['lo']=[min(x,y) for x,y in zip(d.effect_allele,d.other_allele)];pairs['hi']=[max(x,y) for x,y in zip(d.effect_allele,d.other_allele)]
up=pairs[['chrom','pos','lo','hi']].drop_duplicates();multi=up.duplicated(['chrom','pos'],keep=False)
a['multiallelic_positions']=int(up.loc[multi,['chrom','pos']].drop_duplicates().shape[0]);write(up.loc[multi],'multiallelic_positions.tsv.gz')
a['unordered_identity_duplicate_extra_rows']=N-len(up)
is_snp=(d.effect_allele.str.len()==1)&(d.other_allele.str.len()==1)&d.effect_allele.str.fullmatch('[ACGT]')&d.other_allele.str.fullmatch('[ACGT]')
pal=(d.effect_allele+d.other_allele).isin(['AT','TA','CG','GC']);d['variant_type']=np.where(is_snp,np.where(pal,'palindromic_SNP','nonpalindromic_SNP'),'indel_or_complex');a['variant_types']={k:int(v) for k,v in d.variant_type.value_counts().items()};a['identical_alleles']=int((d.effect_allele==d.other_allele).sum())
print('loci',flush=True)
lengths={l.split()[0].replace('chr',''):int(l.split()[1]) for l in (S/'hg19.chrom.sizes').read_text().splitlines()}
auto=d.chromosome.isin(list(map(str,range(1,23))));sig=d.loc[(p<5e-8)&valid&auto].copy();sig['chrom_numeric']=sig.chromosome.astype(int);sig['allele_key']=sig.effect_allele+':'+sig.other_allele
sig.sort_values(['chrom_numeric','base_pair_location_numeric','p_value_numeric','allele_key','source_row'],inplace=True)
loci=[]
for chrom,g in sig.groupby('chromosome',sort=False):
 windows=[]
 for r in g.to_dict('records'):
  st=max(1,int(r['base_pair_location_numeric'])-1500000);en=min(lengths[chrom],int(r['base_pair_location_numeric'])+1500000)
  if windows and st<=windows[-1]['end']:windows[-1]['end']=max(windows[-1]['end'],en);windows[-1]['rows'].append(r)
  else:windows.append({'start':st,'end':en,'rows':[r]})
 for w in windows:
  lead=min(w['rows'],key=lambda r:(r['p_value_numeric'],r['chrom_numeric'],r['base_pair_location_numeric'],r['allele_key'],r['source_row']))
  status='DEFERRED_MHC' if chrom=='6' and w['start']<=36000000 and w['end']>=25000000 else 'PROSPECTIVE'
  loci.append({'track':'A','locus_id':f'A_chr{chrom}_{w["start"]}_{w["end"]}','build':'GRCh37','chrom':chrom,'start':w['start'],'end':w['end'],'lead_pos':int(lead['base_pair_location_numeric']),'lead_p':lead['p_value'],'lead_effect_allele':lead['effect_allele'],'lead_other_allele':lead['other_allele'],'lead_source_row':lead['source_row'],'lead_rsid':lead['variant_id'],'significant_variants':len(w['rows']),'status':status,'edge_clipped':w['start']==1 or w['end']==lengths[chrom],'boundary_rule':'inclusive; all windows +/-1.5Mb transitively merged','lead_rounding_consistency_pass':lead['wald_p_rounding_pass']})
locusdf=pd.DataFrame(loci);write(locusdf,'loci.tsv');write(sig,'significant_variants.tsv.gz');a['significant_autosomal_rows']=len(sig);a['prospective_loci_including_mhc']=len(loci);a['prospective_loci_non_mhc']=sum(l['status']!='DEFERRED_MHC' for l in loci);a['deferred_mhc_loci']=sum(l['status']=='DEFERRED_MHC' for l in loci)
print(locusdf.to_string(index=False),flush=True)
parts=[]
for loc in loci:
 q=d.loc[(d.chromosome==loc['chrom'])&(pos>=loc['start'])&(pos<=loc['end'])].copy();q['locus_id']=loc['locus_id'];q['locus_status']=loc['status'];q['native_duplicate_identity']=dups.loc[q.index];parts.append(q)
v=pd.concat(parts,ignore_index=True);v['track']='A';v['chrom']=v.chromosome;v['pos']=v.base_pair_location;v['ref']='UNRESOLVED';v['alt']='UNRESOLVED';v['beta']=v.beta_logOR;v['se']=v.standard_error;v['p']=v.p_value;v['neglog10_p']=-np.log10(v.p_value_numeric);v['rsid']=v.variant_id;v['af']='NOT_REPORTED';v['n']='NOT_REPORTED';v['qc_status']=np.where(v.wald_p_rounding_pass&~v.native_duplicate_identity,'SIGNED_STATISTIC_PASS','SIGNED_STATISTIC_OR_DUPLICATE_FAIL');v['harmonization_status']='NATIVE_GRCH37_REFERENCE_NOT_YET_VERIFIED';v['variant_identity']=v.chrom+':'+v.pos+':EA='+v.effect_allele+':OA='+v.other_allele
write(v,'locus_variants.tsv.gz')
risk=v[['track','locus_id','variant_identity','effect_allele','other_allele','beta','odds_ratio','se','z','p','source_row','variant_type','harmonization_status']].copy();risk['disease_increasing_allele']='UNRESOLVED';risk['direction_status']=np.where(v.wald_p_rounding_pass&~v.rounding_direction_unresolved,'SIGNED_SOURCE_CONTRAST_VERIFIED_IDENTITY_PENDING','UNRESOLVED_SIGN_OR_STATISTIC');write(risk,'risk_direction.tsv.gz')
a['prospective_locus_variant_rows']=len(v);a['prospective_locus_statistic_fail_rows']=int((~v.wald_p_rounding_pass).sum());a['prospective_locus_duplicate_rows']=int(v.native_duplicate_identity.sum());a['all_summary_gates_loci']=0;a['all_summary_gates_note']='reference identity/normalization, per-variant N contract and LD matching remain pending at native audit';a['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
(O/'native_audit.json').write_text(json.dumps(a,indent=2)+'\n');print(json.dumps(a,indent=2),flush=True)
