from pathlib import Path
import pandas as pd,numpy as np,json,collections,math,datetime,itertools
from scipy.special import log_ndtr
R=Path(__file__).resolve().parents[1];S=R/'sources';O=R/'results';C=json.loads((R/'config/effect_contract.json').read_text());T=json.loads((R/'config/formatted_transport_contract.json').read_text())
path=S/'harmonised/33106845-GCST90016588-EFO_0006527-Build37.f.tsv.gz';m=json.loads(Path(str(path)+'.access.json').read_text());assert m['provider_md5_pass'] and m['complete_http']
def half(x):
 a,*b=x.lower().split('e');return .5*10.**((int(b[0]) if b else 0)-(len(a.split('.')[1]) if '.' in a else 0))
kwargs={'sep':'\t','dtype':str,'keep_default_na':False,'chunksize':200000};nr=pd.read_csv(S/'GCST90016588_buildGRCh37.tsv',**kwargs);fr=pd.read_csv(path,**kwargs)
c=collections.Counter();missing=collections.Counter();fails=[];diff_examples=[];locusrows=set(pd.read_csv(O/'locus_variants.tsv.gz',sep='\t',usecols=['source_row']).source_row);locus=[];total=0
for d,f in itertools.zip_longest(nr,fr):
 assert d is not None and f is not None and len(d)==len(f),'row count mismatch';n=len(d);f['source_row']=np.arange(total+1,total+n+1);total+=n
 for col in f.columns[:-1]:missing[col]+=int(f[col].str.strip('"').isin(['','NA','NaN','.','nan']).sum())
 ids=['variant_id','chromosome','base_pair_location','effect_allele','other_allele'];idgood=np.ones(n,bool)
 for col in ids:
  eq=d[col].str.strip('"').to_numpy()==f[col].str.strip('"').to_numpy();c['identity_'+col+'_differences']+=int((~eq).sum());idgood&=eq
 c['row_identity_exact']+=int(idgood.sum());transport=idgood.copy()
 for col in ['p_value','odds_ratio','standard_error']:
  x=d[col].astype(float).to_numpy();y=f[col].astype(float).to_numpy();eq=x==y;close=abs(x-y)<=T['maximum_ulp_for_serialization_equivalence']*np.spacing(np.maximum(abs(x),abs(y)));c[col+'_string_differences']+=int((d[col].to_numpy()!=f[col].to_numpy()).sum());c[col+'_float64_exact']+=int(eq.sum());c[col+'_within_4ulp']+=int(close.sum());c[col+'_larger_than_4ulp']+=int((~close).sum());transport&=close
  diff=np.where(d[col].to_numpy()!=f[col].to_numpy())[0]
  for k in diff[:5]:diff_examples.append({'source_row':int(f.source_row.iloc[k]),'field':col,'native_string':d[col].iloc[k],'formatted_string':f[col].iloc[k],'float64_difference':float(x[k]-y[k]),'within_4ulp':bool(close[k])})
 c['transport_pass']+=int(transport.sum());p=f.p_value.astype(float).to_numpy();o=f.odds_ratio.astype(float).to_numpy();s=f.standard_error.astype(float).to_numpy();oh=np.array([half(x) for x in f.odds_ratio]);sh=np.array([half(x) for x in f.standard_error]);ph=np.array([half(x) for x in f.p_value]);bl=np.log(o-oh);bh=np.log(o+oh);amin=np.where((bl<=0)&(bh>=0),0,np.minimum(abs(bl),abs(bh)));amax=np.maximum(abs(bl),abs(bh));ml=math.log(2)+log_ndtr(-amax/(s-sh));mh=math.log(2)+log_ndtr(-amin/(s+sh));rl=np.log(np.maximum(p-ph,np.finfo(float).tiny));rh=np.log(np.minimum(p+ph,1));good=(ml<=rh+C['log_p_numerical_slack'])&(rl<=mh+C['log_p_numerical_slack'])&np.isfinite(o)&(o>0)&np.isfinite(s)&(s>0)&np.isfinite(p)&(p>0)&(p<=1)
 f['formatted_rounding_pass']=good;f['native_formatted_transport_pass']=transport;f['beta_logOR']=np.log(o);f['z']=np.log(o)/s;c['effect_pass']+=int(good.sum());c['effect_fail']+=int((~good).sum());c['significant_effect_fail']+=int(((p<5e-8)&~good).sum());fails.append(f.loc[~good]);locus.append(f.loc[f.source_row.isin(locusrows)]);print('formatted rows',total,flush=True)
pd.concat(fails,ignore_index=True).to_csv(O/'formatted_effect_failures.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0});pd.concat(locus,ignore_index=True).to_csv(O/'formatted_locus_source_fields.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0});pd.DataFrame(diff_examples).to_csv(O/'native_formatted_serialization_examples.tsv',sep='\t',index=False)
a={'track':'A','formatted_rows':total,'gzip_full_decompression_crc':'PASS','source_md5_pass':True,'sha256':m['sha256'],'counts':dict(c),'missingness':dict(missing),'gate_unchanged':True,'native_failures_retained':True,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()};(O/'formatted_audit.json').write_text(json.dumps(a,indent=2)+'\n');print(json.dumps(a,indent=2),flush=True)
