from pathlib import Path
import pandas as pd,numpy as np,json,gzip,time,datetime
B=Path(__file__).resolve().parents[1];R=B/'results';F=['chr','pos','ref','alt','af_cases_EUR','af_controls_EUR','beta_EUR','se_EUR','neglog10_pval_EUR','low_confidence_EUR'];N=F[4:9];r={c:{'defined':0,'missing':0,'nonfinite_nonmissing':0,'min':float('inf'),'max':-float('inf'),'invalid_range':0} for c in N};row=0;badrows=[]
for d in pd.read_csv(B/'data/icd10-J44-both_sexes.tsv.bgz',compression='gzip',sep='\t',usecols=F,dtype=str,keep_default_na=False,chunksize=250000):
 bad=np.zeros(len(d),dtype=bool)
 for c in N:
  v=pd.to_numeric(d[c],errors='coerce');missing=d[c].isin(['NA','','nan','NaN']);defined=np.isfinite(v);invalid=(~defined&~missing)|((v<0)|(v>1) if c.startswith('af_') else (v<=0) if c=='se_EUR' else (v<0) if c=='neglog10_pval_EUR' else False);bad|=invalid.to_numpy();q=r[c];q['defined']+=int(defined.sum());q['missing']+=int(missing.sum());q['nonfinite_nonmissing']+=int((~defined&~missing).sum());q['invalid_range']+=int(invalid.sum());q['min']=min(q['min'],float(v.min()));q['max']=max(q['max'],float(v.max()))
 if bad.any():
  d['source_row']=np.arange(row+1,row+len(d)+1);badrows+=d.loc[bad].to_dict('records')
 row+=len(d)
for q in r.values():
 for k,v in list(q.items()):
  if isinstance(v,float) and not np.isfinite(v):q[k]='Infinity' if v>0 else '-Infinity' if v<0 else 'NaN'
pd.DataFrame(badrows).to_csv(R/'whole_file_invalid_numeric_rows.tsv',sep='\t',index=False);(R/'numeric_range_audit.json').write_text(json.dumps({'rows':row,'full_gzip_decode_crc_pass':True,'fields':r,'nonmissing_invalid_rows':len(badrows),'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2)+'\n');print('COMPLETE',row,'INVALID',len(badrows),flush=True)
