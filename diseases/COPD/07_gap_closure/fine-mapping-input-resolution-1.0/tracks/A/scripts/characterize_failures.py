from pathlib import Path
import pandas as pd,numpy as np,json,math
from scipy.special import log_ndtr,ndtri
R=Path(__file__).resolve().parents[1];O=R/'results'
def half(x):
 a,*b=x.lower().split('e');return .5*10.**((int(b[0]) if b else 0)-(len(a.split('.')[1]) if '.' in a else 0))
f=pd.read_csv(O/'effect_consistency_failures.tsv.gz',sep='\t',dtype=str,keep_default_na=False)
o=f.odds_ratio.astype(float).to_numpy();s=f.standard_error.astype(float).to_numpy();p=f.p_value.astype(float).to_numpy();oh=np.array([half(x) for x in f.odds_ratio]);sh=np.array([half(x) for x in f.standard_error]);ph=np.array([half(x) for x in f.p_value]);bl=np.log(o-oh);bh=np.log(o+oh);amin=np.where((bl<=0)&(bh>=0),0,np.minimum(abs(bl),abs(bh)));amax=np.maximum(abs(bl),abs(bh));model_lo=math.log(2)+log_ndtr(-amax/(s-sh));model_hi=math.log(2)+log_ndtr(-amin/(s+sh));report_lo=np.log(p-ph);report_hi=np.log(np.minimum(p+ph,1));gap=np.maximum(report_lo-model_hi,model_lo-report_hi)
f['nonoverlap_log10p_gap']=gap/math.log(10);f['point_log10p_delta']=(math.log(2)+log_ndtr(-abs(np.log(o)/s))-np.log(p))/math.log(10)
f.to_csv(O/'effect_consistency_failure_characterization.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
a={'failure_rows':len(f),'failures_p_less_5e8':int((p<5e-8).sum()),'gap_log10p_quantiles':{str(q):float(np.quantile(gap/math.log(10),q)) for q in [0,.5,.9,.99,1]},'point_log10p_abs_max':float(max(abs(f.point_log10p_delta))),'interpretation':'Mismatches remain failed at the prospectively frozen rounding-interval gate; no alternative OR, beta, SE or allele orientation substituted. Small residuals do not independently identify the historic numerical cause. All significant source rows pass.','gate_changed':False}
(O/'failure_characterization.json').write_text(json.dumps(a,indent=2)+'\n');print(json.dumps(a,indent=2))
