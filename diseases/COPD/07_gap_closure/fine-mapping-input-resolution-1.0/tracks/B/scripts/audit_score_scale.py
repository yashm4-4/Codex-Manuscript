"""Test contemporaneous SAIGE score-column branch lineage; does not alter source effects."""
from pathlib import Path
import pandas as pd,numpy as np,json,datetime
r=Path(__file__).resolve().parents[1];d=pd.read_csv(r/'results/effect_consistency_sample.tsv.gz',sep='\t');normal=d['p.value.NA']>0.05
pred=d.Tstat/d.varT/np.where(normal,1,np.sqrt(d.MAC));err=abs(d.BETA-pred);tol=1e-9+1e-6*abs(pred)
x={'method':'Source-code-derived branch expression; p.value.NA>0.05 fast branch beta=Tstat/varT; otherwise SPA branch beta=Tstat/(varT*sqrt(MAC)). No effects or rows altered.','source_snapshot':'sources/saige_step2_20180921.R','sample_rows':len(d),'normal_branch_rows':int(normal.sum()),'SPA_branch_rows':int((~normal).sum()),'failing_tolerance_rows':int((err>tol).sum()),'maximum_absolute_error':float(err.max()),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'status':'UNRESOLVED_SCORE_COLUMN_LINEAGE','interpretation':'The simple contemporaneous-source branch reconstruction is a diagnostic, not asserted exact provider implementation. Failures include variant-specific-N chr12 rows; no P/beta/SE or allele edits. Exact score-column/software branch provenance remains required for Gaussian/RSS likelihood readiness.'}
(r/'results/score_scale_sample_audit.json').write_text(json.dumps(x,indent=2)+'\n')
