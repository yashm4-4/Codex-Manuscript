"""Exact native allele-key joins for diagnostic inputs; no posterior inference."""
from pathlib import Path
import json
import numpy as np
import pandas as pd

s=Path(__file__).resolve().parents[1]
track='C'
gwas=pd.read_csv(s/'tracks/C/results/locus_variants.tsv.gz',sep='\t',dtype={'chrom':str},keep_default_na=False)
index_summary=json.loads((s/'ld/C/index_export_summary.json').read_text())
results=[]
for locus in index_summary:
    lid=locus['locus_id'];out=s/'ld/C'/lid;out.mkdir(exist_ok=True)
    g=gwas[gwas.locus_id==lid].copy()
    idx=pd.read_csv(s/'ld/C'/(lid+'.index.tsv.gz'),sep='\t',dtype={'chrom':str},keep_default_na=False)
    assert not idx.duplicated(['chrom','pos','ref','alt']).any()
    idx=idx.rename(columns={'af':'ld_af','rsid':'ld_rsid'})
    m=g.merge(idx[['chrom','pos','ref','alt','idx','ld_af','ld_rsid']],on=['chrom','pos','ref','alt'],how='left',validate='many_to_one')
    m['summary_row_eligible']=(m.qc_status=='SOURCE_STATISTIC_QC_PASS') & (m.reference_match.astype(str)=='True') & (m.normalized_identity_duplicate.astype(str)=='False')
    m['ld_present']=m.idx.notna()
    m['retained_for_diagnostics']=m.summary_row_eligible & m.ld_present
    m['exclusion_reason']=np.where(~m.summary_row_eligible,'SOURCE_QC_OR_IDENTITY_UNRESOLVED',np.where(~m.ld_present,'EXACT_SIGNED_LD_IDENTITY_ABSENT',''))
    m['significant']=pd.to_numeric(m.neglog10_p,errors='coerce') > -np.log10(5e-8)
    m[['track','locus_id','source_row','chrom','pos','ref','alt','source_identity','normalized_variant_id','qc_status','harmonization_status','variant_qc_annotation','summary_row_eligible','ld_present','idx','retained_for_diagnostics','exclusion_reason','significant']].to_csv(out/'variant_overlap_ledger.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
    retained=m[m.retained_for_diagnostics].copy().sort_values('idx')
    assert (retained.effect_allele==retained.alt).all() and (retained.other_allele==retained.ref).all()
    assert not retained.idx.duplicated().any()
    retained['idx']=retained.idx.astype('int64');retained['orientation_multiplier']=1
    for c in ['beta','se','z']:retained[c]=pd.to_numeric(retained[c],errors='raise')
    assert np.isfinite(retained[['beta','se','z']].values).all() and (retained.se>0).all()
    retained['z_ld']=retained.z;retained['beta_ld']=retained.beta
    retained.to_csv(out/'ordered_diagnostic_inputs.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
    absent=m[m.significant & ~m.retained_for_diagnostics]
    record=dict(locus,total_gwas_locus_rows=len(m),summary_row_eligible=int(m.summary_row_eligible.sum()),
                retained_diagnostic_variants=len(retained),missing_eligible_variants=int((m.summary_row_eligible&~m.ld_present).sum()),
                significant_rows=int(m.significant.sum()),missing_or_excluded_significant_rows=len(absent),
                complete_variant_order_verified=True,summary_likelihood_contract='SAIGE_0.36.3_SPA_CALIBRATED_SE; reference-LD likelihood not established',
                density_status='BLOCKED_MISSING_OR_EXCLUDED_SIGNIFICANT_VARIANTS' if len(absent) else 'NO_STRONG_SIGNAL_LOSS; full coverage and scientific compatibility still require review')
    (out/'overlap_summary.json').write_text(json.dumps(record,indent=2)+'\n');results.append(record)
    print(lid,'retained',len(retained),'missing/excluded significant',len(absent),flush=True)
(s/'ld/C/gwas_ld_overlap_summary.json').write_text(json.dumps(results,indent=2)+'\n')
