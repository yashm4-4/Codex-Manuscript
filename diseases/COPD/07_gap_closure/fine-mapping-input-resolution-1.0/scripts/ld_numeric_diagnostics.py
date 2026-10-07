"""Prespecified matrix and descriptive z/LD diagnostics; no fine-mapping calls."""
from pathlib import Path
import datetime
import hashlib
import json
import sys
import time
import os
# Future processes only: respect the measured Slurm allocation (8 physical cores).
# Existing solvers are unaffected. This changes resource use, not the diagnostic.
if __name__ == '__main__':
    for key in ['MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS']:
        os.environ[key] = '4'
import numpy as np
import pandas as pd
import scipy
from scipy.linalg import eigh
from scipy.optimize import minimize_scalar

s=Path(__file__).resolve().parents[1]


def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()


def diagnostic(track,lid):
    start=time.monotonic();out=s/'ld'/track/lid
    t=np.load(out/'source_upper_triangle.npy',mmap_mode='r',allow_pickle=False)
    n=t.shape[0];eps=np.finfo(np.float64).eps
    assert t.shape==(n,n) and t.dtype==np.float64
    rec={'track':track,'locus_id':lid,'variants':n,'numpy':np.__version__,'scipy':scipy.__version__,
         'matrix_precision':'float64','posterior_inference':False,'matrix_repair':False,
         'thread_environment':{k:os.environ.get(k) for k in ['MKL_NUM_THREADS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS']},
         'script_sha256':sha(Path(__file__))}
    finite=True;lower_max=0.;upper_min=float('inf');upper_max=-float('inf')
    for a in range(0,n,256):
        block=np.asarray(t[a:a+256]);finite &= bool(np.isfinite(block).all())
        upper_min=min(upper_min,float(np.min(block)));upper_max=max(upper_max,float(np.max(block)))
        for off in range(len(block)):
            if a+off:lower_max=max(lower_max,float(np.max(np.abs(block[off,:a+off]))))
    diag=np.diag(t).copy()
    rec.update(all_finite=finite,original_lower_triangle_max_abs=lower_max,
               raw_entry_min=upper_min,raw_entry_max=upper_max,
               raw_diagonal_quantiles=np.quantile(diag,[0,.01,.25,.5,.75,.99,1]).tolist(),
               raw_diagonal_max_abs_deviation_from_one=float(np.max(np.abs(diag-1))))
    if not finite or lower_max!=0 or not np.all(diag>0):
        rec['numerical_status']='FAIL_ORIGINAL_REPRESENTATION_OR_DIAGONAL'
        (out/'numeric_diagnostics.json').write_text(json.dumps(rec,indent=2)+'\n');return rec
    # This documented representation conversion preserves all original source values.
    r=np.array(t,dtype=np.float64,order='F',copy=True)
    r+=t.T
    np.fill_diagonal(r,diag)
    inv=1/np.sqrt(diag)
    r*=inv[:,None];r*=inv[None,:]
    del t
    asym=0.;maxabs=0.;rmin=float('inf');rmax=-float('inf')
    for a in range(0,n,256):
        b=r[a:a+256];asym=max(asym,float(np.max(np.abs(b-r[:,a:a+256].T))))
        maxabs=max(maxabs,float(np.max(np.abs(b))));rmin=min(rmin,float(np.min(b)));rmax=max(rmax,float(np.max(b)))
    tol=100*eps*max(1,n)
    rec.update(derived_correlation_min=rmin,derived_correlation_max=rmax,
               symmetry_max_abs_error=asym,symmetry_tolerance=tol*max(1,maxabs),
               correlation_range_tolerance=tol,derived_diagonal_max_abs_error=float(np.max(np.abs(np.diag(r)-1))),
               scaling_action='EXPLICIT_Gij_DIVIDED_BY_SQRT_Gii_Gjj; candidate Pearson correlation of residualized dosages',
               release_and_statistic_scaling_compatibility='NOT_ESTABLISHED_BY_NUMERICAL_CONVERSION')
    np.save(out/'diagnostic_signed_correlation.npy',r,allow_pickle=False)
    rec['correlation_sha256']=sha(out/'diagnostic_signed_correlation.npy')
    d=pd.read_csv(out/'ordered_diagnostic_inputs.tsv.gz',sep='\t')
    z=d.z_ld.to_numpy(dtype=np.float64)
    lead=int(np.nanargmax(pd.to_numeric(d.neglog10_p,errors='coerce').to_numpy()))
    corr=r[:,lead].copy();residual=z-corr*z[lead];variance=1-corr*corr
    usable=variance>tol
    conditional=np.full(n,np.nan);conditional[usable]=residual[usable]/np.sqrt(variance[usable])
    pair=pd.DataFrame({'source_row':d.source_row,'ld_idx':d.idx,'z_ld':z,'r_with_lead':corr,
                       'conditional_residual':conditional,'near_perfect_ld':~usable,
                       'unscaled_contrast':residual})
    pair.to_csv(out/'lead_conditional_residuals.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
    rec.update(lead_source_row=int(d.iloc[lead].source_row),lead_z=float(z[lead]),
               pairwise_residual_abs_quantiles=np.nanquantile(np.abs(conditional),[0,.5,.9,.99,1]).tolist(),
               near_perfect_ld_rows=int((~usable).sum()),
               near_perfect_ld_max_abs_contrast=float(np.max(np.abs(residual[~usable]))))
    print(f'{track}/{lid}: eigendecomposition n={n}',flush=True)
    values,vectors=eigh(r,overwrite_a=True,check_finite=False,driver='evr')
    spectral_tol=100*eps*max(1,n)*max(1,float(np.max(np.abs(values))))
    rank=int((values>spectral_tol).sum());theta=vectors.T@z
    np.save(out/'eigenvalues.npy',values,allow_pickle=False)
    cond=float(values[-1]/values[values>spectral_tol][0]) if rank else None
    psd=bool(values[0]>=-spectral_tol)
    rec.update(min_eigenvalue=float(values[0]),max_eigenvalue=float(values[-1]),
               psd_tolerance=spectral_tol,psd_pass=psd,numerical_rank=rank,
               condition_positive_subspace=cond,condition_warning_threshold=1/np.sqrt(eps),
               numerical_null_space_z_energy=float(np.sum(theta[values<=spectral_tol]**2)),
               band_radius_bp=10000000,off_band_imputation=False)
    if psd:
        low=max(0.,float(-values[0]/(1-values[0])))+10*eps
        def objective(x):
            denom=(1-x)*values+x
            if np.any(denom<=0):return float('inf')
            return float(np.sum(np.log(denom)+theta*theta/denom))
        opt=minimize_scalar(objective,bounds=(low,1),method='bounded',options={'xatol':1e-8,'maxiter':1000})
        candidates=[(float(opt.x),float(opt.fun)),(1.,objective(1.)),(low,objective(low))]
        if values[0]>0:candidates.append((0.,objective(0.)))
        chosen=min(candidates,key=lambda x:x[1])
        rec['null_covariance_mismatch_diagnostic']={'s':chosen[0],'objective':chosen[1],
            'optimizer_success':bool(opt.success),'positive_covariance_lower_bound':low,
            'clearance_threshold':None,'interpretation':'Descriptive; not an LD repair or posterior inference'}
    rec['numerical_status']='PASS_NUMERICS_ONLY' if psd and asym<=tol*max(1,maxabs) and maxabs<=1+tol else 'FAIL_NUMERICAL_GATE'
    rec['scientific_LD_gate_pass']=False
    rec['scientific_LD_gate_note']='Numerical checks alone do not establish release/statistic/sample compatibility; see source contracts and track readiness.'
    rec['seconds']=time.monotonic()-start
    rec['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
    (out/'numeric_diagnostics.json').write_text(json.dumps(rec,indent=2,allow_nan=False)+'\n')
    print(json.dumps({k:rec[k] for k in ['track','locus_id','variants','numerical_status','min_eigenvalue','numerical_rank','seconds']}),flush=True)
    return rec


if __name__=='__main__':
    diagnostic(sys.argv[1],sys.argv[2])
