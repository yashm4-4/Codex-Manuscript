"""Read-only, incremental extraction and numeric-output validation. No eigensolver/inference.

Only writes ld/extraction_validation* and provenance/diagnostic_runtime.json.
Use --final after all planned extracts/diagnostics finish. --rehash bypasses
unchanged-file hash cache; source blocks are size/receipt checked, not rehashed.
"""
from pathlib import Path
import argparse,collections,contextlib,datetime,hashlib,io,json,math,os,sys
import numpy as np
import pandas as pd
import scipy

S=Path(__file__).resolve().parents[1];LD=S/'ld';P=argparse.ArgumentParser();P.add_argument('--final',action='store_true');P.add_argument('--rehash',action='store_true');args=P.parse_args();started=datetime.datetime.now(datetime.timezone.utc).isoformat();eps=np.finfo(np.float64).eps
cachepath=LD/'extraction_validation_hash_cache.json';cache=json.loads(cachepath.read_text()) if cachepath.exists() else {};checks=[];loci=[];block_rows=[];seed_boundaries=[]
def identity(path):
 x=path.stat();return {'bytes':x.st_size,'mtime_ns':x.st_mtime_ns,'inode':x.st_ino}
def digest(path):
 k=str(path.relative_to(S));sig=identity(path);previous=cache.get(k,{})
 if not args.rehash and previous.get('stat')==sig:return previous['sha256'],'CACHED_UNCHANGED_FILE_AFTER_PREVIOUS_INDEPENDENT_HASH'
 before=sig;h=hashlib.sha256()
 with path.open('rb') as f:
  for b in iter(lambda:f.read(16*1024*1024),b''):h.update(b)
 if identity(path)!=before:raise RuntimeError('Artifact changed during independent hash: '+k)
 cache[k]={'stat':before,'sha256':h.hexdigest(),'hashed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()};return h.hexdigest(),'RECOMPUTED'
def add(scope,name,value,detail=''):
 status='PENDING' if value is None else 'PASS' if bool(value) else 'FAIL';checks.append({'scope':scope,'check':name,'status':status,'detail':str(detail)});return value

def capture_runtime():
 target=S/'provenance/diagnostic_runtime.json';prior=json.loads(target.read_text()) if target.exists() else {'snapshots':[]};out=io.StringIO()
 with contextlib.redirect_stdout(out):np.show_config();scipy.show_config()
 proc=[];allowed={'OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS','BLIS_NUM_THREADS','VECLIB_MAXIMUM_THREADS','NUMEXPR_NUM_THREADS'}
 for p in Path('/proc').glob('[0-9]*/cmdline'):
  try:
   argv=p.read_bytes().decode().strip('\0').split('\0')
   if not any(Path(x).name=='ld_numeric_diagnostics.py' for x in argv):continue
   env=p.with_name('environ').read_bytes().decode().split('\0');thread_env={v.split('=',1)[0]:v.split('=',1)[1] for v in env if '=' in v and v.split('=',1)[0] in allowed};libs=sorted({line.split()[-1] for line in p.with_name('maps').read_text().splitlines() if any(x in line.lower() for x in ['openblas','libblas','liblapack','libmkl','libblis'])});proc.append({'pid':int(p.parent.name),'argv':argv,'executable':str(p.with_name('exe').resolve()),'thread_environment_whitelist':thread_env,'loaded_BLAS_LAPACK_library_paths':libs})
  except (OSError,UnicodeError):continue
 prior['snapshots'].append({'captured_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'validator_python':sys.version,'validator_executable':sys.executable,'validator_numpy':np.__version__,'validator_scipy':scipy.__version__,'validator_configuration':out.getvalue(),'running_diagnostic_processes':proc,'scope_note':'/proc records actual running numeric processes and whitelisted thread settings/library paths; validator package configuration is identified separately and is not assumed to equal every process runtime. Finished diagnostic reports record their NumPy/SciPy versions. No unrelated environment variables collected.'});target.write_text(json.dumps(prior,indent=2)+'\n')
capture_runtime()
contract=S/'config/ld_diagnostic_contract.json';f=json.loads(contract.with_name('ld_diagnostic_contract.freeze.json').read_text());contract_hash,_=digest(contract);add('GLOBAL','frozen_ld_contract_hash',contract_hash==f['sha256']);residual=S/'config/ld_residual_contract.json';rf=json.loads(residual.with_name('ld_residual_contract.freeze.json').read_text());add('GLOBAL','frozen_residual_contract_hash',digest(residual)[0]==rf['sha256'])
source=LD/'source/UKBB.EUR.ldadj.bm';source_meta=json.loads((source/'metadata.json').read_text());meta_hash,_=digest(source/'metadata.json');B=int(source_meta['blockSize']);grid=math.ceil(source_meta['nRows']/B);part_names=source_meta['partFiles'];block_ids=source_meta['maybeFiltered']
# Provider metadata binds each sparse grid cell to a physical part file.
block_map={int(k):v for k,v in zip(block_ids,part_names)};seen_blocks={}
for track in ['C','A']:
 planpath=LD/track/'matrix_block_plan.json'
 if not planpath.exists():add(track,'matrix_plan_available',None);continue
 plan=json.loads(planpath.read_text());add(track,'matrix_plan_source_metadata_hash',plan['source_metadata_sha256']==meta_hash);add(track,'matrix_plan_block_size',int(plan['source_grid_block_size'])==B);add(track,'only_locus_blocks_planned',plan['entire_release_requested'] is False)
 frozen=pd.read_csv(S/'tracks'/track/'results/loci.tsv',sep='\t',dtype=str,keep_default_na=False).set_index('locus_id')
 seedpath=S/'tracks'/track/'results'/('gwas_significant_seeds.tsv' if track=='C' else 'significant_variants.tsv.gz');seeds=pd.read_csv(seedpath,sep='\t',dtype=str,keep_default_na=False);seedchr=seeds['chrom'] if track=='C' else seeds['chromosome'];seedpos=pd.to_numeric(seeds['pos'] if track=='C' else seeds['base_pair_location'],errors='raise');seedhash=digest(seedpath)[0]
 for q in plan['loci']:
  lid=q['locus_id'];scope=track+'/'+lid;out=LD/track/lid;rec={'track':track,'locus_id':lid,'extraction_status':'PENDING','numeric_output_status':'PENDING'};loci.append(rec)
  left=int(frozen.loc[lid,'start']);right=int(frozen.loc[lid,'end']);sel=(seedchr.astype(str)==str(q['chrom']))&seedpos.between(left,right);ss=seeds.loc[sel];sp=seedpos.loc[sel];expected_seed_n=int(frozen.loc[lid,'n_significant_seeds' if track=='C' else 'significant_variants']);add(scope,'all_frozen_significant_seeds_accounted',len(ss)==expected_seed_n and len(ss)>0);rec.update(significant_seed_table_sha256=seedhash,significant_seed_count=len(ss),min_significant_seed_position=int(sp.min()),max_significant_seed_position=int(sp.max()),closest_significant_seed_distance_left_bp=int(sp.min()-left),closest_significant_seed_distance_right_bp=int(right-sp.max()))
  for ii,sr in ss.iterrows():seed_boundaries.append({'track':track,'locus_id':lid,'source_row':int(sr['source_row']),'chrom':str(q['chrom']),'source_position':int(seedpos.loc[ii]),'source_p_encoding':'negative_log10' if track=='C' else 'ordinary_P','source_p_value':sr['source_neglog10_p'] if track=='C' else sr['p_value'],'frozen_start':left,'frozen_end':right,'distance_left_bp':int(seedpos.loc[ii]-left),'distance_right_bp':int(right-seedpos.loc[ii])})
  for b in q['parts']:
   key=b['part'];bid=int(b['block_id']);part=source/'parts'/key;receipt=part.with_name(key+'.access.json');mapped=(block_map.get(bid)==key and bid==int(b['block_row'])+int(b['block_col'])*grid and int(b['block_row'])<=int(b['block_col']))
   if key not in seen_blocks:
    br={'part':key,'block_id':bid,'metadata_mapping_valid':mapped,'file_exists':part.exists(),'receipt_exists':receipt.exists(),'used_by':[]}
    if part.exists() and receipt.exists():
     rr=json.loads(receipt.read_text());headers=rr.get('response_headers',rr.get('headers',{}));br.update(actual_bytes=part.stat().st_size,receipt_bytes=rr.get('bytes'),http_content_length=headers.get('Content-Length'),receipt_complete=rr.get('complete'),http_status=rr.get('http_status',rr.get('status')),receipt_sha256_present=bool(rr.get('sha256')),receipt_sha256=rr.get('sha256'));br['size_receipt_pass']=(br['actual_bytes']==br['receipt_bytes'] and (br['http_content_length'] is None or br['actual_bytes']==int(br['http_content_length'])) and br['receipt_complete'] is True and br['http_status']==200 and br['receipt_sha256_present']);add('BLOCK/'+key,'complete_source_size_receipt',br['size_receipt_pass'])
    else:br['size_receipt_pass']=None;add('BLOCK/'+key,'complete_source_size_receipt',False if args.final else None)
    add('BLOCK/'+key,'source_metadata_grid_mapping',mapped);seen_blocks[key]=br
   seen_blocks[key]['used_by'].append(scope)
  ip=LD/track/(lid+'.index.tsv.gz');op=out/'ordered_diagnostic_inputs.tsv.gz';ep=out/'matrix_extraction.json';raw=out/'source_upper_triangle.npy';npth=out/'numeric_diagnostics.json'
  if not ip.exists() or not op.exists():add(scope,'index_and_ordered_table_available',False if args.final else None);continue
  ih,im=digest(ip);oh,om=digest(op);rec.update(index_sha256=ih,ordered_inputs_sha256=oh);idx=pd.read_csv(ip,sep='\t',dtype={'chrom':str});d=pd.read_csv(op,sep='\t',dtype=str,keep_default_na=False);ix=pd.to_numeric(d.idx,errors='raise').to_numpy(dtype=np.int64);n=len(d);rec['variants']=n
  add(scope,'ordered_idx_unique_strict',len(np.unique(ix))==n and bool(np.all(np.diff(ix)>0)));add(scope,'ordered_idx_bounds',n>0 and bool(np.all((ix>=0)&(ix<source_meta['nRows']))));add(scope,'source_index_unique_idx',not idx.idx.duplicated().any());indexed=idx.set_index('idx').reindex(ix);add(scope,'every_ordered_idx_in_source_index',not indexed.pos.isna().any());chrom=d['ld_chrom'] if track=='A' else d['chrom'];pos=pd.to_numeric(d['ld_pos'] if track=='A' else d['pos'],errors='raise').to_numpy(dtype=np.int64);ref=d['ld_ref'] if track=='A' else d['ref'];alt=d['ld_alt'] if track=='A' else d['alt'];exact=(np.array_equal(chrom.astype(str).to_numpy(),indexed.chrom.astype(str).to_numpy()) and np.array_equal(pos,indexed.pos.to_numpy()) and np.array_equal(ref.to_numpy(),indexed.ref.to_numpy()) and np.array_equal(alt.to_numpy(),indexed.alt.to_numpy()));add(scope,'ordered_exact_index_alleles_coordinates',exact)
  chroms=chrom.unique().tolist();span=int(pos.max()-pos.min());rec.update(ld_chromosomes=chroms,ld_min_position=int(pos.min()),ld_max_position=int(pos.max()),ld_position_span_bp=span,same_chromosome=len(chroms)==1,within_10Mb_band=span<=10000000,frozen_start=int(frozen.loc[lid,'start']),frozen_end=int(frozen.loc[lid,'end']));add(scope,'single_chromosome',len(chroms)==1 and chroms[0]==str(q['chrom']));add(scope,'all_retained_pairs_inside_10Mb_band',span<=10000000);add(scope,'complete_frozen_locus_inside_10Mb_band',rec['frozen_end']-rec['frozen_start']<=10000000);add(scope,'matrix_plan_equals_frozen_locus',str(q['chrom'])==str(frozen.loc[lid,'chrom']) and int(q['start'])==rec['frozen_start'] and int(q['end'])==rec['frozen_end']);add(scope,'ordered_positions_inside_frozen_locus',bool(np.all((pos>=rec['frozen_start'])&(pos<=rec['frozen_end']))));add(scope,'positions_nondecreasing_in_LD_order',bool(np.all(np.diff(pos)>=0)))
  # Sign transformations are only independently supported source-allele multipliers.
  mult=pd.to_numeric(d.orientation_multiplier,errors='raise').to_numpy();z=pd.to_numeric(d.z,errors='raise').to_numpy();beta=pd.to_numeric(d.beta,errors='raise').to_numpy();zld=pd.to_numeric(d.z_ld,errors='raise').to_numpy();bld=pd.to_numeric(d.beta_ld,errors='raise').to_numpy();add(scope,'signed_orientation_multipliers',bool(np.isin(mult,[-1,1]).all()));add(scope,'signed_z_and_beta_transform',np.allclose(zld,z*mult,rtol=1e-14,atol=1e-14) and np.allclose(bld,beta*mult,rtol=1e-14,atol=1e-14))
  blocks=set((ix//B).tolist());planned={(int(v['block_row']),int(v['block_col'])) for v in q['parts']};needed={(a,b) for a in blocks for b in blocks if a<=b};add(scope,'every_required_grid_block_planned',needed.issubset(planned));rec.update(required_grid_blocks=len(needed),planned_grid_blocks=len(planned))
  overlap=out/'overlap_summary.json'
  if overlap.exists():
   ov=json.loads(overlap.read_text())
   if 'provider_index_sha256' in ov:add(scope,'index_hash_matches_join_receipt',ih==ov['provider_index_sha256'])
   if 'ordered_inputs_sha256' in ov:add(scope,'ordered_hash_matches_join_receipt',oh==ov['ordered_inputs_sha256'])
  if not ep.exists() or not raw.exists():add(scope,'matrix_extraction_completed',False if args.final else None);continue
  e=json.loads(ep.read_text());add(scope,'ordered_hash_matches_extraction_receipt',oh==e['ordered_inputs_sha256']);rawhash,hm=digest(raw);add(scope,'raw_hash_matches_extraction_receipt',rawhash==e['matrix_sha256']);rec.update(raw_matrix_sha256=rawhash,raw_hash_validation=hm);t=np.load(raw,mmap_mode='r',allow_pickle=False);shapeok=t.shape==(n,n) and t.dtype==np.float64;add(scope,'raw_shape_dtype_matches_order',shapeok and e['shape']==[n,n] and e['variants']==n and e['dtype']=='float64');add(scope,'extraction_is_original_triangle_no_genotype_recomputation',e['original_upper_triangle'] is True and e['ld_recomputed_from_genotypes'] is False and e['posterior_inference'] is False);rec['extraction_status']='VALIDATED' if shapeok else 'FAIL';diag=np.diag(t).copy();rec['raw_diagonal_min']=float(diag.min());rec['raw_diagonal_max']=float(diag.max());positive=bool(np.isfinite(diag).all() and (diag>0).all());add(scope,'raw_diagonal_finite_positive',positive)
  if not npth.exists():add(scope,'numeric_diagnostics_completed',False if args.final else None);del t;continue
  num=json.loads(npth.read_text());rec.update(numeric_output_status='VALIDATED',reported_numerical_status=num.get('numerical_status'));rp=out/'diagnostic_signed_correlation.npy'
  if not rp.exists():
   expected_failure=num.get('numerical_status')=='FAIL_ORIGINAL_REPRESENTATION_OR_DIAGONAL';add(scope,'missing_derived_matrix_explained_by_original_failure',expected_failure);del t;continue
  rh,rhm=digest(rp);add(scope,'derived_hash_matches_numeric_report',rh==num.get('correlation_sha256'));rec.update(derived_matrix_sha256=rh,derived_hash_validation=rhm);r=np.load(rp,mmap_mode='r',allow_pickle=False);add(scope,'derived_shape_dtype_matches_order',r.shape==(n,n) and r.dtype==np.float64)
  # Verify the entire deterministic reconstruction, never rewrite matrices.
  tol=100*eps*max(1,n);maxerr=0.;lower=0.;finite=True;scale=1/np.sqrt(diag)
  for a in range(0,n,128):
   stop=min(n,a+128);tt=np.array(t[a:stop],copy=True);finite &= bool(np.isfinite(tt).all())
   for j in range(stop-a):
    if a+j:lower=max(lower,float(np.max(np.abs(tt[j,:a+j]))))
   exp=(tt+np.asarray(t[:,a:stop]).T);exp[np.arange(stop-a),np.arange(a,stop)]=diag[a:stop];exp*=scale[a:stop,None];exp*=scale[None,:];maxerr=max(maxerr,float(np.max(np.abs(exp-np.asarray(r[a:stop])))))
  rec.update(raw_lower_triangle_max_abs=lower,raw_finite=finite,reconstruction_formula='(U+U^T-diag(diag(U)))_ij / sqrt(U_ii U_jj)',reconstruction_max_abs_error=maxerr,reconstruction_absolute_tolerance=tol);add(scope,'raw_upper_triangle_representation',finite and lower==0.);add(scope,'whole_matrix_reconstruction_formula',maxerr<=tol);add(scope,'numeric_report_bound_to_same_n',num.get('variants')==n)
  epth=out/'eigenvalues.npy'
  if epth.exists():
   values=np.load(epth,allow_pickle=False);add(scope,'saved_eigenvalue_shape_and_finiteness',values.shape==(n,) and bool(np.isfinite(values).all()));spectral_tol=100*eps*max(1,n)*max(1,float(np.max(np.abs(values))));rank=int((values>spectral_tol).sum());cond=float(values[-1]/values[values>spectral_tol][0]) if rank else None;threshold=1/np.sqrt(eps);rec.update(numerical_rank=rank,rank_deficient=rank<n,positive_subspace_condition=cond,condition_warning_threshold=threshold,condition_exceeds_warning=cond is not None and cond>threshold,numerical_null_space_z_energy=num.get('numerical_null_space_z_energy'),singular_support_status='DESCRIPTIVE_NULLSPACE_ENERGY_REQUIRES_REVIEW' if rank<n else 'NO_NUMERICAL_NULL_SPACE_AT_FROZEN_THRESHOLD');add(scope,'reported_rank_matches_saved_eigenvalues',rank==num.get('numerical_rank'));add(scope,'reported_PSD_matches_saved_eigenvalues',bool(values[0]>=-spectral_tol)==num.get('psd_pass'));add(scope,'reported_spectral_tolerance_matches_contract',np.isclose(spectral_tol,num.get('psd_tolerance',float('nan')),rtol=1e-14,atol=0));add(scope,'reported_condition_matches_saved_eigenvalues',cond is None and num.get('condition_positive_subspace') is None or cond is not None and np.isclose(cond,num.get('condition_positive_subspace',float('nan')),rtol=1e-12,atol=0));rec['eigenvalues_sha256']=digest(epth)[0]
  else:add(scope,'eigenvalues_available_for_finished_numeric_report',False)
  add(scope,'no_scientific_clearance_from_numerics',num.get('scientific_LD_gate_pass') is False);del t,r
  print(scope,'validated',flush=True)

block_rows=list(seen_blocks.values());counts=collections.Counter(x['status'] for x in checks);final_status='FAIL' if counts['FAIL'] else 'INCOMPLETE_PENDING' if counts['PENDING'] else 'PASS_EXTRACTION_AND_OUTPUT_VALIDATION_ONLY';summary={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'started_utc':started,'mode':'final' if args.final else 'incremental','status':final_status,'check_counts':dict(counts),'planned_loci':len(loci),'source_unique_blocks_checked':len(block_rows),'source_blocks_rehashed':False,'raw_and_derived_matrices_independently_hashed_once_with_stat_bound_cache':True,'cache_note':'Unchanged size/mtime/inode reuses prior independently computed SHA256; --rehash forces recomputation. Root final payload freeze independently hashes all files.','frozen_ld_contract_sha256':contract_hash,'loci':loci,'checks':checks,'scientific_execution_clearance':False,'no_eigendecomposition_or_posterior_inference_performed':True}
(LD/'extraction_validation.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n');pd.DataFrame(loci).to_csv(LD/'extraction_validation_loci.tsv',sep='\t',index=False);pd.DataFrame(checks).to_csv(LD/'extraction_validation_checks.tsv',sep='\t',index=False);pd.DataFrame([{**x,'used_by':';'.join(sorted(x['used_by']))} for x in block_rows]).to_csv(LD/'extraction_validation_source_blocks.tsv',sep='\t',index=False);pd.DataFrame(seed_boundaries).to_csv(LD/'extraction_validation_seed_boundaries.tsv',sep='\t',index=False);cachepath.write_text(json.dumps(cache,indent=2)+'\n');print(json.dumps({'status':final_status,'checks':dict(counts),'loci':len(loci),'source_blocks':len(block_rows)}),flush=True)
if counts['FAIL']:sys.exit(1)
