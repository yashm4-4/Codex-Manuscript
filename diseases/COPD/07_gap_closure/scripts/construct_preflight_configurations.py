#!/usr/bin/env python3
"""Construct fixed A/B/C interval sets, without model inference or training."""
from __future__ import annotations
import os
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
import hashlib
import json
import time
from datetime import datetime, timezone
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT=Path(__file__).resolve().parents[4]
BASE=ROOT/'diseases/COPD/07_gap_closure'
DATA=BASE/'data/COPD-V2-PREFLIGHT'
PROV=BASE/'provenance'
OUT=BASE/'results'
P='COPD-V2-PREFLIGHT'
MODELS=['enhancer','h3k27me3']
SEEDS=[104729,130363,155921]
FEATURES=['gc_fraction','atac_signal_percentile_max','repeat_fraction_2001']
CALIPERS=np.array([.05,.20,.20])
RULES={
 'version':'pretraining-1.0','training_seeds':SEEDS,'matching_seed':271828,
 'bootstrap_seed':314159,'requested_controls_per_B_positive':1,
 'ratio_rationale':'Enhancer eligible pool286943 is less than2*181352; use a prespecified common1:1 target for both models, without replacement, without dropping positives or fabricating controls.',
 'matching':'deterministic greedy nearest Chebyshev distance on GC/.05, within-source ATAC signal percentile/.20, repeat fraction/.20; require all three absolute calipers',
 'exact_strata':['chrom','validation_role','atac_lobes','blacklist_input_any','non_acgt_any'],
 'positive_processing_order':'SHA256(271828|model|interval_id), ascending; candidate coordinate-ID order resolves ties',
 'without_replacement':True,'positive_subsampling':False,
 'unmatched_positive_policy':'retain all V1 positives; no outside-caliper rescue; report every unmatched positive and realized ratio',
 'B_C_controls':'identical selected control IDs; C differs only by same-lobe positive restriction; audit balance against both positive sets',
 'chr7_roles':'union full2001bp-overlap components with encoded-sequence/RC duplicate equivalence; SHA256(chr7-role-v1|component_id) modulo4:0,1 checkpoint;2 selection;3 calibration',
 'bootstrap_unit':'global genomic overlap/encoded-sequence component; matching links are design links, not biological shared-sample units; condition on frozen matching design',
 'hard_balance_limits':{'absolute_SMD':.10,'minimum_B_matching_fraction':.90,'absolute_chromosome_proportion_gap':.02,'absolute_lobe_signature_proportion_gap':.05},
 'no_training':True,'no_model_inference':True,'no_benchmark_outcome_access':True,
}

def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 return h.hexdigest()

def write(path,frame):
 path.parent.mkdir(parents=True,exist_ok=True)
 if path.exists():raise RuntimeError(f'Refusing overwrite: {path}')
 frame.to_csv(path,sep='\t',index=False,header=not path.name.endswith('.bed.gz'),compression={'method':'gzip','mtime':0} if path.suffix=='.gz' else None,float_format='%.10g')

def json_out(path,value):
 with path.open('x') as f:json.dump(value,f,indent=2,sort_keys=True);f.write('\n')

def group_roles(frame):
 # All eligible and historical windows enter this model-independent grouping.
 order=frame.sort_values(['chrom','input_start','input_end']).index
 groups=np.empty(len(frame),dtype=np.int64)
 names=[];last_chrom=None;last_end=-1;group=-1
 for i,c,s,e in zip(order,frame.loc[order,'chrom'],frame.loc[order,'input_start'],frame.loc[order,'input_end']):
  if c!=last_chrom or s>=last_end:
   group+=1;names.append(f'{c}:{s}');last_end=e
  else:last_end=max(last_end,e)
  groups[i]=group;last_chrom=c
 parent=np.arange(group+1,dtype=np.int64)
 def find(x):
  while parent[x]!=x:parent[x]=parent[parent[x]];x=parent[x]
  return int(x)
 def union(a,b):
  a,b=find(a),find(b)
  if a!=b:parent[max(a,b)]=min(a,b)
 dup=frame[frame.canonical_rc_sequence_sha256.duplicated(keep=False)]
 for _,x in dup.groupby('canonical_rc_sequence_sha256',sort=False):
  gs=groups[x.index]
  for g in gs[1:]:union(int(gs[0]),int(g))
 resolved=np.array([find(g) for g in groups])
 frame['component_id']=[names[g] for g in resolved]
 roles={}
 for component in frame.loc[frame.chrom.eq('chr7'),'component_id'].unique():
  bucket=int(hashlib.sha256(('chr7-role-v1|'+component).encode()).hexdigest()[:16],16)%4
  roles[component]='checkpoint' if bucket<2 else 'selection' if bucket==2 else 'calibration'
 frame['validation_role']=[roles[g] if c=='chr7' else 'train' if p=='train' else 'test' for c,p,g in zip(frame.chrom,frame.partition,frame.component_id)]
 return frame

def matched_controls(frame,model):
 positive=frame[frame[f'v1_{model}_positive'].eq(1)]
 controls=frame[frame[f'{model}_accessible_control_eligible'].eq(1)]
 strata=RULES['exact_strata']
 control_groups={k:part for k,part in controls.groupby(strata,sort=True,dropna=False)}
 matched=[];unmatched=[]
 for key,pos in positive.groupby(strata,sort=True,dropna=False):
  pool=control_groups.get(key)
  if pool is None or len(pool)==0:
   unmatched.extend({'model':model,'positive_id':r,'reason':'no_control_in_exact_stratum'} for r in pos.interval_id)
   continue
  pool=pool.sort_values('interval_id')
  ci=pool.index.to_numpy(); values=pool[FEATURES].to_numpy(float)/CALIPERS
  tree=cKDTree(values);used=np.zeros(len(pool),dtype=bool)
  pos=pos.assign(_order=[hashlib.sha256(f'271828|{model}|{x}'.encode()).hexdigest() for x in pos.interval_id]).sort_values('_order')
  for pi,pid,point in zip(pos.index,pos.interval_id,pos[FEATURES].to_numpy(float)/CALIPERS):
   k=min(64,len(pool)); chosen=None; distance=None
   while True:
    distances,idx=tree.query(point,k=k,p=np.inf,distance_upper_bound=1.+1e-12)
    distances=np.atleast_1d(distances);idx=np.atleast_1d(idx)
    valid=np.isfinite(distances)&(idx<len(pool))
    candidates=idx[valid]; dd=distances[valid]
    available=~used[candidates]; candidates=candidates[available];dd=dd[available]
    if len(candidates):
     best=np.lexsort((candidates,dd))[0];chosen=int(candidates[best]);distance=float(dd[best]);break
    if k==len(pool) or not valid.all():break
    if k>=1024:
     candidates=np.asarray(tree.query_ball_point(point,r=1.+1e-12,p=np.inf),dtype=int)
     candidates=candidates[~used[candidates]]
     if len(candidates):
      dd=np.max(np.abs(values[candidates]-point),axis=1);best=np.lexsort((candidates,dd))[0];chosen=int(candidates[best]);distance=float(dd[best])
     break
    k=min(k*4,len(pool))
   if chosen is None:
    unmatched.append({'model':model,'positive_id':pid,'reason':'no_unused_control_within_all_calipers'})
   else:
    used[chosen]=True;row=frame.loc[ci[chosen]];p=frame.loc[pi]
    matched.append({'model':model,'positive_id':pid,'control_id':row.interval_id,'positive_index':int(pi),'control_index':int(ci[chosen]),'chrom':row.chrom,'validation_role':row.validation_role,'atac_lobes':row.atac_lobes,'normalized_chebyshev_distance':distance,**{f+'_absolute_difference':abs(float(p[f])-float(row[f])) for f in FEATURES}})
  print(f'{model}: stratum {key}, targets{len(pos)}, available{len(pool)}, matched cumulative{len(matched)}',flush=True)
 return pd.DataFrame(matched),pd.DataFrame(unmatched,columns=['model','positive_id','reason'])

def diagnostic(frame,model,config,indices):
 rows=[];selected=frame.loc[indices].copy()
 posflag=f'v1_{model}_positive'
 for part in ['train','validation','test']:
  x=selected[selected.partition.eq(part)];pos=x[x[posflag].eq(1)];neg=x[x[posflag].eq(0)]
  for feature in FEATURES+['blacklist_input_any','non_acgt_any']:
   a=pos[feature].to_numpy(float);b=neg[feature].to_numpy(float)
   scale=np.sqrt((np.nanvar(a,ddof=1)+np.nanvar(b,ddof=1))/2) if len(a)>1 and len(b)>1 else np.nan
   diff=np.nanmean(a)-np.nanmean(b) if len(a) and len(b) else np.nan
   smd=diff/scale if scale>0 else 0. if diff==0 else np.nan
   rows.append({'model':model,'configuration':config,'partition':part,'variable':feature,'positive_n':len(pos),'control_n':len(neg),'positive_mean':np.nanmean(a) if len(a) else np.nan,'control_mean':np.nanmean(b) if len(b) else np.nan,'SMD':smd,'status':'PASS' if np.isfinite(smd) and abs(smd)<=.10 else 'FAIL'})
  for feature,limit in [('chrom',.02),('atac_lobes',.05)]:
   a=pos[feature].value_counts(normalize=True);b=neg[feature].value_counts(normalize=True)
   for value in sorted(set(a.index)|set(b.index)):
    delta=float(a.get(value,0)-b.get(value,0))
    rows.append({'model':model,'configuration':config,'partition':part,'variable':feature+'='+str(value),'positive_n':len(pos),'control_n':len(neg),'positive_mean':float(a.get(value,0)),'control_mean':float(b.get(value,0)),'SMD':np.nan,'proportion_gap':delta,'status':'PASS' if abs(delta)<=limit else 'FAIL'})
 return rows

def main():
 started=time.time()
 if (PROV/f'{P}_freeze.json').exists() or (PROV/f'{P}_configuration_manifest.json').exists():raise RuntimeError('Version already constructed/frozen')
 json_out(PROV/f'{P}_matching_specification.json',{**RULES,'recorded_before_matching_utc':datetime.now(timezone.utc).isoformat(),'input_feature_sha256':sha(DATA/'interval_features.tsv.gz'),'implementation_sha256':sha(__file__)})
 frame=pd.read_csv(DATA/'interval_features.tsv.gz',sep='\t',keep_default_na=False)
 assert frame.interval_id.is_unique
 for f in FEATURES:frame[f]=pd.to_numeric(frame[f],errors='coerce')
 frame['blacklist_input_any']=(frame.blacklist_bp_2001>0).astype(int)
 frame['non_acgt_any']=(frame.non_acgt_fraction>0).astype(int)
 frame=group_roles(frame)
 write(DATA/'interval_role_assignment.tsv.gz',frame[['interval_id','chrom','partition','component_id','validation_role','canonical_rc_sequence_sha256']])
 sets={};pair_tables=[];unmatched_tables=[];diagnostics=[];counts=[];panel_counts=[];registry=[]
 cols=['interval_id','chrom','core_start','core_end','input_start','input_end','partition','validation_role','component_id','label','atac_lobes','atac_file_ids','atac_peak_ids','atac_anchor_count','atac_signal_percentile_max','gc_fraction','repeat_fraction_2001','blacklist_bp_1kb','blacklist_bp_2001','promoter_bp_1kb','promoter_bp_2001','sequence_available','sequence_length','non_acgt_fraction','sequence_sha256','canonical_rc_sequence_sha256']
 for model in MODELS:
  pairs,unmatched=matched_controls(frame,model);pair_tables.append(pairs);unmatched_tables.append(unmatched)
  posA=frame.index[frame[f'v1_{model}_positive'].eq(1)].to_numpy()
  negA=frame.index[frame[f'v1_{model}_control'].eq(1)].to_numpy()
  posC=frame.index[frame[f'v1_{model}_positive'].eq(1)&frame[f'{model}_same_lobe_peak_support'].eq(1)].to_numpy()
  negB=pairs.control_index.to_numpy(dtype=int)
  for config,positive,negative in [('V2-A',posA,negA),('V2-B',posA,negB),('V2-C',posC,negB)]:
   ids=np.r_[positive,negative];sets[(model,config)]=set(ids)
   x=frame.loc[ids].copy();x['label']=np.r_[np.ones(len(positive),int),np.zeros(len(negative),int)]
   x=x.sort_values(['chrom','core_start','label'])
   extra=[f'{model}_all_lobe_peak_support',f'{model}_same_lobe_peak_support',f'{model}_all_lobe_histone_file_ids',f'{model}_same_lobe_histone_file_ids',f'{model}_same_lobe_support_lobes',f'{model}_mark_overlap_1kb',f'{model}_mark_overlap_2001',f'{model}_v1_positive_input_overlap',f'{model}_accessible_control_eligible']
   path=DATA/'configurations'/f'{config}_{model}_interval_manifest.tsv.gz';write(path,x[cols+extra])
   registry.append({'configuration':config,'model':model,'manifest':str(path.relative_to(ROOT)),'n_intervals':len(x),'positive':len(positive),'control':len(negative),'sha256':sha(path)})
   for part in ['train','validation','test']:
    subset=x[x.partition.eq(part)]
    for label,word in [(1,'positive'),(0,'control')]:
     y=subset[subset.label.eq(label)]
     bed=DATA/'intervals'/f'{config}_{model}_{part}_{word}.bed.gz'
     write(bed,y[['chrom','core_start','core_end','interval_id']])
     counts.append({'configuration':config,'model':model,'partition':part,'class':word,'n_intervals':len(y),'n_components':y.component_id.nunique(),'bed':str(bed.relative_to(ROOT)),'bed_sha256':sha(bed)})
    for role,z in subset.groupby('validation_role'):
     panel_counts.append({'configuration':config,'model':model,'partition':part,'role':role,'positive':int(z.label.sum()),'control':int((z.label==0).sum()),'components':z.component_id.nunique(),'control_to_positive_ratio':int((z.label==0).sum())/int(z.label.sum()) if z.label.sum() else np.nan})
   diagnostics.extend(diagnostic(frame,model,config,ids))
  # Shared C-label challenge task used unchanged for every configuration.
  common=frame.loc[np.r_[posC,negB]].copy();common['label']=np.r_[np.ones(len(posC),int),np.zeros(len(negB),int)]
  common=common[common.partition.isin(['validation','test'])]
  write(DATA/'evaluation'/f'{model}_common_challenge_panel.tsv.gz',common[cols+extra])
  for role,y in common.groupby('validation_role'):
   panel_counts.append({'configuration':'COMMON','model':model,'partition':str(y.partition.iloc[0]),'role':role,'positive':int(y.label.sum()),'control':int((y.label==0).sum()),'components':y.component_id.nunique(),'control_to_positive_ratio':int((y.label==0).sum())/int(y.label.sum()) if y.label.sum() else np.nan})
 allpairs=pd.concat(pair_tables,ignore_index=True);allunmatched=pd.concat(unmatched_tables,ignore_index=True)
 write(DATA/'control_matching_pairs.tsv.gz',allpairs)
 write(DATA/'unmatched_positive_targets.tsv.gz',allunmatched)
 write(OUT/f'{P}_class_counts.tsv',pd.DataFrame(counts))
 write(OUT/f'{P}_role_class_counts.tsv',pd.DataFrame(panel_counts))
 diag=pd.DataFrame(diagnostics);write(OUT/f'{P}_matching_diagnostics.tsv',diag)
 overlaps=[]
 for model in MODELS:
  for a,b in [('V2-A','V2-B'),('V2-B','V2-C'),('V2-A','V2-C')]:
   x,y=sets[(model,a)],sets[(model,b)]
   overlaps.append({'model':model,'configuration_1':a,'configuration_2':b,'n_1':len(x),'n_2':len(y),'shared':len(x&y),'only_1':len(x-y),'only_2':len(y-x),'Jaccard':len(x&y)/len(x|y)})
 write(OUT/f'{P}_configuration_overlap.tsv',pd.DataFrame(overlaps))
 write(DATA/'configuration_registry.tsv',pd.DataFrame(registry))
 runs=[]
 for r in registry:
  for seed in SEEDS:
   rows=[z for z in panel_counts if z['configuration']==r['configuration'] and z['model']==r['model'] and z['role']=='train']
   n=rows[0]['positive']+rows[0]['control']
   runs.append({'run_id':f"{r['configuration']}_{r['model']}_seed{seed}",'configuration':r['configuration'],'model':r['model'],'seed':seed,'training_intervals':n,'manifest':r['manifest'],'manifest_sha256':r['sha256'],'epochs_max':50,'batch_size':256,'patience':15,'learning_rate':.001,'optimizer':'Adadelta','gpu_type':'a100','gpus':1,'cpus':8,'memory_GiB':48,'walltime_cap_hours':4,'status':'NOT_AUTHORIZED_NOT_STARTED'})
 write(DATA/'run_matrix.tsv',pd.DataFrame(runs))
 active=set().union(*sets.values());active_frame=frame.loc[sorted(active)]
 duplicate=active_frame[active_frame.canonical_rc_sequence_sha256.duplicated(keep=False)].copy()
 cross=duplicate.groupby('canonical_rc_sequence_sha256').partition.nunique()
 write(DATA/'selected_sequence_duplicate_audit.tsv.gz',duplicate[['interval_id','partition','component_id','validation_role','canonical_rc_sequence_sha256']])
 forcol=['sequence_available','blacklist_bp_2001','promoter_bp_2001','non_acgt_fraction']
 sequence_qc=[]
 for (model,config),ids in sets.items():
  x=frame.loc[sorted(ids)]
  sequence_qc.append({'configuration':config,'model':model,'n_intervals':len(x),'unavailable':int((x.sequence_available!=1).sum()),'wrong_length':int((x.sequence_length!=2001).sum()),'full_input_blacklist_overlap':int((x.blacklist_bp_2001>0).sum()),'full_input_promoter_overlap':int((x.promoter_bp_2001>0).sum()),'non_ACGT_sequences':int((x.non_acgt_fraction>0).sum()),'within_set_exact_or_RC_duplicate_rows':int(x.canonical_rc_sequence_sha256.duplicated(keep=False).sum())})
 write(OUT/f'{P}_sequence_context_qc.tsv',pd.DataFrame(sequence_qc))
 summary={'elapsed_seconds':time.time()-started,'all_master_intervals':len(frame),'selected_unique_intervals':len(active),'unique_selected_encoded_sequences':active_frame.canonical_rc_sequence_sha256.nunique(),'future_two_orientation_feature_cache_bytes':int(active_frame.canonical_rc_sequence_sha256.nunique())*36480,'cross_partition_selected_duplicate_groups':int((cross>1).sum()),'matched_controls':allpairs.groupby('model').size().to_dict(),'unmatched_positive_targets':allunmatched.groupby('model').size().to_dict(),'B_C_balance_failures':int(((diag.configuration!='V2-A')&(diag.status=='FAIL')).sum()),'configuration_registry':registry,'training_started':False,'model_inference_performed':False,'benchmark_outcomes_read':False,'matching_rules_sha256':sha(PROV/f'{P}_matching_specification.json')}
 json_out(PROV/f'{P}_configuration_manifest.json',summary)
 print(json.dumps(summary,indent=2),flush=True)

if __name__=='__main__':main()
