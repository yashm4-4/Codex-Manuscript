#!/usr/bin/env python3
"""Summarize constructed interval covariates and costs; no models are loaded."""
from pathlib import Path
import json
import pandas as pd
import numpy as np

ROOT=Path(__file__).resolve().parents[4]
BASE=ROOT/'diseases/COPD/07_gap_closure'
DATA=BASE/'data/COPD-V2-PREFLIGHT'
OUT=BASE/'results'
PROV=BASE/'provenance'
P='COPD-V2-PREFLIGHT'

def write(name,frame):
 path=OUT/f'{P}_{name}.tsv'
 if path.exists():raise RuntimeError(f'Refusing overwrite: {path}')
 frame.to_csv(path,sep='\t',index=False,float_format='%.10g')

def main():
 if (PROV/f'{P}_freeze.json').exists():raise RuntimeError('Version frozen')
 registry=pd.read_csv(DATA/'configuration_registry.tsv',sep='\t')
 covariates=[];contexts=[]
 variables=['gc_fraction','atac_signal_percentile_max','repeat_fraction_2001','blacklist_bp_2001','promoter_bp_2001','non_acgt_fraction']
 for r in registry.itertuples():
  x=pd.read_csv(ROOT/r.manifest,sep='\t')
  for (part,label),y in x.groupby(['partition','label']):
   for variable in variables:
    values=pd.to_numeric(y[variable],errors='coerce').dropna()
    row={'configuration':r.configuration,'model':r.model,'partition':part,'class':'positive' if label else 'control','variable':variable,'total_n':len(y),'observed_n':len(values),'missing_n':len(y)-len(values),'mean':values.mean(),'sd_sample':values.std(),'min':values.min(),'q05':values.quantile(.05),'q25':values.quantile(.25),'median':values.median(),'q75':values.quantile(.75),'q95':values.quantile(.95),'max':values.max()}
    covariates.append(row)
   for col in ['atac_lobes',f'{r.model}_same_lobe_support_lobes']:
    for value,n in y[col].fillna('NO_EXACT_ANCHOR_OR_SUPPORT').value_counts().items():
     contexts.append({'configuration':r.configuration,'model':r.model,'partition':part,'class':'positive' if label else 'control','variable':col,'value':value,'n_intervals':n,'class_n':len(y)})
 write('covariate_distributions',pd.DataFrame(covariates))
 write('source_context_distribution',pd.DataFrame(contexts))
 a=pd.DataFrame(covariates)
 write('ATAC_strength_available_denominators',a[a.variable.eq('atac_signal_percentile_max')])
 pairs=pd.read_csv(DATA/'control_matching_pairs.tsv.gz',sep='\t')
 write('matched_control_counts',pairs.groupby(['model','chrom','validation_role','atac_lobes']).size().rename('matched_controls').reset_index())
 u=pd.read_csv(DATA/'unmatched_positive_targets.tsv.gz',sep='\t')
 u['chrom']=u.positive_id.str.split(':').str[0]
 write('unmatched_target_counts',u.groupby(['model','chrom','reason']).size().rename('unmatched_positives').reset_index())
 matrix=pd.read_csv(DATA/'run_matrix.tsv',sep='\t')
 seconds={'enhancer':2659,'h3k27me3':454};original={'enhancer':464262,'h3k27me3':78165}
 matrix['baseline_scaled_fit_GPU_hours']=[seconds[m]*(n/original[m])/3600 for m,n in zip(matrix.model,matrix.training_intervals)]
 write('compute_projection_by_run',matrix[['run_id','training_intervals','baseline_scaled_fit_GPU_hours','gpus','cpus','memory_GiB','walltime_cap_hours','status']])
 cfg=json.loads((PROV/f'{P}_configuration_manifest.json').read_text())
 value={'n_fits':18,'one_A100_per_fit':True,'scaled_fit_only_GPU_hours':float(matrix.baseline_scaled_fit_GPU_hours.sum()),'conservative_total_GPU_hours':[12,24],'sum_requested_walltime_caps_GPU_hours':72,'max_concurrent_fits':4,'future_per_fit_CPUs':8,'future_per_fit_RAM_GiB':48,'future_feature_cache_bytes_for_all_unique_selected_sequences_both_orientations':cfg['future_two_orientation_feature_cache_bytes'],'future_feature_cache_GiB':cfg['future_two_orientation_feature_cache_bytes']/2**30,'local_scratch_budget_GiB':120,'estimated_execution_walltime_hours_excluding_queue':[4,8],'CPU_preflight_reproduction_request':{'CPUs':8,'RAM_GiB':16,'walltime_hours':1},'actual_CPU_feature_stage_seconds':160.1689863204956,'actual_CPU_feature_peak_RSS_KiB':1941532,'actual_CPU_configuration_stage_seconds':cfg['elapsed_seconds'],'cost_uncertainty':'Linear count scaling of original logged epoch times is an estimate, not a measured V2 fit; symmetric validation, I/O, cache extraction and implementation overhead are budgeted conservatively. No historical MaxRSS was recovered.','cache_release_rule':'Do not compute test model scores before retained-model freeze; construct/cache train/chr7 first, defer final-test execution. Cache shape known without extracting any representation.','training_authorized':False,'training_started':False,'readiness':'NOT_READY_MATCHING_GATES_FAILED'}
 feature_meta=json.loads((PROV/f'{P}_feature_manifest.json').read_text())
 value['actual_CPU_feature_stage_seconds']=feature_meta['elapsed_seconds']
 value['actual_CPU_feature_peak_RSS_KiB']=feature_meta['peak_rss_kib']
 path=PROV/f'{P}_compute_plan.json'
 with path.open('x') as f:json.dump(value,f,indent=2,sort_keys=True);f.write('\n')
 print(json.dumps(value,indent=2))

if __name__=='__main__':main()
