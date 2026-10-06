#!/usr/bin/env python3
"""Prespecified descriptive evaluation of frozen external-stage predictions only."""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path
import resource
import time
from datetime import datetime, timezone
from collections import defaultdict
import numpy as np

STAGE=Path(__file__).resolve().parents[1]
CONTEXTS=('enhancer','h3k27me3')
THRESHOLDS={'enhancer':float.fromhex('0x1.7f39710000001p-1'), 'h3k27me3':float.fromhex('0x1.8a0a12aaaaaacp-1')}

def truth(value):
    if isinstance(value,bool): return value
    if str(value).lower() in ('true','1'): return True
    if str(value).lower() in ('false','0',''): return False
    raise ValueError(f'Unrecognized Boolean {value!r}')

def rows(path):
    with path.open() as handle: return list(csv.DictReader(handle,delimiter='\t'))

def record(path):
    data=path.read_bytes()
    return {'path':str(path.relative_to(STAGE)),'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}

def check_records(records):
    for r in records:
        assert record(STAGE/r['path'])==r,r['path']

def write_json(path,value):
    with path.open('x') as handle:
        json.dump(value,handle,indent=2,sort_keys=True,allow_nan=False)
        handle.write('\n')

def write_tsv(path,data,columns=None):
    if columns is None: columns=list(dict.fromkeys(k for row in data for k in row))
    with path.open('x',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=columns,delimiter='\t',lineterminator='\n')
        writer.writeheader()
        writer.writerows(data)

def sign(x): return (x>0)-(x<0)

def fraction(n,d): return n/d if d else None

def v1_calls(row,orientation):
    available=truth(row['v1_model_evaluable' if orientation=='forward' else 'rc_model_evaluable'])
    if not available: return None
    calls={}
    for context,prefix in [('enhancer','enhancer'),('h3k27me3','silencer')]:
        threshold=float(row[f'{prefix}_region_threshold'])
        ref=float(row[f'{prefix}_{orientation}_ref_score'])
        alt=float(row[f'{prefix}_{orientation}_alt_score'])
        assert all(np.isfinite(v) for v in (threshold,ref,alt))
        region=max(ref,alt)
        assert region==float(row[f'{prefix}_{orientation}_region_score'])
        flag=truth(row[f'{prefix}_{orientation}_region_gate'])
        assert flag==(region>=threshold)
        calls[context]=flag
    calls['union']=calls['enhancer'] or calls['h3k27me3']
    return calls

def direction_summary(data,unit,subset,context):
    n=len(data)
    concordant=sum(r['concordance']=='concordant' for r in data)
    discordant=sum(r['concordance']=='discordant' for r in data)
    tied=sum(r['concordance']=='tie' for r in data)
    assert concordant+discordant+tied==n
    return {'unit':unit,'subset':subset,'context':context,'denominator':n,'concordant':concordant,
            'discordant':discordant,'tied':tied,'strict_concordance_fraction':fraction(concordant,n),
            'non_tied_denominator':n-tied,'secondary_non_tied_fraction':fraction(concordant,n-tied),
            'interpretation':'descriptive; contexts and loci are correlated; no independent-sample inference'}

def self_test():
    assert sign(0.0)==0 and sign(-1e-16)==-1 and sign(1e-16)==1
    assert truth('False') is False and truth('True') is True
    for c,t in THRESHOLDS.items():
        assert max(t,0)>=t and not max(np.nextafter(t,-np.inf),0)>=t
    assert direction_summary([{'concordance':'tie'}],'test','test','enhancer')['strict_concordance_fraction']==0
    assert direction_summary([],'test','test','enhancer')['strict_concordance_fraction'] is None
    print(json.dumps({'status':'PASS','synthetic_checks':9,'benchmark_opened':False}))

def main():
    started=time.monotonic()
    annotation_freeze=json.loads((STAGE/'provenance/annotation_freeze.json').read_text())
    assert annotation_freeze['status']=='PASS'
    check_records(annotation_freeze['files'])
    prospective=json.loads((STAGE/'provenance/prospective_freeze.json').read_text())
    check_records(prospective['files'])
    scoring=json.loads((STAGE/'predictions/inference_manifest.json').read_text())
    assert scoring['status']=='PASS' and scoring['completed'] is True
    check_records(scoring['files'])
    invariance=json.loads((STAGE/'predictions/real_network_invariance.json').read_text())
    assert invariance['status']=='PASS' and invariance['n_seed_audits']==6 and invariance['n_ensemble_audits']==2
    out=STAGE/'results'
    out.mkdir(exist_ok=False)
    write_json(out/'evaluation_started.json',{'started_utc':datetime.now(timezone.utc).isoformat(),
        'script':record(Path(__file__).resolve()),'annotation_freeze':record(STAGE/'provenance/annotation_freeze.json'),
        'prediction_manifest':record(STAGE/'predictions/inference_manifest.json')})
    observations=rows(STAGE/'inputs/observation_annotations.tsv')
    variant_annotations=rows(STAGE/'inputs/variant_annotations.tsv')
    context_observations=rows(STAGE/'inputs/contextual_annotations.tsv')
    input_rows=rows(STAGE/'inputs/scorable_variants.tsv')
    input_ids={r['variant_id'] for r in input_rows}
    assert len(input_ids)==len(input_rows)
    with gzip.open(STAGE/'inputs/benchmark_snapshot.json.gz','rt') as handle: snapshot=json.load(handle)
    source_variants=snapshot['tables']['COPD-V2-BENCH-R015_unique_variant_summary.tsv']['rows']
    source_by_key={r['benchmark_variant_key']:r for r in source_variants}
    assert len(source_by_key)==len(source_variants)==len(variant_annotations)
    score_rows=rows(STAGE/'predictions/variant_context_scores.tsv')
    scores=defaultdict(dict)
    for r in score_rows:
        v,c=r['variant_id'],r['context']
        assert v in input_ids and c in CONTEXTS and c not in scores[v]
        ref,alt=float(r['ref_score']),float(r['alt_score'])
        delta=alt-ref
        assert np.isfinite(ref) and np.isfinite(alt) and 0<=ref<=1 and 0<=alt<=1
        assert float(r['delta'])==delta and float(r['abs_delta'])==abs(delta)
        assert int(r['delta_sign'])==sign(delta)
        assert float(r['max_allele_score'])==max(ref,alt)
        assert truth(r['either_allele_region_call'])==(max(ref,alt)>=THRESHOLDS[c])
        assert truth(r['ref_region_call'])==(ref>=THRESHOLDS[c])
        assert truth(r['alt_region_call'])==(alt>=THRESHOLDS[c])
        assert float(r['threshold'])==THRESHOLDS[c]
        scores[v][c]=r
    assert set(scores)==input_ids and all(set(v)==set(CONTEXTS) for v in scores.values())
    assert len(score_rows)==2*len(input_ids)
    all_variants=[]
    methods={'V1_forward':{},'V1_RC':{},'V2_C_symmetric':{}}
    for annotation in variant_annotations:
        r=dict(annotation)
        key=r['benchmark_variant_key']; v=r['variant_id']
        source=source_by_key[key]
        assert truth(r['positive_in_scope'])==truth(source['any_in_scope_positive_assay'])
        r['v2_scorable']=v in scores
        for name,orientation in [('V1_forward','forward'),('V1_RC','rc')]:
            calls=v1_calls(source,orientation)
            r[name+'_available']=calls is not None
            if calls is not None:
                methods[name][key]=calls
            for context in (*CONTEXTS,'union'):
                r[name+'_'+context+'_region_call']='' if calls is None else calls[context]
        if v in scores:
            calls={context:truth(scores[v][context]['either_allele_region_call']) for context in CONTEXTS}
            calls['union']=any(calls.values())
            methods['V2_C_symmetric'][key]=calls
            for context in CONTEXTS:
                for field in ['ref_score','alt_score','delta','abs_delta','delta_sign','max_allele_score','either_allele_region_call']:
                    r['V2_'+context+'_'+field]=scores[v][context][field]
            r['V2_union_region_call']=calls['union']
            r['V2_missing_reason']=''
        else:
            r['V2_missing_reason']='unresolved_exact_identity' if not v else 'sequence_QC_unevaluable'
            r['V2_union_region_call']=''
        r['variant_type']='unresolved' if not v else ('SNV' if len(v.split(':')[-2])==len(v.split(':')[-1])==1 else 'indel_or_multibase')
        r['V1_presence']='exact_forward_and_RC' if key in methods['V1_forward'] and key in methods['V1_RC'] else ('forward_only' if key in methods['V1_forward'] else ('RC_only' if key in methods['V1_RC'] else 'absent_unavailable'))
        all_variants.append(r)
    by_key={r['benchmark_variant_key']:r for r in all_variants}
    by_id={r['variant_id']:r for r in all_variants if r['variant_id']}
    positive={r['benchmark_variant_key'] for r in all_variants if truth(r['positive_in_scope'])}
    primary=positive & set(methods['V2_C_symmetric'])
    reporter_all={r['benchmark_variant_key'] for r in all_variants if truth(r['reporter_positive'])}
    endogenous_all={r['benchmark_variant_key'] for r in all_variants if truth(r['endogenous_positive'])}
    reporter=reporter_all & primary
    endogenous=endogenous_all & primary
    write_tsv(out/'all_variant_evaluation.tsv',all_variants)
    evaluated_observations=[]
    for annotation in observations:
        r=dict(annotation)
        variant=by_key[r['benchmark_variant_key']]
        for name,value in variant.items():
            if name.startswith(('V1_','V2_')) or name=='v2_scorable': r[name]=value
        evaluated_observations.append(r)
    write_tsv(out/'all_observation_evaluation.tsv',evaluated_observations)
    write_tsv(out/'region_recovery_per_variant.tsv',[by_key[k] for k in sorted(positive)])
    missing=[r for r in all_variants if not r['v2_scorable']]
    write_tsv(out/'unevaluable_variants.tsv',missing,list(all_variants[0]))
    write_tsv(out/'contextual_and_excluded_evidence.tsv',context_observations)
    write_tsv(out/'frozen_denominator_rules.tsv',snapshot['tables']['COPD-V2-BENCH-R005_assayed_denominator_audit.tsv']['rows'])

    recovery=[]
    for subset,keys in [('all_in_scope_positive',positive),('reporter_positive',reporter_all),('endogenous_positive',endogenous_all)]:
        for method,mapping in methods.items():
            usable=keys & set(mapping)
            for context in (*CONTEXTS,'union'):
                n=sum(mapping[k][context] for k in usable)
                recovery.append({'subset':subset,'method':method,'context':context,'positive_keys_before_method_coverage':len(keys),
                    'evaluable_positive_denominator':len(usable),'recovered':n,'recovery_fraction':fraction(n,len(usable)),
                    'unevaluable_count':len(keys)-len(usable),'interpretation':'descriptive region case-series recovery; not sensitivity'})
    write_tsv(out/'region_recovery_summary.tsv',recovery)
    common=positive & set.intersection(*(set(x) for x in methods.values()))
    common_rows=[]
    for method,mapping in methods.items():
        for context in (*CONTEXTS,'union'):
            n=sum(mapping[k][context] for k in common)
            common_rows.append({'population':'three_method_exact_positive_intersection','method':method,'context':context,
                'denominator':len(common),'recovered':n,'fraction':fraction(n,len(common))})
    write_tsv(out/'three_method_common_identity_recovery.tsv',common_rows)
    paired=[]; transitions=[]
    for baseline in ('V1_forward','V1_RC'):
        shared=positive & set(methods[baseline]) & set(methods['V2_C_symmetric'])
        for context in (*CONTEXTS,'union'):
            counts={k:0 for k in ('retained_recovery','newly_recovered','lost','still_missed')}
            for key in sorted(shared):
                old=methods[baseline][key][context]; new=methods['V2_C_symmetric'][key][context]
                category=('retained_recovery' if new else 'lost') if old else ('newly_recovered' if new else 'still_missed')
                counts[category]+=1
                transitions.append({'baseline':baseline,'context':context,'benchmark_variant_key':key,
                    'rsids':by_key[key]['rsids'],'loci':by_key[key]['loci'],'V1_region_call':old,'V2_region_call':new,'transition':category,
                    'preexisting_known_case':by_key[key]['preexisting_known_case']})
            assert sum(counts.values())==len(shared)
            paired.append(dict(baseline=baseline,context=context,shared_positive_denominator=len(shared),**counts))
    write_tsv(out/'paired_V1_V2_region_comparison.tsv',paired)
    write_tsv(out/'paired_V1_V2_variant_transitions.tsv',transitions)
    expansion=[dict(baseline='V1_forward',benchmark_variant_key=k,rsids=by_key[k]['rsids'],loci=by_key[k]['loci'],
                enhancer_recovered=methods['V2_C_symmetric'][k]['enhancer'],h3k27me3_recovered=methods['V2_C_symmetric'][k]['h3k27me3'],
                union_recovered=methods['V2_C_symmetric'][k]['union'],interpretation='V2 sequence coverage expansion, not a paired V1 loss')
                for k in sorted(primary-set(methods['V1_forward']))]
    write_tsv(out/'newly_V2_scorable_positive_variants.tsv',expansion)
    rc_expansion=[dict(baseline='V1_RC',benchmark_variant_key=k,rsids=by_key[k]['rsids'],loci=by_key[k]['loci'],
                enhancer_recovered=methods['V2_C_symmetric'][k]['enhancer'],h3k27me3_recovered=methods['V2_C_symmetric'][k]['h3k27me3'],
                union_recovered=methods['V2_C_symmetric'][k]['union'],interpretation='V2 sequence coverage expansion, not a paired V1 loss')
                for k in sorted(primary-set(methods['V1_RC']))]
    write_tsv(out/'newly_V2_scorable_positive_variants_vs_RC.tsv',rc_expansion)

    direction=[]; direction_missing=[]
    for r in observations:
        for context in CONTEXTS:
            eligible=truth(r[context+'_direction_eligible'])
            v=r['variant_id']
            if not eligible or v not in scores:
                if truth(r['positive_in_scope']) or r.get('assay_class')=='splicing':
                    direction_missing.append({'assay_id':r['assay_id'],'variant_id':v,'rsid':r['rsid'],'context':context,
                        'experimental_state':r['experimental_state'],'assay_class':r['assay_class'],
                        'reason':r.get(context+'_direction_exclusion_reason',r.get('direction_exclusion_reason','strict experimental/mechanism direction unavailable')) if not eligible else 'sequence_or_prediction_unavailable'})
                continue
            expected=int(r[context+'_expected_sign']); assert expected in (-1,1)
            pred=scores[v][context]; predicted=int(pred['delta_sign'])
            category='tie' if predicted==0 else ('concordant' if predicted==expected else 'discordant')
            direction.append({'assay_id':r['assay_id'],'variant_id':v,'rsid':r['rsid'],'study_id':r['study_id'],
                'locus':r['locus'],'cell_context':r['cell_context'],'assay_class':r['assay_class'],'context':context,
                'expected_ALT_minus_REF_sign':expected,'ref_score':pred['ref_score'],'alt_score':pred['alt_score'],
                'delta':pred['delta'],'abs_delta':pred['abs_delta'],'predicted_sign':predicted,'concordance':category,
                'higher_activity_grch38_allele':r.get('higher_activity_grch38_allele',''),
                'preexisting_known_case':r['preexisting_known_case']})
    write_tsv(out/'direction_resolved_allelic_effects.tsv',direction)
    write_tsv(out/'direction_unevaluable_observations.tsv',direction_missing)
    consensus=[]; conflicts=[]
    for v,annotation in by_id.items():
        for context in CONTEXTS:
            if truth(annotation[context+'_direction_conflict']):
                conflicts.append({'variant_id':v,'rsids':annotation['rsids'],'context':context,'reason':'conflicting frozen experimental directions; no majority vote'})
            elif truth(annotation[context+'_direction_eligible']) and v in scores:
                expected=int(annotation[context+'_expected_sign']); predicted=int(scores[v][context]['delta_sign'])
                consensus.append({'variant_id':v,'rsids':annotation['rsids'],'loci':annotation['loci'],'context':context,
                    'expected_ALT_minus_REF_sign':expected,'predicted_sign':predicted,'delta':scores[v][context]['delta'],
                    'concordance':'tie' if predicted==0 else ('concordant' if predicted==expected else 'discordant'),
                    'preexisting_known_case':annotation['preexisting_known_case']})
    write_tsv(out/'direction_unique_variant_consensus.tsv',consensus)
    write_tsv(out/'direction_conflicts.tsv',conflicts,['variant_id','rsids','context','reason'])
    direction_summaries=[]
    for context in CONTEXTS:
        context_rows=[r for r in direction if r['context']==context]
        for subset,predicate in [('all_strict_observations',lambda r:True),('reporter',lambda r:r['assay_class'] in ('conventional_reporter','MPRA_allele_effect')),('endogenous_editing',lambda r:r['assay_class']=='endogenous_allele_editing')]:
            direction_summaries.append(direction_summary([r for r in context_rows if predicate(r)],'observation',subset,context))
        direction_summaries.append(direction_summary([r for r in consensus if r['context']==context],'unique_variant','nonconflicting_consensus',context))
    write_tsv(out/'direction_concordance_summary.tsv',direction_summaries)
    locus_direction=[]
    for context in CONTEXTS:
        for locus in sorted({r['locus'] for r in direction if r['context']==context}):
            group=[r for r in direction if r['context']==context and r['locus']==locus]
            summary=direction_summary(group,'observation',locus,context)
            summary['locus']=locus
            summary['unique_variants']=len({r['variant_id'] for r in group})
            locus_direction.append(summary)
    write_tsv(out/'direction_locus_summary.tsv',locus_direction)

    strata=defaultdict(set)
    for r in observations:
        v=r['variant_id']
        if not truth(r['positive_in_scope']) or v not in primary: continue
        for field in ('study_id','cell_context','assay_class','mechanism_in_model_scope','locus'):
            strata[(field,r.get(field,'') or 'unspecified')].add(v)
        category='reporter' if r['assay_class'] in ('MPRA_allele_effect','conventional_reporter') else ('endogenous_editing' if r['assay_class']=='endogenous_allele_editing' else 'other_in_scope_evidence')
        strata[('evidence_unit',category)].add(v)
        strata[('regulatory_mechanism', 'transcriptional_reporter_not_H3_specific' if category=='reporter' else ('endogenous_expression_partial_scope' if category=='endogenous_editing' else 'TF_binding_partial_scope'))].add(v)
    for v in primary:
        for field in ('variant_type','V1_presence'): strata[(field,by_key[v][field])].add(v)
    stratified=[]
    for (axis,value),keys in sorted(strata.items()):
        for context in (*CONTEXTS,'union'):
            n=sum(methods['V2_C_symmetric'][k][context] for k in keys)
            stratified.append({'axis':axis,'stratum':value,'context':context,'unique_positive_scorable_variants':len(keys),
                'recovered':n,'fraction':fraction(n,len(keys)),'interpretation':'small correlated descriptive strata; overlapping sets'})
    write_tsv(out/'source_mechanism_stratified_recovery.tsv',stratified)
    distributions=[]
    for subset,keys in [('all_scorable_exact',set(scores)),('in_scope_positive',primary),('reporter_positive',reporter),('endogenous_positive',endogenous)]:
        for context in CONTEXTS:
            for field in ('ref_score','alt_score','max_allele_score','delta','abs_delta'):
                values=np.asarray([float(scores[k][context][field]) for k in sorted(keys)],dtype=np.float64)
                q=np.quantile(values,[0,.25,.5,.75,1],method='linear').tolist() if len(values) else [None]*5
                distributions.append(dict(subset=subset,context=context,score=field,n=len(values),**dict(zip(('minimum','q25','median','q75','maximum'),q))))
    write_tsv(out/'continuous_score_distributions.tsv',distributions)
    states=defaultdict(set)
    for r in observations:
        states[(r['experimental_state'],r['mechanism_in_model_scope'])].add(r['benchmark_variant_key'])
    context_states=[{'experimental_state':state,'mechanism_in_model_scope':scope,'unique_benchmark_keys':len(keys),
                    'exact_V2_scorable_keys':len(keys & set(scores)),'interpretation':'overlapping context-specific evidence; null does not mean universal negative'}
                   for (state,scope),keys in sorted(states.items())]
    write_tsv(out/'context_state_populations.tsv',context_states)
    flow=[
        {'population':'original_master_observations','count':len(observations)},
        {'population':'contextual_excluded_source_rows','count':len(context_observations)},
        {'population':'all_benchmark_identity_keys','count':len(all_variants)},
        {'population':'exact_benchmark_identity_keys','count':len(by_id)},
        {'population':'exact_V2_scorable_identity_keys_all_states_mechanisms','count':len(scores)},
        {'population':'experimentally_positive_keys_all_mechanisms_including_unresolved','count':sum(truth(r['experimental_positive']) for r in all_variants)},
        {'population':'exact_experimentally_positive_keys_all_mechanisms','count':sum(truth(r['experimental_positive']) and bool(r['variant_id']) for r in all_variants)},
        {'population':'splice_only_positive_keys_out_of_model_not_false_negatives','count':sum(truth(r['splice_only']) for r in all_variants)},
        {'population':'positive_keys_with_out_of_model_observations','count':sum(truth(r['out_of_model_positive']) for r in all_variants)},
        {'population':'in_scope_positive_identity_keys_including_unresolved','count':len(positive)},
        {'population':'exact_in_scope_positive_identity_keys','count':len(positive & set(by_id))},
        {'population':'exact_V2_scorable_in_scope_positive_primary','count':len(primary)},
        {'population':'exact_V2_scorable_reporter_positive','count':len(reporter)},
        {'population':'exact_V2_scorable_endogenous_positive','count':len(endogenous)},
        {'population':'V1_forward_RC_V2_common_positive','count':len(common)},
        {'population':'strict_enhancer_direction_observations','count':sum(r['context']=='enhancer' for r in direction)},
        {'population':'strict_enhancer_direction_unique_variants_nonconflicting','count':sum(r['context']=='enhancer' for r in consensus)},
        {'population':'strict_H3K27me3_direction_observations','count':sum(r['context']=='h3k27me3' for r in direction)},
    ]
    write_tsv(out/'population_accounting.tsv',flow)
    summary={'status':'PASS','completed_utc':datetime.now(timezone.utc).isoformat(),'stage':'external-benchmark-v2-1.0',
        'population_accounting':flow,'region_recovery':recovery,'common_identity_recovery':common_rows,
        'paired_comparisons':paired,'direction_summaries':direction_summaries,
        'unresolved_in_scope_positive_keys':sorted(positive-set(by_id)),
        'new_V2_coverage_positive_count':len(expansion),'source_denominator_complete':False,
        'classification_metrics_calculated':False,'significance_tests_performed':False,
        'internal_reference_distribution_constructed':False,'candidate_universe_accessed':False,
        'thresholds_changed':False,'models_or_seeds_changed':False,
        'resources':{'wall_seconds':time.monotonic()-started,'peak_RSS_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},
        'files':[record(p) for p in sorted(out.iterdir()) if p.is_file()]}
    write_json(out/'evaluation_manifest.json',summary)
    print(json.dumps({k:summary[k] for k in ('status','population_accounting','common_identity_recovery','paired_comparisons','direction_summaries')},indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--self-test',action='store_true')
    args=parser.parse_args()
    self_test() if args.self_test else main()
