#!/usr/bin/env python3
"""Freeze outcome annotation, exact inputs and code before any V2 prediction."""
import csv
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path

STAGE=Path(__file__).resolve().parents[1]
def record(path):
    data=path.read_bytes()
    return {'path':str(path.relative_to(STAGE)),'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}
def main():
    assert not (STAGE/'predictions').exists(),'Prediction output must not exist before annotation freeze'
    target=STAGE/'provenance/annotation_freeze.json'
    assert not target.exists(),'Never overwrite annotation freeze'
    for name in ['preopen_independent_validation','independent_input_validation','sequence_construction_validation','annotation_rules']:
        data=json.loads((STAGE/f'provenance/{name}.json').read_text())
        assert data['status']=='PASS',name
    with (STAGE/'inputs/variant_annotations.tsv').open() as h:variants=list(csv.DictReader(h,delimiter='\t'))
    with (STAGE/'inputs/scorable_variants.tsv').open() as h:inputs=list(csv.DictReader(h,delimiter='\t'))
    exact={r['variant_id'] for r in variants if r['variant_id']}
    assert exact=={r['variant_id'] for r in inputs}
    paths=['PROSPECTIVE_SPECIFICATION.md','specification/evaluation_specification.json',
        'provenance/prospective_freeze.json','provenance/preopen_source_manifest.json',
        'provenance/preopen_independent_validation.json','provenance/benchmark_opening_event.json',
        'provenance/benchmark_ingestion_complete.json','provenance/independent_input_validation.json',
        'provenance/annotation_rules.json','provenance/sequence_construction_validation.json',
        'inputs/benchmark_snapshot.json.gz','inputs/identity_map.tsv','inputs/scorable_variants.tsv',
        'inputs/allele_sequences.npy','inputs/input_manifest.json','inputs/sequence_qc.tsv',
        'inputs/observation_annotations.tsv','inputs/variant_annotations.tsv','inputs/contextual_annotations.tsv',
        'scripts/prepare_external_inputs.py','scripts/annotate_benchmark.py','scripts/score_frozen_v2.py',
        'scripts/runtime.sh','scripts/evaluate_outputs.py','scripts/run_logged.py','scripts/freeze_annotations.py',
        'scripts/validate_external_inputs_independently.py']
    freeze={'status':'PASS','created_utc':datetime.now(timezone.utc).isoformat(),
        'stage':'external-benchmark-v2-1.0','before_any_V2_external_prediction':True,
        'no_post_outcome_tuning':True,'no_V2_predictions_viewed_for_annotation':True,
        'n_exact_variants':len(exact),'n_frozen_benchmark_identity_keys':len(variants),
        'files':[record(STAGE/p) for p in paths],
        'independent_static_inference_review':{
            'reviewer':'external_v2_geometry_contract','status':'PASS',
            'reviewed_script_sha256':'4ba004c5018d4dcf8eec5816599726c60fb58068e09d063e2e56ff2ba8139ae4',
            'scope':'Full scorer/architecture/helper review; numeric chromosome interface repaired before execution, no input change; no checkpoint, arithmetic or threshold change'},
        'independent_evaluation_design_review':{'reviewer':'external_v2_model_contract','status':'PASS_WITH_IMPLEMENTATION_COMPLETIONS',
            'completions_before_execution':['prediction hash/completion/invariance gates','locus direction denominators','all-observation result join','explicit forward/RC coverage-expansion baselines'],
            'scientific_specification_changed':False},
        'model_pre_execution_interface_repair':'Accepted frozen numeric chromosome convention rather than requiring chr prefix; no model/inference had run and no input identities changed'}
    assert record(STAGE/'scripts/score_frozen_v2.py')['sha256']==freeze['independent_static_inference_review']['reviewed_script_sha256']
    with target.open('x') as h:json.dump(freeze,h,indent=2,sort_keys=True);h.write('\n')
    print(json.dumps({'status':'PASS','annotation_freeze':record(target),'before_V2_inference':True},indent=2))
if __name__=='__main__':main()
