#!/usr/bin/env python3
"""Seal completed external evaluation without rerunning any scientific calculation."""
import csv
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import subprocess

STAGE=Path(__file__).resolve().parents[1]
REPO=STAGE.parents[3]
BASELINE='c9515770db4dc38064d4c1fd7fa323788040ec37'

def record(path,base=STAGE):
    digest=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(8*1024*1024),b''):digest.update(block)
    return {'path':str(path.relative_to(base)),'bytes':path.stat().st_size,'sha256':digest.hexdigest()}

def main():
    ledger=STAGE/'provenance/artifact_checksums.tsv'
    freeze=STAGE/'provenance/freeze.json'
    assert not ledger.exists() and not freeze.exists(),'Never overwrite final scientific freeze'
    assert subprocess.check_output(['git','-C',str(REPO),'rev-parse','HEAD'],text=True).strip()==BASELINE
    for name in ['provenance/prospective_freeze.json','provenance/annotation_freeze.json']:
        data=json.loads((STAGE/name).read_text())
        assert data['status']=='PASS'
        for item in data['files']:
            assert record(STAGE/item['path'])==item,item['path']
    for name in ['predictions/inference_manifest.json','results/evaluation_manifest.json',
                 'provenance/independent_final_validation.json','provenance/shared_register_updates.json']:
        assert json.loads((STAGE/name).read_text())['status']=='PASS',name
    for name in ['EXTERNAL_BENCHMARK_V2_REPORT.md','LIMITATIONS.md']:
        assert (STAGE/name).stat().st_size>0
    register_updates=json.loads((STAGE/'provenance/shared_register_updates.json').read_text())
    for item in register_updates['registers']:
        actual=record(REPO/item['path'],REPO)
        assert actual['bytes']==item['bytes'] and actual['sha256']==item['sha256']
    entries=[record(path) for path in sorted(STAGE.rglob('*')) if path.is_file() and '__pycache__' not in path.parts and path.suffix!='.pyc']
    with ledger.open('x',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=['path','bytes','sha256'],delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(entries)
    data={
        'status':'FROZEN_EXTERNAL_EVALUATION_COMPLETE','version':'external-benchmark-v2-1.0',
        'created_utc':datetime.now(timezone.utc).isoformat(),'baseline_commit_unchanged':BASELINE,
        'n_payloads':len(entries),'total_payload_bytes':sum(r['bytes'] for r in entries),
        'checksum_scope':'All regular new-stage artifacts except pycache/pyc and checksum-ledger/freeze self-references; original source benchmarks and model archives remain external immutable dependencies, not copies.',
        'artifact_checksums':record(ledger),'prospective_freeze':record(STAGE/'provenance/prospective_freeze.json'),
        'annotation_freeze':record(STAGE/'provenance/annotation_freeze.json'),
        'single_benchmark_opening_event':record(STAGE/'provenance/benchmark_opening_event.json'),
        'report':record(STAGE/'EXTERNAL_BENCHMARK_V2_REPORT.md'),
        'independent_final_validation':record(STAGE/'provenance/independent_final_validation.json'),
        'shared_register_updates':record(STAGE/'provenance/shared_register_updates.json'),
        'frozen_C_checkpoint_count':6,'all_original_checkpoints_preserved':True,
        'single_frozen_model_evaluation':True,'retraining_or_retuning_or_recalibration':False,
        'threshold_or_seed_or_checkpoint_change':False,'broader_COPD_candidate_scoring':False,
        'commit_or_push_performed':False,'next_action':'STOP_NO_OTHER_MODULE_AUTHORIZED',
        'external_benchmark_is_no_longer_untouched_for_future_redesign':True}
    with freeze.open('x') as handle:json.dump(data,handle,indent=2,sort_keys=True);handle.write('\n')
    print(json.dumps({'status':'PASS','frozen_payloads':len(entries),'payload_bytes':data['total_payload_bytes'],
                      'freeze':record(freeze),'checksum_ledger':record(ledger)},indent=2))

if __name__=='__main__':main()
