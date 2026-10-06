#!/usr/bin/env python3
"""Freeze the prospective evaluation using opaque hashes, never benchmark labels."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
from datetime import datetime, timezone

REPO = Path('/vf/users/Dcode/yash/Codex-Manuscript/disease-regulatory-genomics-workflow')
GAP = REPO / 'diseases/COPD/07_gap_closure'
STAGE = GAP / 'external-benchmark-v2-1.0'
BASELINE = 'c9515770db4dc38064d4c1fd7fa323788040ec37'
REQUEST = Path('/home/maheshwarany2/.codex/attachments/24f34038-8f60-4753-8ba3-cb6b82b56dda/pasted-text.txt')
EXPECTED = {
    ('enhancer',104729): '1cb670f778a2d148d0ac2b71bea5c97f5a421bb64f744af1dc5b373a60ea3114',
    ('enhancer',130363): '0e42ecc9364f967220ecc6c3bb0007dc25aced5d9c360355f6e48027c39f0795',
    ('enhancer',155921): '733ca1845061508c820240d04d4b7c30d8ba0fa7c1086419eaec0b5c2bb041d5',
    ('h3k27me3',104729): '503683295ebeb48587f42748181bf6c5815f6a8c79b8dfdd4b669c931109be81',
    ('h3k27me3',130363): '56a2d331b6ed27e49f103dce228dda3c26a80ff296e91c9c740e52b9c36841c6',
    ('h3k27me3',155921): '995dc295d163a703b80bb3e451b8fc3bd7af01fb6ce041099fc5ce3d9cc82560',
}

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(8*1024*1024),b''):
            h.update(block)
    return h.hexdigest()

def record(path, base=REPO):
    return {'path':str(path.relative_to(base)), 'bytes':path.stat().st_size, 'sha256':sha(path)}

def write_json(path, data):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x') as handle:
        json.dump(data,handle,indent=2,sort_keys=True)
        handle.write('\n')

def main():
    if (STAGE/'provenance/prospective_freeze.json').exists():
        raise RuntimeError('Prospective freeze already exists; do not overwrite')
    assert subprocess.check_output(['git','-C',str(REPO),'rev-parse','HEAD'],text=True).strip()==BASELINE
    assert not subprocess.check_output(['git','-C',str(REPO),'diff','--name-only']).strip()
    assert not subprocess.check_output(['git','-C',str(REPO),'diff','--cached','--name-only']).strip()
    ledger=GAP/'provenance/COPD-V2-BENCH_artifact_checksums.tsv'
    assert sha(ledger)=='5a177150e5414451f5490c49b8042d48520fc424be138735414eb55f0959e7e0'
    with ledger.open() as handle:
        metadata=list(csv.DictReader(handle,delimiter='\t'))
    benchmark=[]
    for row in metadata:
        if 'COPD-V2-BENCH' not in row['path']:
            continue
        current=record(REPO/row['path'])
        assert current['bytes']==int(row['bytes']) and current['sha256']==row['sha256'],row['path']
        benchmark.append(current)
    assert len(benchmark)==217
    benchmark.append(record(ledger))
    checkpoints=[]
    for (context,seed),expected in EXPECTED.items():
        path=GAP/f'internal-training-1.0/runs/V2-C_{context}_seed{seed}/attempt-001/selected_checkpoint.keras'
        rec=record(path)
        assert rec['sha256']==expected and rec['bytes']==350064314
        checkpoints.append(dict(rec,context=context,seed=seed,configuration='V2-C'))
    dependencies=[]
    for path in [
        GAP/'internal-training-1.0/scripts/phase_two_contract.py',
        GAP/'internal-training-1.0/scripts/extract_phase_one.py',
        GAP/'internal-training-1.0/results/chr7_calibration/C_region_thresholds.json',
        GAP/'internal-training-1.0/provenance/checkpoint_freeze.json',
        GAP/'internal-test-1.0/scripts/runtime.sh',
        REPO/'diseases/COPD/04_modeling/scripts/02_prepare_candidate_alleles.py',
        REPO/'diseases/COPD/04_modeling/scripts/03_score_candidate_alleles.py',
        GAP/'scripts/run_reverse_complement_scoring.py',
        REPO/'diseases/COPD/04_modeling/trednet/model_phase_I/phase_one_weights.h5',
        REPO/'diseases/COPD/04_modeling/trednet/fasta/hg38.fa',
        REPO/'diseases/COPD/04_modeling/trednet/fasta/hg38.fa.fai',
    ]:
        rec=record(path)
        rec['resolved_managed_path']=str(path.resolve())
        dependencies.append(rec)
    registers=[record(GAP/name) for name in ['activity_log.tsv','gap_closure_decision_register.tsv','gap_closure_result_register.tsv']]
    source_manifest={
        'status':'PASS','created_utc':datetime.now(timezone.utc).isoformat(),
        'benchmark_version':'COPD-V2-BENCH1.0','benchmark_payloads_verified_against_original_ledger':217,
        'benchmark_records':benchmark,'checkpoints':checkpoints,'dependencies':dependencies,
        'register_baselines':registers,'baseline_commit':BASELINE,
        'benchmark_outcome_fields_opened':False,
        'preopen_access':'Only checksum-ledger path/bytes/sha256 metadata and opaque byte hashing of benchmark artifacts; no benchmark outcome/label parsing',
    }
    write_json(STAGE/'provenance/preopen_source_manifest.json',source_manifest)
    (STAGE/'provenance/user_authorization.txt').write_bytes(REQUEST.read_bytes())
    specification={
        'stage':'external-benchmark-v2-1.0','version':'1.0','baseline_commit':BASELINE,
        'benchmark_version':'COPD-V2-BENCH1.0','benchmark_records':benchmark,
        'checkpoints':checkpoints,'dependencies':dependencies,'seeds':[104729,130363,155921],
        'contexts':['enhancer','h3k27me3'],'checkpoint_configuration':'V2-C',
        'thresholds':{
            'enhancer':{'decimal17g':'0.74848511815071117','float64_hex':'0x1.7f39710000001p-1','operator':'>='},
            'h3k27me3':{'decimal17g':'0.76960810025533055','float64_hex':'0x1.8a0a12aaaaaacp-1','operator':'>='},
        },
        'ensemble':'float64 mean in ordered seeds of (float64(p_forward)+float64(p_nucleotide_RC))/2',
        'phase_I_batch_size':32,'phase_II_batch_size':256,'invariance_atol':1e-6,'invariance_rtol':1e-6,
        'sequence_length':2001,'allele_start_zero_based':1000,'maximum_allele_length':1001,
        'identity_key':['GRCh38','chrom','pos1','REF','ALT'],'new_identity_rescue':False,
        'internal_reference_percentiles':False,'inferential_statistics':False,
        'classification_metrics':['NOT_PLANNED_NO_SENSITIVITY_SPECIFICITY_FPR_AUROC_AUPRC_ENRICHMENT'],
        'quantiles':[0,.25,.5,.75,1],'quantile_method':'linear',
        'require_all_six_checkpoints_for_primary_recovery':True,
        'single_outcome_opening_event':True,'annotation_map_frozen_before_V2_predictions':True,
        'no_post_outcome_tuning':True,'no_candidate_universe_access':True,
        'full_rules':record(STAGE/'PROSPECTIVE_SPECIFICATION.md',STAGE),
        'preopen_sources':record(STAGE/'provenance/preopen_source_manifest.json',STAGE),
        'user_authorization':record(STAGE/'provenance/user_authorization.txt',STAGE),
    }
    write_json(STAGE/'specification/evaluation_specification.json',specification)
    for item in specification['thresholds'].values():
        assert format(float.fromhex(item['float64_hex']),'.17g')==item['decimal17g']
    frozen_files=[record(STAGE/p,STAGE) for p in ['PROSPECTIVE_SPECIFICATION.md','specification/evaluation_specification.json','provenance/preopen_source_manifest.json','provenance/user_authorization.txt','scripts/initialize_stage.py']]
    freeze={'status':'PASS','created_utc':datetime.now(timezone.utc).isoformat(),'stage':'external-benchmark-v2-1.0','baseline_commit':BASELINE,
            'before_any_benchmark_outcome_or_label_opening':True,'before_any_V2_external_prediction':True,
            'files':frozen_files,'checkpoint_count':6,'seeds':[104729,130363,155921],
            'evaluation_only':True,'no_post_outcome_tuning':True}
    write_json(STAGE/'provenance/prospective_freeze.json',freeze)
    print(json.dumps({'status':'PASS','prospective_freeze':record(STAGE/'provenance/prospective_freeze.json',STAGE),
                      'specification':record(STAGE/'specification/evaluation_specification.json',STAGE),
                      'benchmark_files_verified':217,'benchmark_outcomes_opened':False},indent=2))

if __name__=='__main__':
    main()
