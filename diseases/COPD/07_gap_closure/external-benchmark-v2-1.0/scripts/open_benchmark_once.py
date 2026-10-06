#!/usr/bin/env python3
"""One logged, prospective-gated ingestion of frozen benchmark outcomes."""
import csv
import gzip
import hashlib
import io
import json
from pathlib import Path
from datetime import datetime, timezone

REPO=Path('/vf/users/Dcode/yash/Codex-Manuscript/disease-regulatory-genomics-workflow')
STAGE=REPO/'diseases/COPD/07_gap_closure/external-benchmark-v2-1.0'
EXPECTED_FREEZE='d4725233251785d90b3a768a115452e125fe5323a58768df5047043f6f020830'

def sha(data):
    return hashlib.sha256(data).hexdigest()

def record(path):
    data=path.read_bytes()
    return {'path':str(path.relative_to(STAGE)),'bytes':len(data),'sha256':sha(data)}

def write_json(path,data):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('x') as handle:
        json.dump(data,handle,indent=2,sort_keys=True)
        handle.write('\n')

def main():
    opening=STAGE/'provenance/benchmark_opening_event.json'
    if opening.exists():
        raise RuntimeError('Benchmark has already been opened. Never start a second evaluation event.')
    frozen_path=STAGE/'provenance/prospective_freeze.json'
    raw=frozen_path.read_bytes()
    assert sha(raw)==EXPECTED_FREEZE
    frozen=json.loads(raw)
    for rec in frozen['files']:
        got=record(STAGE/rec['path'])
        assert got==rec,rec['path']
    audit_path=STAGE/'provenance/preopen_independent_validation.json'
    audit=json.loads(audit_path.read_text())
    assert audit['status']=='PASS'
    sources=json.loads((STAGE/'provenance/preopen_source_manifest.json').read_text())
    records=sources['benchmark_records']
    table_records=[r for r in records if '/results/COPD-V2-BENCH-' in r['path'] and r['path'].endswith('.tsv')]
    evidence_names={
        'castaldi_evidence.json','gong_evidence.json','mechanisms_evidence.json',
        'curated_extraction.json','Table1_source_rows.json','merged_alias_resolution.json',
        'GSE109452_count_design_audit.json','COPD-V2-BENCH_benchmark_freeze.json',
        'COPD-V2-BENCH_protocol_lock.json','COPD-V2-BENCH_comparison_manifest.json',
        'COPD-V2-BENCH_source_identity_review.json','COPD-V2-BENCH_experimental_label_review.json',
        'COPD-V2-BENCH_interpretation_review.json','COPD-V2-BENCH_analysis_specification.md',
        'COPD-V2-BENCH_external_functional_benchmark_report.md',
    }
    evidence_records=[r for r in records if Path(r['path']).name in evidence_names]
    chosen=table_records+evidence_records
    assert len({r['path'] for r in chosen})==len(chosen)
    write_json(opening,{
        'status':'OPENING_AUTHORIZED','opened_utc':datetime.now(timezone.utc).isoformat(),
        'stage':'external-benchmark-v2-1.0','benchmark_version':'COPD-V2-BENCH1.0',
        'single_open_event':True,'prospective_freeze':record(frozen_path),
        'independent_preopen_validation':record(audit_path),'source_files':chosen,
        'before_first_outcome_field_read':True,'before_V2_predictions':True,
        'no_post_outcome_tuning':True,'opener_script':record(Path(__file__).resolve()),
        'interpretation':'One scientific evaluation opening. Subsequent reads use this immutable ingested snapshot; no tuning or repeat benchmark selection is permitted.'})
    snapshot={'tables':{},'evidence':{}}
    for rec in chosen:
        content=(REPO/rec['path']).read_bytes()
        assert len(content)==rec['bytes'] and sha(content)==rec['sha256'],rec['path']
        name=Path(rec['path']).name
        if rec in table_records:
            reader=csv.DictReader(io.StringIO(content.decode()),delimiter='\t')
            rows=list(reader)
            assert all(None not in row for row in rows),name
            snapshot['tables'][name]={'source':rec,'columns':reader.fieldnames,'rows':rows}
        else:
            snapshot['evidence'][name]={'source':rec,'content':json.loads(content) if name.endswith('.json') else content.decode()}
    target=STAGE/'inputs/benchmark_snapshot.json.gz'
    target.parent.mkdir(parents=True,exist_ok=True)
    with target.open('xb') as handle:
        with gzip.GzipFile(fileobj=handle,mode='wb',mtime=0) as zipped:
            zipped.write(json.dumps(snapshot,sort_keys=True,separators=(',',':')).encode())
    complete={'status':'PASS','completed_utc':datetime.now(timezone.utc).isoformat(),
        'single_open_event':record(opening),'snapshot':record(target),'table_count':len(table_records),
        'evidence_file_count':len(evidence_records),'tables':{name:{'row_count':len(v['rows']),'columns':v['columns']} for name,v in snapshot['tables'].items()},
        'original_benchmark_modified':False,'inference_started':False,'no_post_outcome_tuning':True}
    write_json(STAGE/'provenance/benchmark_ingestion_complete.json',complete)
    print(json.dumps(complete,indent=2))

if __name__=='__main__':
    main()
