#!/usr/bin/env python3
"""Normalize public acquisition receipts; verify saved bytes only, no scientific analysis."""
import csv,datetime,hashlib,json,pathlib,platform,sys
STAGE=pathlib.Path(__file__).resolve().parents[1]
ROOT=STAGE.parents[3]
def main():
    records=[]
    for receipt in sorted((STAGE/'audits').rglob('*.json')):
        try:x=json.loads(receipt.read_text())
        except Exception:continue
        if not isinstance(x,dict):continue
        url=x.get('requested_url',x.get('url'))
        if not url:continue
        # Accept only first-party acquisition receipts, not derived source reviews.
        if not any(k in x for k in ('request_headers','maximum_response_bytes','downloaded_bytes','acquired_bytes')):continue
        raw_path=x.get('saved_path',x.get('body_path',x.get('stored_path')))
        path=None
        if raw_path:
            candidates=[pathlib.Path(raw_path),STAGE/raw_path,ROOT/raw_path]
            path=next((p.resolve() for p in candidates if p.is_file()),None)
            assert path is not None,(receipt,raw_path)
        expected=x.get('sha256',x.get('acquired_sha256',x.get('stored_sha256',x.get('saved_sha256'))))
        size=x.get('downloaded_bytes',x.get('acquired_bytes',x.get('stored_bytes',x.get('saved_bytes',0))))
        if path:
            body=path.read_bytes();actual=hashlib.sha256(body).hexdigest()
            assert expected==actual and int(size)==len(body),(receipt,expected,actual,size,len(body))
        response={k.lower():v for k,v in x.get('response_headers',{}).items()}
        request_headers={k.lower():v for k,v in x.get('request_headers',{}).items()}
        interval=x.get('requested_range',request_headers.get('range',''))
        records.append({'receipt_path':str(receipt.relative_to(STAGE)),
          'source_url':url,'effective_url':x.get('final_url',x.get('effective_url',x.get('response_url',''))),
          'accessed_utc':x.get('accessed_utc',x.get('requested_utc',x.get('started_utc',''))),
          'http_status':x.get('http_status',x.get('status','')),'saved_path':str(path.relative_to(STAGE)) if path else '',
          'saved_bytes':size,'saved_sha256':expected or '',
          'requested_byte_range':interval or '',
          'response_content_range':response.get('content-range',''),
          'response_content_length':response.get('content-length',''),
          'bounded_or_partial':bool(interval or x.get('bounded_read') or x.get('read_truncated') or x.get('truncated') or x.get('retained_prefix_only')),
          'hash_scope':'exact saved bytes; partial/range hash is not a complete remote-object checksum',
          'error':x.get('error',x.get('transport_error','')) or ''})
    assert records
    dest=STAGE/'provenance/source_access_ledger.tsv'
    with dest.open('w') as f:
        w=csv.DictWriter(f,list(records[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(records)
    unique={r['saved_path']:int(r['saved_bytes']) for r in records if r['saved_path']}
    manifest={'stage':'fine-mapping-preflight-1.0','generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
      'access_attempt_receipts':len(records),'distinct_saved_response_files':len(unique),'saved_response_bytes':sum(unique.values()),
      'saved_response_hashes_verified':True,'all_hashes_are_saved_byte_hashes_not_assumed_full_downloads':True,
      'source_ledger':str(dest.relative_to(STAGE)),'source_ledger_sha256':hashlib.sha256(dest.read_bytes()).hexdigest(),
      'public_source_audit_only':True,'no_statistics_or_LD_matrix_analysis':True,
      'python':sys.version,'platform':platform.platform()}
    (STAGE/'provenance/source_access_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest,indent=2))
if __name__=='__main__':main()
