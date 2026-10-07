"""Independent full-byte verification of a completed input-resolution freeze.

This standalone verifier does not import the freeze implementation, never reads
its hash cache for reuse, and changes only the explicit verification receipt.
Redirect stdout/stderr outside the stage. It performs no scientific analysis.
"""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import argparse,concurrent.futures,csv,datetime,hashlib,json,os,stat,subprocess,traceback

S=Path(__file__).resolve().parents[1];REPO=S.parents[3]
LEDGER='provenance/artifact_checksums.tsv';FREEZE='provenance/freeze.json';RECEIPT='provenance/freeze_verification.json'
CLOSURE={LEDGER,FREEZE,RECEIPT}

def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def digest_bytes(b):return hashlib.sha256(b).hexdigest()
def is_bytecode(rel):return Path(rel).parent.name=='__pycache__' and Path(rel).name.endswith('.pyc')
def sig(p):
    x=p.lstat();return (x.st_size,x.st_mtime_ns,x.st_ctime_ns,x.st_ino,x.st_dev,x.st_mode)
def walk():
    result={}
    for base,dirs,files in os.walk(S,followlinks=False):
        for name in list(dirs):
            p=Path(base)/name
            if p.is_symlink():result[str(p.relative_to(S))]=sig(p);dirs.remove(name)
        for name in files:
            p=Path(base)/name;x=sig(p)
            if not (stat.S_ISREG(x[-1]) or stat.S_ISLNK(x[-1])):raise RuntimeError('Non-regular artifact: '+str(p))
            result[str(p.relative_to(S))]=x
    return result

def full_hash(p,expected_sig=None):
    before=sig(p)
    if expected_sig is not None and before!=expected_sig:raise RuntimeError('File changed before independent hashing: '+str(p))
    h=hashlib.sha256();n=0
    if stat.S_ISLNK(before[-1]):
        b=os.readlink(p).encode('utf-8');h.update(b);n=len(b);kind='symlink_target_utf8'
    else:
        kind='regular_file'
        with p.open('rb') as f:
            for b in iter(lambda:f.read(16*1024*1024),b''):h.update(b);n+=len(b)
    if sig(p)!=before:raise RuntimeError('File changed during independent hashing: '+str(p))
    return n,h.hexdigest(),kind

def require(condition,message):
    if not condition:raise AssertionError(message)
def check_registers(freeze,init):
    baseline={r['path']:r for r in init['baseline_registers']};records=[]
    require({x['path'] for x in freeze['registers']}==set(baseline),'Register inventory differs from initialization')
    for r in freeze['registers']:
        p=REPO/r['path'];before=sig(p);b=p.read_bytes();require(sig(p)==before,'Register changed while reading')
        require(len(b)==r['bytes'] and digest_bytes(b)==r['sha256'],'Register full hash/size mismatch: '+r['path'])
        q=baseline[r['path']];require(r['baseline_prefix_bytes']==q['bytes'] and r['baseline_prefix_sha256']==q['sha256'],'Register baseline declaration mismatch')
        prefix=b[:q['bytes']];require(digest_bytes(prefix)==q['sha256'],'Register original prefix changed: '+r['path'])
        require((S/r['baseline_snapshot']).read_bytes()==prefix,'Saved register baseline prefix differs')
        require(len(b)>q['bytes'] and b'fine-mapping-input-resolution-1.0' in b[q['bytes']:],'New stage register append missing')
        rows=list(csv.reader(b.decode().splitlines(),delimiter='\t'))
        require(all(len(x)==len(rows[0]) for x in rows),'Register field count mismatch')
        require(len({x[0] for x in rows[1:]})==len(rows)-1,'Duplicate register identifier')
        require(len(b)-q['bytes']==r['appended_bytes'] and b[q['bytes']:].count(b'\n')==r['appended_rows'],'Register append count mismatch')
        records.append({'path':r['path'],'sha256':digest_bytes(b),'baseline_prefix_sha256':digest_bytes(prefix),'appended_rows':r['appended_rows'],'append_only_pass':True})
    return records

def verify_required(freeze):
    req=freeze['required_final_artifacts'];require(set(req)=={'integrity','extraction','report','review'},'Missing final prerequisite binding')
    for q in req.values():
        p=S/q['path'];n,h,_=full_hash(p);require(n==q['bytes'] and h==q['sha256'],'Required final artifact mismatch: '+q['path'])
    integrity=json.loads((S/req['integrity']['path']).read_text());require(integrity['status']=='PASS' and integrity['checks_failed']==0 and integrity['checks_passed']>0,'Stage integrity not PASS')
    require(all(x['pass'] for x in integrity['checks']),'Failed stage integrity check')
    extraction=json.loads((S/req['extraction']['path']).read_text());require(extraction['mode']=='final' and extraction['status']=='PASS_EXTRACTION_AND_OUTPUT_VALIDATION_ONLY' and extraction['planned_loci']==16,'Extraction validation not final PASS')
    require(extraction['check_counts'].get('FAIL',0)==0 and extraction['check_counts'].get('PENDING',0)==0,'Extraction checks fail or remain pending')
    p=S/req['review']['path']
    if p.suffix.lower()=='.json':
        r=json.loads(p.read_text());require(r.get('status') in ['PASS','PASS_FINAL_PRE_FREEZE_REVIEW'] and r.get('final') is True,'Review is not final PASS')
    else:require('FREEZE_REVIEW_STATUS: PASS_FINAL' in [x.strip() for x in p.read_text().splitlines()],'Review final PASS marker missing')
    return {k:q['sha256'] for k,q in req.items()}

def check_source_crosslinks(expected,verified):
    """Bind acquisition-time SHA256/bytes to already independently hashed payload.

    This adds no large-payload reads: verified contains the single full rehash.
    Failed/incomplete receipts remain evidence; they are not completed sources.
    The two explicit aliases describe existing receipt/body filename differences,
    never a hash-selected substitute or a transformation.
    """
    aliases={
      'sources/software/susieR_0.16.6_github.access.json':'sources/software/susieR_0.16.6_github.tar.gz',
      'tracks/C/sources/phenotype_manifest_correct.tsv.access.json':'tracks/C/sources/phenotype_manifest_correct.tsv.bgz'}
    completed={};negative=[];error_bodies=[];mapping_records=[]
    def valid_sha(h):return isinstance(h,str) and len(h)==64 and all(c in '0123456789abcdef' for c in h)
    def completed_claim(q):
        flags=[q[k] for k in ['complete','complete_http','integrity_pass','full_size_pass'] if k in q]
        code=q.get('http_status',q.get('status_code',q.get('status')))
        try:good_http=200<=int(code)<300
        except (TypeError,ValueError):good_http=False
        if any(x is True for x in flags):
            require(not q.get('error') and not any(x is False for x in flags),'Conflicting completion/error receipt: '+str(q))
            require(code is None or good_http,'Completed receipt records HTTP failure')
            return True
        return good_http and not q.get('error') and not any(x is False for x in flags)
    def bind(rel,h,n,label):
        require(rel in expected and rel in verified,'Completed source/evidence missing from frozen payload: '+label+' -> '+rel)
        require(valid_sha(h),'Missing/malformed acquisition SHA256: '+label)
        require(verified[rel][1]==h and expected[rel]['sha256']==h,'Original acquisition SHA256 differs from frozen source bytes: '+label+' -> '+rel)
        if n is not None:require(verified[rel][0]==int(n) and int(expected[rel]['bytes'])==int(n),'Acquisition byte count differs from source: '+label)
    def mapped_payload(receipt,q):
        found=set();body=receipt[:-len('.access.json')]
        if body in expected:found.add(body)
        if receipt in aliases:
            # An explicit alias is mandatory when the successful receipt exists.
            found.add(aliases[receipt])
        for field in ['payload_path','path']:
            value=q.get(field)
            if not isinstance(value,str) or not value:continue
            v=Path(value);candidates=[v] if v.is_absolute() else [S/v,REPO/v]
            parts=Path(receipt).parts
            if len(parts)>2 and parts[0]=='tracks':candidates.append(S/parts[0]/parts[1]/v)
            for candidate in candidates:
                candidate=candidate.resolve()
                if S in candidate.parents:
                    rel=str(candidate.relative_to(S))
                    if rel in expected:found.add(rel)
        require(len(found)<=1,'Ambiguous acquisition receipt path mapping: '+receipt+' -> '+str(sorted(found)))
        return next(iter(found)) if found else None
    receipts=sorted(r for r in expected if r.endswith('.access.json'))
    for receipt in receipts:
        q=json.loads((S/receipt).read_text());complete=completed_claim(q);rel=mapped_payload(receipt,q)
        h=q.get('sha256',q.get('acquired_sha256'));n=q.get('bytes',q.get('actual_bytes'))
        if complete:
            require(rel is not None,'Completed acquisition receipt has no frozen source mapping: '+receipt)
            bind(rel,h,n,receipt);completed[receipt]=rel
            if receipt in aliases:mapping_records.append({'receipt':receipt,'payload':rel,'reason':'Explicit existing receipt/body filename difference; exact acquisition hash verified'})
        else:
            negative.append({'receipt':receipt,'status':q.get('http_status',q.get('status_code',q.get('status'))),'error':q.get('error'),'source_completion_claim':False,'retained_body':rel})
            if rel is not None and h:
                bind(rel,h,n,receipt+' (retained negative/error response only)');error_bodies.append(rel)
    completed_payloads=set(completed.values())
    root_groups={}
    for name,prefix in [('ld_matrix_blocks','ld/source/UKBB.EUR.ldadj.bm/parts/'),('ld_variant_index_row_parts','ld/source/UKBB.EUR.ldadj.variant.ht/rows/parts/'),('ld_variant_index_binary_metadata','ld/source/UKBB.EUR.ldadj.variant.ht/index/')]:
        payloads={r for r in expected if r.startswith(prefix) and not r.endswith('.access.json')}
        require(payloads.issubset(completed_payloads),'Root LD source missing a completed acquisition SHA256 crosslink: '+str(sorted(payloads-completed_payloads)))
        root_groups[name]=len(payloads)
    if root_groups['ld_matrix_blocks']:
        require(root_groups['ld_matrix_blocks']==485,'Final completed root LD source block count is not 485')
    table_names=['tables/full_file_acquisition.tsv','tables/complete_file_integrity.tsv']
    expected_ids={'A_NATIVE','A_FORMATTED37','A_HARMONIZED38','B_ARCHIVE','B_COMBINED_AUTOSOMAL','C_J44_EUR','C_VARIANT_MANIFEST'}
    real_source_present=any(r.startswith(('tracks/A/sources/GCST90016588','tracks/B/acquisition/','tracks/C/data/')) for r in expected)
    central={};table_rows={};evidence_links=0;original_receipt_links=0
    if real_source_present or any(t in expected for t in table_names):
        require(all(t in expected for t in table_names),'Both completed central source-acquisition/integrity tables are required')
        for name in table_names:
            rows=list(csv.DictReader((S/name).open(),delimiter='\t'))
            require(len(rows)==7 and {q['source_id'] for q in rows}==expected_ids,'Central completed source inventory must contain the exact seven distinct records: '+name)
            table_rows[name]={q['source_id']:q for q in rows}
            for q in rows:
                label=name+':'+q['source_id'];bind(q['payload_path'],q['recorded_sha256'],q['bytes'],label)
                bind(q['receipt_path'],q['receipt_sha256'],None,label+' receipt');evidence_links+=1
                receipt=json.loads((S/q['receipt_path']).read_text())
                require(completed_claim(receipt),'Central source links to incomplete/failed original receipt: '+label)
                original=receipt
                member=q.get('archive_member','NOT_APPLICABLE')
                if member not in ['', 'NOT_APPLICABLE']:
                    original=receipt.get('analyzed_member',{})
                    require(original.get('name')==member,'Archive member has no exact original acquisition receipt mapping: '+label)
                bind(q['payload_path'],original.get('sha256',original.get('acquired_sha256')),original.get('bytes',original.get('actual_bytes')),label+' original receipt')
                original_receipt_links+=1
                if 'complete_read_status' in q:require(q['complete_read_status'].startswith('PASS_') and q['current_size_matches_recorded']=='true','Central completed full-file audit is not PASS: '+label)
                if q.get('audit_path'):
                    bind(q['audit_path'],q['audit_sha256'],None,label+' audit');evidence_links+=1
                for field in ['all_receipt_paths_sha256','provider_digest_evidence']:
                    if q.get(field):
                        for path,h in json.loads(q[field]).items():bind(path,h,None,label+' '+field);evidence_links+=1
            central[name]={'records':len(rows),'source_ids':sorted(q['source_id'] for q in rows),'payload_sha256_and_bytes_match_original_receipts':True}
        for source_id in expected_ids:
            a=table_rows[table_names[0]][source_id];b=table_rows[table_names[1]][source_id]
            for field in set(a)&set(b):require(a[field]==b[field],'Central source tables disagree: '+source_id+' '+field)
    return {'status':'PASS','hash_basis':'Acquisition-time source digests compared to final ledger and the same independently rehashed payload digests; no extra source-file hash pass',
      'access_receipts_seen':len(receipts),'completed_access_receipt_sha256_links':len(completed),'unique_completed_access_payloads':len(completed_payloads),
      'root_ld_source_crosslinks':root_groups,'central_completed_source_tables':central,'central_unique_completed_source_records':len(expected_ids) if central else 0,
      'central_original_receipt_sha256_links':original_receipt_links,'central_evidence_receipt_audit_sha256_links':evidence_links,
      'negative_or_incomplete_receipts_not_counted_as_completed':negative,'retained_negative_response_bodies_sha256_checked':len(error_bodies),
      'explicit_receipt_body_aliases':mapping_records,'completed_receipts_missing_source_mapping':0}

def run(workers):
    frozen_bytes=(S/FREEZE).read_bytes();freeze=json.loads(frozen_bytes);ledger_bytes=(S/LEDGER).read_bytes()
    require(freeze['checksum_ledger']==LEDGER and freeze['checksum_ledger_sha256']==digest_bytes(ledger_bytes),'Freeze/ledger SHA256 binding mismatch')
    require(set(freeze['exclusions']['closure_paths'])==CLOSURE and freeze['exclusions']['all_other_files_included'] is True,'Unexpected closure exclusions')
    require(freeze['status']=='FROZEN_INPUT_RESOLUTION_STOP_FOR_INVESTIGATOR_REVIEW','Unexpected freeze state')
    require(freeze['scientific_execution_authorized'] is False and freeze['statistical_fine_mapping_executed'] is False,'Execution firewall declaration violated')
    rows=list(csv.DictReader(ledger_bytes.decode().splitlines(),delimiter='\t'));expected={r['path']:r for r in rows}
    require(len(rows)==len(expected)==freeze['payload_files'],'Duplicate or incorrect payload count')
    require(sum(int(r['bytes']) for r in rows)==freeze['payload_bytes'],'Payload byte total mismatch')
    for rel in expected:
        q=Path(rel);require(not q.is_absolute() and '..' not in q.parts and rel not in CLOSURE and not is_bytecode(rel),'Unsafe or excluded ledger path: '+rel)
    before=walk();actual={r:x for r,x in before.items() if r not in CLOSURE and not is_bytecode(r)}
    require(set(actual)==set(expected),'Payload inventory mismatch; missing='+str(sorted(set(expected)-set(actual)))+' extra='+str(sorted(set(actual)-set(expected))))
    ep=S/freeze['exclusions']['ephemeral_inventory'];er=list(csv.DictReader(ep.open(),delimiter='\t'))
    require(all(is_bytecode(x['path']) for x in er),'Non-bytecode artifact claimed as ephemeral')
    require({x['path'] for x in er}=={r for r in before if is_bytecode(r)},'Excluded bytecode inventory set differs from frozen inventory')
    print(json.dumps({'mode':'independent_full_rehash','files':len(rows),'bytes':freeze['payload_bytes'],'workers':workers,'cache_used':False}),flush=True)
    verified={};failures=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        pending={ex.submit(full_hash,S/r,actual[r]):r for r in expected}
        for i,f in enumerate(concurrent.futures.as_completed(pending),1):
            rel=pending[f]
            try:
                n,h,k=f.result();e=expected[rel]
                if n!=int(e['bytes']) or h!=e['sha256'] or k!=e['kind']:failures.append({'path':rel,'expected':e,'observed':{'bytes':n,'sha256':h,'kind':k}})
                verified[rel]=(n,h,k)
            except Exception as e:failures.append({'path':rel,'error':str(e)})
            if i%100==0:print('Independently rehashed',i,'of',len(rows),flush=True)
    require(not failures,'Payload verification failures: '+json.dumps(failures))
    source_crosslinks=check_source_crosslinks(expected,verified)
    init=json.loads((S/'provenance/initialization.json').read_text());registers=check_registers(freeze,init);required=verify_required(freeze)
    head=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip();require(head==init['baseline_commit']==freeze['baseline_commit'],'HEAD changed after stage initialization')
    require(subprocess.check_output(['git','ls-tree','-r',head],cwd=REPO)==(S/'provenance/baseline_git_tree.tsv').read_bytes(),'Historical Git tree changed')
    changed=subprocess.check_output(['git','diff','--name-only',head,'--'],cwd=REPO,text=True).splitlines();require(set(changed).issubset({x['path'] for x in init['baseline_registers']}),'Tracked historical content changed outside registers')
    pre=REPO/init['governing_preflight'];require(full_hash(pre/'provenance/freeze.json')[1]==init['preflight_freeze_sha256'],'Governing preflight freeze changed');require(full_hash(pre/'provenance/artifact_checksums.tsv')[1]==init['preflight_ledger_sha256'],'Governing preflight ledger changed')
    old=list(csv.DictReader((pre/'provenance/artifact_checksums.tsv').open(),delimiter='\t'))
    for q in old:require(full_hash(pre/q['path'])[1]==q['sha256'],'Prior preflight scientific payload changed: '+q['path'])
    after=walk();after_payload={r:x for r,x in after.items() if r not in CLOSURE and not is_bytecode(r)}
    require(after_payload==actual,'Payload changed during independent verification')
    require((S/FREEZE).read_bytes()==frozen_bytes and (S/LEDGER).read_bytes()==ledger_bytes,'Checksum closure changed during verification')
    require(check_registers(freeze,init)==registers,'Registers changed during independent verification')
    return {'status':'PASS','verification_mode':'INDEPENDENT_FULL_REHASH_NO_CACHE','completed_utc':utc(),'stage':'fine-mapping-input-resolution-1.0',
      'freeze_path':FREEZE,'freeze_sha256':digest_bytes(frozen_bytes),'checksum_ledger':LEDGER,'checksum_ledger_sha256':digest_bytes(ledger_bytes),
      'payload_files_verified':len(rows),'payload_bytes_verified':sum(x[0] for x in verified.values()),'exact_payload_inventory_pass':True,'excluded_ephemeral_inventory_pass':True,
      'all_payload_sha256_size_and_kind_matches':True,'source_acquisition_receipt_sha256_checks':source_crosslinks['completed_access_receipt_sha256_links'],'source_byte_provenance_crosslinks':source_crosslinks,'registers':registers,'required_final_artifact_sha256':required,
      'HEAD':head,'baseline_tree_unchanged':True,'tracked_changes_only_authorized_registers':True,'prior_preflight_payloads_fully_rehashed':len(old),
      'no_payload_mutations_during_verification':True,'hash_cache_used':False,'scientific_execution_authorized':False,
      'scope':'Byte/provenance/register/history verification only; no scientific computations or posterior inference. This receipt binds both closure files and does not hash itself.'}

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--workers',choices=[4,8],type=int,default=4);args=parser.parse_args()
    for fd in [1,2]:
        try:t=os.readlink('/proc/self/fd/'+str(fd))
        except OSError:continue
        if t.startswith('/'):
            p=Path(t.removesuffix(' (deleted)')).resolve()
            if p==S or S in p.parents:raise RuntimeError('Verification stdout/stderr must be outside the stage')
    prior=(S/RECEIPT).read_bytes() if (S/RECEIPT).exists() else None
    try:result=run(args.workers);code=0
    except Exception as e:result={'status':'FAIL','verification_mode':'INDEPENDENT_FULL_REHASH_NO_CACHE','completed_utc':utc(),'error':str(e),'traceback':traceback.format_exc(),'hash_cache_used':False};code=1
    if prior is not None:
        old=json.loads(prior);history=old.pop('receipt_history',[]);history.append(old);result['receipt_history']=history;result['previous_receipt_sha256']=digest_bytes(prior)
    # This explicit closure path is excluded from payload to avoid a self-hash
    # cycle. Preserve prior receipt content on retries instead of discarding it.
    (S/RECEIPT).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','completed_utc','verification_mode']},indent=2),flush=True)
    if code:print(result['error'],file=sys.stderr)
    raise SystemExit(code)
