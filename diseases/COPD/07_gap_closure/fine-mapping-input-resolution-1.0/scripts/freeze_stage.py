"""Freeze the completed stage; --prehash never creates a freeze.

Hierarchy: every scientific/source/code/provenance payload (including the static
hash cache) -> SHA256 ledger -> freeze record -> independent verification receipt.
Only the three explicit closure paths and inventoried __pycache__/*.pyc are
excluded. Final freeze never rewrites the cache. Stdout/stderr logs must be
outside this stage. No statistical analysis or storage transformation occurs.
"""
import sys
sys.dont_write_bytecode=True
from pathlib import Path
import argparse,concurrent.futures,csv,datetime,hashlib,json,os,stat,subprocess

S=Path(__file__).resolve().parents[1];REPO=S.parents[3]
LEDGER='provenance/artifact_checksums.tsv';FREEZE='provenance/freeze.json';VERIFY='provenance/freeze_verification.json'
CLOSURE={LEDGER,FREEZE,VERIFY};CACHE='provenance/payload_hash_cache.json'
OVERSIZE='provenance/oversized_artifact_inventory.tsv';EPHEMERAL='provenance/excluded_ephemeral_inventory.tsv'
MIN_BIG=100*1024*1024

def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def sha_bytes(b):return hashlib.sha256(b).hexdigest()
def ephemeral(path):
    p=Path(path);return p.name.endswith('.pyc') and p.parent.name=='__pycache__'
def guard_streams():
    for fd in [1,2]:
        try:t=os.readlink('/proc/self/fd/'+str(fd))
        except OSError:continue
        if t.startswith('/'):
            p=Path(t.removesuffix(' (deleted)')).resolve()
            if p==S or S in p.parents:raise RuntimeError('Redirect stdout/stderr outside the stage; logs cannot mutate bound payload: '+str(p))
def signature(p):
    x=p.lstat();return {'bytes':x.st_size,'mtime_ns':x.st_mtime_ns,'ctime_ns':x.st_ctime_ns,'inode':x.st_ino,'device':x.st_dev,'mode':x.st_mode}
def inventory():
    result={}
    for base,dirs,files in os.walk(S,followlinks=False):
        for name in list(dirs):
            p=Path(base)/name
            if p.is_symlink():result[str(p.relative_to(S))]=signature(p);dirs.remove(name)
        for name in files:
            p=Path(base)/name;q=signature(p)
            if not (stat.S_ISREG(q['mode']) or stat.S_ISLNK(q['mode'])):raise RuntimeError('Unsupported non-regular filesystem artifact: '+str(p))
            result[str(p.relative_to(S))]=q
    return result

def hash_artifact(rel,expected=None):
    p=S/rel;before=signature(p)
    if expected is not None and before!=expected:raise RuntimeError('Changed before hashing: '+rel)
    h=hashlib.sha256();n=0
    if stat.S_ISLNK(before['mode']):
        b=os.readlink(p).encode('utf-8');h.update(b);n=len(b);kind='symlink_target_utf8'
    else:
        kind='regular_file'
        with p.open('rb') as f:
            for b in iter(lambda:f.read(16*1024*1024),b''):h.update(b);n+=len(b)
    if signature(p)!=before:raise RuntimeError('Changed during hashing: '+rel)
    if kind=='regular_file' and n!=before['bytes']:raise RuntimeError('Byte count mismatch: '+rel)
    return {'path':rel,'bytes':n,'sha256':h.hexdigest(),'kind':kind,'stat':before,'hashed_utc':utc()}

def atomic_json(p,obj):
    # The temporary file is removed before final payload collection. Never invoke
    # this after final snapshots except for an explicitly excluded closure path.
    tmp=p.with_name(p.name+'.tmp');tmp.write_text(json.dumps(obj,indent=2)+'\n');os.replace(tmp,p)
def write_tsv(path,rows,columns):
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,delimiter='\t',fieldnames=columns);w.writeheader();w.writerows(rows)
def classification(rel):
    low=rel.lower()
    if low.endswith(('.npy','.npz')) and '/ld/' in '/'+low:return 'diagnostic_matrix'
    if low.startswith('ld/source/') or '/data/' in '/'+low or '/sources/' in '/'+low or '/raw/' in '/'+low:return 'unchanged_local_original_source'
    if any(x in low for x in ['toolchain','software','jdk','r-4.','susier','finemap']) and low.endswith(('.gz','.xz','.zip','.tar','.bz2')):return 'toolchain_or_future_software_archive_context'
    if low.startswith('references/'):return 'reference_asset_original_or_explicit_derived'
    if any(x in low for x in ['attempt','failed','failure']):return 'preserved_failed_attempt_or_context'
    if low.startswith('tracks/'):return 'derived_gwas_audit_or_context'
    return 'archive_or_other_context'
def write_inventories(phase):
    inv=inventory();big=[];excluded=[]
    for rel,sig in sorted(inv.items()):
        if ephemeral(rel):
            q=hash_artifact(rel,sig);excluded.append({'path':rel,'bytes':q['bytes'],'sha256_at_inventory':q['sha256'],'reason':'Interpreter bytecode only: immediate parent __pycache__, suffix .pyc','inventory_phase':phase})
        elif rel not in CLOSURE and sig['bytes']>=MIN_BIG:
            big.append({'path':rel,'bytes':sig['bytes'],'MiB':format(sig['bytes']/1024**2,'.6f'),'classification':classification(rel),'storage_disposition':'RETAINED_UNCHANGED_LOCALLY; no publication/storage transform/deletion performed','checksum_location':LEDGER,'inventory_phase':phase})
    write_tsv(S/OVERSIZE,big,['path','bytes','MiB','classification','storage_disposition','checksum_location','inventory_phase'])
    write_tsv(S/EPHEMERAL,excluded,['path','bytes','sha256_at_inventory','reason','inventory_phase'])
    return {'oversized_files':len(big),'oversized_bytes':sum(x['bytes'] for x in big),'excluded_ephemeral_files':len(excluded)}

def cached_entries():
    p=S/CACHE
    if not p.exists():return {}
    obj=json.loads(p.read_text());return obj.get('entries',{})
def cache_matches(q,sig):
    return isinstance(q,dict) and q.get('stat')==sig and len(q.get('sha256',''))==64 and q.get('kind') in ['regular_file','symlink_target_utf8']
def prehash(workers):
    if (S/FREEZE).exists():raise RuntimeError('Stage is already frozen; prehash would mutate bound cache/inventories')
    info=write_inventories('PRELIMINARY');before=inventory();old=cached_entries();new={};errors=[];todo=[]
    for rel,sig in before.items():
        if rel in CLOSURE or rel==CACHE or ephemeral(rel):continue
        if cache_matches(old.get(rel),sig):new[rel]=old[rel]
        else:todo.append((rel,sig))
    print(json.dumps({'mode':'prehash','files_to_hash':len(todo),'reused_stable_entries':len(new),'bytes_to_hash':sum(x[1]['bytes'] for x in todo),'workers':workers,**info}),flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as ex:
        futures={ex.submit(hash_artifact,rel,sig):rel for rel,sig in todo}
        for i,fu in enumerate(concurrent.futures.as_completed(futures),1):
            rel=futures[fu]
            try:new[rel]=fu.result()
            except Exception as e:errors.append({'path':rel,'error':str(e)})
            if i%100==0:print('Prehashed',i,'of',len(todo),flush=True)
    after=inventory();stable={rel:q for rel,q in new.items() if after.get(rel)==q['stat']}
    for rel in sorted(set(new)-set(stable)):errors.append({'path':rel,'error':'Changed or disappeared after hashing; excluded from cache'})
    atomic_json(S/CACHE,{'format':'stage_payload_hash_cache_v1','authority':'Performance hint only; final freeze stat-checks entries and independent verification ignores all caches','self_inclusion':'Cache never caches itself. Final ledger hashes this static cache file as ordinary payload; final freeze never rewrites it.','created_utc':utc(),'entries':dict(sorted(stable.items())),'unstable_or_failed_entries':errors,'new_paths_after_snapshot':sorted(set(after)-set(before))})
    result={'mode':'prehash','cached_stable_files':len(stable),'cached_stable_bytes':sum(q['bytes'] for q in stable.values()),'unstable_or_failed_entries':len(errors),'cache_path':CACHE,'freeze_created':False}
    print(json.dumps(result,indent=2),flush=True)

def stage_path(value):
    p=(S/value).resolve()
    if p==S or S not in p.parents:raise ValueError('Required artifact must be inside current stage: '+value)
    if not p.is_file() or p.stat().st_size==0:raise RuntimeError('Required nonempty artifact missing: '+value)
    rel=str(p.relative_to(S))
    if rel in CLOSURE or ephemeral(rel):raise RuntimeError('Required scientific/validation artifact cannot be excluded: '+rel)
    return p,rel

def prerequisites(args):
    required={}
    for key in ['integrity','extraction','report','review']:
        p,rel=stage_path(getattr(args,key));required[key]={'path':rel,'sha256':hash_artifact(rel)['sha256']};required[key]['bytes']=p.stat().st_size
    integrity=json.loads((S/required['integrity']['path']).read_text())
    assert integrity['status']=='PASS' and integrity.get('checks_failed',1)==0
    assert all(x.get('pass') is True for x in integrity.get('checks',[])) and integrity.get('checks_passed',0)>0
    extraction=json.loads((S/required['extraction']['path']).read_text())
    assert extraction['mode']=='final' and extraction['status']=='PASS_EXTRACTION_AND_OUTPUT_VALIDATION_ONLY'
    assert extraction['planned_loci']==16 and extraction['check_counts'].get('FAIL',0)==0 and extraction['check_counts'].get('PENDING',0)==0
    rp=S/required['review']['path']
    if rp.suffix.lower()=='.json':
        review=json.loads(rp.read_text());assert review.get('status') in ['PASS','PASS_FINAL_PRE_FREEZE_REVIEW'] and review.get('final') is True
    else:assert 'FREEZE_REVIEW_STATUS: PASS_FINAL' in [x.strip() for x in rp.read_text().splitlines()]
    assert rp.stat().st_mtime_ns>=max((S/required[k]['path']).stat().st_mtime_ns for k in ['integrity','extraction','report']), 'Independent final review must follow final report/integrity/extraction artifacts'
    return required

def register_records(init):
    result=[]
    for r in init['baseline_registers']:
        p=REPO/r['path'];sig=signature(p);b=p.read_bytes();assert signature(p)==sig
        base=r['bytes'];prefix=b[:base]
        assert sha_bytes(prefix)==r['sha256'] and len(b)>base and b'fine-mapping-input-resolution-1.0' in b[base:]
        snapshot=S/'provenance'/('baseline_'+p.name);assert snapshot.read_bytes()==prefix
        rows=list(csv.reader(b.decode().splitlines(),delimiter='\t'))
        assert all(len(x)==len(rows[0]) for x in rows) and len({x[0] for x in rows[1:]})==len(rows)-1
        result.append({'path':r['path'],'bytes':len(b),'sha256':sha_bytes(b),'baseline_prefix_bytes':base,'baseline_prefix_sha256':r['sha256'],'baseline_snapshot':str(snapshot.relative_to(S)),'appended_bytes':len(b)-base,'appended_rows':b[base:].count(b'\n'),'append_only':True,'stat_at_freeze':sig})
    return result

def do_freeze(args):
    if (S/FREEZE).exists() or (S/LEDGER).exists():raise RuntimeError('Closure already exists; do not overwrite an existing freeze or incomplete ledger')
    required=prerequisites(args);init=json.loads((S/'provenance/initialization.json').read_text())
    assert subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip()==init['baseline_commit']
    regs=register_records(init);info=write_inventories('FINAL')
    before=inventory();payload={r:q for r,q in before.items() if r not in CLOSURE and not ephemeral(r)}
    cache=cached_entries();entries={};todo=[];reused=0
    for rel,sig in payload.items():
        if rel!=CACHE and cache_matches(cache.get(rel),sig):entries[rel]=cache[rel];reused+=1
        else:todo.append((rel,sig))
    print(json.dumps({'mode':'freeze','payload_files':len(payload),'cached_stable_files':reused,'files_to_hash':len(todo),'workers':args.workers}),flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
        futures={ex.submit(hash_artifact,rel,sig):rel for rel,sig in todo}
        for fu in concurrent.futures.as_completed(futures):entries[futures[fu]]=fu.result()
    after=inventory()
    assert before==after,'Stage inventory or file metadata changed during final freeze; no closure written'
    assert register_records(init)==regs,'Registers changed during final freeze'
    for q in required.values():assert entries[q['path']]['sha256']==q['sha256'],'Required final artifact changed'
    # The cache itself is now an ordinary payload entry; no further cache writes.
    ledger_rows=[{k:entries[r][k] for k in ['path','bytes','sha256','kind']} for r in sorted(entries)]
    ledger_tmp=S/(LEDGER+'.tmp')
    write_tsv(ledger_tmp,ledger_rows,['path','bytes','sha256','kind'])
    os.replace(ledger_tmp,S/LEDGER)
    ledger_sha=sha_bytes((S/LEDGER).read_bytes())
    freeze={'stage':'fine-mapping-input-resolution-1.0','status':'FROZEN_INPUT_RESOLUTION_STOP_FOR_INVESTIGATOR_REVIEW','frozen_utc':utc(),
      'baseline_commit':init['baseline_commit'],'checksum_ledger':LEDGER,'checksum_ledger_sha256':ledger_sha,
      'payload_files':len(entries),'payload_bytes':sum(x['bytes'] for x in entries.values()),'required_final_artifacts':required,'registers':regs,
      'exclusions':{'closure_paths':sorted(CLOSURE),'ephemeral_rule':'Immediate parent directory __pycache__ AND filename suffix .pyc only','ephemeral_inventory':EPHEMERAL,'all_other_files_included':True},
      'closure_hierarchy':'Scientific/source/code/provenance payload including static hash cache -> SHA256 ledger -> freeze record -> full independent verification receipt. Ledger cannot hash itself; freeze binds ledger; receipt binds both ledger and freeze.',
      'hash_cache_policy':'Size/mtime/ctime/inode/device/mode-bound stable entries may accelerate closure; static cache is itself hashed in payload. Independent verification must rehash every payload without cache.',
      'oversized_artifact_inventory':OVERSIZE,**info,'scientific_execution_authorized':False,'statistical_fine_mapping_executed':False,'publication_authorized':False}
    atomic_json(S/FREEZE,freeze)
    final=inventory();expected=dict(before)
    for rel in CLOSURE:final.pop(rel,None);expected.pop(rel,None)
    assert final==expected,'Payload mutation detected immediately after closure; independent verification must fail'
    assert register_records(init)==regs,'Register mutation detected immediately after closure'
    print(json.dumps({'status':freeze['status'],'payload_files':len(entries),'payload_bytes':freeze['payload_bytes'],'ledger_sha256':ledger_sha,'freeze_sha256':sha_bytes((S/FREEZE).read_bytes()),'independent_full_rehash_required':True},indent=2),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);mode=p.add_mutually_exclusive_group(required=True);mode.add_argument('--prehash',action='store_true');mode.add_argument('--freeze',action='store_true');p.add_argument('--workers',type=int,choices=[4,8],default=4)
    p.add_argument('--integrity',default='provenance/stage_integrity_validation.json');p.add_argument('--extraction',default='ld/extraction_validation.json');p.add_argument('--report',default='FINE_MAPPING_INPUT_RESOLUTION_REPORT.md');p.add_argument('--review',default='provenance/independent_pre_freeze_review.md')
    args=p.parse_args();guard_streams()
    if args.prehash:prehash(args.workers)
    else:do_freeze(args)
