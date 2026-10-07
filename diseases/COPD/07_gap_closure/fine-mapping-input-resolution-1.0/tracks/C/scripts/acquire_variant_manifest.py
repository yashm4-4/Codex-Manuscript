import requests,json,hashlib,datetime,pathlib,time
BASE=pathlib.Path(__file__).resolve().parents[1]
u='https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/full_variant_qc_metrics.txt.bgz'
p=BASE/'data/full_variant_qc_metrics.txt.bgz';a={'url':u,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'purpose':'Variant identity/INFO/QC dependency only; gene annotations will not be parsed or used','expected_bytes':2701503051,'expected_md5_from_single_part_s3_etag':'e70ebc8289f762dd8d5086f54e766654'}
(BASE/'sources/variant_manifest.started.json').write_text(json.dumps(a,indent=2)+'\n')
r=requests.get(u,stream=True,timeout=(30,180));a.update(status=r.status_code,headers=dict(r.headers));r.raise_for_status();h=hashlib.sha256();m=hashlib.md5();n=0;t=time.time()
with p.open('wb') as f:
 for c in r.iter_content(4*1024*1024):
  f.write(c);h.update(c);m.update(c);n+=len(c)
  if n//(200*1024**2)!=(n-len(c))//(200*1024**2):print(n,round(time.time()-t,1),flush=True)
a.update(bytes=n,sha256=h.hexdigest(),md5=m.hexdigest(),finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());a['integrity_pass']=n==a['expected_bytes'] and a['md5']==a['expected_md5_from_single_part_s3_etag'];(BASE/'sources/variant_manifest.acquisition.json').write_text(json.dumps(a,indent=2)+'\n')
print(json.dumps(a,indent=2),flush=True)
