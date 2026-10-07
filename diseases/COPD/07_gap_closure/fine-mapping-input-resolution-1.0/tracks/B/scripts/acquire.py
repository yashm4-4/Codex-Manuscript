from pathlib import Path
import requests, json, hashlib, datetime, zipfile, shutil, gzip
root=Path(__file__).resolve().parents[1]
u='https://humandbs.dbcls.jp/files/hum0014/hum0014.v17.COPD.v1.zip'
p=root/'acquisition/hum0014.v17.COPD.v1.zip'
r={'url':u,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Original distribution archive; only combined-sex autosomal member scientifically inspected. Other members CRC inventory only.'}
try:
 with requests.get(u,stream=True,timeout=(60,180)) as resp:
  r.update(status=resp.status_code,headers=dict(resp.headers),final_url=resp.url); resp.raise_for_status()
  assert resp.status_code==200, 'Full response required; do not accept partial range'
  h=hashlib.sha256(); n=0
  with p.open('wb') as f:
   for b in resp.iter_content(4*1024*1024):
    f.write(b); h.update(b); n+=len(b)
    if n%(128*1024*1024)<4*1024*1024: print(n,flush=True)
  r.update(bytes=n,sha256=h.hexdigest())
  if resp.headers.get('Content-Length'): assert n==int(resp.headers['Content-Length'])
 r['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
 with zipfile.ZipFile(p) as z:
  r['archive_members']=[{'name':i.filename,'bytes':i.file_size,'compressed_bytes':i.compress_size,'crc32':f'{i.CRC:08x}'} for i in z.infolist()]
  bad=z.testzip(); r['all_zip_members_crc32_pass']=bad is None; assert bad is None
  member='hum0014.v17.COPD.v1/COPD.auto.rsq07.mac10.txt.gz'; out=root/'acquisition/COPD.auto.rsq07.mac10.txt.gz'
  with z.open(member) as fi,out.open('wb') as fo: shutil.copyfileobj(fi,fo,4*1024*1024)
  h=hashlib.sha256()
  with out.open('rb') as f:
   for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
  r['analyzed_member']={'member':member,'file':str(out.relative_to(root)),'bytes':out.stat().st_size,'sha256':h.hexdigest()}
except Exception as e:
 r['error']=repr(e)
finally:
 (root/'acquisition/receipt.json').write_text(json.dumps(r,indent=2)+'\n')
 print(json.dumps(r,indent=2),flush=True)
