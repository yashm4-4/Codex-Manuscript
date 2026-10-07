from pathlib import Path
import requests,json,datetime,hashlib,zipfile,shutil
r=Path(__file__).resolve().parents[1];p=r/'acquisition/hum0014.v17.COPD.v1.zip';u='https://humandbs.dbcls.jp/files/hum0014/hum0014.v17.COPD.v1.zip';expected=2244311450
# Preserve failed receipt; partial bytes remain as the exact prefix in the completed artifact.
old=json.loads((r/'acquisition/receipt.json').read_text());(r/'acquisition/initial_incomplete_receipt.json').write_text(json.dumps(old,indent=2)+'\n')
receipts=[]
while p.stat().st_size<expected:
 start=p.stat().st_size;end=min(start+64*1024*1024-1,expected-1);o={'url':u,'range':[start,end],'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 for attempt in range(3):
  try:
   resp=requests.get(u,headers={'Range':f'bytes={start}-{end}','Accept-Encoding':'identity'},timeout=(60,180));o.update(status=resp.status_code,headers=dict(resp.headers),bytes=len(resp.content),sha256=hashlib.sha256(resp.content).hexdigest());assert resp.status_code==206;assert resp.headers['Content-Range']==f'bytes {start}-{end}/{expected}';assert len(resp.content)==end-start+1;assert resp.headers.get('Etag')==old['headers']['Etag'];break
  except Exception as e:
   o.setdefault('errors',[]).append(repr(e))
   if attempt==2:raise
 with p.open('ab') as f:f.write(resp.content)
 o['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();receipts.append(o);(r/'acquisition/range_resume_receipts.json').write_text(json.dumps(receipts,indent=2)+'\n');print('resumed',p.stat().st_size,flush=True)
h=hashlib.sha256()
with p.open('rb') as f:
 for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
o={'url':u,'bytes':p.stat().st_size,'sha256':h.hexdigest(),'full_size_pass':p.stat().st_size==expected,'acquisition':'initial full-response prefix + byte-exact HTTP206 resume; failures preserved','initial_incomplete_receipt':'initial_incomplete_receipt.json','range_receipts':'range_resume_receipts.json'}
with zipfile.ZipFile(p) as z:
 o['archive_members']=[{'name':i.filename,'bytes':i.file_size,'compressed_bytes':i.compress_size,'crc32':f'{i.CRC:08x}'} for i in z.infolist()];bad=z.testzip();o['all_zip_members_crc32_pass']=bad is None;assert bad is None
 member='hum0014.v17.COPD.v1/COPD.auto.rsq07.mac10.txt.gz';dest=r/'acquisition/COPD.auto.rsq07.mac10.txt.gz'
 with z.open(member) as fi,dest.open('wb') as fo:shutil.copyfileobj(fi,fo,4*1024*1024)
 h=hashlib.sha256()
 with dest.open('rb') as f:
  for b in iter(lambda:f.read(4*1024*1024),b''):h.update(b)
 o['analyzed_member']={'name':member,'bytes':dest.stat().st_size,'sha256':h.hexdigest(),'scope':'Combined autosomal only; sex-specific and chrX contents not scientifically analyzed'}
o['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();(r/'acquisition/receipt.json').write_text(json.dumps(o,indent=2)+'\n');print('complete',json.dumps(o),flush=True)
