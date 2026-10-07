from pathlib import Path
import requests,hashlib,json,datetime,concurrent.futures
ROOT=Path(__file__).resolve().parents[1]
BASE='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90016001-GCST90017000/GCST90016588/'
objects=['GCST90016588_buildGRCh37.tsv','GCST90016588_buildGRCh37.tsv-meta.yaml','md5sum.txt','harmonised/33106845-GCST90016588-EFO_0006527.h.tsv.gz','harmonised/33106845-GCST90016588-EFO_0006527.h.tsv.gz-meta.yaml','harmonised/md5sum.txt']
def get(obj):
 p=ROOT/'sources'/obj;p.parent.mkdir(parents=True,exist_ok=True)
 rec={'url':BASE+obj,'accessed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'complete HTTP body','path':str(p.relative_to(ROOT))}
 try:
  r=requests.get(rec['url'],stream=True,timeout=(30,180));rec.update(status=r.status_code,final_url=r.url,response_headers=dict(r.headers));r.raise_for_status()
  h=hashlib.sha256();m=hashlib.md5();n=0
  with p.open('wb') as f:
   for b in r.iter_content(1024*1024):f.write(b);h.update(b);m.update(b);n+=len(b)
  rec.update(bytes=n,sha256=h.hexdigest(),md5=m.hexdigest(),completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),complete_http=True)
 except Exception as e:rec.update(error=repr(e),complete_http=False)
 p.with_name(p.name+'.access.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec),flush=True);return rec
with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:out=list(ex.map(get,objects))
(ROOT/'results/acquisition.json').write_text(json.dumps(out,indent=2)+'\n')
