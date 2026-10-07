from pathlib import Path
import requests,json,hashlib,datetime
R=Path(__file__).resolve().parents[1]
name='33106845-GCST90016588-EFO_0006527-Build37.f.tsv.gz';url='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90016001-GCST90017000/GCST90016588/harmonised/'+name;p=R/'sources/harmonised'/name
rec={'url':url,'accessed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'purpose':'Independent source precision/serialization lineage audit; frozen numeric rule unchanged; native failures retained','scope':'complete HTTP body','path':str(p.relative_to(R))}
try:
 r=requests.get(url,stream=True,timeout=(30,180));rec.update(status=r.status_code,headers=dict(r.headers));r.raise_for_status();h=hashlib.sha256();m=hashlib.md5();n=0
 with p.open('wb') as f:
  for b in r.iter_content(1024*1024):f.write(b);h.update(b);m.update(b);n+=len(b)
 expected=[l.split()[0] for l in (R/'sources/harmonised/md5sum.txt').read_text().splitlines() if l.rstrip().endswith(name)][0]
 rec.update(bytes=n,sha256=h.hexdigest(),md5=m.hexdigest(),expected_md5=expected,provider_md5_pass=m.hexdigest()==expected,complete_http=True,completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
except Exception as e:rec.update(error=repr(e),complete_http=False)
Path(str(p)+'.access.json').write_text(json.dumps(rec,indent=2)+'\n');print(json.dumps(rec,indent=2))
