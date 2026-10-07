from pathlib import Path
import requests,hashlib,json,datetime,concurrent.futures,tarfile
r=Path(__file__).resolve().parent
s={
'R-4.4.3.tar.gz':'https://cran.r-project.org/src/base/R-4/R-4.4.3.tar.gz',
'R443_release.html':'https://stat.ethz.ch/pipermail/r-announce/2025/000708.html',
'susieR_0.16.6.tar.gz':'https://cran.r-project.org/src/contrib/Archive/susieR/susieR_0.16.6.tar.gz',
'finemap_v1.4.2_x86_64.tgz':'https://www.christianbenner.com/finemap_v1.4.2_x86_64.tgz',
'finemap_documentation.html':'https://www.christianbenner.com/',
'lapack-3.12.1.tar.gz':'https://github.com/Reference-LAPACK/lapack/archive/refs/tags/v3.12.1.tar.gz',
'lapack_release.html':'https://www.netlib.org/lapack/lapack-3.12.1.html'
}
def f(item):
 k,u=item;o={'url':u,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Future-software source/archive/docs only; no program fitting or inference execution'}
 try:
  x=requests.get(u,timeout=(60,180));(r/k).write_bytes(x.content);o.update(status=x.status_code,final_url=x.url,headers=dict(x.headers),bytes=len(x.content),sha256=hashlib.sha256(x.content).hexdigest());x.raise_for_status()
  if k.endswith(('.tar.gz','.tgz')):
   with tarfile.open(r/k) as t:o['archive_members']=len(t.getmembers());o['tar_readable']=True
 except Exception as e:o['error']=repr(e)
 o['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();(r/(k+'.access.json')).write_text(json.dumps(o,indent=2)+'\n');return k,o.get('status'),o.get('bytes'),o.get('error')
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as e:
 for z in e.map(f,s.items()):print(z,flush=True)
