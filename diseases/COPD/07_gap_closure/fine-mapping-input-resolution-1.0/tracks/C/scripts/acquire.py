import requests,json,hashlib,datetime,pathlib,time
BASE=pathlib.Path(__file__).resolve().parents[1]
URL='https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_flat_files/icd10-J44-both_sexes.tsv.bgz'
p=BASE/'data/icd10-J44-both_sexes.tsv.bgz'
a={'url':URL,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'expected_bytes':1827545082,'expected_md5':'787f3a8ce3502a9a8f597a17144e07ae'}
(BASE/'sources/acquisition_started.json').write_text(json.dumps(a,indent=2)+'\n')
r=requests.get(URL,stream=True,timeout=(30,180),headers={'Accept-Encoding':'identity','User-Agent':'COPD-input-resolution-1.0'})
a.update(status_code=r.status_code,headers=dict(r.headers),final_url=r.url);r.raise_for_status()
sha=hashlib.sha256(); md5=hashlib.md5();n=0;t=time.time()
with p.open('wb') as f:
 for chunk in r.iter_content(4*1024*1024):
  f.write(chunk);sha.update(chunk);md5.update(chunk);n+=len(chunk)
  if n//(100*1024**2)!=(n-len(chunk))//(100*1024**2):print(n,'bytes',round(time.time()-t,1),'seconds',flush=True)
a.update(actual_bytes=n,sha256=sha.hexdigest(),md5=md5.hexdigest(),finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),path=str(p.relative_to(BASE)))
a['integrity_pass']=n==a['expected_bytes'] and a['md5']==a['expected_md5']
(BASE/'sources/acquisition.json').write_text(json.dumps(a,indent=2)+'\n')
print(json.dumps(a,indent=2),flush=True)
