from pathlib import Path
import requests,hashlib,json,datetime,concurrent.futures,shutil
root=Path(__file__).resolve().parents[1]
urls={'catalog_harmonization.html':'https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics/','plink2_formats.html':'https://www.cog-genomics.org/plink/2.0/formats','hg19.chrom.sizes':'https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes','kim2021_pmc.html':'https://pmc.ncbi.nlm.nih.gov/articles/PMC8096488/'}
def run(item):
 name,url=item;rec={'url':url,'accessed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:
  r=requests.get(url,timeout=60);p=root/'sources'/name;p.write_bytes(r.content);rec.update(status=r.status_code,bytes=len(r.content),sha256=hashlib.sha256(r.content).hexdigest(),path=str(p.relative_to(root)),final_url=r.url)
 except Exception as e:rec['error']=repr(e)
 (root/'sources'/(name+'.access.json')).write_text(json.dumps(rec,indent=2)+'\n');print(rec)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:list(ex.map(run,urls.items()))
