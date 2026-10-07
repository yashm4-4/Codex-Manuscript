from pathlib import Path
import requests,json,hashlib,datetime,concurrent.futures
r=Path(__file__).resolve().parents[1]/'sources'
sources={
'bbj_readme':'https://humandbs.dbcls.jp/files/hum0014/hum0014_v17_v18_v21_README.txt',
'bbj_sample_size.xlsx':'https://humandbs.dbcls.jp/files/hum0014/hum0014_v17_v18_v21_sample_size.xlsx',
'bbj_original_paper':'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/PMC7968075/unicode',
'bbj_available_genomes':'https://biobankjp.org/researchers/730',
'nbdc_reference_panels':'https://humandbs-production.ddbj.nig.ac.jp/en/processed-data-imputation',
'nbdc_jgad000873':'https://humandbs-production.ddbj.nig.ac.jp/dataset/JGAD000873',
'bbj_pheweb_downloads':'https://pheweb.jp/downloads',
'jmorp_genomics':'https://www.megabank.tohoku.ac.jp/english/research/genomics/',
'jmorp_ld_map':'https://www.megabank.tohoku.ac.jp/english/the-linkage-disequilibrium-map-tommo_ld_map_192v1-has-been-released/',
'saige_wiki':'https://github.com/weizhouUMICH/SAIGE/wiki/Genetic-association-tests-using-SAIGE',
'saige_original_paper':'https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/PMC6119127/unicode',
'saige_tags':'https://api.github.com/repos/weizhouUMICH/SAIGE/tags?per_page=100',
'saige_step2_v02942':'https://raw.githubusercontent.com/weizhouUMICH/SAIGE/0.29.4.2/R/SAIGE_SPATest.R'
}
def get(item):
 k,u=item;o={'url':u,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:
  x=requests.get(u,timeout=90);o.update(status=x.status_code,final_url=x.url,headers=dict(x.headers),bytes=len(x.content),sha256=hashlib.sha256(x.content).hexdigest());(r/(k+'.body')).write_bytes(x.content)
 except Exception as e:o['error']=repr(e)
 o['completed_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat();(r/(k+'.json')).write_text(json.dumps(o,indent=2)+'\n');return k,o.get('status'),o.get('bytes')
with concurrent.futures.ThreadPoolExecutor(max_workers=5) as ex:
 for result in ex.map(get,sources.items()):print(result,flush=True)
