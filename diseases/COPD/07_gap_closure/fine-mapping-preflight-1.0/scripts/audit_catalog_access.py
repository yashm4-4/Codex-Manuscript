#!/usr/bin/env python3
"""Metadata/directory/header-only audit for every frozen accession. No association analysis."""
import concurrent.futures, csv, html, io, json, pathlib, re, threading, time, urllib.parse, zlib
from access_probe import probe
STAGE=pathlib.Path(__file__).resolve().parents[1]
OUT=STAGE/'audits/catalog'
lock=threading.Lock();last_api=0.0
def api_probe(url,dest):
    global last_api
    if dest.with_name(dest.name+'.access.json').exists():return probe(url,dest)
    with lock:
        pause=max(0,0.8-(time.monotonic()-last_api))
        if pause:time.sleep(pause)
        last_api=time.monotonic()
    return probe(url,dest)
def links(path,base):
    if not path.exists():return []
    found=re.findall(r'href="([^"]+)"',path.read_text(errors='replace'))
    return [(html.unescape(x),urllib.parse.urljoin(base,html.unescape(x))) for x in found
            if not x.startswith(('/','?')) and x not in ('../','./')]
def sample_header(path):
    b=path.read_bytes()
    if b.startswith(b'\x1f\x8b'):
        try:
            remaining=b;chunks=[];total=0
            while remaining and total<1048576:
                decoder=zlib.decompressobj(16+zlib.MAX_WBITS)
                chunk=decoder.decompress(remaining,1048576-total)
                chunks.append(chunk);total+=len(chunk)
                remaining=decoder.unused_data
                if not decoder.eof:break
            b=b''.join(chunks)
        except zlib.error as e:return {'parse_error':str(e)}
    text=b.decode('utf-8-sig',errors='replace');lines=text.splitlines()
    if b'<' in b[:20] and ('<html' in text.lower() or '<!doctype' in text.lower()):return {'parse_error':'HTML response, not statistics'}
    while lines and lines[0].startswith('##'):lines.pop(0)
    if not lines:return {'parse_error':'no decoded lines'}
    sep='\t' if '\t' in lines[0] else ',' if ',' in lines[0] else None
    header=lines[0].lstrip('#').split(sep)
    records=[x.split(sep) for x in lines[1:9]]
    return {'header':header,'representative_records':records,'sample_records':len(records),
      'records_width_matches':all(len(x)==len(header) for x in records),
      'inspection_scope':'first header and up to eight records only; no locus enumeration, effect analysis or whole-file QC'}
def audit(row):
    accession=row['study_accession'];d=OUT/accession;d.mkdir(parents=True,exist_ok=True)
    summary={'accession':accession,'phenotype_class_frozen':row['phenotype_class'],'publication_pmid':row['pubmed_id'],'file_inventory':[],'headers':[],'metadata_files':[]}
    m=api_probe('https://www.ebi.ac.uk/gwas/rest/api/v2/studies/'+accession,d/'study.json');summary['catalog_api_access']=m
    current={}
    if m.get('status')==200:
        try:current=json.loads((d/'study.json').read_text())
        except Exception:pass
    summary['current_catalog_metadata']=current
    loc=current.get('full_summary_stats')
    if not isinstance(loc,str) or not loc.startswith(('http','ftp')):
        loc=row['catalog_full_stats_location']
    n=int(accession[4:]);lo=((n-1)//1000)*1000+1;hi=lo+999
    width=max(6,len(str(n)))
    guessed='https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST'+str(lo).zfill(width)+'-GCST'+str(hi).zfill(width)+'/'+accession+'/'
    if not loc.startswith(('http','ftp')):loc=guessed;summary['directory_location_basis']='standard accession directory probe, not proof of release'
    else:summary['directory_location_basis']='current or frozen Catalog advertised source'
    loc=loc.replace('http://ftp.ebi.ac.uk','https://ftp.ebi.ac.uk').replace('ftp://ftp.ebi.ac.uk','https://ftp.ebi.ac.uk').rstrip('/')+'/'
    summary['directory_url']=loc
    response=probe(loc,d/'directory.html');summary['directory_access']=response
    directories=[('',loc,d/'directory.html')] if response.get('status')==200 else []
    first_links=links(d/'directory.html',loc) if directories else []
    for name,url in first_links:
        if name.lower()=='harmonised/':
            dest=d/'harmonised_directory.html';receipt=probe(url,dest)
            if receipt.get('status')==200:directories.append(('harmonised_',url,dest))
    for prefix,base,listing in directories:
        for name,url in links(listing,base):
            if name.endswith('/'):continue
            summary['file_inventory'].append({'name':name,'url':url,'directory_evidence':str(listing.relative_to(STAGE))})
            if re.search(r'(meta\.ya?ml|\.ya?ml|md5sum\.txt|README.*|\.readme)$',name,re.I):
                dest=d/(prefix+urllib.parse.unquote(name));receipt=probe(url,dest)
                summary['metadata_files'].append({'path':str(dest.relative_to(STAGE)),'url':url,'status':receipt.get('status')})
            elif re.search(r'\.(tsv|txt|csv)(\.(gz|bgz))?$|\.h\.tsv\.gz$',name,re.I) and name.lower()!='md5sum.txt':
                dest=d/(prefix+urllib.parse.unquote(name)+'.prefix');receipt=probe(url,dest,limit=131072,ranged=True)
                parsed=sample_header(dest) if receipt.get('status') in (200,206) and dest.exists() else {'parse_error':'download failed'}
                parsed.update({'url':url,'prefix_path':str(dest.relative_to(STAGE)),'access_status':receipt.get('status'),
                               'prefix_sha256':receipt.get('sha256'),'whole_source_content_range':receipt.get('response_headers',{}).get('content-range'),
                               'whole_file_downloaded_and_validated':False})
                summary['headers'].append(parsed)
    (d/'audit.json').write_text(json.dumps(summary,indent=2)+'\n')
    return summary
def main():
    rows=list(csv.DictReader((STAGE/'inputs/frozen_study_metadata.tsv').open(),delimiter='\t'))
    OUT.mkdir(parents=True,exist_ok=True)
    results=[]
    with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
        futures={pool.submit(audit,r):r['study_accession'] for r in rows}
        for future in concurrent.futures.as_completed(futures):
            result=future.result();results.append(result)
            print(json.dumps({'accession':result['accession'],'api':result['catalog_api_access'].get('status'),'directory':result['directory_access'].get('status'),'headers':len(result['headers'])}),flush=True)
    results.sort(key=lambda x:x['accession'])
    (OUT/'catalog_audit.json').write_text(json.dumps(results,indent=2)+'\n')
    print(json.dumps({'completed':len(results),'header_audited_accessions':sum(bool(x['headers']) for x in results)}))
if __name__=='__main__':main()
