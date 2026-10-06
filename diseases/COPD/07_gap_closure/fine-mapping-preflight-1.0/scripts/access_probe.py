#!/usr/bin/env python3
"""Bounded public HTTP acquisition. Prefix hashes are never whole-object hashes."""
import datetime, hashlib, json, pathlib, urllib.error, urllib.request
def probe(url, destination, limit=262144, ranged=False, timeout=30):
    p=pathlib.Path(destination);p.parent.mkdir(parents=True,exist_ok=True)
    receipt=p.with_name(p.name+'.access.json')
    if receipt.exists():return json.loads(receipt.read_text())
    meta={'requested_url':url,'accessed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
          'method':'GET','anonymous':True,'read_limit_bytes':limit,'requested_range':f'bytes=0-{limit-1}' if ranged else None}
    headers={'User-Agent':'COPD-fine-mapping-preflight/1.0 (metadata and bounded schema audit)','Accept-Encoding':'identity'}
    if ranged:headers['Range']=f'bytes=0-{limit-1}'
    try:
        response=urllib.request.urlopen(urllib.request.Request(url,headers=headers),timeout=timeout)
    except urllib.error.HTTPError as e:response=e
    except Exception as e:
        meta.update({'status':None,'error_type':type(e).__name__,'error':str(e),'downloaded_bytes':0})
        receipt.write_text(json.dumps(meta,indent=2)+'\n');return meta
    with response:
        body=response.read(limit)
        meta.update({'status':response.status,'final_url':response.geturl(),
          'response_headers':{k.lower():v for k,v in response.headers.items() if k.lower() in
              {'content-length','content-range','content-type','etag','last-modified','accept-ranges','location'}},
          'downloaded_bytes':len(body),'saved_path':str(p),'sha256':hashlib.sha256(body).hexdigest(),
          'bounded_read':len(body)==limit,'hash_scope':'saved response bytes only; range/prefix is not a whole-file checksum'})
    p.write_bytes(body);receipt.write_text(json.dumps(meta,indent=2)+'\n');return meta
if __name__=='__main__':
    import argparse
    a=argparse.ArgumentParser();a.add_argument('url');a.add_argument('destination');a.add_argument('--limit',type=int,default=262144);a.add_argument('--range',action='store_true');x=a.parse_args()
    print(json.dumps(probe(x.url,x.destination,x.limit,x.range),indent=2))
