"""Acquire only Pan-UKB matrix blocks intersecting this track's non-MHC loci."""
from pathlib import Path
import concurrent.futures
import hashlib
import json
import math
import sys
from public_fetch import fetch, utc

s=Path(__file__).resolve().parents[1]
track=sys.argv[1]
base='https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.EUR.ldadj.bm/'
bm=s/'ld/source/UKBB.EUR.ldadj.bm'
meta=json.loads((bm/'metadata.json').read_text())
bs=meta['blockSize'];nb=math.ceil(meta['nRows']/bs)
part_by_id=dict(zip(meta['maybeFiltered'],meta['partFiles']))
index_summary=json.loads((s/'ld'/track/'index_export_summary.json').read_text())
plans=[];needed={}
for locus in index_summary:
    block_range=range(locus['min_idx']//bs,locus['max_idx']//bs+1)
    parts=[];absent=[]
    for bc in block_range:
        for br in block_range:
            if br>bc:continue
            bid=br+bc*nb
            if bid not in part_by_id:
                absent.append({'block_row':br,'block_col':bc,'block_id':bid});continue
            part=part_by_id[bid];needed[part]=bid
            parts.append({'block_row':br,'block_col':bc,'block_id':bid,'part':part})
    plans.append(dict(locus,parts=parts,absent_parts=absent))
plan={'created_utc':utc(),'track':track,'source_metadata_sha256':hashlib.sha256((bm/'metadata.json').read_bytes()).hexdigest(),
      'source_grid_block_size':bs,'matrix_rows':meta['nRows'],'source_upper_triangle':True,
      'loci':plans,'unique_blocks_requested':len(needed),'entire_release_requested':False}
(s/'ld'/track/'matrix_block_plan.json').write_text(json.dumps(plan,indent=2)+'\n')
fetch(base+'_SUCCESS',bm/'_SUCCESS')
print(f'{track}: acquiring {len(needed)} original compressed blocks for {len(plans)} GWAS loci only',flush=True)
errors=[];total=0
with concurrent.futures.ThreadPoolExecutor(max_workers=12) as pool:
    jobs={pool.submit(fetch,base+'parts/'+part,bm/'parts'/part):part for part in needed}
    for i,f in enumerate(concurrent.futures.as_completed(jobs),1):
        try:r=f.result();total+=r['bytes']
        except Exception as e:errors.append({'part':jobs[f],'error':str(e)})
        if i%10==0 or i==len(jobs):print(f'blocks {i}/{len(jobs)}, total completed bytes {total}, errors {len(errors)}',flush=True)
(s/'ld'/track/'matrix_acquisition_result.json').write_text(json.dumps({'completed_utc':utc(),'errors':errors,'objects':len(needed),'bytes':total},indent=2)+'\n')
if errors:raise RuntimeError(errors)
