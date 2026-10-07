"""Acquire only index partitions overlapping the authorized GWAS-defined loci."""
from pathlib import Path
import concurrent.futures
import csv
import gzip
import hashlib
import json
import sys
from public_fetch import fetch, utc

STAGE = Path(__file__).resolve().parents[1]
BASE = 'https://pan-ukb-us-east-1.s3.amazonaws.com/ld_release/UKBB.EUR.ldadj.variant.ht/'
HT = STAGE/'ld/source/UKBB.EUR.ldadj.variant.ht'


def loci_for(track):
    path = STAGE/'tracks'/track/'results/loci.tsv'
    rows = list(csv.DictReader(path.open(), delimiter='\t'))
    out = []
    for r in rows:
        chrom = r.get('chrom', r.get('chromosome'))
        start = r.get('start', r.get('start1', r.get('start_bp')))
        end = r.get('end', r.get('end1', r.get('end_bp')))
        if any('MHC' in str(v).upper() and ('DEFER' in str(v).upper() or 'EXCLUD' in str(v).upper()) for v in r.values()):
            continue
        out.append({'track': track, 'locus_id': r['locus_id'], 'chrom': str(chrom),
                    'start': int(start), 'end': int(end)})
    return out


def pos(x):
    chrom = x['locus']['contig']
    return (int(chrom) if chrom.isdigit() else {'X':23,'Y':24,'MT':25}.get(chrom,100),
            x['locus']['position'])


if __name__ == '__main__':
    track = sys.argv[1]
    loci = loci_for(track)
    meta = json.loads(gzip.decompress((HT/'rows/metadata.json.gz').read_bytes()))
    selected = []
    for i, (part, bounds) in enumerate(zip(meta['_partFiles'], meta['_jRangeBounds'])):
        lo, hi = pos(bounds['start']), pos(bounds['end'])
        if any(lo <= (int(l['chrom']), l['end']) and hi >= (int(l['chrom']), l['start']) for l in loci):
            selected.append((i, part))
    plan = {'created_utc': utc(), 'track': track, 'loci': loci, 'partitions': selected,
            'input_loci_sha256': hashlib.sha256((STAGE/'tracks'/track/'results/loci.tsv').read_bytes()).hexdigest(),
            'entire_ld_release_requested': False}
    out = STAGE/'ld'/track; out.mkdir(exist_ok=True)
    (out/'index_partition_plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    jobs = []
    for _, part in selected:
        for rel in ['rows/parts/'+part, 'index/'+part+'.idx/index', 'index/'+part+'.idx/metadata.json.gz']:
            jobs.append((BASE+rel, HT/rel))
    print(f'{track}: {len(loci)} non-MHC loci, {len(selected)} index partitions, {len(jobs)} objects',flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=12) as pool:
        futs = [pool.submit(fetch,*item) for item in jobs]
        errors = []
        for i, f in enumerate(concurrent.futures.as_completed(futs),1):
            try: f.result()
            except Exception as e: errors.append(str(e))
            if i % 100 == 0: print(f'index objects completed {i}/{len(jobs)}',flush=True)
    (out/'index_acquisition_result.json').write_text(json.dumps({'errors':errors,'objects':len(jobs)},indent=2)+'\n')
    if errors: raise RuntimeError(errors)
