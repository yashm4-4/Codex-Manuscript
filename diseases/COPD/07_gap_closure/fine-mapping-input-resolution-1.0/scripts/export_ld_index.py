"""Use Hail's native decoder on the original downloaded locus index partitions."""
from pathlib import Path
import json
import sys
import hail as hl

s = Path(__file__).resolve().parents[1]
track = sys.argv[1]
plan = json.loads((s/'ld'/track/'index_partition_plan.json').read_text())
hl.init(master='local[4]',log=str(s/'provenance'/f'hail_index_{track}.log'),quiet=True,
        spark_conf={'spark.driver.host':'127.0.0.1','spark.driver.bindAddress':'127.0.0.1',
                    'spark.driver.memory':'12g'})
ht = hl.read_table(str(s/'ld/source/UKBB.EUR.ldadj.variant.ht'))
intervals = [hl.Interval(hl.Locus(l['chrom'],l['start'],'GRCh37'),
                         hl.Locus(l['chrom'],l['end'],'GRCh37'),True,True) for l in plan['loci']]
ht = hl.filter_intervals(ht, intervals)
ht = ht.key_by()
ht = ht.select(chrom=ht.locus.contig,pos=ht.locus.position,
                        ref=ht.alleles[0],alt=ht.alleles[1],idx=ht.idx,af=ht.AF,rsid=ht.rsid)
df = ht.to_pandas()
assert not df.idx.duplicated().any()
out=[]
for l in plan['loci']:
    part=df[(df.chrom==l['chrom']) & (df.pos>=l['start']) & (df.pos<=l['end'])].copy()
    part['locus_id']=l['locus_id'];part['track']=track
    part.sort_values('idx',inplace=True)
    part.to_csv(s/'ld'/track/(l['locus_id']+'.index.tsv.gz'),sep='\t',index=False,
                compression={'method':'gzip','mtime':0})
    out.append(dict(l,index_variants=len(part),min_idx=int(part.idx.min()) if len(part) else None,
                    max_idx=int(part.idx.max()) if len(part) else None))
(s/'ld'/track/'index_export_summary.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2));hl.stop()
