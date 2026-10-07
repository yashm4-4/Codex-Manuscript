from pathlib import Path
import pandas as pd,gzip,csv,json,datetime
r=Path(__file__).resolve().parents[1];d=pd.read_csv(r/'results/duplicate_identities.tsv',sep='\t',dtype=str);ids=set(d.SNPID);found=[]
if ids:
 with gzip.open(r/'acquisition/COPD.auto.rsq07.mac10.txt.gz','rt') as f:
  header=next(f).split()
  for sr,line in enumerate(f,1):
   arr=line.split()
   if arr[2] in ids:found.append(dict(zip(header,arr))|{'source_row':sr})
x=pd.DataFrame(found,columns=header+['source_row']);x.to_csv(r/'results/duplicate_source_rows.tsv',sep='\t',index=False)
summary=[]
for k,g in x.groupby('SNPID'):
 unique=g.drop(columns='source_row').drop_duplicates();summary.append({'variant_identity':k,'source_rows':';'.join(map(str,g.source_row)),'rows':len(g),'exact_duplicate_all_fields':len(unique)==1,'unique_row_patterns':len(unique),'in_prospective_locus':False,'disposition':'Outside prospectively selected loci; all original rows retained; any future locus would defer conflicting duplicate identity or retain first identical row with exclusion ledger'})
(r/'results/duplicate_conflicts.json').write_text(json.dumps({'groups':summary,'scan_completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2)+'\n');print(json.dumps(summary,indent=2))
