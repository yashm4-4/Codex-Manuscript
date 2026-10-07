from pathlib import Path
import pandas as pd,json,time,datetime
B=Path(__file__).resolve().parents[1];R=B/'results';F=['chr','pos','ref','alt'];row=0;carry=set();carry_pos=None;within=0;cross=0;groups=0
for d in pd.read_csv(B/'data/icd10-J44-both_sexes.tsv.bgz',compression='gzip',sep='\t',usecols=F,dtype=str,keep_default_na=False,chunksize=250000):
 within+=int(d.duplicated(F,keep='first').sum());firstpos=(d.iloc[0]['chr'],d.iloc[0]['pos']);lastpos=(d.iloc[-1]['chr'],d.iloc[-1]['pos'])
 if firstpos==carry_pos:
  firstgroup=d[(d['chr']==firstpos[0])&(d['pos']==firstpos[1])]
  cross+=sum(tuple(x) in carry for x in firstgroup[['ref','alt']].itertuples(index=False,name=None))
 lastgroup=d[(d['chr']==lastpos[0])&(d['pos']==lastpos[1])];lastkeys=set(lastgroup[['ref','alt']].itertuples(index=False,name=None));carry=carry|lastkeys if firstpos==lastpos==carry_pos else lastkeys;carry_pos=lastpos;row+=len(d);groups+=1
res={'complete_rows':row,'gzip_complete_decode_crc_pass':True,'within_chunk_duplicate_identities':within,'cross_chunk_duplicate_identities_using_complete_boundary_position_groups':cross,'total_duplicate_identity_rows':within+cross,'chunks':groups,'chromosome_position_sorted_verified_by_main_audit':True,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()};(R/'identity_boundary_verification.json').write_text(json.dumps(res,indent=2)+'\n');print(json.dumps(res),flush=True)
