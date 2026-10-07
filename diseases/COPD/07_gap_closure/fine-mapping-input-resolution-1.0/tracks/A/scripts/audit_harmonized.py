from pathlib import Path
import pandas as pd,numpy as np,json,collections,datetime,gzip,csv,math
R=Path(__file__).resolve().parents[1];S=R/'sources';O=R/'results';C=json.loads((R/'config/effect_contract.json').read_text())
path=S/'harmonised/33106845-GCST90016588-EFO_0006527.h.tsv.gz';rec=json.loads(Path(str(path)+'.access.json').read_text());assert rec['complete_http'] and rec['status']==200
provider=[l.split()[0] for l in (S/'harmonised/md5sum.txt').read_text().splitlines() if l.rstrip().endswith(path.name)][0];assert provider==rec['md5']
cols=['variant_id','p_value','effect_allele','other_allele','odds_ratio','standard_error']
d=pd.read_csv(S/'GCST90016588_buildGRCh37.tsv',sep='\t',dtype=str,keep_default_na=False);d.variant_id=d.variant_id.str.strip('"');mapping={};ambiguous={}
print('build exact transport key dictionary',len(d),flush=True)
for i,key in enumerate(d[cols].itertuples(index=False,name=None)):
 if key in mapping:ambiguous.setdefault(key,[mapping[key]]).append(i)
 else:mapping[key]=i
native_pos=d.base_pair_location.to_numpy();native_chr=d.chromosome.to_numpy();seen=np.zeros(len(d),dtype=np.int32);missing=collections.Counter();codes=collections.Counter();counts=collections.Counter();locusrows=set(pd.read_csv(O/'locus_variants.tsv.gz',sep='\t',usecols=['source_row']).source_row)
ledger=[];fails=[];minpos={};maxpos={};harmkeys=set();dupkeys=[];row=0
comp=str.maketrans('ACGT','TGCA')
for h in pd.read_csv(path,sep='\t',dtype=str,keep_default_na=False,chunksize=200000):
 h.variant_id=h.variant_id.str.strip('"');headers=list(h.columns)
 for col in headers:missing[col]+=int(h[col].str.strip('"').isin(['','NA','NaN','.','nan']).sum())
 codes.update(h.hm_code)
 for r in h.itertuples(index=False,name=None):
  row+=1;q=dict(zip(headers,r));key=tuple(q[k] for k in cols);idx=mapping.get(key)
  if idx is None:counts['native_transport_unmatched']+=1
  elif key in ambiguous:counts['native_transport_ambiguous']+=1
  else:seen[idx]+=1;counts['native_transport_exact']+=1
  code=q['hm_code'];flip=code in ['2','4','6','8','11','13'];reverse=code in ['3','4','7','8','12','13'];valid=code in ['1','2','3','4','5','6','7','8','10','11','12','13']
  status='UNORIENTED_CODE';opass=False;apass=False
  if valid:
   try:
    orig=float(q['odds_ratio']);hor=float(q['hm_odds_ratio']);expected=1/orig if flip else orig
    opass=math.isclose(hor,expected,rel_tol=C['reciprocal_or_relative_tolerance'],abs_tol=0) and abs(math.log(hor)-(-1 if flip else 1)*math.log(orig))<=C['reciprocal_log_or_absolute_tolerance']
    ea=q['effect_allele'];oa=q['other_allele']
    if reverse:ea=ea.translate(comp)[::-1];oa=oa.translate(comp)[::-1]
    if flip:ea,oa=oa,ea
    apass=ea==q['hm_effect_allele'] and oa==q['hm_other_allele']
    status='PASS' if opass and apass else 'FAIL'
   except Exception:status='INVALID_NUMERIC'
   counts['or_transform_'+str(opass)]+=1;counts['allele_transform_'+str(apass)]+=1
   if not (opass and apass):fails.append({'harmonized_row':row,'source_row':None if idx is None else idx+1,**q,'status':status})
  else:counts['unoriented_harmonization']+=1
  if idx is not None and key not in ambiguous:
   oldpos=native_pos[idx];oldchr=native_chr[idx];counts['coordinate_change' if oldpos!=q['base_pair_location'] or oldchr!=q['chromosome'] else 'coordinate_unchanged']+=1
   if idx+1 in locusrows:ledger.append({'source_row':idx+1,'harmonized_row':row,'native_chrom':oldchr,'native_pos':oldpos,'harmonized_build':'GRCh38',**q,'transformation_status':status,'source_transport_match':'EXACT_STATISTIC_AND_ALLELE_SIGNATURE_NOT_RSID_ONLY'})
  hk=(q['hm_chrom'],q['hm_pos'],q['hm_other_allele'],q['hm_effect_allele'])
  if any(v in ['NA','','.'] for v in hk):counts['missing_harmonized_identity']+=1
  elif hk in harmkeys:dupkeys.append({'harmonized_row':row,**q})
  else:harmkeys.add(hk)
 print('harmonized rows',row,flush=True)
for name,frame in [('harmonized_locus_ledger.tsv.gz',pd.DataFrame(ledger)),('harmonized_transform_failures.tsv.gz',pd.DataFrame(fails,columns=['harmonized_row','source_row',*headers,'status'])),('harmonized_duplicate_extra_rows.tsv.gz',pd.DataFrame(dupkeys)),('native_absent_or_ambiguous_in_harmonized.tsv.gz',d.loc[seen==0].assign(source_row=np.where(seen==0)[0]+1))]:frame.to_csv(O/name,sep='\t',index=False,compression={'method':'gzip','mtime':0})
a={'track':'A','rows':row,'harmonized_build':'GRCh38','source_md5_pass':True,'sha256':rec['sha256'],'gzip_full_decompression_crc':'PASS','columns':headers,'missingness':dict(missing),'harmonization_codes':dict(codes),'counts':dict(counts),'native_rows':len(d),'native_seen_once':int((seen==1).sum()),'native_seen_multiple':int((seen>1).sum()),'native_not_uniquely_seen':int((seen==0).sum()),'ambiguous_native_signature_count':len(ambiguous),'harmonized_duplicate_extra_rows':len(dupkeys),'native_locus_rows_harmonized':len(ledger),'transport_key':cols,'transport_note':'Exact tuple of all original allele/effect/test fields and source rsID; coordinate omitted ONLY for provenance crosswalk because Catalog updates positions to GRCh38. Not an analysis identity join. No allele recoding chosen by fit.','completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(O/'harmonized_audit.json').write_text(json.dumps(a,indent=2)+'\n');print(json.dumps(a,indent=2),flush=True)
