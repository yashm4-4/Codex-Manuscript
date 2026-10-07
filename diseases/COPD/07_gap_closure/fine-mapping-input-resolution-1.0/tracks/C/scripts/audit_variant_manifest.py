from pathlib import Path
import pandas as pd,numpy as np,gzip,json,time,collections,datetime
B=Path(__file__).resolve().parents[1];R=B/'results'
while not (R/'loci.tsv').exists() or not (B/'sources/variant_manifest.acquisition.json').exists():time.sleep(5)
a=json.loads((B/'sources/variant_manifest.acquisition.json').read_text());assert a['integrity_pass']
lo=pd.read_csv(R/'loci.tsv',sep='\t',dtype={'chrom':str});F=['chrom','pos','ref','alt','rsid','varid','pass_gnomad_genomes','n_passing_populations','high_quality','info','ac_EUR','af_EUR','an_EUR','gnomad_genomes_ac_EUR','gnomad_genomes_af_EUR','gnomad_genomes_an_EUR'];rows=0;out=gzip.open(R/'variant_manifest_locus_rows.tsv.gz','wt');first=True;counts=collections.Counter();missing=collections.Counter();hq=collections.Counter();amin=mathinf=float('inf');amax=-float('inf');badinfo=0;infole=0;infomin=float('inf');infomax=-float('inf');t=time.time()
for d in pd.read_csv(B/'data/full_variant_qc_metrics.txt.bgz',compression='gzip',sep='\t',usecols=F,dtype=str,keep_default_na=False,chunksize=250000,on_bad_lines='error'):
 p=pd.to_numeric(d['pos'],errors='coerce');sr=np.arange(rows+1,rows+len(d)+1);d['manifest_source_row']=sr
 for f in F:missing[f]+=int(d[f].isin(['NA','','NaN','nan']).sum())
 hq.update(d['high_quality'].value_counts().to_dict());info=pd.to_numeric(d['info'],errors='coerce');an=pd.to_numeric(d['an_EUR'],errors='coerce');badinfo+=int((~np.isfinite(info)).sum());infole+=int((info<=.8).sum());infomin=min(infomin,float(info.min()));infomax=max(infomax,float(info.max()));amin=min(amin,float(an.min()));amax=max(amax,float(an.max()))
 for l in lo.to_dict('records'):
  sel=(d['chrom']==l['chrom'])&(p>=l['start'])&(p<=l['end']);s=d.loc[sel].copy()
  if s.empty:continue
  s.insert(0,'locus_id',l['locus_id']);s.to_csv(out,sep='\t',index=False,header=first);first=False;counts[l['locus_id']]+=len(s)
 rows+=len(d)
 if rows%2000000==0:print('ROWS',rows,'SECONDS',round(time.time()-t,1),'EXPORTED',sum(counts.values()),flush=True)
out.close();res={'rows':rows,'gzip_complete_decode_crc_pass':True,'row_count_matches_provider':rows==28987534,'parsed_columns':F,'gene_columns_parsed_or_used':False,'missingness':dict(missing),'high_quality':dict(hq),'nonfinite_info':badinfo,'info_le_0_8':infole,'info_min':infomin,'info_max':infomax,'EUR_allele_number_min':amin,'EUR_allele_number_max':amax,'locus_rows':dict(counts),'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()};(R/'variant_manifest_audit.json').write_text(json.dumps(res,indent=2)+'\n');print('COMPLETE',rows,flush=True)
