from pathlib import Path
import pandas as pd,numpy as np,json,time,collections,sys,pysam,datetime,hashlib
B=Path(__file__).resolve().parents[1];R=B/'results';fasta=Path(sys.argv[1])
for required in [R/'locus_export.json',R/'variant_manifest_audit.json',fasta,fasta.with_name(fasta.name+'.fai')]:
 while not required.exists():time.sleep(5)
d=pd.read_csv(R/'locus_variants.source.tsv.gz',sep='\t',dtype={'chrom':str},keep_default_na=False);m=pd.read_csv(R/'variant_manifest_locus_rows.tsv.gz',sep='\t',dtype=str,keep_default_na=False);m['pos']=m['pos'].astype(int);d['pos']=d['pos'].astype(int)
keys=['locus_id','chrom','pos','ref','alt'];assert not m.duplicated(keys).any(), 'manifest duplicate keys unresolved'
d=d.merge(m.rename(columns={c:(c if c.startswith('manifest_') else 'manifest_'+c) for c in m.columns if c not in keys}),on=keys,how='left',validate='many_to_one')
d['rsid']=d['manifest_rsid'].fillna('');fa=pysam.FastaFile(str(fasta));ctg=set(fa.references)
def refseq(c,p,n):return fa.fetch(c if c in ctg else 'chr'+c,p-1,p-1+n).upper()
def normalize(c,p,r,a):
 op=p
 # Left normalize indels without allele exchange or complement; retain source identity separately.
 if len(r)!=len(a):
  while r[-1]==a[-1]:
   if len(r)==1 or len(a)==1:
    if p==1:break
    prev=refseq(c,p-1,1);r=prev+r;a=prev+a;p-=1
   r=r[:-1];a=a[:-1]
 while len(r)>1 and len(a)>1 and r[0]==a[0]:r=r[1:];a=a[1:];p+=1
 return p,r,a
norm=[];refpass=[];vtypes=[];pal=[]
for i,r in enumerate(d.itertuples(index=False)):
 c=str(r.chrom);p=int(r.pos);ref=str(r.ref);alt=str(r.alt);ok=refseq(c,p,len(ref))==ref;refpass.append(ok);pal.append(len(ref)==len(alt)==1 and ref+alt in ['AT','TA','CG','GC']);vtypes.append('SNP' if len(ref)==len(alt)==1 else 'INDEL' if len(ref)!=len(alt) else 'MNV')
 np_,nr,na=normalize(c,p,ref,alt) if ok else (p,ref,alt);norm.append((np_,nr,na))
 if i and i%100000==0:print('HARMONIZED',i,flush=True)
d['normalized_pos']=[x[0] for x in norm];d['normalized_ref']=[x[1] for x in norm];d['normalized_alt']=[x[2] for x in norm];d['normalized_variant_id']='GRCh37:'+d['chrom']+':'+d['normalized_pos'].astype(str)+':'+d['normalized_ref']+':'+d['normalized_alt'];d['variant_type']=vtypes;d['palindromic_snp']=pal;d['reference_match']=refpass;d['allele_action']='NONE_FORWARD_ALT_EFFECT';d['normalization_changed']=(d['pos']!=d['normalized_pos'])|(d['ref']!=d['normalized_ref'])|(d['alt']!=d['normalized_alt']);d['source_identity_duplicate']=d.duplicated(keys,keep=False);d['normalized_identity_duplicate']=d.duplicated(['locus_id','normalized_variant_id'],keep=False)
d['harmonization_status']=np.select([~d['reference_match'],d['source_identity_duplicate'],d['normalized_identity_duplicate'],d['normalization_changed']],['UNRESOLVED_REF_MISMATCH','UNRESOLVED_DUPLICATE_SOURCE_IDENTITY','UNRESOLVED_NORMALIZATION_COLLISION','VERIFIED_FORWARD_ALT_CONTRAST_NORMALIZED_ID_RECORDED'],default='VERIFIED_FORWARD_ALT_CONTRAST_REFERENCE_MATCH')
# INFO and provider high_quality are retained, not used to redefine any already frozen locus boundary.
info=pd.to_numeric(d['manifest_info'],errors='coerce');d['variant_qc_annotation']=np.select([d['manifest_source_row'].isna(),~np.isfinite(info),info<=.8,d['manifest_high_quality']!='true'],['MISSING_VARIANT_MANIFEST_KEY','MISSING_INFO','INFO_NOT_GT_0_8','PROVIDER_HIGH_QUALITY_FALSE_OR_MISSING'],default='PROVIDER_VARIANT_QC_PASS')
beta=pd.to_numeric(d['beta'],errors='coerce');se=pd.to_numeric(d['se'],errors='coerce');bwidth=np.where(beta!=0,.0005*np.power(10,np.floor(np.log10(np.abs(beta)))),.0005)
verified=(d['qc_status']=='SOURCE_STATISTIC_QC_PASS')&d['reference_match']&~d['source_identity_duplicate']&~d['normalized_identity_duplicate']&np.isfinite(info)&(info>.8)
pos=verified&((beta-bwidth)>0);neg=verified&((beta+bwidth)<0);d['disease_increasing_allele']=np.select([pos,neg],[d['alt'],d['ref']],default='');d['direction_status']=np.select([pos|neg,verified],['VERIFIED_SIGNED_GWAS_CONTRAST','UNRESOLVED_ZERO_OR_ROUNDING_AMBIGUITY'],default='UNRESOLVED_GWAS_QC_OR_IDENTITY');d['direction_interpretation']='OBSERVED_DISEASE_LOG_ODDS_DIRECTION_NOT_CAUSAL_PROBABILITY'
d['native_signed_ld_join_eligible']=(d['qc_status']=='SOURCE_STATISTIC_QC_PASS')&d['harmonization_status'].str.startswith('VERIFIED')
d['native_signed_ld_join_scope']='DIAGNOSTIC_JOIN_ONLY_NOT_EXECUTION_CLEARANCE; PROVIDER_HIGH_QUALITY_IS_ANNOTATION'
d.to_csv(R/'locus_variants.tsv.gz',sep='\t',index=False,compression='gzip',na_rep='NA')
hcols=['track','locus_id','source_row','variant_id','chrom','pos','ref','alt','effect_allele','other_allele','normalized_pos','normalized_ref','normalized_alt','normalized_variant_id','variant_type','palindromic_snp','reference_match','allele_action','normalization_changed','source_identity_duplicate','normalized_identity_duplicate','harmonization_status','qc_status','variant_qc_annotation'];d[hcols].to_csv(R/'harmonization.tsv.gz',sep='\t',index=False,compression='gzip')
rcols=['track','locus_id','variant_id','effect_allele','other_allele','beta','se','z','p','neglog10_p','disease_increasing_allele','direction_status','harmonization_status','qc_status','direction_interpretation'];d[rcols].to_csv(R/'risk_direction.tsv.gz',sep='\t',index=False,compression='gzip',na_rep='NA')
ex=(d['qc_status']!='SOURCE_STATISTIC_QC_PASS')|~d['harmonization_status'].str.startswith('VERIFIED')|(d['variant_qc_annotation']!='PROVIDER_VARIANT_QC_PASS');d.loc[ex,hcols+['beta','se','neglog10_p']].to_csv(R/'excluded_or_flagged_variants.tsv.gz',sep='\t',index=False,compression='gzip')
res={'total_locus_rows':len(d),'reference_path':str(fasta),'reference_match_counts':d['reference_match'].value_counts().to_dict(),'harmonization_status':d['harmonization_status'].value_counts().to_dict(),'direction_status':d['direction_status'].value_counts().to_dict(),'qc_status':d['qc_status'].value_counts().to_dict(),'variant_qc_annotation':d['variant_qc_annotation'].value_counts().to_dict(),'normalization_changed_rows':int(d['normalization_changed'].sum()),'source_duplicate_rows':int(d['source_identity_duplicate'].sum()),'normalized_duplicate_rows':int(d['normalized_identity_duplicate'].sum()),'variant_type_counts':d['variant_type'].value_counts().to_dict(),'palindromic_snp_rows':int(d['palindromic_snp'].sum()),'per_locus':[],'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
for lid,g in d.groupby('locus_id',sort=False):
 sig=pd.to_numeric(g['neglog10_p'],errors='coerce')>-np.log10(5e-8);source_pass=(g['qc_status']=='SOURCE_STATISTIC_QC_PASS');res['per_locus'].append({'locus_id':lid,'all_rows':len(g),'source_statistic_qc_pass_rows':int(source_pass.sum()),'valid_signed_direction_rows':int((g['direction_status']=='VERIFIED_SIGNED_GWAS_CONTRAST').sum()),'source_qc_pass_significant_rows':int((sig&source_pass).sum()),'reference_mismatch_significant_qc_pass_rows':int((sig&source_pass&~g['reference_match']).sum()),'manifest_missing_rows':int(g['manifest_source_row'].isna().sum()),'variant_qc_pass_rows':int((g['variant_qc_annotation']=='PROVIDER_VARIANT_QC_PASS').sum()),'normalization_collision_rows':int(g['normalized_identity_duplicate'].sum())})
(R/'harmonization_audit.json').write_text(json.dumps(res,indent=2)+'\n');print('COMPLETE',len(d),flush=True)
