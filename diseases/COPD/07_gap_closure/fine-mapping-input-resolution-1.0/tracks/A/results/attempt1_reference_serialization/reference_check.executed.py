from pathlib import Path
import pandas as pd,numpy as np,pysam,pysam.bcftools,json,collections,hashlib,sys,datetime
R=Path(__file__).resolve().parents[1];O=R/'results';stage=R.parents[1]
fa=Path(sys.argv[1]);f=pysam.FastaFile(str(fa));d=pd.read_csv(O/'locus_variants.tsv.gz',sep='\t',dtype={'chrom':str,'pos':int,'effect_allele':str,'other_allele':str},keep_default_na=False)
# Provider REF/ALT is independent of Kim effect-allele designation. Exact positional/allelic
# key matching precedes normalization; rsID never defines the join.
provider={};provider_sources=[]
for pth in sorted((stage/'ld/A').glob('*.index.tsv.gz')):
 ix=pd.read_csv(pth,sep='\t',dtype={'chrom':str,'ref':str,'alt':str},keep_default_na=False)
 provider_sources.append({'path':str(pth.relative_to(stage)),'sha256':hashlib.sha256(pth.read_bytes()).hexdigest(),'rows':len(ix)})
 for row in ix.itertuples(index=False):
  key=(row.chrom,int(row.pos),min(row.ref,row.alt),max(row.ref,row.alt));val=(row.ref,row.alt,int(row.idx),row.af)
  if val not in provider.setdefault(key,[]):provider[key].append(val)
unordered=pd.DataFrame({'chrom':d.chrom,'pos':d.pos,'lo':[min(x,y) for x,y in zip(d.effect_allele,d.other_allele)],'hi':[max(x,y) for x,y in zip(d.effect_allele,d.other_allele)]});d['native_unordered_duplicate_identity']=unordered.duplicated(keep=False)
comp=str.maketrans('ACGT','TGCA');rows=[];norm_records=[]
for r in d.itertuples(index=False):
 ea=r.effect_allele;oa=r.other_allele;c='chr'+r.chrom if 'chr'+r.chrom in f.references else r.chrom;p=r.pos;status='UNRESOLVED';ref='UNRESOLVED';alt='UNRESOLVED';ae=ea;ao=oa;strand='UNRESOLVED';idx='';af='';match=provider.get((r.chrom,p,min(ea,oa),max(ea,oa)),[])
 if ea==oa:status='INVALID_IDENTICAL_ALLELES'
 elif len(ea)==len(oa)==1 and ea in 'ACGT' and oa in 'ACGT':
  base=f.fetch(c,p-1,p).upper();direct=base in [ea,oa];rev=base in [ea.translate(comp),oa.translate(comp)]
  if direct and rev:status='UNRESOLVED_PALINDROMIC_STRAND'
  elif direct or rev:
   strand='FORWARD' if direct else 'REVERSE_COMPLEMENT_LABELS'
   if rev:ae=ea.translate(comp);ao=oa.translate(comp)
   ref=base;alt=ao if ae==base else ae;status='REFERENCE_VERIFIED_SNP'
   match=provider.get((r.chrom,p,min(ae,ao),max(ae,ao)),[])
   if len(match)==1 and match[0][:2]==(ref,alt):idx=match[0][2];af=match[0][3]
   elif len(match)>0:status='PROVIDER_REFERENCE_OR_CONTRAST_CONFLICT'
  else:status='REFERENCE_MISMATCH'
 elif len(match)==1:
  ref,alt,idx,af=match[0]
  if not all(x in 'ACGT' for x in ref+alt):status='UNRESOLVED_NON_SEQUENCE_ALLELE'
  elif f.fetch(c,p-1,p-1+len(ref)).upper()!=ref:status='PROVIDER_REFERENCE_MISMATCH'
  else:
   status='PROVIDER_AND_REFERENCE_VERIFIED_INDEL_PENDING_NORMALIZATION';strand='FORWARD';norm_records.append((c,p,r.source_row,ref,alt))
 elif len(match)>1:status='AMBIGUOUS_PROVIDER_CONTRAST'
 else:status='INDEL_OR_COMPLEX_REQUIRES_INDEPENDENT_REF_ALT'
 if bool(r.native_duplicate_identity or r.native_unordered_duplicate_identity):status+=';DUPLICATE_IDENTITY'
 rows.append((ref,alt,ae,ao,status,strand,idx,af,p))
x=pd.DataFrame(rows,columns=['ref','alt','analysis_effect_allele','analysis_other_allele','harmonization_status','strand_action','provider_ld_idx','provider_af','analysis_pos'])
for k in x:d[k]=x[k]
vcf=O/'reference_verified_indels.pre_normalization.vcf';normalized=O/'reference_verified_indels.normalized.vcf'
with vcf.open('w') as out:
 out.write('##fileformat=VCFv4.2\n')
 for c,n in zip(f.references,f.lengths):
  if c.removeprefix('chr') in set(map(str,range(1,23))):out.write(f'##contig=<ID={c},length={n}>\n')
 out.write('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
 for c,p,sr,ref,alt in norm_records:out.write(f'{c}\t{p}\t{sr}\t{ref}\t{alt}\t.\tPASS\t.\n')
if norm_records:
 result=pysam.bcftools.norm('-f',str(fa),'-c','e','-Ov',str(vcf));normalized.write_text(result)
 lookup={int(r.source_row):i for i,r in enumerate(d.itertuples(index=False))}
 for line in result.splitlines():
  if line.startswith('#'):continue
  c,p,sr,ref,alt,*_=line.split('\t');i=lookup[int(sr)];effect_is_ref=d.at[i,'analysis_effect_allele']==d.at[i,'ref'];d.at[i,'ref']=ref;d.at[i,'alt']=alt;d.at[i,'analysis_effect_allele']=ref if effect_is_ref else alt;d.at[i,'analysis_other_allele']=alt if effect_is_ref else ref;d.at[i,'analysis_pos']=int(p);d.at[i,'harmonization_status']=d.at[i,'harmonization_status'].replace('PROVIDER_AND_REFERENCE_VERIFIED_INDEL_PENDING_NORMALIZATION','PROVIDER_REFERENCE_NORMALIZED_INDEL')
else:normalized.write_text(vcf.read_text())
d['analysis_identity']=d.chrom+':'+d.analysis_pos.astype(str)+':'+d.ref+':'+d.alt;d['source_variant_identity']=d.variant_identity;d['analysis_beta']=d.beta;d['analysis_z']=d.z;d['normalization_shift_bp']=d.analysis_pos-d.pos
resolved=d.harmonization_status.isin(['REFERENCE_VERIFIED_SNP','PROVIDER_REFERENCE_NORMALIZED_INDEL']);dup=resolved&d.duplicated(['chrom','analysis_pos','ref','alt'],keep=False);d.loc[dup,'harmonization_status']+=';NORMALIZED_DUPLICATE_IDENTITY';d['normalized_duplicate_identity']=dup
# A distinct result keeps the original source audit intact.
d.to_csv(O/'locus_variants_reference_checked.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
ledger=d[['track','locus_id','source_row','chrom','pos','analysis_pos','rsid','effect_allele','other_allele','ref','alt','analysis_effect_allele','analysis_other_allele','source_variant_identity','analysis_identity','strand_action','harmonization_status','variant_type','provider_ld_idx','provider_af','normalization_shift_bp']];ledger.to_csv(O/'variant_harmonization.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
risk=d[['track','locus_id','analysis_identity','source_variant_identity','source_row','effect_allele','other_allele','analysis_effect_allele','analysis_other_allele','beta','odds_ratio','se','z','p','harmonization_status','variant_type']].copy();identitypass=d.harmonization_status.isin(['REFERENCE_VERIFIED_SNP','PROVIDER_REFERENCE_NORMALIZED_INDEL']);ready=identitypass&d.wald_p_rounding_pass&~d.rounding_direction_unresolved
risk['disease_increasing_allele']=np.where(ready,np.where(d.beta>0,d.analysis_effect_allele,d.analysis_other_allele),'UNRESOLVED');risk['direction_status']=np.where(ready,'VERIFIED_GWAS_SIGNED_CONTRAST_NOT_CAUSALITY','UNRESOLVED');risk.to_csv(O/'risk_direction_reference_checked.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
excluded=d.loc[~ready].copy();excluded['exclusion_scope']='FUTURE_ANALYSIS_IDENTITY_OR_DIRECTION_ONLY; source record retained';excluded.to_csv(O/'unresolved_variant_ledger.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
l=[]
for loc,g in d.groupby('locus_id',sort=False):
 sig=pd.to_numeric(g.p)<5e-8;ident=g.harmonization_status.isin(['REFERENCE_VERIFIED_SNP','PROVIDER_REFERENCE_NORMALIZED_INDEL'])
 l.append({'track':'A','locus_id':loc,'total_variants':len(g),'significant_variants':int(sig.sum()),'statistic_pass_variants':int(g.wald_p_rounding_pass.sum()),'reference_verified_variants':int(ident.sum()),'significant_unresolved_identity':int((sig&~ident).sum()),'statistic_fail_variants':int((~g.wald_p_rounding_pass).sum()),'provider_index_exact_matches':int((g.provider_ld_idx.astype(str)!='').sum()),'all_summary_gates':'NOT_CLEARED_PENDING_FULL_IDENTITY_AND_SAMPLE_CONTRACT','ld_gates':'ROOT_LD_AUDIT_PENDING','execution_state':'NOT CLEARED'})
pd.DataFrame(l).to_csv(O/'locus_summary_readiness.tsv',sep='\t',index=False)
a={'reference_fasta':str(fa),'reference_index':str(fa)+'.fai','provider_index_sources':provider_sources,'counts':dict(collections.Counter(d.harmonization_status)),'risk_direction_assigned':int(ready.sum()),'risk_direction_unresolved':int((~ready).sum()),'indel_normalization_software':{'pysam':pysam.__version__,'bcftools':pysam.__samtools_version__,'command':'bcftools norm -f verified_reference -c e -Ov'},'indels_normalized':len(norm_records),'positions_shifted_by_normalization':int((d.normalization_shift_bp!=0).sum()),'indels_not_guessed':True,'source_alleles_preserved':True,'unordered_duplicate_locus_rows':int(d.native_unordered_duplicate_identity.sum()),'normalized_duplicate_locus_rows':int(dup.sum()),'beta_unchanged_under_pure_complement':True,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(O/'reference_audit.json').write_text(json.dumps(a,indent=2)+'\n');print(json.dumps(a,indent=2))
