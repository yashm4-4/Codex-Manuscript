"""Reference validation and source-signed direction only; no regulatory annotations or posterior calculations."""
from pathlib import Path
import pandas as pd,numpy as np,json,datetime,hashlib,sys,pysam
from collections import Counter
r=Path(__file__).resolve().parents[1];o=r/'results';ref=Path(sys.argv[1]).resolve();fa=pysam.FastaFile(str(ref));d=pd.read_csv(o/'locus_variants.tsv.gz',sep='\t',dtype=str,keep_default_na=False);ledger=[]
comp=str.maketrans('ACGT','TGCA')
for row in d.itertuples(index=False):
 c='chr'+row.chrom;p=int(row.pos);a1=row.ref;a2=row.alt
 obs=fa.fetch(c,p-1,p-1+len(a1)).upper();isref=obs==a1;pal=len(a1)==len(a2)==1 and a1.translate(comp)==a2
 np1=p;nr=a1;na=a2;changed=False
 sequence_alleles=all(b in 'ACGT' for b in a1+a2)
 if isref and sequence_alleles:
  while nr[-1]==na[-1]:
   if len(nr)==1 or len(na)==1:
    if np1<=1:break
    prev=fa.fetch(c,np1-2,np1-1).upper();nr=prev+nr;na=prev+na;np1-=1
   nr=nr[:-1];na=na[:-1]
  while len(nr)>1 and len(na)>1 and nr[0]==na[0]:nr=nr[1:];na=na[1:];np1+=1
  changed=(np1,nr,na)!=(p,a1,a2)
  assert fa.fetch(c,np1-1,np1-1+len(nr)).upper()==nr
 kind='SYMBOLIC_OR_NONSEQUENCE' if not sequence_alleles else ('SNP' if len(a1)==len(a2)==1 else ('MNV' if len(a1)==len(a2) else 'INDEL'))
 status=('PASS_REF_AND_NORMALIZED_IDENTITY' if not changed else 'PASS_REF_NORMALIZED_IDENTITY_CHANGED') if isref else 'UNRESOLVED_SOURCE_REF_MISMATCH'
 if not sequence_alleles:status='UNRESOLVED_SYMBOLIC_OR_NONSEQUENCE_ALLELE'
 ledger.append({'track':'B','locus_id':row.locus_id,'source_row':row.source_row,'source_variant_identity':row.variant_identity,'chrom':row.chrom,'pos':p,'ref':a1,'alt':a2,'reference_sequence':obs,'reference_match':isref,'variant_type':kind,'palindromic_SNP':pal,'allele_swap_applied':False,'strand_complement_applied':False,'normalized_pos':np1 if isref else '', 'normalized_ref':nr if isref else '', 'normalized_alt':na if isref else '', 'normalized_identity':f'GRCh37:{row.chrom}:{np1}:{nr}:{na}' if isref and sequence_alleles else '', 'normalization_changed':changed,'harmonization_status':status})
h=pd.DataFrame(ledger);duplicate=h.normalized_identity.ne('') & h.normalized_identity.duplicated(keep=False);h['normalized_duplicate']=duplicate
conflict_ids=set(); duplicate_ids=set(h.loc[duplicate,'normalized_identity'])
for identity in sorted(duplicate_ids):
 inds=h.index[h.normalized_identity==identity];sub=d.loc[inds,['beta','se','p','effect_allele','other_allele']]
 if len(sub.drop_duplicates())>1:conflict_ids.add(identity)
h['duplicate_conflict']=h.normalized_identity.isin(conflict_ids);h.loc[h.duplicate_conflict,'harmonization_status']='UNRESOLVED_NORMALIZED_DUPLICATE_CONFLICT'
d['harmonization_status']=h.harmonization_status;d['normalized_identity']=h.normalized_identity;d['normalized_pos']=h.normalized_pos;d['normalized_ref']=h.normalized_ref;d['normalized_alt']=h.normalized_alt
h.to_csv(o/'variant_harmonization.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0});d.to_csv(o/'locus_variants_harmonized.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
beta=pd.to_numeric(d.beta,errors='coerce');signed=(d.qc_status=='PASS') & h.reference_match & ~h.duplicate_conflict & np.isfinite(beta)
risk=d[['track','locus_id','variant_identity','normalized_identity','effect_allele','other_allele','beta','se','z','p','source_row','harmonization_status']].copy();risk['disease_increasing_allele']=np.where(signed&beta.gt(0),d.effect_allele,np.where(signed&beta.lt(0),d.other_allele,''));risk['direction_status']=np.where(~signed,'UNRESOLVED',np.where(beta.eq(0),'ZERO_EFFECT_NO_DIRECTION','VERIFIED_SIGNED_GWAS_CONTRAST_ESTIMATED_DIRECTION'))
risk['interpretation']='Estimated signed COPD case/control contrast; not causal or proof of statistically nonzero risk';risk.to_csv(o/'risk_direction.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
excluded=d.loc[~signed].copy();excluded['reason']='Source-QC or reference-identity unresolved; retained in full locus table';excluded.to_csv(o/'excluded_variant_ledger.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
l=pd.read_csv(o/'loci.tsv',sep='\t');summ=[]
for row in l.itertuples(index=False):
 ix=d.locus_id==row.locus_id;ss=ix&signed;leadok=bool(signed[d.source_row==str(row.lead_source_row)].all())
 summ.append({'track':'B','locus_id':row.locus_id,'chrom':row.chrom,'start':row.start,'end':row.end,'mhc_deferred':row.mhc_deferred,'source_locus_rows':int(ix.sum()),'reference_matched_rows':int((ix&h.reference_match).sum()),'source_qc_pass_rows':int((ix&d.qc_status.eq('PASS')).sum()),'signed_direction_rows':int(ss.sum()),'unresolved_rows':int((ix&~signed).sum()),'lead_source_qc_ref_pass':leadok,'source_schema_signed_identity_gate':'PASS' if bool((signed[ix]).all()) else 'UNRESOLVED_ROWS_EXPLICIT','summary_likelihood_gate':'UNRESOLVED_SAIGE_SPA_CALIBRATED_Z_GAUSSIAN_RSS_SUITABILITY','all_summary_statistic_gates':'NOT_CLEARED','ld_gate':'NOT_RUN_NO_DEFENSIBLE_DENSE_SIGNED_JAPANESE_LD','summary_ld_gate':'NOT_RUN_NO_LD','readiness':'NOT CLEARED'})
pd.DataFrame(summ).to_csv(o/'locus_readiness.tsv',sep='\t',index=False)
missing=d[['track','locus_id','variant_identity','source_row']].copy();missing['ld_status']='UNASSESSED_LD_SOURCE_UNAVAILABLE';missing.to_csv(o/'missing_ld_variant_ledger.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
receipt={'track':'B','reference_path':str(ref),'reference_count':len(fa.references),'locus_rows':len(d),'reference_match_rows':int(h.reference_match.sum()),'reference_mismatch_rows':int((~h.reference_match).sum()),'normalization_changed_rows':int(h.normalization_changed.sum()),'palindromic_SNP_rows':int(h.palindromic_SNP.sum()),'normalized_duplicate_rows':int(duplicate.sum()),'normalized_conflict_identities':len(conflict_ids),'variant_types':dict(Counter(h.variant_type)),'risk_direction_status':dict(Counter(risk.direction_status)),'all_summary_gate_pass_loci':0,'loci_with_verified_signed_rows':sum(x['signed_direction_rows']>0 and not x['mhc_deferred'] for x in summ),'source_schema_signed_identity_pass_nonMHC_loci':sum(x['source_schema_signed_identity_gate']=='PASS' and not x['mhc_deferred'] for x in summ),'ld_gate_pass_loci':0,'execution_cleared_loci':0,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(o/'harmonization_audit.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2));fa.close()
