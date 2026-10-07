"""Provenance-bound Kim/Pan-UKB allele matching for diagnostics only; no LD or inference.

Source strings are retained. Each sign multiplier follows independently established
physical allele coding. Normalized indels are cross-checked against the archived
pre/post-bcftools VCF records and native provider index, never joined by rsID.
"""
from pathlib import Path
import hashlib,json,collections,datetime
import numpy as np
import pandas as pd

s=Path(__file__).resolve().parents[1]
track='A';root=s/'tracks/A';p=root/'results/locus_variants_reference_checked.tsv.gz'
gwas=pd.read_csv(p,sep='\t',dtype=str,keep_default_na=False)
index_summary=json.loads((s/'ld/A/index_export_summary.json').read_text())
reference_audit=json.loads((root/'results/reference_audit.json').read_text())
expected_hash={x['path']:x['sha256'] for x in reference_audit['provider_index_sources']}
input_hash=hashlib.sha256(p.read_bytes()).hexdigest()

def read_vcf(path):
    result={}
    for line in path.read_text().splitlines():
        if line.startswith('#'):continue
        chrom,pos,source_row,ref,alt,*_=line.split('\t')
        assert source_row not in result
        result[source_row]=(chrom.removeprefix('chr'),int(pos),ref,alt)
    return result
prepath=root/'results/reference_verified_indels.pre_normalization.vcf'
postpath=root/'results/reference_verified_indels.normalized.vcf'
pre=read_vcf(prepath);post=read_vcf(postpath)
assert pre.keys()==post.keys()
normalization_provenance={str(q.relative_to(s)):hashlib.sha256(q.read_bytes()).hexdigest() for q in [prepath,postpath]}
results=[]
allowed={'REFERENCE_VERIFIED_SNP','PROVIDER_REFERENCE_NORMALIZED_INDEL'}
for locus in index_summary:
    lid=locus['locus_id'];out=s/'ld/A'/lid;out.mkdir(exist_ok=True)
    g=gwas[gwas.locus_id==lid].copy().reset_index(drop=True)
    assert len(g) and (g.locus_status=='PROSPECTIVE').all()
    indexpath=s/'ld/A'/(lid+'.index.tsv.gz')
    digest=hashlib.sha256(indexpath.read_bytes()).hexdigest()
    assert digest==expected_hash[str(indexpath.relative_to(s))], 'Provider index changed after reference normalization audit'
    index=pd.read_csv(indexpath,sep='\t',dtype=str,keep_default_na=False)
    index_numbers=pd.to_numeric(index.idx,errors='raise').to_numpy(dtype=np.int64)
    assert len(index)==locus['index_variants'] and not index.idx.duplicated().any()
    assert np.all(np.diff(index_numbers)>0)
    assert int(index_numbers[0])==locus['min_idx'] and int(index_numbers[-1])==locus['max_idx']
    assert (index.chrom==str(locus['chrom'])).all()
    assert pd.to_numeric(index.pos).between(locus['start'],locus['end']).all()
    provider={int(r.idx):r for r in index.itertuples(index=False)}
    rows=[]
    for r in g.itertuples(index=False):
        reasons=[]
        if r.qc_status!='SIGNED_STATISTIC_PASS' or r.wald_p_rounding_pass!='True':reasons.append('FROZEN_SIGNED_STATISTIC_GATE_FAIL')
        if r.harmonization_status not in allowed:reasons.append('IDENTITY_UNRESOLVED:'+r.harmonization_status)
        if any(getattr(r,f)!='False' for f in ['native_duplicate_identity','native_unordered_duplicate_identity','normalized_duplicate_identity']):reasons.append('DUPLICATE_OR_CONFLICTING_IDENTITY')
        eligible=not reasons
        q=None;orientation=0;mapping='NOT_PROVEN';source_idx=r.provider_ld_idx
        if source_idx:
            if not source_idx.isdigit():reasons.append('INVALID_PROVIDER_INDEX_ENCODING')
            else:q=provider.get(int(source_idx))
            if q is None:reasons.append('RECORDED_PROVIDER_INDEX_NOT_IN_FROZEN_LOCUS')
        ld_present=q is not None
        if q is not None and eligible:
            analysis=(r.chrom,int(r.analysis_pos),r.ref,r.alt)
            indexed=(q.chrom,int(q.pos),q.ref,q.alt)
            if r.harmonization_status=='REFERENCE_VERIFIED_SNP':
                if analysis!=indexed:reasons.append('PROVIDER_ANALYSIS_SNP_IDENTITY_CONFLICT')
                elif (r.analysis_effect_allele,r.analysis_other_allele)==(q.alt,q.ref):orientation=1;mapping='VERIFIED_ANALYSIS_EFFECT_IS_LD_ALT'
                elif (r.analysis_effect_allele,r.analysis_other_allele)==(q.ref,q.alt):orientation=-1;mapping='VERIFIED_ANALYSIS_EFFECT_IS_LD_REF'
                else:reasons.append('PROVIDER_SNP_ALLELE_CONTRAST_CONFLICT')
            else:
                # Index coding uses the pre-normalization native REF/ALT. The
                # archived normalized record establishes the same physical variant.
                if pre.get(r.source_row)!=indexed or post.get(r.source_row)!=analysis:
                    reasons.append('INDEL_NORMALIZATION_PROVENANCE_CONFLICT')
                elif r.strand_action!='FORWARD':reasons.append('UNSUPPORTED_INDEL_STRAND_ACTION')
                elif (r.effect_allele,r.other_allele)==(q.alt,q.ref):orientation=1;mapping='VERIFIED_NATIVE_EFFECT_IS_LD_ALT_AFTER_NORMALIZATION'
                elif (r.effect_allele,r.other_allele)==(q.ref,q.alt):orientation=-1;mapping='VERIFIED_NATIVE_EFFECT_IS_LD_REF_AFTER_NORMALIZATION'
                else:reasons.append('PROVIDER_NATIVE_INDEL_ALLELE_CONTRAST_CONFLICT')
                if orientation:
                    expect=(r.alt,r.ref) if orientation==1 else (r.ref,r.alt)
                    if (r.analysis_effect_allele,r.analysis_other_allele)!=expect:
                        reasons.append('NORMALIZED_EFFECT_ORIENTATION_CONFLICT');orientation=0
        if eligible and not ld_present:reasons.append('EXACT_SIGNED_LD_IDENTITY_ABSENT')
        retained=eligible and ld_present and not reasons and orientation in [-1,1]
        rows.append({'summary_row_eligible':eligible,'ld_present':ld_present,
          'idx':'' if q is None else int(q.idx),'ld_chrom':'' if q is None else q.chrom,
          'ld_pos':'' if q is None else q.pos,'ld_ref':'' if q is None else q.ref,
          'ld_alt':'' if q is None else q.alt,'ld_af':'' if q is None else q.af,
          'ld_rsid':'' if q is None else q.rsid,'retained_for_diagnostics':retained,
          'exclusion_reason':';'.join(reasons),'orientation_multiplier':orientation if retained else '',
          'allele_mapping_status':mapping,'significant':float(r.p)<5e-8})
    m=pd.concat([g,pd.DataFrame(rows)],axis=1)
    # One physical LD row must not occur twice, including any residual collision
    # across native allele representations. Such a collision is a failure, never
    # resolved by selecting the more significant source row.
    collisions=m.retained_for_diagnostics&m.idx.duplicated(keep=False)
    m.loc[collisions,'retained_for_diagnostics']=False
    m.loc[collisions,'exclusion_reason']='DUPLICATE_SIGNED_LD_INDEX_AFTER_JOIN'
    m.loc[collisions,'orientation_multiplier']=''
    ledger_columns=['track','locus_id','source_row','chrom','pos','analysis_pos','ref','alt',
      'source_variant_identity','analysis_identity','effect_allele','other_allele',
      'analysis_effect_allele','analysis_other_allele','qc_status','harmonization_status',
      'p','se','beta','z','summary_row_eligible','ld_present','idx','ld_pos','ld_ref','ld_alt',
      'retained_for_diagnostics','exclusion_reason','significant','orientation_multiplier','allele_mapping_status']
    m[ledger_columns].to_csv(out/'variant_overlap_ledger.tsv.gz',sep='\t',index=False,compression={'method':'gzip','mtime':0})
    retained=m[m.retained_for_diagnostics].copy()
    retained['idx']=pd.to_numeric(retained.idx,errors='raise').astype('int64')
    retained=retained.sort_values('idx',kind='stable')
    assert not retained.idx.duplicated().any() and np.all(np.diff(retained.idx.to_numpy())>0)
    values=retained[['beta','se','z']].apply(pd.to_numeric,errors='raise').to_numpy(dtype=np.float64)
    assert np.isfinite(values).all() and (values[:,1]>0).all()
    signs=pd.to_numeric(retained.orientation_multiplier,errors='raise').to_numpy(dtype=np.int64)
    assert np.isin(signs,[-1,1]).all()
    retained['orientation_multiplier']=signs
    retained['z_ld']=values[:,2]*signs;retained['beta_ld']=values[:,0]*signs
    retained['n_gwas_nominal']=71591;retained['n_cases_nominal']=12446;retained['n_controls_nominal']=59145
    retained['n_per_variant_status']='NOT_REPORTED; nominal N is not verified observed/effective per-variant N'
    retained['ld_role']='EXTERNAL_FULL_EUR_UKB_REFERENCE_NOT_IN_SAMPLE_SMOKING_STRATUM'
    retained['input_status']='DIAGNOSTICS_ONLY_NOT_EXECUTION_CLEARED'
    # No source scalar is coerced or rewritten in the output table. Only the
    # explicitly named derived columns above introduce numerical values.
    ordered=out/'ordered_diagnostic_inputs.tsv.gz'
    retained.to_csv(ordered,sep='\t',index=False,compression={'method':'gzip','mtime':0})
    check=pd.read_csv(ordered,sep='\t',dtype=str,keep_default_na=False)
    source=g.set_index('source_row',drop=False)
    for col in g.columns:
        assert check[col].tolist()==source.loc[check.source_row,col].tolist(), 'Source scalar string changed: '+col
    absent=m[m.significant&~m.retained_for_diagnostics]
    leadrow=str(pd.read_csv(root/'results/loci.tsv',sep='\t',dtype=str,keep_default_na=False).set_index('locus_id').loc[lid,'lead_source_row'])
    lead=m[m.source_row==leadrow];assert len(lead)==1
    exclusions=collections.Counter(x for x in m.exclusion_reason if x)
    record=dict(locus,total_gwas_locus_rows=len(m),summary_row_eligible=int(m.summary_row_eligible.sum()),
      retained_diagnostic_variants=len(retained),missing_eligible_variants=int((m.summary_row_eligible&~m.ld_present).sum()),
      significant_rows=int(m.significant.sum()),missing_or_excluded_significant_rows=len(absent),
      lead_source_row=int(leadrow),lead_retained_for_diagnostics=bool(lead.retained_for_diagnostics.iloc[0]),
      lead_exclusion_reason=lead.exclusion_reason.iloc[0],orientation_plus=int((signs==1).sum()),orientation_minus=int((signs==-1).sum()),
      ld_index_join_collisions=int(collisions.sum()),exclusion_reasons=dict(exclusions),
      complete_variant_order_verified=True,source_scalar_strings_preserved=True,
      provider_index_sha256=digest,gwas_input_sha256=input_hash,ordered_inputs_sha256=hashlib.sha256(ordered.read_bytes()).hexdigest(),
      provider_index_order_verified=True,provider_index_exact_identity_duplicate_rows=int(index.duplicated(['chrom','pos','ref','alt'],keep=False).sum()),
      indel_normalization_provenance=normalization_provenance,
      summary_likelihood_contract='KIM_MARGINAL_LOGISTIC_LOG_OR_AND_SUPPLIED_SE; per-variant/effective N and smoking-subset/reference-LD compatibility unresolved',
      n_gwas_nominal=71591,n_cases_nominal=12446,n_controls_nominal=59145,n_per_variant='NOT_REPORTED',
      ld_source_role='Pan-UKB EUR full-European UKB external reference; not matched in-sample Kim ever-smoker LD',
      density_status='BLOCKED_MISSING_OR_EXCLUDED_SIGNIFICANT_VARIANTS' if len(absent) else 'NO_STRONG_SIGNAL_LOSS; full coverage and scientific compatibility still require review',
      input_status='DIAGNOSTICS_ONLY_NOT_EXECUTION_CLEARED',posterior_inference=False,
      completed_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
    (out/'overlap_summary.json').write_text(json.dumps(record,indent=2)+'\n');results.append(record)
    print(lid,'retained',len(retained),'missing/excluded significant',len(absent),'minus signs',int((signs==-1).sum()),flush=True)
(s/'ld/A/gwas_ld_overlap_summary.json').write_text(json.dumps(results,indent=2)+'\n')
