"""Assemble existing GWAS-only variant decisions verbatim; no scientific recomputation.

All source scalar strings are read as text. Risk assignments/statuses are copied
from the authoritative per-track risk tables, including every unresolved row.
"""
from pathlib import Path
import collections,datetime,gzip,hashlib,json
import pandas as pd

s=Path(__file__).resolve().parents[1];out=s/'tables';out.mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
specs={
 'A':dict(risk='risk_direction_reference_checked.tsv.gz',variants='locus_variants_reference_checked.tsv.gz',harm='variant_harmonization.tsv.gz',expected_rows=170904,expected_assigned=140886,assigned_status='VERIFIED_GWAS_SIGNED_CONTRAST_NOT_CAUSALITY',accession='GCST90016588',phenotype='Direct ever-smoker spirometric COPD',risk_label_basis='NORMALIZED_ANALYSIS_ALLELE_LABELS',source_id='source_variant_identity',analysis_id='analysis_identity'),
 'B':dict(risk='risk_direction.tsv.gz',variants='locus_variants_harmonized.tsv.gz',harm='variant_harmonization.tsv.gz',expected_rows=52013,expected_assigned=51989,assigned_status='VERIFIED_SIGNED_GWAS_CONTRAST_ESTIMATED_DIRECTION',accession='GCST90013709',phenotype='Direct combined-sex Japanese physician-diagnosed COPD',risk_label_basis='NATIVE_SOURCE_EFFECT_OTHER_ALLELE_LABELS',source_id='variant_identity',analysis_id='normalized_identity'),
 'C':dict(risk='risk_direction.tsv.gz',variants='locus_variants.tsv.gz',harm='harmonization.tsv.gz',expected_rows=259252,expected_assigned=175749,assigned_status='VERIFIED_SIGNED_GWAS_CONTRAST',accession='Pan-UKB_icd10_J44_both_sexes_EUR',phenotype='Secondary EHR/ICD J44 EUR COPD sensitivity',risk_label_basis='NATIVE_FORWARD_ALT_REF_ALLELE_LABELS',source_id='source_identity',analysis_id='normalized_variant_id')}
frames=[];inventory=[];records={};input_hashes={};native_fields={}
canonical=['track','accession_or_release_id','phenotype_role','locus_id','source_row','canonical_row_id','source_genome_build','analysis_genome_build','source_identity','analysis_identity','effect_allele','other_allele','analysis_effect_allele','analysis_other_allele','beta','odds_ratio','se','z','p','neglog10_p','disease_increasing_allele','direction_status','harmonization_status','qc_status','risk_allele_label_basis','native_beta_string','native_odds_ratio_string','native_se_string','native_p_string','native_neglog10_p_string','source_p_encoding','source_af','source_n','source_risk_path','source_harmonization_path','source_variant_path','table_scope']
for track,spec in specs.items():
    root=s/'tracks'/track/'results'
    paths={k:root/spec[k] for k in ['risk','variants','harm']}
    for p in paths.values():input_hashes[str(p.relative_to(s))]=sha(p)
    r=pd.read_csv(paths['risk'],sep='\t',dtype=str,keep_default_na=False)
    v=pd.read_csv(paths['variants'],sep='\t',dtype=str,keep_default_na=False)
    h=pd.read_csv(paths['harm'],sep='\t',dtype=str,keep_default_na=False)
    assert len(r)==len(v)==len(h)==spec['expected_rows']
    keys=['track','locus_id','source_row']
    assert not v.duplicated(keys).any() and not h.duplicated(keys).any()
    assert (v.track==track).all() and (r.track==track).all() and (h.track==track).all()
    # The C risk writer omits source_row. Its unique native variant identity,
    # locus, alleles, all shared scalar strings and row order are verified
    # against the final variant table before copying that existing source_row.
    if 'source_row' not in r:
        assert not v.duplicated(['track','locus_id','variant_id']).any()
        for c in ['track','locus_id','variant_id','effect_allele','other_allele']:
            assert r[c].tolist()==v[c].tolist(), (track,'risk row linkage',c)
        linkage='Exact unique native variant/locus keys, alleles, row order and all shared scalar fields; source_row copied from final variant table'
        r=r.copy();r['source_row']=v.source_row
    else:
        linkage='Exact track/locus/source_row keys; no rsID-only linkage'
        assert not r.duplicated(keys).any()
        v=v.set_index(keys,drop=False).loc[pd.MultiIndex.from_frame(r[keys])].reset_index(drop=True)
    h=h.set_index(keys,drop=False).loc[pd.MultiIndex.from_frame(r[keys])].reset_index(drop=True)
    for c in [x for x in r.columns if x in v.columns]:
        assert r[c].tolist()==v[c].tolist(), (track,'risk/variant scalar mismatch',c)
    assert r.harmonization_status.tolist()==h.harmonization_status.tolist()
    assigned=int((r.direction_status==spec['assigned_status']).sum())
    assert assigned==spec['expected_assigned']
    n=len(r);d=pd.DataFrame('',index=range(n),columns=canonical,dtype=str)
    d['track']=track;d['accession_or_release_id']=spec['accession'];d['phenotype_role']=spec['phenotype']
    d['locus_id']=r.locus_id;d['source_row']=r.source_row;d['canonical_row_id']=track+':'+r.locus_id+':source_row='+r.source_row
    d['source_genome_build']='GRCh37';d['analysis_genome_build']='GRCh37'
    d['source_identity']=v[spec['source_id']];d['analysis_identity']=v[spec['analysis_id']]
    for c in ['effect_allele','other_allele','beta','se','z','p','disease_increasing_allele','direction_status','harmonization_status']:
        d[c]=r[c]
    for c in ['analysis_effect_allele','analysis_other_allele','odds_ratio','neglog10_p']:
        if c in r:d[c]=r[c]
        elif c=='neglog10_p' and c in v:d[c]=v[c]
    d['qc_status']=v.qc_status;d['risk_allele_label_basis']=spec['risk_label_basis']
    d['source_af']=v.af;d['source_n']=v.n
    d['source_risk_path']=str(paths['risk'].relative_to(s));d['source_harmonization_path']=str(paths['harm'].relative_to(s));d['source_variant_path']=str(paths['variants'].relative_to(s));d['table_scope']='COPY_OF_GWAS_ONLY_ESTIMATED_DIRECTION; NO_EXECUTION_CLEARANCE'
    if track=='A':
        originals=['variant_id','p_value','chromosome','base_pair_location','effect_allele','other_allele','odds_ratio','standard_error']
        mapping={c:'native_A_'+c for c in originals}
        d['native_odds_ratio_string']=v.odds_ratio;d['native_se_string']=v.standard_error;d['native_p_string']=v.p_value;d['source_p_encoding']='ORDINARY_P'
    elif track=='B':
        originals=[c for c in v if c.startswith('original_')];mapping={c:'native_B_'+c.removeprefix('original_') for c in originals}
        d['native_beta_string']=v.original_BETA;d['native_se_string']=v.original_SE;d['native_p_string']=v['original_p.value'];d['source_p_encoding']='ORDINARY_SOURCE_SAIGE_P; p.value.NA also preserved separately'
    else:
        originals=[c for c in v if c.startswith('original_')];mapping={c:'native_C_'+c.removeprefix('original_') for c in originals}
        d['native_beta_string']=v.original_beta_EUR;d['native_se_string']=v.original_se_EUR;d['native_neglog10_p_string']=v.original_neglog10_pval_EUR;d['source_p_encoding']='SOURCE_NEGATIVE_LOG10_P; ordinary P copied from existing derivative'
    for c,name in mapping.items():d[name]=v[c]
    native_fields[track]=mapping
    # Copy-only invariants, evaluated for every row, include all unassigned rows.
    for c in ['disease_increasing_allele','direction_status','harmonization_status','effect_allele','other_allele','beta','se','z','p']:
        assert d[c].tolist()==r[c].tolist()
    inventory.append({'track':track,'harmonization_path':str(paths['harm'].relative_to(s)),'harmonization_rows':len(h),'harmonization_bytes':paths['harm'].stat().st_size,'harmonization_sha256':sha(paths['harm']),'authoritative_risk_path':str(paths['risk'].relative_to(s)),'risk_rows':len(r),'risk_bytes':paths['risk'].stat().st_size,'risk_sha256':sha(paths['risk']),'authoritative_variant_path':str(paths['variants'].relative_to(s)),'variant_rows':len(v),'variant_sha256':sha(paths['variants']),'assigned_directions':assigned,'unresolved_directions':n-assigned,'loci_including_deferred':r.locus_id.nunique(),'risk_allele_label_basis':spec['risk_label_basis']})
    records[track]={'rows':n,'assigned':assigned,'unresolved':n-assigned,'assigned_status':spec['assigned_status'],'direction_status_counts':dict(collections.Counter(r.direction_status)),'risk_allele_label_basis':spec['risk_label_basis'],'row_linkage':linkage,'all_source_risk_scalar_strings_copied_verbatim':True,'all_native_source_strings_copied_verbatim':True}
    frames.append(d)
    print(track,'rows',n,'assigned',assigned,'unresolved',n-assigned,flush=True)
all_native=sorted({c for d in frames for c in d if c not in canonical})
combined=pd.concat([d.reindex(columns=canonical+all_native,fill_value='') for d in frames],ignore_index=True)
assert len(combined)==482169 and not combined.canonical_row_id.duplicated().any()
riskpath=out/'gwas_only_risk_direction.tsv.gz'
combined.to_csv(riskpath,sep='\t',index=False,compression={'method':'gzip','mtime':0})
pd.DataFrame(inventory).to_csv(out/'variant_harmonization_inventory.tsv',sep='\t',index=False)
# Reading back as strings validates the assembled payload, not merely the
# in-memory object or its schema. No direction is recalculated.
written=pd.read_csv(riskpath,sep='\t',dtype=str,keep_default_na=False)
assert written.columns.tolist()==combined.columns.tolist() and written.equals(combined)
for path,expected in input_hashes.items():assert sha(s/path)==expected,'Input changed during assembly: '+path
metadata={
 'title':'README metadata for central GWAS-only direction and harmonization ledgers',
 'scope':'Verbatim assembly of final existing A/B/C results; no new allele, statistic, harmonization, selection, QC, risk-direction or execution-readiness decisions',
 'total_rows':len(combined),'assigned_rows':sum(x['assigned'] for x in records.values()),'unresolved_rows':sum(x['unresolved'] for x in records.values()),
 'all_prospective_locus_rows_preserved_including_deferred_MHC_and_unresolved':True,'tracks':records,
 'column_order':combined.columns.tolist(),'native_source_column_mapping':native_fields,
 'column_semantics':{
  'canonical_row_id':'Unique assembly key constructed from existing track/locus/source-row identifiers. Biological variant identities may duplicate or remain unresolved.',
  'source_identity':'Exact source identity string already present in final track variants. A uses chromosome:position:EA/OA; B/C include GRCh37 and REF/ALT. No identity is inferred or joined by rsID.',
  'analysis_identity':'Exact existing normalized identity string. A lacks an inline assembly prefix; the adjacent analysis_genome_build supplies GRCh37. Unresolved components/blank values remain unchanged. This column alone does not authorize analysis.',
  'effect_allele_and_other_allele':'Copied from authoritative risk tables. These are original source allele labels. A supplies separate normalized analysis allele labels. B/C normalized alleles are available in their harmonization tables; no new central relabeling is made.',
  'analysis_effect_allele_and_analysis_other_allele':'Populated only when explicitly present in the authoritative risk table (A); blank for B/C rather than inferring a new representation.',
  'disease_increasing_allele':'Copied exactly, never recalculated. A uses normalized analysis labels; B/C use native source labels. Estimated direction is not causality, statistical significance, or fine-mapping clearance.',
  'beta_odds_ratio_se_z_p':'Exact strings from authoritative risk outputs. A beta is previously computed log(OR). B/C use source-supplied effect/test outputs with their documented scale limitations. OR remains blank for B/C because it was not supplied; nothing is exponentiated.',
  'native_source_strings':'native_beta_string/native_odds_ratio_string/native_se_string/native_p_string/native_neglog10_p_string and track-prefixed native_* fields copy the source strings preserved by final track files. All original per-track source columns are retained; absent track-specific fields are blank.',
  'source_p_encoding':'A/B retain ordinary source P. C source stores negative-log10-P; its ordinary P is a pre-existing derivative copied without recalculation. BBJ p.value.NA is retained independently of p.value.',
  'source_af_and_source_n':'Copied from final variants with their original limitations: A AF/N absent; B original variant AF/N; C weighted AF approximation and phenotype-level nominal N, not measured per-variant N. This assembly performs no AF/N derivation.',
  'status_fields':'Original track-specific status strings retained. They have different scopes and must not be interpreted as a common inference-readiness gate.',
  'missing_values':'Original NA, UNRESOLVED and empty-string conventions are preserved. Blank canonical fields can mean not separately recorded. Track-prefixed fields are blank for other tracks.'},
 'statistic_cautions':{
  'A':'Frozen literal decimal-precision failures remain unresolved; complete sample/LD compatibility remains unclosed.',
  'B':'Approximate signed effect direction can be valid while source score/software lineage and SAIGE SPA/RSS likelihood remain unresolved.',
  'C':'ALT-signed beta and SPA-calibrated SE/P; algebraic beta/SE-P agreement is not independent validation of Gaussian/RSS covariance.'},
 'input_sha256':input_hashes,'output_risk_sha256':sha(riskpath),'output_inventory_sha256':sha(out/'variant_harmonization_inventory.tsv'),
 'validation':{'expected_track_counts_pass':True,'all_shared_risk_variant_strings_match':True,'harmonization_status_alignment_pass':True,'all_output_strings_round_trip_exact':True,'all_input_hashes_unchanged':True,'all_rows_retained':True,'new_scientific_decisions':False},
 'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
(out/'variant_ledgers_metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
print('central rows',len(combined),'assigned',metadata['assigned_rows'],'unresolved',metadata['unresolved_rows'],'all string-preservation checks passed',flush=True)
