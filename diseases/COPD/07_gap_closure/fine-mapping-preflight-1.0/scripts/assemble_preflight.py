#!/usr/bin/env python3
"""Join evidence audits into decision tables. Does not analyze GWAS associations or LD."""
import csv,datetime,hashlib,json,pathlib,re
STAGE=pathlib.Path(__file__).resolve().parents[1]
def read(path):return json.loads((STAGE/path).read_text())
def write_table(name,rows,fields=None):
    fields=fields or list(dict.fromkeys(k for r in rows for k in r))
    with (STAGE/'tables'/name).open('w') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for row in rows:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v for k,v in row.items()})
def main():
    catalog=list(csv.DictReader((STAGE/'tables/study_phenotype_summary_stat_eligibility.tsv').open(),delimiter='\t'))
    direct=read('audits/direct/direct_source_audit.json');alternatives=read('audits/alternatives/alternatives_audit.json');methods=read('audits/methods_ld/methods_ld_audit.json')
    direct_by={r['study_accession']:r for r in direct['studies']}
    external_by={r['accession']:r for r in alternatives['datasets'] if r.get('accession')}
    provider_by={r['dataset_id']:r for r in alternatives['datasets']}
    # Exact frozen phenotype/count and current accession headers match these provider strata.
    crosswalk={'GCST90691934':'PanUKB_J44_EUR','GCST90692407':'PanUKB_J44_CSA'}
    for accession,dataset in crosswalk.items():external_by[accession]=provider_by[dataset]
    ld_by={r['accession']:r for r in methods['dataset_ld_recommendations'] if r.get('accession')}
    kim={'GCST90016588','GCST90016589','GCST90016593','GCST90016594'}
    bbj={'GCST90013709','GCST90013746','GCST90013781'}
    matrix=[]
    for c in catalog:
        a=c['accession'];d=direct_by.get(a,{});external=external_by.get(a,{})
        access=c['public_access_status'];signed=c['signed_effect_nonmissing_in_sample']=='True' and c['effect_allele_fields_present']=='True'
        if a in bbj:
            assert external.get('actual_header'),a
            access='PUBLIC_ORIGINAL_BBJ_SIGNED_AUTOSOMAL_AND_SEX_SPECIFIC_SCHEMA_VERIFIED_BY_BOUNDED_ZIP_RANGES'
            signed=True
        elif a in ('GCST007692','GCST004147'):
            access='PUBLIC_SMALL_DBGAP_EXPORT_VERIFIED_BUT_INCOMPLETE_AND_UNSIGNED; complete signed file requires source request'
            signed=False
        elif a=='GCST90246122':access='CATALOG_SUSPENDED; provider encrypted statistics require application for decryption key';signed=False
        if c['frozen_phenotype_class']=='nonstandard_gene_burden_cnv':signed=False
        role=c['phenotype_eligibility']
        if a=='GCST90016588':role='CONDITIONAL_PRIMARY_DIRECT_COPD_EVER_SMOKERS; preferred public European starting stratum if LD/scale gates resolved'
        elif a=='GCST90016589':role='SEPARATE_DIRECT_COPD_NEVER_SMOKER_STRATUM; not pooled with ever smokers without a new approved model'
        elif a in ('GCST90016593','GCST90016594'):role='OVERLAPPING_DIRECT_COPD_SMOKING_STRATUM_SENSITIVITY; not independent replication'
        elif a=='GCST90013709':role='CONDITIONAL_JAPANESE_DIRECT_COPD_PRIMARY_TRACK; matched Japanese LD access unresolved'
        elif a in bbj:role='DIRECT_COPD_SEX_STRATUM_SENSITIVITY; not independent of combined BBJ cohort'
        if a in crosswalk:
            ld_status='PUBLIC_SAME_PROJECT_ANCESTRY_LD_ROUTE_CONDITIONAL_NOT_LOCUS_VALIDATED'
            ld_preferred='PanUKB '+external['ancestry']+' signed ldadj dosage matrix and exact variant index'
            ld_reason='Exact release/N/scaling, phenotype sample and association-test compatibility and locus QC pending; CSA is one provider ancestry stratum, not a pooled EUR/CSA meta-analysis.'
        elif a in ld_by:ld_status=ld_by[a]['status'];ld_preferred=ld_by[a]['preferred_ld'];ld_reason=ld_by[a]['reason']
        elif a=='GCST007692':ld_status='MATCHED_META_ANALYSIS_LD_NOT_VERIFIED';ld_preferred='Original cohort/meta-analysis-compatible signed LD';ld_reason='Predominantly European is not European-only; original full signed data and compatible ancestry contribution are unresolved.'
        elif external:ld_status='SOURCE_SPECIFIC_LD_OR_LIKELIHOOD_GATES_UNRESOLVED';ld_preferred=external.get('LD_candidate','unverified');ld_reason=';'.join(external.get('LD_blockers',[]))
        else:ld_status='NO_STUDY_SPECIFIC_SIGNED_LD_VALIDATED';ld_preferred=c['LD_candidate_and_caveat'];ld_reason='No complete matching matrix/sample/coding contract established.'
        fields=external.get('actual_header',[]) if a in bbj else json.loads(c['inspected_header'])
        ancestry_caveat=d.get('ancestry_meta_analysis_status',c['ancestry_specific_vs_meta'])
        if a in crosswalk:ancestry_caveat='Verified PanUKB ancestry-specific '+external['ancestry']+' stratum; source-specific analysis, not pooled multi-ancestry meta-analysis'
        elif a=='GCST90692991':ancestry_caveat='PanUKB J44 EUR plus CSA meta-analysis; current accession cases/controls identify both constituent strata; no single-EUR matrix substitution'
        elif a in bbj:ancestry_caveat='Original BBJ Japanese ancestry, '+external['dataset_id'].rsplit('_',1)[-1]+' sex stratum; no trans-ancestry pooling'
        matrix.append({'accession':a,'publication':c['publication'],'pubmed_id':c['pubmed_id'],'doi':c['doi'],
          'frozen_phenotype_class':c['frozen_phenotype_class'],'case_definition':c['case_definition'],'control_definition':c['control_definition'],
          'cohorts':c['cohorts'],'ancestry_and_analysis_specific_caveat':ancestry_caveat,
          'current_catalog_discovery_ancestry':c['discovery_ancestry_current'],'current_catalog_N_description':c['initial_sample_size_current'],
          'external_provider_cases':external.get('cases'),'external_provider_controls':external.get('controls'),
          'overall_public_access_status':access,'signed_variant_effect_schema_verified':signed,
          'native_build_verified':external.get('build') if a in bbj else c['native_genome_build'],
          'preferred_source_header':fields,'signed_effect_details':external.get('effect_field') if a in bbj else d.get('signed_effect_fields',c['signed_effect_fields']),
          'effect_allele_semantics':external.get('effect_allele') if a in bbj else ('Use the actual effect-allele and effect fields from the same source representation; harmonized names may be effect_allele or hm_effect_allele and must remain paired with that file transformed effect, never mixed with raw fields' if signed else 'UNRESOLVED_OR_NOT_VARIANT_SIGNED'),
          'recommended_role':role,'LD_status':ld_status,'preferred_LD':ld_preferred,'LD_blocker_or_caveat':ld_reason,
          'whole_file_integrity_density_QC_passed':False,'locus_QC_passed':False,'inference_executable_now':False,
          'eligible_locus_count':None,'eligible_locus_count_status':'NOT_DETERMINED; zero loci execution-cleared is not evidence of zero GWAS loci',
          'catalog_detail_table':'tables/study_phenotype_summary_stat_eligibility.tsv',
          'direct_review':f'audits/direct/direct_source_audit.json#{a}' if d else '',
          'provider_review':f'audits/alternatives/alternatives_audit.json#{external["dataset_id"]}' if external else '',
          'study_limitation_reference':'tables/study_specific_limitations.tsv;audits/methods_ld/METHOD_AND_ELIGIBILITY_DESIGN.md'})
    assert len(matrix)==104
    accessible_direct=[r for r in matrix if r['frozen_phenotype_class']=='direct_copd_susceptibility' and r['signed_variant_effect_schema_verified']]
    assert {r['accession'] for r in accessible_direct}==kim|bbj
    write_table('study_final_eligibility_matrix.tsv',matrix)
    write_table('provider_accession_crosswalk.tsv',[{'accession':a,'provider_dataset_id':d,'evidence':'Frozen J44 phenotype and current accession cases/controls match exact provider manifest stratum; source/header audit retained','source_review':f'audits/alternatives/alternatives_audit.json#{d}'} for a,d in crosswalk.items()])
    catalog_by={r['accession']:r for r in catalog}
    risks=[];limitations=[];ldrows=[];ancestryrows=[]
    for row in matrix:
        a=row['accession'];c=catalog_by[a];external=external_by.get(a,{})
        fieldnames=row['preferred_source_header']
        def matching(pattern):return ';'.join(f for f in fieldnames if re.search(pattern,f,re.I)) or 'not_present_in_inspected_header'
        signed=row['signed_variant_effect_schema_verified']
        allele=external.get('effect_allele') if a in bbj else c['effect_other_allele_fields']
        effect=external.get('effect_field') if a in bbj else c['signed_effect_fields']
        se=external.get('se_field') if a in bbj else c['SE_fields']
        p=external.get('p_field') if a in bbj else c['P_fields']
        evidence=row['provider_review'] if a in bbj else c['exact_evidence_record']
        risk_status='SIGNED_SCHEMA_AVAILABLE; source outcome/test-scale, exact allele harmonization and variant-level QC still required' if signed else 'UNVERIFIED_OR_NOT_VALID_SIGNED_VARIANT_SOURCE'
        if row['frozen_phenotype_class']=='nonstandard_gene_burden_cnv':risk_status='NOT_A_VARIANT_RISK_ALLELE; aggregate tests do not define individual-variant direction'
        risks.append({'accession':a,'frozen_phenotype_class':row['frozen_phenotype_class'],'evidence_scope':'INTEGRATED_CATALOG_AND_VERIFIED_PROVIDER',
          'source_evidence':evidence,'source_url':'https://humandbs.dbcls.jp/files/hum0014/hum0014.v17.COPD.v1.zip' if a in bbj else c['selected_header_source_url'],
          'effect_allele_fields_and_semantics':allele,'signed_effect_fields':effect,'SE_fields':se,'P_fields':p,
          'frequency_fields':matching(r'frequency|af\.|^af_|AF_'),'per_variant_N_fields':matching(r'^n$|^n_|^N\.'),
          'SE_values_preview_status':'provider representative records and schema verified; no whole-file SE QC' if a in bbj else c['SE_nonmissing_in_sample'],
          'risk_direction_feasibility':risk_status,'risk_alleles_assigned':False,'build':row['native_build_verified'],
          'source_representation_rule':row['effect_allele_semantics']})
        blockers=[]
        if not signed:blockers.append('Valid signed variant source not established; aggregate/unsigned/missing sources are not usable individual-variant directions')
        if a not in bbj and c['SE_nonmissing_in_sample']!='True':blockers.append('Nonmissing SE not demonstrated in usable preview; do not invent from incompatible/rounded P or CI')
        if a in bbj:blockers.append('Exact Japanese cohort/sex/dosage signed LD not verified; complete source QC pending')
        if a in kim:blockers.append('Smoking-stratum LD and Catalog log-OR/SE conversion lineage need validation; AF/N/quality absent from native header')
        if a in ('GCST007692','GCST004147'):blockers.append('Public dbGaP export is incomplete and explicitly absolute/unsigned; request original full signed statistics')
        if a=='GCST90246122':blockers.append('Provider full data encrypted; decryption-key application required')
        blockers += [row['LD_blocker_or_caveat'],'Complete-file integrity, exact variant density/overlap and all scientific locus gates NOT RUN']
        limitations.append({'accession':a,'evidence_scope':'INTEGRATED_CATALOG_AND_VERIFIED_PROVIDER','public_access_status':row['overall_public_access_status'],
          'phenotype_and_design_role':row['recommended_role'],'current_blockers':blockers,'catalog_missing_source_access_resolved_by_provider':a in bbj,
          'frozen_source_uncertainties':c['frozen_uncertainties'],'authoritative_final_matrix':'tables/study_final_eligibility_matrix.tsv'})
        ldrows.append({'accession':a,'analysis_ancestry':row['ancestry_and_analysis_specific_caveat'],'cohorts':row['cohorts'],'preferred_LD':row['preferred_LD'],
          'LD_status':row['LD_status'],'compatibility_caveat':row['LD_blocker_or_caveat'],'matrix_values_acquired':False,'locus_matrix_QC_executed':False})
        ancestryrows.append({'accession':a,'cohorts':row['cohorts'],'analysis_specific_ancestry':row['ancestry_and_analysis_specific_caveat'],
          'discovery_ancestry_current':row['current_catalog_discovery_ancestry'],'sample_description_current':row['current_catalog_N_description'],
          'provider_cases_if_verified':row['external_provider_cases'],'provider_controls_if_verified':row['external_provider_controls'],
          'smoking_design':c['smoking_design'],'sex_design':c['sex_design'],'overlap_design':'STUDY_OVERLAP_CONSIDERATIONS.md; original frozen phenotype cohort and pair ledgers hash-bound by initialization',
          'sample_size_caveat':'Do not sum overlapping strata, duplicate YAML entries or selected follow-up sample counts'})
    write_table('signed_effect_risk_direction_feasibility.tsv',risks)
    write_table('study_specific_limitations.tsv',limitations)
    write_table('study_LD_requirements.tsv',ldrows)
    write_table('ancestry_cohort_sample_size.tsv',ancestryrows)
    write_table('accessible_signed_direct_COPD_sources.tsv',accessible_direct)
    write_table('provider_dataset_ancestry_inventory.tsv',alternatives['datasets'])
    write_table('direct_COPD_source_review.tsv',direct['studies'])
    write_table('LD_resource_inventory.tsv',methods['resources'])
    write_table('LD_source_compatibility_matrix.tsv',methods['dataset_ld_recommendations'])
    write_table('proposed_fine_mapping_eligibility_gates.tsv',methods['gates'])
    write_table('numerical_design_rule_evidence.tsv',methods['numeric_rule_provenance'])
    write_table('proposed_locus_definition_rules.tsv',[{'rule':k,'proposal':v,'executed':False} for k,v in methods['locus_design'].items()])
    write_table('immediately_executable_fine_mapping_datasets.tsv',[],['dataset','status','reason'])
    blocked=[{'accession':r['accession'],'current_access_status':r['overall_public_access_status'],
              'required_resolution':'Obtain original complete signed source and/or resolve documented access restriction; lack of a verified endpoint is not proof of global unavailability',
              'prior_contact_sent':False} for r in matrix if not r['signed_variant_effect_schema_verified'] and not r['recommended_role'].startswith('UNSUITABLE')]
    write_table('datasets_requiring_access_or_signed_source_resolution.tsv',blocked)
    unsuitable=[{'accession':r['accession'],'reason':r['recommended_role'],'not_a_claim_of_no_biological_association':True}
        for r in matrix if r['recommended_role'].startswith(('UNSUITABLE','RESTRICTED','BLOCKED','SEPARATE_ESTIMAND','SEPARATE_ML'))]
    write_table('datasets_unsuitable_or_separate_estimand.tsv',unsuitable)
    readiness=[{'dataset':r['accession'],'permitted_next_stage_design':'Full source ingestion/harmonization design only after investigator approval','signed_source_schema_available':True,
        'LD_resolution':r['LD_status'],'not_execution_cleared':True} for r in accessible_direct]
    readiness.append({'dataset':'PanUKB_J44_EUR','permitted_next_stage_design':'Secondary EHR source/LD-contract validation after investigator approval','signed_source_schema_available':True,
        'LD_resolution':'PUBLIC_SAME_PROJECT_LD_METADATA; release/scaling/N/case-control/subset/variant QC unresolved','not_execution_cleared':True})
    write_table('priority_datasets_for_future_input_validation.tsv',readiness)
    recommendations={'stage':'fine-mapping-preflight-1.0','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
      'core_accessions_audited':104,'direct_core_accessions':27,'direct_accessions_with_verified_public_signed_schema':7,
      'direct_signed_accessions':sorted(kim|bbj),'provider_dataset_ancestry_records':len(alternatives['datasets']),
      'immediately_executable_fine_mapping_dataset_count':0,'loci_execution_cleared':0,'expected_eligible_locus_count':None,
      'locus_count_interpretation':'Unknown until full GWAS-only locus definition and matched-LD QC; neither zero cleared nor published locus totals estimates the eligible count.',
      'best_primary_fully_ready_dataset':None,
      'conditional_primary_recommendation':'Kim GCST90016588 for a specifically ever-smoker direct spirometric COPD primary track, GCST90016589 as separate never-smoker stratum; requires stratum-matched LD or justified QCed reference, full-file integrity and effect/test-scale contract. Not a general-population or Sakornsakolpat replacement.',
      'conditional_east_asian_primary_track':'Ishigaki2020 GCST90013709 original physician-diagnosed Japanese COPD, subject to matched Japanese signed LD and full-file QC; sex-stratified releases are sensitivities, not independent replications.',
      'preferred_direct_anchor_if_access_resolved':'Sakornsakolpat2019 GCST007692; request complete signed meta-analysis/ancestry-specific source and compatible LD. Current public export unsuitable.',
      'best_secondary_recommendation':'PanUKB ICD10 J44 EUR, 11536 cases/408995 controls; best documented public same-project signed-LD route, still pending release/scaling/test/subset/variant and locus QC. Not promoted to direct COPD primary.',
      'proposed_primary_method':methods['method_design']['primary_method'],
      'proposed_sensitivity_method':methods['method_design']['sensitivity_method'],
      'fine_mapping_executed':False,'risk_alleles_assigned':False,'candidate_membership_used':False,'model_scoring_executed':False,
      'complete_genomewide_GWAS_or_LD_matrix_values_downloaded':False,'author_contact_sent':False,'commit_or_push_performed':False,
      'design_only_proposals_require_investigator_review':True}
    (STAGE/'provenance/recommendations.json').write_text(json.dumps(recommendations,indent=2)+'\n')
    (STAGE/'provenance/compute_storage_plan.json').write_text(json.dumps(methods['compute_storage_plan'],indent=2)+'\n')
    bindings=[]
    for path in ['audits/direct/direct_source_audit.json','audits/alternatives/alternatives_audit.json','audits/methods_ld/methods_ld_audit.json']:
        b=(STAGE/path).read_bytes();bindings.append({'path':path,'bytes':len(b),'sha256':hashlib.sha256(b).hexdigest()})
    (STAGE/'provenance/assembly_manifest.json').write_text(json.dumps({'created_utc':recommendations['created_utc'],'bound_source_audits':bindings,'study_matrix_rows':104,'provider_rows':25,'direct_signed_rows':7,'scientific_inference_performed':False},indent=2)+'\n')
    print(json.dumps(recommendations,indent=2))
if __name__=='__main__':main()
