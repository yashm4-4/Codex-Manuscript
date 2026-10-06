#!/usr/bin/env python3
"""Compile metadata/header feasibility, never association results or locus calls."""
import collections,csv,datetime,json,pathlib,re
STAGE=pathlib.Path(__file__).resolve().parents[1]
ROOT=STAGE.parents[3]
def table(name,rows,fields=None):
    target=STAGE/'tables'/name;target.parent.mkdir(exist_ok=True)
    fields=fields or list(rows[0])
    with target.open('w') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader()
        for row in rows:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v for k,v in row.items()})
def nonmissing(h,field):
    if field not in h.get('header',[]):return False
    i=h['header'].index(field)
    return any(len(row)>i and row[i].strip().lower() not in ('','na','nan','null','none','.') for row in h.get('representative_records',[]))
def fields_matching(fields,pattern):return ';'.join(x for x in fields if re.search(pattern,x,re.I)) or 'not_present_in_inspected_header'
def main():
    rows=list(csv.DictReader((STAGE/'inputs/frozen_study_metadata.tsv').open(),delimiter='\t'))
    audits={r['accession']:r for r in json.loads((STAGE/'audits/catalog/catalog_audit.json').read_text())}
    headers=[];inventory=[];matrix=[];risk=[];ancestry=[];limitations=[];ld=[]
    for row in rows:
        accession=row['study_accession'];a=audits[accession];m=a['current_catalog_metadata']
        good=[h for h in a['headers'] if h.get('header') and h.get('sample_records',0)>0 and h.get('records_width_matches')]
        # Prefer native author fields; formatted/harmonized releases are separately represented.
        good.sort(key=lambda h:('harmonised/' in h['url'],'.h.tsv' in h['url']))
        selected=good[0] if good else {}
        fields=selected.get('header',[])
        signed=any(nonmissing(selected,f) for f in ('beta','odds_ratio','hm_beta','hm_odds_ratio'))
        alleles=('effect_allele' in fields and 'other_allele' in fields) or ('hm_effect_allele' in fields and 'hm_other_allele' in fields)
        se=any(nonmissing(selected,f) for f in ('standard_error','se','SE'))
        p=any(nonmissing(selected,f) for f in ('p_value','pval','P','neg_log_10_p_value','neglog10_pval'))
        tech=';'.join(m.get('genotyping_technologies',[]))
        phenotype=row['phenotype_class']
        direct=phenotype=='direct_copd_susceptibility'
        if phenotype=='nonstandard_gene_burden_cnv':scope='UNSUITABLE_FOR_VARIANT_LEVEL_REGULATORY_FINE_MAPPING: aggregate gene/burden/CNV tests'
        elif 'Exome' in tech or 'exome' in tech:scope='RESTRICTED_VARIANT_DENSITY: exome ascertainment is not locus-complete regulatory variation'
        elif direct:scope='PRIMARY_PHENOTYPE_ELIGIBLE_CONDITIONAL_ON_STATISTICS_AND_LD'
        elif phenotype=='ehr_copd_susceptibility':scope='SECONDARY_EHR_ONLY_NOT_PRIMARY_DIRECT_COPD'
        elif phenotype=='mixed_ehr_self_report_susceptibility':scope='SEPARATE_COMPOSITE_SENSITIVITY_NOT_PRIMARY'
        elif phenotype=='ml_surrogate_liability':scope='SEPARATE_ML_LIABILITY_ESTIMAND_NOT_PRIMARY_COPD_RISK'
        elif phenotype=='unresolved':scope='BLOCKED_PHENOTYPE_DEFINITION_UNRESOLVED'
        else:scope='SEPARATE_ESTIMAND_NOT_PRIMARY_COPD_SUSCEPTIBILITY'
        if accession=='GCST90246122':access='RESTRICTED_LICENSE: Catalog suspension readme, no statistics exposed'
        elif good:access='PUBLIC_ANONYMOUS_HEADER_AND_RECORDS_VERIFIED_NOT_FULL_FILE_QC'
        elif a['directory_access'].get('status')==404:access='NO_CATALOG_FILE_FOUND: official alternatives/request assessed separately'
        else:access='DIRECTORY_OR_HEADER_ONLY: no representative variant records verified'
        native_meta=[]
        for item in a['metadata_files']:
            file=STAGE/item['path']
            if file.suffix.lower() not in ('.yaml','.yml') or not file.exists():continue
            text=file.read_text()
            native_meta.append({'path':item['path'],'genome_assembly':';'.join(re.findall(r'^genome_assembly:\s*(.+)$',text,re.M)) or 'not_reported',
                                'declared_sample_sizes':re.findall(r'^\s*sample_size:\s*(\d+)',text,re.M),
                                'data_file_name':';'.join(re.findall(r'^data_file_name:\s*(.+)$',text,re.M))})
        native_builds=[x['genome_assembly'] for x in native_meta if '/harmonised_' not in x['path']]
        builds=';'.join(sorted(set(native_builds))) or 'unverified_no_file_metadata'
        # Catalog ancestry can combine discovery and replication; never pretend it is an LD matching proof.
        ancestry_text=row['catalog_ancestry_categories']
        if any(x in ancestry_text for x in (';',',','NR')):ancestry_type='MIXED_OR_MULTIPLE_CATALOG_LABELS: analysis-specific ancestry weights require verification'
        else:ancestry_type='SINGLE_REPORTED_ANCESTRY: cohort/stratum-specific matching still required'
        cohort=row['cohorts']
        if 'UKB' in cohort and ancestry_text=='European':ld_candidate='Pan-UKB EUR cohort-source LD; compare exact UKB subset, variant set, covariates and build; not automatically study-matched'
        elif 'BBJ' in cohort or 'BioBank' in cohort:ld_candidate='BBJ ancestry/cohort-source signed LD if retrievable; exact release and disease/control sample mismatch require review'
        elif 'MVP' in cohort:ld_candidate='MVP ancestry-specific signed LD preferred; UKB or 1000G is not cohort-matched'
        else:ld_candidate='Study/cohort/ancestry-matched signed LD not established by Catalog; see source-specific LD audit'
        if ancestry_type.startswith('MIXED'):ld_candidate+='; do not substitute single-European LD for pooled mixed-ancestry statistics'
        blockers=[]
        if not good:blockers.append('public representative variant statistics unverified at Catalog')
        if not signed or not alleles:blockers.append('signed effect/effect-allele mapping not demonstrated in inspected usable file')
        if not se:blockers.append('nonmissing SE not demonstrated in inspected usable file')
        if not p:blockers.append('nonmissing P not demonstrated in inspected usable file')
        if scope.startswith(('UNSUITABLE','RESTRICTED','BLOCKED')):blockers.append(scope)
        blockers+=['whole-file completeness and integrity not tested by bounded prefix','per-locus variant density/overlap and signed LD consistency not tested','source-specific LD compatibility/exact access must pass before execution']
        rec={'accession':accession,'pubmed_id':row['pubmed_id'],'publication':row['publication'],'doi':row['doi'],
          'frozen_phenotype_class':phenotype,'case_definition':row['case_definition'],'control_definition':row['control_definition'],
          'smoking_design':row['smoking_design'],'sex_design':row['sex_design'],'cohorts':cohort,'ancestry_labels_frozen':ancestry_text,
          'discovery_ancestry_current':m.get('discovery_ancestry',[]),'replication_ancestry_current':m.get('replication_ancestry',[]),
          'ancestry_specific_vs_meta':ancestry_type,'initial_sample_size_current':m.get('initial_sample_size','unverified'),
          'replication_sample_size_current':m.get('replication_sample_size','unverified'),'genotyping_technology_current':tech,
          'tested_variant_count_catalog_metadata_only':m.get('snp_count','not_reported'),'imputed_catalog':m.get('imputed','not_reported'),
          'native_genome_build':builds,'file_build_metadata':native_meta,'catalog_flag_current_not_usability_proof':m.get('full_summary_stats_available'),
          'catalog_api_http_status':a['catalog_api_access'].get('status'),'directory_http_status':a['directory_access'].get('status'),
          'public_access_status':access,'source_directory':a['directory_url'],'selected_header_source_url':selected.get('url','none_verified'),
          'selected_prefix_path':selected.get('prefix_path','none_verified'),'selected_prefix_sha256':selected.get('prefix_sha256','not_applicable'),
          'inspected_header':fields,'variant_identifier_fields':fields_matching(fields,r'variant|rsid|^rs_id$|^id$|^Name$'),
          'position_fields':fields_matching(fields,r'chrom|pos|location'),'ref_alt_fields':fields_matching(fields,r'^ref$|^alt$'),
          'effect_other_allele_fields':fields_matching(fields,r'allele'),'signed_effect_fields':fields_matching(fields,r'^beta$|odds_ratio|^hm_beta$'),
          'SE_fields':fields_matching(fields,r'^standard_error$|^se$'),'P_fields':fields_matching(fields,r'^p_value$|^pval$|^P$|neg_log_10_p_value|neglog10_pval'),
          'frequency_fields':fields_matching(fields,r'frequency|(^|_)af$'),'per_variant_N_fields':fields_matching(fields,r'^n$|^n_|_n$|num_cases|num_controls'),
          'variant_quality_fields':fields_matching(fields,r'info|r2|qc|quality|is_|errcode|firth|converged'),
          'signed_effect_nonmissing_in_sample':signed,'effect_allele_fields_present':alleles,'SE_nonmissing_in_sample':se,
          'P_nonmissing_in_sample':p,'phenotype_eligibility':scope,
          'complete_locus_statistics_status':'not validated genome-wide/locus-wide; metadata and prefix evidence only',
          'LD_candidate_and_caveat':ld_candidate,'fine_mapping_execution_status':'NOT_EXECUTION_CLEARED_BY_THIS_PREFLIGHT',
          'blockers':blockers,'exact_evidence_record':f'audits/catalog/{accession}/audit.json',
          'phenotype_source_urls':row['source_urls'],'frozen_uncertainties':row['uncertainties']}
        matrix.append(rec)
        risk.append({k:rec[k] for k in ['accession','frozen_phenotype_class','selected_header_source_url','effect_other_allele_fields','signed_effect_fields','SE_fields','P_fields','frequency_fields','signed_effect_nonmissing_in_sample','effect_allele_fields_present']})
        risk_status='ENCODED_EFFECT_SIGN_AVAILABLE; source outcome coding/scale and allele harmonization still mandatory' if signed and alleles else 'UNVERIFIED_OR_NOT_VARIANT_SIGNED'
        if phenotype=='nonstandard_gene_burden_cnv':risk_status='NOT_A_VARIANT_RISK_ALLELE: signed aggregate-test statistics do not define an individual variant effect'
        risk[-1].update({'risk_direction_feasibility':risk_status,
                        'risk_alleles_assigned':False,'not_disease_risk_if_non_susceptibility_estimand':not (direct or phenotype=='ehr_copd_susceptibility')})
        ancestry.append({k:rec[k] for k in ['accession','cohorts','ancestry_labels_frozen','discovery_ancestry_current','replication_ancestry_current','ancestry_specific_vs_meta','initial_sample_size_current','replication_sample_size_current','smoking_design','sex_design']})
        ancestry[-1].update({'overlap_notes':row['overlap_notes'],'sample_size_caveat':'Never sum repeated YAML sample entries, overlapping strata, or discovery/replication indiscriminately; row-level N and effective N remain future checks.'})
        limitations.append({'accession':accession,'phenotype_eligibility':scope,'access_status':access,'blockers':blockers,'frozen_source_limitations':row['uncertainties']})
        ld.append({'accession':accession,'ancestry':ancestry_text,'cohort':cohort,'provisional_LD_candidate':ld_candidate,
                   'matched_signed_matrix_acquired':False,'locus_matrix_QC_executed':False,'compatibility_status':'SOURCE_SPECIFIC_AUDIT_REQUIRED; no automatic eligibility from ancestry label'})
        for i,h in enumerate(a['headers'],1):
            headers.append({'accession':accession,'file_number':i,**h})
        for item in a['file_inventory']:inventory.append({'accession':accession,**item})
    table('study_phenotype_summary_stat_eligibility.tsv',matrix)
    table('header_schema_audit.tsv',headers,sorted(set().union(*(x.keys() for x in headers))))
    table('source_accession_file_inventory.tsv',inventory)
    table('catalog_signed_effect_risk_direction_feasibility.tsv',risk)
    table('catalog_ancestry_cohort_sample_size.tsv',ancestry)
    table('catalog_study_specific_limitations.tsv',limitations)
    table('catalog_study_LD_requirements.tsv',ld)
    stats={'stage':'fine-mapping-preflight-1.0','generated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
      'core_accessions':len(rows),'phenotype_classes':dict(collections.Counter(r['phenotype_class'] for r in rows)),
      'api_200':sum(a['catalog_api_access'].get('status')==200 for a in audits.values()),
      'directory_200':sum(a['directory_access'].get('status')==200 for a in audits.values()),
      'directory_404':sum(a['directory_access'].get('status')==404 for a in audits.values()),
      'accessions_with_inspected_headers':sum(bool(a['headers']) for a in audits.values()),'file_headers':len(headers),
      'accessions_with_selected_header_and_records':sum(r['selected_prefix_path']!='none_verified' for r in matrix),
      'direct_accessions_with_Catalog_signed_header':sum(r['frozen_phenotype_class']=='direct_copd_susceptibility' and r['signed_effect_nonmissing_in_sample'] and r['effect_allele_fields_present'] for r in matrix),
      'candidate_selection_inputs':[],'association_or_locus_analysis_performed':False,'fine_mapping_run':False}
    (STAGE/'provenance/catalog_table_manifest.json').write_text(json.dumps(stats,indent=2)+'\n');print(json.dumps(stats,indent=2))
if __name__=='__main__':main()
