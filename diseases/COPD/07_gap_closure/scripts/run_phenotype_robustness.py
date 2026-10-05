#!/usr/bin/env python3
"""Offline V2 phenotype provenance sensitivity; never writes outside V2.

No external packages, statistical fitting, allele scoring, or LD calculations.
Publication adjudication is an explicit, separately reviewed JSON input.
"""
from __future__ import annotations

import csv
import gzip
import hashlib
import io
import itertools
import json
import platform
import resource
import subprocess
import sys
import time
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

V2 = Path(__file__).resolve().parents[1]
COPD = V2.parent
ROOT = COPD.parent.parent
RESULTS = V2 / 'results'
PREFIX = 'COPD-V2-PHENO'
ML = 'GCST90244098'
INPUTS = {}
OUTPUTS = []
CHECKS = []


def digest(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read_table(rel):
    path = COPD / rel
    INPUTS[str(path.relative_to(ROOT))] = digest(path)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt', newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def tok(value):
    return {x for x in str(value or '').split(';') if x and x not in {'nan', 'not_applicable'}}


def yes(value):
    return str(value).lower() in {'true', 'yes', '1'}


def joined(values):
    return ';'.join(sorted(set(values)))


def scalar(value):
    if isinstance(value, (list, set, tuple)):
        return joined(str(x) for x in value)
    if isinstance(value, dict):
        return json.dumps(value, sort_keys=True, ensure_ascii=False)
    if value is None:
        return 'not_reported'
    return value


def save(name, rows, fields=None):
    path = RESULTS / f'{PREFIX}-{name}'
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(rows)
    fields = fields or list(rows[0])
    buf = io.StringIO(newline='')
    writer = csv.DictWriter(buf, fieldnames=fields, delimiter='\t', lineterminator='\n', extrasaction='raise')
    writer.writeheader()
    for row in rows:
        writer.writerow({k: scalar(v) for k, v in row.items()})
    content = buf.getvalue().encode()
    if path.suffix == '.gz':
        with open(path, 'wb') as f:
            with gzip.GzipFile(fileobj=f, mode='wb', filename='', mtime=0) as z:
                z.write(content)
    else:
        path.write_bytes(content)
    OUTPUTS.append({'path': str(path.relative_to(ROOT)), 'sha256': digest(path), 'rows': len(rows)})
    return path


def check(name, condition, detail):
    CHECKS.append({'check_id': f'{PREFIX}-QC-{len(CHECKS)+1:03d}', 'check': name,
                   'status': 'PASS' if condition else 'FAIL', 'detail': detail})
    if not condition:
        raise AssertionError(f'{name}: {detail}')


def main():
    started=time.perf_counter()
    base = read_table('02_gwas/results/COPD-S2-R001_studies_core.tsv')
    gws = read_table('02_gwas/results/COPD-S2-R002_gws_associations.tsv')
    tags = read_table('02_gwas/results/COPD-S2-R002_gws_unique_tag_variants.tsv')
    ancestry = read_table('02_gwas/data/copd_ancestries.tsv')
    focal = read_table('02_gwas/results/COPD-S2-R006A_focal_match_audit.tsv')
    panels = read_table('02_gwas/results/COPD-S2-R006B_focal_panel_audit.tsv')
    pairs = read_table('02_gwas/results/COPD-S2-R006C_ld_pairs.tsv.gz')
    candidates = read_table('02_gwas/results/COPD-S2-R006E_candidate_variants_grch38.tsv.gz')
    frozen = read_table('04_modeling/results/COPD-S4-R010_THE_LIST.tsv')
    annotation = read_table('03_regulatory_landscape/results/COPD-S3-R003_candidate_classification.tsv.gz')
    locus_definitions = read_table('03_regulatory_landscape/results/COPD-S3-R002_gene_loci_grch38.tsv')
    selection = read_table('06_experimental_validation/results/COPD-S6-R002_selection_audit.tsv.gz')
    shortlist = read_table('06_experimental_validation/results/COPD-S6-R001_candidate_shortlist.tsv')
    for rel in ['manuscript/COPD_regulatory_genomics_manuscript.md', 'manuscript/results/COPD-workflow-validation.tsv']:
        path = COPD / rel
        INPUTS[str(path.relative_to(ROOT))] = digest(path)
    adjudications = []
    for name in ['adjudication_pre2020.json', 'adjudication_2020_2022.json', 'adjudication_2023plus.json']:
        path = V2 / 'data' / name
        INPUTS[str(path.relative_to(ROOT))] = digest(path)
        adjudications.extend(json.loads(path.read_text()))
    # Optional researcher-reviewed amendments are explicit input, never silent overrides.
    amendment = V2 / 'data' / 'phenotype_adjudication_amendments.json'
    if amendment.exists():
        INPUTS[str(amendment.relative_to(ROOT))] = digest(amendment)
        amend = json.loads(amendment.read_text())
        for a in adjudications:
            if a['study_accession'] in amend:
                a.update(amend[a['study_accession']])
    studies = {a['study_accession']: a for a in adjudications}
    check('104 unique adjudicated accessions equal frozen core', len(studies) == len(adjudications) == 104 and set(studies) == {r['STUDY ACCESSION'] for r in base}, str(len(studies)))
    required = ['doi', 'phenotype_class', 'primary_eligible', 'secondary_eligible', 'primary_question', 'case_definition', 'control_definition', 'copd_modeled_outcome', 'ascertainment', 'spirometry', 'clinical_diagnosis', 'ehr_definition', 'surrogate_definition', 'severity_selected', 'within_case_design', 'subtype_secondary_phenotype', 'smoking_design', 'sex_design', 'discovery_replication', 'cohorts', 'overlap_notes', 'full_stats_status', 'full_stats_url', 'rationale', 'confidence', 'source_urls', 'source_evidence', 'uncertainties']
    missing = [(s, k) for s, a in studies.items() for k in required if k not in a]
    check('All requested adjudication fields supplied', not missing, str(missing))
    check('Primary and secondary eligibility mutually exclusive', all(not (yes(a['primary_eligible']) and yes(a['secondary_eligible'])) for a in studies.values()), 'No combined primary/EHR coding')
    check('Primary eligibility uses direct phenotype class', all(not yes(a['primary_eligible']) or a['phenotype_class'] == 'direct_copd_susceptibility' for a in studies.values()), 'Case-control definition is a manual literature input')
    legal_classes = {'direct_copd_susceptibility','ehr_copd_susceptibility','ml_surrogate_liability','lung_function','progression_decline','within_copd_severity','emphysema','chronic_bronchitis','pulmonary_vascular_within_copd','exacerbation_treatment','materially_different_smoking_interaction','shared_cross_disease','nonstandard_gene_burden_cnv','unresolved','mixed_direct_ehr_susceptibility'}
    legal_classes.add('mixed_ehr_self_report_susceptibility')
    check('Only declared phenotype classes used', all(a['phenotype_class'] in legal_classes for a in studies.values()), str(sorted({a['phenotype_class'] for a in studies.values()})))
    check('Secondary eligibility uses EHR phenotype class', all(not yes(a['secondary_eligible']) or a['phenotype_class']=='ehr_copd_susceptibility' for a in studies.values()), 'Mixed and unresolved definitions not forced into EHR')
    check('Sources and cohorts have list representation', all(isinstance(a['cohorts'],list) and isinstance(a['source_urls'],list) and len(a['source_urls'])>0 for a in studies.values()), 'Auditable source links for every accession')
    for a in studies.values():
        a['cohorts_original_review_labels']=list(a['cohorts'])
        a['cohorts']=sorted({{'BEOCOPD':'EOCOPD','EB':'EstBB'}.get(c,c) for c in a['cohorts']})
        a['full_stats_status_original_review']=a['full_stats_status']
        a['full_stats_status']=a['full_stats_status'].replace('not_available_in_frozen_Catalog','not_reported_in_frozen_Catalog').replace('not_available_in_Catalog','not_reported_in_Catalog')
        if a['full_stats_url']=='not_available': a['full_stats_url']='not_reported'
        if a['full_stats_status']=='public_full_statistics_reported_in_Catalog; not_downloaded':
            a['full_stats_status']='Catalog reports full statistics for tested units; accession endpoint/header not verified; complete common-SNP coverage not established'
        if a['study_accession'] in {'GCST90013709','GCST90013746','GCST90013781'}:
            a['full_stats_status']+='; publication-wide 42-disease release established, accession/sex-specific file not verified'
        if a['study_accession']=='GCST007692':
            a['full_stats_status']='not_reported_in_Catalog; paper cites controlled/application dbGaP and UKB routes; availability of a complete meta-analysis summary file not verified; component access does not reconstruct meta-analysis'
        if a['full_stats_status']=='available: GWAS Catalog full summary statistics':
            a['full_stats_status']='Catalog reports full statistics for tested units; complete accession file not downloaded or validated in this module'
        for k in ['severity_selected','within_case_design','subtype_secondary_phenotype','smoking_design','sex_design','uncertainties']:
            if str(a[k]).startswith('negative:'):
                a[k]=a[k][len('negative:'):].strip()
                if 'reported' in a[k] and k in {'severity_selected','smoking_design'}:
                    a[k]='not_reported: '+a[k]
                elif k in {'within_case_design','subtype_secondary_phenotype'}:
                    a[k]='not_applicable: '+a[k]
    sig = defaultdict(set)
    best = {}
    gws_count = Counter()
    for row in gws:
        key = row['normalized_variant_id'], row['STUDY ACCESSION']
        sig[key[0]].add(key[1])
        best[key] = min(best.get(key, 1), float(row['P']))
        gws_count[key[1]] += 1
    agg = {r['normalized_variant_id']: tok(r['study_accessions']) for r in tags}
    check('Frozen universe dimensions', (len(tags), len(candidates), len(frozen), len(shortlist)) == (660,15389,337,12), f'{len(tags)}/{len(candidates)}/{len(frozen)}/{len(shortlist)}')
    check('All GWS tags preserved in association rows', set(sig) == set(agg), str(len(sig)))
    check('GWS support uses unchanged threshold', all(float(r['P']) <= 5e-8 for r in gws), 'P <= 5e-8')
    groups = {'broad': set(studies), 'primary_direct': {s for s,a in studies.items() if yes(a['primary_eligible'])},
              'secondary_ehr': {s for s,a in studies.items() if yes(a['secondary_eligible'])},
              'ml_surrogate': {s for s,a in studies.items() if a['phenotype_class']=='ml_surrogate_liability'},
              'unresolved': {s for s,a in studies.items() if a['phenotype_class']=='unresolved'}}
    groups['other_nonprimary'] = set(studies) - groups['primary_direct'] - groups['secondary_ehr'] - groups['ml_surrogate'] - groups['unresolved']
    groups['primary_or_ehr'] = groups['primary_direct'] | groups['secondary_ehr']
    groups['without_GCST90244098'] = set(studies) - {ML}
    groups['GCST90244098'] = {ML}
    groups['sensitivity_primary_high_confidence_only'] = {s for s in groups['primary_direct'] if studies[s]['confidence']=='high'}
    groups['sensitivity_primary_plus_conditional_smoking_coefficients'] = groups['primary_direct'] | {s for s,a in studies.items() if a['phenotype_class']=='materially_different_smoking_interaction'}
    groups['sensitivity_primary_plus_mixed_direct_ehr'] = groups['primary_direct'] | {s for s,a in studies.items() if a['phenotype_class']=='mixed_direct_ehr_susceptibility'}
    groups['sensitivity_ehr_plus_GCST007996'] = groups['secondary_ehr'] | {'GCST007996'}
    groups['sensitivity_ehr_plus_documented_selfreport_composites'] = groups['secondary_ehr'] | {s for s,a in studies.items() if a['phenotype_class']=='mixed_ehr_self_report_susceptibility'}
    classes = sorted({a['phenotype_class'] for a in studies.values()})
    for cl in classes:
        groups[f'class:{cl}'] = {s for s,a in studies.items() if a['phenotype_class']==cl}
    panelmap = {}
    audit_rows = []
    paircounts = Counter((r['focal_tag_variant'], r['panel']) for r in pairs)
    ancestry_by_study = defaultdict(set)
    ancestry_stages = defaultdict(set)
    for r in ancestry:
        ancestry_by_study[r['STUDY ACCESSION']].add(r['BROAD ANCESTRAL CATEGORY'])
        ancestry_stages[r['STUDY ACCESSION']].add(f"{r['STAGE']}:{r['BROAD ANCESTRAL CATEGORY']}")
    for r in panels:
        key = r['normalized_variant_id'], r['panel']
        supported = tok(r['supporting_studies'])
        if not supported <= sig[key[0]]:
            raise AssertionError(f'Panel audit promoted nonsignificant support: {key}')
        panelmap[key] = supported
        for s in sorted(supported):
            a=studies[s]
            audit_rows.append({'study_accession':s, 'focal_tag_id':key[0], 'ancestry_panel':key[1],
                'study_tag_best_p':best[(key[0],s)], 'phenotype_class':a['phenotype_class'],
                'primary_eligible':yes(a['primary_eligible']), 'secondary_eligible':yes(a['secondary_eligible']),
                'matched':r['matched'], 'match_status':r['match_status'], 'panel_evaluation_status':r['panel_evaluation_status'],
                'panel_variant_id':r['panel_variant_id'], 'n_existing_proxy_links':paircounts[key],
                'supporting_ancestry_categories':r['supporting_ancestry_categories'],
                'cohort_overlap_groups':a['cohorts'], 'source_result':'COPD-S2-R006B;COPD-S2-R006C'})
    check('Every LD pair has a valid frozen tag-panel audit', all((r['focal_tag_variant'],r['panel']) in panelmap for r in pairs), f'{len(pairs)} LD links')
    check('Tag-panel audit keys unique',len(panelmap)==len(panels),f'{len(panelmap)} unique tag-panel keys')
    panel_status={(r['normalized_variant_id'],r['panel']):r['panel_evaluation_status'] for r in panels}
    check('Every LD link arises from a matched polymorphic tag-panel',all(panel_status[(r['focal_tag_variant'],r['panel'])]=='matched_polymorphic_in_panel' for r in pairs),'Unmatched and monomorphic routes never LD-expanded')
    support = defaultdict(set)
    routes = defaultdict(set)
    legacy = {}
    candidate_ids = {r['candidate_record_id'] for r in candidates}
    for row in candidates:
        cid=row['candidate_record_id']; legacy[cid]=set()
        for t in tok(row['gws_tag_ids']):
            support[cid].update(sig[t])
            for s in sig[t]: routes[cid].add((s,t,'DIRECT_TAG'))
        for t in tok(row['gws_tag_ids']) | tok(row['source_focal_tags']):
            legacy[cid].update(agg[t])
    for row in pairs:
        cid=row['proxy_panel_variant_id']; t=row['focal_tag_variant']; panel=row['panel']
        for s in panelmap[(t,panel)]:
            support[cid].add(s); routes[cid].add((s,t,panel))
    check('Broad study/tag/panel reconstruction recovers all 15389 records', set(support)==candidate_ids and all(support[c] for c in candidate_ids), str(len(support)))
    check('Exact reconstructed support is contained in V1 aggregate provenance', all(support[c] <= legacy[c] for c in candidate_ids), 'Extra nonsignificant or wrong-panel links only remain in legacy columns')
    check('Frozen R010 legacy provenance reproduced', all(legacy[r['candidate_record_id']]==tok(r['linked_study_accessions']) for r in frozen), '337 candidate linked_study_accessions sets reproduced')
    # Independent inverse construction: select studies first, then tag-panel pairs,
    # then retained candidates. This differs from candidate-first annotation above.
    primary=groups['primary_direct']; retained_tags={t for t in sig if sig[t] & primary}
    retained_tp={key for key,ss in panelmap.items() if ss & primary}
    independently_retained={r['candidate_record_id'] for r in candidates if tok(r['gws_tag_ids']) & retained_tags}
    independently_retained.update(r['proxy_panel_variant_id'] for r in pairs if (r['focal_tag_variant'],r['panel']) in retained_tp)
    check('Independent study-first reconstruction agrees', independently_retained=={c for c in candidate_ids if support[c]&primary}, 'Study-first tag-panel filter equals candidate-first route projection')
    for group,eligible in groups.items():
        selected_tags={t for t,ss in sig.items() if ss&eligible}
        selected_panels={key for key,ss in panelmap.items() if ss&eligible}
        inverse={r['candidate_record_id'] for r in candidates if tok(r['gws_tag_ids'])&selected_tags}
        inverse.update(r['proxy_panel_variant_id'] for r in pairs if (r['focal_tag_variant'],r['panel']) in selected_panels)
        check(f'Independent reconstruction: {group}',inverse=={c for c in candidate_ids if support[c]&eligible},f'{len(inverse)} records')
    # Audit a tempting but incorrect tag-only ancestry filter, quantitatively.
    naive={r['candidate_record_id'] for r in candidates if tok(r['gws_tag_ids']) & retained_tags}
    naive.update(r['proxy_panel_variant_id'] for r in pairs if r['focal_tag_variant'] in retained_tags)
    naive_excess=naive-independently_retained
    focalmap={r['normalized_variant_id']:r for r in focal}
    frozenmap={r['candidate_record_id']:r for r in frozen}
    ann={r['candidate_record_id']:r for r in annotation}
    selmap={r['candidate_record_id']:r for r in selection}

    def fields(ss, original=None):
        original=ss if original is None else original
        cl={studies[s]['phenotype_class'] for s in ss}
        p=ss & groups['primary_direct']; e=ss & groups['secondary_ehr']; m=ss & groups['ml_surrogate']; u=ss & groups['unresolved']; o=ss & groups['other_nonprimary']
        if p: state='retained_by_primary_direct_copd'
        elif len(cl)>1: state='supported_by_multiple_nonprimary_classes'
        elif e: state='secondary_ehr_only'
        elif m: state='ml_surrogate_only'
        elif o: state='other_nonprimary_only'
        else: state='unresolved'
        co=set().union(*(set(studies[s]['cohorts']) for s in ss)) if ss else set()
        return {'support_state':state,'primary_retained':bool(p),'secondary_ehr_supported':bool(e),
                'ml_supported':bool(m),'other_nonprimary_supported':bool(o),'unresolved_study_supported':bool(u),
                'multiple_phenotype_classes':len(cl)>1,'phenotype_classes':joined(cl),
                'gws_supporting_studies':joined(ss),'n_supporting_accessions_not_independent':len(ss),
                'primary_supporting_studies':joined(p),'ehr_supporting_studies':joined(e),
                'ml_supporting_studies':joined(m),'other_supporting_studies':joined(o),
                'unresolved_supporting_studies':joined(u),'cohort_overlap_groups':joined(co),
                'GCST90244098_supported':ML in ss,'GCST90244098_only_study':ss=={ML},
                'v1_aggregate_linked_studies':joined(original),
                'aggregate_studies_not_supported_by_valid_gws_route':joined(original-ss)}

    allrows=[]
    for i,row in enumerate(candidates,1):
        cid=row['candidate_record_id']; f=frozenmap.get(cid,{})
        out={'v2_record_id':f'{PREFIX}-C{i:05d}', 'candidate_record_id':cid,
             'v1_priority_rank':f.get('predicted_causal_priority_rank','not_applicable'),
             'in_frozen_337':cid in frozenmap, **{k:row[k] for k in ['chromosome_grch38','position_grch38','ref','alt','variant_class','candidate_origin','reference_matched','gws_tag_ids','source_focal_tags','ld_panels']},
             **fields(support[cid],legacy[cid]),
             'primary_valid_study_tag_panel_routes':joined('|'.join(x) for x in routes[cid] if x[0] in primary),
             'ehr_valid_study_tag_panel_routes':joined('|'.join(x) for x in routes[cid] if x[0] in groups['secondary_ehr']),
             'all_valid_study_tag_panel_routes':joined('|'.join(x) for x in routes[cid]),
             'would_be_retained_by_incorrect_tag_only_panel_filter':cid in naive,
             'false_retention_from_ignoring_panel':cid in naive_excess,
             'v1_component_id':selmap.get(cid,{}).get('ld_component_id','not_applicable'),
             'v1_selected_loci':f.get('assigned_selected_loci','not_applicable'),
             'v1_linked_gwas_genes':f.get('linked_gwas_genes','not_applicable'),
             'v1_model_context':f.get('model_context','not_applicable'),
             'v1_exclusive_annotation':ann[cid]['exclusive_class_comprehensive'],
             'v1_coding_CDS':ann[cid]['coding_CDS'], 'v1_donor_refined_any':ann[cid]['donor_refined_any'],
             'v1_known_regulatory_any':ann[cid]['known_regulatory_any']}
        allrows.append(out)
    rowmap={r['candidate_record_id']:r for r in allrows}
    focused=[rowmap[r['candidate_record_id']] for r in frozen]
    check('Original 337 rank/order preserved', [int(r['v1_priority_rank']) for r in focused]==list(range(1,338)), 'No reranking')
    check('Candidate and selection identities unique and complete',len(candidate_ids)==len(candidates) and len(frozenmap)==len(frozen) and len(selmap)==len(selection) and set(selmap)==set(frozenmap),'No duplicate candidate or selection keys')
    shortrows=[]
    for r in shortlist:
        shortrows.append({'v1_shortlist_rank':r['experimental_shortlist_rank'],
                          'v1_selection_anchor':r['selection_anchor_criterion'], **rowmap[r['candidate_record_id']]})
    comp=defaultdict(list)
    for r in selection: comp[r['ld_component_id']].append(r['candidate_record_id'])
    check('Frozen components preserved', len(comp)==153 and sum(map(len,comp.values()))==337, f'{len(comp)} components, {sum(map(len,comp.values()))} members')
    check('Frozen component membership strings reproduced',all(tok(r['ld_component_members'])==set(comp[r['ld_component_id']]) for r in selection),'Members retained exactly; components not rebuilt')
    components=[]
    for key,ids in sorted(comp.items()):
        ss=set().union(*(support[c] for c in ids)); ls=set().union(*(legacy[c] for c in ids))
        nr=sum(bool(support[c]&primary) for c in ids)
        components.append({'v1_component_id':key,'n_frozen_members':len(ids),
            'v1_member_ids_original_order':';'.join(sorted(ids,key=lambda c:int(frozenmap[c]['predicted_causal_priority_rank']))),
            'v1_member_ranks':';'.join(str(frozenmap[c]['predicted_causal_priority_rank']) for c in sorted(ids,key=lambda c:int(frozenmap[c]['predicted_causal_priority_rank']))),
            'n_primary_retained_members':nr,'primary_component_state':'all_members' if nr==len(ids) else ('some_members' if nr else 'no_members'),
            'n_ML_linked_members':sum(ML in support[c] for c in ids),
            'n_GCST90244098_sole_accession_members':sum(support[c]=={ML} for c in ids),
            'n_ML_class_only_members':sum(bool(support[c]) and support[c]<=groups['ml_surrogate'] for c in ids), **fields(ss,ls)})
    tagrows=[]
    for r in tags:
        t=r['normalized_variant_id']; a=focalmap[t]
        tagrows.append({'tag_id':t,'best_p':r['best_p'],'v1_consequence_class':r['consequence_class'],
                        'v1_mapped_genes':r['mapped_genes'],'v1_matched':a['matched'],
                        'v1_match_status':a['match_status'],'v1_panel_variant_id':a['panel_variant_id'],
                        'v1_supported_panels':a['supported_panels'],**fields(sig[t],agg[t])})
    # Per-accession provenance and counts are joined only after classifications exist.
    register=[]; source_rows=[]; cohort_rows=[]
    for i,b in enumerate(base,1):
        s=b['STUDY ACCESSION']; a=studies[s]
        rec={'v2_adjudication_id':f'{PREFIX}-ADJ-{i:03d}','study_accession':s,
             'pubmed_id':b['PUBMED ID'],'publication':b['STUDY'],'publication_date':b['DATE'],
             'reported_phenotype':b['DISEASE/TRAIT'], 'initial_sample_description':b['INITIAL SAMPLE SIZE'],
             'replication_sample_description':b['REPLICATION SAMPLE SIZE'] or 'not_reported',
             **{k:a[k] for k in required}, 'cohorts_original_review_labels':a['cohorts_original_review_labels'],
             'full_stats_status_original_review':a['full_stats_status_original_review'],
             'availability_normalization_rule':'Catalog no means not reported, not globally unavailable; full statistics restricted to actual tested units; unverified endpoints explicit',
             'cohort_alias_rule':'BEOCOPD -> EOCOPD (Boston Early-Onset COPD); EB -> EstBB (GBMI roster/Catalog correspondence, not downloaded supplement); no participant deduplication',
             'catalog_cohorts':b['COHORT'] or 'not_reported',
             'catalog_ancestry_categories':joined(ancestry_by_study[s]) or 'not_reported',
             'catalog_ancestry_by_stage':joined(ancestry_stages[s]) or 'not_reported',
             'catalog_full_summary_statistics_flag':b['FULL SUMMARY STATISTICS'],
             'catalog_full_stats_location':b['SUMMARY STATS LOCATION'] or 'not_reported',
             'v1_gws_association_rows':gws_count[s],'v1_gws_tags':sum(s in ss for ss in sig.values()),
             'v1_candidate_records_supported':sum(s in support[c] for c in candidate_ids),
             'v1_frozen_337_supported':sum(s in support[c] for c in frozenmap),
             'association_coverage_caveat':'0 means no curated V1 GWS contribution, not a negative full GWAS',
             'adjudicated_date':'2026-10-05'}
        register.append(rec)
        for j,url in enumerate(a['source_urls'],1):
            source_rows.append({'source_id':f'{PREFIX}-SRC-{i:03d}-{j:02d}','study_accession':s,'pubmed_id':b['PUBMED ID'],
                'doi':a['doi'],'url':url,'checked_date':'2026-10-05','evidence_paraphrase':a['source_evidence'],
                'phenotype_class':a['phenotype_class'],'confidence':a['confidence'],'uncertainties':a['uncertainties']})
        for co in a['cohorts'] or ['not_reported']:
            cohort_rows.append({'study_accession':s,'cohort_overlap_group':co,'publication_pmid':b['PUBMED ID'],
                'phenotype_class':a['phenotype_class'],'ancestry':joined(ancestry_by_study[s]),'sex_design':a['sex_design'],
                'group_kind':'consortium_umbrella_not_independent_cohort' if co in {'ICGC','CHARGE','UKECC'} else 'named_cohort_or_biobank',
                'overlap_evidence':a['overlap_notes'],'overlap_interpretation':'same cohort provenance; exact participant overlap not identifiable',
                'source_urls':a['source_urls']})
    overlap=[]
    for left,right in itertools.combinations(sorted(studies),2):
        shared=set(studies[left]['cohorts']) & set(studies[right]['cohorts']) - {'not_reported','unclear'}
        if shared:
            overlap.append({'study_accession_a':left,'study_accession_b':right,'shared_cohort_groups':joined(shared),
                            'interpretation':'known shared cohort provenance; participant overlap probable, extent unavailable; not independent replication'})
    universes={'gws_tags':{t:sig[t] for t in sig},'all_15389':dict(support),
               'frozen_337':{c:support[c] for c in frozenmap},
               'components_153':{k:set().union(*(support[c] for c in ids)) for k,ids in comp.items()},
               'shortlist_12':{r['candidate_record_id']:support[r['candidate_record_id']] for r in shortlist}}
    summary=[]
    for name,universe in universes.items():
        for group,eligible in groups.items():
            n=sum(bool(ss & eligible) for ss in universe.values())
            summary.append({'universe':name,'stratum':group,'supported':n,'denominator':len(universe),
                            'percent':round(n/len(universe)*100,4),'counting_unit':'existing component any-member' if name=='components_153' else name,
                            'interpretation':'nonexclusive descriptive support; no independence or causal claim'})
    mlrows=[]
    legacy_universes={'gws_tags':agg,'all_15389':legacy,'frozen_337':{c:legacy[c] for c in frozenmap},
                     'components_153':{k:set().union(*(legacy[c] for c in ids)) for k,ids in comp.items()},
                     'shortlist_12':{r['candidate_record_id']:legacy[r['candidate_record_id']] for r in shortlist}}
    for name,universe in universes.items():
        for mode,umap in [('exact_GWS_tag_panel',universe),('frozen_aggregate_provenance',legacy_universes[name])]:
            linked={k for k,ss in umap.items() if ML in ss}; only={k for k,ss in umap.items() if ss=={ML}}
            primary_keys={k for k,ss in universe.items() if ss&primary}
            mlrows.append({'universe':name,'provenance_mode':mode,'denominator':len(umap),'ML_linked':len(linked),'ML_only_study':len(only),
                           'ML_linked_primary_retained':len(linked&primary_keys),'ML_linked_not_primary_retained':len(linked-primary_keys),
                           'ML_only_primary_retained':len(only&primary_keys)})
    locus_members=defaultdict(set)
    for r in frozen:
        for locus in tok(r['assigned_selected_loci']): locus_members[locus].add(r['candidate_record_id'])
    locusrows=[]
    locusmeta={r['gene']:r for r in locus_definitions}
    check('Frozen locus universes identified',len(locusmeta)==152 and sum(yes(r['is_replicated']) for r in locus_definitions)==140,'152 all-selected; 140 repeated-Catalog-mapping labels, not independent replication')
    for locus,ids in sorted(locus_members.items()):
        ss=set().union(*(support[c] for c in ids))
        locusrows.append({'v1_selected_locus_label':locus,'n_frozen_337_members':len(ids),
            'v1_is_replicated_catalog_mapping':yes(locusmeta[locus]['is_replicated']),
            'v1_locus_selection_reason':locusmeta[locus]['selection_reason'],
            'all_selected_locus_universe_denominator':152,
            'repeated_mapping_locus_universe_denominator':140,
            'n_primary_retained':sum(bool(support[c]&primary) for c in ids),
            'n_EHR_supported':sum(bool(support[c]&groups['secondary_ehr']) for c in ids),
            'n_ML_linked':sum(ML in support[c] for c in ids),'n_ML_only':sum(support[c]=={ML} for c in ids),
            'v1_candidate_ids':joined(ids),'v1_priority_ranks':';'.join(map(str, sorted(int(frozenmap[c]['predicted_causal_priority_rank']) for c in ids))),
            **fields(ss),'boundary':'frozen tag-mapped gene OR coordinate-window assignment; not new target mapping or independent genetic locus'})
    deltas=[]
    for t in agg:
        if agg[t]-sig[t]: deltas.append({'tag_id':t,'v1_aggregate_studies':joined(agg[t]),'GWS_studies':joined(sig[t]),
            'extra_nonsignificant_studies':joined(agg[t]-sig[t]),'interpretation':'reported association below GWS significance for this study-tag pair'})
    bio=[]
    for group in ['broad','primary_direct','secondary_ehr','ml_surrogate']:
        eligible=groups[group]; tt=[r for r in tags if sig[r['normalized_variant_id']] & eligible]
        counts=Counter(r['consequence_class'] for r in tt)
        for label,n in sorted(counts.items()):
            bio.append({'stratum':group,'universe':'gws_tags','feature':f'consequence:{label}','n':n,'denominator':len(tt)})
        for name,ids in [('all_15389',candidate_ids),('frozen_337',set(frozenmap))]:
            cc=[c for c in ids if support[c]&eligible]
            for key in ['coding_CDS','donor_refined_any','known_regulatory_any']:
                bio.append({'stratum':group,'universe':name,'feature':key,'n':sum(yes(ann[c][key]) for c in cc),'denominator':len(cc)})
            if name=='frozen_337':
                for key in ['predicted_causal_enhancer','predicted_causal_silencer']:
                    bio.append({'stratum':group,'universe':name,'feature':key,'n':sum(yes(frozenmap[c][key]) for c in cc),'denominator':len(cc)})
    # Transparently bound all classification decisions: per-study additions/removals.
    influence=[]
    for s,a in sorted(studies.items()):
        alt=primary-{s} if s in primary else primary|{s}
        row={'study_accession':s,'phenotype_class':a['phenotype_class'],'primary_eligible':s in primary,
             'confidence':a['confidence'],'hypothetical_change':'exclude_from_primary' if s in primary else 'include_in_primary',
             'meaning':'classification influence only, not recommended reclassification'}
        for name,umap in universes.items():
            row[f'{name}_primary_change']=sum(bool(ss&alt) for ss in umap.values())-sum(bool(ss&primary) for ss in umap.values())
        influence.append(row)
    publication_influence=[]
    by_pub=defaultdict(set)
    for b in base: by_pub[b['PUBMED ID']].add(b['STUDY ACCESSION'])
    for pmid,accessions in sorted(by_pub.items()):
        alt=primary-accessions
        row={'pubmed_id':pmid,'study_accessions':joined(accessions),
             'primary_accessions':joined(accessions&primary),
             'interpretation':'leave-one-publication provenance sensitivity; publications may share cohorts'}
        for name,umap in universes.items():
            row[f'{name}_primary_lost']=sum(bool(ss&primary) and not bool(ss&alt) for ss in umap.values())
        publication_influence.append(row)
    check('Unmatched tags explicitly retained', sum(not yes(r['matched']) for r in focal)==85, 'All 85 frozen match failures retained')
    unmatched={r['normalized_variant_id'] for r in focal if not yes(r['matched'])}
    check('All 85 unmatched tags retain direct candidate representation',unmatched<=set().union(*(tok(r['gws_tag_ids']) for r in candidates)) and {r['tag_id'] for r in tagrows if not yes(r['v1_matched'])}==unmatched,'Exact ID sets, not only counts')
    check('Frozen shortlist identity/order preserved', [r['candidate_record_id'] for r in shortrows]==[r['candidate_record_id'] for r in shortlist], '12 unchanged shortlist records')
    check('No naive-panel over-retention in main result', not (naive_excess & independently_retained), f'{len(naive_excess)} all-record false retentions avoided; {len(naive_excess & set(frozenmap))} in frozen 337')
    save('R001_accession_adjudication.tsv',register)
    save('R002_cohort_overlap_provenance.tsv',cohort_rows)
    save('R002B_cohort_overlap_pairs.tsv',overlap,['study_accession_a','study_accession_b','shared_cohort_groups','interpretation'])
    save('R003_study_tag_ancestry_support.tsv',audit_rows)
    save('R004_all_candidate_phenotype_support.tsv.gz',allrows)
    save('R005_frozen_337_phenotype_support.tsv',focused)
    save('R006_component_phenotype_support.tsv',components)
    save('R007_shortlist_phenotype_support.tsv',shortrows)
    save('R008_stratum_summary.tsv',summary)
    save('R009_GCST90244098_dependence.tsv',mlrows)
    save('R010_gws_tag_phenotype_support.tsv',tagrows)
    save('R011_frozen_locus_support.tsv',locusrows)
    save('R012_aggregate_nonsignificant_provenance_audit.tsv',deltas,['tag_id','v1_aggregate_studies','GWS_studies','extra_nonsignificant_studies','interpretation'])
    save('R013_unresolved_tag_support.tsv',[r for r in tagrows if not yes(r['v1_matched'])])
    save('R014_frozen_biological_feature_summary.tsv',bio)
    save('R015_classification_influence.tsv',influence)
    save('R016_publication_influence.tsv',publication_influence)
    save('sources.tsv',source_rows)
    # Hash boundary covers every consumed V1 artifact and reviewed input.
    check('All inputs unchanged during execution', all(digest(ROOT/p)==h for p,h in INPUTS.items()), f'{len(INPUTS)} inputs')
    rectangular=True
    for output in OUTPUTS:
        output_path=ROOT/output['path']
        opener=gzip.open if output_path.suffix=='.gz' else open
        with opener(output_path,'rt',newline='') as handle:
            values=list(csv.reader(handle,delimiter='\t'))
        rectangular &= len(values)==output['rows']+1 and len(set(map(len,values)))==1
    check('All emitted TSVs rectangular and row counts verified', rectangular, f'{len(OUTPUTS)} output tables; read-back checked')
    save('validation_checks.tsv',CHECKS)
    manifest={'analysis_id':f'{PREFIX}-RUN-001','generated_utc':datetime.now(timezone.utc).isoformat(),
              'python':sys.version,'platform':platform.platform(),'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
              'command':'python3 diseases/COPD/07_gap_closure/scripts/run_phenotype_robustness.py',
              'script_sha256':digest(Path(__file__)), 'specification_sha256':digest(V2/'provenance'/f'{PREFIX}_analysis_specification.md'),
              'inputs':[{'path':p,'sha256':h} for p,h in sorted(INPUTS.items())], 'outputs':OUTPUTS,
              'gws_accessions':len({r['STUDY ACCESSION'] for r in gws}), 'classification_counts':dict(Counter(a['phenotype_class'] for a in studies.values())),
              'incorrect_tag_only_filter_excess_all':len(naive_excess),'incorrect_tag_only_filter_excess_337':len(naive_excess & set(frozenmap)),
              'elapsed_seconds':round(time.perf_counter()-started,3),
              'max_rss_KiB':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
              'v1_modified':False,'other_modules_executed':False}
    path=V2/'provenance'/f'{PREFIX}_run_manifest.json'
    path.write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n')
    log=json.dumps({'classification_counts':manifest['classification_counts'],'summary':[r for r in summary if r['stratum'] in ['broad','primary_direct','secondary_ehr','ml_surrogate','primary_or_ehr']], 'ML':mlrows,'checks':len(CHECKS),'elapsed_seconds':manifest['elapsed_seconds'],'max_rss_KiB':manifest['max_rss_KiB']},indent=2)
    (V2/'logs'/f'{PREFIX}_run_log.json').write_text(log+'\n')
    print(log)


if __name__ == '__main__':
    main()
