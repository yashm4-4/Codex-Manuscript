#!/usr/bin/env python3
"""Independent preflight artifact/provenance validator; no scientific inference.

Reads source metadata, bounded acquired bytes, reports and Git state. Historical
inputs are opaque-hashed only. Never reads model/candidate/benchmark contents.
Outputs validation receipts only; does not alter source tables or source bytes.
"""
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import subprocess
import zlib

STAGE = Path(__file__).resolve().parents[1]
ROOT = STAGE.parents[3]
CHECKS = []
KIM = {'GCST90016588', 'GCST90016589', 'GCST90016593', 'GCST90016594'}
BBJ = {'GCST90013709', 'GCST90013746', 'GCST90013781'}


def check(name, passed, detail=''):
    CHECKS.append({'check': name, 'passed': bool(passed), 'detail': str(detail)})


def readj(path):
    return json.loads((STAGE / path).read_text())


def table(path):
    with (STAGE / path).open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def git(*args):
    env = dict(os.environ, GIT_OPTIONAL_LOCKS='0')
    return subprocess.check_output(['git', *args], cwd=ROOT, env=env)


def local(path):
    p = Path(path)
    if p.is_absolute():
        return p
    q = STAGE / p
    return q if q.exists() else ROOT / p


def decode_prefix(data):
    if data[:2] != b'\x1f\x8b':
        return data.decode('utf-8-sig', errors='replace')
    out = bytearray()
    pending = data
    while pending and len(out) < 1048576:
        decoder = zlib.decompressobj(31)
        out.extend(decoder.decompress(pending, 1048576-len(out)))
        pending = decoder.unused_data if decoder.eof else b''
    return bytes(out).decode('utf-8-sig', errors='replace')


def main():
    started = datetime.now(timezone.utc).isoformat()
    init = readj('provenance/initialization.json')
    frozen = table('inputs/frozen_study_metadata.tsv')
    core = {r['study_accession']: r for r in frozen}
    check('frozen_metadata_unique104', len(frozen) == len(core) == 104)
    check('candidate_support_columns_removed', not any(k in frozen[0] for k in init['candidate_support_columns_excluded']))
    for binding in init['frozen_inputs']:
        p = ROOT / binding['path']
        check('historical_input_opaque_hash:'+binding['path'], p.stat().st_size == binding['bytes'] and digest(p) == binding['sha256'])
    matrix = table('tables/study_final_eligibility_matrix.tsv')
    catmatrix = table('tables/study_phenotype_summary_stat_eligibility.tsv')
    check('final_matrix_exact104_accessions', len(matrix) == 104 and {r['accession'] for r in matrix} == set(core))
    check('catalog_matrix_exact104_accessions', len(catmatrix) == 104 and {r['accession'] for r in catmatrix} == set(core))
    for row in matrix:
        a = row['accession']
        c = readj(f'audits/catalog/{a}/audit.json')
        m = c['current_catalog_metadata']
        check('frozen_class:'+a, row['frozen_phenotype_class'] == core[a]['phenotype_class'])
        check('publication_identity:'+a, str(m['pubmed_id']) == row['pubmed_id'] == core[a]['pubmed_id'])
        check('current_accession_identity:'+a, m['accession_id'] == a and c['catalog_api_access']['status'] == 200)
        check('no_execution_clearance:'+a, all(row[k] == 'False' for k in ['whole_file_integrity_density_QC_passed', 'locus_QC_passed', 'inference_executable_now']))
        check('no_fabricated_eligible_locus_count:'+a, row['eligible_locus_count'] == '' and 'NOT_DETERMINED' in row['eligible_locus_count_status'])
    check('direct_phenotype27', sum(r['frozen_phenotype_class'] == 'direct_copd_susceptibility' for r in matrix) == 27)
    direct_signed = table('tables/accessible_signed_direct_COPD_sources.tsv')
    check('direct_signed_exact7', len(direct_signed) == 7 and {r['accession'] for r in direct_signed} == KIM | BBJ)
    check('integrated_matrix_direct_signed_exact7', {r['accession'] for r in matrix if r['frozen_phenotype_class'] == 'direct_copd_susceptibility' and r['signed_variant_effect_schema_verified'] == 'True'} == KIM | BBJ)
    check('zero_immediate_execution_datasets', table('tables/immediately_executable_fine_mapping_datasets.tsv') == [])

    headers = table('tables/header_schema_audit.tsv')
    catmanifest = readj('provenance/catalog_table_manifest.json')
    check('catalog_header_count138', len(headers) == catmanifest['file_headers'] == 138)
    check('catalog_header_accessions68', len({r['accession'] for r in headers}) == catmanifest['accessions_with_inspected_headers'] == 68)
    for index, row in enumerate(headers):
        key = row['accession']+':'+str(index)
        p = STAGE / row['prefix_path']
        check('header_prefix_hash:'+key, digest(p) == row['prefix_sha256'])
        lines = decode_prefix(p.read_bytes()).splitlines()
        while lines and lines[0].startswith('##'):
            lines.pop(0)
        sep = '\t' if '\t' in lines[0] else ',' if ',' in lines[0] else None
        parsed_header = lines[0].lstrip('#').split(sep)
        parsed_records = [x.split(sep) for x in lines[1:9]]
        check('header_independent_reparse:'+key, parsed_header == json.loads(row['header']))
        check('representative_records_exact:'+key, parsed_records == json.loads(row['representative_records']))
        check('representative_widths:'+key, all(len(x) == len(parsed_header) for x in parsed_records))
        check('no_prefix_full_file_claim:'+key, row['whole_file_downloaded_and_validated'] == 'False')

    ledger = table('provenance/source_access_ledger.tsv')
    ledger_manifest = readj('provenance/source_access_manifest.json')
    check('source_ledger_count', len(ledger) == ledger_manifest['access_attempt_receipts'])
    check('source_ledger_binding', digest(STAGE/'provenance/source_access_ledger.tsv') == ledger_manifest['source_ledger_sha256'])
    unique = {r['saved_path'] for r in ledger}
    check('source_ledger_unique_paths', len(unique) == ledger_manifest['distinct_saved_response_files'])
    saved_total = 0
    for index, row in enumerate(ledger):
        p = local(row['saved_path'])
        n = int(row['saved_bytes'])
        saved_total += n
        check('source_bytes_and_hash:'+str(index), p.stat().st_size == n and digest(p) == row['saved_sha256'], row['saved_path'])
        check('source_receipt_exists:'+str(index), (STAGE/row['receipt_path']).is_file())
        receipt = readj(row['receipt_path'])
        receipt_hashes = [receipt.get(k) for k in ['sha256','body_sha256','acquired_sha256','saved_sha256','stored_sha256']]
        check('source_receipt_hash_binding:'+str(index), row['saved_sha256'] in receipt_hashes)
        check('source_access_has_url_time:'+str(index), bool(row['source_url']) and bool(row['accessed_utc']))
        if row['http_status'] == '206':
            check('partial_source_not_full_hash:'+str(index), row['bounded_or_partial'] == 'True' and 'partial' in row['hash_scope'])
    check('source_saved_bytes_total', saved_total == ledger_manifest['saved_response_bytes'])

    for a in sorted(KIM):
        c = readj(f'audits/catalog/{a}/audit.json')
        raw = next(h for h in c['headers'] if h['url'].endswith(a+'_buildGRCh37.tsv'))
        check('Kim_raw_explicit_fields:'+a, raw['header'] == ['variant_id','p_value','chromosome','base_pair_location','effect_allele','other_allele','odds_ratio','standard_error'])
        check('Kim_raw_access206:'+a, raw['access_status'] == 206)
        check('Kim_missing_AF_N_QC_not_invented:'+a, not any(f in raw['header'] for f in ['effect_allele_frequency','N','info','Rsq']))
    bbj = readj('audits/alternatives/bbj_archive_schema_evidence.json')
    check('BBJ_four_payload_headers', set(bbj['headers']) == {'autosome','autosome_sex_stratified','chrX','chrX_sex_stratified'})
    check('BBJ_archive_not_fully_validated', all(not h['whole_archive_or_member_CRC_verified'] for h in bbj['headers'].values()))
    alts = readj('audits/alternatives/alternatives_audit.json')
    check('provider_dataset_rows25', len(alts['datasets']) == len(table('tables/provider_dataset_ancestry_inventory.tsv')) == 25)
    expected = {'GCST90013709': ('BETA',3315,201592), 'GCST90013746': ('BETA.x',2855,103089), 'GCST90013781': ('BETA.y',460,98503)}
    for a, (effect, cases, controls) in expected.items():
        r = next(x for x in alts['datasets'] if x.get('accession') == a)
        check('BBJ_effect_contract:'+a, r['effect_allele'] == 'Allele2 (ALT)' and r['effect_field'] == effect and effect in r['actual_header'])
        check('BBJ_sample_release:'+a, r['cases'] == cases and r['controls'] == controls and r['build'] == 'GRCh37/hg19')
        check('BBJ_not_execution_cleared:'+a, r['fine_mapping_ready_now'] is False and r['actual_loci_enumerated'] is False)
    risk = {r['accession']:r for r in table('tables/signed_effect_risk_direction_feasibility.tsv')}
    limits = {r['accession']:r for r in table('tables/study_specific_limitations.tsv')}
    ancestry = {r['accession']:r for r in table('tables/ancestry_cohort_sample_size.tsv')}
    ld_requirements = {r['accession']:r for r in table('tables/study_LD_requirements.tsv')}
    for name, mapping in [('risk',risk),('limitations',limits),('ancestry',ancestry),('LD_requirements',ld_requirements)]:
        check('integrated_'+name+'_all104', set(mapping) == set(core))
    for a in sorted(KIM | BBJ):
        rr = risk[a]
        check('integrated_risk_signed_source:'+a, rr['risk_direction_feasibility'].startswith('SIGNED_SCHEMA_AVAILABLE') and rr['risk_alleles_assigned'] == 'False')
        check('integrated_risk_current_evidence:'+a, rr['evidence_scope'] == 'INTEGRATED_CATALOG_AND_VERIFIED_PROVIDER' and rr['source_url'].startswith('https://'))
        current_limit = limits[a]
        stale = ['public representative variant statistics unverified', 'signed effect/effect-allele mapping not demonstrated', 'nonmissing SE not demonstrated', 'nonmissing P not demonstrated']
        check('integrated_limitations_not_stale:'+a, current_limit['public_access_status'].startswith('PUBLIC_') and not any(s in current_limit['current_blockers'] for s in stale))
        check('provider_override_flag_exact_scope:'+a, current_limit['catalog_missing_source_access_resolved_by_provider'] == ('True' if a in BBJ else 'False'))
    for a, (effect,cases,controls) in expected.items():
        check('integrated_BBJ_risk_mapping:'+a, risk[a]['signed_effect_fields'] == effect and risk[a]['effect_allele_fields_and_semantics'] == 'Allele2 (ALT)')
        check('integrated_BBJ_sample_mapping:'+a, int(ancestry[a]['provider_cases_if_verified']) == cases and int(ancestry[a]['provider_controls_if_verified']) == controls)
    for a in ['GCST90691934','GCST90692407','GCST90692991']:
        check('integrated_PanUKB_negative_log_P_recognized:'+a, 'neg_log_10_p_value' in risk[a]['P_fields'])
    for a, row in ld_requirements.items():
        check('LD_table_no_matrix_or_QC_claim:'+a, row['matrix_values_acquired'] == row['locus_matrix_QC_executed'] == 'False')
    exports = readj('audits/direct/dbgap_public_export_schema_audit.json')
    check('dbGaP_incomplete_export_scope', [r['record_count'] for r in exports] == [19373,25412,25338,25350])
    for r in exports:
        check('dbGaP_absolute_not_signed:'+r['analysis_accession'], 'Absolute value' in r['effect_definition_verbatim'] and not r['usable_as_full_signed_locus_statistics'] and not r['whole_original_gwas_acquired'])

    recommendations = readj('provenance/recommendations.json')
    check('recommendation_counts', recommendations['core_accessions_audited'] == 104 and recommendations['direct_core_accessions'] == 27 and recommendations['direct_accessions_with_verified_public_signed_schema'] == 7)
    check('recommendation_unknown_loci', recommendations['expected_eligible_locus_count'] is None and recommendations['loci_execution_cleared'] == 0)
    check('recommendation_no_ready_primary', recommendations['best_primary_fully_ready_dataset'] is None and recommendations['immediately_executable_fine_mapping_dataset_count'] == 0)
    for flag in ['fine_mapping_executed','risk_alleles_assigned','candidate_membership_used','model_scoring_executed','complete_genomewide_GWAS_or_LD_matrix_values_downloaded','author_contact_sent','commit_or_push_performed']:
        check('declared_prohibition:'+flag, recommendations[flag] is False)
    for binding in readj('provenance/assembly_manifest.json')['bound_source_audits']:
        p = STAGE / binding['path']
        check('assembled_source_binding:'+binding['path'], p.stat().st_size == binding['bytes'] and digest(p) == binding['sha256'])

    # Historical content is only compared/hashes checked, never parsed for science.
    baseline = init['baseline_commit']
    check('HEAD_unchanged_no_commit', git('rev-parse','HEAD').decode().strip() == baseline)
    allowed = {r['path'] for r in init['baseline_registers']}
    stage_rel = str(STAGE.relative_to(ROOT))+'/'
    changed = git('diff','--name-only','--no-ext-diff',baseline,'--').decode().splitlines()
    check('historical_tracked_files_unchanged_except_registers', all(p in allowed or p.startswith(stage_rel) for p in changed), changed)
    check('register_update_receipt_exists', (STAGE/'provenance/shared_register_updates.json').is_file())
    register_receipt = readj('provenance/shared_register_updates.json')
    register_by = {r['path']: r for r in register_receipt['registers']}
    check('register_receipt_no_prior_changes', register_receipt['prior_bytes_changed'] is False and set(register_by) == allowed)
    for row in init['baseline_registers']:
        p = ROOT / row['path']
        old = (STAGE/'provenance'/('baseline_'+p.name)).read_bytes()
        current = p.read_bytes()
        check('register_baseline_hash:'+p.name, hashlib.sha256(old).hexdigest() == row['sha256'])
        check('register_append_only:'+p.name, current.startswith(old) and len(current) > len(old))
        check('register_new_rows_stage_only:'+p.name, all(b'fine-mapping-preflight-1.0' in line for line in current[len(old):].splitlines() if line.strip()))
        rr = register_by[row['path']]
        check('register_receipt_current_binding:'+p.name, len(current) == rr['bytes'] and hashlib.sha256(current).hexdigest() == rr['sha256'])
        check('register_receipt_appended_counts:'+p.name, len(current)-len(old) == rr['appended_bytes'] and len(current[len(old):].splitlines()) == rr['appended_rows'])
    report = STAGE / 'FINE_MAPPING_PREFLIGHT_REPORT.md'
    check('final_report_exists_nonempty', report.is_file() and report.stat().st_size > 5000)
    report_text = report.read_text()
    check('report_source_response_count_matches', f'{len(ledger):,} saved responses' in report_text)
    check('report_source_bytes_match', f'{saved_total:,} bytes' in report_text)
    for term in ['GCST007692','GCST90016588','GCST90013709','SuSiE','FINEMAP','Pan','104','27']:
        check('report_required_anchor:'+term, term in report_text)
    check('no_scientific_array_or_checkpoint_created', not any(p.suffix in {'.keras','.npy','.npz','.bed','.bim','.fam'} for p in STAGE.rglob('*') if p.is_file()))

    outputs = ['FINE_MAPPING_PREFLIGHT_REPORT.md','provenance/recommendations.json','provenance/assembly_manifest.json','tables/study_final_eligibility_matrix.tsv','tables/accessible_signed_direct_COPD_sources.tsv','tables/signed_effect_risk_direction_feasibility.tsv','tables/study_specific_limitations.tsv','tables/ancestry_cohort_sample_size.tsv','tables/study_LD_requirements.tsv','provenance/source_access_ledger.tsv','provenance/shared_register_updates.json']
    bindings = [{'path': p,'bytes':(STAGE/p).stat().st_size,'sha256':digest(STAGE/p)} for p in outputs]
    first_attempt = {
        'provenance/validation_attempt_01/validate_preflight.py':'7ec74a00f3aa72e7e662310c892934953169e650e462ac33b1b1213e1df4a8ef',
        'provenance/validation_attempt_01/independent_validation.json':'7067031852d2ad52c470beb81938c3bdac0e71c251fb18271cb29dcc4cb5c868',
        'provenance/validation_attempt_01/independent_validation.tsv':'4582dd3ef412024122261268609e265fcf53088fdd0bb55e9d824fb71d59f7e8',
    }
    for path, expected_sha in first_attempt.items():
        check('preserved_first_attempt_exact:'+path, digest(STAGE/path) == expected_sha)
    failed = [r for r in CHECKS if not r['passed']]
    result = {
        'stage':'fine-mapping-preflight-1.0','validator':'Independent metadata/schema/provenance audit',
        'started_utc':started,'completed_utc':datetime.now(timezone.utc).isoformat(),
        'status':'PASS' if not failed else 'FAIL','checks':len(CHECKS),'failed_checks':len(failed),'failures':failed,
        'source_receipts':len(ledger),'saved_source_bytes':saved_total,'catalog_headers':len(headers),
        'core_accessions':104,'direct_accessions':27,'verified_signed_direct_sources':7,
        'whole_file_and_LD_scientific_QC_not_performed':True,'fine_mapping_not_executed_by_validator':True,
        'candidate_and_model_contents_not_read_by_validator':True,
        'historical_change_check_scope':'Git tracked content comparison to baseline and opaque frozen-input hashes; historical data not parsed. Append-only register prefixes checked.',
        'prohibition_audit_limit':'File/receipt/static-contract and Git-state checks corroborate declared scope; this is not a system-wide process-forensics audit.',
        'bound_artifacts':bindings,'future_freeze_not_self_referenced':True,
        'prior_attempt': {'files_and_sha256':first_attempt,'status':'FAIL','checks':4805,'failures':4,
                         'resolution':'Prefreeze schema/test semantics mismatch: first validator incorrectly required a provider-override flag for all seven accessible direct sources; the flag applied only to three BBJ releases missing from Catalog. Root renamed it catalog_missing_source_access_resolved_by_provider. Revised checks still require public access and no stale schema blockers for all seven, and verify exact override-flag scope. No source bytes or scientific analysis changed.'},
        'substantive_review': 'Independent reviewer read the report and risk contract against source audits before validation: no execution-ready claim, unsigned dbGaP exports rejected, Kim joint-test locus counts not substituted, BBJ allele/sex schema separated, Pan-UKB remains secondary, eligible loci unknown, candidate/model/benchmark selection firewall explicit.',
    }
    destination = STAGE/'provenance/independent_validation.json'
    with destination.open('x') as handle:
        json.dump(result,handle,indent=2,sort_keys=True);handle.write('\n')
    with (STAGE/'provenance/independent_validation.tsv').open('x') as handle:
        writer=csv.DictWriter(handle,['check','passed','detail'],delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(CHECKS)
    print(json.dumps(result,indent=2))
    raise SystemExit(0 if not failed else 1)


if __name__ == '__main__':
    main()
