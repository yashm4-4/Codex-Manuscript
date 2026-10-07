"""Assemble central inventories from existing completed audits, without rereading GWAS.

Only small receipts/audits are read or hashed. Large source sizes are stat-checked;
their original acquisition SHA-256 values are copied, not newly verified here.
The final stage freeze performs independent full-payload hash verification.
"""
from pathlib import Path
import csv
import datetime
import hashlib
import json

S = Path(__file__).resolve().parents[1]
T = S / 'tables'
T.mkdir(exist_ok=True)
INPUTS = {}
CHECKS = []
NA = 'NOT_APPLICABLE'


def small(rel):
    p = S / rel
    assert p.stat().st_size < 50_000_000, 'Inventory must not scan large source files'
    b = p.read_bytes()
    INPUTS[rel] = {'bytes': len(b), 'sha256': hashlib.sha256(b).hexdigest()}
    return b


def read_json(rel):
    return json.loads(small(rel))


def read_tsv(rel):
    return list(csv.DictReader(small(rel).decode().splitlines(), delimiter='\t'))


def digest(rel):
    if rel not in INPUTS:
        small(rel)
    return INPUTS[rel]['sha256']


def refs(paths):
    return {p: digest(p) for p in paths}


def check(name, passed):
    CHECKS.append({'check': name, 'pass': bool(passed)})
    assert passed, name


def cell(v):
    if isinstance(v, (dict, list)):
        return json.dumps(v, sort_keys=True, separators=(',', ':'))
    if isinstance(v, bool):
        return str(v).lower()
    return v


def write_tsv(name, records):
    cols = list(dict.fromkeys(k for row in records for k in row))
    with (T / name).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t', lineterminator='\n')
        w.writeheader()
        for row in records:
            w.writerow({k: cell(row.get(k, 'NOT_RECORDED')) for k in cols})


acquisition = []
integrity = []
inventory = []


def add_source(*, track, source_id, role, payload, receipt, rec, nrows, build,
               compression, audit, schema, missingness, qc_paths, scope,
               provider_digest_type, provider_digest, provider_digest_evidence,
               provider_digest_status, crc_status, complete_read_status,
               byte_count=None, sha256=None, md5=None, note='', extra_receipts=(),
               headers=None, started=None, completed=None, archive_member=NA):
    headers = headers or rec.get('response_headers', rec.get('headers', {}))
    byte_count = byte_count if byte_count is not None else rec.get('bytes', rec.get('actual_bytes'))
    sha256 = sha256 or rec['sha256']
    md5 = md5 if md5 is not None else rec.get('md5', 'NOT_RECORDED')
    check(source_id + '_payload_exists', (S / payload).is_file())
    check(source_id + '_recorded_size_matches_current_stat', (S / payload).stat().st_size == int(byte_count))
    check(source_id + '_recorded_sha256_well_formed', len(sha256) == 64 and all(x in '0123456789abcdef' for x in sha256))
    receipts = refs([receipt, *extra_receipts])
    a = {
        'track': track, 'source_id': source_id, 'role': role,
        'payload_path': payload, 'archive_member': archive_member,
        'url': rec.get('url', 'NOT_RECORDED'), 'final_url': rec.get('final_url', rec.get('url', 'NOT_RECORDED')),
        'access_started_utc': started or rec.get('accessed_utc', rec.get('started_utc', 'NOT_RECORDED')),
        'access_completed_utc': completed or rec.get('completed_utc', rec.get('finished_utc', 'NOT_RECORDED')),
        'http_status': rec.get('status', rec.get('status_code', 'SEE_RANGE_RECEIPTS')),
        'last_modified': headers.get('Last-Modified', 'NOT_RECORDED'),
        'etag': headers.get('ETag', headers.get('Etag', 'NOT_RECORDED')),
        'object_version': headers.get('x-amz-version-id', 'NOT_RECORDED'),
        'bytes': byte_count, 'recorded_sha256': sha256, 'recorded_md5': md5,
        'provider_digest_type': provider_digest_type, 'provider_digest': provider_digest,
        'provider_digest_status': provider_digest_status,
        'provider_digest_evidence': refs(provider_digest_evidence) if provider_digest_evidence else {},
        'data_rows': nrows, 'genome_build': build, 'compression': compression,
        'scientific_scope': scope, 'receipt_path': receipt,
        'receipt_sha256': digest(receipt), 'all_receipt_paths_sha256': receipts,
        'large_payload_hash_recomputed_by_this_script': False,
        'recorded_hash_basis': 'Existing acquisition/full-file audit; final stage freeze verifies payload',
        'notes': note,
    }
    acquisition.append(a)
    integrity.append({k: a[k] for k in ['track', 'source_id', 'role', 'payload_path', 'archive_member',
        'bytes', 'recorded_sha256', 'recorded_md5', 'provider_digest_type', 'provider_digest',
        'provider_digest_status', 'data_rows', 'genome_build', 'compression', 'scientific_scope',
        'receipt_path', 'receipt_sha256', 'large_payload_hash_recomputed_by_this_script']} | {
        'current_size_matches_recorded': True, 'compression_crc_status': crc_status,
        'complete_read_status': complete_read_status, 'audit_path': audit,
        'audit_sha256': digest(audit), 'notes': note})
    inventory.append({
        'track': track, 'source_id': source_id, 'role': role, 'payload_path': payload,
        'data_rows': nrows, 'genome_build': build, 'scientific_scope': scope,
        'whole_file_audit_path': audit, 'whole_file_audit_sha256': digest(audit),
        'schema_fields': schema, 'schema_evidence_path': audit, 'schema_evidence_sha256': digest(audit),
        'missingness_evidence': refs(missingness), 'QC_and_effect_evidence': refs(qc_paths),
        'QC_interpretation': note,
    })


# Track A primary native file plus two explicitly contextual Catalog transports.
an = read_json('tracks/A/results/native_audit.json')
af = read_json('tracks/A/results/formatted_audit.json')
ah = read_json('tracks/A/results/harmonized_audit.json')
afinal = read_json('tracks/A/audit.json')
for i, r in enumerate(read_tsv('tracks/A/results/full_file_integrity.tsv')):
    rel = 'tracks/A/' + r['path']
    receipt = rel + '.access.json'
    rec = read_json(receipt)
    if i == 0:
        sid, role, audit, doc, schema = 'A_NATIVE', 'PRIMARY_GWAS_INPUT', 'tracks/A/results/native_audit.json', an, an['source_columns']
        note = 'Native is authoritative GRCh37. Mandatory statistic/allele fields complete; 518903 missing rsIDs. AF/N/INFO/explicit REF/ALT absent. 25089 nonsignificant frozen-precision failures retained; no repair. Provider metadata SNP count 9886854 differs from data rows 9886853.'
        qc = [audit, 'tracks/A/results/failure_characterization.json', 'tracks/A/config/effect_contract.json', 'tracks/A/results/reference_audit.json', 'tracks/A/results/artifact_validation.json']
        provider_path = 'tracks/A/sources/md5sum.txt'
    elif i == 1:
        sid, role, audit, doc = 'A_FORMATTED37', 'CONTEXTUAL_FORMATTED_TRANSPORT_NOT_REPLACEMENT', 'tracks/A/results/formatted_audit.json', af
        schema = list(af['missingness'])
        note = 'Context only: complete formatted GRCh37 file agrees with native row identities and OR/SE/P strings after missing-rsID quote handling. No native precision failures replaced or repaired.'
        qc = [audit, 'tracks/A/config/formatted_transport_contract.json']
        provider_path = 'tracks/A/sources/harmonised/md5sum.txt'
    else:
        sid, role, audit, doc = 'A_HARMONIZED38', 'CONTEXTUAL_CATALOG_HARMONIZATION_NOT_ANALYSIS_BUILD', 'tracks/A/results/harmonized_audit.json', ah
        schema = ah['columns']
        note = 'GRCh38 auxiliary only: 9594025 oriented allele/OR transforms pass; 94994 unoriented. 9344885 exact native-transport matches; 344134 unmatched exact keys remain unproven; no rsID-only analysis join.'
        qc = [audit, 'tracks/A/config/effect_contract.json']
        provider_path = 'tracks/A/sources/harmonised/md5sum.txt'
    entries = {x.split()[-1].lstrip('*'): x.split()[0] for x in small(provider_path).decode().splitlines() if x.strip()}
    expected = entries[Path(rel).name]
    check(sid + '_provider_MD5_receipt_match', expected == r['md5'] == rec['md5'] and r['provider_md5_pass'] == 'True')
    check(sid + '_audit_row_count_consistency', int(r['data_rows']) == doc.get('rows', doc.get('formatted_rows')))
    add_source(track='A', source_id=sid, role=role, payload=rel, receipt=receipt, rec=rec,
        nrows=int(r['data_rows']), build=r['genome_build'], compression=r['compression'],
        audit=audit, schema=schema, missingness=[audit], qc_paths=qc,
        scope='Ever-smoker spirometric COPD; ' + ('primary native analysis' if i == 0 else 'contextual transport/allele audit'),
        provider_digest_type='PROVIDER_PUBLISHED_MD5', provider_digest=expected,
        provider_digest_evidence=[provider_path, provider_path + '.access.json'], provider_digest_status='PASS_RECORDED_EXACT_MATCH',
        crc_status='NOT_COMPRESSED' if i == 0 else 'PASS_RECORDED_FULL_GZIP_DECODE_CRC',
        complete_read_status='PASS_RECORDED_FULL_FILE_AUDIT', note=note)


# Track B: preserve original archive and distinguish its analyzed member.
br = read_json('tracks/B/acquisition/receipt.json')
bi = read_json('tracks/B/acquisition/initial_incomplete_receipt.json')
bn = read_json('tracks/B/results/whole_file_audit.json')
bfinal = read_json('tracks/B/results/audit.json')
bscore = read_json('tracks/B/results/score_scale_sample_audit.json')
check('B_archive_recorded_size_and_CRC_pass', br['full_size_pass'] and br['all_zip_members_crc32_pass'])
common_b = dict(track='B', rec=br, receipt='tracks/B/acquisition/receipt.json',
    provider_digest_type='NO_INDEPENDENT_PROVIDER_DIGEST_RECORDED', provider_digest='NOT_AVAILABLE',
    provider_digest_evidence=[], provider_digest_status='NOT_ASSESSED_NO_PROVIDER_DIGEST; recorded local SHA256 and container CRC available',
    headers=bi['headers'], started=bi['started_utc'],
    extra_receipts=['tracks/B/acquisition/initial_incomplete_receipt.json', 'tracks/B/acquisition/range_resume_receipts.json'])
add_source(**common_b, source_id='B_ARCHIVE', role='ORIGINAL_DISTRIBUTION_CONTAINER',
    payload='tracks/B/acquisition/hum0014.v17.COPD.v1.zip', nrows=NA, build='hg19/GRCh37 members', compression='ZIP containing GZIP members',
    audit='tracks/B/acquisition/receipt.json', schema='Archive member inventory only; not a GWAS table', missingness=[],
    qc_paths=['tracks/B/results/complete_file_integrity.tsv'], scope='All members container CRC only; only combined-sex autosomal COPD analyzed',
    crc_status='PASS_RECORDED_ALL_ZIP_MEMBERS_CRC32', complete_read_status='PASS_RECORDED_CONTAINER_AND_SIZE_CHECK; other members not scientifically parsed',
    note='Initial HTTP200 prefix was incomplete and is documented; completed by 18 validated HTTP206 byte ranges. Multipart ETag is not a full-file MD5. Archive retains chrX/sex-specific members without scientific analysis.')
add_source(**common_b, source_id='B_COMBINED_AUTOSOMAL', role='PRIMARY_GWAS_INPUT',
    payload='tracks/B/acquisition/COPD.auto.rsq07.mac10.txt.gz', nrows=bn['rows'], build='hg19/GRCh37', compression='GZIP; space-delimited text',
    audit='tracks/B/results/whole_file_audit.json', schema=bn['schema'],
    missingness=['tracks/B/results/whole_file_schema_missingness.tsv', 'tracks/B/results/whole_file_audit.json'],
    qc_paths=['tracks/B/results/whole_file_qc.tsv', 'tracks/B/results/effect_test_consistency.tsv', 'tracks/B/results/duplicate_conflicts.json', 'tracks/B/results/score_scale_sample_audit.json', 'tracks/B/results/implementation_correction.json', 'tracks/B/config/tolerances.json'],
    scope='Original BBJ combined-sex autosomal COPD only', byte_count=br['analyzed_member']['bytes'], sha256=br['analyzed_member']['sha256'],
    archive_member=br['analyzed_member']['name'], crc_status='PASS_RECORDED_ZIP_MEMBER_CRC_AND_FULL_GZIP_DECODE_CRC',
    complete_read_status='PASS_RECORDED_FULL_FILE_AUDIT',
    note='All mandatory fields present; 4681 symbolic/non-ACGT rows; 3 excess conflicting duplicate identities outside selected loci. N204905-204907; 11542 nominal-N differences are diagnostic, not failures under frozen positive-N rule. Preserved implementation correction restores 5 loci/462 significant rows. No independent archive/member provider MD5 recorded.')


# Track C: analyze EUR statistics only; provider variant manifest is contextual.
cr = read_json('tracks/C/sources/acquisition.json')
cmr = read_json('tracks/C/sources/variant_manifest.acquisition.json')
cn = read_json('tracks/C/results/whole_file_audit.json')
cm = read_json('tracks/C/results/variant_manifest_audit.json')
cs = read_json('tracks/C/results/source_contract.json')
cfinal = read_json('tracks/C/results/audit.json')
manifest = read_json('tracks/C/sources/selected_manifest.json')[0]
check('C_GWAS_provider_manifest_MD5_receipt_match', manifest['md5_hex'] == cr['expected_md5'] == cr['md5'])
check('C_GWAS_recorded_source_hashes_agree', cn['source_sha256'] == cr['sha256'] and cn['source_bytes'] == cr['actual_bytes'])
check('C_GWAS_recorded_full_integrity_pass', cr['integrity_pass'] and cn['gzip_complete_decode_crc_pass'] and cn['row_count_matches_provider'])
add_source(track='C', source_id='C_J44_EUR', role='PRIMARY_GWAS_INPUT',
    payload='tracks/C/' + cr['path'], receipt='tracks/C/sources/acquisition.json', rec=cr,
    nrows=cn['rows'], build=cs['build'], compression='BGZF/gzip complete multiblock stream',
    audit='tracks/C/results/whole_file_audit.json', schema=cn['header'],
    missingness=['tracks/C/results/whole_file_schema_missingness.tsv', 'tracks/C/results/whole_file_audit.json'],
    qc_paths=['tracks/C/results/numeric_range_audit.json', 'tracks/C/results/identity_boundary_verification.json', 'tracks/C/results/whole_file_invalid_numeric_rows.tsv', 'tracks/C/results/source_contract.json', 'tracks/C/results/statistic_semantics_addendum.json'],
    scope='Complete source union file integrity/schema; only EUR statistical columns analyzed; other ancestry/meta columns missingness only',
    provider_digest_type='PROVIDER_PHENOTYPE_MANIFEST_MD5', provider_digest=manifest['md5_hex'],
    provider_digest_evidence=['tracks/C/sources/selected_manifest.json', 'tracks/C/sources/phenotype_manifest_correct.tsv.access.json'],
    provider_digest_status='PASS_RECORDED_EXACT_MATCH', crc_status='PASS_RECORDED_FULL_GZIP_DECODE_CRC', complete_read_status='PASS_RECORDED_FULL_FILE_AUDIT',
    note='28987534 union rows, not all valid EUR. 23861813 valid signed EUR; 18170886 source-confidence-passing autosomal rows. 5125720 EUR missing; one nonmissing invalid SE0/infinite-logP row excluded without repair. Multipart object ETag is not provider MD5; expected MD5 comes from phenotype manifest.')
check('C_manifest_single_part_ETag_MD5_receipt_match', cmr['md5'] == cmr['expected_md5_from_single_part_s3_etag'] == cmr['headers']['ETag'].strip('"'))
check('C_manifest_recorded_full_integrity_pass', cmr['integrity_pass'] and cm['gzip_complete_decode_crc_pass'] and cm['row_count_matches_provider'])
add_source(track='C', source_id='C_VARIANT_MANIFEST', role='CONTEXTUAL_VARIANT_QC_DEPENDENCY',
    payload='tracks/C/data/full_variant_qc_metrics.txt.bgz', receipt='tracks/C/sources/variant_manifest.acquisition.json', rec=cmr,
    nrows=cm['rows'], build=cs['build'], compression='BGZF/gzip complete multiblock stream',
    audit='tracks/C/results/variant_manifest_audit.json', schema=cm['parsed_columns'],
    missingness=['tracks/C/results/variant_manifest_audit.json'], qc_paths=['tracks/C/results/variant_manifest_audit.json', 'tracks/C/results/source_contract.json'],
    scope='Identity/rsID/INFO/QC and EUR frequency/allele-number fields only; gene columns not parsed or used',
    provider_digest_type='SINGLE_PART_S3_ETAG_COMPARED_AS_MD5_NOT_SEPARATE_CHECKSUM_SIDECAR', provider_digest=cmr['expected_md5_from_single_part_s3_etag'],
    provider_digest_evidence=['tracks/C/sources/variant_manifest.head.json', 'tracks/C/sources/variant_manifest.acquisition.json'],
    provider_digest_status='PASS_RECORDED_LOCAL_MD5_EQUALS_SINGLE_PART_ETAG', crc_status='PASS_RECORDED_FULL_GZIP_DECODE_CRC', complete_read_status='PASS_RECORDED_FULL_FILE_AUDIT_OF_SELECTED_COLUMNS',
    note='Contextual manifest, not GWAS statistics. Selected schema only; gene annotations unused. AN/2 describes broader cohort and is not GWAS per-variant N. 397 INFO values print as0.8 despite documented upstream>0.8 filter; flagged, not repaired. Provider high_quality includes ancestry/gnomAD criteria and is not automatic GWAS exclusion.')


effects = [
    dict(track='A', source_id='A_NATIVE', signed_effect_contract='Source effect_allele is coded OR allele; source other_allele is comparison; neither implies REF/ALT',
         beta='natural_log(source odds_ratio)', se='supplied standard_error on log-odds scale; supported by source methods/PLINK2 schema', z='beta / supplied standard_error',
         p='source p_value; two-sided normal Wald consistency evaluated under printed-precision intervals',
         full_file_tested_rows=an['rows'], arithmetic_pass_rows=an['rounding_consistency_pass'], arithmetic_fail_rows=an['rounding_consistency_fail'],
         significant_seed_rows=an['significant_autosomal_rows'], significant_seed_arithmetic_fail_rows=0,
         interpretation='Signed source contract supported; 25089 nonsignificant precision-gate failures quarantined; all 1623 significant rows pass. Sample/N/identity/LD and final method readiness separate.',
         af_n_info='Native AF, per-variant N and INFO absent; nominal12446cases+59145controls=71591 is not observed/effective per-variant N',
         sample_or_lineage='Kim ever-smoker subset; full Pan-UKB EUR LD is external reference, not matched in-sample smoking-stratum LD',
         source_evidence=refs(['tracks/A/config/effect_contract.json', 'tracks/A/results/native_audit.json', 'tracks/A/results/failure_characterization.json', 'tracks/A/sources/plink2_formats.html.access.json', 'tracks/A/sources/kim2021_pmc.html.access.json', 'tracks/A/audit.json'])),
    dict(track='B', source_id='B_COMBINED_AUTOSOMAL', signed_effect_contract='Allele1=REF; Allele2=ALT and signed effect allele under provider README',
         beta='supplied BETA approximate signed ALT log-OR', se='supplied SAIGE SPA-calibrated SE; contemporaneous code constructs abs(logOR/qnorm(SPA_P/2)); exact score branch lineage unresolved',
         z='source BETA / source SE; no P-derived replacement', p='source p.value SPA; p.value.NA separately records unadjusted normal approximation',
         full_file_tested_rows=bn['p_consistency']['p.value']['n'], arithmetic_pass_rows=bn['p_consistency']['p.value']['n']-bn['p_consistency']['p.value']['fail'], arithmetic_fail_rows=bn['p_consistency']['p.value']['fail'],
         significant_seed_rows=bn['significant_valid_rows'], significant_seed_arithmetic_fail_rows=0,
         interpretation='Exact BETA/SE-to-SPA-P arithmetic is built into source SE construction, not independent Wald/RSS likelihood validation; normal-P discrepancies399244. Source-derived score branch reconstruction fails338/9141; no repairs.',
         af_n_info='AF_Allele2,AF.Cases,AF.Controls,N,Rsq,MAC supplied; N204905-204907; original positive-N gate preserved',
         sample_or_lineage='SAIGE0.29.4.2 provider version; score/software branch remains unresolved; no defensible accessible dense signed Japanese LD source obtained',
         source_evidence=refs(['tracks/B/results/whole_file_audit.json', 'tracks/B/results/effect_test_consistency.tsv', 'tracks/B/results/score_scale_sample_audit.json', 'tracks/B/config/tolerances.json', 'tracks/B/sources/bbj_readme.json', 'tracks/B/sources/SAIGE_0.29.4.2_archive.json', 'tracks/B/results/audit.json'])),
    dict(track='C', source_id='C_J44_EUR', signed_effect_contract='Native forward GRCh37 ALT effect allele, REF comparison allele; EUR only',
         beta='supplied beta_EUR signed ALT log-OR', se='supplied se_EUR; source SAIGE constructs abs(logOR/qnorm(SPA_P/2)); not independently a raw information/Wald SE',
         z='source beta_EUR / source se_EUR; no P-derived replacement', p='source neglog10_pval_EUR is -log10(SPA P); original scale preserved',
         full_file_tested_rows=cn['counts']['valid_signed_rows'], arithmetic_pass_rows=cn['counts']['wald_rounding_compatible_valid_rows'], arithmetic_fail_rows=cn['counts']['wald_rounding_discrepant_valid_rows'],
         significant_seed_rows=cn['counts']['genomewide_significant_seed_rows'], significant_seed_arithmetic_fail_rows=0,
         interpretation='Printed-precision arithmetic consistent for valid EUR rows; source SE derives from SPA P, so consistency is not independent evidence of ordinary Gaussian/Wald or RSS covariance suitability.',
         af_n_info='Case/control EUR AF supplied; per-variant GWAS N absent. PhenotypeN420531; LDglobalsN420542; contextual manifestAN/2 broader and not substituted.',
         sample_or_lineage='Exact J44 EUR current Hail metadata SAIGE0.36.3; current metadata container2024 versus GWAS flatfile2023 release relationship unresolved; same-project LD not established exact phenotype sample',
         source_evidence=refs(['tracks/C/results/whole_file_audit.json', 'tracks/C/results/source_contract.json', 'tracks/C/results/statistic_semantics_addendum.json', 'tracks/C/config/gwas_contract.json', 'tracks/C/sources/saige_0.36.3_DESCRIPTION.json', 'tracks/C/sources/saige_0.36.3_R_SAIGE_SPATest.R.json', 'tracks/C/results/audit.json'])),
]

check('exactly_seven_distinct_payloads', len(acquisition) == 7 and len({r['payload_path'] for r in acquisition}) == 7)
check('exactly_one_primary_GWAS_per_track', all(sum(r['track'] == t and r['role'] == 'PRIMARY_GWAS_INPUT' for r in acquisition) == 1 for t in 'ABC'))
check('A_primary_rows_agree_final_audit', an['rows'] == afinal['full_gwas_rows'] == 9886853)
check('B_primary_rows_agree_final_audit', bn['rows'] == bfinal['complete_GWAS_rows'] == 8678470)
check('C_primary_rows_agree_final_audit', cn['rows'] == cfinal['complete_gwas_rows'] == 28987534)
check('B_corrected_loci_and_significant_rows_used', bn['prospective_loci'] == bfinal['prospective_loci'] == 5 and bn['significant_valid_rows'] == 462)
check('B_score_branch_limitation_retained', bscore['failing_tolerance_rows'] == 338 and bscore['sample_rows'] == 9141)
check('all_inventory_provenance_paths_exist', all((S / p).is_file() for p in INPUTS))
check('all_large_hashes_reused_not_recomputed', all(not a['large_payload_hash_recomputed_by_this_script'] for a in acquisition))
for name, rows in [('full_file_acquisition.tsv', acquisition), ('complete_file_integrity.tsv', integrity),
                   ('whole_file_audit_inventory.tsv', inventory), ('effect_scale_summary.tsv', effects)]:
    write_tsv(name, rows)

outputs = {name: {'rows': len(rows), 'sha256': hashlib.sha256((T/name).read_bytes()).hexdigest()} for name, rows in
    [('full_file_acquisition.tsv', acquisition), ('complete_file_integrity.tsv', integrity),
     ('whole_file_audit_inventory.tsv', inventory), ('effect_scale_summary.tsv', effects)]}
result = {
    'status': 'PASS_INVENTORY_ASSEMBLY_AND_RECORDED_AUDIT_CONSISTENCY',
    'completed_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'scope': 'Existing receipt/audit assembly only; no source science rerun, no GWAS download, no inference, no scientific clearance assigned',
    'large_file_hash_verification': 'Not repeated here; recorded acquisition SHA256 copied with present-file size checks. Final freeze must independently hash full payloads.',
    'provider_digest_caveats': ['B has no independent provider full-file digest recorded; local SHA256/range/size/CRC evidence is distinct.',
                               'C manifest expected MD5 comes from single-part S3 ETag, not a separate checksum publication.',
                               'Multipart ETags for C GWAS and B archive are not treated as MD5.'],
    'checks': CHECKS, 'input_small_artifacts': INPUTS, 'outputs': outputs,
    'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
}
(T/'gwas_audit_inventory_validation.json').write_text(json.dumps(result, indent=2, sort_keys=True)+'\n')
print(json.dumps({'status': result['status'], 'checks_passed': len(CHECKS), 'output_rows': {k:v['rows'] for k,v in outputs.items()}}, indent=2))
