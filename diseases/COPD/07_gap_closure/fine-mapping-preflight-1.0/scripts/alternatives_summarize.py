#!/usr/bin/env python3
"""Assemble provider access/schema evidence; never enumerate GWAS loci or run analysis."""
import csv
import datetime
import gzip
import hashlib
import json
from pathlib import Path
import re
import struct
import tarfile
import xml.etree.ElementTree as ET
import zipfile
import zlib

STAGE = Path(__file__).resolve().parents[1]
OUT = STAGE / 'audits/alternatives'
D = OUT / 'downloads'


def save(name, obj):
    with (OUT / name).open('w') as f:
        json.dump(obj, f, indent=2, sort_keys=True)
        f.write('\n')


def xlsx_rows(label, sheet):
    with zipfile.ZipFile(D / (label + '.bin')) as z:
        assert z.testzip() is None, 'Complete workbook CRC failed'
        ns = {'m': 'http://schemas.openxmlformats.org/spreadsheetml/2006/main'}
        ss = [''.join(t.itertext()) for t in ET.fromstring(z.read('xl/sharedStrings.xml')).findall('m:si', ns)]
        rows = []
        for row in ET.fromstring(z.read(sheet)).findall('.//m:row', ns):
            values = []
            for c in row.findall('m:c', ns):
                v = c.find('m:v', ns)
                val = v.text if v is not None else ''
                values.append(ss[int(val)] if c.get('t') == 's' and val else val)
            if any(values):
                rows.append(values)
        return rows


def zip_prefix(label):
    """Decode only first non-directory member of a bounded ZIP range, not whole-file QC."""
    b = (D / (label + '.bin')).read_bytes()
    offset = 0
    while b[offset:offset+4] == b'PK\x03\x04':
        h = struct.unpack_from('<IHHHHHIIIHH', b, offset)
        nlen, xlen = h[-2:]
        name = b[offset+30:offset+30+nlen].decode()
        start = offset + 30 + nlen + xlen
        if name.endswith('/') and h[7] == 0:
            offset = start
            continue
        data = b[start:start+h[7]]
        assert h[3] in (0, 8)
        if h[3] == 8:
            data = zlib.decompressobj(-15).decompress(data, 1048576)
        assert data.startswith(b'\x1f\x8b')
        text = zlib.decompressobj(31).decompress(data, 1048576).decode()
        lines = text.splitlines()
        return {'archive_member': name, 'zip_method': h[3], 'compressed_member_bytes': h[7],
                'gzip_member_bytes': h[8], 'header': lines[0].split(),
                'representative_records': [line.split() for line in lines[1:6]],
                'whole_archive_or_member_CRC_verified': False,
                'inspection_scope': 'header and five records only; no locus or effect analysis'}
    raise ValueError('No local ZIP member in range')


def main():
    receipts = {}
    for file in sorted(D.glob('*.request.json')):
        rec = json.loads(file.read_text())
        label = file.name.removesuffix('.request.json')
        body = D / (label + '.bin')
        assert body.stat().st_size == rec['stored_bytes']
        assert hashlib.sha256(body.read_bytes()).hexdigest() == rec['stored_sha256']
        rec['receipt_path'] = str(file.relative_to(STAGE))
        rec['hash_scope'] = 'saved response bytes only; never a whole remote-file hash for partial responses'
        if label == 'gbmi_manifest_xlsx':
            with zipfile.ZipFile(body) as z:
                assert z.testzip() is None
            rec['audit_completeness_override'] = 'Complete ZIP central directory and all member CRCs pass; original old-helper receipt marked prefix conservatively because response lacked Content-Length.'
        receipts[label] = rec

    def evidence(*labels):
        return [{'label': label, 'url': receipts[label]['requested_url'],
                 'path': receipts[label]['stored_path'], 'receipt': receipts[label]['receipt_path'],
                 'status': receipts[label].get('http_status'), 'bytes': receipts[label]['stored_bytes'],
                 'sha256': receipts[label]['stored_sha256']} for label in labels]

    def record(dataset_id, **kwargs):
        return dict(dataset_id=dataset_id, actual_loci_enumerated=False,
                    published_independent_COPD_locus_count_verified=None,
                    eligible_locus_count='UNKNOWN: future whole-GWAS QC and candidate-blind locus definition required',
                    fine_mapping_ready_now=False, **kwargs)

    rows = []
    with gzip.open(D / 'pan_manifest.bin', 'rt') as f:
        pan = [r for r in csv.DictReader(f, delimiter='\t')
               if r['phenocode'] in ('J44', '496', '496.21') and r['pheno_sex'] == 'both_sexes']
    save('pan_selected_manifest_records.json', pan)
    for r in pan:
        label = {'J44': 'pan_J44_header', '496': 'pan_phecode496_header', '496.21': 'pan_phecode49621_header'}[r['phenocode']]
        for pop in r['pops'].split(','):
            rows.append(record('PanUKB_' + r['phenocode'] + '_' + pop,
                provider='Pan-UK Biobank', phenotype=r['description'], phenotype_role='secondary_EHR_ICD_or_PheCode',
                phenotype_caveat=('PheCode496 also includes J47 bronchiectasis, emphysema and chronic bronchitis; not interchangeable with J44.' if r['phenocode'] == '496' else
                                  'Obstructive chronic bronchitis is a distinct phenotype stratum.' if r['phenocode'] == '496.21' else
                                  'ICD10 J44 endpoint; not a spirometrically adjudicated COPD primary.'),
                ancestry=pop, ancestry_scope='ancestry-specific columns, not cross-ancestry meta columns',
                cases=int(r['n_cases_' + pop]), controls=int(r['n_controls_' + pop]),
                N_source='complete provider phenotype manifest; case/control counts, not effective N',
                build='GRCh37', variant_coding='chr,pos,ref,alt; positive-strand REF/ALT',
                effect_allele='ALT', effect_field='beta_' + pop, se_field='se_' + pop,
                p_field='neglog10_pval_' + pop, af_fields=['af_cases_' + pop, 'af_controls_' + pop],
                phenotype_qc=r['phenotype_qc_' + pop],
                qc_note='not_EUR_plus_1 is the cross-ancestry eligibility criterion; it alone does not disqualify an EUR-only GWAS. low_confidence_POP and variant-manifest INFO/quality still require filtering.',
                stats_access='anonymous HTTP206 actual BGZF header/records verified',
                stats_scope='provider complete per-phenotype release, not GWS-only; full payload not downloaded or row-count/density checked',
                whole_file_bytes=int(r['size_in_bytes']), provider_whole_file_md5=r['md5_hex'],
                whole_file_hash_locally_verified=False, actual_header=receipts[label]['header_and_first_records'][0].split('\t'),
                LD_candidate='public same-ancestry PanUKB covariate-adjusted ALT-dosage LD',
                LD_blockers=['case/control/sample and variant intersection contract must be established', 'residualization/scaling and diagonal QC pending', 'signed-R extraction and per-locus consistency not run'],
                overlap='UK Biobank overlaps other UKB COPD phenotypes, Sakornsakolpat discovery and GBMI UKBB component; not independent replication.',
                recommendation='J44 EUR is the strongest technically promising secondary backup; do not promote EHR to direct-COPD primary.' if r['phenocode'] == 'J44' and pop == 'EUR' else 'Separate secondary sensitivity only; preserve phenotype and ancestry boundaries.',
                evidence=evidence('pan_manifest', 'pan_schema_doc', 'pan_qc_doc', 'pan_phecode_map', label)))

    gbmi_schema = xlsx_rows('gbmi_manifest_xlsx', 'xl/worksheets/sheet1.xml')
    gbmi_files = [r for r in xlsx_rows('gbmi_manifest_xlsx', 'xl/worksheets/sheet2.xml') if any('COPD_Bothsex' in v for v in r)]
    save('gbmi_workbook_extracted_metadata.json', {'schema': gbmi_schema, 'COPD_bothsex_rows': gbmi_files, 'complete_archive_CRC_check': 'PASS'})
    accession_map = {'afr': 'GCST90399691', 'amr': 'GCST90399692', 'eas': 'GCST90399693', 'eur': 'GCST90399694', 'mixed': 'GCST90399695'}
    for pop, accession in accession_map.items():
        label = 'gbmi_' + pop + '_header'
        r = next(r for r in gbmi_files if r[4] == 'all biobanks' and (r[3] == pop or pop == 'mixed' and ',' in r[3]))
        catalog = json.loads((STAGE / 'audits/catalog' / accession / 'audit.json').read_text())
        rows.append(record('GBMI_COPD_' + pop, provider='GBMI May2021 release', accession=accession,
            phenotype='Harmonized biobank COPD endpoint', phenotype_role='secondary_biobank_EHR_expanded', ancestry=r[3],
            ancestry_scope='multi-ancestry combined' if pop == 'mixed' else 'ancestry-specific multi-biobank meta-analysis',
            participating_biobanks=r[5].split(','), N_source=catalog['current_catalog_metadata']['initial_sample_size'],
            per_variant_N_fields=['N_case', 'N_ctrl', 'n_dataset', 'n_bbk'],
            build='GRCh38', variant_coding='#CHR POS REF ALT', effect_allele='ALT', effect_field='inv_var_meta_beta',
            se_field='inv_var_meta_sebeta', p_field='inv_var_meta_p', af_field='all_meta_AF',
            stats_access='anonymous HTTP206 actual gzip header/records verified',
            stats_scope='genome-wide release restricted to >=2 contributing biobanks; not GWS-only; full file QC not performed',
            actual_header=receipts[label]['header_and_first_records'][0].split('\t'),
            LD_candidate='No matching meta-analysis signed-R resource established',
            LD_blockers=['single EUR LD prohibited for mixed-ancestry meta-analysis', 'within-ancestry biobank/sample weights and variant missingness must be matched', 'GBMI SLALOM gnomAD weighted r-squared diagnostics are not a signed-R replacement'],
            overlap='Biobank components overlap PanUKB, BBJ, FinnGen and CKB as explicitly listed; leave-UKBB-out mixed file is still multi-ancestry.',
            recommendation='Secondary ancestry-specific sensitivity only after explicit LD/missingness contract; mixed meta-analysis not approved for single-reference fine-mapping.',
            evidence=evidence('gbmi_manifest_xlsx', label), catalog_evidence='audits/catalog/' + accession + '/audit.json'))

    bbj_headers = {'chrX': zip_prefix('bbj_copd_zip_prefix'), 'autosome': zip_prefix('bbj_copd_autosome_prefix'),
                   'autosome_sex_stratified': zip_prefix('bbj_copd_autosome_sex_prefix'),
                   'chrX_sex_stratified': zip_prefix('bbj_copd_chrx_sex_prefix')}
    b = (D / 'bbj_copd_zip_tail.bin').read_bytes()
    offset = b.find(b'PK\x01\x02')
    members = []
    while offset >= 0 and b[offset:offset+4] == b'PK\x01\x02':
        h = struct.unpack_from('<I6H3I5H2I', b, offset)
        nlen, xlen, clen = h[10:13]
        members.append({'member': b[offset+46:offset+46+nlen].decode(), 'local_header_offset': h[-1],
                        'compressed_bytes': h[8], 'uncompressed_gzip_bytes': h[9]})
        offset += 46+nlen+xlen+clen
    bbj_samples = [r for r in xlsx_rows('bbj_sample_size', 'xl/worksheets/sheet1.xml') if 'COPD' in r]
    save('bbj_archive_schema_evidence.json', {'headers': bbj_headers, 'central_directory_members': members,
        'sample_sheet_COPD_row': bbj_samples, 'archive_total_bytes': 2244311450,
        'whole_archive_CRC_or_SHA_not_verified': True, 'archive_discovery': 'COPD accession hum0014.v17.COPD.v1 in official sample-size sheet, combined with distribution filename pattern; URL verified by HTTP206.'})
    for sex, acc, cases, controls in [('both', 'GCST90013709', 3315, 201592), ('male', 'GCST90013746', 2855, 103089), ('female', 'GCST90013781', 460, 98503)]:
        rows.append(record('BBJ_Ishigaki2020_COPD_' + sex, provider='NBDC hum0014 original release', accession=acc,
            phenotype='Physician-diagnosed COPD; disease indexed as J449', phenotype_role='primary_eligible_direct_clinical_COPD',
            ancestry='Japanese', ancestry_scope='single ancestry', sex=sex, cases=cases, controls=controls,
            build='GRCh37/hg19', variant_coding='CHR POS SNPID Allele1=REF Allele2=ALT', effect_allele='Allele2 (ALT)',
            effect_field={'both': 'BETA', 'male': 'BETA.x', 'female': 'BETA.y'}[sex],
            se_field={'both': 'SE', 'male': 'SE.x', 'female': 'SE.y'}[sex],
            p_field={'both': 'p.value', 'male': 'p.value.x', 'female': 'p.value.y'}[sex],
            sample_size_contract='Allsex autosomal N column verified. Sex-stratified payloads have no N field; source XLSX supplies sex-specific case/control N, requiring later missingness checks.',
            quality='SAIGE0.29.4.2; 1000Gphase3v5 imputation; Rsq>=0.7 and MAC>=10 per source docs',
            stats_access='Original COPD ZIP publicly reachable by HTTP206. JENGER HTTP/HTTPS timed out; Catalog mirrored file absent; OpenGWAS VCF login/JWT required.',
            stats_scope='ZIP central directory and all four payload headers verified: allsex autosomal, allsex chrX, sex-stratified autosomal, sex-stratified chrX; five representative records each. Male.x and female.y share the sex-stratified file and are not independent source files.',
            actual_header=bbj_headers['autosome' if sex == 'both' else 'autosome_sex_stratified']['header'],
            whole_file_bytes=2244311450, whole_file_hash_locally_verified=False,
            LD_candidate='Japanese cohort-compatible dense signed LD not publicly verified',
            LD_blockers=['public fine-mapping outputs do not establish access to original dense dosage LD', 'generic EAS or 1000G JPT is not assumed interchangeable with BBJ LD'],
            overlap='Bothsex contains male/female analyses; all overlap BBJ component in expanded2021/GBMI. Sex results are sensitivity, not independent replication.',
            recommendation='Promising direct-COPD primary phenotype with verified signed source access; blocked on matched dense LD and subsequent whole-file/per-locus QC, not on summary-statistic availability.',
            evidence=evidence('bbj_readme', 'bbj_sample_size', 'bbj_copd_zip_prefix', 'bbj_copd_zip_tail', 'bbj_copd_autosome_prefix', 'bbj_copd_autosome_sex_prefix', 'bbj_copd_chrx_sex_prefix', 'jenger_index', 'jenger_http', 'bbj_opengwas_metadata'),
            primary_methods_evidence='audits/direct/acquired/ishigaki2020_bioc/response.bin'))

    with (D / 'finngen_r11_manifest_correct.bin').open() as f:
        fg = next(r for r in csv.DictReader(f, delimiter='\t') if r['phenocode'] == 'J10_COPD')
    save('finngen_selected_manifest_record.json', fg)
    rows.append(record('FinnGen_R11_J10_COPD', provider='FinnGen R11', phenotype='COPD J10_COPD',
        phenotype_role='secondary_registry_ICD', ancestry='Finnish', ancestry_scope='Finnish founder population; not interchangeable with generic EUR',
        cases=int(fg['num_cases']), controls=int(fg['num_controls']), build='GRCh38',
        variant_coding='#chrom pos ref alt', effect_allele='ALT', effect_field='beta (log OR, REGENIE)',
        se_field='sebeta', p_field='pval', af_fields=['af_alt', 'af_alt_cases', 'af_alt_controls'],
        stats_access='anonymous public GCS HTTP206 verified; official documentation separately asks users to register for download instructions',
        stats_scope='full endpoint summary-statistic release advertised and non-GWS prefix verified; fullfile/density QC not run',
        actual_header=receipts['finngen_r11_header']['header_and_first_records'][0].split('\t'),
        whole_file_bytes=815727196, whole_file_hash_locally_verified=False,
        LD_candidate='public FinnGen LD API based on SISu reference, not verified dense R11 association-sample matrix',
        LD_blockers=['release-matched Finnish signed-R source and index not established', 'endpoint-specific definition page TLS failure; general registry derivation verified; exact endpoint code list still needs freezing'],
        overlap='FinnGen appears in GBMI and GEMINI; repeated releases are not independent replication.',
        recommendation='Separate Finnish registry sensitivity, subject to matched dense LD; no generic-EUR substitution.',
        evidence=evidence('finngen_r11_manifest_correct', 'finngen_r11_header', 'finngen_r11_schema', 'finngen_r11_endpoint_docs', 'finngen_r11_docs_markdown', 'finngen_r11_endpoint')))

    for pop, acc in [('AFR', 'GCST90476026'), ('EUR', 'GCST90476027'), ('AMR', 'GCST90478138'), ('EAS', 'GCST90482079'), ('mixed', 'GCST90480254')]:
        a = json.loads((STAGE / 'audits/catalog' / acc / 'audit.json').read_text())
        headers = a.get('headers', [])
        rows.append(record('MVP_PheCode496_' + pop, provider='MVP 2024 GIA gwPheWAS', accession=acc,
            phenotype='Chronic airway obstruction PheCode496', phenotype_role='secondary_EHR_PheCode',
            ancestry=pop, ancestry_scope='multi-ancestry meta-analysis' if pop == 'mixed' else 'ancestry-specific',
            N_source=a['current_catalog_metadata']['initial_sample_size'], build='GRCh38 (Catalog metadata)',
            effect_allele='explicit effect_allele; do not substitute stale raw alt column after harmonization',
            effect_field='odds_ratio (signed log effect recoverable as log(OR)); standard_error sample values missing',
            se_caveat='Raw and harmonized EUR preview standard_error is NA despite populated header. CI-derived SE would require source CI-level/rounding verification, Wald-vs-SPA consistency and prespecified validation; not certified ready.',
            stats_access='Catalog full-file prefix verified where headers exist; provider dbGaP route separately requires application.' if headers else 'No Catalog payload verified; official policy withholds ancestry-specific groups below500cases (EAS PheCode496 has423cases).',
            actual_header=headers[0]['header'] if headers else None,
            stats_scope='Genome-wide source release, not hits-only; bounded prefix checks only where available.',
            LD_candidate='MVP ancestry-specific dense signed LD access not established',
            LD_blockers=['dense ancestry/sample-matched LD not verified', 'missing SE and rounded OR/CI semantics must be resolved', 'GIA vs HARE releases must not be mixed'],
            overlap='MVP-derived COPD studies can share participants; no independent-replication claim. Chronic bronchitis496.21 is a separate endpoint.',
            recommendation='Secondary sensitivity only after SE and LD blockers; mixed cannot use one ancestry reference.',
            evidence=evidence('mvp_access_doc', 'mvp_dbgap_instructions'), catalog_evidence='audits/catalog/' + acc + '/audit.json'))

    acc = 'GCST90473707'
    a = json.loads((STAGE / 'audits/catalog' / acc / 'audit.json').read_text())
    rows.append(record('UKB_WGS_J44_GCST90473707', provider='UKB WGS via Catalog', accession=acc,
        phenotype='ICD10 J44 other chronic obstructive pulmonary disease', phenotype_role='secondary_EHR_ICD',
        ancestry=a['current_catalog_metadata']['discovery_ancestry'], ancestry_scope='as Catalog accession metadata',
        N_source=a['current_catalog_metadata']['initial_sample_size'], build='see raw versus harmonized Catalog YAML; no build mixing',
        effect_allele='explicit effect_allele', effect_field='beta', se_field='standard_error',
        stats_access='Catalog HTTP206 raw and harmonized actualheaders/records verified',
        actual_header=a['headers'][0]['header'], stats_scope='genome-wide single-variant additive file; not gene-burden file, no wholefile QC',
        LD_candidate='PanUKB imputed-array LD is not an automatic match for UKB WGS rare variants',
        LD_blockers=['WGS-specific sample/build/rare-variant density and dense signed LD contract unverified', 'older PanUKB index may exclude rare WGS alleles'],
        overlap='Same UK Biobank population, overlapping PanUKB and GBMI.',
        recommendation='Separate secondary possibility; no automatic preference for latest/largest over matched LD.',
        catalog_evidence='audits/catalog/' + acc + '/audit.json', evidence=[]))

    with tarfile.open(D / 'gemini_readmes.bin', 'r:gz') as t:
        member = next(m for m in t.getmembers() if m.name == 'GEMINI.sumstats.COPD.README')
        gemini_readme = t.extractfile(member).read().decode()
    gemini_condition = [line for line in (D / 'gemini_conditions.bin').read_text().splitlines() if '\tCOPD\t' in line]
    save('gemini_COPD_source_metadata.json', {'README': gemini_readme, 'condition_header': (D / 'gemini_conditions.bin').read_text().splitlines()[0], 'condition_row': gemini_condition})
    rows.append(record('GEMINI_v1_COPD', provider='GEMINI Zenodo14284047', phenotype='UKB+FinnGen COPD meta-analysis',
        phenotype_role='secondary_mixed_registry_EHR', ancestry='European UKB plus Finnish', ancestry_scope='two-cohort meta-analysis',
        N_source='COPD-specific source README effective N=161709; do not confuse with actual total sample count; per-variant n_samples changes',
        build='GRCh37', variant_coding='chromosome base_pair_location effect_allele other_allele',
        effect_allele='explicit effect_allele', effect_field='beta', se_field='standard_error',
        stats_access='HTTP206 actualgzipheader/records verified', stats_scope='genome-wide GWAMA release, non-GWS records verified; no wholefileQC',
        actual_header=receipts['gemini_copd_prefix']['header_and_first_records'][0].split('\t'), whole_file_bytes=491567826,
        LD_candidate='No matching cohort-weighted signed-R resource verified',
        LD_blockers=['heterogeneous UKB+FinnGen contribution varies by variant', 'single UKB or generic EUR reference not automatically matched'],
        overlap='Explicit UKB+FinnGen; neither independent of PanUKB nor FinnGen.',
        recommendation='Newly located accessible secondary alternative; not direct-COPD consortium data despite general repository mentioning consortia.',
        evidence=evidence('gemini_zenodo_metadata', 'gemini_readmes', 'gemini_conditions', 'gemini_copd_prefix')))
    rows.append(record('CKB_GCST90246122_provider_access', provider='CKB PheWeb', accession='GCST90246122',
        phenotype_role='secondary_ICD_biobank', ancestry='Chinese', ancestry_scope='single biobank',
        stats_access='Landing/about HTTP200. Full downloads encrypted; keys require application. No keys requested or restrictions bypassed.',
        build='Browser GRCh38; Walters2023 downloadable statistics remain GRCh37 according to /about',
        effect_allele='not verified from plaintext payload in this audit', actual_header=None,
        license='CC BY4.0 display policy plus separate key-sharing/decrypted-data redistribution restrictions; must honor both',
        LD_candidate='not established', LD_blockers=['lawful plaintext summary-statistic access and matched LD not established'],
        recommendation='Access-contingent secondary source; landing-page public-access claim is not evidence of unrestricted decrypted bulk files.',
        evidence=evidence('ckb_pheweb_landing', 'ckb_pheweb_about'), catalog_evidence='audits/catalog/GCST90246122/audit.json'))

    report = {'stage': 'fine-mapping-preflight-1.0', 'created_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'scope': 'provider access/phenotype/schema audit only; no candidate/model/benchmark access, no fine-mapping, no locus enumeration',
        'datasets': rows, 'download_receipts': list(receipts.values()),
        'bounded_download_bytes': sum(r['stored_bytes'] for r in receipts.values()),
        'readiness_definition': 'Source access/schema evidence is not proof of whole-file density, locus eligibility, signed-LD consistency or fine-mapping readiness.',
        'signed_effect_rule': 'Positive beta/log(OR) increases endpoint liability for the documented effect allele; negative decreases. P values alone never assign direction. ALT only where provider explicitly says so; harmonized effect_allele controls harmonized rows.',
        'primary_recommendation': 'Ishigaki2020 Japanese direct COPD has verified original signed summary-statistic access but no verified matched dense LD. Root direct-anchor audit also evaluates European direct COPD; do not replace primary with EHR by convenience.',
        'secondary_recommendation': 'PanUKB J44 EUR is a promising same-project signed-LD backup, subject to case/control/scaling/per-locus gates. Others remain sensitivity or access-contingent.',
        'ML_rule': 'No machine-learning COPD phenotype is pooled or promoted here; such frozen104 accessions remain their separate class in root matrix.',
        'forbidden_operations_performed': []}
    save('alternatives_audit.json', report)
    print(json.dumps({'datasets': len(rows), 'retained_responses': len(receipts), 'bytes': report['bounded_download_bytes'], 'receipt_hash_checks': 'PASS'}))


if __name__ == '__main__':
    main()
