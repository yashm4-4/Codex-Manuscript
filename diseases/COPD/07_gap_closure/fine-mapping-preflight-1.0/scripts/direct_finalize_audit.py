#!/usr/bin/env python3
"""Assemble candidate-blind direct-COPD source/access/schema preflight evidence.

No association statistics, LD, fine mapping, or model scores are computed.
The only full GWAS-related files read are four small, already acquired dbGaP
public display exports; counting their rows/field presence audits release scope.
"""
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
from pathlib import Path

STAGE = Path(__file__).resolve().parents[1]
OUT = STAGE / 'audits/direct'
ACQUIRED = OUT / 'acquired'


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(name, value):
    path = OUT / name
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write('\n')


def evidence(label):
    p = ACQUIRED / label / 'access.json'
    item = json.loads(p.read_text())
    return {'label': label, 'access_record': str(p.relative_to(STAGE)),
            'url': item['url'], 'body_path': item['body_path'],
            'acquired_sha256': item['acquired_sha256'],
            'http_status': item['http_status']}


def dbgap_audit(pha):
    label = f'dbgap_{pha}_complete_small_export'
    p = ACQUIRED / label / 'response.bin'
    text = p.read_text()
    data = [x for x in text.splitlines() if x and not x.startswith('#')]
    rows = list(csv.DictReader(data, delimiter='\t'))
    beta = next(k for k in rows[0] if 'beta' in k)
    meta = json.loads((p.parent / 'access.json').read_text())
    headers = meta['response_headers']
    content_length = next(v for k, v in headers.items() if k.lower() == 'content-length')
    assert meta['http_status'] == 200 and not meta['read_truncated']
    assert len(p.read_bytes()) == int(content_length)
    return {
        'analysis_accession': pha,
        'source': evidence(label),
        'whole_public_export_acquired': True,
        'whole_original_gwas_acquired': False,
        'body_size_equals_http_content_length': True,
        'public_export_bytes': len(p.read_bytes()),
        'public_export_sha256': sha(p),
        'header': data[0].split('\t'),
        'record_count': len(rows),
        'effect_column': beta,
        'effect_definition_verbatim': next(x for x in text.splitlines() if x.startswith('# '+beta+':')),
        'negative_effect_strings_present': any(row[beta].startswith('-') for row in rows),
        'missing_coordinate_records': sum(not row['Chr Position'] for row in rows),
        'missing_sample_size_records': sum(not row['Sample size'] for row in rows),
        'first_representative_record': rows[0],
        'first_effect_allele_in_displayed_allele_pair': rows[0]['Coded Allele'] in [rows[0]['Allele1'], rows[0]['Allele2']],
        'metadata_lines': [x for x in text.splitlines() if x.startswith('#') and
                           any(k in x for k in ['Name:', 'Description:', 'Method:', 'genome build:', 'dbSNP build:'])],
        'usable_as_full_signed_locus_statistics': False,
        'reason': 'Small public result/display export, not complete tested-variant release; effect dictionary explicitly absolute, so valid signed effects unavailable. No direction is reconstructed from coded/minor/display alleles.',
        'inspection_scope': 'Release/schema audit only: row counts, lexical sign presence, field missingness and first record. No P-value ranking, locus definition, effect estimation or statistical testing performed.'
    }


def main():
    metadata_path = STAGE / 'inputs/frozen_study_metadata.tsv'
    rows = list(csv.DictReader(metadata_path.open(), delimiter='\t'))
    direct = [r for r in rows if r['phenotype_class'] == 'direct_copd_susceptibility']
    assert len(rows) == 104 and len(direct) == 27
    exports = [dbgap_audit(p) for p in ['pha004766', 'pha004496', 'pha004497', 'pha004498']]
    write_json('dbgap_public_export_schema_audit.json', exports)

    sources = []
    for r in direct:
        a = r['study_accession']
        catalog_path = STAGE / 'audits/catalog' / a / 'audit.json'
        catalog = json.loads(catalog_path.read_text())
        current = catalog.get('current_catalog_metadata', {})
        row = {
            'study_accession': a, 'publication': r['publication'], 'pubmed_id': r['pubmed_id'],
            'doi': r['doi'], 'phenotype_class': r['phenotype_class'],
            'phenotype_definition': {'case': r['case_definition'], 'control': r['control_definition'],
                                     'spirometry': r['spirometry'], 'ascertainment': r['ascertainment'],
                                     'smoking_design': r['smoking_design'], 'sex_design': r['sex_design']},
            'cohorts': r['cohorts'], 'overlap_considerations': r['overlap_notes'],
            'frozen_initial_N_description': r['initial_sample_description'],
            'current_catalog_initial_N_description': current.get('initial_sample_size'),
            'current_catalog_replication_N_description': current.get('replication_sample_size'),
            'current_catalog_discovery_ancestry': current.get('discovery_ancestry'),
            'current_catalog_replication_ancestry': current.get('replication_ancestry'),
            'ancestry_meta_analysis_status': 'See source-specific override or frozen cohort description; aggregate Catalog labels alone do not verify pooled-analysis ancestry.',
            'catalog_api_and_ftp_evidence': {'audit_path': str(catalog_path.relative_to(STAGE)),
                                          'audit_sha256': sha(catalog_path),
                                          'catalog_flag_not_usability_proof': current.get('full_summary_stats_available'),
                                          'directory_url': catalog.get('directory_url'),
                                          'directory_http_status': catalog.get('directory_access', {}).get('status')},
            'source_references_frozen': r['source_urls'].split(';'),
            'new_primary_evidence': [],
            'genome_build': 'Not verified from an actual summary-file header',
            'signed_effect_fields': {'variant_id': 'unverified', 'effect_allele': 'unverified',
                                     'other_allele': 'unverified', 'effect': 'unverified', 'SE': 'unverified',
                                     'P': 'unverified', 'frequency': 'unverified', 'per_variant_N': 'unverified',
                                     'variant_quality': 'unverified'},
            'public_access_status': 'No actual complete signed file verified in direct-source audit; consult current Catalog audit and source-specific alternative audit.',
            'unrestricted_download_verified': False,
            'whole_file_validated': False,
            'locus_completeness_status': 'Not established; Catalog association hits must not be used to reconstruct locus statistics.',
            'risk_direction_status': 'Not executable without signed file and verified effect-allele semantics.',
            'published_locus_count': None,
            'expected_eligible_locus_count': None,
            'LD_requirement': 'Study/ancestry-compatible signed LD with independently verified allele coding, variant coverage and sample provenance; none established by this direct-source artifact.',
            'access_or_execution_blockers': ['Full signed locus-complete statistics not verified', 'Matched signed LD and future locus-level QC not established'],
            'review_scope': 'Frozen phenotype adjudication plus current accession API/FTP checks; no additional full-text source-specific review claimed for this row.',
            'immediately_executable': False,
        }
        if a == 'GCST007692':
            row.update(
                review_scope='Independent full primary article, consortium announcement, dbGaP directories and complete small public export reviewed.',
                new_primary_evidence=[evidence(x) for x in ['sakornsakolpat2019_bioc', 'sakornsakolpat2019_nature', 'consortium_sakorn_announcement', 'dbgap_pha004766_complete_small_export']],
                ancestry_meta_analysis_status='Main fixed-effects UK Biobank + ICGC meta-analysis includes European and non-European cohorts. dbGaP methods list CHS African/European, MESA African/European/Hispanic, KARE and TCGS Korea. Article explicitly says predominantly European, not European-only. Per-variant ancestry weights unavailable.',
                primary_publication_N={'cases': 35735, 'controls': 222076, 'cohort_studies': 25, 'UKB_cases': 21081, 'UKB_controls': 179711, 'ICGC_cases': 14654, 'ICGC_controls': 42365},
                genome_build='Submitted association methods hg19/GRCh37; public dbGaP display-export coordinates labelled build38. Cannot mix these representations.',
                public_access_status='Public dbGaP pha004766 display export retrieved completely; Catalog FTP 404. No unrestricted complete signed meta-analysis file obtained.',
                unrestricted_download_verified=True,
                unrestricted_download_scope='Only incomplete, unsigned public display export; NOT original full summary statistics.',
                locus_completeness_status='19,373 public export records versus 6,224,355 tested variants reported by paper/current Catalog. Not locus-complete.',
                signed_effect_fields={'variant_id': 'SNP ID rs identifier; Submitted SNP ID is placeholder ss in first record', 'effect_allele': 'Coded Allele exists but first record T is not in displayed C/G pair; cannot use without source resolution', 'other_allele': 'Allele1/Allele2 are displayed dbSNP alleles, not a validated association-oriented pair', 'effect': '&beta; dictionary explicitly Absolute value of regression coefficient; all export effect strings nonnegative', 'SE': 'SE present', 'P': 'P-value present', 'frequency': 'No frequency column; Minor allele is an allele label, not AF', 'per_variant_N': 'Sample size column entirely blank in public export', 'variant_quality': 'pHWE and Call Rate columns exist but no original locus-complete QC release verified'},
                risk_direction_status='FAIL for public export: absolute coefficient and inconsistent effect-allele example. Neither minor allele nor coded allele recovers missing sign.',
                published_locus_count={'all_discovery_loci': 82, 'previously_known': 47, 'novel': 35, 'definition': 'Paper used ±1Mb windows; not prospective eligible-locus count.'},
                LD_requirement='Paper used 10,000 unrelated UKB participants for GCTA-COJO LD; that does not establish public signed LD, or adequate matching to the multi-ancestry combined meta-analysis. Request ancestry-specific stats and matched LD or justify a meta-analysis-aware alternative.',
                access_or_execution_blockers=['Request original complete signed meta-analysis statistics and preferably ancestry-specific results from authors/ICGC or authorized dbGaP/UKB route', 'Clarify submission build and effect-allele coding rather than use dbGaP display fields', 'Obtain defensible study/ancestry-matched signed LD and per-variant sample contribution', 'No pseudo-fine-mapping from public display export or Catalog hits'],
                data_availability_statement='Article names dbGaP phs000179.v5.p2 and UK Biobank; consortium announcement links pha4766 under phs000179.v6.p2. Public parent analyses directory is accessible, but its export is not the required full signed file.'
            )
        elif a == 'GCST004147':
            row.update(
                review_scope='Independent full primary article and three complete small dbGaP public exports reviewed.',
                new_primary_evidence=[evidence(x) for x in ['hobbs2017_bioc', 'hobbs2017_nature', 'dbgap_pha004496_complete_small_export', 'dbgap_pha004497_complete_small_export', 'dbgap_pha004498_complete_small_export']],
                ancestry_meta_analysis_status='Stage1 fixed-effects meta-analysis of22 genome-wide +4 custom-content cohorts explicitly included non-European populations; European-only subgroup additionally reported. Not inferred from Catalog labels.',
                primary_publication_N={'stage1_cases': 15256, 'stage1_controls': 47936, 'stage2_cases': 9498, 'stage2_controls': 9748, 'note': 'Stage sizes from paper abstract; Catalog discovery/replication descriptors differ. Do not silently interchange or sum them.'},
                genome_build='Original upload methods hg19 with +strand effect/other alleles; public dbGaP display exports labelled build38 with some coordinates absent.',
                public_access_status='Public pha004496/497/498 exports retrieved; original full signed Stage1 file not obtained; Catalog FTP404.',
                unrestricted_download_verified=True,
                unrestricted_download_scope='Only three incomplete unsigned public display exports, not full signed Stage1 statistics.',
                locus_completeness_status='25,412/25,338/25,350 export rows. Whole-genome Stage1 cannot be reconstructed from these selected display exports. Stage2 only selected top Stage1 associations (P<5e-6), not a dense final genome-wide meta-analysis.',
                signed_effect_fields={'variant_id': 'SNP ID and submitted ss identifiers', 'effect_allele': 'Coded Allele present; not sufficient to recover missing sign', 'other_allele': 'Allele1/Allele2 present but original association alignment not verified', 'effect': '|&beta;| explicitly absolute coefficient', 'SE': 'present', 'P': 'P-value present', 'frequency': 'No AF field', 'per_variant_N': 'Sample size blank throughout exports', 'variant_quality': 'Original QC documented, not reconstructed from public export'},
                risk_direction_status='FAIL for public exports: absolute coefficients. Original signed Stage1 data must be requested.',
                published_locus_count={'stage1_genome_wide_loci': 13, 'combined_after_selected_replication_loci': 22, 'note': 'Do not assign22 as dense Stage1/fine-mapping-eligible count.'},
                LD_requirement='Paper used COPDGene non-Hispanic white reference; not proof of public LD or matching to mixed-ancestry Stage1. Prefer ancestry-specific signed statistics + corresponding LD.',
                access_or_execution_blockers=['Obtain original complete signed Stage1 file with ancestry/cohort scope, excluding inference of completeness from selected Stage2 results', 'Obtain ancestry-appropriate signed LD', 'Resolve build and effect allele without display-export assumptions']
            )
        elif r['pubmed_id'] == '33106845':
            strata = {'GCST90016588': ('ever smokers', 12446, 59145), 'GCST90016589': ('never smokers', 8631, 120544), 'GCST90016593': ('current smokers', 4589, 10001), 'GCST90016594': ('noncurrent smokers', 16488, 169688)}
            stratum, cases, controls = strata[a]
            raw = next(h for h in catalog['headers'] if h['url'].endswith(a+'_buildGRCh37.tsv'))
            row.update(
                review_scope='Primary methods, original/formatted/harmonized Catalog prefixes and metadata, official PLINK output-field definitions reviewed. No full-file statistical checks.',
                new_primary_evidence=[evidence('kim2021_pmc'), evidence('plink2_formats')],
                ancestry_meta_analysis_status='UK Biobank European ancestry; single-cohort smoking-stratum marginal logistic association, not pooled multi-ancestry meta-analysis.',
                primary_publication_N={'stratum': stratum, 'cases': cases, 'controls': controls, 'source': 'Current Catalog accession metadata; paper whole analysis21,077cases/179,689controls.'},
                genome_build='Original raw GRCh37; harmonized hm_* coordinates GRCh38 by corresponding YAML. Do not combine raw allele/effect/position with hm_* counterparts.',
                public_access_status='Anonymous HTTP206 header+representative-record access successful for actual large raw file and formatted/harmonized equivalents.',
                unrestricted_download_verified=True,
                unrestricted_download_scope='Actual file prefixes verified; complete large file deliberately not downloaded in preflight.',
                verified_raw_file=raw,
                locus_completeness_status='Genome-wide release supported by study design, current Catalog9,886,854 tested variants and actual raw prefix containing non-significant records. Whole-file completeness and locus density remain future checks, not proved by header.',
                signed_effect_fields={'variant_id': 'variant_id rsID or quoted NA; preserve chromosome/position/alleles for exact keys', 'effect_allele': 'effect_allele explicit', 'other_allele': 'other_allele explicit', 'effect': 'odds_ratio; tested effect allele associated with COPD outcome', 'SE': 'standard_error; log-odds interpretation supported by logistic PLINK methods and official schema, but Catalog transformation lineage requires future confirmation', 'P': 'p_value', 'frequency': 'absent from raw; formatted/harmonized AF placeholder columns NA in inspected records', 'per_variant_N': 'absent from raw', 'variant_quality': 'absent from raw; paper inclusionMAF>=0.01 and imputation r2>=0.5 are study-level filters, not per-record QC fields'},
                risk_direction_status='Signed disease-risk direction feasible: OR relative to explicit effect allele. Future full-file checks must confirm finite positiveOR, valid paired alleles, scale and build; no risk-allele table constructed here.',
                published_locus_count={'this_stratum': None, 'not_substitutable': '48 ever-interaction and55 current-interaction loci belong to2-df joint tests, NOT the four marginal stratum GWAS.'},
                LD_requirement='UKB European cohort-matched or comparable signed LD, verified exact allele/build coverage; smoking-stratum sample composition and no per-variant N/AF require specific justification. Do not treat all four overlapping smoking strata as independent meta-analysis cohorts.',
                access_or_execution_blockers=['Whole-file download, publisher checksum and genome/locus completeness validation deferred', 'Verify logOR SE lineage and finite valid fields before converting OR to logOR', 'Per-variant N/AF/quality absent; establish acceptable source assumptions and independent allele checks', 'Bind compatible signed UKB LD and perform future locus-level QC', 'Marginal stratum-specific locus counts not yet enumerated; no eligibility count invented'],
                distinction='Direct COPD susceptibility within a smoking-defined stratum. Do not substitute interaction main-effect coefficients GCST90016586/591, interaction terms, or2-df omnibus P values.'
            )
        elif r['pubmed_id'] == '32514122':
            row.update(
                review_scope='Primary Ishigaki methods/phenotype reviewed; actual NBDC/BBJ download verification delegated to audits/alternatives, whose source-specific audit supersedes any unverified access field here.',
                new_primary_evidence=[evidence('ishigaki2020_bioc')],
                ancestry_meta_analysis_status='BioBank Japan Japanese ancestry; physician-diagnosed direct COPD. All-sex versus male/female subgroup releases remain distinct.',
                genome_build='Primary methods GRCh37/hg19; actual provider schema/build corroboration in alternatives audit.',
                public_access_status='Primary article declares unrestricted JENGER/NBDC hum0014 disease summary statistics. CatalogFTP404 is not global unavailability. See independent alternatives audit for verified NHA dataset/file retrieval.',
                locus_completeness_status='Primary study is genome-wide; exact release scope, alleles and actual header verification are recorded in alternatives audit, not inferred solely from article.',
                risk_direction_status='Provider-specific allele2/ALT beta convention must be taken from alternatives file/schema audit; do not infer effect allele from generic REF/ALT labels.',
                published_locus_count={'COPD_specific': None, 'not_substitutable': 'All42-disease study totals are not COPD locus counts.'},
                LD_requirement='BBJ or demonstrably comparable Japanese signed LD matched to released statistics. Genotype controls are controlled-access; methods using1KGphase3v5 imputation do not establish a public cohort-matched LD release.',
                access_or_execution_blockers=['See actual provider file and allele-schema audit in audits/alternatives', 'BBJ/Japanese signed LD provenance and release matching must be established', 'Whole-file/locus missingness and per-variant N/QC require future validation'],
                source_method_QC='SHAPEIT2 v2.778 and minimac3 v2.0.1 with1KGphase3v5; excludes imputationRsq<0.7; actual signed-summary file field coverage is a separate issue.'
            )
        elif r['pubmed_id'] == '35308900':
            row.update(
                review_scope='Primary paper reviewed; no public complete signed summary-file link located; generic CatalogFTP404 retained.',
                new_primary_evidence=[evidence('joo2022_pmc')],
                ancestry_meta_analysis_status='Sex-specific white British UK Biobank analyses, not meta-analysis.',
                genome_build='No actual full-file header verified; do not infer from UKB genotype build alone.',
                primary_publication_N={'male_cases': 12958, 'male_controls': 95631, 'female_cases': 11311, 'female_controls_abstract_and_catalog': 123714, 'female_controls_methods_text': 123741, 'note': 'Paper internal female-control count discrepancy retained, not silently reconciled.'},
                public_access_status='Public primary article with lead/result tables; no unrestricted complete signed genome-wide file obtained.',
                locus_completeness_status='Paper describes ~8.8million GWAS variants per sex, but publication lead tables are not locus-complete files.',
                published_locus_count={'male': 17, 'female': 14, 'shared': 9, 'note': 'Publication FUMA loci; not prospective eligible count.'},
                source_method_QC='SAIGE v0.43.3; imputationINFO>=0.7 andMAF>=0.01; age, age2, height, ever-smoking and4PCs.',
                LD_requirement='Study used UKB white-British FUMA reference for locus identification; not proof of accessible signed full locus LD. Overlaps UKB Kim/Sakornsakolpat participants.',
                access_or_execution_blockers=['Request original complete signed sex-specific summary statistics with exact sample counts and allele/build fields', 'Bind matched UKB signed LD', 'Do not reconstruct genome-wide data from17/14lead loci or gene-level results']
            )
        elif a == 'GCST011766':
            row.update(
                review_scope='Primary Moll paper and figshare API inventory reviewed.',
                new_primary_evidence=[evidence('moll2021_pmc'), evidence('moll2021_figshare_api')],
                ancestry_meta_analysis_status='Mixed-ancestry twelve-cohort meta-analysis, UKB+ICGC; not whole-genome dense coverage.',
                public_access_status='Public supplement DOI10.6084/m9.figshare.14538222 lists only a6,177,157-byte DOCX appendix; no complete dense genome-wide association file verified.',
                locus_completeness_status='Analysis deliberately selected protein-altering exonic variants:109,036initial functional array variants,20,536retained; inappropriate as dense regional full GWAS input even if every retained exonic statistic were obtained.',
                published_locus_count={'exome_significant_variants': 80, 'LD_clumped_lead_variants': 35, 'note': 'Not35dense fine-mapping loci; functional preselection violates region completeness.'},
                LD_requirement='No LD reference can restore association statistics for unassayed/nonselected noncoding variants.',
                access_or_execution_blockers=['Unsuitable as primary dense locus fine-mapping input owing to exonic ascertainment', 'Do not mislabel appendix or20,536selected-variant release as genome-wide locus-complete data']
            )
        elif r['pubmed_id'] == '24621683':
            row.update(
                review_scope='Primary Cho2014 full article reviewed; no unrestricted full signed file verified.',
                new_primary_evidence=[evidence('cho2014_bioc')],
                ancestry_meta_analysis_status='European-ancestry and African-American COPDGene plus ECLIPSE/NETT-NAS/GenKOLS fixed-effects meta-analysis. ICGN follow-up only selected top loci.',
                primary_publication_N={'moderate_to_severe_cases': 6633, 'severe_cases': 3497, 'controls': 5704},
                locus_completeness_status='Discovery genome-wide, but selected ICGN follow-up jointly meta-analyzed only at top loci; public paper/result tables not complete signed files.',
                published_locus_count={'across_moderate_and_severe_analyses': 6, 'known': 3, 'additional': 3, 'note': 'Paper headline across related analyses; do not equate with either accession locus eligibility.'},
                source_method_QC='PLINK1.07 logistic per cohort/race adjusted age/pack-years/PCs; METAL fixed effects; imputation1000GphaseIv3 European or cosmopolitan; marker required pass in all genome-wide cohorts.',
                access_or_execution_blockers=['Request original complete signed discovery summary statistics, preferably ancestry-specific', 'Do not infer final whole-genome completeness from selected ICGN follow-up', 'Matched ancestry-specific signed LD required']
            )
        elif a in ['GCST000603', 'GCST001321']:
            label = 'cho2010_bioc' if a == 'GCST000603' else 'cho2011_pmc'
            row.update(
                review_scope='Primary article reviewed in addition to frozen phenotype/current Catalog audits; no complete signed public file verified.',
                new_primary_evidence=[evidence(label)],
                genome_build='Article coordinates explicitly hg18/NCBI36; no actual full-summary file header obtained.',
                ancestry_meta_analysis_status='European-ancestry discovery cohorts; selected-locus family/case-control replication is not independent dense genome-wide release.',
                published_locus_count={'headline': 'Cho2010 identifies FAM13A additional susceptibility locus; Cho2011 increases then-known COPD GWAS loci to4. These are publication-history claims, not executable-locus counts.'},
                access_or_execution_blockers=['Request original complete signed discovery file and exact release schema', 'Resolve hg18 build and harmonization without reusing lead tables', 'Obtain compatible signed LD and account for shared cohorts with later COPD meta-analyses']
            )
        sources.append(row)

    ledger = []
    for path in sorted(ACQUIRED.glob('*/access.json')):
        item = json.loads(path.read_text())
        body = STAGE / item['body_path']
        assert sha(body) == item['acquired_sha256']
        h = {k.lower(): v for k, v in item['response_headers'].items()}
        ledger.append({
            'label': path.parent.name, 'access_record': str(path.relative_to(STAGE)),
            'url': item['url'], 'effective_url': item['effective_url'],
            'requested_utc': item['requested_utc'], 'completed_utc': item['completed_utc'],
            'http_status': item['http_status'], 'transport_error': item['transport_error'],
            'acquired_bytes': item['acquired_bytes'], 'acquired_sha256': item['acquired_sha256'],
            'body_path': item['body_path'], 'read_truncated_by_client': item['read_truncated'],
            'remote_object_partial': bool(item['read_truncated'] or item['http_status'] == 206),
            'content_range': h.get('content-range'),
            'payload_semantics': 'HTTP200 alone does not imply requested scientific payload; BioC6148-byte text/html responses are error/landing documents, not articleXML.',
            'hash_scope': 'Exact acquired response bytes only. For4 explicitly audited complete small dbGaP exports, schema audit separately proves Content-Length match. These exports are still not complete originalGWAS data.'
        })
    write_json('source_access_ledger.json', ledger)
    result = {
        'stage': 'fine-mapping-preflight-1.0', 'audit': 'direct_COPD_access_and_schema',
        'generated_utc': datetime.now(timezone.utc).isoformat(),
        'selection_firewall': 'Selected all27direct-COPD accessions solely by frozen phenotype class from104-study candidate-blind metadata projection. No candidate membership, model outputs, benchmark outcomes or candidate-support counts used.',
        'frozen_metadata_path': str(metadata_path.relative_to(STAGE)),
        'frozen_metadata_sha256': sha(metadata_path),
        'study_count': len(sources), 'studies': sources,
        'dbgap_schema_audit_path': 'audits/direct/dbgap_public_export_schema_audit.json',
        'access_ledger_path': 'audits/direct/source_access_ledger.json',
        'analysis_prohibition_respected': 'No fine mapping, association recomputation, LD computation, candidate/locus selection, model scoring, colocalization, or retraining performed. Small public-export row/schema checks are release-scope audit only.',
        'global_recommendation': 'Prefer direct phenotype. Kim smoking-stratum files are technically retrievable with signed OR but not execution-cleared until matchedLD and full-file QC. Ishigaki BBJ direct phenotype needs source-specific alternatives/LD audit integration. Sakornsakolpat and Hobbs public exports are unsuitable; request complete signed source files. No prospective eligible-locus count can be claimed from this audit.'
    }
    write_json('direct_source_audit.json', result)
    print(json.dumps({'studies': len(sources), 'direct_source_audit_sha256': sha(OUT/'direct_source_audit.json'),
                      'access_records': len(ledger), 'dbgap_rows': [r['record_count'] for r in exports]}, indent=2))


if __name__ == '__main__':
    main()
