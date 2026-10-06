#!/usr/bin/env python3
"""Source-only benchmark assembly. This script never reads model outputs.

Stages: fetch authoritative rsID mappings; assemble auditable tables; freeze.
Frozen sources/tables cannot be overwritten by this script after freezing.
"""
import argparse
import csv
import gzip
import hashlib
import json
import re
import time
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import pysam
import requests

ROOT = Path(__file__).resolve().parents[4]
BASE = ROOT / 'diseases/COPD/07_gap_closure'
DATA = BASE / 'data/COPD-V2-BENCH'
PROV = BASE / 'provenance'
RES = BASE / 'results'
PREFIX = 'COPD-V2-BENCH'
SPEC_SHA = '46d3aac26b1f980d7c2dcd5034d6cc64566524ca7d99cc635f7d979e04f2de65'
FASTA = ROOT / 'models/TREDNET_v2/fasta/hg38.fa'
SOURCE_FILES = [DATA / (name + '_evidence.json') for name in ['castaldi', 'gong', 'mechanisms']]
FREEZE = PROV / f'{PREFIX}_benchmark_freeze.json'
STATES = {'positive', 'null', 'ambiguous', 'conflicting', 'unavailable', 'unevaluable'}
CLASSES = {'MPRA_allele_effect', 'conventional_reporter', 'endogenous_allele_editing',
           'CRISPRi_region', 'CRISPR_deletion_region', 'chromatin_contact', 'TF_binding',
           'splicing', 'alternative_polyadenylation', 'expression_QTL', 'other_or_unresolved'}


def now():
    return datetime.now(timezone.utc).isoformat()


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def dump(path, obj):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(obj, indent=2, ensure_ascii=False, sort_keys=True) + '\n')


def tsv(path, rows, first=()):
    fields = list(first) + sorted(set().union(*(set(row) for row in rows)) - set(first))
    with Path(path).open('w', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        for row in rows:
            w.writerow({key: json.dumps(value, ensure_ascii=False, sort_keys=True) if isinstance(value, (dict, list))
                        else '' if value is None else value for key, value in row.items()})


def read_sources(require_all=True):
    out = []
    for path in SOURCE_FILES:
        if path.exists():
            out.append((path, json.loads(path.read_text())))
        elif require_all:
            raise ValueError(f'Missing source adjudication: {path}')
    return out


def mapping_cache(build):
    combined = {}
    for path in sorted((DATA / 'allele_mappings').glob(f'{build}_*.json.gz')):
        with gzip.open(path, 'rt') as f:
            obj = json.load(f)
        if isinstance(obj, dict):
            for key, value in obj.items():
                if re.fullmatch(r'rs\d+', key) and isinstance(value, dict):
                    combined[key] = value
    return combined


def fetch(rsids):
    target = DATA / 'allele_mappings'
    target.mkdir(parents=True, exist_ok=True)
    for build, host in [('GRCh38', 'https://rest.ensembl.org'),
                        ('GRCh37', 'https://grch37.rest.ensembl.org')]:
        cached = mapping_cache(build)
        missing = sorted(set(rsids) - set(cached), key=lambda x: int(x[2:]))
        for start in range(0, len(missing), 100):
            ids = missing[start:start + 100]
            key = hashlib.sha256('\n'.join(ids).encode()).hexdigest()[:16]
            path = target / f'{build}_{key}.json.gz'
            url = host + '/variation/homo_sapiens'
            stamp = now()
            attempts = []
            for retry in range(3):
                response = requests.post(url, headers={'Content-Type': 'application/json', 'Accept': 'application/json'},
                                         json={'ids': ids}, timeout=120)
                attempts.append({'time_utc': now(), 'status_code': response.status_code})
                if response.ok:
                    break
                time.sleep(2 * (retry + 1))
            if not response.ok:
                dump(target / f'{build}_{key}.failure.json',
                     {'url': url, 'requested_ids': ids, 'retrieved_utc': stamp,
                      'status_code': response.status_code, 'body': response.text[:4000]})
                raise RuntimeError(f'Authoritative mapping retrieval failed: {response.status_code} {url}')
            obj = response.json()
            if not isinstance(obj, dict):
                raise ValueError('Invalid mapping response')
            with path.open('wb') as raw:
                with gzip.GzipFile(fileobj=raw, mode='wb', mtime=0) as zipped:
                    zipped.write(response.content)
            dump(target / f'{build}_{key}.retrieval.json',
                 {'url': url, 'method': 'POST', 'requested_ids': ids, 'retrieved_utc': stamp,
                  'completed_utc': now(), 'status_code': response.status_code,
                  'content_type': response.headers.get('Content-Type'), 'response_bytes': len(response.content),
                  'attempts': attempts,
                  'response_sha256': hashlib.sha256(response.content).hexdigest(),
                  'local_path': str(path.relative_to(ROOT)), 'sha256': sha(path),
                  'response_keys': sorted(obj)})
            print(f'{build}: retrieved {len(ids)} identifiers ({start + len(ids)}/{len(missing)})', flush=True)
            time.sleep(.2)


def merged_aliases(rsids):
    target = DATA / 'allele_mappings'
    current = mapping_cache('GRCh38')
    missing = [rs for rs in sorted(set(rsids)) if not canonical_mappings(current.get(rs, {}), 'GRCh38')]
    alias = {}
    for rsid in missing:
        url = 'https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/' + rsid[2:]
        response = requests.get(url, timeout=45)
        path = target / ('NCBI_' + rsid + '.json.gz')
        with path.open('wb') as raw:
            with gzip.GzipFile(fileobj=raw, mode='wb', mtime=0) as zipped:
                zipped.write(response.content)
        dump(target / ('NCBI_' + rsid + '.retrieval.json'),
             {'url': url, 'status_code': response.status_code, 'retrieved_utc': now(),
              'response_sha256': hashlib.sha256(response.content).hexdigest(), 'local_sha256': sha(path)})
        if response.ok:
            merged = response.json().get('merged_snapshot_data', {}).get('merged_into', [])
            if len(merged) == 1:
                alias[rsid] = 'rs' + str(merged[0])
        time.sleep(.4)
    if alias:
        fetch(alias.values())
        for build in ['GRCh38', 'GRCh37']:
            cached = mapping_cache(build)
            resolved = {rs: dict(cached[to], source_merged_rsid=rs,
                                merge_evidence_url='https://api.ncbi.nlm.nih.gov/variation/v0/refsnp/' + rs[2:])
                        for rs, to in alias.items() if to in cached}
            path = target / (build + '_resolved_merged_aliases.json.gz')
            if path.exists():
                with gzip.open(path, 'rt') as f:
                    previous = json.load(f)
                previous.update(resolved)
                resolved = previous
            with path.open('wb') as raw:
                with gzip.GzipFile(fileobj=raw, mode='wb', mtime=0) as zipped:
                    zipped.write(json.dumps(resolved, sort_keys=True).encode())
        dump(target / 'merged_alias_resolution.json', {'resolved_utc': now(), 'aliases': alias,
                                                     'method': 'NCBI explicit rsID merge, never LD substitution'})


def complement(allele):
    return allele.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def clean_allele(value):
    text = str(value or '').strip().upper()
    return text if text and not (set(text) - set('ACGT')) else ''


def canonical_mappings(obj, build):
    return [x for x in obj.get('mappings', []) if x.get('assembly_name') == build and
            str(x.get('seq_region_name')) in set(map(str, range(1, 23))) | {'X', 'Y'} and
            x.get('coord_system') == 'chromosome' and int(x.get('strand', 0)) == 1]


def normalize_vcf(chrom, pos, ref, alt, genome):
    """Left-normalize a reference-backed biallelic replacement without model inputs."""
    if not ref or not alt:
        anchor = genome.fetch('chr' + chrom, pos - 2, pos - 1).upper()
        pos, ref, alt = pos - 1, anchor + ref, anchor + alt
    while ref[-1] == alt[-1] and ref != alt:
        if min(len(ref), len(alt)) == 1:
            if pos <= 1:
                break
            anchor = genome.fetch('chr' + chrom, pos - 2, pos - 1).upper()
            pos, ref, alt = pos - 1, anchor + ref[:-1], anchor + alt[:-1]
        else:
            ref, alt = ref[:-1], alt[:-1]
    while min(len(ref), len(alt)) > 1 and ref[0] == alt[0]:
        pos, ref, alt = pos + 1, ref[1:], alt[1:]
    if genome.fetch('chr' + chrom, pos - 1, pos - 1 + len(ref)).upper() != ref:
        raise ValueError('Normalized reference verification failed')
    return pos, ref, alt


def occurrences(sequence, query):
    found, start = set(), 0
    while True:
        at = sequence.find(query, start)
        if at < 0:
            return found
        found.add(at)
        start = at + 1


def construct_mapping(row, mapping, genome):
    """Require the deposited allele inserts to match one reference/alternate contrast.

    This creates only short reference-check windows, never new model sequences.
    Insert orientation/assay design remains separate from the biological allele.
    """
    first = clean_allele(row.get('allele1_forward_insert'))
    second = clean_allele(row.get('allele2_forward_insert'))
    if not first or not second:
        return None
    chrom, pos = str(mapping['seq_region_name']), int(mapping['start'])
    db = mapping['allele_string'].split('/')
    ref = '' if db[0] == '-' else clean_allele(db[0])
    if not ref and db[0] != '-':
        return {'status': 'invalid_database_reference'}
    if genome.fetch('chr' + chrom, pos - 1, pos - 1 + len(ref)).upper() != ref:
        return {'status': 'database_reference_mismatch'}
    start, end = max(0, pos - 201), pos + len(ref) + 200
    wild = genome.fetch('chr' + chrom, start, end).upper()
    offset = pos - 1 - start
    matches = []
    for db_alt in db[1:]:
        alt = '' if db_alt == '-' else clean_allele(db_alt)
        if not alt and db_alt != '-':
            continue
        mutant = wild[:offset] + alt + wild[offset + len(ref):]
        for ref_role, ref_insert, alt_insert in [(1, first, second), (2, second, first)]:
            common = occurrences(wild, ref_insert) & occurrences(mutant, alt_insert)
            common = {x for x in common if x <= offset < x + len(ref_insert)}
            if common:
                normalized = normalize_vcf(chrom, pos, ref, alt, genome)
                matches.append((ref_role, normalized, sorted(common)))
    choices = {(r, v) for r, v, _ in matches}
    if len(choices) != 1:
        # An isolated SNP may have a shared nonreference background in both oligos.
        # Require a unique exact terminal anchor at the reported genomic nucleotide.
        # Never use this fallback for an indel or a multi-base assay contrast.
        differences = [i for i, (a, b) in enumerate(zip(first, second)) if a != b]
        if len(ref) == 1 and int(mapping['end']) == pos and len(first) == len(second) and len(differences) == 1:
            d = differences[0]
            bases = {first[d], second[d]}
            source_pair = {clean_allele(row.get('tested_allele1')), clean_allele(row.get('tested_allele2'))}
            if bases == source_pair and ref in bases and bases <= set(db):
                implied_starts = set()
                for at in (0, len(first) - 20):
                    if at <= d < at + 20 or first[at:at + 20] != second[at:at + 20]:
                        continue
                    hits = occurrences(wild, first[at:at + 20])
                    if len(hits) == 1:
                        implied_starts.add(next(iter(hits)) - at)
                if (len(implied_starts) == 1 and next(iter(implied_starts)) + d == offset and
                    0 <= next(iter(implied_starts)) <= len(wild) - len(first)):
                    sequence_start = start + next(iter(implied_starts))
                    reference_context = genome.fetch('chr' + chrom, sequence_start, sequence_start + len(first)).upper()
                    shared = [{'pos1': sequence_start + i + 1, 'reference': reference_context[i], 'assayed': a}
                              for i, a in enumerate(first) if i != d and a != reference_context[i]]
                    return {'status': 'both_allele_constructs_verified', 'chrom': chrom, 'pos': pos,
                            'ref': ref, 'alt': next(iter(bases - {ref})),
                            'source_reference_allele_role': 1 if first[d] == ref else 2,
                            'validation_detail': 'Isolated SNP; unique exact terminal20bp anchor confirms rsID coordinate; shared background preserved',
                            'shared_nonreference_background': shared}
        return {'status': 'constructs_do_not_identify_one_authoritative_reference_alt_contrast',
                'matching_contrasts': len(choices)}
    ref_role, (normpos, normref, normalt) = next(iter(choices))
    return {'status': 'both_allele_constructs_verified', 'chrom': chrom, 'pos': normpos,
            'ref': normref, 'alt': normalt, 'source_reference_allele_role': ref_role,
            'matching_window_offsets': matches[0][2]}


def harmonize(row, grch38, grch37, genome):
    rsid = str(row.get('rsid', '')).strip()
    a, b = clean_allele(row.get('tested_allele1')), clean_allele(row.get('tested_allele2'))
    pair = {a, b}
    result = {key: row.get(key, '') for key in ['assay_id', 'study_id', 'rsid', 'source_build',
              'source_chrom', 'source_pos', 'tested_allele1', 'tested_allele2', 'reported_strand']}
    result.update({'canonical_variant_id': '', 'grch38_chrom': '', 'grch38_pos': '', 'grch38_ref': '',
                   'grch38_alt': '', 'identity_status': 'unevaluable', 'identity_reason': '',
                   'harmonization_method': 'Ensembl dbSNP GRCh38 primary-chromosome mapping plus existing hg38 FASTA',
                   'liftover_used': False, 'liftover_roundtrip': 'not_applicable_no_liftover',
                   'reference_verified': False, 'source_coordinate_status': 'not_provided',
                   'palindromic': bool(pair in ({'A', 'T'}, {'C', 'G'})),
                   'allele_transform': '', 'direction_identity_resolved': False,
                   'higher_activity_grch38_allele': '', 'reported_direction_alt_minus_ref': '',
                   'mapping_source_url': f'https://rest.ensembl.org/variation/human/{rsid}?content-type=application/json'})
    if not re.fullmatch(r'rs\d+', rsid):
        result['identity_reason'] = 'No exact rsID; coordinate-only source requires separate verified mapping'
        return result
    mapping = canonical_mappings(grch38.get(rsid, {}), 'GRCh38')
    if len(mapping) != 1:
        result['identity_reason'] = f'Expected one GRCh38 primary mapping, observed {len(mapping)}'
        return result
    m = mapping[0]
    chrom, pos = str(m['seq_region_name']), int(m['start'])
    alleles = m.get('allele_string', '').split('/')
    ref = clean_allele(alleles[0])
    result.update({'grch38_chrom': chrom, 'grch38_pos': pos, 'grch38_ref': ref,
                   'database_alleles': '/'.join(alleles), 'resolved_rsid': grch38[rsid].get('name', rsid)})
    if row.get('construct_exact_allele_identity_valid') is False:
        result['identity_reason'] = 'Deposited assay-orientation construct does not preserve the reported biological allele contrast'
        result['construct_reference_check'] = 'source_orientation_construct_invalid'
        return result
    construct = construct_mapping(row, m, genome)
    if construct is not None:
        result['construct_reference_check'] = construct['status']
        result['construct_reference_details'] = construct
        if construct['status'] != 'both_allele_constructs_verified':
            result['identity_reason'] = 'Deposited forward allele inserts fail unique authoritative/reference contrast validation'
            return result
        chrom, pos, ref, alt = (construct[k] for k in ('chrom', 'pos', 'ref', 'alt'))
        result.update({'canonical_variant_id': f'{chrom}:{pos}:{ref}:{alt}',
                       'grch38_chrom': chrom, 'grch38_pos': pos, 'grch38_ref': ref, 'grch38_alt': alt,
                       'reference_verified': True, 'reference_observed': ref,
                       'identity_status': 'exact', 'identity_reason': 'Both deposited allele inserts match one reference-backed database contrast',
                       'allele_transform': 'forward_construct_reference_alignment_and_left_normalization',
                       'direction_identity_resolved': True,
                       'harmonization_method': 'Authoritative GRCh38 rsID mapping, both deposited FOR allele inserts, hg38 FASTA and VCF left normalization'})
        return result
    if not ref or pos != int(m['end']) or len(ref) != 1:
        result['identity_reason'] = 'Non-SNV mapping requires separately validated normalized indel contrast'
        return result
    observed = genome.fetch('chr' + chrom, pos - 1, pos).upper()
    result['reference_observed'] = observed
    result['reference_verified'] = observed == ref
    if observed != ref:
        result['identity_reason'] = 'Authoritative reference and local GRCh38 sequence disagree'
        return result
    source_build_text = str(row.get('source_build', '')).lower()
    build = 'grch38' if 'grch38' in source_build_text or 'hg38' in source_build_text else \
            'grch37' if 'grch37' in source_build_text or 'hg19' in source_build_text else source_build_text
    source_pos = row.get('source_pos', '')
    source_chrom = str(row.get('source_chrom', '')).replace('chr', '')
    if source_pos not in ('', None) and source_chrom:
        source_mappings = mapping if build in ('grch38', 'hg38') else canonical_mappings(grch37.get(rsid, {}), 'GRCh37') if build in ('grch37', 'hg19') else []
        try:
            source_matches = any(str(x['seq_region_name']) == source_chrom and int(x['start']) == int(float(source_pos)) for x in source_mappings)
        except (ValueError, TypeError):
            source_matches = False
        result['source_coordinate_status'] = 'verified_database_assembly_mapping' if source_matches else 'unresolved_source_build_or_coordinate'
        if not source_mappings:
            coordinate_matches = []
            for assembly, candidates in [('GRCh38', mapping), ('GRCh37', canonical_mappings(grch37.get(rsid, {}), 'GRCh37'))]:
                try:
                    if any(str(x['seq_region_name']) == source_chrom and int(x['start']) == int(float(source_pos)) for x in candidates):
                        coordinate_matches.append(assembly)
                except (ValueError, TypeError):
                    pass
            result['source_coordinate_matching_database_assemblies'] = ';'.join(coordinate_matches)
            if coordinate_matches:
                result['source_coordinate_status'] = 'database_position_match_but_source_build_not_explicitly_reported_or_verified'
        if source_mappings and not source_matches:
            result['identity_reason'] = 'Reported source coordinate disagrees with named assembly rsID mapping'
            return result
    if not a or not b or a == b:
        result['identity_reason'] = 'Experimentally tested allele pair not reported/resolved'
        return result
    strand = str(row.get('reported_strand', '')).strip().lower()
    plus = strand in ('+', '+1', '1', 'forward', 'plus', 'genomic_forward')
    minus = strand in ('-', '-1', 'reverse', 'minus', 'genomic_reverse')
    possibilities = []
    for transform, x, y in [('identity', a, b), ('reverse_complement', complement(a), complement(b))]:
        if (plus and transform != 'identity') or (minus and transform != 'reverse_complement'):
            continue
        if len(x) == len(y) == 1 and ref in (x, y) and x in alleles and y in alleles:
            possibilities.append((transform, x, y, y if x == ref else x))
    alt_set = {entry[3] for entry in possibilities}
    if len(alt_set) != 1:
        result['identity_reason'] = 'Tested contrast cannot be uniquely reconciled to reference-backed GRCh38 REF/ALT'
        return result
    alt = next(iter(alt_set))
    identity_resolved = len(possibilities) == 1
    result.update({'canonical_variant_id': f'{chrom}:{pos}:{ref}:{alt}', 'grch38_alt': alt,
                   'identity_status': 'exact', 'identity_reason': 'Tested allele contrast and unique primary mapping agree with GRCh38 reference',
                   'allele_transform': possibilities[0][0] if identity_resolved else 'palindromic_strand_unresolved',
                   'direction_identity_resolved': identity_resolved})
    higher = clean_allele(row.get('higher_activity_allele') or row.get('effect_allele', ''))
    if identity_resolved and higher in pair:
        high = higher if possibilities[0][0] == 'identity' else complement(higher)
        result['higher_activity_grch38_allele'] = high
        result['reported_direction_alt_minus_ref'] = 1 if high == alt else -1
    return result


def assemble():
    sources = read_sources()
    assays, contexts, inventory, searches, denominator = [], [], [], [], []
    for path, obj in sources:
        for key, target in [('assays', assays), ('contexts', contexts), ('sources', inventory), ('searches', searches)]:
            for item in obj.get(key, []):
                if isinstance(item, str):
                    item = {'description': item}
                target.append(dict(item, adjudication_file=str(path.relative_to(ROOT))))
        audit = obj.get('denominator_audit', {})
        for entry in (audit if isinstance(audit, list) else [audit]):
            if 'study_id' in entry:
                denominator.append(dict(entry, adjudication_file=str(path.relative_to(ROOT))))
            else:
                for study in sorted({r['study_id'] for r in obj.get('assays', [])}):
                    denominator.append(dict(entry, study_id=study, complete_assayed_denominator_valid=False,
                                            adjudication_file=str(path.relative_to(ROOT))))
    ids = [r['assay_id'] for r in assays]
    if len(ids) != len(set(ids)):
        raise ValueError('Duplicate assay identifiers')
    for r in assays:
        if r['experimental_state'] not in STATES:
            raise ValueError(f'Unknown experimental state {r}')
        if r['assay_class'] not in CLASSES:
            raise ValueError(f'Unknown assay class {r}')
        if r.get('mechanism_in_model_scope') not in {'yes', 'partial', 'no', 'unevaluable'}:
            raise ValueError(f'Missing typed mechanism scope {r}')
        if not r.get('source_locator') or not r.get('source_url'):
            raise ValueError(f'Missing experimental label provenance {r}')
    grch38, grch37 = mapping_cache('GRCh38'), mapping_cache('GRCh37')
    source_regions = next(obj['denominator_audit'].get('locus_source_regions', [])
                          for path, obj in sources if path.name == 'gong_evidence.json')
    locus_rows = []
    for rsid in sorted({r['rsid'] for r in assays if r['study_id'] == 'GONG2026'}):
        placements = canonical_mappings(grch37.get(rsid, {}), 'GRCh37')
        matching = [(x, reg) for x in placements for reg in source_regions
                    if str(x['seq_region_name']) == reg['chrom'] and reg['start'] <= int(x['start']) <= reg['end']]
        labels = {reg['locus'] for _, reg in matching}
        label = next(iter(labels)) if len(labels) == 1 else 'unresolved_source_locus'
        locus_rows.append({'study_id': 'GONG2026', 'rsid': rsid, 'locus': label,
                           'method': 'Authoritative GRCh37 rsID position within published Figure S1B hg19 MPRA display ranges; not target-gene inference or new LD',
                           'matching_source_regions': [reg for _, reg in matching]})
    locus_by_rsid = {r['rsid']: r for r in locus_rows}
    for assay in assays:
        if assay['study_id'] == 'GONG2026':
            loc = locus_by_rsid[assay['rsid']]
            assay['reported_locus_annotation'] = assay.get('locus', '')
            assay['locus'] = loc['locus']
            assay['locus_assignment_method'] = loc['method']
    with pysam.FastaFile(str(FASTA)) as genome:
        harmonized = [harmonize(r, grch38, grch37, genome) for r in assays]
    master = []
    for assay, mapping in zip(assays, harmonized):
        entry = dict(assay, **{k: v for k, v in mapping.items() if k not in assay})
        entry['benchmark_version'] = '1.0'
        entry['external_evaluation_only'] = True
        entry['model_outcomes_used_for_inclusion'] = False
        master.append(entry)
    paths = {}
    classification = [{k: r.get(k, '') for k in ['assay_id', 'study_id', 'rsid', 'canonical_variant_id',
                      'assay_class', 'cell_context', 'assay_orientation', 'experimental_state',
                      'mechanism_in_model_scope', 'mechanism_scope_rationale', 'evidence_role',
                      'source_locator', 'source_url']} for r in master]
    for name, rows, first in [
        ('R001_source_study_inventory.tsv', inventory, ('study_id',)),
        ('R002_exact_allele_harmonization.tsv', harmonized, ('assay_id', 'study_id', 'rsid', 'canonical_variant_id')),
        ('R003_frozen_benchmark_master.tsv', master, ('assay_id', 'study_id', 'rsid', 'canonical_variant_id', 'assay_class', 'cell_context', 'experimental_state')),
        ('R004_assay_mechanism_classification.tsv', classification, ('assay_id', 'study_id', 'rsid')),
        ('R004B_contextual_and_excluded_evidence.tsv', contexts, ('study_id',)),
        ('R005_assayed_denominator_audit.tsv', denominator, ('study_id',)),
        ('R006_source_search_audit.tsv', searches, ()),
        ('R006B_source_locus_assignment.tsv', locus_rows, ('study_id', 'rsid', 'locus'))]:
        path = RES / f'{PREFIX}-{name}'
        tsv(path, rows, first)
        paths[str(path.relative_to(ROOT))] = {'sha256': sha(path), 'rows': len(rows), 'bytes': path.stat().st_size}
    manifest = {'assembled_utc': now(), 'protocol_sha256': SPEC_SHA, 'assay_rows': len(master),
                'unique_reported_rsids': len({r['rsid'] for r in master if r.get('rsid')}),
                'exact_assay_rows': sum(r['identity_status'] == 'exact' for r in master),
                'unique_exact_variants': len({r['canonical_variant_id'] for r in master if r['canonical_variant_id']}),
                'experimental_states': dict(Counter(r['experimental_state'] for r in master)),
                'Gong_source_locus_counts': dict(Counter(r['locus'] for r in locus_rows)),
                'Gong_figure_stated_locus_counts': {r['locus']: r['designed_variant_count'] for r in source_regions},
                'Gong_locus_count_caveat': 'Source Figure S1A counts and rsID-position assignments differ slightly; retained explicitly without adjusting boundaries or labels',
                'source_files': {str(p.relative_to(ROOT)): sha(p) for p, _ in sources},
                'reference': {'path': str(FASTA.relative_to(ROOT)), 'resolved_path': str(FASTA.resolve()),
                              'bytes': FASTA.stat().st_size, 'sha256': sha(FASTA), 'index_sha256': sha(str(FASTA) + '.fai')},
                'outputs': paths, 'model_outputs_read': False}
    dump(PROV / f'{PREFIX}_assembly_manifest.json', manifest)
    print(json.dumps({k: v for k, v in manifest.items() if k not in ('outputs', 'source_files', 'reference')}, indent=2))


def freeze():
    assembly = json.loads((PROV / f'{PREFIX}_assembly_manifest.json').read_text())
    for path, details in assembly['outputs'].items():
        if sha(ROOT / path) != details['sha256']:
            raise ValueError('Assembly output changed before freeze: ' + path)
    for path, expected in assembly['source_files'].items():
        if sha(ROOT / path) != expected:
            raise ValueError('Source adjudication changed: rerun assembly before freeze')
    # Freeze all raw evidence, extraction/assembly code and derived source tables.
    locked = {}
    for path in sorted(DATA.rglob('*')):
        if path.is_file():
            locked[str(path.relative_to(ROOT))] = sha(path)
    for path in (BASE / 'scripts').glob('*benchmark*.py'):
        if path.name.startswith(('collect_', 'assemble_', 'initialize_')):
            locked[str(path.relative_to(ROOT))] = sha(path)
    for path, details in assembly['outputs'].items():
        locked[path] = details['sha256']
    for suffix in ['analysis_specification.md', 'protocol_lock.json', 'frozen_input_hashes_before.tsv',
                   'assembly_manifest.json', 'source_identity_review.json', 'experimental_label_review.json']:
        path = PROV / f'{PREFIX}_{suffix}'
        if not path.exists():
            raise ValueError('Required pre-freeze review/provenance missing: ' + str(path))
        locked[str(path.relative_to(ROOT))] = sha(path)
    manifest = {'benchmark_version': '1.0', 'status': 'FROZEN_EXTERNAL_EVALUATION_DATA',
                'frozen_utc': now(), 'protocol_sha256': SPEC_SHA,
                'assembly_manifest_sha256': sha(PROV / f'{PREFIX}_assembly_manifest.json'),
                'master_path': f'diseases/COPD/07_gap_closure/results/{PREFIX}-R003_frozen_benchmark_master.tsv',
                'frozen_files': locked, 'new_model_comparison_started': False,
                'no_training_or_model_selection_use': True,
                'known_prior_selected_case_observations_are_not_blinded': True}
    dump(FREEZE, manifest)
    print(json.dumps({'frozen_utc': manifest['frozen_utc'], 'files': len(locked),
                      'master_sha256': locked[manifest['master_path']]}, indent=2))


def main():
    p = argparse.ArgumentParser()
    p.add_argument('stage', choices=['fetch', 'resolve-aliases', 'assemble', 'freeze'])
    p.add_argument('--rsids', nargs='*', default=[])
    p.add_argument('--rsid-file', type=Path)
    args = p.parse_args()
    if sha(PROV / f'{PREFIX}_analysis_specification.md') != SPEC_SHA:
        raise SystemExit('Locked protocol hash mismatch')
    if FREEZE.exists():
        raise SystemExit('Benchmark already frozen; source/label mutation refused')
    if args.stage in ('fetch', 'resolve-aliases'):
        ids = set(args.rsids)
        if args.rsid_file:
            ids.update(x for line in args.rsid_file.read_text().splitlines() for x in line.split('\t') if re.fullmatch(r'rs\d+', x))
        for _, obj in read_sources(require_all=False):
            ids.update(r['rsid'] for key in ['assays', 'contexts'] for r in obj.get(key, [])
                       if re.fullmatch(r'rs\d+', str(r.get('rsid', ''))))
        (fetch if args.stage == 'fetch' else merged_aliases)(ids)
    elif args.stage == 'assemble':
        assemble()
    else:
        freeze()


if __name__ == '__main__':
    main()
