"""Integrate completed audits without recomputing association statistics or eligibility."""
from pathlib import Path
import csv
import datetime
import json
import sys

s = Path(__file__).resolve().parents[1]
final = '--final' in sys.argv


def read(path):
    return json.loads((s/path).read_text())


def tsv(path, rows):
    assert rows
    with (s/path).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)


audit = {'A': read('tracks/A/audit.json'), 'B': read('tracks/B/results/audit.json'),
         'C': read('tracks/C/results/audit.json')}
overlap = {t: {v['locus_id']: v for v in read(f'ld/{t}/gwas_ld_overlap_summary.json')}
           for t in ['A', 'C']}
counts = {'A': (9886853, 11, 10, 1, 140886), 'B': (8678470, 5, 5, 0, 51989),
          'C': (28987534, 7, 6, 1, 175749)}
assert counts['A']==tuple(audit['A'][k] for k in ['full_gwas_rows','prospective_loci_total','prospective_non_mhc_loci','deferred_mhc_loci','risk_direction_assigned_all_including_deferred_mhc'])
assert counts['B']==tuple(audit['B'][k] for k in ['complete_GWAS_rows','prospective_loci_before_MHC','prospective_loci','MHC_deferred_loci','valid_signed_risk_direction_rows'])
assert counts['C']==tuple(audit['C'][k] for k in ['complete_gwas_rows','prospective_merged_loci','prospective_non_MHC_loci','MHC_deferred_loci','signed_direction_assigned_locus_rows'])
blockers = {
 'A': 'Unresolved strong variant identities/LD loss in 9 loci; per-variant/effective N absent; smoking-subset versus full-EUR LD and summary likelihood compatibility unresolved; source precision failures quarantined',
 'B': 'No accessible defensible dense signed Japanese LD; exact SAIGE score/software lineage unresolved; SPA-calibrated z and compatible covariance/effective N unresolved',
 'C': 'SPA-derived SE does not establish Gaussian/RSS likelihood; phenotype N420531 versus LD N420542 and per-variant/effective N unresolved; exact executed LD release configuration unbound; summary-LD covariance compatibility unresolved'}
roles = {'A': 'PRIMARY_PRACTICAL_DIRECT_COPD', 'B': 'INDEPENDENT_ANCESTRY_DIRECT_COPD',
         'C': 'SECONDARY_EHR_SENSITIVITY'}
sources = {'A': 'Pan-UKB UKBB.EUR.ldadj.bm; EXTERNAL_FULL_EUR_REFERENCE_NOT_KIM_SUBGROUP',
           'B': 'NONE; controlled Japanese genotype/reference access required',
           'C': 'Pan-UKB UKBB.EUR.ldadj.bm; SAME_PROJECT_NEAR_OVERLAP_NOT_PROVEN_IDENTICAL_SAMPLE'}
all_loci = []; readiness = []; numeric = []; compatibility = []; gates = []
completed = 0
for t in 'ABC':
    with (s/f'tracks/{t}/results/loci.tsv').open() as f:
        loci = list(csv.DictReader(f, delimiter='\t'))
    assert len(loci) == counts[t][1]
    for l in loci:
        lid = l['locus_id']; start = int(l['start']); end = int(l['end']); chrom = int(l['chrom'])
        mhc = chrom == 6 and start <= 36000000 and end >= 25000000
        lead_pos = int(l.get('lead_pos') or l['lead_variant'].split('_')[1])
        n_sig = int(l.get('significant_variants') or l.get('significant_rows') or l['n_significant_seeds'])
        all_loci.append(dict(track=t, locus_id=lid, build='GRCh37', chrom=chrom, start=start,
            end=end, interval='1_BASED_INCLUSIVE', width_bp=end-start+1,
            significant_seed_rows=n_sig, lead_source_row=l['lead_source_row'], lead_pos=lead_pos,
            lead_distance_left_bp=lead_pos-start, lead_distance_right_bp=end-lead_pos,
            mhc_deferred=mhc, chromosome_edge_clipped=l.get('edge_clipped', l.get('chromosome_edge_clipped')),
            rule='AUTOSOMAL_P_LT_5E-8_PLUS_MINUS_1.5MB_TRANSITIVE_MERGE',
            source_loci_table=f'tracks/{t}/results/loci.tsv'))
        o = overlap.get(t, {}).get(lid, {})
        p = s/f'ld/{t}/{lid}/numeric_diagnostics.json'
        n = json.loads(p.read_text()) if p.exists() else {}
        if n: completed += 1
        status = n.get('numerical_status', 'DEFERRED_MHC' if mhc else 'UNAVAILABLE_NO_LD' if t == 'B' else 'PENDING_DIAGNOSTICS')
        sig_missing = o.get('missing_or_excluded_significant_rows', '')
        numerical = {k: n.get(k, '') for k in ['variants', 'all_finite', 'original_lower_triangle_max_abs',
            'raw_diagonal_max_abs_deviation_from_one', 'derived_diagonal_max_abs_error',
            'symmetry_max_abs_error', 'symmetry_tolerance', 'derived_correlation_min',
            'derived_correlation_max', 'correlation_range_tolerance', 'min_eigenvalue', 'max_eigenvalue',
            'psd_tolerance', 'psd_pass', 'numerical_rank', 'condition_positive_subspace',
            'condition_warning_threshold', 'numerical_null_space_z_energy']}
        numerical.update(rank_deficient=(n['numerical_rank'] < n['variants']) if 'numerical_rank' in n else '',
            condition_warning=(n['condition_positive_subspace'] > n['condition_warning_threshold']) if n.get('condition_positive_subspace') else '',
            null_mismatch_s=n.get('null_covariance_mismatch_diagnostic', {}).get('s', ''),
            scientific_LD_gate_pass=False, source_diagnostic=str(p.relative_to(s)) if p.exists() else '')
        numeric.append(dict(track=t, locus_id=lid, numerical_status=status, **numerical))
        compatibility.append(dict(track=t, locus_id=lid, mhc_deferred=mhc,
            total_gwas_locus_rows=o.get('total_gwas_locus_rows', ''),
            eligible_signed_identity_rows=o.get('summary_row_eligible', ''),
            ordered_diagnostic_variants=o.get('retained_diagnostic_variants', ''),
            missing_eligible_variants=o.get('missing_eligible_variants', ''),
            significant_rows=n_sig, missing_or_excluded_significant_rows=sig_missing,
            lead_retained=o.get('lead_retained_for_diagnostics', bool(o) and sig_missing == 0),
            source_role=sources[t], sample_likelihood_contract='UNRESOLVED',
            diagnostic_path=f'ld/{t}/{lid}' if o else '',
            diagnostic_interpretation='DESCRIPTIVE_ONLY_NO_POSTERIOR_OR_COVARIANCE_REPAIR'))
        reasons = ('FROZEN_COMPLEX_MHC_DEFERRAL; ' if mhc else '') + blockers[t]
        if n and n.get('numerical_status') != 'PASS_NUMERICS_ONLY': reasons += '; NUMERICAL_MATRIX_GATE_FAILED'
        readiness.append(dict(track=t, locus_id=lid, final_state='NOT CLEARED',
            mhc_deferred=mhc, all_summary_statistic_gates_pass=False, all_LD_gates_pass=False,
            numerical_status=status, execution_cleared=False, blockers=reasons,
            future_input_package='NONE_NO_EXECUTION_CLEARANCE'))
        for gate, gate_status, evidence in [
            ('phenotype', 'PASS_WITH_TRACK_ROLE', roles[t]),
            ('full_file_integrity', 'PASS', f'tracks/{t}/'+('audit.json' if t=='A' else 'results/audit.json')),
            ('signed_statistic_and_method_contract', 'UNRESOLVED_METHOD_N_CONTRACT', 'Signed directions verified separately; see effect_scale_summary.tsv and source audit'),
            ('variant_identity', 'DEFERRED_MHC' if mhc else 'UNRESOLVED_STRONG_IDENTITIES' if t=='A' and sig_missing else 'VERIFIED_ROWS_WITH_EXPLICIT_EXCLUSIONS', f'tracks/{t}/results'),
            ('LD_origin', 'UNAVAILABLE' if t=='B' else 'SOURCE_IDENTIFIED_EXACT_EXECUTED_CONFIGURATION_UNBOUND', sources[t]),
            ('locus_density', 'DEFERRED_MHC' if mhc else 'UNAVAILABLE' if t=='B' else 'FAIL_STRONG_SIGNAL_LOSS' if sig_missing else 'NO_SIGNIFICANT_SIGNAL_LOSS; COMPLETE_SCIENTIFIC_DENSITY_REVIEW_OPEN', f'ld/{t}/gwas_ld_overlap_summary.json' if t!='B' else 'tracks/B/results/ld_source_contract.json'),
            ('LD_numerics', status, str(p.relative_to(s)) if p.exists() else f'tracks/{t}/results/loci.tsv' if mhc else 'tracks/B/results/ld_source_contract.json' if t=='B' else 'PENDING_DIAGNOSTICS'),
            ('sample_compatibility', 'UNRESOLVED', blockers[t]),
            ('summary_LD_consistency', 'UNRESOLVED' if not mhc else 'DEFERRED_MHC', 'Descriptive diagnostics do not override sampling/test contract'),
            ('cohort_overlap_documentation', 'PASS_DOCUMENTED', 'COHORT_OVERLAP_AND_SCOPE.md'),
            ('locus_boundary', 'DEFERRED_MHC' if mhc else 'FROZEN_WINDOWS_PASS; LONG_RANGE_STRUCTURE_NOT_EXCLUDED', 'tables/prospective_loci.tsv; ld/extraction_validation.json'),
            ('execution_input_freeze', 'NO_CLEARED_PACKAGE', 'execution_inputs/manifest.json')]:
            gates.append(dict(track=t,locus_id=lid,gate=gate,status=gate_status,evidence=evidence))

if final:
    assert completed == 16, ('Expected all 16 obtainable non-MHC diagnostic matrices', completed)
    assert not any(x['numerical_status']=='PENDING_DIAGNOSTICS' for x in numeric)
    validated=read('ld/extraction_validation.json')
    assert validated['mode']=='final' and validated['status']=='PASS_EXTRACTION_AND_OUTPUT_VALIDATION_ONLY'
    assert validated['planned_loci']==16 and validated['check_counts'].get('FAIL',0)==0
    assert validated['check_counts'].get('PENDING',0)==0
tsv('tables/prospective_loci.tsv', all_loci)
tsv('tables/locus_readiness.tsv', readiness)
tsv('tables/LD_numerical_QC.tsv', numeric)
tsv('tables/GWAS_LD_compatibility.tsv', compatibility)
tsv('tables/readiness_gates.tsv', gates)
tracks=[]
for t in 'ABC':
    full, total, nonmhc, mhc, risk = counts[t]
    tracks.append(dict(track=t,role=roles[t],full_GWAS_rows=full,prospective_loci_total=total,
        prospective_non_MHC_loci=nonmhc,deferred_MHC_loci=mhc,all_summary_statistic_gates_pass_loci=0,
        all_LD_gates_pass_loci=0,execution_cleared_loci=0,final_state='NOT CLEARED',
        diagnostic_matrices_completed=sum(x['track']==t and x['numerical_status'] not in ['DEFERRED_MHC','UNAVAILABLE_NO_LD','PENDING_DIAGNOSTICS'] for x in numeric),
        basic_numerics_pass_loci=sum(x['track']==t and x['numerical_status']=='PASS_NUMERICS_ONLY' for x in numeric),
        risk_direction_assigned_rows_including_deferred_MHC=risk,
        LD_source_or_blocker=sources[t],blockers=blockers[t],
        future_statistical_execution_authorized=False))
tsv('tables/track_readiness.tsv', tracks)
(s/'execution_inputs/manifest.json').write_text(json.dumps({
    'stage':s.name,'execution_cleared_loci':0,'packages':[],
    'reason':'No track closes all frozen summary-statistic/LD/sample gates.',
    'diagnostic_inputs_are_execution_inputs':False,
    'diagnostic_material_location':'ld/{A,C}/{locus_id}/ ordered signed inputs, raw triangles, diagnostic correlations, exclusion ledgers and diagnostics',
    'future_statistical_execution_authorized':False,'new_investigator_authorization_required':True},indent=2)+'\n')
(s/'tables/readiness_integration.json').write_text(json.dumps({
    'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'final':final,'prospective_loci_total':23,'non_MHC_loci_total':21,'deferred_MHC_loci':2,
    'diagnostic_matrices_completed':completed,'execution_cleared_loci':0,
    'track_counts':tracks,'no_new_gate_thresholds':True,'posterior_inference':False,
    'summary_gate_definition':'All method-ready statistic/identity/N contracts; valid signed direction alone is not the complete summary gate.',
    'LD_gate_definition':'All origin/density/numerics/sample/summary-covariance gates; basic numerical pass alone is not the complete LD gate.'},indent=2)+'\n')
print(json.dumps({'completed_diagnostics':completed,'final':final,'cleared_loci':0}))
