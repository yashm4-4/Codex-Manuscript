#!/usr/bin/env python3
"""Freeze reviewed local pretraining-1.1 artifacts; never start an experiment."""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import subprocess
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[5]
BASE = ROOT/'diseases/COPD/07_gap_closure/pretraining-1.1'


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(4*1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(path.read_text())


def tsv(path, rows):
    with path.open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader(); writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--validation-label', required=True)
    args = parser.parse_args()
    if not args.validation_label.replace('_','').replace('-','').isalnum():
        raise ValueError('Invalid validation label')
    freeze = BASE/'provenance/freeze.json'
    checksums = BASE/'provenance/artifact_checksums.tsv'
    inventory = BASE/'provenance/artifact_inventory.tsv'
    if any(path.exists() for path in [freeze, checksums, inventory]):
        raise RuntimeError('Freeze already exists or incomplete prior freeze; refusing overwrite')
    validation_path = BASE/f'provenance/validation_{args.validation_label}_manifest.json'
    validation = read(validation_path)
    if validation['artifact_failures'] or validation['readiness_failures']:
        raise RuntimeError('Cannot declare ready: independent validation failures remain')
    for item in validation['outputs']:
        assert sha(ROOT/item['path']) == item['sha256'], item['path']
    assert sha(BASE/'scripts/validate_preflight11.py') == validation['validator_sha256']
    ancillary = read(BASE/'provenance/final_packaging_validation.json')
    assert ancillary['status'] == 'PASS' and ancillary['failures'] == 0
    assert sha(BASE/'scripts/validate_packaging.py') == ancillary['script_sha256']
    for path, item in ancillary['input_hashes'].items():
        assert sha(ROOT/path) == item['sha256'], path
    pair_audit = read(BASE/'provenance/supplementary_pair_metadata_validation.json')
    assert pair_audit['status'] == 'PASS' and pair_audit['failures'] == 0
    assert sha(BASE/'scripts/validate_pair_metadata.py') == pair_audit['script_sha256']
    for path, item in pair_audit['input_hashes'].items():
        assert sha(ROOT/path) == item['sha256'], path
    cfg = read(BASE/'provenance/configuration_manifest.json')
    for item in cfg['inputs']:
        assert sha(ROOT/item['path']) == item['sha256'], item['path']
    report_review = read(BASE/'provenance/report_final_review.json')
    assert report_review['status'] == 'PASS' and report_review['failures'] == 0
    assert sha(BASE/'MODEL_REDESIGN_PREFLIGHT.md') == report_review['report_sha256']
    attempt = read(BASE/'attempts/001_full_population/attempt_manifest.json')
    assert attempt['balance_failures'] == 0 and attempt['positive_retention'] == 1.
    assert sha(BASE/'scripts/construct_matched_design.py') == attempt['implementation_sha256']
    assert sha(BASE/'specification/matching_attempt_001.json') == attempt['rule_sha256']
    for key in ['training_started','model_inference_performed','benchmark_outcomes_read']:
        assert attempt[key] is False
    assert validation['training_authorized'] is False
    assert validation['GPU_training_or_model_inference_performed'] is False
    assert validation['phase_I_extraction_performed'] is False
    assert validation['benchmark_outcomes_parsed'] is False
    with (BASE/'provenance/historical_tracked_state.tsv').open() as handle:
        protected = list(csv.DictReader(handle, delimiter='\t'))
    for item in protected:
        assert sha(ROOT/item['path']) == item['sha256'], item['path']
    with (BASE.parent/'provenance/COPD-V2-PREFLIGHT_consumed_input_hashes.tsv').open() as handle:
        consumed = list(csv.DictReader(handle, delimiter='\t'))
    for item in consumed:
        assert sha(ROOT/item['path']) == item['sha256'], item['path']
    assert not subprocess.check_output(['git','diff','--name-only','HEAD'],cwd=ROOT).strip()
    readiness = [
        {'condition': 'B_prespecified_balance_gates', 'scope': 'construction_review', 'status': 'SATISFIED', 'observed': '0 failures', 'required': '|SMD|<=.10; chromosome gap<=.02; signature and each-lobe gap<=.05'},
        {'condition': 'C_prespecified_balance_gates', 'scope': 'construction_review', 'status': 'SATISFIED', 'observed': '0 failures with independent C controls', 'required': '|SMD|<=.10; chromosome gap<=.02; signature and each-lobe gap<=.05'},
        {'condition': 'B_positive_retention', 'scope': 'construction_review', 'status': 'SATISFIED', 'observed': '100% all V1 positives; unique 1:1 controls; 0 exclusions', 'required': 'prefer all; >=90% if restriction necessary'},
        {'condition': 'C_positive_retention', 'scope': 'construction_review', 'status': 'SATISFIED', 'observed': '100% same-lobe positives; unique 1:1 controls; 0 exclusions', 'required': 'prefer all; >=90% if restriction necessary'},
        {'condition': 'contamination_and_source_eligibility', 'scope': 'construction_review', 'status': 'SATISFIED', 'observed': 'independent raw-source audit passes', 'required': 'donor ATAC; no relevant full-input mark or known same-model positive; preserved exclusions'},
        {'condition': 'sequence_and_duplicate_leakage', 'scope': 'construction_review', 'status': 'SATISFIED', 'observed': '792681 selected intervals verified; 0 cross-partition/role leakage', 'required': 'available 2001bp sequences; fixed chromosome/role boundary integrity'},
        {'condition': 'C_common_selection_panels', 'scope': 'construction_review', 'status': 'SATISFIED', 'observed': 'enhancer 1755+/1755-/1295 components; H3K27me3 271+/271-/244 components', 'required': '>=100 positives; >=200 controls; >=30 components per context'},
        {'condition': 'C_calibration_control_panels', 'scope': 'construction_review', 'status': 'SATISFIED', 'observed': 'enhancer 1576 controls/809 control components; H3K27me3 263/176', 'required': '>=200 controls and >=30 control components per context'},
        {'condition': 'training_only_ATAC_normalization', 'scope': 'construction_review', 'status': 'SATISFIED', 'observed': 'source ECDFs independently reconstructed from training chromosomes only', 'required': 'frozen ties/outside-range mapping unchanged on held-out rows'},
        {'condition': 'global_tie_and_exchange_contract', 'scope': 'construction_review', 'status': 'SATISFIED', 'observed': '377118 initial assignments and 1371 exchanges independently replayed', 'required': 'complete deterministic global ties; exact final mapping'},
        {'condition': 'independent_validation', 'scope': 'construction_review', 'status': 'SATISFIED', 'observed': f"{validation['artifact_checks']} core checks plus {len(ancillary['checks'])} packaging checks; 0 failures", 'required': 'all scientific and artifact checks pass'},
        {'condition': 'historical_1_0_and_benchmark_preservation', 'scope': 'construction_review', 'status': 'SATISFIED', 'observed': f'{len(protected)} tracked historical byte hashes unchanged; benchmark opaque hashing only', 'required': 'no edits or benchmark outcome access'},
        {'condition': 'all_specifications_and_payloads_frozen', 'scope': 'construction_review', 'status': 'SATISFIED_AT_THIS_FREEZE', 'observed': 'covered by artifact_checksums.tsv and freeze.json', 'required': 'immutable local version before V2 inference'},
        {'condition': 'no_model_or_GPU_execution', 'scope': 'construction_review', 'status': 'SATISFIED', 'observed': 'none', 'required': 'no training, phase-I extraction, inference or scoring'},
        {'condition': 'investigator_training_authorization', 'scope': 'execution', 'status': 'NOT_GRANTED_STOP', 'observed': 'review-only CPU repair; no GPU permission', 'required': 'separate explicit investigator authorization'},
        {'condition': 'future_C_model_absolute_adequacy', 'scope': 'post_training_release', 'status': 'NOT_EVALUATED', 'observed': 'no predictions generated', 'required': 'both contexts pass all fixed adequacy/stability/invariance/calibration gates; no A/B fallback'}
    ]
    tsv(BASE/'provenance/training_readiness.tsv', readiness)
    paths = sorted(path for path in BASE.rglob('*') if path.is_file() and '__pycache__' not in path.parts and path.suffix != '.pyc')
    rows = [{'path': str(path.relative_to(ROOT)), 'bytes': path.stat().st_size, 'sha256': sha(path)} for path in paths]
    tsv(inventory, rows)
    rows.append({'path': str(inventory.relative_to(ROOT)), 'bytes': inventory.stat().st_size, 'sha256': sha(inventory)})
    rows.sort(key=lambda row: row['path'])
    tsv(checksums, rows)
    for item in rows:
        assert sha(ROOT/item['path']) == item['sha256'], item['path']
    summary = {
        'version': 'pretraining-1.1', 'module': 'COPD-V2-PREFLIGHT',
        'status': 'READY_FOR_INVESTIGATOR_TRAINING_REVIEW',
        'status_text': 'READY FOR INVESTIGATOR TRAINING REVIEW',
        'training_authorized': False, 'training_started': False,
        'GPU_training_or_model_inference_performed': False, 'phase_I_extraction_performed': False,
        'benchmark_outcomes_parsed': False, 'historical_pretraining_1_0_changed': False,
        'baseline_commit': '71f7b3105184d7e4425aa39a6d1d441b20bdc473',
        'frozen_utc': datetime.now(timezone.utc).isoformat(),
        'validation_manifest': str(validation_path.relative_to(ROOT)),
        'validation_manifest_sha256': sha(validation_path),
        'artifact_checks': validation['artifact_checks'], 'artifact_failures': 0,
        'packaging_checks': len(ancillary['checks']), 'packaging_failures': 0,
        'pair_metadata_checks': len(pair_audit['checks']), 'pair_metadata_failures': 0,
        'readiness_failures': 0, 'B_C_balance_failures': 0,
        'B_C_positive_retention': 1.0, 'positive_common_support_exclusions': 0,
        'protected_historical_tracked_files_verified': len(protected),
        'protected_original_consumed_inputs_verified': len(consumed),
        'artifact_checksum_ledger': str(checksums.relative_to(ROOT)),
        'artifact_checksum_ledger_sha256': sha(checksums),
        'payload_count': len(rows), 'payload_bytes': sum(row['bytes'] for row in rows),
        'checksum_self_exclusions': ['artifact checksum ledger itself', 'freeze.json itself'],
        'report': str((BASE/'MODEL_REDESIGN_PREFLIGHT.md').relative_to(ROOT)),
        'report_sha256': sha(BASE/'MODEL_REDESIGN_PREFLIGHT.md'),
        'scientific_interpretation': 'B/C construction gates pass using full-positive distributional matching and independent C controls; close-to-boundary balance and removed pair calipers remain explicit limitations. No model adequacy or performance is established.',
        'final_configuration_rule': 'C is the intended corrected model subject to absolute adequacy/stability/invariance/calibration in both contexts after separately authorized training; A/B are diagnostic ablations, never automatic fallback.',
        'next_authority_required': 'Investigator review of this freeze and separate explicit GPU execution authorization. Stop now; do not launch 18 fits.',
        'publication': 'Local artifacts only; no commit or push authorized or performed in this repair.'
    }
    with freeze.open('x') as handle:
        json.dump(summary, handle, indent=2, sort_keys=True); handle.write('\n')
    print(json.dumps(summary, indent=2))
    print('Freeze SHA256:', sha(freeze))


if __name__ == '__main__':
    main()
