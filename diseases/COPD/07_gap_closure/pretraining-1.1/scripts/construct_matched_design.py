#!/usr/bin/env python3
"""Outcome-blind, CPU-only full-population distribution matching for pretraining-1.1.

No models, scores, benchmark tables, or accelerator frameworks are loaded.
Every output is exclusive-create; failed attempts are retained, never overwritten.
"""
from __future__ import annotations

import os
os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
import csv
import hashlib
import json
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[5]
BASE = ROOT / 'diseases/COPD/07_gap_closure'
VERSION = BASE / 'pretraining-1.1'
OLD = BASE / 'data/COPD-V2-PREFLIGHT'
OUT = VERSION / 'attempts/001_full_population'
FEATURES = ['gc_fraction', 'atac_signal_percentile_max_train_only', 'repeat_fraction_2001', 'blacklist_input_any', 'non_acgt_any']
DISTANCE_SCALES = np.array([.05, .20, .20])
SEEDS = [104729, 130363, 155921]
CHROMS = {f'chr{i}': i for i in range(1, 23)} | {'chrX': 23, 'chrY': 24}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def write(path, frame):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        raise RuntimeError(f'Refusing overwrite: {path}')
    frame.to_csv(path, sep='\t', index=False, header=not path.name.endswith('.bed.gz'),
                 compression={'method': 'gzip', 'mtime': 0} if path.suffix == '.gz' else None,
                 float_format='%.17g')


def json_out(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')


def specification():
    return {
        'version': 'pretraining-1.1', 'attempt_id': '001_full_population',
        'recorded_before_execution_utc': datetime.now(timezone.utc).isoformat(),
        'purpose': 'Prefer all original B positives and all same-lobe C positives with unique 1:1 controls; distributional balance, not arbitrary close pairs, is the scientific requirement.',
        'positive_population': {'V2-B': 'all exact V1 positives', 'V2-C': 'all V1 positives with same-lobe original ATAC-peak/relevant-mark support'},
        'control_eligibility': 'Unchanged frozen 1.0 donor-ATAC/full-input-mark/known-positive/core-blacklist/core-TSS/sequence eligibility. No duplicated controls within a configuration/model.',
        'partition_role_contract': 'Unchanged full master 1.0 genomic-overlap/encoded-RC components and chr7 role assignments. Controls never cross a target partition or chr7 role.',
        'positive_order': 'Ascending SHA256(271828|model|configuration|interval_id), with canonical coordinate rank as secondary tie.',
        'canonical_coordinate_order': 'chromosome ordinal chr1..chr22,chrX,chrY; numeric core_start; numeric core_end; interval_id lexicographic',
        'initial_candidate_priority': ['same validation_role + ATAC lobe signature + chromosome', 'same validation_role + ATAC lobe signature, any chromosome', 'same validation_role, any chromosome/signature'],
        'initial_distance': 'float64 Chebyshev on GC/.05, training-only ATAC midpoint-ECDF/.20, repeat/.20; these are distance scales, NOT hard pair calipers',
        'global_nearest_tie': 'Find an unused nearest candidate, then query the COMPLETE radius d+1e-12 (nextafter outward), recompute exact float64 distances for ALL unused candidates in that radius, and minimize (distance, canonical coordinate rank). Never choose from a bounded-k tie subset.',
        'old_pair_calipers': 'Removed in favor of the unchanged prespecified distributional balance gates. No inference of pairwise biological exchangeability; distances and cross-chromosome/signature matches remain auditable.',
        'B_C': 'Independently constructed controls. B-to-C is same-lobe label correction plus induced control rematching, not a pure label-only experiment.',
        'exchange': {
            'scope': 'Each validation_role independently; swaps preserve control counts within exact chromosome x ATAC lobe signature.',
            'variables': FEATURES,
            'scale': 'positive sample SD / sqrt(2); if zero, use same-role eligible-pool sample SD / sqrt(2). If both zero and equal means, assign unit scale to the identically-zero coordinate; unequal constants reject the attempt. This numerical loss scale never replaces actual SMD gates.',
            'loss': 'sum of squared standardized control-minus-positive mean residuals',
            'round': 'Freeze residual gradient. Within every exchange group, globally order selected controls by descending gradient projection (coordinate tie) and unused controls by ascending projection (coordinate tie). Pair these two complete ordered lists by rank up to the shorter length. Globally order all resulting disjoint proposals by (linear loss derivative, removed coordinate rank, added coordinate rank). This explicit proposal family is a heuristic, not a claim of globally optimal Cartesian swaps.',
            'accept': 'Scan the globally ordered proposals and accept a swap only if actual recomputed quadratic loss strictly decreases by >1e-15; accept at most 256 per round, then refresh the gradient. Stop immediately when all five actual role-specific |SMD|<=.10.',
            'max_rounds': 200, 'max_accepted_per_round': 256,
            'failure': 'Keep the failed attempt. No silent parameter/gate modification or model-based tuning.',
            'tie_rule': 'Every projection/derivative tie uses the full global canonical coordinate order; candidate ranking is never truncated before tie resolution.'
        },
        'balance_gates': {'absolute_SMD': .10, 'absolute_chromosome_proportion_gap': .02, 'absolute_ATAC_lobe_signature_proportion_gap': .05, 'absolute_each_lobe_support_proportion_gap': .05},
        'gating_scopes': 'train, whole chr7 validation, chr8-9 test; chr7 roles additionally reported and optimized, not substitutes for whole-partition gates',
        'support': 'This attempt retains 100% of B/C target positives. A separately specified later common-support attempt is permitted only after rejection of full-population construction; minimum retention remains .90.',
        'test_history': 'chr8-9 is a fixed internal V2 training/selection holdout, not historically untouched; only construction/QC covariates may be accessed here.',
        'firewall': {'model_training': False, 'phase_I_extraction': False, 'model_inference': False, 'benchmark_outcome_access': False, 'GWAS_or_functional_priority': False},
        'input_sha256': {'old_features': sha(OLD/'interval_features.tsv.gz'), 'old_roles': sha(OLD/'interval_role_assignment.tsv.gz'), 'train_only_ATAC_sidecar': sha(VERSION/'data/atac_normalized_features.tsv.gz')},
        'implementation_sha256': sha(__file__)
    }


class AvailableNearest:
    def __init__(self, indices, points, used):
        self.indices = np.sort(np.asarray(indices, dtype=np.int64))
        self.points = points[self.indices]
        self.tree = cKDTree(self.points)
        self.used = used

    def choose(self, point):
        n = len(self.indices)
        if not n or self.used[self.indices].all():
            return None
        k = min(16, n)
        while True:
            distances, local = self.tree.query(point, k=k, p=np.inf)
            local = np.atleast_1d(local)
            available = local[~self.used[self.indices[local]]]
            if len(available):
                exact = np.max(np.abs(self.points[available] - point), axis=1)
                best_distance = float(exact.min())
                radius = np.nextafter(best_distance + 1e-12, np.inf)
                all_local = np.asarray(self.tree.query_ball_point(point, radius, p=np.inf), dtype=np.int64)
                all_local = all_local[~self.used[self.indices[all_local]]]
                dd = np.max(np.abs(self.points[all_local] - point), axis=1)
                ranks = self.indices[all_local]
                winner = np.lexsort((ranks, dd))[0]
                return int(ranks[winner]), float(dd[winner]), int(np.sum(dd == dd[winner]))
            if k == n:
                raise AssertionError('Unused candidate disappeared')
            k = min(n, 4 * k)


def initial_match(frame, model, config, pos, controls):
    points = frame[FEATURES[:3]].to_numpy(float) / DISTANCE_SCALES
    assert np.isfinite(points[np.r_[pos, controls]]).all()
    used = np.zeros(len(frame), dtype=bool)
    groups = [defaultdict(list), defaultdict(list), defaultdict(list)]
    role = frame.validation_role.to_numpy()
    lobes = frame.atac_lobes.to_numpy()
    chrom = frame.chrom.to_numpy()
    for i in controls:
        groups[0][(role[i], lobes[i], chrom[i])].append(i)
        groups[1][(role[i], lobes[i])].append(i)
        groups[2][(role[i],)].append(i)
    trees = [{}, {}, {}]
    ids = frame.interval_id.to_numpy()
    order = sorted(pos, key=lambda i: (hashlib.sha256(f'271828|{model}|{config}|{ids[i]}'.encode()).hexdigest(), i))
    rows = []
    for number, pi in enumerate(order):
        keys = [(role[pi], lobes[pi], chrom[pi]), (role[pi], lobes[pi]), (role[pi],)]
        for level, key in enumerate(keys):
            if key not in groups[level]:
                continue
            if key not in trees[level]:
                trees[level][key] = AvailableNearest(groups[level][key], points, used)
            result = trees[level][key].choose(points[pi])
            if result is not None:
                ci, distance, ties = result
                used[ci] = True
                rows.append({'configuration': config, 'model': model, 'processing_order': number,
                             'positive_id': ids[pi], 'control_id': ids[ci], 'positive_index': int(pi),
                             'control_index': ci, 'validation_role': role[pi], 'candidate_priority_level': level,
                             'initial_normalized_chebyshev_distance': distance, 'global_min_distance_tie_count': ties})
                break
        else:
            raise RuntimeError(f'Insufficient role-compatible eligible control pool: {model} {config} {ids[pi]}')
        if (number + 1) % 20000 == 0:
            print(f'{model} {config}: initial matches {number+1}/{len(order)}', flush=True)
    return pd.DataFrame(rows)


def smds(a, b):
    delta = a.mean(axis=0) - b.mean(axis=0)
    denom = np.sqrt((a.var(axis=0, ddof=1) + b.var(axis=0, ddof=1)) / 2)
    result = np.zeros_like(delta)
    np.divide(delta, denom, out=result, where=denom > 0)
    result[(denom == 0) & (delta != 0)] = np.inf
    return result


def exchange_balance(frame, model, config, pairs, eligible):
    values = frame[FEATURES].to_numpy(float)
    ids = frame.interval_id.to_numpy()
    role = frame.validation_role.to_numpy()
    selected = np.zeros(len(frame), dtype=bool)
    selected[pairs.control_index.to_numpy()] = True
    control_to_pair = {int(c): int(k) for k, c in zip(pairs.index, pairs.control_index)}
    swap_rows, round_rows = [], []
    for r in sorted(pairs.validation_role.unique()):
        pr = pairs[pairs.validation_role.eq(r)]
        positive = pr.positive_index.to_numpy(int)
        candidates = np.sort(eligible[role[eligible] == r])
        n = len(positive)
        pv = values[positive]
        scale = pv.std(axis=0, ddof=1) / np.sqrt(2)
        pool_sd = values[candidates].std(axis=0, ddof=1) / np.sqrt(2)
        scale = np.where(scale == 0, pool_sd, scale)
        constant = scale == 0
        if np.any(constant & (pv.mean(axis=0) != values[candidates].mean(axis=0))):
            raise RuntimeError(f'Unsatisfiable constant-feature distribution: {model} {config} {r}')
        scale[constant] = 1.
        normalized = values[candidates] / scale
        target = pv.mean(axis=0) / scale
        group_map = defaultdict(list)
        for local, ci in enumerate(candidates):
            group_map[(frame.at[ci, 'chrom'], frame.at[ci, 'atac_lobes'])].append(local)
        groups = [np.array(group_map[key], dtype=np.int64) for key in sorted(group_map)]
        for round_number in range(201):
            active = selected[candidates]
            residual = normalized[active].mean(axis=0) - target
            loss = float(residual @ residual)
            observed_smd = smds(pv, values[candidates[active]])
            row = {'configuration': config, 'model': model, 'validation_role': r, 'round': round_number,
                   'loss_before': loss, 'maximum_absolute_SMD_before': float(np.max(np.abs(observed_smd))),
                   'proposals': 0, 'accepted': 0, 'rejected_nonimproving': 0, 'status': ''}
            if np.all(np.abs(observed_smd) <= .10):
                row['status'] = 'ALL_ROLE_SMD_GATES_PASS'; round_rows.append(row); break
            if round_number == 200:
                row['status'] = 'STOP_MAX_PRESPECIFIED_ROUNDS'; round_rows.append(row); break
            # Full projections and rankings, without bounded candidate truncation.
            projection = normalized @ (2 * residual / n)
            removes, adds = [], []
            for group in groups:
                old = group[active[group]]; new = group[~active[group]]
                count = min(len(old), len(new))
                if not count:
                    continue
                old = old[np.lexsort((candidates[old], -projection[old]))]
                new = new[np.lexsort((candidates[new], projection[new]))]
                removes.extend(old[:count]); adds.extend(new[:count])
            if not removes:
                row['status'] = 'STOP_NO_EXCHANGE'; round_rows.append(row); break
            removes, adds = np.asarray(removes), np.asarray(adds)
            linear = projection[adds] - projection[removes]
            proposal_order = np.lexsort((candidates[adds], candidates[removes], linear))
            row['proposals'] = len(proposal_order)
            accepted = 0
            selected_values = values[candidates[active]]
            control_sum = selected_values.sum(axis=0)
            control_sum_squares = np.square(selected_values).sum(axis=0)
            positive_mean = pv.mean(axis=0)
            positive_variance = pv.var(axis=0, ddof=1)
            for proposal_rank, ii in enumerate(proposal_order):
                old_local, new_local = removes[ii], adds[ii]
                old, new = int(candidates[old_local]), int(candidates[new_local])
                delta = (normalized[new_local] - normalized[old_local]) / n
                after = residual + delta
                next_loss = float(after @ after)
                if not next_loss < loss - 1e-15:
                    row['rejected_nonimproving'] += 1
                    continue
                pair_index = control_to_pair.pop(old)
                control_to_pair[new] = pair_index
                pairs.at[pair_index, 'control_index'] = new
                pairs.at[pair_index, 'control_id'] = ids[new]
                selected[old] = False; selected[new] = True
                swap_rows.append({'configuration': config, 'model': model, 'validation_role': r,
                                  'round': round_number, 'proposal_rank': proposal_rank, 'swap_in_round': accepted,
                                  'positive_id': pairs.at[pair_index, 'positive_id'], 'removed_control_id': ids[old],
                                  'added_control_id': ids[new], 'removed_control_index': old, 'added_control_index': new,
                                  'frozen_gradient_derivative': float(linear[ii]), 'loss_before': loss, 'loss_after': next_loss})
                residual, loss = after, next_loss
                control_sum += values[new] - values[old]
                control_sum_squares += np.square(values[new]) - np.square(values[old])
                accepted += 1
                # Check using actual selected sample variance, not the surrogate objective.
                control_mean = control_sum / n
                control_variance = np.maximum(0, (control_sum_squares - control_sum**2 / n) / (n-1))
                denominator = np.sqrt((positive_variance + control_variance) / 2)
                provisional = np.zeros(len(FEATURES))
                np.divide(positive_mean-control_mean, denominator, out=provisional, where=denominator>0)
                provisional[(denominator==0)&(positive_mean!=control_mean)] = np.inf
                if np.all(np.abs(provisional) <= .10 + 1e-12):
                    if np.all(np.abs(smds(pv, values[candidates[selected[candidates]]])) <= .10):
                        break
                if accepted >= 256:
                    break
            row['accepted'] = accepted
            row['status'] = 'REFRESH_GRADIENT' if accepted else 'STOP_NO_IMPROVING_PROPOSAL'
            round_rows.append(row)
            print(f'{model} {config} {r}: exchange round {round_number}, {accepted} swaps, loss {loss:.8g}', flush=True)
            if not accepted:
                break
    return pairs, swap_rows, round_rows


def diagnostics(frame, model, config, positives, controls, stage):
    rows = []
    for scope_type, scopes in [('partition', ['train', 'validation', 'test']), ('validation_role', ['checkpoint', 'selection', 'calibration'])]:
        for scope in scopes:
            pi = positives[frame.loc[positives, scope_type].to_numpy() == scope]
            ci = controls[frame.loc[controls, scope_type].to_numpy() == scope]
            a, b = frame.loc[pi], frame.loc[ci]
            if not len(a) or not len(b):
                raise AssertionError(f'Empty diagnostic population: {model} {config} {scope}')
            for variable, smd in zip(FEATURES, smds(a[FEATURES].to_numpy(float), b[FEATURES].to_numpy(float))):
                rows.append({'stage': stage, 'configuration': config, 'model': model, 'scope_type': scope_type,
                             'scope': scope, 'variable': variable, 'kind': 'SMD', 'positive_n': len(a), 'control_n': len(b),
                             'positive_mean': a[variable].mean(), 'control_mean': b[variable].mean(), 'value': smd,
                             'limit': .10, 'status': 'PASS' if np.isfinite(smd) and abs(smd) <= .10 else 'FAIL'})
            for variable, limit in [('chrom', .02), ('atac_lobes', .05)]:
                ap, bp = a[variable].value_counts(normalize=True), b[variable].value_counts(normalize=True)
                for level in sorted(set(ap.index) | set(bp.index)):
                    value = float(ap.get(level, 0) - bp.get(level, 0))
                    rows.append({'stage': stage, 'configuration': config, 'model': model, 'scope_type': scope_type,
                                 'scope': scope, 'variable': variable+'='+str(level), 'kind': 'proportion_gap',
                                 'positive_n': len(a), 'control_n': len(b), 'positive_mean': ap.get(level, 0),
                                 'control_mean': bp.get(level, 0), 'value': value, 'limit': limit,
                                 'status': 'PASS' if abs(value) <= limit else 'FAIL'})
            for lobe in ['lower_left', 'lower_right', 'upper_right']:
                ap = a.atac_lobes.str.split(';').apply(lambda z: lobe in z).mean()
                bp = b.atac_lobes.str.split(';').apply(lambda z: lobe in z).mean()
                rows.append({'stage': stage, 'configuration': config, 'model': model, 'scope_type': scope_type,
                             'scope': scope, 'variable': 'lobe_support='+lobe, 'kind': 'proportion_gap',
                             'positive_n': len(a), 'control_n': len(b), 'positive_mean': ap, 'control_mean': bp,
                             'value': ap-bp, 'limit': .05, 'status': 'PASS' if abs(ap-bp) <= .05 else 'FAIL'})
    return rows


def main():
    started = time.time()
    if OUT.exists() or (VERSION/'provenance/freeze.json').exists():
        raise RuntimeError('Attempt already exists or version frozen')
    rule = specification()
    json_out(VERSION/'specification/matching_attempt_001.json', rule)
    OUT.mkdir(parents=True)
    frame = pd.read_csv(OLD/'interval_features.tsv.gz', sep='\t', keep_default_na=False)
    roles = pd.read_csv(OLD/'interval_role_assignment.tsv.gz', sep='\t', usecols=['interval_id', 'component_id', 'validation_role'])
    rank = pd.read_csv(VERSION/'data/atac_normalized_features.tsv.gz', sep='\t', usecols=['interval_id', FEATURES[1]])
    frame = frame.merge(roles, on='interval_id', validate='one_to_one').merge(rank, on='interval_id', validate='one_to_one')
    for feature in FEATURES[:3]:
        frame[feature] = pd.to_numeric(frame[feature], errors='coerce')
    frame['blacklist_input_any'] = (frame.blacklist_bp_2001 > 0).astype(int)
    frame['non_acgt_any'] = (frame.non_acgt_fraction > 0).astype(int)
    frame['_chrom_order'] = frame.chrom.map(CHROMS)
    frame = frame.sort_values(['_chrom_order', 'core_start', 'core_end', 'interval_id']).reset_index(drop=True)
    initial_all, final_all, swap_all, round_all, diagnostic_all, support = [], [], [], [], [], []
    for model in ['enhancer', 'h3k27me3']:
        eligible = frame.index[frame[f'{model}_accessible_control_eligible'].eq(1)].to_numpy(int)
        for config in ['V2-B', 'V2-C']:
            mask = frame[f'v1_{model}_positive'].eq(1)
            if config == 'V2-C':
                mask &= frame[f'{model}_same_lobe_peak_support'].eq(1)
            positives = frame.index[mask].to_numpy(int)
            initial = initial_match(frame, model, config, positives, eligible)
            initial_all.append(initial.copy())
            diagnostic_all.extend(diagnostics(frame, model, config, positives, initial.control_index.to_numpy(int), 'initial'))
            final, swaps, rounds = exchange_balance(frame, model, config, initial.copy(), eligible)
            final['final_normalized_chebyshev_distance'] = np.max(np.abs(frame.loc[final.positive_index, FEATURES[:3]].to_numpy(float) / DISTANCE_SCALES - frame.loc[final.control_index, FEATURES[:3]].to_numpy(float) / DISTANCE_SCALES), axis=1)
            final_all.append(final); swap_all.extend(swaps); round_all.extend(rounds)
            diagnostic_all.extend(diagnostics(frame, model, config, positives, final.control_index.to_numpy(int), 'final'))
            for partition in ['train', 'validation', 'test']:
                n = int((frame.loc[positives, 'partition'] == partition).sum())
                support.append({'configuration': config, 'model': model, 'partition': partition, 'original_positive_n': n, 'retained_positive_n': n, 'unique_control_n': n, 'excluded_positive_n': 0, 'retention': 1., 'B_represents_complete_V1_positive_population': config == 'V2-B'})
            write(OUT/f'{config}_{model}_initial_pairs.tsv.gz', initial)
            write(OUT/f'{config}_{model}_final_pairs.tsv.gz', final)
            print(f'Completed {model} {config}: {len(positives)} positives, {len(final)} unique controls, {len(swaps)} exchanges', flush=True)
    write(OUT/'initial_pairs.tsv.gz', pd.concat(initial_all, ignore_index=True))
    write(OUT/'control_matching_pairs.tsv.gz', pd.concat(final_all, ignore_index=True))
    write(OUT/'exchange_ledger.tsv.gz', pd.DataFrame(swap_all, columns=['configuration','model','validation_role','round','proposal_rank','swap_in_round','positive_id','removed_control_id','added_control_id','removed_control_index','added_control_index','frozen_gradient_derivative','loss_before','loss_after']))
    write(OUT/'exchange_rounds.tsv', pd.DataFrame(round_all))
    diag = pd.DataFrame(diagnostic_all)
    write(OUT/'covariate_balance.tsv', diag)
    write(OUT/'positive_retention.tsv', pd.DataFrame(support))
    write(OUT/'positive_common_support_exclusions.tsv.gz', pd.DataFrame(columns=['configuration','model','positive_id','partition','reason']))
    gates = diag[(diag.stage == 'final') & (diag.scope_type == 'partition')]
    summary = {'attempt_id': '001_full_population', 'positive_retention': 1., 'unique_one_to_one_controls': True,
               'balance_failures': int((gates.status == 'FAIL').sum()), 'additional_chr7_role_failures': int(((diag.stage == 'final') & (diag.scope_type == 'validation_role') & (diag.status == 'FAIL')).sum()),
               'status': 'CANDIDATE_CONSTRUCTION_PASS_PENDING_INDEPENDENT_VALIDATION' if (gates.status == 'PASS').all() else 'REJECTED_BALANCE_GATES',
               'elapsed_seconds': time.time()-started, 'rule_sha256': sha(VERSION/'specification/matching_attempt_001.json'),
               'implementation_sha256': sha(__file__), 'training_started': False, 'model_inference_performed': False, 'benchmark_outcomes_read': False}
    json_out(OUT/'attempt_manifest.json', summary)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == '__main__':
    main()
