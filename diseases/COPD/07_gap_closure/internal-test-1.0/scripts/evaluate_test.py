#!/usr/bin/env python3
"""One-time, hash-bound frozen chr8/9 common-C-task evaluation; no inference.

All evaluation inputs must already be frozen. This script cannot select models,
change cutpoints, fit calibration, or read candidate/external benchmark data.
Metric arithmetic is reused from the immutable, SHA-bound chr7 implementation.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import numpy as np

sys.dont_write_bytecode = True
FROZEN_EVALUATOR_SHA256 = "55b916dfd94c432b8463b0b6232792222f2fc1b020f060fec2edc4ff11d43e83"
FROZEN_EVALUATOR = Path(__file__).resolve().parents[2] / "internal-training-1.0/scripts/evaluate_chr7.py"
if hashlib.sha256(FROZEN_EVALUATOR.read_bytes()).hexdigest() != FROZEN_EVALUATOR_SHA256:
    raise RuntimeError("Immutable metric implementation hash mismatch")
_module_spec = importlib.util.spec_from_file_location("frozen_chr7_metrics_read_only", FROZEN_EVALUATOR)
ev = importlib.util.module_from_spec(_module_spec)
_module_spec.loader.exec_module(ev)

SEEDS, CONFIGURATIONS, MODELS = ev.SEEDS, ev.CONFIGURATIONS, ev.MODELS
METRICS, COVARIATES = ev.METRICS, ev.COVARIATES
EXPECTED = {"enhancer": {"positive": 12695, "control": 12695, "components": 8754},
            "h3k27me3": {"positive": 1754, "control": 1754, "components": 1809}}
THRESHOLDS = {"enhancer": "0x1.7f39710000001p-1", "h3k27me3": "0x1.8a0a12aaaaaacp-1"}
RATE_NAMES = ("sensitivity", "specificity", "precision", "NPV", "accuracy", "negative_row_FPR")


def validate_test_rows(rows, model, expected=None):
    expected = EXPECTED[model] if expected is None else expected
    ids, components = [], []
    for row in rows:
        if row["chrom"] not in ("chr8", "chr9") or row["partition"] != "test":
            raise ValueError("Forbidden chromosome or partition in common test panel")
        if row["validation_role"] != "test":
            raise ValueError("Unexpected common test role")
        if row["configuration"] != "V2-C" or row["model"] != model:
            raise ValueError("Common panel must retain frozen C-task identity")
        if row["label"] not in ("0", "1", 0, 1):
            raise ValueError("Nonbinary panel label")
        if not row["component_id"] or not row["interval_id"]:
            raise ValueError("Missing component/interval identity")
        if str(row["sequence_available"]) != "1" or int(row["sequence_length"]) != 2001:
            raise ValueError("Frozen sequence eligibility mismatch")
        if int(row["input_end"]) - int(row["input_start"]) != 2001:
            raise ValueError("Frozen input geometry mismatch")
        for key in ("sequence_sha256", "canonical_rc_sequence_sha256"):
            if len(row[key]) != 64 or any(c not in "0123456789abcdef" for c in row[key]):
                raise ValueError("Missing or invalid frozen sequence hash")
        ids.append(row["interval_id"])
        components.append(row["component_id"])
    if not ids or len(ids) != len(set(ids)):
        raise ValueError("Empty or duplicate common-panel membership")
    observed = {"positive": sum(int(r["label"]) == 1 for r in rows),
                "control": sum(int(r["label"]) == 0 for r in rows),
                "components": len(set(components))}
    if observed != expected:
        raise ValueError(f"Exact frozen test population mismatch: {model}: {observed} != {expected}")
    return ids


def load_test_predictions(path, metadata):
    ids = [row["interval_id"] for row in metadata]
    wanted = {row["interval_id"]: row["chrom"] for row in metadata}
    found = {}
    for row in ev.read_rows(path):
        key = row["interval_id"]
        if key not in wanted or row["chrom"] != wanted[key]:
            raise ValueError("Prediction membership/chromosome differs from exact test panel")
        if key in found:
            raise ValueError("Duplicate prediction interval")
        found[key] = row
    if set(found) != set(wanted):
        raise ValueError("Missing frozen test predictions")
    predictions = {}
    for config in CONFIGURATIONS:
        fwd = [[float(found[key][f"{config}_p_forward_{seed}"]) for key in ids] for seed in SEEDS]
        rc = [[float(found[key][f"{config}_p_rc_{seed}"]) for key in ids] for seed in SEEDS]
        symmetric, ensemble = ev.symmetric_ensemble(fwd, rc)
        for i, seed in enumerate(SEEDS):
            stored = ev.probability_array([float(found[key][f"{config}_q_{seed}"]) for key in ids], "stored seed q")
            if not np.array_equal(stored, symmetric[i]):
                raise ValueError("Stored seed score differs from exact float64 symmetric rule")
        stored = ev.probability_array([float(found[key][f"{config}_ensemble_q"]) for key in ids], "stored ensemble")
        if not np.array_equal(stored, ensemble):
            raise ValueError("Stored ensemble differs from equal-three-seed float64 rule")
        predictions[config] = {"seeds": symmetric, "ensemble": ensemble}
    return predictions


def predictor_map(predictions):
    return {(config, predictor): scores for config in CONFIGURATIONS
            for predictor, scores in [(f"seed:{seed}", predictions[config]["seeds"][i])
                                      for i, seed in enumerate(SEEDS)]
            + [("ensemble", predictions[config]["ensemble"])]}


def operating_metrics(labels, scores, threshold, weights=None):
    y = np.asarray(labels, dtype=np.int64)
    p = ev.probability_array(scores, "fixed-threshold probabilities")
    if y.ndim != 1 or p.shape != y.shape or not np.isin(y, [0, 1]).all():
        raise ValueError("Invalid fixed-threshold labels")
    w = np.ones(len(y), dtype=np.int64) if weights is None else np.asarray(weights)
    if w.shape != y.shape or not np.isfinite(w).all() or np.any(w < 0):
        raise ValueError("Invalid fixed-threshold weights")
    call = p >= threshold
    tp = float(w[(y == 1) & call].sum())
    fp = float(w[(y == 0) & call].sum())
    tn = float(w[(y == 0) & ~call].sum())
    fn = float(w[(y == 1) & ~call].sum())
    ratio = lambda a, b: a / b if b else None
    return {"TP": tp, "FP": fp, "TN": tn, "FN": fn,
            "sensitivity": ratio(tp, tp + fn), "specificity": ratio(tn, tn + fp),
            "precision": ratio(tp, tp + fp), "NPV": ratio(tn, tn + fn),
            "accuracy": ratio(tp + tn, tp + fp + tn + fn),
            "negative_row_FPR": ratio(fp, fp + tn)}


def component_bootstrap(labels, predictors, components, threshold, valid_replicates=2000,
                        max_attempts=20000, seed=314159):
    """Reuse exact inherited draw stream for all 12 predictors and C operating rates.

    Undefined precision/NPV are recorded as NA, not a reason to discard a draw
    which is valid for the prespecified discrimination/skill metrics.
    """
    names, inverse = np.unique(np.asarray(components, dtype=str), return_inverse=True)
    if len(inverse) != len(labels) or not len(names):
        raise ValueError("Invalid component vector")
    prepared = {key: ev.PreparedMetrics(labels, value) for key, value in predictors.items()}
    if ("V2-C", "ensemble") not in prepared:
        raise ValueError("Required fixed C ensemble absent")
    rng = np.random.default_rng(seed)
    attempts = invalid = valid = 0
    records, operations = [], []
    while valid < valid_replicates and attempts < max_attempts:
        attempts += 1
        multiplicities = np.bincount(rng.integers(0, len(names), size=len(names)), minlength=len(names))
        weights = multiplicities[inverse]
        try:
            metrics = {key: item.calculate(weights) for key, item in prepared.items()}
        except ValueError:
            invalid += 1
            continue
        valid += 1
        for (config, predictor), values in metrics.items():
            records.append({"valid_replicate": valid, "attempt": attempts,
                            "configuration": config, "predictor": predictor, **values})
        operations.append({"valid_replicate": valid, "attempt": attempts,
                           **operating_metrics(labels, predictors[("V2-C", "ensemble")], threshold, weights)})
    fraction = invalid / attempts if attempts else 1.0
    audit = {"seed": seed, "represented_components": len(names), "attempted": attempts,
             "invalid": invalid, "valid": valid, "invalid_fraction": fraction,
             "required_valid": valid_replicates, "max_attempts": max_attempts,
             "status": "PASS" if valid == valid_replicates and fraction < .1 else "INCONCLUSIVE",
             "all_seed_and_ensemble_predictors_share_draws": True,
             "undefined_operating_ratio_does_not_invalidate_metric_draw": True,
             "bootstrap_unit": "frozen global overlap/encoded-sequence component"}
    return records, operations, audit


def stratification(metadata, cuts, model):
    if cuts["model"] != model or cuts["quantile_method"] != "linear":
        raise ValueError("Wrong inherited train-derived cutpoints")
    strata = {}
    for variable in COVARIATES:
        values = np.asarray([float(row[variable]) for row in metadata])
        boundaries = np.asarray(cuts["cutpoints"][variable], dtype=np.float64)
        if boundaries.shape != (4,) or not np.isfinite(boundaries).all() or np.any(np.diff(boundaries) < 0):
            raise ValueError("Invalid inherited quintile boundaries")
        if not np.isfinite(values).all():
            raise ValueError("Nonfinite frozen test covariate")
        bins = np.searchsorted(boundaries, values, side="right") + 1
        strata[variable] = np.asarray([f"quintile_{i}" for i in bins])
    strata["atac_lobe_signature"] = np.asarray([r["atac_lobes"] or "none" for r in metadata])
    strata["same_lobe_mark_support"] = np.asarray([r[f"{model}_same_lobe_support_lobes"] or "none" for r in metadata])
    strata["chromosome"] = np.asarray([r["chrom"] for r in metadata])
    ambiguous = np.asarray([float(r["non_acgt_fraction"]) for r in metadata])
    if not np.isfinite(ambiguous).all() or np.any((ambiguous < 0) | (ambiguous > 1)):
        raise ValueError("Invalid ambiguous-base fraction")
    strata["ambiguous_base_status"] = np.where(ambiguous > 0, "present", "absent")
    return strata


def holdout_gates(model, labels, components, seed_aps, c_intervals, bootstrap, expected=None):
    expected = EXPECTED[model] if expected is None else expected
    y = np.asarray(labels)
    rows = []

    def gate(name, passed, observed, criterion, inconclusive=False):
        rows.append({"model": model, "gate": name,
                     "status": "PASS" if passed else ("INCONCLUSIVE" if inconclusive else "FAIL"),
                     "observed": observed, "criterion": criterion})

    seed_complete = len(seed_aps) == 3 and np.isfinite(seed_aps).all()
    gate("all_three_C_seeds_evaluable", seed_complete, len(seed_aps), "3 fixed seeds; no seed selection")
    gate("real_network_invariance", True, "PASS", "hash-bound actual-network all seeds/ensembles PASS")
    gate("finite_valid_probabilities", True, "PASS", "all orientations/seeds/ensembles finite in [0,1]")
    for kind, value in (("positive", int((y == 1).sum())), ("control", int((y == 0).sum())),
                        ("components", len(set(components)))):
        gate(f"exact_test_{kind}", value == expected[kind], value, f"exactly {expected[kind]}", True)
    gate("bootstrap_valid_replicates", bootstrap["valid"] == 2000, bootstrap["valid"], "exactly 2000 valid; at most 20000 attempts", True)
    gate("bootstrap_invalid_fraction", bootstrap["invalid_fraction"] < .1, bootstrap["invalid_fraction"], "strictly <0.10", True)
    for metric, boundary in (("AUROC", .5), ("AP_gain", 0.), ("BrierSkill", 0.)):
        lower = c_intervals.get(metric, {}).get("lower95")
        gate(f"lower95_{metric}", lower is not None and lower > boundary,
             lower, f"strictly >{boundary}", lower is None)
    sd = float(np.std(seed_aps, ddof=1)) if seed_complete else None
    spread = float(np.ptp(seed_aps)) if seed_complete else None
    gate("seed_AP_sample_SD", sd is not None and sd <= .03, sd, "<=0.03; ddof=1", sd is None)
    gate("seed_AP_range", spread is not None and spread <= .10, spread, "<=0.10", spread is None)
    statuses = {row["status"] for row in rows}
    status = "INCONCLUSIVE" if "INCONCLUSIVE" in statuses else ("FAIL" if "FAIL" in statuses else "PASS")
    return rows, status


def validate_prerequisites(manifest):
    if manifest["version"] != "internal-test-1.0":
        raise ValueError("Unexpected evaluation version")
    spec_path = ev.checked_reference(manifest["specification"])
    freeze_path = ev.checked_reference(manifest["prospective_freeze"])
    specification, freeze = ev.read_json(spec_path), ev.read_json(freeze_path)
    design = ev.read_json(ev.checked_reference(manifest["frozen_design"]))
    if design["experimental_version"] != "pretraining-1.1" or tuple(design["seeds"]) != SEEDS:
        raise ValueError("Immutable frozen design mismatch")
    bootstrap = design["selection"]["bootstrap"]
    if (bootstrap["seed"], bootstrap["valid_replicates"], bootstrap["max_attempts"], bootstrap["max_invalid_fraction_exclusive"]) != (314159, 2000, 20000, .1):
        raise ValueError("Inherited bootstrap procedure changed")
    # The root preparation audit additionally establishes pre-inference ordering.
    # Accept both path/hash records and an explicit specification_sha256 field.
    frozen_spec_hash = freeze.get("specification_sha256", freeze.get("specification", {}).get("sha256"))
    if frozen_spec_hash != ev.sha256(spec_path):
        raise ValueError("Prospective freeze does not bind this specification")
    if (freeze.get("status") != "PASS" or freeze.get("before_any_test_model_inference") is not True
            or freeze.get("predictions_present") is not False):
        raise ValueError("Specification not prospectively frozen before test inference")
    if specification.get("version", specification.get("experimental_version")) != "internal-test-1.0":
        raise ValueError("Prospective specification version mismatch")
    pre = manifest["prerequisites"]
    for key in ("input_integrity", "checkpoint_verification"):
        record = pre[key]
        evidence = ev.read_json(ev.checked_reference(record["evidence"]))
        if record["status"] != "PASS" or evidence.get("status") != "PASS":
            raise ValueError(f"Mandatory {key} not PASS")
    invariance = pre["invariance"]
    evidence = ev.read_json(ev.checked_reference(invariance["evidence"]))
    if not (invariance["status"] == "PASS" and invariance["real_network"] is True
            and invariance["all_configurations_all_seeds_and_ensembles"] is True
            and evidence.get("status") == "PASS" and evidence.get("real_network") is True
            and evidence.get("all_configurations_all_seeds_and_ensembles") is True):
        raise ValueError("Mandatory real-network RC invariance not PASS")
    expected_units = {(c, m, "seed", str(s)) for c in CONFIGURATIONS for m in MODELS for s in SEEDS}
    expected_units |= {(c, m, "ensemble", "all3") for c in CONFIGURATIONS for m in MODELS}
    seen_units = set()
    for row in evidence["audit_rows"]:
        key = (row["configuration"], row["model"], row["unit"], str(row["seed"]))
        if key in seen_units or key not in expected_units:
            raise ValueError("Duplicate/unexpected actual-network audit unit")
        seen_units.add(key)
        count = EXPECTED[row["model"]]["positive"] + EXPECTED[row["model"]]["control"]
        if not (row["status"] == "PASS" and row["real_network"] is True and row["n_sequences"] == count
                and row["n_failed"] == 0 and row["atol"] == 1e-6 and row["rtol"] == 1e-6
                and np.isfinite(row["max_tolerance_normalized_difference"])
                and 0 <= row["max_tolerance_normalized_difference"] <= 1):
            raise ValueError("Actual-network test invariance coverage/tolerance failure")
    if seen_units != expected_units or evidence["n_seed_audits"] != 18 or evidence["n_ensemble_audits"] != 6:
        raise ValueError("Incomplete actual-network test invariance audit")
    seen = set()
    for checkpoint in pre["checkpoints"]:
        key = (checkpoint["configuration"], checkpoint["model"], int(checkpoint["seed"]))
        if key in seen or checkpoint["completed"] is not True:
            raise ValueError("Repeated/incomplete checkpoint record")
        ev.checked_reference(checkpoint)
        seen.add(key)
    if seen != {(c, m, s) for c in CONFIGURATIONS for m in MODELS for s in SEEDS}:
        raise ValueError("All 18 original frozen checkpoints required")
    if set(manifest["models"]) != set(MODELS):
        raise ValueError("Both contexts required")
    evidence_file_hashes = {row["sha256"] for row in evidence["files"]}
    if any(manifest["models"][m]["predictions"]["sha256"] not in evidence_file_hashes for m in MODELS):
        raise ValueError("Actual-network evidence does not bind evaluation prediction files")
    checkpoint_hashes = {(r["configuration"], r["model"], int(r["seed"])): r["sha256"] for r in pre["checkpoints"]}
    for row in evidence["audit_rows"]:
        if row["unit"] == "seed" and row["checkpoint"]["sha256"] != checkpoint_hashes[(row["configuration"], row["model"], int(row["seed"]))]:
            raise ValueError("Actual-network evidence used a different checkpoint")
    cutpoints = ev.read_json(ev.checked_reference(manifest["frozen_cutpoints"]))
    thresholds = ev.read_json(ev.checked_reference(manifest["frozen_thresholds"]))
    if thresholds["status"] != "PASS":
        raise ValueError("Frozen chr7 threshold source not PASS")
    for model in MODELS:
        row = thresholds["models"][model]
        value = float.fromhex(THRESHOLDS[model])
        if not (row["threshold_hex"] == THRESHOLDS[model]
                and float(row["threshold_decimal_17g"]) == value
                and row["threshold_float64"] == value and row["call_rule"] == "score>=threshold"):
            raise ValueError("Frozen threshold exact-encoding drift")
    return cutpoints, thresholds


def evaluate(manifest_path, output):
    manifest = ev.read_json(manifest_path)
    cutpoints, thresholds = validate_prerequisites(manifest)
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise ValueError("Refusing rerun/overwrite of nonempty test evaluation output")
    output.mkdir(parents=True, exist_ok=True)
    results = {key: [] for key in ("metrics", "reliability_bins", "stratified_metrics", "seed_stability",
                                   "bootstrap_replicates", "bootstrap_summary", "paired_ablation",
                                   "C_holdout_generalization", "fixed_threshold_operating",
                                   "fixed_threshold_intervals", "fixed_threshold_bootstrap", "common_test_scores")}
    decisions, audits = {}, {}
    for model in MODELS:
        item = manifest["models"][model]
        metadata = list(ev.read_rows(ev.checked_reference(item["test"])))
        ids = validate_test_rows(metadata, model)
        predictions = load_test_predictions(ev.checked_reference(item["predictions"]), metadata)
        predictors = predictor_map(predictions)
        labels = np.asarray([int(row["label"]) for row in metadata], dtype=np.int64)
        components = np.asarray([row["component_id"] for row in metadata])
        strata = stratification(metadata, cutpoints[model], model)
        points = {}
        for (config, predictor), scores in predictors.items():
            prefix = {"model": model, "configuration": config, "predictor": predictor}
            points[(config, predictor)] = metrics = ev.PreparedMetrics(labels, scores).calculate()
            results["metrics"].append({**prefix, "n": len(labels), "n_positive": int(labels.sum()),
                                       "n_control": int((1 - labels).sum()), "n_components": len(set(components)), **metrics})
            results["reliability_bins"].extend({**prefix, **row} for row in ev.reliability(labels, scores))
            for variable, bins in strata.items():
                categories = [f"quintile_{i}" for i in range(1, 6)] if variable in COVARIATES else sorted(set(bins))
                for category in categories:
                    mask = bins == category
                    n = int(mask.sum())
                    try:
                        submetrics = ev.PreparedMetrics(labels[mask], scores[mask]).calculate()
                        defined = True
                    except ValueError:
                        submetrics = {metric: None for metric in METRICS}
                        submetrics["prevalence"] = float(labels[mask].mean()) if n else None
                        defined = False
                    results["stratified_metrics"].append({**prefix, "stratifier": variable, "stratum": category,
                        "n": n, "n_positive": int(labels[mask].sum()), "n_control": int((1-labels[mask]).sum()),
                        "n_components": len(set(components[mask])), "cells_lt100_descriptive_only": n < 100,
                        "all_stratified_results_descriptive": True, "metrics_defined": defined, **submetrics})
        for config in CONFIGURATIONS:
            aps = [points[(config, f"seed:{s}")]["AP"] for s in SEEDS]
            results["seed_stability"].append({"model": model, "configuration": config, "seeds_complete": 3,
                "AP_sample_SD": float(np.std(aps, ddof=1)), "AP_range": float(np.ptp(aps)),
                **{f"AP_seed_{s}": aps[i] for i, s in enumerate(SEEDS)}})
            for i, key in enumerate(ids):
                results["common_test_scores"].append({"model": model, "configuration": config,
                    "interval_id": key, "chrom": metadata[i]["chrom"], "partition": "test", "label": int(labels[i]),
                    "component_id": components[i], **{f"q_{s}": float(predictions[config]["seeds"][j, i]) for j, s in enumerate(SEEDS)},
                    "ensemble_q": float(predictions[config]["ensemble"][i])})
        threshold = float.fromhex(THRESHOLDS[model])
        draws, operating_draws, audit = component_bootstrap(labels, predictors, components, threshold)
        audits[model] = audit
        results["bootstrap_replicates"].extend({"model": model, **r} for r in draws)
        results["fixed_threshold_bootstrap"].extend({"model": model, **r} for r in operating_draws)
        grouped = {key: [] for key in predictors}
        for row in draws:
            grouped[(row["configuration"], row["predictor"])].append(row)
        intervals = {key: {} for key in predictors}
        for (config, predictor), rows in grouped.items():
            if not rows:
                continue
            for metric in METRICS:
                limits = ev.interval([r[metric] for r in rows])
                intervals[(config, predictor)][metric] = limits
                results["bootstrap_summary"].append({"model": model, "configuration": config,
                    "predictor": predictor, "metric": metric, "point": points[(config, predictor)][metric], **limits,
                    "valid_replicates": len(rows), "bootstrap_status": audit["status"]})
        if draws:
            for lhs, rhs in (("V2-B", "V2-A"), ("V2-C", "V2-B"), ("V2-C", "V2-A")):
                for metric in ("AP", "AUROC", "Brier"):
                    diffs = [a[metric] - b[metric] for a, b in zip(grouped[(lhs, "ensemble")], grouped[(rhs, "ensemble")])]
                    results["paired_ablation"].append({"model": model, "comparison": f"{lhs} minus {rhs}",
                        "metric": metric, "point": points[(lhs, "ensemble")][metric] - points[(rhs, "ensemble")][metric],
                        **ev.interval(diffs), "diagnostic_only_not_release_gate": True,
                        "bootstrap_status": audit["status"]})
        seed_aps = [points[("V2-C", f"seed:{s}")]["AP"] for s in SEEDS]
        gates, status = holdout_gates(model, labels, components, seed_aps, intervals[("V2-C", "ensemble")], audit)
        results["C_holdout_generalization"].extend(gates)
        decisions[model] = {"status": status, "failed_or_inconclusive_gates": [r for r in gates if r["status"] != "PASS"]}
        operating = operating_metrics(labels, predictors[("V2-C", "ensemble")], threshold)
        row = thresholds["models"][model]
        results["fixed_threshold_operating"].append({"model": model, "configuration": "V2-C", "predictor": "ensemble",
            "threshold_decimal_17g": row["threshold_decimal_17g"], "threshold_hex": row["threshold_hex"],
            "call_rule": "score>=threshold", **{k: int(v) if k in ("TP", "FP", "TN", "FN") else v for k, v in operating.items()},
            "observed_chr7_calibration_FPR": row["observed_calibration_row_FPR"],
            "test_minus_chr7_calibration_FPR": operating["negative_row_FPR"] - row["observed_calibration_row_FPR"],
            "test_FPR_above_0_05": operating["negative_row_FPR"] > .05,
            "FPR_excess_is_not_protocol_failure": True, "threshold_recalibrated": False})
        for metric in RATE_NAMES:
            valid_values = [r[metric] for r in operating_draws if r[metric] is not None]
            confidence = ev.interval(valid_values) if valid_values else {key: None for key in ("lower95", "upper95", "lower90", "upper90")}
            results["fixed_threshold_intervals"].append({"model": model, "metric": metric,
                "point": operating[metric], **confidence, "valid_ratio_replicates": len(valid_values),
                "undefined_ratio_replicates": len(operating_draws) - len(valid_values),
                "total_valid_metric_draws": len(operating_draws), "bootstrap_status": audit["status"],
                "CI_scope": "conditional component-bootstrap percentile interval on defined ratios; not binomial independence"})
    both_pass = all(decisions[m]["status"] == "PASS" for m in MODELS)
    decision = {"module": "COPD-V2-INTERNAL-TEST", "version": "internal-test-1.0",
        "input_manifest_sha256": ev.sha256(manifest_path), "model_contexts": decisions,
        "both_C_contexts_pass": both_pass, "ready_for_investigator_review_before_external_evaluation": both_pass,
        "next_action": "STOP_FOR_INVESTIGATOR_REVIEW_NO_EXTERNAL_OR_CANDIDATE_ACCESS" if both_pass else "STOP_FAILED_OR_INCONCLUSIVE_HOLDOUT_NO_REDESIGN_RETRAINING_OR_EXTERNAL_ACCESS",
        "A_B_ranking_used_for_release": False, "thresholds_changed": False,
        "external_access_authorized": False, "candidate_scoring_authorized": False,
        "interpretation": "fixed internal V2 holdout; not historically pristine or independent external validation"}
    for name, rows in results.items():
        suffix = ".tsv.gz" if name in ("bootstrap_replicates", "fixed_threshold_bootstrap", "common_test_scores") else ".tsv"
        ev.write_rows(output / (name + suffix), rows)
    ev.write_json(output / "bootstrap_audit.json", audits)
    ev.write_json(output / "holdout_decision.json", decision)
    ev.write_json(output / "evaluation_manifest.json", {"input_manifest": str(manifest_path),
        "input_manifest_sha256": ev.sha256(manifest_path), "evaluation_script_sha256": ev.sha256(__file__),
        "immutable_metric_script": {"path": str(FROZEN_EVALUATOR), "sha256": FROZEN_EVALUATOR_SHA256},
        "numpy_version": np.__version__, "frozen_cutpoints_reference": manifest["frozen_cutpoints"],
        "frozen_thresholds_reference": manifest["frozen_thresholds"],
        "artifacts": [{"path": str(p), "bytes": p.stat().st_size, "sha256": ev.sha256(p)} for p in sorted(output.iterdir())],
        "prediction_scope": "only exact frozen C-task chr8/9 region-label common panels",
        "bootstrap_generalization": "conditional fixed-matching genomic-component uncertainty; not donor/population generalization",
        "no_external_or_candidate_access": True, "no_recalibration": True})
    return decision


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-manifest", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    try:
        decision = evaluate(args.input_manifest, args.output)
    except Exception as error:
        print(json.dumps({"status": "HARD_STOP", "error": str(error), "external_access_permitted": False}), file=sys.stderr)
        return 2
    print(json.dumps(decision, sort_keys=True))
    return 0 if decision["both_C_contexts_pass"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
