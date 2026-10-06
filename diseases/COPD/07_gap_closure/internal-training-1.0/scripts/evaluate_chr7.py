#!/usr/bin/env python3
"""Frozen COPD V2 common-C-task chr7 evaluation; no model execution or test access.

The input JSON schema and pre-outcome implementation decisions are documented in
../specification/internal_evaluation_implementation.json. Outputs are append-only:
an existing nonempty output directory is rejected. All path references are
{path, sha256} objects, resolved relative to the current working directory.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import sys

import numpy as np

SEEDS = (104729, 130363, 155921)
CONFIGURATIONS = ("V2-A", "V2-B", "V2-C")
MODELS = ("enhancer", "h3k27me3")
METRICS = ("AP", "AUROC", "Brier", "BrierSkill", "AP_gain", "prevalence")
COVARIATES = ("gc_fraction", "atac_signal_percentile_max_train_only", "repeat_fraction_2001")
TRAIN_CHROMS = {f"chr{x}" for x in (1, 2, 3, 4, 5, 6, *range(10, 23))} | {"chrX", "chrY"}


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checked_reference(reference):
    path = Path(reference["path"])
    if sha256(path) != reference["sha256"]:
        raise ValueError(f"Hash mismatch: {path}")
    return path


def read_json(path):
    with open(path, encoding="utf-8") as stream:
        return json.load(stream)


def write_json(path, obj):
    with open(path, "x", encoding="utf-8") as stream:
        json.dump(obj, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def read_rows(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", newline="") as stream:
        yield from csv.DictReader(stream, delimiter="\t")


def write_rows(path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = list(rows[0]) if rows else []
    if str(path).endswith(".gz"):
        binary = open(path, "xb")
        stream = io.TextIOWrapper(gzip.GzipFile(filename="", mode="wb", fileobj=binary, mtime=0), encoding="utf-8", newline="")
    else:
        binary = None
        stream = open(path, "x", encoding="utf-8", newline="")
    try:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    finally:
        stream.close()
        if binary:
            binary.close()


def probability_array(values, name):
    values = np.asarray(values, dtype=np.float64)
    if not np.all(np.isfinite(values)) or np.any(values < 0) or np.any(values > 1):
        raise ValueError(f"Invalid probabilities: {name}")
    return values


def symmetric_ensemble(forward, reverse):
    forward = probability_array(forward, "forward")
    reverse = probability_array(reverse, "reverse")
    if forward.shape != reverse.shape or forward.ndim != 2 or forward.shape[0] != 3:
        raise ValueError("Expected matching arrays with shape (three fixed seeds, intervals)")
    symmetric = (forward + reverse) / np.float64(2)
    ensemble = np.mean(symmetric, axis=0, dtype=np.float64)
    return symmetric, ensemble


class PreparedMetrics:
    """Weighted tie-aware AP/AUROC; component multiplicity equals row repetition."""

    def __init__(self, labels, scores):
        self.y = np.asarray(labels, dtype=np.int64)
        self.p = probability_array(scores, "metric scores")
        if self.y.shape != self.p.shape or self.y.ndim != 1 or not len(self.y):
            raise ValueError("Labels and probabilities must be nonempty matching vectors")
        if not np.all(np.isin(self.y, [0, 1])):
            raise ValueError("Labels must be binary")
        self.order = np.argsort(-self.p, kind="mergesort")
        ordered_p = self.p[self.order]
        self.ends = np.r_[np.flatnonzero(np.diff(ordered_p)), len(ordered_p) - 1]
        self.error = (self.p - self.y) ** 2

    def calculate(self, weights=None):
        if weights is None:
            weights = np.ones(len(self.y), dtype=np.float64)
        weights = np.asarray(weights, dtype=np.float64)
        if weights.shape != self.y.shape or not np.isfinite(weights).all() or np.any(weights < 0):
            raise ValueError("Invalid metric weights")
        n = float(np.sum(weights))
        pos = float(np.dot(weights, self.y))
        neg = n - pos
        if n <= 0 or pos <= 0 or neg <= 0:
            raise ValueError("Undefined metric in empty or one-class draw")
        y = self.y[self.order]
        w = weights[self.order]
        tp = np.cumsum(w * y)[self.ends]
        fp = np.cumsum(w * (1 - y))[self.ends]
        delta_tp = np.diff(np.r_[0.0, tp])
        delta_fp = np.diff(np.r_[0.0, fp])
        precision = np.divide(tp, tp + fp, out=np.zeros_like(tp), where=(tp + fp) > 0)
        ap = float(np.dot(delta_tp, precision) / pos)
        auc = float(np.dot(delta_fp, tp - 0.5 * delta_tp) / (pos * neg))
        prevalence = pos / n
        brier = float(np.dot(weights, self.error) / n)
        skill = 1 - brier / (prevalence * (1 - prevalence))
        result = {"AP": ap, "AUROC": auc, "Brier": brier, "BrierSkill": skill,
                  "AP_gain": ap - prevalence, "prevalence": prevalence}
        if not all(math.isfinite(x) for x in result.values()):
            raise ValueError("Undefined metric")
        return result


def component_bootstrap(labels, ensembles, components, valid_replicates=2000, max_attempts=20000, seed=314159):
    """Same component draws and multiplicities for all A/B/C predictions."""
    names, inverse = np.unique(np.asarray(components, dtype=str), return_inverse=True)
    if not len(names):
        raise ValueError("No represented components")
    rng = np.random.default_rng(seed)
    prepared = {config: PreparedMetrics(labels, ensembles[config]) for config in CONFIGURATIONS}
    attempts = invalid = valid = 0
    records = []
    while valid < valid_replicates and attempts < max_attempts:
        attempts += 1
        multiplicities = np.bincount(rng.integers(0, len(names), size=len(names)), minlength=len(names))
        weights = multiplicities[inverse]
        try:
            metrics = {config: prepared[config].calculate(weights) for config in CONFIGURATIONS}
        except ValueError:
            invalid += 1
            continue
        valid += 1
        for config in CONFIGURATIONS:
            records.append({"valid_replicate": valid, "attempt": attempts, "configuration": config, **metrics[config]})
    invalid_fraction = invalid / attempts if attempts else 1.0
    return records, {"seed": seed, "represented_components": len(names), "attempted": attempts,
                     "invalid": invalid, "valid": valid, "invalid_fraction": invalid_fraction,
                     "required_valid": valid_replicates, "max_attempts": max_attempts,
                     "status": "PASS" if valid == valid_replicates and invalid_fraction < 0.1 else "INCONCLUSIVE"}


def interval(values):
    q = np.percentile(np.asarray(values, dtype=np.float64), [2.5, 97.5, 5, 95], method="linear")
    return dict(zip(("lower95", "upper95", "lower90", "upper90"), map(float, q)))


def reliability(labels, scores):
    y = np.asarray(labels, dtype=np.int64)
    p = probability_array(scores, "reliability")
    assignment = np.minimum((p * 10).astype(np.int64), 9)
    result = []
    for index in range(10):
        mask = assignment == index
        n = int(mask.sum())
        result.append({"bin": index, "lower": index / 10, "upper": (index + 1) / 10,
                       "upper_inclusive": index == 9, "n": n,
                       "mean_probability": float(p[mask].mean()) if n else None,
                       "observed_positive_fraction": float(y[mask].mean()) if n else None})
    return result


def frozen_cutpoints(train_rows, model):
    if not train_rows:
        raise ValueError("Missing retained C training positives")
    seen = set()
    for row in train_rows:
        if row["chrom"] not in TRAIN_CHROMS or row["partition"] != "train" or int(row["label"]) != 1:
            raise ValueError("Cutpoints require training-chromosome retained C positives only")
        if row.get("configuration") != "V2-C" or row.get("model") != model:
            raise ValueError("Cutpoint configuration/model drift")
        if row["interval_id"] in seen:
            raise ValueError("Duplicate training covariate interval")
        seen.add(row["interval_id"])
    output = {"model": model, "n_retained_C_train_positives": len(train_rows), "quantile_method": "linear", "boundary_assignment": "searchsorted(side=right); repeated boundaries may leave empty bins", "cutpoints": {}}
    for variable in COVARIATES:
        values = np.asarray([float(row[variable]) for row in train_rows], dtype=np.float64)
        if not np.isfinite(values).all():
            raise ValueError(f"Nonfinite training covariates: {variable}")
        output["cutpoints"][variable] = np.quantile(values, [0.2, 0.4, 0.6, 0.8], method="linear").tolist()
    return output


def validate_common_rows(rows, model, role, negatives_only=False):
    if not rows:
        raise ValueError("Empty common panel")
    ids = []
    component_roles = {}
    for row in rows:
        if row["chrom"] != "chr7" or row["partition"] != "validation" or row["validation_role"] != role:
            raise ValueError(f"Forbidden chromosome/role in {role} panel")
        if row.get("configuration") != "V2-C" or row.get("model") != model:
            raise ValueError("Common panel must have frozen C-task labels")
        label = int(row["label"])
        if label not in (0, 1) or (negatives_only and label != 0):
            raise ValueError("Invalid label or forbidden positive in calibration input")
        if not row["component_id"]:
            raise ValueError("Missing component ID")
        component_roles[row["component_id"]] = role
        ids.append(row["interval_id"])
    if len(ids) != len(set(ids)):
        raise ValueError("Duplicate common-panel interval IDs")
    return ids


def load_predictions(path, ids, configurations=CONFIGURATIONS):
    wanted = set(ids)
    found = {}
    for row in read_rows(path):
        # Never convert probabilities for rows outside the authorized role subset.
        if row["interval_id"] not in wanted:
            continue
        if row["chrom"] != "chr7":
            raise ValueError("Prediction input contains forbidden chromosome for requested ID")
        if row["interval_id"] in found:
            raise ValueError("Duplicate prediction interval")
        found[row["interval_id"]] = row
    if set(found) != wanted:
        raise ValueError("Missing requested common-panel predictions")
    output = {}
    for config in configurations:
        fwd = [[float(found[key][f"{config}_p_forward_{seed}"]) for key in ids] for seed in SEEDS]
        rc = [[float(found[key][f"{config}_p_rc_{seed}"]) for key in ids] for seed in SEEDS]
        symmetric, ensemble = symmetric_ensemble(fwd, rc)
        for index, seed in enumerate(SEEDS):
            column = f"{config}_q_{seed}"
            if column in found[ids[0]]:
                stored = probability_array([float(found[key][column]) for key in ids], column)
                if not np.array_equal(stored, symmetric[index]):
                    raise ValueError(f"Stored/recomputed symmetric score differs: {config}, {seed}")
        column = f"{config}_ensemble_q"
        if column in found[ids[0]]:
            stored = probability_array([float(found[key][column]) for key in ids], column)
            if not np.array_equal(stored, ensemble):
                raise ValueError(f"Stored/recomputed ensemble differs: {config}")
        output[config] = {"seeds": symmetric, "ensemble": ensemble}
    return output


def validate_prerequisites(manifest):
    spec = read_json(checked_reference(manifest["specification"]))
    if spec["experimental_version"] != "pretraining-1.1" or tuple(spec["seeds"]) != SEEDS:
        raise ValueError("Wrong immutable design")
    if spec["selection"]["bootstrap"]["valid_replicates"] != 2000 or spec["selection"]["bootstrap"]["seed"] != 314159:
        raise ValueError("Bootstrap contract drift")
    pre = manifest["prerequisites"]
    if pre["construction_gates_pass"] is not True:
        raise ValueError("Construction gates not PASS")
    checked_reference(pre["construction_evidence"])
    checked_reference(pre["checkpoint_freeze"])
    invariance = pre["invariance"]
    checked_reference(invariance["evidence"])
    if not (invariance["status"] == "PASS" and invariance["real_network"] is True
            and invariance["all_configurations_all_seeds_and_ensembles"] is True):
        raise ValueError("Mandatory real-network invariance not PASS")
    observed = set()
    for checkpoint in pre["checkpoints"]:
        key = (checkpoint["configuration"], checkpoint["model"], int(checkpoint["seed"]))
        if key in observed or checkpoint["completed"] is not True:
            raise ValueError("Duplicate or incomplete selected checkpoint")
        checked_reference(checkpoint)
        observed.add(key)
    expected = {(config, model, seed) for config in CONFIGURATIONS for model in MODELS for seed in SEEDS}
    if observed != expected:
        raise ValueError("All 18 fixed-seed selected checkpoints must be frozen before evaluation")
    if set(manifest["models"]) != set(MODELS):
        raise ValueError("Both model contexts are required")
    return spec


def adequacy_gates(model, labels, components, seed_aps, c_intervals, bootstrap):
    y = np.asarray(labels)
    stability_sd = float(np.std(seed_aps, ddof=1))
    stability_range = float(np.ptp(seed_aps))
    rows = []

    def gate(name, passed, observed, criterion, inconclusive=False):
        rows.append({"model": model, "gate": name, "status": "PASS" if passed else ("INCONCLUSIVE" if inconclusive else "FAIL"), "observed": observed, "criterion": criterion})

    gate("frozen_construction", True, "PASS", "B and C construction gates PASS")
    gate("all_three_C_seeds_completed", len(seed_aps) == 3, len(seed_aps), "3 fixed seeds; all checkpoints frozen")
    gate("real_network_invariance", True, "PASS", "actual-network all seeds and ensembles PASS")
    gate("finite_valid_probabilities", True, "PASS", "all orientation/seed/ensemble probabilities finite in [0,1]")
    gate("selection_positives", int(np.sum(y == 1)) >= 100, int(np.sum(y == 1)), ">=100", True)
    gate("selection_controls", int(np.sum(y == 0)) >= 200, int(np.sum(y == 0)), ">=200", True)
    gate("selection_components", len(set(components)) >= 30, len(set(components)), ">=30", True)
    gate("bootstrap_valid_replicates", bootstrap["valid"] == 2000, bootstrap["valid"], "exactly2000valid; at most20000attempts", True)
    gate("bootstrap_invalid_fraction", bootstrap["invalid_fraction"] < 0.1, bootstrap["invalid_fraction"], "strictly<0.10", True)
    for metric, boundary in (("AUROC", 0.5), ("AP_gain", 0.0), ("BrierSkill", 0.0)):
        low = c_intervals.get(metric, {}).get("lower95")
        gate(f"lower95_{metric}", low is not None and low > boundary, low, f"strictly>{boundary}", low is None)
    gate("seed_AP_sample_SD", stability_sd <= 0.03, stability_sd, "<=0.03; ddof=1")
    gate("seed_AP_range", stability_range <= 0.10, stability_range, "<=0.10")
    statuses = [row["status"] for row in rows]
    status = "INCONCLUSIVE" if "INCONCLUSIVE" in statuses else ("FAIL" if "FAIL" in statuses else "PASS")
    return rows, status


def evaluate(manifest_path, output):
    manifest = read_json(manifest_path)
    validate_prerequisites(manifest)
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise ValueError("Refusing outcome-based rerun/overwrite of nonempty output directory")
    output.mkdir(parents=True, exist_ok=True)
    metric_rows, reliability_rows, strata_rows, stability_rows = [], [], [], []
    bootstrap_rows, summary_rows, comparison_rows, gate_rows = [], [], [], []
    decisions, audits, cutpoint_records, score_rows = {}, {}, {}, []
    for model in MODELS:
        item = manifest["models"][model]
        metadata = list(read_rows(checked_reference(item["selection"])))
        ids = validate_common_rows(metadata, model, "selection")
        cuts = frozen_cutpoints(list(read_rows(checked_reference(item["train_covariates"]))), model)
        cutpoint_records[model] = cuts
        predictions = load_predictions(checked_reference(item["predictions"]), ids)
        labels = np.asarray([int(row["label"]) for row in metadata], dtype=np.int64)
        components = [row["component_id"] for row in metadata]
        ensembles = {config: predictions[config]["ensemble"] for config in CONFIGURATIONS}
        point = {}
        seed_aps = {}
        strata = {}
        for variable in COVARIATES:
            values = np.asarray([float(row[variable]) for row in metadata])
            if not np.isfinite(values).all():
                raise ValueError(f"Nonfinite selection covariate: {variable}")
            bins = np.searchsorted(cuts["cutpoints"][variable], values, side="right") + 1
            strata[variable] = np.asarray([f"quintile_{index}" for index in bins])
        strata["atac_lobe_signature"] = np.asarray([row["atac_lobes"] or "none" for row in metadata])
        strata["same_lobe_mark_support"] = np.asarray([row[f"{model}_same_lobe_support_lobes"] or "none" for row in metadata])
        strata["chromosome"] = np.asarray([row["chrom"] for row in metadata])
        strata["ambiguous_base_status"] = np.asarray(["present" if float(row["non_acgt_fraction"]) > 0 else "absent" for row in metadata])
        for config in CONFIGURATIONS:
            seed_aps[config] = []
            for index, key in enumerate(ids):
                score_rows.append({"model": model, "configuration": config, "interval_id": key, "chrom": "chr7", "validation_role": "selection", "label": int(labels[index]), "component_id": components[index], **{f"q_{seed}": float(predictions[config]["seeds"][j, index]) for j, seed in enumerate(SEEDS)}, "ensemble_q": float(ensembles[config][index])})
            predictor_items = [(f"seed:{seed}", predictions[config]["seeds"][index]) for index, seed in enumerate(SEEDS)] + [("ensemble", ensembles[config])]
            for predictor, scores in predictor_items:
                metrics = PreparedMetrics(labels, scores).calculate()
                prefix = {"model": model, "configuration": config, "predictor": predictor}
                metric_rows.append({**prefix, "n": len(labels), "n_positive": int(labels.sum()), "n_control": int((1 - labels).sum()), "n_components": len(set(components)), **metrics})
                if predictor == "ensemble":
                    point[config] = metrics
                else:
                    seed_aps[config].append(metrics["AP"])
                reliability_rows.extend({**prefix, **record} for record in reliability(labels, scores))
                for variable, bins in strata.items():
                    for category in sorted(set(bins)):
                        mask = bins == category
                        n = int(mask.sum())
                        try:
                            submetrics = PreparedMetrics(labels[mask], scores[mask]).calculate()
                            defined = True
                        except ValueError:
                            submetrics = {metric: None for metric in METRICS}
                            submetrics["prevalence"] = float(labels[mask].mean())
                            defined = False
                        strata_rows.append({**prefix, "stratifier": variable, "stratum": category, "n": n, "n_positive": int(labels[mask].sum()), "n_control": int((1 - labels[mask]).sum()), "n_components": len(set(np.asarray(components)[mask])), "cells_lt100_descriptive_only": n < 100, "all_stratified_results_descriptive": True, "metrics_defined": defined, **submetrics})
            stability_rows.append({"model": model, "configuration": config, "seeds_complete": 3, "AP_sample_SD": float(np.std(seed_aps[config], ddof=1)), "AP_range": float(np.ptp(seed_aps[config])), **{f"AP_seed_{seed}": seed_aps[config][index] for index, seed in enumerate(SEEDS)}})
        draws, audit = component_bootstrap(labels, ensembles, components)
        audits[model] = audit
        bootstrap_rows.extend({"model": model, **row} for row in draws)
        intervals = {config: {} for config in CONFIGURATIONS}
        if draws:
            draws_by_config = {config: [row for row in draws if row["configuration"] == config] for config in CONFIGURATIONS}
            for config in CONFIGURATIONS:
                for metric in METRICS:
                    confidence = interval([row[metric] for row in draws_by_config[config]])
                    intervals[config][metric] = confidence
                    summary_rows.append({"model": model, "configuration": config, "metric": metric, "point": point[config][metric], **confidence})
            for lhs, rhs in (("V2-B", "V2-A"), ("V2-C", "V2-B"), ("V2-C", "V2-A")):
                for metric in ("AP", "AUROC", "Brier"):
                    differences = [a[metric] - b[metric] for a, b in zip(draws_by_config[lhs], draws_by_config[rhs])]
                    comparison_rows.append({"model": model, "comparison": f"{lhs} minus {rhs}", "metric": metric, "point": point[lhs][metric] - point[rhs][metric], **interval(differences), "diagnostic_only_not_retention_gate": True})
        gates, status = adequacy_gates(model, labels, components, seed_aps["V2-C"], intervals["V2-C"], audit)
        gate_rows.extend(gates)
        decisions[model] = {"status": status, "failed_or_inconclusive_gates": [row for row in gates if row["status"] != "PASS"]}
    both_pass = all(decisions[model]["status"] == "PASS" for model in MODELS)
    decision = {"module": "COPD-V2-INTERNAL", "version": "internal-training-1.0", "input_manifest_sha256": sha256(manifest_path), "model_contexts": decisions, "both_C_contexts_pass": both_pass, "next_action": "FREEZE_C_CHECKPOINTS_ENSEMBLE_AND_ADEQUACY_THEN_CALIBRATE_CHR7_NEGATIVES_ONLY" if both_pass else "STOP_NO_CALIBRATION_NO_AB_FALLBACK", "A_B_ranking_used_for_retention": False, "test_access_authorized": False, "external_access_authorized": False, "candidate_scoring_authorized": False}
    outputs = {"metrics.tsv": metric_rows, "reliability_bins.tsv": reliability_rows, "stratified_metrics.tsv": strata_rows, "seed_stability.tsv": stability_rows, "bootstrap_replicates.tsv.gz": bootstrap_rows, "bootstrap_summary.tsv": summary_rows, "paired_ablation.tsv": comparison_rows, "C_absolute_adequacy.tsv": gate_rows, "common_selection_scores.tsv.gz": score_rows}
    for filename, rows in outputs.items():
        write_rows(output / filename, rows)
    write_json(output / "bootstrap_audit.json", audits)
    write_json(output / "stratification_cutpoints.json", cutpoint_records)
    write_json(output / "adequacy_decision.json", decision)
    write_json(output / "evaluation_manifest.json", {"input_manifest": str(manifest_path), "input_manifest_sha256": sha256(manifest_path), "evaluation_script_sha256": sha256(__file__), "numpy_version": np.__version__, "artifacts": [{"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)} for path in sorted(output.iterdir())], "prediction_scope": "common C-task chr7 selection only; nonselected prediction rows not converted to numbers", "bootstrap_generalization": "conditional internal genomic-component uncertainty; not donor/population generalization", "no_test_external_or_candidate_access": True})
    return decision


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-manifest", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    try:
        decision = evaluate(args.input_manifest, args.output)
    except Exception as error:
        print(json.dumps({"status": "HARD_STOP", "error": str(error), "calibration_permitted": False}), file=sys.stderr)
        return 2
    print(json.dumps(decision, sort_keys=True))
    return 0 if decision["both_C_contexts_pass"] else 3


if __name__ == "__main__":
    raise SystemExit(main())
