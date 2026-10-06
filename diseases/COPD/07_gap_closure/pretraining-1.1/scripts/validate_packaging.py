#!/usr/bin/env python3
"""Independent CPU-only assembly, manifest, count, distance and compute audit.

Does not import construction/assembly code, rerun matching, load models, or
read benchmark outcomes. Writes one exclusive-create validation artifact.
"""
from __future__ import annotations

from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import itertools
import json
import math
from pathlib import Path
import resource
import time
import traceback

import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve()
VERSION = SCRIPT.parents[1]
ROOT = SCRIPT.parents[5]
OLD = VERSION.parent / "data/COPD-V2-PREFLIGHT"
ATTEMPT = VERSION / "attempts/001_full_population"
OUTPUT = VERSION / "provenance/final_packaging_validation.json"
SEEDS = [104729, 130363, 155921]
MODELS = ["enhancer", "h3k27me3"]
CONFIGS = ["V2-A", "V2-B", "V2-C"]
FEATURES = ["gc_fraction", "atac_signal_percentile_max_train_only", "repeat_fraction_2001", "blacklist_input_any", "non_acgt_any"]
SCALES = np.array([.05, .20, .20])
checks, inputs, observations = [], {}, {}


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def record(path):
    path = Path(path)
    key = str(path.relative_to(ROOT))
    if key not in inputs:
        inputs[key] = {"bytes": path.stat().st_size, "sha256": sha(path)}
    return path


def table(path, **kwargs):
    return pd.read_csv(record(path), sep="\t", dtype=str, keep_default_na=False, **kwargs)


def obj(path):
    return json.loads(record(path).read_text())


def check(name, condition, detail=None):
    item = {"check": name, "status": "PASS" if bool(condition) else "FAIL"}
    if detail is not None:
        item["detail"] = detail
    checks.append(item)


def num(values):
    return pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)


def close(a, b, atol=2e-12):
    return math.isclose(float(a), float(b), rel_tol=2e-12, abs_tol=atol)


def one(frame, **kwargs):
    mask = np.ones(len(frame), dtype=bool)
    for column, value in kwargs.items():
        mask &= frame[column].eq(str(value)).to_numpy()
    result = frame.loc[mask]
    if len(result) != 1:
        raise AssertionError(f"Expected exactly one row {kwargs}; found {len(result)}")
    return result.iloc[0]


def select_scope(frame, kind, scope):
    return frame if kind == "overall" else frame[frame[kind].eq(scope)]


def audit():
    packaging = obj(VERSION / "provenance/configuration_manifest.json")
    compute = obj(VERSION / "provenance/compute_plan.json")
    registry = table(VERSION / "data/configuration_registry.tsv")
    count_table = table(VERSION / "results/class_counts.tsv")
    role_table = table(VERSION / "results/role_class_counts.tsv")
    runs = table(VERSION / "data/run_matrix.tsv")
    projection = table(VERSION / "results/compute_projection_by_run.tsv")
    check("six_configuration_registry_rows", len(registry) == 6 and len(registry.drop_duplicates(["configuration", "model"])) == 6)
    check("36_class_count_rows", len(count_table) == 36)
    check("six_rich_manifests_present", len(list((VERSION / "data/configurations").glob("*_interval_manifest.tsv.gz"))) == 6)
    check("36_BED_files_present", len(list((VERSION / "data/intervals").glob("*.bed.gz"))) == 36)
    check("assembly_no_model_execution", all(packaging[k] is False for k in ["training_started", "model_inference_performed", "phase_I_feature_extraction_performed", "benchmark_outcomes_read"]))

    old_columns = ["interval_id", "chrom", "core_start", "core_end", "partition", "sequence_sha256", "canonical_rc_sequence_sha256", "atac_lobes", "gc_fraction", "repeat_fraction_2001", "blacklist_bp_2001", "non_acgt_fraction"]
    old_columns += [f"v1_{m}_{label}" for m in MODELS for label in ["positive", "control"]]
    old_columns += [f"{m}_same_lobe_peak_support" for m in MODELS]
    master = table(OLD / "interval_features.tsv.gz", usecols=old_columns)
    roles = table(OLD / "interval_role_assignment.tsv.gz", usecols=["interval_id", "component_id", "validation_role"])
    atac = table(VERSION / "data/atac_normalized_features.tsv.gz", usecols=["interval_id", "atac_signal_percentile_max_train_only"])
    check("master_role_ATAC_ID_order", master.interval_id.tolist() == roles.interval_id.tolist() == atac.interval_id.tolist())
    master["component_id"] = roles.component_id
    master["validation_role"] = roles.validation_role
    master[FEATURES[1]] = atac[FEATURES[1]]
    master["blacklist_input_any"] = (num(master.blacklist_bp_2001) > 0).astype(int).astype(str)
    master["non_acgt_any"] = (num(master.non_acgt_fraction) > 0).astype(int).astype(str)
    master = master.set_index("interval_id", drop=False)
    del roles, atac
    manifests, sets, train_counts = {}, {}, {}
    membership = defaultdict(list)
    common_counts = []
    for model, config in itertools.product(MODELS, CONFIGS):
        key = f"{config}:{model}"
        reg = one(registry, configuration=config, model=model)
        path = ROOT / reg.manifest
        frame = table(path, usecols=list(master.columns) + ["configuration", "model", "label"])
        frame = frame.set_index("interval_id", drop=False)
        check(f"{key}:hash", sha(path) == reg.sha256)
        check(f"{key}:registry_authorization", reg.training_authorized == "False")
        check(f"{key}:registry_counts", len(frame) == int(reg.n_intervals) and int(reg.positive) == int(frame.label.eq("1").sum()) and int(reg.control) == int(frame.label.eq("0").sum()))
        check(f"{key}:unique_IDs_labels_identity", frame.index.is_unique and set(frame.label) == {"0", "1"} and frame.configuration.eq(config).all() and frame.model.eq(model).all())
        check(f"{key}:historical_fields_and_roles_preserved", frame[list(master.columns)].equals(master.loc[frame.index, list(master.columns)]))
        positives = set(frame.index[frame.label.eq("1")])
        controls = set(frame.index[frame.label.eq("0")])
        expected_mask = master[f"v1_{model}_positive"].eq("1")
        if config == "V2-C":
            expected_mask &= master[f"{model}_same_lobe_peak_support"].eq("1")
        check(f"{key}:complete_expected_positive_population", positives == set(master.index[expected_mask]))
        if config == "V2-A":
            check(f"{key}:original_controls", controls == set(master.index[master[f"v1_{model}_control"].eq("1")]))
        else:
            check(f"{key}:one_to_one_counts", len(positives) == len(controls) and not positives & controls)
        sets[(model, config, "all")] = set(frame.index)
        sets[(model, config, "positive")] = positives
        sets[(model, config, "control")] = controls
        for interval in frame.index:
            membership[interval].append(key)
        manifests[(model, config)] = frame
        train_counts[(model, config)] = int(frame.partition.eq("train").sum())
        for part, label in itertools.product(["train", "validation", "test"], ["positive", "control"]):
            label_num = "1" if label == "positive" else "0"
            sub = frame[frame.partition.eq(part) & frame.label.eq(label_num)]
            row = one(count_table, configuration=config, model=model, partition=part, **{"class": label})
            bed = table(ROOT / row.bed, header=None, names=["chrom", "core_start", "core_end", "interval_id"])
            check(f"{key}:{part}:{label}:BED_exact", bed.reset_index(drop=True).equals(sub[bed.columns].reset_index(drop=True)))
            check(f"{key}:{part}:{label}:BED_hash_counts", sha(ROOT / row.bed) == row.bed_sha256 and len(sub) == int(row.n_intervals) and sub.component_id.nunique() == int(row.n_components) and sub.canonical_rc_sequence_sha256.nunique() == int(row.n_unique_encoded_sequences))
        if config == "V2-C":
            panel_path = VERSION / f"data/evaluation/{model}_common_challenge_panel.tsv.gz"
            panel = table(panel_path, usecols=list(frame.columns)).set_index("interval_id", drop=False)
            expected = frame[frame.validation_role.isin(["selection", "calibration", "test"])]
            check(f"{key}:common_panel_exact_C_rows", panel.equals(expected))
            check(f"{key}:common_panel_roles", set(panel.validation_role) == {"selection", "calibration", "test"})
            for role, sub in panel.groupby("validation_role"):
                pos = sub[sub.label.eq("1")]; neg = sub[sub.label.eq("0")]
                common_counts.append({"model": model, "role": role, "positive": len(pos), "control": len(neg), "components": sub.component_id.nunique(), "control_components": neg.component_id.nunique()})
                if role == "selection":
                    check(f"{model}:selection_minima", len(pos) >= 100 and len(neg) >= 200 and sub.component_id.nunique() >= 30)
                if role == "calibration":
                    check(f"{model}:calibration_minima", len(neg) >= 200 and neg.component_id.nunique() >= 30)
                row = one(role_table, configuration="COMMON_C_TASK", model=model, role=role)
                check(f"{model}:{role}:common_count_table", int(row.positive) == len(pos) and int(row.control) == len(neg) and int(row.components) == sub.component_id.nunique() and int(row.control_components) == neg.component_id.nunique())
        print(f"Audited {key} manifests and BED files", flush=True)

    for row in role_table.itertuples(index=False):
        if row.configuration == "COMMON_C_TASK":
            frame = manifests[(row.model, "V2-C")]
        else:
            frame = manifests[(row.model, row.configuration)]
        sub = frame[frame.validation_role.eq(row.role)]
        pos = sub[sub.label.eq("1")]; neg = sub[sub.label.eq("0")]
        correct = int(row.positive) == len(pos) and int(row.control) == len(neg)
        correct &= int(row.components) == sub.component_id.nunique() and int(row.positive_components) == pos.component_id.nunique() and int(row.control_components) == neg.component_id.nunique()
        correct &= int(row.positive_unique_encoded_sequences) == pos.canonical_rc_sequence_sha256.nunique() and int(row.control_unique_encoded_sequences) == neg.canonical_rc_sequence_sha256.nunique()
        correct &= close(row.control_to_positive_ratio, len(neg) / len(pos))
        check(f"role_counts:{row.configuration}:{row.model}:{row.role}", correct)

    overlap = table(VERSION / "results/configuration_overlap.tsv")
    for row in overlap.itertuples(index=False):
        allowed = set(master.index if row.scope_type == "overall" else master.index[master.partition.eq(row.scope)])
        a = sets[(row.model, row.configuration_1, row._5)] & allowed
        b = sets[(row.model, row.configuration_2, row._5)] & allowed
        correct = [len(a), len(b), len(a & b), len(a - b), len(b - a)] == [int(row.n_1), int(row.n_2), int(row.shared), int(row.only_1), int(row.only_2)]
        check(f"overlap:{row.model}:{row.configuration_1}:{row.configuration_2}:{row.scope}:{row._5}", correct and close(row.Jaccard, len(a & b) / len(a | b)))

    selected_ids = set().union(*(sets[(m, c, "all")] for m, c in itertools.product(MODELS, CONFIGS)))
    selected = master.loc[sorted(selected_ids)]
    encoded_count = selected.canonical_rc_sequence_sha256.nunique()
    check("selected_union_counts", len(selected) == packaging["selected_unique_intervals"] == compute["selected_unique_intervals"] and encoded_count == packaging["unique_selected_encoded_sequences"] == compute["unique_selected_encoded_sequences"])
    duplicate_table = table(VERSION / "results/selected_sequence_duplicate_audit.tsv.gz")
    for kind, column in [("raw_exact", "sequence_sha256"), ("encoded_forward_RC", "canonical_rc_sequence_sha256")]:
        dup = selected[selected[column].duplicated(keep=False)]
        groups = dup.groupby(column)
        actual = duplicate_table[duplicate_table.identity_kind.eq(kind)]
        check(f"{kind}:duplicate_rows_complete", set(actual.interval_id) == set(dup.interval_id) and len(actual) == len(dup))
        check(f"{kind}:duplicate_counts_and_no_leakage", groups.ngroups == packaging["duplicate_groups"][kind] and int((groups.partition.nunique() > 1).sum()) == packaging["cross_partition_selected_duplicate_groups"][kind] == 0 and int((groups.validation_role.nunique() > 1).sum()) == packaging["cross_validation_role_selected_duplicate_groups"][kind] == 0)
        exact = True
        for row in actual.itertuples(index=False):
            sub = dup[dup[column].eq(row.identity_sha256)]
            exact &= int(row.group_rows) == len(sub) and row.interval_id in set(sub.interval_id)
            exact &= row.group_partitions == ";".join(sorted(sub.partition.unique())) and row.group_validation_roles == ";".join(sorted(sub.validation_role.unique()))
            exact &= set(row.configuration_model_membership.split(";")) == set(membership[row.interval_id])
        check(f"{kind}:duplicate_membership_metadata", exact)

    wanted_runs = {(c, m, str(seed)) for c, m, seed in itertools.product(CONFIGS, MODELS, SEEDS)}
    check("18_exact_runs_seeds_unique", len(runs) == 18 and set(zip(runs.configuration, runs.model, runs.seed)) == wanted_runs and runs.run_id.is_unique)
    check("18_runs_not_authorized_or_started", runs.status.eq("NOT_AUTHORIZED_NOT_STARTED").all() and projection.status.eq("NOT_AUTHORIZED_NOT_STARTED").all())
    totals = 0.0
    references = {"enhancer": (464262, 2659), "h3k27me3": (78165, 454)}
    # These reference times are immutable original 50-epoch resource observations,
    # not any new model execution; verify A counts and prior frozen cost table.
    old_cost = table(VERSION.parent / "results/COPD-V2-PREFLIGHT_compute_projection_by_run.tsv")
    for model, (baseline_n, seconds) in references.items():
        check(f"{model}:historical_reference_count", train_counts[(model, "V2-A")] == baseline_n)
        row = one(old_cost, run_id=f"V2-A_{model}_seed104729")
        check(f"{model}:historical_reference_seconds", close(float(row.baseline_scaled_fit_GPU_hours) * 3600, seconds, atol=1e-6))
    for row in runs.itertuples(index=False):
        reg = one(registry, configuration=row.configuration, model=row.model)
        n = train_counts[(row.model, row.configuration)]
        baseline_n, seconds = references[row.model]
        expected_cost = seconds * n / baseline_n / 3600
        cost = one(projection, run_id=row.run_id)
        check(f"{row.run_id}:manifest_counts_parameters", row.manifest == reg.manifest and row.manifest_sha256 == reg.sha256 and int(row.training_intervals) == n and int(row.epochs_max) == 50 and int(row.batch_size) == 256 and int(row.patience) == 15 and close(row.learning_rate, .001) and row.optimizer == "Adadelta")
        check(f"{row.run_id}:independent_compute", int(cost.training_intervals) == n and int(cost.historical_reference_training_intervals) == baseline_n and int(cost.historical_reference_fit_seconds) == seconds and close(cost.baseline_scaled_fit_GPU_hours, expected_cost))
        totals += expected_cost
    cache = int(encoded_count) * 2 * 4560 * 4
    check("compute_18_fit_total", compute["n_fits"] == 18 and close(compute["scaled_fit_only_GPU_hours"], totals))
    check("compute_dual_orientation_float32_cache", compute["future_two_orientation_float32_cache_bytes_per_encoded_sequence"] == 2 * 4560 * 4 and compute["future_feature_cache_bytes_for_all_unique_selected_sequences_both_orientations"] == packaging["future_two_orientation_feature_cache_bytes"] == cache and close(compute["future_feature_cache_GiB"], cache / 2**30))
    check("compute_planning_not_authorization", compute["training_authorized"] is False and compute["training_started"] is False and compute["model_inference_performed"] is False and compute["training_seeds"] == SEEDS and compute["sum_requested_walltime_caps_GPU_hours"] == 18 * 4 and compute["conservative_total_GPU_hours"] == [12, 24])
    observations.update(selected_unique_intervals=len(selected), unique_selected_encoded_sequences=int(encoded_count), scaled_fit_only_GPU_hours=totals, dual_orientation_cache_bytes=cache, common_C_panel_counts=common_counts)

    retained = table(VERSION / "results/positive_retention.tsv")
    attempt_retained = table(ATTEMPT / "positive_retention.tsv")
    check("retention_table_preserved", retained[list(attempt_retained.columns)].equals(attempt_retained))
    for row in retained.itertuples(index=False):
        ids = sets[(row.model, row.configuration, "positive")]
        sub = master.loc[sorted(ids)]
        n = int(sub.partition.eq(row.partition).sum())
        check(f"retention:{row.configuration}:{row.model}:{row.partition}", int(row.original_positive_n) == int(row.retained_positive_n) == int(row.unique_control_n) == n and int(row.excluded_positive_n) == 0 and float(row.retention) == 1)
    exclusions = table(ATTEMPT / "positive_common_support_exclusions.tsv.gz")
    check("no_positive_common_support_exclusions", len(exclusions) == 0)
    changes = table(VERSION / "results/pure_label_membership_change.tsv")
    for row in changes.itertuples(index=False):
        allowed = set(master.index if row.scope_type == "overall" else master.index[master.partition.eq(row.scope)])
        b = sets[(row.model, "V2-B", "positive")] & allowed
        c = sets[(row.model, "V2-C", "positive")] & allowed
        check(f"pure_label_change:{row.model}:{row.scope}", c <= b and int(row.complete_V1_positive_n) == len(b) and int(row.same_lobe_positive_n) == len(c) and int(row.cross_lobe_only_removed_n) == len(b - c) and close(row.fraction_removed, len(b - c) / len(b)) and int(row.common_support_restriction_exclusions_B) == int(row.common_support_restriction_exclusions_C) == 0)

    distance_summary = table(VERSION / "results/matched_pair_distance_summary.tsv")
    for stage, filename in [("initial", "initial_pairs.tsv.gz"), ("final", "control_matching_pairs.tsv.gz")]:
        pairs = table(ATTEMPT / filename)
        for (config, model), subset in pairs.groupby(["configuration", "model"]):
            a = master.loc[subset.positive_id].reset_index(drop=True)
            b = master.loc[subset.control_id].reset_index(drop=True)
            aa = np.column_stack([num(a[f]) for f in FEATURES[:3]])
            bb = np.column_stack([num(b[f]) for f in FEATURES[:3]])
            distances = np.max(np.abs(aa / SCALES - bb / SCALES), axis=1)
            differences = np.abs(aa - bb)
            check(f"{stage}:{config}:{model}:complete_distance_table", np.allclose(distances, num(subset[f"{stage}_normalized_chebyshev_distance"]), rtol=0, atol=2e-12))
            if stage == "final":
                check(f"{config}:{model}:pair_manifest_membership", set(a.interval_id) == sets[(model, config, "positive")] and set(b.interval_id) == sets[(model, config, "control")] and b.interval_id.is_unique)
            metric_values = {"normalized_chebyshev_distance": (distances, 1.0)}
            metric_values.update({f"{f}_absolute_difference": (differences[:, j], SCALES[j]) for j, f in enumerate(FEATURES[:3])})
            rows = distance_summary[distance_summary.stage.eq(stage) & distance_summary.configuration.eq(config) & distance_summary.model.eq(model)]
            for row in rows.itertuples(index=False):
                mask = np.ones(len(a), dtype=bool) if row.scope_type == "overall" else a[row.scope_type].eq(row.scope).to_numpy()
                values, scale = metric_values[row.metric]
                values = values[mask]
                quantiles = np.quantile(values, [0, .01, .05, .25, .5, .75, .95, .99, 1])
                observed = [row.minimum, row.q01, row.q05, row.q25, row.median, row.q75, row.q95, row.q99, row.maximum]
                correct = int(row.n_total) == int(row.n_finite) == len(values) and int(row.n_missing) == 0
                correct &= all(close(x, y) for x, y in zip(observed, quantiles)) and close(row.mean, np.mean(values)) and close(row.standard_deviation, np.std(values, ddof=1))
                correct &= int(row.above_legacy_scale_n) == int((values > scale).sum()) and close(row.fraction_above_legacy_scale, (values > scale).mean()) and close(row.legacy_scale_not_current_caliper, scale)
                cross_chrom = int((a.chrom.to_numpy()[mask] != b.chrom.to_numpy()[mask]).sum())
                cross_lobe = int((a.atac_lobes.to_numpy()[mask] != b.atac_lobes.to_numpy()[mask]).sum())
                correct &= int(row.cross_chromosome_n) == cross_chrom and int(row.cross_ATAC_signature_n) == cross_lobe and close(row.cross_chromosome_fraction, cross_chrom / len(values)) and close(row.cross_ATAC_signature_fraction, cross_lobe / len(values))
                check(f"distance_summary:{stage}:{config}:{model}:{row.scope}:{row.metric}", correct)
    check("all_224_distance_summary_rows_present", len(distance_summary) == 2 * 2 * 2 * 7 * 4)

    balance = table(ATTEMPT / "covariate_balance.tsv")
    for (config, model, scope_type, scope), sub in balance[balance.stage.eq("final") & balance.kind.eq("SMD")].groupby(["configuration", "model", "scope_type", "scope"]):
        frame = select_scope(manifests[(model, config)], scope_type, scope)
        a = frame[frame.label.eq("1")]; b = frame[frame.label.eq("0")]
        aa = np.column_stack([num(a[f]) for f in FEATURES]); bb = np.column_stack([num(b[f]) for f in FEATURES])
        valid = len(a) >= 2 and len(b) >= 2 and np.isfinite(aa).all() and np.isfinite(bb).all()
        variances = (np.var(aa, axis=0, ddof=1) + np.var(bb, axis=0, ddof=1)) / 2
        valid &= np.isfinite(variances).all()
        check(f"SMD_validity_prerequisites:{config}:{model}:{scope_type}:{scope}", valid, {"positive_n": len(a), "control_n": len(b)})
        delta = aa.mean(axis=0) - bb.mean(axis=0)
        actual = np.zeros(5)
        np.divide(delta, np.sqrt(variances), out=actual, where=variances > 0)
        actual[(variances == 0) & (delta != 0)] = np.inf
        check(f"independent_SMD:{config}:{model}:{scope_type}:{scope}", np.isfinite(actual).all() and (np.abs(actual) <= .10).all() and all(close(one(sub, variable=f).value, value) for f, value in zip(FEATURES, actual)))
    method = obj(VERSION / "provenance/matching_method_review.json")
    check("method_review_hashes_current", method["implementation"]["sha256"] == sha(ROOT / method["implementation"]["path"]) and method["specification"]["sha256"] == sha(ROOT / method["specification"]["path"]))
    check("method_tie_fixtures_pass", method["synthetic_tie_audit"]["total_cases"] == 142 and method["synthetic_tie_audit"]["discrepancies"] == 0)
    observations["method_review_SMD_caveat"] = "All final diagnostic scopes independently checked for >=2 observations per class, finite five-feature arrays and finite pooled variance; undefined-SMD counterexample does not apply to these artifacts."

    selection = obj(VERSION / "provenance/selection_adequacy_review.json")
    check("selection_95_invariant_checks", selection["check_count"] == selection["pass_count"] == 95 and selection["failure_count"] == 0 and all(x["status"] == "PASS" for x in selection["checks"]))
    check("selection_review_current_spec_and_script", selection["specification_sha256"] == sha(ROOT / selection["specification_path"]) and selection["script_sha256"] == sha(ROOT / selection["script_path"]))
    normalization = obj(VERSION / "provenance/atac_normalization_manifest.json")
    norm_qc = table(VERSION / "results/atac_normalization_qc.tsv")
    norm_independent = obj(VERSION / "results/atac_independent_validation.json")
    check("normalization_23_checks", normalization["qc_pass"] == len(norm_qc) == 23 and normalization["qc_fail"] == 0 and norm_qc.status.eq("PASS").all())
    check("normalization_independent_PASS", norm_independent["status"] == "PASS" and norm_independent["heldout_rows_in_fit"] == 0 and norm_independent["interval_rows_validated"] == len(master) and norm_independent["model_or_benchmark_information_read"] is False)
    check("normalization_review_current_outputs", norm_independent["cdf_sha256"] == sha(record(VERSION / "data/atac_training_cdf.tsv.gz")) and norm_independent["sidecar_sha256"] == sha(VERSION / "data/atac_normalized_features.tsv.gz"))
    for entry in packaging["inputs"] + packaging["common_C_panels"] + normalization["outputs"]:
        path = record(ROOT / entry["path"])
        check(f"registered_hash:{entry['path']}", sha(path) == entry["sha256"] and path.stat().st_size == entry["bytes"])


def main():
    if OUTPUT.exists():
        raise RuntimeError("Refusing to overwrite an existing packaging validation, including failures")
    started = time.time()
    try:
        audit()
    except Exception as exc:
        check("unexpected_validation_exception", False, {"type": type(exc).__name__, "message": str(exc), "traceback": traceback.format_exc()})
    failures = sum(x["status"] != "PASS" for x in checks)
    result = {"status": "PASS" if not failures else "FAIL", "failures": failures, "check_count": len(checks), "checks": checks, "input_hashes": inputs, "observations": observations, "script_sha256": sha(SCRIPT), "completed_utc": datetime.now(timezone.utc).isoformat(), "wall_seconds": time.time() - started, "peak_RSS_KiB": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss, "training_authorized": False, "training_or_model_inference_run": False, "benchmark_outcomes_read": False, "limitation": "Packaging/count/compute audit only; complete deterministic matching replay and final historical freeze preservation have independent audits. This pass does not authorize GPU work or establish model adequacy."}
    with OUTPUT.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    print(json.dumps({"status": result["status"], "failures": failures, "checks": len(checks), "output": str(OUTPUT.relative_to(ROOT))}), flush=True)
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
