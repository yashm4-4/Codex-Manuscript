#!/usr/bin/env python3
"""Independently verify ancillary CPU-only preflight summaries, not models."""
from __future__ import annotations
import collections
import csv
import datetime as dt
import gzip
import hashlib
import json
from pathlib import Path
import re
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
BASE = ROOT / "diseases/COPD/07_gap_closure"
DATA = BASE / "data/COPD-V2-PREFLIGHT"
PROV = BASE / "provenance"
OUT = BASE / "results"
P = "COPD-V2-PREFLIGHT"
CHECKS = []
VARIABLES = ["gc_fraction", "atac_signal_percentile_max", "repeat_fraction_2001", "blacklist_bp_2001", "promoter_bp_2001", "non_acgt_fraction"]


def rows(path):
    assert "COPD-V2-BENCH" not in str(path)
    with (gzip.open(path, "rt") if str(path).endswith(".gz") else Path(path).open()) as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def check(name, passed, details=""):
    CHECKS.append({"check": name, "status": "PASS" if passed else "FAIL", "details": details})


def approximately(actual, expected):
    if not np.isfinite(expected):
        return actual in ("", None, "nan")
    return actual not in ("", None) and abs(float(actual)-float(expected)) <= 1e-8*max(1.0, abs(float(expected)))


def main():
    target = PROV / f"{P}_ancillary_validation.json"
    if target.exists():
        raise RuntimeError("Refusing to overwrite an existing independent ancillary audit")
    covariates = {}
    contexts = collections.Counter()
    class_sizes = {}
    registry = list(rows(DATA / "configuration_registry.tsv"))
    for entry in registry:
        config, model = entry["configuration"], entry["model"]
        samples = collections.defaultdict(list)
        totals = collections.Counter()
        for row in rows(ROOT / entry["manifest"]):
            part = row["partition"]
            label = "positive" if row["label"] == "1" else "control"
            group = (config, model, part, label)
            totals[group] += 1
            for variable in VARIABLES:
                if row[variable] not in ("", "nan"):
                    samples[(*group, variable)].append(float(row[variable]))
            for variable in ("atac_lobes", f"{model}_same_lobe_support_lobes"):
                value = row[variable] or "NO_EXACT_ANCHOR_OR_SUPPORT"
                contexts[(*group, variable, value)] += 1
        class_sizes.update(totals)
        for group, total in totals.items():
            for variable in VARIABLES:
                values = np.array(samples[(*group, variable)], dtype=float)
                n = len(values)
                q = np.quantile(values, [.05, .25, .50, .75, .95]) if n else [float("nan")]*5
                covariates[(*group, variable)] = {
                    "total_n": total, "observed_n": n, "missing_n": total-n,
                    "mean": float(np.mean(values)) if n else float("nan"),
                    "sd_sample": float(np.std(values, ddof=1)) if n > 1 else float("nan"),
                    "min": float(np.min(values)) if n else float("nan"),
                    "q05": q[0], "q25": q[1], "median": q[2], "q75": q[3], "q95": q[4],
                    "max": float(np.max(values)) if n else float("nan")
                }
        print(f"Independently summarized {config} {model}", flush=True)
    cov_path = OUT / f"{P}_covariate_distributions.tsv"
    reported_rows = list(rows(cov_path))
    key = lambda r: tuple(r[c] for c in ("configuration", "model", "partition", "class", "variable"))
    reported = {key(r): r for r in reported_rows}
    check("covariates:complete_unique_rows", len(reported) == len(reported_rows) and set(reported) == set(covariates))
    for group, expected in covariates.items():
        actual = reported.get(group, {})
        incorrect = [field for field, value in expected.items() if not approximately(actual.get(field), value)]
        check("covariates:" + ":".join(group), not incorrect, {"incorrect_fields": incorrect})
    available_rows = list(rows(OUT / f"{P}_ATAC_strength_available_denominators.tsv"))
    available = {key(r): r for r in available_rows}
    expected_available = {k: r for k, r in reported.items() if k[-1] == "atac_signal_percentile_max"}
    check("ATAC_denominators:exact_verified_covariate_subset", len(available_rows) == len(available) and available == expected_available)
    context_rows = list(rows(OUT / f"{P}_source_context_distribution.tsv"))
    context_map = {tuple(r[c] for c in ("configuration", "model", "partition", "class", "variable", "value")): r for r in context_rows}
    check("source_context:complete_unique_categories", len(context_map) == len(context_rows) and set(context_map) == set(contexts))
    for group, count in contexts.items():
        row = context_map.get(group, {})
        check("source_context:" + ":".join(group), row.get("n_intervals") == str(count) and row.get("class_n") == str(class_sizes[group[:4]]))
    matched = collections.Counter((r["model"], r["chrom"], r["validation_role"], r["atac_lobes"]) for r in rows(DATA / "control_matching_pairs.tsv.gz"))
    matched_rows = list(rows(OUT / f"{P}_matched_control_counts.tsv"))
    matched_actual = {tuple(r[c] for c in ("model", "chrom", "validation_role", "atac_lobes")): int(r["matched_controls"]) for r in matched_rows}
    check("matched_counts:all_strata_and_counts", len(matched_rows) == len(matched_actual) and matched_actual == dict(matched))
    unmatched = collections.Counter((r["model"], r["positive_id"].split(":")[0], r["reason"]) for r in rows(DATA / "unmatched_positive_targets.tsv.gz"))
    unmatched_rows = list(rows(OUT / f"{P}_unmatched_target_counts.tsv"))
    unmatched_actual = {tuple(r[c] for c in ("model", "chrom", "reason")): int(r["unmatched_positives"]) for r in unmatched_rows}
    check("unmatched_counts:all_strata_and_counts", len(unmatched_rows) == len(unmatched_actual) and unmatched_actual == dict(unmatched))
    historical_seconds = {}
    for model, word in (("enhancer", "Enhancer"), ("h3k27me3", "Silencer")):
        path = ROOT / f"diseases/COPD/04_modeling/trednet/logs/COPD_SevereEmphysema_Lung_{word}_DHS_x2.train.log"
        times = [int(t) for t in re.findall(r"\d+/\d+ - (\d+)s - ", path.read_text())]
        check(f"compute:{model}:fifty_original_epochs", len(times) == 50)
        historical_seconds[model] = sum(times)
    baseline_n = {"enhancer": 464262, "h3k27me3": 78165}
    run_matrix = {r["run_id"]: r for r in rows(DATA / "run_matrix.tsv")}
    projection_rows = list(rows(OUT / f"{P}_compute_projection_by_run.tsv"))
    projection = {r["run_id"]: r for r in projection_rows}
    check("compute:exact_eighteen_run_set", len(projection_rows) == len(projection) == len(run_matrix) == 18 and set(projection) == set(run_matrix))
    expected_total = 0.0
    for run_id, run in run_matrix.items():
        expected = historical_seconds[run["model"]] * (int(run["training_intervals"])/baseline_n[run["model"]])/3600
        expected_total += expected
        actual = projection[run_id]
        check("compute:scaled_fit:" + run_id, approximately(actual["baseline_scaled_fit_GPU_hours"], expected))
        fields = ("training_intervals", "gpus", "cpus", "memory_GiB", "walltime_cap_hours", "status")
        check("compute:request_metadata:" + run_id, all(actual[f] == run[f] for f in fields))
    plan = json.loads((PROV / f"{P}_compute_plan.json").read_text())
    config = json.loads((PROV / f"{P}_configuration_manifest.json").read_text())
    feature = json.loads((PROV / f"{P}_feature_manifest.json").read_text())
    check("compute:total_scaled_fit_hours", approximately(plan["scaled_fit_only_GPU_hours"], expected_total))
    check("compute:cache_bytes", plan["future_feature_cache_bytes_for_all_unique_selected_sequences_both_orientations"] == config["future_two_orientation_feature_cache_bytes"])
    check("compute:cache_GiB", approximately(plan["future_feature_cache_GiB"], config["future_two_orientation_feature_cache_bytes"]/2**30))
    check("compute:actual_feature_seconds", plan["actual_CPU_feature_stage_seconds"] == feature["elapsed_seconds"])
    check("compute:actual_feature_RSS", plan["actual_CPU_feature_peak_RSS_KiB"] == feature["peak_rss_kib"])
    check("compute:actual_configuration_seconds", plan["actual_CPU_configuration_stage_seconds"] == config["elapsed_seconds"])
    check("compute:walltime_cap_sum", plan["sum_requested_walltime_caps_GPU_hours"] == sum(int(r["walltime_cap_hours"])*int(r["gpus"]) for r in run_matrix.values()))
    check("compute:not_authorized_not_started", plan["training_authorized"] is False and plan["training_started"] is False)
    check("runtime:no_TensorFlow_or_Keras", not any(m == "tensorflow" or m.startswith("tensorflow.") or m == "keras" or m.startswith("keras.") for m in sys.modules))
    inputs = [cov_path, OUT / f"{P}_ATAC_strength_available_denominators.tsv", OUT / f"{P}_source_context_distribution.tsv", OUT / f"{P}_matched_control_counts.tsv", OUT / f"{P}_unmatched_target_counts.tsv", OUT / f"{P}_compute_projection_by_run.tsv", PROV / f"{P}_compute_plan.json", DATA / "configuration_registry.tsv", DATA / "run_matrix.tsv", DATA / "control_matching_pairs.tsv.gz", DATA / "unmatched_positive_targets.tsv.gz"] + [ROOT / r["manifest"] for r in registry]
    failures = [r for r in CHECKS if r["status"] == "FAIL"]
    result = {"module": P, "completed_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
              "status": "PASS" if not failures else "FAIL", "n_checks": len(CHECKS), "n_failures": len(failures),
              "checks": CHECKS, "failure_details": failures, "scaled_fit_only_GPU_hours": expected_total,
              "no_model_execution": True, "no_benchmark_contents_read": True,
              "validator": {"path": str(Path(__file__).relative_to(ROOT)), "sha256": sha(__file__)},
              "inputs": [{"path": str(p.relative_to(ROOT)), "sha256": sha(p)} for p in inputs]}
    with target.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    print(json.dumps({k: result[k] for k in ("status", "n_checks", "n_failures", "scaled_fit_only_GPU_hours")}), flush=True)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
