#!/usr/bin/env python3
"""CPU-only complete rich-manifest pair-link metadata closure.

Does not rematch, run a model, or modify any existing output. Complements the
independent 625-check packaging audit without replacing its frozen result.
"""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import itertools
import json
from pathlib import Path
import time
import traceback

import numpy as np
import pandas as pd

SCRIPT = Path(__file__).resolve()
VERSION = SCRIPT.parents[1]
ROOT = SCRIPT.parents[5]
OUTPUT = VERSION / "provenance/supplementary_pair_metadata_validation.json"
LINKS = ["matched_positive_id", "matched_control_id", "matched_final_normalized_chebyshev_distance"]
FEATURES = ["gc_fraction", "atac_signal_percentile_max_train_only", "repeat_fraction_2001"]
SCALES = np.array([.05, .20, .20])
checks, inputs, populations = [], {}, []


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read(path, **kwargs):
    inputs[str(path.relative_to(ROOT))] = {"sha256": sha(path), "bytes": path.stat().st_size}
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, **kwargs)


def check(name, passed, details=None):
    row = {"check": name, "status": "PASS" if bool(passed) else "FAIL"}
    if details is not None:
        row["details"] = details
    checks.append(row)


def audit():
    ledger = read(VERSION / "attempts/001_full_population/control_matching_pairs.tsv.gz")
    for config, model in itertools.product(["V2-A", "V2-B", "V2-C"], ["enhancer", "h3k27me3"]):
        key = f"{config}:{model}"
        path = VERSION / f"data/configurations/{config}_{model}_interval_manifest.tsv.gz"
        frame = read(path, usecols=["interval_id", "configuration", "model", "label"] + LINKS + FEATURES).set_index("interval_id", drop=False)
        check(f"{key}:unique_manifest_identity", frame.index.is_unique and frame.configuration.eq(config).all() and frame.model.eq(model).all())
        if config == "V2-A":
            check(f"{key}:all_historical_link_fields_empty", frame[LINKS].eq("").to_numpy().all())
            populations.append({"configuration": config, "model": model, "manifest_rows": len(frame), "historical_empty_link_rows": len(frame), "matched_pairs": 0})
            continue
        pairs = ledger[ledger.configuration.eq(config) & ledger.model.eq(model)]
        check(f"{key}:unique_one_to_one_pair_endpoints", pairs.positive_id.is_unique and pairs.control_id.is_unique and not set(pairs.positive_id) & set(pairs.control_id))
        check(f"{key}:complete_both_endpoints", len(frame) == 2 * len(pairs) and set(frame.index) == set(pairs.positive_id) | set(pairs.control_id))
        expected = pairs[["positive_id", "control_id", "final_normalized_chebyshev_distance"]].rename(columns=dict(zip(["positive_id", "control_id", "final_normalized_chebyshev_distance"], LINKS)))
        for endpoint, label in [("positive_id", "1"), ("control_id", "0")]:
            actual = frame.loc[pairs[endpoint]]
            check(f"{key}:{endpoint}:correct_class", actual.label.eq(label).all())
            for column in LINKS:
                check(f"{key}:{endpoint}:{column}:exact_ledger_strings", actual[column].tolist() == expected[column].tolist())
        positive = frame.loc[pairs.positive_id]
        control = frame.loc[pairs.control_id]
        a = np.column_stack([pd.to_numeric(positive[c], errors="coerce").to_numpy(float) for c in FEATURES])
        b = np.column_stack([pd.to_numeric(control[c], errors="coerce").to_numpy(float) for c in FEATURES])
        check(f"{key}:finite_distance_covariates", np.isfinite(a).all() and np.isfinite(b).all())
        distances = np.max(np.abs(a / SCALES - b / SCALES), axis=1)
        reported = pd.to_numeric(pairs.final_normalized_chebyshev_distance, errors="coerce").to_numpy(float)
        error = float(np.max(np.abs(distances - reported)))
        check(f"{key}:distances_independently_recomputed", np.isfinite(reported).all() and np.allclose(distances, reported, rtol=0, atol=1e-12), {"maximum_absolute_error": error})
        populations.append({"configuration": config, "model": model, "manifest_rows": len(frame), "matched_pairs": len(pairs), "both_endpoint_link_fields_verified": True, "distance_maximum_absolute_error": error})
    check("all_complete_pairs_covered_once", sum(x["matched_pairs"] for x in populations) == len(ledger))


def main():
    if OUTPUT.exists():
        raise RuntimeError("Refusing overwrite of prior metadata validation, including failures")
    started = time.time()
    try:
        audit()
    except Exception as exc:
        check("unexpected_validation_exception", False, {"type": type(exc).__name__, "message": str(exc), "traceback": traceback.format_exc()})
    failures = sum(x["status"] != "PASS" for x in checks)
    result = {"status": "PASS" if failures == 0 else "FAIL", "failures": failures, "check_count": len(checks), "checks": checks, "input_hashes": inputs, "populations": populations, "script_sha256": sha(SCRIPT), "completed_utc": datetime.now(timezone.utc).isoformat(), "wall_seconds": time.time() - started, "scope": "Complete rich-manifest pair-link strings and independently recomputed Chebyshev distances; historical A rows have no invented pair links.", "previous_packaging_validation_modified": False, "matching_rerun": False, "model_or_benchmark_execution": False, "training_authorized": False}
    with OUTPUT.open("x") as handle:
        json.dump(result, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")
    print(json.dumps({"status": result["status"], "checks": len(checks), "failures": failures, "paired_rows": sum(x["manifest_rows"] for x in populations if x["matched_pairs"]), "matched_pairs": sum(x["matched_pairs"] for x in populations)}), flush=True)
    return bool(failures)


if __name__ == "__main__":
    raise SystemExit(main())
