#!/usr/bin/env python3
"""Calibrate frozen C ensembles only AFTER both contexts pass chr7 adequacy.

No model execution. No positive prediction is loaded or used. The exact
conservative threshold may exceed one, and is never clipped or optimized.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np

from evaluate_chr7 import (MODELS, SEEDS, checked_reference, load_predictions,
                           probability_array, read_json, read_rows, sha256,
                           validate_common_rows, validate_prerequisites,
                           write_json, write_rows)


def conservative_threshold(scores, components, require_minima=True):
    values = probability_array(scores, "calibration negatives")
    if values.ndim != 1 or len(values) != len(components) or not len(values):
        raise ValueError("Calibration arrays must be nonempty and aligned")
    n = len(values)
    n_components = len(set(components))
    if require_minima and (n < 200 or n_components < 30):
        raise ValueError("Calibration safeguard failed: require>=200controls and>=30components")
    k = math.floor(0.05 * n)
    descending = np.sort(values)[::-1]
    boundary = float(descending[k])
    threshold = float(np.nextafter(np.float64(boundary), np.float64(np.inf)))
    hex_value = threshold.hex()
    decimal = format(threshold, ".17g")
    if float.fromhex(hex_value) != threshold or float(decimal) != threshold:
        raise ValueError("Threshold serialization failed exact float64 round-trip")
    calls = int(np.sum(values >= threshold))
    ties = int(np.sum(values == boundary))
    above = int(np.sum(values > boundary))
    if calls > k or calls != above:
        raise ValueError("Frozen threshold failed conservative row-FPR/tie contract")
    return {"n_controls": n, "n_control_components": n_components,
            "target_empirical_FPR": 0.05, "maximum_called_controls": k,
            "zero_based_boundary_index": k, "sorted_boundary_score": boundary,
            "boundary_hex": boundary.hex(), "boundary_decimal_17g": format(boundary, ".17g"),
            "threshold_hex": hex_value, "threshold_decimal_17g": decimal,
            "threshold_float64": threshold, "threshold_above_one": threshold > 1,
            "n_at_boundary": ties, "n_strictly_above_boundary": above,
            "n_called_controls": calls, "observed_calibration_row_FPR": calls / n,
            "call_rule": "score>=threshold", "tie_policy": "all boundary ties excluded; no tie splitting",
            "serialization_roundtrip": "PASS", "empirical_FPR_gate": "PASS",
            "population_FPR_guarantee": False, "variant_FPR_guarantee": False,
            "allele_delta_threshold": None}


def validate_calibration_prerequisites(manifest):
    decision = read_json(checked_reference(manifest["adequacy_decision"]))
    if decision.get("both_C_contexts_pass") is not True or set(decision.get("model_contexts", {})) != set(MODELS):
        raise ValueError("STOP: both C contexts did not pass; no calibration or A/B fallback")
    if any(decision["model_contexts"][model]["status"] != "PASS" or decision["model_contexts"][model]["failed_or_inconclusive_gates"] for model in MODELS):
        raise ValueError("Adequacy decision contains failed or inconclusive conditions")
    evaluation_path = checked_reference(manifest["evaluation_input"])
    if sha256(evaluation_path) != decision["input_manifest_sha256"]:
        raise ValueError("Adequacy/evaluation-input binding failed")
    evaluation = read_json(evaluation_path)
    validate_prerequisites(evaluation)
    if manifest["specification"] != evaluation["specification"]:
        raise ValueError("Calibration/evaluation frozen specification differs")
    release = read_json(checked_reference(manifest["release_freeze"]))
    if release["status"] != "C_READY_FOR_CALIBRATION_FROZEN":
        raise ValueError("C ensemble/checkpoint/adequacy freeze missing")
    if release["adequacy_decision"] != manifest["adequacy_decision"]:
        raise ValueError("Release freeze not bound to the passed adequacy decision")
    checked_reference(release["symmetric_inference_implementation"])
    for supporting in release.get("symmetric_inference_supporting_files", []):
        checked_reference(supporting)
    ensemble = release["ensemble"]
    if (tuple(ensemble["seeds"]) != SEEDS or ensemble["reduction_dtype"] != "float64"
            or ensemble["seed_aggregation"] != "equal-weight"
            or ensemble["symmetric_inference"] != "(p_forward+p_reverse_complement)/2"):
        raise ValueError("Frozen ensemble rule differs")
    expected = {(r["model"], int(r["seed"])): r for r in evaluation["prerequisites"]["checkpoints"] if r["configuration"] == "V2-C"}
    observed = set()
    for record in release["checkpoints"]:
        key = (record["model"], int(record["seed"]))
        if key in observed or record["configuration"] != "V2-C" or key not in expected:
            raise ValueError("C release checkpoint membership differs")
        if record["path"] != expected[key]["path"] or record["sha256"] != expected[key]["sha256"]:
            raise ValueError("C release checkpoint hash differs from evaluated checkpoint")
        checked_reference(record)
        observed.add(key)
    if observed != set(expected) or len(observed) != 6:
        raise ValueError("All six C checkpoints must be frozen")
    if set(manifest["models"]) != set(MODELS):
        raise ValueError("Both C contexts required for calibration")
    for model in MODELS:
        if manifest["models"][model]["predictions"] != evaluation["models"][model]["predictions"]:
            raise ValueError("Calibration must use the same frozen prediction file as selection")
    return release


def calibrate(manifest_path, output):
    manifest = read_json(manifest_path)
    validate_calibration_prerequisites(manifest)
    output = Path(output)
    if output.exists() and any(output.iterdir()):
        raise ValueError("Refusing threshold rerun/overwrite of nonempty output directory")
    output.mkdir(parents=True, exist_ok=True)
    thresholds, negative_rows, descriptive_counts = {}, [], []
    for model in MODELS:
        item = manifest["models"][model]
        metadata = list(read_rows(checked_reference(item["calibration"])))
        validate_common_rows(metadata, model, "calibration")
        negative = [row for row in metadata if int(row["label"]) == 0]
        ids = validate_common_rows(negative, model, "calibration", negatives_only=True)
        predictions = load_predictions(checked_reference(item["predictions"]), ids, configurations=("V2-C",))
        scores = predictions["V2-C"]["ensemble"]
        components = [row["component_id"] for row in negative]
        record = conservative_threshold(scores, components)
        thresholds[model] = record
        descriptive_counts.append({"model": model, "calibration_metadata_rows": len(metadata), "positive_metadata_rows_not_scored_or_used": sum(int(row["label"]) == 1 for row in metadata), "negative_rows_used": len(negative), "negative_components_used": len(set(components))})
        for index, row in enumerate(negative):
            negative_rows.append({"model": model, "configuration": "V2-C", "interval_id": row["interval_id"], "chrom": "chr7", "validation_role": "calibration", "label": 0, "component_id": row["component_id"], **{f"q_{seed}": float(predictions["V2-C"]["seeds"][j, index]) for j, seed in enumerate(SEEDS)}, "ensemble_q": float(scores[index]), "region_score_above_threshold": bool(scores[index] >= record["threshold_float64"])})
    write_json(output / "C_region_thresholds.json", {"status": "PASS", "input_manifest_sha256": sha256(manifest_path), "models": thresholds, "no_allele_delta_threshold": True, "calibration_scope": "dedicated common chr7 calibration negative controls only", "claim_limit": "region-label thresholds; not allele effects, causal probabilities, variant FPR or population-FPR guarantees", "next_action": "STOP_FOR_INVESTIGATOR_REVIEW_BEFORE_CHR8_9_OR_EXTERNAL_ACCESS"})
    write_rows(output / "C_region_thresholds.tsv", [{"model": model, **thresholds[model]} for model in MODELS])
    write_rows(output / "calibration_negative_scores.tsv.gz", negative_rows)
    write_rows(output / "calibration_counts.tsv", descriptive_counts)
    write_json(output / "calibration_manifest.json", {"input_manifest": str(manifest_path), "input_manifest_sha256": sha256(manifest_path), "calibration_script_sha256": sha256(__file__), "evaluation_helper_sha256": sha256(Path(__file__).with_name("evaluate_chr7.py")), "artifacts": [{"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)} for path in sorted(output.iterdir())], "positive_prediction_access": False, "test_external_candidate_access": False})
    return thresholds


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-manifest", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    try:
        thresholds = calibrate(args.input_manifest, args.output)
    except Exception as error:
        print(json.dumps({"status": "HARD_STOP", "error": str(error), "no_calibration_released": True}), file=sys.stderr)
        return 2
    print(json.dumps({"status": "PASS", "thresholds": thresholds}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
