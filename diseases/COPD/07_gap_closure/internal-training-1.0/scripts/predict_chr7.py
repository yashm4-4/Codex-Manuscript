#!/usr/bin/env python3
"""Post-checkpoint chr7 inference plus independent actual-network RC audit.

All 18 selected checkpoints must be frozen before this entry point can open
the union of frozen internal chr7 sequences. It never reads chr8/9 or external
benchmark data, never computes selection metrics, and never chooses thresholds.
The reverse-wrapper pass runs the actual saved network again with swapped
nucleotide-orientation cache lookup and the same fixed order/batch size 256.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import io
from pathlib import Path
import sys
import time
import traceback

import numpy as np

from phase_two_contract import (
    ATOL, CONFIGURATIONS, MODELS, RTOL, SEEDS, FeatureCache, array_digest,
    file_record, inference_function, initialize_runtime, load_input_table,
    predict_orientation, read_json, resource_record, runtime_environment,
    seed_ensemble, symmetric_probability, utcnow, verify_record, write_json,
)


def write_table(path, fields, rows):
    with Path(path).open("xb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            with io.TextIOWrapper(compressed, encoding="utf-8", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)


def numeric_gate(original, reversed_input):
    a = np.asarray(original, dtype=np.float64)
    b = np.asarray(reversed_input, dtype=np.float64)
    residual = np.abs(a - b)
    bound = ATOL + RTOL * np.abs(b)
    failed = residual > bound
    return {"status": "PASS" if not failed.any() else "FAIL", "n_sequences": len(a),
            "n_failed": int(failed.sum()), "atol": ATOL, "rtol": RTOL,
            "max_absolute_difference": float(residual.max(initial=0.0)),
            "max_tolerance_normalized_difference": float(np.max(residual / bound, initial=0.0)),
            "original_wrapper_array_sha256": array_digest(a),
            "reverse_wrapper_array_sha256": array_digest(b)}


def verify_prerequisites(stage):
    input_gate_path = stage / "provenance/execution_input_gate.json"
    inputs = read_json(input_gate_path)
    if inputs.get("status") != "PASS" or inputs.get("input_hashes_verified") is not True:
        raise RuntimeError("Execution input integrity gate is not PASS")
    records = {record["path"]: record for record in inputs["files"]}
    for path in ["cache/cache_manifest.json", "scripts/predict_chr7.py", "scripts/phase_two_contract.py",
                 *[f"inputs/validation/{model}_chr7.tsv.gz" for model in MODELS]]:
        verify_record(stage, records[path])
    gate_path = stage / "provenance/checkpoint_freeze.json"
    gate = read_json(gate_path)
    if (gate.get("status") != "PASS" or gate.get("all_18_selected_checkpoints_frozen") is not True
            or gate.get("selection_prediction_authorized") is not True):
        raise RuntimeError("All-checkpoint freeze and post-checkpoint chr7 inference gate not PASS")
    expected = {(configuration, model, seed) for configuration in CONFIGURATIONS for model in MODELS for seed in SEEDS}
    checkpoints = {}
    for record in gate["checkpoints"]:
        key = record["configuration"], record["model"], int(record["seed"])
        if key not in expected or key in checkpoints or record.get("completed") is not True:
            raise RuntimeError("Checkpoint freeze has unexpected, duplicated, or incomplete fit")
        checkpoints[key] = verify_record(stage, record)
    if set(checkpoints) != expected:
        raise RuntimeError("Not all 18 prescribed selected checkpoints are frozen")
    return checkpoints, gate_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", type=Path, required=True)
    parser.add_argument("--output-relative", default="predictions",
                        help="New stage-relative output directory; alternative only for a separately logged infrastructure retry")
    args = parser.parse_args()
    stage = args.stage.resolve()
    output = (stage / args.output_relative).resolve()
    if not output.is_relative_to(stage):
        raise RuntimeError("Output path must be within new training stage")
    output.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    write_json(output / "started.json", {"started_utc": utcnow(), "argv": sys.argv,
               "scope": "frozen internal chr7 union only; no metric evaluation or threshold selection"})
    try:
        checkpoints, checkpoint_gate = verify_prerequisites(stage)
        cache = FeatureCache(stage)
        tf, keras = initialize_runtime(SEEDS[0])
        write_json(output / "environment.json", runtime_environment(tf, keras))
        audit_rows, output_files = [], []
        for context in MODELS:
            table = load_input_table(stage / f"inputs/validation/{context}_chr7.tsv.gz", "validation")
            cache.check_indices(table)
            n = len(table["interval_id"])
            predictions, repeated = {}, {}
            for configuration in CONFIGURATIONS:
                ordered_original, ordered_reversed = [], []
                for seed in SEEDS:
                    keras.backend.clear_session()
                    model = keras.models.load_model(checkpoints[configuration, context, seed], compile=False)
                    if model.input_shape != (None, 4560, 1) or model.output_shape != (None, 1):
                        raise RuntimeError("Loaded checkpoint has unexpected geometry")
                    infer = inference_function(tf, model)
                    # First complete wrapper evaluation: Q(x), original genomic order.
                    pf = predict_orientation(infer, cache, table, reverse=False, batch_size=256)
                    pr = predict_orientation(infer, cache, table, reverse=True, batch_size=256)
                    q = symmetric_probability(pf, pr)
                    # Independent REAL model calls for Q(RC(x)); do not reuse pf/pr.
                    swapped_first = predict_orientation(infer, cache, table, reverse=True, batch_size=256)
                    swapped_second = predict_orientation(infer, cache, table, reverse=False, batch_size=256)
                    q_rc_input = symmetric_probability(swapped_first, swapped_second)
                    row = {"configuration": configuration, "model": context, "unit": "seed", "seed": seed,
                           **numeric_gate(q, q_rc_input), "real_network": True,
                           "original_batch_size": 256, "reverse_batch_size": 256,
                           "reverse_pass_row_order": "same frozen genomic row order and fixed batching as original wrapper",
                           "reverse_pass_forward_orientation": "nucleotide RC via independently extracted opposite phase-I cache",
                           "original_forward_array_sha256": array_digest(pf),
                           "original_rc_array_sha256": array_digest(pr),
                           "reverse_pass_first_array_sha256": array_digest(swapped_first),
                           "reverse_pass_second_array_sha256": array_digest(swapped_second),
                           "max_same_orientation_repeated_forward_difference": float(np.max(np.abs(pf - swapped_second))),
                           "max_same_orientation_repeated_rc_difference": float(np.max(np.abs(pr - swapped_first))),
                           "checkpoint": file_record(stage, checkpoints[configuration, context, seed])}
                    audit_rows.append(row)
                    write_json(output / f"{configuration}_{context}_seed{seed}_invariance.json", row)
                    if row["status"] != "PASS":
                        raise RuntimeError(f"Actual-network RC invariance HARD STOP: {configuration}/{context}/{seed}")
                    prefix = f"{configuration}"
                    predictions[f"{prefix}_p_forward_{seed}"] = pf
                    predictions[f"{prefix}_p_rc_{seed}"] = pr
                    predictions[f"{prefix}_q_{seed}"] = q
                    repeated[f"{prefix}_rc_input_first_{seed}"] = swapped_first
                    repeated[f"{prefix}_rc_input_second_{seed}"] = swapped_second
                    repeated[f"{prefix}_q_rc_input_{seed}"] = q_rc_input
                    ordered_original.append(q)
                    ordered_reversed.append(q_rc_input)
                    print(f"RC audit PASS: {configuration} {context} seed {seed}; {n} chr7 intervals", flush=True)
                    del infer, model
                ensemble = seed_ensemble(ordered_original)
                ensemble_reversed = seed_ensemble(ordered_reversed)
                row = {"configuration": configuration, "model": context, "unit": "ensemble", "seed": "all3",
                       **numeric_gate(ensemble, ensemble_reversed), "real_network": True,
                       "fixed_seed_order": list(SEEDS), "reduction_dtype": "float64"}
                audit_rows.append(row)
                write_json(output / f"{configuration}_{context}_ensemble_invariance.json", row)
                if row["status"] != "PASS":
                    raise RuntimeError(f"Actual-network ensemble RC invariance HARD STOP: {configuration}/{context}")
                predictions[f"{configuration}_ensemble_q"] = ensemble
                repeated[f"{configuration}_ensemble_q_rc_input"] = ensemble_reversed
            def rows_for(columns):
                for i in range(n):
                    yield {"interval_id": str(table["interval_id"][i]), "chrom": "chr7",
                           **{name: format(float(values[i]), ".17g") for name, values in columns.items()}}
            path = output / f"{context}_chr7.tsv.gz"
            write_table(path, ["interval_id", "chrom", *predictions], rows_for(predictions))
            output_files.append(file_record(stage, path))
            path = output / f"{context}_independent_rc_pass.tsv.gz"
            write_table(path, ["interval_id", "chrom", *repeated], rows_for(repeated))
            output_files.append(file_record(stage, path))
        record = {"status": "PASS", "real_network": True,
                  "all_configurations_all_seeds_and_ensembles": True,
                  "n_seed_audits": 18, "n_ensemble_audits": 6, "audit_rows": audit_rows,
                  "finite_range_gate": "PASS on every raw orientation, repeated orientation, seed symmetric and ensemble probability",
                  "validation_scope": "union of ALL frozen configuration-native chr7 intervals and both common chr7 roles per context",
                  "frozen_input_tables": [file_record(stage, stage / f"inputs/validation/{context}_chr7.tsv.gz") for context in MODELS],
                  "checkpoint_freeze": file_record(stage, checkpoint_gate),
                  "cache_manifest": file_record(stage, cache.manifest_path),
                  "script_files": [file_record(stage, Path(__file__)), file_record(stage, Path(__file__).with_name("phase_two_contract.py"))],
                  "files": output_files, "selection_metrics_computed": False,
                  "calibration_performed": False, "chr8_chr9_access": False,
                  "external_benchmark_access": False, **resource_record(started)}
        write_json(output / "real_network_invariance.json", record)
        write_json(output / "completed.json", {"status": "COMPLETED", "invariance_status": "PASS",
                   "invariance_record": file_record(stage, output / "real_network_invariance.json"),
                   "files": output_files, **resource_record(started)})
    except BaseException as exc:
        write_json(output / "failure.json", {"status": "FAILED", "error": str(exc),
                   "exception_type": type(exc).__name__, "traceback": traceback.format_exc(),
                   "downstream_evaluation_authorized": False, **resource_record(started)})
        raise


if __name__ == "__main__":
    main()
