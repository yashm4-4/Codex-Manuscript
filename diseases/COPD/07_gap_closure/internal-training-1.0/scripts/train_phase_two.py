#!/usr/bin/env python3
"""One authorized, exact frozen phase-II fit; no selection/calibration reads.

Run with PYTHONHASHSEED matching --seed. Only a same-seed retry explicitly
classified as infrastructure failure may create attempt 002 or later. This
entry point never examines common-panel outcome metrics or changes the design.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys
import time
import traceback

import numpy as np

from phase_two_contract import (
    BATCH_SIZE, CONFIGURATIONS, MAX_EPOCHS, MODELS, PATIENCE, SEEDS,
    FeatureCache, array_digest, build_phase_two, epoch_plan, file_record,
    inference_function, initialize_runtime, load_input_table, predict_orientation,
    read_json, resource_record, runtime_environment, sha256, stage_file,
    symmetric_bce, symmetric_probability, utcnow, verify_record, write_json,
)


def input_gate(stage, configuration, model):
    path = stage / "provenance/execution_input_gate.json"
    gate = read_json(path)
    if gate.get("status") != "PASS" or gate.get("training_authorized") is not True or gate.get("input_hashes_verified") is not True:
        raise RuntimeError("Training execution integrity/authorization gate is not PASS")
    records = {record["path"]: record for record in gate["files"]}
    required = [f"inputs/configurations/{configuration}_{model}_train.tsv.gz",
                f"inputs/configurations/{configuration}_{model}_checkpoint.tsv.gz",
                "cache/cache_manifest.json", "scripts/train_phase_two.py", "scripts/phase_two_contract.py"]
    checked = []
    for name in required:
        checked.append(file_record(stage, verify_record(stage, records[name])))
    return path, checked


def prepare_attempt(stage, run_id, attempt, retry_authorization):
    run_root = stage / "runs" / run_id
    prior = sorted(run_root.glob("attempt-*")) if run_root.exists() else []
    if any((path / "completed.json").exists() for path in prior):
        raise RuntimeError("This prescribed fit already completed; duplicate training is forbidden")
    if attempt != len(prior) + 1:
        raise RuntimeError("Attempt number must be consecutive and new")
    retry_record = None
    if attempt > 1:
        if retry_authorization is None:
            raise RuntimeError("Explicit infrastructure retry authorization required")
        retry_path = stage_file(stage, retry_authorization)
        retry_record = read_json(retry_path)
        if (retry_record.get("run_id") != run_id or retry_record.get("attempt") != attempt
                or retry_record.get("reason_class") != "INFRASTRUCTURE"
                or retry_record.get("same_seed_identical_procedure") is not True
                or retry_record.get("authorized") is not True
                or not str(retry_record.get("reason", "")).strip()):
            raise RuntimeError("Invalid infrastructure retry authorization")
        retry_record["authorization_file"] = file_record(stage, retry_path)
    elif retry_authorization is not None:
        raise RuntimeError("Attempt 1 must not masquerade as a retry")
    output = run_root / f"attempt-{attempt:03d}"
    output.mkdir(parents=True, exist_ok=False)
    return output, retry_record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", type=Path, required=True)
    parser.add_argument("--configuration", choices=CONFIGURATIONS, required=True)
    parser.add_argument("--model", choices=MODELS, required=True)
    parser.add_argument("--seed", type=int, choices=SEEDS, required=True)
    parser.add_argument("--attempt", type=int, default=1)
    parser.add_argument("--retry-authorization", help="Stage-relative infrastructure-only retry JSON")
    args = parser.parse_args()
    stage = args.stage.resolve()
    run_id = f"{args.configuration}_{args.model}_seed{args.seed}"
    output, retry_record = prepare_attempt(stage, run_id, args.attempt, args.retry_authorization)
    started = time.monotonic()
    start_record = {"run_id": run_id, "configuration": args.configuration, "model": args.model,
                    "seed": args.seed, "attempt": args.attempt, "started_utc": utcnow(),
                    "argv": sys.argv, "retry": retry_record,
                    "status": "STARTED", "selection_role_access": False,
                    "calibration_role_access": False, "test_access": False,
                    "external_benchmark_access": False}
    write_json(output / "started.json", start_record)
    try:
        gate_path, input_records = input_gate(stage, args.configuration, args.model)
        train = load_input_table(stage / f"inputs/configurations/{args.configuration}_{args.model}_train.tsv.gz", "train")
        checkpoint = load_input_table(stage / f"inputs/configurations/{args.configuration}_{args.model}_checkpoint.tsv.gz", "checkpoint")
        if set(train["interval_id"]) & set(checkpoint["interval_id"]):
            raise RuntimeError("Training/checkpoint interval overlap")
        cache = FeatureCache(stage)
        cache.check_indices(train)
        cache.check_indices(checkpoint)
        tf, keras = initialize_runtime(args.seed)
        write_json(output / "environment.json", runtime_environment(tf, keras))
        model = build_phase_two(keras)
        if int(model.input_shape[1]) != 4560 or model.output_shape != (None, 1):
            raise RuntimeError("Unexpected phase-II input/output geometry")
        model_architecture = json.loads(model.to_json())
        write_json(output / "architecture.json", model_architecture)
        write_json(output / "training_contract.json", {
            "run_id": run_id, "n_train": len(train["label"]), "n_checkpoint": len(checkpoint["label"]),
            "input_gate": file_record(stage, gate_path), "inputs": input_records,
            "cache_manifest": file_record(stage, cache.manifest_path),
            "batch_size": BATCH_SIZE, "max_epochs": MAX_EPOCHS, "patience": PATIENCE,
            "min_delta": 0, "unweighted_BCE": True, "optimizer": model.optimizer.get_config(),
            "trainable_parameters": int(sum(np.prod(weight.shape) for weight in model.trainable_weights)),
            "parameters_total": int(model.count_params()),
            "initial_weights_digest": array_digest(np.concatenate([weight.reshape(-1) for weight in model.get_weights()])),
            "epoch_index": "zero-based", "shuffle_substream": "SeedSequence([seed,1,epoch])",
            "orientation_substream": "SeedSequence([seed,2,epoch])",
            "orientation_draw": "random(n)<0.5 in frozen genomic row order before applying permutation",
            "checkpoint_BCE": "float64 symmetric probabilities; clip to [1e-7,1-1e-7]; float64 mean unweighted Bernoulli log loss",
            "checkpoint_ties": "earliest epoch, update only on strict numerical decrease",
            "selected_checkpoint_scope": "inference-only full architecture and all selected-epoch trainable/BatchNorm weights; no resumable optimizer state",
            "mixed_precision": False, "jit_compile": False,
        })
        infer = inference_function(tf, model)
        best_bce, best_epoch, best_weights, wait = math.inf, None, None, 0
        history_path = output / "history.tsv"
        improvements_path = output / "checkpoint_improvements.tsv"
        history_fields = ["epoch_zero_based", "epoch_one_based", "training_BCE", "checkpoint_symmetric_BCE",
                          "selected_new_best", "wait", "n_train", "steps", "n_reverse_views",
                          "permutation_sha256", "orientation_sha256", "wall_seconds"]
        with history_path.open("x", newline="") as history_stream, improvements_path.open("x", newline="") as improvement_stream:
            writer = csv.DictWriter(history_stream, fieldnames=history_fields, delimiter="\t")
            writer.writeheader()
            improvement_writer = csv.DictWriter(improvement_stream, fieldnames=["epoch_one_based", "checkpoint_symmetric_BCE"], delimiter="\t")
            improvement_writer.writeheader()
            for epoch in range(MAX_EPOCHS):
                epoch_start = time.monotonic()
                order, reverse = epoch_plan(args.seed, epoch, len(train["label"]))
                model.reset_metrics()
                training_bce = None
                steps = 0
                for offset in range(0, len(order), BATCH_SIZE):
                    indices = order[offset:offset + BATCH_SIZE]
                    metrics = model.train_on_batch(cache.batch(train, indices, reverse[indices]),
                                                   train["label"][indices, None], return_dict=True)
                    training_bce = float(metrics["loss"])
                    if not math.isfinite(training_bce):
                        raise RuntimeError("Nonfinite training loss is a hard scientific QC stop")
                    steps += 1
                pf = predict_orientation(infer, cache, checkpoint)
                pr = predict_orientation(infer, cache, checkpoint, reverse=True)
                checkpoint_bce = symmetric_bce(checkpoint["label"], symmetric_probability(pf, pr))
                improved = checkpoint_bce < best_bce
                if improved:
                    best_bce, best_epoch = checkpoint_bce, epoch + 1
                    best_weights = model.get_weights()
                    wait = 0
                    improvement_writer.writerow({"epoch_one_based": best_epoch, "checkpoint_symmetric_BCE": format(best_bce, ".17g")})
                    improvement_stream.flush()
                else:
                    wait += 1
                writer.writerow({"epoch_zero_based": epoch, "epoch_one_based": epoch + 1,
                                 "training_BCE": format(training_bce, ".17g"),
                                 "checkpoint_symmetric_BCE": format(checkpoint_bce, ".17g"),
                                 "selected_new_best": int(improved), "wait": wait,
                                 "n_train": len(order), "steps": steps, "n_reverse_views": int(reverse.sum()),
                                 "permutation_sha256": array_digest(order), "orientation_sha256": array_digest(reverse),
                                 "wall_seconds": time.monotonic() - epoch_start})
                history_stream.flush()
                print(json.dumps({"run_id": run_id, "epoch": epoch + 1, "checkpoint_BCE": checkpoint_bce,
                                  "best_epoch": best_epoch, "best_BCE": best_bce, "wait": wait}), flush=True)
                if wait >= PATIENCE:
                    break
        if best_weights is None:
            raise RuntimeError("No finite selected checkpoint")
        model.set_weights(best_weights)
        selected_pf = predict_orientation(infer, cache, checkpoint)
        selected_pr = predict_orientation(infer, cache, checkpoint, reverse=True)
        restored_bce = symmetric_bce(checkpoint["label"], symmetric_probability(selected_pf, selected_pr))
        if restored_bce != best_bce:
            raise RuntimeError(f"Restored selected weights changed native checkpoint BCE: {restored_bce} vs {best_bce}")
        # JSON contains architecture only, so this is uncompiled and cannot
        # accidentally serialize final-epoch or freshly reset optimizer slots.
        selected = keras.models.model_from_json(model.to_json())
        selected.set_weights(best_weights)
        checkpoint_path = output / "selected_checkpoint.keras"
        if checkpoint_path.exists():
            raise RuntimeError("Refusing checkpoint overwrite")
        selected.save(checkpoint_path)
        reloaded = keras.models.load_model(checkpoint_path, compile=False)
        reloaded_weights = reloaded.get_weights()
        if len(reloaded_weights) != len(best_weights) or any(
                not np.array_equal(left, right) for left, right in zip(reloaded_weights, best_weights)):
            raise RuntimeError("Serialized selected checkpoint weight round-trip failed")
        record = {"run_id": run_id, "configuration": args.configuration, "model": args.model,
                  "seed": args.seed, "attempt": args.attempt, "completed": True,
                  "selected_epoch_one_based": best_epoch, "epochs_completed": epoch + 1,
                  "early_stopped": wait >= PATIENCE, "checkpoint_symmetric_BCE": best_bce,
                  "restored_checkpoint_symmetric_BCE": restored_bce,
                  "checkpoint": file_record(stage, checkpoint_path),
                  "history": file_record(stage, history_path),
                  "native_checkpoint_orientation_hashes": {"forward": array_digest(selected_pf), "rc": array_digest(selected_pr)},
                  "selection_role_access": False, "calibration_role_access": False,
                  "chr8_chr9_access": False, "external_benchmark_access": False,
                  "checkpoint_is_inference_only": True}
        record["serialized_selected_weights_bitwise_roundtrip"] = True
        write_json(output / "selected_checkpoint.json", record)
        write_json(output / "completed.json", {**record, "status": "COMPLETED", **resource_record(started)})
        print(json.dumps({"run_id": run_id, "status": "COMPLETED", "selected_epoch": best_epoch,
                          "checkpoint_sha256": record["checkpoint"]["sha256"]}), flush=True)
    except BaseException as exc:
        write_json(output / "failure.json", {**start_record, "status": "FAILED", "exception_type": type(exc).__name__,
                    "error": str(exc), "traceback": traceback.format_exc(), **resource_record(started),
                    "retry_automatically_authorized": False})
        raise


if __name__ == "__main__":
    main()
