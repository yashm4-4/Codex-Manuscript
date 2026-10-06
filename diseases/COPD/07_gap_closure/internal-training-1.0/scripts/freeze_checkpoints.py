#!/usr/bin/env python3
"""Freeze exactly 18 selected checkpoints before any common chr7 inference.

Read-only validation of existing fit outputs, followed by exclusively created
training outcomes and checkpoint release gate. Does not import a model runtime,
load model weights into a network, or examine selection/calibration performance.
"""
from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import sys

import numpy as np

from phase_two_contract import (
    BATCH_SIZE, CONFIGURATIONS, MAX_EPOCHS, MODELS, PATIENCE, SEEDS,
    array_digest, epoch_plan, file_record, read_json, read_rows,
    stage_file, utcnow, verify_record, write_json,
)


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def inspect_history(history_path, selected, contract):
    rows = read_rows(history_path)
    check(0 < len(rows) <= MAX_EPOCHS, "Empty/overlong training history")
    n_train = int(contract["n_train"])
    best, best_epoch, wait = math.inf, None, 0
    for i, row in enumerate(rows):
        check(int(row["epoch_zero_based"]) == i and int(row["epoch_one_based"]) == i + 1,
              "Nonconsecutive training history epoch")
        train_bce, bce = float(row["training_BCE"]), float(row["checkpoint_symmetric_BCE"])
        check(math.isfinite(train_bce) and math.isfinite(bce), "Nonfinite training/checkpoint loss")
        check(int(row["n_train"]) == n_train and int(row["steps"]) == math.ceil(n_train / BATCH_SIZE),
              "Training exposure differs from frozen one-view-per-region rule")
        permutation, orientation = epoch_plan(int(selected["seed"]), i, n_train)
        check(row["permutation_sha256"] == array_digest(permutation)
              and row["orientation_sha256"] == array_digest(orientation)
              and int(row["n_reverse_views"]) == int(orientation.sum()),
              "Epoch shuffle/orientation differs from frozen RNG substreams")
        improved = bce < best
        if improved:
            best, best_epoch, wait = bce, i + 1, 0
        else:
            wait += 1
        check(int(row["selected_new_best"]) == int(improved) and int(row["wait"]) == wait,
              "History checkpoint tie/early-stopping state is inconsistent")
        check(i == len(rows) - 1 or wait < PATIENCE, "Training continued beyond frozen early-stopping patience")
    check(len(rows) == MAX_EPOCHS or wait >= PATIENCE, "Fit stopped before max epochs without reaching patience")
    check(selected["selected_epoch_one_based"] == best_epoch,
          "Selected epoch is not earliest exact minimum native symmetric BCE")
    check(float(selected["checkpoint_symmetric_BCE"]) == best
          and float(selected["restored_checkpoint_symmetric_BCE"]) == best,
          "Selected/restored checkpoint BCE differs from history minimum")
    check(selected["epochs_completed"] == len(rows) and selected["early_stopped"] == (wait >= PATIENCE),
          "Selected checkpoint epoch/stopping metadata mismatch")
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", type=Path, required=True)
    args = parser.parse_args()
    stage = args.stage.resolve()
    gate_path = stage / "provenance/checkpoint_freeze.json"
    outcomes_path = stage / "results/training_outcomes.tsv"
    check(not gate_path.exists() and not outcomes_path.exists(), "Checkpoint freeze/outcomes already exist; refusing overwrite")
    input_gate_path = stage / "provenance/execution_input_gate.json"
    input_gate = read_json(input_gate_path)
    check(input_gate.get("status") == "PASS" and input_gate.get("training_authorized") is True
          and input_gate.get("input_hashes_verified") is True, "Execution input gate is not PASS")
    verified_inputs = {record["path"]: file_record(stage, verify_record(stage, record)) for record in input_gate["files"]}
    script_names = ("phase_two_contract.py", "train_phase_two.py", "predict_chr7.py", "freeze_checkpoints.py", "train_array.sh", "run_logged.py", "runtime.sh")
    for name in script_names:
        check(f"scripts/{name}" in verified_inputs, f"Prospective execution gate did not freeze scripts/{name}")
    check("inputs/input_manifest.json" in verified_inputs, "Prospective execution gate lacks input manifest hash")
    expected = {(configuration, model, seed) for configuration in CONFIGURATIONS for model in MODELS for seed in SEEDS}
    expected_ids = {f"{configuration}_{model}_seed{seed}" for configuration, model, seed in expected}
    runs = stage / "runs"
    actual_dirs = {path.name for path in runs.iterdir() if path.is_dir()}
    check(actual_dirs == expected_ids, f"Unexpected/missing run directories: {actual_dirs ^ expected_ids}")
    all_success_paths = {path.resolve() for path in runs.glob("*/attempt-*/completed.json")}
    check(len(all_success_paths) == 18, "Exactly 18 successful fit attempts are required")
    checkpoints, outcomes, failures = [], [], []
    for configuration in CONFIGURATIONS:
        for model in MODELS:
            for seed in SEEDS:
                run_id = f"{configuration}_{model}_seed{seed}"
                run_root = runs / run_id
                attempts = sorted(path for path in run_root.glob("attempt-*") if path.is_dir())
                check([path.name for path in attempts] == [f"attempt-{i:03d}" for i in range(1, len(attempts) + 1)],
                      f"Attempt numbers not consecutive: {run_id}")
                successful = [path for path in attempts if (path / "completed.json").exists()]
                check(len(successful) == 1 and successful[0] == attempts[-1], f"Nonunique/nonfinal successful fit: {run_id}")
                output = successful[0]
                for i, attempt_path in enumerate(attempts):
                    started = read_json(attempt_path / "started.json")
                    check(started["run_id"] == run_id and started["seed"] == seed and started["attempt"] == i + 1,
                          f"Attempt identity drift: {run_id}")
                    if i:
                        retry = started.get("retry")
                        check(isinstance(retry, dict) and retry.get("authorized") is True
                              and retry.get("reason_class") == "INFRASTRUCTURE"
                              and retry.get("same_seed_identical_procedure") is True
                              and retry.get("run_id") == run_id and retry.get("attempt") == i + 1,
                              f"Unauthorized same-seed rerun: {run_id}")
                        verify_record(stage, retry["authorization_file"])
                    if attempt_path != output:
                        failure_path = attempt_path / "failure.json"
                        failures.append({"run_id": run_id, "attempt": i + 1,
                                         "failure_record": file_record(stage, failure_path) if failure_path.exists() else None,
                                         "retry_authorization": read_json(attempts[i + 1] / "started.json")["retry"]})
                selected_path = output / "selected_checkpoint.json"
                completed_path = output / "completed.json"
                contract_path = output / "training_contract.json"
                selected, completed, contract = read_json(selected_path), read_json(completed_path), read_json(contract_path)
                for record in (selected, completed):
                    check(record["configuration"] == configuration and record["model"] == model
                          and record["seed"] == seed and record.get("completed") is True
                          and record["run_id"] == run_id, f"Selected/completed identity mismatch: {run_id}")
                    check(all(record.get(name) is False for name in (
                        "selection_role_access", "calibration_role_access", "chr8_chr9_access", "external_benchmark_access")),
                        f"Checkpoint data firewall violation: {run_id}")
                check(completed.get("status") == "COMPLETED", f"Invalid completed status: {run_id}")
                check(selected.get("serialized_selected_weights_bitwise_roundtrip") is True
                      and selected.get("checkpoint_is_inference_only") is True, f"Selected weight roundtrip not PASS: {run_id}")
                check(contract["input_gate"] == file_record(stage, input_gate_path), f"Input gate changed after training: {run_id}")
                for record in contract["inputs"]:
                    verify_record(stage, record)
                verify_record(stage, contract["cache_manifest"])
                check(contract["batch_size"] == BATCH_SIZE and contract["max_epochs"] == MAX_EPOCHS
                      and contract["patience"] == PATIENCE and contract["min_delta"] == 0
                      and contract["unweighted_BCE"] is True and contract["mixed_precision"] is False
                      and contract["jit_compile"] is False, f"Frozen training parameters differ: {run_id}")
                optimizer = contract["optimizer"]
                check(str(optimizer["name"]).lower() == "adadelta"
                      and abs(float(optimizer["learning_rate"]) - .001) < 1e-9
                      and float(optimizer["rho"]) == .95 and float(optimizer["epsilon"]) == 1e-7,
                      f"Frozen Adadelta parameters differ: {run_id}")
                checkpoint_path = verify_record(stage, selected["checkpoint"])
                history_path = verify_record(stage, selected["history"])
                check(completed["checkpoint"] == selected["checkpoint"] and completed["history"] == selected["history"],
                      f"Completion and selected checkpoint provenance disagree: {run_id}")
                history = inspect_history(history_path, selected, contract)
                improvements_path = output / "checkpoint_improvements.tsv"
                improvements = read_rows(improvements_path)
                expected_improvements = [row for row in history if int(row["selected_new_best"]) == 1]
                check(len(improvements) == len(expected_improvements) and all(
                    int(actual["epoch_one_based"]) == int(expected_row["epoch_one_based"])
                    and float(actual["checkpoint_symmetric_BCE"]) == float(expected_row["checkpoint_symmetric_BCE"])
                    for actual, expected_row in zip(improvements, expected_improvements)),
                    f"Checkpoint improvement ledger differs from history: {run_id}")
                command_label = f"{run_id}_attempt{selected['attempt']:03d}"
                command_path = stage / "provenance/commands" / f"{command_label}.completed.json"
                command = read_json(command_path)
                check(command["returncode"] == 0 and command["label"] == command_label
                      and command["environment"]["PYTHONHASHSEED"] == str(seed), f"Missing/suspect exact-command completion: {run_id}")
                flat = {"configuration": configuration, "model": model, "seed": seed, "run_id": run_id,
                        "completed": True, **file_record(stage, checkpoint_path),
                        "selected_checkpoint_record": file_record(stage, selected_path),
                        "completed_record": file_record(stage, completed_path),
                        "history": file_record(stage, history_path),
                        "checkpoint_improvements": file_record(stage, improvements_path),
                        "training_contract": file_record(stage, contract_path),
                        "command_record": file_record(stage, command_path),
                        "selected_epoch_one_based": selected["selected_epoch_one_based"],
                        "checkpoint_symmetric_BCE": selected["checkpoint_symmetric_BCE"]}
                checkpoints.append(flat)
                outcomes.append({"run_id": run_id, "configuration": configuration, "model": model, "seed": seed,
                                 "status": "COMPLETED", "successful_attempt": selected["attempt"],
                                 "infrastructure_retries": len(attempts) - 1, "epochs_completed": len(history),
                                 "selected_epoch_one_based": selected["selected_epoch_one_based"],
                                 "checkpoint_symmetric_BCE": format(float(selected["checkpoint_symmetric_BCE"]), ".17g"),
                                 "early_stopped": selected["early_stopped"],
                                 "checkpoint_sha256": flat["sha256"], "wall_seconds": completed["wall_seconds"]})
    check(len(checkpoints) == 18, "Checkpoint count mismatch")
    outcomes_path.parent.mkdir(parents=True, exist_ok=True)
    with outcomes_path.open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(outcomes[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(outcomes)
    record = {"status": "PASS", "created_utc": utcnow(), "argv": sys.argv,
              "all_18_selected_checkpoints_frozen": True, "selection_prediction_authorized": True,
              "checkpoint_scope": "selected architecture and all selected-epoch trainable/BatchNorm weights; optimizer resumption prohibited",
              "checkpoints": checkpoints, "execution_input_gate": file_record(stage, input_gate_path),
              "input_manifest": verified_inputs["inputs/input_manifest.json"],
              "training_outcomes": file_record(stage, outcomes_path),
              "scripts": [verified_inputs[f"scripts/{name}"] for name in script_names],
              "infrastructure_failures_and_retries": failures,
              "selection_predictions_previously_examined": False, "calibration_previously_performed": False,
              "real_network_invariance": "NEXT_REQUIRED_GATE_NOT_YET_RUN",
              "chr8_chr9_access": False, "external_benchmark_access": False,
              "freeze_allows_only": "chr7 actual-network RC gate and authorized internal common-panel evaluation; no test/external/candidate scoring"}
    write_json(gate_path, record)
    print(f"PASS: froze {len(checkpoints)} selected checkpoints; {len(failures)} recorded infrastructure retries")


if __name__ == "__main__":
    main()
