#!/usr/bin/env python3
"""Independent, non-model-executing audit of COPD V2 internal training.

No TensorFlow imports, model loading, scheduler calls, candidate scoring, or
benchmark/test outcome parsing are permitted in this validator. Input hashing
is deliberately opaque. Scientific outcomes never change the frozen contract.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from datetime import datetime, timezone
from pathlib import Path
import sys
import subprocess

ROOT = Path(__file__).resolve().parents[5]
STAGE = ROOT / "diseases/COPD/07_gap_closure/internal-training-1.0"
BASE = ROOT / "diseases/COPD/07_gap_closure"
PREFLIGHT = BASE / "pretraining-1.1"
SEEDS = (104729, 130363, 155921)
MODELS = ("enhancer", "h3k27me3")
CONFIGS = ("V2-A", "V2-B", "V2-C")
WEIGHT_SHA = "483b6c0cafd750ea8e97932cbc8a949276eaf1e77e7d8766f9a22eeddc2916d6"
FREEZE_SHA = {
    "pretraining-1.0": "2d022d19b4fbdd485d0762c0946519334e654b1de0581c5075758b7550d15fcb",
    "pretraining-1.1": "6933f98bd6a1a7b7afe5adddc8a6bf9d7b1201df4b25f06d0cd8c35e183941a0",
}


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class Audit:
    def __init__(self):
        self.rows = []

    def check(self, name, passed, observed="", expected=""):
        self.rows.append({"check": name, "status": "PASS" if bool(passed) else "FAIL",
                          "observed": str(observed), "expected": str(expected)})

    def opaque_hash(self, name, path, expected, nbytes=None):
        path = Path(path)
        try:
            actual = sha256(path)
            self.check(name + ":sha256", actual == expected, actual, expected)
            if nbytes is not None:
                self.check(name + ":bytes", path.stat().st_size == int(nbytes),
                           path.stat().st_size, nbytes)
        except (FileNotFoundError, PermissionError) as exc:
            self.check(name + ":accessible", False, str(exc), "readable frozen file")

    def write(self, prefix, extra=None):
        prefix = Path(prefix).resolve()
        paths = [prefix.with_suffix(".tsv"), prefix.with_suffix(".json")]
        if any(p.exists() for p in paths):
            raise FileExistsError("Audit outputs are append-only; choose a new --output-prefix")
        prefix.parent.mkdir(parents=True, exist_ok=True)
        with paths[0].open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["check", "status", "observed", "expected"], delimiter="\t")
            writer.writeheader()
            writer.writerows(self.rows)
        failures = sum(row["status"] == "FAIL" for row in self.rows)
        result = {"module": "COPD-V2-INTERNAL-TRAINING-1.0", "validator": "independent-execution-audit",
                  "created_utc": datetime.now(timezone.utc).isoformat(), "checks": len(self.rows),
                  "failures": failures, "status": "PASS" if not failures else "FAIL",
                  "model_execution": False, "test_or_benchmark_outcomes_parsed": False,
                  "checks_path": str(paths[0].relative_to(ROOT)), "checks_sha256": sha256(paths[0]),
                  "script_sha256": sha256(__file__), **(extra or {})}
        paths[1].write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
        print(json.dumps(result, indent=2, sort_keys=True))
        return failures


def validate_inputs(audit):
    freezes = {
        "pretraining-1.0": BASE / "provenance/COPD-V2-PREFLIGHT_freeze.json",
        "pretraining-1.1": PREFLIGHT / "provenance/freeze.json",
    }
    expected_counts = {"pretraining-1.0": 100, "pretraining-1.1": 124}
    for version, path in freezes.items():
        audit.opaque_hash(version + ":freeze", path, FREEZE_SHA[version])
        freeze = json.loads(path.read_text())
        ledger_key = "artifact_checksum_ledger"
        if ledger_key in freeze:
            ledger = ROOT / freeze[ledger_key]
            ledger_sha = freeze[ledger_key + "_sha256"]
        else:
            ledger = BASE / "provenance/COPD-V2-PREFLIGHT_artifact_checksums.tsv"
            ledger_sha = "7a6440c17be3374e5b4732c0bd6a5450fa61b74e4849e99a19b856a02e9e5055"
        audit.opaque_hash(version + ":ledger", ledger, ledger_sha)
        entries = read_tsv(ledger)
        audit.check(version + ":payload_count", len(entries) == expected_counts[version], len(entries), expected_counts[version])
        audit.check(version + ":unique_payload_paths", len({row["path"] for row in entries}) == len(entries))
        for row in entries:
            audit.opaque_hash(version + ":" + row["path"], ROOT / row["path"], row["sha256"], row["bytes"])
    inputs = read_tsv(BASE / "provenance/COPD-V2-PREFLIGHT_consumed_input_hashes.tsv")
    audit.check("consumed_input_count", len(inputs) == 28, len(inputs), 28)
    for row in inputs:
        audit.opaque_hash("consumed_input:" + row["path"], ROOT / row["path"], row["sha256"], row["bytes"])
    audit.opaque_hash("phase_I_weights", ROOT / "diseases/COPD/04_modeling/trednet/model_phase_I/phase_one_weights.h5", WEIGHT_SHA)
    specification = json.loads((PREFLIGHT / "specification/selection_adequacy.json").read_text())
    audit.check("exact_three_seeds", specification["seeds"] == list(SEEDS), specification["seeds"], SEEDS)
    audit.check("eighteen_fits", specification["n_planned_fits"] == 18)
    audit.check("phaseI_frozen_weight_contract", specification["phase_I_weights_sha256"] == WEIGHT_SHA)
    audit.check("fixed_common_C_task", specification["configurations"]["V2-C"]["eligible_final_configuration"] and not specification["configurations"]["V2-A"]["eligible_final_configuration"] and not specification["configurations"]["V2-B"]["eligible_final_configuration"])


def tracked_baseline(audit, prefix):
    """Opaque preservation ledger; benchmark contents are never interpreted."""
    prefix = Path(prefix).resolve()
    ledger = prefix.with_name(prefix.name + "_files.tsv")
    rows = []
    paths = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT).decode().split("\0")
    for name in filter(None, paths):
        path = ROOT / name
        if path.is_symlink():
            data = os.fsencode(os.readlink(path))
            row = {"path": name, "kind": "symlink_target_text", "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
        else:
            row = {"path": name, "kind": "file_bytes", "bytes": path.stat().st_size, "sha256": sha256(path)}
        rows.append(row)
        audit.check("baseline:opaque_hash:" + name, True, row["sha256"])
    ledger.parent.mkdir(parents=True, exist_ok=True)
    with ledger.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["path", "kind", "bytes", "sha256"], delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return {"tracked_ledger": str(ledger.relative_to(ROOT)), "tracked_ledger_sha256": sha256(ledger),
            "tracked_count": len(rows), "baseline_commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT).decode().strip(),
            "benchmark_access": "opaque integrity hashing only"}


def validate_preservation(audit, baseline_manifest):
    baseline = json.loads(Path(baseline_manifest).read_text())
    ledger = ROOT / baseline["tracked_ledger"]
    audit.opaque_hash("baseline_ledger", ledger, baseline["tracked_ledger_sha256"])
    allowed = {"diseases/COPD/07_gap_closure/" + name for name in ("activity_log.tsv", "gap_closure_decision_register.tsv", "gap_closure_result_register.tsv")}
    for row in read_tsv(ledger):
        path = ROOT / row["path"]
        if row["kind"] == "symlink_target_text":
            actual = hashlib.sha256(os.fsencode(os.readlink(path))).hexdigest() if path.is_symlink() else "NOT_SYMLINK"
        else:
            actual = sha256(path)
        if row["path"] not in allowed:
            audit.check("preserved:" + row["path"], actual == row["sha256"], actual, row["sha256"])
        else:
            original = subprocess.check_output(["git", "show", baseline["baseline_commit"] + ":" + row["path"]], cwd=ROOT)
            current = path.read_bytes()
            audit.check("register_append_only:" + row["path"], current.startswith(original))
    audit.check("baseline_commit_preserved", subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT).decode().strip() == baseline["baseline_commit"])


def validate_prepared(audit):
    import numpy as np
    import pandas as pd
    manifest = json.loads((STAGE / "inputs/input_manifest.json").read_text())
    audit.check("prepared:status", manifest["status"] == "PASS" and manifest["completed"] is True)
    for row in manifest["files"]:
        audit.opaque_hash("prepared_artifact:" + row["path"], STAGE / row["path"], row["sha256"], row["bytes"])
    columns = ["interval_id", "chrom", "core_start", "core_end", "input_start", "input_end", "partition", "validation_role", "component_id", "sequence_available", "sequence_length", "sequence_sha256", "canonical_rc_sequence_sha256", "gc_fraction", "non_acgt_fraction"]
    def frame(path, usecols=None):
        return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, usecols=usecols)
    index = frame(STAGE / "inputs/cache_index.tsv.gz")
    allowed = {f"chr{i}" for i in range(1,23) if i not in (8,9)} | {"chrX", "chrY"}
    audit.check("prepared:sequence_scope", set(index.chrom) <= allowed and not index.chrom.isin(["chr8", "chr9"]).any())
    audit.check("prepared:unique_IDs", not index.interval_id.duplicated().any())
    audit.check("prepared:exact_geometry", ((index.core_end.astype(int)-index.core_start.astype(int) == 1000) & (index.input_start.astype(int) == index.core_start.astype(int)-501) & (index.input_end.astype(int) == index.core_end.astype(int)+500) & index.sequence_length.eq("2001") & index.sequence_available.eq("1")).all())
    sequences = np.load(STAGE / "inputs/canonical_sequences.npy", mmap_mode="r", allow_pickle=False)
    n = manifest["encoded_identity_count"]
    audit.check("prepared:array_shape_dtype", sequences.shape == (n,2001) and sequences.dtype == np.uint8)
    by_cache = index.drop_duplicates("cache_row").assign(_row=lambda f:f.cache_row.astype(int)).sort_values("_row")
    audit.check("prepared:dense_cache_rows", np.array_equal(by_cache._row.to_numpy(), np.arange(n)))
    audit.check("prepared:cache_counts", len(index) == manifest["sequence_intervals"] and len(by_cache) == n)
    audit.check("prepared:canonical_identity_bijection", index.groupby("cache_row").canonical_rc_sequence_sha256.nunique().eq(1).all() and by_cache.canonical_rc_sequence_sha256.nunique() == n)
    hashed = [hashlib.sha256(sequences[i].tobytes()).hexdigest() for i in range(n)]
    audit.check("prepared:every_encoded_sequence_SHA", hashed == by_cache.canonical_rc_sequence_sha256.tolist())
    alphabet = np.frombuffer(b"ACGTN", dtype=np.uint8)
    audit.check("prepared:normalized_nucleotides", all(np.isin(sequences[i:i+4096], alphabet).all() for i in range(0,n,4096)))
    audit.check("prepared:forward_orientation_boolean", index.forward_is_canonical.isin(["0", "1"]).all())
    no_ambiguity = index.non_acgt_fraction.astype(float).eq(0)
    audit.check("prepared:unambiguous_forward_mapping", (index.loc[no_ambiguity,"sequence_sha256"].eq(index.loc[no_ambiguity,"canonical_rc_sequence_sha256"]) == index.loc[no_ambiguity,"forward_is_canonical"].eq("1")).all())
    all_frozen = []
    for model in MODELS:
        validation_ids = set()
        for config in CONFIGS:
            original = frame(PREFLIGHT / f"data/configurations/{config}_{model}_interval_manifest.tsv.gz", columns+["label"])
            original = original[original.partition.isin(["train", "validation"])]
            all_frozen.append(original[columns])
            for role in ("train", "checkpoint"):
                actual = frame(STAGE / f"inputs/configurations/{config}_{model}_{role}.tsv.gz", columns+["label"])
                expected = original[original.validation_role.eq(role)]
                audit.check(f"prepared:{config}:{model}:{role}:exact_frozen_rows", actual.reset_index(drop=True).equals(expected.reset_index(drop=True)))
            validation_ids.update(original.loc[original.chrom.eq("chr7"),"interval_id"])
        validation = frame(STAGE / f"inputs/validation/{model}_chr7.tsv.gz", ["interval_id", "chrom"])
        audit.check("prepared:"+model+":validation_union", set(validation.interval_id) == validation_ids and not validation.interval_id.duplicated().any() and validation.chrom.eq("chr7").all())
        panel = frame(PREFLIGHT / f"data/evaluation/{model}_common_challenge_panel.tsv.gz")
        for role in ("selection", "calibration"):
            actual = frame(STAGE / f"inputs/common/{model}_{role}.tsv.gz")
            expected = panel[panel.validation_role.eq(role)]
            audit.check(f"prepared:{model}:{role}:exact_frozen_common", actual.reset_index(drop=True).equals(expected.reset_index(drop=True)))
    frozen = pd.concat(all_frozen).drop_duplicates().set_index("interval_id").sort_index()
    actual = index[columns].set_index("interval_id").sort_index()
    audit.check("prepared:all_sequence_metadata_exact_frozen", frozen.equals(actual))
    for key in ("forbidden_sequence_fetches", "missing_sequences", "hash_failures"):
        audit.check("prepared:"+key, manifest[key] == 0)


def validate_cache(audit):
    import numpy as np
    cache = STAGE / "cache"
    manifest = json.loads((cache / "cache_manifest.json").read_text())
    inputs = json.loads((STAGE / "inputs/input_manifest.json").read_text())
    audit.check("cache:complete", manifest["completed"] is True and manifest["status"] == "PASS")
    audit.check("cache:weight_identity", manifest["phase_I_weights_sha256"] == WEIGHT_SHA)
    audit.check("cache:independent_nucleotide_RC", manifest["nucleotide_rc_independently_evaluated"] is True and manifest["representation"] == "final4560sigmoid")
    n = inputs["encoded_identity_count"]
    audit.check("cache:exact_identity_counts", manifest["encoded_identity_count"] == n and manifest["sequence_intervals"] == inputs["sequence_intervals"] and manifest["orientation_network_evaluations"] == n*2)
    for row in manifest["files"]:
        audit.opaque_hash("cache_artifact:"+row["path"], STAGE / row["path"], row["sha256"], row["bytes"])
    for name in ("features_canonical.npy", "features_rc.npy"):
        value = np.load(cache / name, mmap_mode="r", allow_pickle=False)
        audit.check("cache:"+name+":shape_dtype", value.shape == (n,4560) and value.dtype == np.float32)
        passed = True
        for start in range(0,n,4096):
            block = value[start:start+4096]
            if not (np.isfinite(block).all() and (block >= 0).all() and (block <= 1).all()):
                passed = False
                break
        audit.check("cache:"+name+":all_finite_sigmoid_representations", passed)
    for key in ("test_sequence_extractions", "failures", "missing_sequences"):
        audit.check("cache:"+key, manifest[key] == 0)


def expected_checkpoint(losses):
    """Earliest exact minimum; epochs are one-based."""
    import numpy as np
    losses = np.asarray(losses, dtype=np.float64)
    if not len(losses) or len(losses) > 50 or not np.isfinite(losses).all():
        raise ValueError("Invalid training history")
    return int(np.argmin(losses)) + 1


def symmetric_ensemble(forward, reverse):
    """Independent fixed-order float64 mean, rows x exactly three seeds."""
    import numpy as np
    forward = np.asarray(forward, dtype=np.float64)
    reverse = np.asarray(reverse, dtype=np.float64)
    if forward.shape != reverse.shape or forward.ndim != 2 or forward.shape[1] != 3:
        raise ValueError("Expected same N x 3 orientation matrices")
    if not all(np.isfinite(x).all() and (x >= 0).all() and (x <= 1).all() for x in (forward, reverse)):
        raise ValueError("Invalid probability")
    return ((forward + reverse) / 2).mean(axis=1, dtype=np.float64)


def calibration_threshold(scores):
    import numpy as np
    scores = np.asarray(scores, dtype=np.float64)
    if not len(scores) or not np.isfinite(scores).all() or not ((scores >= 0) & (scores <= 1)).all():
        raise ValueError("Invalid controls")
    k = int(np.floor(0.05 * len(scores)))
    boundary = np.sort(scores)[::-1][k]
    threshold = float(np.nextafter(boundary, np.inf))
    return {"n": len(scores), "k": k, "boundary": float(boundary), "threshold": threshold,
            "threshold_hex": threshold.hex(), "threshold_decimal": format(threshold, ".17g"),
            "calls": int((scores >= threshold).sum())}


def validate_runs(audit):
    """Reconstruct checkpoint/RNG decisions from raw histories, without models."""
    import numpy as np
    expected = {f"{config}_{model}_seed{seed}": (config, model, seed)
                for config in CONFIGS for model in MODELS for seed in SEEDS}
    observed = {path.name for path in (STAGE / "runs").iterdir() if path.is_dir()}
    audit.check("exact_eighteen_run_ids", set(expected) == observed, sorted(observed), sorted(expected))
    initial_digests = {seed: set() for seed in SEEDS}
    architecture_contract = json.loads((BASE / "provenance/COPD-V2-PREFLIGHT_architecture_compute_audit.json").read_text())["phase_two"]
    for run_id, (config, model, seed) in expected.items():
        prefix = "run:" + run_id
        attempts = sorted((STAGE / "runs" / run_id).glob("attempt-*"))
        completed = [path for path in attempts if (path / "completed.json").exists()]
        audit.check(prefix + ":completed_exactly_once", len(completed) == 1, len(completed), 1)
        for index, attempt in enumerate(attempts):
            audit.check(prefix + f":attempt_{index+1}_sequence", attempt.name == f"attempt-{index+1:03d}")
            started = json.loads((attempt / "started.json").read_text())
            audit.check(prefix + f":attempt_{index+1}_seed", int(started["seed"]) == seed)
            if index:
                retry = started.get("retry") or {}
                audit.check(prefix + f":attempt_{index+1}_infrastructure_only",
                            retry.get("reason_class") == "INFRASTRUCTURE" and retry.get("same_seed_identical_procedure") is True and retry.get("authorized") is True)
            for key in ("selection_role_access", "calibration_role_access", "test_access", "external_benchmark_access"):
                audit.check(prefix + f":attempt_{index+1}_firewall:" + key, started.get(key) is False)
        if len(completed) != 1:
            continue
        output = completed[0]
        record = json.loads((output / "selected_checkpoint.json").read_text())
        history = read_tsv(output / "history.tsv")
        contract = json.loads((output / "training_contract.json").read_text())
        initial_digests[seed].add(contract["initial_weights_digest"])
        architecture = json.loads((output / "architecture.json").read_text())
        layers = architecture["config"]["layers"]
        input_layers = [row for row in layers if row["class_name"] == "InputLayer"]
        audit.check(prefix + ":actual_input_geometry", len(input_layers) == 1 and input_layers[0]["config"]["batch_shape"] == [None,4560,1])
        actual_layers = [row for row in layers if row["class_name"] != "InputLayer"]
        audit.check(prefix + ":actual_layer_sequence", [row["class_name"] for row in actual_layers] == [row["type"] for row in architecture_contract["layers"]])
        for layer_number, (actual_layer, expected_layer) in enumerate(zip(actual_layers, architecture_contract["layers"])):
            properties = actual_layer["config"]
            for key, value in expected_layer.items():
                if key == "type":
                    continue
                observed_value = properties.get(key)
                if isinstance(observed_value, list) and len(observed_value) == 1:
                    observed_value = observed_value[0]
                audit.check(prefix + f":actual_layer{layer_number}:" + key, observed_value == value, observed_value, value)
        compiled = architecture["compile_config"]
        audit.check(prefix + ":actual_optimizer_class", compiled["optimizer"]["class_name"] == "Adadelta")
        loss = compiled["loss"]
        audit.check(prefix + ":actual_unweighted_BCE", loss["class_name"] == "BinaryCrossentropy" and loss["config"]["from_logits"] is False and loss["config"]["label_smoothing"] == 0 and loss["config"]["reduction"] == "sum_over_batch_size")
        environment = json.loads((output / "environment.json").read_text())
        audit.check(prefix + ":frozen_python", environment.get("python") == "3.13.0", environment.get("python"), "3.13.0")
        for package, version in (("tensorflow", "2.20.0"), ("keras", "3.14.1"), ("numpy", "2.5.0")):
            audit.check(prefix + ":software:" + package, environment["packages"][package] == version, environment["packages"][package], version)
        audit.check(prefix + ":history_epoch_sequence", [int(row["epoch_one_based"]) for row in history] == list(range(1, len(history)+1)) and [int(row["epoch_zero_based"]) for row in history] == list(range(len(history))))
        losses = [float(row["checkpoint_symmetric_BCE"]) for row in history]
        try:
            best_epoch = expected_checkpoint(losses)
        except ValueError as exc:
            audit.check(prefix + ":valid_history", False, exc)
            continue
        audit.check(prefix + ":earliest_minimum_checkpoint", record["selected_epoch_one_based"] == best_epoch, record["selected_epoch_one_based"], best_epoch)
        audit.check(prefix + ":selected_minimum_BCE", record["checkpoint_symmetric_BCE"] == min(losses), record["checkpoint_symmetric_BCE"], min(losses))
        audit.check(prefix + ":restored_checkpoint_BCE", record["restored_checkpoint_symmetric_BCE"] == min(losses))
        audit.check(prefix + ":serialized_selected_weights", record["serialized_selected_weights_bitwise_roundtrip"] is True and record["checkpoint_is_inference_only"] is True)
        audit.check(prefix + ":epochs_completed", record["epochs_completed"] == len(history))
        best, wait = math.inf, 0
        improvements = []
        n = int(contract["n_train"])
        for index, row in enumerate(history):
            loss = losses[index]
            is_improvement = loss < best
            if is_improvement:
                best, wait = loss, 0
                improvements.append(index+1)
            else:
                wait += 1
            audit.check(prefix + f":epoch{index+1}:improvement", int(row["selected_new_best"]) == int(is_improvement))
            audit.check(prefix + f":epoch{index+1}:patience", int(row["wait"]) == wait and (wait < 15 or index == len(history)-1))
            audit.check(prefix + f":epoch{index+1}:unweighted_exposure", int(row["n_train"]) == n and int(row["steps"]) == math.ceil(n / 256))
            order = np.random.default_rng(np.random.SeedSequence([seed, 1, index])).permutation(n)
            reverse = np.random.default_rng(np.random.SeedSequence([seed, 2, index])).random(n) < .5
            audit.check(prefix + f":epoch{index+1}:shuffle_rng", hashlib.sha256(order.tobytes()).hexdigest() == row["permutation_sha256"])
            audit.check(prefix + f":epoch{index+1}:orientation_rng", hashlib.sha256(reverse.tobytes()).hexdigest() == row["orientation_sha256"] and int(reverse.sum()) == int(row["n_reverse_views"]))
            audit.check(prefix + f":epoch{index+1}:finite_training_loss", math.isfinite(float(row["training_BCE"])))
        audit.check(prefix + ":stopped_only_frozen_rule", len(history) == 50 or wait == 15, {"epochs": len(history), "wait": wait}, "50 epochs or patience15")
        improvement_rows = read_tsv(output / "checkpoint_improvements.tsv")
        audit.check(prefix + ":improvement_ledger", [int(row["epoch_one_based"]) for row in improvement_rows] == improvements)
        for key, expected_value in (("batch_size", 256), ("max_epochs", 50), ("patience", 15), ("min_delta", 0), ("unweighted_BCE", True), ("mixed_precision", False), ("jit_compile", False)):
            audit.check(prefix + ":contract:" + key, contract[key] == expected_value, contract[key], expected_value)
        optimizer = contract["optimizer"]
        audit.check(prefix + ":optimizer_lr", abs(float(optimizer["learning_rate"]) - .001) < 1e-10)
        audit.check(prefix + ":optimizer_rho_epsilon", optimizer["rho"] == .95 and optimizer["epsilon"] == 1e-7)
        for key in ("selection_role_access", "calibration_role_access", "chr8_chr9_access", "external_benchmark_access"):
            audit.check(prefix + ":checkpoint_firewall:" + key, record[key] is False)
        for key in ("checkpoint", "history"):
            item = record[key]
            audit.opaque_hash(prefix + ":" + key, STAGE / item["path"], item["sha256"], item["bytes"])
    for seed, digests in initial_digests.items():
        audit.check(f"initialization:seed{seed}:identical_across_six_contexts", len(digests) == 1, len(digests), 1)


def validate_predictions(audit, directory):
    import numpy as np
    directory = Path(directory)
    evidence = json.loads((directory / "real_network_invariance.json").read_text())
    audit.check("invariance:real_network", evidence["status"] == "PASS" and evidence["real_network"] is True)
    audit.check("invariance:all_audits", evidence["n_seed_audits"] == 18 and evidence["n_ensemble_audits"] == 6 and len(evidence["audit_rows"]) == 24)
    for item in evidence["files"]:
        audit.opaque_hash("prediction_artifact:" + item["path"], STAGE / item["path"], item["sha256"], item["bytes"])
    lookup = {(row["configuration"], row["model"], str(row["seed"])): row for row in evidence["audit_rows"]}
    audit.check("invariance:unique_audits", len(lookup) == 24)
    for model in MODELS:
        original = read_tsv(directory / f"{model}_chr7.tsv.gz")
        repeated = read_tsv(directory / f"{model}_independent_rc_pass.tsv.gz")
        frozen = read_tsv(STAGE / f"inputs/validation/{model}_chr7.tsv.gz")
        keys = [row["interval_id"] for row in original]
        audit.check(model + ":all_validation_IDs", keys == [row["interval_id"] for row in frozen] == [row["interval_id"] for row in repeated])
        audit.check(model + ":chr7_only_predictions", all(row["chrom"] == "chr7" for row in original + repeated))
        audit.check(model + ":unique_prediction_IDs", len(set(keys)) == len(keys))
        for config in CONFIGS:
            f = np.array([[float(row[f"{config}_p_forward_{seed}"]) for seed in SEEDS] for row in original])
            r = np.array([[float(row[f"{config}_p_rc_{seed}"]) for seed in SEEDS] for row in original])
            rf = np.array([[float(row[f"{config}_rc_input_first_{seed}"]) for seed in SEEDS] for row in repeated])
            rr = np.array([[float(row[f"{config}_rc_input_second_{seed}"]) for seed in SEEDS] for row in repeated])
            q, rq = (f + r) / 2, (rf + rr) / 2
            ensemble = symmetric_ensemble(f, r)
            repeated_ensemble = symmetric_ensemble(rf, rr)
            prefix = model + ":" + config
            for index, seed in enumerate(SEEDS):
                audit.check(prefix + f":q{seed}:float64_reduction", np.array_equal(q[:, index], [float(row[f"{config}_q_{seed}"]) for row in original]))
                audit.check(prefix + f":qrc{seed}:float64_reduction", np.array_equal(rq[:, index], [float(row[f"{config}_q_rc_input_{seed}"]) for row in repeated]))
            audit.check(prefix + ":ensemble_float64_fixed_order", np.array_equal(ensemble, [float(row[f"{config}_ensemble_q"]) for row in original]))
            audit.check(prefix + ":reverse_ensemble_float64_fixed_order", np.array_equal(repeated_ensemble, [float(row[f"{config}_ensemble_q_rc_input"]) for row in repeated]))
            pairs = [(str(seed), q[:, index], rq[:, index]) for index, seed in enumerate(SEEDS)] + [("all3", ensemble, repeated_ensemble)]
            for seed, a, b in pairs:
                row = lookup[config, model, seed]
                residual = np.abs(a-b)
                limits = 1e-6 + 1e-6*np.abs(b)
                audit.check(prefix + ":invariance:" + seed, np.all(residual <= limits))
                audit.check(prefix + ":invariance_record:" + seed, row["status"] == "PASS" and row["real_network"] is True and row["n_sequences"] == len(keys) and row["n_failed"] == 0 and row["atol"] == 1e-6 and row["rtol"] == 1e-6 and row["max_absolute_difference"] == float(residual.max()))
                audit.check(prefix + ":original_array_hash:" + seed, row["original_wrapper_array_sha256"] == hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest())
                audit.check(prefix + ":reverse_array_hash:" + seed, row["reverse_wrapper_array_sha256"] == hashlib.sha256(np.ascontiguousarray(b).tobytes()).hexdigest())
    for key in ("selection_metrics_computed", "calibration_performed", "chr8_chr9_access", "external_benchmark_access"):
        audit.check("invariance:firewall:" + key, evidence[key] is False)


def independent_metrics(labels, scores, weights=None):
    """Independent reference uses sklearn, not the deployed custom evaluator."""
    import numpy as np
    from sklearn.metrics import average_precision_score, roc_auc_score
    y, p = np.asarray(labels, dtype=np.int64), np.asarray(scores, dtype=np.float64)
    if weights is None:
        weights = np.ones(len(y))
    weights = np.asarray(weights, dtype=np.float64)
    prevalence = float(np.average(y, weights=weights))
    if not 0 < prevalence < 1:
        raise ValueError("One-class draw")
    brier = float(np.average((p-y)**2, weights=weights))
    ap = float(average_precision_score(y, p, sample_weight=weights))
    return {"AP": ap, "AUROC": float(roc_auc_score(y, p, sample_weight=weights)),
            "Brier": brier, "BrierSkill": 1-brier/(prevalence*(1-prevalence)),
            "AP_gain": ap-prevalence, "prevalence": prevalence}


def validate_evaluation(audit, directory):
    import numpy as np
    directory = Path(directory)
    metrics = read_tsv(directory / "metrics.tsv")
    scores = read_tsv(directory / "common_selection_scores.tsv.gz")
    replicates = read_tsv(directory / "bootstrap_replicates.tsv.gz")
    summaries = read_tsv(directory / "bootstrap_summary.tsv")
    comparisons = read_tsv(directory / "paired_ablation.tsv")
    stability = read_tsv(directory / "seed_stability.tsv")
    bootstrap_audit = json.loads((directory / "bootstrap_audit.json").read_text())
    decision = json.loads((directory / "adequacy_decision.json").read_text())
    audit.check("eval:all_metric_rows", len(metrics) == 24)
    audit.check("eval:chr7_selection_only", all(row["chrom"] == "chr7" and row["validation_role"] == "selection" for row in scores))
    close = lambda a, b: bool(np.isclose(float(a), float(b), rtol=1e-11, atol=1e-12))
    for model in MODELS:
        prefix = "evaluation:" + model
        frozen_rows = [row for row in read_tsv(PREFLIGHT / f"data/evaluation/{model}_common_challenge_panel.tsv.gz") if row["chrom"] == "chr7" and row["validation_role"] == "selection"]
        frozen = {row["interval_id"]: row for row in frozen_rows}
        per_config = {config: [row for row in scores if row["model"] == model and row["configuration"] == config] for config in CONFIGS}
        ids = [row["interval_id"] for row in per_config["V2-C"]]
        audit.check(prefix + ":frozen_selection_IDs", set(ids) == set(frozen) and len(ids) == len(frozen))
        y = np.array([int(frozen[key]["label"]) for key in ids])
        components = np.array([frozen[key]["component_id"] for key in ids])
        all_points, ens = {}, {}
        c_seed_aps = []
        for config in CONFIGS:
            rows = per_config[config]
            audit.check(prefix + ":" + config + ":same_ID_order_and_labels", [row["interval_id"] for row in rows] == ids and all(int(row["label"]) == int(frozen[row["interval_id"]]["label"]) and row["component_id"] == frozen[row["interval_id"]]["component_id"] for row in rows))
            seed_scores = np.array([[float(row[f"q_{seed}"]) for seed in SEEDS] for row in rows])
            ensemble = seed_scores.mean(axis=1, dtype=np.float64)
            audit.check(prefix + ":" + config + ":ensemble_exact", np.array_equal(ensemble, [float(row["ensemble_q"]) for row in rows]))
            ens[config] = ensemble
            seed_aps = []
            for predictor, p in [(f"seed:{seed}", seed_scores[:,i]) for i, seed in enumerate(SEEDS)] + [("ensemble", ensemble)]:
                recomputed = independent_metrics(y, p)
                stored = [r for r in metrics if r["model"] == model and r["configuration"] == config and r["predictor"] == predictor]
                audit.check(prefix + ":" + config + ":metric_membership:" + predictor, len(stored) == 1)
                if len(stored) == 1:
                    audit.check(prefix + ":" + config + ":sklearn_metrics:" + predictor, all(close(stored[0][k], v) for k, v in recomputed.items()))
                if predictor == "ensemble":
                    all_points[config] = recomputed
                else:
                    seed_aps.append(recomputed["AP"])
            st = next(r for r in stability if r["model"] == model and r["configuration"] == config)
            audit.check(prefix + ":" + config + ":sample_seed_SD", close(st["AP_sample_SD"], np.std(seed_aps, ddof=1)) and close(st["AP_range"], np.ptp(seed_aps)) and int(st["seeds_complete"]) == 3)
            if config == "V2-C":
                c_seed_aps = seed_aps
        represented = sorted(set(components))
        inverse = np.searchsorted(represented, components)
        rng = np.random.default_rng(314159)
        attempts, valid, invalid = 0, 0, 0
        replay = {config: [] for config in CONFIGS}
        stored_lookup = {(r["configuration"], int(r["valid_replicate"])): r for r in replicates if r["model"] == model}
        while valid < 2000 and attempts < 20000:
            attempts += 1
            draw = rng.choice(len(represented), size=len(represented), replace=True)
            counts = np.bincount(draw, minlength=len(represented))
            weights = counts[inverse]
            try:
                values = {config: independent_metrics(y, ens[config], weights) for config in CONFIGS}
            except ValueError:
                invalid += 1
                continue
            valid += 1
            for config in CONFIGS:
                row = stored_lookup.get((config, valid))
                audit.check(prefix + f":bootstrap:{config}:{valid}", row is not None and int(row["attempt"]) == attempts and all(close(row[k], v) for k, v in values[config].items()))
                replay[config].append(values[config])
        ba = bootstrap_audit[model]
        audit.check(prefix + ":bootstrap_accounting", ba["attempted"] == attempts and ba["valid"] == valid and ba["invalid"] == invalid and ba["seed"] == 314159 and close(ba["invalid_fraction"], invalid/attempts))
        audit.check(prefix + ":all_bootstrap_records", len(stored_lookup) == 3*valid)
        limits = {}
        for config in CONFIGS:
            for metric in all_points[config]:
                quantiles = np.percentile([r[metric] for r in replay[config]], [2.5,97.5,5,95], method="linear")
                limits[config, metric] = quantiles
                stored = next(r for r in summaries if r["model"] == model and r["configuration"] == config and r["metric"] == metric)
                audit.check(prefix + ":CI:" + config + ":" + metric, all(close(stored[key], value) for key, value in zip(("lower95", "upper95", "lower90", "upper90"), quantiles)))
        for lhs, rhs in (("V2-B", "V2-A"), ("V2-C", "V2-B"), ("V2-C", "V2-A")):
            for metric in ("AP", "AUROC", "Brier"):
                delta = [a[metric]-b[metric] for a,b in zip(replay[lhs], replay[rhs])]
                quantiles = np.percentile(delta, [2.5,97.5,5,95], method="linear")
                stored = next(r for r in comparisons if r["model"] == model and r["comparison"] == f"{lhs} minus {rhs}" and r["metric"] == metric)
                audit.check(prefix + ":paired_CI:" + lhs + ":" + rhs + ":" + metric, close(stored["point"], all_points[lhs][metric]-all_points[rhs][metric]) and all(close(stored[key], value) for key,value in zip(("lower95", "upper95", "lower90", "upper90"), quantiles)))
        count_ok = int(y.sum()) >= 100 and int((1-y).sum()) >= 200 and len(represented) >= 30
        bootstrap_ok = valid == 2000 and invalid/attempts < .1
        performance_ok = all(limits["V2-C", metric][0] > boundary for metric,boundary in (("AUROC", .5), ("AP_gain", 0), ("BrierSkill", 0)))
        seed_ok = np.std(c_seed_aps, ddof=1) <= .03 and np.ptp(c_seed_aps) <= .1
        expected_status = "INCONCLUSIVE" if not (count_ok and bootstrap_ok) else ("PASS" if performance_ok and seed_ok else "FAIL")
        audit.check(prefix + ":absolute_C_decision", decision["model_contexts"][model]["status"] == expected_status, decision["model_contexts"][model]["status"], expected_status)
    audit.check("eval:both_context_gate", decision["both_C_contexts_pass"] == all(decision["model_contexts"][m]["status"] == "PASS" for m in MODELS))
    for key in ("A_B_ranking_used_for_retention", "test_access_authorized", "external_access_authorized", "candidate_scoring_authorized"):
        audit.check("eval:firewall:" + key, decision[key] is False)
    validate_descriptive(audit, directory)


def validate_descriptive(audit, directory):
    """Recompute descriptive bins and train-derived strata independently."""
    import numpy as np
    import pandas as pd
    directory = Path(directory)
    scores = read_tsv(directory / "common_selection_scores.tsv.gz")
    bins = read_tsv(directory / "reliability_bins.tsv")
    strata = read_tsv(directory / "stratified_metrics.tsv")
    cutpoints = json.loads((directory / "stratification_cutpoints.json").read_text())
    covariates = ("gc_fraction", "atac_signal_percentile_max_train_only", "repeat_fraction_2001")
    close = lambda a,b: bool(np.isclose(float(a),float(b),rtol=1e-11,atol=1e-12))
    bin_lookup = {(r["model"],r["configuration"],r["predictor"],int(r["bin"])):r for r in bins}
    stratum_lookup = {(r["model"],r["configuration"],r["predictor"],r["stratifier"],r["stratum"]):r for r in strata}
    audit.check("descriptive:unique_reliability_rows", len(bin_lookup) == len(bins) == 240)
    audit.check("descriptive:unique_stratum_rows", len(stratum_lookup) == len(strata))
    for model in MODELS:
        train = pd.read_csv(PREFLIGHT / f"data/configurations/V2-C_{model}_interval_manifest.tsv.gz", sep="\t", dtype=str, keep_default_na=False, usecols=["chrom","partition","label",*covariates])
        train = train[train.partition.eq("train") & train.label.eq("1")]
        audit.check("descriptive:"+model+":train_only_cutpoint_scope", not train.chrom.isin(["chr7","chr8","chr9"]).any())
        independent_cuts = {name:np.quantile(train[name].astype(float).to_numpy(),[.2,.4,.6,.8],method="linear") for name in covariates}
        for name in covariates:
            audit.check("descriptive:"+model+":training_cutpoints:"+name, np.allclose(independent_cuts[name],cutpoints[model]["cutpoints"][name],rtol=0,atol=1e-14))
        meta = {r["interval_id"]:r for r in read_tsv(STAGE / f"inputs/common/{model}_selection.tsv.gz")}
        for config in CONFIGS:
            rows = [r for r in scores if r["model"] == model and r["configuration"] == config]
            ids = [r["interval_id"] for r in rows]
            y = np.array([int(meta[key]["label"]) for key in ids])
            component = np.array([meta[key]["component_id"] for key in ids])
            categories = {}
            for name in covariates:
                values = np.array([float(meta[key][name]) for key in ids])
                group = np.searchsorted(independent_cuts[name],values,side="right")+1
                categories[name] = np.array([f"quintile_{i}" for i in group])
            categories["atac_lobe_signature"] = np.array([meta[key]["atac_lobes"] or "none" for key in ids])
            categories["same_lobe_mark_support"] = np.array([meta[key][f"{model}_same_lobe_support_lobes"] or "none" for key in ids])
            categories["chromosome"] = np.array([meta[key]["chrom"] for key in ids])
            categories["ambiguous_base_status"] = np.array(["present" if float(meta[key]["non_acgt_fraction"]) > 0 else "absent" for key in ids])
            for predictor, field in [(f"seed:{seed}",f"q_{seed}") for seed in SEEDS]+[("ensemble","ensemble_q")]:
                p = np.array([float(r[field]) for r in rows])
                prefix = f"descriptive:{model}:{config}:{predictor}"
                for i in range(10):
                    mask = (p >= i/10) & ((p <= 1) if i == 9 else (p < (i+1)/10))
                    stored = bin_lookup[model,config,predictor,i]
                    n = int(mask.sum())
                    valid = int(stored["n"]) == n
                    if n:
                        valid = valid and close(stored["mean_probability"],p[mask].mean()) and close(stored["observed_positive_fraction"],y[mask].mean())
                    else:
                        valid = valid and stored["mean_probability"] == "" and stored["observed_positive_fraction"] == ""
                    audit.check(prefix+f":reliability_bin{i}",valid)
                for name, values in categories.items():
                    for value in sorted(set(values)):
                        mask = values == value
                        stored = stratum_lookup[model,config,predictor,name,value]
                        n = int(mask.sum())
                        audit.check(prefix+":"+name+":"+value+":counts", int(stored["n"]) == n and int(stored["n_positive"]) == int(y[mask].sum()) and int(stored["n_control"]) == int((1-y[mask]).sum()) and int(stored["n_components"]) == len(set(component[mask])))
                        audit.check(prefix+":"+name+":"+value+":descriptive_only", stored["all_stratified_results_descriptive"] == "True" and stored["cells_lt100_descriptive_only"] == str(n < 100))
                        try:
                            reference = independent_metrics(y[mask],p[mask])
                            audit.check(prefix+":"+name+":"+value+":sklearn_metrics", stored["metrics_defined"] == "True" and all(close(stored[k],v) for k,v in reference.items()))
                        except ValueError:
                            audit.check(prefix+":"+name+":"+value+":undefined_one_class", stored["metrics_defined"] == "False" and all(stored[k] == "" for k in ("AP","AUROC","Brier","BrierSkill","AP_gain")) and close(stored["prevalence"],y[mask].mean()))


def validate_calibration(audit, directory, evaluation_directory):
    import numpy as np
    directory = Path(directory)
    decision = json.loads((Path(evaluation_directory) / "adequacy_decision.json").read_text())
    both_pass = decision["both_C_contexts_pass"] is True
    if not both_pass:
        audit.check("calibration:prohibited_after_C_failure", not directory.exists() or not any(directory.iterdir()))
        audit.check("calibration:no_release_after_C_failure", not (STAGE / "provenance/C_release_precalibration_freeze.json").exists())
        return
    release = json.loads((STAGE / "provenance/C_release_precalibration_freeze.json").read_text())
    audit.check("calibration:release_status", release["status"] == "C_READY_FOR_CALIBRATION_FROZEN")
    thresholds = json.loads((directory / "C_region_thresholds.json").read_text())
    audit.check("calibration:both_thresholds_PASS", thresholds["status"] == "PASS" and set(thresholds["models"]) == set(MODELS))
    rows = read_tsv(directory / "calibration_negative_scores.tsv.gz")
    audit.check("calibration:negative_chr7_only", all(r["chrom"] == "chr7" and r["validation_role"] == "calibration" and int(r["label"]) == 0 and r["configuration"] == "V2-C" for r in rows))
    for model in MODELS:
        prefix = "calibration:" + model
        selected = [r for r in rows if r["model"] == model]
        frozen = [r for r in read_tsv(PREFLIGHT / f"data/evaluation/{model}_common_challenge_panel.tsv.gz") if r["chrom"] == "chr7" and r["validation_role"] == "calibration" and int(r["label"]) == 0]
        audit.check(prefix + ":exact_frozen_negatives", {r["interval_id"] for r in selected} == {r["interval_id"] for r in frozen} and len(selected) == len(frozen))
        scores = np.array([float(r["ensemble_q"]) for r in selected])
        ensembles = np.mean([[float(r[f"q_{seed}"]) for seed in SEEDS] for r in selected], axis=1, dtype=np.float64)
        audit.check(prefix + ":equal_seed_mean", np.array_equal(ensembles, scores))
        independent = calibration_threshold(scores)
        stored = thresholds["models"][model]
        components = len({r["component_id"] for r in selected})
        audit.check(prefix + ":support_minima", len(selected) >= 200 and components >= 30)
        audit.check(prefix + ":control_counts", stored["n_controls"] == len(selected) and stored["n_control_components"] == components)
        audit.check(prefix + ":boundary_exact", stored["sorted_boundary_score"] == independent["boundary"] and stored["zero_based_boundary_index"] == independent["k"])
        audit.check(prefix + ":threshold_exact", stored["threshold_float64"] == independent["threshold"] and stored["threshold_hex"] == independent["threshold_hex"] and stored["threshold_decimal_17g"] == independent["threshold_decimal"])
        audit.check(prefix + ":FPR_ties", stored["n_called_controls"] == independent["calls"] and independent["calls"] <= independent["k"] and stored["observed_calibration_row_FPR"] == independent["calls"]/len(scores) and stored["n_at_boundary"] == int((scores == independent["boundary"]).sum()))
        audit.check(prefix + ":no_allele_delta", stored["allele_delta_threshold"] is None and stored["population_FPR_guarantee"] is False and stored["variant_FPR_guarantee"] is False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("phase", choices=["inputs", "baseline", "preservation", "prepared", "cache", "runs", "predictions", "evaluation", "calibration"])
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--prediction-dir", type=Path, default=STAGE / "predictions")
    parser.add_argument("--evaluation-dir", type=Path, default=STAGE / "results/chr7_evaluation")
    parser.add_argument("--calibration-dir", type=Path, default=STAGE / "results/chr7_calibration")
    parser.add_argument("--baseline-manifest", type=Path, default=STAGE / "provenance/tracked_baseline_before.json")
    args = parser.parse_args()
    audit = Audit()
    extra = {}
    if args.phase == "inputs":
        validate_inputs(audit)
    elif args.phase == "baseline":
        extra = tracked_baseline(audit, args.output_prefix)
    elif args.phase == "preservation":
        validate_preservation(audit, args.baseline_manifest)
    elif args.phase == "prepared":
        validate_prepared(audit)
    elif args.phase == "cache":
        validate_cache(audit)
    elif args.phase == "runs":
        validate_runs(audit)
    elif args.phase == "predictions":
        validate_predictions(audit, args.prediction_dir)
    elif args.phase == "evaluation":
        validate_evaluation(audit, args.evaluation_dir)
    elif args.phase == "calibration":
        validate_calibration(audit, args.calibration_dir, args.evaluation_dir)
    failures = audit.write(args.output_prefix, {"phase": args.phase, **extra})
    return int(bool(failures))


if __name__ == "__main__":
    sys.exit(main())
