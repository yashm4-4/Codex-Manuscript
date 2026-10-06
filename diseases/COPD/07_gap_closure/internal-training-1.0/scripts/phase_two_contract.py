#!/usr/bin/env python3
"""Frozen pretraining-1.1 phase-II contract; imports no model runtime itself.

Only stage-local train/chr7 tables are accepted. The two cache matrices are
independent frozen phase-I computations on canonical nucleotide strings and
their nucleotide reverse complements, NEVER reversed representation axes.
"""
from __future__ import annotations

import csv
import datetime as dt
import gzip
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import resource
import subprocess
import sys
import time

import numpy as np

SEEDS = (104729, 130363, 155921)
CONFIGURATIONS = ("V2-A", "V2-B", "V2-C")
MODELS = ("enhancer", "h3k27me3")
TRAIN_CHROMS = tuple(f"chr{i}" for i in [1, 2, 3, 4, 5, 6, *range(10, 23)]) + ("chrX", "chrY")
PHASE_I_SHA256 = "483b6c0cafd750ea8e97932cbc8a949276eaf1e77e7d8766f9a22eeddc2916d6"
ATOL = RTOL = 1e-6
BATCH_SIZE = 256
MAX_EPOCHS = 50
PATIENCE = 15


def utcnow():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    """Exclusive creation: finalized provenance is never silently replaced."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def read_json(path):
    with open(path) as stream:
        return json.load(stream)


def stage_file(stage, relative):
    """Reject manifest path traversal and paths outside this execution bundle."""
    stage = Path(stage).resolve()
    path = (stage / relative).resolve()
    if not path.is_relative_to(stage):
        raise RuntimeError(f"Non-stage input path rejected: {relative}")
    return path


def file_record(stage, path):
    path = Path(path)
    return {"path": str(path.resolve().relative_to(Path(stage).resolve())),
            "bytes": path.stat().st_size, "sha256": sha256(path)}


def verify_record(stage, record):
    path = stage_file(stage, record["path"])
    if path.stat().st_size != int(record["bytes"]) or sha256(path) != record["sha256"]:
        raise RuntimeError(f"Frozen input hash/size mismatch: {path}")
    return path


def runtime_environment(tf=None, keras=None):
    packages = {}
    for name in ("numpy", "tensorflow", "keras", "h5py", "scipy", "scikit-learn"):
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    try:
        gpu = subprocess.run(["nvidia-smi", "--query-gpu=name,uuid,driver_version,memory.total", "--format=csv,noheader"],
                             capture_output=True, text=True, timeout=15, check=False)
        gpu_detail = {"returncode": gpu.returncode, "stdout": gpu.stdout, "stderr": gpu.stderr}
    except (OSError, subprocess.TimeoutExpired) as exc:
        gpu_detail = {"error": str(exc)}
    record = {"python": platform.python_version(), "executable": sys.executable,
              "platform": platform.platform(), "packages": packages,
              "nvidia_smi": gpu_detail,
              "environment": {key: os.environ.get(key) for key in (
                  "SLURM_JOB_ID", "SLURM_ARRAY_TASK_ID", "SLURM_CPUS_PER_TASK",
                  "SLURM_JOB_GPUS", "CUDA_VISIBLE_DEVICES", "PYTHONHASHSEED",
                  "TF_DETERMINISTIC_OPS", "TF_CUDNN_DETERMINISTIC", "TF_XLA_FLAGS",
                  "TF_ENABLE_ONEDNN_OPTS")}}
    if tf is not None:
        record["tensorflow_build"] = tf.sysconfig.get_build_info()
        record["logical_devices"] = [str(device) for device in tf.config.list_logical_devices()]
    if keras is not None:
        record["keras_floatx"] = keras.backend.floatx()
        record["keras_epsilon"] = float(keras.backend.epsilon())
    return record


def initialize_runtime(seed):
    """Exact frozen software and deterministic setup; fail, never substitute."""
    if int(seed) not in SEEDS:
        raise RuntimeError("Training seed outside frozen list")
    for key, value in {"TF_DETERMINISTIC_OPS": "1", "TF_CUDNN_DETERMINISTIC": "1",
                       "TF_XLA_FLAGS": "--tf_xla_enable_xla_devices=false"}.items():
        if key in os.environ and os.environ[key] != value:
            raise RuntimeError(f"Conflicting deterministic runtime setting: {key}")
        os.environ[key] = value
    # The launcher sets this before interpreter initialization. Refuse to pretend
    # that assigning it after startup would make Python hashing deterministic.
    if os.environ.get("PYTHONHASHSEED") != str(seed):
        raise RuntimeError("Launch with PYTHONHASHSEED equal to the frozen training seed")
    import tensorflow as tf
    import keras
    if (platform.python_version(), tf.__version__, keras.__version__, np.__version__) != (
            "3.13.0", "2.20.0", "3.14.1", "2.5.0"):
        raise RuntimeError("Software differs from frozen Python3.13.0/TF2.20.0/Keras3.14.1/NumPy2.5.0")
    keras.backend.set_floatx("float32")
    keras.mixed_precision.set_global_policy("float32")
    keras.utils.set_random_seed(seed)
    tf.config.experimental.enable_op_determinism()
    devices = tf.config.list_physical_devices("GPU")
    if len(devices) != 1:
        raise RuntimeError(f"Exactly one allocated GPU is required, saw {len(devices)}")
    tf.config.experimental.set_memory_growth(devices[0], True)
    return tf, keras


def build_phase_two(keras):
    """Layer-for-layer V1 architecture, new seed-controlled initialization."""
    model = keras.Sequential([
        keras.layers.Input(shape=(4560, 1)),
        keras.layers.Conv1D(64, 4, activation="relu", padding="valid", name="conv1d_1"),
        keras.layers.BatchNormalization(name="batch_normalization_1"),
        keras.layers.MaxPooling1D(pool_size=2, name="max_pooling1d_1"),
        keras.layers.Dropout(0.4, name="dropout_1"),
        keras.layers.Conv1D(128, 2, activation="relu", padding="valid", name="conv1d_2"),
        keras.layers.Dropout(0.4, name="dropout_2"),
        keras.layers.Flatten(name="flatten_1"),
        keras.layers.Dense(100, activation="relu", name="dense_1"),
        keras.layers.Dense(50, activation="relu", name="dense_2"),
        keras.layers.Dense(1, activation="sigmoid", name="dense_3"),
    ])
    model.compile(optimizer=keras.optimizers.Adadelta(learning_rate=0.001, rho=0.95, epsilon=1e-7),
                  loss=keras.losses.BinaryCrossentropy(), jit_compile=False)
    return model


def read_rows(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def bool_value(value):
    if str(value).lower() in ("1", "true"):
        return True
    if str(value).lower() in ("0", "false"):
        return False
    raise ValueError(f"Invalid Boolean: {value!r}")


def load_input_table(path, mode):
    rows = read_rows(path)
    if not rows:
        raise RuntimeError(f"Empty input table: {path}")
    allowed_chroms = set(TRAIN_CHROMS) if mode == "train" else {"chr7"}
    for row in rows:
        if row["chrom"] not in allowed_chroms:
            raise RuntimeError(f"Chromosome firewall breach in {path}: {row['chrom']}")
        if mode == "checkpoint" and row.get("validation_role", row.get("role")) != "checkpoint":
            raise RuntimeError("Checkpoint input contains non-checkpoint role")
    identifiers = [row["interval_id"] for row in rows]
    if len(set(identifiers)) != len(identifiers):
        raise RuntimeError(f"Duplicate interval IDs in {path}")
    result = {"rows": rows, "interval_id": np.array(identifiers),
              "cache_row": np.array([int(row["cache_row"]) for row in rows], dtype=np.int64),
              "forward_is_canonical": np.array([bool_value(row["forward_is_canonical"]) for row in rows])}
    if mode in ("train", "checkpoint"):
        labels = np.array([int(row["label"]) for row in rows], dtype=np.float32)
        if set(labels.tolist()) != {0.0, 1.0}:
            raise RuntimeError("Training/checkpoint panel must contain both binary classes")
        result["label"] = labels
    return result


class FeatureCache:
    def __init__(self, stage, manifest_path="cache/cache_manifest.json"):
        self.stage = Path(stage).resolve()
        self.manifest_path = stage_file(self.stage, manifest_path)
        self.manifest = read_json(self.manifest_path)
        if self.manifest.get("completed") is not True or self.manifest.get("phase_I_weights_sha256") != PHASE_I_SHA256:
            raise RuntimeError("Cache is not complete or uses wrong phase-I weights")
        paths = {Path(record["path"]).name: verify_record(self.stage, record)
                 for record in self.manifest["files"]}
        self.canonical = np.load(paths["features_canonical.npy"], mmap_mode="r", allow_pickle=False)
        self.rc = np.load(paths["features_rc.npy"], mmap_mode="r", allow_pickle=False)
        if self.canonical.dtype != np.float32 or self.rc.dtype != np.float32:
            raise RuntimeError("Features are not frozen float32 representations")
        if self.canonical.shape != self.rc.shape or self.canonical.ndim != 2 or self.canonical.shape[1] != 4560:
            raise RuntimeError("Feature geometry failure")
        if tuple(self.manifest["shape"]) != self.canonical.shape:
            raise RuntimeError("Feature manifest/matrix shape mismatch")

    def check_indices(self, table):
        index = table["cache_row"]
        if np.any(index < 0) or np.any(index >= len(self.canonical)):
            raise RuntimeError("Cache row index outside verified matrices")

    def batch(self, table, indices, reverse=False):
        cache_rows = table["cache_row"][indices]
        choose_canonical = table["forward_is_canonical"][indices] ^ np.asarray(reverse, dtype=bool)
        x = np.empty((len(indices), 4560), dtype=np.float32)
        x[choose_canonical] = self.canonical[cache_rows[choose_canonical]]
        x[~choose_canonical] = self.rc[cache_rows[~choose_canonical]]
        if not np.isfinite(x).all() or np.any(x < 0) or np.any(x > 1):
            raise RuntimeError("Nonfinite/out-of-range phase-I sigmoid representation")
        return x[..., None]


def epoch_plan(seed, epoch, n):
    """Orientation assigned to frozen genomic rows, independently of shuffle."""
    order = np.random.default_rng(np.random.SeedSequence([seed, 1, epoch])).permutation(n)
    reverse = np.random.default_rng(np.random.SeedSequence([seed, 2, epoch])).random(n) < 0.5
    return order, reverse


def valid_probabilities(value):
    value = np.asarray(value)
    if not np.isfinite(value).all() or np.any(value < 0) or np.any(value > 1):
        raise RuntimeError("Nonfinite/out-of-range phase-II probability")


def symmetric_probability(forward, reverse):
    valid_probabilities(forward)
    valid_probabilities(reverse)
    return (np.asarray(forward, dtype=np.float64) + np.asarray(reverse, dtype=np.float64)) / 2.0


def seed_ensemble(ordered_symmetric):
    if len(ordered_symmetric) != 3:
        raise RuntimeError("Exactly three seed probabilities required")
    result = np.mean(np.stack(ordered_symmetric, axis=0).astype(np.float64), axis=0, dtype=np.float64)
    valid_probabilities(result)
    return result


def symmetric_bce(labels, probabilities):
    valid_probabilities(probabilities)
    p = np.clip(np.asarray(probabilities, dtype=np.float64), 1e-7, 1.0 - 1e-7)
    y = np.asarray(labels, dtype=np.float64)
    return float(-np.mean(y * np.log(p) + (1.0 - y) * np.log1p(-p), dtype=np.float64))


def inference_function(tf, model):
    @tf.function(input_signature=[tf.TensorSpec([None, 4560, 1], tf.float32)], jit_compile=False)
    def infer(x):
        return model(x, training=False)
    return infer


def predict_orientation(infer, cache, table, reverse=False, reverse_order=False, batch_size=256):
    n = len(table["cache_row"])
    result = np.empty(n, dtype=np.float64)
    order = np.arange(n, dtype=np.int64)
    if reverse_order:
        order = order[::-1]
    for start in range(0, n, batch_size):
        indices = order[start:start + batch_size]
        p = np.asarray(infer(cache.batch(table, indices, reverse=reverse))).reshape(-1)
        if p.shape != (len(indices),):
            raise RuntimeError("Invalid phase-II prediction shape")
        valid_probabilities(p)
        result[indices] = p.astype(np.float64)
    return result


def array_digest(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def resource_record(started):
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return {"finished_utc": utcnow(), "wall_seconds": time.monotonic() - started,
            "user_cpu_seconds": usage.ru_utime, "system_cpu_seconds": usage.ru_stime,
            "peak_RSS_KiB": usage.ru_maxrss}
