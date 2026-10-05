#!/usr/bin/env python3
"""Evaluate the two held-out COPD lung TREDNet classifiers with uncertainty."""

from __future__ import annotations

import hashlib
import json
import os
import platform
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
os.environ.setdefault("TF_DETERMINISTIC_OPS", "1")

import h5py
import keras
import numpy as np
import pandas as pd
import tensorflow as tf
from sklearn import metrics


ROOT = Path(__file__).resolve().parents[4]
SECTION = ROOT / "diseases" / "COPD" / "04_modeling"
MODEL_ROOT = SECTION / "trednet" / "models_output"
RESULTS = SECTION / "results"
PERFORMANCE = RESULTS / "COPD-S4-R001_model_performance.tsv"
THRESHOLDS = RESULTS / "COPD-S4-R001_test_threshold_performance.tsv"
PREDICTIONS = RESULTS / "COPD-S4-R001_test_predictions.tsv.gz"
MANIFEST = RESULTS / "COPD-S4-R001_evaluation_manifest.json"

MODELS = {
    "enhancer": "COPD_SevereEmphysema_Lung_Enhancer_DHS_x2",
    "silencer": "COPD_SevereEmphysema_Lung_Silencer_DHS_x2",
}
SEED = 20261001
BOOTSTRAPS = 500


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pr_auc(labels: np.ndarray, scores: np.ndarray) -> float:
    precision, recall, _ = metrics.precision_recall_curve(labels, scores)
    return float(metrics.auc(recall, precision))


def expected_calibration_error(labels: np.ndarray, scores: np.ndarray, bins: int = 10) -> float:
    edges = np.linspace(0.0, 1.0, bins + 1)
    total = len(labels)
    value = 0.0
    for index in range(bins):
        if index == bins - 1:
            selected = (scores >= edges[index]) & (scores <= edges[index + 1])
        else:
            selected = (scores >= edges[index]) & (scores < edges[index + 1])
        if selected.any():
            value += (
                selected.sum()
                / total
                * abs(float(labels[selected].mean()) - float(scores[selected].mean()))
            )
    return float(value)


def bootstrap_intervals(labels: np.ndarray, scores: np.ndarray, seed: int) -> dict[str, float]:
    rng = np.random.default_rng(seed)
    positive = np.flatnonzero(labels == 1)
    negative = np.flatnonzero(labels == 0)
    roc_values = np.empty(BOOTSTRAPS, dtype=float)
    pr_values = np.empty(BOOTSTRAPS, dtype=float)
    for index in range(BOOTSTRAPS):
        sample = np.concatenate(
            (
                rng.choice(positive, size=len(positive), replace=True),
                rng.choice(negative, size=len(negative), replace=True),
            )
        )
        roc_values[index] = metrics.roc_auc_score(labels[sample], scores[sample])
        pr_values[index] = pr_auc(labels[sample], scores[sample])
    return {
        "roc_auc_ci95_low": float(np.percentile(roc_values, 2.5)),
        "roc_auc_ci95_high": float(np.percentile(roc_values, 97.5)),
        "pr_auc_ci95_low": float(np.percentile(pr_values, 2.5)),
        "pr_auc_ci95_high": float(np.percentile(pr_values, 97.5)),
    }


def read_thresholds(path: Path) -> dict[int, float]:
    values = {}
    with path.open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) == 2:
                values[int(fields[0])] = float(fields[1])
    if set(values) != {1, 3, 5, 10}:
        raise ValueError(f"unexpected FPR thresholds in {path}: {values}")
    return values


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    keras.utils.set_random_seed(SEED)
    try:
        tf.config.experimental.enable_op_determinism()
    except Exception as exc:
        print(f"WARNING: deterministic operations unavailable: {exc}")

    performance_rows = []
    threshold_rows = []
    prediction_frames = []
    input_files: dict[str, Path] = {}

    for model_index, (model_type, eid) in enumerate(MODELS.items()):
        directory = MODEL_ROOT / eid
        dataset_path = directory / f"{eid}_phase_two_dataset.hdf5"
        model_path = directory / "phase_two_model.keras"
        weights_path = directory / f"{eid}_phase_two_weights.weights.h5"
        history_path = directory / "training_history.csv"
        thresholds_path = directory / "fpr_threshold_scores.txt"
        input_files.update(
            {
                f"{model_type}_dataset": dataset_path,
                f"{model_type}_model": model_path,
                f"{model_type}_weights": weights_path,
                f"{model_type}_history": history_path,
                f"{model_type}_thresholds": thresholds_path,
            }
        )

        with h5py.File(dataset_path, "r") as dataset:
            test_data = dataset["test_data"][:]
            labels = dataset["test_labels"][:].astype(np.int8)
        model = keras.models.load_model(model_path)
        model.load_weights(weights_path)
        scores = model.predict(
            test_data[..., np.newaxis], batch_size=256, verbose=0
        ).ravel()
        del test_data, model

        roc_auc = float(metrics.roc_auc_score(labels, scores))
        area_pr = pr_auc(labels, scores)
        average_precision = float(metrics.average_precision_score(labels, scores))
        intervals = bootstrap_intervals(labels, scores, SEED + model_index)
        clipped = np.clip(scores, 1e-7, 1 - 1e-7)

        history = pd.read_csv(history_path)
        best_index = int(history["val_loss"].idxmin())
        stored_roc = float((directory / "auc.txt").read_text().split()[1])
        stored_pr = float((directory / "prc.txt").read_text().split()[1])
        if abs(stored_roc - roc_auc) > 5e-6 or abs(stored_pr - area_pr) > 5e-6:
            raise ValueError(
                f"{model_type}: recomputed metrics disagree with stored values: "
                f"ROC {roc_auc} versus {stored_roc}; PR {area_pr} versus {stored_pr}"
            )

        performance_rows.append(
            {
                "model_type": model_type,
                "eid": eid,
                "regulatory_context": "ENCDO520EJG severe-emphysema lung, three-lobe union",
                "genome_build": "GRCh38",
                "test_chromosomes": "chr8;chr9",
                "n_test": len(labels),
                "n_test_positive": int(labels.sum()),
                "n_test_control": int((labels == 0).sum()),
                "positive_fraction": float(labels.mean()),
                "roc_auc": roc_auc,
                **intervals,
                "pr_auc_trapezoid": area_pr,
                "average_precision": average_precision,
                "brier_score": float(metrics.brier_score_loss(labels, scores)),
                "log_loss": float(metrics.log_loss(labels, clipped, labels=[0, 1])),
                "ece_10_equal_width_bins": expected_calibration_error(labels, scores),
                "best_epoch_one_based": best_index + 1,
                "epochs_run": len(history),
                "best_validation_loss": float(history.loc[best_index, "val_loss"]),
                "bootstrap_replicates": BOOTSTRAPS,
                "bootstrap_seed": SEED + model_index,
                "bootstrap_method": "stratified percentile, resampling positive and control test examples separately",
            }
        )

        for target_fpr, threshold in sorted(read_thresholds(thresholds_path).items(), reverse=True):
            calls = scores >= threshold
            observed_fpr = float(calls[labels == 0].mean())
            observed_tpr = float(calls[labels == 1].mean())
            threshold_rows.append(
                {
                    "model_type": model_type,
                    "target_fpr_percent": target_fpr,
                    "score_threshold": threshold,
                    "observed_fpr": observed_fpr,
                    "observed_tpr": observed_tpr,
                    "n_control_called": int(calls[labels == 0].sum()),
                    "n_positive_called": int(calls[labels == 1].sum()),
                }
            )
        prediction_frames.append(
            pd.DataFrame(
                {
                    "model_type": model_type,
                    "test_index": np.arange(len(labels)),
                    "label": labels,
                    "score": scores,
                }
            )
        )
        print(
            f"{model_type}: n={len(labels)}, ROC AUC={roc_auc:.6f}, "
            f"PR AUC={area_pr:.6f}"
        )

    performance = pd.DataFrame(performance_rows)
    thresholds = pd.DataFrame(threshold_rows)
    predictions = pd.concat(prediction_frames, ignore_index=True)
    performance.to_csv(PERFORMANCE, sep="\t", index=False)
    thresholds.to_csv(THRESHOLDS, sep="\t", index=False)
    predictions.to_csv(
        PREDICTIONS,
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )

    manifest = {
        "result_id": "COPD-S4-R001",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "held_out_unit": "chromosomes chr8 and chr9",
        "bootstrap_replicates": BOOTSTRAPS,
        "software": {
            "python": platform.python_version(),
            "tensorflow": tf.__version__,
            "keras": keras.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "inputs": {},
        "outputs": {},
    }
    for label, path in input_files.items():
        manifest["inputs"][label] = {
            "path": str(path),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
    for label, path in {
        "performance": PERFORMANCE,
        "threshold_performance": THRESHOLDS,
        "test_predictions": PREDICTIONS,
        "script": Path(__file__).resolve(),
    }.items():
        manifest["outputs"][label] = {
            "path": str(path),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
