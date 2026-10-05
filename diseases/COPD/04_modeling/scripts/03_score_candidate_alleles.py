#!/usr/bin/env python3
"""Score paired COPD candidate alleles with the trained lung TREDNet models.

The phase-I embedding is computed once per allele batch and reused for both
the enhancer and silencer classifiers.  This is numerically equivalent to
running the standalone TREDNet inference entry point four times, but avoids
duplicated phase-I computation and bounds host memory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import json
import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

os.environ.setdefault("TF_CPP_MIN_LOG_LEVEL", "2")
os.environ.setdefault("TF_XLA_FLAGS", "--tf_xla_enable_xla_devices=false")
os.environ.setdefault("TF_DETERMINISTIC_OPS", "1")

import keras
import numpy as np
import tensorflow as tf
from Bio import SeqIO


ROOT = Path(__file__).resolve().parents[4]
SECTION = ROOT / "diseases" / "COPD" / "04_modeling"
WORKSPACE = SECTION / "trednet"
MODEL_ROOT = WORKSPACE / "models_output"
PHASE_ONE_DIR = WORKSPACE / "model_phase_I"
REF_FASTA = SECTION / "data" / "COPD_candidate_variants_ref_2001bp.fa"
ALT_FASTA = SECTION / "data" / "COPD_candidate_variants_alt_2001bp.fa"
OUTPUT = SECTION / "results" / "COPD-S4-R003_candidate_allele_scores.tsv.gz"
MANIFEST = SECTION / "results" / "COPD-S4-R003_scoring_manifest.json"

MODELS = {
    "enhancer": "COPD_SevereEmphysema_Lung_Enhancer_DHS_x2",
    "silencer": "COPD_SevereEmphysema_Lung_Silencer_DHS_x2",
}
INPUT_LENGTH = 2001
SEED = 20261001


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_inference_module():
    path = ROOT / "models" / "TREDNET_v2" / "TREDNet_v2_inference.py"
    spec = importlib.util.spec_from_file_location("trednet_inference", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module, path


def paired_batches(ref_path: Path, alt_path: Path, batch_size: int, encode):
    ref_iter = SeqIO.parse(ref_path, "fasta")
    alt_iter = SeqIO.parse(alt_path, "fasta")
    names: list[str] = []
    ref_vectors: list[np.ndarray] = []
    alt_vectors: list[np.ndarray] = []
    n_total = 0
    while True:
        ref_record = next(ref_iter, None)
        alt_record = next(alt_iter, None)
        if ref_record is None and alt_record is None:
            break
        if ref_record is None or alt_record is None:
            raise ValueError("REF and ALT FASTA record counts differ")
        if ref_record.id != alt_record.id:
            raise ValueError(
                f"paired FASTA order mismatch: {ref_record.id} != {alt_record.id}"
            )
        ref_sequence = str(ref_record.seq).upper()
        alt_sequence = str(alt_record.seq).upper()
        if len(ref_sequence) != INPUT_LENGTH or len(alt_sequence) != INPUT_LENGTH:
            raise ValueError(f"{ref_record.id}: sequence length is not {INPUT_LENGTH}")
        names.append(ref_record.id)
        ref_vectors.append(encode(ref_sequence))
        alt_vectors.append(encode(alt_sequence))
        n_total += 1
        if len(names) == batch_size:
            yield names, np.stack(ref_vectors), np.stack(alt_vectors)
            names, ref_vectors, alt_vectors = [], [], []
    if names:
        yield names, np.stack(ref_vectors), np.stack(alt_vectors)
    return n_total


def read_thresholds(path: Path) -> dict[str, float]:
    thresholds = {}
    with path.open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) == 2:
                thresholds[f"fpr{fields[0]}"] = float(fields[1])
    required = {"fpr10", "fpr5", "fpr3", "fpr1"}
    if set(thresholds) != required:
        raise ValueError(f"incomplete threshold file {path}: {thresholds}")
    return thresholds


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--batch-size", type=int, default=256)
    parser.add_argument("--prediction-batch-size", type=int, default=32)
    args = parser.parse_args()

    keras.utils.set_random_seed(SEED)
    try:
        tf.config.experimental.enable_op_determinism()
    except Exception as exc:
        print(f"WARNING: deterministic operations unavailable: {exc}")

    inference, inference_path = load_inference_module()
    phase_one = inference.load_phase_one(str(PHASE_ONE_DIR))
    phase_two = {}
    thresholds = {}
    for model_type, eid in MODELS.items():
        save_dir = MODEL_ROOT / eid
        phase_two[model_type] = inference.load_phase_two(str(save_dir), eid)
        thresholds[model_type] = read_thresholds(save_dir / "fpr_threshold_scores.txt")

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    import gzip

    n_scored = 0
    fieldnames = [
        "candidate_record_id",
        "enhancer_ref_score",
        "enhancer_alt_score",
        "enhancer_delta_alt_minus_ref",
        "silencer_ref_score",
        "silencer_alt_score",
        "silencer_delta_alt_minus_ref",
    ]
    with OUTPUT.open("wb") as raw_handle, gzip.GzipFile(
        filename="", fileobj=raw_handle, mode="wb", compresslevel=6, mtime=0
    ) as gzip_handle, io.TextIOWrapper(gzip_handle) as output_handle:
        writer = csv.DictWriter(output_handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for names, ref_batch, alt_batch in paired_batches(
            REF_FASTA, ALT_FASTA, args.batch_size, inference.seq2one_hot
        ):
            both = np.concatenate((ref_batch, alt_batch), axis=0)
            embeddings = phase_one.predict(
                both,
                batch_size=args.prediction_batch_size,
                verbose=0,
            )
            split = len(names)
            ref_embeddings = embeddings[:split, ..., np.newaxis]
            alt_embeddings = embeddings[split:, ..., np.newaxis]
            scores = {}
            for model_type, model in phase_two.items():
                scores[(model_type, "ref")] = model.predict(
                    ref_embeddings,
                    batch_size=args.prediction_batch_size,
                    verbose=0,
                ).ravel()
                scores[(model_type, "alt")] = model.predict(
                    alt_embeddings,
                    batch_size=args.prediction_batch_size,
                    verbose=0,
                ).ravel()
            for index, name in enumerate(names):
                enhancer_ref = float(scores[("enhancer", "ref")][index])
                enhancer_alt = float(scores[("enhancer", "alt")][index])
                silencer_ref = float(scores[("silencer", "ref")][index])
                silencer_alt = float(scores[("silencer", "alt")][index])
                writer.writerow(
                    {
                        "candidate_record_id": name,
                        "enhancer_ref_score": f"{enhancer_ref:.9g}",
                        "enhancer_alt_score": f"{enhancer_alt:.9g}",
                        "enhancer_delta_alt_minus_ref": f"{enhancer_alt - enhancer_ref:.9g}",
                        "silencer_ref_score": f"{silencer_ref:.9g}",
                        "silencer_alt_score": f"{silencer_alt:.9g}",
                        "silencer_delta_alt_minus_ref": f"{silencer_alt - silencer_ref:.9g}",
                    }
                )
            n_scored += len(names)
            print(f"Scored {n_scored} paired candidates", flush=True)
            del both, embeddings, ref_embeddings, alt_embeddings, scores

    manifest = {
        "result_id": "COPD-S4-R003",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "seed": SEED,
        "input_length_bp": INPUT_LENGTH,
        "paired_candidates_scored": n_scored,
        "outer_batch_size": args.batch_size,
        "prediction_batch_size": args.prediction_batch_size,
        "models": MODELS,
        "thresholds": thresholds,
        "software": {
            "python": platform.python_version(),
            "tensorflow": tf.__version__,
            "keras": keras.__version__,
            "numpy": np.__version__,
        },
        "files": {},
    }
    paths = {
        "ref_fasta": REF_FASTA,
        "alt_fasta": ALT_FASTA,
        "scores": OUTPUT,
        "phase_one_weights": PHASE_ONE_DIR / "phase_one_weights.h5",
        "inference_source": inference_path,
        "script": Path(__file__).resolve(),
    }
    for model_type, eid in MODELS.items():
        paths[f"{model_type}_best_weights"] = (
            MODEL_ROOT / eid / f"{eid}_phase_two_weights.weights.h5"
        )
    for label, path in paths.items():
        manifest["files"][label] = {
            "path": str(path),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        }
    MANIFEST.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(f"Wrote {n_scored} paired scores to {OUTPUT}")


if __name__ == "__main__":
    main()
