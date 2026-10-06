#!/usr/bin/env python3
"""Inference only: six immutable V2-C archives on frozen benchmark allele inputs.

No experimental label/outcome table is parsed. Every input and annotation gate is
verified before any model is imported or evaluated. Original model archives are
never saved, altered, or copied. The separate real-network RC audit repeats both
phase-I and phase-II calls, not merely the arithmetic on stored probabilities.
"""
from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import importlib.util
import json
from pathlib import Path
import resource
import sys
import time
import traceback

import numpy as np

sys.dont_write_bytecode = True
STAGE_NAME = "external-benchmark-v2-1.0"
TRAINING = "diseases/COPD/07_gap_closure/internal-training-1.0"
SEEDS = (104729, 130363, 155921)
CONTEXTS = ("enhancer", "h3k27me3")
ATOL = RTOL = 1e-6
PROSPECTIVE_SHA = "d4725233251785d90b3a768a115452e125fe5323a58768df5047043f6f020830"
SPECIFICATION_SHA = "32432be4a7698256b8ae3cb6fe4139eac8a1ad6c83c5b8de41b880ce771e11bc"
PHASE_ONE_SHA = "483b6c0cafd750ea8e97932cbc8a949276eaf1e77e7d8766f9a22eeddc2916d6"
HELPER_SHA = "f50f226b901df6bbd66455462e4da922d8229a4d030f09613a8483a74f636e2f"
ARCHITECTURE_SHA = "d1b77cf5a8468ee3516d6a89bcfd4f1e7994087642c851f409425242408621dd"
THRESHOLD_SHA = "48f9a7a742a32ebb39d21ecfcde71c340cb2702ef124467d84ce40d0e33037a8"
CHECKPOINT_SHA = {
    ("enhancer", 104729): "1cb670f778a2d148d0ac2b71bea5c97f5a421bb64f744af1dc5b373a60ea3114",
    ("enhancer", 130363): "0e42ecc9364f967220ecc6c3bb0007dc25aced5d9c360355f6e48027c39f0795",
    ("enhancer", 155921): "733ca1845061508c820240d04d4b7c30d8ba0fa7c1086419eaec0b5c2bb041d5",
    ("h3k27me3", 104729): "503683295ebeb48587f42748181bf6c5815f6a8c79b8dfdd4b669c931109be81",
    ("h3k27me3", 130363): "56a2d331b6ed27e49f103dce228dda3c26a80ff296e91c9c740e52b9c36841c6",
    ("h3k27me3", 155921): "995dc295d163a703b80bb3e451b8fc3bd7af01fb6ce041099fc5ce3d9cc82560",
}
THRESHOLDS = {
    "enhancer": ("0.74848511815071117", "0x1.7f39710000001p-1"),
    "h3k27me3": ("0.76960810025533055", "0x1.8a0a12aaaaaacp-1"),
}


def utcnow():
    return dt.datetime.now(dt.timezone.utc).isoformat()


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    with Path(path).open() as stream:
        return json.load(stream)


def write_json(path, data):
    with Path(path).open("x") as stream:
        json.dump(data, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def record(base, path):
    path = Path(path)
    return {"path": str(path.relative_to(base)), "bytes": path.stat().st_size,
            "sha256": sha256(path)}


def verify_record(base, item, *, allow_managed_symlink=False):
    base = Path(base).resolve()
    relative = Path(item["path"])
    if relative.is_absolute() or ".." in relative.parts:
        raise RuntimeError(f"Unscoped manifest path: {relative}")
    path = base / relative
    if not allow_managed_symlink and not path.resolve().is_relative_to(base):
        raise RuntimeError(f"Manifest symlink escapes scope: {relative}")
    if path.stat().st_size != int(item["bytes"]) or sha256(path) != item["sha256"]:
        raise RuntimeError(f"Frozen byte contract failed: {relative}")
    return path


def verified_records(base, items):
    if not isinstance(items, list) or not items:
        raise RuntimeError("Gate requires a nonempty file-record list")
    result = {}
    for item in items:
        if item["path"] in result:
            raise RuntimeError("Duplicate file in gate")
        result[item["path"]] = verify_record(base, item)
    return result


def verify_prospective(stage):
    path = stage / "provenance/prospective_freeze.json"
    if sha256(path) != PROSPECTIVE_SHA:
        raise RuntimeError("Original prospective freeze changed")
    gate = read_json(path)
    if (gate.get("status") != "PASS" or gate.get("stage") != STAGE_NAME
            or gate.get("before_any_V2_external_prediction") is not True
            or gate.get("before_any_benchmark_outcome_or_label_opening") is not True):
        raise RuntimeError("Prospective benchmark freeze is not PASS")
    verified_records(stage, gate["files"])
    path = stage / "specification/evaluation_specification.json"
    if sha256(path) != SPECIFICATION_SHA:
        raise RuntimeError("Prospective scientific specification changed")
    specification = read_json(path)
    if (specification["stage"] != STAGE_NAME or specification["seeds"] != list(SEEDS)
            or specification["contexts"] != list(CONTEXTS)
            or specification["checkpoint_configuration"] != "V2-C"
            or specification["phase_I_batch_size"] != 32
            or specification["phase_II_batch_size"] != 256
            or specification["sequence_length"] != 2001
            or specification["invariance_atol"] != ATOL
            or specification["invariance_rtol"] != RTOL):
        raise RuntimeError("Unrecognized frozen model arithmetic contract")
    for key in ("full_rules", "preopen_sources", "user_authorization"):
        verify_record(stage, specification[key])
    return specification, gate


def verify_inputs(stage, prospective):
    annotation_path = stage / "provenance/annotation_freeze.json"
    annotation = read_json(annotation_path)
    if (annotation.get("status") != "PASS"
            or annotation.get("before_any_V2_external_prediction") is not True):
        raise RuntimeError("Annotation map was not frozen before V2 predictions")
    annotation_files = verified_records(stage, annotation["files"])
    required = {"inputs/input_manifest.json", "inputs/scorable_variants.tsv",
                "inputs/allele_sequences.npy", "provenance/prospective_freeze.json"}
    if not required <= set(annotation_files):
        raise RuntimeError("Annotation freeze must bind specification and exact scoring inputs")
    for gate in (annotation, prospective):
        when = dt.datetime.fromisoformat(gate["created_utc"])
        if when.tzinfo is None or when > dt.datetime.now(dt.timezone.utc):
            raise RuntimeError("Invalid or future freeze timestamp")
    if (dt.datetime.fromisoformat(annotation["created_utc"])
            < dt.datetime.fromisoformat(prospective["created_utc"])):
        raise RuntimeError("Annotation freeze predates prospective plan")
    manifest = read_json(stage / "inputs/input_manifest.json")
    if manifest.get("status") != "PASS":
        raise RuntimeError("Sequence preparation was not PASS")
    input_files = verified_records(stage, manifest["files"])
    if not {"inputs/scorable_variants.tsv", "inputs/allele_sequences.npy"} <= set(input_files):
        raise RuntimeError("Input manifest lacks scoring table or sequence array")
    # Parse identity-only table, never experimental labels or outcomes.
    with (stage / "inputs/scorable_variants.tsv").open(newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if not rows or any(not {"variant_id", "chrom", "pos1", "ref", "alt"} <= set(r) for r in rows):
        raise RuntimeError("Missing identity-only scoring schema or empty population")
    if len({r["variant_id"] for r in rows}) != len(rows):
        raise RuntimeError("Duplicate exact variant identifiers")
    keys = [(r["chrom"], int(r["pos1"]), r["ref"], r["alt"]) for r in rows]
    if len(set(keys)) != len(rows):
        raise RuntimeError("Exact identities were not deduplicated")
    if "n_variants" in manifest and manifest["n_variants"] != len(rows):
        raise RuntimeError("Input population differs from manifest")
    sequences = np.load(stage / "inputs/allele_sequences.npy", allow_pickle=False)
    if sequences.dtype != np.uint8 or sequences.shape != (len(rows), 2, 2001):
        raise RuntimeError("Allele sequence array must be uint8[Nvariants,REF_then_ALT,2001]")
    if not np.isin(sequences, np.frombuffer(b"ACGTN", dtype=np.uint8)).all():
        raise RuntimeError("Unsupported nucleotide in fixed input array")
    sequence_records = []
    for i, row in enumerate(rows):
        exact_key = ":".join(row[k] for k in ("chrom", "pos1", "ref", "alt"))
        allowed_chromosomes = {str(c) for c in range(1, 23)} | {"X", "Y"}
        if (row["variant_id"] != exact_key or row["chrom"] not in allowed_chromosomes
                or int(row["pos1"]) <= 1000):
            raise RuntimeError("Invalid exact variant identity")
        for allele_index, name in enumerate(("ref", "alt")):
            allele = row[name]
            if not 1 <= len(allele) <= 1001 or not set(allele) <= set("ACGT"):
                raise RuntimeError("Ambiguous or unsupported frozen allele")
            sequence = sequences[i, allele_index].tobytes().decode("ascii")
            if sequence[1000:1000 + len(allele)] != allele:
                raise RuntimeError("Frozen allele is not anchored at position 1000")
            rc = sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]
            if rc.translate(str.maketrans("ACGTN", "TGCAN"))[::-1] != sequence:
                raise RuntimeError("Nucleotide RC involution failed")
            if rc[1001 - len(allele):1001] != allele.translate(str.maketrans("ACGT", "TGCA"))[::-1]:
                raise RuntimeError("RC allele-span geometry failed")
            sequence_records.append({"variant_id": row["variant_id"], "allele": name.upper(),
                                     "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
                                     "nucleotide_RC_sha256": hashlib.sha256(rc.encode()).hexdigest()})
        ref_seq, alt_seq = (sequences[i, a].tobytes().decode("ascii") for a in (0, 1))
        if ref_seq[:1000] != alt_seq[:1000]:
            raise RuntimeError("REF/ALT left flanks differ")
        ref_right, alt_right = ref_seq[1000 + len(row["ref"]):], alt_seq[1000 + len(row["alt"]):]
        shared = min(len(ref_right), len(alt_right))
        if ref_right[:shared] != alt_right[:shared]:
            raise RuntimeError("Indel downstream context differs from fixed anchored geometry")
    return rows, sequences, sequence_records


def verify_models(repo, specification):
    checkpoints = {}
    for item in specification["checkpoints"]:
        key = item["context"], int(item["seed"])
        expected_path = f"{TRAINING}/runs/V2-C_{key[0]}_seed{key[1]}/attempt-001/selected_checkpoint.keras"
        if (key not in CHECKPOINT_SHA or key in checkpoints or item["path"] != expected_path
                or item["configuration"] != "V2-C" or item["bytes"] != 350064314
                or item["sha256"] != CHECKPOINT_SHA[key]):
            raise RuntimeError("Unexpected selected checkpoint identity or hash")
        checkpoints[key] = verify_record(repo, item)
    if set(checkpoints) != set(CHECKPOINT_SHA):
        raise RuntimeError("All six frozen V2-C archives are required")
    # Only the model/runtime dependencies actually used here are opened. No prior
    # stage evaluation entrypoint, candidate data or benchmark outcome is read.
    dependencies = {item["path"]: item for item in specification["dependencies"]}
    required = {
        f"{TRAINING}/scripts/phase_two_contract.py": HELPER_SHA,
        f"{TRAINING}/scripts/extract_phase_one.py": ARCHITECTURE_SHA,
        f"{TRAINING}/results/chr7_calibration/C_region_thresholds.json": THRESHOLD_SHA,
        f"{TRAINING}/provenance/checkpoint_freeze.json": "d2330e4881c477cb7d8e5b2b621940674761133b779c89d0db1ad88635eee379",
        "diseases/COPD/04_modeling/trednet/model_phase_I/phase_one_weights.h5": PHASE_ONE_SHA,
        "diseases/COPD/07_gap_closure/internal-test-1.0/scripts/runtime.sh": "bd4ebbfcd3443dece86d4d535f313f3676ead7fb5dd8319e7ce0ec2f843eb0fc",
    }
    verified_dependencies = []
    for name, digest in required.items():
        item = dependencies[name]
        if item["sha256"] != digest:
            raise RuntimeError("Frozen numerical dependency differs from implementation contract")
        verify_record(repo, item, allow_managed_symlink=name.endswith("phase_one_weights.h5"))
        verified_dependencies.append(item)
    calibration = read_json(repo / f"{TRAINING}/results/chr7_calibration/C_region_thresholds.json")
    thresholds = {}
    for context, (decimal, hexadecimal) in THRESHOLDS.items():
        specified = specification["thresholds"][context]
        original = calibration["models"][context]
        value = float.fromhex(hexadecimal)
        if (format(value, ".17g") != decimal or specified != {
                "decimal17g": decimal, "float64_hex": hexadecimal, "operator": ">="}
                or original["threshold_hex"] != hexadecimal
                or original["threshold_decimal_17g"] != decimal
                or original["call_rule"] != "score>=threshold"
                or original["allele_delta_threshold"] is not None):
            raise RuntimeError("Exact frozen region threshold or operator changed")
        thresholds[context] = value
    return checkpoints, thresholds, verified_dependencies


def import_frozen_numerical_helpers(repo):
    directory = repo / TRAINING / "scripts"
    modules = []
    for name, filename, digest in (
            ("phase_two_contract", "phase_two_contract.py", HELPER_SHA),
            ("frozen_external_phase_one_architecture", "extract_phase_one.py", ARCHITECTURE_SHA)):
        path = directory / filename
        if sha256(path) != digest:
            raise RuntimeError("Immutable helper changed before import")
        if name in sys.modules:
            raise RuntimeError("Refuse unverified pre-imported numerical helper")
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
        modules.append(module)
    return modules


def array_sha(array):
    return hashlib.sha256(np.ascontiguousarray(array).tobytes()).hexdigest()


def save_array(path, array):
    with Path(path).open("xb") as stream:
        np.save(stream, array, allow_pickle=False)


def probabilities(value):
    value = np.asarray(value)
    if not np.isfinite(value).all() or np.any((value < 0) | (value > 1)):
        raise RuntimeError("Nonfinite or out-of-range probability/phase-I activation")
    return value


def invariance(original, independent_rc):
    a, b = (np.asarray(x, dtype=np.float64) for x in (original, independent_rc))
    if a.shape != b.shape or not a.size:
        raise RuntimeError("Invariance audit requires aligned nonempty arrays")
    probabilities(a)
    probabilities(b)
    residual = np.abs(a - b)
    bound = ATOL + RTOL * np.abs(b)
    return {"status": "PASS" if np.all(residual <= bound) else "FAIL",
            "n_probabilities": int(a.size), "n_failed": int(np.sum(residual > bound)),
            "atol": ATOL, "rtol": RTOL,
            "max_absolute_difference": float(residual.max()),
            "max_tolerance_normalized_difference": float((residual / bound).max()),
            "original_array_sha256": array_sha(a), "independent_RC_array_sha256": array_sha(b)}


def phase_one_passes(tf, keras, build_phase_one, weights, sequences, output):
    started = time.monotonic()
    chars = sequences.reshape(-1, 2001)
    model = build_phase_one(keras)
    model.load_weights(weights)
    model.trainable = False
    if model.input_shape != (None, 2001, 4) or model.output_shape != (None, 4560):
        raise RuntimeError("Phase-I architecture geometry differs from frozen contract")
    if any(not np.isfinite(variable.numpy()).all() for variable in model.weights):
        raise RuntimeError("Nonfinite frozen phase-I weight")
    write_json(output / "phase_I_architecture.json", json.loads(model.to_json()))
    infer = tf.function(lambda batch: model(batch, training=False),
                        input_signature=[tf.TensorSpec((None, 2001, 4), tf.float32)], jit_compile=False)
    features = np.empty((4, len(chars), 4560), dtype=np.float32)
    alphabet = np.frombuffer(b"ACGT", dtype=np.uint8)
    complement = np.arange(256, dtype=np.uint8)
    complement[np.frombuffer(b"ACGTN", dtype=np.uint8)] = np.frombuffer(b"TGCAN", dtype=np.uint8)
    # Each loop performs independent model calls. Audit phase I repeats nucleotide
    # conversion and network evaluation in RC-first, forward-second order.
    for audit in (False, True):
        for start in range(0, len(chars), 32):
            batch = chars[start:start + 32]
            forward = (batch[:, :, None] == alphabet[None, None, :]).astype(np.float32)
            rc_chars = complement[batch[:, ::-1]]
            reverse = (rc_chars[:, :, None] == alphabet[None, None, :]).astype(np.float32)
            if not np.array_equal(reverse, forward[:, ::-1, ::-1]):
                raise RuntimeError("Nucleotide/one-hot RC orientation mismatch")
            passes = ((2, reverse), (3, forward)) if audit else ((0, forward), (1, reverse))
            for index, data in passes:
                value = np.asarray(infer(np.ascontiguousarray(data)), dtype=np.float32)
                if value.shape != (len(batch), 4560):
                    raise RuntimeError("Invalid phase-I feature geometry")
                probabilities(value)
                features[index, start:start + len(batch)] = value
    names = ("phase_I_forward.npy", "phase_I_nucleotide_RC.npy", "phase_I_independent_RC_first.npy",
             "phase_I_independent_forward_second.npy")
    for name, values in zip(names, features):
        save_array(output / name, values)
    provenance = {
        "status": "PASS", "shape_per_orientation": [len(chars), 4560], "dtype": "float32",
        "flattening": "variant table row order, REF then ALT for each variant",
        "feature_names_in_orientation_order": list(names), "representation": "final4560sigmoid",
        "phase_I_weights_sha256": PHASE_ONE_SHA, "architecture_helper_sha256": ARCHITECTURE_SHA,
        "phase_I_batch_size": 32, "independent_network_calls_per_sequence": 4,
        "orientation_sequence_evaluations": 4 * len(chars), "training": False, "optimizer": None,
        "audit_recomputes_phase_I_features": True, "representation_axis_reversal": False,
        "repeated_forward_max_abs_difference": float(np.max(np.abs(features[0] - features[3]))),
        "repeated_RC_max_abs_difference": float(np.max(np.abs(features[1] - features[2]))),
        "array_data_sha256": [array_sha(values) for values in features],
        "wall_seconds": time.monotonic() - started,
    }
    if sha256(weights) != PHASE_ONE_SHA:
        raise RuntimeError("Original phase-I weights changed during evaluation")
    del infer, model
    keras.backend.clear_session()
    return features, provenance


def phase_two_predict(infer, features):
    result = np.empty(len(features), dtype=np.float64)
    for start in range(0, len(features), 256):
        batch = features[start:start + 256]
        value = np.asarray(infer(batch[..., None])).reshape(-1)
        if value.dtype != np.float32 or value.shape != (len(batch),):
            raise RuntimeError("Phase-II output must be a float32 probability per allele")
        probabilities(value)
        result[start:start + len(batch)] = value.astype(np.float64)
    return result


def numerical_fields(ref_score, alt_score, threshold):
    delta = float(alt_score) - float(ref_score)
    values = {"ref_score": float(ref_score), "alt_score": float(alt_score), "delta": delta,
              "abs_delta": abs(delta), "max_allele_score": max(float(ref_score), float(alt_score)),
              "threshold": threshold}
    result = {}
    for name, value in values.items():
        result[name] = format(value, ".17g")
        result[name + "_hex"] = value.hex()
    result.update({"delta_sign": int(np.sign(delta)), "ref_region_call": int(ref_score >= threshold),
                   "alt_region_call": int(alt_score >= threshold),
                   "either_allele_region_call": int(max(ref_score, alt_score) >= threshold),
                   "region_call_operator": ">="})
    return result


def write_table(path, rows):
    if not rows:
        raise RuntimeError("Refuse empty prediction table")
    with Path(path).open("x", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def resources(started):
    usage = resource.getrusage(resource.RUSAGE_SELF)
    return {"finished_utc": utcnow(), "wall_seconds": time.monotonic() - started,
            "user_cpu_seconds": usage.ru_utime, "system_cpu_seconds": usage.ru_stime,
            "peak_RSS_KiB": usage.ru_maxrss}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stage", required=True, type=Path)
    parser.add_argument("--repo", required=True, type=Path)
    parser.add_argument("--output-relative", default="predictions",
                        help="New stage directory only; retries require a separately logged justification")
    args = parser.parse_args()
    stage, repo = args.stage.resolve(), args.repo.resolve()
    if stage != repo / "diseases/COPD/07_gap_closure" / STAGE_NAME:
        raise RuntimeError("Inference scope must be the exact new external benchmark stage")
    output = (stage / args.output_relative).resolve()
    if output == stage or not output.is_relative_to(stage):
        raise RuntimeError("Output directory must be strictly inside the new stage")
    output.mkdir(parents=False, exist_ok=False)
    started = time.monotonic()
    own_records = [record(stage, Path(__file__).resolve()), record(stage, Path(__file__).with_name("runtime.sh").resolve())]
    write_json(output / "started.json", {
        "started_utc": utcnow(), "argv": sys.argv, "stage": STAGE_NAME,
        "scope": "Only exact sequence-eligible frozen external benchmark variants; six unchanged V2-C models",
        "script_files": own_records, "experimental_outcome_tables_parsed": False,
        "training_or_calibration": False, "candidate_universe_access": False})
    try:
        specification, prospective = verify_prospective(stage)
        rows, sequences, sequence_records = verify_inputs(stage, prospective)
        checkpoints, thresholds, dependencies = verify_models(repo, specification)
        write_json(output / "verified_before_inference.json", {
            "status": "PASS", "verified_utc": utcnow(), "n_variants": len(rows),
            "six_checkpoint_records": specification["checkpoints"], "numerical_dependencies": dependencies,
            "prospective_freeze": record(stage, stage / "provenance/prospective_freeze.json"),
            "annotation_freeze": record(stage, stage / "provenance/annotation_freeze.json"),
            "input_manifest": record(stage, stage / "inputs/input_manifest.json"),
            "script_files": own_records, "sequence_records": sequence_records,
            "annotation_files_verified_as_bytes_only": True})
        helper, architecture = import_frozen_numerical_helpers(repo)
        if helper.SEEDS != SEEDS or helper.MODELS != CONTEXTS or helper.ATOL != ATOL or helper.RTOL != RTOL:
            raise RuntimeError("Imported frozen helper constants differ")
        tf, keras = helper.initialize_runtime(SEEDS[0])
        write_json(output / "environment.json", helper.runtime_environment(tf, keras))
        weights = repo / "diseases/COPD/04_modeling/trednet/model_phase_I/phase_one_weights.h5"
        features, phase_one_provenance = phase_one_passes(
            tf, keras, architecture.build_phase_one, weights, sequences, output)
        phase_one_provenance["files"] = [record(stage, output / name)
                                          for name in phase_one_provenance["feature_names_in_orientation_order"]]
        write_json(output / "phase_I_provenance.json", phase_one_provenance)
        # axes: context, frozen seed, variant row, REF/ALT, original F/RC and
        # independent RC-first/F-second audit orientation. Stored as float64.
        raw = np.empty((2, 3, len(rows), 2, 4), dtype=np.float64)
        symmetric = np.empty((2, 3, len(rows), 2), dtype=np.float64)
        symmetric_rc = np.empty_like(symmetric)
        ensemble = np.empty((2, len(rows), 2), dtype=np.float64)
        ensemble_rc = np.empty_like(ensemble)
        audits, seed_rows, ensemble_rows = [], [], []
        for context_index, context in enumerate(CONTEXTS):
            for seed_index, seed in enumerate(SEEDS):
                keras.backend.clear_session()
                model = keras.models.load_model(checkpoints[context, seed], compile=False)
                if model.input_shape != (None, 4560, 1) or model.output_shape != (None, 1):
                    raise RuntimeError("Selected checkpoint has unexpected inference geometry")
                infer = helper.inference_function(tf, model)
                passes = [phase_two_predict(infer, features[o]) for o in range(4)]
                q = helper.symmetric_probability(passes[0], passes[1]).reshape(len(rows), 2)
                q_rc = helper.symmetric_probability(passes[2], passes[3]).reshape(len(rows), 2)
                raw[context_index, seed_index] = np.stack(passes, axis=-1).reshape(len(rows), 2, 4)
                symmetric[context_index, seed_index] = q
                symmetric_rc[context_index, seed_index] = q_rc
                audit = {"context": context, "unit": "seed", "seed": seed,
                         "phase_I_and_phase_II_independently_repeated": True,
                         "checkpoint_sha256": CHECKPOINT_SHA[context, seed], **invariance(q, q_rc)}
                audits.append(audit)
                write_json(output / f"{context}_seed{seed}_RC_invariance.json", audit)
                if audit["status"] != "PASS":
                    raise RuntimeError(f"Real-network RC invariance failed for {context}/{seed}")
                for index, identity in enumerate(rows):
                    row = {k: identity[k] for k in ("variant_id", "chrom", "pos1", "ref", "alt")}
                    row.update({"context": context, "seed": seed,
                                "seed_calls_are_descriptive_against_ensemble_threshold": True})
                    row.update(numerical_fields(q[index, 0], q[index, 1], thresholds[context]))
                    for a, allele in enumerate(("ref", "alt")):
                        for o, orientation in enumerate(("forward", "nucleotide_RC", "independent_RC_first", "independent_forward_second")):
                            value = float(raw[context_index, seed_index, index, a, o])
                            name = f"{allele}_p_{orientation}"
                            row[name], row[name + "_hex"] = format(value, ".17g"), value.hex()
                        value = float(q_rc[index, a])
                        row[f"{allele}_q_independent_RC"] = format(value, ".17g")
                        row[f"{allele}_q_independent_RC_hex"] = value.hex()
                    seed_rows.append(row)
                print(f"PASS real-network RC: V2-C {context} seed {seed}; {len(rows)} frozen variants", flush=True)
                del infer, model
            ensemble[context_index] = helper.seed_ensemble(list(symmetric[context_index]))
            ensemble_rc[context_index] = helper.seed_ensemble(list(symmetric_rc[context_index]))
            audit = {"context": context, "unit": "ensemble", "seed": "all3",
                     "ordered_seeds": list(SEEDS), "reduction_dtype": "float64",
                     "phase_I_and_phase_II_independently_repeated": True,
                     **invariance(ensemble[context_index], ensemble_rc[context_index])}
            audits.append(audit)
            write_json(output / f"{context}_ensemble_RC_invariance.json", audit)
            if audit["status"] != "PASS":
                raise RuntimeError(f"Real-network ensemble RC invariance failed for {context}")
            for index, identity in enumerate(rows):
                row = {k: identity[k] for k in ("variant_id", "chrom", "pos1", "ref", "alt")}
                row["context"] = context
                row.update(numerical_fields(*ensemble[context_index, index], thresholds[context]))
                for a, allele in enumerate(("ref", "alt")):
                    value = float(ensemble_rc[context_index, index, a])
                    row[f"{allele}_ensemble_independent_RC"] = format(value, ".17g")
                    row[f"{allele}_ensemble_independent_RC_hex"] = value.hex()
                ensemble_rows.append(row)
        for name, array in (("raw_orientation_probabilities.npy", raw),
                            ("seed_symmetric_probabilities.npy", symmetric),
                            ("seed_independent_RC_probabilities.npy", symmetric_rc),
                            ("ensemble_probabilities.npy", ensemble),
                            ("ensemble_independent_RC_probabilities.npy", ensemble_rc)):
            probabilities(array)
            save_array(output / name, array)
        write_table(output / "seed_allele_scores.tsv", seed_rows)
        write_table(output / "variant_context_scores.tsv", ensemble_rows)
        write_json(output / "real_network_invariance.json", {
            "status": "PASS", "n_seed_audits": 6, "n_ensemble_audits": 2,
            "real_network": True, "phase_I_and_phase_II_independently_repeated": True,
            "same_frozen_variant_allele_order_and_batching": True,
            "maximum_absolute_residual": max(row["max_absolute_difference"] for row in audits),
            "atol": ATOL, "rtol": RTOL, "audit_rows": audits})
        # Reverify immutable archives, numerical dependencies, annotation and input
        # contracts after execution. No model save or overwrite occurs anywhere.
        verify_models(repo, specification)
        verify_prospective(stage)
        verify_inputs(stage, prospective)
        for item in own_records:
            verify_record(stage, item)
        files = [record(stage, path) for path in sorted(output.iterdir()) if path.is_file()]
        write_json(output / "inference_manifest.json", {
            "status": "PASS", "completed": True, "stage": STAGE_NAME,
            "n_variants": len(rows), "n_allele_sequences": 2 * len(rows), "n_original_checkpoints": 6,
            "n_context_variant_scores": len(ensemble_rows), "n_seed_variant_scores": len(seed_rows),
            "contexts": list(CONTEXTS), "ordered_seeds": list(SEEDS), "allele_order": ["REF", "ALT"],
            "raw_array_axes": ["context", "seed", "variant_table_row", "allele", "orientation"],
            "raw_orientation_order": ["forward", "nucleotide_RC", "independent_RC_first", "independent_forward_second"],
            "symmetric_array_axes": ["context", "seed", "variant_table_row", "allele"],
            "ensemble_array_axes": ["context", "variant_table_row", "allele"],
            "probability_array_dtype": "float64; raw network float32 values promoted before arithmetic",
            "thresholds": specification["thresholds"], "allele_delta_threshold": None,
            "region_thresholds_are_not_variant_FPR_or_allele_effect_thresholds": True,
            "all_original_archives_unchanged_after_inference": True,
            "experimental_outcome_tables_parsed": False, "candidate_universe_access": False,
            "prior_training_or_internal_test_entrypoint_execution": False,
            "training_or_recalibration": False, "all_six_seeds_contexts_complete": True,
            "files": files, "self_reference_exclusions": ["inference_manifest.json", "completed.json"],
            "resources": resources(started)})
        write_json(output / "completed.json", {
            "status": "COMPLETED", "invariance_status": "PASS",
            "inference_manifest": record(stage, output / "inference_manifest.json"),
            "resources": resources(started)})
        print(json.dumps({"status": "PASS", "n_variants": len(rows), "n_checkpoints": 6,
                          "output": str(output)}), flush=True)
    except BaseException as exc:
        write_json(output / "failure.json", {"status": "FAILED", "error": str(exc),
                   "exception_type": type(exc).__name__, "traceback": traceback.format_exc(),
                   "downstream_evaluation_authorized": False, "resources": resources(started)})
        raise


if __name__ == "__main__":
    main()
