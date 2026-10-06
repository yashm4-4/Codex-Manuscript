#!/usr/bin/env python3
"""CPU-only synthetic contract tests; never load a model or biological data.

These helpers specify the V2 inference arithmetic. Passing these tests does NOT
validate the behavior of any trained network, model loader, or sequence extractor.
Actual trained-network validation remains a mandatory later release gate.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
from datetime import datetime, timezone
from fractions import Fraction
from pathlib import Path
from typing import Callable, Mapping

import numpy as np


TRAINING_SEEDS = (104729, 130363, 155921)
BOOTSTRAP_SEED = 314159
FIXTURE_SEED = 271828
ATOL = 1e-6
RTOL = 1e-6
ALPHABET = "ACGTN"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")
Predictor = Callable[[np.ndarray], np.ndarray]


def reverse_complement_string(sequence: str) -> str:
    sequence = sequence.upper()
    if set(sequence) - set(ALPHABET):
        raise ValueError("Only A, C, G, T and N are allowed in this contract")
    return sequence.translate(COMPLEMENT)[::-1]


def one_hot(sequence: str) -> np.ndarray:
    sequence = sequence.upper()
    if set(sequence) - set(ALPHABET):
        raise ValueError("Only A, C, G, T and N are allowed in this contract")
    letters = np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)
    return (letters[:, None] == np.frombuffer(b"ACGT", dtype=np.uint8)).astype(np.uint8)


def reverse_complement_one_hot(batch: np.ndarray) -> np.ndarray:
    batch = np.asarray(batch)
    if batch.ndim not in (2, 3) or batch.shape[-1] != 4:
        raise ValueError("Expected (length, 4) or (batch, length, 4) in A,C,G,T order")
    if not np.isfinite(batch).all():
        raise ValueError("Nonfinite sequence encoding")
    return batch[..., ::-1, ::-1].copy()


def _probabilities(predictor: Predictor, batch: np.ndarray) -> np.ndarray:
    # Predictors must already be in inference mode. This helper cannot enforce
    # dropout/BatchNorm state or determine whether a callback loads a real model.
    output = np.asarray(predictor(batch), dtype=np.float64)
    if output.shape == (len(batch), 1):
        output = output[:, 0]
    if output.shape != (len(batch),):
        raise ValueError("Predictor must return one scalar probability per sequence")
    if not np.isfinite(output).all() or np.any((output < 0) | (output > 1)):
        raise ValueError("Probabilities must be finite and between zero and one")
    return output


def symmetric_predictions(batch: np.ndarray, predictor: Predictor) -> dict[str, np.ndarray]:
    batch = np.asarray(batch)
    if batch.ndim != 3 or batch.shape[-1] != 4:
        raise ValueError("Expected (batch, length, 4)")
    forward = _probabilities(predictor, batch)
    reverse = _probabilities(predictor, reverse_complement_one_hot(batch))
    return {
        "forward": forward,
        "reverse_complement": reverse,
        "symmetric": (forward + reverse) / np.float64(2.0),
    }


def three_seed_predictions(
    batch: np.ndarray, predictors: Mapping[int, Predictor]
) -> dict[str, object]:
    if set(predictors) != set(TRAINING_SEEDS):
        raise ValueError("All and only the three prespecified training seeds are required")
    by_seed = {str(seed): symmetric_predictions(batch, predictors[seed]) for seed in TRAINING_SEEDS}
    symmetric = np.stack([by_seed[str(seed)]["symmetric"] for seed in TRAINING_SEEDS])
    return {
        "seed_order": TRAINING_SEEDS,
        "by_seed": by_seed,
        "ensemble_symmetric": symmetric.mean(axis=0, dtype=np.float64),
        "seed_sd_sample": symmetric.std(axis=0, ddof=1, dtype=np.float64),
        "seed_min": symmetric.min(axis=0),
        "seed_max": symmetric.max(axis=0),
    }


def validation_control_cutoff(
    control_scores: np.ndarray, target_fpr: float = 0.05, minimum_controls: int = 200
) -> dict[str, object]:
    """Tie-conservative empirical FPR cutoff, for calibration controls ONLY.

    Calls use score >= cutoff. For n controls and k=floor(alpha*n), use the
    next representable float64 above the (k+1)-th descending control score.
    No labels/test scores enter this routine. The empirical bound is not a
    guarantee for future/test population FPR. A cutoff >1 is valid and means
    no probability-valued example is called positive at this operating point.
    """
    scores = np.asarray(control_scores, dtype=np.float64)
    if scores.ndim != 1 or len(scores) < minimum_controls:
        raise ValueError("Insufficient one-dimensional calibration controls")
    if not np.isfinite(scores).all() or np.any((scores < 0) | (scores > 1)):
        raise ValueError("Calibration probabilities must be finite and in [0,1]")
    alpha = Fraction(str(target_fpr))
    if not 0 < alpha < 1:
        raise ValueError("Target FPR must be strictly between zero and one")
    allowed = (len(scores) * alpha.numerator) // alpha.denominator
    boundary = np.sort(scores)[::-1][allowed]
    cutoff = np.nextafter(np.float64(boundary), np.float64(np.inf))
    called = int(np.count_nonzero(scores >= cutoff))
    if called > allowed:
        raise AssertionError("Tie-conservative cutoff contract violated")
    return {
        "cutoff": float(cutoff),
        "cutoff_hex": float(cutoff).hex(),
        "cutoff_decimal_17_significant": format(float(cutoff), ".17g"),
        "target_fpr": float(alpha),
        "n_calibration_controls": len(scores),
        "allowed_false_positives": allowed,
        "called_calibration_controls": called,
        "empirical_calibration_fpr": called / len(scores),
        "comparison": "score >= cutoff",
        "population_fpr_guaranteed": False,
    }


def _mock_predictor(seed: int) -> Predictor:
    # Explicitly orientation-sensitive, deterministic, non-neural fixtures.
    offset = TRAINING_SEEDS.index(seed)

    def predict(batch: np.ndarray) -> np.ndarray:
        x = np.asarray(batch, dtype=np.float64)
        logit = (
            (0.5 + 0.15 * offset) * (x[:, 0, 0] - x[:, -1, 3])
            + 0.4 * (x[:, 17, 1] - x[:, -23, 2])
            + 0.3 * (x[:, 1000, 0] - x[:, 1000, 1])
            + 0.1 * offset
        )
        return 1.0 / (1.0 + np.exp(-logit))

    return predict


def run_synthetic_qc() -> dict[str, object]:
    checks: list[dict[str, object]] = []

    def check(name: str, passed: bool, **detail: object) -> None:
        checks.append({"check": name, "status": "PASS" if passed else "FAIL", **detail})

    def rejects(name: str, callback: Callable[[], object]) -> None:
        try:
            callback()
        except ValueError:
            check(name, True)
        else:
            check(name, False)

    rng = np.random.default_rng(FIXTURE_SEED)
    codes = rng.integers(0, 4, size=(1001, 2001), dtype=np.uint8)
    codes[rng.random(codes.shape) < 0.01] = 4
    alphabet = np.frombuffer(b"ACGTN", dtype=np.uint8)
    references = [alphabet[row].tobytes().decode("ascii") for row in codes]
    half = references[0][:1000]
    references[:5] = ["A" * 2001, "T" * 2001, "N" * 2001, "C" * 2001,
                      half + "N" + reverse_complement_string(half)]
    alternates: list[str] = []
    positions: list[int] = []
    for index, reference in enumerate(references):
        position = (index * 193) % 2001
        alternate = next(base for base in "ACGT" if base != reference[position])
        alternates.append(reference[:position] + alternate + reference[position + 1:])
        positions.append(position)
    check("1001_independent_synthetic_REF_ALT_pairs", len(references) == 1001 and len(alternates) == 1001)
    check("string_RC_involution_REF_and_ALT", all(
        reverse_complement_string(reverse_complement_string(s)) == s for s in references + alternates
    ))
    check("ALT_is_complemented_not_swapped_with_REF", all(
        reverse_complement_string(alt)[2000 - pos] == alt[pos].translate(COMPLEMENT)
        and reverse_complement_string(ref)[2000 - pos] == ref[pos].translate(COMPLEMENT)
        for ref, alt, pos in zip(references, alternates, positions)
    ))
    x_ref = np.stack([one_hot(s) for s in references])
    x_alt = np.stack([one_hot(s) for s in alternates])
    for label, strings, batch in (("REF", references, x_ref), ("ALT", alternates, x_alt)):
        rc = reverse_complement_one_hot(batch)
        check(f"onehot_string_RC_agreement_{label}", np.array_equal(rc, np.stack([
            one_hot(reverse_complement_string(s)) for s in strings
        ])))
        check(f"onehot_RC_involution_{label}", np.array_equal(reverse_complement_one_hot(rc), batch))
    check("N_encoded_as_zero_and_preserved", not one_hot("N").any() and not reverse_complement_one_hot(one_hot("N")).any())
    predictors = {seed: _mock_predictor(seed) for seed in TRAINING_SEEDS}
    ref = three_seed_predictions(x_ref, predictors)
    alt = three_seed_predictions(x_alt, predictors)
    ref_rc = three_seed_predictions(reverse_complement_one_hot(x_ref), predictors)
    alt_rc = three_seed_predictions(reverse_complement_one_hot(x_alt), predictors)
    for label, original, transformed in (("REF", ref, ref_rc), ("ALT", alt, alt_rc)):
        error = float(np.max(np.abs(original["ensemble_symmetric"] - transformed["ensemble_symmetric"])))
        check(f"ensemble_invariance_{label}", np.allclose(original["ensemble_symmetric"], transformed["ensemble_symmetric"], atol=ATOL, rtol=RTOL), maximum_absolute_error=error)
        for seed in TRAINING_SEEDS:
            first, second = original["by_seed"][str(seed)], transformed["by_seed"][str(seed)]
            check(f"per_seed_invariance_{label}_{seed}", np.allclose(first["symmetric"], second["symmetric"], atol=ATOL, rtol=RTOL))
            check(f"orientation_columns_swap_only_{label}_{seed}", np.array_equal(first["forward"], second["reverse_complement"]) and np.array_equal(first["reverse_complement"], second["forward"]))
    raw_gap = max(float(np.max(np.abs(ref["by_seed"][str(seed)]["forward"] - ref["by_seed"][str(seed)]["reverse_complement"]))) for seed in TRAINING_SEEDS)
    check("mock_predictors_are_deliberately_orientation_sensitive", raw_gap > 0.1, maximum_forward_RC_gap=raw_gap)
    delta = alt["ensemble_symmetric"] - ref["ensemble_symmetric"]
    delta_rc = alt_rc["ensemble_symmetric"] - ref_rc["ensemble_symmetric"]
    check("ALT_minus_REF_delta_invariant", np.allclose(delta, delta_rc, atol=ATOL, rtol=RTOL))
    check("max_REF_ALT_region_score_invariant", np.allclose(np.maximum(ref["ensemble_symmetric"], alt["ensemble_symmetric"]), np.maximum(ref_rc["ensemble_symmetric"], alt_rc["ensemble_symmetric"]), atol=ATOL, rtol=RTOL))
    mean_deltas = np.mean(np.stack([alt["by_seed"][str(seed)]["symmetric"] - ref["by_seed"][str(seed)]["symmetric"] for seed in TRAINING_SEEDS]), axis=0, dtype=np.float64)
    check("ensemble_delta_equals_equal_weight_mean_seed_deltas", np.allclose(delta, mean_deltas, atol=1e-15, rtol=0))
    seed_stack = np.stack([ref["by_seed"][str(seed)]["symmetric"] for seed in TRAINING_SEEDS])
    check("seed_SD_is_sample_SD_ddof1", np.array_equal(ref["seed_sd_sample"], np.std(seed_stack, axis=0, ddof=1)))
    reversed_seeds = {seed: predictors[seed] for seed in reversed(TRAINING_SEEDS)}
    check("seed_mapping_insertion_order_does_not_change_mean", np.array_equal(ref["ensemble_symmetric"], three_seed_predictions(x_ref, reversed_seeds)["ensemble_symmetric"]))
    order = rng.permutation(len(x_ref))
    for chunk_size in (1, 17, 256):
        recovered = np.empty(len(x_ref), dtype=np.float64)
        for start in range(0, len(order), chunk_size):
            selected = order[start:start + chunk_size]
            recovered[selected] = three_seed_predictions(x_ref[selected], predictors)["ensemble_symmetric"]
        check(f"shuffled_chunk_order_invariance_{chunk_size}", np.array_equal(recovered, ref["ensemble_symmetric"]))
    threshold_fixtures = {
        "all_zero": np.zeros(200),
        "all_one": np.ones(200),
        "all_tied_half": np.full(200, 0.5),
        "untied_200": np.linspace(0, 1, 200),
        "boundary_ties": np.array([0.9] * 8 + [0.8] * 5 + [0.1] * 187),
        "nonmultiple_n_211": np.linspace(0, 1, 211),
    }
    thresholds: dict[str, object] = {}
    for name, scores in threshold_fixtures.items():
        result = validation_control_cutoff(scores)
        thresholds[name] = result
        check(f"cutoff_empirical_FPR_bound_{name}", result["called_calibration_controls"] <= len(scores) // 20)
        roundtrip = float(result["cutoff_decimal_17_significant"])
        check(f"cutoff_serialization_roundtrip_{name}", roundtrip == result["cutoff"] == float.fromhex(result["cutoff_hex"]))
    check("boundary_ties_are_not_randomly_split", thresholds["boundary_ties"]["called_calibration_controls"] == 8)
    check("all_one_cutoff_can_exceed_one", thresholds["all_one"]["cutoff"] > 1 and thresholds["all_one"]["called_calibration_controls"] == 0)
    rejects("reject_missing_seed", lambda: three_seed_predictions(x_ref[:1], {TRAINING_SEEDS[0]: predictors[TRAINING_SEEDS[0]]}))
    rejects("reject_nonfinite_prediction", lambda: symmetric_predictions(x_ref[:1], lambda _: np.array([np.nan])))
    rejects("reject_out_of_range_prediction", lambda: symmetric_predictions(x_ref[:1], lambda _: np.array([1.1])))
    rejects("reject_incorrect_prediction_shape", lambda: symmetric_predictions(x_ref[:1], lambda _: np.ones((1, 2))))
    rejects("reject_nonfinite_cutoff_input", lambda: validation_control_cutoff(np.full(200, np.nan)))
    rejects("reject_insufficient_calibration_controls", lambda: validation_control_cutoff(np.zeros(199)))
    rejects("reject_unsupported_base", lambda: reverse_complement_string("ACGTR"))
    rejects("reject_incorrect_onehot_shape", lambda: reverse_complement_one_hot(np.zeros((7, 5))))
    script = Path(__file__).resolve()
    failed = sum(row["status"] == "FAIL" for row in checks)
    return {
        "module": "COPD-V2-PREFLIGHT",
        "scope": "CPU-only synthetic mathematical/wrapper contract; NO trained-network validation",
        "status": "PASS" if not failed else "FAIL",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "training_or_model_loading_performed": False,
        "external_benchmark_or_biological_data_read": False,
        "actual_trained_network_invariance_validation": "NOT_RUN; required after training before release",
        "training_seeds": TRAINING_SEEDS,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "fixture_seed": FIXTURE_SEED,
        "synthetic_ref_alt_pairs": len(references),
        "sequence_length": 2001,
        "absolute_tolerance": ATOL,
        "relative_tolerance": RTOL,
        "reduction_dtype": "float64",
        "python": platform.python_version(),
        "numpy": np.__version__,
        "script_path": str(script),
        "script_sha256": hashlib.sha256(script.read_bytes()).hexdigest(),
        "n_checks": len(checks),
        "n_failures": failed,
        "checks": checks,
        "threshold_fixtures": thresholds,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parents[1] / "provenance" / "COPD-V2-PREFLIGHT_orientation_contract_qc.json")
    args = parser.parse_args()
    result = run_synthetic_qc()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"], "n_checks": result["n_checks"], "n_failures": result["n_failures"], "output": str(args.output)}))
    if result["n_failures"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
