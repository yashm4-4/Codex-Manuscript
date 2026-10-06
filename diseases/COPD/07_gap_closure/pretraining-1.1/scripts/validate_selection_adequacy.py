#!/usr/bin/env python3
"""CPU-only specification audit and synthetic decision-boundary checks.

No model, biological outcome table, or external benchmark is loaded. The
synthetic numbers below are contract fixtures, not estimated performance.
"""

from __future__ import annotations

import hashlib
import inspect
import json
import math
from pathlib import Path
import random
import statistics
import sys


HERE = Path(__file__).resolve()
VERSION = HERE.parent.parent
REPO = HERE.parents[5]
NEW = VERSION / "specification/selection_adequacy.json"
OLD = VERSION.parent / "provenance/COPD-V2-PREFLIGHT_selection_specification.json"
OUTPUT = VERSION / "provenance/selection_adequacy_review.json"


def hash_file(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def c_head_passes(c, specification):
    """Synthetic contract predicate; no B/A input can affect the decision."""
    s = specification["selection"]
    a = s["adequacy"]
    minimum = s["minimum_selection"]
    return (
        all(c[key] for key in ("construction", "all_seeds", "invariance", "finite", "bootstrap_valid"))
        and c["positive"] >= minimum["positive"]
        and c["control"] >= minimum["control"]
        and c["components"] >= minimum["components"]
        and c["lower_AUROC"] > a["lower95_AUROC_strictly_gt"]
        and c["lower_AP_gain"] > a["lower95_AP_minus_replicate_prevalence_strictly_gt"]
        and c["lower_BrierSkill"] > a["lower95_BrierSkill_strictly_gt"]
        and c["AP_SD"] <= s["seed_stability_hard_stop"]["AP_sample_SD_gt"]
        and c["AP_range"] <= s["seed_stability_hard_stop"]["AP_range_gt"]
    )


def cutoff(scores, target=0.05):
    if len(scores) < 200 or any(not math.isfinite(x) or not 0 <= x <= 1 for x in scores):
        raise ValueError("invalid synthetic calibration fixture")
    descending = sorted(scores, reverse=True)
    return math.nextafter(descending[math.floor(target * len(scores))], math.inf)


def main():
    s = json.loads(NEW.read_text())
    old = json.loads(OLD.read_text())
    checks = []

    def check(name, passed, detail=None):
        row = {"check": name, "status": "PASS" if passed else "FAIL"}
        if detail is not None:
            row["detail"] = detail
        checks.append(row)

    check("historical_specification_hash_unchanged", hash_file(OLD) == s["pretraining_1_0"]["specification_sha256"])
    for key in ("models", "seeds", "n_planned_fits", "phase_I_weights_sha256", "training", "augmentation", "checkpoint"):
        check(f"unchanged_{key}", s[key] == old[key])
    for key, value in old["inference"].items():
        check(f"unchanged_inference_{key}", s["inference"][key] == value)
    for key in ("chr", "unit", "rule", "assignment"):
        check(f"unchanged_validation_roles_{key}", s["validation_roles"][key] == old["validation_roles"][key])
    for key in ("seed", "valid_replicates", "max_attempts", "max_invalid_fraction_exclusive", "unit", "paired", "resample", "CI95_percentiles", "CI90_percentiles"):
        check(f"unchanged_bootstrap_{key}", s["selection"]["bootstrap"][key] == old["selection"]["bootstrap"][key])
    for key in ("minimum_selection", "seed_stability_hard_stop", "reliability", "test_or_external_tiebreaker"):
        check(f"unchanged_selection_{key}", s["selection"][key] == old["selection"][key])
    for key in ("target_empirical_FPR", "minimum_controls", "formula", "call", "tie_policy", "storage", "test_recalibration", "variant_FPR_guarantee"):
        check(f"unchanged_threshold_{key}", s["threshold"][key] == old["threshold"][key])

    check("training_not_authorized_or_started", s["training_authorized"] is False and s["training_started"] is False)
    check("exact_fixed_seed_list", s["seeds"] == [104729, 130363, 155921])
    check("exact_18_prospective_fits", s["n_planned_fits"] == 3 * 2 * len(s["seeds"]))
    check("C_only_eligible_configuration", [k for k, v in s["configurations"].items() if v["eligible_final_configuration"]] == ["V2-C"])
    check("no_comparative_margins_for_C_retention", s["selection"]["paired_ablation_report"]["comparative_margins_used_for_retention"] is False)
    check("no_comparative_hypothesis_gate", s["selection"]["paired_ablation_report"]["superiority_equivalence_noninferiority_required"] is False)
    check("common_panels_exclude_checkpoint", s["common_panels"]["roles"] == ["selection", "calibration", "test"] and s["common_panels"]["checkpoint_included"] is False)
    check("panels_frozen_before_training", s["common_panels"]["frozen_before_training"] is True)
    check("no_new_allele_delta_cutoff", s["candidate_plan"]["allele_delta_cutoff"] is None)
    check("calibration_component_floor_explicit_conservative_clarification", s["threshold"]["minimum_control_components"] == 30 and "new conservative" in s["threshold"]["component_minimum_provenance"])
    check("absolute_null_margins_preserved", s["numerical_margin_review"]["absolute_null_margins_preserved"] == {"AUROC": 0.5, "AP_gain": 0.0, "BrierSkill": 0.0})
    check("no_outcome_based_tuning", s["numerical_margin_review"]["outcome_based_tuning"] is False)
    check("real_network_gate_not_claimed_run", s["inference"]["real_network_numeric_gate"] == "NOT_RUN_REQUIRED_AFTER_TRAINING")
    check("decision_predicate_has_no_A_or_B_input", list(inspect.signature(c_head_passes).parameters) == ["c", "specification"])

    fixture = {
        "construction": True, "all_seeds": True, "invariance": True,
        "finite": True, "bootstrap_valid": True, "positive": 100,
        "control": 200, "components": 30, "lower_AUROC": 0.501,
        "lower_AP_gain": 0.001, "lower_BrierSkill": 0.001,
        "AP_SD": 0.03, "AP_range": 0.06,
    }
    check("synthetic_adequacy_above_null_passes_without_0_02_improvement", c_head_passes(fixture, s))
    for key, boundary in (("lower_AUROC", 0.5), ("lower_AP_gain", 0.0), ("lower_BrierSkill", 0.0)):
        changed = dict(fixture, **{key: boundary})
        check(f"synthetic_{key}_null_equality_does_not_pass", not c_head_passes(changed, s))
    for key in ("construction", "all_seeds", "invariance", "finite", "bootstrap_valid"):
        check(f"synthetic_{key}_failure_stops", not c_head_passes(dict(fixture, **{key: False}), s))
    for key, count in (("positive", 99), ("control", 199), ("components", 29)):
        check(f"synthetic_insufficient_{key}_stops", not c_head_passes(dict(fixture, **{key: count}), s))
    check("synthetic_seed_SD_above_boundary_stops", not c_head_passes(dict(fixture, AP_SD=math.nextafter(0.03, math.inf)), s))
    check("synthetic_both_C_heads_required", not all(c_head_passes(x, s) for x in (fixture, dict(fixture, lower_AUROC=0.5))))
    check("synthetic_two_adequate_C_heads_pass", all(c_head_passes(x, s) for x in (fixture, fixture)))

    # The extreme pair contributes at least r^2/2 to the squared deviations;
    # with n=3 and ddof=1, sample variance is at least r^2/4.
    rng = random.Random(314159)
    triples = [[0.0, 0.5, 1.0], [0.0, 0.0, 1.0]] + [[rng.random() for _ in range(3)] for _ in range(1000)]
    check("synthetic_three_seed_SD_range_bound_1002_cases", all(statistics.stdev(x) + 1e-15 >= (max(x) - min(x)) / 2 for x in triples))
    check("range_limit_redundant_not_contradictory", 2 * s["selection"]["seed_stability_hard_stop"]["AP_sample_SD_gt"] < s["selection"]["seed_stability_hard_stop"]["AP_range_gt"])

    calibration = {
        "all_zero": [0.0] * 200,
        "all_one": [1.0] * 200,
        "all_half": [0.5] * 200,
        "untied": [i / 199 for i in range(200)],
        "nonmultiple": [i / 210 for i in range(211)],
        "boundary_ties": [0.9] * 11 + [0.1] * 189,
    }
    threshold_fixture_results = {}
    for name, scores in calibration.items():
        t = cutoff(scores)
        calls = sum(x >= t for x in scores)
        check(f"synthetic_threshold_FPR_bound_{name}", calls <= math.floor(0.05 * len(scores)))
        check(f"synthetic_threshold_hex_roundtrip_{name}", float.fromhex(t.hex()) == t)
        check(f"synthetic_threshold_decimal_roundtrip_{name}", float(format(t, ".17g")) == t)
        threshold_fixture_results[name] = {"n": len(scores), "called": calls, "threshold_hex": t.hex()}
    check("synthetic_all_one_threshold_not_clipped", cutoff([1.0] * 200) > 1)
    for name, scores in (("too_few", [0.5] * 199), ("nonfinite", [float("nan")] * 200), ("out_of_range", [1.1] * 200)):
        try:
            cutoff(scores)
        except ValueError:
            rejected = True
        else:
            rejected = False
        check(f"synthetic_threshold_rejects_{name}", rejected)

    review = {
        "experimental_version": "pretraining-1.1",
        "audit": "CPU specification invariants and synthetic decision-boundary contracts only",
        "specification_path": str(NEW.relative_to(REPO)),
        "specification_sha256": hash_file(NEW),
        "historical_specification_sha256": hash_file(OLD),
        "script_path": str(HERE.relative_to(REPO)),
        "script_sha256": hash_file(HERE),
        "checks": checks,
        "check_count": len(checks),
        "pass_count": sum(c["status"] == "PASS" for c in checks),
        "failure_count": sum(c["status"] != "PASS" for c in checks),
        "synthetic_threshold_fixtures": threshold_fixture_results,
        "numerical_review_findings": [
            "Comparative 0.02 AP/AUROC/Brier margins are not absolute adequacy gates and no longer determine C retention.",
            "The three-seed sample-SD limit already implies a stricter range ceiling than the preserved 0.10 range stop; this is redundancy, not a conflicting requirement.",
            "Strict lower95 AUROC>0.5, AP-over-replicate-prevalence>0 and Brier-skill>0 are preserved without tuning.",
            "The new calibration-control component floor of 30 is an explicit conservative clarification; the existing 200-control floor is unchanged.",
            "The observed calibration-row FPR bound is not a population, donor-generalization or variant-level FPR guarantee."
        ],
        "model_training_run": False,
        "model_inference_run": False,
        "phase_I_extraction_run": False,
        "external_benchmark_outcomes_parsed": False,
        "real_network_invariance": "NOT_RUN_REQUIRED_AFTER_SEPARATELY_AUTHORIZED_TRAINING",
        "C_absolute_adequacy": "NOT_EVALUATED_NO_MODEL_PREDICTIONS",
        "construction_panel_minima": "NOT_ASSESSED_BY_THIS_SPECIFICATION_ONLY_REVIEW; require independent final interval validation",
        "training_authorized": False,
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(review, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"checks": review["check_count"], "passed": review["pass_count"], "failed": review["failure_count"], "output": str(OUTPUT.relative_to(REPO))}))
    return bool(review["failure_count"])


if __name__ == "__main__":
    sys.exit(main())
