#!/usr/bin/env python3
"""Synthetic-only tests: never open real chr8/9 scores or protected outcomes."""
import copy
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import evaluate_test as subject


def synthetic_rows():
    return [{"interval_id": f"fake-{i}", "chrom": "chr8" if i < 3 else "chr9", "partition": "test",
             "validation_role": "test", "configuration": "V2-C", "model": "enhancer", "label": str(i % 2),
             "component_id": f"component-{i // 2}", "sequence_available": "1", "sequence_length": "2001",
             "input_start": str(i * 3000), "input_end": str(i * 3000 + 2001),
             "sequence_sha256": "a" * 64, "canonical_rc_sequence_sha256": "b" * 64,
             "gc_fraction": str(i / 10), "atac_signal_percentile_max_train_only": str(i / 10),
             "repeat_fraction_2001": str(i / 10), "atac_lobes": "upper_right",
             "enhancer_same_lobe_support_lobes": "upper_right" if i % 2 else "",
             "non_acgt_fraction": "0" if i < 5 else ".01"} for i in range(6)]


def predictor_set(y):
    p = np.where(y, .8, .2)
    return {(config, predictor): p.copy() for config in subject.CONFIGURATIONS
            for predictor in [*(f"seed:{s}" for s in subject.SEEDS), "ensemble"]}


class FrozenArithmeticTests(unittest.TestCase):
    def test_weighted_tied_metrics_independent_sklearn_reference(self):
        rng = np.random.default_rng(777)
        for _ in range(12):
            y, p, w = rng.integers(0, 2, 100), rng.integers(0, 11, 100) / 10, rng.integers(0, 5, 100)
            result = subject.ev.PreparedMetrics(y, p).calculate(w)
            self.assertAlmostEqual(result["AP"], average_precision_score(y, p, sample_weight=w), places=13)
            self.assertAlmostEqual(result["AUROC"], roc_auc_score(y, p, sample_weight=w), places=13)
            self.assertAlmostEqual(result["Brier"], brier_score_loss(y, p, sample_weight=w), places=13)

    def test_ensemble_uses_all_seeds_float64(self):
        f = np.array([[.1, .9], [.2, .8], [.4, .6]], dtype=np.float32)
        r = np.array([[.3, .7], [.5, .5], [.6, .4]], dtype=np.float32)
        q, ensemble = subject.ev.symmetric_ensemble(f, r)
        self.assertEqual(q.dtype, np.float64)
        np.testing.assert_array_equal(ensemble, ((f.astype(float) + r.astype(float)) / 2).mean(axis=0))

    def test_reliability_boundaries_keep_empty_bins(self):
        rows = subject.ev.reliability([0, 1, 1], [0, .1, 1])
        self.assertEqual([r["n"] for r in rows], [1, 1, 0, 0, 0, 0, 0, 0, 0, 1])
        self.assertIsNone(rows[2]["observed_positive_fraction"])

    def test_bootstrap_exactly_reuses_inherited_ensemble_stream(self):
        y = np.array([0, 1, 1] * 12)
        predictors = predictor_set(y)
        components = [f"g-{i // 2}" for i in range(len(y))]
        rows, operations, audit = subject.component_bootstrap(y, predictors, components, .7, valid_replicates=30, max_attempts=300)
        old_rows, old_audit = subject.ev.component_bootstrap(y, {c: predictors[(c, "ensemble")] for c in subject.CONFIGURATIONS}, components, valid_replicates=30, max_attempts=300)
        new_rows = [{k: v for k, v in row.items() if k != "predictor"} for row in rows if row["predictor"] == "ensemble"]
        self.assertEqual(new_rows, old_rows)
        for field in old_audit:
            self.assertEqual(audit[field], old_audit[field])
        self.assertEqual(len(rows), 30 * 12)
        self.assertEqual(len(operations), 30)
        self.assertGreater(len({r["prevalence"] for r in rows}), 1)
        for row in rows:
            self.assertEqual(row["AP_gain"], row["AP"] - row["prevalence"])
            self.assertAlmostEqual(row["BrierSkill"], 1 - row["Brier"] / (row["prevalence"] * (1 - row["prevalence"])))

    def test_bootstrap_invalid_draw_safeguard(self):
        y = np.array([0, 1])
        _, _, audit = subject.component_bootstrap(y, predictor_set(y), ["a", "b"], .5, valid_replicates=100, max_attempts=1000)
        self.assertEqual(audit["valid"], 100)
        self.assertGreaterEqual(audit["invalid_fraction"], .1)
        self.assertEqual(audit["status"], "INCONCLUSIVE")

    def test_bootstrap_attempt_limit_no_imputation(self):
        y = np.array([0, 1])
        rows, _, audit = subject.component_bootstrap(y, predictor_set(y), ["a", "b"], .5, valid_replicates=100, max_attempts=5)
        self.assertEqual(audit["attempted"], 5)
        self.assertLess(audit["valid"], 100)
        self.assertEqual(len(rows), 12 * audit["valid"])
        self.assertEqual(audit["status"], "INCONCLUSIVE")

    def test_component_multiplicity_is_row_weight_not_matching_pair(self):
        y = np.array([1, 0, 1, 0, 0, 1])
        components = ["a", "a", "b", "c", "c", "c"]
        p = predictor_set(y)
        rows, _, _ = subject.component_bootstrap(y, p, components, .7, valid_replicates=1)
        draw = np.random.default_rng(314159).integers(0, 3, size=3)
        multiplicity = np.bincount(draw, minlength=3)
        weights = multiplicity[np.array([0, 0, 1, 2, 2, 2])]
        expected = subject.ev.PreparedMetrics(y, p[("V2-C", "ensemble")]).calculate(weights)
        self.assertEqual(rows[0]["prevalence"], expected["prevalence"])


class ThresholdTests(unittest.TestCase):
    def test_exact_frozen_threshold_encodings(self):
        self.assertEqual(float.fromhex(subject.THRESHOLDS["enhancer"]), float("0.74848511815071117"))
        self.assertEqual(float.fromhex(subject.THRESHOLDS["h3k27me3"]), float("0.76960810025533055"))

    def test_ties_use_greater_equal_and_no_recalibration(self):
        threshold = float.fromhex(subject.THRESHOLDS["enhancer"])
        scores = [np.nextafter(threshold, -np.inf), threshold, threshold, np.nextafter(threshold, np.inf)]
        result = subject.operating_metrics([0, 0, 1, 1], scores, threshold)
        self.assertEqual({k: result[k] for k in ("TP", "FP", "TN", "FN")}, {"TP": 2, "FP": 1, "TN": 1, "FN": 0})
        self.assertEqual(result["negative_row_FPR"], .5)
        self.assertEqual(result["sensitivity"], 1)
        self.assertEqual(result["precision"], 2/3)

    def test_undefined_ratios_explicit_none(self):
        result = subject.operating_metrics([0, 1], [.1, .2], .8)
        self.assertIsNone(result["precision"])
        self.assertEqual(result["NPV"], .5)
        self.assertIsNone(subject.operating_metrics([0, 1], [.1, .2], 0)["NPV"])

    def test_undefined_operating_ratio_does_not_change_draw_stream(self):
        y = np.array([0, 1] * 20)
        p = predictor_set(y)
        rows, ops, audit = subject.component_bootstrap(y, p, [str(i // 2) for i in range(40)], 1.1, valid_replicates=10)
        self.assertEqual(audit["valid"], 10)
        self.assertEqual(audit["invalid"], 0)
        self.assertTrue(all(r["precision"] is None for r in ops))
        self.assertEqual(len(rows), 120)

    def test_weighted_operating_counts_match_expansion(self):
        y, p, w = np.array([0, 0, 1, 1]), np.array([.1, .8, .4, .9]), np.array([3, 2, 4, 5])
        self.assertEqual(subject.operating_metrics(y, p, .5, w), subject.operating_metrics(np.repeat(y, w), np.repeat(p, w), .5))


class GateAndPanelTests(unittest.TestCase):
    def test_exact_panel_counts_components(self):
        rows = synthetic_rows()
        self.assertEqual(subject.validate_test_rows(rows, "enhancer", {"positive": 3, "control": 3, "components": 3}), [r["interval_id"] for r in rows])
        with self.assertRaises(ValueError):
            subject.validate_test_rows(rows, "enhancer")

    def test_panel_rejects_chromosome_role_sequence_and_duplicates(self):
        for field, value in (("chrom", "chr7"), ("partition", "validation"), ("validation_role", "selection"),
                             ("configuration", "V2-B"), ("sequence_available", "0"), ("sequence_length", "2000"),
                             ("input_end", "2000"), ("sequence_sha256", "unknown"), ("label", "2")):
            rows = synthetic_rows()
            rows[0][field] = value
            with self.subTest(field=field), self.assertRaises(ValueError):
                subject.validate_test_rows(rows, "enhancer", {"positive": 3, "control": 3, "components": 3})
        rows = synthetic_rows()
        rows[-1]["interval_id"] = rows[0]["interval_id"]
        with self.assertRaises(ValueError):
            subject.validate_test_rows(rows, "enhancer", {"positive": 3, "control": 3, "components": 3})

    def test_stratification_uses_frozen_right_boundary(self):
        rows = synthetic_rows()
        cuts = {"model": "enhancer", "quantile_method": "linear", "cutpoints": {c: [.1, .2, .2, .4] for c in subject.COVARIATES}}
        strata = subject.stratification(rows, cuts, "enhancer")
        self.assertEqual(strata["gc_fraction"].tolist(), ["quintile_1", "quintile_2", "quintile_4", "quintile_4", "quintile_5", "quintile_5"])
        self.assertEqual(strata["ambiguous_base_status"][-1], "present")
        self.assertEqual(strata["same_lobe_mark_support"][0], "none")

    def test_strict_C_gate_boundaries_and_seed_sample_SD(self):
        expected = {"positive": 100, "control": 200, "components": 300}
        labels, components = [1]*100 + [0]*200, list(map(str, range(300)))
        cis = {"AUROC": {"lower95": .5}, "AP_gain": {"lower95": 0}, "BrierSkill": {"lower95": 0}}
        rows, status = subject.holdout_gates("enhancer", labels, components, [.4, .46, .52], cis,
                                           {"valid": 2000, "invalid_fraction": 0}, expected)
        gates = {r["gate"]: r for r in rows}
        self.assertEqual(status, "FAIL")
        for metric in cis:
            self.assertEqual(gates[f"lower95_{metric}"]["status"], "FAIL")
        self.assertAlmostEqual(gates["seed_AP_sample_SD"]["observed"], .06)
        self.assertEqual(gates["seed_AP_range"]["status"], "FAIL")

    def test_positive_gates_do_not_require_C_superiority_or_low_test_FPR(self):
        expected = {"positive": 3, "control": 3, "components": 3}
        cis = {"AUROC": {"lower95": .51}, "AP_gain": {"lower95": .001}, "BrierSkill": {"lower95": .001}}
        rows, status = subject.holdout_gates("enhancer", [0, 1]*3, ["a", "a", "b", "b", "c", "c"],
                                          [.7, .7, .7], cis, {"valid": 2000, "invalid_fraction": .099}, expected)
        self.assertEqual(status, "PASS")
        self.assertFalse(any("superiority" in r["gate"] or "FPR" in r["gate"] for r in rows))

    def test_invalid_fraction_equality_is_inconclusive(self):
        cis = {"AUROC": {"lower95": .6}, "AP_gain": {"lower95": .1}, "BrierSkill": {"lower95": .1}}
        _, status = subject.holdout_gates("enhancer", [0, 1], ["a", "b"], [.7]*3, cis,
                                         {"valid": 2000, "invalid_fraction": .1}, {"positive": 1, "control": 1, "components": 2})
        self.assertEqual(status, "INCONCLUSIVE")

    def test_prediction_exactness_and_extra_row_firewall(self):
        metadata = synthetic_rows()
        rows = []
        for i, record in enumerate(metadata):
            row = {"interval_id": record["interval_id"], "chrom": record["chrom"]}
            for config in subject.CONFIGURATIONS:
                for seed in subject.SEEDS:
                    row[f"{config}_p_forward_{seed}"] = .25
                    row[f"{config}_p_rc_{seed}"] = .75
                    row[f"{config}_q_{seed}"] = .5
                row[f"{config}_ensemble_q"] = .5
            rows.append(row)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "synthetic.tsv"
            subject.ev.write_rows(path, rows)
            parsed = subject.load_test_predictions(path, metadata)
            np.testing.assert_array_equal(parsed["V2-C"]["ensemble"], np.full(6, .5))
            bad = copy.deepcopy(rows)
            bad[0]["V2-C_q_104729"] = .5000000001
            bad_path = Path(directory) / "bad.tsv"
            subject.ev.write_rows(bad_path, bad)
            with self.assertRaises(ValueError):
                subject.load_test_predictions(bad_path, metadata)
            extra = copy.deepcopy(rows)
            extra.append({**rows[0], "interval_id": "not-authorized", "chrom": "chr1"})
            extra_path = Path(directory) / "extra.tsv"
            subject.ev.write_rows(extra_path, extra)
            with self.assertRaises(ValueError):
                subject.load_test_predictions(extra_path, metadata)


if __name__ == "__main__":
    unittest.main()
