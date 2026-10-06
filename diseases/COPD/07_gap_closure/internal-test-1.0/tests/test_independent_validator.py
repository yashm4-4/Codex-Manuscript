"""Synthetic checks for independent QC arithmetic; no model/data access."""
import importlib.util
from pathlib import Path
import sys
import unittest

import numpy as np
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

sys.dont_write_bytecode = True
path = Path(__file__).resolve().parents[1] / "scripts/validate_test_stage.py"
spec = importlib.util.spec_from_file_location("test_independent_validator_module", path)
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)


class IndependentValidatorTests(unittest.TestCase):
    def test_weighted_ties_and_zero_weight_rows_against_sklearn(self):
        rng = np.random.default_rng(42)
        for n in (2, 7, 100):
            for _ in range(20):
                y = rng.integers(0, 2, n)
                y[:2] = (0, 1)
                p = rng.integers(0, 5, n) / 4
                w = rng.integers(0, 4, n)
                w[:2] = 1
                result = v.IndependentMetrics(y, p).calculate(w)
                self.assertAlmostEqual(result["AP"], average_precision_score(y, p, sample_weight=w), places=12)
                self.assertAlmostEqual(result["AUROC"], roc_auc_score(y, p, sample_weight=w), places=12)
                self.assertAlmostEqual(result["Brier"], brier_score_loss(y, p, sample_weight=w), places=12)

    def test_bootstrap_null_uses_replicate_prevalence(self):
        y, p, w = np.array([0, 1, 1]), np.array([.1, .3, .8]), np.array([4, 1, 2])
        result = v.IndependentMetrics(y, p).calculate(w)
        self.assertAlmostEqual(result["prevalence"], 3 / 7)
        self.assertAlmostEqual(result["AP_gain"], result["AP"] - 3 / 7)
        self.assertAlmostEqual(result["BrierSkill"], 1 - result["Brier"] / ((3 / 7) * (4 / 7)))

    def test_one_class_component_draw_invalid(self):
        with self.assertRaises(ValueError):
            v.IndependentMetrics([0, 1], [.1, .9]).calculate([1, 0])

    def test_threshold_equality_and_undefined_ratios(self):
        threshold = float.fromhex(v.THRESHOLDS["enhancer"])
        y = np.array([1, 0])
        p = np.array([threshold, np.nextafter(threshold, -np.inf)])
        result = v.independent_operating(y, p, threshold, np.ones(2))
        self.assertEqual((result["TP"], result["FP"], result["TN"], result["FN"]), (1, 0, 1, 0))
        result = v.independent_operating(y, np.zeros(2), threshold, np.ones(2))
        self.assertIsNone(result["precision"])
        self.assertEqual(result["negative_row_FPR"], 0)


if __name__ == "__main__":
    unittest.main()
