"""Independent synthetic audit-helper tests; no actual model execution."""
import sys
from pathlib import Path
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from validate_internal_stage import calibration_threshold, expected_checkpoint, independent_metrics, symmetric_ensemble


class ExecutionContractTests(unittest.TestCase):
    def test_earliest_minimum(self):
        self.assertEqual(expected_checkpoint([.6, .5, .5, .7]), 2)

    def test_single_epoch_checkpoint(self):
        self.assertEqual(expected_checkpoint([.4]), 1)

    def test_nan_history_rejected(self):
        with self.assertRaises(ValueError):
            expected_checkpoint([.6, np.nan])

    def test_long_history_rejected(self):
        with self.assertRaises(ValueError):
            expected_checkpoint([.5] * 51)

    def test_symmetric_ensemble_invariant(self):
        f = np.array([[.1, .2, .3], [.01, .31, .51]], dtype=np.float32)
        r = np.array([[.3, .7, .6], [.91, .73, .21]], dtype=np.float32)
        np.testing.assert_array_equal(symmetric_ensemble(f, r), symmetric_ensemble(r, f))
        self.assertEqual(symmetric_ensemble(f, r).dtype, np.float64)

    def test_all_seeds_required(self):
        with self.assertRaises(ValueError):
            symmetric_ensemble(np.zeros((3, 2)), np.zeros((3, 2)))

    def test_prediction_bounds(self):
        with self.assertRaises(ValueError):
            symmetric_ensemble(np.full((3, 3), 1.01), np.zeros((3, 3)))

    def test_conservative_ties(self):
        result = calibration_threshold(np.full(201, .5))
        self.assertEqual(result["calls"], 0)
        self.assertGreater(result["threshold"], .5)
        self.assertEqual(float.fromhex(result["threshold_hex"]), float(result["threshold_decimal"]))

    def test_exact_fpr_ceiling(self):
        result = calibration_threshold(np.linspace(0, 1, 201))
        self.assertEqual(result["calls"], 10)
        self.assertEqual(result["k"], 10)

    def test_unclipped_one_boundary(self):
        result = calibration_threshold(np.ones(201))
        self.assertGreater(result["threshold"], 1)
        self.assertEqual(result["calls"], 0)

    def test_independent_perfect_metrics(self):
        metric = independent_metrics([0,0,1,1], [0,0,1,1])
        self.assertEqual(metric, {"AP": 1.0, "AUROC": 1.0, "Brier": 0.0,
                                  "BrierSkill": 1.0, "AP_gain": .5, "prevalence": .5})

    def test_independent_tied_null_metrics(self):
        metric = independent_metrics([0,1], [.5,.5])
        self.assertEqual(metric["AP"], .5)
        self.assertEqual(metric["AUROC"], .5)
        self.assertEqual(metric["BrierSkill"], 0)
        self.assertEqual(metric["AP_gain"], 0)

    def test_independent_component_weights_equal_row_repetition(self):
        y, p, w = np.array([0,1,0,1]), np.array([.2,.7,.5,.6]), np.array([2,1,3,2])
        metric = independent_metrics(y,p,w)
        repeated = independent_metrics(np.repeat(y,w),np.repeat(p,w))
        for key in metric:
            self.assertAlmostEqual(metric[key], repeated[key], places=14)

    def test_independent_choice_rng_matches_spec_integer_draw(self):
        left, right = np.random.default_rng(314159), np.random.default_rng(314159)
        for _ in range(20):
            np.testing.assert_array_equal(left.choice(43,size=43,replace=True), right.integers(0,43,size=43))


if __name__ == "__main__":
    unittest.main()
