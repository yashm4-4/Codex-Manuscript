"""Synthetic tests for prespecified training-only ECDF contract."""
import importlib.util
from pathlib import Path
import unittest

import numpy as np

PATH = Path(__file__).resolve().parents[1] / "scripts/normalize_atac.py"
SPEC = importlib.util.spec_from_file_location("normalize_atac", PATH)
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class AtacNormalizationTests(unittest.TestCase):
    def test_ties_boundaries_and_steps(self):
        mapping = MODULE.fit([5, 1, 3, 1])
        actual = MODULE.apply_mapping(mapping, [-1, 1, 2, 3, 4, 5, 6])
        np.testing.assert_array_equal(actual, [0, .25, .5, .625, .75, .875, 1])

    def test_identical_values_share_midpoint(self):
        mapping = MODULE.fit([7, 7, 7])
        np.testing.assert_array_equal(MODULE.apply_mapping(mapping, [6, 7, 8]), [0, .5, 1])

    def test_heldout_replacement_cannot_change_training_mapping(self):
        train = [1, 2, 3]
        mapping = MODULE.fit(train)
        before = mapping.copy()
        MODULE.apply_mapping(mapping, [-1e20, 1e20, 2.5])
        np.testing.assert_array_equal(mapping, before)
        np.testing.assert_array_equal(mapping, MODULE.fit(train))

    def test_partition_framework(self):
        self.assertEqual(MODULE.TRAIN & {"chr7", "chr8", "chr9"}, set())
        self.assertEqual(MODULE.partition("chr7"), "validation")
        self.assertEqual(MODULE.partition("chr8"), "test")
        self.assertEqual(MODULE.partition("chr9"), "test")
        self.assertTrue(all(MODULE.partition(chrom) == "train" for chrom in MODULE.TRAIN))
        with self.assertRaises(ValueError):
            MODULE.partition("chrM")

    def test_invalid_values_stop(self):
        for values in ([], [float("nan")], [float("inf")]):
            with self.assertRaises(ValueError):
                MODULE.fit(values)
        with self.assertRaises(ValueError):
            MODULE.apply_mapping(MODULE.fit([1]), [float("nan")])

    def test_order_independence(self):
        np.testing.assert_array_equal(MODULE.fit([3, 1, 1, 2]), MODULE.fit([1, 3, 2, 1]))


if __name__ == "__main__":
    unittest.main()
