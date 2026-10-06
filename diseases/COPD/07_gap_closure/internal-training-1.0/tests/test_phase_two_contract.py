#!/usr/bin/env python3
"""CPU-only implementation tests; deliberately do not import TensorFlow/Keras."""
import csv
import gzip
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import phase_two_contract as contract
from predict_chr7 import numeric_gate
from train_phase_two import prepare_attempt


class ContractTests(unittest.TestCase):
    def test_rng_substreams_and_genomic_assignment(self):
        seed, epoch, n = 104729, 7, 1031
        order, reverse = contract.epoch_plan(seed, epoch, n)
        self.assertTrue(np.array_equal(order, np.random.default_rng(np.random.SeedSequence([seed, 1, epoch])).permutation(n)))
        self.assertTrue(np.array_equal(reverse, np.random.default_rng(np.random.SeedSequence([seed, 2, epoch])).random(n) < .5))
        self.assertEqual(sorted(order.tolist()), list(range(n)))
        self.assertTrue(np.array_equal(reverse, contract.epoch_plan(seed, epoch, n)[1]))
        self.assertFalse(np.array_equal(reverse, contract.epoch_plan(seed, epoch + 1, n)[1]))

    def test_orientation_lookup_never_reverses_features(self):
        cache = contract.FeatureCache.__new__(contract.FeatureCache)
        cache.canonical = np.arange(4 * 4560, dtype=np.float32).reshape(4, 4560) / 20000
        cache.rc = np.ones((4, 4560), dtype=np.float32) - cache.canonical
        table = {"cache_row": np.array([2, 0, 3, 1]), "forward_is_canonical": np.array([True, False, True, False])}
        indices = np.array([3, 0, 2, 1])
        forward = cache.batch(table, indices)[..., 0]
        reverse = cache.batch(table, indices, reverse=True)[..., 0]
        chosen = np.array([False, True, True, False])
        for i, index in enumerate(indices):
            r = table["cache_row"][index]
            self.assertTrue(np.array_equal(forward[i], cache.canonical[r] if chosen[i] else cache.rc[r]))
            self.assertTrue(np.array_equal(reverse[i], cache.rc[r] if chosen[i] else cache.canonical[r]))
        mixed = cache.batch(table, indices, np.array([True, False, True, False]))[..., 0]
        self.assertTrue(np.array_equal(mixed[[0, 2]], reverse[[0, 2]]))
        self.assertTrue(np.array_equal(mixed[[1, 3]], forward[[1, 3]]))

    def test_float64_symmetry_and_ordered_three_seed_mean(self):
        f = np.array([.1, .25, 1], dtype=np.float32)
        r = np.array([.7, .5, 0], dtype=np.float32)
        q = contract.symmetric_probability(f, r)
        self.assertEqual(q.dtype, np.float64)
        self.assertTrue(np.array_equal(q, contract.symmetric_probability(r, f)))
        qs = [q, q / 2, q / 3]
        self.assertTrue(np.array_equal(contract.seed_ensemble(qs), np.mean(np.stack(qs), axis=0, dtype=np.float64)))
        with self.assertRaises(RuntimeError):
            contract.seed_ensemble(qs[:2])

    def test_invalid_probability_stops(self):
        for p in ([np.nan], [np.inf], [-1e-12], [1 + 1e-12]):
            with self.assertRaises(RuntimeError):
                contract.valid_probabilities(p)

    def test_symmetric_bce_finite_boundaries(self):
        got = contract.symmetric_bce([0, 1], [0, 1])
        self.assertTrue(np.isfinite(got))
        self.assertAlmostEqual(got, -np.log(1 - 1e-7), places=15)
        self.assertAlmostEqual(contract.symmetric_bce([0, 1], [.5, .5]), np.log(2), places=15)

    def test_frozen_invariance_bound_is_not_an_exact_equality_test(self):
        self.assertEqual(numeric_gate(np.array([.5]), np.array([.500001]))["status"], "PASS")
        self.assertEqual(numeric_gate(np.array([.5]), np.array([.500002]))["status"], "FAIL")

    def test_cpu_fake_inference_order_restored(self):
        cache = contract.FeatureCache.__new__(contract.FeatureCache)
        cache.canonical = np.zeros((17, 4560), dtype=np.float32)
        cache.rc = np.zeros((17, 4560), dtype=np.float32)
        cache.canonical[:, 0] = np.arange(17) / 20
        cache.rc[:, 0] = 1 - np.arange(17) / 20
        table = {"cache_row": np.arange(17), "forward_is_canonical": np.arange(17) % 2 == 0}
        def fake_infer(x):
            return x[:, 0, :]
        p = contract.predict_orientation(fake_infer, cache, table, batch_size=5)
        p2 = contract.predict_orientation(fake_infer, cache, table, reverse_order=True, batch_size=3)
        self.assertTrue(np.array_equal(p, p2))

    def test_input_firewalls(self):
        with tempfile.TemporaryDirectory(prefix="phase-two-contract-test-") as temporary:
            path = Path(temporary) / "rows.tsv"
            rows = [{"interval_id": "chr8:1-1001", "chrom": "chr8", "label": "0", "cache_row": "0", "forward_is_canonical": "1", "validation_role": "test"}]
            with path.open("w") as stream:
                writer = csv.DictWriter(stream, delimiter="\t", fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            for mode in ("train", "checkpoint", "validation"):
                with self.assertRaises(RuntimeError):
                    contract.load_input_table(path, mode)
            rows[0]["chrom"] = "chr7"
            rows[0]["validation_role"] = "selection"
            with path.open("w") as stream:
                writer = csv.DictWriter(stream, delimiter="\t", fieldnames=list(rows[0]))
                writer.writeheader()
                writer.writerows(rows)
            with self.assertRaises(RuntimeError):
                contract.load_input_table(path, "checkpoint")

    def test_no_completed_seed_can_repeat_and_retry_needs_authority(self):
        with tempfile.TemporaryDirectory(prefix="phase-two-attempt-test-") as temporary:
            stage = Path(temporary)
            output, _ = prepare_attempt(stage, "V2-A_enhancer_seed104729", 1, None)
            with self.assertRaises(RuntimeError):
                prepare_attempt(stage, "V2-A_enhancer_seed104729", 2, None)
            contract.write_json(output / "completed.json", {"status": "COMPLETED"})
            with self.assertRaises(RuntimeError):
                prepare_attempt(stage, "V2-A_enhancer_seed104729", 2, None)

    def test_no_model_runtime_imported(self):
        self.assertNotIn("tensorflow", sys.modules)
        self.assertNotIn("keras", sys.modules)


if __name__ == "__main__":
    unittest.main()
