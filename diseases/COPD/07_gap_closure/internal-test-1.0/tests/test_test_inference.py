#!/usr/bin/env python3
"""Synthetic inference-contract tests; no TensorFlow, checkpoint or test data read."""
import ast
import contextlib
import copy
import io
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest
from unittest import mock

import numpy as np

sys.dont_write_bytecode = True
SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
import predict_test as inference
import test_contract as contract
import phase_two_contract as frozen


def panel_rows():
    return [{"interval_id": f"synthetic-{i}", "chrom": "chr8" if i < 2 else "chr9",
             "partition": "test", "validation_role": "test", "configuration": "V2-C", "model": "enhancer",
             "label": str(i % 2), "component_id": f"component-{i // 2}",
             "cache_row": str(i), "forward_is_canonical": str(i % 2)} for i in range(4)]


def table(rows=None):
    rows = panel_rows() if rows is None else rows
    return {"rows": rows, "interval_id": np.array([r["interval_id"] for r in rows]),
            "cache_row": np.array([int(r["cache_row"]) for r in rows]),
            "forward_is_canonical": np.array([contract.bool_value(r["forward_is_canonical"]) for r in rows])}


def synthetic_cache(stage):
    cache = contract.FeatureCache.__new__(contract.FeatureCache)
    cache.manifest_path = Path(stage) / "cache/cache_manifest.json"
    cache.canonical = np.full((4, 4560), .2, dtype=np.float32)
    cache.rc = np.full((4, 4560), .8, dtype=np.float32)
    return cache


class AdapterFirewallTests(unittest.TestCase):
    def test_adapter_preserves_original_helper_identities(self):
        for name in ("FeatureCache", "predict_orientation", "symmetric_probability", "seed_ensemble", "inference_function"):
            self.assertIs(getattr(contract, name), getattr(frozen, name))
        self.assertEqual(contract.SEEDS, (104729, 130363, 155921))
        self.assertEqual(contract.ATOL, 1e-6)
        self.assertEqual(contract.RTOL, 1e-6)

    def test_exact_small_synthetic_test_table(self):
        with mock.patch.object(contract, "read_rows", return_value=panel_rows()), mock.patch.dict(contract.EXPECTED_COUNTS, {"enhancer": (4, 2, 2)}):
            actual = contract.load_test_table(Path("synthetic.tsv"), "enhancer")
        self.assertEqual(actual["interval_id"].tolist(), [f"synthetic-{i}" for i in range(4)])
        self.assertEqual(actual["forward_is_canonical"].dtype, bool)

    def test_test_table_rejects_wrong_task_role_chromosome_and_counts(self):
        cases = [("chrom", "chr7"), ("partition", "train"), ("validation_role", "selection"),
                 ("configuration", "V2-B"), ("model", "h3k27me3"), ("label", "2"),
                 ("component_id", "new-component"), ("forward_is_canonical", "unknown")]
        for field, value in cases:
            rows = panel_rows()
            rows[0][field] = value
            with self.subTest(field=field), mock.patch.object(contract, "read_rows", return_value=rows), mock.patch.dict(contract.EXPECTED_COUNTS, {"enhancer": (4, 2, 2)}):
                with self.assertRaises((RuntimeError, ValueError)):
                    contract.load_test_table(Path("synthetic.tsv"), "enhancer")

    def test_duplicate_test_id_rejected(self):
        rows = panel_rows()
        rows[1]["interval_id"] = rows[0]["interval_id"]
        with mock.patch.object(contract, "read_rows", return_value=rows), mock.patch.dict(contract.EXPECTED_COUNTS, {"enhancer": (4, 2, 2)}):
            with self.assertRaises(RuntimeError):
                contract.load_test_table(Path("synthetic.tsv"), "enhancer")

    def test_cache_uses_independent_orientation_rows_not_feature_reversal(self):
        cache = synthetic_cache(Path("synthetic"))
        data = table()
        forward = cache.batch(data, np.arange(4), reverse=False)
        reverse = cache.batch(data, np.arange(4), reverse=True)
        self.assertEqual(forward.shape, (4, 4560, 1))
        np.testing.assert_array_equal(forward[:, 0, 0], np.array([.8, .2, .8, .2], dtype=np.float32))
        np.testing.assert_array_equal(reverse[:, 0, 0], np.array([.2, .8, .2, .8], dtype=np.float32))
        self.assertFalse(np.array_equal(reverse, forward[:, ::-1, :]))

    def test_cache_rejects_bad_index_or_feature(self):
        cache = synthetic_cache(Path("synthetic"))
        data = table()
        data["cache_row"][0] = -1
        with self.assertRaises(RuntimeError):
            cache.check_indices(data)
        data = table()
        cache.rc[0, 0] = np.nan
        with self.assertRaises(RuntimeError):
            cache.batch(data, np.arange(4))

    def test_predict_orientation_preserves_row_order_and_float64(self):
        cache = synthetic_cache(Path("synthetic"))
        cache.canonical[:] = np.arange(4, dtype=np.float32)[:, None] / 4
        cache.rc[:] = 1 - cache.canonical
        fn = lambda x: x[:, 0, 0, None]
        normal = contract.predict_orientation(fn, cache, table(), batch_size=2)
        reversed_order = contract.predict_orientation(fn, cache, table(), reverse_order=True, batch_size=3)
        np.testing.assert_array_equal(normal, reversed_order)
        self.assertEqual(normal.dtype, np.float64)

    def test_invalid_orientation_probability_and_seed_count_hard_stop(self):
        for value in (np.nan, np.inf, -.001, 1.001):
            with self.subTest(value=value), self.assertRaises(RuntimeError):
                contract.predict_orientation(lambda x: np.full((len(x), 1), value), synthetic_cache(Path("synthetic")), table())
        with self.assertRaises(RuntimeError):
            contract.seed_ensemble([np.array([.4]), np.array([.6])])
        with self.assertRaises(RuntimeError):
            contract.symmetric_probability(np.array([np.nan]), np.array([.4]))

    def test_numeric_gate_strict_tolerance_and_finite_inputs(self):
        row = inference.numeric_gate([.5, .2], [.5, .2])
        self.assertEqual(row["status"], "PASS")
        self.assertEqual(row["n_failed"], 0)
        self.assertEqual(inference.numeric_gate([.50001], [.5])["status"], "FAIL")
        for a, b in (([np.nan], [.5]), ([np.inf], [.5]), ([-.1], [.5]), ([1.1], [.5]), ([.5, .2], [.5])):
            with self.subTest(a=a), self.assertRaises((RuntimeError, ValueError)):
                inference.numeric_gate(a, b)

    def test_entrypoint_ast_has_four_real_orientation_calls_no_fit_or_save(self):
        tree = ast.parse((SCRIPTS / "predict_test.py").read_text())
        calls = [node for node in ast.walk(tree) if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == "predict_orientation"]
        self.assertEqual(len(calls), 4)
        self.assertEqual([next(k.value.value for k in node.keywords if k.arg == "reverse") for node in calls], [False, True, True, False])
        forbidden = [node.func.attr for node in ast.walk(tree) if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr in ("fit", "train_on_batch", "save", "save_weights")]
        self.assertEqual(forbidden, [])
        loads = [node for node in ast.walk(tree) if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute) and node.func.attr == "load_model"]
        self.assertEqual(len(loads), 1)
        self.assertIs(next(k.value.value for k in loads[0].keywords if k.arg == "compile"), False)


class CheckpointIdentityTests(unittest.TestCase):
    def fixture(self, directory):
        repo = Path(directory)
        training = repo / "frozen-training"
        checkpoints = []
        for config in contract.CONFIGURATIONS:
            for model in contract.MODELS:
                for seed in contract.SEEDS:
                    path = training / f"runs/{config}_{model}_seed{seed}/attempt-001/selected_checkpoint.keras"
                    contract.write_json(path, {"synthetic_fixture": True, "configuration": config, "model": model, "seed": seed})
                    checkpoints.append({"configuration": config, "model": model, "seed": seed,
                                        "path": str(path.relative_to(repo)), "bytes": path.stat().st_size,
                                        "sha256": contract.sha256(path), "status": "PASS", "preserved_original_archive": True})
        gate = {"status": "PASS", "n_selected_original_checkpoints": 18,
                "prospective_specification_freeze": {}, "checkpoints": checkpoints}
        return repo, training, gate

    def checked(self, repo, training, gate):
        with mock.patch.object(contract, "TRAINING_STAGE", training), mock.patch.object(contract, "verify_prospective"), mock.patch.object(contract, "verify_record"), mock.patch.object(contract, "read_json", return_value=gate):
            return contract.verify_test_checkpoints(repo / "new-stage", repo)

    def test_all_18_exact_paths_and_frozen_hashes(self):
        with tempfile.TemporaryDirectory() as directory:
            repo, training, gate = self.fixture(directory)
            result, _ = self.checked(repo, training, gate)
            self.assertEqual(len(result), 18)
            for (config, model, seed), path in result.items():
                self.assertEqual(path, training / f"runs/{config}_{model}_seed{seed}/attempt-001/selected_checkpoint.keras")

    def test_checkpoint_hash_identity_duplicate_and_archive_flag_rejections(self):
        with tempfile.TemporaryDirectory() as directory:
            repo, training, gate = self.fixture(directory)
            cases = []
            changed = copy.deepcopy(gate)
            changed["checkpoints"][0]["sha256"] = "f" * 64
            cases.append(changed)
            changed = copy.deepcopy(gate)
            for key in ("path", "bytes", "sha256"):
                changed["checkpoints"][0][key] = changed["checkpoints"][1][key]
            cases.append(changed)
            changed = copy.deepcopy(gate)
            changed["checkpoints"][1] = copy.deepcopy(changed["checkpoints"][0])
            cases.append(changed)
            changed = copy.deepcopy(gate)
            changed["checkpoints"][0]["preserved_original_archive"] = False
            cases.append(changed)
            changed = copy.deepcopy(gate)
            changed["checkpoints"].pop()
            cases.append(changed)
            for i, bad in enumerate(cases):
                with self.subTest(case=i), self.assertRaises(RuntimeError):
                    self.checked(repo, training, bad)


class StubbedExecutionTests(unittest.TestCase):
    def test_complete_loop_loads_18_original_archives_calls_actual_wrapper_72_times(self):
        class FakeTF:
            float32 = np.float32
            TensorSpec = staticmethod(lambda *a, **k: (a, k))
            function = staticmethod(lambda **k: lambda f: f)

        class Network:
            input_shape, output_shape = (None, 4560, 1), (None, 1)

            def __call__(self, x, training):
                calls.append({"shape": x.shape, "training": training})
                return x[:, 0, 0, None]

        loads, calls, records, tables = [], [], {}, {}

        def load_model(path, compile):
            loads.append((path, compile))
            return Network()

        def record(stage, path):
            return {"path": str(path), "bytes": 1, "sha256": "a" * 64}

        def store_json(path, value):
            records[Path(path).name] = value

        def store_table(path, fields, rows):
            tables[Path(path).name] = list(rows)

        def load_table(path, context):
            rows = panel_rows()
            for row in rows:
                row["model"] = context
            return table(rows)

        checkpoints = {(c, m, s): Path(f"synthetic/{c}_{m}_{s}.keras") for c in contract.CONFIGURATIONS for m in contract.MODELS for s in contract.SEEDS}
        keras = SimpleNamespace(backend=SimpleNamespace(clear_session=mock.Mock()), models=SimpleNamespace(load_model=load_model))
        with tempfile.TemporaryDirectory() as directory, contextlib.ExitStack() as stack:
            stage = Path(directory) / "new-stage"
            stage.mkdir()
            stack.enter_context(mock.patch.object(sys, "argv", ["predict_test.py", "--stage", str(stage), "--repo", directory]))
            for name, value in {"verify_prerequisites": lambda *a: (checkpoints, stage / "gate.json"),
                                "verify_test_checkpoints": mock.Mock(), "FeatureCache": synthetic_cache,
                                "initialize_runtime": lambda *a: (FakeTF(), keras), "runtime_environment": lambda *a: {"synthetic": True},
                                "load_test_table": load_table, "file_record": record,
                                "external_checkpoint_record": record, "write_json": store_json, "write_table": store_table}.items():
                stack.enter_context(mock.patch.object(inference, name, value))
            with contextlib.redirect_stdout(io.StringIO()):
                inference.main()
        self.assertEqual(len(loads), 18)
        self.assertEqual({path for path, _ in loads}, set(checkpoints.values()))
        self.assertTrue(all(compile is False for _, compile in loads))
        self.assertEqual(len(calls), 18 * 4)
        self.assertTrue(all(row["training"] is False for row in calls))
        self.assertEqual(keras.backend.clear_session.call_count, 18)
        audit = records["real_network_invariance.json"]
        self.assertEqual(audit["status"], "PASS")
        self.assertEqual(len(audit["audit_rows"]), 24)
        self.assertEqual((audit["n_seed_audits"], audit["n_ensemble_audits"]), (18, 6))
        self.assertEqual(len(tables), 4)
        for name, rows in tables.items():
            self.assertEqual(len(rows), 4)
            self.assertEqual([r["chrom"] for r in rows], ["chr8", "chr8", "chr9", "chr9"])
        self.assertNotIn("failure.json", records)


if __name__ == "__main__":
    unittest.main()
