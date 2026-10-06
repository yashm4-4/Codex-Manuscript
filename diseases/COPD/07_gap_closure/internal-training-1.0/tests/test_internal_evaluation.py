#!/usr/bin/env python3
"""Synthetic-only metric, bootstrap, firewall, and calibration contract tests."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

import numpy as np
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import evaluate_chr7 as ev
import calibrate_chr7 as cal


class SyntheticMetricTests(unittest.TestCase):
    def test_weighted_tied_metrics_match_sklearn_and_expansion(self):
        rng = np.random.default_rng(9041)
        for _ in range(30):
            y = rng.integers(0, 2, 101)
            p = rng.integers(0, 11, 101) / 10
            w = rng.integers(0, 5, 101)
            actual = ev.PreparedMetrics(y, p).calculate(w)
            self.assertAlmostEqual(actual["AP"], average_precision_score(y, p, sample_weight=w), places=13)
            self.assertAlmostEqual(actual["AUROC"], roc_auc_score(y, p, sample_weight=w), places=13)
            self.assertAlmostEqual(actual["Brier"], brier_score_loss(y, p, sample_weight=w), places=13)
            expanded = ev.PreparedMetrics(np.repeat(y, w), np.repeat(p, w)).calculate()
            for key in ev.METRICS:
                self.assertAlmostEqual(actual[key], expanded[key], places=13)

    def test_constant_predictor_null(self):
        actual = ev.PreparedMetrics([1, 0, 0, 0], [.25] * 4).calculate()
        self.assertEqual(actual["AP"], .25)
        self.assertEqual(actual["AUROC"], .5)
        self.assertEqual(actual["AP_gain"], 0)
        self.assertEqual(actual["BrierSkill"], 0)

    def test_undefined_one_class_draw(self):
        with self.assertRaises(ValueError):
            ev.PreparedMetrics([1, 1], [.4, .5]).calculate()

    def test_invalid_probability_rejected(self):
        for value in (np.nan, np.inf, -1e-9, 1.0001):
            with self.assertRaises(ValueError):
                ev.probability_array([value], "synthetic")

    def test_fixed_float64_seed_ensemble(self):
        f = np.asarray([[.1, .9], [.2, .8], [.3, .7]], dtype=np.float32)
        r = np.asarray([[.2, .8], [.3, .7], [.4, .6]], dtype=np.float32)
        seed, ensemble = ev.symmetric_ensemble(f, r)
        reversed_seed, reversed_ensemble = ev.symmetric_ensemble(r, f)
        self.assertEqual(seed.dtype, np.dtype("float64"))
        np.testing.assert_array_equal(seed, reversed_seed)
        np.testing.assert_array_equal(ensemble, reversed_ensemble)
        np.testing.assert_array_equal(ensemble, np.mean((f.astype(np.float64) + r.astype(np.float64)) / 2, axis=0))

    def test_reliability_bin_boundaries(self):
        bins = ev.reliability([0, 1, 1], [0., .1, 1.])
        self.assertEqual([b["n"] for b in bins], [1, 1, 0, 0, 0, 0, 0, 0, 0, 1])
        self.assertIsNone(bins[2]["mean_probability"])

    def test_bootstrap_pairing_prevalence_and_reproducibility(self):
        y = np.asarray([1, 0, 0] * 30)
        p = np.where(y, .8, .2)
        components = [str(i // 2) for i in range(len(y))]
        predictions = {config: p for config in ev.CONFIGURATIONS}
        rows, audit = ev.component_bootstrap(y, predictions, components, valid_replicates=50, max_attempts=500)
        rows2, audit2 = ev.component_bootstrap(y, predictions, components, valid_replicates=50, max_attempts=500)
        self.assertEqual(rows, rows2)
        self.assertEqual(audit, audit2)
        self.assertEqual(audit["status"], "PASS")
        for index in range(0, len(rows), 3):
            a, b, c = rows[index:index + 3]
            for metric in ev.METRICS:
                self.assertEqual(a[metric], b[metric])
                self.assertEqual(b[metric], c[metric])
            self.assertAlmostEqual(a["AP_gain"], a["AP"] - a["prevalence"])
            self.assertAlmostEqual(a["BrierSkill"], 1 - a["Brier"] / (a["prevalence"] * (1 - a["prevalence"])))
        self.assertGreater(len({row["prevalence"] for row in rows}), 1)

    def test_bootstrap_invalid_fraction_is_a_stop(self):
        scores = {config: np.asarray([.9, .1]) for config in ev.CONFIGURATIONS}
        _, audit = ev.component_bootstrap([1, 0], scores, ["a", "b"], valid_replicates=100, max_attempts=1000)
        self.assertEqual(audit["valid"], 100)
        self.assertGreaterEqual(audit["invalid_fraction"], .1)
        self.assertEqual(audit["status"], "INCONCLUSIVE")

    def test_exact_seed_SD_and_strict_null_gates(self):
        y = np.asarray([1] * 100 + [0] * 200)
        intervals = {"AUROC": {"lower95": .5}, "AP_gain": {"lower95": 0}, "BrierSkill": {"lower95": 0}}
        rows, status = ev.adequacy_gates("enhancer", y, [str(i) for i in range(300)], [.4, .46, .52], intervals, {"valid": 2000, "invalid_fraction": 0})
        self.assertEqual(status, "FAIL")
        gates = {row["gate"]: row for row in rows}
        self.assertEqual(gates["lower95_AUROC"]["status"], "FAIL")
        self.assertEqual(gates["lower95_AP_gain"]["status"], "FAIL")
        self.assertEqual(gates["lower95_BrierSkill"]["status"], "FAIL")
        self.assertAlmostEqual(gates["seed_AP_sample_SD"]["observed"], .06)
        self.assertEqual(gates["seed_AP_sample_SD"]["status"], "FAIL")

    def test_forbidden_panel_chromosome_rejected(self):
        row = {"interval_id": "synthetic", "chrom": "chr8", "partition": "test", "validation_role": "selection", "configuration": "V2-C", "model": "enhancer", "label": "0", "component_id": "g1"}
        with self.assertRaises(ValueError):
            ev.validate_common_rows([row], "enhancer", "selection")

    def test_cutpoints_training_C_only(self):
        rows = [{"interval_id": str(i), "chrom": "chr1", "partition": "train", "label": "1", "configuration": "V2-C", "model": "enhancer", **{name: str(i / 10) for name in ev.COVARIATES}} for i in range(10)]
        cuts = ev.frozen_cutpoints(rows, "enhancer")
        np.testing.assert_allclose(cuts["cutpoints"]["gc_fraction"], [.18, .36, .54, .72])
        rows[0]["chrom"] = "chr7"
        with self.assertRaises(ValueError):
            ev.frozen_cutpoints(rows, "enhancer")


class SyntheticCalibrationTests(unittest.TestCase):
    def test_exact_threshold_untied(self):
        values = np.linspace(0, 1, 400)
        record = cal.conservative_threshold(values, [str(i // 2) for i in range(400)])
        self.assertEqual(record["n_called_controls"], 20)
        self.assertEqual(record["observed_calibration_row_FPR"], .05)
        self.assertEqual(record["threshold_float64"], float.fromhex(record["threshold_hex"]))
        self.assertEqual(record["threshold_float64"], float(record["threshold_decimal_17g"]))
        self.assertEqual(record["sorted_boundary_score"], np.sort(values)[::-1][20])

    def test_ties_conservative(self):
        values = np.r_[np.full(30, .9), np.full(370, .1)]
        record = cal.conservative_threshold(values, [str(i) for i in range(400)])
        self.assertEqual(record["n_at_boundary"], 30)
        self.assertEqual(record["n_called_controls"], 0)

    def test_boundary_one_not_clipped(self):
        record = cal.conservative_threshold(np.ones(400), [str(i) for i in range(400)])
        self.assertGreater(record["threshold_float64"], 1)
        self.assertTrue(record["threshold_above_one"])
        self.assertEqual(record["n_called_controls"], 0)

    def test_row_and_component_minima(self):
        with self.assertRaises(ValueError):
            cal.conservative_threshold([.1] * 199, [str(i) for i in range(199)])
        with self.assertRaises(ValueError):
            cal.conservative_threshold([.1] * 400, ["one"] * 400)

    def test_failed_adequacy_never_touches_predictions(self):
        with tempfile.TemporaryDirectory(prefix="copd-synthetic-eval-") as temp:
            path = Path(temp) / "decision.json"
            ev.write_json(path, {"both_C_contexts_pass": False})
            with self.assertRaisesRegex(ValueError, "both C contexts"):
                cal.validate_calibration_prerequisites({"adequacy_decision": {"path": str(path), "sha256": ev.sha256(path)}})


class SyntheticEndToEndTests(unittest.TestCase):
    def test_full_evaluation_then_frozen_negative_only_calibration(self):
        with tempfile.TemporaryDirectory(prefix="copd-synthetic-eval-") as temp:
            root = Path(temp)

            def ref(path):
                return {"path": str(path), "sha256": ev.sha256(path)}

            spec_path = root / "spec.json"
            ev.write_json(spec_path, {"experimental_version": "pretraining-1.1", "seeds": list(ev.SEEDS), "selection": {"bootstrap": {"valid_replicates": 2000, "seed": 314159}}})
            evidence_path = root / "synthetic_evidence.json"
            ev.write_json(evidence_path, {"synthetic_fixture_only": True})
            checkpoints = []
            for config in ev.CONFIGURATIONS:
                for model in ev.MODELS:
                    for seed in ev.SEEDS:
                        path = root / f"{config}_{model}_{seed}.synthetic.json"
                        ev.write_json(path, {"synthetic_fixture_only": True, "seed": seed})
                        checkpoints.append({"configuration": config, "model": model, "seed": seed, "completed": True, **ref(path)})
            manifest = {"specification": ref(spec_path), "prerequisites": {"construction_gates_pass": True, "construction_evidence": ref(evidence_path), "checkpoint_freeze": ref(evidence_path), "invariance": {"status": "PASS", "real_network": True, "all_configurations_all_seeds_and_ensembles": True, "evidence": ref(evidence_path)}, "checkpoints": checkpoints}, "models": {}}
            calibration_refs = {}
            for model in ev.MODELS:
                panel_paths = {}
                prediction_rows = []
                for role in ("selection", "calibration"):
                    rows = []
                    for i in range(400):
                        label = int(i % 4 == 0)
                        identifier = f"synthetic_{role}_{i:04d}"
                        row = {"interval_id": identifier, "chrom": "chr7", "partition": "validation", "validation_role": role, "configuration": "V2-C", "model": model, "label": label, "component_id": f"{role}_g{i // 2:04d}", "gc_fraction": .2 + i / 1000, "atac_signal_percentile_max_train_only": .1 + i / 800, "repeat_fraction_2001": i / 800, "atac_lobes": "left", f"{model}_same_lobe_support_lobes": "left" if label else "", "non_acgt_fraction": 0}
                        rows.append(row)
                        pred = {"interval_id": identifier, "chrom": "chr7"}
                        for config in ev.CONFIGURATIONS:
                            forward = np.asarray([[.8 if label else .2], [.801 if label else .201], [.802 if label else .202]], dtype=np.float64)
                            reverse = forward + .002
                            sym, ensemble = ev.symmetric_ensemble(forward, reverse)
                            for j, seed in enumerate(ev.SEEDS):
                                pred[f"{config}_p_forward_{seed}"] = float(forward[j, 0])
                                pred[f"{config}_p_rc_{seed}"] = float(reverse[j, 0])
                                pred[f"{config}_q_{seed}"] = float(sym[j, 0])
                            pred[f"{config}_ensemble_q"] = float(ensemble[0])
                        prediction_rows.append(pred)
                    path = root / f"{model}_{role}.tsv.gz"
                    ev.write_rows(path, rows)
                    panel_paths[role] = ref(path)
                train_rows = [{"interval_id": f"synthetic_train_{i}", "chrom": "chr1", "partition": "train", "label": 1, "configuration": "V2-C", "model": model, **{field: .1 + i / 200 for field in ev.COVARIATES}} for i in range(100)]
                train_path = root / f"{model}_train_covariates.tsv.gz"
                ev.write_rows(train_path, train_rows)
                pred_path = root / f"{model}_predictions.tsv.gz"
                ev.write_rows(pred_path, prediction_rows)
                manifest["models"][model] = {"selection": panel_paths["selection"], "train_covariates": ref(train_path), "predictions": ref(pred_path)}
                calibration_refs[model] = {"calibration": panel_paths["calibration"], "predictions": ref(pred_path)}
            manifest_path = root / "evaluation_input.json"
            ev.write_json(manifest_path, manifest)
            output = root / "evaluation"
            decision = ev.evaluate(manifest_path, output)
            self.assertTrue(decision["both_C_contexts_pass"])
            audit = ev.read_json(output / "bootstrap_audit.json")
            self.assertEqual(audit["enhancer"]["valid"], 2000)
            self.assertEqual(audit["h3k27me3"]["valid"], 2000)
            with self.assertRaisesRegex(ValueError, "Refusing"):
                ev.evaluate(manifest_path, output)
            release_path = root / "release_freeze.json"
            decision_ref = ref(output / "adequacy_decision.json")
            ev.write_json(release_path, {"status": "C_READY_FOR_CALIBRATION_FROZEN", "adequacy_decision": decision_ref, "checkpoints": [r for r in checkpoints if r["configuration"] == "V2-C"], "ensemble": {"seeds": list(ev.SEEDS), "reduction_dtype": "float64", "seed_aggregation": "equal-weight", "symmetric_inference": "(p_forward+p_reverse_complement)/2"}, "symmetric_inference_implementation": ref(evidence_path)})
            calibration_path = root / "calibration_input.json"
            ev.write_json(calibration_path, {"specification": ref(spec_path), "adequacy_decision": decision_ref, "evaluation_input": ref(manifest_path), "release_freeze": ref(release_path), "models": calibration_refs})
            thresholds = cal.calibrate(calibration_path, root / "calibration")
            for model in ev.MODELS:
                self.assertEqual(thresholds[model]["n_controls"], 300)
                self.assertEqual(thresholds[model]["n_called_controls"], 0)
            negatives = list(ev.read_rows(root / "calibration/calibration_negative_scores.tsv.gz"))
            self.assertEqual(len(negatives), 600)
            self.assertTrue(all(row["label"] == "0" for row in negatives))


if __name__ == "__main__":
    unittest.main()
