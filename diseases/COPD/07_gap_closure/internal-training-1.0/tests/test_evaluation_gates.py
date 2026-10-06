#!/usr/bin/env python3
"""Synthetic-only tests: gates bind original files and never execute a model."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import evaluate_chr7 as ev
import calibrate_chr7 as cal
import prepare_evaluation_gates as gates


def fixture(root):
    stage = root / "stage"
    for directory in ("provenance", "inputs/common", "inputs/validation", "cache", "scripts", "runs", "predictions", "results/chr7_evaluation"):
        (stage / directory).mkdir(parents=True)

    def j(relative, obj):
        path = stage / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        ev.write_json(path, obj)
        return path

    def r(path, base=stage):
        return {"path": str(path.relative_to(base)), "bytes": path.stat().st_size, "sha256": ev.sha256(path)}

    spec = j("synthetic_spec.json", {"experimental_version": "pretraining-1.1", "seeds": list(ev.SEEDS), "selection": {"bootstrap": {"valid_replicates": 2000, "seed": 314159}}})
    construction = j("synthetic_construction.json", {"synthetic_only": True})
    source = j("synthetic_frozen_source.json", {"synthetic_only": True})
    source_record = r(source, root)
    prepared = []
    for model in ev.MODELS:
        for suffix in ("selection", "calibration", "train_C_positive_covariates"):
            path = stage / f"inputs/common/{model}_{suffix}.tsv.gz"
            ev.write_rows(path, [{"interval_id": "synthetic", "chrom": "chr7"}])
            prepared.append(r(path))
        path = stage / f"inputs/validation/{model}_chr7.tsv.gz"
        ev.write_rows(path, [{"interval_id": "synthetic", "chrom": "chr7"}])
        prepared.append(r(path))
    input_manifest = j("inputs/input_manifest.json", {"status": "PASS", "completed": True, "files": prepared, "frozen_sources": [source_record]})
    cache_manifest = j("cache/cache_manifest.json", {"synthetic_only": True})
    scripts = [j(f"scripts/{name}", {"synthetic_script_never_executed": name}) for name in ("predict_chr7.py", "phase_two_contract.py", "evaluate_chr7.py", "calibrate_chr7.py", "prepare_evaluation_gates.py")]
    implementation = j("specification/internal_evaluation_implementation.json", {"synthetic_prospective_conventions": True})
    input_gate = j("provenance/execution_input_gate.json", {"status": "PASS", "input_hashes_verified": True, "training_authorized": True, "files": [r(p) for p in [input_manifest, cache_manifest, *scripts, implementation]]})
    checkpoint_records, audit_rows = [], []
    for config in ev.CONFIGURATIONS:
        for model in ev.MODELS:
            for seed in ev.SEEDS:
                prefix = f"runs/{config}_{model}_{seed}"
                checkpoint = j(prefix + "/checkpoint.json", {"synthetic_weights_never_loaded": True})
                selected = j(prefix + "/selected.json", {"configuration": config, "model": model, "seed": seed, "completed": True, "checkpoint": r(checkpoint)})
                completed = j(prefix + "/completed.json", {"status": "COMPLETED"})
                history = j(prefix + "/history.json", {"synthetic_history_never_evaluated": True})
                contract = j(prefix + "/contract.json", {"input_gate": r(input_gate)})
                checkpoint_records.append({"configuration": config, "model": model, "seed": seed, "completed": True, **r(checkpoint), "selected_checkpoint_record": r(selected), "completed_record": r(completed), "history": r(history), "training_contract": r(contract)})
                audit_rows.append({"configuration": config, "model": model, "unit": "seed", "seed": seed, "checkpoint": r(checkpoint), "status": "PASS", "real_network": True, "n_failed": 0, "n_sequences": 1, "atol": 1e-6, "rtol": 1e-6, "max_absolute_difference": 0., "max_tolerance_normalized_difference": 0.})
            audit_rows.append({"configuration": config, "model": model, "unit": "ensemble", "seed": "all3", "status": "PASS", "real_network": True, "n_failed": 0, "n_sequences": 1, "atol": 1e-6, "rtol": 1e-6, "max_absolute_difference": 0., "max_tolerance_normalized_difference": 0., "fixed_seed_order": list(ev.SEEDS), "reduction_dtype": "float64"})
    freeze = j("provenance/checkpoint_freeze.json", {"status": "PASS", "all_18_selected_checkpoints_frozen": True, "selection_prediction_authorized": True, "checkpoints": checkpoint_records, "execution_input_gate": r(input_gate), "input_manifest": r(input_manifest)})
    predictions = [j(f"predictions/{model}_{suffix}.tsv.gz", {"opaque_non_score_fixture": True}) for model in ev.MODELS for suffix in ("chr7", "independent_rc_pass")]
    audit = j("predictions/real_network_invariance.json", {"status": "PASS", "real_network": True, "all_configurations_all_seeds_and_ensembles": True, "n_seed_audits": 18, "n_ensemble_audits": 6, "audit_rows": audit_rows, "checkpoint_freeze": r(freeze), "cache_manifest": r(cache_manifest), "script_files": [r(p) for p in scripts[:2]], "frozen_input_tables": [r(stage / f"inputs/validation/{model}_chr7.tsv.gz") for model in ev.MODELS], "selection_metrics_computed": False, "calibration_performed": False, "chr8_chr9_access": False, "external_benchmark_access": False, "files": [r(p) for p in predictions]})
    j("predictions/completed.json", {"status": "COMPLETED", "invariance_status": "PASS", "invariance_record": r(audit), "files": [r(p) for p in predictions]})
    original = (gates.absolute_record(spec), gates.absolute_record(construction), {str(source): source_record}, {"synthetic_anchor": True})
    return stage, original


def adequate_outputs(stage, both_pass=True):
    directory = stage / "results/chr7_evaluation"
    inp = stage / "provenance/evaluation_input.json"
    decision = {"both_C_contexts_pass": both_pass, "model_contexts": {model: {"status": "PASS" if both_pass else "FAIL", "failed_or_inconclusive_gates": [] if both_pass else ["synthetic failure"]} for model in ev.MODELS}, "input_manifest_sha256": ev.sha256(inp)}
    ev.write_json(directory / "adequacy_decision.json", decision)
    observed = {"frozen_construction": "PASS", "real_network_invariance": "PASS", "finite_valid_probabilities": "PASS", "all_three_C_seeds_completed": 3, "selection_positives": 100, "selection_controls": 200, "selection_components": 30, "bootstrap_valid_replicates": 2000, "bootstrap_invalid_fraction": 0., "lower95_AUROC": .7, "lower95_AP_gain": .2, "lower95_BrierSkill": .1, "seed_AP_sample_SD": .01, "seed_AP_range": .02}
    ev.write_rows(directory / "C_absolute_adequacy.tsv", [{"model": model, "gate": name, "status": "PASS", "observed": value} for model in ev.MODELS for name, value in observed.items()])
    ev.write_json(directory / "bootstrap_audit.json", {model: {"status": "PASS", "seed": 314159, "valid": 2000, "required_valid": 2000, "max_attempts": 20000, "attempted": 2000, "invalid": 0, "invalid_fraction": 0.} for model in ev.MODELS})
    artifacts = [gates.absolute_record(path) for path in sorted(directory.iterdir())]
    ev.write_json(directory / "evaluation_manifest.json", {"input_manifest_sha256": ev.sha256(inp), "evaluation_script_sha256": ev.sha256(stage / "scripts/evaluate_chr7.py"), "artifacts": artifacts})


class SyntheticGateTests(unittest.TestCase):
    def test_evaluation_then_calibration_source_binding_chain(self):
        with tempfile.TemporaryDirectory(prefix="copd-synthetic-gates-") as tmp:
            repo = Path(tmp)
            stage, original = fixture(repo)
            with patch.object(gates, "frozen_construction_sources", return_value=original):
                result = gates.prepare_evaluation(stage, repo)
                self.assertEqual(result["status"], "PASS")
                self.assertFalse(result["evaluation_executed"])
                adequate_outputs(stage)
                result = gates.prepare_calibration(stage, repo)
                self.assertFalse(result["calibration_executed"])
                calibration = ev.read_json(stage / "provenance/calibration_input.json")
                release = cal.validate_calibration_prerequisites(calibration)
                self.assertEqual(len(release["checkpoints"]), 6)
                self.assertEqual(len(release["symmetric_inference_supporting_files"]), 1)
                with self.assertRaisesRegex(ValueError, "replace"):
                    gates.prepare_evaluation(stage, repo)

    def test_failed_C_does_not_create_release_or_calibration(self):
        with tempfile.TemporaryDirectory(prefix="copd-synthetic-gates-") as tmp:
            repo = Path(tmp)
            stage, original = fixture(repo)
            with patch.object(gates, "frozen_construction_sources", return_value=original):
                gates.prepare_evaluation(stage, repo)
                adequate_outputs(stage, both_pass=False)
                with self.assertRaisesRegex(ValueError, "BOTH C"):
                    gates.prepare_calibration(stage, repo)
                self.assertFalse((stage / "provenance/calibration_input.json").exists())
                self.assertFalse((stage / "provenance/C_release_precalibration_freeze.json").exists())

    def test_inference_script_change_rejected(self):
        with tempfile.TemporaryDirectory(prefix="copd-synthetic-gates-") as tmp:
            repo = Path(tmp)
            stage, original = fixture(repo)
            (stage / "scripts/predict_chr7.py").write_text("synthetic tamper\n")
            with patch.object(gates, "frozen_construction_sources", return_value=original):
                with self.assertRaisesRegex(ValueError, "hash/size"):
                    gates.prepare_evaluation(stage, repo)
                self.assertFalse((stage / "provenance/evaluation_input.json").exists())

    def test_metadata_change_rejected(self):
        with tempfile.TemporaryDirectory(prefix="copd-synthetic-gates-") as tmp:
            repo = Path(tmp)
            stage, original = fixture(repo)
            (stage / "inputs/common/enhancer_calibration.tsv.gz").write_text("synthetic tamper\n")
            with patch.object(gates, "frozen_construction_sources", return_value=original):
                with self.assertRaisesRegex(ValueError, "hash/size"):
                    gates.prepare_evaluation(stage, repo)

    def test_post_freeze_evaluation_implementation_change_rejected(self):
        with tempfile.TemporaryDirectory(prefix="copd-synthetic-gates-") as tmp:
            repo = Path(tmp)
            stage, original = fixture(repo)
            (stage / "scripts/evaluate_chr7.py").write_text("synthetic post-freeze tamper\n")
            with patch.object(gates, "frozen_construction_sources", return_value=original):
                with self.assertRaisesRegex(ValueError, "hash/size"):
                    gates.prepare_evaluation(stage, repo)

    def test_audit_missing_unit_rejected(self):
        with tempfile.TemporaryDirectory(prefix="copd-synthetic-gates-") as tmp:
            stage, _ = fixture(Path(tmp))
            audit = ev.read_json(stage / "predictions/real_network_invariance.json")
            freeze = ev.read_json(stage / "provenance/checkpoint_freeze.json")
            checkpoints = {(r["configuration"], r["model"], r["seed"]): r for r in freeze["checkpoints"]}
            audit["audit_rows"].pop()
            with self.assertRaisesRegex(ValueError, "Missing"):
                gates.validate_audit_rows(audit, checkpoints, {model: 1 for model in ev.MODELS})

    def test_numeric_adequacy_value_cannot_hide_under_PASS_label(self):
        with tempfile.TemporaryDirectory(prefix="copd-synthetic-gates-") as tmp:
            repo = Path(tmp)
            stage, original = fixture(repo)
            with patch.object(gates, "frozen_construction_sources", return_value=original):
                gates.prepare_evaluation(stage, repo)
            adequate_outputs(stage)
            directory = stage / "results/chr7_evaluation"
            rows = list(ev.read_rows(directory / "C_absolute_adequacy.tsv"))
            next(row for row in rows if row["gate"] == "lower95_AUROC")["observed"] = ".5"
            with self.assertRaisesRegex(ValueError, "numeric value"):
                gates.validate_passed_adequacy(ev.read_json(directory / "adequacy_decision.json"), rows, ev.read_json(directory / "bootstrap_audit.json"))

    def test_path_escape_rejected(self):
        with tempfile.TemporaryDirectory(prefix="copd-synthetic-gates-") as tmp:
            root = Path(tmp)
            (root / "stage").mkdir()
            outside = root / "outside.json"
            ev.write_json(outside, {"synthetic": True})
            with self.assertRaisesRegex(ValueError, "outside"):
                gates.verified_record(root / "stage", {"path": "../outside.json", "sha256": ev.sha256(outside)})


if __name__ == "__main__":
    unittest.main()
