#!/usr/bin/env python3
"""Synthetic path/gate-equivalence fixtures only; never run the real adapter."""
import csv
import hashlib
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import prepare_calibration_path_adapter as adapter


class AdapterTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="calibration-path-adapter-test-")
        self.repo = Path(self.temporary.name) / "repo"
        self.stage = self.repo / "diseases/COPD/07_gap_closure/internal-training-1.0"
        self.directory = self.stage / "results/chr7_evaluation"
        self.directory.mkdir(parents=True)
        self.input = self.stage / "provenance/evaluation_input.json"
        self.input.parent.mkdir()
        adapter.write_json(self.input, {"synthetic_fixture_only": True})
        decision = {"both_C_contexts_pass": True, "input_manifest_sha256": adapter.sha256(self.input),
                    "model_contexts": {m: {"status": "PASS", "failed_or_inconclusive_gates": []} for m in adapter.MODELS}}
        adapter.write_json(self.directory / "adequacy_decision.json", decision)
        simple = ("frozen_construction", "real_network_invariance", "finite_valid_probabilities")
        numeric = {"all_three_C_seeds_completed": 3, "selection_positives": 100, "selection_controls": 200,
                   "selection_components": 30, "bootstrap_valid_replicates": 2000, "bootstrap_invalid_fraction": 0,
                   "lower95_AUROC": .6, "lower95_AP_gain": .1, "lower95_BrierSkill": .1,
                   "seed_AP_sample_SD": .01, "seed_AP_range": .02}
        rows = [{"model": m, "gate": name, "status": "PASS", "observed": "PASS"} for m in adapter.MODELS for name in simple]
        rows += [{"model": m, "gate": name, "status": "PASS", "observed": str(value)} for m in adapter.MODELS for name, value in numeric.items()]
        with (self.directory / "C_absolute_adequacy.tsv").open("x", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=("model", "gate", "status", "observed"), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        bootstrap = {m: {"status": "PASS", "seed": 314159, "valid": 2000, "required_valid": 2000,
                         "max_attempts": 20000, "attempted": 2000, "invalid": 0, "invalid_fraction": 0} for m in adapter.MODELS}
        adapter.write_json(self.directory / "bootstrap_audit.json", bootstrap)
        self.files = sorted(self.directory.iterdir())

    def tearDown(self):
        self.temporary.cleanup()

    def records(self, spelling):
        result = []
        for path in self.files:
            row = adapter.absolute_record(path)
            if spelling == "stage":
                row["path"] = str(path.relative_to(self.stage))
            elif spelling == "repo":
                row["path"] = str(path.relative_to(self.repo))
            result.append(row)
        return result

    def test_absolute_stage_repo_metadata_are_gate_equivalent(self):
        outputs = []
        before = {p: p.read_bytes() for p in self.files}
        for spelling in ("absolute", "stage", "repo"):
            outputs.append(adapter.verify_passed_evaluation_artifacts(
                self.stage, self.repo, {"artifacts": self.records(spelling)}, self.input)[:3])
        self.assertEqual(outputs[0], outputs[1])
        self.assertEqual(outputs[1], outputs[2])
        self.assertTrue(all(p.read_bytes() == content for p, content in before.items()))
        self.assertFalse((self.stage / "provenance/C_release_precalibration_freeze.json").exists())
        self.assertFalse((self.stage / "provenance/calibration_input.json").exists())

    def test_alias_duplicates_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            adapter.normalize_evaluation_artifacts(self.stage, self.repo, self.records("repo") + [self.records("absolute")[0]])

    def test_hash_and_size_drift_rejected(self):
        for key, value in (("sha256", "0" * 64), ("bytes", 999999)):
            records = self.records("repo")
            records[0][key] = value
            with self.assertRaisesRegex(ValueError, "bytes/hash mismatch"):
                adapter.normalize_evaluation_artifacts(self.stage, self.repo, records)

    def test_missing_file_and_missing_byte_binding_rejected(self):
        records = self.records("repo")
        records[0]["path"] = str((self.directory / "missing.json").relative_to(self.repo))
        with self.assertRaisesRegex(ValueError, "Missing"):
            adapter.normalize_evaluation_artifacts(self.stage, self.repo, records)
        records = self.records("repo")
        records[0].pop("bytes")
        with self.assertRaisesRegex(ValueError, "both SHA-256 and byte count"):
            adapter.normalize_evaluation_artifacts(self.stage, self.repo, records)

    def test_parent_escape_subdirectory_and_symlink_escape_rejected(self):
        other = self.stage / "results/other.json"
        adapter.write_json(other, {"outside": True})
        nested = self.directory / "nested"
        nested.mkdir()
        nested_path = nested / "artifact.json"
        adapter.write_json(nested_path, {"nested": True})
        link = self.directory / "escaped_link.json"
        link.symlink_to(other)
        for path in (other, nested_path, link):
            record = {"path": str(path), "sha256": adapter.sha256(path), "bytes": path.stat().st_size}
            with self.assertRaisesRegex(ValueError, "exact evaluation directory"):
                adapter.normalize_evaluation_artifacts(self.stage, self.repo, [record])

    def test_traversal_rejected_even_if_final_target_would_be_inside(self):
        records = self.records("stage")
        records[0]["path"] = "results/chr7_evaluation/../chr7_evaluation/" + self.files[0].name
        with self.assertRaisesRegex(ValueError, "Traversal"):
            adapter.normalize_evaluation_artifacts(self.stage, self.repo, records)

    def test_frozen_fail_gate_still_blocks_every_path_spelling(self):
        decision_path = self.directory / "adequacy_decision.json"
        decision = json.loads(decision_path.read_text())
        decision["both_C_contexts_pass"] = False
        with decision_path.open("w") as stream:
            json.dump(decision, stream)
        for spelling in ("absolute", "stage", "repo"):
            with self.assertRaisesRegex(ValueError, "BOTH C contexts must PASS"):
                adapter.verify_passed_evaluation_artifacts(self.stage, self.repo, {"artifacts": self.records(spelling)}, self.input)

    def test_json_reference_prediction_matches_frozen_serializer(self):
        path = self.stage / "provenance/synthetic_serialization.json"
        value = {"z": [1, "alpha"], "a": False}
        predicted = adapter.serialized_record(path, value)
        adapter.write_json(path, value)
        self.assertEqual(predicted, adapter.absolute_record(path))

    def test_command_flag_ambiguity_rejected(self):
        self.assertEqual(adapter.option_value(["script.py", "--output", "one"], "--output"), "one")
        self.assertEqual(adapter.option_value(["script.py", "--output=one"], "--output"), "one")
        for command in (["--output"], ["--output", "one", "--output=two"], ["script.py"]):
            with self.assertRaises(ValueError):
                adapter.option_value(command, "--output")


if __name__ == "__main__":
    unittest.main()
