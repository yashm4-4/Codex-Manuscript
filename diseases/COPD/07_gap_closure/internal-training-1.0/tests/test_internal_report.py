#!/usr/bin/env python3
"""Formatting/gate tests only; no real model output or metadata is consumed."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
PRELOADED_MODULES = set(sys.modules)
import write_internal_report as report
REPORT_IMPORTED_MODULES = set(sys.modules) - PRELOADED_MODULES


class ReportTests(unittest.TestCase):
    def test_markdown_escaping_and_missing_values(self):
        rendered = report.table(["A", "B"], [("x|y", None), ("new\nline", .5)])
        self.assertIn("x\\|y", rendered)
        self.assertIn("—", rendered)
        self.assertIn("new line", rendered)
        self.assertTrue(rendered.endswith("\n"))

    def test_display_precision_does_not_change_exact_threshold_strings(self):
        self.assertEqual(report.number(.123456789, 7), "0.1234568")
        exact = "0.12345678901234568"
        rendered = report.table(["exact threshold"], [(exact,)])
        self.assertIn(exact, rendered)

    def test_failed_c_stops_no_fallback(self):
        text = report.derive_disposition({"both_C_contexts_pass": False}, None, True)
        self.assertIn("STOP", text)
        self.assertIn("No calibration", text)
        self.assertIn("no A/B fallback", text)

    def test_passed_adequacy_without_calibration_not_released(self):
        text = report.derive_disposition({"both_C_contexts_pass": True}, None, True)
        self.assertIn("release remains blocked", text)

    def test_calibration_without_final_qc_not_released(self):
        text = report.derive_disposition({"both_C_contexts_pass": True}, {"status": "PASS"}, False)
        self.assertIn("QC remains incomplete", text)

    def test_pass_is_review_only_not_test_authorization(self):
        text = report.derive_disposition({"both_C_contexts_pass": True}, {"status": "PASS"}, True)
        self.assertIn("READY FOR INVESTIGATOR REVIEW", text)
        self.assertIn("Separate test/external authorization remains required", text)

    def test_calibration_qc_not_applicable_when_c_does_not_pass(self):
        other_phases = ("inputs", "prepared", "cache", "runs", "packaging", "predictions", "evaluation", "preservation")
        latest = {phase: {"status": "PASS"} for phase in other_phases}
        statuses, complete = report.qc_phase_statuses(False, latest)
        self.assertEqual(statuses["calibration"], "NOT_APPLICABLE (scientifically gated off)")
        self.assertTrue(complete)
        self.assertIn("STOP", report.derive_disposition({"both_C_contexts_pass": False}, None, complete))

    def test_calibration_qc_required_when_both_c_pass(self):
        other_phases = ("inputs", "prepared", "cache", "runs", "packaging", "predictions", "evaluation", "preservation")
        latest = {phase: {"status": "PASS"} for phase in other_phases}
        statuses, complete = report.qc_phase_statuses(True, latest)
        self.assertEqual(statuses["calibration"], "NOT PRESENT")
        self.assertFalse(complete)
        latest["calibration"] = {"status": "PASS"}
        statuses, complete = report.qc_phase_statuses(True, latest)
        self.assertEqual(statuses["calibration"], "PASS")
        self.assertTrue(complete)

    def test_gate_prevents_any_pre_freeze_outcome_access(self):
        with tempfile.TemporaryDirectory(prefix="report-gate-test-") as directory:
            with mock.patch.object(report, "read_json", return_value={"status": "NOT_READY", "checkpoints": []}) as reader:
                with self.assertRaisesRegex(RuntimeError, "before the complete all-18 checkpoint freeze"):
                    report.report(Path(directory))
                self.assertEqual(reader.call_count, 1)
                self.assertEqual(reader.call_args.args[0].name, "checkpoint_freeze.json")

    def test_no_scientific_or_model_runtime_imported(self):
        for module in ("tensorflow", "keras", "numpy", "sklearn"):
            self.assertNotIn(module, REPORT_IMPORTED_MODULES)

    def test_packaging_qualification_preserves_inference_only_meaning(self):
        source = (SCRIPTS / "write_internal_report.py").read_text()
        self.assertIn("not literally uncompiled or optimizer-free", source)
        self.assertIn("Actual prediction loads use `compile=False`", source)
        self.assertIn("must not be used to resume training", source)
        self.assertIn("checkpoints and scientific procedure were not rewritten or rerun", source)


if __name__ == "__main__":
    unittest.main()
