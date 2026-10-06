"""Synthetic administrative seal tests; never operate on the scientific stage."""
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.dont_write_bytecode = True
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import finalize_test_stage as finalizer


class AdministrativeSealTests(unittest.TestCase):
    def test_stage_reference_cannot_escape(self):
        with tempfile.TemporaryDirectory() as directory:
            stage = Path(directory).resolve()
            self.assertEqual(finalizer.local_path(stage, "x.json"), stage / "x.json")
            with self.assertRaises(RuntimeError):
                finalizer.local_path(stage, "../outside.json")

    def test_exclusive_write_and_hash_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "record.json"
            finalizer.write(path, {"synthetic": True})
            with self.assertRaises(FileExistsError):
                finalizer.write(path, {"synthetic": False})
            with self.assertRaises(RuntimeError):
                finalizer.verify(path, {"sha256": "0" * 64})

    def test_unfinished_command_blocks_freeze_but_is_explicit_in_snapshot(self):
        with tempfile.TemporaryDirectory() as directory:
            stage = Path(directory)
            finalizer.write(stage / "provenance/commands/synthetic.started.json", {
                "label": "synthetic", "argv": ["no-model"], "cwd": str(stage), "started_utc": "test", "environment": {}})
            rows, ledger = finalizer.collect_commands(stage)
            self.assertEqual(rows[0]["status"], "RUNNING_AT_PROVENANCE_SNAPSHOT")
            self.assertEqual(ledger["unfinished_commands"], 1)
            with self.assertRaises(RuntimeError):
                finalizer.collect_commands(stage, require_complete=True)

    def test_payload_scope_excludes_only_pycache_and_seal_selfreferences(self):
        with tempfile.TemporaryDirectory() as directory:
            stage = Path(directory)
            for name in ("results/synthetic.json", "scripts/keep.json", "__pycache__/omit.json", "provenance/artifact_checksums.tsv", "provenance/freeze.json"):
                finalizer.write(stage / name, {"synthetic": True})
            self.assertEqual([str(p.relative_to(stage)) for p in finalizer.payload_files(stage)], ["results/synthetic.json", "scripts/keep.json"])

    def test_freeze_ledger_and_seal_are_last_writes_and_hashes_match(self):
        with tempfile.TemporaryDirectory() as directory:
            repo = Path(directory)
            stage = repo / "internal-test-1.0"
            stage.mkdir()
            finalizer.write(stage / "provenance/prospective_specification_freeze.json", {"baseline_commit": "frozen-head"})
            finalizer.write(stage / finalizer.RESULTS / "holdout_decision.json", {
                "version": "internal-test-1.0", "external_access_authorized": False, "candidate_scoring_authorized": False,
                "model_contexts": {"enhancer": {"status": "FAIL"}, "h3k27me3": {"status": "PASS"}},
                "both_C_contexts_pass": False, "next_action": "STOP"})
            finalizer.write(stage / "INTERNAL_TEST_REPORT.md", {"synthetic_report": True})
            script = stage / "scripts/finalize_test_stage.py"
            finalizer.write(script, {"synthetic_script": True})
            writes = []
            original_write, original_tsv = finalizer.write, finalizer.tsv

            def track_write(path, value):
                writes.append(str(Path(path).relative_to(stage)))
                return original_write(path, value)

            def track_tsv(path, rows, fields):
                writes.append(str(Path(path).relative_to(stage)))
                return original_tsv(path, rows, fields)

            rows = [{"label": "synthetic", "returncode": 0, "elapsed_seconds": 0}]
            checks = [{"check": "synthetic", "status": "PASS", "detail": "synthetic only"}]
            with mock.patch.object(finalizer, "collect_commands", return_value=(rows, {})), \
                 mock.patch.object(finalizer, "check_source_chain", return_value={"status": "PASS"}), \
                 mock.patch.object(finalizer, "check_validation_and_registers", return_value=checks), \
                 mock.patch.object(finalizer, "resource_summaries", return_value={}), \
                 mock.patch.object(finalizer.subprocess, "check_output", return_value="frozen-head\n"), \
                 mock.patch.object(finalizer, "__file__", str(script)), \
                 mock.patch.object(finalizer, "write", side_effect=track_write), \
                 mock.patch.object(finalizer, "tsv", side_effect=track_tsv):
                result = finalizer.freeze(stage, repo)
            self.assertEqual(writes[-2:], ["provenance/artifact_checksums.tsv", "provenance/freeze.json"])
            self.assertLess(writes.index("provenance/finalization_record.json"), writes.index("provenance/artifact_checksums.tsv"))
            with (stage / "provenance/artifact_checksums.tsv").open() as stream:
                ledger = list(csv.DictReader(stream, delimiter="\t"))
            self.assertEqual(result["n_payloads"], len(ledger))
            for row in ledger:
                self.assertEqual(finalizer.sha(stage / row["path"]), row["sha256"])
                self.assertEqual((stage / row["path"]).stat().st_size, int(row["bytes"]))
            frozen = finalizer.read(stage / "provenance/freeze.json")
            self.assertFalse(frozen["both_C_contexts_pass"])
            self.assertEqual(frozen["next_action"], "STOP")
            self.assertEqual(frozen["artifact_checksums"]["sha256"], finalizer.sha(stage / "provenance/artifact_checksums.tsv"))
            with self.assertRaises(RuntimeError):
                finalizer.freeze(stage, repo)


if __name__ == "__main__":
    unittest.main()
