#!/usr/bin/env python3
"""Synthetic collector tests: never call sacct or read actual stage outputs."""
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import collect_execution_provenance as collector


def synthetic_rows():
    rows = []
    for index in range(18):
        base = dict(zip(collector.FIELDS, ["32021359_" + str(index), "COMPLETED", "0:0",
                        "2026-10-06T10:00:00", "2026-10-06T10:10:00", "600", "00:02:00",
                        "8", "48Gn", "", "cpu=8,gres/gpu=1,gres/gpu:a100=1,mem=48G,node=1", "synthetic-node"]))
        rows.append(base)
        for step in ("batch", "extern"):
            rows.append({**base, "JobID": base["JobID"] + "." + step, "MaxRSS": "123456K"})
    return rows


class ProvenanceTests(unittest.TestCase):
    def test_steps_not_double_counted(self):
        summary = collector.allocation_summary(synthetic_rows(), "32021359")
        self.assertEqual(summary["n_task_allocations"], 18)
        self.assertEqual(summary["n_job_step_rows"], 36)
        self.assertEqual(summary["n_batch_step_rows"], 18)
        self.assertEqual(summary["n_extern_step_rows"], 18)
        self.assertEqual(summary["allocated_GPU_seconds"], 18 * 600)
        self.assertEqual(summary["allocated_GPU_hours"], 3)
        self.assertEqual(summary["allocated_CPU_core_hours"], 24)
        self.assertIsNone(summary["measured_GPU_utilization"])

    def test_frozen_array_mapping(self):
        summary = collector.allocation_summary(synthetic_rows(), "32021359")
        runs = [row["run_id"] for row in summary["task_allocations"]]
        expected = [f"{config}_{model}_seed{seed}" for config in collector.CONFIGURATIONS for model in collector.MODELS for seed in collector.SEEDS]
        self.assertEqual(runs, expected)

    def test_failure_missing_duplicate_or_outside_array_rejected(self):
        for mutate in (lambda rows: rows.pop(0),
                       lambda rows: rows.append(rows[0]),
                       lambda rows: rows[0].update(State="FAILED"),
                       lambda rows: rows[0].update(ExitCode="1:0"),
                       lambda rows: rows[0].update(JobID="32021359_18")):
            rows = synthetic_rows()
            mutate(rows)
            with self.assertRaises(RuntimeError):
                collector.allocation_summary(rows, "32021359")

    def test_gpu_tres_generic_and_typed_not_double_counted(self):
        self.assertEqual(collector.gpu_count("gres/gpu=1,gres/gpu:a100=1"), 1)
        self.assertEqual(collector.gpu_count("gres/gpu:a100=1"), 1)
        for value in ("cpu=8", "gres/gpu=2,gres/gpu:a100=2", "gres/gpu=1,gres/gpu:a100=2"):
            with self.assertRaises(RuntimeError):
                collector.gpu_count(value)

    def test_pipe_parse_preserves_all_requested_fields(self):
        rows = synthetic_rows()
        raw = ("\n".join("|".join(row[field] for field in collector.FIELDS) for row in rows) + "\n").encode()
        self.assertEqual(collector.parse_sacct(raw), rows)

    def test_no_failures_is_explicit_header_only_table(self):
        payload = collector.tsv_bytes([], collector.FAILURE_FIELDS).decode()
        self.assertEqual(payload, "\t".join(collector.FAILURE_FIELDS) + "\n")


if __name__ == "__main__":
    unittest.main()
