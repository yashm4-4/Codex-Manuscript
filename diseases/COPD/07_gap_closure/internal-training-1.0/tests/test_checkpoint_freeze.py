#!/usr/bin/env python3
"""Synthetic CPU-only checkpoint history/release tests; no actual fit data."""
import csv
import os
from pathlib import Path
import re
import sys
import subprocess
import tempfile
import unittest

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
from freeze_checkpoints import inspect_history
from phase_two_contract import array_digest, epoch_plan


def history_rows(values, n=19, seed=104729):
    rows = []
    best, wait = float("inf"), 0
    for epoch, value in enumerate(values):
        order, reverse = epoch_plan(seed, epoch, n)
        improved = value < best
        if improved:
            best, wait = value, 0
        else:
            wait += 1
        rows.append({"epoch_zero_based": epoch, "epoch_one_based": epoch + 1,
                     "training_BCE": .7, "checkpoint_symmetric_BCE": value,
                     "n_train": n, "steps": 1, "n_reverse_views": int(reverse.sum()),
                     "permutation_sha256": array_digest(order), "orientation_sha256": array_digest(reverse),
                     "selected_new_best": int(improved), "wait": wait})
    return rows


class FreezeTests(unittest.TestCase):
    def run_case(self, rows, selected):
        with tempfile.TemporaryDirectory(prefix="checkpoint-freeze-test-") as directory:
            path = Path(directory) / "history.tsv"
            with path.open("x", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
                writer.writeheader()
                writer.writerows(rows)
            return inspect_history(path, selected, {"n_train": 19})

    def test_earliest_tie_and_patience(self):
        rows = history_rows([.5] * 16)
        selected = {"seed": 104729, "selected_epoch_one_based": 1, "checkpoint_symmetric_BCE": .5,
                    "restored_checkpoint_symmetric_BCE": .5, "epochs_completed": 16, "early_stopped": True}
        self.assertEqual(len(self.run_case(rows, selected)), 16)
        with self.assertRaises(RuntimeError):
            self.run_case(rows, {**selected, "selected_epoch_one_based": 2})

    def test_fifty_epochs_valid_without_earlystop(self):
        values = [1 - i / 100 for i in range(50)]
        rows = history_rows(values)
        selected = {"seed": 104729, "selected_epoch_one_based": 50, "checkpoint_symmetric_BCE": values[-1],
                    "restored_checkpoint_symmetric_BCE": values[-1], "epochs_completed": 50, "early_stopped": False}
        self.assertEqual(len(self.run_case(rows, selected)), 50)

    def test_rng_drift_rejected(self):
        rows = history_rows([.5] * 16)
        rows[3]["orientation_sha256"] = "0" * 64
        selected = {"seed": 104729, "selected_epoch_one_based": 1, "checkpoint_symmetric_BCE": .5,
                    "restored_checkpoint_symmetric_BCE": .5, "epochs_completed": 16, "early_stopped": True}
        with self.assertRaises(RuntimeError):
            self.run_case(rows, selected)

    def test_incomplete_or_late_stopping_rejected(self):
        selected = {"seed": 104729, "selected_epoch_one_based": 1, "checkpoint_symmetric_BCE": .5,
                    "restored_checkpoint_symmetric_BCE": .5, "epochs_completed": 15, "early_stopped": False}
        with self.assertRaises(RuntimeError):
            self.run_case(history_rows([.5] * 15), selected)
        with self.assertRaises(RuntimeError):
            self.run_case(history_rows([.5] * 17), {**selected, "epochs_completed": 17, "early_stopped": True})

    def test_array_declares_exact_frozen_cartesian_product(self):
        script = (SCRIPTS / "train_array.sh").read_text()
        def array(name):
            return re.search(rf"^{name}=\(([^)]*)\)$", script, re.MULTILINE).group(1).split()
        configs, models, seeds = array("TREDNET_CONFIGURATIONS"), array("TREDNET_MODELS"), array("TREDNET_SEEDS")
        actual = [(configs[i // 6], models[(i % 6) // 3], int(seeds[i % 3])) for i in range(18)]
        expected = [(c, m, s) for c in ("V2-A", "V2-B", "V2-C") for m in ("enhancer", "h3k27me3") for s in (104729, 130363, 155921)]
        self.assertEqual(actual, expected)
        self.assertEqual(len(set(actual)), 18)
        for text in ("#SBATCH --array=0-17%4", "#SBATCH --gres=gpu:a100:1", "#SBATCH --cpus-per-task=8",
                     "#SBATCH --mem=48G", "#SBATCH --time=04:00:00", 'export PYTHONHASHSEED="$TREDNET_SEED"'):
            self.assertIn(text, script)

    def test_no_model_runtime_imported(self):
        self.assertNotIn("tensorflow", sys.modules)
        self.assertNotIn("keras", sys.modules)

    def test_actual_shell_mapping_ignores_submission_or_spool_location(self):
        script = (SCRIPTS / "train_array.sh").read_text()
        prefix = script.split("\nexec bash ", 1)[0]
        prefix += '\nprintf "%s\\t%s\\t%s\\t%s\\t%s\\n" "$TREDNET_STAGE" "$TREDNET_CONFIGURATION" "$TREDNET_MODEL" "$TREDNET_SEED" "$PYTHONHASHSEED"\n'
        expected_stage = str(SCRIPTS.parent)
        with tempfile.TemporaryDirectory(prefix="slurm-spool-simulation-") as directory:
            for index in range(18):
                result = subprocess.run(["bash", "-c", prefix], cwd=directory,
                                        env={**os.environ, "SLURM_ARRAY_TASK_ID": str(index)},
                                        text=True, capture_output=True, check=True)
                stage, configuration, model, seed, pythonhashseed = result.stdout.strip().split("\t")
                self.assertEqual(stage, expected_stage)
                self.assertEqual(configuration, ("V2-A", "V2-B", "V2-C")[index // 6])
                self.assertEqual(model, ("enhancer", "h3k27me3")[(index % 6) // 3])
                self.assertEqual(int(seed), (104729, 130363, 155921)[index % 3])
                self.assertEqual(seed, pythonhashseed)


if __name__ == "__main__":
    unittest.main()
