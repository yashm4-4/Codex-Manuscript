"""CPU-only independent phase-I architecture and nucleotide input contracts."""
import json
from pathlib import Path
import sys
import tempfile
from types import SimpleNamespace
import unittest

import numpy as np

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
sys.path.insert(0, str(SCRIPTS))
from extract_phase_one import build_phase_one
from prepare_internal_inputs import IndexedFasta, NORMALIZE, COMPLEMENT


class FakeLayers:
    def __getattr__(self, name):
        return lambda *args, **kwargs: {"type": name, "args": args, **kwargs}


class PhaseOneInputTests(unittest.TestCase):
    def test_architecture_matches_independent_frozen_audit(self):
        fake = SimpleNamespace(layers=FakeLayers(), constraints=FakeLayers(),
                               initializers=FakeLayers(), Sequential=lambda layers: layers)
        layers = build_phase_one(fake)
        audit_path = SCRIPTS.parents[1] / "provenance/COPD-V2-PREFLIGHT_architecture_compute_audit.json"
        frozen = json.loads(audit_path.read_text())["phase_one"]
        self.assertEqual(layers[0], {"type": "Input", "args": (), "shape": (2001, 4)})
        self.assertEqual([x["type"] for x in layers[1:]], [x["type"] for x in frozen["layers"]])
        for actual, expected in zip(layers[1:], frozen["layers"]):
            kind = expected["type"]
            if kind == "Conv1D":
                self.assertEqual(actual["args"], (expected["filters"], expected["kernel_size"]))
                self.assertEqual(actual["padding"], expected["padding"])
                self.assertEqual(actual["activation"], expected["activation"])
                self.assertEqual(actual["kernel_constraint"], {"type": "MaxNorm", "args": (), "max_value": .9, "axis": 0})
            elif kind in ("Dropout", "Dense", "Activation"):
                key = {"Dropout": "rate", "Dense": "units", "Activation": "activation"}[kind]
                self.assertEqual(actual["args"], (expected[key],))
                if kind == "Dense":
                    self.assertEqual(actual["activation"], expected["activation"])
            elif kind == "MaxPooling1D":
                self.assertEqual(actual["pool_size"], expected["pool_size"])

    def test_normalize_ambiguity_matches_all_zero_encoding(self):
        raw = b"ACGTNRYKBDHVSW?"
        normalized = raw.translate(NORMALIZE)
        self.assertEqual(normalized, b"ACGT" + b"N"*(len(raw)-4))
        alphabet = np.frombuffer(b"ACGT", dtype=np.uint8)
        original_onehot = np.frombuffer(raw, dtype=np.uint8)[:,None] == alphabet
        normalized_onehot = np.frombuffer(normalized, dtype=np.uint8)[:,None] == alphabet
        np.testing.assert_array_equal(original_onehot, normalized_onehot)

    def test_nucleotide_axis_channel_reverse_equals_string_RC(self):
        normalized = (b"ACGTNNAGCT"*201)[:2001]
        reverse = normalized.translate(COMPLEMENT)[::-1]
        alphabet = np.frombuffer(b"ACGT", dtype=np.uint8)
        onehot = (np.frombuffer(normalized, dtype=np.uint8)[:,None] == alphabet).astype(np.float32)
        expected = (np.frombuffer(reverse, dtype=np.uint8)[:,None] == alphabet).astype(np.float32)
        np.testing.assert_array_equal(onehot[::-1,::-1], expected)
        self.assertEqual(expected.shape, (2001,4))

    def test_RC_involution_and_canonical_identity(self):
        for seq in (b"ACGTN", b"N"*2001, b"AG"*1000+b"C"):
            rc = seq.translate(COMPLEMENT)[::-1]
            self.assertEqual(rc.translate(COMPLEMENT)[::-1], seq)
            self.assertEqual(min(seq, rc), min(rc, rc.translate(COMPLEMENT)[::-1]))

    def fixture(self, directory):
        path = Path(directory) / "fixture.fa"
        sequence = (b"ACGTN" * 430)[:2101]
        chunks = [sequence[i:i+60] for i in range(0,len(sequence),60)]
        path.write_bytes(b">chr1\n" + b"\n".join(chunks) + b"\n")
        Path(str(path)+".fai").write_text("chr1\t2101\t6\t60\t61\n")
        return IndexedFasta(path), sequence

    def test_indexed_fasta_exact_2001_geometry(self):
        with tempfile.TemporaryDirectory() as directory:
            fasta, sequence = self.fixture(directory)
            try:
                self.assertEqual(fasta.fetch("chr1", 50, 2051), sequence[50:2051])
                self.assertEqual(len(fasta.fetch("chr1", 50, 2051)), 2001)
            finally:
                fasta.handle.close()

    def test_indexed_fasta_test_chromosome_guard(self):
        with tempfile.TemporaryDirectory() as directory:
            fasta, _ = self.fixture(directory)
            try:
                for chrom in ("chr8", "chr9", "chrM"):
                    with self.assertRaisesRegex(ValueError, "Forbidden chromosome"):
                        fasta.fetch(chrom, 0, 2001)
            finally:
                fasta.handle.close()

    def test_indexed_fasta_bounds_guard(self):
        with tempfile.TemporaryDirectory() as directory:
            fasta, _ = self.fixture(directory)
            try:
                for start, end in ((-1,2000),(101,2102),(100,100)):
                    with self.assertRaises(ValueError):
                        fasta.fetch("chr1", start, end)
            finally:
                fasta.handle.close()


if __name__ == "__main__":
    unittest.main()
