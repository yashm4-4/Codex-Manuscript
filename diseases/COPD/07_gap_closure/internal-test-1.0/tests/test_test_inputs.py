#!/usr/bin/env python3
"""Synthetic-only tests: never load real sequences, panels, or a model runtime."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.dont_write_bytecode = True
SPEC = importlib.util.spec_from_file_location("prepare_test_inputs", Path(__file__).parents[1] / "scripts/prepare_test_inputs.py")
prep = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(prep)


def sequence_row(raw=b"A" * 2001, start=1000, chrom="chr8", label="1", model="enhancer"):
    normalized = raw.translate(prep.NORMALIZE)
    canonical = min(normalized, normalized.translate(prep.COMPLEMENT)[::-1])
    return {"interval_id": f"{chrom}:{start}-{start + 1000}", "chrom": chrom,
            "core_start": str(start), "core_end": str(start + 1000),
            "input_start": str(start - 501), "input_end": str(start + 1500),
            "partition": "test", "validation_role": "test", "component_id": f"{chrom}:{start - 501}",
            "sequence_available": "1", "sequence_length": "2001",
            "sequence_sha256": hashlib.sha256(raw).hexdigest(),
            "canonical_rc_sequence_sha256": hashlib.sha256(canonical).hexdigest(),
            "gc_fraction": str((raw.count(b"G") + raw.count(b"C")) / 2001),
            "non_acgt_fraction": str(normalized.count(b"N") / 2001),
            "configuration": "V2-C", "model": model, "label": label,
            model + "_same_lobe_peak_support": "1" if label == "1" else "0",
            model + "_accessible_control_eligible": "0" if label == "1" else "1",
            "non_acgt_any": str(int(b"N" in normalized)),
            "atac_signal_percentile_max_train_only": "0.1250000000",
            "preserved_empty_field": ""}


def panel():
    rows = [sequence_row(), sequence_row(b"C" * 2001, start=4000, chrom="chr9", label="0")]
    expected = {"rows": 2, "positive": 1, "negative": 1, "components": 2,
                "encoded_identities": 2, "chromosome_counts": {"chr8": 1, "chr9": 1}}
    return list(rows[0]), rows, expected


class GeometryAndSequenceTests(unittest.TestCase):
    def test_exact_2001_geometry_passes(self):
        prep.verify_geometry(sequence_row())

    def test_reject_non_test_chromosome(self):
        with self.assertRaises(RuntimeError):
            prep.verify_geometry(sequence_row(chrom="chr7"))

    def test_reject_wrong_partition_or_role(self):
        for column in ("partition", "validation_role"):
            with self.subTest(column=column):
                row = sequence_row()
                row[column] = "train"
                with self.assertRaises(RuntimeError):
                    prep.verify_geometry(row)

    def test_reject_wrong_geometry(self):
        for column in ("input_start", "input_end", "core_start", "core_end"):
            with self.subTest(column=column):
                row = sequence_row()
                row[column] = str(int(row[column]) + 1)
                with self.assertRaises(RuntimeError):
                    prep.verify_geometry(row)

    def test_reject_unavailable_sequence(self):
        row = sequence_row()
        row["sequence_available"] = "0"
        with self.assertRaises(RuntimeError):
            prep.verify_geometry(row)

    def test_reject_malformed_hash(self):
        row = sequence_row()
        row["sequence_sha256"] = "a" * 63
        with self.assertRaises(RuntimeError):
            prep.verify_geometry(row)

    def test_forward_orientation(self):
        canonical, flag = prep.verify_sequence(sequence_row(), b"A" * 2001)
        self.assertEqual(canonical, b"A" * 2001)
        self.assertEqual(flag, 1)

    def test_reverse_orientation(self):
        canonical, flag = prep.verify_sequence(sequence_row(b"T" * 2001), b"T" * 2001)
        self.assertEqual(canonical, b"A" * 2001)
        self.assertEqual(flag, 0)

    def test_ambiguous_symbols_normalize_before_canonicalization(self):
        raw = b"R" + b"A" * 2000
        canonical, flag = prep.verify_sequence(sequence_row(raw), raw)
        self.assertEqual(canonical, b"N" + b"A" * 2000)
        self.assertEqual(flag, 1)

    def test_reject_wrong_raw_sequence_hash(self):
        with self.assertRaises(RuntimeError):
            prep.verify_sequence(sequence_row(), b"C" * 2001)

    def test_reject_wrong_gc_or_non_acgt(self):
        for column in ("gc_fraction", "non_acgt_fraction"):
            with self.subTest(column=column):
                row = sequence_row()
                row[column] = "0.1"
                with self.assertRaises(RuntimeError):
                    prep.verify_sequence(row, b"A" * 2001)

    def test_reject_wrong_canonical_hash(self):
        row = sequence_row()
        row["canonical_rc_sequence_sha256"] = "b" * 64
        with self.assertRaises(RuntimeError):
            prep.verify_sequence(row, b"A" * 2001)


class FrozenPanelTests(unittest.TestCase):
    def test_panel_exact_metadata_passes(self):
        columns, rows, expected = panel()
        audit = prep.validate_panel("enhancer", columns, rows, columns, copy.deepcopy(rows), expected)
        self.assertEqual(audit["status"], "PASS")
        self.assertFalse(audit["matching_rebuilt"])

    def test_any_frozen_metadata_change_rejected(self):
        for column, replacement in (("label", "0"), ("component_id", "chr8:0"),
                                    ("atac_signal_percentile_max_train_only", "0.125")):
            with self.subTest(column=column):
                columns, rows, expected = panel()
                original = copy.deepcopy(rows)
                rows[0][column] = replacement
                with self.assertRaises(RuntimeError):
                    prep.validate_panel("enhancer", columns, rows, columns, original, expected)

    def test_duplicate_membership_rejected(self):
        columns, rows, expected = panel()
        rows[1] = copy.deepcopy(rows[0])
        with self.assertRaises(RuntimeError):
            prep.validate_panel("enhancer", columns, rows, columns, copy.deepcopy(rows), expected)

    def test_missing_rows_rejected(self):
        columns, rows, expected = panel()
        with self.assertRaises(RuntimeError):
            prep.validate_panel("enhancer", columns, rows[:1], columns, rows, expected)

    def test_wrong_component_count_rejected(self):
        columns, rows, expected = panel()
        expected["components"] = 1
        with self.assertRaises(RuntimeError):
            prep.validate_panel("enhancer", columns, rows, columns, copy.deepcopy(rows), expected)

    def test_invalid_source_eligibility_rejected_even_if_shared(self):
        for row_index, column, value in ((0, "enhancer_same_lobe_peak_support", "0"),
                                         (1, "enhancer_accessible_control_eligible", "0"),
                                         (0, "atac_signal_percentile_max_train_only", "nan"),
                                         (0, "non_acgt_any", "1")):
            with self.subTest(column=column):
                columns, rows, expected = panel()
                rows[row_index][column] = value
                with self.assertRaises(RuntimeError):
                    prep.validate_panel("enhancer", columns, rows, columns, copy.deepcopy(rows), expected)

    def test_digest_preserves_row_order_and_string_encodings(self):
        columns, rows, _ = panel()
        digest = prep.table_digest(columns, rows)
        self.assertNotEqual(digest, prep.table_digest(columns, list(reversed(rows))))
        rows[0]["atac_signal_percentile_max_train_only"] = "0.125"
        self.assertNotEqual(digest, prep.table_digest(columns, rows))

    def test_table_roundtrip_preserves_frozen_strings_and_empty_fields(self):
        columns, rows, _ = panel()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "panel.tsv.gz"
            prep.write_table(path, columns, rows)
            new_columns, new_rows = prep.read_table(path)
            self.assertEqual((columns, rows), (new_columns, new_rows))
            with self.assertRaises(FileExistsError):
                prep.write_table(path, columns, rows)

    def test_deterministic_gzip_has_no_filename_or_timestamp(self):
        columns, rows, _ = panel()
        with tempfile.TemporaryDirectory() as directory:
            first, second = Path(directory) / "first.tsv.gz", Path(directory) / "second.tsv.gz"
            prep.write_table(first, columns, rows)
            prep.write_table(second, columns, rows)
            self.assertEqual(first.read_bytes(), second.read_bytes())

    def test_cross_model_shared_interval_reuses_one_cache_row(self):
        _, rows, _ = panel()
        union = prep.unique_intervals({"enhancer": rows, "h3k27me3": [copy.deepcopy(rows[0])]})
        self.assertEqual(len(union), 2)

    def test_cross_model_component_conflict_is_hard_stop(self):
        _, rows, _ = panel()
        second = copy.deepcopy(rows[0])
        second["component_id"] = "chr8:0"
        with self.assertRaises(RuntimeError):
            prep.unique_intervals({"enhancer": rows, "h3k27me3": [second]})


class FastaAndGateTests(unittest.TestCase):
    def test_reader_exact_panel_allowlist_and_bounds(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "synthetic.fa"
            path.write_bytes(b">chr8\nACGT\nTGCA\nAAAA\n")
            Path(str(path) + ".fai").write_text("chr8\t12\t6\t4\t5\n")
            fasta = prep.TestIndexedFasta(path, [("chr8", 2, 10)])
            try:
                self.assertEqual(fasta.fetch("chr8", 2, 10), b"GTTGCAAA")
                for request in (("chr7", 2, 10), ("chr8", 1, 10), ("chr9", 2, 10)):
                    with self.assertRaises(RuntimeError):
                        fasta.fetch(*request)
                self.assertEqual(fasta.fetches, 1)
            finally:
                fasta.close()

    def test_reader_rejects_broad_or_non_test_allowlist(self):
        with self.assertRaises(RuntimeError):
            prep.TestIndexedFasta(Path("not_opened.fa"), [("chr1", 0, 2001)])
        with self.assertRaises(RuntimeError):
            prep.TestIndexedFasta(Path("not_opened.fa"), [])

    def test_test_only_reader_keeps_inconsistent_rows_for_rejection(self):
        columns, rows, _ = panel()
        rows[0]["partition"] = "train"
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "panel.tsv.gz"
            prep.write_table(path, columns, rows)
            _, selected = prep.read_table(path, test_only=True)
            self.assertEqual(len(selected), 2)
            with self.assertRaises(RuntimeError):
                prep.verify_geometry(selected[0])

    def test_prospective_hash_and_independent_gate_required(self):
        with tempfile.TemporaryDirectory() as directory:
            stage = Path(directory)
            (stage / "provenance").mkdir()
            (stage / "specification").mkdir()
            spec_path = stage / "specification/test_specification.json"
            spec_path.write_text('{"synthetic_only": true}\n')
            integrity_path = stage / "provenance/input_integrity_verified_before.json"
            freeze_path = stage / "provenance/prospective_specification_freeze.json"
            freeze_path.write_text(json.dumps({"status": "PASS", "before_any_test_model_inference": True,
                                              "specification": prep.record(spec_path, stage)}))
            integrity_path.write_text(json.dumps({"status": "PASS", "prospective_specification_freeze": prep.record(freeze_path, stage)}))
            self.assertEqual(len(prep.verify_preparation_gate(stage)), 3)
            spec_path.write_text('{}\n')
            with self.assertRaises(RuntimeError):
                prep.verify_preparation_gate(stage)

    def test_freeze_chain_substitution_rejected_even_at_same_size(self):
        with tempfile.TemporaryDirectory() as directory:
            stage = Path(directory)
            (stage / "provenance").mkdir()
            (stage / "specification").mkdir()
            spec_path = stage / "specification/test_specification.json"
            spec_path.write_text('{"synthetic_only": true}\n')
            freeze_path = stage / "provenance/prospective_specification_freeze.json"
            freeze = {"status": "PASS", "before_any_test_model_inference": True,
                      "specification": prep.record(spec_path, stage), "synthetic_revision": "original"}
            freeze_path.write_text(json.dumps(freeze))
            integrity_path = stage / "provenance/input_integrity_verified_before.json"
            frozen_record = prep.record(freeze_path, stage)
            integrity_path.write_text(json.dumps({"status": "PASS", "prospective_specification_freeze": frozen_record}))
            prep.verify_preparation_gate(stage)
            freeze["synthetic_revision"] = "replaced"
            freeze_path.write_text(json.dumps(freeze))
            self.assertEqual(freeze_path.stat().st_size, frozen_record["bytes"])
            with self.assertRaisesRegex(RuntimeError, "freeze hash/size mismatch"):
                prep.verify_preparation_gate(stage)

    def test_preinference_attestation_required_even_if_hash_chain_valid(self):
        with tempfile.TemporaryDirectory() as directory:
            stage = Path(directory)
            (stage / "provenance").mkdir()
            (stage / "specification").mkdir()
            spec_path = stage / "specification/test_specification.json"
            spec_path.write_text('{"synthetic_only": true}\n')
            freeze_path = stage / "provenance/prospective_specification_freeze.json"
            freeze_path.write_text(json.dumps({"status": "PASS", "before_any_test_model_inference": False,
                                              "specification": prep.record(spec_path, stage)}))
            integrity_path = stage / "provenance/input_integrity_verified_before.json"
            integrity_path.write_text(json.dumps({"status": "PASS", "prospective_specification_freeze": prep.record(freeze_path, stage)}))
            with self.assertRaisesRegex(RuntimeError, "specification freeze is not PASS"):
                prep.verify_preparation_gate(stage)
            integrity_path.write_text('{"status":"FAIL"}\n')
            with self.assertRaises(RuntimeError):
                prep.verify_preparation_gate(stage)


if __name__ == "__main__":
    unittest.main(verbosity=2)
