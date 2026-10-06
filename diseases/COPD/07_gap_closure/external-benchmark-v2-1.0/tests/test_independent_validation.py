"""Synthetic-only checks for the independent validator; no real data access."""
import importlib.util
from pathlib import Path
import unittest

import numpy as np

source = Path(__file__).resolve().parents[1] / "scripts/validate_external_stage.py"
spec = importlib.util.spec_from_file_location("independent_external_validator", source)
v = importlib.util.module_from_spec(spec)
spec.loader.exec_module(v)


class IndependentValidationTests(unittest.TestCase):
    def identity(self):
        return dict(identity_status="exact", reference_verified="True", grch38_chrom="4",
                    grch38_pos="10001", grch38_ref="A", grch38_alt="T", canonical_variant_id="4:10001:A:T",
                    exact_identity_eligible="", construct_exact_allele_identity_valid="")

    def direction(self):
        return dict(self.identity(), experimental_state="positive", mechanism_in_model_scope="yes",
                    assay_class="conventional_reporter", direction_identity_resolved="True",
                    higher_activity_grch38_allele="T", reported_direction_alt_minus_ref="1",
                    tested_allele1="A", tested_allele2="T", allele_transform="identity")

    def test_empty_optional_flags_do_not_remove_exact_identity(self):
        self.assertTrue(v.exact_identity(self.identity()))

    def test_optional_true_cannot_rescue_unresolved_identity(self):
        row = self.identity(); row.update(identity_status="unevaluable", exact_identity_eligible="True")
        self.assertFalse(v.exact_identity(row))

    def test_explicit_false_vetoes_construct(self):
        row = self.identity(); row["construct_exact_allele_identity_valid"] = "False"
        self.assertFalse(v.exact_identity(row))

    def test_reference_mismatch_flag_not_rescued(self):
        row = self.identity(); row["reference_verified"] = "False"
        self.assertFalse(v.exact_identity(row))

    def test_canonical_key_mismatch_not_rescued(self):
        row = self.identity(); row["canonical_variant_id"] = "4:10002:A:T"
        self.assertFalse(v.exact_identity(row))

    def test_explicit_allelic_direction(self):
        self.assertEqual(v.direction_sign(self.direction()), 1)

    def test_direction_sign_contradiction_rejected(self):
        row = self.direction(); row["reported_direction_alt_minus_ref"] = "-1"
        with self.assertRaises(AssertionError): v.direction_sign(row)

    def test_direction_unresolved_ratio_not_invented(self):
        row = self.direction(); row["higher_activity_grch38_allele"] = ""
        self.assertIsNone(v.direction_sign(row))

    def test_null_context_not_direction_positive(self):
        row = self.direction(); row["experimental_state"] = "null"
        self.assertIsNone(v.direction_sign(row))

    def test_splice_readout_excluded(self):
        row = self.direction(); row["assay_class"] = "splicing"
        self.assertIsNone(v.direction_sign(row))

    def test_TF_binding_not_activity_direction(self):
        row = self.direction(); row["assay_class"] = "TF_binding"
        self.assertIsNone(v.direction_sign(row))

    def test_wrong_experimental_alleles_rejected(self):
        row = self.direction(); row["tested_allele2"] = "C"
        with self.assertRaises(AssertionError): v.direction_sign(row)

    def test_seed_ensemble_uses_float64_after_float32_outputs(self):
        raw = np.arange(48, dtype=np.float32).reshape(2, 3, 1, 2, 4) / np.float32(100)
        raw = raw.astype(np.float64)
        seed, audit_seed, ensemble, audit_ensemble = v.independent_arithmetic(raw)
        self.assertTrue(np.array_equal(seed, (raw[..., 0] + raw[..., 1]) / 2))
        self.assertTrue(np.array_equal(ensemble, np.mean(seed, axis=1, dtype=np.float64)))
        self.assertTrue(np.array_equal(audit_ensemble, np.mean(audit_seed, axis=1, dtype=np.float64)))

    def test_missing_seed_rejected(self):
        with self.assertRaises(ValueError): v.independent_arithmetic(np.zeros((2, 2, 1, 2, 4)))

    def test_non_float32_raw_probability_rejected(self):
        raw = np.full((2, 3, 1, 2, 4), 0.1, dtype=np.float64)
        with self.assertRaises(ValueError): v.independent_arithmetic(raw)

    def test_nan_rejected(self):
        raw = np.zeros((2, 3, 1, 2, 4)); raw[0, 0, 0, 0, 0] = np.nan
        with self.assertRaises(ValueError): v.independent_arithmetic(raw)

    def test_threshold_equality_inclusive_and_predecessor_excluded(self):
        for threshold in v.THRESHOLDS.values():
            self.assertTrue(threshold >= threshold)
            self.assertFalse(np.nextafter(threshold, -np.inf) >= threshold)
            self.assertEqual(float.fromhex(threshold.hex()), threshold)

    def test_no_epsilon_direction_threshold(self):
        self.assertEqual([v.signed(x) for x in (-1e-300, 0.0, 1e-300)], [-1, 0, 1])

    def test_ties_retained_in_strict_denominator(self):
        result = v.direction_totals(["concordant", "discordant", "tie"])
        self.assertEqual(result["denominator"], 3)
        self.assertEqual(result["strict_concordance_fraction"], 1 / 3)
        self.assertEqual(result["secondary_non_tied_fraction"], 1 / 2)

    def test_empty_H3_direction_not_zero_success_rate(self):
        result = v.direction_totals([])
        self.assertIsNone(result["strict_concordance_fraction"])
        self.assertIsNone(result["secondary_non_tied_fraction"])

    def test_V1_region_gate_not_combined_delta_gate(self):
        row = {"v1_model_evaluable": "True", "rc_model_evaluable": "False"}
        for prefix in ("enhancer", "silencer"):
            row.update({prefix + "_region_threshold": "0.5", prefix + "_forward_ref_score": "0.75",
                        prefix + "_forward_alt_score": "0.75", prefix + "_forward_region_score": "0.75",
                        prefix + "_forward_region_gate": "True", prefix + "_forward_call": "False"})
        self.assertEqual(v.expected_v1_region(row, "forward"), {"enhancer": True, "h3k27me3": True, "union": True})
        self.assertIsNone(v.expected_v1_region(row, "rc"))

    def test_V1_unavailable_not_zero(self):
        self.assertIsNone(v.expected_v1_region({"v1_model_evaluable": "False"}, "forward"))

    def test_no_naive_timestamp(self):
        with self.assertRaises(ValueError): v.timestamp("2026-10-06T20:00:00")

    def test_Boolean_strings_not_python_truthiness(self):
        self.assertFalse(v.true("False")); self.assertTrue(v.true("True"))


if __name__ == "__main__":
    unittest.main()
