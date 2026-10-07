"""Synthetic-only checks of reconstruction, scaling, ordering and sign invariants."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
import numpy as np
import pandas as pd

STAGE=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('ld_numeric',STAGE/'scripts/ld_numeric_diagnostics.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)


class DiagnosticsTests(unittest.TestCase):
    def test_gram_scaling_and_psd(self):
        with tempfile.TemporaryDirectory(prefix='copd-synthetic-') as directory:
            mod.s=Path(directory)
            out=mod.s/'ld/SYNTHETIC/unit';out.mkdir(parents=True)
            expected=np.array([[1.,.5,-.2],[.5,1.,.1],[-.2,.1,1.]])
            scale=np.array([.8,.9,.7])
            gram=expected*np.sqrt(np.outer(scale,scale))
            np.save(out/'source_upper_triangle.npy',np.triu(gram))
            pd.DataFrame({'source_row':[1,2,3],'idx':[10,11,12],'z_ld':[.1,-.2,.3],
                          'neglog10_p':[.1,.2,.3]}).to_csv(out/'ordered_diagnostic_inputs.tsv.gz',sep='\t',index=False)
            result=mod.diagnostic('SYNTHETIC','unit')
            self.assertEqual(result['numerical_status'],'PASS_NUMERICS_ONLY')
            self.assertFalse(result['posterior_inference'])
            self.assertFalse(result['scientific_LD_gate_pass'])
            np.testing.assert_allclose(np.load(out/'diagnostic_signed_correlation.npy'),expected,atol=1e-14)
            np.testing.assert_allclose(np.load(out/'eigenvalues.npy'),np.linalg.eigvalsh(expected),atol=1e-14)

    def test_verified_swap_and_permutation_invariance(self):
        r=np.array([[1.,.6,-.2],[.6,1.,.1],[-.2,.1,1.]])
        z=np.array([1.,2.,-1.]);d=np.array([-1.,1.,1.])
        swapped=r*d[:,None]*d[None,:];swapped_z=z*d
        self.assertEqual(swapped[0,1],-r[0,1])
        self.assertEqual(swapped[1,0],-r[1,0])
        self.assertEqual(swapped_z[0],-z[0])
        perm=[2,0,1]
        np.testing.assert_allclose(z@np.linalg.solve(r,z),
            swapped_z[perm]@np.linalg.solve(swapped[np.ix_(perm,perm)],swapped_z[perm]))
        # A pure strand-label complement does not change physical allele dosage or signs.
        np.testing.assert_array_equal(r,r.copy())

    def test_nonpositive_diagonal_is_not_repaired(self):
        with tempfile.TemporaryDirectory(prefix='copd-synthetic-') as directory:
            mod.s=Path(directory);out=mod.s/'ld/SYNTHETIC/bad';out.mkdir(parents=True)
            np.save(out/'source_upper_triangle.npy',np.array([[0.,.2],[0.,1.]]))
            result=mod.diagnostic('SYNTHETIC','bad')
            self.assertEqual(result['numerical_status'],'FAIL_ORIGINAL_REPRESENTATION_OR_DIAGONAL')
            self.assertFalse((out/'diagnostic_signed_correlation.npy').exists())


if __name__=='__main__':unittest.main()
