"""Bounded engineering-only comparison; deterministic synthetic matrix, no GWAS/LD."""
import contextlib
import datetime
import hashlib
import io
import json
import os
from pathlib import Path
import sys
import time

for name in ['MKL_NUM_THREADS', 'OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS']:
    assert os.environ.get(name) == '2', (name, os.environ.get(name))

import numpy as np
import scipy
from scipy.linalg import eigh

n = 1200
bits = np.random.PCG64(161803398).random_raw(n * n)
a = ((bits >> 11).astype(np.float64) * 2.0**-53 - 0.5).reshape(n, n)
a = (a + a.T) * 0.5
assert np.array_equal(a, a.T)
matrix_sha = hashlib.sha256(a.tobytes(order='C')).hexdigest()
work = np.array(a, dtype=np.float64, order='F', copy=True)
start_cpu, start_wall = time.process_time(), time.perf_counter()
values, vectors = eigh(work, overwrite_a=True, check_finite=False, driver='evr')
elapsed_wall, elapsed_cpu = time.perf_counter()-start_wall, time.process_time()-start_cpu
# Same dimension/spectral scaling as the frozen LD diagnostic tolerance.
tolerance = 100 * np.finfo(np.float64).eps * n * max(1., float(np.max(np.abs(values))))
# Eight deterministic eigenpairs verify the equation without another eigensolver.
indices = np.linspace(0, n-1, 8, dtype=int)
residual = a @ vectors[:, indices] - vectors[:, indices] * values[indices]
residual_max = float(np.max(np.abs(residual)))
config = io.StringIO()
with contextlib.redirect_stdout(config):
    np.show_config()
    scipy.show_config()
libraries = sorted({line.split()[-1] for line in Path('/proc/self/maps').read_text().splitlines()
                    if any(k in line.lower() for k in ['openblas', 'mkl', 'libblas', 'lapack', 'libiomp'])})
result = {
    'scope': 'ENGINEERING_SYNTHETIC_ONLY_NO_GWAS_OR_LD_INPUTS',
    'utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'executable': sys.executable, 'python': sys.version,
    'numpy': np.__version__, 'scipy': scipy.__version__,
    'algorithm': 'scipy.linalg.eigh(float64 F-order, driver=evr, overwrite_a=True, check_finite=False)',
    'dimension': n, 'seed': 161803398, 'generator': 'PCG64.random_raw; top53bit exact dyadic uniforms; symmetric arithmetic average',
    'synthetic_input_sha256_C_order': matrix_sha,
    'eigensolver_wall_seconds': elapsed_wall, 'eigensolver_cpu_seconds': elapsed_cpu,
    'thread_environment': {k: os.environ.get(k) for k in ['MKL_NUM_THREADS','OMP_NUM_THREADS','OPENBLAS_NUM_THREADS']},
    'affinity': sorted(os.sched_getaffinity(0)), 'nice': os.getpriority(os.PRIO_PROCESS, 0),
    'loaded_blas_lapack': libraries, 'build_configuration': config.getvalue(),
    'frozen_style_backward_error_tolerance': tolerance,
    'sampled_eigenpair_indices': indices.tolist(), 'sampled_eigenpair_residual_max_abs': residual_max,
    'sampled_eigenpair_residual_pass': bool(residual_max <= tolerance),
    'eigenvalues': values.tolist(),
    'limitation': 'One small synthetic run under shared-node load; no guarantee of large-matrix speedup or identical real-input rounding. No actual GWAS/LD diagnostics run or changed.'
}
Path(sys.argv[1]).write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps({k:v for k,v in result.items() if k in ['executable','numpy','scipy','dimension','eigensolver_wall_seconds','eigensolver_cpu_seconds','sampled_eigenpair_residual_pass','synthetic_input_sha256_C_order']}))
