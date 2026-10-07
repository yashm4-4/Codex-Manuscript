"""Run the remaining frozen diagnostics; no inference or eligibility tuning."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
import subprocess
import sys
import time

s = Path(__file__).resolve().parents[1]
jobs = [('C', x) for x in ['C_L004', 'C_L005', 'C_L006', 'C_L007']]
jobs += [('A', x['locus_id']) for x in json.loads((s/'ld/A/gwas_ld_overlap_summary.json').read_text())]


def run(job):
    track, locus = job
    target = s/'ld'/track/locus
    while not (target/'matrix_extraction.json').exists():
        time.sleep(5)
    if (target/'numeric_diagnostics.json').exists():
        raise RuntimeError('Refusing to overwrite completed diagnostic: '+locus)
    env = dict(os.environ, OPENBLAS_NUM_THREADS='16', OMP_NUM_THREADS='16')
    log = s/'provenance'/f'ld_numeric_{locus}.log'
    print('START', track, locus, flush=True)
    with log.open('x') as f:
        p = subprocess.run([sys.executable, str(s/'scripts/ld_numeric_diagnostics.py'), track, locus],
                           stdout=f, stderr=subprocess.STDOUT, env=env)
    print('FINISH', track, locus, p.returncode, flush=True)
    return {'track': track, 'locus_id': locus, 'returncode': p.returncode,
            'log': str(log.relative_to(s)), 'posterior_inference': False}


with ThreadPoolExecutor(max_workers=2) as pool:
    results = [f.result() for f in as_completed([pool.submit(run, j) for j in jobs])]
(s/'provenance/remaining_ld_diagnostic_jobs.json').write_text(json.dumps(results, indent=2)+'\n')
if any(x['returncode'] != 0 for x in results):
    sys.exit(1)
