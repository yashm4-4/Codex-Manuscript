"""Decode only ordered, explicitly matched diagnostic submatrices with Hail."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
import pandas as pd
import hail as hl
from hail.linalg import BlockMatrix

s=Path(__file__).resolve().parents[1];track=sys.argv[1]
hl.init(master='local[4]',log=str(s/'provenance'/f'hail_matrices_{track}.log'),quiet=True,
        tmp_dir=str(s/'ld/runtime_tmp'),
        spark_conf={'spark.driver.host':'127.0.0.1','spark.driver.bindAddress':'127.0.0.1',
                    'spark.driver.memory':'16g','spark.executor.memory':'8g'})
bm=BlockMatrix.read(str(s/'ld/source/UKBB.EUR.ldadj.bm'))
summaries=json.loads((s/'ld'/track/'gwas_ld_overlap_summary.json').read_text())
for l in summaries:
    lid=l['locus_id'];out=s/'ld'/track/lid
    data=pd.read_csv(out/'ordered_diagnostic_inputs.tsv.gz',sep='\t')
    indices=data.idx.astype('int64').tolist()
    if not indices:continue
    assert indices==sorted(set(indices))
    path=out/'source_upper_triangle.npy'
    if path.exists():raise RuntimeError('Refusing to overwrite prior matrix '+str(path))
    print(f'{track}/{lid}: decoding {len(indices)} x {len(indices)} original signed matrix',flush=True)
    a=bm.filter(indices,indices).to_numpy()
    assert a.shape==(len(indices),len(indices))
    np.save(path,a,allow_pickle=False)
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    (out/'matrix_extraction.json').write_text(json.dumps({'locus_id':lid,'track':track,'variants':len(indices),
        'ordered_inputs_sha256':hashlib.sha256((out/'ordered_diagnostic_inputs.tsv.gz').read_bytes()).hexdigest(),
        'matrix_sha256':h.hexdigest(),'shape':list(a.shape),'dtype':str(a.dtype),
        'original_upper_triangle':True,'ld_recomputed_from_genotypes':False,'posterior_inference':False},indent=2)+'\n')
    del a
    print(f'{track}/{lid}: saved raw matrix',flush=True)
hl.stop()
