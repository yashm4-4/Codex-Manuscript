"""Check local Hail compatibility and public LD index globals; no inference."""
from pathlib import Path
import json
import os
import subprocess
import hail as hl

STAGE = Path(__file__).resolve().parents[1]
hl.init(master='local[4]', default_reference='GRCh37',
        log=str(STAGE/'provenance/hail_initialization.log'),
        tmp_dir='/tmp/copd-input-resolution-hail',
        spark_conf={'spark.driver.memory': '12g', 'spark.executor.memory': '8g',
                    'spark.driver.host': '127.0.0.1', 'spark.driver.bindAddress': '127.0.0.1'},
        quiet=True)
ht = hl.read_table(str(STAGE/'ld/source/UKBB.EUR.ldadj.variant.ht'))
result = dict(hl.eval(ht.index_globals()))
result.update(hail_version=hl.__version__, java_home=os.environ.get('JAVA_HOME'),
              table_schema=str(ht.row.dtype), matrix_values_acquired=False)
(STAGE/'ld/panukb_index_globals.json').write_text(json.dumps(result, indent=2)+'\n')
(STAGE/'provenance/diagnostic_python_packages.txt').write_bytes(
    subprocess.check_output([os.sys.executable, '-m', 'pip', 'freeze']))
print(json.dumps(result, indent=2))
hl.stop()
