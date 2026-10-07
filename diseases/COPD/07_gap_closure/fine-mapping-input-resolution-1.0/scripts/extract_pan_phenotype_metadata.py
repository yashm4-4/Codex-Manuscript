"""Extract only the authorized J44 EUR metadata from the provider column container."""
from pathlib import Path
import json
import hail as hl

s = Path(__file__).resolve().parents[1]
hl.init(master='local[2]', log=str(s/'provenance/hail_phenotype_metadata_retry.log'), quiet=True,
        spark_conf={'spark.driver.host':'127.0.0.1','spark.driver.bindAddress':'127.0.0.1'})
ht = hl.read_table(str(s/'sources/results_full.mt/cols'))
ht = ht.filter((ht.trait_type == 'icd10') & (ht.phenocode == 'J44') & (ht.pheno_sex == 'both_sexes'))
ht = ht.annotate(pheno_data=ht.pheno_data.filter(lambda x: x.pop == 'EUR'))
rows = ht.select(value=hl.json(ht.row)).collect()
assert len(rows) == 1
result = json.loads(rows[0].value)
result['provenance_caveat'] = 'Current provider Hail column container last modified 2024; flat-file manifest/checksum separately validated. Only EUR J44 association-specific metadata used.'
(s/'ld/panukb_J44_EUR_exact_phenotype_metadata.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
hl.stop()
