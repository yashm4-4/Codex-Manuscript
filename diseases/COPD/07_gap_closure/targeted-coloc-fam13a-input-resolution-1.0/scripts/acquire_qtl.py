#!/usr/bin/env python3
"""Retrieve only FAM13A records from a broad chr4 tabix slice; no inference."""
import csv
import gzip
import hashlib
import json
from pathlib import Path
from datetime import datetime, timezone
import requests
import pysam

STAGE = Path(__file__).resolve().parents[1]
SOURCE = STAGE / 'sources'
SOURCE.mkdir(exist_ok=True)
URL = 'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000271/QTD000271.all.tsv.gz'
CS_URL = 'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/susie/QTS000015/QTD000271/QTD000271.credible_sets.tsv.gz'
PERM_URL = 'https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/QTD000271/QTD000271.permuted.tsv.gz'
META_URL = 'https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/data_tables/dataset_metadata_r7.tsv'
FIELDS = 'molecular_trait_id chromosome position ref alt variant ma_samples maf pvalue beta se type ac an r2 molecular_trait_object_id gene_id median_tpm rsid'.split()
REGION = ('4', 87000000, 94000000)
GENE = 'ENSG00000138640'

def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

meta=SOURCE/'dataset_metadata_r7.tsv'
if not meta.exists():meta.write_bytes(requests.get(META_URL,timeout=30).content)
cs=SOURCE/'QTD000271.credible_sets.tsv.gz'
if not cs.exists():cs.write_bytes(requests.get(CS_URL,timeout=60).content)
perm=SOURCE/'QTD000271.permuted.tsv.gz'
if not perm.exists():perm.write_bytes(requests.get(PERM_URL,timeout=60).content)
index=SOURCE/'QTD000271.all.tsv.gz.tbi'
if not index.exists():index.write_bytes(requests.get(URL+'.tbi',timeout=60).content)
out=SOURCE/'fam13a_full_cis_source.tsv.gz'
if not out.exists():
    tb=pysam.TabixFile(URL)
    count=0
    with gzip.open(out,'wt',newline='') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(FIELDS)
        for line in tb.fetch(*REGION):
            fields=line.rstrip('\r\n').split('\t')
            if fields[0]==GENE and fields[16]==GENE:
                assert len(fields)==len(FIELDS)
                w.writerow(fields);count+=1
    print('FAM13A rows:',count,flush=True)
receipt={
 'retrieved_utc':datetime.now(timezone.utc).isoformat(),
 'source_url':URL,'source_study':'QTS000015','dataset':'QTD000271',
 'source_release':'2023-04-06 archived file, eQTL Catalogue r6/r7 GTEx v8 uniformly reprocessed',
 'query':'4:87000001-94000000 (pysam zero-based half-open 87000000-94000000); filter both gene columns ENSG00000138640',
 'full_file_downloaded':False,
 'source_head':dict(requests.head(URL,timeout=30).headers),
 'artifacts':{p.name:{'bytes':p.stat().st_size,'sha256':sha(p)} for p in [meta,cs,perm,index,out]},
 'metadata_url':META_URL,'credible_sets_url':CS_URL,'permuted_url':PERM_URL,
 'columns':FIELDS}
(STAGE/'acquisition_provenance.json').write_text(json.dumps(receipt,indent=2)+'\n')
