#!/usr/bin/env python3
"""Reconstruct ignored small public inputs; compare hashes to frozen receipts."""
import csv
import hashlib
import io
import json
import zipfile
from pathlib import Path

import requests
from pypdf import PdfReader

STAGE=Path(__file__).resolve().parents[1]
DEST=STAGE/'sources'
URLS={
 'saferali_fulltext.xml':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12481885/fullTextXML',
 'saferali_supplementary.zip':'https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12481885/supplementaryFiles',
 'ENCFF324XYW.bed.gz':'https://www.encodeproject.org/files/ENCFF324XYW/@@download/ENCFF324XYW.bed.gz',
 'encode_ENCSR528UQX.json':'https://www.encodeproject.org/annotations/ENCSR528UQX/?format=json&frame=embedded',
 'ensembl_rs2013701.json':'https://rest.ensembl.org/variation/human/rs2013701?content-type=application/json',
}
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 DEST.mkdir(exist_ok=True)
 expected={x['source_id']:x['sha256'] for x in csv.DictReader(open(STAGE/'source_acquisition_receipts.tsv'),delimiter='\t')}
 for name,url in URLS.items():
  path=DEST/name
  if not path.exists():
   r=requests.get(url,timeout=120);r.raise_for_status();path.write_bytes(r.content)
  assert digest(path)==expected[name],f'PUBLIC SOURCE HASH CHANGED: {name}'
 with zipfile.ZipFile(DEST/'saferali_supplementary.zip') as z:
  pdf=z.read('mmc1.pdf')
 pdfpath=DEST/'saferali_mmc1.pdf';pdfpath.write_bytes(pdf)
 assert digest(pdfpath)==expected['saferali_mmc1.pdf']
 lines=[]
 for i,page in enumerate(PdfReader(io.BytesIO(pdf)).pages,1):
  lines.append(f'\n===PAGE {i}===\n'+(page.extract_text() or ''))
 (DEST/'saferali_mmc1.txt').write_text(''.join(lines))
 assert digest(DEST/'saferali_mmc1.txt')==expected['saferali_mmc1.txt']
 print('PASS: public inputs match frozen receipts')
if __name__=='__main__':main()
