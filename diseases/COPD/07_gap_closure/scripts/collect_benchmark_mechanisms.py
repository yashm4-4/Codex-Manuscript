#!/usr/bin/env python3
"""Collect public COPD benchmark evidence independently of all model scores.

Authorized by attachment31c0f2f1-803e-4c5e-81ad-bd39763c2f0a. Downloaded raw
responses are immutable; content failures are preserved, not treated as sources.
Manual primary-source extraction is in mechanisms/curated_extraction.json.
"""
import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from urllib.request import Request, urlopen
import xml.etree.ElementTree as ET

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'data/COPD-V2-BENCH/mechanisms'

def retrieve(url,name):
    OUT.mkdir(parents=True,exist_ok=True)
    path=OUT/name
    record=dict(url=url,local_path=str(path.relative_to(ROOT)),retrieval_utc=datetime.now(timezone.utc).isoformat())
    try:
        if path.exists():
            data=path.read_bytes();record['status']='cached_existing_download'
        else:
            with urlopen(Request(url,headers={'User-Agent':'COPD-public-literature-audit/1.0'}),timeout=60) as response:
                data=response.read()
                record.update(status='downloaded',content_type=response.headers.get('Content-Type'),final_url=response.url)
            path.write_bytes(data)
        record.update(bytes=len(data),sha256=hashlib.sha256(data).hexdigest())
    except Exception as exc:record.update(status='failed',error=str(exc))
    ledger_path=OUT/'retrievals.json'
    ledger=json.loads(ledger_path.read_text()) if ledger_path.exists() else []
    ledger.append(record);ledger_path.write_text(json.dumps(ledger,indent=2)+'\n')
    return record

def xmltext(name):
    for elem in ET.parse(OUT/name).iter():
        if elem.tag in ('p','title','caption','table','supplementary-material'):
            print(elem.tag,elem.get('id',''),' '.join(''.join(elem.itertext()).split()))

def assemble():
    result=json.loads((OUT/'curated_extraction.json').read_text())
    result['sources'].append(dict(study_id='LIN2020_FAM13A',doi='10.1016/j.isci.2020.100928',
        pmcid='PMC7063182',title='FAM13A Represses AMPK Activity and Regulates Hepatic Glucose and Lipid Metabolism',
        ascertainment='HDL/adiposity/metabolic traits; shared FAM13A locus previously screened by Castaldi COPD MPRA',
        role='shared_locus_metabolic_reporter_followup',files=['Lin2020.html','Lin2020_supplements.zip'],
        notes='Reuses45Castaldi MPRA hits, not independent MPRA replication. New HepG2 reporter assays are eligible through COPD FAM13A locus connection; metabolic ascertainment is not independent COPD disease validation.'))
    result['searches'].append(dict(query='PMC7063182 FAM13A Represses AMPK Activity directly linked Castaldi follow-up',
        outcome='Reused MPRA distinguished from new hepatic reporter/deletion; all6 reporter contrasts retained with exact figure labels, then bounded search closed.'))
    for rs,pos,orientation,a1,a2,state in [
        ('rs2276936',89726283,'forward','A','C','positive'),('rs2276936',89726283,'reverse','T','G','ambiguous'),
        ('rs2167750',89730074,'forward','C','T','ambiguous'),('rs2167750',89730074,'reverse','G','A','ambiguous'),
        ('rs7695177',89737558,'forward','C','G','ambiguous'),('rs7695177',89737558,'reverse','G','C','ambiguous')]:
        result['assays'].append(dict(assay_id='LIN2020_'+rs+'_'+orientation,study_id='LIN2020_FAM13A',
            locus='FAM13A_4q22',rsid=rs,source_build='GRCh37',source_chrom='4',source_pos=pos,
            tested_allele1=a1,tested_allele2=a2,reported_strand='construct orientation shown; genomic forward strand not explicitly stated',
            assay_class='conventional_reporter',cell_context='HepG2',experimental_state=state,assay_orientation=orientation,
            reported_effect_direction='A higher reporter activity than C' if state=='positive' else 'not_established',
            reported_effect_allele='A' if state=='positive' else '',higher_activity_allele='A' if state=='positive' else '',
            p_value='<0.05' if state=='positive' else 'not_reported',
            source_rule='Figure1B *P<0.05 unpaired Student t test; unstarred contrasts without NS/numerical P remain ambiguous, not inferred null.',
            mechanism_in_model_scope='yes',source_locator='Figure1A GRCh37 identity; Figure1B '+orientation+' reporter pair',
            notes='Three biological replicates in triplicate. Metabolic follow-up of COPD FAM13A locus, not independent COPD disease validation. Figure explicitly labels reverse-insert complementary bases; tested letters preserved. rs7695177 is palindromic. Reused Castaldi MPRA is separate context, not a new experiment.'))
    result['contexts'].append(dict(context_id='LIN2020_reused_MPRA',study_id='LIN2020_FAM13A',
        locus='FAM13A_4q22',rsid='rs2276936;rs2167750;rs7695177',assay_class='MPRA_allele_effect',
        cell_context='Original Castaldi MPRA; not new Lin experiment',experimental_state='positive',
        source_locator='Results first paragraph',evidence_unit='reused_source_evidence',
        notes='Three members of45published Castaldi hits selected for HDL association. Do not count as independent assays/replication or duplicate benchmark positives.'))
    result['contexts'].append(dict(context_id='LIN2020_regional_deletion',study_id='LIN2020_FAM13A',
        locus='FAM13A_4q22',rsid='rs2276936',assay_class='CRISPR_deletion_region',cell_context='HepG2',
        experimental_state='positive',source_locator='Figure1C-D',evidence_unit='region_non_COPD_followup',
        notes='Approximately100bp deletion or arbitrary indels lower FAM13A by approximately30%; not exact C/A nucleotide editing.'))
    sources={s['study_id']:s for s in result['sources']}
    for s in sources.values():
        s['source_url']=f"https://pmc.ncbi.nlm.nih.gov/articles/{s['pmcid']}/"
        s['complete_assay_denominator']=False
    # Explicit regional positive/null readouts; not natural-allele labels.
    def ctx(id,rs,cell,gene,state,panel,pvalue,note):
        result['contexts'].append(dict(context_id=id,study_id='STUART2020',locus='TGFB1_19q13',
            rsid=rs,assay_class='CRISPRi_region',cell_context=cell,target_gene=gene,
            experimental_state=state,source_locator=panel,p_value=pvalue,source_rule='Student t test P<0.05',
            evidence_unit='region_gene_readout',notes=note))
    for cell in ['A549','H292','H441','primary lung fibroblasts']:
        for gene in ['TGFB1','B9D2','TMEM91']:
            ctx('STUART2020_rs1800469_'+cell+'_'+gene,'rs1800469',cell,gene,
                'positive' if gene=='TGFB1' else 'null','Figure1C' if cell=='A549' else 'Figure2',
                {'TGFB1':'0.002','B9D2':'0.37','TMEM91':'0.999'}[gene] if cell=='A549' else 'see Figure2',
                'Lentiviral regional CRISPRi: TGFB1 decreased, other2genes explicitly nonsignificant. Cells carry different genotypes but no isogenic natural allele comparison was performed.')
    for gene,pvalue in [('TGFB1','0.00002'),('B9D2','0.00026'),('TMEM91','0.00061')]:
        ctx('STUART2020_rs2241712_'+gene,'rs2241712','A549',gene,'positive','Figure6C',pvalue,
            'Synthetic guide regional perturbation lowers expression; no natural allele comparison.')
    for gene,pvalue,state in [('TGFB1','0.0045','positive'),('B9D2','>0.05','null'),('TMEM91','>0.05','null')]:
        ctx('STUART2020_rs1800469_synthetic_'+gene,'rs1800469','A549',gene,state,'Figure5',pvalue,
            'Synthetic guide regional perturbation; kept distinct from lentiviral guide experiments.')
    # Preserve NPNT null endpoints explicitly without presenting a splice-positive
    # variant as enhancer-null or treating repeated endpoint tests as variants.
    for iso,pval in [(1,'0.06'),(2,'0.11'),(3,'0.30'),(6,'0.25')]:
        result['contexts'].append(dict(context_id='SAFERALI2025_isoform'+str(iso)+'_null',
            study_id='SAFERALI2025_NPNT',locus='NPNT_4q24',rsid='rs34712979',assay_class='splicing',
            cell_context='Human LTRC lung; fiveAA/fiveGG donor set',experimental_state='null',
            source_locator='Table2 isoform'+str(iso),p_value=pval,evidence_unit='within_assay_endpoint',
            mechanism_in_model_scope='no',mechanism_scope_rationale='Splice isoform endpoint; not enhancer/H3K27me3 function.',
            notes='Nonsignificant isoform-proportion contrast within an otherwise splice-positive assay. Same donor set and SNP as master assay; not an additional independent variant.'))
    rationales={
        'yes':'Allele-specific transcriptional reporter is in broad regulatory-sequence scope; plasmid and cell context are not equivalent to frozen bulk-lung models.',
        'partial':'Isolated TF-binding or endogenous expression effect supports regulation but is not a direct enhancer/H3K27me3 classification assay.',
        'no':'Splice acceptor/exon inclusion mechanism; no independent enhancer or H3K27me3 allele mechanism established.'}
    for a in result['assays']:
        for k,v in dict(source_build='not_reported',source_chrom='',source_pos='',reported_strand='not_explicitly_reported',fdr='',assay_orientation='not_reported').items():a.setdefault(k,v)
        a['mechanism_scope_rationale']=rationales[a['mechanism_in_model_scope']]
        a['source_url']=sources[a['study_id']]['source_url'];a['ascertainment']=sources[a['study_id']]['ascertainment']
    for c in result['contexts']:
        c['source_url']=sources[c['study_id']]['source_url']
        c.setdefault('mechanism_in_model_scope','unevaluable')
        c.setdefault('mechanism_scope_rationale','Region/contact/gene/haplotype evidence is not an exact assayed natural nucleotide contrast.')
    retrievals=json.loads((OUT/'retrievals.json').read_text())
    for r in retrievals:
        if r.get('bytes'):
            p=ROOT/r['local_path'];data=p.read_bytes();extension=p.suffix.lower()
            valid=(extension=='.pdf' and data.startswith(b'%PDF')) or (extension in ('.zip','.docx','.xlsx') and data.startswith(b'PK')) or (extension=='.jpg' and data.startswith(b'\xff\xd8')) or (extension=='.xml' and b'<article' in data[:5000]) or (extension=='.html' and b'<article' in data)
            r['content_validation']='expected_content' if valid else 'non_source_response_preserved'
            if not valid:r['status']='unavailable_non_source_response'
    result['retrievals']=retrievals
    for q in result['searches']:q.update(search_date='2026-10-05',system='web search; primary sources used for labels')
    result['denominator_audit'].update(exact_assay_rows=len(result['assays']),
        unique_assayed_rsids=len({a['rsid'] for a in result['assays']}),context_rows=len(result['contexts']))
    ids=[a['assay_id'] for a in result['assays']]
    assert len(ids)==len(set(ids))
    assert len(result['assays'])==16
    assert len({a['rsid'] for a in result['assays']})==8
    assert all(a['experimental_state'] in ('positive','null','ambiguous') for a in result['assays'])
    dest=ROOT/'data/COPD-V2-BENCH/mechanisms_evidence.json'
    dest.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(output=str(dest),assay_rows=len(result['assays']),context_rows=len(result['contexts']),sources=len(sources)),indent=2))

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--url');p.add_argument('--name');p.add_argument('--text');p.add_argument('--assemble',action='store_true')
    args=p.parse_args()
    if args.url:print(json.dumps(retrieve(args.url,args.name),indent=2))
    if args.text:xmltext(args.text)
    if args.assemble:assemble()

if __name__=='__main__':main()
