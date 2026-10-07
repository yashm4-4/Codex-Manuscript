#!/usr/bin/env python3
"""Reference-check and enumerate GWAS/QTL identities. Never run coloc."""
import collections
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import requests
import pysam
from pyliftover import LiftOver

STAGE=Path(__file__).resolve().parents[1]
GAP=STAGE.parent
A=GAP/'fine-mapping-input-resolution-1.0/tracks/A/results'
FA37=GAP/'fine-mapping-input-resolution-1.0/references/GRCh37.p13.genome.fa'
SOURCE=STAGE/'sources'
START,END=88300892,91560531
GENE='ENSG00000138640'
FIELDS=['source_row','native_variant','native_pos','native_effect','native_other','native_or','native_se','native_p','native_beta_log_or','native_statistic_pass','native_identity_state','provider_harmonized_row','provider_code','chain_status','b38_position','b38_ref','b38_alt','b38_effect','b38_other','b38_beta','b38_varbeta','qtl_variant','qtl_beta_alt','qtl_se','qtl_p','qtl_maf','qtl_n','qtl_source_rows','status','reason']

def read(path):
    op=gzip.open if str(path).endswith('.gz') else open
    with op(path,'rt',newline='') as f:
        yield from csv.DictReader(f,delimiter='\t')
def write(path,rows,fields):
    with open(path,'w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore',lineterminator='\n')
        w.writeheader();w.writerows(rows)
def sha(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''):h.update(b)
    return h.hexdigest()
def parseid(s):
    parts=s.split('_',3)
    if len(parts)!=4:return None
    try:return (parts[0] if parts[0].startswith('chr') else 'chr'+parts[0],int(parts[1]),parts[2].upper(),parts[3].upper())
    except ValueError:return None
def pal(a,b):return len(a)==len(b)==1 and {a,b} in ({'A','T'},{'C','G'})
def finite(s):
    try:return math.isfinite(float(s))
    except:return False

# A frozen source row is joined to the provider harmonized row only by the
# frozen transport audit's complete source statistic/allele signature.
g=[x for x in read(A/'locus_variants_reference_checked.tsv.gz') if x['locus_id']=='A_chr4_88300892_91560531']
assert len(g)==12094
h=collections.defaultdict(list)
for x in read(A/'harmonized_locus_ledger.tsv.gz'):
    if x['native_chrom']=='4' and START<=int(x['native_pos'])<=END:h[x['source_row']].append(x)
q=list(read(SOURCE/'fam13a_full_cis_source.tsv.gz'))
qby=collections.defaultdict(list)
qpositions=collections.defaultdict(set)
for x in q:
    assert x['gene_id']==GENE and x['molecular_trait_id']==GENE
    qby[x['variant']].append(x)
    qpositions[x['position']].add(x['variant'])
assert len(q)==7264

# One public hg38 regional reference query, independent of the GWAS and QTL
# allele records. Coordinates in the API are zero-based, half-open.
refpath=SOURCE/'ucsc_hg38_chr4_87000000_92000000.fa'
refurl='https://api.genome.ucsc.edu/getData/sequence?genome=hg38;chrom=chr4;start=87000000;end=92000000'
if not refpath.exists():
    z=requests.get(refurl,timeout=90);z.raise_for_status();j=z.json()
    assert j['chrom']=='chr4' and j['start']==87000000 and j['end']==92000000 and len(j['dna'])==5000000
    with open(refpath,'w') as f:
        f.write('>chr4:87000000-92000000 UCSC hg38 zero-based half-open\n')
        for i in range(0,len(j['dna']),80):f.write(j['dna'][i:i+80].upper()+'\n')
fa38=pysam.FastaFile(str(refpath))
fa37=pysam.FastaFile(str(FA37))
lift=LiftOver(str(SOURCE/'hg19ToHg38.over.chain.gz'))
back=LiftOver(str(SOURCE/'hg38ToHg19.over.chain.gz'))
def ref38(pos,length):
    if not (87000001<=pos and pos+length-1<=92000000):return ''
    return fa38.fetch('chr4:87000000-92000000',pos-1-87000000,pos-1-87000000+length).upper()

qrows=[]
for v,group in sorted(qby.items(),key=lambda z: (int(z[1][0]['position']),z[0])):
    first=group[0]
    same=all(all(x[k]==first[k] for k in ['chromosome','position','ref','alt','maf','pvalue','beta','se','ac','an']) for x in group)
    parsed=parseid(v)
    refok=bool(parsed and parsed[0]=='chr4' and parsed[1]==int(first['position']) and parsed[2]==first['ref'] and parsed[3]==first['alt'] and ref38(parsed[1],len(parsed[2]))==parsed[2])
    qrows.append({**first,'source_rows':len(group),'duplicate_consistent':str(same),'b38_reference_pass':str(refok)})
write(STAGE/'fam13a_gtex_v8_lung_full_cis.tsv',qrows,list(qrows[0]))

cross=[];reject=[];gout=[];usedq=set()
for x in g:
    row=x['source_row']; hist=h.get(row,[])
    rec={k:'' for k in FIELDS}
    rec.update(source_row=row,native_variant=x['variant_id'],native_pos=x['base_pair_location'],native_effect=x['effect_allele'],native_other=x['other_allele'],native_or=x['odds_ratio'],native_se=x['standard_error'],native_p=x['p_value'],native_beta_log_or=x['beta_logOR'],native_statistic_pass=x['wald_p_rounding_pass'],native_identity_state=x['harmonization_status'])
    reasons=[]
    if not x['wald_p_rounding_pass']=='True':reasons.append('GWAS_PRECISION_EXCEPTION')
    if not x['harmonization_status'].startswith('REFERENCE_VERIFIED') and not x['harmonization_status'].startswith('PROVIDER_REFERENCE_NORMALIZED'):reasons.append('NATIVE_IDENTITY_UNRESOLVED')
    if len(hist)!=1:reasons.append('MISSING_OR_AMBIGUOUS_PROVIDER_MAPPING')
    native37=fa37.fetch('chr4',int(x['base_pair_location'])-1,int(x['base_pair_location'])).upper()
    if len(x['effect_allele'])==len(x['other_allele'])==1 and native37 not in (x['effect_allele'],x['other_allele']):reasons.append('GRCH37_REF_CONFLICT')
    if len(hist)==1:
        a=hist[0];parsed=parseid(a['hm_variant_id'])
        rec.update(provider_harmonized_row=a['harmonized_row'],provider_code=a['hm_code'])
        if a['transformation_status']!='PASS' or a['source_transport_match']!='EXACT_STATISTIC_AND_ALLELE_SIGNATURE_NOT_RSID_ONLY':reasons.append('PROVIDER_TRANSPORT_UNRESOLVED')
        if not parsed or parsed[0]!='chr4':reasons.append('BAD_GRCH38_ID')
        else:
            _,pos,ref,alt=parsed
            rec.update(b38_position=pos,b38_ref=ref,b38_alt=alt,b38_effect=a['hm_effect_allele'],b38_other=a['hm_other_allele'])
            fw=lift.convert_coordinate('chr4',int(x['base_pair_location'])-1)
            bw=back.convert_coordinate('chr4',pos-1)
            chain_ok=(len(fw)==1 and fw[0][0]=='chr4' and fw[0][1]==pos-1 and len(bw)==1 and bw[0][0]=='chr4' and bw[0][1]==int(x['base_pair_location'])-1 and fw[0][2]==bw[0][2])
            rec['chain_status']='UNIQUE_ROUNDTRIP_EXACT' if chain_ok else 'CHAIN_ROUNDTRIP_OR_PROVIDER_POSITION_CONFLICT'
            if not chain_ok:reasons.append('CHAIN_ROUNDTRIP_OR_PROVIDER_POSITION_CONFLICT')
            if ref38(pos,len(ref))!=ref:reasons.append('GRCH38_REF_CONFLICT')
            if {a['hm_effect_allele'],a['hm_other_allele']}!={ref,alt}:reasons.append('HARMONIZED_ALLELE_CONFLICT')
            if pal(ref,alt):reasons.append('PALINDROMIC_UNRESOLVED_NO_GWAS_AF')
            qv='chr4_'+str(pos)+'_'+ref+'_'+alt;rec['qtl_variant']=qv
            group=qby.get(qv,[])
            if not group:
                reasons.append('NO_EXACT_QTL_VARIANT')
                if qpositions.get(str(pos)):reasons.append('SAME_POSITION_DIFFERENT_ALLELE')
            else:
                z=group[0]; rec.update(qtl_beta_alt=z['beta'],qtl_se=z['se'],qtl_p=z['pvalue'],qtl_maf=z['maf'],qtl_n=str(int(z['an'])//2),qtl_source_rows=len(group))
                if not all(all(t[k]==z[k] for k in ['position','ref','alt','maf','pvalue','beta','se','ac','an']) for t in group):reasons.append('QTL_DUPLICATE_CONFLICT')
                if ref38(pos,len(ref))!=z['ref']:reasons.append('QTL_REF_CONFLICT')
                if not all(finite(z[k]) for k in ['beta','se','pvalue','maf']) or float(z['se'])<=0 or not 0<float(z['maf'])<=.5:reasons.append('QTL_STATISTIC_INVALID')
                if finite(a['hm_odds_ratio']) and finite(x['odds_ratio']) and finite(x['beta_logOR']):
                    hb=math.log(float(a['hm_odds_ratio']));nb=float(x['beta_logOR'])
                    sign=1 if abs(hb-nb)<abs(hb+nb) else -1
                    if abs(hb-sign*nb)<.00002:
                        rec['b38_beta']=str(sign*nb)
                        rec['b38_varbeta']=str(float(x['standard_error'])**2) if finite(x['standard_error']) else ''
                if not rec['b38_beta']:reasons.append('GWAS_EFFECT_ORIENTATION_UNRESOLVED')
                if not reasons:usedq.add(qv)
    rec['status']='EXACT_PASS' if not reasons else 'REJECTED'
    rec['reason']=';'.join(reasons) if reasons else '.'
    cross.append(rec);gout.append(rec.copy())
    if reasons:reject.append({'source':'GWAS','source_id':row,'variant':rec['qtl_variant'],'reasons':rec['reason']})
write(STAGE/'fam13a_gwas_region.tsv',gout,FIELDS)
write(STAGE/'fam13a_exact_gwas_eqtl_crosswalk.tsv',cross,FIELDS)
for v,group in qby.items():
    if v not in usedq:reject.append({'source':'QTL','source_id':v,'variant':v,'reasons':'NO_ELIGIBLE_GWAS_MATCH'})
write(STAGE/'harmonization_rejections.tsv',reject,['source','source_id','variant','reasons'])
counts=collections.Counter(x['status'] for x in cross)
rc=collections.Counter(reason for x in cross for reason in x['reason'].split(';') if reason and reason!='.')
coverage=[{'metric':'gwas_rows_native_locus','count':len(g)}, {'metric':'qtl_tested_rows_full_cis','count':len(q)}, {'metric':'qtl_distinct_variants_full_cis','count':len(qby)}, {'metric':'gwas_provider_mapping_rows','count':sum(bool(h.get(x['source_row'])) for x in g)}, {'metric':'exact_coordinate_ref_alt_candidate_rows','count':sum(bool(x['qtl_beta_alt']) for x in cross)}, {'metric':'eligible_exact_intersection','count':counts['EXACT_PASS']}, {'metric':'gwas_without_eligible_qtl','count':len(g)-counts['EXACT_PASS']}, {'metric':'qtl_without_eligible_gwas','count':len(qby)-len(usedq)}]
coverage += [{'metric':'direct_reference_or_harmonized_allele_conflicts','count':sum(rc.get(k,0) for k in ['GRCH37_REF_CONFLICT','GRCH38_REF_CONFLICT','HARMONIZED_ALLELE_CONFLICT','QTL_REF_CONFLICT'])}, {'metric':'chain_roundtrip_conflicts','count':rc.get('CHAIN_ROUNDTRIP_OR_PROVIDER_POSITION_CONFLICT',0)}, {'metric':'unresolved_palindromic_gwas_rows','count':rc.get('PALINDROMIC_UNRESOLVED_NO_GWAS_AF',0)}]
coverage += [{'metric':'rejection_'+k.lower(),'count':v} for k,v in sorted(rc.items())]
write(STAGE/'regional_coverage_summary.tsv',coverage,['metric','count'])
audit={'gwas_rows':len(g),'qtl_rows':len(q),'qtl_distinct':len(qby),'candidate_exact_rows':sum(bool(x['qtl_beta_alt']) for x in cross),'eligible_exact_intersection':counts['EXACT_PASS'],'gwas_reason_counts':dict(rc),'qtl_duplicate_extra':len(q)-len(qby),'qtl_ref_mismatch':sum(x['b38_reference_pass']=='False' for x in qrows),'qtl_duplicate_conflicts':sum(x['duplicate_consistent']=='False' for x in qrows),'gwas_region_sha256':sha(STAGE/'fam13a_gwas_region.tsv'),'qtl_source_sha256':sha(SOURCE/'fam13a_full_cis_source.tsv.gz'),'hg38_reference_sha256':sha(refpath),'hg19_to_hg38_chain_sha256':sha(SOURCE/'hg19ToHg38.over.chain.gz'),'hg38_to_hg19_chain_sha256':sha(SOURCE/'hg38ToHg19.over.chain.gz'),'ucsc_reference_url':refurl,'no_coloc_executed':True}
(STAGE/'validation_records.json').write_text(json.dumps(audit,indent=2)+'\n')
print(json.dumps(audit,indent=2))
