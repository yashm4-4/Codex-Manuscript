#!/usr/bin/env python3
"""Independent COPD external-benchmark validation, after source freeze only.

Does not import assembly/comparison code. Reconstructs scores, eligibility, calls,
rank fidelity, missingness categories and case-series summaries from raw frozen
inputs. No fitting, sequence inference, model execution or orientation averaging.
Run only after investigator workflow has completed source freeze and comparison.
"""
import argparse
import ast
import csv
import gzip
import hashlib
import json
import math
import re
import subprocess
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[4]
BASE=ROOT/'diseases/COPD/07_gap_closure'
PROV=BASE/'provenance'
RES=BASE/'results'
DATA=BASE/'data/COPD-V2-BENCH'
P='COPD-V2-BENCH'
LOCKED_SPEC='46d3aac26b1f980d7c2dcd5034d6cc64566524ca7d99cc635f7d979e04f2de65'
CHECKS=[]
SCRIPTS={'initialize_external_benchmark.py','assemble_external_benchmark.py',
         'compare_external_benchmark.py','validate_external_benchmark.py',
         'collect_benchmark_castaldi.py','collect_benchmark_gong.py','collect_benchmark_mechanisms.py'}
SHARED={'README.md','.gitignore','activity_log.tsv','gap_closure_decision_register.tsv',
        'gap_closure_result_register.tsv','data/README.md','logs/README.md',
        'provenance/README.md','results/README.md','scripts/README.md','manuscript/README.md'}

def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1048576),b''):h.update(block)
    return h.hexdigest()

def read(path):
    op=gzip.open if str(path).endswith('.gz') else open
    with op(path,'rt',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))

def write(path,rows,fields):
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)

def check(name,ok,detail=''):
    CHECKS.append(dict(check=name,status='PASS' if ok else 'FAIL',detail=str(detail)))

def batch(name,errors,count):
    check(name,not errors,f'n={count}; mismatches={len(errors)}; examples={errors[:5]}')

def boolean(value):
    if value in (True,False):return value
    if str(value).lower() in ('true','1'):return True
    if str(value).lower() in ('false','0'):return False
    raise ValueError('Invalid boolean '+repr(value))

def serial(value):
    if value is None:return ''
    if isinstance(value,(dict,list)):return json.dumps(value,ensure_ascii=False,sort_keys=True)
    return str(value)

def equal(a,b):
    if a is None:return b==''
    if isinstance(a,bool):return str(b).lower()==str(a).lower()
    if isinstance(a,(int,float)):
        try:return math.isclose(float(a),float(b),rel_tol=1e-12,abs_tol=2e-15)
        except (ValueError,TypeError):return False
    if isinstance(a,(dict,list)):
        try:return a==json.loads(b)
        except (ValueError,TypeError):return False
    return str(a)==str(b)

def complement(x):return x.translate(str.maketrans('ACGT','TGCA'))[::-1]

def variant_key(r):
    return r.get('canonical_variant_id') or 'unresolved:'+r.get('rsid','')+':'+ '/'.join(sorted([r.get('tested_allele1',''),r.get('tested_allele2','')]))

def category(a,b):
    if a is None or b is None:return 'unevaluable'
    return {(True,True):'stable_both',(True,False):'forward_only',(False,True):'RC_only',(False,False):'negative_both'}[(a,b)]

def model_context(e,s):return 'both' if e and s else 'enhancer_only' if e else 'H3K27me3_only' if s else 'neither'

def tstamp(x):return datetime.fromisoformat(x.replace('Z','+00:00'))

class Reference:
    """Independent faidx access using only Python standard library."""
    def __init__(self,path):
        self.f=Path(path).open('rb');self.index={}
        with Path(str(path)+'.fai').open() as f:
            for line in f:
                name,length,offset,width,bytes_per_line,*_=line.rstrip().split('\t')
                self.index[name]=tuple(map(int,(length,offset,width,bytes_per_line)))
    def fetch(self,chrom,start,end):
        name=chrom if chrom in self.index else 'chr'+chrom
        length,offset,width,stride=self.index[name]
        if start<0 or end>length or start>end:raise ValueError('Reference interval invalid')
        chunks=[]
        while start<end:
            count=min(end-start,width-start%width)
            self.f.seek(offset+start//width*stride+start%width)
            chunks.append(self.f.read(count));start+=count
        return b''.join(chunks).decode().upper()

def normalize(chrom,pos,ref,alt,genome):
    ref='' if ref=='-' else ref;alt='' if alt=='-' else alt
    if not ref or not alt:
        anchor=genome.fetch(chrom,pos-2,pos-1);pos-=1;ref=anchor+ref;alt=anchor+alt
    while ref[-1]==alt[-1] and ref!=alt:
        if min(len(ref),len(alt))==1:
            if pos==1:break
            anchor=genome.fetch(chrom,pos-2,pos-1);pos-=1;ref=anchor+ref[:-1];alt=anchor+alt[:-1]
        else:ref=ref[:-1];alt=alt[:-1]
    while min(len(ref),len(alt))>1 and ref[0]==alt[0]:pos+=1;ref=ref[1:];alt=alt[1:]
    return f'{chrom}:{pos}:{ref}:{alt}'

def protect_inputs(lock,freeze,run,assembly):
    spec=ROOT/lock['protocol_path']
    check('Locked protocol exact SHA256',digest(spec)==lock['protocol_sha256']==LOCKED_SPEC)
    ledger=ROOT/lock['frozen_input_ledger']
    check('Initial immutable ledger SHA256',digest(ledger)==lock['frozen_input_ledger_sha256'])
    rows=read(ledger)
    check('All316 protected baseline files enumerated',len(rows)==lock['frozen_input_count']==316)
    for r in rows:
        p=ROOT/r['path']
        check('Protected baseline: '+r['path'],p.exists() and digest(p)==r['sha256'] and p.stat().st_size==int(r['bytes']))
    check('No unauthorized commit',subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()==lock['baseline_git_head'])
    check('Freeze external evaluation status',freeze['status']=='FROZEN_EXTERNAL_EVALUATION_DATA' and freeze['no_training_or_model_selection_use'] is True)
    check('Protocol hash consistent across lifecycle',all(x.get('protocol_sha256')==LOCKED_SPEC for x in [freeze,run,assembly]))
    check('Lock precedes assembly precedes freeze precedes comparison',tstamp(lock['protocol_locked_utc'])<=tstamp(assembly['assembled_utc'])<=tstamp(freeze['frozen_utc'])<tstamp(run['comparison_started_utc'])<=tstamp(run['completed_utc']))
    check('Comparison links exact freeze',run['benchmark_freeze_sha256']==digest(PROV/f'{P}_benchmark_freeze.json'))
    check('Comparison implementation hash matches executed code',run['script_sha256']==digest(BASE/'scripts/compare_external_benchmark.py'))
    check('Every disclosed comparison run follows source freeze',all(tstamp(x['comparison_started_utc'])>tstamp(freeze['frozen_utc']) for x in run.get('previous_runs',[])))
    check('Freeze links exact assembly',freeze['assembly_manifest_sha256']==digest(PROV/f'{P}_assembly_manifest.json'))
    for path,expected in freeze['frozen_files'].items():
        check('Frozen evidence/table: '+path,(ROOT/path).exists() and digest(ROOT/path)==expected)
    for path,expected in run['inputs'].items():check('Comparison input: '+path,digest(ROOT/path)==expected)
    for path,expected in run['outputs'].items():check('Comparison output: '+path,digest(ROOT/path)==expected)
    check('No inference/new sequences/threshold changes in run',all(run.get(k) is False for k in ['new_model_inference','new_sequences','thresholds_changed']))
    check('Source assembler did not read model outcomes',assembly.get('model_outputs_read') is False)
    check('Prespecified prior-knowledge disclosure retained',freeze.get('known_prior_selected_case_observations_are_not_blinded') is True)

def source_checks(master):
    source_rows=[];contexts=[]
    for name in ['castaldi','gong','mechanisms']:
        obj=json.loads((DATA/(name+'_evidence.json')).read_text());source_rows.extend(obj['assays']);contexts.extend(obj.get('contexts',[]))
    source_by={r['assay_id']:r for r in source_rows};by={r['assay_id']:r for r in master}
    check('Source/master exact assay identity multiset',len(source_rows)==len(source_by)==len(master)==len(by) and set(source_by)==set(by))
    errors=[]
    for aid,source in source_by.items():
        m=by.get(aid,{})
        for key,value in source.items():
            target='reported_locus_annotation' if source['study_id']=='GONG2026' and key=='locus' else key
            if not equal(value,m.get(target,'')):errors.append((aid,key,str(value)[:60],m.get(target,'')))
    batch('Every source allele/label/context/value preserved in master',errors,len(source_rows))
    harm={r['assay_id']:r for r in read(RES/f'{P}-R002_exact_allele_harmonization.tsv')}
    classrows={r['assay_id']:r for r in read(RES/f'{P}-R004_assay_mechanism_classification.tsv')}
    for name,table in [('harmonization',harm),('classification',classrows)]:
        errors=[]
        for aid,row in table.items():
            for k,v in row.items():
                if k in by.get(aid,{}) and v!=by[aid][k]:errors.append((aid,k))
        if set(table)!=set(by):errors.append('Assay sets differ')
        batch('Source '+name+' agrees with master',errors,len(table))
    ctx=read(RES/f'{P}-R004B_contextual_and_excluded_evidence.tsv')
    errors=[]
    if len(ctx)!=len(contexts):errors.append('Context counts differ')
    for i,(original,observed) in enumerate(zip(contexts,ctx)):
        for k,v in original.items():
            if not equal(v,observed.get(k,'')):errors.append((i,k))
    batch('Region/haplotype/excluded contexts preserved separately',errors,len(contexts))
    check('All source labels use prespecified states',all(r['experimental_state'] in {'positive','null','ambiguous','conflicting','unavailable','unevaluable'} for r in master))
    check('All rows retain external-evaluation firewall',all(boolean(r['external_evaluation_only']) and not boolean(r['model_outcomes_used_for_inclusion']) for r in master))

def allele_checks(master):
    mapping={}
    for p in sorted((DATA/'allele_mappings').glob('GRCh38_*.json.gz')):
        with gzip.open(p,'rt') as f:obj=json.load(f)
        mapping.update({k:v for k,v in obj.items() if re.fullmatch('rs[0-9]+',k) and isinstance(v,dict)})
    genome=Reference(ROOT/'models/TREDNET_v2/fasta/hg38.fa')
    ref_errors=[];map_errors=[];pair_errors=[];gate_errors=[];cache={}
    for r in master:
        canonical=r['canonical_variant_id']
        if not canonical:
            if r['identity_status']=='exact':gate_errors.append(r['assay_id'])
            continue
        if not re.fullmatch(r'(?:[0-9]+|X|Y):[0-9]+:[ACGT]+:[ACGT]+',canonical):ref_errors.append(canonical);continue
        chrom,pos,ref,alt=canonical.split(':');pos=int(pos)
        if canonical not in cache:cache[canonical]=genome.fetch(chrom,pos-1,pos-1+len(ref))==ref
        if not cache[canonical] or ref==alt:ref_errors.append(canonical)
        if r['identity_status']!='exact' or not boolean(r['reference_verified']):gate_errors.append(r['assay_id'])
        if canonical!=':'.join([r['grch38_chrom'],r['grch38_pos'],r['grch38_ref'],r['grch38_alt']]):gate_errors.append(r['assay_id']+' fields')
        mapped=set()
        for m in mapping.get(r['rsid'],{}).get('mappings',[]):
            if m.get('assembly_name')!='GRCh38' or m.get('coord_system')!='chromosome' or int(m.get('strand',0))!=1:continue
            if str(m.get('seq_region_name'))!=chrom:continue
            alleles=m.get('allele_string','').split('/')
            for other in alleles[1:]:
                if re.fullmatch('[ACGT-]+',alleles[0]) and re.fullmatch('[ACGT-]+',other):
                    mapped.add(normalize(chrom,int(m['start']),alleles[0],other,genome))
        if canonical not in mapped:map_errors.append((r['rsid'],canonical))
        a,b=r.get('tested_allele1',''),r.get('tested_allele2','')
        if len(ref)==len(alt)==len(a)==len(b)==1 and set(a+b)<=set('ACGT'):
            if {a,b}!={ref,alt} and {complement(a),complement(b)}!={ref,alt}:pair_errors.append(r['assay_id'])
        if boolean(r.get('liftover_used','False')):gate_errors.append('Unexpected liftover:'+r['assay_id'])
        if boolean(r.get('palindromic','False')) and not boolean(r.get('direction_identity_resolved','False')) and r.get('reported_direction_alt_minus_ref'):gate_errors.append('Palindromic direction:'+r['assay_id'])
    for name,errors in [('GRCh38 references independently faidx-verified',ref_errors),('Exact canonical identities supported by authoritative cached mappings',map_errors),('Source SNP contrasts match canonical alleles or documented complements',pair_errors),('Identity gates and palindromic direction safeguards',gate_errors)]:batch(name,errors,len(master))
    genome.f.close()

def reconstruct(master):
    model=ROOT/'diseases/COPD/04_modeling/results'
    fwd={r['candidate_record_id']:r for r in read(model/'COPD-S4-R003_candidate_allele_scores.tsv.gz')}
    meta={r['candidate_record_id']:r for r in read(model/'COPD-S4-R004_prioritized_candidates.tsv.gz')}
    ranks={r['candidate_record_id']:int(r['predicted_causal_priority_rank']) for r in read(model/'COPD-S4-R010_THE_LIST.tsv')}
    rc={r['candidate_record_id']:r for r in read(RES/'COPD-V2-RC-R001_reverse_complement_scores.tsv.gz')}
    check('Independent frozen universe counts',(len(fwd),len(meta),len(ranks),len(rc))==(15303,15389,337,15303))
    check('Raw forward and RC exact pair sets equal',set(fwd)==set(rc))
    result={}
    # Numeric constants are independently pinned to original V1 thresholds,
    # not read from benchmark predictions or chosen using experimental labels.
    thresholds={'enhancer':(.643623,.05706318769999998,.04936093850000001),
                'silencer':(.58505,.028802613899999996,.026183359500000003)}
    for v in {r['canonical_variant_id'] for r in master if r['canonical_variant_id']}:
        row=dict(v1_model_evaluable=v in fwd,rc_model_evaluable=v in rc,
                 v1_candidate_rank=ranks.get(v),v1_candidate_status=(v in ranks) if v in meta else None,
                 v1_match_status='exact_existing_scored_pair' if v in fwd else 'exact_existing_unscorable_record' if v in meta else 'no_frozen_V1_record')
        if v not in fwd:
            row.update({name+'_orientation_category':'unevaluable' for name in ['enhancer','silencer','union']});result[v]=row;continue
        m=meta[v];blacklist=boolean(m['blacklisted']);eligible=boolean(m['causal_call_eligible']);kind=m['variant_class_group']
        check('Scored exact eligibility consistent with blacklist '+v,eligible==(not blacklist))
        row.update(v1_blacklisted=blacklist,v1_call_eligible=eligible,v1_class_group=kind)
        for name,(region_cut,snv_cut,indel_cut) in thresholds.items():
            delta_cut=snv_cut if kind=='SNV' else indel_cut
            row[name+'_region_threshold']=region_cut;row[name+'_abs_delta_threshold']=delta_cut
            for orient,raws in [('forward',fwd),('rc',rc)]:
                a=float(raws[v][name+'_ref_score']);b=float(raws[v][name+'_alt_score']);d=b-a
                pre=name+'_'+orient
                rg=max(a,b)>=region_cut;dg=abs(d)>=delta_cut
                row.update({pre+'_ref_score':a,pre+'_alt_score':b,pre+'_delta':d,pre+'_abs_delta':abs(d),
                    pre+'_serialized_delta':float(raws[v][name+'_delta_alt_minus_ref']),pre+'_region_score':max(a,b),
                    pre+'_region_gate':rg,pre+'_delta_gate':dg,pre+'_call':eligible and rg and dg})
                if orient=='forward':check('Original authoritative call '+v+' '+name,row[pre+'_call']==boolean(m['predicted_causal_'+name]))
            row[name+'_orientation_category']=category(row[name+'_forward_call'],row[name+'_rc_call'])
        for orient in ['forward','rc']:
            e,s=row['enhancer_'+orient+'_call'],row['silencer_'+orient+'_call']
            row[orient+'_union_call']=e or s;row[orient+'_model_context']=model_context(e,s)
        row['union_orientation_category']=category(row['forward_union_call'],row['rc_union_call'])
        check('Original337membership '+v,row['forward_union_call']==(v in ranks))
        result[v]=row
    sentinel=result.get('4:88963935:G:T',{})
    check('rs2013701 frozen rank210 and forward enhancer recovery',sentinel.get('v1_candidate_rank')==210 and sentinel.get('enhancer_forward_call') is True)
    check('rs2013701 RC loss solely fixed delta gate',sentinel.get('enhancer_rc_call') is False and sentinel.get('enhancer_rc_region_gate') is True and sentinel.get('enhancer_rc_delta_gate') is False and sentinel.get('v1_call_eligible') is True)
    return result

def compare_checks(master,computed):
    tables={name:read(RES/f'{P}-{name}') for name in ['R007_V1_forward_comparison.tsv','R008_forward_RC_comparison.tsv']}
    masterby={r['assay_id']:r for r in master};independent=[]
    for a in master:
        r=dict(a);r.update(computed.get(a['canonical_variant_id'],dict(v1_model_evaluable=False,rc_model_evaluable=False,
            v1_match_status='unresolved_assayed_identity',v1_candidate_status=None,v1_candidate_rank=None,
            union_orientation_category='unevaluable',enhancer_orientation_category='unevaluable',silencer_orientation_category='unevaluable')))
        state,scope=a['experimental_state'],a['mechanism_in_model_scope']
        if scope=='no':con='outside_model_scope_not_a_sequence_model_false_negative'
        elif not r['v1_model_evaluable']:con='model_unevaluable_not_a_false_negative'
        elif scope!='yes':con='partial_or_unresolved_mechanism_descriptive_only'
        elif state=='positive':con='reported_active_forward_recovered' if r['forward_union_call'] else 'reported_active_forward_not_recovered'
        elif state=='null':con='assay_null_model_positive' if r['forward_union_call'] else 'assay_null_model_negative'
        else:con='experimental_label_not_evaluable_for_binary_performance'
        r['concordance_category']=con;r['enhancer_reporter_direction_agreement']='unevaluable'
        if r['v1_model_evaluable'] and state=='positive' and a['assay_class'] in ['conventional_reporter','MPRA_allele_effect'] and a.get('reported_direction_alt_minus_ref') in ['1','-1']:
            predicted=(r['enhancer_forward_delta']>0)-(r['enhancer_forward_delta']<0)
            r['enhancer_reporter_direction_agreement']='same' if predicted==int(a['reported_direction_alt_minus_ref']) else 'opposite' if predicted else 'zero'
        r['H3K27me3_activity_direction_interpretation']='not_a_validated_reporter_activation_or_repression_direction'
        independent.append(r)
    expected={r['assay_id']:r for r in independent}
    for name,rows in tables.items():
        by={r['assay_id']:r for r in rows};errors=[];source_errors=[]
        check(name+' every source assay once',len(rows)==len(master)==len(by) and set(by)==set(masterby))
        for aid,r in by.items():
            for k,v in masterby.get(aid,{}).items():
                if r.get(k)!=v:source_errors.append((aid,k))
            for k,v in expected.get(aid,{}).items():
                if name.startswith('R007') and ('_rc_' in k or k.startswith('rc_') or 'orientation_category' in k):continue
                if not equal(v,r.get(k,'')):errors.append((aid,k,str(v)[:40],r.get(k,'')))
        batch(name+' all frozen source fields preserved',source_errors,len(rows))
        batch(name+' independently reconstructed scores/calls/directions/categories',errors,len(rows))
    npnt=[r for r in independent if r['rsid']=='rs34712979' and r['assay_class']=='splicing']
    check('NPNT splice mechanism retained out-of-scope, never enhancer/H3K27me3 false negative',bool(npnt) and all(r['mechanism_in_model_scope']=='no' and r['concordance_category']=='outside_model_scope_not_a_sequence_model_false_negative' for r in npnt))
    check('Missing model pairs never negative predictions',all(r['union_orientation_category']=='unevaluable' and ('not_a_false_negative' in r['concordance_category'] or r['mechanism_in_model_scope']=='no') for r in independent if not r['v1_model_evaluable']))
    return independent

def summary_checks(rows,run):
    groups=defaultdict(list)
    for r in rows:groups[variant_key(r)].append(r)
    actual={r['benchmark_variant_key']:r for r in read(RES/f'{P}-R015_unique_variant_summary.tsv')}
    check('Unique variant summary contains every exact/unresolved key',set(groups)==set(actual))
    errors=[];positive=[]
    for key,group in groups.items():
        inscope=[r for r in group if r['mechanism_in_model_scope']=='yes'];labels={r['experimental_state'] for r in inscope}
        directions={r['reported_direction_alt_minus_ref'] for r in inscope if r['experimental_state']=='positive' and r.get('reported_direction_alt_minus_ref') in ['1','-1']}
        state='conflicting' if 'conflicting' in labels or {'positive','null'}<=labels or len(directions)>1 else 'positive' if 'positive' in labels else 'null' if labels=={'null'} else 'unevaluable'
        values=dict(evidence_summary_state=state,assay_rows=len(group),any_in_scope_positive_assay='positive' in labels,
            in_scope_positive_assay_rows=sum(r['experimental_state']=='positive' for r in inscope),
            in_scope_null_assay_rows=sum(r['experimental_state']=='null' for r in inscope),opposite_activity_directions_across_contexts=len(directions)>1,
            rsids=';'.join(sorted({r['rsid'] for r in group})),study_ids=';'.join(sorted({r['study_id'] for r in group})))
        values['positive_reporter_direction_assay_counts']=dict(Counter(r['enhancer_reporter_direction_agreement'] for r in group
            if r['experimental_state']=='positive' and r['assay_class'] in ['MPRA_allele_effect','conventional_reporter']))
        for k,v in values.items():
            if not equal(v,actual.get(key,{}).get(k,'')):errors.append((key,k))
        for k,v in group[0].items():
            if k in actual.get(key,{}) and (k.startswith(('v1_','enhancer_','silencer_','forward_','rc_')) or k in ['canonical_variant_id','union_orientation_category']):
                if not equal(v,actual[key][k]):errors.append((key,k))
        if 'positive' in labels:positive.append((group[0],state))
    batch('R015 independent evidence states, counts, model values',errors,len(groups))
    covered=[(r,state) for r,state in positive if r['v1_model_evaluable']]
    expected=dict(any_positive_unique_variant_keys=len(positive),model_covered=len(covered),model_unevaluable=len(positive)-len(covered),
        forward_union_recovered=sum(r['forward_union_call'] for r,_ in covered),rc_union_recovered=sum(r['rc_union_call'] for r,_ in covered),
        orientation_categories=dict(Counter(r['union_orientation_category'] for r,_ in covered)),context_heterogeneous_among_covered=sum(s=='conflicting' for _,s in covered))
    check('Independent aggregate unique-positive recovery',all(run['unique_any_positive_in_scope_summary'].get(k)==v for k,v in expected.items()),expected)
    for suffix,fields in [('R009_study_specific_summary.tsv',('study_id','assay_class','cell_context','assay_orientation')),
                         ('R010_mechanism_specific_summary.tsv',('assay_class','mechanism_in_model_scope')),
                         ('R011_locus_specific_summary.tsv',('study_id','locus')),
                         ('R014_benchmark_coverage_flow.tsv',('study_id',))]:
        grouped=defaultdict(list)
        for r in rows:grouped[tuple(r.get(k,'') for k in fields)].append(r)
        actual={tuple(r.get(k,'') for k in fields):r for r in read(RES/f'{P}-{suffix}')};errors=[]
        if set(actual)!=set(grouped):errors.append('Grouping key mismatch')
        for key,group in grouped.items():
            values=dict(assay_rows=len(group),unique_variant_keys=len({variant_key(r) for r in group}),
                unique_reported_rsids=len({r['rsid'] for r in group if r['rsid']}),unique_exact_variants=len({r['canonical_variant_id'] for r in group if r['canonical_variant_id']}),
                model_evaluable_assay_rows=sum(r['v1_model_evaluable'] for r in group),model_evaluable_unique_variants=len({r['canonical_variant_id'] for r in group if r['v1_model_evaluable']}))
            for state in ['positive','null','ambiguous','conflicting','unavailable','unevaluable']:
                sub=[r for r in group if r['experimental_state']==state];values[state+'_assay_rows']=len(sub);values[state+'_unique_variant_keys']=len({variant_key(r) for r in sub})
                covered={r['canonical_variant_id']:r for r in sub if r['v1_model_evaluable']}
                values[state+'_model_covered_unique_variants']=len(covered)
                values[state+'_forward_union_call_unique_variants']=sum(r['forward_union_call'] for r in covered.values())
                values[state+'_rc_union_call_unique_variants']=sum(r['rc_union_call'] for r in covered.values())
            active={r['canonical_variant_id']:r for r in group if r['experimental_state']=='positive' and r['mechanism_in_model_scope']=='yes' and r['v1_model_evaluable']}
            values['in_scope_positive_model_evaluable_unique']=len(active)
            for name,field in [('enhancer','enhancer_forward_call'),('H3K27me3','silencer_forward_call'),('union','forward_union_call')]:
                count=sum(r[field] for r in active.values());values['in_scope_positive_forward_'+name+'_count']=count;values['in_scope_positive_forward_'+name+'_case_recovery_fraction']=count/len(active) if active else None
            for cat in ['stable_both','forward_only','RC_only','negative_both']:values['in_scope_positive_union_'+cat]=sum(r['union_orientation_category']==cat for r in active.values())
            for k,v in values.items():
                if not equal(v,actual.get(key,{}).get(k,'')):errors.append((key,k))
        batch(suffix+' independent strata and recovery counts',errors,len(grouped))

def statistical_gates(rows):
    approved={r['study_id'] for r in read(RES/f'{P}-R005_assayed_denominator_audit.tsv') if r.get('complete_assayed_denominator_valid','').lower()=='true'}
    groups=defaultdict(list)
    for r in rows:groups[(r['study_id'],r['assay_class'],r['cell_context'],r.get('assay_orientation',''))].append(r)
    actual={(r['study_id'],r['assay_class'],r['cell_context'],r.get('assay_orientation','')):r for r in read(RES/f'{P}-R012_quantitative_metric_eligibility.tsv')}
    errors=[]
    for key,group in groups.items():
        pos={variant_key(r) for r in group if r['experimental_state']=='positive'};null={variant_key(r) for r in group if r['experimental_state']=='null'}
        valid=key[0] in approved and all(r['experimental_state'] in ['positive','null'] and r['mechanism_in_model_scope']=='yes' and r['v1_model_evaluable'] for r in group) and bool(pos) and bool(null) and not(pos&null)
        out=actual.get(key,{})
        if valid or out.get('full_panel_metrics_valid','').lower()!='false':errors.append((key,'valid panel requires prespecified metrics'))
        if out.get('positive_variant_keys')!=str(len(pos)) or out.get('null_variant_keys')!=str(len(null)):errors.append((key,'class counts'))
        if not out.get('reason'):errors.append((key,'missing omission rationale'))
        if any(out.get(k)!='not_calculated' for k in ['AUROC','average_precision','sensitivity_specificity_FPR','confidence_interval']):errors.append((key,'invalid quantitative metrics'))
    if set(actual)!=set(groups):errors.append('Missing/extra metric strata')
    batch('All performance gates independently fail; no unsupported AUROC/AP/FPR/CI',errors,len(groups))

def distribution_checks(rows):
    groups=defaultdict(list)
    for r in rows:
        if r['v1_model_evaluable']:groups[(r['study_id'],r['assay_class'],r['cell_context'],r.get('assay_orientation',''),r['experimental_state'])].append(r)
    table=read(RES/f'{P}-R013_score_distributions.tsv')
    actual={(r['study_id'],r['assay_class'],r['cell_context'],r.get('assay_orientation',''),r['experimental_state'],r['model'],r['score']):r for r in table}
    errors=[];expected_keys=set()
    for groupkey,group in groups.items():
        for model in ['enhancer','silencer']:
            for metric in ['region_score','abs_delta']:
                key=groupkey+(model,metric);expected_keys.add(key)
                values=sorted({r['canonical_variant_id']:r[model+'_forward_'+metric] for r in group}.values())
                def percentile(frac):
                    at=(len(values)-1)*frac;lo=math.floor(at);hi=math.ceil(at)
                    return values[lo]*(hi-at)+values[hi]*(at-lo) if hi!=lo else values[lo]
                expected=dict(n=len(values),min=values[0],max=values[-1],mean=sum(values)/len(values),
                              p25=percentile(.25),median=percentile(.5),p75=percentile(.75))
                for field,value in expected.items():
                    if not equal(value,actual.get(key,{}).get(field,'')):errors.append((key,field))
    if set(actual)!=expected_keys or len(actual)!=len(table):errors.append('Distribution grouping mismatch/duplicates')
    batch('R013 independently recomputed model-covered score distributions',errors,len(expected_keys))

def additional_summary_checks(rows):
    actual={(r['study_id'],r['model'],r['orientation']):r for r in read(RES/f'{P}-R017_positive_case_gate_failures.tsv')}
    errors=[];keys=set()
    for study in ['ALL_DEDUPLICATED']+sorted({r['study_id'] for r in rows}):
        selected={r['canonical_variant_id']:r for r in rows if (study=='ALL_DEDUPLICATED' or r['study_id']==study)
            and r['experimental_state']=='positive' and r['mechanism_in_model_scope']=='yes' and r['v1_model_evaluable']}
        for model in ['enhancer','silencer']:
            for orientation in ['forward','rc']:
                key=(study,model,orientation);keys.add(key);counts=Counter()
                for r in selected.values():
                    prefix=model+'_'+orientation
                    if r[prefix+'_call']:state='recovered'
                    elif not r['v1_call_eligible']:state='original_eligibility_exclusion'
                    elif not r[prefix+'_region_gate'] and not r[prefix+'_delta_gate']:state='both_gates_fail'
                    elif not r[prefix+'_region_gate']:state='region_only_fail'
                    else:state='delta_only_fail'
                    counts[state]+=1
                values={'in_scope_positive_scored_unique_variants':len(selected),**{s:counts[s] for s in
                    ['recovered','original_eligibility_exclusion','both_gates_fail','region_only_fail','delta_only_fail']}}
                for field,value in values.items():
                    if not equal(value,actual.get(key,{}).get(field,'')):errors.append((key,field))
    if set(actual)!=keys:errors.append('Gate-failure grouping mismatch')
    batch('R017 independently reconstructed positive-case fixed-gate failures',errors,len(keys))
    exact={r['canonical_variant_id']:r for r in rows if r['v1_model_evaluable']}
    positive={r['canonical_variant_id']:r for r in rows if r['v1_model_evaluable'] and r['mechanism_in_model_scope']=='yes' and r['experimental_state']=='positive'}
    actual={(r['subset'],r['model']):r for r in read(RES/f'{P}-R018_orientation_recovery_summary.tsv')}
    errors=[];keys=set()
    for subset,selected in [('all_exact_scored_variants',exact),('any_in_scope_positive_scored_variants',positive)]:
        for model in ['enhancer','silencer','union']:
            key=(subset,model);keys.add(key)
            counts=Counter(r[model+'_orientation_category'] for r in selected.values())
            values={'unique_variants':len(selected),**{s:counts[s] for s in ['stable_both','forward_only','RC_only','negative_both']}}
            for field,value in values.items():
                if not equal(value,actual.get(key,{}).get(field,'')):errors.append((key,field))
    if set(actual)!=keys:errors.append('Orientation grouping mismatch')
    batch('R018 independently reconstructed orientation recovery counts',errors,len(keys))

def boundary_checks(require_handoff=True):
    changed=set(filter(None,subprocess.check_output(['git','diff','--name-only','-z','HEAD'],cwd=ROOT).decode().split('\0')))
    changed.update(filter(None,subprocess.check_output(['git','ls-files','--others','--exclude-standard','-z'],cwd=ROOT).decode().split('\0')))
    allowed=[];bad=[]
    for path in sorted(changed):
        prefix='diseases/COPD/07_gap_closure/'
        if not path.startswith(prefix):bad.append(path);continue
        relative=path[len(prefix):]
        okay=relative in SHARED or relative.startswith('data/COPD-V2-BENCH/') or Path(relative).name.startswith(P) or (relative.startswith('scripts/') and Path(relative).name in SCRIPTS)
        if 'COPD-V2-PHENO' in relative or 'COPD-V2-RC' in relative:okay=False
        (allowed if okay else bad).append(path)
    check('Git write boundary only benchmark and named shared V2 files',not bad,bad)
    # Parse rather than execute all new code. Training/inference libraries are absent.
    dangerous=[]
    for name in SCRIPTS:
        path=BASE/'scripts'/name
        tree=ast.parse(path.read_text())
        imported=[]
        for node in ast.walk(tree):
            if isinstance(node,ast.Import):imported.extend(n.name for n in node.names)
            elif isinstance(node,ast.ImportFrom):imported.append(node.module or '')
        if any(x.split('.')[0] in {'tensorflow','keras','torch','TREDNet_v2_inference'} for x in imported):dangerous.append(name)
    check('Benchmark scripts contain no model-training/inference imports',not dangerous,dangerous)
    if not require_handoff:return allowed
    reports=list(RES.glob(P+'*report.md'))
    check('Detailed Markdown report present',len(reports)==1)
    if reports:
        text=reports[0].read_text().lower()
        check('Report explicitly marks frozen external evaluation data','frozen external evaluation data' in text)
        check('Report retains rs2013701 and NPNT interpretation','rs2013701' in text and 'rs34712979' in text)
        check('Report discusses future-data firewall','training' in text and 'hyperparameter' in text and 'threshold' in text and 'internal' in text)
    for name in ['activity_log.tsv','gap_closure_decision_register.tsv','gap_closure_result_register.tsv']:
        check('Updated V2 register '+name,P in (BASE/name).read_text())
    return allowed

def finish(started,allowed,preliminary=False):
    if preliminary:
        out=PROV/f'{P}_preliminary_validation.tsv'
        write(out,CHECKS,['check','status','detail'])
        failures=[r for r in CHECKS if r['status']=='FAIL']
        print(json.dumps(dict(mode='preliminary_no_final_checksums',checks=len(CHECKS),failed=len(failures),failures=failures),indent=2))
        return 0 if not failures else 1
    out=PROV/f'{P}_final_validation.tsv'
    write(out,CHECKS,['check','status','detail'])
    manifest=PROV/f'{P}_final_validation_manifest.json'
    failures=[r for r in CHECKS if r['status']=='FAIL']
    meta=dict(started_utc=started,completed_utc=datetime.now(timezone.utc).isoformat(),
        validator_sha256=digest(__file__),protocol_sha256=LOCKED_SPEC,
        checks=len(CHECKS),passed=len(CHECKS)-len(failures),failed=len(failures),
        status='PASS' if not failures else 'FAIL',failures=failures,
        independent_of_assembly_comparison_imports=True,new_model_inference=False,
        output_path=str(out.relative_to(ROOT)),output_sha256=digest(out),
        artifact_checksum_path=f'diseases/COPD/07_gap_closure/provenance/{P}_artifact_checksums.tsv',
        artifact_checksum_self_excluded=True)
    manifest.write_text(json.dumps(meta,indent=2)+'\n')
    paths=set()
    for p in DATA.rglob('*'):
        if p.is_file() and '__pycache__' not in p.parts:paths.add(p)
    for folder in [PROV,RES,BASE/'logs']:
        for p in folder.glob(P+'*'):
            if p.is_file():paths.add(p)
    for name in SCRIPTS:paths.add(BASE/'scripts'/name)
    for path in allowed:
        p=ROOT/path
        if p.is_file() and p.relative_to(BASE).as_posix() in SHARED:paths.add(p)
    ledger=PROV/f'{P}_artifact_checksums.tsv';paths.discard(ledger)
    entries=[dict(path=str(p.relative_to(ROOT)),bytes=p.stat().st_size,sha256=digest(p)) for p in sorted(paths)]
    write(ledger,entries,['path','bytes','sha256'])
    print(json.dumps(dict(status=meta['status'],checks=meta['checks'],failed=meta['failed'],artifacts=len(entries)),indent=2))
    return 0 if not failures else 1

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--comparison-complete',action='store_true');parser.add_argument('--preliminary',action='store_true');args=parser.parse_args()
    if not args.comparison_complete:raise SystemExit('Refusing model-value access without explicit completed-freeze/comparison workflow signal.')
    started=datetime.now(timezone.utc).isoformat();allowed=[]
    try:
        freeze=json.loads((PROV/f'{P}_benchmark_freeze.json').read_text())
        run=json.loads((PROV/f'{P}_comparison_manifest.json').read_text())
        if not run.get('completed_utc'):raise ValueError('Comparison not complete; refusing model-value access')
        lock=json.loads((PROV/f'{P}_protocol_lock.json').read_text());assembly=json.loads((PROV/f'{P}_assembly_manifest.json').read_text())
        protect_inputs(lock,freeze,run,assembly)
        if any(r['status']=='FAIL' for r in CHECKS):raise ValueError('Immutable-input or lifecycle validation failed before model reconstruction')
        master=read(ROOT/freeze['master_path'])
        source_checks(master);allele_checks(master)
        computed=reconstruct(master);independent=compare_checks(master,computed)
        summary_checks(independent,run);statistical_gates(independent);distribution_checks(independent);additional_summary_checks(independent)
        analysis=read(RES/f'{P}-R016_analysis_validation.tsv')
        check('Analysis QC all PASS',bool(analysis) and all(r['status']=='PASS' for r in analysis),len(analysis))
        allowed=boundary_checks(require_handoff=not args.preliminary)
    except Exception as exc:check('Validator completed all checks without exception',False,repr(exc))
    raise SystemExit(finish(started,allowed,args.preliminary))

if __name__=='__main__':main()
