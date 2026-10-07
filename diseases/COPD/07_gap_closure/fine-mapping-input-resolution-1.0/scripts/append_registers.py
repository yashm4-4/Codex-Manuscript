"""Append this completed stage once; preserve byte-exact previous register prefixes."""
from pathlib import Path
import csv
import datetime
import hashlib
import io
import json

s=Path(__file__).resolve().parents[1];repo=s.parents[3];gap=s.parent
i=json.loads((s/'tables/readiness_integration.json').read_text())
assert i['final'] and i['diagnostic_matrices_completed']==16 and i['execution_cleared_loci']==0
assert (s/'FINE_MAPPING_INPUT_RESOLUTION_REPORT.md').is_file()
init=json.loads((s/'provenance/initialization.json').read_text())
date=datetime.datetime.now(datetime.timezone.utc).strftime('%Y-%m-%d')
stage=s.name
rows={
'activity_log.tsv':[
 ['COPD-V2-ACT-056',date,'complete_input_resolution','Acquire and audit complete authorized A/B/C GWAS inputs',
  'Frozen fine-mapping-preflight-1.0 and new investigator authorization',
  stage+'/tables/full_file_acquisition.tsv;'+stage+'/tables/complete_file_integrity.tsv',
  'A9886853 B8678470 C28987534 complete rows; native schemas signed contrasts missingness duplicates provider transformations and failed attempts preserved.'],
 ['COPD-V2-ACT-057',date,'complete_diagnostics_only','Define own-GWAS loci and resolve signed variant and LD inputs',
  'A Kim ever-smoker;B combined BBJ;C PanUKB J44 EUR;GRCh37 reference;selected PanUKB signed LD blocks',
  stage+'/tables/prospective_loci.tsv;'+stage+'/tables/gwas_only_risk_direction.tsv.gz;'+stage+'/ld/extraction_validation.json',
  '23 merged loci;21 non-MHC;2 MHC deferred. 482169 locus rows;368624 verified directions.16 A/C matrices diagnosed;no Japanese matrix substituted. No candidates/models/functional evidence used.'],
 ['COPD-V2-ACT-058',date,'complete_input_resolution_stop','Integrate readiness and freeze for investigator review',
  'Complete GWAS harmonization signed-LD provenance numerical diagnostics and method/sample contracts',
  stage+'/FINE_MAPPING_INPUT_RESOLUTION_REPORT.md;'+stage+'/provenance/freeze.json',
  'All A/B/C NOT CLEARED;0 complete summary gates;0 complete LD gates;0 execution-cleared loci. No SuSiE/FINEMAP fitting PIPs credible sets posterior inference model/candidate scoring target genes colocalization manuscript edit author contact commit or push.']
],
'gap_closure_decision_register.tsv':[
 ['COPD-V2-DEC-059',date,'preserve_phenotype_and_ancestry_roles','Retain A primary practical direct COPD B independent Japanese ancestry C secondary EHR sensitivity',
  'Availability does not promote J44 to primary COPD;A/C overlap in UKB and are not independent replication.Only three authorized tracks analyzed.',stage+'/COHORT_OVERLAP_AND_SCOPE.md'],
 ['COPD-V2-DEC-060',date,'retain_unresolved_contracts','Distinguish signed direction arithmetic agreement and matrix numerics from method-ready summary-LD covariance',
  'A smoking-subgroup versus external full-EUR LD N and significant identity loss remain unresolved;B/C SPA-derived SE cannot independently establish ordinary Wald/RSS covariance;B dense signed Japanese LD unavailable.',stage+'/tables/readiness_gates.tsv;'+stage+'/ld/source_contracts.json'],
 ['COPD-V2-DEC-061',date,'no_gate_relaxation','Retain precision identity density numerical and provenance failures under frozen tolerances',
  'No LD repair sign optimization missing-band zero fill biological selection or inferred risk direction. Corrected implementation deviations retain prior attempts and do not redefine gates.',stage+'/config/ld_diagnostic_contract.json;'+stage+'/tables/GWAS_LD_compatibility.tsv'],
 ['COPD-V2-DEC-062',date,'stop_for_investigator_review','Do not authorize statistical fine-mapping or publication from completion of input resolution',
  'A/B/C all NOT CLEARED;no complete execution package. Future uniform-prior multi-signal SuSiE-RSS/FINEMAP design pinned but software/sample contracts and explicit new investigator authorization remain required. Sakornsakolpat request draft unsent.',stage+'/execution_inputs/manifest.json;'+stage+'/provenance/freeze.json']
],
'gap_closure_result_register.tsv':[
 ['COPD-V2-FMR-GWAS-001',date,'complete','whole_file_signed_input_audit',
  'Complete GWAS rows A9886853 B8678470 C28987534. Full compression transport/checksum and schema audits distinguish provider digest from local SHA/CRC evidence.',stage+'/tables/complete_file_integrity.tsv','no'],
 ['COPD-V2-FMR-LOCI-001',date,'complete','prospective_GWAS_only_loci',
  'Own-GWAS P<5e-8 autosomal +/-1.5Mb transitive merge: A11/10nonMHC B5/5 C7/6;23total21nonMHC. Entire intersecting MHC loci deferred;chr20 boundary clipped. No candidate or gene selection.',stage+'/tables/prospective_loci.tsv','no'],
 ['COPD-V2-FMR-DIRECTION-001',date,'complete_with_explicit_unresolved_rows','GWAS_only_signed_risk_direction',
  '482169 rows retained;368624 signed directions assigned A140886 B51989 C175749;113545 unresolved. Estimates of association direction only;no causality or regulatory-model claims.',stage+'/tables/gwas_only_risk_direction.tsv.gz','no'],
 ['COPD-V2-FMR-LD-001',date,'complete_diagnostics_only','signed_LD_input_and_compatibility_audit',
  'Selected signed PanUKB EUR reference used externally for A and as same-project near-overlap for C.16 matrices with source block order allele numeric residual rank provenance checks. B dense signed Japanese LD unavailable. Complete scientific LD gates0.',stage+'/tables/LD_numerical_QC.tsv','no'],
 ['COPD-V2-FMR-READINESS-001',date,'NOT CLEARED','input_execution_readiness',
  'A/B/C all NOT CLEARED;0 loci passing all method-ready summary gates;0 all-LD-gate loci;0 execution-cleared;no future input packages. Signed effects and basic matrix passes do not close sampling/likelihood contracts.',stage+'/tables/track_readiness.tsv','no'],
 ['COPD-V2-FMR-FREEZE-001',date,'complete_input_resolution_stop','input_resolution_freeze',
  'Versioned local input/diagnostic bundle and uniform-prior future software design frozen for review;no statistical fine-mapping PIPs credible sets posterior effects model/candidate scores target genes colocalization manuscript change external author contact commit or push.',stage+'/provenance/freeze.json','no']
]}
receipts=[]
for baseline in init['baseline_registers']:
    p=repo/baseline['path']; original=p.read_bytes()
    assert len(original)==baseline['bytes'] and hashlib.sha256(original).hexdigest()==baseline['sha256'], 'Refusing duplicate append or changed baseline: '+str(p)
    newrows=rows[p.name]
    header=next(csv.reader(io.StringIO(original.decode()),delimiter='\t'))
    assert all(len(r)==len(header) for r in newrows)
    buf=io.StringIO();csv.writer(buf,delimiter='\t',lineterminator='\n').writerows(newrows)
    appended=buf.getvalue().encode()
    with p.open('ab') as f:f.write(appended)
    current=p.read_bytes();assert current[:len(original)]==original
    receipts.append(dict(path=baseline['path'],before_bytes=len(original),before_sha256=baseline['sha256'],
        appended_rows=len(newrows),appended_bytes=len(appended),bytes=len(current),sha256=hashlib.sha256(current).hexdigest(),append_only=True))
(s/'provenance/register_updates.json').write_text(json.dumps(receipts,indent=2)+'\n')
print(json.dumps({'registers':len(receipts),'appended_rows':sum(x['appended_rows'] for x in receipts)}))
