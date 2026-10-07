from pathlib import Path
import pandas as pd,json,datetime,shutil,hashlib
R=Path(__file__).resolve().parents[1];O=R/'results';backup=O/'attempt1_reference_serialization';backup.mkdir(exist_ok=True)
raw=pd.read_csv(O/'locus_variants.tsv.gz',sep='\t',dtype=str,keep_default_na=False).set_index('source_row',drop=False)
changes=[]
for name in ['locus_variants_reference_checked.tsv.gz','risk_direction_reference_checked.tsv.gz','unresolved_variant_ledger.tsv.gz']:
 p=O/name;shutil.copy2(p,backup/name);d=pd.read_csv(p,sep='\t',dtype=str,keep_default_na=False);source=raw.loc[d.source_row].reset_index(drop=True);assert list(d.source_row)==list(source.source_row)
 updated={}
 for col in ['variant_id','p_value','chromosome','base_pair_location','effect_allele','other_allele','odds_ratio','standard_error','p','se','beta','beta_logOR','z','neglog10_p']:
  if col in d and col in source:updated[col]=int((d[col]!=source[col]).sum());d[col]=source[col]
 for col,alias in [('analysis_beta','beta'),('analysis_z','z')]:
  if col in d:updated[col]=int((d[col]!=source[alias]).sum());d[col]=source[alias]
 d.to_csv(p,sep='\t',index=False,compression={'method':'gzip','mtime':0});changes.append({'path':str(p.relative_to(R)),'changed_string_counts':updated,'before_sha256':hashlib.sha256((backup/name).read_bytes()).hexdigest(),'after_sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
for name in ['artifact_validation.json','finalize_audit.log','reference_audit.log','reference_audit.json']:
 if (O/name).exists():shutil.copy2(O/name,backup/name)
(O/'source_string_restoration.json').write_text(json.dumps({'scope':'Output serialization only; restored exact original/previously computed source strings by source_row, with all original first-attempt derivatives and failed validation retained','scientific_gate_changes':False,'allele_or_normalization_changes':False,'risk_direction_changes':False,'files':changes,'completed_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()},indent=2)+'\n');print(json.dumps(changes,indent=2))
