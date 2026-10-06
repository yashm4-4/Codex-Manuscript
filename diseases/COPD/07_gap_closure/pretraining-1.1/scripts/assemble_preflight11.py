#!/usr/bin/env python3
"""Assemble the completed outcome-blind attempt001; never rematch or run models.

Original string-valued feature fields are copied without numeric reformatting.
All output files use exclusive creation and reproducible gzip containers.
"""
from __future__ import annotations

import collections
import contextlib
import datetime as dt
import gzip
import hashlib
import io
import itertools
import json
from pathlib import Path
import resource
import time

import numpy as np
import pandas as pd

VERSION = Path(__file__).resolve().parents[1]
V2 = VERSION.parent
ROOT = V2.parents[2]
OLD = V2 / "data/COPD-V2-PREFLIGHT"
ATTEMPT = VERSION / "attempts/001_full_population"
MODELS = ["enhancer", "h3k27me3"]
CONFIGS = ["V2-A", "V2-B", "V2-C"]
PARTS = ["train", "validation", "test"]
SEEDS = [104729, 130363, 155921]
DISTANCE_VARS = ["gc_fraction", "atac_signal_percentile_max_train_only", "repeat_fraction_2001"]
SCALES = [.05, .20, .20]
COVARIATES = ["gc_fraction", "atac_signal_percentile_max", "atac_signal_percentile_max_train_only", "repeat_fraction_2001", "blacklist_input_any", "non_acgt_any", "non_acgt_fraction", "blacklist_bp_2001", "promoter_bp_2001", "atac_anchor_count"]
SCOPES = [("overall", "all")] + [("partition", p) for p in PARTS] + [("validation_role", r) for r in ["checkpoint", "selection", "calibration"]]


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8388608), b""):
            digest.update(block)
    return digest.hexdigest()


def info(path):
    path = Path(path)
    return {"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size, "sha256": sha(path)}


@contextlib.contextmanager
def output_text(path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.suffix == ".gz":
        with path.open("xb") as raw:
            with gzip.GzipFile(filename="", mode="wb", fileobj=raw, compresslevel=6, mtime=0) as zipped:
                with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as handle:
                    yield handle
    else:
        with path.open("x", newline="") as handle:
            yield handle


def write(path, frame, header=True):
    with output_text(path) as handle:
        frame.to_csv(handle, sep="\t", index=False, header=header, lineterminator="\n", float_format="%.17g")


def write_json(path, data):
    with output_text(path) as handle:
        json.dump(data, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write("\n")


def read(path, **kwargs):
    return pd.read_csv(path, sep="\t", dtype=str, keep_default_na=False, **kwargs)


def numeric(frame, column):
    return pd.to_numeric(frame[column], errors="coerce").to_numpy(dtype=float)


def subset(frame, scope_type, scope):
    return frame if scope_type == "overall" else frame[frame[scope_type].eq(scope)]


def distribution(values):
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    result = {"n_total": len(values), "n_finite": len(finite), "n_missing": len(values)-len(finite)}
    names = ["minimum", "q01", "q05", "q25", "median", "q75", "q95", "q99", "maximum"]
    quantiles = np.quantile(finite, [0, .01, .05, .25, .5, .75, .95, .99, 1]) if len(finite) else [np.nan]*9
    result.update(zip(names, quantiles))
    result.update(mean=float(finite.mean()) if len(finite) else np.nan, standard_deviation=float(finite.std(ddof=1)) if len(finite)>1 else np.nan)
    return result


def main():
    started = time.time()
    if (VERSION / "provenance/configuration_manifest.json").exists() or (VERSION / "provenance/freeze.json").exists():
        raise RuntimeError("Version is assembled/frozen; refusing overwrite")
    attempt = json.loads((ATTEMPT / "attempt_manifest.json").read_text())
    rule = json.loads((VERSION / "specification/matching_attempt_001.json").read_text())
    assert attempt["balance_failures"] == attempt["additional_chr7_role_failures"] == 0
    assert attempt["positive_retention"] == 1 and attempt["unique_one_to_one_controls"]
    inputs = [info(path) for path in [OLD / "interval_features.tsv.gz", OLD / "interval_role_assignment.tsv.gz", VERSION / "data/atac_normalized_features.tsv.gz", VERSION / "specification/matching_attempt_001.json", ATTEMPT / "attempt_manifest.json", ATTEMPT / "control_matching_pairs.tsv.gz", ATTEMPT / "initial_pairs.tsv.gz", ATTEMPT / "positive_retention.tsv", ATTEMPT / "covariate_balance.tsv", VERSION / "specification/selection_adequacy.json", Path(__file__)]]
    frame = read(OLD / "interval_features.tsv.gz")
    original_columns = frame.columns.tolist()
    roles = read(OLD / "interval_role_assignment.tsv.gz", usecols=["interval_id", "component_id", "validation_role"])
    ranks = read(VERSION / "data/atac_normalized_features.tsv.gz", usecols=["interval_id", "atac_signal_percentile_max_train_only"])
    assert frame.interval_id.is_unique and frame.interval_id.tolist() == roles.interval_id.tolist() == ranks.interval_id.tolist()
    frame["component_id"] = roles.component_id
    frame["validation_role"] = roles.validation_role
    frame["atac_signal_percentile_max_train_only"] = ranks.atac_signal_percentile_max_train_only
    frame["blacklist_input_any"] = (numeric(frame, "blacklist_bp_2001")>0).astype(int)
    frame["non_acgt_any"] = (numeric(frame, "non_acgt_fraction")>0).astype(int)
    del roles, ranks
    positions = pd.Series(frame.index.to_numpy(), index=frame.interval_id)
    assert len(frame) == 1143370
    # Frozen master uses canonical chromosome/numeric coordinate order.
    chrom_ord = frame.chrom.map({f"chr{i}": i for i in range(1,23)} | {"chrX":23,"chrY":24}).to_numpy()
    order = np.lexsort((frame.interval_id.to_numpy(), numeric(frame,"core_end"), numeric(frame,"core_start"), chrom_ord))
    assert np.array_equal(order, frame.index.to_numpy())
    final_pairs = read(ATTEMPT / "control_matching_pairs.tsv.gz")
    initial_pairs = read(ATTEMPT / "initial_pairs.tsv.gz")
    registry, counts, role_counts, distributions, lobes, sequence_rows = [], [], [], [], [], []
    sets, class_sets = {}, {}
    membership = collections.defaultdict(list)
    panel_infos = []
    rich_columns = original_columns + ["component_id", "validation_role", "atac_signal_percentile_max_train_only", "blacklist_input_any", "non_acgt_any", "configuration", "model", "label", "matched_positive_id", "matched_control_id", "matched_final_normalized_chebyshev_distance"]
    for model in MODELS:
        positive_a = frame.index[frame[f"v1_{model}_positive"].eq("1")].to_numpy()
        control_a = frame.index[frame[f"v1_{model}_control"].eq("1")].to_numpy()
        positive_c = frame.index[frame[f"v1_{model}_positive"].eq("1") & frame[f"{model}_same_lobe_peak_support"].eq("1")].to_numpy()
        for config in CONFIGS:
            pairs = final_pairs[final_pairs.configuration.eq(config) & final_pairs.model.eq(model)]
            positives = positive_c if config == "V2-C" else positive_a
            controls = control_a if config == "V2-A" else positions.loc[pairs.control_id].to_numpy(dtype=int)
            assert len(set(positives) & set(controls)) == 0 and len(set(controls)) == len(controls)
            if config != "V2-A":
                assert set(pairs.positive_id) == set(frame.loc[positives,"interval_id"])
                assert len(pairs) == len(positives) == len(controls)
            indices = np.sort(np.r_[positives, controls])
            x = frame.loc[indices].copy()
            x["configuration"], x["model"] = config, model
            x["label"] = np.isin(indices, positives).astype(int)
            for column in ["matched_positive_id", "matched_control_id", "matched_final_normalized_chebyshev_distance"]:
                x[column] = ""
            if config != "V2-A":
                pidx = positions.loc[pairs.positive_id].to_numpy(dtype=int)
                cidx = positions.loc[pairs.control_id].to_numpy(dtype=int)
                for target in (pidx, cidx):
                    x.loc[target,"matched_positive_id"] = pairs.positive_id.to_numpy()
                    x.loc[target,"matched_control_id"] = pairs.control_id.to_numpy()
                    x.loc[target,"matched_final_normalized_chebyshev_distance"] = pairs.final_normalized_chebyshev_distance.to_numpy()
            manifest_path = VERSION / "data/configurations" / f"{config}_{model}_interval_manifest.tsv.gz"
            write(manifest_path, x[rich_columns])
            registry.append({"configuration":config,"model":model,"manifest":str(manifest_path.relative_to(ROOT)),"n_intervals":len(x),"positive":len(positives),"control":len(controls),"sha256":sha(manifest_path),"configuration_role":"intended_corrected_model" if config=="V2-C" else "diagnostic_ablation","training_authorized":False})
            sets[(model,config)] = set(indices)
            class_sets[(model,config,"positive")] = set(positives)
            class_sets[(model,config,"control")] = set(controls)
            for idx in indices:
                membership[int(idx)].append(f"{config}:{model}")
            for part in PARTS:
                y = x[x.partition.eq(part)]
                for label, word in [(1,"positive"),(0,"control")]:
                    z = y[y.label.eq(label)]
                    bed_path = VERSION / "data/intervals" / f"{config}_{model}_{part}_{word}.bed.gz"
                    write(bed_path,z[["chrom","core_start","core_end","interval_id"]],header=False)
                    counts.append({"configuration":config,"model":model,"partition":part,"class":word,"n_intervals":len(z),"n_components":z.component_id.nunique(),"n_unique_encoded_sequences":z.canonical_rc_sequence_sha256.nunique(),"bed":str(bed_path.relative_to(ROOT)),"bed_sha256":sha(bed_path)})
            for role,y in x.groupby("validation_role",sort=True):
                positive=y[y.label.eq(1)];control=y[y.label.eq(0)]
                role_counts.append({"configuration":config,"model":model,"partition":y.partition.iloc[0],"role":role,"positive":len(positive),"control":len(control),"components":y.component_id.nunique(),"positive_components":positive.component_id.nunique(),"control_components":control.component_id.nunique(),"positive_unique_encoded_sequences":positive.canonical_rc_sequence_sha256.nunique(),"control_unique_encoded_sequences":control.canonical_rc_sequence_sha256.nunique(),"control_to_positive_ratio":len(control)/len(positive)})
            for scope_type,scope in SCOPES:
                y=subset(x,scope_type,scope)
                for label,word in [(1,"positive"),(0,"control")]:
                    z=y[y.label.eq(label)]
                    for variable in COVARIATES:
                        distributions.append({"configuration":config,"model":model,"scope_type":scope_type,"scope":scope,"class":word,"variable":variable,**distribution(numeric(z,variable))})
                    sequence_rows.append({"configuration":config,"model":model,"scope_type":scope_type,"scope":scope,"class":word,"n_intervals":len(z),"n_components":z.component_id.nunique(),"sequence_unavailable":int((z.sequence_available!="1").sum()),"wrong_sequence_length":int((numeric(z,"sequence_length")!=2001).sum()),"non_ACGT_sequences":int((numeric(z,"non_acgt_fraction")>0).sum()),"core_blacklist_overlap":int((numeric(z,"blacklist_bp_1kb")>0).sum()),"full_input_blacklist_overlap":int((numeric(z,"blacklist_bp_2001")>0).sum()),"core_promoter_overlap":int((numeric(z,"promoter_bp_1kb")>0).sum()),"full_input_promoter_overlap":int((numeric(z,"promoter_bp_2001")>0).sum()),"ATAC_anchored":int((numeric(z,"atac_anchor_count")>0).sum()),"ATAC_full_input_overlap":int((z.atac_overlap_2001=="1").sum()),"relevant_mark_full_input_overlap":int((z[f"{model}_mark_overlap_2001"]=="1").sum()),"known_positive_input_overlap":int((z[f"{model}_v1_positive_input_overlap"]=="1").sum()),"accessible_control_eligible":int((z[f"{model}_accessible_control_eligible"]=="1").sum()),"same_lobe_positive_support":int((z[f"{model}_same_lobe_peak_support"]=="1").sum()),"within_class_encoded_duplicate_rows":int(z.canonical_rc_sequence_sha256.duplicated(keep=False).sum()),"historical_or_positive_context_counts_are_descriptive":config=="V2-A" or label==1})
                a=y[y.label.eq(1)];b=y[y.label.eq(0)]
                for lobe in ["lower_left","lower_right","upper_right"]:
                    ap=int(a.atac_lobes.str.split(";").map(lambda values:lobe in values).sum());bp=int(b.atac_lobes.str.split(";").map(lambda values:lobe in values).sum())
                    gap=ap/len(a)-bp/len(b)
                    lobes.append({"configuration":config,"model":model,"scope_type":scope_type,"scope":scope,"provenance_type":"each_ATAC_lobe_support","level":lobe,"positive_n":len(a),"control_n":len(b),"positive_support_n":ap,"control_support_n":bp,"positive_proportion":ap/len(a),"control_proportion":bp/len(b),"proportion_gap":gap,"absolute_gap":abs(gap),"limit":.05,"status":"DESCRIPTIVE_HISTORICAL" if config=="V2-A" else "PASS" if abs(gap)<=.05 else "FAIL"})
                ac=a.atac_lobes.value_counts();bc=b.atac_lobes.value_counts()
                for signature in sorted(set(ac.index)|set(bc.index)):
                    ap=int(ac.get(signature,0));bp=int(bc.get(signature,0));gap=ap/len(a)-bp/len(b)
                    lobes.append({"configuration":config,"model":model,"scope_type":scope_type,"scope":scope,"provenance_type":"ATAC_lobe_signature","level":signature or "UNANCHORED","positive_n":len(a),"control_n":len(b),"positive_support_n":ap,"control_support_n":bp,"positive_proportion":ap/len(a),"control_proportion":bp/len(b),"proportion_gap":gap,"absolute_gap":abs(gap),"limit":.05,"status":"DESCRIPTIVE_HISTORICAL" if config=="V2-A" else "PASS" if abs(gap)<=.05 else "FAIL"})
            if config=="V2-C":
                panel=x[x.validation_role.isin(["selection","calibration","test"])].copy()
                assert not panel.validation_role.isin(["train","checkpoint"]).any()
                panel_path=VERSION/"data/evaluation"/f"{model}_common_challenge_panel.tsv.gz"
                write(panel_path,panel[rich_columns])
                panel_infos.append(info(panel_path))
                for role,y in panel.groupby("validation_role",sort=True):
                    p=y[y.label.eq(1)];c=y[y.label.eq(0)]
                    role_counts.append({"configuration":"COMMON_C_TASK","model":model,"partition":y.partition.iloc[0],"role":role,"positive":len(p),"control":len(c),"components":y.component_id.nunique(),"positive_components":p.component_id.nunique(),"control_components":c.component_id.nunique(),"positive_unique_encoded_sequences":p.canonical_rc_sequence_sha256.nunique(),"control_unique_encoded_sequences":c.canonical_rc_sequence_sha256.nunique(),"control_to_positive_ratio":len(c)/len(p)})
            print(f"Assembled {config} {model}: {len(positives)} positives / {len(controls)} controls",flush=True)
            del x
    # Quantify pair distances rather than imply reinstated local calipers.
    pair_summaries=[]
    for stage,table in [("initial",initial_pairs),("final",final_pairs)]:
        for (config,model),pairs in table.groupby(["configuration","model"],sort=True):
            pi=positions.loc[pairs.positive_id].to_numpy(dtype=int);ci=positions.loc[pairs.control_id].to_numpy(dtype=int)
            p=frame.loc[pi].reset_index(drop=True);c=frame.loc[ci].reset_index(drop=True)
            deltas=np.stack([np.abs(numeric(p,var)-numeric(c,var)) for var in DISTANCE_VARS],axis=1)
            cheby=np.max(deltas/np.asarray(SCALES),axis=1)
            exact_metric=np.max(np.abs(np.stack([numeric(p,var) for var in DISTANCE_VARS],axis=1)/np.asarray(SCALES)-np.stack([numeric(c,var) for var in DISTANCE_VARS],axis=1)/np.asarray(SCALES)),axis=1)
            stored=pd.to_numeric(pairs[f"{stage}_normalized_chebyshev_distance"]).to_numpy()
            assert np.allclose(exact_metric,stored,rtol=0,atol=1e-12)
            assert np.array_equal(p.partition.to_numpy(),c.partition.to_numpy()) and np.array_equal(p.validation_role.to_numpy(),c.validation_role.to_numpy())
            for scope_type,scope in SCOPES:
                mask=np.ones(len(p),dtype=bool) if scope_type=="overall" else p[scope_type].eq(scope).to_numpy()
                if not mask.any():continue
                for variable,values,scale in [("normalized_chebyshev_distance",exact_metric,1.)]+[(var+"_absolute_difference",deltas[:,j],SCALES[j]) for j,var in enumerate(DISTANCE_VARS)]:
                    chosen=values[mask]
                    pair_summaries.append({"stage":stage,"configuration":config,"model":model,"scope_type":scope_type,"scope":scope,"metric":variable,**distribution(chosen),"legacy_scale_not_current_caliper":scale,"above_legacy_scale_n":int((chosen>scale).sum()),"fraction_above_legacy_scale":float((chosen>scale).mean()),"cross_chromosome_n":int((p.chrom.to_numpy()[mask]!=c.chrom.to_numpy()[mask]).sum()),"cross_chromosome_fraction":float((p.chrom.to_numpy()[mask]!=c.chrom.to_numpy()[mask]).mean()),"cross_ATAC_signature_n":int((p.atac_lobes.to_numpy()[mask]!=c.atac_lobes.to_numpy()[mask]).sum()),"cross_ATAC_signature_fraction":float((p.atac_lobes.to_numpy()[mask]!=c.atac_lobes.to_numpy()[mask]).mean()),"interpretation":"distributional balance; no claim of local pairwise common support or biological exchangeability"})
    overlaps=[]
    label_changes=[]
    for model in MODELS:
        for scope_type,scope in [("overall","all")]+[("partition",part) for part in PARTS]:
            allowed=set(frame.index) if scope_type=="overall" else set(frame.index[frame.partition.eq(scope)])
            for a,b in itertools.combinations(CONFIGS,2):
                for word in ["all","positive","control"]:
                    aset=(sets[(model,a)] if word=="all" else class_sets[(model,a,word)]) & allowed
                    bset=(sets[(model,b)] if word=="all" else class_sets[(model,b,word)]) & allowed
                    overlaps.append({"model":model,"configuration_1":a,"configuration_2":b,"scope_type":scope_type,"scope":scope,"class":word,"n_1":len(aset),"n_2":len(bset),"shared":len(aset&bset),"only_1":len(aset-bset),"only_2":len(bset-aset),"Jaccard":len(aset&bset)/len(aset|bset),"B_C_interpretation":"same-lobe label correction plus induced control rematching" if (a,b)==("V2-B","V2-C") else "staged diagnostic ablation"})
            bpos=class_sets[(model,"V2-B","positive")]&allowed;cpos=class_sets[(model,"V2-C","positive")]&allowed
            assert cpos<=bpos
            label_changes.append({"model":model,"scope_type":scope_type,"scope":scope,"complete_V1_positive_n":len(bpos),"same_lobe_positive_n":len(cpos),"cross_lobe_only_removed_n":len(bpos-cpos),"fraction_removed":len(bpos-cpos)/len(bpos),"common_support_restriction_exclusions_B":0,"common_support_restriction_exclusions_C":0,"pure_membership_change_is_directly_measured":True,"trained_B_C_contrast_is_label_plus_control_rematching":True})
    active=sorted(set().union(*sets.values()))
    selected=frame.loc[active]
    duplicates=[]
    cross_partition_groups={};cross_role_groups={};duplicate_groups={}
    for kind,column in [("raw_exact","sequence_sha256"),("encoded_forward_RC","canonical_rc_sequence_sha256")]:
        duplicated=selected[selected[column].duplicated(keep=False)]
        groups=duplicated.groupby(column,sort=True)
        cross_partition_groups[kind]=int((groups.partition.nunique()>1).sum())
        cross_role_groups[kind]=int((groups.validation_role.nunique()>1).sum())
        duplicate_groups[kind]=groups.ngroups
        for identity,y in groups:
            for idx,row in y.iterrows():
                duplicates.append({"identity_kind":kind,"identity_sha256":identity,"group_rows":len(y),"group_partitions":";".join(sorted(y.partition.unique())),"group_validation_roles":";".join(sorted(y.validation_role.unique())),"interval_id":row.interval_id,"partition":row.partition,"component_id":row.component_id,"validation_role":row.validation_role,"configuration_model_membership":";".join(membership[int(idx)])})
    assert max(cross_partition_groups.values())==max(cross_role_groups.values())==0
    dup_columns=["identity_kind","identity_sha256","group_rows","group_partitions","group_validation_roles","interval_id","partition","component_id","validation_role","configuration_model_membership"]
    write(VERSION/"results/selected_sequence_duplicate_audit.tsv.gz",pd.DataFrame(duplicates,columns=dup_columns))
    write(VERSION/"data/configuration_registry.tsv",pd.DataFrame(registry))
    results={"class_counts":counts,"role_class_counts":role_counts,"configuration_overlap":overlaps,"sequence_context_qc":sequence_rows,"covariate_distributions":distributions,"lobe_provenance_balance":lobes,"matched_pair_distance_summary":pair_summaries,"pure_label_membership_change":label_changes}
    for name,rows in results.items():
        write(VERSION/"results"/f"{name}.tsv",pd.DataFrame(rows))
    retention=read(ATTEMPT/"positive_retention.tsv")
    retention["interpretation"]="full target population retained; distributional matching, not a local-caliper common-support claim"
    write(VERSION/"results/positive_retention.tsv",retention)
    runs=[];projection=[]
    for reg in registry:
        train=[row for row in role_counts if row["configuration"]==reg["configuration"] and row["model"]==reg["model"] and row["role"]=="train"][0]
        train_n=train["positive"]+train["control"]
        baseline_n,baseline_seconds=(464262,2659) if reg["model"]=="enhancer" else (78165,454)
        for seed in SEEDS:
            run={"run_id":f"{reg['configuration']}_{reg['model']}_seed{seed}","configuration":reg["configuration"],"model":reg["model"],"seed":seed,"training_intervals":train_n,"manifest":reg["manifest"],"manifest_sha256":reg["sha256"],"epochs_max":50,"batch_size":256,"patience":15,"learning_rate":.001,"optimizer":"Adadelta","gpu_type":"a100","gpus":1,"cpus":8,"memory_GiB":48,"walltime_cap_hours":4,"status":"NOT_AUTHORIZED_NOT_STARTED"}
            runs.append(run)
            projection.append({"run_id":run["run_id"],"configuration":reg["configuration"],"model":reg["model"],"seed":seed,"training_intervals":train_n,"historical_reference_training_intervals":baseline_n,"historical_reference_fit_seconds":baseline_seconds,"baseline_scaled_fit_GPU_hours":baseline_seconds*train_n/baseline_n/3600,"gpus":1,"cpus":8,"memory_GiB":48,"walltime_cap_hours":4,"status":"NOT_AUTHORIZED_NOT_STARTED"})
    write(VERSION/"data/run_matrix.tsv",pd.DataFrame(runs))
    write(VERSION/"results/compute_projection_by_run.tsv",pd.DataFrame(projection))
    normalization=json.loads((VERSION/"provenance/atac_normalization_manifest.json").read_text())
    unique_sequences=int(selected.canonical_rc_sequence_sha256.nunique())
    cache_bytes=unique_sequences*36480
    compute={"version":"pretraining-1.1","training_authorized":False,"training_started":False,"model_inference_performed":False,"readiness":"AWAITING_INDEPENDENT_FINAL_VALIDATION_AND_INVESTIGATOR_REVIEW","n_fits":18,"training_seeds":SEEDS,"scaled_fit_only_GPU_hours":sum(row["baseline_scaled_fit_GPU_hours"] for row in projection),"scaling_method":"historical fit seconds * current training interval count / historical training interval count; enhancer 2659 s / 464262 intervals, H3K27me3-associated 454 s / 78165 intervals","conservative_total_GPU_hours":[12,24],"cost_uncertainty":"Count scaling is not a measured V2 fit; frozen phase-I two-orientation cache extraction, symmetric validation, I/O and implementation overhead remain unmeasured and are covered by a conservative unchanged planning range, not a guarantee.","estimated_execution_walltime_hours_excluding_queue":[4,8],"one_A100_per_fit":True,"future_per_fit_CPUs":8,"future_per_fit_RAM_GiB":48,"max_concurrent_fits":4,"sum_requested_walltime_caps_GPU_hours":72,"selected_unique_intervals":len(active),"unique_selected_encoded_sequences":unique_sequences,"future_two_orientation_float32_cache_bytes_per_encoded_sequence":36480,"future_feature_cache_bytes_for_all_unique_selected_sequences_both_orientations":cache_bytes,"future_feature_cache_GiB":cache_bytes/2**30,"cache_release_rule":"Construct/cache train and chr7 only after separate authorization; defer chr8-9 model execution until final retained-model freeze. No representation has been extracted in this preflight.","local_scratch_budget_GiB":120,"CPU_preflight_reproduction_request":{"CPUs":8,"RAM_GiB":24,"walltime_hours":1},"actual_CPU_ATAC_normalization_seconds":normalization["wall_seconds"],"actual_CPU_matching_attempt001_seconds":attempt["elapsed_seconds"],"actual_CPU_assembly_seconds":time.time()-started,"actual_CPU_assembly_peak_RSS_KiB":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    write_json(VERSION/"provenance/compute_plan.json",compute)
    manifest={"version":"pretraining-1.1","completed_utc":dt.datetime.now(dt.timezone.utc).isoformat(),"assembly_only_no_new_matching":True,"all_original_feature_strings_preserved":True,"rich_manifest_columns":rich_columns,"all_master_intervals":len(frame),"selected_unique_intervals":len(active),"unique_selected_encoded_sequences":unique_sequences,"future_two_orientation_feature_cache_bytes":cache_bytes,"duplicate_groups":duplicate_groups,"cross_partition_selected_duplicate_groups":cross_partition_groups,"cross_validation_role_selected_duplicate_groups":cross_role_groups,"configuration_registry":registry,"common_C_panels":panel_infos,"common_panel_roles":["selection","calibration","test"],"checkpoint_role_excluded_from_common_panels":True,"V2_B_complete_V1_positive_population":True,"B_C_unique_one_to_one_controls":True,"B_C_independent_rematching":True,"B_C_interpretation":"same-lobe label correction plus induced control rematching","local_pair_calipers_removed":True,"pair_distance_summary":str((VERSION/"results/matched_pair_distance_summary.tsv").relative_to(ROOT)),"no_local_common_support_or_pairwise_exchangeability_claim":True,"attempt_rule_sha256":sha(VERSION/"specification/matching_attempt_001.json"),"inputs":inputs,"training_started":False,"model_inference_performed":False,"phase_I_feature_extraction_performed":False,"benchmark_outcomes_read":False,"status":"ASSEMBLED_PENDING_INDEPENDENT_FINAL_VALIDATION","elapsed_seconds":time.time()-started,"peak_RSS_KiB":resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
    write_json(VERSION/"provenance/configuration_manifest.json",manifest)
    print(json.dumps({"status":manifest["status"],"selected_unique_intervals":len(active),"unique_encoded_sequences":unique_sequences,"scaled_fit_GPU_hours":compute["scaled_fit_only_GPU_hours"],"cache_GiB":compute["future_feature_cache_GiB"],"elapsed_seconds":manifest["elapsed_seconds"]}),flush=True)


if __name__=="__main__":
    main()
