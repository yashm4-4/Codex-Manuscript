#!/usr/bin/env python3
"""Plot the prespecified frozen COPD reverse-complement diagnostic only."""
from __future__ import annotations
import os, sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[4]
V2=ROOT/"diseases/COPD/07_gap_closure"
os.environ["MPLCONFIGDIR"]=str(V2/"logs/matplotlib")
os.environ.setdefault("PYTHONDONTWRITEBYTECODE","1")
sys.dont_write_bytecode=True
import hashlib,json,platform
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.stats import pearsonr,spearmanr
EXPECTED_SPEC="09d3656f3e06c24fa5402c28b976d6be06fff91604ad0b7b29c3150120d76ac1"
R=V2/"results"
S4=ROOT/"diseases/COPD/04_modeling/results"
MODELS=("enhancer","silencer")
LABELS={"enhancer":"Enhancer","silencer":"H3K27me3-associated"}
COLORS={"SNV":"#0072B2","indel_or_complex":"#D55E00"}
plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,"pdf.fonttype":42,"ps.fonttype":42})
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return pd.read_csv(p,sep="\t",float_precision="round_trip")
def bools(s): return s.astype(str).str.lower().isin(["true","1"])
def main():
    spec=V2/"provenance/COPD-V2-RC_analysis_specification.md"
    assert sha(spec)==EXPECTED_SPEC
    fpath=S4/"COPD-S4-R003_candidate_allele_scores.tsv.gz"
    rpath=R/"COPD-V2-RC-R001_reverse_complement_scores.tsv.gz"
    mpath=S4/"COPD-S4-R004_prioritized_candidates.tsv.gz"
    tpath=S4/"COPD-S4-R004_delta_thresholds.tsv"
    f=read(fpath).set_index("candidate_record_id")
    rc=read(rpath).set_index("candidate_record_id")
    assert len(f)==15303 and f.index.equals(rc.index) and f.index.is_unique
    meta=read(mpath).set_index("candidate_record_id").loc[f.index]
    thresholds=read(tpath)
    cls=meta.variant_class_group
    eligible=~bools(meta.encode_blacklist)
    data={}
    for model in MODELS:
        t=thresholds[thresholds.model_type.eq(model)]
        T=float(t.region_score_cutoff.iloc[0])
        D=cls.map(t.set_index("variant_class_group").abs_delta_cutoff).to_numpy()
        ff=f[[model+"_ref_score",model+"_alt_score"]].to_numpy()
        rr=rc[[model+"_ref_score",model+"_alt_score"]].to_numpy()
        fs=ff.max(axis=1);rs=rr.max(axis=1)
        fd=ff[:,1]-ff[:,0];rd=rr[:,1]-rr[:,0]
        fc=eligible.to_numpy()&(fs>=T)&(abs(fd)>=D)
        cr=eligible.to_numpy()&(rs>=T)&(abs(rd)>=D)
        data[model]=dict(T=T,D=D,ff=ff,rr=rr,fs=fs,rs=rs,fd=fd,rd=rd,fc=fc,rc=cr)
    outputs=[]
    def save(fig,name):
        for ext in ("png","pdf"):
            p=R/(name+"."+ext)
            fig.savefig(p,dpi=220,bbox_inches="tight")
            outputs.append({"path":str(p.relative_to(ROOT)),"bytes":p.stat().st_size,"sha256":sha(p)})
        plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(14,8),layout="constrained")
    for row,model in enumerate(MODELS):
        a=data[model]
        for col,kind in enumerate(("REF","ALT","Region = max(REF, ALT)")):
            x=a["ff"][:,col] if col<2 else a["fs"]
            y=a["rr"][:,col] if col<2 else a["rs"]
            ax=axes[row,col];ax.hexbin(x,y,gridsize=55,extent=(0,1,0,1),mincnt=1,bins="log",cmap="viridis")
            ax.plot([0,1],[0,1],"--",color="#777777",lw=.8)
            ax.axvline(a["T"],color="#D55E00",lw=.7,alpha=.7);ax.axhline(a["T"],color="#D55E00",lw=.7,alpha=.7)
            ax.set(xlim=(0,1),ylim=(0,1),xlabel="Frozen forward score",ylabel="Reverse-complement score",title=f"{LABELS[model]}: {kind}")
            ax.text(.04,.95,f"r = {pearsonr(x,y).statistic:.3f}\nρ = {spearmanr(x,y).statistic:.3f}\nMAE = {np.mean(abs(y-x)):.4f}",transform=ax.transAxes,va="top",bbox=dict(facecolor="white",alpha=.9,edgecolor="none"))
    fig.suptitle("COPD-V2-RC-F001 | All 15,303 frozen allele pairs\nDensity on log scale; orange lines mark the frozen region-score threshold",fontsize=13)
    save(fig,"COPD-V2-RC-F001_score_concordance")
    fig,axes=plt.subplots(2,2,figsize=(12,9),layout="constrained")
    for row,model in enumerate(MODELS):
        a=data[model];x=a["fd"];y=a["rd"];ax=axes[row,0]
        limit=max(abs(x).max(),abs(y).max())*1.03
        ax.hexbin(x,y,gridsize=65,mincnt=1,bins="log",cmap="viridis",extent=(-limit,limit,-limit,limit))
        ax.plot([-limit,limit],[-limit,limit],"--",color="grey",lw=.8)
        ax.axhline(0,color="grey",lw=.5);ax.axvline(0,color="grey",lw=.5)
        ax.set(xlabel="Frozen forward ALT − REF",ylabel="RC ALT − REF",title=LABELS[model],xlim=(-limit,limit),ylim=(-limit,limit))
        ax.text(.04,.96,f"r = {pearsonr(x,y).statistic:.3f}\nρ = {spearmanr(x,y).statistic:.3f}",transform=ax.transAxes,va="top",bbox=dict(facecolor="white",alpha=.9,edgecolor="none"))
        ax=axes[row,1]
        for group,color in COLORS.items():
            vals=np.sort((abs(y-x)/a["D"])[cls.eq(group)])
            ax.plot(vals,np.arange(1,len(vals)+1)/len(vals),label=f"{group} (n={len(vals):,})",color=color)
        ax.set_xscale("symlog",linthresh=.01)
        ax.axvline(1,color="grey",ls="--",lw=.8)
        ax.set(xlabel="|RC delta − forward delta| / frozen class cutoff",ylabel="Cumulative fraction",ylim=(0,1.02),title="Delta disagreement by variant class")
        ax.legend(loc="lower right",fontsize=9)
    fig.suptitle("COPD-V2-RC-F002 | Allele effects and class-normalized disagreement\nDelta sign describes model activity, not disease-risk direction",fontsize=13)
    save(fig,"COPD-V2-RC-F002_delta_concordance")
    fig,axes=plt.subplots(2,3,figsize=(16,9),layout="constrained")
    for row,model in enumerate(MODELS):
        a=data[model];err=abs(a["rs"]-a["fs"]);T=a["T"]
        ax=axes[row,0];ax.hexbin(a["fs"],err,gridsize=50,mincnt=1,bins="log",cmap="viridis")
        ax.axvspan(T-.02,T+.02,color="#E69F00",alpha=.25)
        ax.set(xlabel="Frozen forward region score",ylabel="|RC − forward region score|",title=LABELS[model])
        ax=axes[row,1]
        bins=np.array([0,T-.10,T-.02,T+.02,T+.10,1.000001])
        means=[];p95=[];counts=[]
        for low,high in zip(bins[:-1],bins[1:]):
            vals=err[(a["fs"]>=low)&(a["fs"]<high)]
            means.append(np.mean(vals) if len(vals) else np.nan)
            p95.append(np.quantile(vals,.95) if len(vals) else np.nan)
            counts.append(len(vals))
        x=np.arange(5);ax.bar(x,means,color="#0072B2",alpha=.8,label="Mean")
        ax.plot(x,p95,"o-",color="#D55E00",label="95th percentile")
        ax.set_xticks(x,[f"{lo:.3f}–{min(hi,1):.3f}\nn={n:,}" for lo,hi,n in zip(bins[:-1],bins[1:],counts)],rotation=35,ha="right",fontsize=8)
        ax.set(ylabel="Absolute region-score difference",title="Disagreement across score bins");ax.legend(fontsize=9)
        ax=axes[row,2];fcall=a["fc"];lost=fcall&~a["rc"];stable=fcall&a["rc"]
        xm=a["fs"]-T;ym=abs(a["fd"])/a["D"]
        ax.scatter(xm[stable],ym[stable],s=22,alpha=.7,label=f"Retained (n={stable.sum()})",color="#009E73")
        ax.scatter(xm[lost],ym[lost],s=32,alpha=.85,marker="x",label=f"Lost (n={lost.sum()})",color="#D55E00")
        ax.axvspan(0,.02,color="#E69F00",alpha=.13);ax.axhspan(1,1.1,color="#E69F00",alpha=.13)
        ax.axvline(.05,color="grey",ls=":",lw=.8);ax.axhline(1.5,color="grey",ls=":",lw=.8)
        ax.set(xlabel="Forward region score − threshold",ylabel="|Forward delta| / frozen cutoff",title="Original positive calls: fixed margins")
        ax.legend(fontsize=8)
    fig.suptitle("COPD-V2-RC-F003 | Disagreement and losses relative to prespecified thresholds\nOrange: borderline bands; dotted margins: well-separated calls",fontsize=13)
    save(fig,"COPD-V2-RC-F003_threshold_proximity")
    fig,axes=plt.subplots(2,2,figsize=(11,10),layout="constrained")
    for ax,label,ff,rr in [
        (axes[0,0],LABELS["enhancer"],data["enhancer"]["fc"],data["enhancer"]["rc"]),
        (axes[0,1],LABELS["silencer"],data["silencer"]["fc"],data["silencer"]["rc"]),
        (axes[1,0],"Union",data["enhancer"]["fc"]|data["silencer"]["fc"],data["enhancer"]["rc"]|data["silencer"]["rc"])]:
        matrix=np.array([[(~ff&~rr).sum(),(~ff&rr).sum()],[(ff&~rr).sum(),(ff&rr).sum()]])
        ax.imshow(matrix,cmap="Blues",norm=LogNorm(vmin=1,vmax=15303))
        for i in range(2):
            for j in range(2): ax.text(j,i,f"{matrix[i,j]:,}",ha="center",va="center",color="white" if matrix[i,j]>500 else "black",fontsize=15)
        both=int((ff&rr).sum());u=int((ff|rr).sum())
        ax.set(xticks=[0,1],yticks=[0,1],xticklabels=["Negative","Positive"],yticklabels=["Negative","Positive"],xlabel="Reverse complement",ylabel="Frozen forward",title=f"{label}: Jaccard={both/u:.3f}")
    ax=axes[1,1]
    fc=data["enhancer"]["fc"].astype(int)+2*data["silencer"]["fc"].astype(int)
    rr=data["enhancer"]["rc"].astype(int)+2*data["silencer"]["rc"].astype(int)
    mat=np.array([[( (fc==i)&(rr==j)).sum() for j in range(4)] for i in range(4)])
    ax.imshow(mat,cmap="Blues",norm=LogNorm(vmin=1,vmax=15303))
    for i in range(4):
        for j in range(4): ax.text(j,i,f"{mat[i,j]:,}",ha="center",va="center",color="white" if mat[i,j]>500 else "black",fontsize=10)
    names=["Neither","Enhancer","H3K27me3","Both"]
    ax.set(xticks=range(4),yticks=range(4),xticklabels=names,yticklabels=names,xlabel="RC model support",ylabel="Forward model support",title="Full model-support contexts")
    plt.setp(ax.get_xticklabels(),rotation=30,ha="right")
    fig.suptitle("COPD-V2-RC-F004 | Frozen-rule candidate transitions\nAll 15,303 pairs; blacklist exclusions preserved; no replacement candidate set",fontsize=13)
    save(fig,"COPD-V2-RC-F004_call_transitions")
    manifest={"result_id":"COPD-V2-RC-FIGURES","spec_sha256":EXPECTED_SPEC,"python":platform.python_version(),"matplotlib":matplotlib.__version__,"inputs":[{"path":str(p.relative_to(ROOT)),"sha256":sha(p)} for p in (fpath,rpath,mpath,tpath,Path(__file__))],"outputs":outputs}
    (V2/"provenance/COPD-V2-RC_figure_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(f"Wrote {len(outputs)} figure files",flush=True)
if __name__=="__main__":main()
