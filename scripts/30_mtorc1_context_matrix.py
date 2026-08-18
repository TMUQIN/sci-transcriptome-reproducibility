#!/usr/bin/env python3
"""Assemble an auditable mTORC1 context/method matrix and data figure."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


TERM = "HALLMARK_MTORC1_SIGNALING"


def sha256(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""): h.update(b)
    return h.hexdigest()


def cls(p: float|None, fdr: float|None, status: str="tested") -> str:
    if status != "tested": return "not_evaluable"
    if fdr is not None and np.isfinite(fdr) and fdr <= .05: return "FDR_supported_association_or_enrichment"
    if p is not None and np.isfinite(p) and p <= .05: return "nominal_only_not_replication"
    return "directional_or_no_statistical_support"


def build(root:Path)->tuple[pd.DataFrame,list[Path]]:
    rows=[]; inputs=[]
    gsea_path=root/"results/whole_lesion_programs/hallmark_gsea_by_study.tsv"; inputs.append(gsea_path)
    gsea=pd.read_csv(gsea_path,sep="\t"); gsea=gsea[gsea.term.eq(TERM)]
    for r in gsea.itertuples():
        rows.append({"source_layer":"primary_three_study_whole_lesion","dataset":r.dataset,
            "tissue_scope":"whole_lesion","cell_type":"all_cells","genotype_context":"source_cohort",
            "contrast_id":r.contrast_id,"scoring_method":"preranked_GSEA","effect_metric":"NES",
            "estimate":r.nes,"se":np.nan,"ci_low":np.nan,"ci_high":np.nan,"nominal_p":r.nominal_p,"fdr":r.fdr,
            "direction":int(np.sign(r.nes)),"evidence_class":cls(r.nominal_p,r.fdr),
            "interpretation_limit":"rank-based transcript enrichment; not pathway activation"})
    eff_path=root/"reports/phase_reproducibility_calibration_2026_07/pathway_effects_by_study.tsv"; inputs.append(eff_path)
    eff=pd.read_csv(eff_path,sep="\t"); eff=eff[eff.program.eq(TERM)]
    for r in eff.itertuples():
        rows.append({"source_layer":"biological_sample_score","dataset":r.dataset,"tissue_scope":"whole_lesion",
            "cell_type":"all_cells","genotype_context":"source_cohort","contrast_id":r.contrast_id,
            "scoring_method":r.method,"effect_metric":"standardized_score_contrast","estimate":r.estimate,"se":r.se,
            "ci_low":r.ci_low,"ci_high":r.ci_high,"nominal_p":r.p,"fdr":r.fdr,"direction":int(np.sign(r.estimate)),
            "evidence_class":cls(r.p,r.fdr),"interpretation_limit":"within-study score association; FDR is within-study Hallmark family"})
    meta_path=root/"reports/phase_reproducibility_calibration_2026_07/pathway_meta_mkh.tsv"; inputs.append(meta_path)
    meta=pd.read_csv(meta_path,sep="\t"); meta=meta[meta.program.eq(TERM)]
    for r in meta.itertuples():
        rows.append({"source_layer":"three_study_mKH_meta","dataset":"three_study_meta","tissue_scope":"whole_lesion",
            "cell_type":"all_cells","genotype_context":"source_cohorts","contrast_id":r.contrast_id,
            "scoring_method":r.method,"effect_metric":"mKH_meta_standardized_effect","estimate":r.estimate,"se":r.se_hksj,
            "ci_low":r.ci_low,"ci_high":r.ci_high,"nominal_p":r.p,"fdr":r.meta_fdr,"direction":int(np.sign(r.estimate)),
            "evidence_class":cls(r.p,r.meta_fdr),
            "interpretation_limit":"k=3 modified Hartung-Knapp; wide small-k uncertainty retained"})
    ge_path=root/"reports/phase_reproducibility_calibration_2026_07/gsema_benchmark/gsema_effects_by_study.tsv"; inputs.append(ge_path)
    ge=pd.read_csv(ge_path,sep="\t"); ge=ge[ge.pathway.eq(TERM)]
    for r in ge.itertuples():
        rows.append({"source_layer":"GSEMA_benchmark_study_effect","dataset":r.dataset,"tissue_scope":"whole_lesion",
            "cell_type":"all_cells","genotype_context":"source_cohort","contrast_id":r.contrast_id,
            "scoring_method":"GSEMA_"+r.method,"effect_metric":"GSEMA_standardized_effect","estimate":r.estimate,
            "se":r.se,"ci_low":r.estimate-1.96*r.se,"ci_high":r.estimate+1.96*r.se,"nominal_p":np.nan,"fdr":np.nan,
            "direction":int(np.sign(r.estimate)),"evidence_class":"effect_only_no_study_level_test",
            "interpretation_limit":"benchmark score effect; inference comes from separate mKH meta row"})
    gm_path=root/"reports/phase_reproducibility_calibration_2026_07/gsema_benchmark/gsema_mkh_random_effects.tsv"; inputs.append(gm_path)
    gm=pd.read_csv(gm_path,sep="\t"); gm=gm[gm.pathway.eq(TERM)]
    for r in gm.itertuples():
        rows.append({"source_layer":"GSEMA_benchmark_mKH_meta","dataset":"three_study_meta","tissue_scope":"whole_lesion",
            "cell_type":"all_cells","genotype_context":"source_cohorts","contrast_id":r.contrast_id,
            "scoring_method":"GSEMA_"+r.method,"effect_metric":"mKH_meta_standardized_effect","estimate":r.estimate,
            "se":r.se_hksj,"ci_low":r.ci_low,"ci_high":r.ci_high,"nominal_p":r.p,"fdr":r.meta_fdr,
            "direction":int(np.sign(r.estimate)),"evidence_class":cls(r.p,r.meta_fdr),
            "interpretation_limit":"GSEMA score with project mKH; k=3 and method dependence retained"})
    rs_path=root/"reports/phase_reproducibility_calibration_2026_07/gsema_benchmark/gsema_run_status.tsv"; inputs.append(rs_path)
    rs=pd.read_csv(rs_path,sep="\t")
    for r in rs[rs.status.ne("completed")].itertuples():
        rows.append({"source_layer":"GSEMA_benchmark_failure","dataset":"three_study_meta","tissue_scope":"whole_lesion",
            "cell_type":"all_cells","genotype_context":"source_cohorts","contrast_id":r.contrast_id,
            "scoring_method":"GSEMA_"+r.method,"effect_metric":"not_evaluable","estimate":np.nan,"se":np.nan,
            "ci_low":np.nan,"ci_high":np.nan,"nominal_p":np.nan,"fdr":np.nan,"direction":0,
            "evidence_class":"not_evaluable","interpretation_limit":str(r.error).replace("\n"," ")[:500]})
    ext_path=root/"results/GSE304361_support/programs/GSE304361_frozen_program_results.tsv"; inputs.append(ext_path)
    ext=pd.read_csv(ext_path,sep="\t"); ext=ext[ext.term.eq(TERM)]
    for r in ext.itertuples():
        status="tested" if r.gsea_status=="tested" else "not_evaluable"
        context=("wild_type_injury" if r.contrast_id=="wt_d7_vs_wt_uninjured" else
                 "Plxnb1_KO_injury" if r.contrast_id=="ko_d7_vs_ko_uninjured" else
                 "Plxnb1_genotype_by_injury" if "interaction" in r.contrast_id else "genotype_contrast")
        rows.append({"source_layer":"GSE304361_endpoint_stress_test","dataset":"GSE304361",
            "tissue_scope":r.analysis_scope,"cell_type":r.cell_state,"genotype_context":context,
            "contrast_id":r.contrast_id,"scoring_method":"preranked_GSEA","effect_metric":"NES",
            "estimate":r.nes,"se":np.nan,"ci_low":np.nan,"ci_high":np.nan,"nominal_p":r.nominal_p,
            "fdr":r.fdr,"direction":int(np.sign(r.nes)) if np.isfinite(r.nes) else 0,
            "evidence_class":cls(r.nominal_p,r.fdr,status),
            "interpretation_limit":"one study family; context stress test, not independent pathway mechanism"})
    comp_path=root/"reports/phase_final_calibration_and_v5_2026_07/composition_program_association.tsv"; inputs.append(comp_path)
    comp=pd.read_csv(comp_path,sep="\t"); comp=comp[(comp.program.eq(TERM))&(comp.analysis_type.eq("composition_PC1_adjusted_group_model"))]
    for r in comp.itertuples():
        rows.append({"source_layer":"composition_PC1_sensitivity","dataset":r.dataset,"tissue_scope":"whole_lesion",
            "cell_type":"all_cells","genotype_context":"source_cohort","contrast_id":r.contrast_id,
            "scoring_method":"centered_rank_mean_plus_CLR_PC1","effect_metric":"adjusted_standardized_score_contrast",
            "estimate":r.estimate,"se":r.se,"ci_low":r.ci_low,"ci_high":r.ci_high,"nominal_p":r.nominal_p,
            "fdr":r.nominal_fdr_within_analysis_type,"direction":int(np.sign(r.estimate)) if np.isfinite(r.estimate) else 0,
            "evidence_class":cls(r.nominal_p,r.nominal_fdr_within_analysis_type,"tested" if r.status=="evaluable_exploratory" else "not_evaluable"),
            "interpretation_limit":"one-PC sensitivity; not proof of composition independence"})
    return pd.DataFrame(rows),inputs


def plot_context(matrix:pd.DataFrame,path:Path)->None:
    contrasts=["injury_d1_vs_uninjured","injury_d7_vs_uninjured","change_d7_minus_d1"]
    fig,axes=plt.subplots(1,4,figsize=(15,7),gridspec_kw={"width_ratios":[1,1,1,1.25]},constrained_layout=True)
    methods={"centered_rank_mean":("o","#0072B2"),"magnitude_mean_z":("s","#D55E00")}
    for ax,contrast in zip(axes[:3],contrasts):
        d=matrix[(matrix.source_layer=="biological_sample_score")&(matrix.contrast_id==contrast)].copy()
        d["label"]=d.dataset+" | "+d.scoring_method.str.replace("_mean","",regex=False)
        d=d.sort_values(["dataset","scoring_method"]); y=np.arange(len(d))[::-1]
        for yi,r in zip(y,d.itertuples()):
            marker,color=methods.get(r.scoring_method,("o","black"))
            ax.errorbar(r.estimate,yi,xerr=[[r.estimate-r.ci_low],[r.ci_high-r.estimate]],fmt=marker,color=color,capsize=2,ms=5)
        m=matrix[(matrix.source_layer=="three_study_mKH_meta")&(matrix.contrast_id==contrast)]
        for j,r in enumerate(m.itertuples()):
            yi=-1-j; marker,color=methods.get(r.scoring_method,("D","black"))
            ax.errorbar(r.estimate,yi,xerr=[[r.estimate-r.ci_low],[r.ci_high-r.estimate]],fmt="D",color=color,capsize=3,ms=6)
            d.loc[len(d),"label"]="mKH | "+r.scoring_method.replace("_mean","")
            y=np.r_[y,yi]
        ax.axvline(0,color="0.5",lw=.8); ax.set_yticks(y); ax.set_yticklabels(d.label,fontsize=7)
        ax.set_title(contrast.replace("injury_","").replace("_vs_uninjured"," vs U").replace("change_",""),fontsize=9)
        ax.set_xlabel("Standardized score effect (95% CI)",fontsize=8); ax.grid(axis="x",alpha=.2)
    ax=axes[3]
    ext=matrix[(matrix.source_layer=="GSE304361_endpoint_stress_test") & matrix.contrast_id.isin(
        ["wt_d7_vs_wt_uninjured","plxnb1_genotype_by_injury_interaction"]) & matrix.cell_type.isin(
        ["whole_lesion_bulk_equivalent","Astrocyte"])].copy()
    ext["label"]=ext.genotype_context+" | "+ext.cell_type
    y=np.arange(len(ext))[::-1]
    colors=np.where(ext.fdr.le(.05),"#CC79A7","#999999")
    for yi,(_,r),c in zip(y,ext.iterrows(),colors): ax.plot(r.estimate,yi,"o",color=c,ms=7)
    ax.axvline(0,color="0.5",lw=.8); ax.set_yticks(y); ax.set_yticklabels(ext.label,fontsize=7)
    ax.set_xlabel("GSE304361 preranked GSEA NES\n(points only; no comparable CI)",fontsize=8)
    ax.set_title("Context stress test",fontsize=9); ax.grid(axis="x",alpha=.2)
    fig.suptitle("mTORC1 transcriptional evidence varies by contrast, score, cell type and genotype",fontsize=11)
    fig.savefig(path,bbox_inches="tight"); plt.close(fig)


def write_md(matrix:pd.DataFrame,path:Path)->None:
    meta=matrix[matrix.source_layer.isin(["three_study_mKH_meta","GSEMA_benchmark_mKH_meta"])]
    n_meta=len(meta); n_meta_fdr=int(meta.fdr.le(.05).sum())
    wt=matrix[(matrix.source_layer=="GSE304361_endpoint_stress_test")&(matrix.contrast_id=="wt_d7_vs_wt_uninjured")]
    astro=wt[wt.cell_type.eq("Astrocyte")].iloc[0]
    inter=matrix[(matrix.source_layer=="GSE304361_endpoint_stress_test")&
                 (matrix.contrast_id=="plxnb1_genotype_by_injury_interaction")&
                 matrix.cell_type.isin(["whole_lesion_bulk_equivalent","Astrocyte"])]
    text=f"""# mTORC1 evidence classification

## Frozen classification

**mTORC1 is a context-dependent and method-sensitive injury-associated transcriptional program.**

This is an enrichment/score classification, not a claim of biochemical activation, pathway causality, therapeutic validation, or direct regulation by Plexin-B1.

## Evidence synthesis

- Preranked GSEA in the three discovery cohorts gave the expected positive day-1/day-7 and negative direct-delta directions, but sample-level scoring did not reproduce the same inferential strength across methods.
- Across {n_meta} mKH rows from centered-rank, magnitude-aware, GSVA, ssGSEA and Zscore implementations, {n_meta_fdr} reached the corresponding Hallmark meta-FDR threshold. Method disagreement is therefore a result, not a nuisance to suppress.
- The GSE304361 wild-type astrocyte endpoint was negative (NES {astro.estimate:.3f}, FDR {astro.fdr:.4g}) while the whole-lesion endpoint was weakly positive; this rejects a cross-lineage universal-positive description.
- The Plexin-B1 genotype-by-injury contrast was negative in whole lesion and astrocytes ({'; '.join(f'{r.cell_type}: NES {r.estimate:.3f}, FDR {r.fdr:.4g}' for r in inter.itertuples())}). These are within-family context modifications and do not establish a Plexin-B1→mTORC1 causal path.
- One-PC composition sensitivity retained most discovery-cohort directions but changed the GSE234774 day-7 estimate to approximately null/slightly negative. This further supports context sensitivity.

## Prohibited wording

- universally activated;
- validated therapeutic target;
- Plexin-B1 directly regulates mTORC1;
- mTORC1 mediates Plexin-B1 effects;
- replicated when only nominal P or direction is available.

## Manuscript use

Use mTORC1 as the principal falsification case demonstrating that a Hallmark label can appear coherent under preranked enrichment yet change with biological-sample scoring, cellular context, genotype background and composition adjustment. Preserve all negative and non-significant rows in the main context matrix.
"""
    path.write_text(text,encoding="utf-8")


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument("--root",type=Path,default=Path(__file__).resolve().parents[1])
    ap.add_argument("--out-dir",type=Path,default=Path("reports/phase_final_calibration_and_v5_2026_07"))
    args=ap.parse_args(); root=args.root.resolve(); out=args.out_dir if args.out_dir.is_absolute() else root/args.out_dir
    matrix,inputs=build(root); matrix=matrix.sort_values(["source_layer","contrast_id","dataset","cell_type","scoring_method"])
    matrix_path=out/"mtorc1_context_matrix.tsv"; matrix.to_csv(matrix_path,sep="\t",index=False,na_rep="NA")
    md_path=out/"mtorc1_evidence_classification.md"; write_md(matrix,md_path)
    fig_path=out/"figure_mtorc1_context_forest.pdf"; plot_context(matrix,fig_path)
    outputs=[matrix_path,md_path,fig_path]
    prov={"analysis":"mTORC1 context matrix","created_at":datetime.now().astimezone().isoformat(),
          "classification":"context-dependent and method-sensitive injury-associated transcriptional program",
          "boundaries":["enrichment is not activation","no causal Plexin-B1-to-mTORC1 claim","nominal P is not replication",
                        "incompatible contexts are displayed, not forced into one meta-analysis"],
          "inputs":[{"path":str(p.relative_to(root)),"bytes":p.stat().st_size,"sha256":sha256(p)} for p in inputs],
          "outputs":[{"path":str(p.relative_to(root)),"bytes":p.stat().st_size,"sha256":sha256(p)} for p in outputs]}
    (out/"mtorc1_context_provenance.json").write_text(json.dumps(prov,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":"ok","matrix_rows":len(matrix),"figure_bytes":fig_path.stat().st_size}))


if __name__=="__main__": main()
