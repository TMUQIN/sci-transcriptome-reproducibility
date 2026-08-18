#!/usr/bin/env python3
"""Composition and ontology sensitivity with frozen small-sample safeguards.

Only GSE162610 and GSE234774 have usable sample-by-state cell counts. GSE304399 is
carried through every requested output as not evaluable. Adjustment uses one CLR-PC,
never more, and must pass rank, residual-df, VIF and condition-number gates.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import rankdata, spearmanr, t as student_t
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor


LABELED = ("GSE162610", "GSE234774")
ALL_STUDIES = (*LABELED, "GSE304399")
GROUPS = ("uninjured", "dpi_1", "dpi_7")
CONTRASTS = {
    "injury_d1_vs_uninjured": np.array([0, 0, 1, 0.]),
    "injury_d7_vs_uninjured": np.array([0, 0, 0, 1.]),
    "change_d7_minus_d1": np.array([0, 0, -1, 1.]),
}
FOCAL = (
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_HYPOXIA",
    "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_MTORC1_SIGNALING",
)
AMBIGUOUS = {"myeloid_broad", "leukocyte_broad", "lymphoid_broad"}
BROAD = {
    "microglia": "myeloid", "macrophage_monocyte": "myeloid", "myeloid_broad": "myeloid",
    "dendritic": "myeloid", "neutrophil": "myeloid", "leukocyte_broad": "myeloid",
    "t_cell": "lymphoid", "b_cell": "lymphoid", "lymphoid_broad": "lymphoid",
    "oligodendrocyte": "oligodendroglial", "opc": "oligodendroglial", "schwann": "oligodendroglial",
    "endothelial": "vascular", "pericyte_vsmc": "vascular",
    "fibroblast": "stromal", "meningeal_stromal": "stromal",
    "astrocyte": "astrocyte", "ependymal": "ependymal", "neuron": "neuron",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_gmt(path: Path) -> dict[str, set[str]]:
    out = {}
    with path.open(encoding="utf-8") as f:
        for line in f:
            x = line.rstrip("\n").split("\t")
            if len(x) >= 3: out[x[0]] = set(x[2:])
    return out


def bh(x: pd.Series) -> pd.Series:
    p = x.to_numpy(float); order = np.argsort(p)
    q = p[order]*len(p)/np.arange(1, len(p)+1); q = np.minimum.accumulate(q[::-1])[::-1]
    z = np.empty_like(q); z[order] = np.clip(q, 0, 1)
    return pd.Series(z, index=x.index)


def load_harmonized(root: Path, study: str) -> tuple[np.ndarray, list[str], pd.DataFrame, list[Path]]:
    base = root/"data_processed"/"harmonized"
    paths = [base/f"{study}_pseudobulk_counts.mtx.gz", base/f"{study}_pseudobulk_genes.tsv",
             base/f"{study}_pseudobulk_coldata.tsv"]
    with gzip.open(paths[0], "rb") as f: matrix = mmread(f)
    matrix = matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)
    genes = pd.read_csv(paths[1], sep="\t").gene.astype(str).tolist()
    meta = pd.read_csv(paths[2], sep="\t")
    keep = meta.group.isin(GROUPS).to_numpy(); matrix = matrix[keep]; meta = meta.loc[keep].reset_index(drop=True)
    if matrix.shape != (len(meta), len(genes)): raise ValueError(f"{study}: matrix/metadata mismatch")
    if set(meta.cell_state) - set(BROAD): raise ValueError(f"{study}: unmapped strict lineages {set(meta.cell_state)-set(BROAD)}")
    return matrix.astype(float), genes, meta, paths


def composition_long(meta: pd.DataFrame, study: str) -> pd.DataFrame:
    samples = meta[["sample", "subject_id", "group"]].drop_duplicates()
    lineages = sorted(meta.cell_state.unique())
    grid = samples.assign(_=1).merge(pd.DataFrame({"lineage": lineages, "_": 1}), on="_").drop(columns="_")
    counts = meta.groupby(["sample", "subject_id", "group", "cell_state"], observed=True).n_cells.sum().rename("n_cells").reset_index().rename(columns={"cell_state":"lineage"})
    grid = grid.merge(counts, how="left", on=["sample", "subject_id", "group", "lineage"])
    grid["n_cells"] = grid.n_cells.fillna(0).astype(int)
    grid["sample_total_cells"] = grid.groupby("sample").n_cells.transform("sum")
    grid["proportion"] = grid.n_cells/grid.sample_total_cells
    grid.insert(0, "dataset", study); grid["composition_status"] = "evaluable_author_labels"
    return grid


def clr_pca(comp: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for study, d in comp[comp.composition_status.eq("evaluable_author_labels")].groupby("dataset", observed=True):
        wide = d.pivot(index=["sample", "subject_id", "group"], columns="lineage", values="n_cells").fillna(0)
        pseudo = wide.to_numpy(float)+.5; prop = pseudo/pseudo.sum(axis=1, keepdims=True)
        clr = np.log(prop)-np.log(prop).mean(axis=1, keepdims=True)
        centered = clr-clr.mean(axis=0, keepdims=True); u,s,vt = np.linalg.svd(centered, full_matrices=False)
        denom = max((s**2).sum(), np.finfo(float).eps); var = s**2/denom
        idx = wide.index.to_frame(index=False)
        for i, row in idx.iterrows():
            rows.append({"dataset":study, **row.to_dict(), "PC1":u[i,0]*s[0],
                         "PC2":u[i,1]*s[1] if len(s)>1 else np.nan,
                         "PC1_variance_ratio":var[0], "PC2_variance_ratio":var[1] if len(var)>1 else np.nan,
                         "n_lineages":wide.shape[1], "zero_replacement":"0.5 cell before closure",
                         "status":"evaluable"})
    return pd.DataFrame(rows)


def aggregate_ontology(matrix: np.ndarray, meta: pd.DataFrame, mode: str) -> tuple[np.ndarray, pd.DataFrame]:
    work = meta.copy(); keep = np.ones(len(work), dtype=bool)
    if mode == "strict": work["ontology_lineage"] = work.cell_state.astype(str)
    elif mode == "broad_family": work["ontology_lineage"] = work.cell_state.map(BROAD)
    elif mode == "exclude_ambiguous":
        keep = ~work.cell_state.isin(AMBIGUOUS).to_numpy(); work = work.loc[keep].reset_index(drop=True)
        matrix = matrix[keep]; work["ontology_lineage"] = work.cell_state.astype(str)
    else: raise ValueError(mode)
    keys = work["sample"].astype(str)+"\x1f"+work.ontology_lineage.astype(str)
    rows, mats = [], []
    for key in pd.unique(keys):
        idx = np.where(keys.to_numpy()==key)[0]; sample, lineage = key.split("\x1f")
        block = work.iloc[idx]; first = block.iloc[0]
        mats.append(matrix[idx].sum(axis=0))
        rows.append({"sample":sample, "subject_id":first.subject_id, "group":first.group,
                     "ontology_lineage":lineage, "n_cells":int(block.n_cells.sum()),
                     "source_lineages":";".join(sorted(block.cell_state.astype(str).unique()))})
    return np.vstack(mats), pd.DataFrame(rows)


def rank_scores(matrix: np.ndarray, genes: list[str], sets: dict[str,set[str]]) -> pd.DataFrame:
    lib = matrix.sum(axis=1); logcpm = np.log2(matrix/lib[:,None]*1e6+.5)
    ranks = np.apply_along_axis(rankdata, 1, logcpm, method="average")
    rankall = 2*(ranks-(matrix.shape[1]+1)/2)/max(matrix.shape[1]-1,1)
    gi = {g:i for i,g in enumerate(genes)}; rows=[]
    for term in FOCAL:
        idx = [gi[g] for g in sets[term] if g in gi]
        if len(idx)<10: continue
        for i,x in enumerate(rankall[:,idx].mean(axis=1)):
            rows.append({"row":i,"program":term,"score":x,"n_genes_observed":len(idx)})
    return pd.DataFrame(rows)


def contrast_stats(values: pd.DataFrame) -> list[dict]:
    group = values.group.astype(str)
    x = pd.DataFrame({"const":1.0, "PC1":0.0, "dpi_1":(group=="dpi_1").astype(float),
                      "dpi_7":(group=="dpi_7").astype(float)})
    # PC1 is a placeholder so the same frozen contrast vectors can be reused.
    x = x.drop(columns="PC1")
    model = sm.OLS(values.value.to_numpy(float), x.to_numpy(float)).fit()
    rows=[]
    vectors={"injury_d1_vs_uninjured":np.array([0,1,0.]),
             "injury_d7_vs_uninjured":np.array([0,0,1.]),
             "change_d7_minus_d1":np.array([0,-1,1.])}
    for name,v in vectors.items():
        est=float(v@model.params); se=float(math.sqrt(v@model.cov_params()@v)); crit=float(student_t.ppf(.975,model.df_resid))
        rows.append({"contrast_id":name,"estimate":est,"se":se,"df_residual":model.df_resid,
                     "ci_low":est-crit*se,"ci_high":est+crit*se,
                     "nominal_p":float(2*student_t.sf(abs(est/se),model.df_resid)) if se>0 else np.nan})
    return rows


def recomposition(root:Path, study:str, matrix:np.ndarray, genes:list[str], meta:pd.DataFrame,
                  sets:dict[str,set[str]], original:pd.DataFrame) -> pd.DataFrame:
    rows=[]
    for mode in ("strict","broad_family","exclude_ambiguous"):
        mat, md = aggregate_ontology(matrix,meta,mode); scored=rank_scores(mat,genes,sets)
        score = scored.merge(md.reset_index(names="row"),on="row",validate="many_to_one")
        samples=sorted(md.subject_id.unique())
        present=md.groupby("ontology_lineage",observed=True).subject_id.nunique()
        common=sorted(present[present.eq(len(samples))].index)
        if len(common)<2:
            for term in FOCAL:
                for contrast in CONTRASTS:
                    rows.append({"dataset":study,"ontology_mode":mode,"program":term,"contrast_id":contrast,
                                 "status":"not_evaluable","reason":"fewer than two lineages present in every eligible sample"})
            continue
        md_common=md[md.ontology_lineage.isin(common)].copy()
        md_common["sample_common_total"]=md_common.groupby("subject_id").n_cells.transform("sum")
        md_common["sample_proportion"]=md_common.n_cells/md_common.sample_common_total
        weights=md_common.groupby("ontology_lineage",observed=True).sample_proportion.mean(); weights=weights/weights.sum()
        score=score.merge(weights.rename("fixed_weight"),left_on="ontology_lineage",right_index=True,how="inner")
        composed=score.groupby(["subject_id","group","program"],observed=True).apply(
            lambda z: np.average(z.score,weights=z.fixed_weight),include_groups=False).rename("value").reset_index()
        composed["value"] = composed.groupby("program",observed=True).value.transform(
            lambda z:(z-z.mean())/z.std(ddof=1))
        for term, d in composed.groupby("program",observed=True):
            for stat in contrast_stats(d):
                orig=original[(original.dataset==study)&(original.method=="centered_rank_mean")&
                              (original.program==term)&(original.contrast_id==stat["contrast_id"])]
                orig_est=float(orig.estimate.iloc[0]) if len(orig) else np.nan
                rows.append({"dataset":study,"ontology_mode":mode,"program":term,**stat,
                             "original_centered_rank_estimate":orig_est,
                             "direction_matches_original":bool(np.sign(stat["estimate"])==np.sign(orig_est)),
                             "n_common_lineages":len(common),"common_lineages":";".join(common),
                             "fixed_weights":";".join(f"{k}:{weights[k]:.8g}" for k in common),
                             "n_biological_samples":d.subject_id.nunique(),"status":"evaluable_descriptive",
                             "reason":"fixed-weight recomposition over lineages observed in every eligible sample"})
    return pd.DataFrame(rows)


def associations(comp:pd.DataFrame,pca:pd.DataFrame,scores:pd.DataFrame,original:pd.DataFrame) -> pd.DataFrame:
    rows=[]
    focal=scores[(scores.program.isin(FOCAL))&(scores.method.eq("centered_rank_mean"))]
    for study in LABELED:
        sc=focal[focal.dataset.eq(study)][["subject_id","group","program","standardized_score"]]
        cp=comp[(comp.dataset.eq(study))&(comp.composition_status.eq("evaluable_author_labels"))]
        wide=cp.pivot(index="subject_id",columns="lineage",values="proportion").fillna(0)
        for term,d in sc.groupby("program",observed=True):
            y=d.set_index("subject_id").standardized_score
            for lineage in wide.columns:
                ids=y.index.intersection(wide.index); rho,p=spearmanr(y.loc[ids],wide.loc[ids,lineage])
                rows.append({"analysis_type":"spearman_program_vs_lineage_proportion","dataset":study,
                             "program":term,"contrast_id":"not_applicable","covariate":lineage,
                             "n_biological_samples":len(ids),"estimate":rho,"nominal_p":p,
                             "status":"evaluable_exploratory","reason":"descriptive association; no causal interpretation"})
        pc=pca[pca.dataset.eq(study)][["subject_id","PC1"]]
        for term,d in sc.groupby("program",observed=True):
            z=d.merge(pc,on="subject_id",validate="one_to_one")
            x=pd.DataFrame({"const":1.0,"PC1":z.PC1,"dpi_1":z.group.eq("dpi_1").astype(float),
                            "dpi_7":z.group.eq("dpi_7").astype(float)})
            rank=int(np.linalg.matrix_rank(x)); df=len(x)-rank; cond=float(np.linalg.cond(x))
            try: vifs={c:float(variance_inflation_factor(x.to_numpy(float),i)) for i,c in enumerate(x.columns) if c!="const"}
            except Exception: vifs={c:np.inf for c in x.columns if c!="const"}
            maxv=max(vifs.values()); ok=rank==x.shape[1] and df>=3 and maxv<=10 and cond<1000
            model=sm.OLS(z.standardized_score.to_numpy(float),x.to_numpy(float)).fit()
            for contrast,v in CONTRASTS.items():
                est=float(v@model.params); se=float(math.sqrt(v@model.cov_params()@v)); crit=float(student_t.ppf(.975,model.df_resid))
                orig=original[(original.dataset==study)&(original.method=="centered_rank_mean")&
                              (original.program==term)&(original.contrast_id==contrast)]
                oe=float(orig.estimate.iloc[0]) if len(orig) else np.nan
                rows.append({"analysis_type":"composition_PC1_adjusted_group_model","dataset":study,
                    "program":term,"contrast_id":contrast,"covariate":"CLR_PC1","n_biological_samples":len(z),
                    "estimate":est,"se":se,"df_residual":model.df_resid,"ci_low":est-crit*se,"ci_high":est+crit*se,
                    "nominal_p":float(2*student_t.sf(abs(est/se),model.df_resid)) if se>0 else np.nan,
                    "original_centered_rank_estimate":oe,"direction_matches_original":bool(np.sign(est)==np.sign(oe)),
                    "design_rank":rank,"design_columns":x.shape[1],"condition_number":cond,"max_vif":maxv,
                    "status":"evaluable_exploratory" if ok else "not_evaluable",
                    "reason":"one-PC frozen adjustment" if ok else "frozen rank/df/VIF/condition-number gate failed"})
    out=pd.DataFrame(rows); out["nominal_fdr_within_analysis_type"] = np.nan
    for _,idx in out[out.nominal_p.notna()].groupby(["analysis_type","dataset"],observed=True).groups.items():
        out.loc[idx,"nominal_fdr_within_analysis_type"]=bh(out.loc[idx,"nominal_p"])
    return out


def ontology_summary(recomp:pd.DataFrame) -> pd.DataFrame:
    rows=[]
    for (study,term,contrast),d in recomp.groupby(["dataset","program","contrast_id"],observed=True):
        ev=d[d.status.eq("evaluable_descriptive")]
        if len(ev)!=3: cls="not_evaluable"
        elif ev.direction_matches_original.all() and len(set(np.sign(ev.estimate)))==1: cls="ontology_robust_descriptive"
        else: cls="ontology_sensitive"
        rows.append({"dataset":study,"program":term,"contrast_id":contrast,"classification":cls,
                     "n_modes_evaluable":len(ev),"modes":";".join(ev.ontology_mode.astype(str)),
                     "estimates":";".join(f"{r.ontology_mode}:{r.estimate:.8g}" for r in ev.itertuples()),
                     "all_modes_match_original":bool(len(ev)==3 and ev.direction_matches_original.all()),
                     "status":"evaluable_descriptive" if len(ev)==3 else "not_evaluable"})
    for term in FOCAL:
        for contrast in CONTRASTS:
            rows.append({"dataset":"GSE304399","program":term,"contrast_id":contrast,
                         "classification":"not_evaluable","n_modes_evaluable":0,"modes":"",
                         "estimates":"","all_modes_match_original":False,"status":"not_evaluable",
                         "reason":"no barcode-level author cell labels"})
    return pd.DataFrame(rows)


def main()->None:
    ap=argparse.ArgumentParser(); ap.add_argument("--root",type=Path,default=Path(__file__).resolve().parents[1])
    ap.add_argument("--out-dir",type=Path,default=Path("reports/phase_final_calibration_and_v5_2026_07"))
    args=ap.parse_args(); root=args.root.resolve(); out=args.out_dir if args.out_dir.is_absolute() else root/args.out_dir
    out.mkdir(parents=True,exist_ok=True)
    gmt_path=root/"references/msigdb_mh.all.v2026.1.Mm.symbols.gmt"; sets=read_gmt(gmt_path)
    score_path=root/"reports/phase_reproducibility_calibration_2026_07/sample_level_hallmark_scores.tsv"
    effect_path=root/"reports/phase_reproducibility_calibration_2026_07/pathway_effects_by_study.tsv"
    scores=pd.read_csv(score_path,sep="\t"); original=pd.read_csv(effect_path,sep="\t")
    comp_frames=[]; recomp_frames=[]; inputs=[gmt_path,score_path,effect_path]
    for study in LABELED:
        matrix,genes,meta,paths=load_harmonized(root,study); inputs.extend(paths)
        comp_frames.append(composition_long(meta,study)); recomp_frames.append(recomposition(root,study,matrix,genes,meta,sets,original))
    unresolved=pd.read_csv(root/"data_processed/whole_lesion/GSE304399_pseudobulk_coldata.tsv",sep="\t")
    unresolved=unresolved[unresolved.group.isin(GROUPS)]
    comp_frames.append(pd.DataFrame({"dataset":"GSE304399","sample":unresolved["sample"],
        "subject_id":unresolved["subject_id"],"group":unresolved["group"],"lineage":"not_evaluable_unresolved",
        "n_cells":unresolved["n_cells"],"sample_total_cells":unresolved["n_cells"],"proportion":np.nan,
        "composition_status":"not_evaluable_no_barcode_level_labels"}))
    comp=pd.concat(comp_frames,ignore_index=True); comp.to_csv(out/"biological_sample_composition.tsv",sep="\t",index=False,na_rep="NA")
    pca=clr_pca(comp)
    pca=pd.concat([pca,pd.DataFrame({"dataset":["GSE304399"],"status":["not_evaluable"],
        "reason":["no barcode-level author cell labels"]})],ignore_index=True,sort=False)
    pca.to_csv(out/"composition_pca.tsv",sep="\t",index=False,na_rep="NA")
    assoc=associations(comp,pca,scores,original)
    missing=[]
    for term in FOCAL:
        for contrast in CONTRASTS:
            missing.append({"analysis_type":"composition_PC1_adjusted_group_model","dataset":"GSE304399",
                            "program":term,"contrast_id":contrast,"covariate":"not_evaluable",
                            "status":"not_evaluable","reason":"no barcode-level author cell labels"})
    assoc=pd.concat([assoc,pd.DataFrame(missing)],ignore_index=True,sort=False)
    assoc.to_csv(out/"composition_program_association.tsv",sep="\t",index=False,na_rep="NA")
    recomp=pd.concat(recomp_frames,ignore_index=True); 
    missing_re=[]
    for mode in ("strict","broad_family","exclude_ambiguous"):
        for term in FOCAL:
            for contrast in CONTRASTS:
                missing_re.append({"dataset":"GSE304399","ontology_mode":mode,"program":term,
                                   "contrast_id":contrast,"status":"not_evaluable",
                                   "reason":"no barcode-level author cell labels"})
    recomp=pd.concat([recomp,pd.DataFrame(missing_re)],ignore_index=True,sort=False)
    recomp.to_csv(out/"recomposed_pseudobulk_results.tsv",sep="\t",index=False,na_rep="NA")
    onto=ontology_summary(recomp); onto.to_csv(out/"ontology_sensitivity_results.tsv",sep="\t",index=False,na_rep="NA")
    outputs=[out/x for x in ["biological_sample_composition.tsv","composition_pca.tsv",
        "composition_program_association.tsv","recomposed_pseudobulk_results.tsv","ontology_sensitivity_results.tsv"]]
    prov={"analysis":"composition and ontology sensitivity","created_at":datetime.now().astimezone().isoformat(),
          "adjustment":"one CLR composition PC only; frozen rank>=columns, residual df>=3, VIF<=10, condition number<1000",
          "ontology_modes":["strict","broad_family","exclude_ambiguous"],
          "recomposition":"fixed average lineage weights among lineages present in every eligible sample; descriptive",
          "boundaries":["GSE304399 is not evaluable","no missing lineage is imputed","nominal P is not replication",
                        "composition association is not causality","enrichment is not activation"],
          "inputs":[{"path":str(p.relative_to(root)),"bytes":p.stat().st_size,"sha256":sha256(p)} for p in sorted(set(inputs))],
          "outputs":[{"path":str(p.relative_to(root)),"bytes":p.stat().st_size,"sha256":sha256(p)} for p in outputs]}
    (out/"composition_ontology_provenance.json").write_text(json.dumps(prov,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":"ok","composition_rows":len(comp),"adjustment_evaluable":int(assoc.status.eq("evaluable_exploratory").sum()),
                      "ontology_robust":int(onto.classification.eq("ontology_robust_descriptive").sum())}))


if __name__=="__main__": main()
