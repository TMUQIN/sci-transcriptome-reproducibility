#!/usr/bin/env python3
"""Frozen Hallmark stability and MYC/cell-cycle audit for the final calibration phase.

The script deliberately separates four objects that are often conflated: fixed gene-set
membership, original GSEA leading-edge membership, biological-sample score direction,
and gene-level direction conditional on an observed leading edge. It never labels a
nominal P value or a conditional bootstrap as pathway replication.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import math
import subprocess
import sys
from datetime import datetime
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import pearsonr, spearmanr, t as student_t
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor


STUDIES = ("GSE162610", "GSE234774", "GSE304399")
GROUPS = ("uninjured", "dpi_1", "dpi_7")
CONTRASTS = {
    "injury_d1_vs_uninjured": ("dpi_1", "uninjured"),
    "injury_d7_vs_uninjured": ("dpi_7", "uninjured"),
    "change_d7_minus_d1": ("dpi_7", "dpi_1"),
}
EXPECTED = {"injury_d1_vs_uninjured": 1, "injury_d7_vs_uninjured": 1,
            "change_d7_minus_d1": -1}
FOCAL = (
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_HYPOXIA",
    "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_MTORC1_SIGNALING",
)
CELL_CYCLE = ("HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT")
MYC = ("HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_gmt(path: Path) -> dict[str, set[str]]:
    out: dict[str, set[str]] = {}
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                out[fields[0]] = set(fields[2:])
    if len(out) != 50:
        raise ValueError(f"expected 50 Hallmarks, observed {len(out)}")
    return out


def bh(values: pd.Series) -> pd.Series:
    x = values.to_numpy(float)
    order = np.argsort(x)
    q = x[order] * len(x) / np.arange(1, len(x) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    ans = np.empty_like(q); ans[order] = np.clip(q, 0, 1)
    return pd.Series(ans, index=values.index)


def membership_redundancy(sets: dict[str, set[str]]) -> pd.DataFrame:
    rows = []
    for a, b in combinations(sorted(sets), 2):
        inter, union = sets[a] & sets[b], sets[a] | sets[b]
        rows.append({"term_1": a, "term_2": b, "n_1": len(sets[a]), "n_2": len(sets[b]),
                     "n_intersection": len(inter), "n_union": len(union),
                     "jaccard": len(inter) / len(union),
                     "overlap_coefficient": len(inter) / min(len(sets[a]), len(sets[b])),
                     "intersection_genes": ";".join(sorted(inter)),
                     "either_focal": a in FOCAL or b in FOCAL})
    return pd.DataFrame(rows)


def score_correlations(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (dataset, method), frame in scores.groupby(["dataset", "method"], observed=True):
        wide = frame.pivot(index="subject_id", columns="program", values="standardized_score")
        for a, b in combinations(sorted(wide.columns), 2):
            pair = wide[[a, b]].dropna()
            r, p = spearmanr(pair[a], pair[b]) if len(pair) >= 4 else (np.nan, np.nan)
            rp, pp = pearsonr(pair[a], pair[b]) if len(pair) >= 4 else (np.nan, np.nan)
            rows.append({"dataset": dataset, "method": method, "term_1": a, "term_2": b,
                         "n_samples": len(pair), "spearman_rho": r, "spearman_p": p,
                         "pearson_r": rp, "pearson_p": pp,
                         "either_focal": a in FOCAL or b in FOCAL})
    for (dataset, program), frame in scores.groupby(["dataset", "program"], observed=True):
        wide = frame.pivot(index="subject_id", columns="method", values="standardized_score").dropna()
        if {"centered_rank_mean", "magnitude_mean_z"} <= set(wide) and len(wide) >= 4:
            r, p = spearmanr(wide["centered_rank_mean"], wide["magnitude_mean_z"])
            rp, pp = pearsonr(wide["centered_rank_mean"], wide["magnitude_mean_z"])
            rows.append({"dataset": dataset, "method": "cross_method", "term_1": program,
                         "term_2": program, "n_samples": len(wide), "spearman_rho": r,
                         "spearman_p": p, "pearson_r": rp, "pearson_p": pp,
                         "either_focal": program in FOCAL})
    result = pd.DataFrame(rows)
    result["spearman_fdr_within_dataset_method"] = np.nan
    for _, idx in result.groupby(["dataset", "method"], observed=True).groups.items():
        ok = result.loc[idx, "spearman_p"].notna()
        result.loc[np.asarray(idx)[ok], "spearman_fdr_within_dataset_method"] = bh(
            result.loc[np.asarray(idx)[ok], "spearman_p"])
    return result


def myc_overlap(sets: dict[str, set[str]]) -> pd.DataFrame:
    rows = []
    terms = (*MYC, *CELL_CYCLE)
    cycle_union = sets[CELL_CYCLE[0]] | sets[CELL_CYCLE[1]]
    for a, b in combinations(terms, 2):
        inter, union = sets[a] & sets[b], sets[a] | sets[b]
        rows.append({"term_1": a, "term_2": b, "n_1": len(sets[a]), "n_2": len(sets[b]),
                     "n_intersection": len(inter), "jaccard": len(inter) / len(union),
                     "intersection_genes": ";".join(sorted(inter))})
    for term in MYC:
        unique = sets[term] - cycle_union
        rows.append({"term_1": term, "term_2": "E2F_OR_G2M_UNION", "n_1": len(sets[term]),
                     "n_2": len(cycle_union), "n_intersection": len(sets[term] & cycle_union),
                     "jaccard": len(sets[term] & cycle_union) / len(sets[term] | cycle_union),
                     "intersection_genes": ";".join(sorted(sets[term] & cycle_union)),
                     "n_myc_specific": len(unique), "myc_specific_genes": ";".join(sorted(unique))})
    return pd.DataFrame(rows)


def load_score_module(root: Path):
    path = root / "scripts" / "22_sample_level_pathway_scores.py"
    spec = importlib.util.spec_from_file_location("pathway_scores", path)
    mod = importlib.util.module_from_spec(spec); assert spec and spec.loader
    spec.loader.exec_module(mod)
    return mod


def build_variant_sets(sets: dict[str, set[str]]) -> tuple[dict[str, set[str]], pd.DataFrame]:
    variants: dict[str, set[str]] = {}
    audit = []
    cycle = sets[CELL_CYCLE[0]] | sets[CELL_CYCLE[1]]
    for term in FOCAL:
        other = set().union(*(sets[x] for x in FOCAL if x != term))
        name = term + "__EXCL_OTHER_FOCAL"
        variants[name] = sets[term] - other
        audit.append({"variant": name, "source_term": term, "exclusion": "other_three_focal_Hallmarks",
                      "n_source": len(sets[term]), "n_removed": len(sets[term] & other),
                      "n_retained": len(variants[name])})
    for term in MYC:
        name = term + "__EXCL_E2F_G2M"
        variants[name] = sets[term] - cycle
        audit.append({"variant": name, "source_term": term, "exclusion": "E2F_or_G2M_members",
                      "n_source": len(sets[term]), "n_removed": len(sets[term] & cycle),
                      "n_retained": len(variants[name])})
    return variants, pd.DataFrame(audit)


def score_variants(root: Path, variants: dict[str, set[str]], out_dir: Path) -> tuple[pd.DataFrame, list[Path]]:
    mod = load_score_module(root)
    frames, inputs = [], []
    for study in STUDIES:
        counts, genes, coldata, paths = mod.load_dataset(root, study)
        frame = mod.score_dataset(counts, genes, coldata, variants, study)
        frames.append(frame); inputs.extend(paths)
    scores = pd.concat(frames, ignore_index=True)
    score_path = out_dir / "deoverlapped_hallmark_scores.tsv"
    scores.to_csv(score_path, sep="\t", index=False, na_rep="NA", float_format="%.17g")
    effect_path = out_dir / "deoverlapped_hallmark_results.tsv"
    proc = subprocess.run(["Rscript", str(root / "scripts" / "22_pathway_score_DE.R"),
                           str(score_path), str(effect_path)], cwd=root, text=True,
                          capture_output=True, timeout=300)
    (out_dir / "hallmark_myc_r_execution.log").write_text(
        f"returncode={proc.returncode}\nstdout:\n{proc.stdout}\nstderr:\n{proc.stderr}\n", encoding="utf-8")
    if proc.returncode:
        raise RuntimeError(f"de-overlap limma failed: {proc.stderr[-1000:]}")
    return pd.read_csv(effect_path, sep="\t"), [*inputs, score_path, effect_path]


def residual_models(scores: pd.DataFrame) -> pd.DataFrame:
    rows = []
    needed = (*MYC, *CELL_CYCLE)
    use = scores[scores.program.isin(needed)].copy()
    for (dataset, method), frame in use.groupby(["dataset", "method"], observed=True):
        wide = frame.pivot(index=["subject_id", "group"], columns="program",
                           values="standardized_score").reset_index()
        if not set(needed) <= set(wide):
            continue
        wide = wide[wide.group.isin(GROUPS)].dropna(subset=list(needed)).copy()
        d1 = (wide.group == "dpi_1").astype(float).to_numpy()
        d7 = (wide.group == "dpi_7").astype(float).to_numpy()
        x = pd.DataFrame({"const": 1.0, "E2F": wide[CELL_CYCLE[0]].to_numpy(float),
                          "G2M": wide[CELL_CYCLE[1]].to_numpy(float), "dpi_1": d1, "dpi_7": d7})
        rank = int(np.linalg.matrix_rank(x.to_numpy(float)))
        df = len(x) - rank
        cond = float(np.linalg.cond(x.to_numpy(float)))
        try:
            vifs = {col: float(variance_inflation_factor(x.to_numpy(float), i))
                    for i, col in enumerate(x.columns) if col != "const"}
        except Exception:
            vifs = {col: np.inf for col in x.columns if col != "const"}
        max_vif = max(vifs.values())
        evaluable = rank == x.shape[1] and df >= 3 and np.isfinite(max_vif) and max_vif <= 10 and cond < 1000
        for myc in MYC:
            model = sm.OLS(wide[myc].to_numpy(float), x.to_numpy(float)).fit()
            for contrast, vector in {
                "injury_d1_vs_uninjured": np.array([0, 0, 0, 1, 0.]),
                "injury_d7_vs_uninjured": np.array([0, 0, 0, 0, 1.]),
                "change_d7_minus_d1": np.array([0, 0, 0, -1, 1.]),
            }.items():
                est = float(vector @ model.params)
                se = float(math.sqrt(vector @ model.cov_params() @ vector))
                stat = est / se if se > 0 else np.nan
                p = float(2 * student_t.sf(abs(stat), model.df_resid)) if np.isfinite(stat) else np.nan
                crit = float(student_t.ppf(.975, model.df_resid))
                rows.append({"dataset": dataset, "method": method, "myc_program": myc,
                             "contrast_id": contrast, "estimate_adjusted": est, "se": se,
                             "df_residual": float(model.df_resid), "ci_low": est-crit*se,
                             "ci_high": est+crit*se, "nominal_p": p,
                             "direction": int(np.sign(est)), "expected_direction": EXPECTED[contrast],
                             "expected_direction_retained": bool(np.sign(est) == EXPECTED[contrast]),
                             "design_rank": rank, "design_columns": x.shape[1],
                             "condition_number": cond, "max_vif": max_vif,
                             "vif_E2F": vifs.get("E2F"), "vif_G2M": vifs.get("G2M"),
                             "vif_dpi_1": vifs.get("dpi_1"), "vif_dpi_7": vifs.get("dpi_7"),
                             "evaluable": evaluable,
                             "status": "evaluable_exploratory" if evaluable else "not_evaluable",
                             "reason": "" if evaluable else "frozen rank/df/VIF/condition-number gate failed"})
    result = pd.DataFrame(rows)
    result["nominal_p_fdr_within_method_contrast"] = np.nan
    for _, idx in result.groupby(["method", "contrast_id"], observed=True).groups.items():
        result.loc[idx, "nominal_p_fdr_within_method_contrast"] = bh(result.loc[idx, "nominal_p"])
    return result


def leading_lists(gsea: pd.DataFrame) -> dict[tuple[str, str, str], list[str]]:
    cache = {}
    for row in gsea.itertuples():
        genes = [] if pd.isna(row.leading_edge_genes) else [x for x in str(row.leading_edge_genes).split(";") if x]
        cache[(row.term, row.contrast_id, row.dataset)] = genes
    return cache


def leading_overlap(gsea: pd.DataFrame, sets: dict[str, set[str]], rng: np.random.Generator,
                    n_null: int) -> tuple[pd.DataFrame, pd.DataFrame]:
    cache = leading_lists(gsea)
    rows, loso = [], []
    for term in FOCAL:
        universe = np.array(sorted(sets[term]), dtype=object)
        for contrast in CONTRASTS:
            study_lists = {s: cache[(term, contrast, s)] for s in STUDIES}
            for a, b in combinations(STUDIES, 2):
                la, lb = study_lists[a], study_lists[b]
                sa, sb = set(la), set(lb); union = sa | sb
                jac = len(sa & sb) / len(union) if union else np.nan
                weights_a = {g: 1 / math.log2(i + 2) for i, g in enumerate(la)}
                weights_b = {g: 1 / math.log2(i + 2) for i, g in enumerate(lb)}
                wj = sum(min(weights_a.get(g, 0), weights_b.get(g, 0)) for g in union) / sum(
                    max(weights_a.get(g, 0), weights_b.get(g, 0)) for g in union) if union else np.nan
                null = np.empty(n_null)
                for i in range(n_null):
                    xa = set(rng.choice(universe, len(sa), replace=False))
                    xb = set(rng.choice(universe, len(sb), replace=False))
                    null[i] = len(xa & xb) / len(xa | xb)
                rows.append({"term": term, "contrast_id": contrast, "dataset_1": a, "dataset_2": b,
                             "n_1": len(sa), "n_2": len(sb), "jaccard": jac,
                             "weighted_jaccard_rank_weighted": wj,
                             "empirical_p_random_subsets_within_fixed_term": (1 + int((null >= jac).sum()))/(n_null+1),
                             "null_B": n_null})
            for heldout in STUDIES:
                train = [s for s in STUDIES if s != heldout]
                train_core = set(study_lists[train[0]]) & set(study_lists[train[1]])
                held = set(study_lists[heldout])
                loso.append({"term": term, "contrast_id": contrast, "heldout_dataset": heldout,
                             "train_datasets": ";".join(train), "n_train_pair_core": len(train_core),
                             "n_heldout": len(held), "n_core_retained_in_heldout": len(train_core & held),
                             "train_core_retention": len(train_core & held)/len(train_core) if train_core else np.nan,
                             "heldout_coverage_by_train_core": len(train_core & held)/len(held) if held else np.nan,
                             "status": "evaluable" if train_core and held else "not_evaluable"})
    pair = pd.DataFrame(rows)
    pair["empirical_fdr_within_contrast"] = np.nan
    for _, idx in pair.groupby("contrast_id", observed=True).groups.items():
        pair.loc[idx, "empirical_fdr_within_contrast"] = bh(pair.loc[idx, "empirical_p_random_subsets_within_fixed_term"])
    return pair, pd.DataFrame(loso)


def load_logcpm(root: Path, study: str, genes_needed: set[str]) -> tuple[np.ndarray, pd.DataFrame, list[str], list[Path]]:
    base = root / "data_processed" / "whole_lesion"
    paths = [base/f"{study}_pseudobulk_counts.mtx.gz", base/f"{study}_pseudobulk_genes.tsv",
             base/f"{study}_pseudobulk_coldata.tsv"]
    with gzip.open(paths[0], "rb") as handle:
        counts = mmread(handle)
    counts = counts.toarray() if hasattr(counts, "toarray") else np.asarray(counts)
    genes = pd.read_csv(paths[1], sep="\t").gene.astype(str)
    meta = pd.read_csv(paths[2], sep="\t")
    keep_sample = meta.group.isin(GROUPS).to_numpy(); meta = meta.loc[keep_sample].reset_index(drop=True)
    counts = counts[keep_sample]
    idx = [i for i, g in enumerate(genes) if g in genes_needed]
    selected = [str(genes.iloc[i]) for i in idx]
    lib = counts.sum(axis=1)
    logcpm = np.log2(counts[:, idx] / lib[:, None] * 1e6 + .5)
    return logcpm, meta, selected, paths


def perturbation_audits(root: Path, gsea: pd.DataFrame, sets: dict[str, set[str]], scores: pd.DataFrame,
                        rng: np.random.Generator, n_boot: int) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, list[Path]]:
    cache = leading_lists(gsea)
    b_rows, j_rows, ps_rows, inputs = [], [], [], []
    union = set().union(*(sets[x] for x in FOCAL))
    for study in STUDIES:
        logcpm, meta, genes, paths = load_logcpm(root, study, union); inputs.extend(paths)
        gi = {g: i for i, g in enumerate(genes)}
        for contrast, (num, den) in CONTRASTS.items():
            ni = np.where(meta.group.to_numpy() == num)[0]; di = np.where(meta.group.to_numpy() == den)[0]
            bn = rng.choice(ni, (n_boot, len(ni)), replace=True)
            bd = rng.choice(di, (n_boot, len(di)), replace=True)
            boot = logcpm[bn].mean(axis=1) - logcpm[bd].mean(axis=1)
            observed = logcpm[ni].mean(axis=0) - logcpm[di].mean(axis=0)
            for term in FOCAL:
                for gene in cache[(term, contrast, study)]:
                    if gene not in gi: continue
                    k = gi[gene]
                    b_rows.append({"term": term, "gene": gene, "dataset": study,
                                   "contrast_id": contrast, "original_leading_edge_member": True,
                                   "selection_recomputed": False,
                                   "metric": "direction_retention_conditional_on_original_leading_edge",
                                   "observed_logcpm_difference": observed[k],
                                   "bootstrap_original_direction_frequency": float((np.sign(boot[:, k]) == np.sign(observed[k])).mean()),
                                   "bootstrap_expected_program_direction_frequency": float((np.sign(boot[:, k]) == EXPECTED[contrast]).mean()),
                                   "bootstrap_B": n_boot})
                    signs = []
                    for drop in np.r_[ni, di]:
                        n2 = ni[ni != drop]; d2 = di[di != drop]
                        if len(n2) == 0 or len(d2) == 0: continue
                        signs.append(np.sign(logcpm[n2, k].mean() - logcpm[d2, k].mean()))
                    j_rows.append({"term": term, "gene": gene, "dataset": study,
                                   "contrast_id": contrast, "selection_recomputed": False,
                                   "metric": "LOSO_direction_retention_conditional_on_original_leading_edge",
                                   "n_leave_one_sample_out": len(signs),
                                   "jackknife_original_direction_frequency": float(np.mean(np.asarray(signs) == np.sign(observed[k]))),
                                   "jackknife_expected_program_direction_frequency": float(np.mean(np.asarray(signs) == EXPECTED[contrast]))})
    for (dataset, method, program), frame in scores[scores.program.isin(FOCAL)].groupby(
            ["dataset", "method", "program"], observed=True):
        for contrast, (num, den) in CONTRASTS.items():
            x = frame.loc[frame.group.eq(num), "standardized_score"].to_numpy(float)
            y = frame.loc[frame.group.eq(den), "standardized_score"].to_numpy(float)
            bs = rng.choice(x, (n_boot, len(x)), replace=True).mean(axis=1) - rng.choice(
                y, (n_boot, len(y)), replace=True).mean(axis=1)
            obs = float(x.mean() - y.mean())
            ps_rows.append({"dataset": dataset, "method": method, "program": program,
                            "contrast_id": contrast, "observed_mean_difference": obs,
                            "bootstrap_median": float(np.median(bs)), "bootstrap_ci_low": float(np.quantile(bs, .025)),
                            "bootstrap_ci_high": float(np.quantile(bs, .975)),
                            "original_direction_frequency": float((np.sign(bs) == np.sign(obs)).mean()),
                            "expected_direction_frequency": float((np.sign(bs) == EXPECTED[contrast]).mean()),
                            "bootstrap_B": n_boot,
                            "interpretation": "sample-score direction stability; not GSEA leading-edge reselection"})
    return pd.DataFrame(b_rows), pd.DataFrame(j_rows), pd.DataFrame(ps_rows), inputs


def rank_overlap(root: Path, sets: dict[str, set[str]]) -> pd.DataFrame:
    frames = []
    for study in STUDIES:
        path = root / "results" / "whole_lesion" / f"effects_{study}.tsv"
        x = pd.read_csv(path, sep="\t", usecols=["gene", "contrast_id", "estimate", "se"])
        x["dataset"] = study; x["rank_statistic"] = x.estimate / x.se
        frames.append(x); 
    allx = pd.concat(frames, ignore_index=True)
    rows = []
    for term in FOCAL:
        for contrast in CONTRASTS:
            sub = allx[allx.gene.isin(sets[term]) & allx.contrast_id.eq(contrast)]
            wide = sub.pivot(index="gene", columns="dataset", values="rank_statistic")
            for a, b in combinations(STUDIES, 2):
                pair = wide[[a, b]].dropna()
                rho, p = spearmanr(pair[a], pair[b]) if len(pair) >= 10 else (np.nan, np.nan)
                q = max(int(math.ceil(.2 * len(pair))), 1)
                ta = set(pair[a].abs().nlargest(q).index); tb = set(pair[b].abs().nlargest(q).index)
                top = ta | tb
                rows.append({"term": term, "contrast_id": contrast, "dataset_1": a,
                             "dataset_2": b, "n_genes": len(pair), "rank_stat_spearman": rho,
                             "rank_stat_spearman_p": p, "top20_abs_union_n": len(top),
                             "top20_abs_jaccard": len(ta & tb)/len(top),
                             "top20_union_sign_agreement": float((np.sign(pair.loc[list(top), a]) == np.sign(pair.loc[list(top), b])).mean())})
    return pd.DataFrame(rows)


def stable_core(root: Path, bootstrap: pd.DataFrame, jackknife: pd.DataFrame, sets: dict[str, set[str]]) -> pd.DataFrame:
    core = pd.read_csv(root / "results/whole_lesion_programs/leading_edge_audit/temporal_gene_core.tsv", sep="\t")
    b = bootstrap.groupby(["term", "gene"], observed=True).agg(
        bootstrap_mean_expected=("bootstrap_expected_program_direction_frequency", "mean"),
        bootstrap_min_expected=("bootstrap_expected_program_direction_frequency", "min"),
        n_bootstrap_instances=("contrast_id", "size")).reset_index()
    j = jackknife.groupby(["term", "gene"], observed=True).agg(
        jackknife_mean_expected=("jackknife_expected_program_direction_frequency", "mean"),
        jackknife_min_expected=("jackknife_expected_program_direction_frequency", "min"),
        n_jackknife_instances=("contrast_id", "size")).reset_index()
    out = core.merge(b, on=["term", "gene"], how="left").merge(j, on=["term", "gene"], how="left")
    membership_count = {g: sum(g in s for s in sets.values()) for g in set().union(*sets.values())}
    out["n_hallmark_memberships"] = out.gene.map(membership_count)
    out["stability_class"] = "unstable"
    out.loc[out.partial_shared_temporal_driver & out.bootstrap_mean_expected.ge(.70) &
            out.jackknife_mean_expected.ge(.80), "stability_class"] = "moderately_stable"
    out.loc[out.strong_shared_temporal_driver & out.bootstrap_min_expected.ge(.70) &
            out.jackknife_min_expected.eq(1), "stability_class"] = "stable_core"
    out.loc[out.bootstrap_mean_expected.isna(), "stability_class"] = "not_evaluable"
    return out


def myc_celltypes(root: Path) -> pd.DataFrame:
    path = root / "results/GSE304361_support/programs/GSE304361_frozen_program_results.tsv"
    x = pd.read_csv(path, sep="\t")
    x = x[x.term.isin(MYC)].copy()
    x["interpretation_class"] = np.where(x.gsea_status.ne("tested"), "not_evaluable",
        np.where(x.fdr <= .05, "FDR_supported_enrichment", np.where(x.nominal_p <= .05,
                 "nominal_enrichment_only", "directional_or_no_support")))
    extra = []
    for stratum in ("cycling_cells", "non_cycling_cells"):
        extra.append({"dataset": "GSE304361", "analysis_scope": "cell_cycle_stratum",
                      "cell_state": stratum, "contrast_id": "wt_d7_vs_wt_uninjured",
                      "effect_type": "endpoint", "term": "HALLMARK_MYC_TARGETS_V1",
                      "gsea_status": "not_evaluable", "interpretation_class": "not_evaluable",
                      "reason": "deposited analysis has no frozen cycling/non-cycling biological-sample stratum"})
    return pd.concat([x, pd.DataFrame(extra)], ignore_index=True, sort=False)


def validation_candidates(core: pd.DataFrame, sets: dict[str, set[str]], root: Path) -> pd.DataFrame:
    cycle = sets[CELL_CYCLE[0]] | sets[CELL_CYCLE[1]]
    myc = core[core.term.eq(MYC[0])].copy()
    myc["not_E2F_G2M_member"] = ~myc.gene.isin(cycle)
    effects_path = root / "results/GSE304361_support/cell_states/effects_GSE304361.tsv"
    effects = pd.read_csv(effects_path, sep="\t", usecols=["gene", "cell_state", "contrast_id", "estimate", "fdr"])
    wt = effects[effects.contrast_id.eq("wt_d7_vs_wt_uninjured")]
    cell_rows = []
    for row in myc.itertuples():
        g = wt[wt.gene.eq(row.gene)]
        positive = sorted(g.loc[g.estimate.gt(0), "cell_state"].astype(str).unique())
        sig = sorted(g.loc[g.estimate.gt(0) & g.fdr.le(.05), "cell_state"].astype(str).unique())
        cell_rows.append({"gene": row.gene, "n_GSE304361_celltypes_tested": g.cell_state.nunique(),
                          "n_positive_celltypes": len(positive), "positive_celltypes": ";".join(positive),
                          "positive_FDR_celltypes": ";".join(sig),
                          "rna_detectability_status": "tested_in_celltype_pseudobulk" if len(g) else "not_evaluable"})
    myc = myc.merge(pd.DataFrame(cell_rows), on="gene", how="left")
    myc["protein_detectability_status"] = "not_evaluated_no_proteomic_input"
    myc["frozen_candidate_eligibility"] = (myc.strong_shared_temporal_driver &
        myc.not_E2F_G2M_member & myc.bootstrap_min_expected.ge(.70) &
        myc.n_positive_celltypes.fillna(0).ge(1))
    myc["candidate_score"] = (4*myc.frozen_candidate_eligibility.astype(int) +
                              myc.bootstrap_min_expected.fillna(0) +
                              .25*myc.bootstrap_mean_expected.fillna(0) +
                              .05*myc.n_positive_celltypes.fillna(0))
    myc = myc.sort_values(["frozen_candidate_eligibility", "bootstrap_min_expected",
                           "n_positive_celltypes"], ascending=False)
    myc["candidate_rank"] = np.arange(1, len(myc)+1)
    myc["candidate_status"] = np.where(myc.frozen_candidate_eligibility & myc.candidate_rank.le(3),
        "RNA_prioritized_protein_feasibility_pending", "not_shortlisted")
    return myc


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument("--out-dir", type=Path, default=Path("reports/phase_final_calibration_and_v5_2026_07"))
    ap.add_argument("--seed", type=int, default=20260729)
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--overlap-null", type=int, default=10000)
    args = ap.parse_args(); root = args.root.resolve()
    out = args.out_dir if args.out_dir.is_absolute() else root / args.out_dir; out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    gmt = root / "references/msigdb_mh.all.v2026.1.Mm.symbols.gmt"
    score_path = root / "reports/phase_reproducibility_calibration_2026_07/sample_level_hallmark_scores.tsv"
    gsea_path = root / "results/whole_lesion_programs/hallmark_gsea_by_study.tsv"
    sets = read_gmt(gmt); scores = pd.read_csv(score_path, sep="\t"); gsea = pd.read_csv(gsea_path, sep="\t")
    membership_redundancy(sets).to_csv(out/"hallmark_redundancy.tsv", sep="\t", index=False)
    score_correlations(scores).to_csv(out/"hallmark_score_correlations.tsv", sep="\t", index=False, na_rep="NA")
    myc_overlap(sets).to_csv(out/"myc_e2f_g2m_overlap.tsv", sep="\t", index=False, na_rep="NA")
    variants, variant_audit = build_variant_sets(sets)
    deoverlap, variant_inputs = score_variants(root, variants, out)
    variant_audit.to_csv(out/"deoverlapped_hallmark_variant_audit.tsv", sep="\t", index=False)
    deoverlap[deoverlap.program.str.startswith("HALLMARK_MYC_TARGETS")].to_csv(
        out/"myc_deoverlap_results.tsv", sep="\t", index=False, na_rep="NA")
    residual = residual_models(scores); residual.to_csv(out/"myc_residual_models.tsv", sep="\t", index=False, na_rep="NA")
    pair, loso = leading_overlap(gsea, sets, rng, args.overlap_null)
    pair.to_csv(out/"leading_edge_pairwise_calibrated.tsv", sep="\t", index=False, na_rep="NA")
    loso.to_csv(out/"leading_edge_leave_one_study_out.tsv", sep="\t", index=False, na_rep="NA")
    boot, jack, sample_boot, perturb_inputs = perturbation_audits(root, gsea, sets, scores, rng, args.bootstrap)
    boot.to_csv(out/"leading_edge_bootstrap.tsv", sep="\t", index=False, na_rep="NA")
    jack.to_csv(out/"leading_edge_jackknife.tsv", sep="\t", index=False, na_rep="NA")
    sample_boot.to_csv(out/"sample_score_bootstrap.tsv", sep="\t", index=False, na_rep="NA")
    ranks = rank_overlap(root, sets); ranks.to_csv(out/"continuous_rank_overlap.tsv", sep="\t", index=False, na_rep="NA")
    core = stable_core(root, boot, jack, sets); core.to_csv(out/"stable_driver_core.tsv", sep="\t", index=False, na_rep="NA")
    cell = myc_celltypes(root); cell.to_csv(out/"myc_celltype_results.tsv", sep="\t", index=False, na_rep="NA")
    cand = validation_candidates(core, sets, root); cand.to_csv(out/"myc_validation_candidates.tsv", sep="\t", index=False, na_rep="NA")
    outputs = [out/name for name in ["hallmark_redundancy.tsv", "hallmark_score_correlations.tsv",
        "myc_e2f_g2m_overlap.tsv", "deoverlapped_hallmark_scores.tsv", "deoverlapped_hallmark_results.tsv",
        "deoverlapped_hallmark_variant_audit.tsv", "myc_deoverlap_results.tsv", "myc_residual_models.tsv",
        "leading_edge_pairwise_calibrated.tsv", "leading_edge_leave_one_study_out.tsv",
        "leading_edge_bootstrap.tsv", "leading_edge_jackknife.tsv", "sample_score_bootstrap.tsv",
        "continuous_rank_overlap.tsv", "stable_driver_core.tsv", "myc_celltype_results.tsv",
        "myc_validation_candidates.tsv", "hallmark_myc_r_execution.log"]]
    inputs = sorted(set([gmt, score_path, gsea_path, root/"scripts/22_sample_level_pathway_scores.py",
                         root/"scripts/22_pathway_score_DE.R", *perturb_inputs, *variant_inputs]))
    prov = {"analysis": "Hallmark stability and MYC/cell-cycle audit", "created_at": datetime.now().astimezone().isoformat(),
            "seed": args.seed, "bootstrap_B": args.bootstrap, "overlap_null_B": args.overlap_null,
            "evidence_boundaries": ["leading-edge bootstrap conditions on original membership and does not rerun GSEA selection",
                "sample-score bootstrap is not pathway replication", "MYC residual models are exploratory and must pass rank/df/VIF gates",
                "nominal P values are not replication", "transcript enrichment is not pathway activation"],
            "inputs": [{"path": str(p.relative_to(root)), "bytes": p.stat().st_size, "sha256": sha256(p)} for p in inputs if p.exists()],
            "outputs": [{"path": str(p.relative_to(root)), "bytes": p.stat().st_size, "sha256": sha256(p)} for p in outputs]}
    (out/"hallmark_myc_stability_provenance.json").write_text(json.dumps(prov, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"status":"ok", "stable_core": int((core.stability_class=="stable_core").sum()),
                      "residual_evaluable": int(residual.evaluable.sum()), "outputs": len(outputs)}))


if __name__ == "__main__":
    main()
