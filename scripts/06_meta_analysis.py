#!/usr/bin/env python3
"""Comparability-aware random-effects meta-analysis of prespecified effects.

Inputs are study-level direct contrasts from ``05_pseudobulk_DE.R``.  The mapping table
decides which contrasts are scientifically comparable and which datasets are independent.
Same-family modalities and lesion-remote tissue can therefore never inflate the number of
replications by accident.

Estimator: REML tau² + modified Hartung-Knapp (q>=1) uncertainty.  With only two
independent studies, estimates are reported as ``limited_two_study_*`` rather than strong
replication.  No result from this module is called a direction reversal; that joint claim
is evaluated transparently by ``06c_reversal_statistic.py``.
"""
from __future__ import annotations

import argparse
import glob
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.stats import chi2, t


REQUIRED_EFFECT_COLUMNS = {
    "gene", "cell_state", "dataset", "contrast_id", "effect_type", "estimate", "se", "p", "fdr"
}


def bh_adjust(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    out = np.full(values.shape, np.nan)
    ok = np.isfinite(values)
    p = values[ok]
    if not len(p):
        return out
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    adjusted = np.empty_like(ranked)
    adjusted[order] = np.clip(ranked, 0, 1)
    out[ok] = adjusted
    return out


def harmonize_cell_states(values: pd.Series, mapping: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    ordered = mapping.sort_values("priority", kind="stable")
    known_lineages = set(ordered["harmonized_lineage"].astype(str))
    harmonized, matched_pattern = [], []
    for value in values.astype(str):
        if value in known_lineages:
            harmonized.append(value); matched_pattern.append("__already_harmonized__")
            continue
        matched = ordered[ordered["pattern"].apply(lambda p: bool(re.search(str(p), value, flags=re.I)))]
        if matched.empty:
            harmonized.append(pd.NA); matched_pattern.append(pd.NA)
        else:
            row = matched.iloc[0]
            harmonized.append(row["harmonized_lineage"]); matched_pattern.append(row["pattern"])
    return pd.Series(harmonized, index=values.index, dtype="string"), pd.Series(matched_pattern, index=values.index, dtype="string")


def reml_tau2(y: np.ndarray, v: np.ndarray) -> float:
    y, v = np.asarray(y, float), np.asarray(v, float)
    if len(y) < 2 or np.allclose(y, y[0]):
        return 0.0
    if len(y) == 2:
        # With an intercept-only two-study model, the REML likelihood is the
        # likelihood of the study difference: Var(y1-y2)=v1+v2+2*tau2.
        return max(0.0, float(((y[0] - y[1]) ** 2 - v[0] - v[1]) / 2.0))
    upper = max(float(np.var(y, ddof=1) * 10), float(np.max(v) * 100), 1.0)

    def objective(tau2: float) -> float:
        w = 1.0 / (v + tau2)
        mu = float(np.sum(w * y) / np.sum(w))
        return 0.5 * (np.sum(np.log(v + tau2)) + np.log(np.sum(w)) + np.sum(w * (y - mu) ** 2))

    result = minimize_scalar(objective, bounds=(0.0, upper), method="bounded", options={"xatol": 1e-12})
    if not result.success:
        return 0.0
    tau2 = float(result.x)
    return 0.0 if tau2 < 1e-10 else tau2


def meta_reml_hksj(y: np.ndarray, se: np.ndarray) -> dict:
    y, se = np.asarray(y, float), np.asarray(se, float)
    if len(y) != len(se) or len(y) == 0 or np.any(~np.isfinite(y)) or np.any(~np.isfinite(se)) or np.any(se <= 0):
        raise ValueError("meta inputs require finite estimates and positive SE")
    k = len(y)
    v = se ** 2
    if k == 1:
        return dict(k=1, estimate=float(y[0]), se_hksj=float(se[0]), ci_low=float("nan"), ci_high=float("nan"),
                    p=float("nan"), tau2=float("nan"), i2=float("nan"), q=float("nan"), q_p=float("nan"),
                    prediction_low=float("nan"), prediction_high=float("nan"), max_weight_fraction=1.0,
                    direction_concordance=1.0)
    tau2 = reml_tau2(y, v)
    w = 1.0 / (v + tau2)
    estimate = float(np.sum(w * y) / np.sum(w))
    q_re = float(np.sum(w * (y - estimate) ** 2))
    hk_scale = max(q_re / (k - 1), 1.0)  # modified HK avoids variance deflation
    se_hk = math.sqrt(hk_scale / np.sum(w))
    crit = float(t.ppf(0.975, df=k - 1))
    statistic = estimate / se_hk
    p_value = float(2 * t.sf(abs(statistic), df=k - 1))
    ci_low, ci_high = estimate - crit * se_hk, estimate + crit * se_hk
    w_fixed = 1.0 / v
    mu_fixed = float(np.sum(w_fixed * y) / np.sum(w_fixed))
    q_fixed = float(np.sum(w_fixed * (y - mu_fixed) ** 2))
    q_p = float(chi2.sf(q_fixed, df=k - 1))
    i2 = max(0.0, (q_fixed - (k - 1)) / q_fixed) * 100 if q_fixed > 0 else 0.0
    if k >= 3:
        pred_crit = float(t.ppf(0.975, df=k - 2))
        pred_half = pred_crit * math.sqrt(tau2 + se_hk ** 2)
        pred_low, pred_high = estimate - pred_half, estimate + pred_half
    else:
        pred_low = pred_high = float("nan")
    nonzero = y[np.abs(y) > 0]
    concordance = float(np.mean(np.sign(nonzero) == np.sign(estimate))) if len(nonzero) and estimate else float("nan")
    return dict(k=k, estimate=estimate, se_hksj=se_hk, ci_low=ci_low, ci_high=ci_high, p=p_value,
                tau2=tau2, i2=i2, q=q_fixed, q_p=q_p, prediction_low=pred_low,
                prediction_high=pred_high, max_weight_fraction=float(np.max(w) / np.sum(w)),
                direction_concordance=concordance)


def leave_one_out(y: np.ndarray, se: np.ndarray, labels: list[str]) -> pd.DataFrame:
    rows = []
    for i, label in enumerate(labels):
        keep = np.arange(len(y)) != i
        if keep.sum() == 1:
            result = dict(estimate=float(y[keep][0]), ci_low=float("nan"), ci_high=float("nan"), p=float("nan"))
        else:
            result = meta_reml_hksj(y[keep], se[keep])
        rows.append(dict(omitted=label, k_remaining=int(keep.sum()), estimate=result["estimate"],
                         ci_low=result["ci_low"], ci_high=result["ci_high"], p=result["p"]))
    return pd.DataFrame(rows)


def classify(row: pd.Series, study_signs: np.ndarray) -> str:
    if row["k"] < 2:
        return "underpowered_single_independent_study"
    mixed = len(set(np.sign(study_signs[np.abs(study_signs) > 0]))) > 1
    significant = row["meta_fdr"] <= 0.05 and (row["ci_low"] > 0 or row["ci_high"] < 0)
    if mixed and (row["i2"] >= 50 or row["q_p"] < 0.10):
        return "cross_study_discordance"
    if significant and row["max_weight_fraction"] > 0.80:
        return "dominant_study_sensitive"
    if significant and row["direction_concordance"] == 1.0:
        return "limited_two_study_change" if row["k"] == 2 else "reproducible_change"
    if row["i2"] >= 50:
        return "heterogeneous_underpowered"
    return "no_reproducible_change"


def load_effects(patterns: list[str]) -> pd.DataFrame:
    paths = []
    for pattern in patterns:
        paths.extend(glob.glob(pattern))
    paths = sorted(set(paths))
    if not paths:
        raise FileNotFoundError(f"no effect files matched {patterns}")
    frames = []
    for path in paths:
        frame = pd.read_csv(path, sep="\t")
        missing = REQUIRED_EFFECT_COLUMNS - set(frame.columns)
        if missing:
            raise ValueError(f"{path} missing columns {sorted(missing)}")
        frame["source_file"] = path
        frames.append(frame)
    effects = pd.concat(frames, ignore_index=True)
    if effects[["dataset", "gene", "cell_state", "contrast_id"]].duplicated().any():
        raise ValueError("duplicate dataset/gene/cell_state/contrast effect rows")
    return effects


def run(effects: pd.DataFrame, contrast_map: pd.DataFrame, cell_map: pd.DataFrame,
        min_studies: int = 2) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    merged = effects.merge(contrast_map, on=["dataset", "contrast_id"], how="left", validate="many_to_one",
                           indicator=True)
    unmapped_contrasts = merged[merged["_merge"].ne("both")].copy()
    merged = merged[merged["_merge"].eq("both")].drop(columns="_merge")
    merged["harmonized_lineage"], merged["cell_mapping_pattern"] = harmonize_cell_states(merged["cell_state"], cell_map)
    unmapped_cells = merged[merged["harmonized_lineage"].isna()].copy()
    merged = merged[merged["harmonized_lineage"].notna()].copy()
    primary = merged[merged["include_primary"].astype(str).str.lower().eq("yes")].copy()
    duplicate_cluster = primary.duplicated(
        ["gene", "harmonized_lineage", "meta_contrast_id", "independence_cluster"], keep=False)
    if duplicate_cluster.any():
        bad = primary.loc[duplicate_cluster, ["dataset", "meta_contrast_id", "independence_cluster"]].drop_duplicates()
        raise ValueError(f"multiple primary effects from one independence cluster:\n{bad.to_string(index=False)}")

    rows, group_indices = [], {}
    group_cols = ["gene", "harmonized_lineage", "meta_contrast_id", "compartment", "effect_type"]
    for key, group in primary.groupby(group_cols, sort=False, observed=True):
        y = group["estimate"].to_numpy(float)
        se = group["se"].to_numpy(float)
        result = meta_reml_hksj(y, se)
        row = dict(zip(group_cols, key))
        row.update(result)
        row.update(datasets=";".join(group["dataset"].astype(str)),
                   independence_clusters=";".join(group["independence_cluster"].astype(str)),
                   min_input_fdr=float(group["fdr"].min()), max_input_fdr=float(group["fdr"].max()),
                   n_input_fdr_005=int((group["fdr"] <= 0.05).sum()),
                   study_estimates=";".join(f"{x:.12g}" for x in y))
        rows.append(row)
        if len(group) >= min_studies:
            group_indices[key] = group.index.to_numpy()
    meta = pd.DataFrame(rows)
    if meta.empty:
        return meta, pd.DataFrame(), unmapped_contrasts, unmapped_cells
    meta["meta_fdr"] = np.nan
    family_cols = ["meta_contrast_id", "harmonized_lineage", "effect_type"]
    for _, index in meta.groupby(family_cols, observed=True).groups.items():
        meta.loc[index, "meta_fdr"] = bh_adjust(meta.loc[index, "p"].to_numpy())
    # Retrieve study signs for the final class without selecting on study P/FDR.
    meta["evidence_class"] = [
        classify(row, np.fromstring(str(row.study_estimates), sep=";")) for _, row in meta.iterrows()
    ]
    # Full leave-one-out rows are retained for k>=3, potentially interesting FDR<=0.10
    # results, and heterogeneity/influence flags. This avoids a huge uninformative k=2 file.
    loo_mask = (meta["k"] >= 3) | (meta["meta_fdr"] <= 0.10)
    loo_rows = []
    for _, row in meta[loo_mask].iterrows():
        key = tuple(row[column] for column in group_cols)
        indices = group_indices.get(key)
        if indices is None:
            continue
        group = primary.loc[indices]
        loo = leave_one_out(group["estimate"].to_numpy(float), group["se"].to_numpy(float),
                            group["dataset"].astype(str).tolist())
        for column in group_cols:
            loo[column] = row[column]
        loo_rows.append(loo)
    meta["leave_one_out_evaluated"] = loo_mask
    loo = pd.concat(loo_rows, ignore_index=True) if loo_rows else pd.DataFrame()
    return meta, loo, unmapped_contrasts, unmapped_cells


def selftest() -> None:
    y = np.array([1.0, 1.1, 0.9])
    result = meta_reml_hksj(y, np.array([0.2, 0.2, 0.2]))
    assert result["ci_low"] > 0 and result["direction_concordance"] == 1.0
    discordant = meta_reml_hksj(np.array([2.0, -2.0, 1.8]), np.array([0.1, 0.1, 0.1]))
    assert discordant["i2"] > 90
    assert np.isclose(reml_tau2(np.array([2.0, -2.0]), np.array([0.01, 0.01])), 7.99)
    assert np.allclose(bh_adjust(np.array([0.01, 0.04, 0.2])), [0.03, 0.06, 0.2])
    try:
        meta_reml_hksj(np.array([1.0, 2.0]), np.array([0.0, 1.0]))
    except ValueError:
        pass
    else:
        raise AssertionError("zero SE accepted")
    print("[selftest] PASS: REML/mHK positive signal, heterogeneity, BH and invalid-SE checks")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--effects", nargs="+", default=["results/effects_GSE*.tsv"])
    ap.add_argument("--contrast-map", default="tables/meta_contrast_map.tsv")
    ap.add_argument("--cell-map", default="tables/cell_state_harmonization.tsv")
    ap.add_argument("--out-dir", default="results/meta")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest(); return
    effects = load_effects(args.effects)
    contrast_map = pd.read_csv(args.contrast_map, sep="\t")
    cell_map = pd.read_csv(args.cell_map, sep="\t")
    meta, loo, unmapped_contrasts, unmapped_cells = run(effects, contrast_map, cell_map)
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    meta.to_csv(out / "meta_effects.tsv", sep="\t", index=False)
    loo.to_csv(out / "leave_one_out.tsv", sep="\t", index=False)
    unmapped_contrasts.to_csv(out / "unmapped_contrasts.tsv", sep="\t", index=False)
    unmapped_cells.to_csv(out / "unmapped_cell_states.tsv", sep="\t", index=False)
    print(f"Wrote {len(meta)} meta rows; {len(unmapped_contrasts)} unmapped contrasts; "
          f"{len(unmapped_cells)} unmapped cell-state rows")


if __name__ == "__main__":
    main()
