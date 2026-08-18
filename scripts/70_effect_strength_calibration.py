#!/usr/bin/env python3
"""Effect-strength calibration for the frozen three-study whole-lesion branch.

This is an exploratory post-lock sensitivity analysis. It does not alter the
frozen transfer metrics or meta-analysis.  It asks whether the endpoint versus
cross-sectional contrast pattern remains visible within comparable absolute
study-score strata and in a deterministic strength-overlap window.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image


STUDIES = ["GSE162610", "GSE234774", "GSE304399"]
CONTRASTS = [
    ("injury_d1_vs_uninjured", "day 1 vs uninjured", "endpoint"),
    ("injury_d7_vs_uninjured", "day 7 vs uninjured", "endpoint"),
    ("change_d7_minus_d1", "day 7 vs day 1", "between_timepoint"),
]


def _sign(x: pd.Series) -> pd.Series:
    return np.sign(pd.to_numeric(x, errors="coerce")).astype(float)


def _heldout_balanced_accuracy(wide: pd.DataFrame) -> float:
    values = []
    for study in STUDIES:
        others = [s for s in STUDIES if s != study]
        pred = np.sign(wide[others].mean(axis=1))
        truth = np.sign(wide[study])
        keep = pred.ne(0) & truth.ne(0) & pred.notna() & truth.notna()
        pred = pred[keep]
        truth = truth[keep]
        if truth.empty or truth.nunique() < 2:
            continue
        tp = int(((truth == 1) & (pred == 1)).sum())
        tn = int(((truth == -1) & (pred == -1)).sum())
        fp = int(((truth == -1) & (pred == 1)).sum())
        fn = int(((truth == 1) & (pred == -1)).sum())
        sens = tp / (tp + fn) if tp + fn else np.nan
        spec = tn / (tn + fp) if tn + fp else np.nan
        if np.isfinite(sens) and np.isfinite(spec):
            values.append((sens + spec) / 2)
    return float(np.mean(values)) if values else np.nan


def _metric(wide: pd.DataFrame) -> dict[str, float]:
    signs = wide.apply(_sign)
    nonzero = (signs != 0).all(axis=1)
    signs = signs.loc[nonzero]
    if signs.empty:
        return {"n_features": 0, "all_study_concordance": np.nan, "heldout_balanced_accuracy": np.nan}
    all_same = signs.nunique(axis=1).eq(1)
    return {
        "n_features": int(len(signs)),
        "all_study_concordance": float(all_same.mean()),
        "heldout_balanced_accuracy": _heldout_balanced_accuracy(wide.loc[signs.index]),
    }


def _load_scale(root: Path, scale: str) -> pd.DataFrame:
    rows = []
    if scale == "gene":
        for study in STUDIES:
            p = root / "results" / "whole_lesion" / f"effects_{study}.tsv"
            x = pd.read_csv(p, sep="\t")
            x = x.loc[x["effect_type"].isin(["endpoint", "temporal_delta"])].copy()
            x["score"] = x["estimate"] / x["se"].replace(0, np.nan)
            x["feature"] = x["gene"]
            x["study"] = study
            rows.append(x[["feature", "study", "contrast_id", "score"]])
    else:
        p = root / "results" / "whole_lesion_programs" / "hallmark_gsea_by_study.tsv"
        x = pd.read_csv(p, sep="\t")
        x["feature"] = x["term"]
        x["score"] = x["nes"]
        x["study"] = x["dataset"]
        rows.append(x[["feature", "study", "contrast_id", "score"]])
    long = pd.concat(rows, ignore_index=True)
    long = long.loc[long["study"].isin(STUDIES)].dropna(subset=["score"])
    # Keep only features observed in all three studies for each contrast.
    counts = long.groupby(["contrast_id", "feature"])["study"].nunique()
    keep = counts[counts == len(STUDIES)].index
    long = long.set_index(["contrast_id", "feature"]).loc[keep].reset_index()
    return long


def _analyse_scale(long: pd.DataFrame, scale: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    wide_rows = []
    for contrast_id, label, contrast_class in CONTRASTS:
        x = long.loc[long["contrast_id"] == contrast_id].pivot(index="feature", columns="study", values="score")
        x = x.reindex(columns=STUDIES).dropna()
        x["strength"] = x[STUDIES].abs().mean(axis=1)
        wide_rows.append(
            x.reset_index().assign(
                scale=scale,
                contrast_id=contrast_id,
                contrast_label=label,
                contrast_class=contrast_class,
            )
        )
    all_wide = pd.concat(wide_rows, ignore_index=True)
    # Common absolute-score cut points within scale; this makes bins comparable
    # across contrasts instead of defining a separate strength scale for each.
    q = np.quantile(all_wide["strength"].to_numpy(), [0, 0.25, 0.5, 0.75, 1])
    q = np.maximum.accumulate(q)
    labels = ["Q1 (lowest)", "Q2", "Q3", "Q4 (highest)"]
    all_wide["strength_bin"] = pd.cut(
        all_wide["strength"], bins=np.unique(q), labels=labels[: len(np.unique(q)) - 1], include_lowest=True,
    )
    metrics = []
    for (contrast_id, strength_bin), g in all_wide.groupby(["contrast_id", "strength_bin"], observed=False):
        x = g.set_index("feature").reindex(columns=STUDIES)
        m = _metric(x)
        metrics.append({
            "scale": scale,
            "contrast_id": contrast_id,
            "contrast_label": g["contrast_label"].iloc[0],
            "contrast_class": g["contrast_class"].iloc[0],
            "strength_bin": str(strength_bin),
            "strength_low": float(g["strength"].min()) if len(g) else np.nan,
            "strength_high": float(g["strength"].max()) if len(g) else np.nan,
            **m,
        })
    # Deterministic overlap-window comparison: features in the same absolute
    # score window shared by each endpoint and the cross-sectional contrast.
    delta = all_wide.loc[all_wide["contrast_id"] == "change_d7_minus_d1", "strength"]
    overlap = []
    if not delta.empty:
        lo, hi = float(delta.quantile(0.1)), float(delta.quantile(0.9))
        for contrast_id, label, contrast_class in CONTRASTS:
            g = all_wide.loc[(all_wide["contrast_id"] == contrast_id) & all_wide["strength"].between(lo, hi)]
            x = g.set_index("feature").reindex(columns=STUDIES)
            m = _metric(x)
            overlap.append({
                "scale": scale,
                "contrast_id": contrast_id,
                "contrast_label": label,
                "contrast_class": contrast_class,
                "window_low": lo,
                "window_high": hi,
                **m,
            })
    return pd.DataFrame(metrics), pd.DataFrame(overlap)


def _plot(metrics: pd.DataFrame, out: Path) -> None:
    colors = {"injury_d1_vs_uninjured": "#0072B2", "injury_d7_vs_uninjured": "#009E73", "change_d7_minus_d1": "#D55E00"}
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)
    for ax, (metric, ylabel) in zip(
        axes,
        [("all_study_concordance", "All-study directional concordance"), ("heldout_balanced_accuracy", "Held-out balanced accuracy")],
    ):
        for scale, marker in [("gene", "o"), ("hallmark", "D")]:
            for contrast_id, label, _ in CONTRASTS:
                g = metrics.loc[(metrics["scale"] == scale) & (metrics["contrast_id"] == contrast_id)].copy()
                if g.empty:
                    continue
                x = np.arange(len(g)) + (0.0 if scale == "gene" else 0.14)
                ax.plot(x, g[metric], marker=marker, color=colors[contrast_id], lw=1.4, ms=4,
                        label=f"{scale} · {label}")
        ax.set_xticks(np.arange(4) + 0.07)
        ax.set_xticklabels(["Q1", "Q2", "Q3", "Q4"], fontsize=8)
        ax.set_xlabel("Pooled absolute-score stratum", fontsize=8)
        ax.set_ylabel(ylabel, fontsize=8)
        ax.set_ylim(0, 1.02)
        ax.grid(axis="y", color="#D9D9D9", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.10), ncol=3, frameon=False, fontsize=7)
    fig.suptitle("Effect-strength-stratified transfer is a sensitivity analysis", fontsize=9, y=1.18)
    fig.tight_layout()
    fig.savefig(out.with_suffix(".png"), dpi=600, bbox_inches="tight")
    fig.savefig(out.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(out.with_suffix(".svg"), bbox_inches="tight")
    tiff_path = out.with_suffix(".tif")
    fig.savefig(tiff_path, dpi=600, bbox_inches="tight", pil_kwargs={"compression": "tiff_lzw"})
    with Image.open(tiff_path) as im:
        im.convert("RGB").save(tiff_path, dpi=(600, 600), compression="tiff_lzw")
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    all_metrics, all_overlap = [], []
    for scale in ["gene", "hallmark"]:
        m, o = _analyse_scale(_load_scale(args.root, scale), scale)
        all_metrics.append(m)
        all_overlap.append(o)
    metrics = pd.concat(all_metrics, ignore_index=True)
    overlap = pd.concat(all_overlap, ignore_index=True)
    metrics.to_csv(args.out / "effect_strength_stratified_metrics.csv", index=False)
    overlap.to_csv(args.out / "effect_strength_overlap_window_metrics.csv", index=False)
    summary = (
        metrics.groupby(["scale", "contrast_id"], as_index=False)
        .agg(n_features=("n_features", "sum"), median_strength=("strength_low", "median"), max_strength=("strength_high", "max"))
    )
    summary.to_csv(args.out / "effect_strength_summary.csv", index=False)
    _plot(metrics, args.out / "Fig3_EffectStrengthSensitivity")
    (args.out / "README.md").write_text(
        "# Effect-strength sensitivity\n\n"
        "Post-lock sensitivity analysis of the frozen three-study whole-lesion score matrices. "
        "Features were stratified by pooled absolute mean study score (moderated t for genes; NES for Hallmarks). "
        "The common 10th-90th percentile window of the between-timepoint contrast was also applied to each contrast. "
        "These outputs do not alter the primary discovery matrices or their frozen metrics.\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
