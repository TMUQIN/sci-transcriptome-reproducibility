#!/usr/bin/env python3
"""Post-lock direction-imbalance and feature-identity calibration.

This script reads only frozen three-study signed-score matrices. It does not
recompute differential expression, enrichment, sample scores or meta-analysis.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
from datetime import datetime
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import rankdata


STUDIES = ("GSE162610", "GSE234774", "GSE304399")
CONTRASTS = (
    "injury_d1_vs_uninjured",
    "injury_d7_vs_uninjured",
    "change_d7_minus_d1",
)
ENDPOINTS = CONTRASTS[:2]
DIRECT_CHANGE = CONTRASTS[2]
PROGRAM_EFFECT_TYPES = {
    "injury_d1_vs_uninjured": "endpoint",
    "injury_d7_vs_uninjured": "endpoint",
    "change_d7_minus_d1": "temporal_delta",
}
METRIC_ORDER = (
    "all_study_sign_concordance_rate",
    "mean_pairwise_sign_concordance",
    "mean_pairwise_spearman",
    "mean_held_out_direction_accuracy",
    "mean_pairwise_cohen_kappa",
    "fleiss_kappa",
    "mean_held_out_balanced_accuracy",
    "mean_held_out_mcc",
    "mean_agreeing_train_coverage",
    "mean_agreeing_train_conditional_accuracy",
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, na_rep="NA", float_format="%.17g", lineterminator="\n")


def sign_array(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, float)
    return np.where(values > 0, 1, np.where(values < 0, -1, 0)).astype(np.int8)


def load_gene_scores(root: Path) -> pd.DataFrame:
    frames = []
    usecols = ["gene", "dataset", "contrast_id", "estimate", "se"]
    for study in STUDIES:
        path = root / "results" / "whole_lesion" / f"effects_{study}.tsv"
        frame = pd.read_csv(path, sep="\t", usecols=usecols)
        frame = frame[frame["contrast_id"].isin(CONTRASTS)].copy()
        frame["score"] = frame["estimate"] / frame["se"]
        frame = frame[np.isfinite(frame["score"])].copy()
        if frame.duplicated(["gene", "dataset", "contrast_id"]).any():
            raise ValueError(f"Duplicate gene/study/contrast rows in {path}")
        frames.append(frame[["gene", "dataset", "contrast_id", "score"]])
    long = pd.concat(frames, ignore_index=True)
    counts = long.groupby(["gene", "contrast_id"], observed=True)["dataset"].nunique()
    complete = counts[counts == len(STUDIES)].index
    index = pd.MultiIndex.from_frame(long[["gene", "contrast_id"]])
    return long[index.isin(complete)].reset_index(drop=True)


def load_program_scores(root: Path) -> pd.DataFrame:
    path = root / "results" / "whole_lesion_programs" / "hallmark_gsea_by_study.tsv"
    frame = pd.read_csv(path, sep="\t")
    frame = frame[
        frame["dataset"].isin(STUDIES)
        & frame["contrast_id"].isin(CONTRASTS)
        & frame.apply(lambda row: str(row["effect_type"]) == PROGRAM_EFFECT_TYPES[str(row["contrast_id"])], axis=1)
    ].copy()
    frame = frame.rename(columns={"term": "program", "nes": "score"})
    if frame.duplicated(["program", "dataset", "contrast_id"]).any():
        raise ValueError("Duplicate program/study/contrast rows")
    counts = frame.groupby(["program", "contrast_id"], observed=True)["dataset"].nunique()
    complete = counts[counts == len(STUDIES)].index
    index = pd.MultiIndex.from_frame(frame[["program", "contrast_id"]])
    return frame[index.isin(complete)][["program", "dataset", "contrast_id", "score"]].reset_index(drop=True)


def score_matrix(long: pd.DataFrame, feature: str, contrast: str) -> pd.DataFrame:
    pivot = long[long["contrast_id"].eq(contrast)].pivot(index=feature, columns="dataset", values="score")
    return pivot.dropna().loc[:, list(STUDIES)].sort_index()


def binary_performance(pred: np.ndarray, obs: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return accuracy, balanced accuracy and MCC for B x N sign matrices."""
    pred = np.asarray(pred, np.int8)
    obs = np.asarray(obs, np.int8)
    valid = (pred != 0) & (obs != 0)
    n = obs.shape[1]
    accuracy = np.sum((pred == obs) & (obs != 0), axis=1) / n
    tp = np.sum(valid & (pred == 1) & (obs == 1), axis=1).astype(float)
    tn = np.sum(valid & (pred == -1) & (obs == -1), axis=1).astype(float)
    fp = np.sum(valid & (pred == 1) & (obs == -1), axis=1).astype(float)
    fn = np.sum(valid & (pred == -1) & (obs == 1), axis=1).astype(float)
    with np.errstate(divide="ignore", invalid="ignore"):
        tpr = tp / (tp + fn)
        tnr = tn / (tn + fp)
        balanced = (tpr + tnr) / 2
        denominator = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        mcc = (tp * tn - fp * fn) / denominator
    balanced[~np.isfinite(balanced)] = np.nan
    mcc[~np.isfinite(mcc)] = 0.0
    return accuracy, balanced, mcc


def pairwise_kappa_batch(signs: list[np.ndarray]) -> np.ndarray:
    kappas = []
    for left, right in combinations(signs, 2):
        valid = (left != 0) & (right != 0)
        denom = valid.sum(axis=1).astype(float)
        observed = np.sum(valid & (left == right), axis=1) / denom
        p_left = np.sum(valid & (left == 1), axis=1) / denom
        p_right = np.sum(valid & (right == 1), axis=1) / denom
        expected = p_left * p_right + (1 - p_left) * (1 - p_right)
        with np.errstate(divide="ignore", invalid="ignore"):
            kappa = (observed - expected) / (1 - expected)
        kappa[~np.isfinite(kappa)] = np.nan
        kappas.append(kappa)
    return np.nanmean(np.column_stack(kappas), axis=1)


def fleiss_kappa_batch(signs: list[np.ndarray]) -> np.ndarray:
    stack = np.stack(signs, axis=2)  # B x N x 3
    valid_items = np.all(stack != 0, axis=2)
    n_valid = valid_items.sum(axis=1).astype(float)
    n_pos = np.sum(stack == 1, axis=2)
    n_neg = np.sum(stack == -1, axis=2)
    item_agreement = (n_pos * (n_pos - 1) + n_neg * (n_neg - 1)) / 6.0
    p_bar = np.sum(item_agreement * valid_items, axis=1) / n_valid
    total_ratings = 3.0 * n_valid
    p_pos = np.sum(n_pos * valid_items, axis=1) / total_ratings
    p_neg = np.sum(n_neg * valid_items, axis=1) / total_ratings
    p_expected = p_pos**2 + p_neg**2
    with np.errstate(divide="ignore", invalid="ignore"):
        kappa = (p_bar - p_expected) / (1 - p_expected)
    kappa[~np.isfinite(kappa)] = np.nan
    return kappa


def metrics_batch(scores: list[np.ndarray], standardized_ranks: list[np.ndarray]) -> dict[str, np.ndarray]:
    """Calculate all metrics for aligned B x N arrays."""
    signs = [sign_array(array) for array in scores]
    n = scores[0].shape[1]
    all_concordant = np.mean(
        (signs[0] == signs[1]) & (signs[0] == signs[2]) & (signs[0] != 0), axis=1
    )
    pair_concordance = np.mean(
        np.column_stack([
            np.mean((left == right) & (left != 0), axis=1)
            for left, right in combinations(signs, 2)
        ]),
        axis=1,
    )
    spearman = np.mean(
        np.column_stack([np.sum(left * right, axis=1) for left, right in combinations(standardized_ranks, 2)]),
        axis=1,
    )
    held_accuracy = []
    held_balanced = []
    held_mcc = []
    held_coverage = []
    held_conditional = []
    for held_idx in range(3):
        train = [idx for idx in range(3) if idx != held_idx]
        predicted = sign_array((scores[train[0]] + scores[train[1]]) / 2)
        observed = signs[held_idx]
        accuracy, balanced, mcc = binary_performance(predicted, observed)
        held_accuracy.append(accuracy)
        held_balanced.append(balanced)
        held_mcc.append(mcc)
        covered = (signs[train[0]] == signs[train[1]]) & (signs[train[0]] != 0)
        coverage_count = covered.sum(axis=1).astype(float)
        held_coverage.append(coverage_count / n)
        with np.errstate(divide="ignore", invalid="ignore"):
            conditional = np.sum(covered & (signs[train[0]] == observed) & (observed != 0), axis=1) / coverage_count
        conditional[~np.isfinite(conditional)] = np.nan
        held_conditional.append(conditional)
    return {
        "all_study_sign_concordance_rate": all_concordant,
        "mean_pairwise_sign_concordance": pair_concordance,
        "mean_pairwise_spearman": spearman,
        "mean_held_out_direction_accuracy": np.nanmean(np.column_stack(held_accuracy), axis=1),
        "mean_pairwise_cohen_kappa": pairwise_kappa_batch(signs),
        "fleiss_kappa": fleiss_kappa_batch(signs),
        "mean_held_out_balanced_accuracy": np.nanmean(np.column_stack(held_balanced), axis=1),
        "mean_held_out_mcc": np.nanmean(np.column_stack(held_mcc), axis=1),
        "mean_agreeing_train_coverage": np.nanmean(np.column_stack(held_coverage), axis=1),
        "mean_agreeing_train_conditional_accuracy": np.nanmean(np.column_stack(held_conditional), axis=1),
    }


def standardized_rank(values: np.ndarray) -> np.ndarray:
    ranks = rankdata(values, method="average")
    ranks = ranks - ranks.mean()
    norm = np.sqrt(np.sum(ranks**2))
    if norm == 0:
        raise ValueError("Cannot standardize a constant rank vector")
    return ranks / norm


def observed_metrics(pivot: pd.DataFrame) -> dict[str, float]:
    scores = [pivot.iloc[:, idx].to_numpy(float)[None, :] for idx in range(3)]
    ranks = [standardized_rank(pivot.iloc[:, idx].to_numpy(float))[None, :] for idx in range(3)]
    return {name: float(values[0]) for name, values in metrics_batch(scores, ranks).items()}


def permutation_metrics(
    pivot: pd.DataFrame,
    n_iter: int,
    rng: np.random.Generator,
    batch_size: int,
) -> dict[str, np.ndarray]:
    vectors = [pivot.iloc[:, idx].to_numpy(float) for idx in range(3)]
    rank_vectors = [standardized_rank(vector) for vector in vectors]
    n = len(pivot)
    outputs = {metric: np.empty(n_iter, float) for metric in METRIC_ORDER}
    start = 0
    while start < n_iter:
        size = min(batch_size, n_iter - start)
        permutations = [np.vstack([rng.permutation(n) for _ in range(size)]) for _ in range(3)]
        scores = [vectors[idx][permutations[idx]] for idx in range(3)]
        ranks = [rank_vectors[idx][permutations[idx]] for idx in range(3)]
        batch = metrics_batch(scores, ranks)
        for metric in METRIC_ORDER:
            outputs[metric][start:start + size] = batch[metric]
        start += size
    return outputs


def study_prevalence(pivot: pd.DataFrame, scale: str, contrast: str) -> list[dict[str, object]]:
    rows = []
    for study in STUDIES:
        signs = sign_array(pivot[study].to_numpy(float))
        rows.append({
            "scale": scale,
            "contrast_id": contrast,
            "dataset": study,
            "n_features": len(signs),
            "n_positive": int(np.sum(signs == 1)),
            "n_negative": int(np.sum(signs == -1)),
            "n_zero": int(np.sum(signs == 0)),
            "positive_proportion": float(np.mean(signs == 1)),
            "negative_proportion": float(np.mean(signs == -1)),
            "zero_proportion": float(np.mean(signs == 0)),
        })
    return rows


def audit_frozen_metric_match(root: Path, observed: pd.DataFrame) -> pd.DataFrame:
    frozen = pd.read_csv(root / "reports/phase_reproducibility_calibration_2026_07/gene_program_global_metrics.tsv", sep="\t")
    names = {
        "all_study_sign_concordance_rate",
        "mean_pairwise_sign_concordance",
        "mean_pairwise_spearman",
        "mean_held_out_direction_accuracy",
    }
    left = observed[observed["metric"].isin(names)][["scale", "contrast_id", "metric", "observed"]]
    right = frozen[frozen["metric"].isin(names)][["scale", "contrast_id", "metric", "estimate"]]
    merged = left.merge(right, on=["scale", "contrast_id", "metric"], how="outer", validate="one_to_one")
    merged["absolute_difference"] = (merged["observed"] - merged["estimate"]).abs()
    merged["matches_within_1e_12"] = merged["absolute_difference"] <= 1e-12
    if len(merged) != 24 or not merged["matches_within_1e_12"].all():
        raise ValueError("Post-lock implementation does not reproduce all frozen global metrics")
    return merged


def classify(program_summary: pd.DataFrame) -> tuple[str, dict[str, object]]:
    table = program_summary.set_index(["contrast_id", "metric"])
    identity_metrics = (
        "mean_pairwise_spearman",
        "mean_held_out_balanced_accuracy",
        "mean_held_out_mcc",
        "mean_agreeing_train_conditional_accuracy",
    )
    imbalance_metrics = (
        "mean_pairwise_cohen_kappa",
        "mean_held_out_balanced_accuracy",
        "mean_held_out_mcc",
        "mean_agreeing_train_conditional_accuracy",
    )

    def p(contrast: str, metric: str) -> float:
        return float(table.loc[(contrast, metric), "empirical_p_null_ge_observed"])

    def value(contrast: str, metric: str) -> float:
        return float(table.loc[(contrast, metric), "observed"])

    direct_identity_count = sum(p(DIRECT_CHANGE, metric) <= 0.05 for metric in identity_metrics)
    endpoint_details = {}
    strengthened = True
    weakened = False
    for endpoint in ENDPOINTS:
        identity_count = sum(p(endpoint, metric) <= 0.05 for metric in identity_metrics)
        higher_count = sum(value(endpoint, metric) > value(DIRECT_CHANGE, metric) for metric in imbalance_metrics)
        raw_higher = value(endpoint, "all_study_sign_concordance_rate") > value(DIRECT_CHANGE, "all_study_sign_concordance_rate")
        all_p = p(endpoint, "all_study_sign_concordance_rate") <= 0.05
        endpoint_details[endpoint] = {
            "identity_metric_p_le_005_count": identity_count,
            "imbalance_metric_higher_than_direct_count": higher_count,
            "raw_all_study_concordance_higher_than_direct": raw_higher,
            "all_study_concordance_empirical_p_le_005": all_p,
        }
        strengthened &= all_p and identity_count >= 2 and higher_count >= 3
        weakened |= (not raw_higher) or higher_count <= 1
    strengthened &= direct_identity_count < 2
    if strengthened:
        decision = "A_strengthened"
    elif weakened:
        decision = "C_weakened"
    else:
        decision = "B_prevalence_qualified"
    return decision, {"direct_identity_metric_p_le_005_count": direct_identity_count, "endpoints": endpoint_details}


def build_report(
    decision: str,
    details: dict[str, object],
    prevalence: pd.DataFrame,
    summary: pd.DataFrame,
) -> str:
    program = summary[summary["scale"].eq("program")]
    lines = [
        "# Post-lock direction-imbalance and feature-identity calibration audit",
        "",
        f"Decision: **{decision}**.",
        "",
        "This sensitivity branch used only the frozen three-study signed-score matrices. It did not recompute or alter any v5.3.2 primary model, threshold, dataset decision or inferential result.",
        "",
        "## Hallmark-scale results",
        "",
    ]
    for contrast in CONTRASTS:
        subset = program[program["contrast_id"].eq(contrast)].set_index("metric")
        lines.append(
            f"For `{contrast}`, all-study sign concordance was {subset.loc['all_study_sign_concordance_rate', 'observed']:.3f} "
            f"(identity-permutation P={subset.loc['all_study_sign_concordance_rate', 'empirical_p_null_ge_observed']:.4g}), "
            f"mean pairwise kappa was {subset.loc['mean_pairwise_cohen_kappa', 'observed']:.3f}, "
            f"held-out balanced accuracy was {subset.loc['mean_held_out_balanced_accuracy', 'observed']:.3f}, "
            f"MCC was {subset.loc['mean_held_out_mcc', 'observed']:.3f}, and agreement-gated conditional accuracy was "
            f"{subset.loc['mean_agreeing_train_conditional_accuracy', 'observed']:.3f}."
        )
        lines.append("")
    lines.extend([
        "## Direction prevalence",
        "",
        "The study-specific positive and negative proportions are reported in `study_direction_prevalence.tsv`. They define the global sign-prevalence structure preserved by every permutation.",
        "",
        "## Prespecified interpretation",
        "",
        json.dumps(details, ensure_ascii=False, indent=2),
        "",
        "The feature-identity null tests whether the same feature transfers across studies beyond marginal sign prevalence. It is not a correlation-matched competitive gene-set null and cannot establish that curated Hallmarks outperform matched synthetic sets.",
        "",
        "All observed and null summaries are retained irrespective of the decision category.",
    ])
    return "\n".join(lines) + "\n"


def selftest() -> None:
    n = 100
    base = np.tile(np.array([-1.0, 1.0]), n // 2)
    pivot = pd.DataFrame({study: base * (1 + idx / 10) for idx, study in enumerate(STUDIES)})
    result = observed_metrics(pivot)
    expected_one = [
        "all_study_sign_concordance_rate",
        "mean_pairwise_sign_concordance",
        "mean_held_out_direction_accuracy",
        "mean_pairwise_cohen_kappa",
        "fleiss_kappa",
        "mean_held_out_balanced_accuracy",
        "mean_held_out_mcc",
        "mean_agreeing_train_coverage",
        "mean_agreeing_train_conditional_accuracy",
    ]
    if not all(abs(result[name] - 1.0) < 1e-12 for name in expected_one):
        raise AssertionError(result)
    print("[selftest] PASS")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--out-dir", type=Path, default=Path("reports/phase_v5_3_3_postlock_upgrade_2026_08/sign_identity"))
    parser.add_argument("--iterations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20260803)
    parser.add_argument("--batch-size", type=int, default=50)
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()
    if args.selftest:
        selftest()
        return
    if args.iterations != 10000 or args.seed != 20260803:
        raise ValueError("The frozen protocol requires 10,000 iterations and seed 20260803")
    root = args.root.resolve()
    out = args.out_dir if args.out_dir.is_absolute() else root / args.out_dir
    out.mkdir(parents=True, exist_ok=True)

    gene = load_gene_scores(root)
    program = load_program_scores(root)
    if gene.groupby("contrast_id")["gene"].nunique().to_dict() != {contrast: 15331 for contrast in CONTRASTS}:
        raise ValueError("Frozen complete gene universe is not 15,331 for every contrast")
    if program.groupby("contrast_id")["program"].nunique().to_dict() != {contrast: 50 for contrast in CONTRASTS}:
        raise ValueError("Frozen complete Hallmark universe is not 50 for every contrast")

    rng = np.random.default_rng(args.seed)
    prevalence_rows = []
    observed_rows = []
    summary_rows = []
    null_archive: dict[str, np.ndarray] = {}
    for scale, long, feature in (("gene", gene, "gene"), ("program", program, "program")):
        for contrast in CONTRASTS:
            pivot = score_matrix(long, feature, contrast)
            prevalence_rows.extend(study_prevalence(pivot, scale, contrast))
            observed = observed_metrics(pivot)
            null = permutation_metrics(pivot, args.iterations, rng, args.batch_size)
            for metric in METRIC_ORDER:
                values = null[metric]
                key = f"{scale}__{contrast}__{metric}"
                null_archive[key] = values
                empirical = float((1 + np.sum(values >= observed[metric])) / (args.iterations + 1))
                observed_rows.append({
                    "scale": scale,
                    "contrast_id": contrast,
                    "metric": metric,
                    "observed": observed[metric],
                    "n_features": len(pivot),
                })
                summary_rows.append({
                    "scale": scale,
                    "contrast_id": contrast,
                    "metric": metric,
                    "observed": observed[metric],
                    "n_features": len(pivot),
                    "n_permutations": args.iterations,
                    "null_median": float(np.nanmedian(values)),
                    "null_q025": float(np.nanquantile(values, 0.025)),
                    "null_q975": float(np.nanquantile(values, 0.975)),
                    "empirical_p_null_ge_observed": empirical,
                    "null_interpretation": "feature identity permuted within study; score and sign-prevalence margins preserved",
                })

    prevalence = pd.DataFrame(prevalence_rows)
    observed_frame = pd.DataFrame(observed_rows)
    summary = pd.DataFrame(summary_rows)
    frozen_match = audit_frozen_metric_match(root, observed_frame)
    decision, decision_details = classify(summary[summary["scale"].eq("program")])

    write_tsv(prevalence, out / "study_direction_prevalence.tsv")
    write_tsv(observed_frame, out / "observed_identity_transfer_metrics.tsv")
    write_tsv(summary, out / "feature_identity_permutation_summary.tsv")
    write_tsv(frozen_match, out / "frozen_metric_reproduction_audit.tsv")
    np.savez_compressed(out / "feature_identity_permutation_null_draws.npz", **null_archive)
    report = build_report(decision, decision_details, prevalence, summary)
    (out / "sign_identity_calibration_audit.md").write_text(report, encoding="utf-8")

    inputs = [
        *(root / "results" / "whole_lesion" / f"effects_{study}.tsv" for study in STUDIES),
        root / "results" / "whole_lesion_programs" / "hallmark_gsea_by_study.tsv",
        root / "reports" / "phase_reproducibility_calibration_2026_07" / "gene_program_global_metrics.tsv",
        root / "reports" / "phase_v5_3_3_postlock_upgrade_2026_08" / "01_sign_identity_calibration_protocol.md",
    ]
    outputs = sorted(path for path in out.iterdir() if path.is_file())
    provenance = {
        "analysis": "post-lock direction-imbalance and feature-identity calibration",
        "created_at": datetime.now().astimezone().isoformat(),
        "status": decision,
        "scientific_scope": "sensitivity branch; v5.3.2 frozen analyses not recomputed",
        "seed": args.seed,
        "iterations": args.iterations,
        "batch_size": args.batch_size,
        "python": sys.version,
        "platform": platform.platform(),
        "packages": {"numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__},
        "inputs": [{"path": path.relative_to(root).as_posix(), "bytes": path.stat().st_size, "sha256": sha256(path)} for path in inputs],
        "outputs": [{"path": path.relative_to(root).as_posix(), "bytes": path.stat().st_size, "sha256": sha256(path)} for path in outputs],
        "warnings": [
            "Gene moderated statistics and Hallmark normalized enrichment scores are evaluated within scale.",
            "The identity-permutation null preserves study-level score and sign-prevalence margins.",
            "This is not a correlation-matched competitive gene-set null and cannot establish curated-set advantage.",
        ],
    }
    (out / "sign_identity_calibration_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": decision,
        "gene_features": 15331,
        "program_features": 50,
        "iterations": args.iterations,
        "out_dir": str(out),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
