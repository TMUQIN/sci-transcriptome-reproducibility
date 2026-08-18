#!/usr/bin/env python3
"""Align descriptive reproducibility metrics for gene and Hallmark scales.

This reviewer-triggered analysis uses the frozen three-study whole-lesion branch.
It does not test matched random gene sets and must not be used alone to claim that
programs are intrinsically more reproducible than genes.
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
from scipy.stats import spearmanr


SEED = 20260729
STUDIES = ("GSE162610", "GSE234774", "GSE304399")
CONTRASTS = (
    "injury_d1_vs_uninjured",
    "injury_d7_vs_uninjured",
    "change_d7_minus_d1",
)
PROGRAM_EFFECT_TYPES = {
    "injury_d1_vs_uninjured": "endpoint",
    "injury_d7_vs_uninjured": "endpoint",
    "change_d7_minus_d1": "temporal_delta",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False, na_rep="NA", float_format="%.17g", lineterminator="\n")


def sign_array(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, float)
    return np.where(values > 0, 1, np.where(values < 0, -1, 0)).astype(int)


def all_nonzero_same_sign(values: np.ndarray) -> bool:
    signs = sign_array(values)
    return bool(np.all(signs == signs[0]) and signs[0] != 0)


def pairwise_sign_fraction(values: np.ndarray) -> float:
    signs = sign_array(values)
    pairs = list(combinations(range(len(signs)), 2))
    return float(np.mean([signs[i] == signs[j] and signs[i] != 0 for i, j in pairs]))


def wilson_interval(successes: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n <= 0:
        return float("nan"), float("nan")
    p = successes / n
    den = 1 + z * z / n
    center = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return max(0.0, center - half), min(1.0, center + half)


def load_gene_scores(root: Path) -> pd.DataFrame:
    frames = []
    columns = ["gene", "dataset", "contrast_id", "estimate", "se", "fdr"]
    for study in STUDIES:
        path = root / "results" / "whole_lesion" / f"effects_{study}.tsv"
        frame = pd.read_csv(path, sep="\t", usecols=columns)
        frame = frame[frame["contrast_id"].isin(CONTRASTS)].copy()
        frame["score"] = frame["estimate"] / frame["se"]
        frame = frame[np.isfinite(frame["score"]) & np.isfinite(frame["fdr"])].copy()
        if frame.duplicated(["gene", "dataset", "contrast_id"]).any():
            raise ValueError(f"duplicate gene/study/contrast rows in {path}")
        frames.append(frame[["gene", "dataset", "contrast_id", "score", "fdr"]])
    long = pd.concat(frames, ignore_index=True)
    counts = long.groupby(["gene", "contrast_id"], observed=True)["dataset"].nunique()
    complete = counts[counts == len(STUDIES)].index
    index = pd.MultiIndex.from_frame(long[["gene", "contrast_id"]])
    return long[index.isin(complete)].reset_index(drop=True)


def parse_gmt_sizes(path: Path) -> dict[str, int]:
    sizes: dict[str, int] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                sizes[fields[0]] = len(set(fields[2:]))
    return sizes


def load_program_scores(root: Path) -> pd.DataFrame:
    path = root / "results" / "whole_lesion_programs" / "hallmark_gsea_by_study.tsv"
    frame = pd.read_csv(path, sep="\t")
    frame = frame[
        frame["dataset"].isin(STUDIES)
        & frame["contrast_id"].isin(CONTRASTS)
        & frame.apply(lambda r: str(r["effect_type"]) == PROGRAM_EFFECT_TYPES[str(r["contrast_id"])], axis=1)
    ].copy()
    frame = frame.rename(columns={"term": "program", "nes": "score"})
    if frame.duplicated(["program", "dataset", "contrast_id"]).any():
        raise ValueError("duplicate program/study/contrast rows")
    counts = frame.groupby(["program", "contrast_id"], observed=True)["dataset"].nunique()
    complete = counts[counts == len(STUDIES)].index
    index = pd.MultiIndex.from_frame(frame[["program", "contrast_id"]])
    frame = frame[index.isin(complete)].copy()
    sizes = parse_gmt_sizes(root / "references" / "msigdb_mh.all.v2026.1.Mm.symbols.gmt")
    frame["set_size"] = frame["program"].map(sizes)
    return frame[["program", "dataset", "contrast_id", "score", "fdr", "set_size"]].reset_index(drop=True)


def temporal_patterns(long: pd.DataFrame, feature: str) -> pd.DataFrame:
    rows = []
    for feature_value, group in long.groupby(feature, observed=True, sort=False):
        pivot = group.pivot(index="dataset", columns="contrast_id", values="score")
        evaluable = set(STUDIES).issubset(pivot.index) and set(CONTRASTS).issubset(pivot.columns)
        patterns: list[str] = []
        if evaluable:
            for study in STUDIES:
                patterns.append(",".join(str(x) for x in sign_array(pivot.loc[study, list(CONTRASTS)].to_numpy(float))))
        same = bool(evaluable and len(set(patterns)) == 1 and "0" not in patterns[0])
        rows.append({feature: feature_value, "complete_temporal_pattern_evaluable": evaluable,
                     "complete_temporal_pattern_concordant": same,
                     "temporal_sign_patterns_by_study": ";".join(f"{s}:{p}" for s, p in zip(STUDIES, patterns))})
    return pd.DataFrame(rows)


def feature_metrics(long: pd.DataFrame, feature: str, scale: str) -> pd.DataFrame:
    temporal = temporal_patterns(long, feature).set_index(feature)
    rows = []
    for (feature_value, contrast), group in long.groupby([feature, "contrast_id"], observed=True, sort=False):
        group = group.set_index("dataset").loc[list(STUDIES)]
        scores = group["score"].to_numpy(float)
        fdr = group["fdr"].to_numpy(float)
        sig = fdr <= 0.05
        row = {
            "scale": scale,
            "feature": feature_value,
            "contrast_id": contrast,
            "n_studies": len(scores),
            "all_study_signs_concordant": all_nonzero_same_sign(scores),
            "pairwise_sign_concordance": pairwise_sign_fraction(scores),
            "n_study_fdr_le_005": int(sig.sum()),
            "same_direction_fdr_replication_2plus": bool(sig.sum() >= 2 and all_nonzero_same_sign(scores)),
            "same_direction_fdr_replication_all": bool(sig.all() and all_nonzero_same_sign(scores)),
            "mean_abs_signed_score": float(np.mean(np.abs(scores))),
            "across_study_score_sd": float(np.std(scores, ddof=1)),
            "across_study_score_iqr": float(np.subtract(*np.percentile(scores, [75, 25]))),
            "signed_scores_by_study": ";".join(f"{s}:{x:.12g}" for s, x in zip(STUDIES, scores)),
            "fdr_by_study": ";".join(f"{s}:{x:.12g}" for s, x in zip(STUDIES, fdr)),
        }
        if "set_size" in group.columns:
            row["set_size"] = int(group["set_size"].iloc[0]) if pd.notna(group["set_size"].iloc[0]) else pd.NA
        t = temporal.loc[feature_value]
        row.update(t.to_dict())
        rows.append(row)
    return pd.DataFrame(rows)


def rank_matrix(long: pd.DataFrame, feature: str, contrast: str) -> pd.DataFrame:
    pivot = long[long["contrast_id"].eq(contrast)].pivot(index=feature, columns="dataset", values="score")
    pivot = pivot.dropna().loc[:, list(STUDIES)]
    return pivot.rank(axis=0, method="average", pct=True)


def global_metrics(long: pd.DataFrame, feature: str, scale: str) -> pd.DataFrame:
    rows = []
    for contrast in CONTRASTS:
        pivot = long[long["contrast_id"].eq(contrast)].pivot(index=feature, columns="dataset", values="score")
        pivot = pivot.dropna().loc[:, list(STUDIES)]
        signs = sign_array(pivot.to_numpy(float))
        all_sign = np.all(signs == signs[:, [0]], axis=1) & (signs[:, 0] != 0)
        pair_sign = np.mean(
            np.column_stack([(signs[:, i] == signs[:, j]) & (signs[:, i] != 0)
                             for i, j in combinations(range(len(STUDIES)), 2)]), axis=1)
        rank = rank_matrix(long, feature, contrast)
        correlations = []
        heldout_accuracy = []
        loo_rank = []
        for i, j in combinations(range(len(STUDIES)), 2):
            correlations.append(float(spearmanr(pivot.iloc[:, i], pivot.iloc[:, j]).statistic))
        for held_idx, held in enumerate(STUDIES):
            train = [s for s in STUDIES if s != held]
            predicted = sign_array(pivot[train].mean(axis=1).to_numpy(float))
            observed = sign_array(pivot[held].to_numpy(float))
            heldout_accuracy.append(float(np.mean((predicted == observed) & (observed != 0))))
            loo_rank.append(float(spearmanr(rank[train].mean(axis=1), rank[held]).statistic))
        fdr_pivot = long[long["contrast_id"].eq(contrast)].pivot(index=feature, columns="dataset", values="fdr").loc[pivot.index]
        replicated = (fdr_pivot.loc[:, list(STUDIES)].le(0.05).sum(axis=1) >= 2) & all_sign

        def add_rate(name: str, values: np.ndarray, comparable: bool = True) -> None:
            successes = int(np.sum(values))
            low, high = wilson_interval(successes, len(values))
            rows.append({"scale": scale, "contrast_id": contrast, "metric": name,
                         "estimate": successes / len(values), "ci_low": low, "ci_high": high,
                         "n_features": len(values), "directly_comparable": comparable})

        add_rate("all_study_sign_concordance_rate", all_sign)
        rows.append({"scale": scale, "contrast_id": contrast, "metric": "mean_pairwise_sign_concordance",
                     "estimate": float(pair_sign.mean()), "ci_low": float("nan"), "ci_high": float("nan"),
                     "n_features": len(pair_sign), "directly_comparable": True})
        add_rate("same_direction_fdr_replication_2plus_rate", replicated.to_numpy(bool))
        rows.append({"scale": scale, "contrast_id": contrast, "metric": "mean_pairwise_spearman",
                     "estimate": float(np.mean(correlations)), "ci_low": float("nan"), "ci_high": float("nan"),
                     "n_features": len(pivot), "directly_comparable": True})
        rows.append({"scale": scale, "contrast_id": contrast, "metric": "mean_leave_one_study_out_rank_retention",
                     "estimate": float(np.mean(loo_rank)), "ci_low": float("nan"), "ci_high": float("nan"),
                     "n_features": len(pivot), "directly_comparable": True})
        rows.append({"scale": scale, "contrast_id": contrast, "metric": "mean_held_out_direction_accuracy",
                     "estimate": float(np.mean(heldout_accuracy)), "ci_low": float("nan"), "ci_high": float("nan"),
                     "n_features": len(pivot), "directly_comparable": True})
        rows.append({"scale": scale, "contrast_id": contrast, "metric": "median_across_study_score_sd",
                     "estimate": float(np.median(np.std(pivot.to_numpy(float), axis=1, ddof=1))),
                     "ci_low": float("nan"), "ci_high": float("nan"), "n_features": len(pivot),
                     "directly_comparable": False})
    patterns = temporal_patterns(long, feature)
    values = patterns.loc[patterns["complete_temporal_pattern_evaluable"], "complete_temporal_pattern_concordant"].to_numpy(bool)
    low, high = wilson_interval(int(values.sum()), len(values))
    rows.append({"scale": scale, "contrast_id": "complete_temporal_pattern", "metric": "complete_temporal_pattern_concordance_rate",
                 "estimate": float(values.mean()), "ci_low": low, "ci_high": high, "n_features": len(values),
                 "directly_comparable": True})
    return pd.DataFrame(rows)


def equal_feature_count_calibration(gene: pd.DataFrame, program: pd.DataFrame, global_table: pd.DataFrame,
                                    n_iter: int, seed: int) -> pd.DataFrame:
    rng = np.random.default_rng(seed)
    rows = []
    metrics = ("all_study_sign_concordance_rate", "mean_pairwise_sign_concordance",
               "mean_held_out_direction_accuracy")
    for contrast in CONTRASTS:
        gp = gene[gene["contrast_id"].eq(contrast)].pivot(index="gene", columns="dataset", values="score").dropna()
        gp = gp.loc[:, list(STUDIES)]
        pp = program[program["contrast_id"].eq(contrast)].pivot(index="program", columns="dataset", values="score").dropna()
        pp = pp.loc[:, list(STUDIES)]
        n_pick = len(pp)
        gsign = sign_array(gp.to_numpy(float))
        all_sign = (np.all(gsign == gsign[:, [0]], axis=1) & (gsign[:, 0] != 0)).astype(float)
        pair_sign = np.mean(np.column_stack([(gsign[:, i] == gsign[:, j]) & (gsign[:, i] != 0)
                                             for i, j in combinations(range(3), 2)]), axis=1)
        held = []
        for h in range(3):
            train = [i for i in range(3) if i != h]
            held.append((sign_array(gp.iloc[:, train].mean(axis=1).to_numpy(float)) == gsign[:, h]).astype(float))
        held_mean = np.mean(np.column_stack(held), axis=1)
        arrays = dict(zip(metrics, (all_sign, pair_sign, held_mean)))
        indices = np.arange(len(gp))
        for metric, values in arrays.items():
            null = np.empty(n_iter, float)
            for b in range(n_iter):
                null[b] = values[rng.choice(indices, size=n_pick, replace=False)].mean()
            observed = float(global_table[
                global_table["scale"].eq("program") & global_table["contrast_id"].eq(contrast)
                & global_table["metric"].eq(metric)
            ]["estimate"].iloc[0])
            rows.append({"contrast_id": contrast, "metric": metric, "program_observed": observed,
                         "n_program_features": n_pick, "n_gene_features": len(gp), "n_random_gene_feature_subsets": n_iter,
                         "random_gene_subset_median": float(np.median(null)),
                         "random_gene_subset_q025": float(np.quantile(null, 0.025)),
                         "random_gene_subset_q975": float(np.quantile(null, 0.975)),
                         "empirical_p_program_ge_random_gene_features": float((1 + np.sum(null >= observed)) / (n_iter + 1)),
                         "interpretation_limit": "feature-count diagnostic only; not a matched random gene-set null"})
    return pd.DataFrame(rows)


def build_estimability(gene: pd.DataFrame, program: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for scale, frame, feature in (("gene", gene, "gene"), ("program", program, "program")):
        for contrast in CONTRASTS:
            subset = frame[frame["contrast_id"].eq(contrast)]
            rows.append({"branch": "post_primary_three_study_whole_lesion", "scale": scale,
                         "study_set": ";".join(STUDIES), "tissue_compartment": "whole_lesion_bulk_equivalent",
                         "lineage": "not_applicable", "contrast_id": contrast, "estimable": "yes",
                         "n_studies": subset["dataset"].nunique(), "n_complete_features": subset[feature].nunique(),
                         "reason": "complete signed score in all three independent studies"})
    for contrast in ("injury_d1_vs_uninjured", "change_d7_minus_d1"):
        rows.append({"branch": "independent_external_day7_endpoint_stress_test", "scale": "program",
                     "study_set": "GSE304361", "tissue_compartment": "whole_lesion_and_author_cell_types",
                     "lineage": "author_cell_types", "contrast_id": contrast, "estimable": "no",
                     "n_studies": 1, "n_complete_features": 0, "reason": "GSE304361 has no day-1 endpoint"})
    rows.append({"branch": "independent_external_day7_endpoint_stress_test", "scale": "program",
                 "study_set": "GSE304361", "tissue_compartment": "whole_lesion_and_author_cell_types",
                 "lineage": "author_cell_types", "contrast_id": "injury_d7_vs_uninjured", "estimable": "yes",
                 "n_studies": 1, "n_complete_features": 50, "reason": "endpoint stress test only; not temporal replication"})
    return pd.DataFrame(rows)


def build_manifest(root: Path, input_paths: list[Path]) -> pd.DataFrame:
    metadata = {
        "GSE162610": ("GSE162610", "GSM library/mouse", "Milich study family", "lesion_site", "whole_lesion_bulk_equivalent"),
        "GSE234774": ("GSE234774", "GSM library/mouse", "Tabulae Paralytica family", "mid_thoracic", "whole_lesion_bulk_equivalent"),
        "GSE304399": ("GSE304399", "biological library/subject_id", "independent 2026 family", "lesion_site", "whole_lesion_bulk_equivalent"),
    }
    rows = []
    for path in input_paths:
        rel = path.relative_to(root).as_posix()
        dataset = next((x for x in STUDIES if x in path.name), "multiple_or_reference")
        md = metadata.get(dataset, (dataset, "not_applicable", "not_applicable", "not_applicable", "not_applicable"))
        role = "reference_gene_sets" if path.suffix == ".gmt" else ("study_level_effects" if "effects_" in path.name else "pseudobulk_input")
        rows.append({"path": rel, "bytes": path.stat().st_size, "sha256": sha256(path), "dataset": md[0],
                     "biological_sample_definition": md[1], "study_family": md[2], "tissue_compartment": md[3],
                     "lineage": md[4], "contrast_scope": "d1_vs_uninjured;d7_vs_uninjured;direct_d7_minus_d1",
                     "gene_universe": "dataset-filtered modeled genes" if role != "reference_gene_sets" else "MSigDB mouse Hallmark 2026.1",
                     "input_role": role})
    return pd.DataFrame(rows)


def selftest() -> None:
    rows = []
    for feature, multiplier in (("a", 1.0), ("b", -1.0), ("c", 0.5)):
        for contrast_idx, contrast in enumerate(CONTRASTS, start=1):
            for study_idx, study in enumerate(STUDIES, start=1):
                rows.append({"gene": feature, "dataset": study, "contrast_id": contrast,
                             "score": multiplier * contrast_idx * (1 + study_idx / 100), "fdr": 0.01})
    frame = pd.DataFrame(rows)
    metrics = feature_metrics(frame, "gene", "gene")
    assert metrics["all_study_signs_concordant"].all()
    global_frame = global_metrics(frame, "gene", "gene")
    rate = global_frame.loc[global_frame["metric"].eq("all_study_sign_concordance_rate"), "estimate"]
    assert np.allclose(rate, 1.0)
    print("[selftest] PASS: aligned signs, temporal patterns and global metrics")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--out-dir", type=Path, default=Path("reports/phase_reproducibility_calibration_2026_07"))
    parser.add_argument("--iterations", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()
    if args.selftest:
        selftest()
        return
    root = args.root.resolve()
    out_dir = args.out_dir if args.out_dir.is_absolute() else root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    gene = load_gene_scores(root)
    program = load_program_scores(root)
    gene_metrics = feature_metrics(gene, "gene", "gene")
    program_metrics = feature_metrics(program, "program", "program")
    global_table = pd.concat([global_metrics(gene, "gene", "gene"),
                              global_metrics(program, "program", "program")], ignore_index=True)
    comparison = global_table.pivot_table(index=["contrast_id", "metric", "directly_comparable"],
                                          columns="scale", values="estimate").reset_index()
    comparison.columns.name = None
    comparison["program_minus_gene"] = comparison.get("program", np.nan) - comparison.get("gene", np.nan)
    count_calibration = equal_feature_count_calibration(gene, program, global_table, args.iterations, args.seed)
    estimability = build_estimability(gene, program)

    input_paths = [
        *(root / "results" / "whole_lesion" / f"effects_{s}.tsv" for s in STUDIES),
        root / "results" / "whole_lesion_programs" / "hallmark_gsea_by_study.tsv",
        root / "references" / "msigdb_mh.all.v2026.1.Mm.symbols.gmt",
    ]
    manifest = build_manifest(root, list(input_paths))
    outputs = {
        "gene_reproducibility_metrics.tsv": gene_metrics,
        "program_reproducibility_metrics.tsv": program_metrics,
        "gene_program_global_metrics.tsv": global_table,
        "gene_program_metric_comparison.tsv": comparison,
        "equal_feature_count_diagnostic.tsv": count_calibration,
        "estimability_matrix.tsv": estimability,
        "frozen_input_manifest.tsv": manifest,
    }
    for name, frame in outputs.items():
        write_tsv(frame, out_dir / name)

    output_records = []
    for name in outputs:
        path = out_dir / name
        output_records.append({"path": path.relative_to(root).as_posix(), "bytes": path.stat().st_size,
                               "sha256": sha256(path)})
    provenance = {
        "analysis": "aligned gene/program reproducibility metrics",
        "scientific_status": "phase_A_descriptive_calibration; not matched random gene-set inference",
        "created_at": datetime.now().astimezone().isoformat(),
        "seed": args.seed,
        "iterations": args.iterations,
        "python": sys.version,
        "platform": platform.platform(),
        "packages": {"numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__},
        "inputs": manifest.to_dict(orient="records"),
        "outputs": output_records,
        "warnings": [
            "Gene estimate/SE and program NES are different signed-score constructions.",
            "The equal-feature-count diagnostic samples single genes, not matched random gene sets.",
            "No Phase A result alone licenses a program-superiority claim.",
        ],
    }
    provenance_path = out_dir / "phase_A_metric_calibration_provenance.json"
    provenance_path.write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": "ok", "gene_rows": len(gene_metrics), "program_rows": len(program_metrics),
                      "out_dir": str(out_dir)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
