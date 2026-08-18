#!/usr/bin/env python3
"""Three-level matched random gene-set pilot for four focal Hallmarks.

This is an exchangeability and compute-budget pilot. It must pass the frozen
diagnostics before expansion to 10,000 sets per Hallmark/contrast/null level.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.util
import json
import math
import platform
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.io import mmread
from scipy.stats import rankdata
from sklearn.cluster import KMeans
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import StandardScaler


SEED = 20260729
STUDIES = ("GSE162610", "GSE234774", "GSE304399")
CONTRASTS = {
    "injury_d1_vs_uninjured": ("dpi_1", "uninjured"),
    "injury_d7_vs_uninjured": ("dpi_7", "uninjured"),
    "change_d7_minus_d1": ("dpi_7", "dpi_1"),
}
FOCAL = (
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "HALLMARK_HYPOXIA",
    "HALLMARK_MYC_TARGETS_V1",
    "HALLMARK_MTORC1_SIGNALING",
)
METHODS = ("magnitude_mean_z", "centered_rank_mean")


def bh_adjust(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, float)
    order = np.argsort(values)
    ranked = values[order] * len(values) / np.arange(1, len(values) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty_like(ranked)
    out[order] = np.clip(ranked, 0, 1)
    return out


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def parse_gmt(path: Path) -> dict[str, set[str]]:
    result = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            result[fields[0]] = set(fields[2:])
    return result


def load_study(root: Path, study: str) -> dict:
    base = root / "data_processed" / "whole_lesion"
    matrix_path = base / f"{study}_pseudobulk_counts.mtx.gz"
    genes_path = base / f"{study}_pseudobulk_genes.tsv"
    coldata_path = base / f"{study}_pseudobulk_coldata.tsv"
    with gzip.open(matrix_path, "rb") as handle:
        counts = mmread(handle)
    counts = counts.toarray() if hasattr(counts, "toarray") else np.asarray(counts)
    genes = pd.read_csv(genes_path, sep="\t")["gene"].astype(str).tolist()
    meta = pd.read_csv(coldata_path, sep="\t")
    keep = meta["group"].isin({"uninjured", "dpi_1", "dpi_7"}).to_numpy()
    counts = counts[keep].astype(float, copy=False)
    meta = meta.loc[keep].reset_index(drop=True)
    if counts.shape != (len(meta), len(genes)) or meta["subject_id"].duplicated().any():
        raise ValueError(f"{study}: invalid biological-sample matrix contract")
    lib = counts.sum(axis=1)
    logcpm = np.log2(counts / lib[:, None] * 1_000_000 + 0.5)
    gene_sd = logcpm.std(axis=0, ddof=1)
    magnitude = np.zeros_like(logcpm)
    variable = gene_sd > 0
    magnitude[:, variable] = (logcpm[:, variable] - logcpm[:, variable].mean(axis=0)) / gene_sd[variable]
    ranks = np.apply_along_axis(rankdata, 1, logcpm, method="average")
    centered_rank = 2 * (ranks - (len(genes) + 1) / 2) / (len(genes) - 1)
    return {
        "counts": counts, "logcpm": logcpm, "magnitude": magnitude,
        "centered_rank": centered_rank, "genes": genes,
        "gene_index": {g: i for i, g in enumerate(genes)}, "meta": meta,
        "input_paths": [matrix_path, genes_path, coldata_path],
    }


def common_gene_universe(root: Path) -> tuple[list[str], dict[str, pd.DataFrame], list[Path]]:
    frames = {}
    inputs = []
    common = None
    usecols = ["gene", "dataset", "contrast_id", "estimate", "se"]
    for study in STUDIES:
        path = root / "results" / "whole_lesion" / f"effects_{study}.tsv"
        frame = pd.read_csv(path, sep="\t", usecols=usecols)
        frame = frame[frame["contrast_id"].isin(CONTRASTS)].copy()
        frame["z"] = frame["estimate"] / frame["se"]
        frames[study] = frame
        genes = set(frame.groupby("gene")["contrast_id"].nunique().loc[lambda x: x == 3].index)
        common = genes if common is None else common & genes
        inputs.append(path)
    return sorted(common), frames, inputs


def residual_profile(study_data: dict, common: list[str]) -> np.ndarray:
    index = np.array([study_data["gene_index"][g] for g in common])
    x = study_data["logcpm"][:, index].copy()
    groups = study_data["meta"]["group"].to_numpy()
    for group in np.unique(groups):
        mask = groups == group
        x[mask] -= x[mask].mean(axis=0)
    sd = x.std(axis=0, ddof=1)
    valid = sd > 1e-12
    x[:, valid] /= sd[valid]
    x[:, ~valid] = 0
    return x


def build_gene_features(common: list[str], studies: dict[str, dict], effects: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, np.ndarray]:
    feature = pd.DataFrame(index=common)
    means, variances, detects, residual_blocks = [], [], [], []
    for study in STUDIES:
        data = studies[study]
        idx = np.array([data["gene_index"][g] for g in common])
        means.append(data["logcpm"][:, idx].mean(axis=0))
        variances.append(np.log1p(data["logcpm"][:, idx].var(axis=0, ddof=1)))
        detects.append((data["counts"][:, idx] > 0).mean(axis=0))
        residual_blocks.append(residual_profile(data, common))
    feature["mean_expression"] = np.mean(means, axis=0)
    feature["log_expression_variance"] = np.mean(variances, axis=0)
    feature["detectability"] = np.mean(detects, axis=0)
    for contrast in CONTRASTS:
        se_mat, z_mat = [], []
        for study in STUDIES:
            sub = effects[study][effects[study]["contrast_id"].eq(contrast)].set_index("gene").loc[common]
            se_mat.append(np.log(sub["se"].to_numpy(float)))
            z_mat.append(sub["z"].to_numpy(float))
        z = np.column_stack(z_mat)
        signs = np.sign(z)
        feature[f"mean_log_se__{contrast}"] = np.mean(se_mat, axis=0)
        feature[f"positive_study_count__{contrast}"] = (signs > 0).sum(axis=1)
        feature[f"sign_concordant__{contrast}"] = (np.all(signs == signs[:, [0]], axis=1) & (signs[:, 0] != 0)).astype(int)
        feature[f"estimable_studies__{contrast}"] = 3
    residual = np.vstack(residual_blocks)
    residual -= residual.mean(axis=0)
    sd = residual.std(axis=0, ddof=1)
    valid = sd > 1e-12
    residual[:, valid] /= sd[valid]
    residual[:, ~valid] = 0
    return feature, residual


def correlation_structure(indices: np.ndarray, residual: np.ndarray, modules: np.ndarray, n_modules: int) -> dict:
    x = residual[:, indices]
    n, m = x.shape
    gram_small = x @ x.T / max(n - 1, 1)
    sum_vec = x.sum(axis=1)
    sum_corr = (float(sum_vec @ sum_vec) / max(n - 1, 1) - m) / 2
    mean_corr = 2 * sum_corr / max(m * (m - 1), 1)
    sum_r2 = float(np.sum(gram_small * gram_small))
    effective_size = m * m / max(sum_r2, 1e-12)
    module_prop = np.bincount(modules[indices], minlength=n_modules) / m
    return {"mean_pairwise_residual_correlation": mean_corr,
            "effective_gene_set_size": effective_size, "module_proportions": module_prop}


def build_neighbors(feature: pd.DataFrame, columns: list[str], targets: np.ndarray, n_neighbors: int) -> tuple[np.ndarray, np.ndarray]:
    scaled = StandardScaler().fit_transform(feature[columns].to_numpy(float))
    model = NearestNeighbors(n_neighbors=min(n_neighbors, len(feature)), algorithm="auto").fit(scaled)
    distances, neighbors = model.kneighbors(scaled[targets])
    return neighbors, distances


def sample_matched_sets(neighbors: np.ndarray, distances: np.ndarray, excluded: set[int], n_sets: int,
                        rng: np.random.Generator, max_attempts_factor: int = 30) -> tuple[list[np.ndarray], int]:
    results, seen = [], set()
    attempts = 0
    target_order = np.arange(len(neighbors))
    max_attempts = n_sets * max_attempts_factor
    while len(results) < n_sets and attempts < max_attempts:
        attempts += 1
        chosen: set[int] = set()
        success = True
        for t in rng.permutation(target_order):
            candidates = neighbors[t]
            d = distances[t]
            ok = np.array([c not in excluded and c not in chosen for c in candidates])
            if not ok.any():
                success = False
                break
            candidates = candidates[ok]
            d = d[ok]
            scale = max(float(np.median(d[d > 0])) if np.any(d > 0) else 1.0, 1e-6)
            weights = np.exp(-0.5 * (d / scale) ** 2)
            weights /= weights.sum()
            chosen.add(int(rng.choice(candidates, p=weights)))
        if not success:
            continue
        key = tuple(sorted(chosen))
        if len(key) == len(neighbors) and key not in seen:
            seen.add(key)
            results.append(np.asarray(key, dtype=int))
    return results, attempts


def set_effects(indices: np.ndarray, studies: dict[str, dict], common: list[str], contrast: str) -> dict[str, dict]:
    exp_group, ref_group = CONTRASTS[contrast]
    result = {}
    for method in METHODS:
        effects = []
        for study in STUDIES:
            data = studies[study]
            dataset_idx = np.array([data["gene_index"][common[i]] for i in indices])
            matrix = data["magnitude"] if method == "magnitude_mean_z" else data["centered_rank"]
            score = matrix[:, dataset_idx].mean(axis=1)
            sd = score.std(ddof=1)
            score = (score - score.mean()) / sd if sd > 0 else np.zeros_like(score)
            group = data["meta"]["group"].to_numpy()
            effects.append(float(score[group == exp_group].mean() - score[group == ref_group].mean()))
        effects = np.asarray(effects)
        signs = np.sign(effects)
        directional = abs(float(effects.mean())) / max(float(np.sqrt(np.mean(effects ** 2))), 1e-12)
        pairwise = np.mean([signs[i] == signs[j] and signs[i] != 0 for i in range(3) for j in range(i + 1, 3)])
        heldout = []
        for h in range(3):
            train = [i for i in range(3) if i != h]
            heldout.append(np.sign(effects[train].mean()) == signs[h] and signs[h] != 0)
        result[method] = {
            "directional_consistency": directional,
            "all_study_sign_concordance": bool(np.all(signs == signs[0]) and signs[0] != 0),
            "pairwise_sign_concordance": float(pairwise),
            "held_out_direction_accuracy": float(np.mean(heldout)),
            "study_effects": effects,
        }
    return result


def row_for_set(program: str, contrast: str, level: str, set_id: int, indices: np.ndarray,
                common: list[str], studies: dict[str, dict], structure: dict) -> dict:
    metrics = set_effects(indices, studies, common, contrast)
    row = {"program": program, "contrast_id": contrast, "null_level": level, "set_id": set_id,
           "set_size": len(indices), "set_hash": hashlib.sha256(";".join(common[i] for i in indices).encode()).hexdigest(),
           "members": ";".join(common[i] for i in indices),
           "mean_pairwise_residual_correlation": structure["mean_pairwise_residual_correlation"],
           "effective_gene_set_size": structure["effective_gene_set_size"]}
    for method, values in metrics.items():
        prefix = method + "__"
        row[prefix + "directional_consistency"] = values["directional_consistency"]
        row[prefix + "all_study_sign_concordance"] = values["all_study_sign_concordance"]
        row[prefix + "pairwise_sign_concordance"] = values["pairwise_sign_concordance"]
        row[prefix + "held_out_direction_accuracy"] = values["held_out_direction_accuracy"]
        row[prefix + "study_effects"] = ";".join(f"{s}:{e:.12g}" for s, e in zip(STUDIES, values["study_effects"]))
    return row


def exchangeability_rows(program: str, contrast: str, level: str, target_idx: np.ndarray, sets: list[np.ndarray],
                         feature: pd.DataFrame, residual: np.ndarray, modules: np.ndarray, n_modules: int) -> list[dict]:
    rows = []
    continuous = ["mean_expression", "log_expression_variance", "detectability"]
    if level in {"level2", "level3"}:
        continuous += [f"mean_log_se__{contrast}", f"positive_study_count__{contrast}", f"sign_concordant__{contrast}"]
    for covariate in continuous:
        obs = feature.iloc[target_idx][covariate].to_numpy(float)
        null_pooled = np.concatenate([feature.iloc[x][covariate].to_numpy(float) for x in sets])
        pooled_sd = math.sqrt((np.var(obs, ddof=1) + np.var(null_pooled, ddof=1)) / 2)
        smd = (float(np.mean(obs)) - float(np.mean(null_pooled))) / pooled_sd if pooled_sd > 0 else 0.0
        rows.append({"program": program, "contrast_id": contrast, "null_level": level, "diagnostic": covariate,
                     "observed": float(np.mean(obs)), "null_median": float(np.median([feature.iloc[x][covariate].mean() for x in sets])),
                     "difference": float(np.mean(obs) - np.mean(null_pooled)), "standardized_difference": smd,
                     "threshold": 0.25, "pass": abs(smd) <= 0.25})
    target_s = correlation_structure(target_idx, residual, modules, n_modules)
    null_s = [correlation_structure(x, residual, modules, n_modules) for x in sets]
    corr_median = float(np.median([x["mean_pairwise_residual_correlation"] for x in null_s]))
    eff_median = float(np.median([x["effective_gene_set_size"] for x in null_s]))
    module_l1 = float(np.mean([np.abs(x["module_proportions"] - target_s["module_proportions"]).sum() for x in null_s]))
    corr_diff = target_s["mean_pairwise_residual_correlation"] - corr_median
    eff_rel = (target_s["effective_gene_set_size"] - eff_median) / max(target_s["effective_gene_set_size"], 1e-12)
    rows.extend([
        {"program": program, "contrast_id": contrast, "null_level": level,
         "diagnostic": "mean_pairwise_residual_correlation", "observed": target_s["mean_pairwise_residual_correlation"],
         "null_median": corr_median, "difference": corr_diff, "standardized_difference": np.nan,
         "threshold": 0.05, "pass": abs(corr_diff) <= 0.05},
        {"program": program, "contrast_id": contrast, "null_level": level,
         "diagnostic": "effective_gene_set_size_relative_difference", "observed": target_s["effective_gene_set_size"],
         "null_median": eff_median, "difference": eff_rel, "standardized_difference": np.nan,
         "threshold": 0.20, "pass": abs(eff_rel) <= 0.20},
        {"program": program, "contrast_id": contrast, "null_level": level,
         "diagnostic": "coexpression_module_composition_l1", "observed": 0.0,
         "null_median": module_l1, "difference": module_l1, "standardized_difference": np.nan,
         "threshold": 0.20, "pass": module_l1 <= 0.20},
    ])
    return rows


def summarize_lift(observed: pd.DataFrame, null_tables: dict[str, pd.DataFrame]) -> pd.DataFrame:
    metrics = ("directional_consistency", "all_study_sign_concordance", "held_out_direction_accuracy")
    rows = []
    for level, null in null_tables.items():
        for program in FOCAL:
            for contrast in CONTRASTS:
                obs = observed[(observed["program"] == program) & (observed["contrast_id"] == contrast)].iloc[0]
                sub = null[(null["program"] == program) & (null["contrast_id"] == contrast)]
                for method in METHODS:
                    for metric in metrics:
                        column = method + "__" + metric
                        observed_value = float(obs[column])
                        values = sub[column].astype(float).to_numpy()
                        rows.append({"program": program, "contrast_id": contrast, "null_level": level,
                                     "method": method, "metric": metric, "observed": observed_value,
                                     "null_median": float(np.median(values)), "null_q025": float(np.quantile(values, 0.025)),
                                     "null_q975": float(np.quantile(values, 0.975)), "n_null": len(values),
                                     "empirical_p": float((1 + np.sum(values >= observed_value)) / (len(values) + 1))})
    result = pd.DataFrame(rows)
    result["fdr"] = np.nan
    for _, idx in result.groupby(["null_level", "method", "metric"], observed=True).groups.items():
        result.loc[idx, "fdr"] = bh_adjust(result.loc[idx, "empirical_p"].to_numpy(float))
    return result


def selftest() -> None:
    e = np.array([1.0, 1.0, 1.0])
    d = abs(e.mean()) / np.sqrt(np.mean(e ** 2))
    assert np.isclose(d, 1.0)
    e2 = np.array([1.0, -1.0, 1.0])
    assert abs(e2.mean()) / np.sqrt(np.mean(e2 ** 2)) < 0.5
    assert np.allclose(bh_adjust(np.array([0.01, 0.04, 0.2])), [0.03, 0.06, 0.2])
    print("[selftest] PASS: directional consistency and BH")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--out-dir", type=Path, default=Path("reports/phase_final_calibration_and_v5_2026_07"))
    parser.add_argument("--pilot-sets", type=int, default=500)
    parser.add_argument("--level3-candidates", type=int, default=2500)
    parser.add_argument("--neighbors", type=int, default=250)
    parser.add_argument("--modules", type=int, default=20)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()
    if args.selftest:
        selftest(); return
    root = args.root.resolve()
    out_dir = args.out_dir if args.out_dir.is_absolute() else root / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(args.seed)
    common, effect_frames, effect_paths = common_gene_universe(root)
    studies = {s: load_study(root, s) for s in STUDIES}
    feature, residual = build_gene_features(common, studies, effect_frames)
    clustering = KMeans(n_clusters=args.modules, random_state=args.seed, n_init=10).fit(residual.T)
    modules = clustering.labels_
    gmt_path = root / "references" / "msigdb_mh.all.v2026.1.Mm.symbols.gmt"
    gene_sets = parse_gmt(gmt_path)
    common_index = {g: i for i, g in enumerate(common)}

    null_rows = {"level1": [], "level2": [], "level3": []}
    diagnostic_rows, observed_rows, failure_lines = [], [], []
    for program in FOCAL:
        target = np.array(sorted(common_index[g] for g in gene_sets[program] if g in common_index), dtype=int)
        excluded = set(target.tolist())
        target_structure = correlation_structure(target, residual, modules, args.modules)
        for contrast in CONTRASTS:
            observed_rows.append(row_for_set(program, contrast, "observed", 0, target, common, studies, target_structure))
            level_specs = {
                "level1": ["mean_expression", "log_expression_variance", "detectability"],
                "level2": ["mean_expression", "log_expression_variance", "detectability",
                           f"mean_log_se__{contrast}", f"positive_study_count__{contrast}",
                           f"sign_concordant__{contrast}"],
            }
            generated = {}
            for level, columns in level_specs.items():
                neighbors, distances = build_neighbors(feature, columns, target, args.neighbors)
                sets, attempts = sample_matched_sets(neighbors, distances, excluded, args.pilot_sets, rng)
                generated[level] = sets
                if len(sets) < args.pilot_sets:
                    failure_lines.append(f"- {program} | {contrast} | {level}: generated {len(sets)}/{args.pilot_sets} after {attempts} attempts.")
                for set_id, indices in enumerate(sets, 1):
                    structure = correlation_structure(indices, residual, modules, args.modules)
                    null_rows[level].append(row_for_set(program, contrast, level, set_id, indices, common, studies, structure))
                diagnostic_rows.extend(exchangeability_rows(program, contrast, level, target, sets, feature,
                                                             residual, modules, args.modules))

            # Level 3: generate a larger Level-2-matched candidate pool and retain the best structure matches.
            columns = level_specs["level2"]
            neighbors, distances = build_neighbors(feature, columns, target, args.neighbors)
            candidates, attempts = sample_matched_sets(neighbors, distances, excluded, args.level3_candidates, rng, 50)
            scored = []
            for indices in candidates:
                s = correlation_structure(indices, residual, modules, args.modules)
                corr_d = abs(s["mean_pairwise_residual_correlation"] - target_structure["mean_pairwise_residual_correlation"]) / 0.05
                eff_d = abs(s["effective_gene_set_size"] - target_structure["effective_gene_set_size"]) / max(0.20 * target_structure["effective_gene_set_size"], 1e-12)
                mod_d = np.abs(s["module_proportions"] - target_structure["module_proportions"]).sum() / 0.20
                scored.append((corr_d + eff_d + mod_d, indices, s))
            scored.sort(key=lambda x: x[0])
            selected = scored[:args.pilot_sets]
            if len(selected) < args.pilot_sets:
                failure_lines.append(f"- {program} | {contrast} | level3: retained {len(selected)}/{args.pilot_sets} from {len(candidates)} candidates.")
            level3_sets = [x[1] for x in selected]
            for set_id, (_, indices, structure) in enumerate(selected, 1):
                null_rows["level3"].append(row_for_set(program, contrast, "level3", set_id, indices, common, studies, structure))
            diagnostic_rows.extend(exchangeability_rows(program, contrast, "level3", target, level3_sets, feature,
                                                         residual, modules, args.modules))

    observed = pd.DataFrame(observed_rows)
    null_tables = {k: pd.DataFrame(v) for k, v in null_rows.items()}
    diagnostics = pd.DataFrame(diagnostic_rows)
    lift = summarize_lift(observed, null_tables)
    # Pilot expansion requires every diagnostic for a program/contrast/level to pass and full unique-set yield.
    pass_table = diagnostics.groupby(["program", "contrast_id", "null_level"], observed=True)["pass"].all().reset_index(name="all_exchangeability_diagnostics_pass")
    counts = pd.concat([v.groupby(["program", "contrast_id", "null_level"], observed=True).size().reset_index(name="n_unique_sets")
                        for v in null_tables.values()], ignore_index=True)
    pass_table = pass_table.merge(counts, on=["program", "contrast_id", "null_level"], how="left")
    pass_table["pilot_set_target"] = args.pilot_sets
    pass_table["approved_for_10000"] = pass_table["all_exchangeability_diagnostics_pass"] & pass_table["n_unique_sets"].eq(args.pilot_sets)

    diagnostics.to_csv(out_dir / "matched_null_exchangeability_diagnostics.tsv", sep="\t", index=False,
                       na_rep="NA", float_format="%.17g", lineterminator="\n")
    observed.to_csv(out_dir / "matched_null_observed_focal.tsv", sep="\t", index=False, na_rep="NA",
                    float_format="%.17g", lineterminator="\n")
    for level, frame in null_tables.items():
        frame.to_csv(out_dir / f"matched_null_{level}.tsv.gz", sep="\t", index=False, compression="gzip",
                     na_rep="NA", float_format="%.17g", lineterminator="\n")
    lift.to_csv(out_dir / "calibrated_program_lift_pilot.tsv", sep="\t", index=False, na_rep="NA",
                float_format="%.17g", lineterminator="\n")
    pass_table.to_csv(out_dir / "matched_null_pilot_expansion_decision.tsv", sep="\t", index=False,
                      na_rep="NA", float_format="%.17g", lineterminator="\n")
    failure_text = "# Matched-null pilot failure and warning log\n\n"
    failure_text += "Outcome-conditioned Level 2 is a conservative sensitivity null, not the primary null.\n\n"
    failure_text += "## Generation failures\n\n" + ("\n".join(failure_lines) if failure_lines else "None. All requested pilot sets were unique and generated.") + "\n"
    failure_text += "\n## Expansion rule\n\nOnly rows with `approved_for_10000 = True` may be expanded without revising the frozen matching design. Failed rows remain not fully calibrated.\n"
    (out_dir / "matched_null_failure_log.md").write_text(failure_text, encoding="utf-8")

    inputs = [gmt_path, *effect_paths]
    for s in STUDIES:
        inputs.extend(studies[s]["input_paths"])
    output_paths = [out_dir / "matched_null_exchangeability_diagnostics.tsv",
                    out_dir / "matched_null_observed_focal.tsv", out_dir / "calibrated_program_lift_pilot.tsv",
                    out_dir / "matched_null_pilot_expansion_decision.tsv", out_dir / "matched_null_failure_log.md",
                    *(out_dir / f"matched_null_{x}.tsv.gz" for x in ("level1", "level2", "level3"))]
    provenance = {
        "analysis": "three-level matched random gene-set pilot",
        "created_at": datetime.now().astimezone().isoformat(), "seed": args.seed,
        "pilot_sets": args.pilot_sets, "level3_candidates": args.level3_candidates,
        "n_common_genes": len(common), "n_coexpression_modules": args.modules,
        "correlation_profile": "within-study group-centered residual logCPM, standardized per gene and concatenated",
        "python": sys.version, "platform": platform.platform(),
        "packages": {"numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__},
        "inputs": [{"path": p.relative_to(root).as_posix(), "bytes": p.stat().st_size, "sha256": sha256(p)} for p in dict.fromkeys(inputs)],
        "outputs": [{"path": p.relative_to(root).as_posix(), "bytes": p.stat().st_size, "sha256": sha256(p)} for p in output_paths],
    }
    (out_dir / "matched_null_pilot_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "ok", "common_genes": len(common), "observed_rows": len(observed),
                      "null_rows": {k: len(v) for k, v in null_tables.items()},
                      "approved_rows": int(pass_table["approved_for_10000"].sum()),
                      "decision_rows": len(pass_table)}))


if __name__ == "__main__":
    main()
