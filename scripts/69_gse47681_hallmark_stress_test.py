#!/usr/bin/env python3
"""Hallmark GSEA and frozen-reference comparison for GSE47681.

This is a post-lock, orthogonal cross-platform context evaluation.  It never writes
to the frozen three-study results directories and never pools GSE47681 into
the discovery meta-analysis.
"""
from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import platform
import sys
from pathlib import Path

import gseapy
import numpy as np
import pandas as pd
import scipy
from scipy.stats import spearmanr


STUDIES = ("GSE162610", "GSE234774", "GSE304399")
CONTRASTS = (
    "injury_d1_vs_uninjured",
    "injury_d7_vs_uninjured",
    "change_d7_minus_d1",
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


def load_gsea_helpers(root: Path):
    script = root / "scripts" / "06e_hallmark_gsea.py"
    spec = importlib.util.spec_from_file_location("frozen_hallmark_gsea", script)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import helpers from {script}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def binary_metrics(reference: np.ndarray, external: np.ndarray) -> dict[str, float | int]:
    reference = np.asarray(reference, float)
    external = np.asarray(external, float)
    ref_sign = np.sign(reference).astype(np.int8)
    ext_sign = np.sign(external).astype(np.int8)
    valid = (ref_sign != 0) & (ext_sign != 0)
    ref = ref_sign[valid]
    ext = ext_sign[valid]
    tp = int(np.sum((ref == 1) & (ext == 1)))
    tn = int(np.sum((ref == -1) & (ext == -1)))
    fp = int(np.sum((ref == -1) & (ext == 1)))
    fn = int(np.sum((ref == 1) & (ext == -1)))
    accuracy = float((tp + tn) / len(ref)) if len(ref) else math.nan
    tpr = tp / (tp + fn) if tp + fn else math.nan
    tnr = tn / (tn + fp) if tn + fp else math.nan
    balanced = float((tpr + tnr) / 2) if np.isfinite(tpr) and np.isfinite(tnr) else math.nan
    denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = float((tp * tn - fp * fn) / denom) if denom else 0.0
    return {
        "n_binary": int(len(ref)),
        "n_zero_excluded": int(len(reference) - len(ref)),
        "raw_direction_accuracy": accuracy,
        "balanced_direction_accuracy": balanced,
        "mcc": mcc,
        "tp": tp,
        "tn": tn,
        "fp": fp,
        "fn": fn,
    }


def comparison_metrics(matrix: pd.DataFrame, external_col: str) -> dict[str, float | int]:
    frozen = matrix.loc[:, list(STUDIES)].to_numpy(float)
    external = matrix[external_col].to_numpy(float)
    reference = frozen.mean(axis=1)
    binary = binary_metrics(reference, external)
    unanimous = np.all(np.sign(frozen) == np.sign(frozen[:, [0]]), axis=1) & (np.sign(frozen[:, 0]) != 0)
    ext_nonzero = np.sign(external) != 0
    gated_valid = unanimous & ext_nonzero
    conditional = float(np.mean(np.sign(external[gated_valid]) == np.sign(frozen[gated_valid, 0]))) if np.any(gated_valid) else math.nan
    rho = float(spearmanr(external, reference).statistic)
    return {
        "n_shared_features": int(len(matrix)),
        "spearman_vs_frozen_mean": rho,
        **binary,
        "consensus_coverage": float(np.mean(unanimous)),
        "n_unanimous_frozen": int(np.sum(unanimous)),
        "n_unanimous_with_nonzero_external": int(np.sum(gated_valid)),
        "conditional_accuracy_given_unanimous_frozen": conditional,
        "frozen_mean_positive_fraction": float(np.mean(reference > 0)),
        "external_positive_fraction": float(np.mean(external > 0)),
    }


def build_gene_matrix(root: Path, external_effects: Path, contrast: str) -> pd.DataFrame:
    columns = []
    for study in STUDIES:
        path = root / "results" / "whole_lesion" / f"effects_{study}.tsv"
        frame = pd.read_csv(path, sep="\t", usecols=["gene", "contrast_id", "estimate", "se"])
        frame = frame[frame["contrast_id"].eq(contrast)].copy()
        frame[study] = frame["estimate"] / frame["se"]
        if frame["gene"].duplicated().any():
            raise ValueError(f"Duplicate genes in {path} for {contrast}")
        columns.append(frame.set_index("gene")[[study]])
    ext = pd.read_csv(external_effects, sep="\t", usecols=["gene", "contrast_id", "estimate", "se"])
    ext = ext[ext["contrast_id"].eq(contrast)].copy()
    ext["GSE47681"] = ext["estimate"] / ext["se"]
    if ext["gene"].duplicated().any():
        raise ValueError(f"Duplicate genes in {external_effects} for {contrast}")
    columns.append(ext.set_index("gene")[["GSE47681"]])
    matrix = pd.concat(columns, axis=1, join="inner").dropna()
    if not np.isfinite(matrix.to_numpy(float)).all():
        raise ValueError("Non-finite gene score in shared matrix")
    return matrix.sort_index()


def build_hallmark_matrix(root: Path, external_gsea: pd.DataFrame, contrast: str) -> pd.DataFrame:
    frozen_path = root / "results" / "whole_lesion_programs" / "hallmark_gsea_by_study.tsv"
    frozen = pd.read_csv(frozen_path, sep="\t", usecols=["dataset", "contrast_id", "term", "nes"])
    frozen = frozen[frozen["contrast_id"].eq(contrast) & frozen["dataset"].isin(STUDIES)]
    matrix = frozen.pivot(index="term", columns="dataset", values="nes")
    ext = external_gsea[external_gsea["contrast_id"].eq(contrast)].set_index("term")["nes"].rename("GSE47681")
    matrix = matrix.join(ext, how="inner").dropna()
    matrix = matrix.loc[:, [*STUDIES, "GSE47681"]].sort_index()
    if matrix.shape[0] != 50:
        raise ValueError(f"Expected 50 shared Hallmarks for {contrast}; found {matrix.shape[0]}")
    if not np.isfinite(matrix.to_numpy(float)).all():
        raise ValueError("Non-finite Hallmark score in shared matrix")
    return matrix


def run_external_gsea(root: Path, effects_path: Path, out_dir: Path) -> pd.DataFrame:
    helper = load_gsea_helpers(root)
    gmt = root / "references" / "msigdb_mh.all.v2026.1.Mm.symbols.gmt"
    effects = pd.read_csv(effects_path, sep="\t")
    outputs = []
    for (dataset, contrast, effect_type), group in effects.groupby(
        ["dataset", "contrast_id", "effect_type"], sort=True, observed=True
    ):
        if contrast not in CONTRASTS:
            continue
        group = group.copy()
        group["rank_statistic"] = group["estimate"] / group["se"]
        ranked = (
            group[["gene", "rank_statistic"]]
            .replace([np.inf, -np.inf], np.nan)
            .dropna()
            .drop_duplicates("gene", keep=False)
            .sort_values(["rank_statistic", "gene"], ascending=[False, True], kind="stable")
        )
        ranked, n_tied = helper.deterministic_break_ties(ranked)
        pre = gseapy.prerank(
            rnk=ranked,
            gene_sets=str(gmt),
            min_size=15,
            max_size=500,
            permutation_num=5000,
            weight=1.0,
            ascending=False,
            threads=1,
            seed=20260803,
            outdir=None,
            no_plot=True,
            verbose=False,
        )
        outputs.append(
            helper.normalise_result(
                pre.res2d, str(dataset), str(contrast), str(effect_type), len(ranked), n_tied
            )
        )
    result = pd.concat(outputs, ignore_index=True)
    if len(result) != 150 or result.groupby("contrast_id")["term"].nunique().ne(50).any():
        raise ValueError(f"Unexpected external GSEA dimensions: {result.shape}")
    write_tsv(result, out_dir / "GSE47681_hallmark_gsea.tsv")
    provenance = {
        "analysis_role": "orthogonal_cross_platform_context_evaluation",
        "dataset": "GSE47681",
        "not_pooled_with_frozen_meta_analysis": True,
        "rank_statistic": "moderated t statistic (estimate/se) from raw-CEL RMA plus robust limma",
        "gene_sets": str(gmt.relative_to(root)),
        "gene_sets_bytes": gmt.stat().st_size,
        "gene_sets_sha256": sha256(gmt),
        "permutations": 5000,
        "seed": 20260803,
        "min_size": 15,
        "max_size": 500,
        "weight": 1.0,
        "gseapy_version": gseapy.__version__,
        "exact_tie_handling": "same frozen deterministic IEEE-754 nextafter rule as script 06e",
        "input_effects": {"path": str(effects_path.relative_to(root)), "bytes": effects_path.stat().st_size, "sha256": sha256(effects_path)},
    }
    (out_dir / "GSE47681_hallmark_gsea_provenance.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8"
    )
    return result


def main() -> None:
    root = Path(__file__).resolve().parents[1]
    out_dir = root / "reports" / "phase_v5_3_3_postlock_upgrade_2026_08" / "GSE47681"
    out_dir.mkdir(parents=True, exist_ok=True)
    effects_path = out_dir / "effects_GSE47681.tsv"
    if not effects_path.exists():
        raise FileNotFoundError(effects_path)

    external_gsea = run_external_gsea(root, effects_path, out_dir)
    metric_rows = []
    pair_rows = []
    shared_rows = []
    for scale in ("gene", "hallmark"):
        for contrast in CONTRASTS:
            matrix = (
                build_gene_matrix(root, effects_path, contrast)
                if scale == "gene"
                else build_hallmark_matrix(root, external_gsea, contrast)
            )
            metrics = comparison_metrics(matrix, "GSE47681")
            metric_rows.append({"scale": scale, "contrast_id": contrast, **metrics})
            for study in STUDIES:
                pair_rows.append(
                    {
                        "scale": scale,
                        "contrast_id": contrast,
                        "reference": study,
                        "n_shared_features": len(matrix),
                        "spearman_rho": float(spearmanr(matrix[study], matrix["GSE47681"]).statistic),
                        "direction_accuracy": float(np.mean(np.sign(matrix[study]) == np.sign(matrix["GSE47681"]))),
                    }
                )
            export = matrix.reset_index().rename(columns={matrix.index.name or "index": "feature"})
            export.insert(0, "contrast_id", contrast)
            export.insert(0, "scale", scale)
            export["frozen_mean_score"] = export.loc[:, list(STUDIES)].mean(axis=1)
            signs = np.sign(export.loc[:, list(STUDIES)].to_numpy(float))
            export["frozen_unanimous_direction"] = np.where(
                np.all(signs == signs[:, [0]], axis=1) & (signs[:, 0] != 0), signs[:, 0], 0
            ).astype(int)
            shared_rows.append(export)

    metrics = pd.DataFrame(metric_rows)
    pairs = pd.DataFrame(pair_rows)
    shared = pd.concat(shared_rows, ignore_index=True)
    write_tsv(metrics, out_dir / "GSE47681_frozen_reference_comparison.tsv")
    write_tsv(pairs, out_dir / "GSE47681_pairwise_reference_diagnostics.tsv")
    shared.to_csv(
        out_dir / "GSE47681_shared_feature_score_matrix.tsv.gz",
        sep="\t", index=False, na_rep="NA", float_format="%.17g", compression="gzip", lineterminator="\n"
    )

    provenance = {
        "analysis_role": "orthogonal_cross_platform_context_evaluation",
        "metric_contract": {
            "reference_score": "arithmetic mean of the three frozen signed scores",
            "raw_binary_metrics": "sign(GSE47681) versus sign(reference mean)",
            "consensus_gate": "all three frozen studies share the same non-zero sign",
            "conditional_accuracy": "GSE47681 sign accuracy among unanimity-gated features",
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
            "gseapy": gseapy.__version__,
        },
        "inputs": [
            {"path": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in [
                effects_path,
                root / "results" / "whole_lesion_programs" / "hallmark_gsea_by_study.tsv",
                *[root / "results" / "whole_lesion" / f"effects_{study}.tsv" for study in STUDIES],
            ]
        ],
    }
    (out_dir / "GSE47681_stress_test_provenance.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8"
    )
    print(metrics.to_string(index=False))


if __name__ == "__main__":
    main()
