#!/usr/bin/env python3
"""Compute magnitude-aware and rank-based Hallmark scores per biological pseudobulk."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import platform
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.io import mmread
from scipy.stats import rankdata


STUDIES = ("GSE162610", "GSE234774", "GSE304399")
GROUPS = ("uninjured", "dpi_1", "dpi_7")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def parse_gmt(path: Path) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                sets[fields[0]] = set(fields[2:])
    if len(sets) != 50:
        raise ValueError(f"expected 50 Hallmarks, observed {len(sets)}")
    return sets


def load_dataset(root: Path, study: str) -> tuple[np.ndarray, pd.Series, pd.DataFrame, list[Path]]:
    base = root / "data_processed" / "whole_lesion"
    matrix_path = base / f"{study}_pseudobulk_counts.mtx.gz"
    genes_path = base / f"{study}_pseudobulk_genes.tsv"
    coldata_path = base / f"{study}_pseudobulk_coldata.tsv"
    with gzip.open(matrix_path, "rb") as handle:
        matrix = mmread(handle)
    matrix = matrix.toarray() if hasattr(matrix, "toarray") else np.asarray(matrix)
    genes = pd.read_csv(genes_path, sep="\t")["gene"].astype(str)
    coldata = pd.read_csv(coldata_path, sep="\t")
    if matrix.shape != (len(coldata), len(genes)):
        raise ValueError(f"{study}: matrix {matrix.shape} does not match samples/genes {len(coldata), len(genes)}")
    if genes.duplicated().any():
        raise ValueError(f"{study}: duplicate gene symbols")
    if coldata["subject_id"].duplicated().any():
        raise ValueError(f"{study}: duplicated biological subject in whole-lesion pseudobulk")
    keep = coldata["group"].isin(GROUPS).to_numpy()
    matrix = matrix[keep, :].astype(float, copy=False)
    coldata = coldata.loc[keep].reset_index(drop=True)
    observed = coldata.groupby("group")["subject_id"].nunique().to_dict()
    if any(observed.get(group, 0) < 2 for group in GROUPS):
        raise ValueError(f"{study}: fewer than two biological replicates in required group: {observed}")
    return matrix, genes, coldata, [matrix_path, genes_path, coldata_path]


def score_dataset(counts: np.ndarray, genes: pd.Series, coldata: pd.DataFrame,
                  gene_sets: dict[str, set[str]], study: str) -> pd.DataFrame:
    libsize = counts.sum(axis=1)
    if np.any(libsize <= 0):
        raise ValueError(f"{study}: zero pseudobulk library size")
    logcpm = np.log2(counts / libsize[:, None] * 1_000_000 + 0.5)
    gene_sd = logcpm.std(axis=0, ddof=1)
    variable = np.isfinite(gene_sd) & (gene_sd > 0)
    z = np.zeros_like(logcpm)
    z[:, variable] = (logcpm[:, variable] - logcpm[:, variable].mean(axis=0)) / gene_sd[variable]
    ranks = np.apply_along_axis(rankdata, 1, logcpm, method="average")
    rank_score_all = 2 * (ranks - (logcpm.shape[1] + 1) / 2) / max(logcpm.shape[1] - 1, 1)
    gene_to_index = {g: i for i, g in enumerate(genes)}
    rows = []
    for program, members in gene_sets.items():
        index = np.array([gene_to_index[g] for g in sorted(members) if g in gene_to_index], dtype=int)
        index_variable = index[variable[index]]
        if len(index_variable) < 10:
            continue
        method_scores = {
            "magnitude_mean_z": z[:, index_variable].mean(axis=1),
            "centered_rank_mean": rank_score_all[:, index].mean(axis=1),
        }
        for method, raw_score in method_scores.items():
            sd = float(np.std(raw_score, ddof=1))
            if not np.isfinite(sd) or sd <= 0:
                continue
            standardized = (raw_score - float(np.mean(raw_score))) / sd
            for i, meta in coldata.iterrows():
                rows.append({
                    "dataset": study,
                    "sample": meta["sample"],
                    "subject_id": meta["subject_id"],
                    "group": meta["group"],
                    "program": program,
                    "method": method,
                    "raw_score": raw_score[i],
                    "standardized_score": standardized[i],
                    "n_genes_in_set_reference": len(members),
                    "n_genes_observed": len(index),
                    "n_variable_genes_scored": len(index_variable),
                    "library_size": libsize[i],
                })
    return pd.DataFrame(rows)


def selftest() -> None:
    rng = np.random.default_rng(7)
    counts = rng.poisson(lam=np.arange(2, 12), size=(6, 10)).astype(float) + 1
    counts[2:4, :4] += 10
    counts[4:6, 4:8] += 10
    coldata = pd.DataFrame({"sample": list("abcdef"), "subject_id": list("abcdef"),
                            "group": ["uninjured"] * 2 + ["dpi_1"] * 2 + ["dpi_7"] * 2})
    sets = {f"SET_{i}": {"A", "B", "C", "D", "E", "F", "G", "H", "I", "J"} for i in range(50)}
    genes = pd.Series([chr(65 + i) for i in range(10)])
    result = score_dataset(counts, genes, coldata, sets, "TEST")
    assert result["subject_id"].nunique() == 6 and result["method"].nunique() == 2
    assert np.allclose(result.groupby(["program", "method"])["standardized_score"].std().to_numpy(), 1.0)
    print("[selftest] PASS: biological-sample scores and within-study standardization")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--out-dir", type=Path, default=Path("reports/phase_reproducibility_calibration_2026_07"))
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()
    if args.selftest:
        selftest()
        return
    root = args.root.resolve()
    out_dir = args.out_dir if args.out_dir.is_absolute() else root / args.out_dir
    score_dir = out_dir / "sample_level_hallmark_scores"
    score_dir.mkdir(parents=True, exist_ok=True)
    gmt_path = root / "references" / "msigdb_mh.all.v2026.1.Mm.symbols.gmt"
    sets = parse_gmt(gmt_path)
    all_scores = []
    input_paths = [gmt_path]
    for study in STUDIES:
        counts, genes, coldata, paths = load_dataset(root, study)
        scores = score_dataset(counts, genes, coldata, sets, study)
        out_path = score_dir / f"{study}_hallmark_scores.tsv"
        scores.to_csv(out_path, sep="\t", index=False, na_rep="NA", float_format="%.17g", lineterminator="\n")
        all_scores.append(scores)
        input_paths.extend(paths)
    combined = pd.concat(all_scores, ignore_index=True)
    combined_path = out_dir / "sample_level_hallmark_scores.tsv"
    combined.to_csv(combined_path, sep="\t", index=False, na_rep="NA", float_format="%.17g", lineterminator="\n")
    outputs = [combined_path, *(score_dir / f"{s}_hallmark_scores.tsv" for s in STUDIES)]
    provenance = {
        "analysis": "sample-level Hallmark scoring",
        "created_at": datetime.now().astimezone().isoformat(),
        "methods": {
            "magnitude_mean_z": "mean of within-study gene-standardized log2(CPM+0.5)",
            "centered_rank_mean": "centered mean within-sample rank; rank sensitivity, not the official singscore package",
            "model_input": "each program score standardized to SD=1 across uninjured/d1/d7 samples within study",
        },
        "biological_unit": "subject_id/library pseudobulk; no cell or technical-partition inference",
        "python": sys.version,
        "platform": platform.platform(),
        "packages": {"numpy": np.__version__, "pandas": pd.__version__, "scipy": scipy.__version__},
        "inputs": [{"path": p.relative_to(root).as_posix(), "bytes": p.stat().st_size, "sha256": sha256(p)} for p in input_paths],
        "outputs": [{"path": p.relative_to(root).as_posix(), "bytes": p.stat().st_size, "sha256": sha256(p)} for p in outputs],
    }
    (out_dir / "sample_level_hallmark_score_provenance.json").write_text(
        json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": "ok", "rows": len(combined), "samples": combined["subject_id"].nunique(),
                      "programs": combined["program"].nunique(), "methods": combined["method"].nunique()}))


if __name__ == "__main__":
    main()
