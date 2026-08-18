#!/usr/bin/env python3
"""REML + modified Hartung-Knapp meta-analysis of sample-level pathway effects."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


def load_meta_module(root: Path):
    path = root / "scripts" / "06_meta_analysis.py"
    spec = importlib.util.spec_from_file_location("sci_meta", path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--out-dir", type=Path, default=Path("reports/phase_reproducibility_calibration_2026_07"))
    args = parser.parse_args()
    root = args.root.resolve()
    out_dir = args.out_dir if args.out_dir.is_absolute() else root / args.out_dir
    input_path = out_dir / "pathway_effects_by_study.tsv"
    effects = pd.read_csv(input_path, sep="\t")
    meta_mod = load_meta_module(root)
    rows = []
    for (method, contrast, effect_type, program), group in effects.groupby(
        ["method", "contrast_id", "effect_type", "program"], observed=True, sort=False
    ):
        group = group.sort_values("dataset")
        result = meta_mod.meta_reml_hksj(group["estimate"].to_numpy(float), group["se"].to_numpy(float))
        rows.append({"method": method, "contrast_id": contrast, "effect_type": effect_type, "program": program,
                     **result, "n_study_fdr_le_005": int((group["fdr"] <= 0.05).sum()),
                     "all_study_effects_same_sign": bool(np.all(np.sign(group["estimate"]) == np.sign(group["estimate"].iloc[0]))),
                     "datasets": ";".join(group["dataset"]),
                     "study_effects": ";".join(f"{d}:{e:.12g}" for d, e in zip(group["dataset"], group["estimate"])),
                     "study_fdr": ";".join(f"{d}:{q:.12g}" for d, q in zip(group["dataset"], group["fdr"]))})
    meta = pd.DataFrame(rows)
    meta["meta_fdr"] = np.nan
    for _, idx in meta.groupby(["method", "contrast_id"], observed=True).groups.items():
        meta.loc[idx, "meta_fdr"] = meta_mod.bh_adjust(meta.loc[idx, "p"].to_numpy(float))
    meta = meta.sort_values(["method", "contrast_id", "meta_fdr", "program"], kind="stable")
    meta_path = out_dir / "pathway_meta_mkh.tsv"
    meta.to_csv(meta_path, sep="\t", index=False, na_rep="NA", float_format="%.17g", lineterminator="\n")

    gsea_path = root / "results" / "whole_lesion_programs" / "hallmark_gsea_by_study.tsv"
    gsea = pd.read_csv(gsea_path, sep="\t").rename(columns={"term": "program", "nes": "gsea_nes"})
    merged = effects.merge(gsea[["dataset", "contrast_id", "program", "gsea_nes", "fdr"]].rename(columns={"fdr": "gsea_fdr"}),
                           on=["dataset", "contrast_id", "program"], how="left", validate="many_to_one")
    merged["score_gsea_same_sign"] = np.sign(merged["estimate"]) == np.sign(merged["gsea_nes"])
    concordance = merged.groupby(["method", "contrast_id"], observed=True).agg(
        n_rows=("program", "size"), direction_concordance=("score_gsea_same_sign", "mean")
    ).reset_index()
    pair = effects.pivot_table(index=["dataset", "contrast_id", "program"], columns="method", values="estimate").dropna().reset_index()
    pair["magnitude_rank_same_sign"] = np.sign(pair["magnitude_mean_z"]) == np.sign(pair["centered_rank_mean"])
    method_agreement = pair.groupby("contrast_id", observed=True).agg(
        n_study_program_rows=("program", "size"), magnitude_rank_direction_concordance=("magnitude_rank_same_sign", "mean")
    ).reset_index()
    concordance = concordance.merge(method_agreement, on="contrast_id", how="left")
    concordance_path = out_dir / "gsea_vs_score_meta_concordance.tsv"
    concordance.to_csv(concordance_path, sep="\t", index=False, na_rep="NA", float_format="%.17g", lineterminator="\n")

    provenance = {
        "analysis": "sample-level Hallmark effect meta-analysis",
        "created_at": datetime.now().astimezone().isoformat(),
        "estimator": "REML tau2 plus modified Hartung-Knapp (q>=1), reused from scripts/06_meta_analysis.py",
        "multiplicity": "BH within scoring method and contrast across 50 Hallmarks",
        "inputs": [{"path": p.relative_to(root).as_posix(), "bytes": p.stat().st_size, "sha256": sha256(p)}
                   for p in (input_path, gsea_path, root / "scripts" / "06_meta_analysis.py")],
        "outputs": [{"path": p.relative_to(root).as_posix(), "bytes": p.stat().st_size, "sha256": sha256(p)}
                    for p in (meta_path, concordance_path)],
        "warnings": ["k=3 heterogeneity statistics and prediction intervals are descriptive.",
                     "Pathway scores are transcriptomic summaries, not direct pathway activation measurements."],
    }
    (out_dir / "pathway_effect_meta_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "ok", "study_effect_rows": len(effects), "meta_rows": len(meta)}))


if __name__ == "__main__":
    main()
