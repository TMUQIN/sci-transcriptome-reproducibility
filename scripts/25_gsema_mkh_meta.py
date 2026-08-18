#!/usr/bin/env python3
"""Apply the project's REML+mKH estimator to GSEMA study-level effects."""
from __future__ import annotations

import hashlib
import importlib.util
import json
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "reports" / "phase_reproducibility_calibration_2026_07" / "gsema_benchmark"


def load_module():
    path = ROOT / "scripts" / "06_meta_analysis.py"
    spec = importlib.util.spec_from_file_location("sci_meta", path)
    module = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    h = hashlib.sha256(path.read_bytes()).hexdigest()
    return h


def main() -> None:
    source = OUT / "gsema_effects_by_study.tsv"
    effects = pd.read_csv(source, sep="\t")
    effects = effects[np.isfinite(effects["estimate"]) & np.isfinite(effects["se"]) & effects["se"].gt(0)].copy()
    mm = load_module()
    rows = []
    for (method, contrast, pathway), group in effects.groupby(["method", "contrast_id", "pathway"], observed=True):
        group = group.sort_values("dataset")
        if group["dataset"].nunique() != 3:
            continue
        result = mm.meta_reml_hksj(group["estimate"].to_numpy(float), group["se"].to_numpy(float))
        rows.append({"method": method, "contrast_id": contrast, "pathway": pathway, **result,
                     "all_study_effects_same_sign": bool(np.all(np.sign(group["estimate"]) == np.sign(group["estimate"].iloc[0]))),
                     "study_effects": ";".join(f"{d}:{e:.12g}" for d, e in zip(group["dataset"], group["estimate"]))})
    meta = pd.DataFrame(rows)
    meta["meta_fdr"] = np.nan
    for _, idx in meta.groupby(["method", "contrast_id"], observed=True).groups.items():
        meta.loc[idx, "meta_fdr"] = mm.bh_adjust(meta.loc[idx, "p"].to_numpy(float))
    meta = meta.sort_values(["method", "contrast_id", "meta_fdr", "pathway"], kind="stable")
    target = OUT / "gsema_mkh_random_effects.tsv"
    meta.to_csv(target, sep="\t", index=False, na_rep="NA", float_format="%.17g", lineterminator="\n")
    status = OUT / "gsema_run_status.tsv"
    native = OUT / "gsema_native_random_effects.tsv"
    provenance = {
        "analysis": "GSEMA benchmark with project REML+mKH uncertainty",
        "created_at": datetime.now().astimezone().isoformat(),
        "inputs": [{"path": str(p.relative_to(ROOT)).replace("\\", "/"), "bytes": p.stat().st_size, "sha256": sha256(p)}
                   for p in (source, status, native, ROOT / "scripts" / "06_meta_analysis.py")],
        "outputs": [{"path": str(target.relative_to(ROOT)).replace("\\", "/"), "bytes": target.stat().st_size, "sha256": sha256(target)}],
        "boundary": "GSEMA native Z inference is retained separately; manuscript inference uses mKH because k=3.",
    }
    (OUT / "gsema_benchmark_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "ok", "rows": len(meta), "methods": meta["method"].nunique()}))


if __name__ == "__main__":
    main()
