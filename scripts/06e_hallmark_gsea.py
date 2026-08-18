#!/usr/bin/env python3
"""Secondary Hallmark GSEA for the composition-sensitive whole-lesion branch.

Each dataset/contrast is ranked independently by moderated t statistic (estimate / SE).
No DEG threshold is applied.  Cross-study output is descriptive: NES signs and within-study
GSEA FDR are compared without treating NES as an effect with a known standard error.
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
from itertools import combinations
from pathlib import Path

import gseapy
import numpy as np
import pandas as pd


USECOLS = ["gene", "dataset", "contrast_id", "effect_type", "estimate", "se"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalise_result(result: pd.DataFrame, dataset: str, contrast: str,
                     effect_type: str, n_ranked: int, n_tied_genes: int) -> pd.DataFrame:
    aliases = {
        "Term": "term", "ES": "es", "NES": "nes", "NOM p-val": "nominal_p",
        "FDR q-val": "fdr", "FWER p-val": "fwer_p", "Tag %": "tag_fraction",
        "Gene %": "gene_fraction", "Lead_genes": "leading_edge_genes",
    }
    out = result.rename(columns=aliases).copy()
    required = {"term", "es", "nes", "nominal_p", "fdr"}
    missing = required - set(out.columns)
    if missing:
        raise ValueError(f"GSEA result missing {sorted(missing)}; columns={list(result.columns)}")
    keep = [column for column in aliases.values() if column in out.columns]
    out = out[keep]
    out.insert(0, "effect_type", effect_type)
    out.insert(0, "contrast_id", contrast)
    out.insert(0, "dataset", dataset)
    out["n_ranked_genes"] = n_ranked
    out["n_exact_tied_genes_before_deterministic_break"] = n_tied_genes
    for column in ("es", "nes", "nominal_p", "fdr", "fwer_p"):
        if column in out:
            out[column] = pd.to_numeric(out[column], errors="coerce")
    return out


def deterministic_break_ties(ranked: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """Resolve exact statistic ties by gene order using adjacent floats.

    Only exact ties are changed.  The first alphabetically ordered gene retains the
    original value and subsequent genes receive successive ``nextafter(..., -inf)``
    values, so downstream sorting is deterministic without adding a material effect.
    """
    ranked = ranked.sort_values(["rank_statistic", "gene"], ascending=[False, True], kind="stable").copy()
    tied = ranked.rank_statistic.duplicated(keep=False)
    n_tied = int(tied.sum())
    if not n_tied:
        return ranked, 0
    for _, index in ranked.loc[tied].groupby("rank_statistic", sort=False).groups.items():
        previous = float(ranked.loc[index[0], "rank_statistic"])
        for row_index in index[1:]:
            previous = float(np.nextafter(previous, -np.inf))
            ranked.loc[row_index, "rank_statistic"] = previous
    if ranked.rank_statistic.duplicated().any():
        raise AssertionError("deterministic tie break left exact duplicate statistics")
    return ranked, n_tied


def summarize_replication(gsea: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (contrast, effect_type, term), group in gsea.groupby(
        ["contrast_id", "effect_type", "term"], sort=True, observed=True
    ):
        group = group.sort_values("dataset")
        significant = group[group.fdr.le(.05)]
        significant_signs = set(np.sign(significant.nes[significant.nes.ne(0)]))
        all_signs = set(np.sign(group.nes[group.nes.ne(0)]))
        rows.append({
            "contrast_id": contrast,
            "effect_type": effect_type,
            "term": term,
            "n_studies": len(group),
            "n_study_fdr_005": int(group.fdr.le(.05).sum()),
            "all_nes_same_sign": len(all_signs) == 1 and len(group) > 0,
            "significant_nes_same_sign": len(significant_signs) <= 1,
            "replicated_fdr_2plus_same_sign": len(significant) >= 2 and len(significant_signs) == 1,
            "replicated_fdr_all_studies_same_sign": len(significant) == len(group) and len(group) >= 2 and len(significant_signs) == 1,
            "datasets": ";".join(group.dataset.astype(str)),
            "nes_by_dataset": ";".join(f"{row.dataset}:{row.nes:.8g}" for row in group.itertuples()),
            "fdr_by_dataset": ";".join(f"{row.dataset}:{row.fdr:.8g}" for row in group.itertuples()),
            "median_nes": float(group.nes.median()),
            "min_fdr": float(group.fdr.min()),
            "max_fdr": float(group.fdr.max()),
        })
    return pd.DataFrame(rows)


def nes_agreement(gsea: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for contrast, group in gsea.groupby("contrast_id", sort=True, observed=True):
        pivot = group.pivot(index="term", columns="dataset", values="nes")
        for left, right in combinations(sorted(pivot.columns), 2):
            pair = pivot[[left, right]].dropna()
            rows.append({
                "contrast_id": contrast,
                "dataset_1": left,
                "dataset_2": right,
                "n_shared_programs": len(pair),
                "pearson_r": pair[left].corr(pair[right], method="pearson"),
                "spearman_rho": pair[left].corr(pair[right], method="spearman"),
                "sign_agreement_fraction": float((np.sign(pair[left]) == np.sign(pair[right])).mean()),
            })
    return pd.DataFrame(rows)


def temporal_program_patterns(gsea: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    wanted = {
        "early": "injury_d1_vs_uninjured",
        "late": "injury_d7_vs_uninjured",
        "delta": "change_d7_minus_d1",
    }
    components = []
    for label, contrast in wanted.items():
        part = gsea[gsea.contrast_id.eq(contrast)][["dataset", "term", "nes", "fdr"]].copy()
        part = part.rename(columns={"nes": f"{label}_nes", "fdr": f"{label}_fdr"})
        components.append(part)
    temporal = components[0]
    for component in components[1:]:
        temporal = temporal.merge(component, on=["dataset", "term"], how="inner", validate="one_to_one")
    temporal["all_three_fdr_005"] = temporal[
        ["early_fdr", "late_fdr", "delta_fdr"]
    ].le(.05).all(axis=1)
    e, l, d = np.sign(temporal.early_nes), np.sign(temporal.late_nes), np.sign(temporal.delta_nes)
    conditions = [
        (e > 0) & (l > 0) & (d < 0),
        (e > 0) & (l > 0) & (d > 0),
        (e < 0) & (l < 0) & (d > 0),
        (e < 0) & (l < 0) & (d < 0),
        (e > 0) & (l < 0) & (d < 0),
        (e < 0) & (l > 0) & (d > 0),
    ]
    labels = [
        "activated_but_attenuating", "activated_and_strengthening",
        "suppressed_but_recovering", "suppressed_and_deepening",
        "positive_to_negative_reversal", "negative_to_positive_reversal",
    ]
    temporal["direction_pattern"] = np.select(conditions, labels, default="direction_inconsistent")
    temporal["fdr_supported_pattern"] = np.where(
        temporal.all_three_fdr_005, temporal.direction_pattern, "incomplete_fdr_support"
    )
    rows = []
    for term, group in temporal.groupby("term", sort=True, observed=True):
        supported = group[group.all_three_fdr_005]
        counts = supported.direction_pattern.value_counts()
        rows.append({
            "term": term,
            "n_studies_complete": len(group),
            "n_studies_all_three_fdr_005": len(supported),
            "n_distinct_supported_patterns": int(supported.direction_pattern.nunique()),
            "replicated_same_pattern_2plus": bool(len(counts) and counts.iloc[0] >= 2),
            "replicated_same_pattern_all_studies": bool(len(counts) == 1 and counts.iloc[0] == len(group)),
            "dominant_supported_pattern": counts.index[0] if len(counts) else "none",
            "datasets_and_patterns": ";".join(
                f"{row.dataset}:{row.fdr_supported_pattern}" for row in group.sort_values("dataset").itertuples()
            ),
        })
    return temporal, pd.DataFrame(rows)


def validate_resumable_result(result: pd.DataFrame, effect_paths: list[Path],
                              contrasts: list[str]) -> None:
    required = {"dataset", "contrast_id", "effect_type", "term", "nes", "fdr"}
    missing = required - set(result.columns)
    if missing:
        raise ValueError(f"resume table missing {sorted(missing)}")
    expected_combinations = set()
    for path in effect_paths:
        keys = pd.read_csv(path, sep="\t", usecols=["dataset", "contrast_id"]).drop_duplicates()
        keys = keys[keys.contrast_id.isin(contrasts)]
        expected_combinations.update(map(tuple, keys[["dataset", "contrast_id"]].itertuples(index=False, name=None)))
    actual_combinations = set(map(tuple, result[["dataset", "contrast_id"]].drop_duplicates()
                                  .itertuples(index=False, name=None)))
    if actual_combinations != expected_combinations:
        raise ValueError(f"resume combinations differ: expected={sorted(expected_combinations)}, "
                         f"actual={sorted(actual_combinations)}")
    term_sets = [frozenset(group.term.astype(str)) for _, group in result.groupby(
        ["dataset", "contrast_id"], observed=True
    )]
    if not term_sets or len(term_sets[0]) < 1 or any(terms != term_sets[0] for terms in term_sets[1:]):
        raise ValueError("resume table does not contain one identical non-empty gene-set collection per combination")
    if result.duplicated(["dataset", "contrast_id", "term"]).any():
        raise ValueError("resume table contains duplicate dataset/contrast/term rows")


def selftest() -> None:
    toy = pd.DataFrame([
        {"dataset": d, "contrast_id": "C", "effect_type": "endpoint", "term": "A", "nes": n, "fdr": q}
        for d, n, q in (("D1", 2.0, .01), ("D2", 1.5, .04), ("D3", 1.0, .2))
    ])
    out = summarize_replication(toy).iloc[0]
    assert out.n_studies == 3 and out.n_study_fdr_005 == 2
    assert bool(out.replicated_fdr_2plus_same_sign)
    assert not bool(out.replicated_fdr_all_studies_same_sign)
    ranked, n_tied = deterministic_break_ties(pd.DataFrame({
        "gene": ["B", "A", "C"], "rank_statistic": [1.0, 1.0, 0.5]
    }))
    assert n_tied == 2 and ranked.iloc[0].gene == "A"
    assert ranked.rank_statistic.is_unique and ranked.rank_statistic.is_monotonic_decreasing
    toy_temporal = pd.DataFrame([
        {"dataset": dataset, "term": "A", "effect_type": "endpoint" if "change" not in contrast else "temporal_delta",
         "contrast_id": contrast, "nes": nes, "fdr": .01}
        for dataset in ("D1", "D2")
        for contrast, nes in (("injury_d1_vs_uninjured", 2.0),
                              ("injury_d7_vs_uninjured", 1.0),
                              ("change_d7_minus_d1", -1.5))
    ])
    _, temporal_summary = temporal_program_patterns(toy_temporal)
    assert bool(temporal_summary.iloc[0].replicated_same_pattern_all_studies)
    assert temporal_summary.iloc[0].dominant_supported_pattern == "activated_but_attenuating"
    print("[selftest] PASS: program replication and deterministic exact-tie handling")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--effects", nargs="+", default=["results/whole_lesion/effects_GSE*.tsv"])
    parser.add_argument("--gmt", default="references/msigdb_mh.all.v2026.1.Mm.symbols.gmt")
    parser.add_argument("--contrasts", nargs="+", default=[
        "injury_d1_vs_uninjured", "injury_d7_vs_uninjured", "change_d7_minus_d1"
    ])
    parser.add_argument("--out-dir", default="results/whole_lesion_programs")
    parser.add_argument("--permutations", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20260728)
    parser.add_argument("--min-size", type=int, default=15)
    parser.add_argument("--max-size", type=int, default=500)
    parser.add_argument("--resume", action="store_true",
                        help="reuse a complete by-study table after strict combination/gene-set validation")
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()
    if args.selftest:
        selftest()
        return
    if args.permutations < 1000:
        raise ValueError("program-level analysis requires at least 1000 permutations")
    paths = [Path(path) for path in sorted(set(path for pattern in args.effects for path in glob.glob(pattern)))]
    if not paths:
        raise FileNotFoundError(f"no effect files matched {args.effects}")
    gmt = Path(args.gmt)
    if not gmt.exists():
        raise FileNotFoundError(gmt)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    by_study_path = out_dir / "hallmark_gsea_by_study.tsv"
    if args.resume:
        if not by_study_path.exists():
            raise FileNotFoundError(f"--resume requested but missing {by_study_path}")
        result = pd.read_csv(by_study_path, sep="\t")
        validate_resumable_result(result, paths, args.contrasts)
    else:
        outputs, failures = [], []
        for path in paths:
            effects = pd.read_csv(path, sep="\t", usecols=USECOLS)
            effects = effects[effects.contrast_id.isin(args.contrasts)].copy()
            for (dataset, contrast, effect_type), group in effects.groupby(
                ["dataset", "contrast_id", "effect_type"], sort=True, observed=True
            ):
                group["rank_statistic"] = group.estimate / group.se
                ranked = (group[["gene", "rank_statistic"]]
                          .replace([np.inf, -np.inf], np.nan).dropna()
                          .drop_duplicates("gene", keep=False)
                          .sort_values(["rank_statistic", "gene"], ascending=[False, True], kind="stable"))
                ranked, n_tied = deterministic_break_ties(ranked)
                try:
                    pre = gseapy.prerank(
                        rnk=ranked, gene_sets=str(gmt), min_size=args.min_size,
                        max_size=args.max_size, permutation_num=args.permutations,
                        weight=1.0, ascending=False, threads=1, seed=args.seed,
                        outdir=None, no_plot=True, verbose=False,
                    )
                    outputs.append(normalise_result(
                        pre.res2d, str(dataset), str(contrast), str(effect_type), len(ranked), n_tied
                    ))
                except Exception as exc:
                    failures.append({"source_file": str(path), "dataset": dataset,
                                     "contrast_id": contrast, "error": repr(exc)})
        if failures:
            raise RuntimeError(f"GSEA failures: {failures}")
        result = pd.concat(outputs, ignore_index=True)
        validate_resumable_result(result, paths, args.contrasts)
        result.to_csv(by_study_path, sep="\t", index=False)
    replication = summarize_replication(result)
    replication.to_csv(out_dir / "hallmark_cross_study_replication.tsv", sep="\t", index=False)
    nes_agreement(result).to_csv(out_dir / "hallmark_nes_agreement.tsv", sep="\t", index=False)
    temporal, temporal_summary = temporal_program_patterns(result)
    temporal.to_csv(out_dir / "hallmark_temporal_patterns_by_study.tsv", sep="\t", index=False)
    temporal_summary.to_csv(out_dir / "hallmark_temporal_pattern_replication.tsv", sep="\t", index=False)
    provenance = {
        "analysis_role": "secondary_exploratory_program_level_whole_lesion",
        "rank_statistic": "estimate/se from sample-level limma-voom contrast",
        "deg_threshold_used": False,
        "gene_sets": str(gmt), "gene_sets_bytes": gmt.stat().st_size,
        "gene_sets_sha256": sha256(gmt), "gseapy_version": gseapy.__version__,
        "permutations": args.permutations, "seed": args.seed,
        "zero_empirical_p_or_fdr_note": ("reported zero means no exceedance at the configured finite "
                                         "permutation count; it is not a mathematical probability of zero"),
        "min_size": args.min_size, "max_size": args.max_size, "weight": 1.0,
        "exact_tie_handling": ("rank statistic descending, gene ascending; subsequent exact ties "
                               "shifted by successive IEEE-754 nextafter toward -infinity"),
        "resumed_from_validated_complete_by_study_table": bool(args.resume),
        "input_effects": [{"path": str(path), "bytes": path.stat().st_size,
                            "sha256": sha256(path)} for path in paths],
        "interpretation_limit": ("NES is compared descriptively across independent studies; "
                                 "it is not meta-analysed as an effect with a known SE"),
    }
    (out_dir / "hallmark_gsea_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"Wrote {len(result)} study-program rows and {len(replication)} cross-study summaries")


if __name__ == "__main__":
    main()
