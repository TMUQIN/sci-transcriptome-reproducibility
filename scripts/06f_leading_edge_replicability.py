#!/usr/bin/env python3
"""Audit cross-study leading-edge reproducibility of post-primary Hallmark hits.

This is an exploratory robustness analysis.  It asks whether programs that showed the
same whole-lesion temporal pattern in all three cohorts are driven by shared genes, and
whether leading genes have broad directionality across author-labelled lineages in the
two cohorts where those labels exist.  It does not remove composition confounding.
"""
from __future__ import annotations

import argparse
import glob
import hashlib
import json
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom


CONTRASTS = {
    "early": "injury_d1_vs_uninjured",
    "late": "injury_d7_vs_uninjured",
    "delta": "change_d7_minus_d1",
}
EXPECTED_SIGN = {
    "injury_d1_vs_uninjured": 1,
    "injury_d7_vs_uninjured": 1,
    "change_d7_minus_d1": -1,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_bool(series: pd.Series) -> pd.Series:
    values = series.astype(str).str.strip().str.lower()
    allowed = {"true": True, "false": False, "1": True, "0": False, "yes": True, "no": False}
    bad = sorted(set(values) - set(allowed))
    if bad:
        raise ValueError(f"invalid boolean values: {bad}")
    return values.map(allowed).astype(bool)


def read_gmt(path: Path) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                sets[fields[0]] = set(fields[2:])
    return sets


def bh_adjust(values: pd.Series) -> pd.Series:
    p = values.to_numpy(float)
    order = np.argsort(p)
    ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    output = np.empty_like(ranked)
    output[order] = np.clip(ranked, 0, 1)
    return pd.Series(output, index=values.index)


def leading_set(value: object) -> set[str]:
    if pd.isna(value) or not str(value).strip():
        return set()
    return {gene for gene in str(value).split(";") if gene}


def build_overlap(gsea: pd.DataFrame, terms: list[str], gene_sets: dict[str, set[str]]) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    pair_rows, consensus_rows = [], []
    cache: dict[tuple[str, str, str], set[str]] = {}
    for term in terms:
        if term not in gene_sets:
            raise ValueError(f"selected term absent from GMT: {term}")
        for contrast in CONTRASTS.values():
            block = gsea[gsea.term.eq(term) & gsea.contrast_id.eq(contrast)].sort_values("dataset")
            if len(block) != 3 or block.dataset.nunique() != 3:
                raise ValueError(f"expected three study rows for {term}/{contrast}, found {len(block)}")
            study_sets = {}
            for row in block.itertuples():
                genes = leading_set(row.leading_edge_genes)
                if not genes <= gene_sets[term]:
                    raise ValueError(f"leading edge contains genes outside GMT for {term}/{row.dataset}")
                cache[(term, contrast, row.dataset)] = genes
                study_sets[row.dataset] = genes
            for left, right in combinations(sorted(study_sets), 2):
                a, b = study_sets[left], study_sets[right]
                intersection, union = a & b, a | b
                pair_rows.append({
                    "term": term, "contrast_id": contrast,
                    "dataset_1": left, "dataset_2": right,
                    "gene_set_size": len(gene_sets[term]),
                    "n_leading_1": len(a), "n_leading_2": len(b),
                    "n_intersection": len(intersection), "n_union": len(union),
                    "jaccard": len(intersection) / len(union) if union else np.nan,
                    "overlap_coefficient": len(intersection) / min(len(a), len(b)) if a and b else np.nan,
                    "hypergeom_p_within_gene_set": hypergeom.sf(
                        len(intersection) - 1, len(gene_sets[term]), len(a), len(b)
                    ),
                    "intersection_genes": ";".join(sorted(intersection)),
                })
            sets = list(study_sets.values())
            all_three = set.intersection(*sets)
            counts: dict[str, int] = {}
            for genes in sets:
                for gene in genes:
                    counts[gene] = counts.get(gene, 0) + 1
            at_least_two = {gene for gene, n in counts.items() if n >= 2}
            union = set.union(*sets)
            consensus_rows.append({
                "term": term, "contrast_id": contrast, "gene_set_size": len(gene_sets[term]),
                "n_union_leading_edge": len(union), "n_at_least_two_studies": len(at_least_two),
                "n_all_three_studies": len(all_three),
                "fraction_gene_set_at_least_two": len(at_least_two) / len(gene_sets[term]),
                "fraction_gene_set_all_three": len(all_three) / len(gene_sets[term]),
                "genes_at_least_two_studies": ";".join(sorted(at_least_two)),
                "genes_all_three_studies": ";".join(sorted(all_three)),
            })
    pair = pd.DataFrame(pair_rows)
    pair["overlap_fdr_across_selected_tests"] = bh_adjust(pair.hypergeom_p_within_gene_set)
    return pair, pd.DataFrame(consensus_rows), cache


def load_whole_effects(patterns: list[str], genes: set[str]) -> pd.DataFrame:
    paths = [Path(path) for pattern in patterns for path in glob.glob(pattern)]
    if not paths:
        raise FileNotFoundError("no whole-lesion effect files")
    frames = []
    usecols = ["gene", "dataset", "contrast_id", "estimate", "se", "fdr"]
    for path in sorted(set(paths)):
        for chunk in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=100_000):
            chunk = chunk[chunk.gene.isin(genes) & chunk.contrast_id.isin(CONTRASTS.values())]
            if len(chunk):
                frames.append(chunk)
    effects = pd.concat(frames, ignore_index=True)
    if effects.duplicated(["gene", "dataset", "contrast_id"]).any():
        raise ValueError("duplicate whole-lesion gene/dataset/contrast effects")
    effects["rank_statistic"] = effects.estimate / effects.se
    return effects


def gene_core_table(effects: pd.DataFrame, terms: list[str], gene_sets: dict[str, set[str]], cache: dict) -> pd.DataFrame:
    datasets = sorted(effects.dataset.unique())
    if len(datasets) != 3:
        raise ValueError(f"expected three whole-lesion datasets, found {datasets}")
    lookup = effects.set_index(["gene", "dataset", "contrast_id"])
    rows = []
    for term in terms:
        for gene in sorted(gene_sets[term]):
            row: dict[str, object] = {"term": term, "gene": gene}
            study_patterns = []
            for dataset in datasets:
                signs = []
                for label, contrast in CONTRASTS.items():
                    key = (gene, dataset, contrast)
                    if key in lookup.index:
                        effect = lookup.loc[key]
                        statistic = float(effect.rank_statistic)
                        row[f"{dataset}_{label}_rank_statistic"] = statistic
                        row[f"{dataset}_{label}_fdr"] = float(effect.fdr)
                        signs.append(np.sign(statistic) == EXPECTED_SIGN[contrast])
                    else:
                        row[f"{dataset}_{label}_rank_statistic"] = np.nan
                        row[f"{dataset}_{label}_fdr"] = np.nan
                        signs.append(False)
                    row[f"{dataset}_{label}_leading_edge"] = gene in cache[(term, contrast, dataset)]
                pattern = all(signs)
                row[f"{dataset}_expected_temporal_sign_pattern"] = pattern
                study_patterns.append(pattern)
            row["n_studies_expected_temporal_sign_pattern"] = int(sum(study_patterns))
            for label, contrast in CONTRASTS.items():
                row[f"n_studies_{label}_leading_edge"] = int(sum(
                    gene in cache[(term, contrast, dataset)] for dataset in datasets
                ))
            row["strong_shared_temporal_driver"] = (
                row["n_studies_expected_temporal_sign_pattern"] == 3 and
                all(row[f"n_studies_{label}_leading_edge"] >= 2 for label in CONTRASTS)
            )
            row["partial_shared_temporal_driver"] = (
                row["n_studies_expected_temporal_sign_pattern"] >= 2 and
                all(row[f"n_studies_{label}_leading_edge"] >= 2 for label in CONTRASTS)
            )
            rows.append(row)
    return pd.DataFrame(rows)


def program_redundancy(terms: list[str], gene_sets: dict[str, set[str]], gene_core: pd.DataFrame) -> pd.DataFrame:
    rows = []
    strong = {term: set(group.loc[group.strong_shared_temporal_driver, "gene"])
              for term, group in gene_core.groupby("term", observed=True)}
    partial = {term: set(group.loc[group.partial_shared_temporal_driver, "gene"])
               for term, group in gene_core.groupby("term", observed=True)}
    for left, right in combinations(sorted(terms), 2):
        a, b = gene_sets[left], gene_sets[right]
        s1, s2 = strong[left], strong[right]
        p1, p2 = partial[left], partial[right]
        rows.append({
            "term_1": left, "term_2": right,
            "gmt_intersection": len(a & b), "gmt_union": len(a | b),
            "gmt_jaccard": len(a & b) / len(a | b),
            "strong_driver_intersection": len(s1 & s2), "strong_driver_union": len(s1 | s2),
            "strong_driver_jaccard": len(s1 & s2) / len(s1 | s2) if s1 | s2 else np.nan,
            "partial_driver_intersection": len(p1 & p2), "partial_driver_union": len(p1 | p2),
            "partial_driver_jaccard": len(p1 & p2) / len(p1 | p2) if p1 | p2 else np.nan,
            "shared_strong_driver_genes": ";".join(sorted(s1 & s2)),
        })
    return pd.DataFrame(rows)


def load_lineage_effects(patterns: list[str], genes: set[str]) -> pd.DataFrame:
    paths = [Path(path) for pattern in patterns for path in glob.glob(pattern)]
    if not paths:
        raise FileNotFoundError("no lineage effect files")
    frames = []
    usecols = ["gene", "cell_state", "dataset", "contrast_id", "estimate", "se", "fdr"]
    for path in sorted(set(paths)):
        for chunk in pd.read_csv(path, sep="\t", usecols=usecols, chunksize=150_000):
            chunk = chunk[chunk.gene.isin(genes) & chunk.contrast_id.isin(CONTRASTS.values())]
            if len(chunk):
                frames.append(chunk)
    effects = pd.concat(frames, ignore_index=True)
    if effects.duplicated(["gene", "cell_state", "dataset", "contrast_id"]).any():
        raise ValueError("duplicate lineage effect key")
    return effects


def lineage_breadth(lineage: pd.DataFrame, terms: list[str], cache: dict) -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    for (gene, dataset, contrast), group in lineage.groupby(
        ["gene", "dataset", "contrast_id"], sort=False, observed=True
    ):
        expected = EXPECTED_SIGN[contrast]
        signs = np.sign(group.estimate.to_numpy(float))
        for term in terms:
            if gene not in cache.get((term, contrast, dataset), set()):
                continue
            same = signs == expected
            rows.append({
                "term": term, "gene": gene, "dataset": dataset, "contrast_id": contrast,
                "expected_program_sign": expected, "n_lineages_tested": len(group),
                "n_lineages_expected_sign": int(same.sum()),
                "fraction_lineages_expected_sign": float(same.mean()),
                "n_lineage_fdr005_expected_sign": int(((group.fdr <= .05).to_numpy() & same).sum()),
                "n_lineage_fdr005_opposite_sign": int(((group.fdr <= .05).to_numpy() & ~same).sum()),
                "lineages_expected_sign": ";".join(sorted(group.loc[same, "cell_state"].astype(str))),
                "lineages_opposite_sign": ";".join(sorted(group.loc[~same, "cell_state"].astype(str))),
                "broad_direction_075": bool(len(group) >= 3 and same.mean() >= .75),
            })
    gene = pd.DataFrame(rows)
    summaries = []
    for key, group in gene.groupby(["term", "dataset", "contrast_id"], observed=True):
        summaries.append({
            "term": key[0], "dataset": key[1], "contrast_id": key[2],
            "n_leading_edge_genes_with_lineage_data": len(group),
            "median_n_lineages_tested": float(group.n_lineages_tested.median()),
            "median_fraction_lineages_expected_sign": float(group.fraction_lineages_expected_sign.median()),
            "fraction_genes_broad_direction_075": float(group.broad_direction_075.mean()),
            "fraction_genes_any_lineage_fdr005_expected": float(group.n_lineage_fdr005_expected_sign.gt(0).mean()),
            "fraction_genes_any_lineage_fdr005_opposite": float(group.n_lineage_fdr005_opposite_sign.gt(0).mean()),
        })
    return gene, pd.DataFrame(summaries)


def audit_summary(terms: list[str], pair: pd.DataFrame, consensus: pd.DataFrame,
                  gene_core: pd.DataFrame, breadth_summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for term in terms:
        core = gene_core[gene_core.term.eq(term)]
        row = {
            "term": term,
            "n_strong_shared_temporal_drivers": int(core.strong_shared_temporal_driver.sum()),
            "n_partial_shared_temporal_drivers": int(core.partial_shared_temporal_driver.sum()),
            "strong_driver_genes": ";".join(sorted(core.loc[core.strong_shared_temporal_driver, "gene"])),
            "partial_driver_genes": ";".join(sorted(core.loc[core.partial_shared_temporal_driver, "gene"])),
        }
        for label, contrast in CONTRASTS.items():
            c = consensus[(consensus.term.eq(term)) & (consensus.contrast_id.eq(contrast))].iloc[0]
            p = pair[(pair.term.eq(term)) & (pair.contrast_id.eq(contrast))]
            row[f"{label}_leading_edge_all3_genes"] = int(c.n_all_three_studies)
            row[f"{label}_leading_edge_atleast2_genes"] = int(c.n_at_least_two_studies)
            row[f"{label}_median_pairwise_jaccard"] = float(p.jaccard.median())
        subset = breadth_summary[breadth_summary.term.eq(term)]
        row["median_lineage_broad_fraction_across_available_study_contrasts"] = (
            float(subset.fraction_genes_broad_direction_075.median()) if len(subset) else np.nan
        )
        rows.append(row)
    return pd.DataFrame(rows)


def selftest() -> None:
    assert leading_set("A;B;A") == {"A", "B"}
    adjusted = bh_adjust(pd.Series([.01, .04, .2]))
    assert np.allclose(adjusted, [.03, .06, .2])
    assert parse_bool(pd.Series(["True", "no", "1"])).tolist() == [True, False, True]
    print("[selftest] PASS: leading-edge parsing, BH and strict booleans")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gsea", default="results/whole_lesion_programs/hallmark_gsea_by_study.tsv")
    parser.add_argument("--temporal-summary", default="results/whole_lesion_programs/hallmark_temporal_pattern_replication.tsv")
    parser.add_argument("--gmt", default="references/msigdb_mh.all.v2026.1.Mm.symbols.gmt")
    parser.add_argument("--whole-effects", nargs="+", default=["results/whole_lesion/effects_GSE*.tsv"])
    parser.add_argument("--lineage-effects", nargs="+", default=[
        "results/harmonized/effects_GSE162610.tsv", "results/harmonized/effects_GSE234774.tsv"
    ])
    parser.add_argument("--out-dir", default="results/whole_lesion_programs/leading_edge_audit")
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()
    if args.selftest:
        selftest(); return
    gsea_path, temporal_path, gmt_path = Path(args.gsea), Path(args.temporal_summary), Path(args.gmt)
    gsea = pd.read_csv(gsea_path, sep="\t")
    temporal = pd.read_csv(temporal_path, sep="\t")
    selected = parse_bool(temporal.replicated_same_pattern_all_studies) & temporal.dominant_supported_pattern.eq(
        "activated_but_attenuating"
    )
    terms = sorted(temporal.loc[selected, "term"].astype(str))
    if not terms:
        raise ValueError("no programs meet the declared all-study temporal selection rule")
    gene_sets = read_gmt(gmt_path)
    pair, consensus, cache = build_overlap(gsea, terms, gene_sets)
    union_genes = set().union(*(gene_sets[term] for term in terms))
    whole = load_whole_effects(args.whole_effects, union_genes)
    gene_core = gene_core_table(whole, terms, gene_sets, cache)
    redundancy = program_redundancy(terms, gene_sets, gene_core)
    lineage = load_lineage_effects(args.lineage_effects, union_genes)
    breadth_gene, breadth_summary = lineage_breadth(lineage, terms, cache)
    summary = audit_summary(terms, pair, consensus, gene_core, breadth_summary)
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    pair.to_csv(out / "leading_edge_pairwise_overlap.tsv", sep="\t", index=False)
    consensus.to_csv(out / "leading_edge_consensus.tsv", sep="\t", index=False)
    gene_core.to_csv(out / "temporal_gene_core.tsv", sep="\t", index=False)
    redundancy.to_csv(out / "program_redundancy.tsv", sep="\t", index=False)
    breadth_gene.to_csv(out / "lineage_breadth_by_gene.tsv", sep="\t", index=False)
    breadth_summary.to_csv(out / "lineage_breadth_summary.tsv", sep="\t", index=False)
    summary.to_csv(out / "leading_edge_audit_summary.tsv", sep="\t", index=False)
    input_paths = [gsea_path, temporal_path, gmt_path]
    input_paths += [Path(path) for pattern in args.whole_effects for path in glob.glob(pattern)]
    input_paths += [Path(path) for pattern in args.lineage_effects for path in glob.glob(pattern)]
    provenance = {
        "analysis_role": "post-primary_exploratory_leading_edge_replicability_audit",
        "program_selection_rule": "all three studies have the same FDR-supported activated_but_attenuating temporal pattern",
        "selected_terms": terms,
        "pairwise_overlap_null": "random leading-edge subsets drawn from the corresponding fixed Hallmark gene set",
        "strong_shared_temporal_driver_rule": "expected +,+,- gene rank signs in all 3 studies and leading-edge membership in >=2 studies for each of d1,d7,delta",
        "partial_shared_temporal_driver_rule": "expected +,+,- gene rank signs in >=2 studies and leading-edge membership in >=2 studies for each contrast",
        "lineage_breadth_rule": ">=3 tested author-labelled broad lineages and >=75% with the expected whole-lesion program direction",
        "interpretation_limits": [
            "programs and genes were selected after observing the primary and GSEA results",
            "Hallmark gene sets overlap and are not independent mechanisms",
            "lineage breadth is available only for GSE162610 and GSE234774",
            "same-direction lineage effects do not eliminate cell-composition confounding",
            "GSE304399 still lacks barcode-level lineage labels",
        ],
        "inputs": [{"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
                   for path in sorted(set(input_paths))],
    }
    (out / "leading_edge_audit_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"Wrote leading-edge audit for {len(terms)} programs; "
          f"{int(gene_core.strong_shared_temporal_driver.sum())} term-gene strong drivers")


if __name__ == "__main__":
    main()
