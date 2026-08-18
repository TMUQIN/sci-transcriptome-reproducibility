#!/usr/bin/env python3
"""Aggregate author-state pseudobulk counts before differential expression.

Cross-study harmonisation must occur on counts within each sample, before differential
expression.  Mapping several fitted substate effects to one lineage after modelling would
double-count one study and is prohibited.

Two explicitly separated modes are supported:

* ``lineage`` maps author states to prespecified broad lineages.
* ``whole_lesion`` sums every retained state within a sample.  This is a
  composition-sensitive, bulk-equivalent sensitivity analysis and is never a substitute
  for lineage-resolved inference.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


INVARIANT = ["dataset", "sample", "subject_id", "condition", "group", "dpi", "phase", "region", "segment"]


def map_states(states: pd.Series, mapping: pd.DataFrame) -> tuple[pd.Series, pd.Series]:
    mapping = mapping.sort_values("priority", kind="stable")
    lineages, patterns = [], []
    for state in states.astype(str):
        hit = mapping[mapping.pattern.apply(lambda p: bool(re.search(str(p), state, flags=re.I)))]
        if hit.empty:
            lineages.append(pd.NA); patterns.append(pd.NA)
        else:
            lineages.append(hit.iloc[0].harmonized_lineage); patterns.append(hit.iloc[0].pattern)
    return (pd.Series(lineages, index=states.index, dtype="string"),
            pd.Series(patterns, index=states.index, dtype="string"))


def group_id(dataset: str, sample: str, lineage: str) -> str:
    digest = hashlib.sha1(lineage.encode()).hexdigest()[:8]
    return f"{dataset}::{sample}::{lineage}::{digest}"


def aggregate_whole_lesion(matrix: sp.spmatrix, coldata: pd.DataFrame) -> tuple[sp.csr_matrix, pd.DataFrame]:
    """Sum retained author-state pseudobulks within each biological sample.

    The input has already passed dataset-specific cell QC and the sample x state minimum
    cell gate.  Consequently this branch inherits any cells excluded by that upstream
    gate; the provenance records that limitation rather than reconstructing missing cells.
    """
    if matrix.shape[0] != len(coldata):
        raise ValueError("count rows and coldata rows differ")
    missing = set(INVARIANT + ["cell_state", "n_cells"]) - set(coldata.columns)
    if missing:
        raise ValueError(f"coldata missing {sorted(missing)}")
    codes, levels = pd.factorize(coldata["sample"].astype(str), sort=True)
    indicator = sp.csr_matrix(
        (np.ones(len(codes), dtype=np.int8), (codes, np.arange(len(codes)))),
        shape=(len(levels), len(codes)),
    )
    aggregated = (indicator @ sp.csr_matrix(matrix)).tocsr()
    rows = []
    label = "whole_lesion_bulk_equivalent"
    for code, sample in enumerate(levels):
        block = coldata.iloc[np.flatnonzero(codes == code)]
        values = {}
        for column in INVARIANT:
            unique = block[column].drop_duplicates()
            if len(unique) != 1:
                raise ValueError(f"non-constant {column} within sample {sample}")
            values[column] = unique.iloc[0]
        rows.append(dict(
            group_id=group_id(values["dataset"], sample, label),
            cell_state=label,
            n_cells=int(block.n_cells.sum()),
            source_states=";".join(sorted(block.cell_state.astype(str).unique())),
            **values,
        ))
    out_coldata = pd.DataFrame(rows).set_index("group_id")
    if out_coldata.index.duplicated().any():
        raise ValueError("whole-lesion group_id collision")
    return aggregated, out_coldata


def aggregate_counts(matrix: sp.spmatrix, coldata: pd.DataFrame,
                     mapping: pd.DataFrame) -> tuple[sp.csr_matrix, pd.DataFrame, pd.DataFrame]:
    if matrix.shape[0] != len(coldata):
        raise ValueError("count rows and coldata rows differ")
    missing = set(INVARIANT + ["cell_state", "n_cells"]) - set(coldata.columns)
    if missing:
        raise ValueError(f"coldata missing {sorted(missing)}")
    lineages, patterns = map_states(coldata.cell_state, mapping)
    audit = coldata[["dataset", "cell_state"]].drop_duplicates().copy()
    audit["harmonized_lineage"], audit["mapping_pattern"] = map_states(audit.cell_state, mapping)
    if lineages.isna().any():
        unmapped = sorted(coldata.loc[lineages.isna(), "cell_state"].unique().tolist())
        raise ValueError(f"unmapped cell states; update mapping explicitly: {unmapped}")
    work = coldata.copy(); work["harmonized_lineage"] = lineages.to_numpy()
    keys = work["sample"].astype(str) + "\x1f" + work["harmonized_lineage"].astype(str)
    codes, levels = pd.factorize(keys, sort=True)
    indicator = sp.csr_matrix((np.ones(len(codes), dtype=np.int8), (codes, np.arange(len(codes)))),
                              shape=(len(levels), len(codes)))
    aggregated = (indicator @ sp.csr_matrix(matrix)).tocsr()
    rows = []
    for code, level in enumerate(levels):
        block = work.iloc[np.flatnonzero(codes == code)]
        sample, lineage = level.split("\x1f", 1)
        values = {}
        for column in INVARIANT:
            unique = block[column].drop_duplicates()
            if len(unique) != 1:
                raise ValueError(f"non-constant {column} within {sample}/{lineage}")
            values[column] = unique.iloc[0]
        rows.append(dict(group_id=group_id(values["dataset"], sample, lineage),
                         cell_state=lineage, n_cells=int(block.n_cells.sum()),
                         source_states=";".join(sorted(block.cell_state.astype(str).unique())), **values))
    out_coldata = pd.DataFrame(rows).set_index("group_id")
    if out_coldata.index.duplicated().any():
        raise ValueError("harmonized group_id collision")
    return aggregated, out_coldata, audit.sort_values(["dataset", "cell_state"])


def read_matrix(path: Path) -> sp.csr_matrix:
    with gzip.open(path, "rb") as handle:
        return sp.csr_matrix(scipy.io.mmread(handle))


def write_matrix(path: Path, matrix: sp.spmatrix) -> None:
    with gzip.open(path, "wb") as handle:
        scipy.io.mmwrite(handle, matrix, field="integer")


def process(dataset: str, in_dir: Path, out_dir: Path, mapping_path: Path, mode: str) -> None:
    matrix = read_matrix(in_dir / f"{dataset}_pseudobulk_counts.mtx.gz")
    genes = pd.read_csv(in_dir / f"{dataset}_pseudobulk_genes.tsv", sep="\t")
    coldata = pd.read_csv(in_dir / f"{dataset}_pseudobulk_coldata.tsv", sep="\t", index_col=0)
    if matrix.shape[1] != len(genes):
        raise ValueError("count columns and gene table differ")
    if mode == "lineage":
        mapping = pd.read_csv(mapping_path, sep="\t")
        aggregated, out_coldata, audit = aggregate_counts(matrix, coldata, mapping)
        rule = "sum counts within sample across author states mapped to one broad lineage"
    elif mode == "whole_lesion":
        aggregated, out_coldata = aggregate_whole_lesion(matrix, coldata)
        audit = coldata[["dataset", "cell_state"]].drop_duplicates().sort_values(["dataset", "cell_state"])
        audit["harmonized_lineage"] = "whole_lesion_bulk_equivalent"
        audit["mapping_pattern"] = "__all_retained_states__"
        rule = ("sum counts within sample across all states retained after upstream QC and "
                "sample-by-state minimum-cell gating; composition-sensitive secondary analysis")
    else:
        raise ValueError(f"unknown mode {mode}")
    out_dir.mkdir(parents=True, exist_ok=True)
    write_matrix(out_dir / f"{dataset}_pseudobulk_counts.mtx.gz", aggregated)
    genes.to_csv(out_dir / f"{dataset}_pseudobulk_genes.tsv", sep="\t", index=False)
    out_coldata.to_csv(out_dir / f"{dataset}_pseudobulk_coldata.tsv", sep="\t")
    audit.to_csv(out_dir / f"{dataset}_cell_state_mapping_audit.tsv", sep="\t", index=False)
    upstream_provenance_path = in_dir / f"{dataset}_processing_provenance.json"
    upstream_provenance = (json.loads(upstream_provenance_path.read_text(encoding="utf-8"))
                           if upstream_provenance_path.exists() else {})
    provenance = dict(dataset=dataset, mode=mode, source_groups=int(matrix.shape[0]),
                      harmonized_groups=int(aggregated.shape[0]), genes=int(aggregated.shape[1]),
                      retained_cells_in_aggregated_groups=int(out_coldata.n_cells.sum()),
                      upstream_cells_qc_pass=upstream_provenance.get("cells_qc_pass"),
                      upstream_cells_retained_after_group_gate=upstream_provenance.get("cells_retained_after_group_gate"),
                      rule=rule, mapping=str(mapping_path) if mode == "lineage" else None)
    (out_dir / f"{dataset}_harmonization_provenance.json").write_text(
        json.dumps(provenance, indent=2), encoding="utf-8")
    unit = "sample-lineage" if mode == "lineage" else "whole-lesion sample"
    print(f"[{dataset}] {matrix.shape[0]} author-state groups -> {aggregated.shape[0]} {unit} groups")


def selftest() -> None:
    matrix = sp.csr_matrix([[1, 2], [3, 4], [5, 6]])
    coldata = pd.DataFrame([
        dict(dataset="GSE", sample="S1", subject_id="S1", condition="injured", group="dpi_1",
             dpi=1, phase="acute", region="lesion", segment="T", cell_state="Macrophage", n_cells=10),
        dict(dataset="GSE", sample="S1", subject_id="S1", condition="injured", group="dpi_1",
             dpi=1, phase="acute", region="lesion", segment="T", cell_state="Monocyte", n_cells=20),
        dict(dataset="GSE", sample="S2", subject_id="S2", condition="uninjured", group="uninjured",
             dpi=0, phase="uninjured", region="lesion", segment="T", cell_state="Microglia", n_cells=30),
    ])
    mapping = pd.DataFrame([
        dict(priority=1, pattern="macrophage|monocyte", harmonized_lineage="macrophage_monocyte"),
        dict(priority=2, pattern="microglia", harmonized_lineage="microglia")])
    out, meta, _ = aggregate_counts(matrix, coldata, mapping)
    assert out.shape == (2, 2)
    row = meta[meta.cell_state.eq("macrophage_monocyte")].index[0]
    pos = meta.index.get_loc(row)
    assert np.array_equal(out[pos].toarray().ravel(), [4, 6])
    assert meta.loc[row, "n_cells"] == 30
    whole, whole_meta = aggregate_whole_lesion(matrix, coldata)
    assert whole.shape == (2, 2)
    s1 = whole_meta.index[whole_meta["sample"].eq("S1")][0]
    s1_pos = whole_meta.index.get_loc(s1)
    assert np.array_equal(whole[s1_pos].toarray().ravel(), [4, 6])
    assert whole_meta.loc[s1, "cell_state"] == "whole_lesion_bulk_equivalent"
    assert whole_meta.loc[s1, "n_cells"] == 30
    print("[selftest] PASS: substates summed within lineage or whole sample before modelling")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", nargs="+")
    ap.add_argument("--in-dir", default="data_processed")
    ap.add_argument("--out-dir", default="data_processed/harmonized")
    ap.add_argument("--mapping", default="tables/cell_state_harmonization.tsv")
    ap.add_argument("--mode", choices=["lineage", "whole_lesion"], default="lineage")
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest(); return
    if not args.dataset:
        ap.error("provide --dataset GSE...")
    for dataset in args.dataset:
        process(dataset, Path(args.in_dir), Path(args.out_dir), Path(args.mapping), args.mode)


if __name__ == "__main__":
    main()
