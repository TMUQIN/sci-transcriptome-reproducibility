#!/usr/bin/env python3
"""Dataset-aware QC, canonical sample crosswalk and sparse pseudobulk.

This module deliberately separates three identities:

* cell barcode (``adata.obs_names``),
* GEO library (canonical ``sample``/GSM), and
* biological/paired-modality unit (``subject_id``).

Inferential mode fails if any cell cannot be mapped to a real GEO library.  There is no
dummy ``S1`` fallback.  Published annotations are preserved as ``cell_state_original``;
an unsupervised cluster is used when no author annotation exists, without forcing the
cluster to the highest-scoring marker label.

Primary pseudobulk output is sparse Matrix Market (groups x genes):

``<GSE>_pseudobulk_counts.mtx.gz``, ``<GSE>_pseudobulk_genes.tsv`` and
``<GSE>_pseudobulk_coldata.tsv``.
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


DATASETS = {
    "GSE234774": dict(
        kind="mtx", matrix="GSE234774_rnaseq_filtered_scRNA.mtx.gz",
        features="GSE234774_rnaseq_features.txt.gz", barcodes="GSE234774_rnaseq_barcodes.txt.gz",
        meta="GSE234774_rnaseq_meta.txt.gz", meta_join="barcode", barcode_col="barcode",
        library_col="library", author_annotation="layer3", modality="snRNA-seq",
        default_mode="inferential",
        stream_mtx=True, stream_qc_mode="filtered_metrics",
    ),
    "GSE230765": dict(
        kind="mtx", matrix="GSE230765_scRNA_filtered_scRNA.mtx.gz",
        features="GSE230765_scRNA_features.txt.gz", barcodes="GSE230765_scRNA_barcodes.txt.gz",
        meta="GSE230765_scRNA_meta.txt.gz", meta_join="barcode", barcode_col="barcode",
        library_col="replicate", author_annotation="layer3", modality="snRNA-seq",
        default_mode="support",
        stream_mtx=True, stream_qc_mode="filtered_no_metrics",
    ),
    "GSE162610": dict(
        kind="mtx", matrix="GSE162610_sci_mat.mtx.gz", features="GSE162610_genes.tsv.gz",
        barcodes="GSE162610_barcodes.tsv.gz", meta="GSE162610_barcode_metadata.tsv.gz",
        meta_join="index", library_col="orig.ident", author_annotation="celltype",
        modality="scRNA-seq", default_mode="inferential",
        stream_mtx=True, stream_qc_mode="author_flags",
    ),
    "GSE172167": dict(
        kind="mtx", matrix="GSE172167_spinalcord_SCI_counts_matrix.mtx.gz",
        features="GSE172167_features.txt.gz", barcodes="GSE172167_barcodes.txt.gz",
        feature_header=0, barcode_header=0, meta="GSE172167_meta.txt.gz", meta_join="position",
        author_annotation="subclass", modality="snRNA-seq", default_mode="inferential",
        stream_mtx=True, stream_qc_mode="filtered_no_metrics",
        stream_format="dense_whitespace_header",
    ),
    "GSE192824": dict(
        kind="csv", matrix="GSE192824_All_cells_UMI_Matrix.csv.gz",
        author_annotation=None, modality="scRNA-seq", default_mode="descriptive",
    ),
    "GSE304399": dict(
        kind="mtx", matrix="GSE304399_rna_matrix.mtx.gz",
        features="GSE304399_rna_features.tsv.gz", barcodes="GSE304399_rna_barcodes.tsv.gz",
        feature_header=0, feature_sep=r"\s+",
        meta="GSE304399_rna_barcode_meta.tsv", meta_join="position", library_col="library_key",
        author_annotation="cell_state", modality="multiome-RNA", default_mode="inferential",
        stream_mtx=True, stream_qc_mode="filtered_no_metrics", stream_format="matrix_market_coordinate",
    ),
}

CANONICAL_COLUMNS = [
    "sample", "subject_id", "subject_verification", "condition", "group", "dpi",
    "phase", "region", "segment", "modality", "primary_arm", "inference_role",
    "exclusion_reason", "verification_status",
]


def read_one_column(path: Path, header=None, sep="\t") -> np.ndarray:
    frame = pd.read_csv(path, sep=sep, header=header, dtype=str)
    return frame.iloc[:, -1].astype(str).to_numpy()


def read_mtx_dataset(cfg: dict, raw_dir: Path):
    import anndata as ad

    features = read_one_column(raw_dir / cfg["features"], header=cfg.get("feature_header"),
                               sep=cfg.get("feature_sep", "\t"))
    barcodes = read_one_column(raw_dir / cfg["barcodes"], header=cfg.get("barcode_header"),
                               sep=cfg.get("barcode_sep", "\t"))
    with gzip.open(raw_dir / cfg["matrix"], "rb") as handle:
        matrix = scipy.io.mmread(handle)
    matrix = sp.csr_matrix(matrix)
    if matrix.shape == (len(features), len(barcodes)):
        matrix = matrix.T.tocsr()
    elif matrix.shape != (len(barcodes), len(features)):
        raise ValueError(
            f"matrix {matrix.shape} does not match {len(barcodes)} barcodes x {len(features)} features"
        )
    adata = ad.AnnData(X=matrix)
    adata.obs_names = pd.Index(barcodes, dtype=str)
    adata.var_names = pd.Index(features, dtype=str)
    adata.var_names_make_unique()
    return adata


def read_csv_dataset(cfg: dict, raw_dir: Path):
    import anndata as ad

    frame = pd.read_csv(raw_dir / cfg["matrix"], index_col=0)
    if frame.index.duplicated().any():
        frame = frame.groupby(level=0, sort=False).sum()
    matrix = sp.csr_matrix(frame.to_numpy().T)
    adata = ad.AnnData(X=matrix)
    adata.obs_names = pd.Index(frame.columns.astype(str))
    adata.var_names = pd.Index(frame.index.astype(str))
    adata.var_names_make_unique()
    return adata


def attach_author_metadata(adata, dataset: str, cfg: dict, raw_dir: Path) -> None:
    if dataset == "GSE192824":
        adata.obs["library_key"] = np.where(
            adata.obs_names.str.startswith("Troy_UI_"), "Troy_UI",
            np.where(adata.obs_names.str.startswith("Troy_SCI_"), "Troy_SCI", pd.NA),
        )
        return
    meta_path = raw_dir / cfg["meta"]
    if not meta_path.exists():
        raise FileNotFoundError(f"author cell metadata is required: {meta_path}")
    join = cfg["meta_join"]
    if join == "index":
        meta = pd.read_csv(meta_path, sep="\t", index_col=0, low_memory=False)
        meta.index = meta.index.astype(str)
        if not adata.obs_names.isin(meta.index).all():
            missing = adata.obs_names[~adata.obs_names.isin(meta.index)][:5].tolist()
            raise ValueError(f"cell metadata missing {len(missing)}+ barcode(s), e.g. {missing}")
        meta = meta.reindex(adata.obs_names)
    elif join == "barcode":
        meta = pd.read_csv(meta_path, sep="\t", low_memory=False)
        key = cfg["barcode_col"]
        if meta[key].duplicated().any():
            raise ValueError(f"duplicate author metadata barcode in {key}")
        meta = meta.set_index(meta[key].astype(str)).reindex(adata.obs_names)
        if meta.isna().all(axis=1).any():
            missing = meta.index[meta.isna().all(axis=1)][:5].tolist()
            raise ValueError(f"cell metadata join failed, e.g. {missing}")
    elif join == "position":
        meta = pd.read_csv(meta_path, sep="\t", low_memory=False)
        if len(meta) != adata.n_obs:
            raise ValueError(f"positional metadata has {len(meta)} rows; expected {adata.n_obs}")
        meta.index = adata.obs_names
    else:
        raise ValueError(f"unsupported meta_join={join}")
    for column in meta.columns:
        name = column if column not in adata.obs.columns else f"author_{column}"
        adata.obs[name] = meta[column].to_numpy()

    if dataset == "GSE172167":
        adata.obs["library_key"] = adata.obs["orig.ident"].astype(str)
    else:
        adata.obs["library_key"] = adata.obs[cfg["library_col"]].astype(str)


def canonical_sample_lookup(metadata: pd.DataFrame, dataset: str, modality: str) -> pd.DataFrame:
    lookup = metadata[metadata["dataset"].eq(dataset)].copy()
    if "modality" in lookup and lookup["modality"].eq(modality).any():
        lookup = lookup[lookup["modality"].eq(modality)].copy()
    if lookup["library_key"].duplicated().any():
        dup = lookup.loc[lookup["library_key"].duplicated(False), "library_key"].unique().tolist()
        raise ValueError(f"non-unique library_key in canonical metadata: {dup[:5]}")
    return lookup.set_index("library_key", drop=False)


def attach_canonical_metadata(adata, dataset: str, cfg: dict, sample_metadata: Path) -> None:
    metadata = pd.read_csv(sample_metadata, sep="\t", dtype={"sample": str, "library_key": str})
    lookup = canonical_sample_lookup(metadata, dataset, cfg["modality"])
    keys = adata.obs["library_key"].astype("string")
    missing = sorted(set(keys.dropna()) - set(lookup.index.astype(str)))
    if missing or keys.isna().any():
        raise ValueError(
            f"{dataset}: cell-to-library crosswalk incomplete; missing keys={missing[:10]}, "
            f"null_cells={int(keys.isna().sum())}"
        )
    for column in CANONICAL_COLUMNS:
        adata.obs[column] = keys.map(lookup[column]).to_numpy()
    adata.obs["dataset"] = dataset
    validate_obs_contract(adata.obs, inferential=False)


def validate_obs_contract(obs: pd.DataFrame, inferential: bool) -> None:
    required = ["dataset", "sample", "subject_id", "condition", "group", "library_key"]
    missing_cols = [c for c in required if c not in obs]
    if missing_cols:
        raise ValueError(f"missing canonical obs columns: {missing_cols}")
    nulls = {c: int(obs[c].isna().sum()) for c in required if obs[c].isna().any()}
    if nulls:
        raise ValueError(f"null canonical sample identities: {nulls}")
    if inferential:
        if (obs["sample"].astype(str).str.upper() == "S1").any():
            raise ValueError("dummy sample S1 is prohibited in inferential mode")
        bad = ~obs["condition"].isin(["injured", "uninjured"])
        if bad.any():
            raise ValueError(f"unknown condition in {int(bad.sum())} inferential cells")
        roles = set(obs["inference_role"].astype(str)) if "inference_role" in obs else set()
        if "within_study_inference" not in roles:
            raise ValueError(f"dataset is not eligible for inferential mode; roles={sorted(roles)}")


def analysis_scope_mask(obs: pd.DataFrame, mode: str) -> pd.Series:
    if mode == "inferential":
        return obs["inference_role"].astype(str).eq("within_study_inference")
    if mode == "support":
        return obs["inference_role"].astype(str).isin(["regulatory_support", "spatial_support"])
    return pd.Series(True, index=obs.index)


def _mad(values: np.ndarray) -> float:
    med = float(np.nanmedian(values))
    return float(np.nanmedian(np.abs(values - med)))


def qc_by_sample(adata, min_genes: int, max_genes: int, max_mito: float,
                 mad_multiplier: float, doublet_method: str) -> tuple[object, pd.DataFrame]:
    adata.var["mt"] = adata.var_names.str.lower().str.startswith("mt-")
    x = sp.csr_matrix(adata.X)
    total = np.asarray(x.sum(axis=1)).ravel()
    n_genes = np.asarray((x > 0).sum(axis=1)).ravel()
    mt_counts = np.asarray(x[:, adata.var["mt"].to_numpy()].sum(axis=1)).ravel()
    adata.obs["total_counts"] = total
    adata.obs["n_genes_by_counts"] = n_genes
    adata.obs["pct_counts_mt"] = np.divide(mt_counts * 100.0, total, out=np.zeros_like(total, dtype=float),
                                            where=total > 0)
    keep = pd.Series(False, index=adata.obs_names)
    rows = []
    for sample, index in adata.obs.groupby("sample", observed=True).groups.items():
        frame = adata.obs.loc[index]
        genes = frame["n_genes_by_counts"].to_numpy(float)
        counts = frame["total_counts"].to_numpy(float)
        mito = frame["pct_counts_mt"].to_numpy(float)
        gmed, cmed, mmed = map(float, (np.median(genes), np.median(counts), np.median(mito)))
        gmad, cmad, mmad = _mad(genes), _mad(counts), _mad(mito)
        gene_low = max(min_genes, gmed - 3 * mad_multiplier * max(gmad, 1.0))
        gene_high = min(max_genes, gmed + mad_multiplier * max(gmad, 1.0))
        count_high = cmed + mad_multiplier * max(cmad, 1.0)
        mito_high = min(max_mito, mmed + mad_multiplier * max(mmad, 0.1))
        sample_keep = ((frame.n_genes_by_counts >= gene_low) &
                       (frame.n_genes_by_counts <= gene_high) &
                       (frame.total_counts <= count_high) &
                       (frame.pct_counts_mt <= mito_high))
        keep.loc[index] = sample_keep
        rows.append(dict(sample=sample, n_before=len(frame), n_qc_pass=int(sample_keep.sum()),
                         gene_low=gene_low, gene_high=gene_high, count_high=count_high,
                         mito_high=mito_high, threshold_method="per_sample_median_MAD_plus_absolute_caps"))
    adata.obs["qc_pass_metrics"] = keep.reindex(adata.obs_names).to_numpy()

    if "is_doublet" in adata.obs:
        raw = adata.obs["is_doublet"]
        if raw.dtype == bool:
            predicted = raw.fillna(False)
        else:
            predicted = raw.astype(str).str.strip().str.lower().isin({"true", "1", "yes"})
        adata.obs["predicted_doublet"] = predicted.to_numpy()
        doublet_source = "author"
    elif doublet_method == "scrublet":
        import scanpy as sc
        adata.obs["predicted_doublet"] = False
        for sample, index in adata.obs.groupby("sample", observed=True).groups.items():
            if len(index) < 100:
                continue
            sub = adata[index].copy()
            try:
                sc.pp.scrublet(sub)
                adata.obs.loc[index, "predicted_doublet"] = sub.obs["predicted_doublet"].to_numpy()
            except Exception as exc:
                print(f"[warn] scrublet failed for {sample}: {type(exc).__name__}: {exc}")
        doublet_source = "scrublet_or_not_estimable_for_small_sample"
    else:
        adata.obs["predicted_doublet"] = False
        doublet_source = "not_run"
    adata.obs["doublet_source"] = doublet_source
    final_keep = adata.obs["qc_pass_metrics"] & ~adata.obs["predicted_doublet"].astype(bool)
    adata.obs["qc_pass"] = final_keep
    thresholds = pd.DataFrame(rows)
    thresholds["doublet_source"] = doublet_source
    thresholds["ambient_rna_status"] = "not_estimable_from_filtered_supplementary_matrix"
    adata = adata[final_keep].copy()
    gene_keep = np.asarray((sp.csr_matrix(adata.X) > 0).sum(axis=0)).ravel() >= 3
    return adata[:, gene_keep].copy(), thresholds


def annotate_and_embed(adata, cfg: dict, compute_embedding: bool) -> object:
    adata.layers["counts"] = adata.X.copy()
    author_col = cfg.get("author_annotation")
    needs_cluster = not author_col or author_col not in adata.obs
    if compute_embedding:
        import scanpy as sc
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        try:
            sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3",
                                        layer="counts", batch_key="sample")
        except Exception:
            sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="cell_ranger", batch_key="sample")
        sc.pp.pca(adata, n_comps=min(50, adata.n_vars - 1), mask_var="highly_variable")
        sc.pp.neighbors(adata, n_neighbors=15)
        sc.tl.leiden(adata, resolution=1.0, key_added="unsupervised_cluster")
        sc.tl.umap(adata)
    if needs_cluster and compute_embedding:
        adata.obs["cell_state_original"] = "cluster_" + adata.obs["unsupervised_cluster"].astype(str)
        adata.obs["cell_state_source"] = "unsupervised_unresolved"
    elif needs_cluster:
        adata.obs["cell_state_original"] = "unresolved_all_cells"
        adata.obs["cell_state_source"] = "unresolved_no_forced_annotation"
    else:
        labels = adata.obs[author_col].astype("string").fillna("author_unresolved")
        adata.obs["cell_state_original"] = labels.to_numpy()
        adata.obs["cell_state_source"] = f"author:{author_col}"
    # Harmonisation is a later, explicit ontology step.  Never silently relabel here.
    adata.obs["cell_state"] = adata.obs["cell_state_original"].astype("category")
    return adata


def stable_group_id(dataset: str, sample: str, cell_state: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9._-]+", "_", cell_state).strip("_")[:48] or "state"
    suffix = hashlib.sha1(cell_state.encode("utf-8")).hexdigest()[:8]
    return f"{dataset}::{sample}::{slug}::{suffix}"


def bool_series(values: pd.Series, field: str) -> pd.Series:
    if values.dtype == bool:
        return values.fillna(False)
    normalized = values.astype("string").str.strip().str.lower()
    known = normalized.isin(["true", "false", "1", "0", "yes", "no"])
    if not known.all():
        bad = normalized[~known].drop_duplicates().tolist()[:5]
        raise ValueError(f"ambiguous boolean values in {field}: {bad}")
    return normalized.isin(["true", "1", "yes"])


def coldata_from_stream_obs(obs: pd.DataFrame, codes: np.ndarray, levels: pd.Index) -> pd.DataFrame:
    rows = []
    for code, level in enumerate(levels):
        frame = obs.iloc[np.flatnonzero(codes == code)]
        sample, state = level.split("\x1f", 1)
        invariant = ["dataset", "subject_id", "condition", "group", "dpi", "phase", "region", "segment"]
        values = {}
        for column in invariant:
            unique = frame[column].drop_duplicates()
            if len(unique) != 1:
                raise ValueError(f"non-constant {column} within {sample}/{state}: {unique.tolist()[:5]}")
            values[column] = unique.iloc[0]
        rows.append(dict(group_id=stable_group_id(values["dataset"], sample, state), sample=sample,
                         cell_state=state, n_cells=len(frame), **values))
    return pd.DataFrame(rows).set_index("group_id")


def stream_author_pseudobulk(dataset: str, cfg: dict, raw_root: Path, out_dir: Path,
                             sample_metadata: Path, mode: str, min_cells: int,
                             chunk_rows: int) -> None:
    """Memory-bounded MatrixMarket aggregation for a dataset with trusted cell metadata.

    This is currently used for GSE162610 (151,791,147 non-zero entries). Author QC flags
    and doublet calls are applied before streaming; the raw integer matrix is aggregated
    directly to GSM x author-celltype groups.
    """
    raw_dir = raw_root / dataset
    matrix_path = raw_dir / cfg["matrix"]
    for name in (cfg["matrix"], cfg["features"], cfg["barcodes"], cfg["meta"]):
        if not (raw_dir / name).exists():
            raise FileNotFoundError(raw_dir / name)
    genes = read_one_column(raw_dir / cfg["features"], header=cfg.get("feature_header"),
                            sep=cfg.get("feature_sep", "\t"))
    barcodes = read_one_column(raw_dir / cfg["barcodes"], header=cfg.get("barcode_header"),
                               sep=cfg.get("barcode_sep", "\t"))
    if pd.Index(genes).duplicated().any():
        raise ValueError("streaming pseudobulk requires unique gene names")
    if cfg["meta_join"] == "index":
        author = pd.read_csv(raw_dir / cfg["meta"], sep="\t", index_col=0, low_memory=False)
        author.index = author.index.astype(str)
        if len(author) != len(barcodes) or not pd.Index(barcodes).isin(author.index).all():
            raise ValueError("author metadata/barcode contract failed for streaming aggregation")
        obs = author.reindex(barcodes).copy()
    elif cfg["meta_join"] == "barcode":
        author = pd.read_csv(raw_dir / cfg["meta"], sep="\t", low_memory=False)
        key = cfg["barcode_col"]
        if author[key].duplicated().any():
            raise ValueError(f"duplicate author barcode in {key}")
        obs = author.set_index(author[key].astype(str)).reindex(barcodes).copy()
        if obs.isna().all(axis=1).any():
            raise ValueError("author barcode join produced empty rows")
    elif cfg["meta_join"] == "position":
        author = pd.read_csv(raw_dir / cfg["meta"], sep="\t", low_memory=False)
        if len(author) != len(barcodes):
            raise ValueError(f"positional metadata has {len(author)} rows; expected {len(barcodes)}")
        author.index = pd.Index(barcodes)
        obs = author.copy()
    else:
        raise ValueError(f"unsupported streaming meta_join={cfg['meta_join']}")
    if dataset == "GSE172167":
        obs["library_key"] = obs["orig.ident"].astype(str)
    else:
        obs["library_key"] = obs[cfg["library_col"]].astype(str)
    canonical = pd.read_csv(sample_metadata, sep="\t", dtype={"sample": str, "library_key": str})
    lookup = canonical_sample_lookup(canonical, dataset, cfg["modality"])
    missing = sorted(set(obs.library_key) - set(lookup.index.astype(str)))
    if missing:
        raise ValueError(f"unmapped streaming library keys: {missing[:10]}")
    for column in CANONICAL_COLUMNS:
        obs[column] = obs.library_key.map(lookup[column]).to_numpy()
    obs["dataset"] = dataset
    obs["cell_state"] = obs[cfg["author_annotation"]].astype("string")
    scope = analysis_scope_mask(obs, mode)
    if not scope.any():
        raise ValueError(f"no cells remain in requested analysis mode={mode}")
    validate_obs_contract(obs.loc[scope], inferential=(mode == "inferential"))
    qc_mode = cfg.get("stream_qc_mode", "filtered_no_metrics")
    if qc_mode == "author_flags":
        required_qc = ["pass_umi", "pass_n_genes", "pass_percent_mt", "is_doublet"]
        missing_qc = [column for column in required_qc if column not in obs]
        if missing_qc:
            raise ValueError(f"streaming path requires author QC columns: {missing_qc}")
        qc_pass = (scope & bool_series(obs.pass_umi, "pass_umi") &
                   bool_series(obs.pass_n_genes, "pass_n_genes") &
                   bool_series(obs.pass_percent_mt, "pass_percent_mt") &
                   ~bool_series(obs.is_doublet, "is_doublet") & obs.cell_state.notna())
        qc_method = "author_pass_umi+pass_n_genes+pass_percent_mt+not_doublet"
        doublet_status = "author_calls_applied"
    elif qc_mode == "filtered_metrics":
        required = ["nCount_RNA", "nFeature_RNA"]
        missing_metrics = [column for column in required if column not in obs]
        if missing_metrics:
            raise ValueError(f"filtered-metrics QC missing {missing_metrics}")
        qc_pass = pd.Series(False, index=obs.index)
        for sample, index in obs.loc[scope].groupby("sample", observed=True).groups.items():
            counts = pd.to_numeric(obs.loc[index, "nCount_RNA"], errors="coerce")
            features = pd.to_numeric(obs.loc[index, "nFeature_RNA"], errors="coerce")
            c_med, f_med = counts.median(), features.median()
            c_mad = (counts - c_med).abs().median()
            f_mad = (features - f_med).abs().median()
            keep = (features >= max(200, f_med - 15 * max(f_mad, 1))) & \
                   (features <= min(7500, f_med + 5 * max(f_mad, 1))) & \
                   (counts <= c_med + 5 * max(c_mad, 1))
            qc_pass.loc[index] = keep.fillna(False)
        qc_pass &= scope & obs.cell_state.notna()
        qc_method = "filtered_matrix_per_sample_count_and_feature_MAD;mitochondrial_fraction_unavailable"
        doublet_status = "not_available_in_supplementary_metadata"
    elif qc_mode == "filtered_no_metrics":
        qc_pass = scope & obs.cell_state.notna()
        qc_method = "published_filtered_matrix_no_additional_metric_filter"
        doublet_status = "not_available_in_supplementary_metadata"
    else:
        raise ValueError(f"unknown stream_qc_mode={qc_mode}")
    keys = obs[["sample", "cell_state"]].astype(str).agg("\x1f".join, axis=1)
    group_cell_counts = keys[qc_pass].value_counts(sort=False)
    eligible_levels = group_cell_counts[group_cell_counts >= min_cells].index
    retained = qc_pass & keys.isin(eligible_levels)
    kept_obs = obs.loc[retained].copy()
    kept_keys = keys.loc[retained]
    codes, levels = pd.factorize(kept_keys, sort=True)
    levels = pd.Index(levels)
    cell_to_group = np.full(len(obs), -1, dtype=np.int32)
    cell_to_group[np.flatnonzero(retained.to_numpy())] = codes.astype(np.int32)
    coldata = coldata_from_stream_obs(kept_obs, codes, levels)
    dropped = group_cell_counts[group_cell_counts < min_cells].rename("n_cells").reset_index()
    dropped = dropped.rename(columns={dropped.columns[0]: "sample_cell_state_key"})

    stream_format = cfg.get("stream_format", "matrix_market_coordinate")
    values_seen = 0
    nonzero_seen = 0
    if stream_format == "matrix_market_coordinate":
        aggregate = sp.csr_matrix((len(levels), len(genes)), dtype=np.int64)
        with gzip.open(matrix_path, "rt", encoding="ascii", newline="") as handle:
            line = handle.readline()
            if not line.startswith("%%MatrixMarket matrix coordinate"):
                raise ValueError("expected coordinate MatrixMarket input")
            line = handle.readline()
            while line.startswith("%"):
                line = handle.readline()
            n_gene, n_cell, expected_entries = map(int, line.split())
            if (n_gene, n_cell) != (len(genes), len(barcodes)):
                raise ValueError(f"matrix header {(n_gene, n_cell)} != genes/barcodes {(len(genes), len(barcodes))}")
            reader = pd.read_csv(handle, sep=" ", header=None, names=["gene", "cell", "count"],
                                 dtype=np.int32, chunksize=chunk_rows, engine="c")
            for chunk_i, chunk in enumerate(reader, start=1):
                values_seen += len(chunk)
                nonzero_seen += int((chunk["count"] != 0).sum())
                gene_idx = chunk.gene.to_numpy(np.int32) - 1
                cell_idx = chunk.cell.to_numpy(np.int32) - 1
                if gene_idx.min() < 0 or gene_idx.max() >= len(genes) or cell_idx.min() < 0 or cell_idx.max() >= len(obs):
                    raise ValueError("out-of-bounds MatrixMarket coordinate")
                group_idx = cell_to_group[cell_idx]
                use = group_idx >= 0
                if use.any():
                    part = sp.coo_matrix((chunk.loc[use, "count"].to_numpy(np.int64),
                                          (group_idx[use], gene_idx[use])),
                                         shape=aggregate.shape, dtype=np.int64).tocsr()
                    aggregate = aggregate + part
                if chunk_i % 10 == 0:
                    print(f"[{dataset}] streamed {values_seen:,}/{expected_entries:,} non-zero entries", flush=True)
        if values_seen != expected_entries:
            raise ValueError(f"stream ended at {values_seen} entries; expected {expected_entries}")
    elif stream_format == "dense_whitespace_header":
        dense_aggregate = np.zeros((len(levels), len(genes)), dtype=np.int64)
        valid_cells = np.flatnonzero(cell_to_group >= 0)
        valid_groups = cell_to_group[valid_cells]
        with gzip.open(matrix_path, "rt", encoding="ascii", newline="") as handle:
            header_barcodes = handle.readline().split()
            if header_barcodes != list(barcodes):
                mismatch = next((i for i, (a, b) in enumerate(zip(header_barcodes, barcodes)) if a != b), None)
                raise ValueError(f"dense matrix header/barcode order mismatch at position {mismatch}")
            gene_idx = 0
            pending = np.empty(0, dtype=np.int64)
            for line in handle:
                physical = np.fromstring(line, dtype=np.int64, sep=" ")
                values = np.concatenate([pending, physical]) if len(pending) else physical
                while len(values) >= len(barcodes):
                    if gene_idx >= len(genes):
                        raise ValueError("dense matrix contains more gene rows than feature table")
                    row_values = values[:len(barcodes)]
                    values = values[len(barcodes):]
                    selected = row_values[valid_cells]
                    dense_aggregate[:, gene_idx] = np.bincount(
                        valid_groups, weights=selected, minlength=len(levels)).astype(np.int64)
                    values_seen += len(row_values)
                    nonzero_seen += int(np.count_nonzero(row_values))
                    gene_idx += 1
                    if gene_idx % 1000 == 0:
                        print(f"[{dataset}] streamed {gene_idx:,}/{len(genes):,} dense gene rows", flush=True)
                pending = values
            if len(pending):
                raise ValueError(f"dense matrix ended with {len(pending)} unassigned values")
            if gene_idx != len(genes):
                raise ValueError(f"dense matrix ended at {gene_idx} genes; expected {len(genes)}")
        aggregate = sp.csr_matrix(dense_aggregate)
    else:
        raise ValueError(f"unknown stream_format={stream_format}")

    out_dir.mkdir(parents=True, exist_ok=True)
    write_sparse_mtx(out_dir / f"{dataset}_pseudobulk_counts.mtx.gz", aggregate)
    pd.DataFrame({"gene": genes}).to_csv(out_dir / f"{dataset}_pseudobulk_genes.tsv", sep="\t", index=False)
    coldata.to_csv(out_dir / f"{dataset}_pseudobulk_coldata.tsv", sep="\t")
    dropped.to_csv(out_dir / f"{dataset}_pseudobulk_dropped_groups.tsv", sep="\t", index=False)
    qc_summary = (obs.assign(in_analysis_scope=scope, qc_pass=qc_pass, retained=retained)
                  .groupby("sample", observed=True)
                  .agg(n_before=("sample", "size"), n_in_analysis_scope=("in_analysis_scope", "sum"),
                       n_qc_pass=("qc_pass", "sum"),
                       n_retained_min_cell_gate=("retained", "sum")).reset_index())
    qc_summary["threshold_method"] = qc_method
    qc_summary["doublet_status"] = doublet_status
    qc_summary["ambient_rna_status"] = "not_estimable_from_filtered_supplementary_matrix"
    qc_summary.to_csv(out_dir / f"{dataset}_qc_thresholds.tsv", sep="\t", index=False)
    provenance = dict(dataset=dataset, mode=mode, processing_mode="streaming_author_annotation_pseudobulk",
                      matrix_format=stream_format, matrix_shape_genes_by_cells=[len(genes), len(barcodes)],
                      matrix_values_streamed=values_seen, matrix_nonzero_entries=nonzero_seen,
                      cells_in_analysis_scope=int(scope.sum()), cells_qc_pass=int(qc_pass.sum()),
                      cells_retained_after_group_gate=int(retained.sum()),
                      pseudobulk_groups=int(aggregate.shape[0]), min_cells_per_group=min_cells,
                      author_annotation=cfg["author_annotation"], qc_method=qc_method,
                      doublet_status=doublet_status, h5ad_written=False,
                      h5ad_reason="source matrix exceeds safe materialization memory in current environment",
                      ambient_rna_status="not_estimable_from_filtered_supplementary_matrix")
    (out_dir / f"{dataset}_processing_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[{dataset}] streamed {values_seen:,} values ({nonzero_seen:,} non-zero) -> "
          f"{aggregate.shape[0]} pseudobulk groups; "
          f"{int(retained.sum()):,}/{len(obs):,} cells retained", flush=True)


def sparse_pseudobulk(adata, min_cells: int) -> tuple[sp.csr_matrix, pd.DataFrame, pd.DataFrame]:
    validate_obs_contract(adata.obs, inferential=False)
    x = adata.layers["counts"] if "counts" in adata.layers else adata.X
    x = sp.csr_matrix(x)
    keys = adata.obs[["sample", "cell_state"]].astype(str).agg("\x1f".join, axis=1)
    counts = keys.value_counts(sort=False)
    retained_keys = counts[counts >= min_cells].index
    retained_mask = keys.isin(retained_keys).to_numpy()
    dropped = counts[counts < min_cells].rename("n_cells").reset_index()
    dropped = dropped.rename(columns={dropped.columns[0]: "sample_cell_state_key"})
    obs = adata.obs.loc[retained_mask].copy()
    x = x[retained_mask]
    keys = keys.loc[retained_mask]
    codes, levels = pd.factorize(keys, sort=True)
    aggregator = sp.csr_matrix((np.ones(len(codes), dtype=np.int8), (codes, np.arange(len(codes)))),
                               shape=(len(levels), len(codes)))
    pb = (aggregator @ x).tocsr()
    rows = []
    for code, level in enumerate(levels):
        idx = np.flatnonzero(codes == code)
        frame = obs.iloc[idx]
        sample, state = level.split("\x1f", 1)
        invariant = ["dataset", "subject_id", "condition", "group", "dpi", "phase", "region", "segment"]
        values = {}
        for column in invariant:
            unique = frame[column].drop_duplicates()
            if len(unique) != 1:
                raise ValueError(f"non-constant {column} within {sample}/{state}: {unique.tolist()[:5]}")
            values[column] = unique.iloc[0]
        rows.append(dict(group_id=stable_group_id(values["dataset"], sample, state), sample=sample,
                         cell_state=state, n_cells=len(frame), **values))
    coldata = pd.DataFrame(rows).set_index("group_id")
    if coldata.index.duplicated().any():
        raise ValueError("stable pseudobulk group_id collision")
    return pb, coldata, dropped


def write_sparse_mtx(path: Path, matrix: sp.spmatrix) -> None:
    with gzip.open(path, "wb") as handle:
        scipy.io.mmwrite(handle, matrix, field="integer" if np.issubdtype(matrix.dtype, np.integer) else "real")


def process(dataset: str, cfg: dict, raw_root: Path, out_dir: Path, sample_metadata: Path,
            mode: str, min_cells: int, compute_embedding: bool, doublet_method: str,
            min_genes: int, max_genes: int, max_mito: float, mad_multiplier: float,
            write_dense_tsv: bool, streaming: str, stream_chunk_rows: int) -> None:
    use_streaming = streaming == "always" or (streaming == "auto" and cfg.get("stream_mtx", False))
    if use_streaming:
        if cfg["kind"] != "mtx" or not cfg.get("author_annotation"):
            raise ValueError("streaming mode requires MatrixMarket input and author annotation")
        if compute_embedding or write_dense_tsv:
            raise ValueError("streaming mode cannot compute an embedding or write a dense count table")
        stream_author_pseudobulk(dataset, cfg, raw_root, out_dir, sample_metadata, mode,
                                 min_cells, stream_chunk_rows)
        return
    raw_dir = raw_root / dataset
    matrix_path = raw_dir / cfg["matrix"]
    if not matrix_path.exists():
        raise FileNotFoundError(f"expression matrix not downloaded: {matrix_path}")
    adata = read_mtx_dataset(cfg, raw_dir) if cfg["kind"] == "mtx" else read_csv_dataset(cfg, raw_dir)
    attach_author_metadata(adata, dataset, cfg, raw_dir)
    attach_canonical_metadata(adata, dataset, cfg, sample_metadata)
    scope = analysis_scope_mask(adata.obs, mode)
    if not scope.any():
        raise ValueError(f"no cells remain in requested analysis mode={mode}")
    adata = adata[scope.to_numpy()].copy()
    validate_obs_contract(adata.obs, inferential=(mode == "inferential"))
    n_loaded = adata.n_obs
    adata, thresholds = qc_by_sample(adata, min_genes, max_genes, max_mito,
                                     mad_multiplier, doublet_method)
    adata = annotate_and_embed(adata, cfg, compute_embedding)
    pb, coldata, dropped = sparse_pseudobulk(adata, min_cells)

    out_dir.mkdir(parents=True, exist_ok=True)
    prefix = out_dir / dataset
    adata.write_h5ad(prefix.with_suffix(".h5ad"), compression="gzip")
    write_sparse_mtx(out_dir / f"{dataset}_pseudobulk_counts.mtx.gz", pb)
    pd.Series(adata.var_names, name="gene").to_csv(
        out_dir / f"{dataset}_pseudobulk_genes.tsv", sep="\t", index=False)
    coldata.to_csv(out_dir / f"{dataset}_pseudobulk_coldata.tsv", sep="\t")
    thresholds.to_csv(out_dir / f"{dataset}_qc_thresholds.tsv", sep="\t", index=False)
    dropped.to_csv(out_dir / f"{dataset}_pseudobulk_dropped_groups.tsv", sep="\t", index=False)
    if write_dense_tsv:
        pd.DataFrame(pb.toarray(), index=coldata.index, columns=adata.var_names).to_csv(
            out_dir / f"{dataset}_pseudobulk_counts.tsv", sep="\t")
    provenance = dict(dataset=dataset, mode=mode, cells_loaded=n_loaded, cells_retained=adata.n_obs,
                      genes=adata.n_vars, pseudobulk_groups=pb.shape[0], min_cells_per_group=min_cells,
                      author_annotation=cfg.get("author_annotation"), embedding_computed=compute_embedding,
                      ambient_rna_status="not_estimable_from_filtered_supplementary_matrix",
                      inference_warning=("filtered matrix lacks empty droplets; ambient correction cannot be estimated"))
    (out_dir / f"{dataset}_processing_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[{dataset}] {n_loaded} -> {adata.n_obs} cells; {pb.shape[0]} pseudobulk groups; mode={mode}")


def selftest() -> None:
    good = pd.DataFrame(dict(dataset=["GSE1"], sample=["GSM1"], subject_id=["GSM1"],
                             condition=["injured"], group=["dpi_1"], library_key=["lib1"],
                             inference_role=["within_study_inference"]))
    validate_obs_contract(good, inferential=True)
    bad = good.copy(); bad["sample"] = "S1"
    try:
        validate_obs_contract(bad, inferential=True)
    except ValueError as exc:
        assert "dummy sample" in str(exc)
    else:
        raise AssertionError("dummy sample was accepted")
    assert stable_group_id("GSE1", "GSM1", "Microglia activated") == stable_group_id(
        "GSE1", "GSM1", "Microglia activated")
    meta = pd.DataFrame(dict(dataset=["GSE1"], modality=["scRNA-seq"], library_key=["lib1"], sample=["GSM1"]))
    assert canonical_sample_lookup(meta, "GSE1", "scRNA-seq").loc["lib1", "sample"] == "GSM1"
    print("[selftest] PASS: canonical library contract enforced; dummy S1 rejected; group IDs stable")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", choices=sorted(DATASETS))
    ap.add_argument("--all-available", action="store_true")
    ap.add_argument("--raw-root", default="data_raw")
    ap.add_argument("--out-dir", default="data_processed")
    ap.add_argument("--sample-metadata", default="tables/sample_metadata.tsv")
    ap.add_argument("--mode", choices=["inferential", "support", "descriptive"])
    ap.add_argument("--min-cells-per-pseudobulk", type=int, default=20)
    ap.add_argument("--min-genes", type=int, default=200)
    ap.add_argument("--max-genes", type=int, default=7500)
    ap.add_argument("--max-mito", type=float, default=15.0)
    ap.add_argument("--mad-multiplier", type=float, default=5.0)
    ap.add_argument("--doublet-method", choices=["author_or_none", "scrublet"], default="author_or_none")
    ap.add_argument("--compute-embedding", action="store_true")
    ap.add_argument("--write-dense-tsv", action="store_true")
    ap.add_argument("--streaming", choices=["auto", "always", "never"], default="auto")
    ap.add_argument("--stream-chunk-rows", type=int, default=2_000_000)
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    if args.selftest:
        selftest(); return
    names = []
    if args.dataset:
        names.append(args.dataset)
    if args.all_available:
        for name, cfg in DATASETS.items():
            if (Path(args.raw_root) / name / cfg["matrix"]).exists() and name not in names:
                names.append(name)
    if not names:
        ap.error("provide --dataset GSE... or --all-available")
    for name in names:
        cfg = DATASETS[name]
        process(name, cfg, Path(args.raw_root), Path(args.out_dir), Path(args.sample_metadata),
                args.mode or cfg["default_mode"], args.min_cells_per_pseudobulk,
                args.compute_embedding, args.doublet_method, args.min_genes, args.max_genes,
                args.max_mito, args.mad_multiplier, args.write_dense_tsv,
                args.streaming, args.stream_chunk_rows)


if __name__ == "__main__":
    main()
