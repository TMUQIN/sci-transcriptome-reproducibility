#!/usr/bin/env Rscript
# Cross-dataset cell-state replicability with MetaNeighbor (Crow et al., Nat Commun 2018).
#
# Answers "which annotated SCI cell states replicate across independent datasets?" via
# unsupervised AUROC. States with high mean AUROC are trustworthy for the reproducible
# cell-state set (draft_v2 Result 2); low-AUROC states are flagged as dataset-specific.
# Runs in YOUR env (needs MetaNeighbor, SummarizedExperiment). Not runnable in the sandbox.
#
# Input : a merged SingleCellExperiment / matrix of shared HVGs across datasets, with
#         colData columns study_id (dataset) and cell_type (harmonized label).
# Output: results/metaneighbor_auroc.tsv (cell_type x replicability AUROC) + a summary.
#
# Usage:
#   Rscript 10_metaneighbor_replicability.R --sce data_processed/merged_for_metaneighbor.rds \
#       --study-col study_id --celltype-col cell_type --out results/metaneighbor_auroc.tsv
suppressMessages({library(MetaNeighbor); library(SummarizedExperiment)})

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- which(args == flag); if (length(i) && i < length(args)) args[i + 1] else default
}
sce_path   <- getarg("--sce")
study_col  <- getarg("--study-col", "study_id")
ct_col     <- getarg("--celltype-col", "cell_type")
out_path   <- getarg("--out", "results/metaneighbor_auroc.tsv")
stopifnot(!is.null(sce_path))
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

sce <- readRDS(sce_path)
# highly variable genes shared across datasets (unsupervised MetaNeighbor)
vg <- variableGenes(dat = sce, exp_labels = colData(sce)[[study_col]])
celltype_NV <- MetaNeighborUS(var_genes = vg, dat = sce,
                              study_id = colData(sce)[[study_col]],
                              cell_type = colData(sce)[[ct_col]], fast_version = TRUE)
# mean cross-study AUROC per harmonized cell type (diagonal blocks vs off-diagonal)
top_hits <- topHits(cell_NV = celltype_NV, dat = sce,
                    study_id = colData(sce)[[study_col]],
                    cell_type = colData(sce)[[ct_col]], threshold = 0.9)
write.table(as.data.frame(as.table(celltype_NV)), out_path, sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(top_hits, sub("\\.tsv$", "_tophits.tsv", out_path), sep = "\t",
            quote = FALSE, row.names = FALSE)
message("Wrote ", out_path, " and *_tophits.tsv. ",
        "High-AUROC (>0.9) states = replicable across datasets -> trust in Result 2; ",
        "low-AUROC states -> flag as dataset-specific.")
