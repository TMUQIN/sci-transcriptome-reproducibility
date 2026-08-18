#!/usr/bin/env Rscript

# Pre-registered orthogonal cross-platform context evaluation for GSE47681.
#
# This script intentionally analyzes only the 13 WT SHAM/day-1/day-7 arrays
# specified in reports/phase_v5_3_3_postlock_upgrade_2026_08/
# provenance/GSE47681_orthogonal_microarray_protocol.md. It does not alter or pool
# results from the frozen three-study branch.

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1L) normalizePath(args[[1L]], mustWork = TRUE) else normalizePath(".", mustWork = TRUE)
raw_dir <- if (length(args) >= 2L) normalizePath(args[[2L]], mustWork = TRUE) else file.path(root, "data_raw", "GSE47681", "raw_cel")
out_dir <- if (length(args) >= 3L) args[[3L]] else file.path(root, "reports", "phase_v5_3_3_postlock_upgrade_2026_08", "GSE47681")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

local_lib <- file.path(root, ".deps", "R")
if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))

suppressPackageStartupMessages({
  library(affy)
  library(mouse4302cdf)
  library(mouse4302.db)
  library(AnnotationDbi)
  library(limma)
})

expected <- data.frame(
  gsm = c(
    "GSM1154538", "GSM1154539", "GSM1154540", "GSM1154542",
    "GSM1154530", "GSM1154541", "GSM1154547", "GSM1154548", "GSM1154549",
    "GSM1154529", "GSM1154531", "GSM1154532", "GSM1154533"
  ),
  group = c(rep("SHAMWT", 4), rep("DAY1WT", 5), rep("DAY7WT", 4)),
  biological_replicate = c(1:4, 1:5, 1:4),
  stringsAsFactors = FALSE
)

all_files <- list.files(raw_dir, pattern = "\\.CEL\\.gz$", full.names = TRUE, ignore.case = TRUE)
if (length(all_files) != 13L) stop("Expected exactly 13 eligible CEL.gz files; found ", length(all_files))
file_gsm <- sub("_.*$", "", basename(all_files))
if (anyDuplicated(file_gsm)) stop("Duplicate GSM accession among CEL files")
if (!setequal(file_gsm, expected$gsm)) {
  stop("CEL accessions do not match the frozen eligibility list. Missing: ",
       paste(setdiff(expected$gsm, file_gsm), collapse = ","), "; unexpected: ",
       paste(setdiff(file_gsm, expected$gsm), collapse = ","))
}
expected$file <- all_files[match(expected$gsm, file_gsm)]
expected$filename <- basename(expected$file)
expected$compressed_bytes <- file.info(expected$file)$size

write.table(expected[, c("gsm", "group", "biological_replicate", "filename", "compressed_bytes")],
            file.path(out_dir, "GSE47681_WT_sample_manifest.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)

message("Reading 13 pre-registered WT CEL arrays")
raw <- ReadAffy(filenames = expected$file)
sampleNames(raw) <- expected$gsm
if (length(sampleNames(raw)) != 13L) stop("AffyBatch did not contain 13 arrays")
if (!identical(sampleNames(raw), expected$gsm)) stop("AffyBatch sample order differs from manifest")
if (!identical(annotation(raw), "mouse4302")) stop("Unexpected Affymetrix annotation: ", annotation(raw))

raw_intensity <- intensity(raw)
raw_quantile_probs <- c(0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)
raw_q <- t(apply(raw_intensity, 2, quantile, probs = raw_quantile_probs, na.rm = TRUE, names = FALSE))
colnames(raw_q) <- paste0("q", gsub("\\.", "_", format(raw_quantile_probs, trim = TRUE)))
raw_q <- data.frame(gsm = rownames(raw_q), raw_q, check.names = FALSE)
rownames(raw_q) <- NULL
write.table(raw_q, file.path(out_dir, "GSE47681_raw_intensity_quantiles.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
rm(raw_intensity)
invisible(gc())

message("Running joint RMA background correction, quantile normalization and median-polish summarization")
rma_eset <- rma(raw, background = TRUE, normalize = TRUE, verbose = TRUE)
probe_expr <- exprs(rma_eset)
if (ncol(probe_expr) != 13L || any(!is.finite(probe_expr))) stop("Invalid RMA expression matrix")
if (!identical(colnames(probe_expr), expected$gsm)) stop("RMA output sample order differs from manifest")

norm_q <- t(apply(probe_expr, 2, quantile, probs = raw_quantile_probs, na.rm = TRUE, names = FALSE))
colnames(norm_q) <- paste0("q", gsub("\\.", "_", format(raw_quantile_probs, trim = TRUE)))
norm_q <- data.frame(gsm = rownames(norm_q), norm_q, check.names = FALSE)
rownames(norm_q) <- NULL
write.table(norm_q, file.path(out_dir, "GSE47681_RMA_expression_quantiles.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

# RLE-like diagnostics on the normalized probe-set matrix.
row_medians <- apply(probe_expr, 1, median, na.rm = TRUE)
rle <- sweep(probe_expr, 1, row_medians, FUN = "-")
rle_q <- t(apply(rle, 2, quantile, probs = c(0.05, 0.25, 0.5, 0.75, 0.95), na.rm = TRUE, names = FALSE))
colnames(rle_q) <- c("rle_q05", "rle_q25", "rle_median", "rle_q75", "rle_q95")
rle_q <- data.frame(gsm = rownames(rle_q), rle_q, rle_iqr = rle_q[, "rle_q75"] - rle_q[, "rle_q25"], check.names = FALSE)
rownames(rle_q) <- NULL
write.table(rle_q, file.path(out_dir, "GSE47681_RLE_diagnostics.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
rm(rle)
invisible(gc())

# Deterministic PCA on the 5,000 most variable probe sets for sample-level QC.
probe_var <- apply(probe_expr, 1, var)
top_n <- min(5000L, length(probe_var))
top_idx <- order(probe_var, rownames(probe_expr), decreasing = c(TRUE, FALSE), method = "radix")[seq_len(top_n)]
pca <- prcomp(t(probe_expr[top_idx, , drop = FALSE]), center = TRUE, scale. = FALSE)
pca_scores <- data.frame(gsm = rownames(pca$x), group = expected$group[match(rownames(pca$x), expected$gsm)],
                         pca$x[, seq_len(min(5L, ncol(pca$x))), drop = FALSE], check.names = FALSE)
pca_var <- data.frame(component = paste0("PC", seq_along(pca$sdev)),
                      variance_explained = pca$sdev^2 / sum(pca$sdev^2))
write.table(pca_scores, file.path(out_dir, "GSE47681_RMA_PCA_scores.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(pca_var, file.path(out_dir, "GSE47681_RMA_PCA_variance.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

sample_cor <- cor(probe_expr, method = "pearson")
sample_cor_out <- data.frame(gsm = rownames(sample_cor), sample_cor, check.names = FALSE)
write.table(sample_cor_out, file.path(out_dir, "GSE47681_RMA_sample_correlations.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

# Fixed GPL1261 probe-to-symbol mapping and deterministic IQR collapse.
probe_ids <- rownames(probe_expr)
map <- AnnotationDbi::select(mouse4302.db, keys = probe_ids, columns = "SYMBOL", keytype = "PROBEID")
map <- map[!is.na(map$SYMBOL) & nzchar(trimws(map$SYMBOL)), c("PROBEID", "SYMBOL"), drop = FALSE]

# Expand the rare delimiter-encoded multi-symbol record if present; duplicated
# rows returned directly by AnnotationDbi are naturally retained.
expanded <- lapply(seq_len(nrow(map)), function(i) {
  symbols <- trimws(unlist(strsplit(as.character(map$SYMBOL[[i]]), "\\s*(///|;|\\|)\\s*", perl = TRUE)))
  symbols <- unique(symbols[nzchar(symbols)])
  data.frame(PROBEID = as.character(map$PROBEID[[i]]), SYMBOL = symbols, stringsAsFactors = FALSE)
})
map <- do.call(rbind, expanded)
map <- unique(map)
map <- map[map$PROBEID %in% probe_ids, , drop = FALSE]

probe_iqr <- apply(probe_expr, 1, IQR, na.rm = TRUE)
map$probe_iqr <- unname(probe_iqr[map$PROBEID])
map <- map[is.finite(map$probe_iqr), , drop = FALSE]
map <- map[order(map$SYMBOL, -map$probe_iqr, map$PROBEID, method = "radix"), , drop = FALSE]
selected <- map[!duplicated(map$SYMBOL), , drop = FALSE]
selected$selection_rule <- "maximum_across_sample_IQR_then_lexicographic_probe_id"
gene_expr <- probe_expr[selected$PROBEID, , drop = FALSE]
rownames(gene_expr) <- selected$SYMBOL
if (anyDuplicated(rownames(gene_expr))) stop("Gene-level collapse left duplicate symbols")
if (any(!is.finite(gene_expr))) stop("Gene-level expression contains non-finite values")
write.table(selected, file.path(out_dir, "GSE47681_selected_probe_by_gene.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

gene_con <- gzfile(file.path(out_dir, "GSE47681_RMA_gene_expression.tsv.gz"), open = "wt")
write.table(data.frame(gene = rownames(gene_expr), gene_expr, check.names = FALSE), gene_con,
            sep = "\t", quote = FALSE, row.names = FALSE)
close(gene_con)

# Frozen no-intercept limma model and three contrasts.
group <- factor(expected$group, levels = c("SHAMWT", "DAY1WT", "DAY7WT"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
rownames(design) <- expected$gsm
if (qr(design)$rank != ncol(design)) stop("Rank-deficient design")
contrast_matrix <- makeContrasts(
  injury_d1_vs_uninjured = DAY1WT - SHAMWT,
  injury_d7_vs_uninjured = DAY7WT - SHAMWT,
  change_d7_minus_d1 = DAY7WT - DAY1WT,
  levels = design
)
fit <- lmFit(gene_expr, design)
fit2 <- eBayes(contrasts.fit(fit, contrast_matrix), robust = TRUE)

contrast_meta <- data.frame(
  contrast_id = colnames(contrast_matrix),
  effect_type = c("endpoint", "endpoint", "temporal_delta"),
  numerator_group = c("DAY1WT", "DAY7WT", "DAY7WT"),
  denominator_group = c("SHAMWT", "SHAMWT", "DAY1WT"),
  n_num = c(5L, 4L, 4L),
  n_den = c(4L, 4L, 5L),
  stringsAsFactors = FALSE
)

effect_list <- vector("list", nrow(contrast_meta))
for (j in seq_len(nrow(contrast_meta))) {
  estimate <- fit2$coefficients[, j]
  se <- sqrt(fit2$s2.post) * fit2$stdev.unscaled[, j]
  p_value <- fit2$p.value[, j]
  df_total <- fit2$df.total
  crit <- qt(0.975, df = df_total)
  effect_list[[j]] <- data.frame(
    gene = rownames(fit2$coefficients), cell_state = "whole_lesion_bulk_microarray",
    dataset = "GSE47681", axis = "time", contrast_id = contrast_meta$contrast_id[[j]],
    effect_type = contrast_meta$effect_type[[j]], numerator_group = contrast_meta$numerator_group[[j]],
    denominator_group = contrast_meta$denominator_group[[j]], primary = "stress_test",
    estimate = estimate, se = se, df_residual = df_total,
    ci_low = estimate - crit * se, ci_high = estimate + crit * se,
    p = p_value, fdr = p.adjust(p_value, method = "BH"),
    n_num = contrast_meta$n_num[[j]], n_den = contrast_meta$n_den[[j]],
    stringsAsFactors = FALSE, check.names = FALSE
  )
}
effects <- do.call(rbind, effect_list)
write.table(effects, file.path(out_dir, "effects_GSE47681.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)

design_out <- data.frame(gsm = rownames(design), group = expected$group, design, check.names = FALSE)
write.table(design_out, file.path(out_dir, "GSE47681_design_matrix.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(data.frame(coefficient = rownames(contrast_matrix), contrast_matrix, check.names = FALSE),
            file.path(out_dir, "GSE47681_contrast_matrix.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)

qc_summary <- data.frame(
  item = c("arrays_read", "eligible_SHAMWT", "eligible_DAY1WT", "eligible_DAY7WT",
           "rma_probe_sets", "mapped_probe_symbol_pairs", "selected_unique_genes",
           "arrays_excluded_after_qc", "model_design_rank", "model_design_columns"),
  value = c(13L, 4L, 5L, 4L, nrow(probe_expr), nrow(map), nrow(gene_expr), 0L,
            qr(design)$rank, ncol(design)),
  note = c(
    "All frozen-eligible CEL files parsed successfully",
    "Independent WT sham biological samples", "Independent WT day-1 biological samples",
    "Independent WT day-7 biological samples", "Affymetrix RMA summaries",
    "After missing-symbol removal and multi-symbol expansion",
    "Maximum IQR probe; lexicographic tie-break", "No technical parse failure; no outcome-based exclusion",
    "No-intercept limma model", "SHAMWT, DAY1WT, DAY7WT"
  ),
  stringsAsFactors = FALSE
)
write.table(qc_summary, file.path(out_dir, "GSE47681_reconstruction_summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

capture.output(sessionInfo(), file = file.path(out_dir, "GSE47681_R_sessionInfo.txt"))
message("Completed GSE47681 RMA/limma reconstruction: ", nrow(gene_expr),
        " unique symbols and ", nrow(effects), " contrast rows")
