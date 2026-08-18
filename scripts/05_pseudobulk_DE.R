#!/usr/bin/env Rscript
# Sparse pseudobulk differential expression with prespecified direct contrasts.
#
# Important estimand distinction:
# Most SCI time-course datasets have one uninjured baseline rather than time-matched
# sham animals at every post-injury time. A full condition*phase factorial model is
# therefore not identifiable. We fit a single group-factor model and estimate:
#   endpoint       = injured_time - uninjured
#   temporal_delta = injured_time_B - injured_time_A
# The latter is algebraically (B-control) - (A-control), so the shared control cancels
# inside one fitted model. It is a direct change in injury-associated expression, not
# automatically a causal treatment-by-time interaction.

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag, default = NULL) {
  i <- which(args == flag)
  if (length(i) && i < length(args)) args[i + 1] else default
}
hasflag <- function(flag) flag %in% args

validate_coldata <- function(coldata, counts_nrow = NULL) {
  required <- c("sample", "cell_state", "dataset", "subject_id", "condition", "group")
  missing <- setdiff(required, names(coldata))
  if (length(missing)) stop("Missing coldata columns: ", paste(missing, collapse = ", "))
  if (anyDuplicated(rownames(coldata))) stop("Duplicated pseudobulk group_id")
  if (anyNA(coldata[, required, drop = FALSE])) stop("NA in canonical coldata fields")
  if (any(toupper(as.character(coldata$sample)) == "S1")) stop("Dummy sample S1 is prohibited")
  if (!is.null(counts_nrow) && nrow(coldata) != counts_nrow) {
    stop("Counts groups (", counts_nrow, ") != coldata rows (", nrow(coldata), ")")
  }
  invisible(TRUE)
}

if (hasflag("--selftest")) {
  x <- data.frame(sample = c("GSM2", "GSM1"), cell_state = c("Micro", "Micro"),
                  dataset = "GSE_TEST", subject_id = c("A2", "A1"),
                  condition = c("injured", "uninjured"), group = c("dpi_7", "uninjured"),
                  row.names = c("group_B", "group_A"), check.names = FALSE)
  validate_coldata(x, 2)
  a <- 2; b <- -1; control <- 0.5
  stopifnot(all.equal((b - control) - (a - control), b - a))
  bad <- x; bad$sample[1] <- "S1"
  rejected <- try(validate_coldata(bad), silent = TRUE)
  stopifnot(inherits(rejected, "try-error"))
  message("[selftest] PASS: canonical order contract; dummy S1 rejection; shared-control contrast algebra")
  quit(status = 0)
}

dataset <- getarg("--dataset")
if (is.null(dataset)) stop("Provide --dataset GSE...")
pb_dir <- getarg("--pb-dir", "data_processed")
out_dir <- getarg("--out-dir", "results")
contrast_path <- getarg("--contrasts", "tables/analysis_contrasts.tsv")
min_reps <- as.integer(getarg("--min-reps", "2"))
min_cells <- as.integer(getarg("--min-cells", "20"))
fdr_scope <- getarg("--fdr-scope", "dataset_cellstate_axis_effecttype")
if (fdr_scope != "dataset_cellstate_axis_effecttype") {
  stop("Unsupported --fdr-scope; use dataset_cellstate_axis_effecttype")
}

required_packages <- c("Matrix", "limma", "edgeR")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "),
       ". Install from Bioconductor before inferential execution.")
}

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
count_path <- file.path(pb_dir, paste0(dataset, "_pseudobulk_counts.mtx.gz"))
gene_path <- file.path(pb_dir, paste0(dataset, "_pseudobulk_genes.tsv"))
coldata_path <- file.path(pb_dir, paste0(dataset, "_pseudobulk_coldata.tsv"))
for (path in c(count_path, gene_path, coldata_path, contrast_path)) {
  if (!file.exists(path)) stop("Required input missing: ", path)
}

con <- gzfile(count_path, open = "rb")
on.exit(close(con), add = TRUE)
pb <- Matrix::readMM(con)  # groups x genes
close(con)
on.exit(NULL, add = FALSE)
genes <- read.delim(gene_path, check.names = FALSE, stringsAsFactors = FALSE)
coldata <- read.delim(coldata_path, row.names = 1, check.names = FALSE,
                      stringsAsFactors = FALSE)
validate_coldata(coldata, nrow(pb))
if (nrow(genes) != ncol(pb)) stop("Gene rows do not match count matrix columns")
if (anyDuplicated(genes$gene)) stop("Duplicated gene names in pseudobulk gene table")

counts <- as(pb, "dgCMatrix")
counts <- Matrix::t(counts)  # genes x groups
rownames(counts) <- genes$gene
colnames(counts) <- rownames(coldata)  # positional contract, validated above
if (!identical(colnames(counts), rownames(coldata))) stop("Internal group order mismatch")
if (any(coldata$n_cells < min_cells)) {
  stop("Pseudobulk group below --min-cells detected; regenerate or lower threshold explicitly")
}

contrast_spec <- read.delim(contrast_path, check.names = FALSE, stringsAsFactors = FALSE)
contrast_spec <- contrast_spec[contrast_spec$dataset == dataset, , drop = FALSE]
if (!nrow(contrast_spec)) stop("No prespecified contrasts for ", dataset)
if (anyDuplicated(contrast_spec$contrast_id)) stop("Duplicated contrast_id for ", dataset)

results <- list()
diagnostics <- list()
result_i <- 0L
diag_i <- 0L

for (cell_state in unique(coldata$cell_state)) {
  in_state <- which(coldata$cell_state == cell_state)
  state_meta <- coldata[in_state, , drop = FALSE]
  available <- table(state_meta$group)
  eligible_spec <- contrast_spec[
    contrast_spec$numerator_group %in% names(available) &
      contrast_spec$denominator_group %in% names(available), , drop = FALSE]
  if (nrow(eligible_spec)) {
    enough <- vapply(seq_len(nrow(eligible_spec)), function(i) {
      available[[eligible_spec$numerator_group[i]]] >= min_reps &&
        available[[eligible_spec$denominator_group[i]]] >= min_reps
    }, logical(1))
    eligible_spec <- eligible_spec[enough, , drop = FALSE]
  }
  if (!nrow(eligible_spec)) {
    diag_i <- diag_i + 1L
    diagnostics[[diag_i]] <- data.frame(
      dataset = dataset, cell_state = cell_state, status = "skipped_no_estimable_contrast",
      n_groups = length(in_state), group_counts = paste(names(available), available, sep = ":", collapse = ";"),
      design_rank = NA_integer_, design_columns = NA_integer_, genes_tested = 0L)
    next
  }

  needed_groups <- unique(c(eligible_spec$numerator_group, eligible_spec$denominator_group))
  use <- in_state[state_meta$group %in% needed_groups]
  meta <- coldata[use, , drop = FALSE]
  y <- counts[, use, drop = FALSE]
  group_factor <- factor(meta$group, levels = unique(needed_groups))
  design <- model.matrix(~ 0 + group_factor)
  colnames(design) <- make.names(levels(group_factor))
  if (qr(design)$rank < ncol(design)) stop("Rank-deficient design in ", cell_state)

  dge <- edgeR::DGEList(counts = y)
  keep_gene <- edgeR::filterByExpr(dge, design = design)
  dge <- dge[keep_gene, , keep.lib.sizes = FALSE]
  if (nrow(dge) < 50) {
    diag_i <- diag_i + 1L
    diagnostics[[diag_i]] <- data.frame(
      dataset = dataset, cell_state = cell_state, status = "skipped_too_few_genes",
      n_groups = nrow(meta), group_counts = paste(names(available), available, sep = ":", collapse = ";"),
      design_rank = qr(design)$rank, design_columns = ncol(design), genes_tested = nrow(dge))
    next
  }
  dge <- edgeR::calcNormFactors(dge, method = "TMM")
  voomed <- limma::voom(dge, design = design, plot = FALSE)
  fit <- limma::lmFit(voomed, design)
  contrast_strings <- paste0(make.names(eligible_spec$numerator_group), "-",
                             make.names(eligible_spec$denominator_group))
  contrast_matrix <- limma::makeContrasts(contrasts = contrast_strings, levels = design)
  colnames(contrast_matrix) <- eligible_spec$contrast_id
  fit2 <- limma::eBayes(limma::contrasts.fit(fit, contrast_matrix), robust = TRUE)

  for (j in seq_len(nrow(eligible_spec))) {
    spec <- eligible_spec[j, , drop = FALSE]
    estimate <- fit2$coefficients[, j]
    se <- sqrt(fit2$s2.post) * fit2$stdev.unscaled[, j]
    p_value <- fit2$p.value[, j]
    df_total <- fit2$df.total
    crit <- qt(0.975, df = df_total)
    result_i <- result_i + 1L
    results[[result_i]] <- data.frame(
      gene = rownames(fit2$coefficients), cell_state = cell_state, dataset = dataset,
      axis = spec$axis, contrast_id = spec$contrast_id, effect_type = spec$effect_type,
      numerator_group = spec$numerator_group, denominator_group = spec$denominator_group,
      primary = spec$primary, estimate = estimate, se = se, df_residual = df_total,
      ci_low = estimate - crit * se, ci_high = estimate + crit * se,
      p = p_value, n_num = unname(available[[spec$numerator_group]]),
      n_den = unname(available[[spec$denominator_group]]), stringsAsFactors = FALSE,
      check.names = FALSE)
  }
  diag_i <- diag_i + 1L
  diagnostics[[diag_i]] <- data.frame(
    dataset = dataset, cell_state = cell_state, status = "modelled",
    n_groups = nrow(meta), group_counts = paste(names(available), available, sep = ":", collapse = ";"),
    design_rank = qr(design)$rank, design_columns = ncol(design), genes_tested = nrow(dge))
}

diagnostic_out <- do.call(rbind, diagnostics)
write.table(diagnostic_out, file.path(out_dir, paste0("model_diagnostics_", dataset, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (!length(results)) {
  message("No cell state passed replicate/estimability gates for ", dataset)
  quit(status = 0)
}

out <- do.call(rbind, results)
family <- interaction(out$dataset, out$cell_state, out$axis, out$effect_type, drop = TRUE, lex.order = TRUE)
out$fdr <- ave(out$p, family, FUN = function(x) p.adjust(x, method = "BH"))
out <- out[, c("gene", "cell_state", "dataset", "axis", "contrast_id", "effect_type",
               "numerator_group", "denominator_group", "primary", "estimate", "se",
               "df_residual", "ci_low", "ci_high", "p", "fdr", "n_num", "n_den")]
outfile <- file.path(out_dir, paste0("effects_", dataset, ".tsv"))
write.table(out, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote ", outfile, " with ", nrow(out), " rows. Temporal deltas are direct within-model contrasts; ",
        "they are not labelled causal condition-by-time interactions.")
