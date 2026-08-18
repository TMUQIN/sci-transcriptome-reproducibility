#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(limma))

args <- commandArgs(trailingOnly=TRUE)
input <- if (length(args) >= 1) args[[1]] else "reports/phase_reproducibility_calibration_2026_07/sample_level_hallmark_scores.tsv"
output <- if (length(args) >= 2) args[[2]] else "reports/phase_reproducibility_calibration_2026_07/pathway_effects_by_study.tsv"

scores <- read.delim(input, check.names=FALSE, stringsAsFactors=FALSE)
required <- c("dataset", "sample", "subject_id", "group", "program", "method", "standardized_score", "n_genes_observed")
if (!all(required %in% colnames(scores))) stop("missing required score columns")
if (any(duplicated(scores[, c("dataset", "subject_id", "program", "method")]))) stop("duplicate subject/program/method rows")

contrast_defs <- list(
  injury_d1_vs_uninjured=c(uninjured=-1, dpi_1=1, dpi_7=0),
  injury_d7_vs_uninjured=c(uninjured=-1, dpi_1=0, dpi_7=1),
  change_d7_minus_d1=c(uninjured=0, dpi_1=-1, dpi_7=1)
)
effect_types <- c(injury_d1_vs_uninjured="endpoint", injury_d7_vs_uninjured="endpoint", change_d7_minus_d1="temporal_delta")
rows <- list()
ri <- 1L

for (dataset in sort(unique(scores$dataset))) {
  for (method in sort(unique(scores$method))) {
    x <- scores[scores$dataset == dataset & scores$method == method, ]
    meta <- unique(x[, c("sample", "subject_id", "group")])
    meta <- meta[match(unique(x$subject_id), meta$subject_id), ]
    if (!all(c("uninjured", "dpi_1", "dpi_7") %in% meta$group)) stop(paste(dataset, method, "missing required group"))
    if (any(table(meta$group) < 2)) stop(paste(dataset, method, "insufficient biological replicates"))
    wide <- reshape(x[, c("program", "subject_id", "standardized_score")], idvar="program", timevar="subject_id", direction="wide")
    rownames(wide) <- wide$program
    wide$program <- NULL
    colnames(wide) <- sub("^standardized_score\\.", "", colnames(wide))
    meta <- meta[match(colnames(wide), meta$subject_id), ]
    if (any(is.na(meta$subject_id))) stop("score/metadata alignment failure")
    mat <- as.matrix(wide)
    storage.mode(mat) <- "double"
    group <- factor(meta$group, levels=c("uninjured", "dpi_1", "dpi_7"))
    design <- model.matrix(~0 + group)
    colnames(design) <- levels(group)
    fit <- lmFit(mat, design)
    cm <- do.call(cbind, contrast_defs)
    fit2 <- contrasts.fit(fit, cm)
    fit2 <- eBayes(fit2, robust=TRUE)
    for (contrast in colnames(cm)) {
      j <- match(contrast, colnames(fit2$coefficients))
      estimate <- fit2$coefficients[, j]
      se <- fit2$stdev.unscaled[, j] * sqrt(fit2$s2.post)
      df <- fit2$df.total
      p <- fit2$p.value[, j]
      fdr <- p.adjust(p, method="BH")
      crit <- qt(0.975, df=df)
      for (program in rownames(mat)) {
        k <- match(program, rownames(mat))
        ng <- unique(x$n_genes_observed[x$program == program])
        rows[[ri]] <- data.frame(
          program=program, dataset=dataset, method=method, contrast_id=contrast,
          effect_type=unname(effect_types[[contrast]]), estimate=estimate[[k]], se=se[[k]],
          df_residual=df[[k]], ci_low=estimate[[k]]-crit[[k]]*se[[k]],
          ci_high=estimate[[k]]+crit[[k]]*se[[k]], p=p[[k]], fdr=fdr[[k]],
          n_num=if (contrast == "injury_d1_vs_uninjured") sum(group == "dpi_1") else sum(group == "dpi_7"),
          n_den=if (contrast == "injury_d7_vs_uninjured") sum(group == "uninjured") else if (contrast == "injury_d1_vs_uninjured") sum(group == "uninjured") else sum(group == "dpi_1"),
          n_genes_observed=ng[[1]], stringsAsFactors=FALSE
        )
        ri <- ri + 1L
      }
    }
  }
}
result <- do.call(rbind, rows)
result <- result[order(result$method, result$contrast_id, result$program, result$dataset), ]
write.table(result, output, sep="\t", row.names=FALSE, quote=FALSE, na="NA", eol="\n")
cat(sprintf("[pathway-score-DE] PASS rows=%d datasets=%d methods=%d programs=%d\n", nrow(result), length(unique(result$dataset)), length(unique(result$method)), length(unique(result$program))))

