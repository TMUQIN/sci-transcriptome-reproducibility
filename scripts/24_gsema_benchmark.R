#!/usr/bin/env Rscript
.libPaths(c(normalizePath(".deps/R"), .libPaths()))
suppressPackageStartupMessages({
  library(Matrix)
  library(GSEMA)
})

args <- commandArgs(trailingOnly=TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else "reports/phase_reproducibility_calibration_2026_07/gsema_benchmark"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

studies <- c("GSE162610", "GSE234774", "GSE304399")
contrasts <- list(
  injury_d1_vs_uninjured=c(exp="dpi_1", ref="uninjured"),
  injury_d7_vs_uninjured=c(exp="dpi_7", ref="uninjured"),
  change_d7_minus_d1=c(exp="dpi_7", ref="dpi_1")
)
methods <- c("Zscore", "Singscore", "GSVA", "ssGSEA")

read_gmt <- function(path) {
  lines <- readLines(path, warn=FALSE)
  fields <- strsplit(lines, "\t", fixed=TRUE)
  result <- lapply(fields, function(x) unique(x[-c(1, 2)]))
  names(result) <- vapply(fields, `[[`, character(1), 1)
  result
}

load_study <- function(study) {
  base <- file.path("data_processed", "whole_lesion")
  con <- gzfile(file.path(base, paste0(study, "_pseudobulk_counts.mtx.gz")), "rb")
  on.exit(close(con), add=TRUE)
  counts <- readMM(con)
  genes <- read.delim(file.path(base, paste0(study, "_pseudobulk_genes.tsv")), stringsAsFactors=FALSE)$gene
  meta <- read.delim(file.path(base, paste0(study, "_pseudobulk_coldata.tsv")), stringsAsFactors=FALSE)
  if (nrow(counts) != nrow(meta) || ncol(counts) != length(genes)) stop(paste(study, "dimension mismatch"))
  if (anyDuplicated(meta$subject_id)) stop(paste(study, "duplicate biological subject"))
  lib <- Matrix::rowSums(counts)
  logcpm <- log2(sweep(as.matrix(counts), 1, lib, "/") * 1e6 + 0.5)
  rownames(logcpm) <- meta$subject_id
  colnames(logcpm) <- genes
  list(logcpm=logcpm, meta=meta)
}

gene_sets <- read_gmt(file.path("references", "msigdb_mh.all.v2026.1.Mm.symbols.gmt"))
if (length(gene_sets) != 50) stop("expected 50 Hallmarks")
loaded <- lapply(studies, load_study)
names(loaded) <- studies

native_rows <- list()
effect_rows <- list()
status_rows <- list()
ni <- 1L; ei <- 1L; si <- 1L

for (contrast_id in names(contrasts)) {
  groups <- contrasts[[contrast_id]]
  list_ex <- list(); list_pheno <- list()
  for (study in studies) {
    x <- loaded[[study]]
    keep <- x$meta$group %in% groups
    meta <- x$meta[keep, , drop=FALSE]
    expr <- t(x$logcpm[keep, , drop=FALSE])
    colnames(expr) <- meta$subject_id
    pheno <- data.frame(Condition=meta$group, row.names=meta$subject_id, stringsAsFactors=FALSE)
    if (any(table(pheno$Condition) < 2)) stop(paste(study, contrast_id, "insufficient biological replicates"))
    list_ex[[study]] <- expr
    list_pheno[[study]] <- pheno
  }
  for (method in methods) {
    started <- Sys.time()
    run <- tryCatch({
      object <- createObjectMApath(
        listEX=list_ex, listPheno=list_pheno,
        namePheno=rep("Condition", length(studies)),
        expGroups=rep(list(unname(groups[["exp"]])), length(studies)),
        refGroups=rep(list(unname(groups[["ref"]])), length(studies)),
        geneSets=gene_sets, pathMethod=method, minSize=10,
        kcdf="Gaussian", normalize=TRUE, n.cores=1, internal.n.cores=1
      )
      effects <- calculateESpath(object, measure="limma", WithinVarCorrect=TRUE, missAllow=0.3)
      native <- metaAnalysisESpath(objectMApath=object, measure="limma", WithinVarCorrect=TRUE,
                                   typeMethod="REM", missAllow=0.3, numData=length(studies))
      native$method <- method
      native$contrast_id <- contrast_id
      native_rows[[ni]] <- native; ni <- ni + 1L
      common <- intersect(rownames(effects$ES), rownames(effects$Var))
      for (pathway in common) {
        for (j in seq_along(studies)) {
          effect_rows[[ei]] <- data.frame(
            method=method, contrast_id=contrast_id, pathway=pathway, dataset=studies[[j]],
            estimate=effects$ES[pathway, j], variance=effects$Var[pathway, j],
            se=sqrt(effects$Var[pathway, j]), stringsAsFactors=FALSE
          )
          ei <- ei + 1L
        }
      }
      list(ok=TRUE, n_pathways=nrow(native), error="")
    }, error=function(e) list(ok=FALSE, n_pathways=0L, error=conditionMessage(e)))
    status_rows[[si]] <- data.frame(
      method=method, contrast_id=contrast_id, status=if (run$ok) "completed" else "failed",
      n_pathways=run$n_pathways, elapsed_seconds=as.numeric(difftime(Sys.time(), started, units="secs")),
      error=run$error, stringsAsFactors=FALSE
    )
    si <- si + 1L
    write.table(do.call(rbind, status_rows), file.path(out_dir, "gsema_run_status.tsv"),
                sep="\t", row.names=FALSE, quote=FALSE, na="NA", eol="\n")
  }
}

if (length(native_rows)) {
  native <- do.call(rbind, native_rows)
  native <- native[order(native$method, native$contrast_id, native$FDR, native$Pathway), ]
  write.table(native, file.path(out_dir, "gsema_native_random_effects.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE, na="NA", eol="\n")
}
if (length(effect_rows)) {
  effects <- do.call(rbind, effect_rows)
  effects <- effects[order(effects$method, effects$contrast_id, effects$pathway, effects$dataset), ]
  write.table(effects, file.path(out_dir, "gsema_effects_by_study.tsv"),
              sep="\t", row.names=FALSE, quote=FALSE, na="NA", eol="\n")
}
session <- capture.output(sessionInfo())
writeLines(session, file.path(out_dir, "sessionInfo.txt"), useBytes=TRUE)
cat(sprintf("[GSEMA benchmark] completed=%d failed=%d\n",
            sum(vapply(status_rows, function(x) x$status == "completed", logical(1))),
            sum(vapply(status_rows, function(x) x$status == "failed", logical(1)))))
