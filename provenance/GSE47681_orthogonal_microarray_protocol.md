# GSE47681 orthogonal cross-platform context evaluation protocol

Protocol frozen: 2026-08-03, before download or inspection of expression values.

## Purpose and role

GSE47681 will be used only as an orthogonal cross-platform context evaluation. It will not be treated as a homogeneous fourth discovery cohort and will not enter the existing three-study meta-analysis. Platform, sham baseline and trkB.T1 study context are retained as explicit boundaries.

## Eligibility contract

The target organism is *Mus musculus*. The platform must be GPL1261, Affymetrix Mouse Genome 430 2.0 Array. Only wild-type samples labelled SHAMWT, DAY1WT and DAY7WT are eligible. The expected group sizes from the GEO record are 4, 5 and 4 independent biological replicates, respectively. DAY3 samples and every knockout sample are excluded before expression analysis.

The analysis proceeds only if the official GEO sample metadata reconcile with the raw CEL archive and each included group retains at least two independent biological samples.

## Frozen preprocessing and mapping

Raw CEL files will be normalized together with robust multi-array average background correction, quantile normalization and median-polish summarization. Quality checks will include raw and normalized intensity distributions, relative log expression or an equivalent sample-level diagnostic, and sample clustering/PCA. No sample will be excluded solely because its inclusion weakens the desired conclusion; any exclusion requires a documented technical failure defined before outcome modeling.

Probe sets will be mapped with the platform-specific Bioconductor annotation for GPL1261. Control probes and probe sets without a mouse gene symbol will be removed. When a probe set maps to multiple gene symbols, each symbol will be expanded before gene-level collapse. Multiple probe sets mapping to the same symbol will be collapsed by the probe set with the highest across-sample interquartile range, with lexicographic probe-set ID used as a deterministic tie-breaker. This rule is fixed before contrast inspection.

## Models and contrasts

A no-intercept limma model will include SHAMWT, DAY1WT and DAY7WT. Three contrasts will be estimated with empirical-Bayes moderation: day 1 versus sham, day 7 versus sham and direct day-7-minus-day-1 change. The signed gene score is the moderated t statistic.

The primary Hallmark score will be study-level preranked GSEA normalized enrichment score using the same mouse Hallmark GMT and deterministic ranking/tie handling used by the frozen three-study branch. This maintains the program-scale construction already used for the cross-study comparison; no additional pathway-scoring method will be selected after inspecting results.

## Comparison with the frozen three-study branch

For genes and Hallmarks available in GSE47681 and the frozen three-study matrices, the GSE47681 signed score will be compared with the mean frozen-study score and the frozen-study consensus direction. Reported quantities will include overlap count, Spearman rank correlation, raw direction accuracy, balanced direction accuracy, Matthews correlation coefficient, consensus coverage and conditional accuracy. Endpoint results will be compared with direct change, without pooling GSE47681 into the meta-analysis.

### Operational metric addendum

Frozen on 2026-08-03 after CEL reconstruction and limma estimation but before inspection of any GSE47681-to-reference comparison or Hallmark result.

For each shared feature, the reference score is the arithmetic mean of the three frozen study-specific signed scores (moderated *t* at gene level and NES at Hallmark level). Raw direction accuracy, balanced accuracy and Matthews correlation coefficient compare the sign of the GSE47681 score with the sign of this reference mean. Exact zero scores, if any, are excluded from binary confusion-matrix quantities and their number is reported.

The higher-specificity consensus gate requires all three frozen studies to have the same non-zero sign for a feature. Consensus coverage is the fraction of shared features passing this unanimity gate; conditional accuracy is the fraction of gated features for which GSE47681 has the same sign. These gated metrics are descriptive and are not used to select features or alter the frozen discovery branch. Pairwise GSE47681-to-study correlations are retained as diagnostics, while the reference-mean comparison is primary.

Hallmark GSEA uses 5,000 phenotype-independent rank permutations, seed 20260803, minimum set size 15 and maximum set size 500, with the same mouse Hallmark GMT and deterministic exact-tie handling as the frozen branch.

## Interpretation

Supportive results may be described only as extension of the endpoint/direct-change distinction to an orthogonal cross-platform context. Discordant results will be reported as a platform, baseline, genotype-design or injury-context boundary. Neither outcome will be called an external validation of the original three-study estimand.
