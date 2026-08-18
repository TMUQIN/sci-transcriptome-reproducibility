# GSE47681 orthogonal cross-platform context evaluation

## Decision

GSE47681 is retained as an orthogonal cross-platform context evaluation. It is not a fourth discovery cohort and must not enter the frozen three-study meta-analysis.

The result strengthens the manuscript’s positive claim at the Hallmark scale: endpoint transfer is high across an independent Affymetrix bulk-microarray context, whereas direct day-7-minus-day-1 transfer is materially weaker. At gene level, all three contrasts show moderate transfer and the direct contrast is not worse than the endpoints. The correct interpretation is therefore a **scale-by-estimand distinction**, not a universal claim that direct temporal change is always less reproducible at every molecular level.

## Raw-data reconstruction

- Official GEO accession: GSE47681; GPL1261 Affymetrix Mouse Genome 430 2.0 Array.
- Context: moderate thoracic contusion; lesion-site tissue; trkB.T1 WT/KO design.
- Eligible samples fixed before expression inspection: WT sham n=4, WT day 1 n=5, WT day 7 n=4.
- Excluded by design: all KO arrays and all day-3 arrays.
- Thirteen sample-level CEL.gz files were downloaded from NCBI GEO and passed gzip integrity and CEL-header checks.
- All 13 arrays were normalized together by RMA background correction, quantile normalization and median-polish summarization.
- All CEL files parsed successfully; no array was excluded after QC.
- RMA yielded 45,101 probe-set summaries. GPL1261 annotation yielded 40,713 valid probe–symbol pairs after multi-symbol expansion. The fixed maximum-IQR/lexicographic rule selected one probe for each of 22,007 unique mouse gene symbols.
- The no-intercept robust limma design had rank 3/3 and generated three direct contrasts (66,021 gene-contrast rows).

## Quality control

- RLE medians ranged from -0.00032 to 0.00904 and RLE IQRs from 0.133 to 0.226.
- The lowest pairwise normalized probe-set correlation was 0.967.
- PCA of the 5,000 most variable probe sets separated sham, day 1 and day 7: PC1 explained 58.8% and PC2 27.7% of variance.
- No technical parse failure, gross array-level distribution failure or isolated correlation failure was observed. All eligible arrays were retained.

## Frozen-reference comparison

The reference score is the mean of the three frozen study scores. Binary metrics compare GSE47681 with the sign of that mean. The high-specificity gate requires all three frozen studies to share a non-zero sign.

| Scale | Contrast | Shared features | Spearman vs frozen mean | Direction accuracy | Balanced accuracy | MCC | Frozen unanimity coverage | Accuracy given unanimity |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| Gene | day 1 | 13,749 | 0.571 | 0.705 | 0.709 | 0.419 | 0.483 | 0.815 |
| Gene | day 7 | 13,749 | 0.550 | 0.679 | 0.685 | 0.371 | 0.485 | 0.787 |
| Gene | direct change | 13,749 | 0.594 | 0.704 | 0.705 | 0.411 | 0.452 | 0.811 |
| Hallmark | day 1 | 50 | 0.879 | 0.880 | 0.821 | 0.475 | 0.820 | 0.927 |
| Hallmark | day 7 | 50 | 0.812 | 0.940 | 0.969 | 0.612 | 0.820 | 0.976 |
| Hallmark | direct change | 50 | 0.649 | 0.640 | 0.644 | 0.294 | 0.360 | 0.778 |

No exact-zero score was excluded from any binary metric.

## Interpretation for the manuscript

Recommended main-text statement:

> In an independently reconstructed WT bulk-microarray cohort, Hallmark transfer remained high for the day-1 and day-7 endpoints (Spearman ρ=0.879 and 0.812; direction accuracy=0.88 and 0.94) but was lower for direct day-7-minus-day-1 change (ρ=0.649; accuracy=0.64). Among Hallmarks with unanimous direction in the three discovery cohorts, GSE47681 retained 92.7% and 97.6% of endpoint directions versus 77.8% of direct-change directions. Gene-level transfer was moderate and similar across the three contrasts, localizing the strongest endpoint/direct-change separation to the program scale.

Required boundary sentence:

> Because GSE47681 differs in platform, sham baseline and trkB.T1 study context, it was used as an orthogonal cross-platform context evaluation and was not pooled with the three discovery cohorts.

Avoid the terms “independent validation cohort,” “replication cohort,” or “four-study meta-analysis.”

## Reproducibility assets

- Frozen protocol and operational addendum: `provenance/GSE47681_orthogonal_microarray_protocol.md`
- Raw RMA/limma script: `scripts/68_gse47681_raw_rma_limma.R`
- GSEA and comparison script: `scripts/69_gse47681_hallmark_stress_test.py`
- Sample manifest, QC tables, design/contrast matrices, selected-probe dictionary, gene-level effects, Hallmark results, comparison metrics and SHA256 provenance are stored in this directory.
