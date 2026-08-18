# Analysis workflow

This document describes the analytical pipeline that produced the frozen
results in `data/` and the figures in the manuscript, and maps each stage to its
authoritative script in `scripts/`.

## Stage 1 — Sample metadata & QC
- `00_build_sample_metadata.py` — auditable sample-level metadata from frozen GEO check records.
- `03b_qc_pseudobulk.py` — dataset-aware QC, canonical sample crosswalk, sparse pseudobulk.
- `03c_harmonize_pseudobulk.py` — aggregate author-state pseudobulk counts before DE.

Inferential unit: **sample-level pseudobulk**. Cells/nuclei are never independent replicates.

## Stage 2 — Pseudobulk differential expression
- `05_pseudobulk_DE.R` — no-intercept group design; filterByExpr; TMM; voom; lmFit;
  robust eBayes (edgeR/limma). Produces `results/whole_lesion/effects_GSE*.tsv`.

## Stage 3 — Gene-level & Hallmark scoring / GSEA
- `22_sample_level_pathway_scores.py` — magnitude-aware and rank-based Hallmark scores.
- `06e_hallmark_gsea.py` — secondary Hallmark GSEA for the whole-lesion branch.
- `22_pathway_score_DE.R` — pathway-score differential expression.

## Stage 4 — Cross-study meta-analysis & calibration
- `06_meta_analysis.py` — comparability-aware REML random-effects meta-analysis.
- `23_pathway_effect_meta.py` — REML + modified Hartung-Knapp meta of pathway effects.
- `25_gsema_mkh_meta.py` — REML+mKH estimator on GSEMA study-level effects.
- `70_effect_strength_calibration.py` — effect-strength calibration.

## Stage 5 — Sensitivity & context evaluations
- `28_hallmark_myc_stability_audit.py` — Hallmark stability & MYC/cell-cycle audit.
- `29_composition_ontology_sensitivity.py` — composition/ontology sensitivity.
- `30_mtorc1_context_matrix.py` — mTORC1 context/method matrix.
- `26_matched_gene_set_null_pilot.py` — matched random gene-set pilot.
- `17_audit_geo_candidates.py` — audit of five recent SCI single-cell GEO candidates.
- `65_audit_gse205029_temporal_eligibility.py` — GSE205029 temporal estimability (no counts read).

## Stage 6 — Feature-identity permutation
- `67_postlock_sign_identity_calibration.py` — 10,000 within-study feature-label
  permutations preserving score/rank/sign distributions; seed 20260803.
  Outputs `data/sign_identity_permutation_null_draws.npz` (+ summary).

## Stage 7 — Orthogonal GSE47681 context
- `68_gse47681_raw_rma_limma.R` — joint RMA + robust limma on 13 WT CEL files.
- `69_gse47681_hallmark_stress_test.py` — Hallmark GSEA & frozen-reference comparison.

## Stage 8 - Figure assembly
- The current manuscript-facing organization is Fig. 1 analytical sample
  reconstruction, Fig. 2 gene/Hallmark transfer, Fig. 3 effect-strength
  sensitivity, Fig. 4 method-specific Hallmark support, Fig. 5 focal-program
  context-dependent transfer profiles, and Fig. S1-S6 as documented in
  `docs/figure_reproduction_map.md`.
- Final TIFF presentation assembly used project-level plotting infrastructure
  that is not included in this minimal public release. The former
  `72_build_v5_3_3_figures.py` is retained under `scripts/legacy/` for
  provenance only and is not a supported command.

## Companion provenance
`provenance/` contains the protocols and JSON provenance records (hashes, seeds,
software) that anchor every stage above. See `metadata/CODE_TO_OUTPUT_MAP.csv`
for the figure/table → data → script mapping.
