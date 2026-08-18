# Current figure and table reproduction map

This map is synchronized to the 2026-08-18 manuscript-facing figure
organization. It replaces the obsolete v5.3.3 map. Source data and upstream
analysis scripts are archived where available; final TIFF assembly is not
claimed clean-clone reproducible from this minimal release.

## Main figures

| Figure | Manuscript role | Archived source data | Upstream analysis | Presentation assembly | Status |
|---|---|---|---|---|---|
| Fig. 1 | Analytical sample reconstruction and estimability | `tables/Supplementary_Table_S1_Analytical_Samples.csv`; `data/estimability_matrix.csv`; `data/GSE304361_program_support_summary.csv` | `00_build_sample_metadata.py`; `03b_qc_pseudobulk.py`; `21_reproducibility_metric_calibration.py`; `17_audit_geo_candidates.py` | Historical project-level plotting infrastructure | Numerical source archived; presentation assembly not clean-clone runnable |
| Fig. 2 | Gene/Hallmark transfer | `data/sign_identity_observed_metrics.csv`; `data/sign_identity_permutation_summary.csv`; `data/GSE47681_frozen_reference_comparison.csv` | `67_postlock_sign_identity_calibration.py`; `69_gse47681_hallmark_stress_test.py` | Historical project-level plotting infrastructure | Numerical source archived; presentation assembly not clean-clone runnable |
| Fig. 3 | Effect-strength sensitivity | `data/effect_strength_summary.csv`; `data/effect_strength_stratified_metrics.csv`; `data/effect_strength_overlap_window_metrics.csv` | `70_effect_strength_calibration.py` | Historical project-level plotting infrastructure | Numerical source archived; presentation assembly not clean-clone runnable |
| Fig. 4 | Method-specific Hallmark support | `data/pathway_meta_mkh.csv`; `data/gsema_mkh_random_effects.csv`; `data/hallmark_gsea_by_study.csv` | `23_pathway_effect_meta.py`; `25_gsema_mkh_meta.py`; `06e_hallmark_gsea.py`; `24_gsema_benchmark.R` | Historical project-level plotting infrastructure | Numerical source archived; presentation assembly not clean-clone runnable |
| Fig. 5 | Focal-program context-dependent transfer profiles | `data/GSE304361_program_support_summary.csv`; `data/composition_program_association.csv`; `data/hallmark_gsea_by_study.csv`; `data/pathway_meta_mkh.csv`; `data/myc_deoverlap_results.csv`; `data/mtorc1_context_matrix.csv` | `29_composition_ontology_sensitivity.py`; `28_hallmark_myc_stability_audit.py`; `30_mtorc1_context_matrix.py`; `06e_hallmark_gsea.py`; `23_pathway_effect_meta.py` | Historical project-level plotting infrastructure | Numerical source archived; presentation assembly not clean-clone runnable |

## Supplementary figures

| Figure | Manuscript role | Archived source data | Upstream analysis | Status |
|---|---|---|---|---|
| Fig. S1 | Composition sensitivity | `data/composition_program_association.csv`; `data/mtorc1_context_matrix.csv` | `29_composition_ontology_sensitivity.py`; `30_mtorc1_context_matrix.py` | Source data archived; presentation assembly not clean-clone runnable |
| Fig. S2 | Estimability | `data/estimability_flow_matrix.csv`; `data/estimability_matrix.csv`; `data/GSE205029_contrast_estimability.csv` | `65_audit_gse205029_temporal_eligibility.py`; `21_reproducibility_metric_calibration.py` | Source data archived; presentation assembly not clean-clone runnable |
| Fig. S3 | Equal-feature-count diagnostic | `data/equal_feature_count_diagnostic.csv` | `21_reproducibility_metric_calibration.py` | Source data archived; presentation assembly not clean-clone runnable |
| Fig. S4 | GSE47681 QC | `data/GSE47681_RMA_PCA_scores.csv`; `data/GSE47681_RMA_PCA_variance.csv`; `data/GSE47681_RLE_diagnostics.csv`; `data/GSE47681_RMA_sample_correlations.csv` | `68_gse47681_raw_rma_limma.R` | Reproducible with public CEL inputs; final assembly not clean-clone runnable |
| Fig. S5 | Conditional member stability | `data/conditional_member_stability.csv`; `data/leading_edge_pairwise_calibrated.csv`; `data/myc_deoverlap_results.csv` | `06f_leading_edge_replicability.py`; `28_hallmark_myc_stability_audit.py` | Source data archived; presentation assembly not clean-clone runnable |
| Fig. S6 | Matched-set exchangeability | `data/matched_null_exchangeability_diagnostics.csv`; `data/matched_null_expansion_decision.csv` | `26_matched_gene_set_null_pilot.py` | Source data archived; presentation assembly not clean-clone runnable |

## Supplementary tables

| Table | Source | Synchronization status |
|---|---|---|
| S1 | `tables/Supplementary_Table_S1_Analytical_Samples.csv` | Byte-synchronized with current ESM2 archive |
| S2 | `tables/Supplementary_Table_S2_Dataset_Disposition.csv` | Byte-synchronized with current ESM2 archive |
| S3 | `tables/Supplementary_Table_S3_Reconstruction_Parameters.csv` | Byte-synchronized with current ESM2 archive |

GSE47681 is consistently described as an orthogonal cross-platform context
evaluation. It is not pooled into the primary three-dataset synthesis and is
not described as formal external validation.
