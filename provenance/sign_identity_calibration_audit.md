# Post-lock direction-imbalance and feature-identity calibration audit

Decision: **B_prevalence_qualified**.

This sensitivity branch used only the frozen three-study signed-score matrices. It did not recompute or alter any v5.3.2 primary model, threshold, dataset decision or inferential result.

## Hallmark-scale results

For `injury_d1_vs_uninjured`, all-study sign concordance was 0.820 (identity-permutation P=0.0002), mean pairwise kappa was 0.354, held-out balanced accuracy was 0.790, MCC was 0.502, and agreement-gated conditional accuracy was 0.934.

For `injury_d7_vs_uninjured`, all-study sign concordance was 0.820 (identity-permutation P=0.009199), mean pairwise kappa was 0.259, held-out balanced accuracy was 0.748, MCC was 0.359, and agreement-gated conditional accuracy was 0.933.

For `change_d7_minus_d1`, all-study sign concordance was 0.360 (identity-permutation P=0.0003), mean pairwise kappa was 0.222, held-out balanced accuracy was 0.662, MCC was 0.295, and agreement-gated conditional accuracy was 0.659.

## Direction prevalence

The study-specific positive and negative proportions are reported in `study_direction_prevalence.tsv`. They define the global sign-prevalence structure preserved by every permutation.

## Prespecified interpretation

{
  "direct_identity_metric_p_le_005_count": 4,
  "endpoints": {
    "injury_d1_vs_uninjured": {
      "identity_metric_p_le_005_count": 4,
      "imbalance_metric_higher_than_direct_count": 4,
      "raw_all_study_concordance_higher_than_direct": true,
      "all_study_concordance_empirical_p_le_005": true
    },
    "injury_d7_vs_uninjured": {
      "identity_metric_p_le_005_count": 4,
      "imbalance_metric_higher_than_direct_count": 4,
      "raw_all_study_concordance_higher_than_direct": true,
      "all_study_concordance_empirical_p_le_005": true
    }
  }
}

The feature-identity null tests whether the same feature transfers across studies beyond marginal sign prevalence. It is not a correlation-matched competitive gene-set null and cannot establish that curated Hallmarks outperform matched synthetic sets.

All observed and null summaries are retained irrespective of the decision category.
