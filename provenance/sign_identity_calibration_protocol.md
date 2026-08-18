# Post-lock direction-imbalance and feature-identity calibration protocol

Protocol frozen: 2026-08-03, before execution of the analyses defined below.

## Scientific status

This is an independent post-lock sensitivity branch. It does not replace or recompute the frozen differential-expression, GSEA, sample-score, meta-analysis, matched-null, composition, leading-edge or context analyses used in manuscript v5.3.2. The v5.3.2 package is retained unchanged until this branch has been audited.

## Inputs

The branch uses the already frozen complete-feature score matrices from GSE162610, GSE234774 and GSE304399. Gene scores are moderated statistics (`estimate / se`) for the 15,331 genes with complete signed effects in all three studies. Program scores are normalized enrichment scores for the 50 mouse Hallmarks with complete results in all three studies. The three contrasts are day 1 versus uninjured, day 7 versus uninjured and direct day-7-minus-day-1 change.

Gene and Hallmark scores are evaluated within scale. Their numerical magnitudes are not treated as directly exchangeable.

## Observed metrics

For every scale and contrast, the following will be reported:

1. Study-specific positive, negative and exact-zero feature counts and proportions.
2. All-study nonzero sign-concordance rate.
3. Mean pairwise raw sign concordance.
4. Mean pairwise Spearman correlation of signed scores.
5. Mean leave-one-study-out direction accuracy, using the sign of the mean signed score in the other two studies.
6. Mean pairwise Cohen's kappa for binary positive/negative directions.
7. Fleiss' kappa across the three studies as raters and positive/negative direction as categories.
8. Mean leave-one-study-out balanced accuracy and Matthews correlation coefficient using the same prediction rule as item 5.
9. Mean leave-one-study-out agreement-gated coverage: the proportion of features for which the two training studies have the same nonzero direction.
10. Mean conditional held-out accuracy among the covered features.

Exact-zero signed scores will be excluded only from a metric when its binary definition requires positive/negative categories. Their counts will be disclosed. When an MCC denominator is zero because a predictor is constant, MCC is recorded as zero rather than missing, following the conventional finite-value definition for a non-informative constant classifier. No significance or false-discovery threshold is used to define direction.

## Feature-identity permutation null

The random seed is 20260803 and the number of permutations is 10,000. Within each scale and contrast, feature labels will be independently permuted within each study while the study-specific signed-score vector itself remains unchanged. This preserves sample size, score distribution, ranks and marginal positive/negative prevalence in every study while breaking the correspondence of gene or Hallmark identity across studies.

All metrics above will be recomputed for every permutation. For metrics where larger values indicate stronger transfer, the one-sided empirical probability is `(1 + count(null >= observed)) / 10001`. Null median and 2.5th and 97.5th percentiles will also be reported. Coverage will be treated primarily as a descriptive availability quantity even though its null distribution is retained for completeness.

This null tests identity-specific transfer beyond global sign-prevalence agreement. It is not a correlation-matched random gene-set null, does not compare curated Hallmarks with matched synthetic gene sets and cannot establish a curated-gene-set advantage.

## Decision rule for manuscript integration

The audit will assign one of three prespecified interpretations. The formal decision is applied to the Hallmark scale because the reviewer concern targets the high Hallmark endpoint rates; gene-scale results are reported as contextual calibration.

- **A, strengthened:** for each endpoint, all-study concordance has empirical P <= 0.05; at least two of mean pairwise Spearman correlation, balanced accuracy, Matthews correlation coefficient and conditional accuracy have empirical P <= 0.05; and at least three of mean pairwise Cohen's kappa, balanced accuracy, Matthews correlation coefficient and conditional accuracy exceed the corresponding direct-change value. Direct change must not meet the same identity-support rule.
- **C, weakened:** either endpoint has raw all-study concordance no greater than direct change, or direct change is at least as high as that endpoint for at least three of mean pairwise Cohen's kappa, balanced accuracy, Matthews correlation coefficient and conditional accuracy.
- **B, prevalence-qualified:** every outcome that is neither A nor C. The manuscript must explicitly distinguish the supported identity-specific components from global sign-prevalence agreement.

All observed and null metrics will be released regardless of outcome. No threshold, metric, seed or reporting family will be changed after result inspection.
