# Final calibration and manuscript v5 protocol — frozen 2026-07-29

## Inputs frozen

Manuscript v4, all historical analyses, and every file already present in `reports/phase_reproducibility_calibration_2026_07/` are frozen inputs. This phase writes only to `reports/phase_final_calibration_and_v5_2026_07/`, `manuscript_v5/`, and new versioned scripts. No historical result may be silently regenerated or replaced.

## Five decision questions

1. Does d1/d7 program concordance exceed expression-, statistical-feature- and correlation-matched random gene sets?
2. Is the absence of direct d7-minus-d1 program lift stable to score method and sample perturbation?
3. Do EMT, hypoxia, MYC V1 and mTORC1 have stable rather than label-level leading-edge support?
4. How much whole-lesion evidence is composition- or ontology-sensitive?
5. Do MYC and mTORC1 support, respectively, partial cell-cycle separation and formal context dependence?

## Matched-null pilot rules

The four focal Hallmarks are tested before extension to all 50. Pilot size is 500 unique random sets per Hallmark, contrast, null level and scoring method. Expansion to 10,000 is allowed only after pilot exchangeability diagnostics pass.

- Level 1 primary null: exact set size; member-level nearest-neighbour matching on mean expression, detectability and expression variance.
- Level 2 outcome-conditioned sensitivity null: Level 1 variables plus gene-level SE, number of estimable studies, positive-study count and all-study sign concordance. Because the last two encode outcome information, Level 2 is explicitly a conservative sensitivity analysis, not the primary null.
- Level 3 correlation-preserving null: Level 2 matching plus selection for mean pairwise correlation, mean squared correlation/effective set size and coexpression-module composition.

Exchangeability targets are frozen as absolute standardized mean difference <=0.25 for continuous member-level covariates, exact set size, module-composition L1 distance <=0.20, absolute mean-correlation difference <=0.05, and effective-size relative difference <=20%. These are feasibility thresholds, not significance thresholds. Duplicate-set rate must be <=5% and every retained set must contain unique genes.

The primary set-level reproducibility statistic is the normalized directional consistency of three standardized study effects:

`abs(mean(effect)) / root_mean_square(effect)`.

It lies in [0,1], reaches 1 when study effects have identical direction and magnitude, and falls as signs or magnitudes disagree. All-study sign concordance and held-out direction accuracy are reported separately. Observed and null sets use identical biological samples, gene universe, score construction and contrasts.

Empirical P uses `(1 + number(null >= observed))/(B + 1)`. BH correction is applied within null level, score method and metric across the four focal programs and three contrasts. A binary sign-concordance P is descriptive because of its coarse support.

## Conclusion classes

- all three nulls support: `calibrated program lift`;
- Levels 1/2 support but Level 3 does not: `aggregation-associated concordance`, not biological-program superiority;
- endpoint support without direct-delta support: `endpoint-specific lift without temporal-change lift`;
- failed exchangeability: `not fully calibrated`.

## Remaining analyses and stop rule

After matched-null feasibility, run Hallmark/leading-edge stability, composition/ontology sensitivity, MYC de-overlap/residual analyses and the mTORC1 context matrix. Preserve all method disagreement. Once these predefined analyses and claim audit are complete, stop adding pathways or datasets and draft v5. Do not lower FDR, replace mKH, or search additional pathway libraries for a positive result.

