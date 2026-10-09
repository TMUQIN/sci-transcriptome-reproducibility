
# Figure Assembly Status

## Current manuscript

This repository accompanies the manuscript:

"Cross-study reproducibility of molecular responses after mouse spinal cord injury varies by temporal contrast"

Prepared for consideration in *Experimental Neurology*.

The current manuscript contains five main figures and six supplementary figures.

## Main figures

- **Figure 1. Contrast-aligned evidence space in mouse SCI.**
  Overview of the three temporal contrasts, eligible sample-level units, and the shared analytical space of 15,331 genes and 50 Hallmark programs. GSE304361 is shown as an additional context dataset and is not included in the three-study temporal synthesis.

- **Figure 2. Contrast-specific molecular correspondence across SCI studies.**
  Gene-level and Hallmark-level directional concordance, leave-one-study-out balanced accuracy, rank correlation, feature-identity permutation-null comparisons, and the independent GSE47681 cross-platform context evaluation.

- **Figure 3. Effect-strength sensitivity analysis.**
  Cross-study concordance and balanced accuracy across pooled absolute-score strata and within a common effect-strength window.

- **Figure 4. Method-specific support for Hallmark programs.**
  Summary of pathway-level support under preranked GSEA and alternative sample-level scoring approaches. The different methods use distinct statistical criteria and are not interpreted as directly interchangeable significance tests.

- **Figure 5. Contrast-specific enrichment of representative SCI-associated programs.**
  Study-specific GSEA normalized enrichment scores for EMT, hypoxia, MYC targets V1, and mTORC1 signaling across the three temporal contrasts. Additional panels show MYC targets V1 overlap-removal sensitivity and mTORC1 responses in the GSE304361 Plexin-B1 context.

## Supplementary figures

- **Figure S1:** Composition sensitivity of four focal programs.
- **Figure S2:** Candidate-dataset eligibility and estimability assessment.
- **Figure S3:** Equal-feature-count diagnostic.
- **Figure S4:** GSE47681 microarray reconstruction quality control.
- **Figure S5:** Conditional leading-edge member stability.
- **Figure S6:** Matched-set comparability and exchangeability assessment.

## Numerical source data

Numerical outputs supporting the figures are preserved in the archived derived-data release and manuscript Supplementary File 2.

The repository provides source tables, analysis scripts, provenance records, and documentation linking analytical outputs to reported figures.

The current figure organization is consistent with the Experimental Neurology manuscript.

## Reproducibility boundary

Numerical result reconstruction and final figure assembly are distinct processes.

The historical figure-generation workflow depended on project-level intermediate files and helper modules that are not all included in this minimal public release.

Accordingly, this repository does not claim that a clean clone can regenerate every final submission TIFF file pixel-for-pixel.

The archived numerical data remain available for inspection and independent comparison with the reported values.

The historical v5.3.3 figure builder is retained under `scripts/legacy/` for provenance only. It is not the current supported figure-generation workflow.

For additional details, see:

- `REPRODUCE.md`
- `docs/figure_reproduction_map.md`
- `metadata/CODE_TO_OUTPUT_MAP.csv`

## Release integrity

This documentation update synchronizes figure descriptions with the current manuscript.

It does not modify source datasets, frozen statistical results, archived numerical outputs, or the Zenodo v1.0.0 scientific release.
