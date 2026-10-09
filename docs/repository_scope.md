
# Repository scope

## What this repository is

This repository is a public reproducibility companion to the manuscript:

"Cross-study reproducibility of molecular responses after mouse spinal cord injury varies by temporal contrast"

The manuscript is prepared for consideration in *Experimental Neurology*.

The repository provides versioned analysis code, author-generated numerical source data, and analytical provenance supporting the reported findings.

It contains:

- Analysis scripts supporting the archived statistical results.
- Frozen derived numerical data and Supplementary Tables S1-S3.
- Provenance records, analytical protocols, hashes, and random seeds.
- Metadata linking code, datasets, and reported numerical outputs.
- Documentation of software environments, data provenance, and reproducibility boundaries.

The versioned scientific release is archived on Zenodo:

DOI: https://doi.org/10.5281/zenodo.22005689
        
        
        
        

Release version: v1.0.0

## What this repository is NOT

- It is not the journal submission package containing the manuscript, title page, cover letter, or final figure files.
- It does not replace manuscript Supplementary File 2 (the frozen source-data ZIP archive) or Supplementary File 3 (the formatted XLSX tables).
- It does not redistribute third-party raw GEO matrices, CEL files, or the MSigDB GMT.
- It does not contain every historical project intermediate or helper module required to rerun all analytical branches from a clean clone.
- It does not claim pixel-identical regeneration of every final publication TIFF figure.

Numerical reconstruction and final figure assembly have distinct reproducibility boundaries, as documented in `REPRODUCE.md` and `docs/FIGURE_ASSEMBLY_STATUS.md`.

## Inclusion criteria for scripts

Scripts are included when they directly generated, processed, audited, or supported numerical results reported in the manuscript, including Figures 1-5, Supplementary Figures S1-S6, Supplementary Tables S1-S3, and the archived source-data outputs.

Historical figure-generation and manuscript-packaging utilities are not necessarily included in the minimal public release. Their exclusion does not imply that the archived numerical results were recalculated or altered.

## Exclusions

The public release does not redistribute:

- Raw third-party GEO expression matrices or CEL files.
- The MSigDB 2026.1 mouse Hallmark GMT file.
- Internal manuscript drafts, reviewer correspondence, or private project records.
- Local absolute paths, credentials, or personal configuration files.

## Licensing

- Repository code and scripts are distributed under the MIT license.
- Author-generated derived data, tables, provenance records, and documentation are covered by the accompanying CC BY 4.0 notice.
- Third-party datasets and reference resources retain their original upstream licensing and usage conditions.

## Scientific data integrity

This documentation update changes manuscript-facing metadata and repository descriptions only.

It does not modify frozen scientific results, statistical thresholds, source datasets, analysis outputs, or the archived Zenodo v1.0.0 scientific release.
