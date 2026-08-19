# Changelog

## v1.0.0 - release (2026-08-19)

- Finalized the canonical GitHub repository URL and reserved Zenodo DOI.
- Finalized v1.0.0 citation and release metadata.
- Clarified reproducibility and reconstruction boundaries.
- Refreshed repository integrity hashes after release-metadata updates.
- Preserved the audited derived data, code, provenance, licensing scope, and historical RC1 record.

## v1.0.0-rc1 - private release candidate (2026-08-18)

- Synchronized Supplementary Tables S1-S3 with the current ESM2 archive.
- Updated manuscript-facing terminology to `orthogonal cross-platform context evaluation` for GSE47681.
- Synchronized the figure map to the current Fig. 1-5 and Fig. S1-S6 organization.
- Retained the former v5.3.3 figure builder under `scripts/legacy/` for provenance and retired it at the top level.
- Added MIT licensing for code and a scoped CC BY 4.0 notice for author-generated derived data, tables and documentation.
- Confirmed the corresponding-author address as Tianjin 300070, China.
- Documented Python 3.14.2 / Matplotlib 3.10.8 as the historical result-generation environment and Python 3.11 as an advisory reconstruction environment.
- Added an explicit `.gitignore` exception for the archived permutation NPZ.
- Created the private GitHub repository at https://github.com/TMUQIN/sci-transcriptome-reproducibility; the Zenodo DOI remains pending until a DOI is reserved.

## Planned v1.0.0 public release

- Initialize Git locally, run a clean-clone smoke test, create the GitHub repository and reserve the Zenodo DOI.
- Replace pending identifiers only after the real GitHub URL and reserved DOI exist.
- Tag the final approved release as `v1.0.0` and upload the GitHub tag ZIP to the Zenodo draft.
