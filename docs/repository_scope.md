# Repository scope

## What this repository is

A **public reproducibility companion** for the *Molecular Neurobiology* manuscript
"Cross-Study Reproducibility in Mouse Spinal Cord Injury Transcriptomes Varies by
Gene-versus-Hallmark Representation and Temporal Contrast."

It contains:
- the **minimal authoritative analysis scripts** that produced the frozen numerical results;
- the **frozen derived source data** and supplementary tables (content identical to
  Online Resource 2, minus the ZIP packaging and its integrity manifest);
- **provenance** (protocols, audit records, JSON provenance with hashes/seeds);
- **metadata** (code manifest, code-to-output map, integrity & GEO manifests);
- **documentation** (README, REPRODUCE, ENVIRONMENT, analysis workflow, license guide).

## What this repository is NOT

- Not the journal submission package (Manuscript, TitlePage, CoverLetter, ESM1, ESM3, Figures 1–5).
- Not a replacement for Online Resource 2 (ESM2) — ESM2 remains the frozen
  source-data archive; this repo holds the executable code and provenance that
  complement it.
- Not a distributor of third-party data — raw GEO matrices, CEL files, and the
  MSigDB GMT are **not** included.
- This repository is the versioned **v1.0.0 release** accompanying the manuscript.
  GitHub URL and Zenodo DOI remain placeholders until the real repository and
  archive identifiers are assigned.

## Inclusion criteria for `scripts/`

A script is included if it **directly generated, processed, audited, or supported**
a numerical result reported in the manuscript (Figs 1–5, S1–S6, Tables S1–S3,
or ESM2). Final presentation assembly is documented separately because the
historical plotting infrastructure is not part of the minimal release.
Manuscript-generation, packaging, validation, and QA utilities (e.g. the
`*_build_molecular_neurobiology_*` / `*_validate_*` / `*_prepare_*` family) are
**excluded** from the minimal set because they do not change scientific results.

## Exclusions (explicitly not redistributed)

- Any raw GEO count matrix or CEL file.
- The MSigDB 2026.1 GMT file.
- Internal redlines, reviewer correspondence, and failed experimental branches.
- Local absolute paths, user-profile paths, and credentials (none present).
