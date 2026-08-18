# Reproducibility instructions

This repository separates numerical reconstruction from final presentation
assembly. It does not claim that a clean clone can regenerate pixel-identical
TIFF figures.

## Level 1 - archived numerical source data

The frozen numerical source tables are available under `data/` and `tables/`.
Their figure/table relationships are listed in
`metadata/CODE_TO_OUTPUT_MAP.csv` and
`docs/figure_reproduction_map.md`.

Supplementary Tables S1-S3 are synchronized byte-for-byte with the current
2026-08-18 Online Resource 2 archive. The repository integrity manifest records
their SHA256 values.

Status: `EXACTLY_REPRODUCIBLE_FROM_ARCHIVED_INPUT` for frozen tables and source
data that are directly archived.

## Level 2 - selected statistical reconstruction

Selected permutation, meta-analysis, sensitivity and GSE47681 branches can be
recomputed where their required archived intermediate files are present:

| Branch | Primary script | Seed | Boundary |
|---|---|---:|---|
| Feature-identity permutation | `scripts/67_postlock_sign_identity_calibration.py` | 20260803 | Complete null draws are archived; rebuilding them requires frozen effects or Level 3 inputs. |
| GSE47681 Hallmark comparison | `scripts/69_gse47681_hallmark_stress_test.py` | 20260803 | Uses archived derived matrices and MSigDB 2026.1 Hallmarks. |
| Effect-strength sensitivity | `scripts/70_effect_strength_calibration.py` | - | Uses archived effect-strength source tables. |
| Pathway meta-analysis | `scripts/23_pathway_effect_meta.py`, `scripts/25_gsema_mkh_meta.py` | - | Small-k inference depends on the archived study-level effects. |

Status: `PARTIAL_ENVIRONMENT_CAPTURE` or
`REPRODUCIBLE_WITH_PUBLIC_THIRD_PARTY_INPUT`, depending on the branch.

## Level 3 - public-source reconstruction

Source reconstruction requires public third-party inputs that are not
redistributed:

- GSE162610, GSE234774 and GSE304399 source matrices from NCBI GEO;
- GSE47681 13 eligible WT CEL files from NCBI GEO, verified against
  `provenance/GSE47681_CEL_SHA256_manifest.csv`;
- MSigDB 2026.1 mouse Hallmark GMT, verified against the recorded hash and used
  under upstream terms.

See `THIRD_PARTY_DATA.md` for accession roles and restrictions.

## Final figure assembly

The current figure hierarchy is Fig. 1-5 plus Fig. S1-S6. The historical
v5.3.3 builder generated an obsolete hierarchy and is retained under
`scripts/legacy/` for provenance only. Final presentation assembly depended on
project-level `reports/`, `results/` and helper modules not present in this
minimal release. Do not use `scripts/72_build_v5_3_3_figures.py` as a supported
reproduction command.

See `docs/FIGURE_ASSEMBLY_STATUS.md` for the explicit boundary.

## Environment

- Historical result-generation environment: Python 3.14.2, NumPy 2.4.0,
  pandas 2.3.3, SciPy 1.17.1, gseapy 1.1.13, Matplotlib 3.10.8, R 4.4.2,
  Bioconductor 3.20, limma 3.62.2, affy 1.84.0 and metafor 5.0-1 where
  applicable.
- Advisory reconstruction environment: `environment/environment.yml`, using
  Python 3.11 for compatibility.

The reconstruction environment is not presented as the historical analysis
environment, and exact transitive dependency capture is not claimed.
