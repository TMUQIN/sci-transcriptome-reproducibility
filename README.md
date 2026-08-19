# sci-transcriptome-reproducibility

**Cross-Study Reproducibility in Mouse Spinal Cord Injury Transcriptomes: Analysis Code and Derived Data**

Status: **PRIVATE RELEASE CANDIDATE RC1 - NOT YET PUBLISHED**  
GitHub repository: https://github.com/TMUQIN/sci-transcriptome-reproducibility  
Archived release: `v1.0.0` (planned)  
Zenodo DOI: `ZENODO_DOI_PENDING`  
Code license: MIT  
Author-generated derived data, tables and documentation: CC BY 4.0

## 1. Purpose

This repository supports the accompanying *Molecular Neurobiology* manuscript.
It assesses how cross-study reproducibility of mouse spinal-cord-injury
transcriptome conclusions depends on gene-level versus Hallmark-level
representation, the temporal contrast, and the biological-sample analysis unit.
It contains frozen derived outputs, analysis scripts and provenance; it is not
the journal submission package.

Cells and nuclei are never treated as independent biological replicates.
GSE47681 is an **orthogonal cross-platform context evaluation** and is not
pooled into the three-dataset synthesis. It is not formal external validation.

## 2. Authors and contact

1. Dexiang Qin
2. Junrui Guo
3. Baicao Li
4. Hebin Hou
5. Qi Zhang
6. Shouchen Li
7. Guangzhi Ning

Corresponding author: Guangzhi Ning  
Department of Orthopedics, Tianjin Medical University General Hospital,
Tianjin 300070, China  
Email: gzning@tmu.edu.cn  
ORCID: 0000-0002-1635-9902

## 3. Repository contents

```text
README.md
REPRODUCE.md
CITATION.cff
LICENSE
DATA_LICENSE_CC_BY_4.0.md
THIRD_PARTY_DATA.md
scripts/                 authoritative analysis scripts and legacy provenance
environment/             historical environment and reconstruction guidance
data/                    frozen derived data
tables/                  synchronized Online Resource 2 S1-S3 source tables
provenance/              protocols, hashes and provenance records
metadata/                code/output map, GEO manifest and integrity manifest
docs/                    workflow, scope, figure map and license boundaries
```

## 4. Current manuscript-facing figures

The public map is synchronized to the final manuscript organization:

- Fig. 1: analytical sample reconstruction and estimability;
- Fig. 2: gene/Hallmark transfer;
- Fig. 3: effect-strength sensitivity;
- Fig. 4: method-specific Hallmark support;
- Fig. 5: focal-program context-dependent transfer profiles;
- Fig. S1: composition sensitivity;
- Fig. S2: estimability;
- Fig. S3: equal-feature-count diagnostic;
- Fig. S4: GSE47681 QC;
- Fig. S5: conditional member stability;
- Fig. S6: matched-set exchangeability.

See `docs/figure_reproduction_map.md` and
`docs/FIGURE_ASSEMBLY_STATUS.md` for the source-data and reproducibility
boundary of each figure.

## 5. Reproducibility boundary

The repository supports inspection and reconstruction of the numerical source
tables. The final presentation assembly used historical project-level plotting
infrastructure that depends on `reports/`, `results/` and helper modules not
included in this minimal public release. Therefore this repository does not
claim clean-clone, pixel-identical regeneration of the published TIFF figures.
The former v5.3.3 builder is retained under `scripts/legacy/` for provenance;
the top-level file with that name is a retired compatibility marker.

Three reconstruction levels are documented in `REPRODUCE.md`:

- Level 1: inspect frozen data and tables from this archive;
- Level 2: recompute selected statistical or permutation analyses from archived
  intermediate data where the required inputs are present;
- Level 3: reconstruct source analyses after downloading public GEO and MSigDB
  inputs under their upstream terms.

## 6. Third-party data

Raw GEO matrices, GSE47681 CEL files and the MSigDB GMT are not redistributed.
The GEO accession roles, the 13-CEL SHA256 manifest and the MSigDB release hash
are documented in `THIRD_PARTY_DATA.md` and `provenance/`.

## 7. Historical and reconstruction environments

The frozen result-generation environment is documented as:

- Python 3.14.2;
- NumPy 2.4.0, pandas 2.3.3, SciPy 1.17.1, gseapy 1.1.13;
- Matplotlib 3.10.8;
- R 4.4.2, Bioconductor 3.20, limma 3.62.2, affy 1.84.0;
- metafor 5.0-1 where applicable.

`environment/environment.yml` is an advisory Python 3.11 reconstruction
environment chosen for compatibility on the current workstation. It is not a
claim that Python 3.11 generated the frozen results. See `environment/ENVIRONMENT.md`.

## 8. Quick start and citation

Do not use the retired v5.3.3 builder to generate current figures. Start with
`REPRODUCE.md`, read the frozen tables under `tables/`, and use
`metadata/CODE_TO_OUTPUT_MAP.csv` to identify the relevant numerical source
files and analysis scripts.

Until a DOI is assigned, cite the manuscript and this repository by title.
After publication, replace the pending GitHub URL and Zenodo DOI in this file
and in `CITATION.cff` with the real identifiers.

## 9. Licensing scope

- `LICENSE` applies to repository code and scripts.
- `DATA_LICENSE_CC_BY_4.0.md` applies to author-generated derived data, tables,
  provenance records and documentation.
- GEO and MSigDB materials retain upstream terms and are not relicensed here.
