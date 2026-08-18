# Third-party data — provenance and redistribution restrictions

This repository **does not redistribute** any third-party raw data. It contains
only **derived, frozen outputs** (normalized matrices, effect tables, scores,
null draws, summaries) produced by the analysis scripts from publicly available
source series. Raw GEO matrices, raw CEL files, and the MSigDB GMT are **not**
included.

## GEO series

| Accession | Role in manuscript | Platform | Raw required to reproduce a branch? | Where to obtain | Repository content |
|-----------|--------------------|----------|--------------------------------------|----------------|--------------------|
| GSE162610 | Primary scRNA-seq (all-state aggregate pseudobulk) | 10x/V3 | Yes (Level 3) | NCBI GEO | Derived effects only |
| GSE234774 | Primary scRNA-seq | 10x/V3 | Yes (Level 3) | NCBI GEO | Derived effects only |
| GSE304399 | Primary scRNA-seq | 10x | Yes (Level 3) | NCBI GEO | Derived effects only |
| GSE304361 | Endpoint/genotype context evaluation | scRNA-seq | No lesion-site result | NCBI GEO | Derived support summary only |
| GSE159638 | Candidate audit | scRNA-seq | Audit only | NCBI GEO | Provenance only |
| GSE172167 | Candidate audit | scRNA-seq | Contributes no lesion-site result | NCBI GEO | Provenance only |
| GSE205029 | Temporal eligibility assessment | scRNA-seq | Excluded before expression modeling | NCBI GEO | Eligibility provenance only (`data/GSE205029_wt_sample_eligibility.csv`, `data/GSE205029_contrast_estimability.csv`) |
| GSE47681 | Orthogonal cross-platform context evaluation | Affymetrix mouse4302 | Yes (Level 3, 13 WT CEL) | NCBI GEO | Derived RMA/gene effects + CEL SHA256 manifest; **CEL files not redistributed** |

### GSE47681 CEL integrity manifest

The 13 eligible WT CEL-file SHA256 hashes are recorded in
`provenance/GSE47681_CEL_SHA256_manifest.csv` for integrity verification. The CEL
files themselves are **not** included. Obtain them from NCBI GEO and verify
against that manifest before running `scripts/68_gse47681_raw_rma_limma.R`.

## MSigDB

- Exact release used: **MSigDB 2026.1** Hallmark gene sets, mouse symbols
  (`msigdb_mh.all.v2026.1.Mm.symbols.gmt`).
- Expected source: MSigDB / Gene Set Enrichment Analysis (GSEA) official distribution.
- Locally recorded file hash (from provenance):
  `3a21be724a87dc0375955e725ca9688b87a26e7b74b62fba0c62da0967b789f7`
  (48,007 bytes) — see `provenance/GSE47681_hallmark_gsea_provenance.json`.
- **Redistribution:** the GMT file is **not** redistributed by this repository.
  Verify MSigDB licensing/redistribution permission before including a GMT in any
  public release.

## Licensing note

Third-party GEO and MSigDB material keeps its own upstream license. This
repository does **not** relicense any third-party content.

Legacy internal identifiers containing `stress_test` refer to the manuscript's
orthogonal cross-platform or genotype-context evaluation. They are retained
for provenance and do not denote formal external validation.
