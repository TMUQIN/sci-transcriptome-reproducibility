# Environment record

This directory distinguishes the historical environment used to generate the
frozen results from the advisory environment supplied for reconstruction.

## Historical result-generation environment

The frozen analysis and figure-related numerical outputs were generated or
documented with:

- Python 3.14.2
- NumPy 2.4.0
- pandas 2.3.3
- SciPy 1.17.1
- gseapy 1.1.13
- Matplotlib 3.10.8
- R 4.4.2
- Bioconductor 3.20
- limma 3.62.2
- affy 1.84.0
- metafor 5.0-1 where applicable

The captured `R_sessionInfo_GSE47681.txt` is authoritative for the GSE47681
RMA/limma branch. The Python package record is documented rather than a full
transitive lock.

## Reconstruction environment

`environment.yml` creates Python 3.11 with bounded package ranges and
Matplotlib 3.10.8. This choice was made for compatibility on the current
workstation. It is an advisory reconstruction environment, not a statement
that Python 3.11 generated the frozen results.

`requirements.txt` records key historical Python pins for reference. It is not
a complete lockfile. Exact transitive versions and platform-specific build
details were not captured at freeze time.

## Figure assembly boundary

The final TIFF presentation assembly used project-level plotting infrastructure
outside this minimal repository. Consequently, exact pixel-identical figure
regeneration from a clean clone is not claimed. See
`docs/FIGURE_ASSEMBLY_STATUS.md`.
