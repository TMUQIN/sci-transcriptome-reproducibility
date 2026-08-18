#!/usr/bin/env python3
"""Build auditable sample-level metadata from the frozen GEO check records.

The GEO JSON records in ``logs/GSE*_geo_check.json`` are the source of accession and
sample-title facts.  This script adds only conservative, title-derived fields.  It does
not invent animal identities: when a subject cannot be verified, the GEO library is the
independence unit and ``subject_verification`` is marked ``library_only``.

Outputs
-------
tables/sample_metadata.tsv
    Canonical machine-readable sample crosswalk used by downstream loaders/models.
tables/sample_metadata_schema.tsv
    Field definitions and allowed values.
tables/dataset_eligibility.tsv
    Dataset/arm-level role and inferential restrictions.
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import pandas as pd


DATASETS = (
    "GSE234774", "GSE230765", "GSE162610", "GSE172167",
    "GSE182803", "GSE192824", "GSE256397", "GSE304399",
)


def phase_from_dpi(dpi: float | None) -> str:
    if dpi is None:
        return "unknown"
    if dpi == 0:
        return "uninjured"
    if dpi <= 3:
        return "acute"
    if dpi <= 14:
        return "subacute"
    if dpi < 60:
        return "intermediate"
    return "chronic"


def base_row(dataset: str, accession: str, title: str) -> dict:
    return dict(
        dataset=dataset,
        sample=accession,
        sample_alias=title,
        library_key=title,
        library_key_source="geo_record.samples.title",
        subject_id=accession,
        subject_verification="library_only",
        condition="unknown",
        group="unknown",
        dpi=pd.NA,
        phase="unknown",
        region="unknown",
        segment="spinal_cord",
        modality="unknown",
        primary_arm="no",
        inference_role="descriptive",
        exclusion_reason="not_yet_classified",
        source_record=f"logs/{dataset}_geo_check.json",
        source_field="geo_record.samples.title",
        verification_status="title_derived",
    )


def parse_gse162610(row: dict) -> dict:
    title = row["sample_alias"].lower()
    row["modality"] = "scRNA-seq"
    row["library_key"] = row["sample_alias"]
    row["region"] = "lesion_site"
    row["primary_arm"] = "yes"
    if title.startswith("uninj"):
        dpi = 0
        condition = "uninjured"
    else:
        m = re.match(r"(\d+)dpi_", title)
        if not m:
            return row
        dpi = int(m.group(1))
        condition = "injured"
    row.update(condition=condition, group="uninjured" if dpi == 0 else f"dpi_{dpi}",
               dpi=dpi, phase=phase_from_dpi(dpi), inference_role="within_study_inference",
               exclusion_reason="", verification_status="geo_title_verified")
    return row


def parse_gse172167(row: dict) -> dict:
    title = row["sample_alias"].lower()
    row.update(modality="snRNA-seq", region="lesion_remote_lumbar", primary_arm="yes")
    if title.startswith("uninjured"):
        dpi, condition = 0, "uninjured"
    elif title.startswith("1dpi"):
        dpi, condition = 1, "injured"
    else:
        m = re.match(r"(\d+)wpi_", title)
        if not m:
            return row
        dpi, condition = 7 * int(m.group(1)), "injured"
    row.update(condition=condition, group="uninjured" if dpi == 0 else f"dpi_{dpi}",
               library_key=title,
               dpi=dpi, phase=phase_from_dpi(dpi), inference_role="within_study_inference",
               exclusion_reason="remote_compartment_analyse_separately",
               verification_status="geo_title_verified")
    return row


def parse_gse230765(row: dict) -> dict:
    title = row["sample_alias"]
    m = re.match(r"(RNA|ATAC)_(2M|D7|UN)_(.+)$", title, flags=re.I)
    if not m:
        return row
    modality, time, animal = m.groups()
    dpi = {"UN": 0, "D7": 7, "2M": 60}[time.upper()]
    condition = "uninjured" if dpi == 0 else "injured"
    # RNA and ATAC entries with the same phase/animal suffix are paired modalities.
    row.update(
        library_key=f"{time.upper()}_{animal}",
        subject_id=f"{row['dataset']}_{time.upper()}_{animal}",
        subject_verification="title_matched_cross_modality",
        condition=condition,
        group="uninjured" if dpi == 0 else f"dpi_{dpi}",
        dpi=dpi,
        phase=phase_from_dpi(dpi),
        region="lesion_site",
        modality="snRNA-seq" if modality.upper() == "RNA" else "snATAC-seq",
        primary_arm="yes",
        inference_role="regulatory_support",
        exclusion_reason="same_study_family_as_GSE234774_not_independent_replication",
        verification_status="geo_title_verified",
    )
    return row


def parse_gse234774(row: dict) -> dict:
    title = row["sample_alias"].lower()
    key_match = re.search(r",\s*([^,]+?)\s*\((?:rnaseq|spatial)\)\s*$", title)
    if key_match:
        row["library_key"] = key_match.group(1)
    row.update(region="mid_thoracic", modality="snRNA-seq" if "(rnaseq)" in title else "spatial")
    if "uninjured_" in title and "(rnaseq)" in title:
        row.update(condition="uninjured", group="uninjured", dpi=0, phase="uninjured",
                   primary_arm="yes", inference_role="within_study_inference",
                   exclusion_reason="", verification_status="geo_title_verified")
        return row
    m = re.search(r"timecourse_(1d|4d|7d|14d|1m|2m)_", title)
    if m:
        dpi = {"1d": 1, "4d": 4, "7d": 7, "14d": 14, "1m": 30, "2m": 60}[m.group(1)]
        row.update(condition="injured", group=f"dpi_{dpi}", dpi=dpi,
                   phase=phase_from_dpi(dpi), primary_arm="yes",
                   inference_role="within_study_inference", exclusion_reason="",
                   verification_status="geo_title_verified")
        return row
    if "timecourse_old_" in title:
        row.update(condition="injured", group="timecourse_old_unresolved",
                   inference_role="excluded", exclusion_reason="OLD_timepoint_not_resolved_from_GEO_title")
    elif "(spatial)" in title:
        row.update(inference_role="spatial_support", exclusion_reason="condition_and_section_map_require_spatial_metadata")
    else:
        row.update(inference_role="secondary_experimental_arm",
                   exclusion_reason="drug_severity_sex_or_mechanism_arm_not_primary_timecourse")
    return row


def parse_gse182803(row: dict) -> dict:
    title = row["sample_alias"].lower()
    row.update(modality="scRNA-seq", region="lesion_site", primary_arm="yes")
    if "healthy" in title:
        row.update(condition="uninjured", group="uninjured_actseq" if "act-seq" in title else "uninjured",
                   dpi=0, phase="uninjured")
    else:
        m = re.search(r"(3|14)dpi", title)
        if m:
            dpi = int(m.group(1))
            row.update(condition="injured", group=f"dpi_{dpi}", dpi=dpi, phase=phase_from_dpi(dpi))
    row.update(inference_role="descriptive", exclusion_reason="one_library_per_group_and_Act-seq_not_a_control_replicate",
               verification_status="geo_title_verified")
    return row


def parse_gse192824(row: dict) -> dict:
    title = row["sample_alias"].lower()
    dpi = 3 if title.startswith("injured") else 0
    row.update(condition="injured" if dpi else "uninjured", group=f"dpi_{dpi}" if dpi else "uninjured",
               library_key="Troy_SCI" if dpi else "Troy_UI",
               dpi=dpi, phase=phase_from_dpi(dpi), region="ependymal_enriched", modality="scRNA-seq",
               primary_arm="yes", inference_role="descriptive",
               exclusion_reason="one_library_per_group", verification_status="geo_title_verified")
    return row


def parse_gse256397(row: dict) -> dict:
    title = row["sample_alias"].lower()
    if title.startswith("uninjure"):
        dpi, condition = 0.0, "uninjured"
    else:
        m = re.match(r"(3|24|72)h after injure", title)
        if not m:
            return row
        dpi, condition = int(m.group(1)) / 24.0, "injured"
    side = "rostral" if "rostral" in title else "caudal"
    distance = "0.5mm" if "0.5mm" in title else "1mm"
    row.update(condition=condition, group="uninjured" if dpi == 0 else f"hour_{int(dpi * 24)}",
               dpi=dpi, phase=phase_from_dpi(dpi), region=f"{side}_{distance}_from_epicenter",
               modality="spatial_transcriptomics", primary_arm="yes", inference_role="spatial_support",
               exclusion_reason="one_section_per_time_by_region_cell; subject_pairing_not_verified",
               verification_status="geo_title_verified")
    return row


def parse_gse304399(row: dict) -> dict:
    """Parse GEX/ATAC titles and map them to barcode prefixes."""
    title = row["sample_alias"]
    match = re.match(
        r"(Unsorted|NeuN-/Sox10-|Foxj1\+),\s*(Uninjured|1dpi|3dpi|7dpi|28dpi),\s*rep(\d+),\s*(ATAC|GEX)$",
        title, flags=re.I,
    )
    if not match:
        row.update(inference_role="excluded", exclusion_reason="unparsed_GEO_title")
        return row
    population, timepoint, rep, assay = match.groups()
    population_key = {"Unsorted": "ALL", "NeuN-/Sox10-": "DEPL", "Foxj1+": "FT"}[population]
    stage_key = "U" if timepoint.lower() == "uninjured" else timepoint.lower()
    # The author barcode prefix uses internal library numbers for the lone 3-dpi
    # and two 28-dpi unsorted libraries (verified from the downloaded barcode list),
    # while U/1/7 dpi retain the GEO replicate number.
    if population_key in {"ALL", "DEPL"} and stage_key == "3dpi":
        internal_suffix = "6"
    elif population_key == "ALL" and stage_key == "28dpi":
        internal_suffix = str(int(rep) + 1)
    else:
        internal_suffix = {"ALL": rep, "DEPL": "5", "FT": "4"}[population_key]
    dpi = 0 if stage_key == "U" else int(stage_key.replace("dpi", ""))
    primary = population_key == "ALL" and assay.upper() == "GEX"
    row.update(
        library_key=f"{population_key}_{stage_key}_{internal_suffix}",
        library_key_source="GEO_title_to_barcode_prefix_crosswalk",
        subject_id=f"GSE304399_{population_key}_{stage_key}_{internal_suffix}",
        subject_verification="title_matched_cross_modality",
        condition="uninjured" if dpi == 0 else "injured",
        group="uninjured" if dpi == 0 else f"dpi_{dpi}",
        dpi=dpi, phase=phase_from_dpi(dpi), region="lesion_site", segment="mid_thoracic",
        modality="multiome-RNA" if assay.upper() == "GEX" else "multiome-ATAC",
        primary_arm="yes" if primary else "no",
        inference_role="within_study_inference" if primary else "regulatory_support",
        exclusion_reason="" if primary else "sorted_or_ATAC_same_cohort_support",
        verification_status="geo_title_verified",
    )
    return row


PARSERS = {
    "GSE234774": parse_gse234774,
    "GSE230765": parse_gse230765,
    "GSE162610": parse_gse162610,
    "GSE172167": parse_gse172167,
    "GSE182803": parse_gse182803,
    "GSE192824": parse_gse192824,
    "GSE256397": parse_gse256397,
    "GSE304399": parse_gse304399,
}


SCHEMA = [
    ("dataset", "GEO series accession", "GSE[0-9]+"),
    ("sample", "Canonical library accession and join key", "GSM[0-9]+; unique within dataset"),
    ("sample_alias", "Original GEO sample title used for crosswalk", "verbatim string"),
    ("library_key", "Key present in cell-level metadata/barcodes", "dataset-specific exact string"),
    ("library_key_source", "Provenance for the cell-to-GSM join key", "source field or override table"),
    ("subject_id", "Biological independence/paired-modality unit", "verified ID or sample accession"),
    ("subject_verification", "Confidence in subject mapping", "library_only|title_matched_cross_modality|verified_metadata"),
    ("condition", "Injury exposure", "uninjured|injured|unknown"),
    ("group", "Model group; one level per sampled condition/time", "uninjured|dpi_N|hour_N|other"),
    ("dpi", "Days post injury; fractional values allowed", ">=0 or NA"),
    ("phase", "Coarse descriptive time bin, not a substitute for group", "uninjured|acute|subacute|intermediate|chronic|unknown"),
    ("region", "Anatomical sampling compartment", "controlled dataset-specific string"),
    ("segment", "Spinal segment", "controlled string"),
    ("modality", "Assay modality", "controlled string"),
    ("primary_arm", "Included in prespecified primary arm", "yes|no"),
    ("inference_role", "Permitted evidence role", "within_study_inference|regulatory_support|spatial_support|descriptive|secondary_experimental_arm|excluded"),
    ("exclusion_reason", "Why a sample/arm is restricted", "empty or explicit reason"),
    ("source_record", "Local provenance record", "relative path"),
    ("source_field", "Field used from source", "JSON field path"),
    ("verification_status", "Level of metadata verification", "geo_title_verified|title_derived|full_metadata_verified"),
]


ELIGIBILITY = [
    ("GSE234774", "primary_timecourse_snRNA", "lesion_site_temporal_discovery", "eligible_after_cell_metadata_crosswalk", "2-3 libraries/group; OLD and secondary arms excluded; same study family as GSE230765"),
    ("GSE162610", "primary_timecourse_scRNA", "independent_acute_lesion_replication", "eligible_with_small_n_caveat", "3 controls, 3 at 1 dpi, 2 at 3/7 dpi; verify pooled-animal structure"),
    ("GSE230765", "paired_RNA_ATAC", "regulatory_support", "not_independent_replication", "paired modalities and same Tabulae study family; use for chromatin support"),
    ("GSE172167", "remote_lumbar_snRNA", "separate_remote_compartment_arm", "eligible_with_compartment_restriction", "3 libraries/time; never pool as lesion-site replication"),
    ("GSE182803", "immune_scRNA", "descriptive_reference", "not_inferential", "one library/group; Act-seq control is not a replicate"),
    ("GSE192824", "ependymal_scRNA", "descriptive_reference", "not_inferential", "one library/group"),
    ("GSE256397", "spatial_time_region", "spatial_support", "not_primary_gene_level_inference", "one section per time-region cell in GEO; subject structure unresolved"),
    ("GSE304399", "unsorted_multiome_RNA", "independent_lesion_site_replication", "eligible_after_barcode_crosswalk", "2026 Nature Neuroscience cohort; unsorted GEX has 3 libraries at uninjured/1/7 dpi and 2 at 28 dpi; ATAC/sorted fractions are same-cohort support"),
]


def build(project_root: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows = []
    for dataset in DATASETS:
        path = project_root / "logs" / f"{dataset}_geo_check.json"
        if not path.exists():
            raise FileNotFoundError(path)
        record = json.loads(path.read_text(encoding="utf-8"))["geo_record"]
        for sample in record["samples"]:
            row = base_row(dataset, sample["accession"], sample["title"])
            rows.append(PARSERS[dataset](row))
    meta = pd.DataFrame(rows).sort_values(["dataset", "sample"], kind="stable")
    override_path = project_root / "tables" / "sample_library_key_overrides.tsv"
    if override_path.exists():
        overrides = pd.read_csv(override_path, sep="\t", dtype=str)
        if overrides.duplicated(["dataset", "sample"]).any():
            raise ValueError("duplicate dataset+sample in sample_library_key_overrides.tsv")
        override_map = overrides.set_index(["dataset", "sample"])
        for index, row in meta.iterrows():
            key = (row["dataset"], row["sample"])
            if key in override_map.index:
                override = override_map.loc[key]
                meta.loc[index, "library_key"] = override["library_key"]
                meta.loc[index, "library_key_source"] = (
                    f"{override['source_record']}::{override['source_field']}")
                meta.loc[index, "verification_status"] = override["verification_status"]
    if meta.duplicated(["dataset", "sample"]).any():
        raise ValueError("duplicate dataset+sample keys")
    schema = pd.DataFrame(SCHEMA, columns=["field", "definition", "allowed_or_format"])
    eligibility = pd.DataFrame(ELIGIBILITY, columns=["dataset", "arm", "analysis_role", "eligibility", "rationale"])
    return meta, schema, eligibility


def selftest(project_root: Path) -> None:
    meta, _, _ = build(project_root)
    assert len(meta) == 77 + 16 + 10 + 15 + 4 + 2 + 16 + 42
    assert not meta.duplicated(["dataset", "sample"]).any()
    g162 = meta[meta.dataset.eq("GSE162610")]
    assert g162.groupby("group").size().to_dict() == {"dpi_1": 3, "dpi_3": 2, "dpi_7": 2, "uninjured": 3}
    g192 = meta[meta.dataset.eq("GSE192824")]
    assert set(g192.inference_role) == {"descriptive"}
    paired = meta[(meta.dataset.eq("GSE230765")) & meta.sample_alias.str.contains("_2M_18")]
    assert paired.subject_id.nunique() == 1 and paired.modality.nunique() == 2
    g304 = meta[meta.dataset.eq("GSE304399")]
    assert g304.groupby(["modality", "inference_role"]).size().to_dict() == {
        ("multiome-ATAC", "regulatory_support"): 21,
        ("multiome-RNA", "regulatory_support"): 9,
        ("multiome-RNA", "within_study_inference"): 12,
    }
    print("[selftest] PASS: 182 unique GEO samples; group counts and modality pairing verified")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", default=str(Path(__file__).resolve().parents[1]))
    ap.add_argument("--selftest", action="store_true")
    args = ap.parse_args()
    root = Path(args.project_root)
    if args.selftest:
        selftest(root)
        return
    meta, schema, eligibility = build(root)
    out = root / "tables"
    out.mkdir(parents=True, exist_ok=True)
    meta.to_csv(out / "sample_metadata.tsv", sep="\t", index=False)
    schema.to_csv(out / "sample_metadata_schema.tsv", sep="\t", index=False)
    eligibility.to_csv(out / "dataset_eligibility.tsv", sep="\t", index=False)
    print(f"Wrote {len(meta)} samples to {out / 'sample_metadata.tsv'}")


if __name__ == "__main__":
    main()
