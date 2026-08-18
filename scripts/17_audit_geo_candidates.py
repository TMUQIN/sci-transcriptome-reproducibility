#!/usr/bin/env python3
"""Audit five recent SCI single-cell GEO candidates from official family SOFT files."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import re
from pathlib import Path

import pandas as pd


CANDIDATES = {
    "GSE240727": {
        "expected_gsm": 6,
        "expression_units": "6 condition-level libraries; one at each of uninjured/d1/d3/d7/d14/d28",
        "exact_d1_d7_replication_eligible": "no_one_library_per_endpoint",
        "recommended_role": "descriptive_immune_timecourse_atlas_only",
        "priority": "exclude_from_inferential_download",
        "blocking_reason": "no biological library replication at any endpoint",
        "innovation_conflict": "broad lesion immune-cell heterogeneity across d1-d28 is already the published study premise",
    },
    "GSE247844": {
        "expected_gsm": 8,
        "expression_units": "d28 injured n=4 and sex/age-matched uninjured n=4",
        "exact_d1_d7_replication_eligible": "no_d1_and_d7_absent",
        "recommended_role": "independent_chronic_d28_astrocyte_endpoint_support",
        "priority": "secondary_download_after_main_rewrite",
        "blocking_reason": "single chronic d28 endpoint cannot test acute-to-subacute attenuation",
        "innovation_conflict": "healthy-versus-d28 astrocyte state mapping is already the published contribution",
    },
    "GSE275982": {
        "expected_gsm": 24,
        "expression_units": "12 GEX biological libraries (4 old-naive, 4 old-d28, 4 young-d28) plus paired TCR modality",
        "exact_d1_d7_replication_eligible": "no_early_timepoints",
        "recommended_role": "independent_chronic_T_cell_age_context",
        "priority": "targeted_chronic_immune_sensitivity_only",
        "blocking_reason": "d28-only injury groups; no young-naive group for a full age-by-injury interaction; GEX/TCR pairs are not independent",
        "innovation_conflict": "age-dependent clonal NK-like T-cell expansion and impaired wound healing are already claimed",
    },
    "GSE298545": {
        "expected_gsm": 12,
        "expression_units": "new scRNA aged d3 n=3/uninjured n=3, young d3 n=1/uninjured n=1; each GSM pools two animals; spatial d3 n=2 per age; published combined young arm reuses GSE162610",
        "exact_d1_d7_replication_eligible": "no_only_d3",
        "recommended_role": "aged_within_cohort_d3_support_only_combined_age_analysis_reuses_GSE162610",
        "priority": "secondary_age_generalizability_download",
        "blocking_reason": "only d3; new young endpoints have n=1; each scRNA library pools two mice and combines unenriched with astrocyte-enriched inputs; published age analysis reuses a current primary cohort",
        "innovation_conflict": "age comparison, acute cell-subpopulation shifts and d3 spatial injury zones are already the published story",
    },
    "GSE304361": {
        "expected_gsm": 8,
        "expression_units": "12 biological mice (WT/Plxnb1-KO x uninjured/7dpi, n=3 per condition) multiplexed into 8 Parse sublibraries",
        "exact_d1_d7_replication_eligible": "no_d1_and_sample_demultiplex_key_missing_from_GEO",
        "recommended_role": "conditional_WT_d7_astrocyte_support_if_author_sample_map_obtained",
        "priority": "highest_conditional_metadata_request",
        "blocking_reason": "GSMs are technical sublibraries mixing genotype and injury; GEO has no supplementary cell-to-mouse assignment",
        "innovation_conflict": "Plexin-B1 astrocyte agility, wound corralling and astrocyte-microglia signaling are already claimed",
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def slug(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.lower()).strip("_")


def parse_fields(block: str, prefix: str) -> dict[str, list[str]]:
    fields: dict[str, list[str]] = {}
    for line in block.splitlines():
        if not line.startswith(prefix) or " = " not in line:
            continue
        key, value = line[1:].split(" = ", 1)
        fields.setdefault(key, []).append(value.strip())
    return fields


def first(fields: dict[str, list[str]], key: str, default: str = "") -> str:
    values = fields.get(key, [])
    return values[0] if values else default


def parse_soft(path: Path, accession: str) -> tuple[dict, list[dict]]:
    with gzip.open(path, "rt", encoding="utf-8", errors="strict") as handle:
        text = handle.read()
    series_block = text.split("^SAMPLE = ", 1)[0]
    series_fields = parse_fields(series_block, "!Series_")
    series = {
        "accession": accession,
        "title": first(series_fields, "Series_title"),
        "pubmed": first(series_fields, "Series_pubmed_id"),
        "summary": " ".join(series_fields.get("Series_summary", [])),
        "overall_design": " ".join(series_fields.get("Series_overall_design", [])),
        "bioproject": next(
            (value.rsplit("/", 1)[-1] for value in series_fields.get("Series_relation", [])
             if "bioproject" in value.lower()), ""
        ),
    }
    samples: list[dict] = []
    for block in re.split(r"(?=\^SAMPLE = )", text):
        match = re.match(r"\^SAMPLE = (GSM\d+)", block)
        if not match:
            continue
        fields = parse_fields(block, "!Sample_")
        characteristics: dict[str, list[str]] = {}
        for value in fields.get("Sample_characteristics_ch1", []):
            if ":" in value:
                key, item = value.split(":", 1)
                characteristics.setdefault(slug(key), []).append(item.strip())
        relations = fields.get("Sample_relation", [])
        supplementary = [value for key, values in fields.items()
                         if key.startswith("Sample_supplementary_file_") for value in values]
        row = {
            "accession": accession,
            "gsm": match.group(1),
            "title": first(fields, "Sample_title"),
            "source_name": first(fields, "Sample_source_name_ch1"),
            "biosample": next((value.rsplit("/", 1)[-1] for value in relations
                               if "biosample" in value.lower()), ""),
            "sra": next((value.rsplit("=", 1)[-1] for value in relations
                         if "sra?term=" in value.lower()), ""),
            "n_supplementary_files": len([value for value in supplementary if value != "NONE"]),
            "supplementary_files": ";".join(supplementary),
        }
        for key, values in characteristics.items():
            row[f"characteristic_{key}"] = ";".join(values)
        samples.append(row)
    if len(samples) != CANDIDATES[accession]["expected_gsm"]:
        raise ValueError(
            f"{accession}: expected {CANDIDATES[accession]['expected_gsm']} GSM, found {len(samples)}"
        )
    if pd.Series([row["gsm"] for row in samples]).duplicated().any():
        raise ValueError(f"{accession}: duplicate GSM")
    return series, samples


def annotate_sample_roles(samples: pd.DataFrame) -> pd.DataFrame:
    out = samples.copy()
    out["assay_role"] = "expression"
    out.loc[out.accession.eq("GSE275982") & out.title.str.contains("_TCR", regex=False), "assay_role"] = "paired_TCR"
    out.loc[out.accession.eq("GSE275982") & out.title.str.contains("_GEX", regex=False), "assay_role"] = "GEX"
    out.loc[out.accession.eq("GSE298545") & out.title.str.startswith("Spatial_"), "assay_role"] = "spatial"
    out.loc[out.accession.eq("GSE298545") & ~out.title.str.startswith("Spatial_"), "assay_role"] = "scRNA_composite"
    out.loc[out.accession.eq("GSE304361"), "assay_role"] = "technical_Parse_sublibrary"
    out["independence_note"] = "one GSM treated as one library unless official design states otherwise"
    out.loc[out.accession.eq("GSE240727"), "independence_note"] = "one condition-level GSM per endpoint; no replicate inference"
    out.loc[out.accession.eq("GSE247844"), "independence_note"] = "official design states four biological mice per endpoint"
    out.loc[out.accession.eq("GSE275982"), "independence_note"] = "GEX and TCR from the same condition/replicate are one independence cluster"
    out.loc[out.accession.eq("GSE298545") & out.assay_role.eq("scRNA_composite"), "independence_note"] = (
        "one library/GSM pools two mice: one unenriched whole-tissue plus one astrocyte-enriched input; tag numbers retained"
    )
    out.loc[out.accession.eq("GSE298545") & out.assay_role.eq("spatial"), "independence_note"] = "one spatial section/GSM; only two per age at d3"
    out.loc[out.accession.eq("GSE304361"), "independence_note"] = (
        "technical sublibrary contains mixed genotype/injury samples; GSM is not a biological replicate"
    )
    return out.sort_values(["accession", "gsm"], kind="stable")


def run(raw_root: Path, out_dir: Path) -> None:
    series_rows, sample_rows, inputs = [], [], []
    for accession, policy in CANDIDATES.items():
        path = raw_root / accession / f"{accession}_family.soft.gz"
        if not path.exists() or path.stat().st_size <= 0:
            raise FileNotFoundError(f"missing/non-positive official SOFT: {path}")
        series, samples = parse_soft(path, accession)
        series_rows.append({**series, **policy})
        sample_rows.extend(samples)
        inputs.append({
            "accession": accession,
            "url": f"https://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:-3]}nnn/{accession}/soft/{accession}_family.soft.gz",
            "path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path),
        })
    samples = annotate_sample_roles(pd.DataFrame(sample_rows))
    series = pd.DataFrame(series_rows)
    observed = samples.groupby("accession", observed=True).size().to_dict()
    expected = {key: value["expected_gsm"] for key, value in CANDIDATES.items()}
    if observed != expected:
        raise ValueError(f"series GSM totals differ: observed={observed}, expected={expected}")
    out_dir.mkdir(parents=True, exist_ok=True)
    sample_path = out_dir / "recent_geo_candidate_sample_audit.tsv"
    series_path = out_dir / "recent_geo_candidate_eligibility.tsv"
    samples.to_csv(sample_path, sep="\t", index=False)
    series.to_csv(series_path, sep="\t", index=False)
    provenance = {
        "audit_date": "2026-07-28",
        "audit_role": "official_GEO_sample_design_eligibility_before_expression_download",
        "n_series": len(series), "n_gsm": len(samples),
        "decision_rule": "biological libraries/mice, timepoint coverage, assay pairing and sample resolvability; never cell count",
        "inputs": inputs,
        "outputs": [
            {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in (sample_path, series_path)
        ],
        "interpretation_limit": "eligibility uses GEO SOFT and must be reconciled with the primary paper before inferential analysis",
    }
    (out_dir / "recent_geo_candidate_audit_provenance.json").write_text(
        json.dumps(provenance, ensure_ascii=False, indent=2), encoding="utf-8"
    )
    print(f"PASS: {len(series)} series and {len(samples)} GSM audited from official family SOFT")


def selftest() -> None:
    text = """^SERIES = GSE1
!Series_title = T
!Series_pubmed_id = 1
!Series_overall_design = design
^SAMPLE = GSM1
!Sample_title = A
!Sample_source_name_ch1 = cord
!Sample_characteristics_ch1 = injury: 7 dpi
!Sample_relation = BioSample: https://www.ncbi.nlm.nih.gov/biosample/SAMN1
!Sample_supplementary_file_1 = NONE
"""
    fields = parse_fields(text, "!Sample_")
    assert first(fields, "Sample_title") == "A"
    assert slug("Days post-SCI") == "days_post_sci"
    print("[selftest] PASS: SOFT repeated fields and characteristic slugs")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--raw-root", default="data_raw")
    parser.add_argument("--out-dir", default="tables")
    parser.add_argument("--selftest", action="store_true")
    args = parser.parse_args()
    if args.selftest:
        selftest()
        return
    run(Path(args.raw_root), Path(args.out_dir))


if __name__ == "__main__":
    main()
