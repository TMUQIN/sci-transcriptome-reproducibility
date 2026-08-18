#!/usr/bin/env python3
"""Audit GSE205029 temporal estimability without reading expression counts."""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime
from pathlib import Path

import pandas as pd


ACCESSION = "GSE205029"
MIN_CELLS = 20
MIN_REPLICATES = 2


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    root = args.root.resolve()
    raw = root / "data_raw/GSE205029"
    out = root / "reports/phase_submission_extension_audit_2026_08"
    out.mkdir(parents=True, exist_ok=True)

    metadata_path = raw / "GSE205029_meta_data.csv.gz"
    barcodes_path = raw / "GSE205029_sample_barcodes.csv.gz"
    features_path = raw / "GSE205029_features.csv.gz"
    metadata = pd.read_csv(metadata_path)
    barcodes = pd.read_csv(barcodes_path, header=None, names=["sample_barcode"])
    features = pd.read_csv(features_path, header=None, names=["gene_symbol"])

    if len(metadata) != 18_203 or len(barcodes) != 18_203:
        raise ValueError(f"Unexpected cell totals: metadata={len(metadata)}, barcodes={len(barcodes)}")
    if len(features) != 55_487:
        raise ValueError(f"Unexpected feature total: {len(features)}")
    if metadata["barcode"].duplicated().any():
        # Cell barcodes can repeat across samples; the composite identifier must not.
        composite = metadata["orig.ident"].astype(str) + ":" + metadata["barcode"].astype(str)
        if composite.duplicated().any():
            raise ValueError("Duplicate composite cell identifier")

    wt = metadata.loc[metadata["genotype"].eq("WT")].copy()
    sample_counts = (
        wt.groupby(["treatment", "sample", "orig.ident"], as_index=False, observed=True)
        .size()
        .rename(columns={"size": "retained_cells"})
        .sort_values(["treatment", "sample"], kind="stable")
    )
    gsm_map = {
        1771: "GSM6204332", 10120: "GSM6204335",
        3821: "GSM6204331", 10115: "GSM6204337",
        3820: "GSM6204330", 10117: "GSM6204339",
        10113: "GSM6204340", 10114: "GSM6204341",
    }
    time_map = {"WT_nsci": "uninjured", "Ryk_WT_1": "day_1", "Ryk_WT_7": "day_7", "Ryk_WT_14": "day_14"}
    sample_counts["gsm"] = sample_counts["sample"].map(gsm_map)
    sample_counts["timepoint"] = sample_counts["treatment"].map(time_map)
    sample_counts["passes_min_cells"] = sample_counts["retained_cells"].ge(MIN_CELLS)
    sample_counts["sample_barcode_labels"] = sample_counts["sample"].map(
        lambda value: ";".join(sorted({
            item.split("-", 1)[1].rsplit("-", 1)[0]
            for item in barcodes.loc[barcodes.sample_barcode.str.startswith(f"{value}-"), "sample_barcode"].astype(str)
        }))
    )
    expected_barcode_label = {
        "WT_nsci": "Ryk_WT_NSCI", "Ryk_KO_nsci": "Ryk_KO_NSCI",
        "Ryk_WT_1": "Ryk_WT_1", "Ryk_KO_1": "Ryk_KO_1",
        "Ryk_WT_7": "Ryk_WT_7", "Ryk_KO_7": "Ryk_KO_7",
        "Ryk_WT_14": "Ryk_WT_14", "Ryk_KO_14": "Ryk_KO_14",
    }
    sample_counts["label_consistent"] = sample_counts.apply(
        lambda row: expected_barcode_label[row["treatment"]].upper() in {
            item.upper() for item in row["sample_barcode_labels"].split(";")
        }, axis=1
    )
    if int(sample_counts.loc[sample_counts["sample"].eq(10115), "retained_cells"].iloc[0]) != 1:
        raise ValueError("The critical WT day-1 sample no longer has one retained cell")
    if bool(sample_counts.loc[sample_counts["sample"].eq(10115), "label_consistent"].iloc[0]):
        raise ValueError("Expected the documented treatment-label conflict for sample 10115")

    eligible_counts = sample_counts.groupby("timepoint", observed=True)["passes_min_cells"].sum().to_dict()
    endpoint_rows = [
        {
            "candidate_dataset": ACCESSION,
            "contrast": "day_1_vs_uninjured",
            "eligible_numerator_samples": int(eligible_counts.get("day_1", 0)),
            "eligible_denominator_samples": int(eligible_counts.get("uninjured", 0)),
        },
        {
            "candidate_dataset": ACCESSION,
            "contrast": "day_7_vs_uninjured",
            "eligible_numerator_samples": int(eligible_counts.get("day_7", 0)),
            "eligible_denominator_samples": int(eligible_counts.get("uninjured", 0)),
        },
        {
            "candidate_dataset": ACCESSION,
            "contrast": "direct_day_7_minus_day_1",
            "eligible_numerator_samples": int(eligible_counts.get("day_7", 0)),
            "eligible_denominator_samples": int(eligible_counts.get("day_1", 0)),
        },
    ]
    estimability = pd.DataFrame(endpoint_rows)
    estimability["estimable"] = (
        estimability["eligible_numerator_samples"].ge(MIN_REPLICATES)
        & estimability["eligible_denominator_samples"].ge(MIN_REPLICATES)
    )
    estimability["decision"] = estimability["estimable"].map({
        True: "eligible as a separate injury-model endpoint context only",
        False: "not modeled; no threshold relaxation or imputation",
    })
    if estimability.set_index("contrast").loc["direct_day_7_minus_day_1", "estimable"]:
        raise ValueError("Direct temporal change unexpectedly became estimable")

    sample_out = out / "GSE205029_wt_sample_eligibility.tsv"
    contrast_out = out / "GSE205029_contrast_estimability.tsv"
    sample_counts.to_csv(sample_out, sep="\t", index=False, lineterminator="\n")
    estimability.to_csv(contrast_out, sep="\t", index=False, lineterminator="\n")
    provenance = {
        "created_at": datetime.now().astimezone().isoformat(),
        "accession": ACCESSION,
        "status": "excluded_before_expression_modeling",
        "expression_counts_accessed": False,
        "eligible_samples_by_timepoint": eligible_counts,
        "direct_change_estimable": False,
        "organism": "Mus musculus",
        "injury_context": "dorsal column lesion; tissue around lesion core",
        "audit_scope": "eligibility before expression-count download or outcome analysis",
        "fixed_rules": {"minimum_retained_cells_per_sample": MIN_CELLS, "minimum_biological_samples_per_group": MIN_REPLICATES},
        "inputs": [
            {"path": str(path.relative_to(root)), "bytes": path.stat().st_size, "sha256": sha256(path)}
            for path in [metadata_path, barcodes_path, features_path, raw / "GSE205029_esearch.json", raw / "GSE205029_esummary.json", raw / "GSM6204337_brief.soft.txt"]
        ],
        "count_reconciliation": {"metadata_cells": len(metadata), "barcode_cells": len(barcodes), "features": len(features)},
        "critical_finding": "WT day-1 sample GSM6204337/#10115 retained one post-QC cell; its barcode file label says Ryk_WT_14 while GEO and meta_data identify WT day 1",
        "wt_sample_barcode_label_conflicts": int((~sample_counts["label_consistent"]).sum()),
        "final_role": "audited candidate; excluded from expression modeling and all inference",
        "outputs": [str(sample_out.relative_to(root)), str(contrast_out.relative_to(root))],
    }
    provenance_path = out / "GSE205029_eligibility_provenance.json"
    provenance_path.write_text(json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": "excluded_before_expression_modeling",
        "eligible_samples_by_timepoint": eligible_counts,
        "direct_change_estimable": False,
        "sample_table": str(sample_out),
        "contrast_table": str(contrast_out),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
