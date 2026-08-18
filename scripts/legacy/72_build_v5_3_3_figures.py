#!/usr/bin/env python3
"""Build reference-guided v5.3.3 main and supplementary figures."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from datetime import datetime
from pathlib import Path


PACKAGE = "submission_Molecular_Neurobiology_v5_3_3_postlock_upgrade_2026_08_03"
PACKAGE_V534 = "submission_Molecular_Neurobiology_v5_3_4_editorial_journal_compliance_2026_08_04"


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def p_label(value: float, exact: bool) -> str:
    if exact:
        return f"P = {value:.4f}"
    return "P < 0.001" if value < 0.001 else f"P = {value:.3f}"


def figure2(m, root: Path, out: Path, inputs: list[Path]) -> list[Path]:
    comp = m.read_tsv(root, "reports/phase_reproducibility_calibration_2026_07/gene_program_metric_comparison.tsv", inputs)
    null = m.read_tsv(root, "reports/phase_v5_3_3_postlock_upgrade_2026_08/sign_identity/feature_identity_permutation_summary.tsv", inputs)
    external = m.read_tsv(root, "reports/phase_v5_3_3_postlock_upgrade_2026_08/GSE47681/GSE47681_frozen_reference_comparison.tsv", inputs)

    fig = m.plt.figure(figsize=(m.FULL_WIDTH_IN, 7.35), layout="constrained")
    gs = fig.add_gridspec(4, 3, height_ratios=[1.00, 1.00, 1.00, 0.17], hspace=0.10)
    a_axes = [fig.add_subplot(gs[0, j]) for j in range(3)]
    b_axes = [fig.add_subplot(gs[1, j]) for j in range(3)]
    c_axes = [fig.add_subplot(gs[2, j]) for j in range(3)]
    legend_ax = fig.add_subplot(gs[3, :]); legend_ax.axis("off")

    metrics_a = ["all_study_sign_concordance_rate", "mean_held_out_direction_accuracy", "mean_pairwise_spearman"]
    labels_a = ["all-study direction", "held-out direction", "rank correlation"]
    for j, (ax, contrast) in enumerate(zip(a_axes, m.CONTRASTS)):
        z = comp[(comp["contrast_id"] == contrast) & comp["metric"].isin(metrics_a)].set_index("metric").reindex(metrics_a)
        y = m.np.arange(len(metrics_a))
        for yi, row in zip(y, z.itertuples()):
            ax.plot([row.gene, row.program], [yi, yi], color="#B5B5B5", lw=1.1, zorder=1)
        ax.scatter(z["gene"], y, s=34, facecolor="white", edgecolor=m.MID, linewidth=1.0, marker="o", zorder=3)
        ax.scatter(z["program"], y, s=38, facecolor=m.BLUE, edgecolor=m.BLUE, marker="D", zorder=4)
        ax.set_xlim(0, 1); ax.set_ylim(2.35, -0.35); ax.set_xticks(m.np.arange(0, 1.01, 0.2))
        ax.set_title(m.CONTRAST_LABEL[contrast], pad=4); ax.set_xlabel("metric value")
        ax.set_yticks(y, labels_a if j == 0 else [])
        if j == 0: m.panel_label(ax, "a", x=-0.30)
        m.clean_axis(ax, "x")

    metrics_b = ["all_study_sign_concordance_rate", "mean_held_out_balanced_accuracy", "mean_held_out_mcc"]
    labels_b = ["all-study direction", "balanced accuracy", "Matthews correlation"]
    for j, (ax, contrast) in enumerate(zip(b_axes, m.CONTRASTS)):
        z = null[(null["scale"] == "program") & (null["contrast_id"] == contrast) & null["metric"].isin(metrics_b)].set_index("metric").reindex(metrics_b)
        y = m.np.arange(len(metrics_b))
        low = z["null_q025"].to_numpy(float); med = z["null_median"].to_numpy(float); high = z["null_q975"].to_numpy(float); observed = z["observed"].to_numpy(float)
        ax.hlines(y, low, high, color=m.MID, lw=2.5, zorder=1)
        ax.scatter(med, y, s=26, facecolor="white", edgecolor=m.MID, marker="o", zorder=3)
        ax.scatter(observed, y, s=40, color=m.BLUE, marker="D", zorder=4)
        for yi, p, obs in zip(y, z["empirical_p_null_ge_observed"], observed):
            label = p_label(float(p), bool(getattr(m, "EXACT_EMPIRICAL_P", False)))
            ax.text(min(0.985, obs + 0.035), yi - 0.18, label, fontsize=6.4, ha="right" if obs > 0.85 else "left", va="bottom")
        ax.set_xlim(-0.15, 1.0); ax.set_ylim(2.35, -0.35); ax.set_title(m.CONTRAST_LABEL[contrast], pad=4); ax.set_xlabel("feature-identity calibration")
        ax.set_yticks(y, labels_b if j == 0 else [])
        if j == 0: m.panel_label(ax, "b", x=-0.30)
        m.clean_axis(ax, "x")

    metrics_c = ["spearman_vs_frozen_mean", "balanced_direction_accuracy", "conditional_accuracy_given_unanimous_frozen"]
    labels_c = ["rank correlation", "balanced accuracy", "unanimity-gated accuracy"]
    for j, (ax, contrast) in enumerate(zip(c_axes, m.CONTRASTS)):
        q = external[external["contrast_id"] == contrast].set_index("scale")
        gene = q.loc["gene", metrics_c].astype(float).to_numpy(); hallmark = q.loc["hallmark", metrics_c].astype(float).to_numpy(); y = m.np.arange(len(metrics_c))
        for yi, left, right in zip(y, gene, hallmark):
            ax.plot([left, right], [yi, yi], color="#B5B5B5", lw=1.1, zorder=1)
        ax.scatter(gene, y, s=34, facecolor="white", edgecolor=m.MID, linewidth=1.0, marker="o", zorder=3)
        ax.scatter(hallmark, y, s=38, facecolor=m.TEAL, edgecolor=m.TEAL, marker="D", zorder=4)
        ax.set_xlim(0.25, 1.0); ax.set_ylim(2.35, -0.35); ax.set_xticks(m.np.arange(0.25, 1.01, 0.25))
        ax.set_title(m.CONTRAST_LABEL[contrast], pad=4); ax.set_xlabel("GSE47681 transfer")
        ax.set_yticks(y, labels_c if j == 0 else [])
        if j == 0: m.panel_label(ax, "c", x=-0.30)
        m.clean_axis(ax, "x")

    panel_refs = "a,c" if bool(getattr(m, "LOWERCASE_PANEL_REFERENCES", False)) else "A,C"
    null_ref = "b" if bool(getattr(m, "LOWERCASE_PANEL_REFERENCES", False)) else "B"
    handles = [
        m.plt.Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor=m.MID, markersize=5.5, label=f"gene observed ({panel_refs}) / null median ({null_ref})"),
        m.plt.Line2D([0], [0], marker="D", color="none", markerfacecolor=m.BLUE, markeredgecolor=m.BLUE, markersize=5.5, label="Hallmark observed"),
        m.plt.Line2D([0], [0], color=m.MID, lw=2.5, label="identity-permutation 95% interval"),
        m.plt.Line2D([0], [0], marker="D", color="none", markerfacecolor=m.TEAL, markeredgecolor=m.TEAL, markersize=5.5, label="GSE47681 Hallmark"),
    ]
    legend_ax.legend(handles=handles, frameon=False, ncol=4, loc="center", handletextpad=0.45, columnspacing=1.0, fontsize=7.0)
    return m.save_figure(fig, out, "Fig2_GeneProgramTransfer")


def build_estimability_data(root: Path):
    m51 = load_module(root / "scripts" / "51_molecular_neurobiology_v5_3_figures.py", "m51_estimability_data")
    m = m51.load_base(root)
    datasets = ["GSE162610", "GSE234774", "GSE304399", "GSE304361", "GSE205029", "GSE159638", "GSE172167", "GSE240727", "GSE247844", "GSE275982", "GSE298545"]
    stages = ["exact_d1_d7_coverage", "independent_replication", "sample_mapping", "retained_cell_eligibility", "lesion_site_compatibility", "direct_change_estimable"]
    status = {
        "GSE162610": ["pass"]*6, "GSE234774": ["pass"]*6, "GSE304399": ["pass"]*6,
        "GSE304361": ["fail", "pass", "pass", "pass", "pass", "fail"],
        "GSE205029": ["pass", "pass", "qualified", "fail", "pass", "fail"],
        "GSE159638": ["fail", "pass", "pass", "pass", "pass", "fail"],
        "GSE172167": ["pass", "pass", "pass", "pass", "fail", "fail"],
        "GSE240727": ["pass", "fail", "qualified", "pass", "pass", "fail"],
        "GSE247844": ["fail", "pass", "pass", "pass", "pass", "fail"],
        "GSE275982": ["fail", "pass", "pass", "pass", "qualified", "fail"],
        "GSE298545": ["fail", "qualified", "qualified", "pass", "pass", "fail"],
    }
    rows = []
    for dataset in datasets:
        rows.append({"dataset": dataset, **dict(zip(stages, status[dataset]))})
    return m, m.pd.DataFrame(rows), stages


def figure_s2(m, frame, stages, root: Path, out: Path, inputs: list[Path]) -> list[Path]:
    data_path = root / "reports" / "phase_v5_3_3_postlock_upgrade_2026_08" / "estimability_flow_matrix.tsv"
    frame.to_csv(data_path, sep="\t", index=False, lineterminator="\n")
    inputs.append(data_path)
    flow_counts = [11, 6, 5, 5, 4, 3, 3]
    flow_labels = ["candidate\ndatasets", "exact d1+d7\ncoverage", "independent\nreplication", "resolvable or\nqualified mapping", "retained-cell\neligibility", "lesion-site\ncompatible", "between-timepoint\nestimable"]
    fig = m.plt.figure(figsize=(m.FULL_WIDTH_IN, 5.10), layout="constrained")
    gs = fig.add_gridspec(2, 1, height_ratios=[0.72, 1.35])
    ax1 = fig.add_subplot(gs[0, 0]); ax2 = fig.add_subplot(gs[1, 0])
    x = m.np.arange(len(flow_counts))
    ax1.plot(x, flow_counts, color=m.BLUE, lw=1.7, marker="o", ms=5.5)
    for xi, count in zip(x, flow_counts): ax1.text(xi, count + 0.45, str(count), ha="center", va="bottom", fontsize=7.5, fontweight="bold")
    ax1.set_xticks(x, flow_labels); ax1.set_ylim(0, 12.5); ax1.set_ylabel("datasets retained"); m.clean_axis(ax1, "y"); m.panel_label(ax1, "a", x=-0.08)
    ax1.text(0.99, 0.05, "Sequential counts retain qualified mapping; later criteria can still render a dataset ineligible.", transform=ax1.transAxes, ha="right", va="bottom", fontsize=6.3, color=m.MID)

    display = ["d1+d7", "replication", "mapping", "retained cells", "lesion site", "between-timepoint"]
    color = {"pass": "#DCEFE4", "qualified": "#F7E6B5", "fail": "#F4D6D6"}
    label = {"pass": "pass", "qualified": "qualified", "fail": "not met"}
    for i, row in frame.iterrows():
        for j, stage in enumerate(stages):
            value = row[stage]
            ax2.add_patch(m.Rectangle((j-0.48, i-0.39), 0.96, 0.78, facecolor=color[value], edgecolor="white", linewidth=0.8))
            ax2.text(j, i, label[value], ha="center", va="center", fontsize=6.1, color=m.DARK, fontweight="bold" if value == "pass" else "normal")
    ax2.set_xlim(-0.52, len(stages)-0.48); ax2.set_ylim(len(frame)-0.5, -0.5); ax2.set_xticks(range(len(stages)), display); ax2.xaxis.tick_top(); ax2.set_yticks(range(len(frame)), frame["dataset"]); ax2.tick_params(length=0)
    for spine in ax2.spines.values(): spine.set_visible(False)
    m.panel_label(ax2, "b", x=-0.10, y=1.10)
    return m.save_figure(fig, out, "FigS2_EstimabilityFlow")


def figure_s3(m, root: Path, out: Path, inputs: list[Path]) -> list[Path]:
    diag = m.read_tsv(root, "reports/phase_reproducibility_calibration_2026_07/equal_feature_count_diagnostic.tsv", inputs)
    metrics = ["all_study_sign_concordance_rate", "mean_pairwise_sign_concordance", "mean_held_out_direction_accuracy"]
    labels = ["all-study direction", "pairwise direction", "held-out direction"]
    fig, axes = m.plt.subplots(1, 3, figsize=(m.FULL_WIDTH_IN, 2.65), layout="constrained")
    for j, (ax, contrast) in enumerate(zip(axes, m.CONTRASTS)):
        z = diag[(diag["contrast_id"] == contrast) & diag["metric"].isin(metrics)].set_index("metric").reindex(metrics)
        x = m.np.arange(len(metrics)); low=z["random_gene_subset_q025"].astype(float); med=z["random_gene_subset_median"].astype(float); high=z["random_gene_subset_q975"].astype(float); obs=z["program_observed"].astype(float)
        ax.vlines(x, low, high, color=m.MID, lw=2.2); ax.scatter(x, med, s=25, facecolor="white", edgecolor=m.MID); ax.scatter(x, obs, s=38, color=m.BLUE, marker="D")
        for xi, p, y, upper in zip(x, z["empirical_p_program_ge_random_gene_features"], obs, high):
            txt=p_label(float(p), bool(getattr(m, "EXACT_EMPIRICAL_P", False)))
            xpos = xi + 0.04 if xi == 0 else xi
            ax.text(xpos, min(1.0, max(float(y), float(upper))+0.055), txt, ha="left" if xi == 0 else "center", va="bottom", fontsize=6.2)
        ax.set_ylim(0.2,1.05); ax.set_xticks(x, labels, rotation=18, ha="right"); ax.set_title(m.CONTRAST_LABEL[contrast]);
        if j==0: ax.set_ylabel("metric value")
        m.panel_label(ax, chr(ord("a") + j), x=-0.28 if j == 0 else -0.18)
        m.clean_axis(ax,"y")
    return m.save_figure(fig, out, "FigS3_EqualFeatureDiagnostic")


def figure_s4(m, root: Path, out: Path, inputs: list[Path]) -> list[Path]:
    pca = m.read_tsv(root, "reports/phase_v5_3_3_postlock_upgrade_2026_08/GSE47681/GSE47681_RMA_PCA_scores.tsv", inputs)
    var = m.read_tsv(root, "reports/phase_v5_3_3_postlock_upgrade_2026_08/GSE47681/GSE47681_RMA_PCA_variance.tsv", inputs)
    rle = m.read_tsv(root, "reports/phase_v5_3_3_postlock_upgrade_2026_08/GSE47681/GSE47681_RLE_diagnostics.tsv", inputs)
    fig, (ax1, ax2) = m.plt.subplots(1, 2, figsize=(m.FULL_WIDTH_IN, 2.85), layout="constrained", gridspec_kw={"width_ratios":[0.9,1.3]})
    colors = {"SHAMWT": m.MID, "DAY1WT": m.ORANGE, "DAY7WT": m.TEAL}; markers={"SHAMWT":"o","DAY1WT":"s","DAY7WT":"D"}
    for group in ["SHAMWT","DAY1WT","DAY7WT"]:
        q=pca[pca["group"]==group]; ax1.scatter(q["PC1"],q["PC2"],s=42,color=colors[group],marker=markers[group],edgecolor="white",linewidth=0.5,label=group.replace("WT"," WT"))
    ax1.set_xlabel(f"PC1 ({100*float(var.iloc[0]['variance_explained']):.1f}%)"); ax1.set_ylabel(f"PC2 ({100*float(var.iloc[1]['variance_explained']):.1f}%)"); ax1.legend(frameon=False,fontsize=6.8); m.clean_axis(ax1,"both"); m.panel_label(ax1,"a",x=-0.16)
    ordered_pca = pca.assign(group_order=m.pd.Categorical(pca["group"], categories=["SHAMWT", "DAY1WT", "DAY7WT"], ordered=True))
    order=ordered_pca.sort_values(["group_order","gsm"])["gsm"].tolist(); q=rle.set_index("gsm").loc[order].reset_index(); x=m.np.arange(len(q)); lower=q["rle_median"]-q["rle_q25"]; upper=q["rle_q75"]-q["rle_median"]
    group_map=pca.set_index("gsm")["group"].to_dict()
    for i,row in q.iterrows():
        group=group_map[row.gsm]; ax2.errorbar(i,row.rle_median,yerr=[[float(lower.iloc[i])],[float(upper.iloc[i])]],fmt=markers[group],ms=4,color=colors[group],ecolor=colors[group],capsize=2,lw=1)
    ax2.axhline(0,color=m.DARK,lw=0.7); ax2.set_xticks(x,[gsm.replace("GSM115","") for gsm in order],rotation=45,ha="right"); ax2.set_ylabel("relative log expression"); ax2.set_xlabel("GSM suffix"); m.clean_axis(ax2,"y"); m.panel_label(ax2,"b",x=-0.10)
    return m.save_figure(fig,out,"FigS4_GSE47681_QC")


def main() -> None:
    parser=argparse.ArgumentParser(); parser.add_argument("--root",type=Path,default=Path(__file__).resolve().parents[1]); parser.add_argument("--out-dir",type=Path,default=None); parser.add_argument("--journal-compliance-v534",action="store_true"); args=parser.parse_args(); root=args.root.resolve(); default_package=PACKAGE_V534 if args.journal_compliance_v534 else PACKAGE; out_arg=args.out_dir if args.out_dir is not None else Path(default_package)/"Figures"; out=out_arg if out_arg.is_absolute() else root/out_arg; out.mkdir(parents=True,exist_ok=True)
    m51=load_module(root/"scripts"/"51_molecular_neurobiology_v5_3_figures.py","m51_v533"); m=m51.load_base(root); m.set_style(); inputs=[]; outputs=[]
    # The legacy branch keeps the prior uppercase presentation. The v5.3.4
    # journal-compliance branch follows Springer multi-panel lowercase grammar.
    original_panel_label = m.panel_label
    if not args.journal_compliance_v534:
        m.panel_label = lambda ax, label, x=-0.12, y=1.08: original_panel_label(ax, str(label).upper(), x=x, y=y)
    else:
        m.plt.rcParams.update({"font.size": 8.5, "axes.labelsize": 8.5, "axes.titlesize": 9.0, "xtick.labelsize": 8.0, "ytick.labelsize": 8.0, "legend.fontsize": 8.0})
        m.SHOW_INLINE_PANEL_TITLES = False
        m.NOTE_FONTSIZE = 7.6
        m.STABILITY_LABELS = ["conditionally direction-stable", "moderately stable", "unstable", "not evaluable"]
        m.FIG6_COLUMNS = ["d1", "d7", "day 7 vs\nday 1", "d1", "direction\nretained", "d7"]
        m.FIG1_CONTRAST_LABELS = ["day 1", "day 7", "day 7 vs day 1"]
        m.LOWERCASE_PANEL_REFERENCES = True
        m.ACCESSIBLE_HATCH = True
        m.EXACT_EMPIRICAL_P = True
        m.FIG4_HEIGHT = 5.0
        m.FIG4_HEIGHT_RATIOS = [0.85, 1.0, 0.30]
        m.FIG4_NOTE_Y = 0.10
    m.CONTRAST_LABEL[m.CONTRASTS[0]] = "day 1 vs uninjured"
    m.CONTRAST_LABEL[m.CONTRASTS[1]] = "day 7 vs uninjured"
    m.CONTRAST_LABEL[m.CONTRASTS[2]] = "day 7 vs day 1" if args.journal_compliance_v534 else "direct day 7 - day 1"
    outputs += m51.figure1(m,root,out,inputs); outputs += figure2(m,root,out,inputs); outputs += m51.figure3(m,root,out,inputs); outputs += m51.figure4(m,root,out,inputs); outputs += m51.figure5(m,root,out,inputs); outputs += m51.figure6(m,root,out,inputs); outputs += m.figure_s1(root,out,inputs)
    _, flow, stages = build_estimability_data(root); outputs += figure_s2(m,flow,stages,root,out,inputs); outputs += figure_s3(m,root,out,inputs); outputs += figure_s4(m,root,out,inputs)
    version = "v5.3.4" if args.journal_compliance_v534 else "v5.3.3"
    state = "presentation-only journal-compliance revision using frozen results; includes post-lock feature-identity calibration and orthogonal stress test" if args.journal_compliance_v534 else "main figures integrate pre-registered post-lock calibration and orthogonal stress test"
    provenance={"analysis":f"Molecular Neurobiology {version} reference-guided figures","created_at":datetime.now().astimezone().isoformat(),"scientific_state":state,"inputs":[{"path":str(p.relative_to(root)).replace('\\','/'),"bytes":p.stat().st_size,"sha256":sha256(p)} for p in sorted(set(inputs))],"outputs":[{"path":str(p.relative_to(root)).replace('\\','/'),"bytes":p.stat().st_size,"sha256":sha256(p)} for p in sorted(outputs)]}
    path=out/"figure_provenance.json"; path.write_text(json.dumps(provenance,ensure_ascii=False,indent=2)+"\n",encoding="utf-8"); print({"figures":10,"files":len(outputs),"provenance":str(path)})


if __name__ == "__main__":
    main()
