#!/usr/bin/env python3
"""Report figures for the M1 (MOp) and V1 (VISp) neuron-type expression analysis.

10x (whole transcriptome, within-region contrasts)
  fig1_10x_subclass_markers_heatmap        top replicated markers per MOp subclass
  fig2_10x_adenosine_dotplot               adenosine pathway genes by subclass, MOp and VISp
  fig3_10x_energy_module_heatmap           OXPHOS complexes and energy sets by subclass
  fig4_10x_pbm_candidates_heatmap          PBM same-direction DEG candidates by MOp subclass
MERFISH Zhuang-ABCA (four animals, MOp versus VISp paired by animal)
  fig5_merfish_paired_deg_counts           called regional DEG per subclass
  fig6_merfish_paired_deg_top_heatmap      top regional DEG per subclass with animal agreement
  fig7_merfish_adora1_penk_animals         per-animal paired differences for Adora1 and Penk

Every figure is written as PNG (180 dpi) and PDF under figures/report_m1v1/.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
WT = ROOT / "results/rnaseq/whole_transcriptome"
MK = WT / "subclass_markers"
ZD = ROOT / "results/zhuang/paired_deg"
FIG = ROOT / "figures/report_m1v1"
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "font.family": "DejaVu Sans"})

GLUT_ORDER = ["L2/3 IT CTX", "L4/5 IT CTX", "L5 IT CTX", "L6 IT CTX", "L5 PT CTX", "L5/6 NP CTX", "L6 CT CTX",
              "L6b CTX", "Car3"]
GABA_ORDER = ["Pvalb", "Sst", "Sst Chodl", "Vip", "Lamp5", "Sncg"]
ORDER_10X = GLUT_ORDER + GABA_ORDER


def save(fig, name):
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{name}.png", dpi=180, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def present(order, available):
    return [s for s in order if s in set(available)]


def load_pseudobulk():
    data = np.load(WT / "pseudobulk_counts.npz", allow_pickle=True)
    groups = pd.read_csv(WT / "pseudobulk_groups.csv").set_index("group").loc[data["groups"]]
    return data["sums"].astype(float), data["detected"].astype(float), pd.Index(data["genes"]), groups


def subclass_log2cpm(sums, genes, groups, region, gene_list):
    """Donor-averaged log2(CPM+1) per subclass for the given genes."""
    g = groups[groups.region == region]
    idx = genes.get_indexer(gene_list)
    out = {}
    for s, f in g.groupby("subclass"):
        vals = []
        for key, row in f.iterrows():
            c = groups.index.get_loc(key)
            vals.append(np.log2(sums[idx, c] / row.library_counts * 1e6 + 1))
        out[s] = np.mean(vals, axis=0)
    return pd.DataFrame(out, index=gene_list)


def fig1(sums, genes, groups):
    mk = pd.read_csv(MK / "MOp_subclass_markers.csv")
    rep = mk[mk.replicated_MOp_and_VISp]
    order = present(ORDER_10X, rep.subclass.unique())
    chosen = []
    for s in order:
        top = rep[rep.subclass == s].sort_values("log2FC", ascending=False)
        for gname in top.gene:
            if gname not in chosen:
                chosen.append(gname)
                if sum(1 for x in chosen if x in set(top.gene)) >= 5:
                    break
    mat = subclass_log2cpm(sums, genes, groups, "MOp", chosen)[order]
    z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, np.nan), axis=0)
    fig, ax = plt.subplots(figsize=(0.42 * len(order) + 2.2, 0.16 * len(chosen) + 1.6))
    im = ax.imshow(z.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    ax.set_xticks(range(len(order)), order, rotation=60, ha="right")
    ax.set_yticks(range(len(chosen)), chosen, fontsize=6.5, style="italic")
    ax.set_title("MOp neuron subclass markers (10x)\nreplicated in VISp donors from separate batches")
    fig.colorbar(im, ax=ax, shrink=0.4, label="row z-score of log2(CPM+1)")
    save(fig, "fig1_10x_subclass_markers_heatmap")
    return chosen


def dotplot(ax, mean, pct, title, vmax):
    genes_, subs = mean.index.tolist(), mean.columns.tolist()
    for i, gname in enumerate(genes_):
        for j, s in enumerate(subs):
            ax.scatter(j, i, s=8 + 260 * pct.loc[gname, s], c=[[mean.loc[gname, s]]], cmap="viridis",
                       vmin=0, vmax=vmax, edgecolors="0.3", linewidths=0.3)
    ax.set_xticks(range(len(subs)), subs, rotation=60, ha="right")
    ax.set_yticks(range(len(genes_)), genes_, style="italic")
    ax.set_xlim(-0.6, len(subs) - 0.4)
    ax.set_ylim(len(genes_) - 0.5, -0.5)
    ax.set_title(title)


def fig2(sums, det, genes, groups):
    sets = pd.read_csv(WT / "focus_gene_sets.csv")
    aden = sets[(sets.theme == "adenosine") & sets.present]
    order_sets = ["adenosine_receptor", "adenosine_production", "adenosine_clearance", "adenosine_transport",
                  "ATP_release"]
    gl = []
    for s in order_sets:
        gl += [x for x in aden[aden.set == s].gene if x not in gl]
    fig, axes = plt.subplots(1, 2, figsize=(11, 0.32 * len(gl) + 2.2), sharey=True)
    vmax = 0
    panels = []
    for region in ("MOp", "VISp"):
        mean = subclass_log2cpm(sums, genes, groups, region, gl)
        subs = present(ORDER_10X, mean.columns)
        mean = mean[subs]
        g = groups[groups.region == region]
        idx = genes.get_indexer(gl)
        pct = pd.DataFrame({s: det[np.ix_(idx, [groups.index.get_loc(k) for k in f.index])].sum(axis=1)
                            / f.n_cells.sum() for s, f in g.groupby("subclass")}, index=gl)[subs]
        panels.append((region, mean, pct))
        vmax = max(vmax, float(mean.to_numpy().max()))
    for ax, (region, mean, pct) in zip(axes, panels):
        dotplot(ax, mean, pct, f"{region} (10x, donor mean)", vmax)
    sm = plt.cm.ScalarMappable(cmap="viridis", norm=plt.Normalize(0, vmax))
    fig.colorbar(sm, ax=axes, shrink=0.5, label="log2(CPM+1)")
    for frac in (0.1, 0.4, 0.8):
        axes[1].scatter([], [], s=8 + 260 * frac, c="0.6", label=f"{int(frac * 100)}% cells")
    axes[1].legend(loc="upper left", bbox_to_anchor=(1.18, 1), frameon=False, title="detected")
    fig.suptitle("Adenosine pathway genes by neuron subclass. Regions are shown side by side only; "
                 "10x region contrasts are confounded with donor and batch", fontsize=8)
    save(fig, "fig2_10x_adenosine_dotplot")


def fig3():
    sc = pd.read_csv(MK / "focus_set_scores_by_donor.csv")
    sc = sc[sc.region == "MOp"]
    sets = ["OXPHOS_CI", "OXPHOS_CII", "OXPHOS_CIII", "OXPHOS_CIV", "OXPHOS_CV", "cytochrome_c",
            "OXPHOS_CIV_mt_encoded", "mito_biogenesis", "creatine_buffer", "glycolysis_lactate"]
    sets = [s for s in sets if s in set(sc.set)]
    subs = present(ORDER_10X, sc.subclass.unique())
    mean = sc.pivot_table(index="set", columns="subclass", values="log2_ratio_to_all_neurons").loc[sets, subs]
    agree = sc.groupby(["set", "subclass"]).log2_ratio_to_all_neurons.apply(
        lambda v: (v > 0).all() or (v < 0).all()).unstack().loc[sets, subs]
    lim = float(np.nanmax(np.abs(mean.to_numpy())))
    fig, ax = plt.subplots(figsize=(0.45 * len(subs) + 2.5, 0.38 * len(sets) + 1.6))
    im = ax.imshow(mean.to_numpy(), aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim))
    for i in range(len(sets)):
        for j in range(len(subs)):
            v = mean.iat[i, j]
            ax.text(j, i, f"{v:+.2f}{'' if agree.iat[i, j] else '*'}", ha="center", va="center", fontsize=5.5,
                    color="white" if abs(v) > 0.6 * lim else "black")
    ax.set_xticks(range(len(subs)), subs, rotation=60, ha="right")
    ax.set_yticks(range(len(sets)), sets)
    ax.set_title("MOp energy metabolism sets by neuron subclass (10x)\nlog2 share of library vs all MOp neurons, "
                 "3-donor mean. * donors disagree in sign")
    fig.colorbar(im, ax=ax, shrink=0.6, label="log2 ratio to all neurons")
    save(fig, "fig3_10x_energy_module_heatmap")


def fig4(sums, genes, groups):
    common = pd.read_csv(ROOT / "results/pbm/same_direction_52.csv", encoding="utf-8-sig")
    dcol = "hippocampus_direction"
    common = common[["gene", dcol]].rename(columns={dcol: "direction"})
    common = common[common.gene.isin(genes)]
    mat = subclass_log2cpm(sums, genes, groups, "MOp", common.gene.tolist())
    subs = present(ORDER_10X, mat.columns)
    mat = mat[subs]
    z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1).replace(0, np.nan), axis=0)
    common["peak"] = [subs[int(np.nanargmax(z.loc[g].to_numpy()))] if z.loc[g].notna().any() else ""
                      for g in common.gene]
    common["order"] = common.peak.map({s: i for i, s in enumerate(subs)})
    common = common.sort_values(["direction", "order"])
    z = z.loc[common.gene]
    fig, ax = plt.subplots(figsize=(0.45 * len(subs) + 2.6, 0.15 * len(z) + 1.6))
    im = ax.imshow(z.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5)
    labels = [f"{g} ({d})" for g, d in zip(common.gene, common.direction)]
    ax.set_yticks(range(len(z)), labels, fontsize=6)
    ax.set_xticks(range(len(subs)), subs, rotation=60, ha="right")
    ax.set_title("PBM same-direction DEG (hippocampus and cortex) across MOp neuron subclasses (10x)\n"
                 "row z-score of donor-mean log2(CPM+1); baseline atlas, not PBM response")
    fig.colorbar(im, ax=ax, shrink=0.4, label="row z-score")
    save(fig, "fig4_10x_pbm_candidates_heatmap")
    common.to_csv(MK / "pbm_candidates_peak_subclass_MOp.csv", index=False)


def fig5():
    s = pd.read_csv(ZD / "paired_region_deg_summary.csv")
    s = s.sort_values("cell_class", ascending=False)
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    y = np.arange(len(s))
    ax.barh(y, s.get("MOp_higher", 0), color="#1f5fa8", label="MOp higher")
    ax.barh(y, -s.get("VISp_higher", 0), color="#c2571a", label="VISp higher")
    ax.set_yticks(y, s.subclass_clean)
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("called genes (1,122-gene panel)")
    ax.set_title("MERFISH MOp vs VISp within subclass, paired by animal (Zhuang-ABCA-1 to 4)\n"
                 "rule: all animals same sign, |diff| >= 0.10 log2, BH q < 0.10")
    ax.legend(frameon=False, fontsize=7)
    xmax = max(int(s.get("MOp_higher", 0).max()), int(s.get("VISp_higher", 0).max())) + 8
    ax.set_xlim(-xmax, xmax)
    save(fig, "fig5_merfish_paired_deg_counts")


def fig6():
    c = pd.read_csv(ZD / "paired_region_deg_called.csv")
    subs = c.groupby("subclass_clean").size()
    subs = subs[subs >= 5].index.tolist()
    order = [s for s in ["L2/3 IT CTX Glut", "L4/5 IT CTX Glut", "L5 IT CTX Glut", "L6 IT CTX Glut",
                         "L5 ET CTX Glut", "L6 CT CTX Glut"] if s in subs]
    chosen = []
    for s in order:
        f = c[c.subclass_clean == s]
        for direction in ("MOp_higher", "VISp_higher"):
            top = f[f.direction == direction].reindex(
                f[f.direction == direction].mean_difference.abs().sort_values(ascending=False).index).gene.head(4)
            chosen += [g for g in top if g not in chosen]
    alld = pd.read_csv(ZD / "paired_region_deg_all.csv.gz")
    mat = alld.pivot_table(index="gene", columns="subclass_clean", values="mean_difference").reindex(chosen)[order]
    called = c.assign(v=1).pivot_table(index="gene", columns="subclass_clean", values="v").reindex(chosen)
    called = called.reindex(columns=order)
    lim = float(np.nanmax(np.abs(mat.to_numpy())))
    fig, ax = plt.subplots(figsize=(0.7 * len(order) + 2.4, 0.17 * len(chosen) + 1.6))
    im = ax.imshow(mat.to_numpy(), aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim))
    for i in range(len(chosen)):
        for j in range(len(order)):
            if called.iat[i, j] == 1:
                ax.text(j, i, "●", ha="center", va="center", fontsize=5, color="black")
    ax.set_xticks(range(len(order)), order, rotation=45, ha="right")
    ax.set_yticks(range(len(chosen)), chosen, fontsize=6.5, style="italic")
    ax.set_title("Top MERFISH regional DEG per subclass (4-animal mean MOp minus VISp)\n"
                 "dot = called in that subclass")
    fig.colorbar(im, ax=ax, shrink=0.5, label="MOp minus VISp (log2 units)")
    save(fig, "fig6_merfish_paired_deg_top_heatmap")


def fig7():
    alld = pd.read_csv(ZD / "paired_region_deg_all.csv.gz")
    animals = [c for c in alld.columns if c.startswith("diff_")]
    order = ["L2/3 IT CTX Glut", "L4/5 IT CTX Glut", "L5 IT CTX Glut", "L6 IT CTX Glut", "L5 ET CTX Glut",
             "L6 CT CTX Glut", "Pvalb Gaba", "Sst Gaba", "Vip Gaba", "Lamp5 Gaba"]
    colors = ["#1f5fa8", "#5aa0d8", "#c2571a", "#e8a06a"]
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.8), sharey=True)
    for ax, gene in zip(axes, ["Adora1", "Penk"]):
        f = alld[alld.gene == gene].set_index("subclass_clean").reindex(order)
        for k, (a, col) in enumerate(zip(animals, colors)):
            ax.scatter(f[a], np.arange(len(order)) + (k - 1.5) * 0.12, color=col, s=16,
                       label=a.replace("diff_", ""))
        ax.scatter(f.mean_difference, np.arange(len(order)), marker="|", s=180, color="black")
        for i, (_, row) in enumerate(f.iterrows()):
            if row.get("called") is True:
                ax.annotate("called", (row.mean_difference, i), xytext=(0, 6), textcoords="offset points",
                            ha="center", fontsize=5.5)
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_yticks(range(len(order)), order)
        ax.invert_yaxis()
        ax.set_title(gene)
        ax.set_xlabel("MOp minus VISp (log2 units), dot = animal, bar = mean")
        ax.grid(axis="x", alpha=0.25)
    axes[1].legend(fontsize=6.5, frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1))
    fig.suptitle("MERFISH per-animal regional differences for Adora1 and Penk ('called' = paired DEG rule met)",
                 fontsize=9)
    save(fig, "fig7_merfish_adora1_penk_animals")


def main():
    sums, det, genes, groups = load_pseudobulk()
    fig1(sums, genes, groups)
    fig2(sums, det, genes, groups)
    fig3()
    fig4(sums, genes, groups)
    fig5()
    fig6()
    fig7()
    print(sorted(p.name for p in FIG.glob("*.png")))


if __name__ == "__main__":
    main()
