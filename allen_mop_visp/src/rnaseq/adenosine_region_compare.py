#!/usr/bin/env python3
"""Adenosine pathway MOp versus VISp difference by neuron subclass, checked against MERFISH.

Problem
-------
In the 10x data MOp and VISp come from different donors and library batches, so a raw
MOp minus VISp difference mixes region with batch. Zhuang MERFISH has both regions in the
same sections of four animals, so its paired difference is free of that confound, but its
panel holds only Adora1 and Adora2a from the adenosine pathway.

Approach
--------
1. Validation. For every gene shared by both datasets and expressed in the subclass, compare
   the 10x difference with the MERFISH four-animal paired difference. High sign agreement
   in a subclass means 10x regional differences in that subclass are mostly biological.
2. Show the 10x difference for all adenosine genes (Fig2b), labelling each subclass with its
   validation agreement and masking genes too low to compare.
3. Show Adora1 in both datasets side by side (Fig2c).
4. Supplementary scatter of the validation per subclass (Fig S6).

Units: 10x = difference of donor-mean log2(CPM+1); MERFISH = difference of mean official log2
expression. The scales differ, so agreement is judged on sign and correlation, not size.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import TwoSlopeNorm  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
WT = ROOT / "results/rnaseq/whole_transcriptome"
OUT = WT / "adenosine_region"
MF = ROOT / "results/zhuang/paired_deg/paired_region_deg_all.csv.gz"
FIG = ROOT / "figures/report_m1v1"
MIN_CELLS = 20
MIN_CPM = 5.0          # mask a gene in a subclass when both regions are below this mean CPM
EXPR_FOR_VALIDATION = 1.0   # log2(CPM+1) mean over regions, genes used in the validation
BIG = 0.10             # MERFISH |difference| counted as a clear regional difference
MAP = {"L2/3 IT CTX": "L2/3 IT CTX Glut", "L4/5 IT CTX": "L4/5 IT CTX Glut", "L5 IT CTX": "L5 IT CTX Glut",
       "L6 IT CTX": "L6 IT CTX Glut", "L5 PT CTX": "L5 ET CTX Glut", "L5/6 NP CTX": "L5 NP CTX Glut",
       "L6 CT CTX": "L6 CT CTX Glut", "L6b CTX": "L6b CTX Glut", "Pvalb": "Pvalb Gaba", "Sst": "Sst Gaba",
       "Vip": "Vip Gaba", "Lamp5": "Lamp5 Gaba", "Sncg": "Sncg Gaba"}
ORDER = ["L2/3 IT CTX", "L4/5 IT CTX", "L5 IT CTX", "L6 IT CTX", "L5 PT CTX", "L5/6 NP CTX", "L6 CT CTX",
         "L6b CTX", "Car3", "Pvalb", "Sst", "Vip", "Lamp5", "Sncg"]
SET_ORDER = ["adenosine_receptor", "adenosine_production", "adenosine_clearance", "adenosine_transport",
             "ATP_release"]
SET_LABEL = {"adenosine_receptor": "receptor", "adenosine_production": "production",
             "adenosine_clearance": "clearance", "adenosine_transport": "transport", "ATP_release": "ATP release"}
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "font.family": "DejaVu Sans"})


def save(fig, name):
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{name}.png", dpi=180, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def load():
    data = np.load(WT / "pseudobulk_counts.npz", allow_pickle=True)
    groups = pd.read_csv(WT / "pseudobulk_groups.csv").set_index("group").loc[data["groups"]]
    return data["sums"].astype(np.float64), pd.Index(data["genes"]), groups


def donor_log2cpm(sums, genes, groups):
    """Long table region, donor, subclass, gene index -> log2(CPM+1) and CPM, eligible groups only."""
    out = {}
    for (region, subclass), f in groups.groupby(["region", "subclass"]):
        if f.donor.nunique() < 3 or f.n_cells.min() < MIN_CELLS:
            continue
        cpm = np.stack([sums[:, groups.index.get_loc(k)] / r.library_counts * 1e6 for k, r in f.iterrows()])
        out[(region, subclass)] = cpm
    return out


def validation(cpm, genes, mf):
    rows, points = [], []
    for s10, smf in MAP.items():
        if ("MOp", s10) not in cpm or ("VISp", s10) not in cpm:
            continue
        mop = np.log2(cpm[("MOp", s10)] + 1).mean(axis=0)
        vis = np.log2(cpm[("VISp", s10)] + 1).mean(axis=0)
        diff = pd.Series(mop - vis, index=genes)
        expr = pd.Series((mop + vis) / 2, index=genes)
        m = mf[mf.subclass_clean == smf].set_index("gene")
        if m.empty:
            continue
        common = m.index.intersection(expr.index[expr > EXPR_FOR_VALIDATION])
        x, y = diff[common], m.loc[common, "mean_difference"]
        big = y.abs() >= BIG
        called = m.loc[common, "called"].astype(bool)
        rows.append({"subclass": s10, "merfish_subclass": smf, "n_genes": len(common),
                     "pearson_r": float(np.corrcoef(x, y)[0, 1]),
                     "sign_agreement_clear": float((np.sign(x[big]) == np.sign(y[big])).mean()),
                     "n_clear": int(big.sum()),
                     "sign_agreement_called": float((np.sign(x[called]) == np.sign(y[called])).mean())
                     if called.any() else np.nan, "n_called": int(called.sum())})
        points.append(pd.DataFrame({"subclass": s10, "gene": common, "tenx_diff": x.to_numpy(),
                                    "merfish_diff": y.to_numpy(), "merfish_called": called.to_numpy()}))
    return pd.DataFrame(rows), pd.concat(points, ignore_index=True)


def adenosine_table(cpm, genes, sets):
    aden = sets[(sets.theme == "adenosine") & sets.present]
    order = []
    for s in SET_ORDER:
        order += [(g, s) for g in aden[aden.set == s].gene if g not in [o[0] for o in order]]
    rows = []
    for gene, s in order:
        j = genes.get_loc(gene)
        for sub in ORDER:
            if ("MOp", sub) not in cpm or ("VISp", sub) not in cpm:
                continue
            m, v = cpm[("MOp", sub)][:, j], cpm[("VISp", sub)][:, j]
            rows.append({"gene": gene, "set": s, "subclass": sub, "cpm_MOp": m.mean(), "cpm_VISp": v.mean(),
                         "diff_log2cpm": np.log2(m + 1).mean() - np.log2(v + 1).mean(),
                         "min_MOp_donor": np.log2(m + 1).min(), "max_MOp_donor": np.log2(m + 1).max(),
                         "min_VISp_donor": np.log2(v + 1).min(), "max_VISp_donor": np.log2(v + 1).max()})
    t = pd.DataFrame(rows)
    # Donor ranges do not overlap: all three MOp donors above (or below) all three VISp donors.
    t["donor_ranges_separate"] = (t.min_MOp_donor > t.max_VISp_donor) | (t.max_MOp_donor < t.min_VISp_donor)
    t["low_expression"] = np.maximum(t.cpm_MOp, t.cpm_VISp) < MIN_CPM
    return t, [g for g, _ in order], dict(order)


def fig2b(table, gene_order, gene_set, val):
    subs = [s for s in ORDER if s in set(table.subclass)]
    diff = table.pivot_table(index="gene", columns="subclass", values="diff_log2cpm").reindex(index=gene_order,
                                                                                            columns=subs)
    low = table.pivot_table(index="gene", columns="subclass", values="low_expression").reindex(index=gene_order,
                                                                                             columns=subs)
    sep = table.pivot_table(index="gene", columns="subclass", values="donor_ranges_separate").reindex(
        index=gene_order, columns=subs)
    shown = diff.where(low != 1)
    lim = 1.5
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#ececec")
    fig, ax = plt.subplots(figsize=(0.62 * len(subs) + 3.2, 0.33 * len(gene_order) + 2.6))
    ax.imshow(np.ma.masked_invalid(shown.to_numpy().clip(-lim, lim)), aspect="auto", cmap=cmap,
              norm=TwoSlopeNorm(0, -lim, lim))
    for i, gene in enumerate(gene_order):
        for j, s in enumerate(subs):
            if low.iat[i, j] == 1:
                ax.text(j, i, "low", ha="center", va="center", fontsize=5, color="#9a9a9a")
                continue
            v = diff.iat[i, j]
            mark = "†" if sep.iat[i, j] == 1 else ""
            ax.text(j, i, f"{v:+.2f}{mark}", ha="center", va="center", fontsize=5.4,
                    color="white" if abs(v) > 0.9 else "#222222")
    # set separators and labels
    prev = None
    for i, gene in enumerate(gene_order):
        if prev is not None and gene_set[gene] != prev:
            ax.axhline(i - 0.5, color="white", linewidth=2.5)
        prev = gene_set[gene]
    for s in SET_ORDER:
        idx = [i for i, g in enumerate(gene_order) if gene_set[g] == s]
        if idx:
            ax.text(len(subs) - 0.3, np.mean(idx), SET_LABEL[s], va="center", ha="left", fontsize=7,
                    color="#444444")
    v = val.set_index("subclass").reindex(subs)
    labels = []
    for s in subs:
        a = v.loc[s, "sign_agreement_clear"] if s in v.index else np.nan
        labels.append(f"{s}\n{a * 100:.0f}%" if np.isfinite(a) else f"{s}\nn/a")
    ax.set_xticks(range(len(subs)), labels, rotation=60, ha="right")
    for tick, s in zip(ax.get_xticklabels(), subs):
        a = v.loc[s, "sign_agreement_clear"] if s in v.index else np.nan
        if not np.isfinite(a) or a < 0.85:
            tick.set_color("#b05a00")
    ax.set_yticks(range(len(gene_order)), gene_order, style="italic")
    gaba_start = subs.index("Pvalb") if "Pvalb" in subs else None
    if gaba_start:
        ax.axvline(gaba_start - 0.5, color="black", linewidth=0.8)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=TwoSlopeNorm(0, -lim, lim))
    fig.colorbar(sm, ax=ax, shrink=0.5, pad=0.12, label="MOp minus VISp, log2(CPM+1)\nred = higher in MOp")
    ax.set_title("Adenosine pathway, MOp minus VISp within each neuron subclass (10x)\n"
                 "% under subclass = sign agreement of 10x regional differences with MERFISH 4-animal paired "
                 "differences (shared genes).\nOrange = under 85%, read with caution. † all 3 MOp donors "
                 f"beyond all 3 VISp donors. low = under {MIN_CPM:g} CPM in both regions", fontsize=7.5)
    save(fig, "fig2b_adenosine_region_difference")


def fig2c(table, mf):
    a = table[table.gene == "Adora1"].set_index("subclass")
    subs = [s for s in ORDER if s in a.index and s in MAP]
    animals = [c for c in mf.columns if c.startswith("diff_")]
    m = mf[mf.gene == "Adora1"].set_index("subclass_clean")
    colors = ["#1f5fa8", "#5aa0d8", "#c2571a", "#e8a06a"]
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2), sharey=True, gridspec_kw={"wspace": 0.08})
    y = np.arange(len(subs))
    # Left: 10x absolute level in each region, donor range and mean.
    ax = axes[0]
    for region, col, off in (("MOp", "#1f5fa8", -0.13), ("VISp", "#c2571a", 0.13)):
        lo = a.loc[subs, f"min_{region}_donor"]
        hi = a.loc[subs, f"max_{region}_donor"]
        mean = np.log2(a.loc[subs, f"cpm_{region}"] + 1)
        ax.hlines(y + off, lo, hi, color=col, linewidth=2)
        ax.scatter(mean, y + off, color=col, s=22, zorder=3, label=f"{region} (3 donors, range)")
    ax.set_yticks(y, subs)
    ax.invert_yaxis()
    ax.set_xlabel("Adora1 log2(CPM+1)")
    ax.set_title("10x level by region\n(different donors and batches)")
    ax.legend(fontsize=6.5, frameon=False, loc="lower left")
    ax.grid(axis="x", alpha=0.25)
    # Right: MOp minus VISp in both datasets.
    ax = axes[1]
    ax.scatter(a.loc[subs, "diff_log2cpm"], y - 0.18, marker="D", s=26, color="#333333", label="10x difference")
    for k, (an, col) in enumerate(zip(animals, colors)):
        vals = [m.loc[MAP[s], an] if MAP[s] in m.index else np.nan for s in subs]
        ax.scatter(vals, y + 0.08 + 0.06 * k, s=12, color=col, label=an.replace("diff_", "MERFISH "))
    mean = [m.loc[MAP[s], "mean_difference"] if MAP[s] in m.index else np.nan for s in subs]
    ax.scatter(mean, y + 0.17, marker="|", s=160, color="black", label="MERFISH mean")
    ax.axvline(0, color="black", linewidth=0.8)
    ax.set_xlabel("MOp minus VISp (log2 units, scales differ)")
    ax.set_title("Regional difference\n10x vs MERFISH paired by animal")
    ax.legend(fontsize=6, frameon=False, loc="upper left", bbox_to_anchor=(1.01, 1))
    ax.grid(axis="x", alpha=0.25)
    fig.suptitle("Adora1 in MOp and VISp. Glutamatergic subclasses agree between datasets; the 10x MOp excess in "
                 "GABAergic subclasses is not seen in MERFISH", fontsize=8.5, y=1.02)
    save(fig, "fig2c_adora1_10x_vs_merfish")


def figS6(points, val):
    subs = [s for s in ORDER if s in set(val.subclass)]
    n = len(subs)
    ncol = 4
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 2.7 * nrow))
    v = val.set_index("subclass")
    for ax, s in zip(axes.flat, subs):
        p = points[points.subclass == s]
        ax.scatter(p.tenx_diff, p.merfish_diff, s=4, color="#9aa5b1", alpha=0.6, linewidths=0)
        c = p[p.merfish_called]
        ax.scatter(c.tenx_diff, c.merfish_diff, s=6, color="#1f5fa8", linewidths=0)
        ad = p[p.gene == "Adora1"]
        ax.scatter(ad.tenx_diff, ad.merfish_diff, s=30, color="#c2571a", edgecolors="black", linewidths=0.5,
                   zorder=4)
        ax.axhline(0, color="black", linewidth=0.5)
        ax.axvline(0, color="black", linewidth=0.5)
        ax.set_title(f"{s}\nr = {v.loc[s, 'pearson_r']:.2f}, sign agree {v.loc[s, 'sign_agreement_clear'] * 100:.0f}% "
                     f"(n = {v.loc[s, 'n_clear']})", fontsize=7)
        ax.tick_params(labelsize=6)
    for ax in list(axes.flat)[n:]:
        ax.axis("off")
    fig.supxlabel("10x MOp minus VISp, log2(CPM+1)  (region confounded with donor and batch)", fontsize=8)
    fig.supylabel("MERFISH MOp minus VISp, 4-animal paired mean", fontsize=8)
    fig.suptitle("Validation of 10x regional differences against MERFISH, genes in both datasets.\nBlue = called "
                 "regional DEG in MERFISH. Orange = Adora1. Sign agreement is over genes with MERFISH |diff| >= 0.10",
                 fontsize=8.5)
    fig.tight_layout()
    save(fig, "figS6_10x_vs_merfish_region_concordance")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sums, genes, groups = load()
    sets = pd.read_csv(WT / "focus_gene_sets.csv")
    mf = pd.read_csv(MF)
    cpm = donor_log2cpm(sums, genes, groups)
    val, points = validation(cpm, genes, mf)
    val.to_csv(OUT / "tenx_vs_merfish_region_concordance.csv", index=False)
    points.to_csv(OUT / "tenx_vs_merfish_region_points.csv.gz", index=False)
    table, order, gene_set = adenosine_table(cpm, genes, sets)
    table.to_csv(OUT / "adenosine_region_difference_10x.csv", index=False)
    fig2b(table, order, gene_set, val)
    fig2c(table, mf)
    figS6(points, val)
    print(val.round(2).to_string(index=False))
    t = table[~table.low_expression & table.donor_ranges_separate & table.diff_log2cpm.abs().ge(0.3)]
    print(t[["gene", "subclass", "diff_log2cpm", "cpm_MOp", "cpm_VISp"]].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
