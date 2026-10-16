#!/usr/bin/env python3
"""Energy metabolism set scores in MOp and VISp, gene-level detail and leave-one-out check.

Purpose
-------
Fig3 summarises each energy set as one score per neuron subclass. A summed score can be
driven by one abundant gene, so this script also writes, for every set, the per-gene values
behind the score, each gene's share of the set's counts, and a leave-one-out (LOO) check
that recomputes the score with each gene removed in turn.

Definitions (10x GSE185862, neurons only, per region)
-----------------------------------------------------
- Set score of subclass S in donor d
      log2( share of S library from set genes / same share in all neurons of donor d )
  then averaged over the 3 donors of that region. Each region is scaled to its own neurons,
  so MOp and VISp panels are independent replications of the same within-region pattern;
  cells must not be compared across panels (region equals donor and batch in this dataset).
- Gene value: log2( (CPM in S + 1) / (CPM in all neurons + 1) ) per donor, donor mean.
- Contribution: share of a gene in the set's summed counts across all neurons of the region.
- OXPHOS rows use nuclear-encoded subunits. The 13 mtDNA-encoded subunits are one separate
  row for all complexes, because their level also tracks mtDNA copy number and cell quality.
- A subclass is shown in a region only when every donor of that region has >= MIN_CELLS cells.

Outputs: results/rnaseq/whole_transcriptome/energy_detail/ and figures/report_m1v1/
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
OUT = WT / "energy_detail"
FIG = ROOT / "figures/report_m1v1"
MIN_CELLS = 20
NEURON_CLASSES = ("Glutamatergic", "GABAergic")
EXCLUDE_SUBCLASS = ("CR", "Meis2", "L5/6 IT TPE-ENT", "L2/3 IT PPP")
ORDER = ["L2/3 IT CTX", "L4/5 IT CTX", "L5 IT CTX", "L6 IT CTX", "L5 PT CTX", "L5/6 NP CTX", "L6 CT CTX",
         "L6b CTX", "Car3", "Pvalb", "Sst", "Sst Chodl", "Vip", "Lamp5", "Sncg"]
ROWS = [  # (row key, readable label, source set, gene filter)
    ("CI", "Complex I, NADH dehydrogenase", "OXPHOS_CI", "nuclear"),
    ("CII", "Complex II, succinate dehydrogenase", "OXPHOS_CII", "nuclear"),
    ("CIII", "Complex III, cytochrome bc1", "OXPHOS_CIII", "nuclear"),
    ("CIV", "Complex IV, cytochrome c oxidase", "OXPHOS_CIV", "nuclear"),
    ("CV", "Complex V, ATP synthase", "OXPHOS_CV", "nuclear"),
    ("mtOX", "mtDNA-encoded OXPHOS (13 genes)", "OXPHOS_*", "mt"),
    ("CYCS", "Cytochrome c (Cycs)", "cytochrome_c", "all"),
    ("BIO", "Mito biogenesis (Ppargc1a, Tfam, ...)", "mito_biogenesis", "all"),
    ("CK", "Creatine / phosphocreatine", "creatine_buffer", "all"),
    ("GLY", "Glycolysis and lactate", "glycolysis_lactate", "all"),
]
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "font.family": "DejaVu Sans"})


def row_genes(sets):
    present = sets[sets.present]
    out = {}
    for key, _, source, filt in ROWS:
        if source == "OXPHOS_*":
            genes = present[present.set.str.startswith("OXPHOS_")].gene
        else:
            genes = present[present.set == source].gene
        genes = sorted(set(genes))
        if filt == "nuclear":
            genes = [g for g in genes if not g.startswith("mt-")]
        elif filt == "mt":
            genes = [g for g in genes if g.startswith("mt-")]
        out[key] = genes
    return out


def load():
    data = np.load(WT / "pseudobulk_counts.npz", allow_pickle=True)
    groups = pd.read_csv(WT / "pseudobulk_groups.csv").set_index("group").loc[data["groups"]]
    return data["sums"].astype(np.float64), pd.Index(data["genes"]), groups


def eligible(groups, region):
    g = groups[(groups.region == region) & groups.cell_class.isin(NEURON_CLASSES)
               & ~groups.subclass.isin(EXCLUDE_SUBCLASS)]
    donors = g.donor.nunique()
    ok = [s for s, f in g.groupby("subclass") if f.donor.nunique() == donors and f.n_cells.min() >= MIN_CELLS]
    return g, ok


def compute(sums, genes, groups, rows):
    set_rows, gene_rows, loo_rows, contrib_rows = [], [], [], []
    for region in ("MOp", "VISp"):
        g, ok = eligible(groups, region)
        col = {k: groups.index.get_loc(k) for k in g.index}
        for key, gl in rows.items():
            ix = genes.get_indexer(gl)
            pool = sums[np.ix_(ix, [col[k] for k in g.index])].sum(axis=1)
            for gene, share in zip(gl, pool / pool.sum()):
                contrib_rows.append({"region": region, "row": key, "gene": gene, "share_of_set_counts": share})
        for donor, f in g.groupby("donor"):
            cols = [col[k] for k in f.index]
            lib_all = f.library_counts.sum()
            all_counts = sums[:, cols].sum(axis=1)
            for key_g, r in f.iterrows():
                if r.subclass not in ok:
                    continue
                c = col[key_g]
                for key, gl in rows.items():
                    ix = genes.get_indexer(gl)
                    sub = sums[ix, c]
                    ref = all_counts[ix]
                    score = np.log2((sub.sum() / r.library_counts) / (ref.sum() / lib_all))
                    set_rows.append({"region": region, "donor": donor, "subclass": r.subclass, "row": key,
                                     "score": score})
                    for j, gene in enumerate(gl):
                        cpm_s = sub[j] / r.library_counts * 1e6
                        cpm_a = ref[j] / lib_all * 1e6
                        gene_rows.append({"region": region, "donor": donor, "subclass": r.subclass, "row": key,
                                          "gene": gene, "log2_ratio": np.log2((cpm_s + 1) / (cpm_a + 1)),
                                          "cpm_all_neurons": cpm_a})
                        if len(gl) > 1:
                            keep = np.arange(len(gl)) != j
                            loo = np.log2((sub[keep].sum() / r.library_counts) / (ref[keep].sum() / lib_all))
                            loo_rows.append({"region": region, "donor": donor, "subclass": r.subclass, "row": key,
                                             "dropped_gene": gene, "loo_score": loo})
    sets_df = pd.DataFrame(set_rows)
    gene_df = pd.DataFrame(gene_rows)
    loo_df = pd.DataFrame(loo_rows)
    contrib = pd.DataFrame(contrib_rows)
    return sets_df, gene_df, loo_df, contrib


def summarise(sets_df, loo_df, contrib):
    mean = sets_df.groupby(["region", "subclass", "row"]).score.agg(
        score="mean", donors_pos=lambda v: int((v > 0).sum()), donors_neg=lambda v: int((v < 0).sum()),
        n_donors="size").reset_index()
    loo_mean = loo_df.groupby(["region", "subclass", "row", "dropped_gene"]).loo_score.mean().reset_index()
    loo_sum = loo_mean.groupby(["region", "subclass", "row"]).loo_score.agg(loo_min="min", loo_max="max").reset_index()
    out = mean.merge(loo_sum, on=["region", "subclass", "row"], how="left")
    out["loo_max_shift"] = np.maximum((out.loo_max - out.score).abs(), (out.loo_min - out.score).abs())
    out["loo_sign_stable"] = np.where(out.loo_min.isna(), np.nan,
                                      ((out.loo_min > 0) & (out.score > 0)) | ((out.loo_max < 0) & (out.score < 0))
                                      | (out.score.abs() < 0.1))
    top = contrib.sort_values("share_of_set_counts", ascending=False).groupby(["region", "row"]).head(1)
    top = top.rename(columns={"gene": "top_contributor", "share_of_set_counts": "top_contributor_share"})
    out = out.merge(top, on=["region", "row"], how="left")
    return out


def subclass_columns(groups):
    _, mop = eligible(groups, "MOp")
    _, vis = eligible(groups, "VISp")
    return [s for s in ORDER if s in set(mop) | set(vis)]


def fig3(summary, cols):
    labels = {k: lab for k, lab, _, _ in ROWS}
    keys = [k for k, *_ in ROWS]
    vals = {}
    for region in ("MOp", "VISp"):
        s = summary[summary.region == region]
        vals[region] = (s.pivot_table(index="row", columns="subclass", values="score").reindex(index=keys, columns=cols),
                        s.assign(agree=(s.donors_pos == s.n_donors) | (s.donors_neg == s.n_donors))
                        .pivot_table(index="row", columns="subclass", values="agree").reindex(index=keys, columns=cols))
    lim = max(float(np.nanmax(np.abs(v[0].to_numpy()))) for v in vals.values())
    norm = TwoSlopeNorm(0, -lim, lim)
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 4.6), sharey=True, gridspec_kw={"wspace": 0.04})
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#e6e6e6")
    for ax, region in zip(axes, ("MOp", "VISp")):
        mean, agree = vals[region]
        ax.imshow(np.ma.masked_invalid(mean.to_numpy()), aspect="auto", cmap=cmap, norm=norm)
        for i in range(len(keys)):
            for j in range(len(cols)):
                v = mean.iat[i, j]
                if np.isnan(v):
                    ax.text(j, i, "n/a", ha="center", va="center", fontsize=5, color="#777777")
                    continue
                mark = "" if agree.iat[i, j] == 1 else "*"
                ax.text(j, i, f"{v:+.2f}{mark}", ha="center", va="center", fontsize=5.4,
                        color="white" if abs(v) > 0.6 * lim else "#222222")
        ax.axhline(5.5, color="white", linewidth=2)
        ax.set_xticks(range(len(cols)), cols, rotation=60, ha="right")
        ax.set_title(f"{region} (3 donors)")
    axes[0].set_yticks(range(len(keys)), [labels[k] for k in keys])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    fig.colorbar(sm, ax=axes, shrink=0.7, pad=0.01, label="log2 ratio to all neurons of the same region")
    fig.suptitle("Energy metabolism gene sets by neuron subclass (10x). Each panel is scaled to its own region's "
                 "neurons; compare patterns, not cells across panels.\n* donors disagree in sign. n/a = fewer than "
                 f"{MIN_CELLS} cells in a donor. Gene-level detail in Fig S1 to S5", fontsize=8, y=1.03)
    save(fig, "fig3_10x_energy_module_heatmap")


def save(fig, name):
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{name}.png", dpi=180, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


SUPP = [("figS1_energy_genes_complexI", ["CI", "mtOX"]),
        ("figS2_energy_genes_complexII_III", ["CII", "CIII"]),
        ("figS3_energy_genes_complexIV_cycs", ["CIV", "CYCS"]),
        ("figS4_energy_genes_complexV", ["CV"]),
        ("figS5_energy_genes_biogenesis_creatine_glycolysis", ["BIO", "CK", "GLY"])]
CLIP = 2.0


SHORT = {"CI": "Complex I", "CII": "Complex II", "CIII": "Complex III", "CIV": "Complex IV",
         "CV": "Complex V", "mtOX": "mtDNA OXPHOS", "CYCS": "Cytochrome c", "BIO": "Mito biogenesis",
         "CK": "Creatine", "GLY": "Glycolysis/lactate"}


def supp_figures(gene_df, contrib, summary, cols):
    """One page per group of sets. Layout, left to right:
    gene names | MOp heatmap | VISp heatmap | bar of each gene's share of the set counts.
    Gene names are the y tick labels of the MOp panel, so they always sit next to the cells."""
    labels = {k: lab for k, lab, _, _ in ROWS}
    gmean = gene_df.groupby(["region", "subclass", "row", "gene"]).agg(
        log2_ratio=("log2_ratio", "mean"), cpm_all=("cpm_all_neurons", "mean")).reset_index()
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#e6e6e6")
    norm = TwoSlopeNorm(0, -CLIP, CLIP)
    for name, keys in SUPP:
        blocks = []
        for key in keys:
            c = contrib[(contrib.region == "MOp") & (contrib.row == key)].sort_values("share_of_set_counts",
                                                                                       ascending=False)
            blocks.append((key, c.gene.tolist(), c.share_of_set_counts.to_numpy()))
        n_rows = sum(len(b[1]) + 1 for b in blocks)
        row_h = 0.24
        height = row_h * (n_rows + 1.2 * len(blocks)) + 2.6
        fig = plt.figure(figsize=(17.5, height))
        top = 1 - 1.05 / height
        bottom = 1.45 / height
        fig.subplots_adjust(left=0.13, right=0.90, top=top, bottom=bottom, wspace=0.05)
        grid = fig.add_gridspec(1, 3, width_ratios=[len(cols), len(cols), 3.2], left=0.13, right=0.90,
                                top=top, bottom=bottom, wspace=0.05)
        ax_m = fig.add_subplot(grid[0])
        ax_v = fig.add_subplot(grid[1], sharey=ax_m)
        ax_bar = fig.add_subplot(grid[2], sharey=ax_m)
        y = 0.0
        ticks, ticklabels, kinds = [], [], []
        for b, (key, gl, share) in enumerate(blocks):
            for ax, region in ((ax_m, "MOp"), (ax_v, "VISp")):
                srow = summary[(summary.region == region) & (summary.row == key)].set_index("subclass").reindex(cols)
                vals = srow.score.to_numpy()
                ax.imshow(np.ma.masked_invalid(vals[None, :]), extent=(-0.5, len(cols) - 0.5, y + 0.5, y - 0.5),
                          aspect="auto", cmap=cmap, norm=norm)
                for j, v in enumerate(vals):
                    if np.isfinite(v):
                        ax.text(j, y, f"{v:+.2f}", ha="center", va="center", fontsize=6.3, fontweight="bold",
                                color="white" if abs(v) > 0.6 * CLIP else "#222222")
                    else:
                        ax.text(j, y, "n/a", ha="center", va="center", fontsize=5.5, color="#888888")
                ax.add_patch(plt.Rectangle((-0.5, y - 0.5), len(cols), 1, fill=False, edgecolor="black",
                                           linewidth=1.0))
            ticks.append(y)
            ticklabels.append(f"{SHORT[key]} SET SCORE")
            kinds.append("set")
            y += 1.25
            for gene, sh in zip(gl, share):
                for ax, region in ((ax_m, "MOp"), (ax_v, "VISp")):
                    row = gmean[(gmean.region == region) & (gmean.row == key) & (gmean.gene == gene)]
                    vals = row.set_index("subclass").reindex(cols).log2_ratio.to_numpy()
                    ax.imshow(np.ma.masked_invalid(np.clip(vals, -CLIP, CLIP)[None, :]),
                              extent=(-0.5, len(cols) - 0.5, y + 0.5, y - 0.5), aspect="auto", cmap=cmap, norm=norm)
                ax_bar.barh(y, sh * 100, height=0.72, color="#5b6b7a")
                ax_bar.text(sh * 100 + 2, y, f"{sh * 100:.1f}%", va="center", fontsize=6, color="#333333")
                cpm = gmean[(gmean.region == "MOp") & (gmean.row == key) & (gmean.gene == gene)].cpm_all.mean()
                ticks.append(y)
                ticklabels.append(gene)
                kinds.append("low" if cpm < 10 else "gene")
                y += 1
            y += 1.0
        for ax in (ax_m, ax_v):
            ax.set_xlim(-0.5, len(cols) - 0.5)
            ax.set_xticks(range(len(cols)), cols, rotation=60, ha="right", fontsize=7.5)
        ax_m.set_ylim(y - 1.5, -0.55)
        ax_m.set_yticks(ticks, ticklabels, fontsize=7.5)
        for lab, kind in zip(ax_m.get_yticklabels(), kinds):
            if kind == "set":
                lab.set_fontweight("bold")
            else:
                lab.set_fontstyle("italic")
                if kind == "low":
                    lab.set_color("#9a9a9a")
        ax_v.tick_params(axis="y", left=False, labelleft=False)
        ax_bar.tick_params(axis="y", left=False, labelleft=False)
        ax_m.set_title("MOp (3 donors)", fontsize=9)
        ax_v.set_title("VISp (3 other donors)", fontsize=9)
        ax_bar.set_title("share of set counts\n(all MOp neurons)", fontsize=7.5)
        ax_bar.set_xlim(0, 125)
        ax_bar.set_xticks([0, 50, 100], ["0", "50", "100%"], fontsize=7)
        ax_bar.spines[["right", "top", "left"]].set_visible(False)
        cax = fig.add_axes([0.925, bottom + (top - bottom) * 0.25, 0.010, (top - bottom) * 0.5])
        fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), cax=cax,
                     label=f"log2 ratio to all neurons\nof the same region (genes clipped at +/-{CLIP:g})")
        fig.text(0.5, 1 - 0.22 / height, f"{' / '.join(labels[k] for k in keys)}", ha="center", va="top",
                 fontsize=10, fontweight="bold")
        fig.text(0.5, 1 - 0.47 / height, "Boxed bold row = set score. Rows below = genes in that set, ordered by "
                 "their share of the set counts. Gray name = under 10 CPM in all MOp neurons. "
                 "Each panel is scaled to its own region.", ha="center", va="top", fontsize=7.5, color="#333333")
        save(fig, name)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sums, genes, groups = load()
    sets = pd.read_csv(WT / "focus_gene_sets.csv")
    rows = row_genes(sets)
    pd.DataFrame([{"row": k, "label": lab, "genes": ";".join(rows[k]), "n_genes": len(rows[k])}
                  for k, lab, _, _ in ROWS]).to_csv(OUT / "energy_row_definitions.csv", index=False)
    sets_df, gene_df, loo_df, contrib = compute(sums, genes, groups, rows)
    summary = summarise(sets_df, loo_df, contrib)
    sets_df.to_csv(OUT / "energy_set_scores_by_donor.csv", index=False)
    gene_df.to_csv(OUT / "energy_gene_values_by_donor.csv.gz", index=False)
    loo_df.to_csv(OUT / "energy_leave_one_out_by_donor.csv.gz", index=False)
    contrib.to_csv(OUT / "energy_gene_contribution.csv", index=False)
    summary.to_csv(OUT / "energy_set_summary.csv", index=False)
    cols = subclass_columns(groups)
    fig3(summary, cols)
    supp_figures(gene_df, contrib, summary, cols)
    s = summary[summary.n_donors > 0]
    focus = s[s.subclass.isin(["Pvalb", "L5 PT CTX"])][["region", "subclass", "row", "score", "loo_min", "loo_max",
                                                         "top_contributor", "top_contributor_share"]]
    print(focus.round(3).to_string(index=False))
    print("max LOO shift overall", round(float(s.loo_max_shift.max()), 3),
          "| rows with LOO sign change:", int((s.loo_sign_stable == False).sum()))  # noqa: E712
    print("MOp-VISp score r:", round(s.pivot_table(index=["subclass", "row"], columns="region", values="score")
                                    .dropna().corr().iloc[0, 1], 3))


if __name__ == "__main__":
    main()
