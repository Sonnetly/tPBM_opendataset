#!/usr/bin/env python3
"""PBM same-direction DEG across neuron subclasses in MOp and VISp, and their MOp minus VISp difference.

Fig4a  MOp and VISp side by side. Row z-score is computed over BOTH regions' columns together,
       so the same color means the same log2(CPM+1) level in either panel.
Fig4b  MOp minus VISp of donor-mean log2(CPM+1) per subclass (10x). Each subclass is labelled with
       the sign agreement of 10x regional differences with MERFISH four-animal paired differences
       (from adenosine_region_compare.py); a dagger marks cells where all three MOp donors lie beyond
       all three VISp donors; genes under MIN_CPM in both regions are masked.

Rows keep the order of the original Fig4 (PBM direction, then the MOp subclass of peak expression).
The 10x region contrast is confounded with donor and batch; read Fig4b only where the MERFISH
agreement is high. No PBM condition is present in the atlas.
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
OUT = WT / "pbm_candidates_region"
FIG = ROOT / "figures/report_m1v1"
MIN_CELLS = 20
MIN_CPM = 5.0
ORDER = ["L2/3 IT CTX", "L4/5 IT CTX", "L5 IT CTX", "L6 IT CTX", "L5 PT CTX", "L5/6 NP CTX", "L6 CT CTX",
         "L6b CTX", "Car3", "Pvalb", "Sst", "Sst Chodl", "Vip", "Lamp5", "Sncg"]
plt.rcParams.update({"font.size": 8, "axes.titlesize": 9, "font.family": "DejaVu Sans"})


def save(fig, name):
    FIG.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIG / f"{name}.png", dpi=180, bbox_inches="tight")
    fig.savefig(FIG / f"{name}.pdf", bbox_inches="tight")
    plt.close(fig)


def donor_cpm(sums, genes, groups, gene_list):
    """dict (region, subclass) -> donors x genes CPM, eligible groups only."""
    idx = genes.get_indexer(gene_list)
    out = {}
    for (region, sub), f in groups.groupby(["region", "subclass"]):
        if sub not in ORDER or f.donor.nunique() < 3 or f.n_cells.min() < MIN_CELLS:
            continue
        out[(region, sub)] = np.stack([sums[idx, groups.index.get_loc(k)] / r.library_counts * 1e6
                                       for k, r in f.iterrows()])
    return out


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    data = np.load(WT / "pseudobulk_counts.npz", allow_pickle=True)
    genes = pd.Index(data["genes"])
    sums = data["sums"].astype(np.float64)
    groups = pd.read_csv(WT / "pseudobulk_groups.csv").set_index("group").loc[data["groups"]]
    rows = pd.read_csv(WT / "subclass_markers/pbm_candidates_peak_subclass_MOp.csv")
    gl = rows.gene.tolist()
    labels = [f"{g} ({d})" for g, d in zip(rows.gene, rows.direction)]
    cpm = donor_cpm(sums, genes, groups, gl)
    subs = [s for s in ORDER if ("MOp", s) in cpm or ("VISp", s) in cpm]

    def mean_log(region):
        return pd.DataFrame({s: np.log2(cpm[(region, s)] + 1).mean(axis=0) if (region, s) in cpm
                             else np.full(len(gl), np.nan) for s in subs}, index=gl)

    mop, vis = mean_log("MOp"), mean_log("VISp")
    both = pd.concat([mop, vis], axis=1)
    mu = both.mean(axis=1)
    sd = both.std(axis=1).replace(0, np.nan)
    zm, zv = mop.sub(mu, axis=0).div(sd, axis=0), vis.sub(mu, axis=0).div(sd, axis=0)

    # Long table with donor ranges for the difference figure.
    recs = []
    for i, g in enumerate(gl):
        for s in subs:
            if ("MOp", s) not in cpm or ("VISp", s) not in cpm:
                continue
            m, v = cpm[("MOp", s)][:, i], cpm[("VISp", s)][:, i]
            lm, lv = np.log2(m + 1), np.log2(v + 1)
            recs.append({"gene": g, "direction": rows.direction.iloc[i], "subclass": s,
                         "cpm_MOp": m.mean(), "cpm_VISp": v.mean(), "diff_log2cpm": lm.mean() - lv.mean(),
                         "donor_ranges_separate": bool(lm.min() > lv.max() or lm.max() < lv.min()),
                         "low_expression": bool(max(m.mean(), v.mean()) < MIN_CPM)})
    t = pd.DataFrame(recs)
    t.to_csv(OUT / "pbm_candidates_region_difference_10x.csv", index=False)
    pd.concat({"MOp": mop, "VISp": vis}, axis=1).to_csv(OUT / "pbm_candidates_log2cpm_by_region.csv")

    # Fig4a, shared z-score scale.
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#e6e6e6")
    h = 0.17 * len(gl) + 2.0
    fig, axes = plt.subplots(1, 2, figsize=(16, h), sharey=True, gridspec_kw={"wspace": 0.04})
    fig.subplots_adjust(top=1 - 0.75 / h)
    for ax, z, region in ((axes[0], zm, "MOp (3 donors)"), (axes[1], zv, "VISp (3 other donors)")):
        ax.imshow(np.ma.masked_invalid(z.to_numpy()), aspect="auto", cmap=cmap, vmin=-2.5, vmax=2.5)
        ax.set_xticks(range(len(subs)), subs, rotation=60, ha="right")
        ax.set_title(region)
        n_down = int((rows.direction == "down").sum())
        ax.axhline(n_down - 0.5, color="black", linewidth=1)
    axes[0].set_yticks(range(len(gl)), labels, fontsize=6.3)
    fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(-2.5, 2.5)), ax=axes, shrink=0.35,
                 pad=0.01, label="row z-score over both regions")
    fig.suptitle("PBM same-direction DEG (hippocampus and cortex) across neuron subclasses, MOp and VISp (10x).\n"
                 "Row z-score uses both regions together, so equal colors mean equal log2(CPM+1). Gray = subclass "
                 f"under {MIN_CELLS} cells in a donor. Baseline atlas, not PBM response", fontsize=8.5,
                 y=1 - 0.12 / h)
    save(fig, "fig4a_10x_pbm_candidates_MOp_VISp")

    # Fig4b, MOp minus VISp.
    val = pd.read_csv(WT / "adenosine_region/tenx_vs_merfish_region_concordance.csv").set_index("subclass")
    dsub = [s for s in subs if ("MOp", s) in cpm and ("VISp", s) in cpm]
    diff = t.pivot_table(index="gene", columns="subclass", values="diff_log2cpm").reindex(index=gl, columns=dsub)
    low = t.pivot_table(index="gene", columns="subclass", values="low_expression").reindex(index=gl, columns=dsub)
    sep = t.pivot_table(index="gene", columns="subclass", values="donor_ranges_separate").reindex(index=gl,
                                                                                               columns=dsub)
    lim = 1.5
    fig, ax = plt.subplots(figsize=(0.62 * len(dsub) + 3.4, 0.2 * len(gl) + 2.6))
    ax.imshow(np.ma.masked_invalid(diff.where(low != 1).to_numpy().clip(-lim, lim)), aspect="auto", cmap=cmap,
              norm=TwoSlopeNorm(0, -lim, lim))
    for i in range(len(gl)):
        for j in range(len(dsub)):
            if low.iat[i, j] == 1:
                ax.text(j, i, "low", ha="center", va="center", fontsize=4.6, color="#9a9a9a")
                continue
            v = diff.iat[i, j]
            ax.text(j, i, f"{v:+.2f}{chr(8224) if sep.iat[i, j] == 1 else ''}", ha="center", va="center",
                    fontsize=4.8, color="white" if abs(v) > 0.9 else "#222222")
    n_down = int((rows.direction == "down").sum())
    ax.axhline(n_down - 0.5, color="black", linewidth=1)
    if "Pvalb" in dsub:
        ax.axvline(dsub.index("Pvalb") - 0.5, color="black", linewidth=0.8)
    xl = []
    for s in dsub:
        a = val.loc[s, "sign_agreement_clear"] if s in val.index else np.nan
        xl.append(f"{s}\n{a * 100:.0f}%" if np.isfinite(a) else f"{s}\nn/a")
    ax.set_xticks(range(len(dsub)), xl, rotation=60, ha="right")
    for tick, s in zip(ax.get_xticklabels(), dsub):
        a = val.loc[s, "sign_agreement_clear"] if s in val.index else np.nan
        if not np.isfinite(a) or a < 0.85:
            tick.set_color("#b05a00")
    ax.set_yticks(range(len(gl)), labels, fontsize=6.3)
    fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=TwoSlopeNorm(0, -lim, lim)), ax=ax, shrink=0.4, pad=0.02,
                 label="MOp minus VISp, log2(CPM+1)\nred = higher in MOp")
    ax.set_title("PBM same-direction DEG, MOp minus VISp within each neuron subclass (10x)\n"
                 "% under subclass = sign agreement of 10x regional differences with MERFISH 4-animal paired "
                 "differences.\nOrange = under 85%, read with caution. " + chr(8224) + " all 3 MOp donors beyond "
                 f"all 3 VISp donors. low = under {MIN_CPM:g} CPM in both regions", fontsize=7.5)
    save(fig, "fig4b_10x_pbm_candidates_region_difference")

    glut = t[t.subclass.isin(["L2/3 IT CTX", "L4/5 IT CTX", "L5 IT CTX", "L6 IT CTX", "L5 PT CTX",
                              "L5/6 NP CTX", "L6 CT CTX"]) & ~t.low_expression]
    strong = glut[glut.donor_ranges_separate & glut.diff_log2cpm.abs().ge(0.5)]
    print(strong.sort_values("diff_log2cpm")[["gene", "direction", "subclass", "diff_log2cpm", "cpm_MOp",
                                              "cpm_VISp"]].round(2).to_string(index=False))


if __name__ == "__main__":
    main()
