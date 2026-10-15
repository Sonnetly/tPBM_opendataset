#!/usr/bin/env python3
"""Neuron-subclass marker DEG within each region, replicated across regions.

Question: which neuron subclasses express which genes more than the other neurons?

Design
------
The 10x MOp and VISp cells come from disjoint donors and library batches, so a region
contrast is confounded and is not attempted here. A subclass-versus-other-neurons contrast
inside one region uses cells from the same donors and batch, so it is not confounded.

For each region separately and each neuronal subclass S with >= MIN_CELLS cells in every
donor of that region:
- pseudobulk log2 CPM of S and of all other neurons pooled, per donor
  (CPM uses the full 31,053-gene library size)
- per-donor log2 fold change S over others (pseudocount 1 CPM), paired by donor
- one-sample t test of the per-donor log2FC across donors, BH within S
- detection fraction in S and in the other neurons
A gene is a marker of S in MOp when log2FC >= 1, every donor agrees in sign,
detection in S >= 0.25 and q < 0.05. It is "replicated" when, in the VISp donors (3 other
donors from other batches), log2FC >= 1 with all donors up and detection >= 0.25. The q
filter is not applied to replication because a 3-donor t test (df = 2) is underpowered.

Also writes a donor x subclass set-level table for the prespecified focus sets
(share of library counts from the set, as log2 ratio to the all-neuron share of that donor).
OXPHOS sets are reported for nuclear-encoded genes; mt-encoded transcripts are reported
separately because they also track mtDNA copy and cell quality.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
WT = ROOT / "results/rnaseq/whole_transcriptome"
OUT = WT / "subclass_markers"
MIN_CELLS = 20
NEURON_CLASSES = ("Glutamatergic", "GABAergic")
EXCLUDE_SUBCLASS = ("CR", "Meis2", "L5/6 IT TPE-ENT", "L2/3 IT PPP")


def bh(p):
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return q
    order = np.argsort(p[ok])
    r = p[ok][order] * ok.sum() / np.arange(1, ok.sum() + 1)
    r = np.minimum.accumulate(r[::-1])[::-1].clip(max=1)
    tmp = np.empty(ok.sum())
    tmp[order] = r
    q[ok] = tmp
    return q


def load():
    data = np.load(WT / "pseudobulk_counts.npz", allow_pickle=True)
    groups = pd.read_csv(WT / "pseudobulk_groups.csv").set_index("group").loc[data["groups"]]
    return data["sums"].astype(np.float64), data["detected"].astype(np.float64), pd.Index(data["genes"]), groups


def markers_for_region(region, sums, det, genes, groups):
    g = groups[(groups.region == region) & groups.cell_class.isin(NEURON_CLASSES)
               & ~groups.subclass.isin(EXCLUDE_SUBCLASS)]
    donors = sorted(g.donor.unique())
    eligible = [s for s, f in g.groupby("subclass") if len(f) == len(donors) and f.n_cells.min() >= MIN_CELLS]
    pos = {k: groups.index.get_loc(k) for k in g.index}
    rows = []
    for s in eligible:
        lfc = np.zeros((len(genes), len(donors)))
        det_in = np.zeros(len(genes))
        det_out = np.zeros(len(genes))
        n_in = n_out = 0
        cpm_in_all = np.zeros(len(genes))
        for j, d in enumerate(donors):
            gd = g[g.donor == d]
            ins = [pos[k] for k in gd.index[gd.subclass == s]]
            outs = [pos[k] for k in gd.index[gd.subclass != s]]
            lib_in = groups.library_counts.iloc[ins].sum()
            lib_out = groups.library_counts.iloc[outs].sum()
            cpm_in = sums[:, ins].sum(axis=1) / lib_in * 1e6
            cpm_out = sums[:, outs].sum(axis=1) / lib_out * 1e6
            lfc[:, j] = np.log2(cpm_in + 1) - np.log2(cpm_out + 1)
            cpm_in_all += cpm_in / len(donors)
            det_in += det[:, ins].sum(axis=1)
            det_out += det[:, outs].sum(axis=1)
            n_in += groups.n_cells.iloc[ins].sum()
            n_out += groups.n_cells.iloc[outs].sum()
        mean = lfc.mean(axis=1)
        sd = lfc.std(axis=1, ddof=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            t = mean / (sd / np.sqrt(len(donors)))
        p = 2 * stats.t.sf(np.abs(t), df=len(donors) - 1)
        p = np.where(np.isfinite(p), p, np.nan)
        frame = pd.DataFrame({"region": region, "subclass": s, "cell_class": g[g.subclass == s].cell_class.iloc[0],
                              "gene": genes, "log2FC": mean, "sd_log2FC": sd, "t": t, "p": p,
                              "donors_up": (lfc > 0).sum(axis=1), "donors_down": (lfc < 0).sum(axis=1),
                              "n_donors": len(donors), "pct_in": det_in / n_in, "pct_out": det_out / n_out,
                              "mean_cpm_in": cpm_in_all})
        frame = frame[(frame.pct_in >= 0.05) | (frame.pct_out >= 0.05)].copy()
        frame["q"] = bh(frame.p.to_numpy())
        frame["marker"] = ((frame.log2FC >= 1) & (frame.donors_up == len(donors)) & (frame.pct_in >= 0.25)
                           & (frame.q < 0.05))
        rows.append(frame)
    return pd.concat(rows, ignore_index=True), eligible, donors


def set_scores(sums, genes, groups, sets):
    neurons = groups[groups.cell_class.isin(NEURON_CLASSES) & ~groups.subclass.isin(EXCLUDE_SUBCLASS)]
    pos = {k: groups.index.get_loc(k) for k in neurons.index}
    present = sets[sets.present]
    defs = {}
    for name, f in present.groupby("set"):
        if name == "PBM_candidate62":
            continue
        members = f.gene.tolist()
        nuclear = [m for m in members if not m.startswith("mt-")]
        defs[name] = nuclear if name.startswith("OXPHOS") else members
        mt = [m for m in members if m.startswith("mt-")]
        if mt:
            defs[name + "_mt_encoded"] = mt
    idx = {k: genes.get_indexer(v) for k, v in defs.items()}
    rows = []
    for (region, donor), f in neurons.groupby(["region", "donor"]):
        cols = [pos[k] for k in f.index]
        lib_all = groups.library_counts.iloc[cols].sum()
        share_all = {k: sums[np.ix_(ix, cols)].sum() / lib_all for k, ix in idx.items()}
        for key, row in f.iterrows():
            c = groups.index.get_loc(key)
            for k, ix in idx.items():
                share = sums[ix, c].sum() / row.library_counts
                rows.append({"region": region, "donor": donor, "subclass": row.subclass,
                             "cell_class": row.cell_class, "n_cells": row.n_cells, "set": k,
                             "n_genes": len(ix), "share_of_library": share,
                             "log2_ratio_to_all_neurons": np.log2(share / share_all[k])})
    return pd.DataFrame(rows)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sums, det, genes, groups = load()
    sets = pd.read_csv(WT / "focus_gene_sets.csv")
    tables, meta = [], {}
    for region in ("MOp", "VISp"):
        table, eligible, donors = markers_for_region(region, sums, det, genes, groups)
        tables.append(table)
        meta[region] = {"donors": [str(d) for d in donors], "subclasses": eligible}
    allm = pd.concat(tables, ignore_index=True)
    allm.to_csv(OUT / "subclass_marker_stats_all.csv.gz", index=False)
    # Replication in the independent VISp donors uses effect size, sign agreement in all three
    # donors and detection, without the q filter: a 3-donor t test (df = 2) has too little power
    # to serve as the replication gate (e.g. VISp Car3 has 314 genes passing effect and sign but 3
    # passing q). Discovery in MOp keeps the q filter.
    visp = allm[allm.region == "VISp"]
    consistent = visp[(visp.log2FC >= 1) & (visp.donors_up == visp.n_donors) & (visp.pct_in >= 0.25)]
    rep = consistent[["subclass", "gene"]].copy()
    rep["replicated_MOp_and_VISp"] = True
    mop = allm[(allm.region == "MOp") & allm.marker].merge(rep, on=["subclass", "gene"], how="left")
    mop["replicated_MOp_and_VISp"] = mop.replicated_MOp_and_VISp.fillna(False)
    vis = allm[allm.region == "VISp"][["subclass", "gene", "log2FC", "q", "pct_in"]].rename(
        columns={"log2FC": "VISp_log2FC", "q": "VISp_q", "pct_in": "VISp_pct_in"})
    mop = mop.merge(vis, on=["subclass", "gene"], how="left").sort_values(["subclass", "log2FC"], ascending=[True, False])
    mop.to_csv(OUT / "MOp_subclass_markers.csv", index=False)
    focus = sets.groupby("gene").set.apply(lambda s: ";".join(sorted(set(s)))).rename("focus_sets")
    foc = allm.merge(focus, left_on="gene", right_index=True)
    foc.to_csv(OUT / "focus_gene_subclass_stats.csv", index=False)
    scores = set_scores(sums, genes, groups, sets)
    scores.to_csv(OUT / "focus_set_scores_by_donor.csv", index=False)
    summary = mop.groupby("subclass").agg(markers_MOp=("gene", "size"),
                                          replicated_in_VISp=("replicated_MOp_and_VISp", "sum")).reset_index()
    summary.to_csv(OUT / "marker_counts.csv", index=False)
    report = {"design": "within-region subclass vs other neurons, paired by donor; replication across disjoint "
                        "MOp and VISp donors and batches",
              "rule": {"log2FC_min": 1, "all_donors_same_sign": True, "pct_in_min": 0.25, "q_max_MOp_discovery": 0.05, "VISp_replication": "log2FC>=1, all 3 donors up, pct_in>=0.25, no q filter",
                       "min_cells_per_donor_subclass": MIN_CELLS},
              "regions": meta, "markers_MOp": int(len(mop)),
              "replicated_in_VISp": int(mop.replicated_MOp_and_VISp.sum()),
              "caution": "Snap25-IRES2-Cre;Ai14 RFP+ sorted neurons, male P57-P59. Not WT, no PBM condition. "
                         "Region contrasts are not computed because region equals donor and batch."}
    (OUT / "subclass_marker_report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))
    print(summary.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
