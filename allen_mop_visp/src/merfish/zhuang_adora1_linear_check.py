#!/usr/bin/env python3
"""Adora1 MOp versus VISp per animal and subclass on the linear scale, from Zhuang MERFISH cells.

Why
---
run_zhuang.py reports the difference of mean official log2(x+1) values. For a gene near
the detection floor (Adora1 in GABAergic neurons, mean log2 about 0.1 to 0.3) a small
log2(x+1) difference can still be a large fold difference, so that scale can understate
regional differences in low-expressing cell types. This script recomputes, per animal,
region and subclass, from every cell:
- mean linear expression  x = 2**v - 1  (v = official log2 value)
- detection fraction       v > 0
- MOp / VISp ratio of linear means and its log2
Only the Adora1 column is read (h5py, column chunks), so memory stays small.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_zhuang as rz  # noqa: E402

ANIMALS = ["Zhuang-ABCA-1", "Zhuang-ABCA-2", "Zhuang-ABCA-3", "Zhuang-ABCA-4"]
GENES = ["Adora1"]
OUT = rz.RESULTS / "adora1_linear"
MIN_CELLS = rz.MIN_CELLS_PER_STRATUM


def gene_columns(path, genes):
    with h5py.File(path, "r") as f:
        var = f["var"]
        key = "gene_symbol" if "gene_symbol" in var else var.attrs["_index"]
        col = var[key]
        if isinstance(col, h5py.Group):
            symbols = col["categories"][:][col["codes"][:]]
        else:
            symbols = col[:]
        symbols = np.array([s.decode() if isinstance(s, bytes) else str(s) for s in symbols])
        obs = f["obs"]
        okey = obs.attrs["_index"]
        cells = np.array([s.decode() if isinstance(s, bytes) else str(s) for s in obs[okey][:]])
        values = {g: f["X"][:, int(np.flatnonzero(symbols == g)[0])] for g in genes}
    return cells, values


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    manifest, _, _ = rz.fetch_manifest(rz.RELEASE, live=True)
    rows = []
    for animal in ANIMALS:
        cells = rz.select_cells(manifest, animal)
        spec = rz.resolve(manifest, animal, "expression_matrices", animal, kind="log2")
        path = rz.download_verified(spec)
        labels, values = gene_columns(path, GENES)
        index = pd.Series(np.arange(len(labels)), index=labels)
        meta = cells.set_index("cell_label")[["region", "subclass_clean", "cell_class"]]
        meta = meta[meta.index.isin(index.index)]
        rowidx = index.loc[meta.index].to_numpy()
        for gene, vec in values.items():
            v = vec[rowidx].astype(np.float64)
            frame = meta.assign(v=v, x=np.exp2(v) - 1, det=(v > 0))
            for (sub, cls, reg), f in frame.groupby(["subclass_clean", "cell_class", "region"]):
                rows.append({"animal": animal, "gene": gene, "subclass": sub, "cell_class": cls, "region": reg,
                             "n_cells": len(f), "mean_log2": f.v.mean(), "mean_linear": f.x.mean(),
                             "detection": f.det.mean()})
        print(animal, "done", flush=True)
    long = pd.DataFrame(rows)
    long.to_csv(OUT / "adora1_by_animal_region_subclass.csv", index=False)
    w = long.pivot_table(index=["gene", "animal", "subclass", "cell_class"], columns="region",
                         values=["n_cells", "mean_log2", "mean_linear", "detection"])
    w.columns = [f"{a}_{b}" for a, b in w.columns]
    w = w.reset_index()
    w = w[(w.n_cells_MOp >= MIN_CELLS) & (w.n_cells_VISp >= MIN_CELLS)].copy()
    w["log2_diff_of_means"] = w.mean_log2_MOp - w.mean_log2_VISp
    w["linear_ratio_MOp_over_VISp"] = w.mean_linear_MOp / w.mean_linear_VISp
    w["log2_linear_ratio"] = np.log2(w.linear_ratio_MOp_over_VISp)
    w["detection_diff"] = w.detection_MOp - w.detection_VISp
    w.to_csv(OUT / "adora1_region_comparison_per_animal.csv", index=False)
    summary = (w[w.cell_class != "Other"].groupby(["gene", "subclass", "cell_class"])
               .agg(n_animals=("animal", "size"), log2_linear_ratio_mean=("log2_linear_ratio", "mean"),
                    log2_linear_ratio_min=("log2_linear_ratio", "min"),
                    log2_linear_ratio_max=("log2_linear_ratio", "max"),
                    animals_MOp_higher=("log2_linear_ratio", lambda s: int((s > 0).sum())),
                    detection_MOp=("detection_MOp", "mean"), detection_VISp=("detection_VISp", "mean"),
                    log2_diff_of_means=("log2_diff_of_means", "mean"))
               .reset_index().sort_values(["cell_class", "subclass"]))
    summary.to_csv(OUT / "adora1_region_summary.csv", index=False)
    (OUT / "adora1_linear_report.json").write_text(json.dumps({
        "purpose": "Linear-scale MOp/VISp ratio of Adora1 per animal; checks whether log2(x+1) differences "
                   "understate regional differences at low expression",
        "min_cells_per_region": MIN_CELLS, "animals": ANIMALS}, indent=2))
    print(summary.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
