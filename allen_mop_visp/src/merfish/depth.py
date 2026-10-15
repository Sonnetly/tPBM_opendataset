#!/usr/bin/env python3
"""Physical cortical depth from the pial surface for B6J MERFISH MOp/VISp cells.

Why this exists
---------------
`depth_group` in run_merfish.py is derived from the CCF anatomical layer label
(L1/L2/3 -> superficial, L4 -> middle, L5/L6 -> deep). Layer labels are NOT
depth-matched between regions: MOp cortex is roughly twice as thick as VISp, so
"MOp superficial" and "VISp superficial" occupy different physical depth ranges.
Photobiomodulation dose falls off with physical depth, not with layer identity,
so any PBM-relevant composition comparison has to be cut in micrometres.

Method
------
Allen CCF is ASL with x = anterior-posterior, y = dorsal-ventral, z = left-right,
in millimetres. Within one coronal section the cortical surface normal lies in the
(y, z) plane, so depth is computed in that plane only.

For each (brain_section_label, region, hemisphere) group the pial reference is the
set of CCF-layer-1 cells in that same group. Each cell's depth is its Euclidean
distance, in micrometres, to the k-th nearest reference cell. k=1 is the primary
definition; a larger k is reported as a robustness column because a single
mislabelled L1 cell can drag a nearest-neighbour distance.

Known limitations, carried into the output as columns and into the audit JSON
------------------------------------------------------------------------------
1. L1 cell somata sit a few tens of micrometres below the true pia, so absolute
   depth is systematically underestimated. The bias is shared by both regions and
   largely cancels in a MOp-minus-VISp comparison. These are NOT depths below the
   skull: scalp and skull thickness and the cortical surface offset are not
   included, so this cannot be converted to an optical dose without them.
2. Distance is measured inside the section plane. Where cortex curves sharply the
   nearest L1 cell can belong to an adjacent bank. MOp and VISp are comparatively
   flat, but `pia_reference_cells` is emitted per cell so sparse groups can be
   filtered.
3. Groups with fewer than `min_reference_cells` L1 cells get depth NaN rather than
   a fabricated value.
4. One animal. No animal-level inference, same as the rest of this pipeline.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

# Shell edges in micrometres. The last bin is open-ended.
DEFAULT_SHELL_EDGES = (0, 100, 200, 300, 400, 600, 800, 1000, 1200)
PIA_LAYER = "1"
MIN_REFERENCE_CELLS = 10
ROBUSTNESS_K = 5

DEPTH_CAUTION = (
    "Depth is distance in the coronal plane to the nearest CCF layer-1 cell of the same "
    "section, region and hemisphere. L1 somata lie below the true pia, so absolute depth is "
    "slightly underestimated; the offset is shared by both regions. Depth is measured from the "
    "cortical surface, not from the skull, and is not an optical dose. One C57BL/6J male donor."
)


def shell_labels(edges=DEFAULT_SHELL_EDGES):
    labels = [f"{a:04d}-{b:04d}" for a, b in zip(edges[:-1], edges[1:])]
    labels.append(f"{edges[-1]:04d}+")
    return labels


def assign_shell(depth_um, edges=DEFAULT_SHELL_EDGES):
    labels = shell_labels(edges)
    bins = list(edges) + [np.inf]
    out = pd.cut(depth_um, bins=bins, labels=labels, right=False, include_lowest=True)
    return out.astype(object).where(pd.notna(out), "unknown")


def cortical_depth(cells, min_reference_cells=MIN_REFERENCE_CELLS,
                   robustness_k=ROBUSTNESS_K, edges=DEFAULT_SHELL_EDGES):
    """Return `cells` with depth columns added. Input is not mutated.

    Added columns
        cortical_depth_um            distance to the nearest L1 reference cell
        cortical_depth_um_k{k}       distance to the k-th nearest, robustness only
        depth_shell                  binned `cortical_depth_um`
        pia_reference_cells          L1 cells available in this cell's group
        depth_reference_ok           group had enough L1 cells and valid coordinates
    """
    from scipy.spatial import cKDTree

    required = {"brain_section_label", "region", "hemisphere", "anatomical_layer", "y_ccf", "z_ccf"}
    missing = required - set(cells.columns)
    if missing:
        raise ValueError(f"cortical_depth needs columns {sorted(missing)}")

    out = cells.copy()
    out["cortical_depth_um"] = np.nan
    out[f"cortical_depth_um_k{robustness_k}"] = np.nan
    out["pia_reference_cells"] = 0

    coordinates = out[["y_ccf", "z_ccf"]].to_numpy(dtype=float)
    finite = np.isfinite(coordinates).all(axis=1)
    group_keys = ["brain_section_label", "region", "hemisphere"]
    is_pia = out["anatomical_layer"].astype(str).eq(PIA_LAYER).to_numpy()

    for _, index in out.groupby(group_keys, dropna=False, observed=True).indices.items():
        index = np.asarray(index)
        usable = index[finite[index]]
        reference = usable[is_pia[usable]]
        out.iloc[index, out.columns.get_loc("pia_reference_cells")] = len(reference)
        if len(reference) < min_reference_cells or usable.size == 0:
            continue
        tree = cKDTree(coordinates[reference])
        k = min(robustness_k, len(reference))
        distances, _ = tree.query(coordinates[usable], k=k)
        distances = np.atleast_2d(distances.T).T if distances.ndim == 1 else distances
        out.iloc[usable, out.columns.get_loc("cortical_depth_um")] = distances[:, 0] * 1000.0
        out.iloc[usable, out.columns.get_loc(f"cortical_depth_um_k{robustness_k}")] = distances[:, k - 1] * 1000.0

    out["depth_shell"] = assign_shell(out["cortical_depth_um"], edges)
    out["depth_reference_ok"] = out["cortical_depth_um"].notna()
    return out


def depth_layer_profile(cells):
    """Median and quartile depth of each CCF layer, the evidence for thickness."""
    usable = cells[cells.get("qc_annotated", True) & cells["cortical_depth_um"].notna()]
    rows = []
    for (region, layer), frame in usable.groupby(["region", "anatomical_layer"], observed=True):
        depth = frame["cortical_depth_um"]
        rows.append({
            "region": region, "anatomical_layer": layer, "cells": len(frame),
            "depth_um_p05": depth.quantile(0.05), "depth_um_p25": depth.quantile(0.25),
            "depth_um_median": depth.median(), "depth_um_p75": depth.quantile(0.75),
            "depth_um_p95": depth.quantile(0.95),
            "depth_um_median_robustness_k": frame.filter(like="cortical_depth_um_k").median().iloc[0],
        })
    columns = ["region", "anatomical_layer", "cells", "depth_um_p05", "depth_um_p25",
               "depth_um_median", "depth_um_p75", "depth_um_p95", "depth_um_median_robustness_k"]
    table = pd.DataFrame(rows, columns=columns).sort_values(["region", "anatomical_layer"])
    table["donor_label"] = "C57BL6J-638850"
    table["independent_animals"] = 1
    table["caution"] = DEPTH_CAUTION
    return table


def depth_cumulative(cells, step_um=50, max_um=1800):
    """Cumulative share of each region's cells shallower than a given depth.

    This is the table that shows an identical surface irradiance does not reach an
    equivalent share of the two regions.
    """
    usable = cells[cells.get("qc_annotated", True) & cells["cortical_depth_um"].notna()]
    rows = []
    for region, frame in usable.groupby("region", observed=True):
        depth = frame["cortical_depth_um"].to_numpy()
        for cut in range(0, max_um + step_um, step_um):
            rows.append({"region": region, "depth_um": cut, "cells_shallower": int((depth < cut).sum()),
                         "region_cells": len(depth), "cumulative_percentage": 100.0 * (depth < cut).mean()})
    table = pd.DataFrame(rows)
    table["donor_label"] = "C57BL6J-638850"
    table["independent_animals"] = 1
    return table


def cortex_depth_limit(cells, quantile=0.95):
    """Deepest cortical depth per region, taken as a quantile of the L6b depths.

    Below this a shell is mostly white matter and subcortical tissue, where an
    MOp-to-VISp composition ratio is not a cortical comparison. VISp reaches it at
    roughly half the depth MOp does, so shells have to be marked rather than
    silently compared.
    """
    usable = cells[cells.get("qc_annotated", True) & cells["cortical_depth_um"].notna()]
    deepest = usable[usable["anatomical_layer"].astype(str).eq("6b")]
    return deepest.groupby("region", observed=True)["cortical_depth_um"].quantile(quantile).to_dict()


def _shell_end(label):
    """Deep edge of a shell. The open-ended last bin and 'unknown' have no edge."""
    if not label or not label[0].isdigit() or label.endswith("+"):
        return np.inf
    return float(label.split("-")[1])


def depth_shell_ratios(cells, proportions_fn, min_shell_cells=200):
    """Neuron-normalised composition per depth shell, with the E to I ratio.

    Percentages of all QC cells share one denominator across every cell class, so a
    change in the neuron share mechanically moves every other class. Expressing glia
    per 100 neurons removes that closure effect; `proportions_fn` output is kept
    alongside for continuity with the existing tables.
    """
    usable = cells[cells.get("qc_annotated", True) & cells["cortical_depth_um"].notna()]
    counts = (usable.groupby(["region", "depth_shell", "major_group"], observed=True)
              .size().unstack("major_group", fill_value=0))
    glutamatergic = counts.get("Glutamatergic neuron", pd.Series(0, index=counts.index))
    gabaergic = counts.get("GABAergic neuron", pd.Series(0, index=counts.index))
    neurons = glutamatergic + gabaergic
    limits = cortex_depth_limit(cells)
    rows = []
    for key, row in counts.iterrows():
        region, shell = key
        total = int(row.sum())
        n_neurons = int(neurons.loc[key])
        if total < min_shell_cells or n_neurons == 0:
            continue
        record = {"region": region, "depth_shell": shell, "shell_cells": total, "shell_neurons": n_neurons,
                  "shell_sufficient": True,
                  "within_cortex_depth": bool(_shell_end(shell) <= limits.get(region, np.inf))}
        for group, value in row.items():
            record[f"{group} per 100 neurons"] = 100.0 * value / n_neurons
        record["E_to_I_ratio"] = (gabaergic.loc[key] / glutamatergic.loc[key]
                                  if glutamatergic.loc[key] else np.nan)
        rows.append(record)
    columns = ["region", "depth_shell", "shell_cells", "shell_neurons", "shell_sufficient",
               "within_cortex_depth", "E_to_I_ratio"]
    table = pd.DataFrame(rows, columns=columns if not rows else None)
    if len(table):
        table = table.sort_values(["depth_shell", "region"])
    table["denominator"] = "neurons_within_depth_shell"
    table["donor_label"] = "C57BL6J-638850"
    table["independent_animals"] = 1
    table["caution"] = DEPTH_CAUTION
    return table


def depth_shell_difference(shell_table):
    """MOp minus VISp within each depth shell, only where both regions qualify."""
    columns = ["depth_shell", "metric", "MOp", "VISp", "MOp_minus_VISp",
               "MOp_to_VISp_ratio", "both_regions_within_cortex"]
    if shell_table.empty:
        empty = pd.DataFrame(columns=columns)
        empty["independent_animals"] = pd.Series(dtype=int)
        empty["inference"] = pd.Series(dtype=object)
        return empty
    metrics = [c for c in shell_table.columns if c.endswith("per 100 neurons")] + ["E_to_I_ratio"]
    wide = shell_table.pivot(index="depth_shell", columns="region", values=metrics)
    cortical = shell_table.groupby("depth_shell", observed=True)["within_cortex_depth"].all()
    rows = []
    for shell in wide.index:
        for metric in metrics:
            mop = wide.loc[shell, (metric, "MOp")] if (metric, "MOp") in wide.columns else np.nan
            visp = wide.loc[shell, (metric, "VISp")] if (metric, "VISp") in wide.columns else np.nan
            if not (np.isfinite(mop) and np.isfinite(visp)):
                continue
            rows.append({"depth_shell": shell, "metric": metric, "MOp": mop, "VISp": visp,
                         "MOp_minus_VISp": mop - visp,
                         "MOp_to_VISp_ratio": mop / visp if visp else np.nan,
                         "both_regions_within_cortex": bool(cortical.get(shell, False))})
    table = pd.DataFrame(rows, columns=columns if not rows else None)
    table["independent_animals"] = 1
    table["inference"] = ("Descriptive single-animal difference. No p value, no FDR, no density. "
                          "Rows with both_regions_within_cortex False compare cortex against mostly "
                          "white matter and are not a cortical comparison.")
    return table
