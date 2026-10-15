#!/usr/bin/env python3
"""Physical depth shell composition for B6J MERFISH MOp/VISp, plus figures.

Runs off the saved cell table that run_merfish.py already produced, so none of the
multi-gigabyte Allen inputs have to be downloaded again.

    python src/merfish/run_depth_shells.py

Outputs
    results/merfish/composition_depth_shell_major.csv
    results/merfish/composition_depth_shell_subclass_all_cells.csv
    results/merfish/composition_depth_shell_subclass_neurons.csv
    results/merfish/composition_depth_shell_gaba_within_gaba.csv
    results/merfish/composition_depth_shell_per100_neurons.csv
    results/merfish/depth_shell_regional_difference.csv
    results/merfish/depth_layer_profile.csv
    results/merfish/depth_cumulative_distribution.csv
    results/merfish/selected_cells_with_depth.csv.gz
    results/merfish/depth_report.json
    figures/merfish/depth_*.png
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from depth import (DEFAULT_SHELL_EDGES, DEPTH_CAUTION, cortical_depth, depth_cumulative,
                   depth_layer_profile, depth_shell_difference, depth_shell_ratios, shell_labels)

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "results" / "merfish"
FIGURES = ROOT / "figures" / "merfish"
CELLS = RESULTS / "selected_cells_raw_labels_and_coordinates.csv.gz"
LAYERS = ["1", "2/3", "4", "5", "6a", "6b"]
REGION_COLOR = {"MOp": "#1f5fa8", "VISp": "#d2691e"}


def proportions(cells, dimensions, category, denominator="all_QC_annotated"):
    """Same contract as run_merfish.proportions, duplicated to keep this runnable alone."""
    if denominator == "neurons":
        cells = cells[cells.is_neuron]
    elif denominator == "GABAergic_neurons":
        cells = cells[cells.is_gaba]
    elif denominator != "all_QC_annotated":
        raise ValueError(denominator)
    totals = cells.groupby(dimensions, observed=True, dropna=False).size().rename("denominator_cells")
    counts = cells.groupby(dimensions + [category], observed=True, dropna=False).size().rename("cell_count")
    categories = sorted(cells[category].dropna().unique())
    expanded = totals.reset_index()[dimensions].merge(pd.DataFrame({category: categories}), how="cross")
    result = expanded.merge(counts.reset_index(), on=dimensions + [category], how="left")
    result["cell_count"] = result.cell_count.fillna(0).astype(int)
    result = result.merge(totals.reset_index(), on=dimensions, validate="many_to_one")
    result["fraction"] = result.cell_count / result.denominator_cells
    result["percentage"] = 100 * result.fraction
    result["denominator"] = denominator
    result["donor_label"] = "C57BL6J-638850"
    result["independent_animals"] = 1
    return result


def setup_matplotlib():
    cache = ROOT / ".cache"
    cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(cache / "matplotlib"))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def plot_depth(retained, layer_table, cumulative, shell_table, figdir):
    plt = setup_matplotlib()
    figdir.mkdir(parents=True, exist_ok=True)

    # 1. Layer depth ranges. The thickness mismatch is the whole argument, so it leads.
    fig, ax = plt.subplots(figsize=(8, 5))
    for offset, region in zip([-0.18, 0.18], ["MOp", "VISp"]):
        sub = layer_table[layer_table.region.eq(region)].set_index("anatomical_layer").reindex(LAYERS)
        y = np.arange(len(LAYERS)) + offset
        ok = sub.depth_um_median.notna().to_numpy()
        ax.hlines(y[ok], sub.depth_um_p25[ok], sub.depth_um_p75[ok],
                  color=REGION_COLOR[region], linewidth=7, alpha=0.35)
        ax.plot(sub.depth_um_median[ok], y[ok], "o", color=REGION_COLOR[region], label=region)
        for yi, (p50, n) in zip(y[ok], zip(sub.depth_um_median[ok], sub.cells[ok])):
            ax.annotate(f"{p50:.0f} (n={int(n):,})", (p50, yi), textcoords="offset points",
                        xytext=(8, -3), fontsize=7, color=REGION_COLOR[region])
    ax.set_yticks(np.arange(len(LAYERS)), [f"L{x}" for x in LAYERS])
    ax.invert_yaxis()
    ax.set(xlabel="depth below pia (um), median with interquartile range",
           title="MOp cortex is about twice as thick as VISp\nCCF layer labels are not depth matched")
    ax.legend(title="region")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout(); fig.savefig(figdir / "depth_layer_profile.png", dpi=180); plt.close(fig)

    # 2. Cumulative reach. Same surface irradiance, different share of tissue covered.
    fig, ax = plt.subplots(figsize=(8, 5))
    for region in ["MOp", "VISp"]:
        sub = cumulative[cumulative.region.eq(region)]
        ax.plot(sub.depth_um, sub.cumulative_percentage, color=REGION_COLOR[region], label=region, linewidth=2)
    for cut in (400, 800):
        ax.axvline(cut, color="grey", linestyle=":", linewidth=1)
        for region in ["MOp", "VISp"]:
            sub = cumulative[cumulative.region.eq(region) & cumulative.depth_um.eq(cut)]
            if len(sub):
                ax.annotate(f"{sub.cumulative_percentage.iloc[0]:.0f}%", (cut, sub.cumulative_percentage.iloc[0]),
                            textcoords="offset points", xytext=(6, -10), fontsize=8, color=REGION_COLOR[region])
    ax.set(xlabel="depth below pia (um)", ylabel="% of region cells shallower than depth",
           ylim=(0, 101), title="Cumulative cell depth distribution\nDepth is from the cortical surface, not the skull")
    ax.legend(title="region"); ax.grid(alpha=0.25)
    fig.tight_layout(); fig.savefig(figdir / "depth_cumulative_distribution.png", dpi=180); plt.close(fig)

    # 3. Depth-matched composition, normalised to neurons to avoid the closure artefact.
    metrics = ["Astrocyte per 100 neurons", "Oligodendrocyte per 100 neurons",
               "Endothelial per 100 neurons", "GABAergic neuron per 100 neurons"]
    metrics = [m for m in metrics if m in shell_table.columns]
    shells = [s for s in shell_labels() if s in set(shell_table.depth_shell)]
    fig, axes = plt.subplots(1, len(metrics), figsize=(3.4 * len(metrics), 4.6), sharey=True)
    axes = np.atleast_1d(axes)
    cortical = shell_table.groupby("depth_shell", observed=True).within_cortex_depth.all().reindex(shells).fillna(False)
    first_outside = int(np.argmax(~cortical.to_numpy())) if (~cortical).any() else None
    for ax, metric in zip(axes, metrics):
        for region in ["MOp", "VISp"]:
            sub = shell_table[shell_table.region.eq(region)].set_index("depth_shell").reindex(shells)
            ax.plot(sub[metric], np.arange(len(shells)), "o-", color=REGION_COLOR[region], label=region)
        if first_outside is not None:
            ax.axhspan(first_outside - 0.5, len(shells) - 0.5, color="lightgrey", alpha=0.45, zorder=0)
        ax.set_title(metric.replace(" per 100 neurons", ""), fontsize=10)
        ax.grid(alpha=0.25)
    axes[0].set_yticks(np.arange(len(shells)), [s.replace("-", " to ") for s in shells], fontsize=8)
    axes[0].invert_yaxis()
    axes[0].set_ylabel("depth shell below pia (um)")
    axes[0].legend(title="region", fontsize=8)
    fig.suptitle("Depth-matched composition, per 100 neurons in the same shell (one animal)")
    fig.supxlabel("cells per 100 neurons within shell", fontsize=9)
    fig.tight_layout(); fig.savefig(figdir / "depth_shell_composition.png", dpi=180); plt.close(fig)

    # 4. E to I ratio. The quantity that tracks the PBM dose gradient.
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    width = 0.36
    x = np.arange(len(shells))
    cortical = shell_table.groupby("depth_shell", observed=True).within_cortex_depth.all().reindex(shells).fillna(False)
    for offset, region in zip([-width / 2, width / 2], ["MOp", "VISp"]):
        sub = shell_table[shell_table.region.eq(region)].set_index("depth_shell").reindex(shells)
        inside = sub.within_cortex_depth.fillna(False).to_numpy(dtype=bool)
        ax.bar(x + offset, sub.E_to_I_ratio, width, color=REGION_COLOR[region], label=region)
        # Shells past the region's own L6b depth are mostly white matter; hatch them.
        ax.bar(x[~inside] + offset, sub.E_to_I_ratio.to_numpy()[~inside], width,
               color="none", edgecolor="white", hatch="///", linewidth=0)
    ratio = (shell_table.pivot(index="depth_shell", columns="region", values="E_to_I_ratio")
             .reindex(shells))
    if {"MOp", "VISp"}.issubset(ratio.columns):
        for xi, value, ok in zip(x, (ratio.MOp / ratio.VISp).to_numpy(), cortical.to_numpy()):
            if np.isfinite(value) and ok:
                ax.annotate(f"x{value:.2f}", (xi, 0.005), ha="center", fontsize=8, color="black")
    ax.annotate("hatched = below that region's L6b, mostly white matter",
                (0.98, 0.72), xycoords="axes fraction", ha="right", fontsize=8, color="dimgrey")
    ax.set_xticks(x, [s.replace("-", "\nto ") for s in shells], fontsize=8)
    ax.set(ylabel="GABAergic / glutamatergic neurons", xlabel="depth shell below pia (um)",
           title="Inhibitory to excitatory ratio by physical depth\nLabel is the MOp to VISp ratio in that shell")
    ax.legend(title="region"); ax.grid(axis="y", alpha=0.25)
    fig.tight_layout(); fig.savefig(figdir / "depth_shell_e_to_i_ratio.png", dpi=180); plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cells", type=Path, default=CELLS)
    parser.add_argument("--min-shell-cells", type=int, default=200)
    parser.add_argument("--no-figures", action="store_true")
    args = parser.parse_args(argv)

    RESULTS.mkdir(parents=True, exist_ok=True)
    report = {"run_time_UTC": datetime.now(timezone.utc).isoformat(), "status": "pending",
              "dataset": "MERFISH-C57BL6J-638850", "caution": DEPTH_CAUTION,
              "shell_edges_um": list(DEFAULT_SHELL_EDGES), "pending": []}
    try:
        if not args.cells.exists():
            raise FileNotFoundError(f"Run src/merfish/run_merfish.py first: {args.cells} is absent")
        cells = pd.read_csv(args.cells, dtype={"cell_label": str, "anatomical_layer": str}, low_memory=False)
        cells = cortical_depth(cells)
        cells.to_csv(RESULTS / "selected_cells_with_depth.csv.gz", index=False)

        retained = cells[cells.qc_annotated.astype(str).isin(["True", "TRUE", "1"])].copy()
        report["cells"] = {
            "selected_cells": int(len(cells)), "QC_annotated": int(len(retained)),
            "with_depth": int(retained.cortical_depth_um.notna().sum()),
            "without_depth_insufficient_pia_reference": int(retained.cortical_depth_um.isna().sum()),
            "pia_reference_groups": int(retained.groupby(
                ["brain_section_label", "region", "hemisphere"], observed=True).ngroups),
        }

        for category, denominator, suffix in [
                ("major_group", "all_QC_annotated", "major"),
                ("raw_subclass", "all_QC_annotated", "subclass_all_cells"),
                ("raw_subclass", "neurons", "subclass_neurons"),
                ("gaba_family", "GABAergic_neurons", "gaba_within_gaba")]:
            table = proportions(retained[retained.cortical_depth_um.notna()],
                                ["region", "depth_shell"], category, denominator)
            table["caution"] = DEPTH_CAUTION
            table.to_csv(RESULTS / f"composition_depth_shell_{suffix}.csv", index=False)

        layer_table = depth_layer_profile(retained)
        layer_table.to_csv(RESULTS / "depth_layer_profile.csv", index=False)
        cumulative = depth_cumulative(retained)
        cumulative.to_csv(RESULTS / "depth_cumulative_distribution.csv", index=False)
        shell_table = depth_shell_ratios(retained, proportions, args.min_shell_cells)
        shell_table.to_csv(RESULTS / "composition_depth_shell_per100_neurons.csv", index=False)
        difference = depth_shell_difference(shell_table)
        difference.to_csv(RESULTS / "depth_shell_regional_difference.csv", index=False)

        thickness = layer_table[layer_table.anatomical_layer.eq("6b")].set_index("region").depth_um_median
        report["cortical_thickness_proxy_L6b_median_depth_um"] = {
            region: float(value) for region, value in thickness.items()}
        report["shells_reported"] = sorted(set(shell_table.depth_shell))
        report["shells_dropped_below_min_cells"] = args.min_shell_cells

        if not args.no_figures:
            plot_depth(retained, layer_table, cumulative, shell_table, FIGURES)
            report["figures"] = ["depth_layer_profile.png", "depth_cumulative_distribution.png",
                                 "depth_shell_composition.png", "depth_shell_e_to_i_ratio.png"]
        report["status"] = "depth_shell_composition_completed"
        return 0
    except (OSError, ValueError, KeyError) as error:
        report["pending"].append(f"{type(error).__name__}: {error}")
        print(report["pending"][-1], file=sys.stderr)
        return 2
    finally:
        (RESULTS / "depth_report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False))
        print(json.dumps({"status": report["status"], "pending": report["pending"]}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
