#!/usr/bin/env python3
"""Plot subclass-level versus supertype-standardized MOp minus VISp for chosen genes.

Reads results/zhuang/supertype/candidate_supertype_decomposition.csv written by
zhuang_supertype_check.py and writes figures/zhuang/zhuang_supertype_standardized.png.
Hollow marker = plain difference over the kept supertypes. Filled marker = supertype
standardized difference. The gap between them is the within-subclass composition part.
Marker size scales with coverage (share of subclass cells in kept supertypes).
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
TABLE = ROOT / "results" / "zhuang" / "supertype" / "candidate_supertype_decomposition.csv"
OUT = ROOT / "figures" / "zhuang" / "zhuang_supertype_standardized.png"
ANIMAL_COLORS = {"Zhuang-ABCA-1": "#1f5fa8", "Zhuang-ABCA-2": "#5aa0d8",
                 "Zhuang-ABCA-3": "#c2571a", "Zhuang-ABCA-4": "#e8a06a"}


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--genes", nargs="+", default=["Adora1", "Zeb2"])
    args = parser.parse_args(argv)
    table = pd.read_csv(TABLE)
    fig, axes = plt.subplots(1, len(args.genes), figsize=(6.2 * len(args.genes), 6.4), sharey=True)
    axes = [axes] if len(args.genes) == 1 else axes
    order = (table[table.cell_class.ne("Other")]
             .sort_values(["cell_class", "subclass_clean"]).subclass_clean.drop_duplicates().tolist())
    for ax, gene in zip(axes, args.genes):
        sub = table[table.gene.eq(gene) & table.subclass_clean.isin(order)]
        for k, (animal, color) in enumerate(ANIMAL_COLORS.items()):
            frame = sub[sub.animal.eq(animal)]
            y = [order.index(s) + (k - 1.5) * 0.17 for s in frame.subclass_clean]
            size = 15 + 55 * frame.coverage
            ax.scatter(frame.raw_restricted, y, s=size, facecolors="none", edgecolors=color, linewidths=1.2)
            ax.scatter(frame.standardized, y, s=size, color=color, label=animal)
            for x0, x1, yy in zip(frame.raw_restricted, frame.standardized, y):
                ax.plot([x0, x1], [yy, yy], color=color, linewidth=0.8, alpha=0.7)
        ax.axvline(0, color="black", linewidth=0.8)
        ax.set_yticks(range(len(order)), order, fontsize=8)
        ax.set_title(gene)
        ax.set_xlabel("MOp minus VISp (official log2 units)")
        ax.grid(axis="x", alpha=0.25)
    axes[0].legend(fontsize=8, loc="lower left")
    fig.suptitle("Hollow = plain difference over kept supertypes, filled = supertype standardized; "
                 "size = coverage", fontsize=9)
    fig.tight_layout()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=180)
    print(OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
