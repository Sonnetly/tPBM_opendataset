#!/usr/bin/env python3
"""Explicit age-unmatched exploratory screen of the P57–P59 10x cohort.

Primary exact-age results remain in analyze_expression.py outputs. These
supplementary tables retain recorded donor ages and date/region confounding;
they contain descriptive averages only and cannot identify a region effect.
"""
import argparse
import json
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

DEFAULT_ROOT = Path(__file__).resolve().parents[2]
KEYS = ["chemistry", "genotype", "sex", "facs", "class", "subclass", "gene"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=DEFAULT_ROOT)
    args = parser.parse_args()
    root = args.project_root.resolve()
    out = root / "results/rnaseq"
    source = root / "input/rnaseq/ctx_hpf_mop_visp_candidate_counts.h5ad"
    data = ad.read_h5ad(source)
    metadata = data.obs.copy().reset_index()
    counts = data.X.toarray()
    totals = metadata.full_library_counts.to_numpy(dtype=float)
    expression = np.log1p(counts / totals[:, None] * 1e4)
    genes = data.var_names.tolist()
    panel = set(pd.read_csv(root / "results/pbm/panel22.csv").gene)
    common = pd.read_csv(root / "results/pbm/same_direction_52.csv")
    direction_column = "hippocampus_direction" if "hippocampus_direction" in common else "A_direction"
    direction = common.set_index("gene")[direction_column].to_dict()
    if set(metadata.age) != {"P57", "P58", "P59"}:
        raise RuntimeError("Exploratory cohort age scope differs from the explicitly declared P57–P59 screen")
    fixed = ["chemistry", "genotype", "sex", "facs"]
    if not metadata[fixed].nunique().eq(1).all():
        raise RuntimeError("Unexpected mixed chemistry/genotype/sex/FACS; refine cohort scope explicitly")
    rows = []
    grouping = ["region"] + fixed + ["class", "subclass"]
    for labels, positions in metadata.groupby(grouping, observed=True).indices.items():
        positions = np.asarray(positions)
        group_metadata = metadata.iloc[positions]
        for j, gene in enumerate(genes):
            rows.append(dict(zip(grouping, labels), gene=gene,
                n_cells=len(positions), n_recorded_donors=group_metadata.donor.nunique(),
                donor_labels=";".join(sorted(group_metadata.donor.astype(str).unique())),
                ages=";".join(sorted(group_metadata.age.astype(str).unique())),
                preparation_dates=";".join(sorted(group_metadata.batch.astype(str).unique())),
                mean_log1p_CP10k=float(expression[positions, j].mean()),
                mean_linear_CP10k=float((counts[positions, j] / totals[positions] * 1e4).mean()),
                detection_fraction=float((counts[positions, j] > 0).mean()),
                panel22=gene in panel, common52_direction=direction.get(gene, "not_common52"),
                scope="Exploratory P57–P59 sampled-cell cohort; ages not exactly matched; region/preparation date confounded"))
    summary = pd.DataFrame(rows)
    summary.to_csv(out / "exploratory_cohort_region_subclass_gene_summary.csv", index=False)
    comparison = summary.loc[summary.region.eq("MOp")].merge(summary.loc[summary.region.eq("VISp")], on=KEYS, how="outer", suffixes=("_MOp", "_VISp"), indicator=True)
    comparison["both_regions_observed"] = comparison["_merge"].eq("both")
    comparison["MOp_minus_VISp_mean_log1p_CP10k"] = comparison.mean_log1p_CP10k_MOp - comparison.mean_log1p_CP10k_VISp
    comparison["MOp_minus_VISp_detection_fraction"] = comparison.detection_fraction_MOp - comparison.detection_fraction_VISp
    comparison["interpretation"] = "Observed descriptive contrast; age unmatched and preparation date confounded; no P/FDR or causal region effect"
    comparison.drop(columns="_merge").to_csv(out / "exploratory_cohort_within_subclass_region_difference.csv", index=False)
    cluster_rows = []
    cluster_grouping = grouping + ["cluster"]
    known_cluster = metadata.cluster.notna() & ~metadata.cluster.astype(str).str.lower().isin(["unknown", "na", "", "nan"])
    for labels, group in metadata.loc[known_cluster].groupby(cluster_grouping, observed=True):
        positions = group.index.to_numpy(dtype=int)
        for j, gene in enumerate(genes):
            cluster_rows.append(dict(zip(cluster_grouping, labels), gene=gene, n_cells=len(positions), mean_log1p_CP10k=float(expression[positions, j].mean())))
    clusters = pd.DataFrame(cluster_rows)
    shared = clusters.loc[clusters.region.eq("MOp")].merge(clusters.loc[clusters.region.eq("VISp")], on=KEYS + ["cluster"], suffixes=("_MOp", "_VISp"))
    shared["MOp_minus_VISp"] = shared.mean_log1p_CP10k_MOp - shared.mean_log1p_CP10k_VISp
    shared.to_csv(out / "exploratory_cohort_shared_cluster_region_difference.csv", index=False)
    standardized = shared.groupby(KEYS, observed=True).agg(
        n_shared_clusters=("cluster", "nunique"),
        equal_cluster_weighted_MOp_mean=("mean_log1p_CP10k_MOp", "mean"),
        equal_cluster_weighted_VISp_mean=("mean_log1p_CP10k_VISp", "mean"),
        equal_cluster_weighted_MOp_minus_VISp=("MOp_minus_VISp", "mean"),
        n_positive_cluster_contrasts=("MOp_minus_VISp", lambda x: int(x.gt(0).sum())),
        n_negative_cluster_contrasts=("MOp_minus_VISp", lambda x: int(x.lt(0).sum())),
    ).reset_index()
    standardized["scope"] = "Equal matched-cluster weights; exploratory age-unmatched P57–P59 cohort; date confounded; no animal inference"
    standardized.to_csv(out / "exploratory_cohort_shared_cluster_standardized_difference.csv", index=False)
    donor_rows = []
    for labels, positions in metadata.groupby(["region", "donor", "age", "batch", "class", "subclass"], observed=True).indices.items():
        positions = np.asarray(positions)
        for j, gene in enumerate(genes):
            donor_rows.append(dict(zip(["region", "donor", "age", "batch", "class", "subclass"], labels),
                gene=gene, n_cells=len(positions), mean_log1p_CP10k=float(expression[positions, j].mean()),
                detection_fraction=float((counts[positions, j] > 0).mean()),
                candidate_raw_count_sum=int(counts[positions, j].sum()),
                full_library_count_sum=int(totals[positions].sum()),
                inference="Recorded donor summaries only; pool independence and confounding preclude animal tests"))
    pd.DataFrame(donor_rows).to_csv(out / "recorded_donor_subclass_candidate_summaries.csv", index=False)
    os.environ.setdefault("MPLCONFIGDIR", str(root / ".cache/matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figdir = root / "figures/rnaseq"
    figdir.mkdir(parents=True, exist_ok=True)
    for name, selection in [("panel22", genes and sorted(panel)), ("common_up34", sorted(g for g, d in direction.items() if d == "up")), ("common_down18", sorted(g for g, d in direction.items() if d == "down"))]:
        frame = summary.loc[summary.gene.isin(selection)]
        labels = sorted(set(zip(frame.region, frame["class"], frame.subclass)))
        row_counts = frame.groupby(["region", "class", "subclass"])["n_cells"].first().to_dict()
        xmap, ymap = {g: i for i, g in enumerate(selection)}, {label: i for i, label in enumerate(labels)}
        fig, ax = plt.subplots(figsize=(max(11, len(selection) * .40), max(5, len(labels) * .29)))
        points = ax.scatter(frame.gene.map(xmap), [ymap[(r, c, s)] for r, c, s in zip(frame.region, frame["class"], frame.subclass)],
            s=4 + frame.detection_fraction.to_numpy() * 120, c=frame.mean_log1p_CP10k, cmap="viridis")
        ax.set_xticks(range(len(selection)), selection, rotation=65, ha="right")
        ax.set_yticks(range(len(labels)), [f"{r} | {c} | {s} (n={int(row_counts[(r, c, s)]):,})" for r, c, s in labels], fontsize=7)
        ax.set_title(f"{name}: exploratory P57–P59 cohort; ages not exactly matched\nMale Snap25;Ai14, 10x v2; region/preparation date confounded", fontsize=10)
        fig.colorbar(points, ax=ax, label="mean log1p(CP10k), full-gene library denominator")
        for fraction in [.25, .5, 1.0]:
            ax.scatter([], [], s=4 + fraction * 120, c="gray", label=f"{fraction:.0%} detected")
        ax.legend(loc="upper left", bbox_to_anchor=(1.2, 1), frameon=False, fontsize=8)
        fig.tight_layout()
        fig.savefig(figdir / f"exploratory_cohort_{name}_dotplot.png", dpi=160)
        fig.savefig(figdir / f"exploratory_cohort_{name}_dotplot.pdf")
        plt.close(fig)
    status = {"status": "exploratory_cohort_complete", "n_cells": len(metadata), "n_genes": len(genes), "n_recorded_donors": int(metadata.donor.nunique()), "ages": sorted(metadata.age.unique().tolist()), "strict_age_comparison": "Separate primary within_subclass_region_difference.csv; exact shared age P58 only", "animal_inference": "not performed; regional library preparation dates are disjoint and independent pools not verified"}
    (out / "exploratory_cohort_status.json").write_text(json.dumps(status, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(status, ensure_ascii=False))


if __name__ == "__main__":
    main()
