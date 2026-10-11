#!/usr/bin/env python3
"""Export auditable core-panel and sensitivity summaries from measured results."""
import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_ROOT = Path(__file__).resolve().parents[2]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=DEFAULT_ROOT)
    args = parser.parse_args()
    root = args.project_root.resolve()
    out = root / "results/rnaseq"
    panel = pd.read_csv(root / "results/pbm/panel22.csv")
    strict = pd.read_csv(out / "within_subclass_region_difference.csv")
    strict = strict.loc[strict.both_regions_observed].copy()
    if set(strict.age) != {"P58"}:
        raise RuntimeError("Exact age-matched report scope changed; review before generating interpretation")
    strict_clusters = pd.read_csv(out / "shared_cluster_standardized_difference.csv")
    strict_clusters = strict_clusters.loc[strict_clusters.age.eq("P58")]
    cohort = pd.read_csv(out / "exploratory_cohort_within_subclass_region_difference.csv")
    cohort = cohort.loc[cohort.both_regions_observed].copy()
    cohort_clusters = pd.read_csv(out / "exploratory_cohort_shared_cluster_standardized_difference.csv")
    keys = ["class", "subclass", "gene"]
    strict_cols = keys + ["age", "n_cells_MOp", "n_cells_VISp", "mean_expression_MOp", "mean_expression_VISp", "detection_fraction_MOp", "detection_fraction_VISp", "MOp_minus_VISp_mean_expression"]
    screening = strict[strict_cols].merge(strict_clusters[keys + ["n_shared_clusters", "equal_cluster_weighted_MOp_minus_VISp"]], on=keys, how="left")
    screening = screening.rename(columns={"MOp_minus_VISp_mean_expression": "P58_cell_weighted_MOp_minus_VISp", "equal_cluster_weighted_MOp_minus_VISp": "P58_shared_cluster_equal_weight_MOp_minus_VISp", "n_shared_clusters": "P58_n_shared_clusters"})
    screening = screening.merge(cohort[keys + ["MOp_minus_VISp_mean_log1p_CP10k"]].rename(columns={"MOp_minus_VISp_mean_log1p_CP10k": "P57_59_cell_weighted_MOp_minus_VISp"}), on=keys, how="left")
    screening = screening.merge(cohort_clusters[keys + ["n_shared_clusters", "equal_cluster_weighted_MOp_minus_VISp"]].rename(columns={"n_shared_clusters": "P57_59_n_shared_clusters", "equal_cluster_weighted_MOp_minus_VISp": "P57_59_shared_cluster_equal_weight_MOp_minus_VISp"}), on=keys, how="left")
    screening["P58_cell_vs_cluster_comparison_available"] = np.isfinite(screening.P58_cell_weighted_MOp_minus_VISp) & np.isfinite(screening.P58_shared_cluster_equal_weight_MOp_minus_VISp)
    screening["P58_vs_P57_59_cluster_comparison_available"] = np.isfinite(screening.P58_shared_cluster_equal_weight_MOp_minus_VISp) & np.isfinite(screening.P57_59_shared_cluster_equal_weight_MOp_minus_VISp)
    screening["P58_cell_vs_cluster_direction_changes"] = (np.sign(screening.P58_cell_weighted_MOp_minus_VISp) != np.sign(screening.P58_shared_cluster_equal_weight_MOp_minus_VISp)).astype("boolean").where(screening.P58_cell_vs_cluster_comparison_available, pd.NA)
    screening["P58_vs_P57_59_cluster_direction_changes"] = (np.sign(screening.P58_shared_cluster_equal_weight_MOp_minus_VISp) != np.sign(screening.P57_59_shared_cluster_equal_weight_MOp_minus_VISp)).astype("boolean").where(screening.P58_vs_P57_59_cluster_comparison_available, pd.NA)
    screening["units"] = "Difference of mean log1p(CP10k); full-source-gene library denominator; not log2FC"
    screening["inference"] = "Descriptive; region/preparation date confounded; P58 MOp recorded donor n=1, VISp n=2; independent pools unverified"
    screening.to_csv(out / "all62_P58_cohort_cluster_sensitivity.csv", index=False)
    screening.loc[screening.gene.isin(panel.gene)].to_csv(out / "panel22_P58_cohort_cluster_sensitivity.csv", index=False)
    rows = []
    for gene in panel.gene:
        selected = screening.loc[screening.gene.eq(gene) & screening["class"].isin(["GABAergic", "Glutamatergic"])]
        rows.append({"gene": gene, "n_age_matched_neuron_subclasses": len(selected),
            "max_detection_fraction_MOp": selected.detection_fraction_MOp.max(),
            "max_detection_fraction_VISp": selected.detection_fraction_VISp.max(),
            "min_detection_fraction_MOp": selected.detection_fraction_MOp.min(),
            "min_detection_fraction_VISp": selected.detection_fraction_VISp.min(),
            "n_positive_P58_shared_cluster_contrasts": int(selected.P58_shared_cluster_equal_weight_MOp_minus_VISp.gt(0).sum()),
            "n_negative_P58_shared_cluster_contrasts": int(selected.P58_shared_cluster_equal_weight_MOp_minus_VISp.lt(0).sum()),
            "n_P58_cell_vs_shared_cluster_sign_changes": int(selected.P58_cell_vs_cluster_direction_changes.sum()),
            "n_P58_vs_cohort_shared_cluster_sign_changes": int(selected.P58_vs_P57_59_cluster_direction_changes.sum()),
            "interpretation": "Counts of descriptive contrast directions only; magnitude/detection must be inspected; no significance or qPCR priority rank"})
    pd.DataFrame(rows).to_csv(out / "panel22_screening_coverage_and_sensitivity.csv", index=False)
    core_pairs = [("Adora1", "L2/3 IT CTX"), ("Adora1", "Pvalb"), ("Adora2a", "L6 CT CTX"),
        ("Nt5e", "Pvalb"), ("Nt5e", "L2/3 IT CTX"), ("Grm2", "L4/5 IT CTX"),
        ("Grm2", "L2/3 IT CTX"), ("Penk", "Vip"), ("Penk", "L5 IT CTX"),
        ("Kcnj6", "L6 CT CTX"), ("Adk", "L4/5 IT CTX"), ("Slc29a2", "L2/3 IT CTX")]
    spotlight = pd.concat([screening.loc[screening.gene.eq(gene) & screening.subclass.eq(subclass)] for gene, subclass in core_pairs])
    spotlight.to_csv(out / "core_gene_observation_examples.csv", index=False)
    status = {"status": "measured_findings_summary_complete", "n_genes": screening.gene.nunique(), "n_core_panel_genes": len(panel), "P58_matched_subclasses": screening.subclass.nunique(), "animal_Pvalues": "not computed", "qPCR_rank": "not assigned; observations remain confounded screening candidates"}
    (out / "findings_summary_status.json").write_text(json.dumps(status, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(status, ensure_ascii=False))


if __name__ == "__main__":
    main()
