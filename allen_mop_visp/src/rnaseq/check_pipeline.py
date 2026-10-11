#!/usr/bin/env python3
"""Verify normalization, exact ROI selection and matched-cluster summaries.

Fixtures live only in a temporary directory. These checks do not validate any
Allen biological conclusion and are never stored as project RNA-seq results.
"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def main():
    script = Path(__file__).with_name("analyze_expression.py")
    with tempfile.TemporaryDirectory(prefix="allen-rnaseq-check-") as temp:
        root = Path(temp)
        (root / "results/pbm").mkdir(parents=True)
        pd.DataFrame({"gene": ["Entpd1", "Ada"]}).to_csv(root / "results/pbm/panel22.csv", index=False)
        pd.DataFrame({"gene": ["Grm2"]}).to_csv(root / "results/pbm/same_direction_52.csv", index=False)
        obs = pd.DataFrame({
            "region_label": ["MOp", "MOp", "VISp", "VISp", "VIS"],
            "subclass_label": ["IT"] * 5,
            "class_label": ["Glutamatergic"] * 5,
            "cluster_label": ["c1", "c2", "c1", "c2", "c1"],
            "chemistry": ["10x_v3"] * 5,
            "genotype": ["B6J reporter"] * 5,
            "sex": ["M"] * 5,
            "age": ["P56"] * 5,
            "facs": ["NeuN+"] * 5,
        }, index=[f"cell{i}" for i in range(5)])
        # The final Background gene is essential to the library denominator.
        counts = np.array([[1, 0, 9], [3, 1, 6], [0, 1, 9], [2, 0, 8], [100, 0, 0]])
        adata = ad.AnnData(sparse.csr_matrix(counts), obs=obs, var=pd.DataFrame(index=["Entpd1", "Grm2", "Background"]))
        path = root / "counts.h5ad"
        adata.write_h5ad(path)
        command = [sys.executable, str(script), "--project-root", str(root), "--h5ad", str(path), "--value-kind", "raw_counts", "--full-transcriptome-counts-confirmed"]
        subprocess.run(command + ["--chunk-size", "2"], check=True, capture_output=True, text=True)
        result = root / "results/rnaseq"
        status = json.loads((result / "status.json").read_text())
        assert status["n_selected_cells"] == 4, "VIS must not substitute for VISp"
        assert status["n_candidates_present"] == 2
        assert "not performed" in status["animal_statistics"]
        coverage = pd.read_csv(result / "candidate_gene_coverage.csv").set_index("gene")
        assert not bool(coverage.loc["Ada", "present_in_matrix"]), "Missing gene is unavailable, not zero"
        summary = pd.read_csv(result / "region_subclass_gene_summary.csv")
        mop_entpd1 = summary.loc[(summary.region == "MOp") & (summary.gene == "Entpd1")].iloc[0]
        assert np.isclose(mop_entpd1.mean_expression, np.mean(np.log1p([1000, 3000])))
        visp_entpd1 = summary.loc[(summary.region == "VISp") & (summary.gene == "Entpd1")].iloc[0]
        assert visp_entpd1.detection_fraction == .5
        standardized = pd.read_csv(result / "shared_cluster_standardized_difference.csv")
        assert standardized.n_shared_clusters.eq(2).all()
        saved = np.load(result / "selected_candidate_expression.npz")
        assert saved["library_total"].tolist() == [10, 10, 10, 10]
        # A verified external full-library total also supports feature subsets.
        adata.obs["full_library_umis"] = 20
        adata.X = sparse.csr_matrix(np.array([[0, 0, 0], *counts[1:]]))
        adata.write_h5ad(root / "external_totals.h5ad")
        external_command = [x for x in command if x != "--full-transcriptome-counts-confirmed"]
        external_command[external_command.index(str(path))] = str(root / "external_totals.h5ad")
        subprocess.run(external_command + ["--library-total-column", "full_library_umis"], check=True, capture_output=True, text=True)
        external_summary = pd.read_csv(result / "region_subclass_gene_summary.csv")
        target_row = external_summary.loc[(external_summary.region == "MOp") & (external_summary.gene == "Entpd1")].iloc[0]
        assert np.isclose(target_row.mean_expression, np.log1p(1500) / 2)
        adata.X = sparse.csr_matrix(counts)
        # Unknown cluster labels must never become a matched subtype.
        adata.obs.loc[["cell1", "cell3"], "cluster_label"] = np.nan
        adata.write_h5ad(root / "unknown_cluster.h5ad")
        cluster_command = command.copy()
        cluster_command[cluster_command.index(str(path))] = str(root / "unknown_cluster.h5ad")
        subprocess.run(cluster_command, check=True, capture_output=True, text=True)
        standardized = pd.read_csv(result / "shared_cluster_standardized_difference.csv")
        assert standardized.n_shared_clusters.eq(1).all()
        cluster_audit = json.loads((result / "shared_cluster_audit.json").read_text())
        assert cluster_audit["n_cells_without_cluster_excluded_from_standardization"] == 2
        # Restore complete taxonomy for the raw-count validation checks.
        adata.obs["cluster_label"] = ["c1", "c2", "c1", "c2", "c1"]
        # Non-integer normalized data cannot be claimed to be raw counts.
        adata.X = adata.X.astype(float) + sparse.csr_matrix(np.full(counts.shape, .5))
        adata.write_h5ad(root / "bad_counts.h5ad")
        bad_command = command.copy()
        bad_command[bad_command.index(str(path))] = str(root / "bad_counts.h5ad")
        bad = subprocess.run(bad_command, capture_output=True, text=True)
        assert bad.returncode != 0 and "non-integer" in bad.stderr
        assert json.loads((result / "status.json").read_text())["status"] == "failed"
        # The same validation must reject negative and all-zero cell rows.
        for name, replacement, expected in [
            ("negative", np.array([[-1, 0, 9], *counts[1:]]), "negative"),
            ("zero", np.array([[0, 0, 0], *counts[1:]]), "zero library"),
        ]:
            adata.X = sparse.csr_matrix(replacement)
            invalid_path = root / f"{name}.h5ad"
            adata.write_h5ad(invalid_path)
            invalid_command = command.copy()
            invalid_command[invalid_command.index(str(path))] = str(invalid_path)
            invalid = subprocess.run(invalid_command, capture_output=True, text=True)
            assert invalid.returncode != 0 and expected in invalid.stderr
        unconfirmed = subprocess.run([x for x in command if x != "--full-transcriptome-counts-confirmed"], capture_output=True, text=True)
        assert unconfirmed.returncode != 0 and "target/HVG-only" in unconfirmed.stderr
        # Unknown subclasses cannot produce a within-type comparison.
        adata.X = sparse.csr_matrix(counts)
        adata.obs.loc["cell0", "subclass_label"] = np.nan
        adata.write_h5ad(root / "unknown_subclass.h5ad")
        bad_command[bad_command.index(str(root / "bad_counts.h5ad"))] = str(root / "unknown_subclass.h5ad")
        unknown = subprocess.run(bad_command, capture_output=True, text=True)
        assert unknown.returncode != 0 and "complete subclass" in unknown.stderr
        no_input = subprocess.run([sys.executable, str(script), "--project-root", str(root / "empty")], capture_output=True, text=True)
        assert no_input.returncode == 2 and "pending_external_data" in no_input.stdout
    print("PASS: exact ROI, full-library normalization, detection, unavailable gene, shared-cluster weighting, raw-count guard, taxonomy guard and explicit pending exit")


if __name__ == "__main__":
    main()
