#!/usr/bin/env python3
"""Descriptive MOp/VISp analysis of an explicitly identified CTX-HPF H5AD.

This script never treats cells as independent animals and does not calculate
animal-level P values. A local subset is retained; the complete matrix stays at
the input location. The matrix value scale must be supplied explicitly.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse

DEFAULT_ROOT = Path(__file__).resolve().parents[2]
ALIASES = {
    "cell_id": ["cell_id", "cell_label", "sample_name", "exp_component_name"],
    "region": ["region", "region_label", "structure_acronym", "brain_region", "roi"],
    "class": ["class", "class_label", "cell_class", "class_name"],
    "subclass": ["subclass", "subclass_label", "cell_subclass", "subclass_name"],
    "cluster": ["cluster", "cluster_label", "cell_cluster", "cluster_name"],
    "donor": ["donor", "donor_id", "donor_name", "external_donor_name"],
    "pool": ["pool", "pool_id", "biological_pool", "pool_name"],
    "chemistry": ["chemistry", "chemistry_label", "library_prep", "platform", "assay"],
    "genotype": ["genotype", "genotype_label", "donor_genotype", "full_genotype"],
    "sex": ["sex", "donor_sex", "donor_sex_label", "sex_label"],
    "age": ["age", "donor_age", "age_days", "age_label"],
    "facs": ["facs", "facs_gate", "facs_population", "selection", "sort_gate"],
    "batch": ["batch", "batch_id", "library_id", "experiment_id"],
    "library_total": ["library_total", "total_umis", "nCount_RNA", "total_counts"],
}
STRATA = ["chemistry", "genotype", "sex", "age", "facs"]


def write_json(path: Path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def load_candidates(root):
    files = [root / "results/pbm/panel22.csv", root / "results/pbm/same_direction_52.csv"]
    missing = [str(p) for p in files if not p.exists()]
    if missing:
        raise ValueError("Candidate tables unavailable: " + ", ".join(missing))
    panel = pd.read_csv(files[0])["gene"].dropna().astype(str).tolist()
    common = pd.read_csv(files[1])["gene"].dropna().astype(str).tolist()
    return sorted(set(panel + common)), set(panel), set(common)


def standardize_metadata(df, overrides=None):
    overrides = overrides or {}
    fields = {}
    source = df.copy()
    if "cell_id" not in source and source.index.name not in (None, "index"):
        source["cell_id"] = source.index.astype(str)
    for name, aliases in ALIASES.items():
        if name in overrides:
            column = overrides[name]
            if column not in source:
                raise ValueError(f"Explicit schema field {name!r} refers to absent {column!r}")
        else:
            available = [a for a in aliases if a in source]
            column = available[0] if available else None
            if name == "region" and len(available) > 1:
                # A structure column may include layers; prefer the explicit
                # column that actually contains both exact requested regions.
                candidates = [a for a in available if {"MOp", "VISp"} <= set(source[a].astype(str))]
                if candidates:
                    column = candidates[0]
        fields[name] = {
            "source_column": column,
            "present": column is not None,
            "n_missing": int(source[column].isna().sum()) if column else len(source),
        }
    if not fields["region"]["present"]:
        raise ValueError("No region field found; use --schema-map to identify an exact MOp/VISp column")
    out = pd.DataFrame(index=source.index)
    for name, info in fields.items():
        out[name] = source[info["source_column"]].astype("string").fillna("unknown").astype(str).str.strip().replace("", "unknown") if info["present"] else "unknown"
        missing_tokens = out[name].str.casefold().isin(["unknown", "unassigned", "none", "nan", "na", "n/a"])
        out.loc[missing_tokens, name] = "unknown"
        if info["present"]:
            info["n_missing"] = int(out[name].eq("unknown").sum())
    if not fields["cell_id"]["present"]:
        out["cell_id"] = source.index.astype(str)
        fields["cell_id"] = {"source_column": "index", "present": True, "n_missing": 0}
    if out["cell_id"].duplicated().any():
        raise ValueError("Cell IDs are duplicated; identify the globally unique cell ID column")
    return out, fields


def audit_metadata(metadata, fields, result_dir):
    result_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame([{"field": name, **info} for name, info in fields.items()]).to_csv(result_dir / "metadata_fields.csv", index=False)
    metadata.groupby("region", dropna=False).size().rename("n_cells").reset_index().to_csv(result_dir / "all_roi_counts.csv", index=False)
    selected = metadata.loc[metadata.region.isin(["MOp", "VISp"])].copy()
    groups = ["region", "class", "subclass"] + STRATA
    selected.groupby(groups, dropna=False).size().rename("n_cells").reset_index().to_csv(result_dir / "selected_cell_counts.csv", index=False)
    selected.to_csv(result_dir / "selected_cell_metadata.csv.gz", index=False)
    for unit in ["donor", "pool", "batch"]:
        selected.groupby([unit, "region"] + STRATA, dropna=False).size().rename("n_cells").reset_index().to_csv(result_dir / f"{unit}_region_counts.csv", index=False)
    regions = sorted(selected.region.unique().tolist())
    unknown = [f for f in ["class", "subclass", "cluster", "chemistry", "genotype", "sex", "age", "facs", "donor", "pool", "batch"] if not fields[f]["present"] or selected[f].eq("unknown").any()]
    status = {
        "scope": "CTX-HPF MOp and VISp; RNA-seq cell counts are sampling counts, not tissue composition",
        "n_all_cells": len(metadata), "n_selected_cells": len(selected),
        "exact_regions_present": regions, "unknown_or_partial_fields": unknown,
        "genotype_policy": "Retain raw annotation; do not infer WT or C57BL/6J from missing values",
        "animal_inference": "withheld; donor/pool independence and region/batch confounding require documented review",
    }
    write_json(result_dir / "metadata_audit.json", status)
    return selected, status


def get_gene_symbols(adata, mapping_path=None):
    if mapping_path:
        mapping = pd.read_csv(mapping_path)
        if not {"gene_id", "gene_symbol"} <= set(mapping):
            raise ValueError("--gene-map must have gene_id and gene_symbol columns")
        if mapping.gene_id.duplicated().any():
            raise ValueError("--gene-map has duplicate gene_id entries")
        lookup = mapping.set_index("gene_id").gene_symbol
        return pd.Index(adata.var_names.astype(str)).map(lookup)
    for name in ["gene_symbol", "gene_name", "feature_name", "symbol"]:
        if name in adata.var:
            return pd.Index(adata.var[name].astype(str))
    return pd.Index(adata.var_names.astype(str))


def read_subset(adata, rows, target_indexes, value_kind, chunk_size, external_library_totals=False):
    out = np.zeros((len(rows), len(target_indexes)), dtype=np.float64)
    libraries = np.zeros(len(rows), dtype=np.float64)
    for start in range(0, len(rows), chunk_size):
        stop = min(start + chunk_size, len(rows))
        # Row access can require all genes for sparse H5AD; chunks cap memory.
        block = adata.X[rows[start:stop], :]
        if sparse.issparse(block):
            values = block.data
            chosen = block[:, target_indexes].toarray()
            totals = np.asarray(block.sum(axis=1)).ravel()
        else:
            block = np.asarray(block)
            values = block.ravel()
            chosen = block[:, target_indexes]
            totals = block.sum(axis=1)
        if not np.isfinite(values).all() or np.any(values < 0):
            raise ValueError("Expression matrix has negative or non-finite entries")
        if value_kind == "raw_counts":
            if np.any(np.abs(values - np.rint(values)) > 1e-6):
                raise ValueError("raw_counts was specified but non-integer values were found")
            if not external_library_totals and np.any(totals <= 0):
                raise ValueError("Selected cells have zero library size; do not normalize these cells silently")
        out[start:stop] = chosen
        libraries[start:stop] = totals
    return out, libraries


def summarize_expression(selected, values, genes, value_kind, libraries, panel, common, result_dir):
    normalized = np.log1p(values / libraries[:, None] * 1e4) if value_kind == "raw_counts" else values
    normalization = "log1p(CP10k), confirmed full-transcriptome library denominator" if value_kind == "raw_counts" else value_kind
    keys = ["region", "class", "subclass"] + STRATA
    records = []
    # Work on the positional index because values use the same cell ordering.
    positional = selected.reset_index(drop=True)
    for labels, indexes in positional.groupby(keys, dropna=False).indices.items():
        indexes = np.asarray(indexes)
        block = normalized[indexes]
        detection = values[indexes] > 0
        for j, gene in enumerate(genes):
            records.append(dict(zip(keys, labels), gene=gene, n_cells=len(indexes),
                mean_expression=float(block[:, j].mean()), detection_fraction=float(detection[:, j].mean()),
                in_panel22=gene in panel, in_common52=gene in common, value_kind=normalization))
    summary = pd.DataFrame(records)
    summary.to_csv(result_dir / "region_subclass_gene_summary.csv", index=False)
    overall = []
    for labels, indexes in positional.groupby(["region"] + STRATA, dropna=False).indices.items():
        indexes = np.asarray(indexes)
        for j, gene in enumerate(genes):
            overall.append(dict(zip(["region"] + STRATA, labels), gene=gene, n_cells=len(indexes),
                mean_expression=float(normalized[indexes, j].mean()), detection_fraction=float((values[indexes, j] > 0).mean()),
                value_kind=normalization, interpretation="FACS/RNA-seq sampled-cell average; not an estimate of tissue-average expression"))
    pd.DataFrame(overall).to_csv(result_dir / "region_all_selected_cells_gene_summary.csv", index=False)
    comparison_keys = ["class", "subclass"] + STRATA + ["gene"]
    mop = summary.loc[summary.region.eq("MOp")].drop(columns="region")
    visp = summary.loc[summary.region.eq("VISp")].drop(columns="region")
    comparison = mop.merge(visp, on=comparison_keys, how="outer", suffixes=("_MOp", "_VISp"), indicator=True)
    comparison["both_regions_observed"] = comparison["_merge"].eq("both")
    comparison["MOp_minus_VISp_mean_expression"] = comparison.mean_expression_MOp - comparison.mean_expression_VISp
    comparison["MOp_minus_VISp_detection_fraction"] = comparison.detection_fraction_MOp - comparison.detection_fraction_VISp
    comparison["analysis_level"] = "descriptive cell-weighted summaries; no animal-level significance test"
    comparison.drop(columns="_merge").to_csv(result_dir / "within_subclass_region_difference.csv", index=False)

    # Cluster means and shared-cluster standardization separate some subtype
    # mixture effects from the gross subclass comparison. No animal test.
    cluster_known = positional.loc[~positional.cluster.eq("unknown")]
    write_json(result_dir / "shared_cluster_audit.json", {
        "n_cells_without_cluster_excluded_from_standardization": int(positional.cluster.eq("unknown").sum()),
        "n_cells_with_cluster": len(cluster_known),
        "interpretation": "Only observed, labeled clusters shared across both exact regions enter equal-cluster weighting",
    })
    if len(cluster_known):
        cluster_keys = keys + ["cluster"]
        cluster_records = []
        for labels, cluster_group in cluster_known.groupby(cluster_keys, dropna=False):
            indexes = cluster_group.index.to_numpy(dtype=int)
            means = normalized[np.asarray(indexes)].mean(axis=0)
            for j, gene in enumerate(genes):
                cluster_records.append(dict(zip(cluster_keys, labels), gene=gene, n_cells=len(indexes), mean_expression=float(means[j])))
        clusters = pd.DataFrame(cluster_records)
        clusters.to_csv(result_dir / "region_cluster_gene_summary.csv", index=False)
        match_keys = ["class", "subclass"] + STRATA + ["cluster", "gene"]
        shared = clusters.loc[clusters.region.eq("MOp")].merge(clusters.loc[clusters.region.eq("VISp")], on=match_keys, suffixes=("_MOp", "_VISp"))
        shared["MOp_minus_VISp"] = shared.mean_expression_MOp - shared.mean_expression_VISp
        shared.to_csv(result_dir / "shared_cluster_region_difference.csv", index=False)
        if len(shared):
            standardized = shared.groupby(comparison_keys, dropna=False).agg(
                n_shared_clusters=("cluster", "nunique"),
                equal_cluster_weighted_MOp_mean=("mean_expression_MOp", "mean"),
                equal_cluster_weighted_VISp_mean=("mean_expression_VISp", "mean"),
                equal_cluster_weighted_MOp_minus_VISp=("MOp_minus_VISp", "mean"),
            ).reset_index()
            standardized["weighting"] = "equal weight for each matched cluster; exploratory, not tissue composition"
            standardized.to_csv(result_dir / "shared_cluster_standardized_difference.csv", index=False)
    np.savez_compressed(result_dir / "selected_candidate_expression.npz", expression=values, cell_id=selected.cell_id.to_numpy(dtype=str), gene_symbol=np.array(genes), library_total=libraries, value_kind=np.array(value_kind))
    return summary, comparison


def dotplots(summary, root):
    os.environ.setdefault("MPLCONFIGDIR", str(root / ".cache/matplotlib"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figure_dir = root / "figures/rnaseq"
    figure_dir.mkdir(parents=True, exist_ok=True)
    panel = summary.loc[summary.in_panel22].copy()
    # Separate each chemistry/genotype/sex/age/FACS stratum instead of silently
    # mixing biological or assay conditions in a single plot.
    index = []
    for number, (labels, data) in enumerate(panel.groupby(STRATA, dropna=False)):
        genes = sorted(data.gene.unique())
        rows = sorted(set(zip(data.region, data["class"], data.subclass)))
        fig, ax = plt.subplots(figsize=(max(9, len(genes) * .44), max(4, len(rows) * .34)))
        coords = {row: i for i, row in enumerate(rows)}
        x = data.gene.map({gene: i for i, gene in enumerate(genes)})
        y = [coords[(r, c, s)] for r, c, s in zip(data.region, data["class"], data.subclass)]
        points = ax.scatter(x, y, s=5 + data.detection_fraction.to_numpy() * 170, c=data.mean_expression, cmap="viridis")
        ax.set_xticks(range(len(genes)), genes, rotation=70, ha="right")
        ax.set_yticks(range(len(rows)), [f"{r} | {c} | {s}" for r, c, s in rows])
        ax.set_title("MOp / VISp panel22: dot size = detection fraction\n" + "; ".join(f"{k}={v}" for k, v in zip(STRATA, labels)), fontsize=9)
        fig.colorbar(points, ax=ax, label=data.value_kind.iloc[0])
        for fraction in [.25, .5, 1.0]:
            ax.scatter([], [], s=5 + fraction * 170, c="gray", label=f"{fraction:.0%} detected")
        ax.legend(loc="upper left", bbox_to_anchor=(1.20, 1), frameon=False, fontsize=8)
        fig.tight_layout()
        filename = f"panel22_stratum_{number:03d}.png"
        fig.savefig(figure_dir / filename, dpi=160)
        plt.close(fig)
        index.append({"filename": filename, **dict(zip(STRATA, labels))})
    pd.DataFrame(index).to_csv(figure_dir / "dotplot_strata.csv", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--h5ad", type=Path)
    parser.add_argument("--schema-map", type=Path, help="JSON canonical field → actual metadata column")
    parser.add_argument("--gene-map", type=Path, help="CSV with gene_id,gene_symbol")
    parser.add_argument("--value-kind", choices=["raw_counts", "log1p_cp10k", "log2_cpm_plus1", "linear_cpm"])
    parser.add_argument("--full-transcriptome-counts-confirmed", action="store_true", help="Explicitly attest the raw-count matrix contains the full normalization gene universe")
    parser.add_argument("--library-total-column", help="Metadata column with externally computed full-transcriptome raw-count totals; alternative to full-matrix confirmation")
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--chunk-size", type=int, default=2048)
    args = parser.parse_args()
    root = args.project_root.resolve()
    result_dir = root / "results/rnaseq"
    result_dir.mkdir(parents=True, exist_ok=True)
    write_json(result_dir / "status.json", {"status": "running_incomplete", "executed": False,
        "note": "Only a descriptive_expression_complete status validates this run; older result files may remain"})
    if args.chunk_size <= 0:
        parser.error("--chunk-size must be positive")
    if not args.metadata and not args.h5ad:
        status = {"status": "pending_external_data", "executed": False,
                  "missing": ["CTX-HPF MOp/VISp cell metadata", "cell-level expression matrix and documented scale"],
                  "not_generated": ["RNA-seq region expression differences", "animal-level P values"],
                  "next_command": "python src/rnaseq/analyze_expression.py --metadata input/rnaseq/cell_metadata.csv --h5ad input/rnaseq/ctx_hpf.h5ad --value-kind raw_counts --full-transcriptome-counts-confirmed"}
        write_json(result_dir / "status.json", status)
        print(json.dumps(status, ensure_ascii=False))
        raise SystemExit(2)
    overrides = json.loads(args.schema_map.read_text(encoding="utf-8")) if args.schema_map else {}
    if args.library_total_column:
        overrides["library_total"] = args.library_total_column
    adata = None
    if args.h5ad:
        import anndata as ad
        adata = ad.read_h5ad(args.h5ad, backed="r")
    raw_metadata = pd.read_csv(args.metadata, low_memory=False) if args.metadata else adata.obs.copy()
    metadata, fields = standardize_metadata(raw_metadata, overrides)
    selected, audit = audit_metadata(metadata, fields, result_dir)
    if args.audit_only or not args.h5ad:
        write_json(result_dir / "status.json", {"status": "metadata_audit_complete_expression_pending", "metadata": audit})
        print("Metadata audit saved; expression analysis awaits an explicitly documented matrix")
        if adata is not None:
            adata.file.close()
        return
    if not args.value_kind:
        parser.error("--value-kind is required; do not infer normalization from matrix values")
    if args.value_kind == "raw_counts" and not (args.full_transcriptome_counts_confirmed or args.library_total_column):
        parser.error("Raw counts require --full-transcriptome-counts-confirmed or a documented --library-total-column; target/HVG-only sums are not full library sizes")
    if {"MOp", "VISp"} != set(selected.region):
        raise ValueError("Both exact MOp and VISp regions must be present; VIS/VIS-PTLp/ALM cannot substitute")
    for field in ["class", "subclass"]:
        if not fields[field]["present"] or selected[field].eq("unknown").any():
            raise ValueError(f"Selected cells lack complete {field} annotation; provide the correct --schema-map before within-type analysis")
    genes, panel, common = load_candidates(root)
    symbols = get_gene_symbols(adata, args.gene_map)
    target = [g for g in genes if g in symbols]
    if not target:
        raise ValueError("No candidate gene symbols matched; supply a validated gene-ID mapping")
    duplicates = [g for g in target if (symbols == g).sum() != 1]
    if duplicates:
        raise ValueError("Candidate symbols map to multiple matrix features: " + ", ".join(duplicates))
    lookup = pd.Series(np.arange(adata.n_obs), index=adata.obs_names.astype(str))
    if lookup.index.duplicated().any():
        raise ValueError("H5AD obs_names are duplicated")
    absent = selected.loc[~selected.cell_id.isin(lookup.index), "cell_id"]
    if len(absent):
        raise ValueError(f"{len(absent)} selected metadata cells are missing from matrix; do not silently drop them")
    selected["matrix_row"] = selected.cell_id.map(lookup)
    selected = selected.sort_values("matrix_row").reset_index(drop=True)
    target_indexes = [int(np.flatnonzero(symbols == g)[0]) for g in target]
    values, libraries = read_subset(adata, selected.matrix_row.to_numpy(dtype=int), target_indexes, args.value_kind, args.chunk_size, external_library_totals=bool(args.library_total_column))
    adata.file.close()
    if args.value_kind == "raw_counts" and args.library_total_column:
        supplied_totals = pd.to_numeric(selected.library_total, errors="coerce").to_numpy(dtype=float)
        if not np.isfinite(supplied_totals).all() or np.any(supplied_totals <= 0) or np.any(supplied_totals + 1e-6 < libraries) or np.any(np.abs(supplied_totals - np.rint(supplied_totals)) > 1e-6):
            raise ValueError("Full-transcriptome library totals must be finite, positive integers, and at least the counts in the supplied matrix")
        libraries = supplied_totals
    elif args.value_kind != "raw_counts":
        libraries[:] = np.nan  # A sum of normalized values is not library size.
    pd.DataFrame({"gene": genes, "present_in_matrix": [g in target for g in genes], "in_panel22": [g in panel for g in genes], "in_common52": [g in common for g in genes], "absence_interpretation": ["measured" if g in target else "unavailable; not zero expression" for g in genes]}).to_csv(result_dir / "candidate_gene_coverage.csv", index=False)
    summary, comparison = summarize_expression(selected, values, target, args.value_kind, libraries, panel, common, result_dir)
    dotplots(summary, root)
    status = {"status": "descriptive_expression_complete", "n_selected_cells": len(selected), "n_candidates_present": len(target), "n_candidates_total": len(genes), "matrix": str(args.h5ad), "value_kind": args.value_kind, "metadata_audit": audit,
              "animal_statistics": "not performed; independent biological-unit mapping and design review required",
              "normalization_denominator_evidence": ("external metadata column " + args.library_total_column) if args.library_total_column else ("user explicitly confirmed full-transcriptome raw-count matrix" if args.value_kind == "raw_counts" else "already-normalized input; no library size inferred"),
              "comparison_interpretation": "descriptive matched annotation labels; matching chemistry/genotype/age/sex/FACS is unverified where unknown" if audit["unknown_or_partial_fields"] else "descriptive matched recorded strata; biological independence/design review still required",
              "reported_difference": "MOp minus VISp in the supplied normalized expression scale; not log2FC unless separately modeled"}
    write_json(result_dir / "status.json", status)
    print(json.dumps(status, ensure_ascii=False))


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        root = DEFAULT_ROOT
        if "--project-root" in sys.argv:
            position = sys.argv.index("--project-root")
            if position + 1 < len(sys.argv):
                root = Path(sys.argv[position + 1]).resolve()
        write_json(root / "results/rnaseq/status.json", {
            "status": "failed", "error": str(error), "executed": False,
            "note": "Older result files may remain; do not interpret them as outputs from this failed run",
        })
        raise
