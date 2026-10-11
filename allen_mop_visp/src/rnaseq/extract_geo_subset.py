#!/usr/bin/env python3
"""Extract MOp/VISp candidate counts and full-gene library totals from GEO HDF5.

The source gene-by-cell matrix is dense, gzip-compressed and chunked. Reading
aligned chunks avoids large fancy-index selections and retains all 31,053 genes
in each selected cell's normalization denominator, while saving only candidates.
"""
import argparse
import concurrent.futures
import hashlib
import json
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

DEFAULT_ROOT = Path(__file__).resolve().parents[2]


def read_gene_block(task):
    path, start, stop, selected_rows, target_rows, column_chunk = task
    totals = np.zeros(len(selected_rows), dtype=np.int64)
    requested = [(gene, row - start) for gene, row in target_rows if start <= row < stop]
    chosen = np.zeros((len(selected_rows), len(requested)), dtype=np.int32)
    with h5py.File(path, "r") as source:
        matrix = source["data/counts"]
        occupied_chunks = np.unique(selected_rows // column_chunk)
        for column_number in occupied_chunks:
            column_start = int(column_number * column_chunk)
            column_stop = min(column_start + column_chunk, matrix.shape[1])
            mask = (selected_rows >= column_start) & (selected_rows < column_stop)
            positions = np.flatnonzero(mask)
            local_columns = selected_rows[mask] - column_start
            block = matrix[start:stop, column_start:column_stop]
            if np.any(block < 0):
                raise RuntimeError("Negative entries in a purported raw-count artifact")
            selected_block = block[:, local_columns]
            totals[positions] = selected_block.sum(axis=0, dtype=np.int64)
            for j, (_, local_row) in enumerate(requested):
                chosen[positions, j] = selected_block[local_row]
    return start, stop, totals, [gene for gene, _ in requested], chosen


def build_metadata(root):
    source = root / "input/rnaseq/raw/GSE185862_metadata_10x.csv.gz"
    columns = ["sample_name", "region_label", "platform_label", "class_label", "subclass_label", "cluster_label", "external_donor_name_label", "full_genotype_label", "donor_sex_label", "facs_population_plan_label", "expected_i1_fastq", "cortical_layer_label"]
    subsets, roi_counts = [], []
    for data in pd.read_csv(source, usecols=columns, chunksize=100000, dtype=str, keep_default_na=False):
        roi_counts.append(data.groupby(["region_label", "platform_label", "class_label"], dropna=False).size().rename("n_cells").reset_index())
        subsets.append(data.loc[data.region_label.isin(["MOp", "VISp"])])
    count_table = pd.concat(roi_counts).groupby(["region_label", "platform_label", "class_label"], dropna=False).n_cells.sum().reset_index()
    count_table.to_csv(root / "results/rnaseq/full_metadata_roi_platform_class_counts.csv", index=False)
    selected = pd.concat(subsets, ignore_index=True)
    selected.to_csv(root / "input/rnaseq/selected_original_metadata.csv.gz", index=False)
    original_donors = pd.read_excel(root / "input/rnaseq/raw/GSE185862_CTX-HPF-73538sample_metadata.xlsx", header=3, dtype=str)
    donors = original_donors.loc[original_donors["characteristics: donor_id"].isin(selected.external_donor_name_label)].copy()
    donors.to_csv(root / "input/rnaseq/selected_donor_original_records.csv.gz", index=False)
    age_counts = donors.groupby("characteristics: donor_id")["characteristics: age"].nunique()
    if not age_counts.eq(1).all():
        raise RuntimeError("Donor ages are missing or inconsistent in the original donor records")
    ages = donors.groupby("characteristics: donor_id")["characteristics: age"].first()
    for sample, region in [("GSM5629156", "MOp"), ("GSM5629158", "VISp")]:
        record = (root / f"input/rnaseq/{sample}.soft.txt").read_text()
        if "Chromium Single Cell 3’ Reagent Kit v2" not in record or "!Sample_title = Mouse 10x " + region not in record:
            raise RuntimeError("The required region-specific official GEO chemistry record is not verified")
    mapping = {
        "sample_name": "cell_id", "region_label": "region", "class_label": "class",
        "subclass_label": "subclass", "cluster_label": "cluster",
        "external_donor_name_label": "donor", "full_genotype_label": "genotype",
        "donor_sex_label": "sex", "facs_population_plan_label": "facs",
    }
    result = selected.rename(columns=mapping)
    result["chemistry"] = "10x_v2"
    result["age"] = result.donor.map(ages)
    result["library_preparation_id"] = result.cell_id.str.split("-", n=1).str[1]
    result["library_preparation_date"] = result.library_preparation_id.str.extract(r"L8TX_([0-9]{6})_")[0]
    result["batch"] = result.library_preparation_date
    # Original protocol states same-ROI samples were pooled before loading.
    # Donor labels are retained, but do not invent an independent pool ID.
    result["pool"] = "unknown"
    result["genetic_background"] = "C57BL/6J; transgenic reporter background, not WT"
    result["age_source"] = "GSE185862_CTX-HPF-73538sample_metadata.xlsx donor_id→unique age"
    result["chemistry_source"] = "GSM5629156/GSM5629158 official SOFT extract protocol"
    result.groupby(["region", "donor", "age", "batch"], dropna=False).size().rename("n_cells").reset_index().to_csv(root / "results/rnaseq/observed_donor_age_batch_region_counts.csv", index=False)
    return result, count_table


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args()
    if not 1 <= args.workers <= 4:
        parser.error("Use 1–4 workers to respect the cloud machine's four CPU quota")
    root = args.project_root.resolve()
    (root / "results/rnaseq").mkdir(parents=True, exist_ok=True)
    metadata, counts = build_metadata(root)
    matrix_path = root / "input/rnaseq/raw/GSE185862_expression_matrix_10x.hdf5"
    panel = pd.read_csv(root / "results/pbm/panel22.csv").gene.astype(str)
    common = pd.read_csv(root / "results/pbm/same_direction_52.csv").gene.astype(str)
    candidates = sorted(set(panel) | set(common))
    with h5py.File(matrix_path, "r") as source:
        matrix = source["data/counts"]
        shape, chunks, dtype = matrix.shape, matrix.chunks, str(matrix.dtype)
        symbols = pd.Index([x.decode() for x in source["data/gene"][:]])
        cell_ids = pd.Index([x.decode() for x in source["data/samples"][:]])
        declared_shape = source["data/shape"][:].tolist()
    # The released counts dataset is gene×cell, but data/shape declares
    # the reverse cell×gene orientation. Identify axes using exact identity
    # lengths instead of blindly reshaping counts according to this field.
    if declared_shape not in [list(shape), list(shape[::-1])] or len(symbols) != shape[0] or len(cell_ids) != shape[1]:
        raise RuntimeError("GEO matrix identity dimensions are inconsistent")
    if cell_ids.duplicated().any() or symbols.duplicated().any():
        raise RuntimeError("GEO matrix cell IDs or gene symbols are duplicated")
    row_lookup = pd.Series(np.arange(len(cell_ids)), index=cell_ids)
    metadata["source_matrix_column"] = metadata.cell_id.map(row_lookup)
    if metadata.source_matrix_column.isna().any() or metadata.cell_id.duplicated().any():
        raise RuntimeError("Selected cells fail a unique exact match to original matrix columns")
    metadata = metadata.sort_values("source_matrix_column").reset_index(drop=True)
    positions = metadata.source_matrix_column.to_numpy(dtype=np.int64)
    available = [gene for gene in candidates if gene in symbols]
    target_rows = [(gene, int(symbols.get_loc(gene))) for gene in available]
    target_lookup = {gene: j for j, gene in enumerate(available)}
    extracted = np.zeros((len(metadata), len(available)), dtype=np.int32)
    totals = np.zeros(len(metadata), dtype=np.int64)
    tasks = [(str(matrix_path), start, min(start + chunks[0], shape[0]), positions, target_rows, chunks[1]) for start in range(0, shape[0], chunks[0])]
    start_time = time.monotonic()
    completed = 0
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = [executor.submit(read_gene_block, task) for task in tasks]
        for future in concurrent.futures.as_completed(futures):
            start, stop, contribution, genes, values = future.result()
            totals += contribution
            for j, gene in enumerate(genes):
                extracted[:, target_lookup[gene]] = values[:, j]
            completed += 1
            print(json.dumps({"gene_blocks_complete": completed, "gene_blocks_total": len(tasks), "source_gene_range": [start, stop], "elapsed_seconds": round(time.monotonic() - start_time, 1)}), flush=True)
    if np.any(totals <= 0) or np.any(extracted.sum(axis=1, dtype=np.int64) > totals):
        raise RuntimeError("Full-gene library denominators are invalid")
    metadata["full_library_counts"] = totals
    metadata.to_csv(root / "input/rnaseq/selected_cell_metadata.csv.gz", index=False)
    obs = metadata.set_index("cell_id")
    artifact = ad.AnnData(sparse.csr_matrix(extracted), obs=obs, var=pd.DataFrame({"gene_symbol": available}, index=available))
    artifact.uns["source"] = "GEO GSE185862 10x_v2"
    artifact.uns["matrix_scope"] = "62-candidate subset; full-library counts separately sum all 31053 source genes"
    artifact.uns["denominator"] = "obs.full_library_counts, exact source raw-count sum over all source genes"
    subset_path = root / "input/rnaseq/ctx_hpf_mop_visp_candidate_counts.h5ad"
    artifact.write_h5ad(subset_path, compression="gzip")
    manifest = {"status": "subset_extracted", "source_matrix": str(matrix_path.relative_to(root)), "source_matrix_shape_gene_cell": list(shape), "source_declared_data_shape": declared_shape, "orientation_resolution": "counts axes identified by lengths of data/gene and data/samples; data/shape is recorded separately", "source_matrix_chunks": list(chunks), "source_matrix_dtype": dtype, "n_source_metadata_cells": int(counts.n_cells.sum()), "n_selected_cells": len(metadata), "n_source_genes_in_library_denominator": shape[0], "n_candidate_genes_requested": len(candidates), "n_candidate_genes_available": len(available), "normalization_denominator": "Exact raw-count sum over all source genes; not candidate/HVG sum", "subset_file": str(subset_path.relative_to(root)), "subset_bytes": subset_path.stat().st_size, "subset_sha256": hashlib.file_digest(subset_path.open("rb"), "sha256").hexdigest(), "selected_metadata_sha256": hashlib.file_digest((root / "input/rnaseq/selected_cell_metadata.csv.gz").open("rb"), "sha256").hexdigest(), "source_matrix_cell_names_absent_from_metadata": shape[1] - int(counts.n_cells.sum()), "workers": args.workers, "pool_independence": "not established; official same-ROI pooling protocol retained; no animal inference", "region_library_preparation_confounding": "No shared library-preparation date between MOp and VISp; cannot isolate a region effect from date"}
    (root / "input/rnaseq/subset_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n")
    print(json.dumps(manifest, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
