#!/usr/bin/env python3
"""Whole-transcriptome donor x subclass pseudobulk from the GEO GSE185862 10x matrix.

Reads the dense, gzip-chunked gene x cell matrix block by block (same access pattern as
extract_geo_subset.py) for the MOp and VISp cells already listed in
input/rnaseq/selected_cell_metadata.csv.gz, and writes for every gene and every
region x donor x subclass group:
- the raw-count sum
- the number of cells with a nonzero count (detection)
plus per-cell raw counts of the prespecified focus genes (focus_gene_sets.py) for
dotplots and module scores. Library-size denominators are the full_library_counts already
stored in the metadata (sum over all 31,053 source genes).

Outputs (results/rnaseq/whole_transcriptome/):
- pseudobulk_counts.npz   sums, detection, genes, group table
- focus_gene_sets.csv     gene membership and presence
- focus_cells.h5ad        per-cell counts of focus genes (sparse, gzip)
"""
from __future__ import annotations

import argparse
import concurrent.futures
import json
import sys
import time
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import focus_gene_sets as fgs  # noqa: E402

MATRIX = ROOT / "input/rnaseq/raw/GSE185862_expression_matrix_10x.hdf5"
META = ROOT / "input/rnaseq/selected_cell_metadata.csv.gz"
OUT = ROOT / "results/rnaseq/whole_transcriptome"
CANDIDATES = ROOT / "results/pbm/exploration_genes62.csv"


def read_block(task):
    path, start, stop, columns, group_index, n_groups, focus_rows, column_chunk = task
    n_genes = stop - start
    sums = np.zeros((n_genes, n_groups), dtype=np.float64)
    detected = np.zeros((n_genes, n_groups), dtype=np.int64)
    local_focus = [(j, row - start) for j, row in focus_rows if start <= row < stop]
    focus = np.zeros((len(columns), len(local_focus)), dtype=np.int32)
    with h5py.File(path, "r") as source:
        matrix = source["data/counts"]
        for chunk_number in np.unique(columns // column_chunk):
            c0 = int(chunk_number * column_chunk)
            c1 = min(c0 + column_chunk, matrix.shape[1])
            mask = (columns >= c0) & (columns < c1)
            positions = np.flatnonzero(mask)
            block = matrix[start:stop, c0:c1][:, columns[mask] - c0]
            indicator = sparse.csr_matrix(
                (np.ones(len(positions)), (np.arange(len(positions)), group_index[positions])),
                shape=(len(positions), n_groups))
            sums += sparse.csr_matrix(block.astype(np.float64)) @ indicator
            detected += (sparse.csr_matrix((block > 0).astype(np.int64)) @ indicator).astype(np.int64)
            for k, (_, local_row) in enumerate(local_focus):
                focus[positions, k] = block[local_row]
    return start, stop, sums, detected, [j for j, _ in local_focus], focus


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workers", type=int, default=2)
    args = parser.parse_args(argv)
    OUT.mkdir(parents=True, exist_ok=True)

    meta = pd.read_csv(META, low_memory=False)
    meta = meta.sort_values("source_matrix_column").reset_index(drop=True)
    meta["group"] = meta.region + "|" + meta.donor.astype(str) + "|" + meta.subclass
    groups = pd.Index(sorted(meta.group.unique()))
    group_index = groups.get_indexer(meta.group).astype(np.int64)
    columns = meta.source_matrix_column.to_numpy(np.int64)

    with h5py.File(MATRIX, "r") as source:
        matrix = source["data/counts"]
        shape, chunks = matrix.shape, matrix.chunks
        genes = pd.Index([x.decode() for x in source["data/gene"][:]])
        cells = pd.Index([x.decode() for x in source["data/samples"][:]])
    if len(genes) != shape[0] or len(cells) != shape[1]:
        raise RuntimeError("Matrix orientation is not gene x cell as expected")
    if not (cells[columns] == pd.Index(meta.cell_id)).all():
        raise RuntimeError("Metadata cell IDs do not match matrix columns")

    pbm = pd.read_csv(CANDIDATES).iloc[:, 0].astype(str).tolist()
    sets = fgs.build(genes, pbm)
    sets.to_csv(OUT / "focus_gene_sets.csv", index=False)
    focus_genes = sorted(sets[sets.present].gene.unique())
    focus_rows = [(j, int(genes.get_loc(g))) for j, g in enumerate(focus_genes)]

    sums = np.zeros((shape[0], len(groups)), dtype=np.float64)
    detected = np.zeros((shape[0], len(groups)), dtype=np.int64)
    focus = np.zeros((len(meta), len(focus_genes)), dtype=np.int32)
    tasks = [(str(MATRIX), s, min(s + chunks[0], shape[0]), columns, group_index, len(groups), focus_rows, chunks[1])
             for s in range(0, shape[0], chunks[0])]
    t0 = time.monotonic()
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as pool:
        for k, future in enumerate(concurrent.futures.as_completed([pool.submit(read_block, t) for t in tasks]), 1):
            start, stop, s, d, idx, f = future.result()
            sums[start:stop] = s
            detected[start:stop] = d
            for c, j in enumerate(idx):
                focus[:, j] = f[:, c]
            print(json.dumps({"blocks_done": k, "blocks": len(tasks), "seconds": round(time.monotonic() - t0)}),
                  flush=True)

    library = meta.groupby("group").full_library_counts.sum().reindex(groups).to_numpy()
    if not np.allclose(sums.sum(axis=0), library):
        raise RuntimeError("Group count sums do not equal the stored full-library totals")
    table = meta.groupby("group").agg(region=("region", "first"), donor=("donor", "first"),
                                      subclass=("subclass", "first"), cell_class=("class", "first"),
                                      n_cells=("cell_id", "size"),
                                      library_counts=("full_library_counts", "sum")).reindex(groups)
    table.index.name = "group"
    table.reset_index().to_csv(OUT / "pseudobulk_groups.csv", index=False)
    np.savez_compressed(OUT / "pseudobulk_counts.npz", sums=sums.astype(np.float32), detected=detected.astype(np.int32),
                        genes=np.array(genes), groups=np.array(groups))
    obs = meta.set_index("cell_id")[["region", "donor", "subclass", "class", "cluster", "full_library_counts"]]
    focus_ad = ad.AnnData(sparse.csr_matrix(focus), obs=obs, var=pd.DataFrame(index=focus_genes))
    focus_ad.write_h5ad(OUT / "focus_cells.h5ad", compression="gzip")
    report = {"source": "GEO GSE185862 10x_v2", "cells": int(len(meta)), "genes": int(shape[0]),
              "groups": int(len(groups)), "focus_genes": len(focus_genes),
              "check": "group sums equal stored full-library totals", "seconds": round(time.monotonic() - t0)}
    (OUT / "pseudobulk_report.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
