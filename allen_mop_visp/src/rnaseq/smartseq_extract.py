#!/usr/bin/env python3
"""Extract exact MOp/VISp SMART-seq candidates and validate full exon libraries.

The GEO source stores both sparse orientations. `exon` is CSC cell×gene and
`t_exon` is CSC gene×cell; floating storage contains integer read counts.
"""
import hashlib,json,re,time
from pathlib import Path
import anndata as ad
import h5py
import numpy as np
import pandas as pd
from scipy import sparse
ROOT=Path(__file__).resolve().parents[2]

def main():
 out=ROOT/'input/rnaseq/smartseq';result=ROOT/'results/rnaseq/smartseq';out.mkdir(parents=True,exist_ok=True);result.mkdir(parents=True,exist_ok=True)
 original=pd.read_csv(ROOT/'input/rnaseq/raw/GSE185862_metadata_ssv4.csv.gz',dtype=str,keep_default_na=False)
 original.groupby(['region_label','class_label','platform_label']).size().rename('n_cells').reset_index().to_csv(result/'all_metadata_region_class_counts.csv',index=False)
 selected=original.loc[original.region_label.isin(['MOp','VISp'])].copy()
 original_samples=pd.read_excel(ROOT/'input/rnaseq/raw/GSE185862_CTX-HPF-73538sample_metadata.xlsx',header=3,dtype=str).fillna('')
 donors=original_samples.loc[original_samples['characteristics: donor_id'].isin(selected.external_donor_name_label)]
 agegroups=donors.groupby('characteristics: donor_id')['characteristics: age'].agg(lambda x:sorted(set(x)))
 if any(len(v)!=1 for v in agegroups): raise RuntimeError('Contradictory donor ages in original GEO sample records')
 ages=agegroups.map(lambda v:v[0])
 donors.to_csv(out/'selected_original_donor_records.csv.gz',index=False)
 ages.rename('age_original').to_csv(out/'donor_age_map.csv')
 mapping={'sample_name':'cell_id','region_label':'region','class_label':'class','subclass_label':'subclass','cluster_label':'cluster','external_donor_name_label':'donor','full_genotype_label':'genotype','donor_sex_label':'sex','facs_population_plan_label':'facs'}
 metadata=selected.rename(columns=mapping)
 metadata['age_original']=metadata.donor.map(ages).fillna('unknown')
 metadata['age']=metadata.age_original.map(lambda x:'P'+re.search(r'\d+',x)[0] if re.fullmatch(r'(P\d+|\d+ days)',x) else 'unknown')
 sample_batch=original_samples.drop_duplicates('Sample name').set_index('Sample name')['characteristics: sequencing_batch']
 metadata['batch']=metadata.cell_id.map(sample_batch).fillna('unknown').replace('','unknown')
 for col in ['genotype','sex','facs','class','subclass','cluster','donor']: metadata[col]=metadata[col].replace('','unknown')
 metadata['chemistry']='SMART_seq_v4';metadata['pool']='single_cell_well; donor independence not assessed';metadata['genetic_background']='C57BL/6J background; reporter/transgenic, not WT';metadata['no_recorded_injection']=metadata.injection_materials_label.eq('')&metadata.injection_roi_label.eq('')&metadata.injection_method_label.eq('')
 metadata['age_source']='GSE185862 original donor_id to unique reported age, unknown retained'
 panel=set(pd.read_csv(ROOT/'results/pbm/panel22.csv').gene);common=set(pd.read_csv(ROOT/'results/pbm/same_direction_52.csv').gene);candidates=sorted(panel|common)
 source_path=ROOT/'input/rnaseq/raw/GSE185862_expression_matrix_SSv4.hdf5'
 start=time.monotonic()
 with h5py.File(source_path,'r') as f:
  genes=pd.Index([x.decode() for x in f['gene_names'][:]]);cells=pd.Index([x.decode() for x in f['sample_names'][:]])
  shape=f['data/exon/dims'][:].tolist();transshape=f['data/t_exon/dims'][:].tolist()
  if shape!=[len(cells),len(genes)] or transshape!=shape[::-1] or cells.duplicated().any() or genes.duplicated().any():raise RuntimeError('Source identity or dimensions inconsistent')
  lookup=pd.Series(np.arange(len(cells)),index=cells);metadata['source_matrix_row']=metadata.cell_id.map(lookup)
  if metadata.source_matrix_row.isna().any() or metadata.cell_id.duplicated().any():raise RuntimeError('Cell identity fails exact unique alignment')
  metadata=metadata.sort_values('source_matrix_row').reset_index(drop=True);positions=metadata.source_matrix_row.to_numpy(dtype=np.int64)
  available=[g for g in candidates if g in genes];counts=np.zeros((len(metadata),len(available)),dtype=np.int64);reverse=np.full(len(cells),-1,dtype=np.int64);reverse[positions]=np.arange(len(positions))
  ptr=f['data/exon/p'][:]
  for j,gene in enumerate(available):
   gi=genes.get_loc(gene);a,b=int(ptr[gi]),int(ptr[gi+1]);ix=f['data/exon/i'][a:b];val=f['data/exon/x'][a:b]
   if not np.isfinite(val).all() or np.any(val<0) or np.any(val!=np.rint(val)) or np.any(ix!=np.rint(ix)):raise RuntimeError('Source values or indices invalid')
   ix=ix.astype(np.int64)
   if np.any(np.diff(ix)<=0) or np.any(ix<0) or np.any(ix>=len(cells)):raise RuntimeError('Sparse gene column indices invalid')
   local=reverse[ix];keep=local>=0;counts[local[keep],j]=val[keep].astype(np.int64)
  declared=f['data/total_exon_counts'][:][positions];tptr=f['data/t_exon/p'][:];totals=np.zeros(len(metadata),dtype=np.int64);validated=0
  for block in np.unique(positions//128):
   lo=int(block*128);hi=min(lo+128,len(cells));a,b=int(tptr[lo]),int(tptr[hi]);val=f['data/t_exon/x'][a:b]
   if not np.isfinite(val).all() or np.any(val<0) or np.any(val!=np.rint(val)):raise RuntimeError('Full source exon counts invalid')
   if np.any(np.diff(tptr[lo:hi+1])<=0):raise RuntimeError('Zero-library source cell requires explicit handling')
   rawtotals=np.add.reduceat(val,tptr[lo:hi].astype(np.int64)-a)
   keep=(positions>=lo)&(positions<hi);totals[keep]=rawtotals[positions[keep]-lo].astype(np.int64);validated+=int(keep.sum())
   if validated//4096!=(validated-int(keep.sum()))//4096:print(json.dumps({'validated_selected_libraries':validated,'total':len(metadata),'elapsed_s':round(time.monotonic()-start,1)}),flush=True)
  if not np.array_equal(totals,declared) or np.any(totals<=0) or np.any(counts.sum(axis=1)>totals):raise RuntimeError('Full exon libraries disagree with source declared totals')
  # Independently verify both sparse orientations for ten selected cells.
  for row in np.linspace(0,len(metadata)-1,10,dtype=int):
   ci=positions[row];a,b=int(tptr[ci]),int(tptr[ci+1]);ix=f['data/t_exon/i'][a:b].astype(np.int64);val=f['data/t_exon/x'][a:b];record=dict(zip(ix,val))
   for j,g in enumerate(available):
    if counts[row,j]!=record.get(genes.get_loc(g),0):raise RuntimeError('Sparse orientation candidate mismatch')
 metadata['full_exon_read_counts']=totals;metadata.to_csv(out/'selected_cell_metadata.csv.gz',index=False)
 h5ad=out/'mop_visp_candidate_exon_counts.h5ad';adata=ad.AnnData(sparse.csr_matrix(counts),obs=metadata.set_index('cell_id'),var=pd.DataFrame({'gene_symbol':available},index=available));adata.uns['value_kind']='integer exon-aligned uniquely mapped read counts, PCR duplicates removed per GEO processing';adata.uns['denominator']='full_exon_read_counts sums all 45768 source genes, independently verified';adata.write_h5ad(h5ad,compression='gzip')
 manifest={'status':'complete','n_source_metadata_cells':len(original),'n_source_matrix_cells':len(cells),'n_selected_cells':len(metadata),'n_selected_MOp':int(metadata.region.eq('MOp').sum()),'n_selected_VISp':int(metadata.region.eq('VISp').sum()),'n_source_genes':len(genes),'n_candidates_requested':len(candidates),'n_candidates_available':len(available),'missing_candidates':sorted(set(candidates)-set(available)),'full_library_validation':'Exact recomputed sum across all 45768 source exon-count features for every selected cell equals data/total_exon_counts','sparse_orientation_validation':'Independent candidate alignment in both CSC orientations for ten deterministically selected cells','missing_age_cells':int(metadata.age.eq('unknown').sum()),'cells_with_recorded_injection':int((~metadata.no_recorded_injection).sum()),'normalization':'log1p(exon read counts / full exon library * 1e6), per-cell logCPM; read counts, not UMI or TPM','source_matrix':'input/rnaseq/raw/GSE185862_expression_matrix_SSv4.hdf5','subset_file':str(h5ad.relative_to(ROOT)),'subset_sha256':hashlib.file_digest(h5ad.open('rb'),'sha256').hexdigest(),'source_matrix_dimensions_cell_gene':shape,'source_sparse_nonzero_entries':int(tptr[-1]),'animal_inference':'withheld; library batch and donor sampling reviewed descriptively','elapsed_s':round(time.monotonic()-start,1)}
 (out/'subset_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n');print(json.dumps(manifest),flush=True)
if __name__=='__main__':main()
