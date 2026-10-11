#!/usr/bin/env python3
"""Validate released-data SMART-seq identities, denominators and reported rows."""
import hashlib,json
from pathlib import Path
import anndata as ad
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
KEYS=['class','subclass','genotype','sex','age','facs']
def main():
 result=ROOT/'results/rnaseq/smartseq';subset=ROOT/'input/rnaseq/smartseq';a=ad.read_h5ad(subset/'mop_visp_candidate_exon_counts.h5ad');m=a.obs.reset_index();counts=a.X.toarray();totals=m.full_exon_read_counts.to_numpy();manifest=json.loads((subset/'subset_manifest.json').read_text());status=json.loads((result/'status.json').read_text())
 assert set(m.region)=={'MOp','VISp'} and len(m)==22160 and len(a.var_names)==62
 assert m.region.value_counts().to_dict()=={'VISp':15572,'MOp':6588}
 assert not m.cell_id.duplicated().any() and np.isfinite(counts).all() and (counts>=0).all() and (counts==np.rint(counts)).all()
 assert np.all(totals>0) and np.all(counts.sum(axis=1)<=totals)
 actualhash=hashlib.file_digest((subset/'mop_visp_candidate_exon_counts.h5ad').open('rb'),'sha256').hexdigest();assert actualhash==manifest['subset_sha256']
 assert manifest['n_source_genes']==45768 and manifest['missing_candidates']==[]
 portable=subset/'GSE185862_CTX_HPF_2021_SS_readme.txt';raw=ROOT/'input/rnaseq/raw/GSE185862_CTX_HPF_2021_SS_readme.txt';assert portable.read_bytes()==raw.read_bytes()
 table=pd.read_csv(result/'shared_strata_region_difference_min10.csv');key=pd.read_csv(result/'displayed_stratum_key.csv');std=pd.read_csv(result/'shared_cluster_standardized_difference_min10.csv')
 checked=[]
 for _,row in table.loc[table.gene.isin(['Adora1','Adora2a','Nt5e','Grm2','Penk'])].iloc[::19].iterrows():
  mask=np.ones(len(m),dtype=bool)
  for k in KEYS:mask&=m[k].astype(str).to_numpy()==str(row[k])
  mask&=m.no_recorded_injection.to_numpy(dtype=bool)
  j=a.var_names.get_loc(row.gene)
  for region in ['MOp','VISp']:
   ix=np.flatnonzero(mask&(m.region.astype(str).to_numpy()==region));expected=float(np.log1p(counts[ix,j]/totals[ix]*1e6).mean());expected10k=float(np.log1p(counts[ix,j]/totals[ix]*1e4).mean());expected_detection=float((counts[ix,j]>0).mean())
   assert len(ix)==row['n_cells_'+region]
   assert np.isclose(expected,row['mean_expression_'+region],rtol=0,atol=1e-12)
   assert np.isclose(expected10k,row['mean_log1p_reads_CP10k_'+region],rtol=0,atol=1e-12)
   assert np.isclose(expected_detection,row['detection_fraction_'+region],rtol=0,atol=1e-12)
  checked.append({k:row[k] for k in KEYS+['gene']})
 assert len(key)==38 and len(table)==38*62
 assert (table.n_cells_MOp>=10).all() and (table.n_cells_VISp>=10).all()
 assert status['n_shared_metadata_strata_unfiltered']==67 and status['n_same_dissection_layer_strata_min10_display']==22
 assert set(key.age)=={'P53','P54','P55','P56','P59'}
 cross=pd.read_csv(result/'exact_metadata_10x_direction_crosscheck.csv');assert len(cross)==0 and status['n_crosscheck_10x_exact_metadata_gene_rows']==0
 assert len(pd.read_csv(result/'gene_direction_concordance_summary.csv'))==0
 assert not any(('pvalue' in col.lower() or col.lower() in ['p_value','pval','fdr']) for col in table.columns)
 readme_record=next(x for x in json.loads((subset/'source_manifest.json').read_text()) if x['source_url'].endswith('GSE185862_CTX_HPF_2021_SS_readme.txt'))
 assert readme_record['sha256']==hashlib.file_digest(portable.open('rb'),'sha256').hexdigest()
 # Recompute the shared-cluster standardized effect for one entire stratum.
 cluster=pd.read_csv(result/'shared_cluster_region_difference_unfiltered.csv.gz');row=std.loc[std.gene.eq('Adora1')].iloc[0];mask=cluster.both_regions_observed.to_numpy(copy=True)
 for k in KEYS+['gene']:mask&=cluster[k].astype(str).to_numpy()==str(row[k])
 block=cluster.loc[mask];w=(block.n_cells_MOp+block.n_cells_VISp);w=w/w.sum();expected=float((w*(block.mean_expression_MOp-block.mean_expression_VISp)).sum());assert np.isclose(expected,row.MOp_minus_VISp_standardized_mean_expression,rtol=0,atol=1e-12)
 payload={'status':'passed','data_kind':'actual public GEO SMART-seq counts, not synthetic','source_all_gene_library_validation':manifest['full_library_validation'],'sparse_orientation_validation':manifest['sparse_orientation_validation'],'n_candidate_genes':62,'n_selected_cells':22160,'n_independently_recomputed_gene_stratum_rows':len(checked),'checks':['exact MOp/VISp identities and counts','integer nonnegative read counts','subset hash and readme source hash','all-source-gene library denominator recorded and positive','reported mean logCPM/logCP10k/detection recomputed from subset','shared-cluster weighted mean difference independently recomputed','all67shared strata retained and38display filter strata','actual age metadata and22same-layer strata','no false exact-condition10x cross-assay results','no P values/FDR invented'],'recomputed_gene_stratum_rows':checked,'strict_direction_cross_assay_match_rows':0,'minimum_cells_note':'display filter only; no biological-replicate independence supplied','animal_inference':'not calculated','portable_result_bytes':sum(p.stat().st_size for p in result.rglob('*') if p.is_file())}
 output=ROOT/'results/validation/smartseq_actual_data_validation.json';output.parent.mkdir(parents=True,exist_ok=True);output.write_text(json.dumps(payload,ensure_ascii=False,indent=2)+'\n');print(json.dumps({k:v for k,v in payload.items() if k!='recomputed_gene_stratum_rows'},ensure_ascii=False),flush=True)
if __name__=='__main__':main()
