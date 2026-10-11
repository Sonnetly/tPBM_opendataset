#!/usr/bin/env python3
"""Descriptive reporter-stratified SMART-seq MOp/VISp cross-check.

No cell or animal P values are computed. A minimum cell number is only a
small-sample display filter; every shared stratum remains in the unfiltered CSV.
"""
import json
import os
from pathlib import Path
import anndata as ad
os.environ.setdefault('MPLCONFIGDIR', '/tmp/allen_smartseq_matplotlib')
os.environ.setdefault('XDG_CACHE_HOME', '/tmp/allen_smartseq_cache')
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parents[2]
RESULT=ROOT/'results/rnaseq/smartseq';FIG=ROOT/'figures/rnaseq/smartseq'
KEYS=['class','subclass','genotype','sex','age','facs']
VALUE_KIND='log1p(exon read CPM); full 45768-gene exon library denominator'

def summarize(meta,raw,genes,libs,keys,panel,common):
 records=[];cp10k=np.log1p(raw/libs[:,None]*1e4);cpm=np.log1p(raw/libs[:,None]*1e6)
 for labels,ix in meta.groupby(keys,dropna=False).indices.items():
  ix=np.asarray(ix);labels=(labels,) if len(keys)==1 else labels
  for j,gene in enumerate(genes):
   records.append(dict(zip(keys,labels),gene=gene,n_cells=len(ix),mean_expression=float(cpm[ix,j].mean()),mean_log1p_reads_CP10k=float(cp10k[ix,j].mean()),detection_fraction=float((raw[ix,j]>0).mean()),raw_read_sum=int(raw[ix,j].sum()),in_panel22=gene in panel,in_common52=gene in common,value_kind=VALUE_KIND))
 return pd.DataFrame(records)

def compare(summary,keys):
 a=summary.loc[summary.region.eq('MOp')].drop(columns='region');b=summary.loc[summary.region.eq('VISp')].drop(columns='region')
 merged=a.merge(b,on=keys+['gene'],how='outer',suffixes=('_MOp','_VISp'),validate='one_to_one')
 merged['both_regions_observed']=merged.n_cells_MOp.notna()&merged.n_cells_VISp.notna()
 for measure in ['mean_expression','mean_log1p_reads_CP10k','detection_fraction']:merged['MOp_minus_VISp_'+measure]=merged[measure+'_MOp']-merged[measure+'_VISp']
 merged['analysis_level']='descriptive cell-weighted summaries; no animal-level significance test'
 return merged

def main():
 RESULT.mkdir(parents=True,exist_ok=True);FIG.mkdir(parents=True,exist_ok=True)
 artifact=ad.read_h5ad(ROOT/'input/rnaseq/smartseq/mop_visp_candidate_exon_counts.h5ad');meta=artifact.obs.reset_index();genes=artifact.var_names.tolist();raw=artifact.X.toarray().astype(np.int64);libs=meta.full_exon_read_counts.to_numpy(dtype=np.int64)
 panel=set(pd.read_csv(ROOT/'results/pbm/panel22.csv').gene);common=set(pd.read_csv(ROOT/'results/pbm/same_direction_52.csv').gene)
 if not set(meta.region)=={'MOp','VISp'} or len(set(meta.cell_id))!=len(meta):raise RuntimeError('Exact region/identity validation failed')
 meta.groupby(['region']+KEYS+['cortical_layer_label'],dropna=False).size().rename('n_cells').reset_index().to_csv(RESULT/'all_selected_metadata_counts.csv',index=False)
 for unit in ['donor','batch']:
  meta.groupby(['region',unit]+KEYS,dropna=False).size().rename('n_cells').reset_index().to_csv(RESULT/f'{unit}_region_stratum_counts.csv',index=False)
 # Keep all source cells, including missing-age cells, in this descriptive audit.
 full=summarize(meta,raw,genes,libs,['region']+KEYS,panel,common);full.to_csv(RESULT/'all_selected_region_subclass_gene_summary.csv.gz',index=False)
 allregion=summarize(meta,raw,genes,libs,['region'],panel,common);allregion.to_csv(RESULT/'all_selected_region_gene_summary_audit_only.csv',index=False)
 missing=meta.age.eq('unknown');missing_summary=summarize(meta.loc[missing].reset_index(drop=True),raw[missing],genes,libs[missing],['region','class','subclass'],panel,common) if missing.any() else pd.DataFrame();missing_summary.to_csv(RESULT/'missing_age_candidate_detection_audit.csv',index=False)
 known=~meta[KEYS].isin(['unknown','ZZ_Missing','']).any(axis=1);uninjected=meta.no_recorded_injection.astype(bool);strict=known&uninjected
 clean=meta.loc[strict].reset_index(drop=True);values=raw[strict];libraries=libs[strict]
 clean.to_csv(ROOT/'input/rnaseq/smartseq/strict_known_metadata_no_recorded_injection.csv.gz',index=False)
 summary=summarize(clean,values,genes,libraries,['region']+KEYS,panel,common);summary.to_csv(RESULT/'strict_region_subclass_gene_summary.csv.gz',index=False)
 diff=compare(summary,KEYS);diff.to_csv(RESULT/'within_subclass_region_difference_unfiltered.csv.gz',index=False)
 shared=diff.loc[diff.both_regions_observed].copy();shared.to_csv(RESULT/'shared_strata_region_difference_unfiltered.csv',index=False)
 displayed=shared.loc[(shared.n_cells_MOp>=10)&(shared.n_cells_VISp>=10)].copy();displayed.to_csv(RESULT/'shared_strata_region_difference_min10.csv',index=False)
 strata=displayed[KEYS+['n_cells_MOp','n_cells_VISp']].drop_duplicates().sort_values(KEYS).reset_index(drop=True);strata['stratum_id']=[f'S{i+1:03d}' for i in range(len(strata))];strata.to_csv(RESULT/'displayed_stratum_key.csv',index=False)
 # Actual anatomical dissection-layer labels remain separate from subclasses.
 layerkeys=KEYS+['cortical_layer_label'];layer_summary=summarize(clean,values,genes,libraries,['region']+layerkeys,panel,common);layerdiff=compare(layer_summary,layerkeys);layerdiff.to_csv(RESULT/'same_dissection_layer_region_difference_unfiltered.csv.gz',index=False)
 layerdisplay=layerdiff.loc[layerdiff.both_regions_observed&(layerdiff.n_cells_MOp>=10)&(layerdiff.n_cells_VISp>=10)];layerdisplay.to_csv(RESULT/'same_dissection_layer_region_difference_min10.csv',index=False)
 # Restrict to clusters shared in each metadata stratum and use identical weights.
 csum=summarize(clean,values,genes,libraries,['region']+KEYS+['cluster'],panel,common);csum.to_csv(RESULT/'strict_region_cluster_gene_summary.csv.gz',index=False)
 cdiff=compare(csum,KEYS+['cluster']);cdiff.to_csv(RESULT/'shared_cluster_region_difference_unfiltered.csv.gz',index=False)
 cshared=cdiff.loc[cdiff.both_regions_observed].copy();std=[]
 for labels,g in cshared.groupby(KEYS+['gene'],dropna=False):
  weights=(g.n_cells_MOp.to_numpy()+g.n_cells_VISp.to_numpy());weights/=weights.sum()
  record=dict(zip(KEYS+['gene'],labels),n_shared_clusters=len(g),n_cells_MOp=int(g.n_cells_MOp.sum()),n_cells_VISp=int(g.n_cells_VISp.sum()),weight_definition='pooled MOp/VISp cell counts over shared clusters, identical weights for both regions')
  for measure in ['mean_expression','mean_log1p_reads_CP10k','detection_fraction']:
   for region in ['MOp','VISp']:record['standardized_'+measure+'_'+region]=float(np.dot(weights,g[measure+'_'+region]))
   record['MOp_minus_VISp_standardized_'+measure]=record['standardized_'+measure+'_MOp']-record['standardized_'+measure+'_VISp']
  std.append(record)
 standard=pd.DataFrame(std);standard.to_csv(RESULT/'shared_cluster_standardized_difference_unfiltered.csv',index=False)
 if len(standard):standard.loc[(standard.n_cells_MOp>=10)&(standard.n_cells_VISp>=10)].to_csv(RESULT/'shared_cluster_standardized_difference_min10.csv',index=False)
 # Cross-assay directions only: same gene/subclass/reporter/sex/age/FACS;
 # relative read CP10k and UMI CP10k are retained in separate columns.
 primary=pd.read_csv(ROOT/'results/rnaseq/within_subclass_region_difference.csv')
 primary=primary.loc[primary.both_regions_observed&(primary.n_cells_MOp>=10)&(primary.n_cells_VISp>=10)]
 harmonized=displayed.merge(primary,on=KEYS+['gene'],how='inner',suffixes=('_SMARTseq','_10x'),validate='one_to_one')
 if len(harmonized):
  harmonized['SMARTseq_MOp_minus_VISp_log1p_reads_CP10k']=harmonized['MOp_minus_VISp_mean_log1p_reads_CP10k'];harmonized['10x_MOp_minus_VISp_log1p_UMI_CP10k']=harmonized['MOp_minus_VISp_mean_expression_10x']
  harmonized['same_direction']=np.sign(harmonized.SMARTseq_MOp_minus_VISp_log1p_reads_CP10k)==np.sign(harmonized['10x_MOp_minus_VISp_log1p_UMI_CP10k'])
  harmonized['comparison_limit']='directional exploratory check; no pooling of read and UMI measurements; library batch and dissection are not jointly matched'
 harmonized.to_csv(RESULT/'exact_metadata_10x_direction_crosscheck.csv',index=False)
 concordance=[]
 if len(harmonized):
  for gene,g in harmonized.groupby('gene'):
   usable=g.loc[(g.SMARTseq_MOp_minus_VISp_log1p_reads_CP10k.abs()>=0.05)&(g['10x_MOp_minus_VISp_log1p_UMI_CP10k'].abs()>=0.05)]
   concordance.append({'gene':gene,'in_panel22':gene in panel,'n_exact_metadata_subclasses':len(g),'n_same_direction':int(g.same_direction.sum()),'same_direction_fraction':float(g.same_direction.mean()),'n_both_abs_difference_ge_0p05_display_sensitivity':len(usable),'n_same_direction_in_display_sensitivity':int(usable.same_direction.sum()),'threshold_note':'0.05 removes tiny descriptive mean-log differences only; not a statistical significance threshold'})
 concordance_columns=['gene','in_panel22','n_exact_metadata_subclasses','n_same_direction','same_direction_fraction','n_both_abs_difference_ge_0p05_display_sensitivity','n_same_direction_in_display_sensitivity','threshold_note']
 pd.DataFrame(concordance,columns=concordance_columns).to_csv(RESULT/'gene_direction_concordance_summary.csv',index=False)
 # Counts of descriptive direction patterns across heterogeneous strata are
 # an audit, not independent replicates or a pooled biological region effect.
 sensitivity=displayed.merge(standard,on=KEYS+['gene'],how='left',validate='one_to_one');pattern_rows=[]
 for gene,g in sensitivity.loc[sensitivity.gene.isin(panel)].groupby('gene'):
  delta=g.MOp_minus_VISp_mean_expression;both=g.dropna(subset=['MOp_minus_VISp_standardized_mean_expression']);cs=both.MOp_minus_VISp_standardized_mean_expression
  pattern_rows.append({'gene':gene,'n_strata':len(g),'MOp_higher_n':int((delta>0).sum()),'VISp_higher_n':int((delta<0).sum()),'equal_n':int((delta==0).sum()),'min_MOp_detection':g.detection_fraction_MOp.min(),'max_MOp_detection':g.detection_fraction_MOp.max(),'min_VISp_detection':g.detection_fraction_VISp.min(),'max_VISp_detection':g.detection_fraction_VISp.max(),'sharedcluster_std_strata':len(cs),'std_MOp_higher_n':int((cs>0).sum()),'std_VISp_higher_n':int((cs<0).sum()),'direction_same_after_std_n':int((np.sign(both.MOp_minus_VISp_mean_expression)==np.sign(cs)).sum())})
 pd.DataFrame(pattern_rows).to_csv(RESULT/'panel22_direction_sensitivity_audit.csv',index=False)
 # A heatmap of the 22 prespecified genes, with the full stratum key beside CSVs.
 if len(strata):
  shown=displayed.loc[displayed.gene.isin(panel)].merge(strata,on=KEYS+['n_cells_MOp','n_cells_VISp'],validate='many_to_one');matrix=shown.pivot(index='stratum_id',columns='gene',values='MOp_minus_VISp_mean_expression').reindex(strata.stratum_id).reindex(columns=sorted(panel));limit=np.nanpercentile(np.abs(matrix.to_numpy()),98) or 1
  fig,ax=plt.subplots(figsize=(14,max(6,len(strata)*0.27)));im=ax.imshow(matrix.to_numpy(),aspect='auto',cmap='RdBu_r',vmin=-limit,vmax=limit)
  ax.set_xticks(range(len(matrix.columns)),matrix.columns,rotation=90);ax.set_yticks(range(len(strata)),[f'{r.stratum_id} | {r.subclass} | {r.sex} {r.age}' for r in strata.itertuples()],fontsize=7);ax.set_title('SMART-seq: MOp minus VISp mean log1p(read CPM)\nGenotype/sex/exact-age/FACS/subclass strata; no recorded injection; n >= 10 display filter');fig.colorbar(im,ax=ax,label='Descriptive mean-log difference (no animal P values)');fig.tight_layout();fig.savefig(FIG/'panel22_strict_metadata_difference_heatmap.png',dpi=180);fig.savefig(FIG/'panel22_strict_metadata_difference_heatmap.svg');plt.close(fig)
 # Portable workbook avoids duplicating large raw matrix/whole-atlas inputs.
 with pd.ExcelWriter(RESULT/'SMARTseq_MOp_VISp_candidate_crosscheck.xlsx') as writer:
  strata.to_excel(writer,sheet_name='stratum_key',index=False);displayed.to_excel(writer,sheet_name='matched62_min10',index=False);layerdisplay.to_excel(writer,sheet_name='same_layer_min10',index=False);harmonized.to_excel(writer,sheet_name='exact_10x_direction',index=False);pd.DataFrame(concordance).to_excel(writer,sheet_name='direction_summary',index=False);allregion.to_excel(writer,sheet_name='allcells_audit_only',index=False)
 shared_strata=shared[KEYS].drop_duplicates();sharedclean=clean.merge(shared_strata,on=KEYS,how='inner',validate='many_to_one');displayclean=clean.merge(strata[KEYS],on=KEYS,how='inner',validate='many_to_one')
 donor_sets={region:set(displayclean.loc[displayclean.region.eq(region)].donor.astype(str)) for region in ['MOp','VISp']};batch_sets={region:set(displayclean.loc[displayclean.region.eq(region)].batch.astype(str)) for region in ['MOp','VISp']}
 for unit in ['donor','batch']:
  displayclean.groupby(['region',unit]+KEYS,dropna=False).size().rename('n_cells').reset_index().to_csv(RESULT/f'displayed_{unit}_region_stratum_counts.csv',index=False)
 status={'status':'complete_descriptive_supplement','source':'GEO GSE185862 SMART-seq v4 exon read counts','exact_regions':['MOp','VISp'],'n_selected_cells':len(meta),'n_candidates_measured':len(genes),'strict_known_metadata_no_recorded_injection_cells':len(clean),'n_shared_metadata_strata_unfiltered':len(shared_strata),'shared_metadata_cells_unfiltered':sharedclean.region.value_counts().to_dict(),'n_shared_metadata_strata_min10_display':len(strata),'displayed_metadata_cells':displayclean.region.value_counts().to_dict(),'displayed_ages':sorted(strata.age.unique()),'displayed_donor_counts':{region:len(v) for region,v in donor_sets.items()},'displayed_shared_donor_ids':sorted(donor_sets['MOp']&donor_sets['VISp']),'displayed_shared_library_batch_ids':sorted(batch_sets['MOp']&batch_sets['VISp']),'n_same_dissection_layer_strata_min10_display':len(layerdisplay[layerkeys].drop_duplicates()),'n_crosscheck_10x_exact_metadata_gene_rows':len(harmonized),'n_crosscheck_10x_exact_metadata_subclasses':len(harmonized[KEYS].drop_duplicates()) if len(harmonized) else 0,'missing_age_cells_excluded_from_strict_comparison':int(missing.sum()),'recorded_injection_cells_excluded_from_strict_comparison':int((~uninjected).sum()),'normalization':VALUE_KIND,'unit':'PCR-duplicate-removed uniquely mapped exon read counts, not UMI or TPM','minimum_cells_policy':'n>=10 is only a display filter to omit very small strata; no inferential threshold','genotype_policy':'B6J background transgenic reporter lines; no WT comparison','anatomical_layer_policy':'Original dissection labels retained separately; same-dissection-layer sensitivity supplied','animal_inference':'not performed; cells are not independent animals; region/library-batch/dissection confounding unresolved','cross_assay_policy':'Separate read and UMI assays; direction-only comparison at same metadata stratum using log1p(relative CP10k) in separate columns'}
 (RESULT/'status.json').write_text(json.dumps(status,ensure_ascii=False,indent=2)+'\n');print(json.dumps(status,ensure_ascii=False),flush=True)
if __name__=='__main__':main()
