"""Independent numerical and provenance checks for the ROS analysis."""
from pathlib import Path
import hashlib
import itertools
import json
import numpy as np
import pandas as pd
from scipy import stats
from analyze_ros import enrichment_score

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
RES=HERE/'results'

def main():
    # Validate the sparse GSEA extrema against a dense independent running sum,
    # including every hit pattern in a small universe and both signs.
    weights=np.array([4.2,3.1,2.9,2.2,1.7,1.4,.8,.2])
    patterns=0
    for k in [1,2,3,5,7]:
        for pos in itertools.combinations(range(8),k):
            mask=np.zeros(8,dtype=bool);mask[list(pos)]=True
            increments=np.where(mask,weights/weights[mask].sum(),-1/(8-k))
            walk=np.r_[0,np.cumsum(increments)]
            expected=walk.max() if walk.max()>=abs(walk.min())-1e-12 else walk.min()
            actual,_=enrichment_score(np.array(pos),weights)
            assert np.isclose(actual,expected,atol=1e-12)
            patterns+=1
    # Cross-check all ORA p values via Fisher's exact contingency-table test.
    ora=pd.read_csv(RES/'focused_GEO_ORA.csv')
    for r in ora.itertuples():
        matrix=[[r.k,r.n-r.k],[r.K-r.k,r.N-r.K-r.n+r.k]]
        assert min(v for row in matrix for v in row)>=0
        p=stats.fisher_exact(matrix,alternative='greater').pvalue
        assert np.isclose(p,r.pvalue,rtol=1e-9,atol=1e-12)
    # Cross-check every sample-level test with SciPy's independent exhaustive engine.
    scores=pd.read_csv(RES/'sample_level_pathway_scores.csv')
    tests=pd.read_csv(RES/'sample_level_exact_permutation.csv')
    for r in tests.itertuples():
        d=scores[(scores.dataset==r.dataset)&(scores.set_name==r.set_name)]
        result=stats.permutation_test((d[d.condition=='PBM'].score.to_numpy(),d[d.condition=='control'].score.to_numpy()),
            lambda a,b:np.mean(a)-np.mean(b),vectorized=False,n_resamples=np.inf,alternative='two-sided')
        assert np.isclose(result.pvalue,r.exact_two_sided_p,atol=1e-12)
    # Validate original-source snapshot integrity and the stored analysis input hashes.
    manifest=json.loads((ROOT/'ora_analysis/input/ros_geo/source_manifest.json').read_text())
    assert all(r['status']=='ok' for r in manifest)
    for r in manifest:
        path=ROOT/'ora_analysis/input/ros_geo'/r['file']
        assert hashlib.sha256(path.read_bytes()).hexdigest()==r['sha256']
    summary=json.loads((RES/'validation_summary.json').read_text())
    for path,digest in summary['input_sha256'].items():
        assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
    # The current supplied-table cross-study conclusion is independently bounded
    # by the pre-existing all-gene overlap and both gene FDRs.
    common=pd.read_csv(RES/'ROS_same_direction_DEG.csv')
    prior=pd.read_csv(ROOT/'ora_analysis/deg_validation/same_direction_52.csv')
    assert set(common.gene_symbol)<=set(prior.gene)
    assert (common.FDR_padj_A<.05).all() and (common.FDR_padj_B<.05).all()
    assert (common.direction_A==common.direction_B).all()
    assert common.gene_symbol.tolist()==['Penk']
    assert not (common.antioxidant_activity | common.hydrogen_peroxide_catabolism | common.superoxide_removal).any()
    panel=pd.read_csv(RES/'canonical_and_support_reference_panel.csv')
    b=pd.read_csv(RES/'GSE325653_GEO_full_statistics.csv').set_index('gene_symbol')
    assert np.isclose(panel.set_index('gene_symbol').loc['Gpx1','GSE325653_GEO_log2FC'],b.loc['Gpx1','log2FC'])
    # Confirm all expected formats are nonempty; Arial is checked externally with pdffonts.
    figures=list((HERE/'figures').glob('*.png'))
    assert len(figures)==8
    for p in figures:
        for suffix in ['.png','.pdf','.svg']:
            assert p.with_suffix(suffix).stat().st_size>1000
    font=json.loads((RES/'font_provenance.json').read_text())
    assert font['family']=='Arial'
    out=dict(sparse_GSEA_vs_dense_patterns=patterns,Fisher_vs_hypergeom_tests=len(ora),
        exact_permutation_vs_scipy_tests=len(tests),source_files_verified=len(manifest),
        figure_count=len(figures),figure_formats=['PNG','PDF','SVG'],font='Arial',status='passed')
    (RES/'independent_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))

if __name__=='__main__':main()
