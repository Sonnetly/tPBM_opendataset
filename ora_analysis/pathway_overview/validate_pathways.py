"""Independent checks of whole-pathway statistics, sources and Arial exports."""
import hashlib
import json
import re
import subprocess
import numpy as np
import pandas as pd
from scipy import stats
from analyze_pathways import HERE,ROOT,OLD,OUT,INP,DATASETS

def main():
    catalog=pd.read_csv(OUT/'gene_set_catalog.csv').set_index('set_id')
    counts={'dense_ES_checks':0,'Fisher_ORA_checks':0,'ROS_sample_exact_checks':0,'focused_ES_comparisons':0}
    for acc in DATASETS:
        g=pd.read_csv(OUT/acc/'GSEA_all_pathways.csv');o=pd.read_csv(OUT/acc/'ORA_all_pathways.csv')
        ranks=pd.read_csv(OLD/f'{acc}_GSEA_ranked_genes.csv');symbols=ranks.gene_symbol.to_numpy();weights=abs(ranks.signed_z.to_numpy());N=len(ranks)
        assert g.set_id.is_unique and (g.nominal_p>0).all()
        assert np.allclose(stats.false_discovery_control(g.nominal_p),g.FDR_KEGG_plus_ROS_per_dataset)
        assert np.allclose(stats.false_discovery_control(o.pvalue),o.FDR_all_sets_directions_per_dataset)
        for r in g.itertuples():
            hit=np.isin(symbols,catalog.loc[r.set_id,'genes'].split(';'))
            increments=np.where(hit,weights/weights[hit].sum(),-1/(N-hit.sum()))
            walk=np.r_[0,np.cumsum(increments)]
            expected=walk.max() if walk.max()>=abs(walk.min())-1e-12 else walk.min()
            assert np.isclose(expected,r.ES,atol=1e-11)
            assert hit.sum()==r.matched_genes
            assert np.sign(r.ES)==np.sign(r.NES)
            counts['dense_ES_checks']+=1
        for r in o.itertuples():
            matrix=[[r.k,r.n-r.k],[r.K-r.k,r.N-r.K-r.n+r.k]]
            p=stats.fisher_exact(matrix,alternative='greater').pvalue
            assert np.isclose(p,r.pvalue,atol=1e-12,rtol=1e-9)
            counts['Fisher_ORA_checks']+=1
        # Independent SciPy exact enumeration for every selected ROS/context set.
        scores=pd.read_csv(OUT/acc/'sample_pathway_scores.csv')
        for r in g[g.ROS_related].itertuples():
            d=scores[scores.set_id==r.set_id]
            p=stats.permutation_test((d[d.condition=='PBM'].score.to_numpy(),d[d.condition=='control'].score.to_numpy()),
                lambda x,y:np.mean(x)-np.mean(y),vectorized=False,n_resamples=np.inf,alternative='two-sided').pvalue
            assert np.isclose(p,r.sample_exact_p,atol=1e-12)
            counts['ROS_sample_exact_checks']+=1
        old=pd.read_csv(OLD/'focused_GSEA.csv');old=old[old.dataset==acc]
        definitions=pd.read_csv(OLD/'gene_set_definitions.csv')
        for r in old.itertuples():
            sid=definitions.loc[definitions.set_name==r.Term,'root_id'].iloc[0]
            now=g[g.set_id==sid].iloc[0]
            assert np.isclose(now.ES,r.ES,atol=1e-12)
            assert now.matched_genes==r.matched_genes
            counts['focused_ES_comparisons']+=1
    provenance=json.loads((OUT/'provenance.json').read_text())
    for path,digest in provenance['SHA256'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest
    for record in json.loads((INP/'source_manifest.json').read_text()):
        assert hashlib.sha256((INP/record['file']).read_bytes()).hexdigest()==record['sha256']
    pdfs=sorted((HERE/'figures').glob('*.pdf'));assert len(pdfs)==10
    font_checks=[]
    for pdf in pdfs:
        data=subprocess.check_output(['pdffonts',str(pdf)],text=True).splitlines()[2:]
        assert data and all('Arial' in r and re.search(r'\byes\s+yes\s+yes\b',r) for r in data),pdf
        for ext in ['.png','.svg']:assert pdf.with_suffix(ext).stat().st_size>1000
        font_checks.append(pdf.name)
    for name in ['Dataset_pathway_report_ko.md','README.md']:
        text=(HERE/name).read_text()
        for target in re.findall(r'\]\(([^)]+)\)',text):
            if not target.startswith(('https:','http:')):assert (HERE/target).exists(),target
    counts.update(status='passed',Arial_embedded_PDFs=len(pdfs),source_hashes_verified=len(provenance['SHA256']))
    (OUT/'independent_checks.json').write_text(json.dumps(counts,indent=2)+'\n')
    print(json.dumps(counts,indent=2))

if __name__=='__main__':main()
