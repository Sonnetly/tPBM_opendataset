"""Download public GEO/QuickGO source snapshots; preserve provenance and hashes."""
from pathlib import Path
import concurrent.futures
import hashlib
import json
import urllib.parse
import urllib.request
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[2]
DEST = ROOT / 'ora_analysis/input/ros_geo'
DEST.mkdir(parents=True, exist_ok=True)
SOURCES = {}
for accession in ['GSE245658', 'GSE325653']:
    base = f'https://ftp.ncbi.nlm.nih.gov/geo/series/{accession[:-3]}nnn/{accession}'
    SOURCES[f'{accession}_series.txt'] = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}&targ=self&form=text&view=full'
    SOURCES[f'{accession}_samples.txt'] = f'https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}&targ=gsm&form=text&view=full'
SOURCES['GSE245658_all.genes.expression.annot.LT.txt.gz'] = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/GSE245658/suppl/GSE245658_all.genes.expression.annot.LT.txt.gz'
SOURCES['GSE325653_raw_counts.xlsx'] = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE325nnn/GSE325653/suppl/GSE325653_raw_counts.xlsx'

# is_a/part_of descendants, not regulatory descendants; NOT annotations excluded.
TERMS = {
    'ros_metabolism': 'GO:0072593',
    'antioxidant_activity': 'GO:0016209',
    'oxidative_stress_response': 'GO:0006979',
    'glutathione_metabolism': 'GO:0006749',
    'ros_biosynthesis': 'GO:1903409',
    'hydrogen_peroxide_catabolism': 'GO:0042744',
    'superoxide_removal': 'GO:0019430',
    'ros_metabolism_regulation': 'GO:2000377',
    'nadph_regeneration': 'GO:0006740',
}
for name, term in TERMS.items():
    params = urllib.parse.urlencode(dict(taxonId='10090', taxonUsage='exact', goId=term,
        goUsage='descendants', goUsageRelationships='is_a,part_of', geneProductType='protein'))
    SOURCES[f'quickgo_{name}.tsv'] = 'https://www.ebi.ac.uk/QuickGO/services/annotation/downloadSearch?' + params
    SOURCES[f'quickgo_{name}_term.json'] = 'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/' + term
SOURCES['kegg_glutathione_mouse.tsv'] = 'https://rest.kegg.jp/link/mmu/mmu00480'
SOURCES['kegg_mouse_gene_names.tsv'] = 'https://rest.kegg.jp/list/mmu'

def download(item):
    name, url = item
    path = DEST / name
    record = {'file': name, 'url': url, 'retrieved_utc': datetime.now(timezone.utc).isoformat()}
    try:
        if not path.exists():
            req = urllib.request.Request(url, headers={'Accept': 'text/tsv' if name.endswith('.tsv') else '*/*', 'User-Agent': 'ROS-research/1.0'})
            for attempt in range(3):
                try:
                    data = urllib.request.urlopen(req, timeout=180).read()
                    break
                except Exception:
                    if attempt == 2: raise
            if not data: raise ValueError('Empty download')
            path.write_bytes(data)
        data = path.read_bytes()
        record.update(bytes=len(data), sha256=hashlib.sha256(data).hexdigest(), status='ok')
    except Exception as exc:
        record.update(status='error', error=str(exc))
    print(name, record['status'], record.get('bytes', record.get('error')), flush=True)
    return record

if __name__ == '__main__':
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as pool:
        records = list(pool.map(download, SOURCES.items()))
    (DEST / 'source_manifest.json').write_text(json.dumps(records, indent=2) + '\n')
