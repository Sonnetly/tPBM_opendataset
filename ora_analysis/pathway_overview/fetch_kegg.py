"""Download/cache the mouse KEGG catalog without changing previous ROS snapshots."""
from pathlib import Path
import concurrent.futures
import datetime
import hashlib
import json
import urllib.request

DEST=Path(__file__).resolve().parents[2]/'ora_analysis/input/pathway_overview'
SOURCES={'kegg_mouse_pathways.tsv':'https://rest.kegg.jp/list/pathway/mmu',
 'kegg_mouse_gene_pathways.tsv':'https://rest.kegg.jp/link/pathway/mmu'}

def fetch(item):
    name,url=item;path=DEST/name
    # Preserve original retrieval timestamps when reusing saved snapshots.
    old=json.loads((DEST/'source_manifest.json').read_text()) if (DEST/'source_manifest.json').exists() else []
    previous=next((r for r in old if r['file']==name),None)
    if not path.exists():path.write_bytes(urllib.request.urlopen(url,timeout=120).read())
    data=path.read_bytes();digest=hashlib.sha256(data).hexdigest()
    if previous and previous['sha256']==digest:return previous
    return dict(file=name,url=url,sha256=digest,bytes=len(data),retrieved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())

if __name__=='__main__':
    DEST.mkdir(parents=True,exist_ok=True)
    with concurrent.futures.ThreadPoolExecutor(2) as pool:records=list(pool.map(fetch,SOURCES.items()))
    (DEST/'source_manifest.json').write_text(json.dumps(records,indent=2)+'\n')
    print('KEGG snapshots verified:',len(records))
