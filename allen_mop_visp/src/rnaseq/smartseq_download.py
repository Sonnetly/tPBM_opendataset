#!/usr/bin/env python3
"""Fetch the authoritative GEO SMART-seq release without bypassing TLS checks."""
import argparse, hashlib, json, time, urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
BASE='https://ftp.ncbi.nlm.nih.gov/geo/series/GSE185nnn/GSE185862/suppl/'
FILES=['GSE185862_CTX_HPF_2021_SS_readme.txt','GSE185862_metadata_ssv4.csv.gz','GSE185862_expression_matrix_SSv4.hdf5']
def main():
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('filename',choices=FILES);args=parser.parse_args()
 path=ROOT/'input/rnaseq/raw'/args.filename;url=BASE+args.filename
 with urllib.request.urlopen(urllib.request.Request(url,method='HEAD'),timeout=60) as response: expected=int(response.headers['Content-Length'])
 records_path=ROOT/'input/rnaseq/smartseq/source_manifest.json'
 records=json.loads(records_path.read_text()) if records_path.exists() else []
 if path.exists():
  if path.stat().st_size!=expected: raise RuntimeError('Existing final file length mismatch; preserve and diagnose it')
  digest=hashlib.file_digest(path.open('rb'),'sha256').hexdigest(); headers={}
 else:
  partial=path.with_name(path.name+'.part')
  if partial.exists(): raise RuntimeError('Partial file exists; diagnose before retry')
  h=hashlib.sha256();n=0;start=time.monotonic()
  with urllib.request.urlopen(url,timeout=60) as response,partial.open('wb') as target:
   headers=dict(response.headers)
   while block:=response.read(8*1024*1024):
    target.write(block);h.update(block);n+=len(block)
    if n//(256*1024*1024)!=(n-len(block))//(256*1024*1024): print(json.dumps({'bytes':n,'expected':expected,'elapsed_s':round(time.monotonic()-start,1)}),flush=True)
  if n!=expected: raise RuntimeError('Download truncated; partial preserved')
  partial.rename(path);digest=h.hexdigest()
 record={'source_url':url,'local_path':str(path.relative_to(ROOT)),'bytes':expected,'sha256':digest,'TLS_verification':'preserved','publisher_checksum_available':False,'response_headers':headers}
 records=[x for x in records if x['source_url']!=url]+[record]
 records_path.write_text(json.dumps(records,ensure_ascii=False,indent=2)+'\n')
 print(json.dumps(record),flush=True)
if __name__=='__main__':main()
