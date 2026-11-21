"""Compress a completed engineering product reproducibly without changing data."""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
from pathlib import Path


def sha(data):return hashlib.sha256(data).hexdigest()


def package(source,output):
    raw=Path(source).read_bytes();product=json.loads(raw)
    canonical=(json.dumps(product,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()
    compressed=gzip.compress(canonical,compresslevel=9,mtime=0)
    output=Path(output);output.parent.mkdir(parents=True,exist_ok=True);output.write_bytes(compressed)
    manifest=dict(schema='terluna.provisioning.compressed-product/1',product_schema=product['schema'],
        path=output.name,compression='gzip; mtime0',sha256=sha(compressed),decompressed_sha256=sha(canonical),
        producer=product['producer'],dependencies=product['dependencies'],
        packager=dict(path='provisioning/enclosed_industry/package.py',sha256=sha(Path(__file__).read_bytes())),
        reading_rule='Verify compressed hash, decompress, verify decompressed hash and parse the unchanged engineering-product schema. Compression changes storage only; every annual row is retained.')
    output.with_suffix('.manifest.json').write_text(json.dumps(manifest,sort_keys=True,indent=2)+'\n')
    return manifest


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('source',type=Path);ap.add_argument('output',type=Path);a=ap.parse_args();print(json.dumps(package(a.source,a.output),indent=2))

if __name__=='__main__':main()
