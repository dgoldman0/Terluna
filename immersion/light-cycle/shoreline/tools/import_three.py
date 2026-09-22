"""Import only the verified r186 ESM modules and MIT license from a release ZIP.

Usage: python tools/import_three.py /path/to/three.js-r186.zip
Existing reference hashes are checked; a mismatch never updates the manifest.
"""
from pathlib import Path
import argparse, hashlib, json
from zipfile import ZipFile
ROOT=Path(__file__).resolve().parents[1]

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('archive',type=Path)
    args=parser.parse_args()
    out=ROOT/'vendor/three-r186'
    manifest=json.loads((out/'manifest.json').read_text())
    digest=hashlib.file_digest(args.archive.open('rb'),'sha256').hexdigest()
    if digest!=manifest['archive_sha256']:
        raise SystemExit('Archive hash mismatch. Reference input has been preserved.')
    pending={}
    with ZipFile(args.archive) as archive:
        prefix='three.js-r186/'
        meta=json.loads(archive.read(prefix+'package.json'))
        if meta.get('version')!='0.186.0':raise ValueError('Expected Three.js r186')
        for name,spec in manifest['files'].items():
            data=archive.read(prefix+'build/'+name)
            if len(data)!=spec['bytes'] or hashlib.sha256(data).hexdigest()!=spec['sha256']:
                raise ValueError('Release-file hash mismatch: '+name)
            pending[name]=data
        pending['LICENSE']=archive.read(prefix+'LICENSE')
    # All verification succeeds before any dependency file is replaced.
    out.mkdir(parents=True,exist_ok=True)
    for name,data in pending.items():(out/name).write_bytes(data)
    (out/'package.json').write_text('{"type":"module"}\n')
    print(f'Imported {len(pending)} files; {sum(map(len,pending.values())):,} bytes. Examples and fonts excluded.')
if __name__=='__main__':main()
