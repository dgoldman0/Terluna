"""Restore exact retained survey bytes; changing live responses are rejected."""
import json,shutil,urllib.request
from pathlib import Path
from . import neo,lunar
HERE=Path(__file__).resolve().parent

def main():
    c=json.loads((HERE/'access_inputs.json').read_text());folder=neo.ROOT/c['archive_directory'];folder.mkdir(parents=True,exist_ok=True)
    for source in c['external_files']:
        dest=folder/source['path']
        if dest.exists():raw=dest.read_bytes()
        elif source.get('archive') and Path(source['archive']).exists():raw=Path(source['archive']).read_bytes()
        elif source.get('url'):raw=urllib.request.urlopen(source['url'],timeout=60).read()
        else:raise FileNotFoundError('Restore historical archive manually: '+source['path'])
        if len(raw)!=source['bytes'] or neo.sha(raw)!=source['sha256']:raise ValueError('Source changed; restore the retained snapshot or explicitly revise the study: '+source['path'])
        if not dest.exists():dest.write_bytes(raw)
    lunar.load_inputs(download=True)
    print('Exact source bytes verified; raw files remain in ignored research/runs.')

if __name__=='__main__':main()
