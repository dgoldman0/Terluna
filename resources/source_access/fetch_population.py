"""Restore the byte-pinned public population grid to ignored local storage."""
import argparse,hashlib,json,ssl,urllib.request
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def main():
 p=argparse.ArgumentParser();p.add_argument('--allow-unverified-tls',action='store_true',help='For the documented SWRI chain failure only; exact SHA256 remains mandatory.');args=p.parse_args()
 cfg=json.loads((HERE/'population_inputs.json').read_text());s=cfg['source_files']['neomod_grid'];path=ROOT/cfg['archive_directory']/s['filename']
 if path.exists() and hashlib.sha256(path.read_bytes()).hexdigest()==s['sha256']:print('Pinned input present:',path);return
 context=ssl._create_unverified_context() if args.allow_unverified_tls else None
 raw=urllib.request.urlopen(s['url'],timeout=120,context=context).read()
 if len(raw)!=s['bytes'] or hashlib.sha256(raw).hexdigest()!=s['sha256']:raise ValueError('Changed live model; retained snapshot is required')
 path.parent.mkdir(parents=True,exist_ok=True);path.write_bytes(raw);print('Restored:',path)
if __name__=='__main__':main()
