"""Package the static calendar from checked domain products; output is ignored."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
sys.path.insert(0,str(ROOT))
from shared.provenance import constants_changed
sys.path.insert(0,str(HERE))
import assets as explorer_assets  # noqa: E402


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_product(path,schema):
    if not path.exists():raise FileNotFoundError(f'Missing calendar data product: {path}')
    product=json.loads(path.read_text())
    if product.get('schema')!=schema:raise ValueError(f'Unsupported schema in {path}')
    for producer in (product.get('producer',{}),product.get('optical_producer',{})):
        for name,expected in producer.get('files',{}).items():
            if digest(ROOT/name)!=expected:raise ValueError(f'Stale calendar producer: {name}')
        changed=constants_changed(producer.get('constants'))
        if changed:raise ValueError(f'Stale calendar constants: {changed}')
        for name,expected in producer.get('external_inputs',{}).items():
            candidate=ROOT/name
            if candidate.exists() and digest(candidate)!=expected:
                raise ValueError(f'Calendar external input differs: {name}')
    for name,expected in product.get('inputs',{}).items():
        # Large transport archives are audit identities, not runtime build inputs.
        if '/' in name and digest(ROOT/name)!=expected:
            raise ValueError(f'Stale calendar input: {name}')
    return product


PAGES=['index.html','explorer.css','explorer.mjs','explorer-worker.mjs','sky-frame.mjs','sky-render.mjs',
       'globe.mjs','words.mjs','almanac.html','almanac.css','almanac.mjs','worker.mjs']


def build(destination=None,require_assets=False):
    destination=Path(destination) if destination else HERE/'dist'
    domain=ROOT/'illumination/calendar'
    products=[domain/'results'/f'{name}.json' for name in ['astronomy','transfer','validation']]
    for p in products:check_product(p,f'terluna.illumination.calendar-{p.stem}/1')
    destination.mkdir(parents=True,exist_ok=True)
    for name in PAGES:shutil.copy2(HERE/name,destination/name)
    (destination/'data').mkdir(exist_ok=True);(destination/'model').mkdir(exist_ok=True)
    for p in products:shutil.copy2(p,destination/'data'/p.name)
    shutil.copy2(domain/'evaluator.mjs',destination/'model/evaluator.mjs')
    source=[HERE/p for p in PAGES+['build.py','assets.py']]+products+[domain/'evaluator.mjs']
    metadata=dict(schema='terluna.visualization.light-calendar-build/1',
                  files={str(p.relative_to(ROOT)):digest(p) for p in source})
    try:
        record=explorer_assets.build(destination)
        metadata['explorer_assets']=record
    except explorer_assets.MissingSource as missing:
        # The almanac still builds; the explorer page needs these display assets.
        if require_assets:raise
        metadata['explorer_assets']=dict(missing=str(missing))
        print('Explorer assets not built; missing',missing)
    (destination/'build.json').write_text(json.dumps(metadata,indent=2)+'\n')
    print('Built static calendar:',destination)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path)
    parser.add_argument('--require-assets',action='store_true',help='fail when an explorer asset source is missing')
    args=parser.parse_args()
    build(args.out,args.require_assets)
