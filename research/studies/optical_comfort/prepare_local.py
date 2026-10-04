"""Copy verified ignored inputs and link an isolated optical checkout to local runs.

python -m research.studies.optical_comfort.prepare_local \
    --source-checkout /path/to/existing/terluna \
    --storage /large-drive/terluna-research/optical-comfort-runs

The source checkout is read only. Existing destination bytes must agree.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import shutil

from climate.crm.cloud_columns import sha256
from .run import HERE, ROOT


def prepare(source, storage, check_inputs=False):
    source, storage = Path(source).resolve(), Path(storage).resolve()
    if source == ROOT or storage.is_relative_to(ROOT):
        raise ValueError('Use a separate source checkout and external bulk-data storage')
    admission = json.loads((HERE / 'local_inputs.json').read_text())
    records = []

    def copy(src, dst, expected=None, size=None):
        actual = sha256(src)
        if (expected is not None and actual != expected) or (size is not None and src.stat().st_size != size):
            raise ValueError(f'Source input differs from its manifest: {src}')
        dst.parent.mkdir(parents=True, exist_ok=True)
        if dst.exists():
            if sha256(dst) != actual:
                raise ValueError(f'Existing destination differs: {dst}')
        else:
            shutil.copy2(src, dst)
        if sha256(dst) != actual:
            raise ValueError(f'Copy verification failed: {dst}')
        records.append(dict(source=str(src), stored=str(dst), sha256=actual, bytes=dst.stat().st_size))

    def link(relative, target):
        p = ROOT / relative
        p.parent.mkdir(parents=True, exist_ok=True)
        if p.is_symlink():
            if p.resolve() != target.resolve():
                raise ValueError(f'Existing link differs: {p}')
        elif p.exists():
            raise ValueError(f'Existing destination requires review: {p}')
        else:
            p.symlink_to(target, target_is_directory=True)

    storage.mkdir(parents=True, exist_ok=True)
    for domain in ('atmosphere/radiative_convective', 'atmosphere/middle_atmosphere', 'protection'):
        manifest = json.loads((ROOT / domain / 'inputs.json').read_text())
        relative = Path(domain) / manifest['input_directory']
        destination = storage / 'inputs' / relative
        for item in manifest['files']:
            copy(source / relative / item['name'], destination / item['name'], item['sha256'], item['bytes'])
        link(relative, destination)
    relative = Path('illumination/sky/data')
    destination = storage / 'inputs' / relative
    for name, item in admission['native_atlases'].items():
        copy(source / relative / name, destination / name, item['sha256'], item['bytes'])
    link(relative, destination)
    for relative, pattern in (('climate/gcm/products', '*'), ('climate/crm/products', '*'),
                              ('atmosphere/column/products', '*.json')):
        destination = storage / 'inputs' / relative
        for p in sorted((source / relative).glob(pattern)):
            if p.is_file():
                copy(p, destination / p.name)
        link(relative, destination)
    for relative, name in (('climate/crm/runs', 'crm-runs'), ('climate/gcm/runs', 'gcm-runs'),
                           ('research/runs/waves', 'wave-runs')):
        target = storage.parent / name
        if not target.is_dir():
            raise FileNotFoundError(f'Missing local simulation store: {target}')
        link(relative, target)
    for relative in ('atmosphere/radiative_convective/cache', 'atmosphere/middle_atmosphere/cache', 'illumination/sky/work'):
        target = storage / 'cache' / relative
        target.mkdir(parents=True, exist_ok=True)
        link(relative, target)
    link('research/runs/optical_comfort', storage)
    if check_inputs:
        lock = 'immersion/package-lock.json'
        if sha256(ROOT / lock) != sha256(source / lock):
            raise ValueError('Source and optical dependency locks differ')
        target = storage / 'environment/node_modules'
        if not target.exists():
            shutil.copytree(source / 'immersion/node_modules', target, symlinks=True)
        link('immersion/node_modules', target)
    report = dict(schema='terluna.local-optical-inputs/1', source_checkout=str(source), destination_checkout=str(ROOT),
                  files=records, mode='Verified copies of ignored inputs; shared raw runs linked on the research drive.')
    (storage / 'local-setup.json').write_text(json.dumps(report, indent=2) + '\n')
    print(f'Verified {len(records)} local input files; storage: {storage}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-checkout', type=Path, required=True)
    parser.add_argument('--storage', type=Path, required=True)
    parser.add_argument('--include-check-inputs', action='store_true', help='Copy matching local Node dependencies for make check')
    args = parser.parse_args()
    prepare(args.source_checkout, args.storage, args.include_check_inputs)


if __name__ == '__main__':
    main()
