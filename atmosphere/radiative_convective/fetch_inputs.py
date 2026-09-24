#!/usr/bin/env python3
"""Restore the exact external spectroscopic inputs; never accept changed bytes.

    python -m atmosphere.radiative_convective.fetch_inputs            # report state
    python -m atmosphere.radiative_convective.fetch_inputs --download # fetch missing files
    python -m atmosphere.radiative_convective.fetch_inputs --source DIR

Each file must match the size and SHA-256 recorded in inputs.json. A mismatch
means the upstream data changed (for example a new HITRAN edition) and needs
review; the script stops instead of substituting anything.
"""
from __future__ import annotations
import argparse
import hashlib
import http.client
import json
from pathlib import Path
import sys
import tempfile
import urllib.request

ROOT = Path(__file__).resolve().parent
MANIFEST = ROOT / 'inputs.json'


def manifest() -> dict:
    return json.loads(MANIFEST.read_text(encoding='utf-8'))


def input_dir() -> Path:
    return ROOT / manifest()['input_directory']


def path(name: str) -> Path:
    """Path of a verified input; raises if it is missing (use --download)."""
    item = next((f for f in manifest()['files'] if f['name'] == name), None)
    if item is None:
        raise KeyError(f'{name} is not a recorded input')
    target = input_dir() / name
    if not target.is_file():
        raise FileNotFoundError(f'Missing input {target}; run fetch_inputs --download')
    return target


def _digest(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _download(url: str, limit: int) -> bytes:
    req = urllib.request.Request(url, headers={'User-Agent': 'Terluna-input-restorer/1'})
    chunks = []
    try:
        with urllib.request.urlopen(req, timeout=300) as response:
            while True:
                chunk = response.read(1 << 20)
                if not chunk:
                    break
                chunks.append(chunk)
                if sum(map(len, chunks)) > limit:
                    break
    except http.client.IncompleteRead as err:
        # hitran.org closes its connection abruptly after the last byte; the
        # size and hash checks below decide whether the transfer is complete.
        chunks.append(err.partial)
    return b''.join(chunks)


def restore(download: bool = False, source: Path | None = None) -> dict:
    info = manifest()
    target = input_dir()
    report = {'present': [], 'restored': [], 'missing': []}
    for item in info['files']:
        name = item['name']
        if Path(name).name != name:
            raise ValueError(f'Unsafe input name: {name}')
        dest = target / name
        if dest.is_file() and dest.stat().st_size == item['bytes'] and _digest(dest.read_bytes()) == item['sha256']:
            report['present'].append(name)
            continue
        if source is not None and (source / name).is_file():
            data = (source / name).read_bytes()
        elif download:
            data = _download(item['url'], item['bytes'])
        else:
            report['missing'].append(name)
            continue
        if len(data) != item['bytes'] or _digest(data) != item['sha256']:
            raise ValueError(f'Size or checksum mismatch for {name} ({len(data)} bytes); '
                             'upstream data changed and need review before use.')
        target.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=target, delete=False) as tmp:
            tmp.write(data)
            tmp_path = Path(tmp.name)
        tmp_path.replace(dest)
        report['restored'].append(name)
    report['status'] = 'PASS' if not report['missing'] else 'BLOCKED'
    report['scope'] = 'exact input bytes only; not scientific source admission'
    return report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--download', action='store_true', help='fetch missing files from their URLs')
    parser.add_argument('--source', type=Path, help='copy missing files from a local directory')
    args = parser.parse_args(argv)
    report = restore(download=args.download, source=args.source)
    print(json.dumps({k: (len(v) if isinstance(v, list) else v) for k, v in report.items()}, indent=2))
    for name in report['missing']:
        print('missing:', name)
    return 0 if report['status'] == 'PASS' else 1


if __name__ == '__main__':
    sys.exit(main())
