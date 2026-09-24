"""Restore exact external input files listed in a domain's inputs.json; never accept changed bytes.

A manifest lists each file's name, URL, size and SHA-256. Files live in the
manifest's input_directory beside it and are ignored by Git. A size or hash
mismatch means the upstream data changed (a new edition, a moved file) and
needs review; the restorer stops instead of substituting anything.
"""
from __future__ import annotations
import hashlib
import http.client
import json
from pathlib import Path
import tempfile
import urllib.request


class Inputs:
    def __init__(self, manifest_path: Path):
        self.manifest_path = Path(manifest_path)
        self.root = self.manifest_path.parent

    def manifest(self) -> dict:
        return json.loads(self.manifest_path.read_text(encoding='utf-8'))

    def input_dir(self) -> Path:
        return self.root / self.manifest()['input_directory']

    def path(self, name: str) -> Path:
        """Path of a recorded input; raises if it is missing."""
        if not any(f['name'] == name for f in self.manifest()['files']):
            raise KeyError(f'{name} is not a recorded input')
        target = self.input_dir() / name
        if not target.is_file():
            raise FileNotFoundError(f'Missing input {target}; run fetch_inputs --download')
        return target

    def restore(self, download: bool = False, source: Path | None = None) -> dict:
        target = self.input_dir()
        report = {'present': [], 'restored': [], 'missing': []}
        for item in self.manifest()['files']:
            name = item['name']
            if Path(name).name != name:
                raise ValueError(f'Unsafe input name: {name}')
            dest = target / name
            if dest.is_file() and dest.stat().st_size == item['bytes'] and _digest(dest.read_bytes()) == item['sha256']:
                report['present'].append(name)
                continue
            if source is not None and (Path(source) / name).is_file():
                data = (Path(source) / name).read_bytes()
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

    def main(self, argv=None, description='Restore external inputs') -> int:
        import argparse
        parser = argparse.ArgumentParser(description=description)
        parser.add_argument('--download', action='store_true', help='fetch missing files from their URLs')
        parser.add_argument('--source', type=Path, help='copy missing files from a local directory')
        args = parser.parse_args(argv)
        report = self.restore(download=args.download, source=args.source)
        print(json.dumps({k: (len(v) if isinstance(v, list) else v) for k, v in report.items()}, indent=2))
        for name in report['missing']:
            print('missing:', name)
        return 0 if report['status'] == 'PASS' else 1


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
        # Some servers (hitran.org) close abruptly after the last byte; the size
        # and hash checks decide whether the transfer is complete.
        chunks.append(err.partial)
    return b''.join(chunks)
