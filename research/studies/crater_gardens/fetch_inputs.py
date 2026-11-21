"""Restore the pinned LUC-GRAS200 archive and its documentation to ignored research storage."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from urllib.request import Request, urlopen

from research.studies.crater_gardens.count import DEFAULT_WORK, HERE, verified


def restore(work_dir=DEFAULT_WORK, download=False):
    manifest = json.loads((HERE / 'inputs.json').read_text())
    folder = Path(work_dir) / 'inputs'
    folder.mkdir(parents=True, exist_ok=True)
    for key, filename in (('catalogue', 'LUC-GRAS200_full.zip'), ('documentation', 'LUC-GRAS200_README.md')):
        record = manifest[key]
        target = folder / filename
        if target.exists():
            verified(target, record)
            continue
        if not download:
            raise FileNotFoundError(f'Restore {target} with --download')
        partial = target.with_name(target.name + '.partial')
        size = record['bytes']
        if partial.exists() and partial.stat().st_size > size:
            raise ValueError(f'Retained partial exceeds the recorded source size: {partial}')
        if key == 'documentation':
            with urlopen(record['url'], timeout=45) as response:
                data = response.read(size + 1)
            if len(data) != size:
                raise ValueError('Documentation response differs in length; inspect the intended source')
            partial.write_bytes(data)
        else:
            while (offset := partial.stat().st_size if partial.exists() else 0) < size:
                end = min(offset + 16 * 1024**2 - 1, size - 1)
                request = Request(record['url'], headers={'User-Agent': 'Mozilla/5.0',
                                  'Range': f'bytes={offset}-{end}', 'If-Range': record['etag']})
                with urlopen(request, timeout=60) as response:
                    if response.status != 206 or response.headers.get('Content-Range') != f'bytes {offset}-{end}/{size}':
                        raise ValueError('The source must return the recorded byte range')
                    data = response.read(end - offset + 2)
                if len(data) != end - offset + 1:
                    raise ValueError('Incomplete byte range; the completed prefix remains available for resume')
                with partial.open('ab') as stream:
                    stream.write(data)
                    stream.flush()
                    os.fsync(stream.fileno())
        verified(partial, record)
        partial.replace(target)
    return folder / 'LUC-GRAS200_full.zip'


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, default=DEFAULT_WORK)
    parser.add_argument('--download', action='store_true')
    args = parser.parse_args()
    print(restore(args.work_dir, args.download))
