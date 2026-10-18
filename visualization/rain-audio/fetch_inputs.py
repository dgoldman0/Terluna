#!/usr/bin/env python3
"""Restore the three CC0 recordings behind the texture models; never accept changed bytes.

    python3 visualization/rain-audio/fetch_inputs.py                     # report state
    python3 visualization/rain-audio/fetch_inputs.py --archive <zip>     # from the rain work's archive
    python3 visualization/rain-audio/fetch_inputs.py --source <folder>   # from a folder holding the files
    python3 visualization/rain-audio/fetch_inputs.py --download          # the Freesound file only; see inputs.json

Each file must match the size and SHA-256 recorded in inputs.json; a mismatch stops the restore instead of
substituting anything. The recordings are needed only to refit or compare the models.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import tempfile
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from shared.external_inputs import Inputs  # noqa: E402

INPUTS = Inputs(HERE / 'inputs.json')


def from_archive(archive: Path) -> dict:
    """Restore each file from its member of the archive, checking size and hash before writing any."""
    manifest = INPUTS.manifest()
    target = INPUTS.input_dir()
    pending, present = [], []
    with zipfile.ZipFile(archive) as source:
        for item in manifest['files']:
            path = target / item['name']
            if path.is_file() and hashlib.sha256(path.read_bytes()).hexdigest() == item['sha256']:
                present.append(item['name'])
                continue
            data = source.read(item['archive_member'])
            if len(data) != item['bytes'] or hashlib.sha256(data).hexdigest() != item['sha256']:
                raise ValueError(f'Size or checksum mismatch for {item["name"]} in {archive.name}')
            pending.append((path, data))
    target.mkdir(parents=True, exist_ok=True)
    for path, data in pending:
        with tempfile.NamedTemporaryFile(dir=target, delete=False) as tmp:
            tmp.write(data)
        Path(tmp.name).replace(path)
    return dict(status='PASS', present=present, restored=[p.name for p, _ in pending])


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--archive', type=Path, help='the rain work\'s archive, open-moon-rain-generator.zip')
    args, rest = parser.parse_known_args(argv)
    if args.archive is None:
        return INPUTS.main(rest, description=__doc__.splitlines()[0])
    try:
        print(json.dumps(from_archive(args.archive), indent=2))
        return 0
    except (OSError, KeyError, ValueError, zipfile.BadZipFile) as exc:
        print(f'INPUTS BLOCKED: {exc}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
