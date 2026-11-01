"""Map every /tmp path the snow records name to its copy on the project-space drive.

The snow work of 21-24 October left working material in the laptop's /tmp,
which a reboot clears.  On 24 October it was copied, with every file's SHA-256
checked, to the drive folder below.  The committed products keep their original
/tmp values, because other records bind those products by hash; this module
writes a separate record that gives, for each /tmp path, the copy's path and
identity.  Paths that were transient run folders and are absent from the copy
are listed as such.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
COPY = Path('/media/projectspace/terluna-research/atmosphere-snow/tmp-copy-2026-10-24')
LINK = ROOT / 'research/runs/snow_tmp_copy'   # ignored link to COPY/'tmp', like research/runs/snow_scale
PRODUCT = HERE / 'results/tmp_relocation.json'
PATTERN = re.compile(r'/tmp/[A-Za-z0-9_./-]*')
SCOPE = ('atmosphere/snow', 'atmosphere/tests/test_snow_scale_column.py')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def manifest():
    rows = {}
    for line in (COPY / 'manifest-sha256.txt').read_text().splitlines():
        digest, name = line.split(None, 1)
        rows[name.strip()] = digest
    return rows


def references():
    """Distinct /tmp paths in the committed snow records, with the files naming them."""
    out = subprocess.run(['git', 'grep', '-o', '-I', PATTERN.pattern.replace('\\', ''), '--', *SCOPE],
                         cwd=ROOT, capture_output=True, text=True, check=False).stdout
    found = {}
    for line in out.splitlines():
        name, _, path = line.partition(':')
        path = path.rstrip('.,)')
        if name == str(PRODUCT.relative_to(ROOT)) or name.endswith('tmp_relocation.py'):
            continue
        found.setdefault(path, set()).add(name)
    return {p: sorted(v) for p, v in sorted(found.items())}


def locate(path, files):
    relative = path[len('/tmp/'):].rstrip('/')
    target = COPY / 'tmp' / relative
    if relative in files:
        return dict(kind='file', copy=str(target), sha256=files[relative],
                    verified=sha(target) == files[relative])
    inside = [n for n in files if n.startswith(relative + '/')]
    if inside:
        return dict(kind='directory', copy=str(target), files=len(inside))
    if target.is_symlink():
        return dict(kind='link', copy=str(target), target=str(target.readlink()))
    return dict(kind='absent', note='Not in the 24 October copy; a transient run or check location.')


def run():
    files = manifest()
    entries = []
    for path, named_by in references().items():
        entries.append(dict(tmp_path=path, named_by=named_by, **locate(path, files)))
    kinds = {}
    for e in entries:
        kinds[e['kind']] = kinds.get(e['kind'], 0) + 1
    record = dict(schema='terluna.atmosphere.snow-tmp-relocation/1',
                  copied='2026-10-24', copy_root=str(COPY), copy_readme_sha256=sha(COPY / 'README.md'),
                  copy_manifest_sha256=sha(COPY / 'manifest-sha256.txt'), copy_files=len(files),
                  repository_link=str(LINK.relative_to(ROOT)) + ' -> ' + str(COPY / 'tmp'),
                  rule=('Committed products keep their original /tmp values; read each through this record. '
                        'New work keeps nothing it needs in /tmp.'),
                  counts=kinds, entries=entries, producer_sha256=sha(Path(__file__)))
    if any(e.get('verified') is False for e in entries):
        raise ValueError('A copied file no longer matches its manifest hash')
    PRODUCT.write_text(json.dumps(record, separators=(',', ':')) + '\n')
    return record


def main():
    argparse.ArgumentParser(description=__doc__).parse_args()
    record = run()
    print(json.dumps(dict(counts=record['counts'], entries=len(record['entries'])), indent=1))


if __name__ == '__main__':
    main()
