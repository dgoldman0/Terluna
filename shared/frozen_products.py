#!/usr/bin/env python3
"""Frozen products: large results of closed studies, kept out of every checkout's tracked files.

A frozen product is a generated result that only its own closed study reads and that no longer changes. Git holds its
bytes in history. From the commit that froze it, its path is ignored, and a checkout holds one of three things there:
a link to the one copy in a store outside the repository, a plain copy written from Git's history, or nothing.
`frozen_products.json` beside this file lists every frozen product with its size, its SHA-256 and the Git blob that
holds its bytes, so each of the three states can be told apart and the second can be made anywhere the history is.

    python3 shared/frozen_products.py                # each product: present, absent or differing
    python3 shared/frozen_products.py --link         # link the absent ones to the store's copies, checked first
    python3 shared/frozen_products.py --restore      # write the absent ones from Git's history as plain files
    python3 shared/frozen_products.py --fill-store   # give the store the copies it lacks, from here or from history
    python3 shared/frozen_products.py --require      # as the report, with exit 2 if any is absent

The store is the folder named by TERLUNA_FROZEN_STORE, or else the manifest's `store`. A product that differs from
its recorded digest is never overwritten, linked over or copied into the store: the report names it and exits 1.
Tests that open an absent frozen product are skipped with the restoring command in the reason (the repository's
conftest.py).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MANIFEST = Path(__file__).resolve().with_suffix('.json')
RESTORE = 'python3 shared/frozen_products.py --link, or --restore'


def manifest() -> dict:
    return json.loads(MANIFEST.read_text(encoding='utf-8'))


def entries() -> list[dict]:
    return manifest()['files']


def store() -> Path:
    return Path(os.environ.get('TERLUNA_FROZEN_STORE') or manifest()['store'])


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as fh:
        for block in iter(lambda: fh.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def state(entry: dict, root: Path = ROOT) -> str:
    """'present' with the recorded bytes, 'absent' (nothing there, or a link to nothing), or 'differs'."""
    path = root / entry['path']
    if not path.exists():
        return 'absent'
    if path.stat().st_size != entry['bytes'] or digest(path) != entry['sha256']:
        return 'differs'
    return 'present'


def absent_paths(root: Path = ROOT) -> set[str]:
    """The frozen products with nothing to read at their path, without hashing anything."""
    return {e['path'] for e in entries() if not (root / e['path']).exists()}


def absent_frozen_product(exc: BaseException | None, root: Path = ROOT) -> str | None:
    """The absent frozen product an exception stems from: a FileNotFoundError naming it, anywhere in its chain."""
    absent = None
    seen = set()
    while exc is not None and id(exc) not in seen:
        seen.add(id(exc))
        name = getattr(exc, 'filename', None)
        if isinstance(exc, FileNotFoundError) and name:
            absent = absent_paths(root) if absent is None else absent
            try:
                rel = os.path.relpath(os.path.abspath(os.fspath(name)), root)
            except ValueError:
                rel = None
            if rel in absent:
                return rel
        exc = exc.__cause__ or exc.__context__
    return None


def from_history(entry: dict) -> bytes:
    """The product's bytes from Git's object store, checked against its digest."""
    done = subprocess.run(['git', '-C', str(ROOT), 'cat-file', 'blob', entry['git_blob']], capture_output=True)
    if done.returncode:
        raise RuntimeError(f"{entry['path']}: blob {entry['git_blob']} is not in this repository's history "
                           '(a shallow clone lacks it)')
    if hashlib.sha256(done.stdout).hexdigest() != entry['sha256']:
        raise RuntimeError(f"{entry['path']}: blob {entry['git_blob']} does not have the recorded digest")
    return done.stdout


def _write(target: Path, data: bytes) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    partial = target.with_name(target.name + '.partial')
    partial.write_bytes(data)
    partial.replace(target)


def _clear_dead_link(path: Path) -> None:
    if path.is_symlink() and not path.exists():
        path.unlink()


def fill_store(target: Path | None = None) -> list[str]:
    """Copy into the store every product it lacks, from this checkout's copy or else from history."""
    target = store() if target is None else target
    made = []
    for entry in entries():
        copy = target / entry['path']
        if copy.exists():
            if copy.stat().st_size != entry['bytes'] or digest(copy) != entry['sha256']:
                raise RuntimeError(f'{copy} differs from the recorded digest; the store is left as it is')
            continue
        here = ROOT / entry['path']
        data = here.read_bytes() if state(entry) == 'present' else from_history(entry)
        if hashlib.sha256(data).hexdigest() != entry['sha256']:
            raise RuntimeError(f"{entry['path']}: the copy here changed while it was read")
        _write(copy, data)
        made.append(entry['path'])
    return made


def link(source: Path | None = None) -> list[str]:
    """Link every absent product to the store's copy, each copy checked against its digest first."""
    source = store() if source is None else source
    made = []
    for entry in entries():
        if state(entry) != 'absent':
            continue
        copy = source / entry['path']
        if not copy.is_file():
            raise RuntimeError(f'{copy} is not in the store; --fill-store writes it, or --restore copies it here')
        if copy.stat().st_size != entry['bytes'] or digest(copy) != entry['sha256']:
            raise RuntimeError(f'{copy} differs from the recorded digest; nothing is linked to it')
        here = ROOT / entry['path']
        _clear_dead_link(here)
        here.parent.mkdir(parents=True, exist_ok=True)
        here.symlink_to(copy)
        made.append(entry['path'])
    return made


def restore() -> list[str]:
    """Write every absent product from Git's history as a plain file."""
    made = []
    for entry in entries():
        if state(entry) != 'absent':
            continue
        here = ROOT / entry['path']
        _clear_dead_link(here)
        _write(here, from_history(entry))
        made.append(entry['path'])
    return made


def report() -> dict[str, list[str]]:
    found: dict[str, list[str]] = {'present': [], 'absent': [], 'differs': []}
    for entry in entries():
        found[state(entry)].append(entry['path'])
    return found


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    action = parser.add_mutually_exclusive_group()
    action.add_argument('--link', action='store_true', help="link the absent products to the store's copies")
    action.add_argument('--restore', action='store_true', help="write the absent products from Git's history")
    action.add_argument('--fill-store', action='store_true', help='copy into the store what it lacks')
    parser.add_argument('--require', action='store_true', help='exit 2 if any product is absent')
    args = parser.parse_args(argv)
    try:
        if args.fill_store:
            for path in fill_store():
                print(f'stored    {path}')
        if args.link:
            for path in link():
                print(f'linked    {path}')
        if args.restore:
            for path in restore():
                print(f'restored  {path}')
    except RuntimeError as error:
        print(f'frozen_products: {error}', file=sys.stderr)
        return 1
    found = report()
    size = {e['path']: e['bytes'] for e in entries()}
    for kind in ('differs', 'absent'):
        for path in found[kind]:
            print(f'{kind:8}  {path}')
    megabytes = {kind: sum(size[p] for p in paths) / 1e6 for kind, paths in found.items()}
    print(f"{len(found['present'])} present ({megabytes['present']:.0f} MB), {len(found['absent'])} absent "
          f"({megabytes['absent']:.0f} MB), {len(found['differs'])} differing; store: {store()}")
    if found['differs']:
        return 1
    if found['absent'] and args.require:
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
