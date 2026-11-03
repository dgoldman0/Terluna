"""The ledger of the Terluna patches in an ExoPlaSim source folder: which version of each patch each file carries
(3 November 2026).

exoplasim_run.ensure_patched applies a patch to ExoPlaSim's Fortran only where its marker is missing, and the twilight
copy's sources are shared by every run, so a patch revised after a file took it never reaches that file: the run builds
from the earlier version while its key records the file's hash. The ledger, a JSON file beside the sources
(LEDGER_NAME), holds for each file and each patch's marker the fingerprint of the edits the file took (the patch's
version); ensure_patched records a patch as it applies it and refuses to build from a file that carries a patch's marker
with no record or with another version's.

A folder patched before the ledger existed is checked once against an untouched copy of ExoPlaSim's sources
(exoplasim_run.UNTOUCHED_SOURCE, the pip wheel's): each file that carries a marker must be exactly the untouched file
with the current version of every patch whose marker it carries, applied in ensure_patched's order (exoplasim_run's
patch_catalogue). `record` checks every file so and writes the ledger for those that match; `rebuild` writes a file
anew from the untouched copy with the current versions of the patches it carries.

    python -m climate.gcm.patch_ledger [--twilight] verify|record   # the installed sources, or the twilight copy's
    python -m climate.gcm.patch_ledger [--twilight] rebuild FILE     # e.g. rainmod.f90; removes the executables
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

LEDGER_NAME = 'terluna_patches.json'


def fingerprint(edits) -> str:
    """A patch's version: the hash of its edits, each old text and the new text that replaces it, in order."""
    return hashlib.sha256(json.dumps([[old, new] for old, new in edits]).encode()).hexdigest()[:16]


def read(src: Path) -> dict:
    """The ledger of the folder src ({} if it has none): file name -> {marker: fingerprint}."""
    path = Path(src) / LEDGER_NAME
    return json.loads(path.read_text()) if path.exists() else {}


def write(src: Path, ledger: dict):
    path = Path(src) / LEDGER_NAME
    tmp = path.with_name(path.name + '.part')
    tmp.write_text(json.dumps(ledger, indent=1, sort_keys=True) + '\n')
    tmp.replace(path)


def note(ledger: dict, name: str, marker: str, edits):
    """Record in the ledger that the file name took this version of the patch marked marker."""
    ledger.setdefault(name, {})[marker] = fingerprint(edits)


def check(ledger: dict, src: Path, name: str, marker: str, edits):
    """Refuse a file that carries the patch's marker unless the ledger holds the patch's current version for it."""
    path = Path(src) / name
    held = ledger.get(name, {}).get(marker)
    if held is None:
        raise RuntimeError(f"{path} carries the Terluna patch marked {marker!r} with no record of its version in "
                           f"{Path(src) / LEDGER_NAME}; check the folder against the untouched sources with "
                           "`python -m climate.gcm.patch_ledger record` (--twilight for the twilight copy), which "
                           'records every file it finds to be the untouched file with the current patches')
    if held != fingerprint(edits):
        raise RuntimeError(f"{path} carries another version of the Terluna patch marked {marker!r} (the ledger holds "
                           f"{held}, the runner's patch is {fingerprint(edits)}); rebuild the file from the untouched "
                           f"sources with `python -m climate.gcm.patch_ledger rebuild {name}` (--twilight for the "
                           'twilight copy) and run again')


def apply(text: str, edits, name: str) -> str:
    """One patch applied as ensure_patched applies it: each old text found exactly once and replaced."""
    for old, new in edits:
        if text.count(old) != 1:
            raise RuntimeError(f'ExoPlaSim {name} is not the version the Terluna patch expects near {old.strip()!r}')
        text = text.replace(old, new)
    return text


def carried(text: str, name: str, catalogue) -> list:
    """The catalogue's patches for the file name whose markers text carries, in the catalogue's order."""
    return [(marker, edits) for n, marker, edits in catalogue if n == name and marker in text]


def expected(untouched: str, text: str, name: str, catalogue) -> str:
    """The untouched file with the current version of every patch whose marker text carries, in catalogue order."""
    for marker, edits in carried(text, name, catalogue):
        untouched = apply(untouched, edits, name)
    return untouched


def verify(src: Path, untouched: Path, catalogue) -> dict:
    """For each file of the catalogue in src that carries a patch's marker: whether it is exactly the untouched file
    with the current versions of the patches it carries, with the markers it carries."""
    out = {}
    for name in sorted({n for n, _, _ in catalogue}):
        path = Path(src) / name
        if not path.exists():
            continue
        text = path.read_text()
        markers = [m for m, _ in carried(text, name, catalogue)]
        if not markers:
            continue
        base = Path(untouched) / name
        try:
            ok = base.exists() and expected(base.read_text(), text, name, catalogue) == text
            reason = None if ok else ('no untouched copy' if not base.exists() else 'differs from the untouched file '
                                      'with the current patches')
        except RuntimeError as err:
            ok, reason = False, str(err)
        out[name] = dict(ok=ok, markers=markers, reason=reason)
    return out


def record(src: Path, untouched: Path, catalogue) -> dict:
    """verify, then write the ledger for every file found to carry the current patches; returns verify's result."""
    result = verify(src, untouched, catalogue)
    ledger = read(src)
    for name, r in result.items():
        if r['ok']:
            ledger[name] = {marker: fingerprint(edits) for marker, edits in carried((Path(src) / name).read_text(), name,
                                                                                  catalogue)}
    write(src, ledger)
    return result


def rebuild(src: Path, untouched: Path, name: str, catalogue) -> str:
    """Write the file name anew in src: the untouched file with the current version of every patch it carries, with
    the ledger's record; returns the new text. The caller removes the executables built from the old text."""
    path = Path(src) / name
    text = expected((Path(untouched) / name).read_text(), path.read_text(), name, catalogue)
    tmp = path.with_name(path.name + '.part')
    tmp.write_text(text)
    tmp.replace(path)
    ledger = read(src)
    ledger[name] = {marker: fingerprint(edits) for marker, edits in carried(text, name, catalogue)}
    write(src, ledger)
    return text


def main(argv=None):
    import importlib.util
    from climate.gcm import exoplasim_run as er
    parser = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    parser.add_argument('--twilight', action='store_true',
                        help="the twilight copy's sources (default: the installed package's)")
    parser.add_argument('action', choices=('verify', 'record', 'rebuild'))
    parser.add_argument('file', nargs='?', help='for rebuild: the file, e.g. rainmod.f90')
    args = parser.parse_args(argv)
    if args.twilight:
        src = er.TWILIGHT_PACKAGE / 'exoplasim' / 'plasim' / 'src'
    else:
        src = Path(importlib.util.find_spec('exoplasim').origin).parent / 'plasim' / 'src'
    catalogue = er.patch_catalogue()
    if args.action == 'rebuild':
        if not args.file:
            parser.error('rebuild needs the file')
        rebuild(src, er.UNTOUCHED_SOURCE, args.file, catalogue)
        for exe in (src.parent / 'run').glob('most_plasim_*.x'):
            exe.unlink()
        print(f'{src / args.file} rebuilt; the executables in {src.parent / "run"} removed')
        return
    result = (record if args.action == 'record' else verify)(src, er.UNTOUCHED_SOURCE, catalogue)
    for name, r in result.items():
        print(f"{name:16s} {'current' if r['ok'] else 'NOT CURRENT'}  {len(r['markers'])} patches"
              + ('' if r['ok'] else f"  ({r['reason']})"))
    if args.action == 'record':
        print(f'ledger: {src / LEDGER_NAME}')


if __name__ == '__main__':
    main()
