"""Files kept byte for byte as the baseline of an earlier state, and how products built from them are checked.

When a file changes after products were built from it, those products go on recording the file's earlier digest
under its live path. Where the change corrects physics that the current standing of the work does not rest on, the
earlier file is kept verbatim beside the live one and listed with its hash in research/provenance.json, the products
built from it say so in their records, and their identity checks read the baseline whenever a product recorded the
baseline's digest.

    from shared.baselines import recorded_file
    assert sha256(recorded_file(path, expected).read_bytes()) == expected
"""
from __future__ import annotations

import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

# The live path, and the baselines kept verbatim for it (research/provenance.json), oldest first.
BASELINES = {
    # The rain solver before its vapour diffusivity took Massman's (1998) measured fit on 25 October 2026; its
    # Chapman-Enskog value stood 17.6% lower.
    'atmosphere/rain/drops.py': ('atmosphere/rain/baseline/drops_2026-10-24.py',),
}


def _sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def recorded_file(path, recorded):
    """The file whose bytes a recorded sha256 digest, in full or its leading characters, names: the file at path
    when its digest matches, or else a baseline of it whose digest does; otherwise the file at path."""
    path = str(path)
    live = ROOT / path
    if live.is_file() and _sha256(live).startswith(recorded):
        return live
    for baseline in BASELINES.get(path, ()):
        kept = ROOT / baseline
        if kept.is_file() and _sha256(kept).startswith(recorded):
            return kept
    return live


def rebased(identity, recorded):
    """A current producer identity ({'files': {path: digest}, 'inputs': ...}) with each baselined path's digest
    replaced by the recorded one when the recorded digest is that path's baseline: the identity a product built
    from the baseline should still match."""
    out = dict(identity)
    for group in ('files', 'inputs'):
        if not isinstance(identity.get(group), dict) or not isinstance(recorded.get(group), dict):
            continue
        files = dict(identity[group])
        for path, digest in recorded[group].items():
            if path in BASELINES and isinstance(digest, str) and recorded_file(path, digest) != ROOT / path:
                files[path] = digest
        out[group] = files
    return out
