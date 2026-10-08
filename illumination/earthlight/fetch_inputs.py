#!/usr/bin/env python3
"""Restore the exact table of the Earth's spectral radiance at the Moon; never accept changed bytes.

    python -m illumination.earthlight.fetch_inputs            # report state
    python -m illumination.earthlight.fetch_inputs --download # fetch missing files

The file must match the size and SHA-256 recorded in inputs.json; a mismatch
stops the restore instead of substituting anything.
"""
from __future__ import annotations
from pathlib import Path
import sys

from shared.external_inputs import Inputs

INPUTS = Inputs(Path(__file__).resolve().parent / 'inputs.json')
manifest, input_dir, path, restore = INPUTS.manifest, INPUTS.input_dir, INPUTS.path, INPUTS.restore

if __name__ == '__main__':
    sys.exit(INPUTS.main(description=__doc__.splitlines()[0]))
