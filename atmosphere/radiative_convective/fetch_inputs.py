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
from pathlib import Path
import sys

from shared.external_inputs import Inputs

INPUTS = Inputs(Path(__file__).resolve().parent / 'inputs.json')
manifest, input_dir, path, restore = INPUTS.manifest, INPUTS.input_dir, INPUTS.path, INPUTS.restore

if __name__ == '__main__':
    sys.exit(INPUTS.main(description=__doc__.splitlines()[0]))
