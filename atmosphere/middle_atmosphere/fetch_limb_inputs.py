#!/usr/bin/env python3
"""Restore the exact X-ray and extreme-ultraviolet absorption tables of limb_heat.py; never accept changed bytes.

    python -m atmosphere.middle_atmosphere.fetch_limb_inputs            # report state
    python -m atmosphere.middle_atmosphere.fetch_limb_inputs --download # fetch missing files

Each file must match the size and SHA-256 recorded in limb_inputs.json; a mismatch
stops the restore instead of substituting anything. The files share the inputs folder
of fetch_inputs.py, whose own manifest other products pin by hash.
"""
from __future__ import annotations
from pathlib import Path
import sys

from shared.external_inputs import Inputs

INPUTS = Inputs(Path(__file__).resolve().parent / 'limb_inputs.json')
manifest, input_dir, path, restore = INPUTS.manifest, INPUTS.input_dir, INPUTS.path, INPUTS.restore

if __name__ == '__main__':
    sys.exit(INPUTS.main(description=__doc__.splitlines()[0]))
