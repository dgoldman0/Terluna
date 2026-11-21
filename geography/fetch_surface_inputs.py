#!/usr/bin/env python3
"""Restore the inputs for the ground's rock (the LROC mare boundaries and Lunar Prospector's oxides); never accept changed bytes.

    python -m geography.fetch_surface_inputs            # report state
    python -m geography.fetch_surface_inputs --download # fetch missing files (about 19 MB)

They have a manifest of their own (surface_inputs.json) so that the topography and gravity inputs, which the atlas and
the drainage chain require whole, do not depend on them.
"""
from __future__ import annotations
from pathlib import Path
import sys

from shared.external_inputs import Inputs

INPUTS = Inputs(Path(__file__).resolve().parent / 'surface_inputs.json')
manifest, input_dir, path, restore = INPUTS.manifest, INPUTS.input_dir, INPUTS.path, INPUTS.restore

if __name__ == '__main__':
    sys.exit(INPUTS.main(description=__doc__.splitlines()[0]))
