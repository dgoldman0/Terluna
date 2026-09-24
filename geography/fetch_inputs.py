#!/usr/bin/env python3
"""Restore the exact LOLA topography and GRAIL gravity inputs; never accept changed bytes.

    python -m geography.fetch_inputs            # report state
    python -m geography.fetch_inputs --download # fetch missing files (about 46 MB)
"""
from __future__ import annotations
from pathlib import Path
import sys

from shared.external_inputs import Inputs

INPUTS = Inputs(Path(__file__).resolve().parent / 'inputs.json')
manifest, input_dir, path, restore = INPUTS.manifest, INPUTS.input_dir, INPUTS.path, INPUTS.restore

if __name__ == '__main__':
    sys.exit(INPUTS.main(description=__doc__.splitlines()[0]))
