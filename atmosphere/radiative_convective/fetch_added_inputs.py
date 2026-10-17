#!/usr/bin/env python3
"""Restore the spectroscopic inputs added after inputs.json was bound: CH4 and N2O lines and partition sums.

    python -m atmosphere.radiative_convective.fetch_added_inputs            # report state
    python -m atmosphere.radiative_convective.fetch_added_inputs --download # fetch missing files

Each file must match the size and SHA-256 recorded in inputs_added.json; a mismatch stops the script instead of
substituting anything. The files share inputs.json's input directory.
"""
from __future__ import annotations
from pathlib import Path
import sys

from shared.external_inputs import Inputs

INPUTS = Inputs(Path(__file__).resolve().parent / 'inputs_added.json')
manifest, input_dir, path, restore = INPUTS.manifest, INPUTS.input_dir, INPUTS.path, INPUTS.restore

if __name__ == '__main__':
    sys.exit(INPUTS.main(description=__doc__.splitlines()[0]))
