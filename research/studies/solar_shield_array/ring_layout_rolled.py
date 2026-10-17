"""The ring stack's layout under the rolled attitude law's push (research/array-industry, 13 October 2026).

The glow screen's answer rolls each ring's tiles 80 degrees about the line through their centres on both halves of
the orbit outside the service arc (attitude_schemes, rolled_80_deg). Along the kept orbits that law leaves 0.481 of
the radial-facing tile's push on the eccentricity, so the forced eccentricity, the strips' swing at perilune and with
them the stack's radii change. This runner lays the stack out again exactly as ring_layout does, with one change: the
calibration and sample flights' sail force is scaled by the law's eccentricity drive, flown as an effective areal
mass of 26 g/m2 over the drive in the flights alone, and each calibration ring starts on frozen_rings' centre scaled
by the drive. The tiles stay 26 g/m2 everywhere else, the photon steering's reach among them.

Departures from the as-built layout (results/ring_layout.json), which stays as it is:
- the whole sail force is scaled, the small out-of-plane push of the tile normal's 1.5-degree tilt included, while the
  rolled law's own out-of-plane forces, which give the attitude study's strip turn and slow plane drift (in this
  product beside the layout), are not flown;
- the photon steering's reach stays attitude_schemes' steering roll at 26 g/m2, as in ring_layout;
- the checkpoints go to research/runs/array_industry/ring_layout_rolled, the product to
  results/ring_layout_rolled.json.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.ring_layout_rolled          # about 25 CPU min
    OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.ring_layout_rolled --smoke  # eight days
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

from . import frozen_rings as fr
from . import ring_layout as rl

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE/'results/ring_layout_rolled.json'
RUN = ROOT/'research/runs/array_industry/ring_layout_rolled'
LAW = 'rolled_80_deg'
SIGMA = rl.DESIGN['sigma_kg_m2']


def drive():
    """The law's eccentricity drive over the as-built tile's, and its plane effects, from attitude_schemes."""
    scheme = json.loads(rl.ATTITUDE.read_text())['schemes'][LAW]
    return (float(scheme['eccentricity_drive']['over_as_built']),
            dict(strip_turn_deg_per_day=scheme['strip_turn_deg_per_day'],
                 crossing_height_rate_km_per_day=scheme['crossing_height_rate_km_per_day']))


def main(argv=None):
    share, plane = drive()

    centre = rl.frozen_centre

    def frozen_centre(frozen, radius_km, tilt_deg):
        """frozen_rings' 26 g/m2 centre, as ring_layout takes it, scaled by the law's drive."""
        return share*centre(frozen, radius_km, tilt_deg)

    def fly(env, states, sigma):
        """frozen_rings' flight with the sail force scaled by the law's drive."""
        return fr.fly(env, states, np.asarray(sigma, float)/share)

    rl.frozen_centre, rl.fly = frozen_centre, fly
    rl.RUN, rl.OUT = RUN, OUT
    rl.FILES = rl.FILES+['research/studies/solar_shield_array/ring_layout_rolled.py']
    rl.main(argv)
    smoke = '--smoke' in (argv or [])
    path = RUN.parent/'ring_layout_smoke'/'ring_layout_smoke.json' if smoke else OUT
    out = json.loads(path.read_text())
    out['variant'] = dict(
        law=LAW, eccentricity_drive_over_as_built=share, tile_areal_mass_kg_m2=SIGMA,
        flown_areal_mass_kg_m2=SIGMA/share, plane_effects_not_flown=plane,
        departures=('The calibration and sample flights scale frozen_rings\' whole sail force by the law\'s '
                    'eccentricity drive, flown as the tile\'s areal mass over the drive in the flights alone, and '
                    'start each calibration ring on frozen_rings\' centre scaled by it. The rolled law\'s own '
                    'out-of-plane forces (the strip turn and crossing-height drift above, from attitude_schemes along '
                    'the kept orbits) are not flown; the photon steering\'s reach is attitude_schemes\' steering roll '
                    'at the tile\'s 26 g/m2, as in ring_layout.'))
    out['evidence'] = ('The as-built layout\'s flights and layouts (ring_layout) with the sail force scaled by the '
                       f'{LAW} law\'s eccentricity drive, {share:.4f} of the radial-facing tile\'s. '+out['evidence'])
    path.write_text(json.dumps(out, indent=1, default=float)+'\n')


if __name__ == '__main__':
    main(sys.argv[1:])
