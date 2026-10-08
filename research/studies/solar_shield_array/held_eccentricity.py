"""How much of each ring's eccentricity drive the stack's radius layout asks keeping to hold back.

ring_layout.py lays the stack out with keeping that holds close rings' eccentricity vectors together, which moves the
edge rings below their forced (Sun-tracking) eccentricity. In the averaged balance that sets the forced value, the
sail's push turns the eccentricity vector across the Sun line as fast as the apse line's motion against the Sun turns
it back. Holding a ring at e_held below e_forced therefore takes a steady control drive of about
(e_forced - e_held)/e_forced of the push's own drive. This runner rebuilds the coupled layouts from ring_layout.json's
calibration, checks that they reproduce its largest departures, and records that share ring by ring, for the trim's
duties in the planned array-industry branch and for the parked keeping study.

    python -m research.studies.solar_shield_array.held_eccentricity
"""
from __future__ import annotations

import hashlib
import json
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np

from . import ring_layout as rl

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE/'results/held_eccentricity.json'
FILES = ['research/studies/solar_shield_array/held_eccentricity.py', 'research/studies/solar_shield_array/ring_layout.py']
SCHEMA = 'terluna.research.held-eccentricity/1'
THRESHOLDS = (.05, .1, .15, .2, .25)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def calibration(product):
    c = product['calibration']
    return rl.Calibration(c['radii_km'], c['tilts_deg'], c['turning_deg_per_day'], c['forced_eccentricity'])


def shares(clearance_m):
    product = json.loads(rl.OUT.read_text())
    cal = calibration(product)
    heights, _ = rl.stack_heights(rl.DESIGN['middle_radius_km'])
    radii, kept, change, need, passes = rl.coupled_matched(cal, heights, clearance_m/1e3, 0.)
    tilts, _, forced, _ = rl.evaluate(cal, radii, heights, 0.)
    share = (forced-kept)/forced
    order = np.argsort(heights)
    worst = int(np.argmax(share))
    tenths = np.array_split(order, 10)
    return dict(clearance_m=clearance_m, rings=int(len(radii)),
                largest_departure=float(np.abs(kept-forced).max()),
                layout_largest_departure=product['layouts'][str(int(clearance_m))]['held']['matched_coupled']
                ['held_eccentricity']['largest_departure'],
                share=dict(max=float(share.max()), mean=float(share.mean()),
                           rings_over={f'{t:g}': int(np.sum(share > t)) for t in THRESHOLDS},
                           mean_by_height_tenth=[float(share[k].mean()) for k in tenths],
                           below_orbit_plane_max=float(share[heights < 0].max()),
                           above_orbit_plane_max=float(share[heights > 0].max())),
                worst=dict(height_km=float(heights[worst]), tilt_deg=float(tilts[worst]), radius_km=float(radii[worst]),
                           forced=float(forced[worst]), held=float(kept[worst])),
                forced_range=[float(forced.min()), float(forced.max())])


def main():
    clearances = [float(c) for c in rl.DESIGN['clearance_m']]
    with ProcessPoolExecutor(max_workers=2) as pool:
        rows = list(pool.map(shares, clearances))
    out = dict(schema=SCHEMA, producer=dict(files={f: digest(ROOT/f) for f in FILES}),
               inputs={'results/ring_layout.json': digest(rl.OUT)},
               evidence=('The coupled matched layouts of ring_layout.py rebuilt from its product\'s calibration; the '
                         'share is the averaged balance\'s estimate of the steady drive that holding takes, with no '
                         'flight of held rings.'),
               reading_rule=('share = (forced - held)/forced eccentricity, the steady part of a ring\'s eccentricity '
                             'drive that keeping has to cancel, as a fraction of the drive its sail gives; the trim '
                             'reaches 25% of the redirected band\'s push, which also holds attitude.'),
               by_clearance={f'{r["clearance_m"]:g}': r for r in rows})
    out['checks'] = dict(reproduces_layout=all(abs(r['largest_departure']-r['layout_largest_departure']) < 1e-5
                                               for r in rows))
    OUT.write_text(json.dumps(out, indent=1)+'\n')
    print(json.dumps({k: dict(share=v['share'], worst=v['worst']) for k, v in out['by_clearance'].items()}, indent=1))
    print(out['checks'])


if __name__ == '__main__':
    main()
