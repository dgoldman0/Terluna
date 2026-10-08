"""The Moon's night with the ring fleet beside it: the fleet's light against the night's own, place by place.

    python -m research.studies.joint_synthesis.night_sky
    # -> results/night_sky.json

The natural night comes from the light calendar's transfer product (illumination/calendar, sea-appearance branch):
the Sun's scattered light on level ground by solar elevation, with the Earth's light on the nearside taken at the
decision register's 2.3-2.9 lux under a high Earth. The fleet's light comes from illumination/fleet_light: diffuse
scatter at a stated reflectance, and mirror reflection for each tilt case. For each latitude the study follows the
night from dusk to dawn, counts the hours below practical dusk (2.98 lux) and below a tenth of a lux, with and
without the fleet, and finds the diffuse reflectance at which the fleet stays at a tenth of the night's own light at
local midnight.
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

from shared.constants import SYNODIC_MONTH_DAYS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.joint-synthesis-night-sky/1'
CALENDAR = ROOT / 'illumination' / 'calendar' / 'results' / 'transfer.json'
FLEET = ROOT / 'illumination' / 'fleet_light' / 'results' / 'night_light.json'
OUT = HERE / 'results' / 'night_sky.json'
DUSK_LUX = 2.98                       # practical dusk (research/decisions.md, Light and time)
DARK_LUX = 0.1
NEARSIDE_EARTHLIGHT_LUX = 2.3         # the register's darkest nearside nights under a high Earth (2.3-2.9 lux)
LATITUDES = (0.0, 15.0, 30.0, 45.0, 60.0)
RHO_D_LEVELS = (1e-3, 1e-4, 1e-5)
DEG_PER_HOUR = 360.0 / (SYNODIC_MONTH_DAYS * 24.0)
EVIDENCE = ('Clear-sky light from the calendar\'s solved column (molecular air, no refraction or haze) beside a '
            'geometric model of the lead ring fleet; the Moon\'s orbit plane stands for its equator, the nearside\'s '
            'Earthlight is held at a high Earth\'s darkest level, and the fleet\'s night-side attitude is a tilt case, '
            'not a flown law.')
READING_RULE = ('Lux on level ground. hours_below[threshold] count the night\'s hours (dusk to dawn, 354 h) below '
                'each level. rho_d_for_a_tenth is the diffuse reflectance of the fleet\'s Moon-facing faces that keeps '
                'its glow at a tenth of the night\'s own light at local midnight. resource gives the light the night '
                'half intercepts and mirrors, in lumens and as watts of sunlight.')


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def twilight(calendar):
    e = np.array(calendar['diffuse_elevation_deg'])
    lux = np.array(calendar['solar']['diffuse_lux'])
    return lambda elevation: np.interp(elevation, e, lux)


def fleet_map(fleet):
    alpha, beta = np.array(fleet['alpha_deg']), np.array(fleet['beta_deg'])
    lux = np.array(fleet['diffuse_lux'])
    def at(h_deg, lat):
        h = np.clip(np.abs(h_deg), alpha[0], alpha[-1])
        j = np.interp(lat, beta, np.arange(beta.size))
        j0 = int(np.clip(np.floor(j), 0, beta.size - 2)); t = j - j0
        col = (1 - t) * lux[:, j0] + t * lux[:, j0 + 1]
        return np.interp(h, alpha, col)
    return at


def specular_map(fleet, case):
    s = fleet['specular'][case]['lux_by_alpha_band']
    edges, lux = np.array(s['edges_deg']), np.array(s['lux'])
    return lambda angle: lux[np.clip(np.searchsorted(edges, angle) - 1, 0, lux.size - 1)]


def night(lat, sky, glow, rho_scale, nearside, specular=None):
    """Light through one night (hour angle from local midnight, -90 to 90 degrees), natural and with the fleet."""
    h = np.linspace(-90.0, 90.0, 721)
    angle = np.degrees(np.arccos(np.cos(np.radians(lat)) * np.cos(np.radians(h))))    # from the anti-solar point
    natural = sky(angle - 90.0) + (NEARSIDE_EARTHLIGHT_LUX if nearside else 0.0)
    fleet = glow(h, lat) * rho_scale + (specular(angle) if specular else 0.0)
    hours = lambda below: float(below.mean() * 180.0 / DEG_PER_HOUR)
    return dict(natural_midnight_lux=round(float(natural[360]), 4), fleet_midnight_lux=round(float(fleet[360]), 4),
                hours_below_dusk=dict(natural=round(hours(natural < DUSK_LUX), 1), with_fleet=round(hours(natural + fleet < DUSK_LUX), 1)),
                hours_below_dark=dict(natural=round(hours(natural < DARK_LUX), 1), with_fleet=round(hours(natural + fleet < DARK_LUX), 1)))


def main(argv=None) -> int:
    calendar, fleet = json.loads(CALENDAR.read_text()), json.loads(FLEET.read_text())
    sky, glow = twilight(calendar), fleet_map(fleet)
    rho0 = fleet['rho_d']
    places = {}
    for side in ('far side', 'nearside'):
        for lat in LATITUDES:
            row = {}
            for rho in RHO_D_LEVELS:
                row[f'diffuse_rho_{rho:g}'] = night(lat, sky, glow, rho / rho0, side == 'nearside')
            for case in ('random_1.5', 'random_3'):
                row[f'diffuse_rho_1e-05_with_specular_{case}'] = night(lat, sky, glow, 1e-5 / rho0, side == 'nearside',
                                                                       specular_map(fleet, case))
            natural0 = row[f'diffuse_rho_{RHO_D_LEVELS[0]:g}']['natural_midnight_lux']
            fleet0 = row[f'diffuse_rho_{RHO_D_LEVELS[0]:g}']['fleet_midnight_lux']
            row['rho_d_for_a_tenth'] = float(f'{rho0 * 0.1 * natural0 / fleet0:.2g}') if fleet0 > 0 else None
            places[f'{side}, {lat:g} deg'] = row
    lm_per_w = 94.0
    resource = dict(night_half_intercepted_lm=fleet['night_half_luminous_flux_lm'],
                    night_half_intercepted_w=fleet['night_half_luminous_flux_lm'] / lm_per_w,
                    mirrored_visible_lm=fleet['specular_reflected_lm'],
                    note='sunlight on the night half of the fleet, which the Moon never receives; the mirrored part is '
                         'the films\' visible reflection, and steering it is the same act as keeping it off the Moon')
    product = dict(schema=SCHEMA, producer=dict(study='joint_synthesis', files={'research/studies/joint_synthesis/night_sky.py': digest(__file__)},
                                                inputs={str(p.relative_to(ROOT)): digest(p) for p in (CALENDAR, FLEET)}),
                   evidence=EVIDENCE, reading_rule=READING_RULE, thresholds_lux=dict(dusk=DUSK_LUX, dark=DARK_LUX),
                   places=places, resource=resource,
                   specular_tilt_away_share={k: v['tiles_that_must_tilt_away_share'] for k, v in fleet['specular'].items() if k.startswith('random')})
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    for name, row in places.items():
        r = row['diffuse_rho_0.001']
        print(f"{name:16s} natural midnight {r['natural_midnight_lux']:8.4f} lux | fleet at 0.1%: {r['fleet_midnight_lux']:7.2f} lux | "
              f"hours below 0.1 lux {r['hours_below_dark']['natural']:5.1f} -> {r['hours_below_dark']['with_fleet']:5.1f} | "
              f"below dusk {r['hours_below_dusk']['natural']:5.1f} -> {r['hours_below_dusk']['with_fleet']:5.1f} | rho_d for a tenth {row['rho_d_for_a_tenth']}")
    print(f"night half intercepts {resource['night_half_intercepted_w']:.3g} W of sunlight")
    print(OUT)
    return 0


if __name__ == '__main__':
    sys.exit(main())
