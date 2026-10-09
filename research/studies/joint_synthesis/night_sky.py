"""The Moon's night with the ring fleet beside it: the fleet's light against the night's own, place by place.

    python -m research.studies.joint_synthesis.night_sky
    # -> results/night_sky.json

The natural night comes from the light calendar's transfer product (illumination/calendar): the Sun's scattered light
on level ground by solar elevation, and the Earth's light through the same air by the Earth's elevation and phase,
its disk integrated as the calendar's Python evaluator does (illumination/calendar/validate.py). The geometry is the
mean one: the Earth over the sub-Earth point at its mean distance and the Sun in the Moon's equatorial plane, with no
libration. The far-side rows lie on the 180° meridian and the nearside rows on the central meridian; places.py takes
each place at its own longitude. The fleet's light comes from illumination/fleet_light: diffuse scatter at a stated
reflectance, and mirror reflection for each tilt case. For each latitude the study follows the night from dusk to
dawn, counts the hours below practical dusk (2.98 lux) and below a tenth of a lux, with and without the fleet, and
finds the diffuse reflectance at which the fleet stays at a tenth of the night's own light at local midnight.
"""
from __future__ import annotations
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

from illumination.calendar.geometry import geometry_from_vectors
from illumination.calendar.validate import Response
from shared.constants import AU, EARTH_MOON_DISTANCE, SYNODIC_MONTH_DAYS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.joint-synthesis-night-sky/1'
CALENDAR = ROOT / 'illumination' / 'calendar' / 'results' / 'transfer.json'
FLEET = ROOT / 'illumination' / 'fleet_light' / 'results' / 'night_light.json'
CODE = tuple(ROOT / 'illumination' / 'calendar' / name for name in ('geometry.py', 'validate.py', 'disk.py'))
OUT = HERE / 'results' / 'night_sky.json'
DUSK_LUX = 2.98                       # practical dusk (research/decisions.md, Light and time)
DARK_LUX = 0.1
LATITUDES = (0.0, 15.0, 30.0, 45.0, 60.0)
MERIDIANS = {'far side': 180.0, 'nearside': 0.0}
HOUR_ANGLE = np.linspace(-90.0, 90.0, 721)   # the Sun's hour angle from local midnight, degrees; midnight in the middle
RHO_D_LEVELS = (1e-3, 1e-4, 1e-5)
DEG_PER_HOUR = 360.0 / (SYNODIC_MONTH_DAYS * 24.0)
EVIDENCE = ('Clear-sky light from the calendar\'s solved column (molecular air, no refraction or haze) beside a '
            'geometric model of the lead ring fleet. The geometry is the mean one: the Moon\'s orbit plane stands for '
            'its equator, the Earth stands over the sub-Earth point at its mean distance, and the libration that moves '
            'it by up to about 8 degrees is left out. The fleet\'s night-side attitude is a tilt case; no attitude law '
            'has been flown.')
READING_RULE = ('Lux on level ground. Far-side rows lie on the 180-degree meridian and nearside rows on the central '
                'meridian; earth_elevation_deg is the Earth\'s mean height there. natural_darkest_lux is the least of '
                'the night\'s own light. hours_below[threshold] count the night\'s hours (dusk to dawn, 354 h) below '
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


def earthlight(response, lat, lon):
    """The Earth's light on level ground through the night at a place, lux on HOUR_ANGLE, and the Earth's height."""
    n = HOUR_ANGLE.size
    sun_lon = np.radians(lon - 180.0 - HOUR_ANGLE)       # the subsolar point moves west as the night passes
    g = dict(earth=np.tile([1.0, 0.0, 0.0], (n, 1)), earth_distance_m=np.full(n, EARTH_MOON_DISTANCE),
             sun=np.stack([np.cos(sun_lon), np.sin(sun_lon), np.zeros(n)], axis=-1), sun_distance_m=np.full(n, AU))
    geo = geometry_from_vectors(g, lon, lat)
    lux = np.array([response.light({k: float(v[i]) for k, v in geo.items()})['earth'] for i in range(n)])
    return lux, float(geo['earth_elevation_deg'][0])


def natural(lat, sky, earth):
    """The night's own light on HOUR_ANGLE (twilight and the Earth's light), and the angle from the anti-solar point."""
    angle = np.degrees(np.arccos(np.cos(np.radians(lat)) * np.cos(np.radians(HOUR_ANGLE))))
    return sky(angle - 90.0) + earth, angle


def hours(below):
    return float(below.mean() * 180.0 / DEG_PER_HOUR)


def night(lat, sky, glow, rho_scale, earth, specular=None):
    """Light through one night (hour angle from local midnight, -90 to 90 degrees), natural and with the fleet."""
    own, angle = natural(lat, sky, earth)
    fleet = glow(HOUR_ANGLE, lat) * rho_scale + (specular(angle) if specular else 0.0)
    mid = HOUR_ANGLE.size // 2
    return dict(natural_midnight_lux=round(float(own[mid]), 4), natural_darkest_lux=round(float(own.min()), 4),
                fleet_midnight_lux=round(float(fleet[mid]), 4),
                hours_below_dusk=dict(natural=round(hours(own < DUSK_LUX), 1), with_fleet=round(hours(own + fleet < DUSK_LUX), 1)),
                hours_below_dark=dict(natural=round(hours(own < DARK_LUX), 1), with_fleet=round(hours(own + fleet < DARK_LUX), 1)))


def main(argv=None) -> int:
    calendar, fleet = json.loads(CALENDAR.read_text()), json.loads(FLEET.read_text())
    sky, glow, response = twilight(calendar), fleet_map(fleet), Response(calendar)
    rho0 = fleet['rho_d']
    places = {}
    for side, lon in MERIDIANS.items():
        for lat in LATITUDES:
            earth, elevation = earthlight(response, lat, lon)
            row = dict(longitude_deg=lon, earth_elevation_deg=round(elevation, 2))
            for rho in RHO_D_LEVELS:
                row[f'diffuse_rho_{rho:g}'] = night(lat, sky, glow, rho / rho0, earth)
            for case in ('random_1.5', 'random_3'):
                row[f'diffuse_rho_1e-05_with_specular_{case}'] = night(lat, sky, glow, 1e-5 / rho0, earth,
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
    files = {'research/studies/joint_synthesis/night_sky.py': digest(__file__)}
    files.update({str(p.relative_to(ROOT)): digest(p) for p in CODE})
    product = dict(schema=SCHEMA, producer=dict(study='joint_synthesis', files=files,
                                                inputs={str(p.relative_to(ROOT)): digest(p) for p in (CALENDAR, FLEET)}),
                   evidence=EVIDENCE, reading_rule=READING_RULE, thresholds_lux=dict(dusk=DUSK_LUX, dark=DARK_LUX),
                   places=places, resource=resource,
                   specular_tilt_away_share={k: v['tiles_that_must_tilt_away_share'] for k, v in fleet['specular'].items() if k.startswith('random')})
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    for name, row in places.items():
        r = row['diffuse_rho_0.001']
        print(f"{name:16s} Earth {row['earth_elevation_deg']:6.1f} deg | natural midnight {r['natural_midnight_lux']:8.4f} lux, "
              f"darkest {r['natural_darkest_lux']:8.4f} | fleet at 0.1%: {r['fleet_midnight_lux']:7.2f} lux | "
              f"hours below 0.1 lux {r['hours_below_dark']['natural']:5.1f} -> {r['hours_below_dark']['with_fleet']:5.1f} | "
              f"below dusk {r['hours_below_dusk']['natural']:5.1f} -> {r['hours_below_dusk']['with_fleet']:5.1f} | rho_d for a tenth {row['rho_d_for_a_tenth']}")
    print(f"night half intercepts {resource['night_half_intercepted_w']:.3g} W of sunlight")
    print(OUT)
    return 0


if __name__ == '__main__':
    sys.exit(main())
