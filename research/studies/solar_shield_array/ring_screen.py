"""Year-long screen of ring planes against the Sun (stage 3 of relative_orbits.md).

Full-disk coverage needs nested rings whose sunward crossings stay at fixed
heights across the protected disk. A ring tilted by beta about the line through
the Moon perpendicular to the Sun crosses the Sun meridian at a height of about
r sin(beta); the height holds while Earth's tide turns the ring's line of nodes
with the Sun. Earth's tide turns ring planes about the pole of the Moon's own
orbit, 5.1 degrees from the ecliptic, while the Sun stays in the ecliptic, so
rings are set up against either plane and every crossing height is measured
from the Sun-Moon axis. This screen follows single tiles on such rings for a year in
DE440s point-mass gravity with the filter's sail force on a radial-facing
normal tilted 1.5 degrees, a mutual-shadow factor for bundle interiors, and
finite-Sun Moon and Earth eclipses (overlap of angular disks). At every sunward
crossing it records the tile's radius and height, and from their drift the
steering that would hold the strip.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.ring_screen
"""
from __future__ import annotations

import hashlib
import json
import resource
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.cycling import CyclingEnvironment
from protection.dynamics.ephemeris import DEFAULT_KERNEL, Ephemeris
from protection.dynamics.optical import disk_overlap, unit
from protection.dynamics.ring_bundle import tile_frames
from .cycling_search import CENTRAL

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE/'results/ring_screen.json'
FILES = ['research/studies/solar_shield_array/ring_screen.py', 'protection/dynamics/ring_bundle.py',
         'protection/dynamics/cycling.py', 'protection/dynamics/ephemeris.py', 'protection/dynamics/optical.py']
RADII_KM = (15000., 16000., 17000., 18000., 19000., 20000., 21000.)
TILTS_DEG = (0., 5., 10., 15., 20.)
SETTINGS = dict(days=365., tile_tilt_deg=1.5, interior_light=0.77, mass_ratio=1.2539, sample_s=600.,
                rtol=1e-10, atol=1e-3, max_step_s=1800.)
DAY = K.JULIAN_DAY


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def environment(days):
    eph = Ephemeris(epoch='2026-10-04')
    samples = eph.sample(np.arange(-3600., days*DAY+7200., 1800.))
    eph.close()
    return CyclingEnvironment(samples, .05, CENTRAL, 'filter_only')


def visible_fraction(q, sample):
    """Fraction of the solar disk visible past the Moon and Earth from each point (largest single overlap)."""
    sun = sample['positions']['sun']-q
    ds = np.linalg.norm(sun, axis=1)
    rs = np.arcsin(np.clip(K.SUN_RADIUS/ds, 0, 1))
    blocked = np.zeros(len(q))
    for body, radius in ((-q, K.MOON_RADIUS), (sample['positions']['earth']-q, K.EARTH_RADIUS)):
        db = np.linalg.norm(body, axis=1)
        rb = np.arcsin(np.clip(radius/db, 0, 1))
        sep = np.arccos(np.clip(np.sum(unit(sun)*unit(body), axis=1), -1, 1))
        behind = db < ds
        blocked = np.maximum(blocked, np.where(behind, disk_overlap(rs, rb, sep)/(np.pi*rs**2), 0.))
    return 1-blocked


def initial_states(env):
    """Retrograde circular rings tilted about a node line perpendicular to the Sun, starting at the sub-solar point.

    The reference plane is the ecliptic or the Moon's orbit about Earth; tilts are
    measured from it, toward its pole, at the sub-solar point.
    """
    sample = env.at(0.)
    sun = sample['positions']['sun']
    ecliptic = np.array([0., -np.sin(np.radians(K.EARTH_OBLIQUITY_DEG)), np.cos(np.radians(K.EARTH_OBLIQUITY_DEG))])
    earth = sample['positions']['earth']
    moon_orbit = unit(np.cross(-earth, sample['moon_v']-sample['earth_v']))
    states, cases = [], []
    for plane, pole in (('ecliptic', ecliptic), ('moon_orbit', moon_orbit)):
        s = unit(sun-np.dot(sun, pole)*pole)
        node = np.cross(pole, s)
        for radius in RADII_KM:
            for tilt in TILTS_DEG:
                b = np.radians(tilt)
                up = np.cos(b)*s+np.sin(b)*pole
                # Retrograde about the pole: at the sub-solar point the motion runs along -node.
                v = -np.sqrt(K.MOON_GM/(radius*1e3))*node
                states.append(np.r_[radius*1e3*up, v])
                height = np.dot(radius*up, ecliptic)
                cases.append(dict(plane=plane, radius_km=radius, tilt_deg=tilt, design_height_km=float(height)))
    return np.array(states), cases, ecliptic


def main():
    wall, cpu0 = time.monotonic(), time.process_time()
    env = environment(SETTINGS['days'])
    states, cases, pole = initial_states(env)
    m = len(states)
    tilt = np.radians(SETTINGS['tile_tilt_deg'])
    per_mass = 1/(env.sigma*SETTINGS['mass_ratio'])
    area_scale = SETTINGS['interior_light']

    def rhs(t, y):
        s = y.reshape(m, 6)
        q = s[:, :3]
        sample = env.at(t)
        n = tile_frames(s, tilt)[:, :, 2]
        sun = sample['positions']['sun']-q
        u = unit(sun)
        cosine = np.sum(u*n, axis=1)
        flux = K.SOLAR_CONSTANT*(K.AU/np.linalg.norm(sun, axis=1))**2
        pressure = -2*env.central_fraction/K.SPEED_OF_LIGHT*flux*np.abs(cosine)*cosine
        light = (pressure*area_scale*visible_fraction(q, sample)*per_mass)[:, None]*n
        return np.hstack([s[:, 3:], env.gravity(t, q)+light]).ravel()
    T = SETTINGS['days']*DAY
    sol = solve_ivp(rhs, (0., T), states.ravel(), method='DOP853', rtol=SETTINGS['rtol'], atol=SETTINGS['atol'],
                    max_step=SETTINGS['max_step_s'],
                    dense_output=True)
    times = np.arange(0., T, SETTINGS['sample_s'])
    all_states = sol.sol(times).T.reshape(len(times), m, 6)
    suns = np.stack([env.at(t)['positions']['sun'] for t in times])
    su = unit(suns-np.sum(suns*pole, axis=1)[:, None]*pole)
    out_cases = []
    for k, case in enumerate(cases):
        traj = all_states[:, k]
        # Angle of the tile from the Sun meridian about the ecliptic pole, and its height above the ecliptic.
        lon = np.arctan2(np.sum(np.cross(su, traj[:, :3])*pole, axis=1), np.sum(su*traj[:, :3], axis=1))
        height = traj[:, :3]@pole
        radius = np.linalg.norm(traj[:, :3], axis=1)
        idx = np.nonzero((np.abs(lon[:-1]) < np.pi/2)&(np.sign(lon[:-1]) != np.sign(lon[1:])))[0]
        crossings = []
        for i in idx:
            f = lon[i]/(lon[i]-lon[i+1])
            crossings.append((times[i]+f*(times[i+1]-times[i]), height[i]+f*(height[i+1]-height[i]),
                              radius[i]+f*(radius[i+1]-radius[i])))
        c = np.array(crossings)
        days = c[:, 0]/DAY
        h0 = case['design_height_km']*1e3
        fit = np.polyfit(days, c[:, 1], 1) if len(c) > 2 else [np.nan, np.nan]
        speed = np.sqrt(K.MOON_GM/(case['radius_km']*1e3))
        out_cases.append(dict(case, crossings=len(c), height_drift_km_per_day=float(fit[0]/1e3),
                              height_range_km=[float(c[:, 1].min()/1e3), float(c[:, 1].max()/1e3)],
                              height_error_max_km=float(np.abs(c[:, 1]-h0).max()/1e3),
                              radius_range_km=[float(c[:, 2].min()/1e3), float(c[:, 2].max()/1e3)],
                              # Holding the height against a steady drift: rotate the plane by the drift
                              # rate over the radius, an out-of-plane velocity change of v times that angle.
                              steering_m_s_per_day=float(speed*abs(fit[0])/(case['radius_km']*1e3)),
                              _days=days.tolist(), _error_m=(c[:, 1]-h0).tolist(),
                              series=dict(day=days.round(3).tolist()[::10],
                                          height_km=(c[:, 1]/1e3).round(2).tolist()[::10],
                                          radius_km=(c[:, 2]/1e3).round(2).tolist()[::10])))
    # The whole strip pattern can shift with the season; coverage needs the strips' relative heights.
    grid = np.arange(2., SETTINGS['days']-2., 1.)
    error = {}
    for case in out_cases:
        days = np.array(case['_days'])
        error[(case['plane'], case['radius_km'], case['tilt_deg'])] = np.interp(grid, days, np.array(case['_error_m']))
    differential = dict(between_tilts=[], between_radii=[])
    for plane in ('ecliptic', 'moon_orbit'):
        for radius in RADII_KM:
            for a, b in zip(TILTS_DEG[:-1], TILTS_DEG[1:]):
                d = error[(plane, radius, b)]-error[(plane, radius, a)]
                differential['between_tilts'].append(dict(plane=plane, radius_km=radius, tilts_deg=[a, b],
                                                          max_km=float(np.abs(d).max()/1e3),
                                                          drift_km_per_day=float(np.polyfit(grid, d, 1)[0]/1e3)))
        for tilt in TILTS_DEG:
            for a, b in zip(RADII_KM[:-1], RADII_KM[1:]):
                d = error[(plane, b, tilt)]-error[(plane, a, tilt)]
                differential['between_radii'].append(dict(plane=plane, tilt_deg=tilt, radii_km=[a, b],
                                                          max_km=float(np.abs(d).max()/1e3),
                                                          drift_km_per_day=float(np.polyfit(grid, d, 1)[0]/1e3)))
    for case in out_cases:
        del case['_days'], case['_error_m']
    out = dict(schema='terluna.research.ring-screen/1',
               producer=dict(files={f: digest(ROOT/f) for f in FILES},
                             constants=constants_used([f for f in FILES if f.endswith('.py')]),
                             inputs={'de440s.bsp': digest(DEFAULT_KERNEL)}),
               evidence=('Single tiles on retrograde rings about the Moon, one year from 2026-10-04 TDB in DE440s '
                         'point-mass gravity, with the filter sail force on a radial-facing normal tilted 1.5 '
                         'degrees, a uniform mutual-shadow factor for bundle interiors and Moon/Earth eclipse '
                         'fractions from angular-disk overlap. Rings start tilted from the ecliptic or from the '
                         'Moon\'s orbit plane; heights are offsets from the Sun-Moon axis (the ecliptic) at each '
                         'crossing of the Sun meridian.'),
               reading_rule=('A ring holds its strip when its crossing height stays within the strip overlap. '
                             'steering_m_s_per_day is the out-of-plane rate that cancels a steady height drift; '
                             'the monthly oscillation around it needs margin or separate trims.'),
               settings=SETTINGS, cases=out_cases, differential=differential,
               resources=dict(cpu_s=time.process_time()-cpu0, wall_s=time.monotonic()-wall, nfev=int(sol.nfev),
                              max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, numerical_threads=1))
    OUT.write_text(json.dumps(out, indent=1)+'\n')
    for c in out_cases:
        print(f"{c['plane']:10s} r {c['radius_km']:7.0f} tilt {c['tilt_deg']:4.0f}  drift {c['height_drift_km_per_day']:8.3f} km/day  "
              f"height {c['height_range_km'][0]:8.1f}..{c['height_range_km'][1]:8.1f}  "
              f"radius {c['radius_range_km'][0]:8.1f}..{c['radius_range_km'][1]:8.1f}  steer {c['steering_m_s_per_day']:.3f} m/s/day")
    for row in differential['between_tilts']+differential['between_radii']:
        print('differential', row)
    print('cpu', round(out['resources']['cpu_s'], 1))


if __name__ == '__main__':
    main()
