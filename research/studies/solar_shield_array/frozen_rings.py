"""Frozen common orbits for the ring screen (step 2 of the narrowing; stage 3 of relative_orbits.md).

A radial-facing tile is pushed toward the Moon on the day side and away from it on the night side. Averaged over
an orbit, that drives the eccentricity vector at a steady rate perpendicular to the Sun line, while Earth's tide
turns the orbit's apse line against the Sun. In the frame that follows the Sun the push therefore balances at a
forced eccentricity along the Sun line, and a ring started on a circle traces a circle around it, which is the
swing ring_screen found. Starting a ring on that equilibrium (a Sun-tracking, heliotropic orbit) should hold its
sunward crossing in place.

The runner finds the equilibrium by iteration. Single tiles fly a year in DE440s point-mass gravity with the
filter's sail force on a radial-facing normal tilted 1.5 degrees, the bundle interior's mutual-shadow factor and
Moon and Earth eclipses (as ring_screen). The orbit-averaged eccentricity vector, in the frame of the Sun's
direction projected into each ring's plane, traces a circle; its centre is the forced eccentricity. Pass 1 starts
on circles; each later pass starts on the centre the previous pass found, with the sunward crossing kept at the
same radius. The forced eccentricity scales inversely with areal mass, so light tiles start from the heavier
tiles' centre scaled to their mass.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.frozen_rings
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
from protection.dynamics.ephemeris import DEFAULT_KERNEL
from protection.dynamics.optical import unit
from protection.dynamics.ring_bundle import tile_frames
from .ring_screen import SETTINGS as SCREEN, environment, visible_fraction

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE/'results/frozen_rings.json'
FILES = ['research/studies/solar_shield_array/frozen_rings.py', 'research/studies/solar_shield_array/ring_screen.py',
         'protection/dynamics/ring_bundle.py', 'protection/dynamics/cycling.py', 'protection/dynamics/ephemeris.py',
         'protection/dynamics/optical.py']
DESIGN = dict(days=365., radii_km=[15000., 19000., 20000.], tilts_deg=[0., 20.], planes=['ecliptic', 'moon_orbit'],
              areal_masses_kg_m2=[.0627, .026, .005], passes=3, sample_s=600., rtol=1e-10, atol=1e-3,
              max_step_s=1800., tile_tilt_deg=SCREEN['tile_tilt_deg'], interior_light=SCREEN['interior_light'])
DAY = K.JULIAN_DAY


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def geometry(env):
    """Sun direction, ecliptic pole and the Moon's orbit pole at the epoch."""
    sample = env.at(0.)
    sun = sample['positions']['sun']
    ecliptic = np.array([0., -np.sin(np.radians(K.EARTH_OBLIQUITY_DEG)), np.cos(np.radians(K.EARTH_OBLIQUITY_DEG))])
    moon_orbit = unit(np.cross(-sample['positions']['earth'], sample['moon_v']-sample['earth_v']))
    return sun, dict(ecliptic=ecliptic, moon_orbit=moon_orbit)


def start(sun, pole, radius_km, tilt_deg, ex, ey):
    """A retrograde ring through the sub-solar point at the given crossing radius with eccentricity vector
    (ex, ey) in the ring plane: ex along the Sun direction, ey along the direction of motion there."""
    s = unit(sun-np.dot(sun, pole)*pole)
    node = np.cross(pole, s)
    b = np.radians(tilt_deg)
    up = np.cos(b)*s+np.sin(b)*pole
    along = -node
    r = radius_km*1e3
    # Orbit with e = ex*up + ey*along, position r*up: from e = ((v^2 - mu/r) r - (r.v) v)/mu, choose v.
    # With v = vr*up + vt*along: e.up = (vt^2 r)/mu - 1 - ... solve for the tangential and radial speeds.
    mu = K.MOON_GM
    vt = np.sqrt(mu*(1+ex)/r)
    vr = -ey*mu/(r*vt)
    return np.r_[r*up, vr*up+vt*along], up, along


def eccentricity(states):
    r, v = states[..., :3], states[..., 3:]
    rn = np.linalg.norm(r, axis=-1)
    return ((np.sum(v*v, axis=-1)-K.MOON_GM/rn)[..., None]*r-np.sum(r*v, axis=-1)[..., None]*v)/K.MOON_GM


def fly(env, states, sigma):
    """One year of single tiles; returns sample times and states."""
    m = len(states)
    tilt = np.radians(DESIGN['tile_tilt_deg'])
    per_mass = 1/np.asarray(sigma)

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
        light = (pressure*DESIGN['interior_light']*visible_fraction(q, sample)*per_mass)[:, None]*n
        return np.hstack([s[:, 3:], env.gravity(t, q)+light]).ravel()
    T = DESIGN['days']*DAY
    sol = solve_ivp(rhs, (0., T), states.ravel(), method='DOP853', rtol=DESIGN['rtol'], atol=DESIGN['atol'],
                    max_step=DESIGN['max_step_s'], dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    times = np.arange(0., T, DESIGN['sample_s'])
    return times, sol.sol(times).T.reshape(len(times), m, 6), int(sol.nfev)


def measure(env, times, traj, case):
    """Orbit-averaged eccentricity vector in the Sun-following frame, its circle, and the sunward crossings."""
    pole = case['pole']
    suns = np.stack([env.at(t)['positions']['sun'] for t in times])
    h = unit(np.cross(traj[:, :3], traj[:, 3:]))
    s = unit(suns-np.sum(suns*h, axis=1)[:, None]*h)
    # h x s is the direction of motion where the ring crosses the Sun line, as in start().
    w = np.cross(h, s)
    e = eccentricity(traj)
    ex, ey = np.sum(e*s, axis=1), np.sum(e*w, axis=1)
    period = 2*np.pi*np.sqrt(np.mean(np.linalg.norm(traj[:, :3], axis=1))**3/K.MOON_GM)
    width = max(1, int(round(period/DESIGN['sample_s'])))
    kernel = np.ones(width)/width
    ax, ay = np.convolve(ex, kernel, 'valid'), np.convolve(ey, kernel, 'valid')
    # Least-squares circle through the averaged points: x^2 + y^2 = 2 cx x + 2 cy y + c.
    A = np.c_[2*ax, 2*ay, np.ones_like(ax)]
    (cx, cy, c), *_ = np.linalg.lstsq(A, ax**2+ay**2, rcond=None)
    radius = float(np.sqrt(max(c+cx**2+cy**2, 0.)))
    spread = float(np.max(np.hypot(ax-ax.mean(), ay-ay.mean())))
    # Sunward crossings, as ring_screen measures them, with heights above the ecliptic.
    ecl = case['ecliptic']
    su = unit(suns-np.sum(suns*ecl, axis=1)[:, None]*ecl)
    lon = np.arctan2(np.sum(np.cross(su, traj[:, :3])*ecl, axis=1), np.sum(su*traj[:, :3], axis=1))
    height = traj[:, :3]@ecl
    rad = np.linalg.norm(traj[:, :3], axis=1)
    idx = np.nonzero((np.abs(lon[:-1]) < np.pi/2) & (np.sign(lon[:-1]) != np.sign(lon[1:])))[0]
    frac = lon[idx]/(lon[idx]-lon[idx+1])
    cross_t = times[idx]+frac*(times[idx+1]-times[idx])
    cross_h = height[idx]+frac*(height[idx+1]-height[idx])
    cross_r = rad[idx]+frac*(rad[idx+1]-rad[idx])
    fit = np.polyfit(cross_t/DAY, cross_h, 1) if len(idx) > 2 else [np.nan, np.nan]
    speed = np.sqrt(K.MOON_GM/(case['radius_km']*1e3))
    return dict(centre=[float(cx), float(cy)], circle_radius=radius, averaged_spread=spread,
                e_mean=float(np.mean(np.hypot(ax, ay))), e_max=float(np.max(np.hypot(ex, ey))),
                periapsis_min_km=float(rad.min()/1e3), crossings=int(len(idx)),
                crossing_radius_km=[float(cross_r.min()/1e3), float(cross_r.max()/1e3)] if len(idx) else None,
                crossing_height_error_max_km=float(np.abs(cross_h-case['design_height_m']).max()/1e3) if len(idx) else None,
                height_drift_km_per_day=float(fit[0]/1e3),
                steering_m_s_per_day=float(speed*abs(fit[0])/(case['radius_km']*1e3)))


def main():
    wall, cpu0 = time.monotonic(), time.process_time()
    env = environment(DESIGN['days'])
    sun, poles = geometry(env)
    base = []
    for plane in DESIGN['planes']:
        for radius in DESIGN['radii_km']:
            for tilt in DESIGN['tilts_deg']:
                for sigma in DESIGN['areal_masses_kg_m2']:
                    up = start(sun, poles[plane], radius, tilt, 0., 0.)[1]
                    base.append(dict(plane=plane, radius_km=radius, tilt_deg=tilt, sigma_kg_m2=sigma,
                                     pole=poles[plane], ecliptic=poles['ecliptic'],
                                     design_height_m=float(np.dot(radius*1e3*up, poles['ecliptic']))))
    heaviest = max(DESIGN['areal_masses_kg_m2'])
    passes, centres, nfev = [], {}, 0
    for k in range(DESIGN['passes']):
        cases = []
        for case in base:
            key = (case['plane'], case['radius_km'], case['tilt_deg'], case['sigma_kg_m2'])
            if k == 0:
                # Lighter tiles start where the heaviest tiles' forced eccentricity, scaled to their mass, lies;
                # from a circle their swing would be too large to fly.
                e0 = [0., 0.] if case['sigma_kg_m2'] == heaviest else None
            else:
                e0 = centres.get(key)
            if e0 is None:
                ref = centres.get((case['plane'], case['radius_km'], case['tilt_deg'], heaviest))
                if ref is None:
                    continue
                e0 = [c*heaviest/case['sigma_kg_m2'] for c in ref]
            if np.hypot(*e0) >= .9:
                continue
            cases.append(dict(case, e0=list(map(float, e0))))
        states = np.array([start(sun, c['pole'], c['radius_km'], c['tilt_deg'], *c['e0'])[0] for c in cases])
        times, traj, used = fly(env, states, [c['sigma_kg_m2'] for c in cases])
        nfev += used
        rows = []
        for j, case in enumerate(cases):
            result = measure(env, times, traj[:, j], case)
            key = (case['plane'], case['radius_km'], case['tilt_deg'], case['sigma_kg_m2'])
            centres[key] = result['centre']
            rows.append(dict({k2: v for k2, v in case.items() if k2 not in ('pole', 'ecliptic', 'design_height_m')},
                             design_height_km=case['design_height_m']/1e3, **result))
        passes.append(rows)
        print(f'pass {k}: {len(rows)} rings')
        for r in rows:
            print(f"  {r['plane']:10s} {r['radius_km']:6.0f} km tilt {r['tilt_deg']:4.0f} {r['sigma_kg_m2']*1e3:5.1f} g/m2 "
                  f"e0 {np.hypot(*r['e0']):.3f} centre {np.hypot(*r['centre']):.3f} circle {r['circle_radius']:.3f} "
                  f"spread {r['averaged_spread']:.3f} periapsis {r['periapsis_min_km']:7.0f} km "
                  f"crossing r {r['crossing_radius_km']} h err {r['crossing_height_error_max_km']}")
    out = dict(schema='terluna.research.frozen-rings/1',
               producer=dict(files={f: digest(ROOT/f) for f in FILES},
                             constants=constants_used([f for f in FILES if f.endswith('.py')]),
                             inputs={'de440s.bsp': digest(DEFAULT_KERNEL)}),
               evidence=('Single tiles on retrograde rings about the Moon, one year from 2026-10-04 TDB in DE440s '
                         'point-mass gravity, with ring_screen\'s sail force (radial-facing normal tilted 1.5 degrees, '
                         'interior shadow factor, Moon and Earth eclipses) at each areal mass. The eccentricity vector '
                         'is averaged over one orbit and expressed in the frame of the Sun\'s direction projected into '
                         'the ring plane (x toward the Sun, y along the motion at the sub-solar point); the least-squares '
                         'circle through it gives the forced eccentricity. Each pass starts on the previous pass\'s '
                         'centre with the sunward crossing at the design radius.'),
               reading_rule=('A ring is frozen when its averaged eccentricity stays near its centre (small spread) and '
                             'its sunward crossing radius and height stay put. The plane\'s drift against the Sun, '
                             'which eccentricity cannot hold, remains in the height drift and its steering.'),
               design=DESIGN, passes=passes,
               resources=dict(cpu_s=time.process_time()-cpu0, wall_s=time.monotonic()-wall, nfev=nfev,
                              max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, numerical_threads=1))
    OUT.write_text(json.dumps(out, indent=1, default=float)+'\n')


if __name__ == '__main__':
    main()
