"""Plane motion of the ring screen's strip pattern (item 2 of the narrowing).

frozen_rings held the sunward crossings' radius on a Sun-tracking eccentricity and left the ring planes free: at
19,000-20,000 km the strip pattern still moved 2,400-3,400 km over a year, and at 15,000 km tilted rings lagged the
Sun by 23-28 km a day. This runner splits that motion into the part common to every strip, which oversizing the
pattern covers, and the part that differs between strips, which changes the overlap between neighbours and needs
steering.

Nine rings span the protected aperture's height at each radius from 15,000 to 22,000 km and each reference plane,
from near its bottom edge to near its top, started on their forced eccentricity (frozen_rings' 26 g/m2 centres at
the nearest radius it flew, interpolated in tilt, then one pass that starts on the centres the first pass found).
They fly a year as frozen_rings' tiles do: 26 g/m2, radial-facing, the interior shadow factor, Moon and Earth
eclipses. Each pass's measurements are checkpointed in research/runs/solar_shield_array/plane_motion. Each sunward passage is measured where the ring
crosses three vertical planes of the aperture, the Sun meridian and 80% of the aperture's half-width at that height
on either side: the height above the ecliptic and the distance from the Moon.

The elevation of each crossing seen from the Moon's centre depends only on the ring's plane; its change times the
starting distance is the plane's share of the strip's motion, and the rest is the radius's share, which the in-plane
keeping holds. A vertical shift and a rotation of the whole pattern in the aperture plane, fitted to the plane's
share of every crossing each day, are the common motion. A pattern of straight full-width strips covers a disk in any
rotation, so oversizing has to cover only the shift. What remains differs between strips: differentiated across the
rings' starting heights and scaled to the bundle's height step, it is the change in the overlap between neighbouring
strips; its steady drift and the oscillation about it, the stack's half-monthly breathing, are given apart. The
plane turns that would remove the residual from each ring, about the line of nodes (the crossing's height) and about
the Sun line (the strip's slope), are set against the plane steering that photon roll gives (attitude_schemes).

    OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.plane_motion
"""
from __future__ import annotations

import hashlib
import json
import resource
import time
from pathlib import Path

import numpy as np

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.ephemeris import DEFAULT_KERNEL
from protection.dynamics.optical import covering_radius, unit
from .frozen_rings import DESIGN as FROZEN, FILES as FROZEN_FILES, OUT as FROZEN_OUT, fly, geometry, measure, start
from .ring_screen import environment
from .run import SCENARIO

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE/'results/plane_motion.json'
ATTITUDE = HERE/'results/attitude_schemes.json'
RUN = ROOT/'research/runs/solar_shield_array/plane_motion'
FILES = ['research/studies/solar_shield_array/plane_motion.py', 'research/studies/solar_shield_array/run.py',
         'research/studies/solar_shield_array/scenario.json'] + FROZEN_FILES
DESIGN = dict(radii_km=[15000., 19000., 20000., 21000., 22000.], planes=['ecliptic', 'moon_orbit'],
              height_fractions=[-.95, -.75, -.5, -.25, 0., .25, .5, .75, .95], width_fraction=.8, sigma_kg_m2=.026,
              passes=2, grid_days=[2., 363.], height_step_m=9286.456, design_overlap_m=600., series_every_days=10)
DAY = K.JULIAN_DAY


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def tilt_for(sun, pole, ecliptic, radius_km, height_km):
    """Tilt (degrees) about the line through the Moon perpendicular to the Sun at which a ring of the given radius in
    the reference plane with this pole crosses the Sun meridian at the given height above the ecliptic."""
    s = unit(sun-np.dot(sun, pole)*pole)
    a, b = radius_km*np.dot(s, ecliptic), radius_km*np.dot(pole, ecliptic)
    return float(np.degrees(np.arcsin(np.clip(height_km/np.hypot(a, b), -1, 1))-np.arctan2(a, b)))


def passages(env, times, traj, ecliptic, widths):
    """Sunward crossings of each ring through the vertical planes x = -X, 0, +X of the aperture (x horizontal and
    perpendicular to the Sun in the ecliptic): times (s), heights above the ecliptic (m) and distances (m)."""
    suns = np.stack([env.at(t)['positions']['sun'] for t in times])
    s = unit(suns-(suns@ecliptic)[:, None]*ecliptic)
    w = np.cross(ecliptic, s)
    out = []
    for k, width in enumerate(widths):
        p = traj[:, k, :3]
        x, y, depth, dist = np.sum(p*w, axis=1), p@ecliptic, np.sum(p*s, axis=1), np.linalg.norm(p, axis=1)
        rows = []
        for target in (-width, 0., width):
            f = x-target
            idx = np.nonzero((depth[:-1] > 0) & (np.sign(f[:-1]) != np.sign(f[1:])))[0]
            frac = f[idx]/(f[idx]-f[idx+1])
            rows.append(dict(t=times[idx]+frac*(times[idx+1]-times[idx]), y=y[idx]+frac*(y[idx+1]-y[idx]),
                             r=dist[idx]+frac*(dist[idx+1]-dist[idx]), x=target))
        out.append(rows)
    return out


def split(heights, xs, y_plane, step):
    """Common motion (shift c0 and rotation theta per day) fitted to the plane's share of every crossing, the
    residual, and the overlap change between neighbouring strips per height step.

    heights: (rings,) starting heights (m); xs: (rings, 3) crossing abscissae (m); y_plane: (days, rings, 3) (m)."""
    days = y_plane.shape[0]
    A = np.c_[np.ones(xs.size), xs.ravel()]
    coef = np.linalg.lstsq(A, y_plane.reshape(days, -1).T, rcond=None)[0]
    residual = y_plane-(coef[0][:, None, None]+coef[1][:, None, None]*xs[None])
    order = np.argsort(heights)
    h = heights[order]
    d = residual[:, order]
    change = (d[:, 1:]-d[:, :-1])/(h[1:]-h[:-1])[None, :, None]*step
    return coef[0], coef[1], residual, change


def trend(days, series):
    """Least-squares slope per day along the first axis."""
    return np.polyfit(days, series.reshape(len(days), -1), 1)[0].reshape(series.shape[1:])


def main():
    wall, cpu0 = time.monotonic(), time.process_time()
    env = environment(FROZEN['days'])
    sun, poles = geometry(env)
    ecliptic = poles['ecliptic']
    frozen = json.loads(FROZEN_OUT.read_text())['passes'][-1]
    sigma = DESIGN['sigma_kg_m2']
    cases = []
    for plane in DESIGN['planes']:
        for radius in DESIGN['radii_km']:
            aperture = (covering_radius(radius*1e3, K.AU, SCENARIO['protected_radii']*K.MOON_RADIUS)
                        +SCENARIO['formation_margin_m'])
            # frozen_rings' 26 g/m2 centres at the nearest radius it flew.
            near = min({row['radius_km'] for row in frozen}, key=lambda x: abs(x-radius))
            ref = {row['tilt_deg']: np.array(row['centre']) for row in frozen
                   if row['plane'] == plane and row['radius_km'] == near and abs(row['sigma_kg_m2']-sigma) < 1e-9}
            for fraction in DESIGN['height_fractions']:
                height = fraction*aperture
                tilt = tilt_for(sun, poles[plane], ecliptic, radius, height/1e3)
                # Forced eccentricity interpolated in the size of the tilt between frozen_rings' 0 and 20 degrees.
                e0 = ref[0.]+(ref[20.]-ref[0.])*abs(tilt)/20.
                up = start(sun, poles[plane], radius, tilt, 0., 0.)[1]
                cases.append(dict(plane=plane, radius_km=radius, tilt_deg=tilt, sigma_kg_m2=sigma, pole=poles[plane],
                                  ecliptic=ecliptic, aperture_m=float(aperture), height_fraction=fraction,
                                  design_height_m=float(np.dot(radius*1e3*up, ecliptic)),
                                  width_m=float(DESIGN['width_fraction']*np.sqrt(max(aperture**2-height**2, 0.))),
                                  e0=[float(e0[0]), float(e0[1])]))
    nfev, centres = 0, []
    RUN.mkdir(parents=True, exist_ok=True)
    code = ''.join(digest(ROOT/f) for f in FROZEN_FILES)
    for k in range(DESIGN['passes']):
        key = hashlib.sha256(json.dumps([[c['plane'], c['radius_km'], c['tilt_deg'], c['e0']] for c in cases]).encode()+
                             json.dumps({n: FROZEN[n] for n in ('days', 'sample_s', 'rtol', 'atol', 'max_step_s',
                                                                'tile_tilt_deg', 'interior_light')}).encode()+
                             code.encode()).hexdigest()[:16]
        path = RUN/f'pass_{k}_{key}.json'
        if path.exists():
            saved = json.loads(path.read_text())
            rows, crossings = saved['rows'], saved['crossings']
        else:
            states = np.array([start(sun, c['pole'], c['radius_km'], c['tilt_deg'], *c['e0'])[0] for c in cases])
            times, traj, used = fly(env, states, [c['sigma_kg_m2'] for c in cases])
            nfev += used
            rows = [measure(env, times, traj[:, j], c) for j, c in enumerate(cases)]
            crossings = [[{n: (v.tolist() if isinstance(v, np.ndarray) else v) for n, v in row.items()} for row in ring]
                         for ring in passages(env, times, traj, ecliptic, [c['width_m'] for c in cases])]
            tmp = path.with_suffix('.tmp')
            tmp.write_text(json.dumps(dict(rows=rows, crossings=crossings, nfev=used), default=float))
            tmp.rename(path)
            del traj
        centres.append([r['centre'] for r in rows])
        print(f'pass {k}', flush=True)
        for c, r in zip(cases, rows):
            print(f"  {c['plane']:10s} {c['radius_km']:6.0f} km h {c['design_height_m']/1e3:8.0f} km tilt "
                  f"{c['tilt_deg']:6.2f} e0 {np.hypot(*c['e0']):.4f} centre {np.hypot(*r['centre']):.4f} "
                  f"circle {r['circle_radius']:.4f} drift {r['height_drift_km_per_day']:+.2f} km/d", flush=True)
        if k < DESIGN['passes']-1:
            for c, r in zip(cases, rows):
                c['e0'] = list(map(float, r['centre']))
    measured = rows
    grid = np.arange(DESIGN['grid_days'][0], DESIGN['grid_days'][1], 1.)
    attitude = json.loads(ATTITUDE.read_text()) if ATTITUDE.exists() else None
    results = []
    for plane in DESIGN['planes']:
        for radius in DESIGN['radii_km']:
            idx = [j for j, c in enumerate(cases) if c['plane'] == plane and c['radius_km'] == radius]
            heights = np.array([cases[j]['design_height_m'] for j in idx])
            y = np.zeros((len(grid), len(idx), 3)); r = np.zeros_like(y); xs = np.zeros((len(idx), 3))
            for a, j in enumerate(idx):
                for b, row in enumerate(crossings[j]):
                    y[:, a, b] = np.interp(grid, np.asarray(row['t'])/DAY, row['y'])
                    r[:, a, b] = np.interp(grid, np.asarray(row['t'])/DAY, row['r'])
                    xs[a, b] = row['x']
            elevation = np.arcsin(np.clip(y/r, -1, 1))
            y_plane = r[0][None]*(np.sin(elevation)-np.sin(elevation[0])[None])
            y_total = y-y[0][None]
            y_radius = y_total-y_plane
            c0, theta, residual, change = split(heights, xs, y_plane, DESIGN['height_step_m'])
            _, _, _, change_total = split(heights, xs, y_total, DESIGN['height_step_m'])
            _, _, _, change_radius = split(heights, xs, y_radius, DESIGN['height_step_m'])
            # The overlap's steady drift and, about it, the oscillation (the stack's half-monthly breathing).
            rate = trend(grid, change)
            fit = np.polyfit(grid, change.reshape(len(grid), -1), 1)
            swing = change.reshape(len(grid), -1)-(np.outer(grid, fit[0])+fit[1])
            # Plane turns that would remove each ring's residual: about the line of nodes (the crossing's height over
            # its distance) and about the Sun line (the strip's slope across the aperture).
            node_turn = residual[:, :, 1]/r[0][None, :, 1]
            sun_turn = (residual[:, :, 2]-residual[:, :, 0])/(xs[None, :, 2]-xs[None, :, 0])
            speed = np.sqrt(K.MOON_GM/(radius*1e3))
            node_rate, sun_rate = trend(grid, node_turn), trend(grid, sun_turn)
            dv = speed*(np.abs(node_rate)+np.abs(sun_rate))
            every = DESIGN['series_every_days']
            aperture = cases[idx[0]]['aperture_m']
            row = dict(plane=plane, radius_km=radius, aperture_km=aperture/1e3,
                       heights_km=(heights/1e3).round(1).tolist(), tilts_deg=[round(cases[j]['tilt_deg'], 3) for j in idx],
                       forced_eccentricity=[round(float(np.hypot(*measured[j]['centre'])), 4) for j in idx],
                       eccentricity_circle=[round(measured[j]['circle_radius'], 4) for j in idx],
                       meridian_drift_km_per_day=[round(measured[j]['height_drift_km_per_day'], 3) for j in idx],
                       common=dict(shift_km=dict(min=float(c0.min()/1e3), max=float(c0.max()/1e3)),
                                   rotation_deg=dict(min=float(np.degrees(theta.min())),
                                                     max=float(np.degrees(theta.max()))),
                                   oversizing_km_per_edge=float(np.abs(c0).max()/1e3),
                                   oversizing_share_of_rings=float(np.abs(c0).max()/aperture),
                                   centred_oversizing_km_per_edge=float(np.ptp(c0)/2e3),
                                   centred_oversizing_share_of_rings=float(np.ptp(c0)/2/aperture),
                                   series=dict(day=grid[::every].tolist(), shift_km=(c0[::every]/1e3).round(1).tolist(),
                                               rotation_deg=np.degrees(theta[::every]).round(3).tolist())),
                       plane_share_km=dict(max=float(np.abs(y_plane).max()/1e3)),
                       radius_share_km=dict(max=float(np.abs(y_radius).max()/1e3)),
                       differential=dict(
                           residual_km_max=float(np.abs(residual).max()/1e3),
                           overlap_change_per_step_m=dict(opening_max=float(change.max()), closing_max=float(-change.min()),
                                                          by_pair_and_place=change.max(axis=0).round(1).tolist()),
                           opening_rate_m_per_day=dict(max=float(rate.max()), by_pair_and_place=rate.round(2).tolist()),
                           closing_rate_m_per_day_max=float(-rate.min()),
                           oscillation_half_range_m=float(np.ptp(swing, axis=0).max()/2),
                           overlap_change_rate_m_per_day=float(np.abs(rate).max()),
                           days_to_consume_overlap=float(DESIGN['design_overlap_m']/max(rate.max(), 1e-12)),
                           total_with_radius_opening_m=float(change_total.max()),
                           radius_only_opening_m=float(change_radius.max()),
                           series=dict(day=grid[::every].tolist(),
                                       opening_m=change.max(axis=(1, 2))[::every].round(1).tolist())),
                       steering=dict(node_turn_deg_per_day=np.degrees(node_rate).round(5).tolist(),
                                     sun_line_turn_deg_per_day=np.degrees(sun_rate).round(5).tolist(),
                                     velocity_m_s_per_day_max=float(dv.max()),
                                     acceleration_m_s2_max=float(dv.max()/DAY)))
            if attitude is not None:
                # attitude_schemes' one-sided roll at its tiles' areal mass and radius; the plane's turn rate scales as
                # the force per unit mass over the orbital speed, the crossing height's with the radius as well.
                steer = attitude['schemes']['steering']
                base = attitude['tile_areal_mass_kg_m2']
                ratio = radius/attitude['radius_km']
                height = abs(steer['crossing_height_rate_km_per_day']['median'])*base/sigma*ratio**1.5
                turn = abs(steer['strip_turn_deg_per_day']['median'])*base/sigma*ratio**.5
                need_height = np.abs(node_rate)*radius*1e3/1e3
                need_turn = np.degrees(np.abs(sun_rate))
                row['photon_steering'] = dict(
                    crossing_height_km_per_day=float(height), strip_turn_deg_per_day=float(turn),
                    needed_km_per_day_max=float(need_height.max()), needed_strip_turn_deg_per_day_max=float(need_turn.max()),
                    needed_over_available_by_ring=np.maximum(need_height/height, need_turn/turn).round(2).tolist(),
                    note='attitude_schemes\' steering roll at 26 g/m2, scaled from its tiles and radius')
            results.append(row)
            print(json.dumps({k: v for k, v in row.items() if k in ('plane', 'radius_km', 'common')}, default=float)[:600],
                  flush=True)
    out = dict(schema='terluna.research.plane-motion/1',
               producer=dict(files={f: digest(ROOT/f) for f in FILES},
                             constants=constants_used([f for f in FILES if f.endswith('.py')]),
                             inputs={'results/frozen_rings.json': digest(FROZEN_OUT), 'de440s.bsp': digest(DEFAULT_KERNEL),
                                     **({'results/attitude_schemes.json': digest(ATTITUDE)} if attitude else {})}),
               evidence=('Single tiles on retrograde rings about the Moon, one year from 2026-10-04 TDB in DE440s '
                         'point-mass gravity with frozen_rings\' sail force at 26 g/m2, nine rings per radius and '
                         'reference plane from 95% below to 95% above the aperture, each started on its forced '
                         'eccentricity after one refinement pass. Crossings of the vertical planes at the Sun meridian '
                         'and at 80% of the aperture\'s half-width, interpolated between ten-minute samples and then '
                         'to whole days.'),
               reading_rule=('Oversizing covers the common shift on each edge. Where the overlap change between '
                             'neighbouring strips opens past the bundle\'s overlap, steering must hold the residual; '
                             'its plane turns per day set the out-of-plane velocity a day. The nine rings resolve '
                             'deformations smooth over about 1,800 km of height.'),
               design=DESIGN, cases=results,
               resources=dict(cpu_s=time.process_time()-cpu0, wall_s=time.monotonic()-wall, nfev=nfev,
                              max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, numerical_threads=1))
    OUT.write_text(json.dumps(out, indent=1, default=float)+'\n')


if __name__ == '__main__':
    main()
