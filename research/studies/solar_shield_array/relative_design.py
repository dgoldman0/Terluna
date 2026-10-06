"""Linear relative-orbit screen of formation layouts for the shield (stage 2 of relative_orbits.md).

Tiles are parallel 10 km squares in one common attitude, as in the coupled
model. Positions follow the linear relative motion about a circular reference
orbit; u = 0 is the sub-solar point. Each layout is screened over one full
orbit for the smallest surface distance between squares, and over the service
arc for overlap between neighbouring rows and along each row.

Layouts:
- equal_energy_patch: rows alternate in depth through the relative eccentricity
  (no drift), with small fold steps, as a finite 2D patch;
- drifting_patch: rows alternate in depth through the semi-major axis;
- ring_bundle: each row is one orbit (a string of tiles at equal energy), rows
  step outward in radius, i.e. a bundle of adjacent rings that never reassemble.
The radial-facing attitude is tilted by a small angle about the orbit normal
so that neighbours on one orbit overlap like shingles; a Sun-facing patch is
screened for comparison. Ring bundles are then replayed with the coupled
model's ephemeris gravity. Each ring starts on a circle about the Moon, either
as two-body circular states or as time-shifted copies of its centre
trajectory, with one common attitude or each tile's own attitude.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.relative_design
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

from shared import constants as K
from protection.dynamics.relative_orbit import hcw_state, rtn_frame
from protection.dynamics.square_distance import closest_squares
from shared.provenance import constants_used
from .relative_diagnosis import PACKAGE, ROOT, SEED, digest, environment, propagate

HERE = Path(__file__).resolve().parent
OUT = HERE/'results/relative_design.json'
FILES = ['research/studies/solar_shield_array/relative_design.py', 'research/studies/solar_shield_array/relative_diagnosis.py',
         'protection/dynamics/relative_orbit.py', 'protection/dynamics/square_distance.py',
         'protection/dynamics/cycling.py', 'protection/dynamics/ephemeris.py']
SIDE, CLEAR, CLEARANCE = 10e3, 9.89e3, 150.
A_REF = 15.0e6
N_REF = np.sqrt(K.MOON_GM/A_REF**3)
SERVICE = np.radians([-22.5, 25.])


def frame(tilt_deg):
    """Common tile axes (columns u, v, n) in RTN: radial-facing, tilted about N."""
    t = np.radians(tilt_deg)
    eu, ev, en = np.array([np.sin(t), np.cos(t), 0.]), np.array([0., 0., 1.]), np.array([np.cos(t), -np.sin(t), 0.])
    return np.c_[eu, ev, en]


def layout(kind, rows=19, cols=21, pitch=8.5e3, depth=5e3, fold_step=0.2e3, tilt_deg=2.):
    """Relative orbital elements (divided by A_REF) at u = 0 for every tile, and the common frame."""
    r, c = np.meshgrid(np.arange(rows)-(rows-1)/2, np.arange(cols)-(cols-1)/2, indexing='ij')
    r, c = r.ravel(), c.ravel()
    e = np.zeros((len(r), 6))
    # Cross-track rows peak at mid-service and return to the nominal pitch at the service ends.
    amplitude = r*pitch/np.cos(np.mean(np.abs(SERVICE)))
    e[:, 5] = -amplitude/A_REF
    e[:, 1] = c*pitch/A_REF
    parity = np.where(np.arange(rows)[np.round(r+(rows-1)/2).astype(int)] % 2 == 0, 1., -1.)
    if kind == 'equal_energy_patch':
        e[:, 2] = -parity*depth/A_REF
        e[:, 3] = r*fold_step/A_REF
        e[:, 1] += 2*r*fold_step/A_REF
    elif kind == 'drifting_patch':
        e[:, 0] = parity*depth/A_REF
    elif kind == 'ring_bundle':
        e[:, 0] = (r-r.min())*depth/A_REF
    else:
        raise ValueError(kind)
    return e, frame(tilt_deg), r, c


def sun_facing(axes, u):
    """Rotate a common frame about the orbit normal so that it keeps facing the Sun."""
    c, s_ = np.cos(-u), np.sin(-u)
    rot = np.array([[c, -s_, 0.], [s_, c, 0.], [0., 0., 1.]])
    return rot@axes


def positions(elements, u):
    return np.stack([hcw_state(elements, A_REF, N_REF, 0., uu/N_REF)[:, :3] for uu in u])


def clearance(pos, axes, u=None, attitude='radial'):
    if attitude == 'sun':
        return np.array([closest_squares(p, sun_facing(axes, uu), SIDE)[0] for p, uu in zip(pos, u)])
    return np.array([closest_squares(p, axes, SIDE)[0] for p in pos])


def service_overlap(pos, axes, u, r, c):
    """Smallest overlap in projection along the Sun between neighbouring rows and along rows."""
    rows_min, along_min = np.inf, np.inf
    for p, uu in zip(pos, u):
        sun = np.array([np.cos(uu), -np.sin(uu), 0.])
        e1 = np.array([np.sin(uu), np.cos(uu), 0.])
        e2 = np.array([0., 0., 1.])
        # Projected clear half-widths of a square along e1 and e2.
        half1 = CLEAR/2*(abs(axes[:, 0]@e1)+abs(axes[:, 1]@e1))
        half2 = CLEAR/2*(abs(axes[:, 0]@e2)+abs(axes[:, 1]@e2))
        x1, x2 = p@e1, p@e2
        for k in np.unique(r)[:-1]:
            lower, upper = r == k, r == k+1
            # Overlap of the cross-track bands of adjacent rows over their common along-track span.
            gap = (x2[upper].min()-half2)-(x2[lower].max()+half2)
            rows_min = min(rows_min, -gap)
        for k in np.unique(r):
            row = np.argsort(x1[r == k])
            xs = x1[r == k][row]
            along_min = min(along_min, float(np.min(2*half1-np.diff(xs))))
    return float(rows_min), float(along_min)


def screen(kind, orbits=1, attitude='radial', **kw):
    e, axes, r, c = layout(kind, **kw)
    u = np.radians(np.linspace(-180., -180.+360.*orbits, 720*orbits+1))
    pos = positions(e, u)
    dist = clearance(pos, axes, u, attitude)
    s = np.radians(np.linspace(np.degrees(SERVICE[0]), np.degrees(SERVICE[1]), 49))
    rows_overlap, along_overlap = service_overlap(positions(e, s), axes, s, r, c)
    k = int(np.argmin(dist))
    drift = 3*np.pi*A_REF*np.ptp(e[:, 0])
    return dict(kind=kind, attitude=attitude, orbits=orbits, params=kw, tiles=len(e), min_surface_distance_m=float(dist[k]),
                at_u_deg=float(np.degrees(u[k])), samples_below_clearance=int(np.count_nonzero(dist < CLEARANCE)),
                service_min_surface_distance_m=float(clearance(positions(e, s), axes, s, attitude).min()),
                adjacent_row_overlap_min_m=rows_overlap, along_row_overlap_min_m=along_overlap,
                semi_major_axis_spread_m=float(A_REF*np.ptp(e[:, 0])), drift_spread_per_orbit_m=float(drift))


def tile_axes(states, tilt):
    """Per-tile radial-facing axes from each tile's own velocity, tilted about the orbit normal."""
    v = states[:, 3:6]/np.linalg.norm(states[:, 3:6], axis=1)[:, None]
    n = np.cross(states[:, :3], states[:, 3:6])
    n /= np.linalg.norm(n, axis=1)[:, None]
    r = np.cross(v, n)
    return np.stack([np.sin(tilt)*r+np.cos(tilt)*v, n, np.cos(tilt)*r-np.sin(tilt)*v], axis=-1)


def own_attitude_distance(states, tilt):
    """Smallest surface distance between nearly parallel squares, each pair in its mean frame."""
    from scipy.spatial import cKDTree
    q = states[:, :3]
    axes = tile_axes(states, tilt)
    pairs = cKDTree(q).query_pairs(np.sqrt(2)*SIDE+500., output_type='ndarray')
    best = np.inf
    for chunk in np.array_split(pairs, max(1, len(pairs)//4000+1)):
        f = axes[chunk[:, 0]]+axes[chunk[:, 1]]
        f /= np.linalg.norm(f, axis=1)[:, None, :]
        d = np.einsum('ki,kij->kj', q[chunk[:, 0]]-q[chunk[:, 1]], f)
        gap = np.linalg.norm(np.maximum(np.abs(d)-[SIDE, SIDE, 0.], 0.), axis=1)
        best = min(best, float(gap.min()))
    return best


def ephemeris_ring_bundle(rings=19, per_ring=21, pitch=8.5e3, step=450., tilt_deg=1.5, hours=93., time_shift=True,
                          attitude='common'):
    """Exact circular rings about the seed's epoch centroid, replayed with ephemeris gravity."""
    env = environment()
    seed = np.array(json.loads(SEED.read_text())['arrays']['initial_epoch_state'])
    ref = seed.mean(axis=0)
    R, T, N = rtn_frame(ref)
    sun = env.at(0.)['positions']['sun']
    sp = sun-np.dot(sun, N)*N
    sp /= np.linalg.norm(sp)
    w = np.cross(N, sp)
    u0 = np.arctan2(np.dot(R, w), np.dot(R, sp))
    r0 = np.linalg.norm(ref[:3])
    half = np.mean(np.abs(SERVICE))
    states, ring_of = [], []
    for k in range(rings):
        radius = r0+(k-(rings-1)/2)*step
        height = (k-(rings-1)/2)*pitch/np.cos(half)
        b = np.arcsin(height/radius)
        s_k = np.cos(b)*sp+np.sin(b)*N
        speed = np.sqrt(K.MOON_GM/radius)
        centre = np.r_[radius*(np.cos(u0)*s_k+np.sin(u0)*w), speed*(-np.sin(u0)*s_k+np.cos(u0)*w)]
        # Tiles on one ring are time-shifted copies of the ring's centre trajectory, so the tide
        # perturbs them alike; two-body circles diverge by hundreds of metres within a day.
        lags = (np.arange(per_ring)-(per_ring-1)/2)*pitch/speed
        if time_shift:
            track = propagate(env, centre[None], lags.min()-1., lags.max()+1.)
            states.extend(track(lag)[0] for lag in lags)
        else:
            for lag in lags:
                phi = u0+lag*speed/radius
                states.append(np.r_[radius*(np.cos(phi)*s_k+np.sin(phi)*w), speed*(-np.sin(phi)*s_k+np.cos(phi)*w)])
        ring_of.extend([k]*per_ring)
    states = np.array(states)
    at = propagate(env, states, 0., hours*3600.)
    times = np.linspace(0., hours*3600., int(hours*16)+1)
    t = np.radians(tilt_deg)
    distances, thickness = [], []
    for tt in times:
        s = at(tt)
        if attitude == 'own':
            distances.append(own_attitude_distance(s, t))
            continue
        Rt, Tt, Nt = rtn_frame(s.mean(axis=0))
        axes = np.c_[np.sin(t)*Rt+np.cos(t)*Tt, Nt, np.cos(t)*Rt-np.sin(t)*Tt]
        distances.append(closest_squares(s[:, :3], axes, SIDE)[0])
    distances = np.array(distances)
    k = int(np.argmin(distances))
    radii = np.linalg.norm(at(times[-1])[:, :3], axis=1)
    return dict(rings=rings, per_ring=per_ring, pitch_m=pitch, radius_step_m=step, tilt_deg=tilt_deg,
                time_shifted_rings=time_shift, attitude=attitude,
                hours=hours, samples=len(times), start_u_deg=float(np.degrees(u0)),
                min_surface_distance_m=float(distances[k]), at_hour=float(times[k]/3600.),
                samples_below_clearance=int(np.count_nonzero(distances < CLEARANCE)),
                final_radius_range_km=[float(radii.min()/1e3), float(radii.max()/1e3)])


def main():
    t0 = time.process_time()
    cases = [screen('equal_energy_patch', depth=d, fold_step=f) for d in (2e3, 5e3, 10e3) for f in (0.2e3, 0.5e3)]
    cases += [screen('drifting_patch', depth=d) for d in (0.3e3, 1e3, 5e3)]
    cases += [screen('equal_energy_patch', attitude='sun', depth=d) for d in (2e3, 5e3)]
    cases += [screen('ring_bundle', depth=d, tilt_deg=t) for d in (0.2e3, 0.3e3, 0.5e3) for t in (1., 2.)]
    # Shingle tilt and radius step chosen from the two clearance conditions, run over three orbits.
    cases += [screen('ring_bundle', orbits=3, depth=d, tilt_deg=t) for d, t in ((0.4e3, 1.2), (0.45e3, 1.5), (0.6e3, 2.))]
    replay = [ephemeris_ring_bundle(step=0.45e3, tilt_deg=1.5, time_shift=False)]
    replay += [ephemeris_ring_bundle(step=d, tilt_deg=t) for d, t in ((0.45e3, 1.5), (0.6e3, 2.))]
    # Each tile keeps its own radial-facing attitude, referenced to its velocity.
    replay += [ephemeris_ring_bundle(step=d, tilt_deg=t, attitude='own') for d, t in
               ((0.6e3, 1.5), (0.8e3, 1.2), (0.8e3, 1.5), (1.0e3, 1.5))]
    out = dict(schema='terluna.research.relative-design/1',
               producer=dict(files={f: digest(ROOT/f) for f in FILES}, constants=constants_used(FILES),
                             inputs={'results/joint_seed.json': digest(SEED),
                                     'local_continuation/ephemeris.npz': digest(PACKAGE/'ephemeris.npz')}),
               evidence=('Layout cases: linear relative motion about a circular 15,000 km reference with parallel '
                         'squares in one common attitude, sampled every 0.5 degrees over the stated orbits and about '
                         'every degree over the service arc; sail force, mutual shadows, tides and nonlinearity are '
                         'absent.'),
               reading_rule=('A layout passes when no sample falls below clearance_m and both service overlaps stay '
                             'positive. The coupled model with sail force and mutual shadows verifies a passing layout '
                             'before any acceptance, and the protection gates decide it.'),
               clearance_m=CLEARANCE, cases=cases, ephemeris_ring_bundle=replay,
               ephemeris_evidence=('Gravity-only replay with the compact DE440 samples from the seed epoch, 93 hours '
                                   'sampled every 3.75 minutes. Rings start on circles about the Moon; with '
                                   'time_shifted_rings each ring is time-shifted copies of its centre trajectory. '
                                   'attitude common uses the bundle centroid frame; own gives each tile its own '
                                   'velocity-referenced radial-facing frame, with pair distances in the pair mean '
                                   'frame. Sail force and mutual shadows are absent.'),
               cpu_s=time.process_time()-t0)
    OUT.write_text(json.dumps(out, indent=1)+'\n')
    for case in cases:
        print(f"{case['kind']:20s} {case['attitude']:6s} {case['orbits']} {json.dumps({k: v for k, v in case['params'].items()})[:60]:60s} "
              f"min {case['min_surface_distance_m']:9.1f} m at u={case['at_u_deg']:7.1f}  "
              f"rows {case['adjacent_row_overlap_min_m']:8.1f}  along {case['along_row_overlap_min_m']:8.1f}  "
              f"drift {case['drift_spread_per_orbit_m']/1e3:7.1f} km")
    for r in replay:
        print('ephemeris', {k: r[k] for k in ('time_shifted_rings', 'attitude', 'radius_step_m', 'tilt_deg', 'min_surface_distance_m', 'at_hour',
                                               'samples_below_clearance', 'final_radius_range_km')})
    print('cpu', round(out['cpu_s'], 1))


if __name__ == '__main__':
    main()
