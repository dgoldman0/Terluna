"""Coupled ring-bundle verification (stage 2b of relative_orbits.md).

Nineteen nested rings fly from the seed epoch under DE440s point-mass gravity
with finite-square quadrature, the filter's sail force on per-tile radial-
facing normals, finite-Sun mutual shadows and Moon/Earth eclipses
(ring_bundle). Three cases are run:

- a finite segment of 31 tiles per ring for two orbits, whose segment ends
  receive less shadow than interior tiles;
- continuous rings for six orbits, about twelve days: one representative
  tile per ring, with its ring and the adjacent rings present as two-body
  shifted copies (ring_bundle.ContinuousRings), so every tile is shadowed as on
  a complete ring;
- a lone string of seven tiles under gravity alone for 120 hours, which shows
  how the time-varying tide moves neighbours on one ring relative to each other
  and to the shifted copies.

The continuous rings' representatives also fly sixteen orbits, a month, under
gravity alone. At each mid-service of both continuous runs the runner records
every ring's elements about its inner neighbour carried to the same phase under
the tide (as the stage 2c keeping measures them), so their orbit-to-orbit
changes separate the tide's share of the drift between rings from the light's.
For the month under gravity it also measures the tilts in the frame of each
pair's node on Earth's orbital plane, which turns with the stack's common
regression. It records the fraction of each ring's light that reaches it past
the neighbouring rings' shadows inside the service arc.

Every five minutes the runner records the smallest surface distance between
tiles and, inside the service arc, the projected overlap between adjacent
rings and along each ring. At each mid-service it records each ring's
osculating semi-major axis, plane and tile pitch; their orbit-to-orbit changes
give a first estimate of the corrections that keep the bundle's geometry. The
estimate compares rings at different phases, whose osculating elements carry
the tide's short-period swing, so stage 2c's matched-phase keeping supersedes
it. Each orbit is a checkpointed segment in
research/runs/solar_shield_array/ring_bundle.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.ring_bundle_run
"""
from __future__ import annotations

import hashlib
import inspect
import json
import resource
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.optimize import brentq

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics import ring_bundle as rb
from protection.dynamics.cycling import CyclingEnvironment
from protection.dynamics.ephemeris import DEFAULT_KERNEL, Ephemeris
from protection.dynamics.natural_pattern import finite_gravity
from protection.dynamics.relative_orbit import rtn_frame, semi_major_axis
from protection.dynamics.ring_keeping import ContinuousKeeping, geometric_elements, gravity_shift
from .cycling_search import CENTRAL

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SEED = HERE/'results/joint_seed.json'
RUN = ROOT/'research/runs/solar_shield_array/ring_bundle'
OUT = HERE/'results/ring_bundle.json'
FILES = ['research/studies/solar_shield_array/ring_bundle_run.py', 'protection/dynamics/ring_bundle.py',
         'protection/dynamics/eclipse_parallel.py', 'protection/dynamics/eclipse_parallel.cpp',
         'protection/dynamics/fast_parallel.py', 'protection/dynamics/fast_parallel.cpp',
         'protection/dynamics/natural_pattern.py', 'protection/dynamics/cycling.py',
         'protection/dynamics/ephemeris.py', 'protection/dynamics/relative_orbit.py']
DESIGN = dict(rings=19, per_ring=31, pitch_m=8.5e3, radius_step_m=600., tilt_deg=1.5,
              segment_orbits=2, continuous_orbits=6, gravity_only_orbits=16, copy_reach=3, service_half_angle_deg=25.,
              sample_s=300., interior_margin=6, rtol=1e-10, atol=1e-3, max_step_s=600.)
DYNAMICS = [f for f in FILES if f != 'research/studies/solar_shield_array/ring_bundle_run.py']
ANALYSIS = ['protection/dynamics/ring_keeping.py']
HOUR = 3600.


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def environment(days):
    eph = Ephemeris(epoch='2026-10-04')
    samples = eph.sample(np.arange(-3600., days*K.JULIAN_DAY+7200., 1800.))
    eph.close()
    return CyclingEnvironment(samples, .05, CENTRAL, 'filter_only')


def initial(env, step=None):
    """Bundle states from the seed epoch; step overrides the designed radius step (m)."""
    step = DESIGN['radius_step_m'] if step is None else step
    seed = json.loads(SEED.read_text())['arrays']
    ref = np.array(seed['initial_epoch_state']).mean(axis=0)
    ratio = float(np.mean(seed['mass_ratio']))
    R, T, N = rtn_frame(ref)
    sun = env.at(0.)['positions']['sun']
    sp = sun-np.dot(sun, N)*N
    sp /= np.linalg.norm(sp)
    w = np.cross(N, sp)
    u0 = np.arctan2(np.dot(R, w), np.dot(R, sp))

    def track(centre, lags):
        at, _ = rb.propagate(lambda t, s: env.gravity(t, s[:, :3]), centre[None], lags.min()-1., lags.max()+1.,
                             max_step=60.)
        return [at(lag)[0] for lag in lags]
    height = DESIGN['pitch_m']/np.cos(np.radians(23.75))
    states, ring, slot = rb.build(np.linalg.norm(ref[:3]), track, DESIGN['rings'], DESIGN['per_ring'],
                                  DESIGN['pitch_m'], step, height, sp, N, u0)
    return states, ring, slot, ratio, dict(start_u_deg=float(np.degrees(u0)), height_pitch_m=height,
                                           reference_radius_m=float(np.linalg.norm(ref[:3])))


def sun_angle(env, t, states):
    """Angle of the bundle centroid from the Sun direction projected into its orbit plane, toward the motion."""
    c = states.mean(axis=0)
    R, T, N = rtn_frame(c)
    sun = env.at(t)['positions']['sun']
    sp = sun-np.dot(sun, N)*N
    sp /= np.linalg.norm(sp)
    return float(np.arctan2(np.dot(R, np.cross(N, sp)), np.dot(R, sp)))


def run(env, states, model, orbits, label):
    """Propagate orbit by orbit with checkpoints; return sample times, states, period and CPU."""
    RUN.mkdir(parents=True, exist_ok=True)
    period = 2*np.pi*np.sqrt(semi_major_axis(states.mean(axis=0)[None])[0]**3/K.MOON_GM)
    settings = {k: DESIGN[k] for k in ('rings', 'per_ring', 'pitch_m', 'radius_step_m', 'tilt_deg', 'copy_reach',
                                       'sample_s', 'rtol', 'atol', 'max_step_s')}
    key = hashlib.sha256(json.dumps(settings, sort_keys=True).encode()+label.encode()+states.tobytes()+
                         ''.join(digest(ROOT/f) for f in DYNAMICS).encode()).hexdigest()[:16]
    times, samples, current, t0, cpu = [], [], states, 0., 0.
    for k in range(orbits):
        path = RUN/f'{label}_orbit_{k}_{key}.npz'
        t1 = (k+1)*period
        if path.exists():
            z = np.load(path)
            times.append(z['t']); samples.append(z['s']); current = z['s'][-1]; cpu += float(z['cpu'])
        else:
            start = time.process_time()
            at, sol = rb.propagate(model.acceleration, current, t0, t1, rtol=DESIGN['rtol'], atol=DESIGN['atol'],
                                   max_step=DESIGN['max_step_s'])
            grid = np.arange(t0, t1, DESIGN['sample_s'])
            if grid[-1] < t1:
                grid = np.r_[grid, t1]
            s = np.stack([at(t) for t in grid])
            used = time.process_time()-start
            tmp = path.with_name(path.stem+'.tmp.npz')
            np.savez(tmp, t=grid, s=s, cpu=used, nfev=sol.nfev)
            tmp.rename(path)
            times.append(grid); samples.append(s); current = s[-1]; cpu += used
        t0 = t1
    t = np.concatenate([x if i == 0 else x[1:] for i, x in enumerate(times)])
    s = np.concatenate([x if i == 0 else x[1:] for i, x in enumerate(samples)])
    return t, s, period, cpu, key


def overlaps(states, ring, slot, sun_dir, tilt, inner=None):
    """Projected overlaps along the Sun: between adjacent rings and along each ring, over the chosen tiles."""
    if inner is None:
        m = DESIGN['interior_margin']
        inner = (slot >= m) & (slot < DESIGN['per_ring']-m)
    frames = rb.tile_frames(states, tilt)
    c = states.mean(axis=0)
    normal = rtn_frame(c)[2]
    e2 = normal-np.dot(normal, sun_dir)*sun_dir
    e2 /= np.linalg.norm(e2)
    e1 = np.cross(e2, sun_dir)
    half1 = rb.CLEAR/2*(np.abs(frames[:, :, 0]@e1)+np.abs(frames[:, :, 1]@e1))
    half2 = rb.CLEAR/2*(np.abs(frames[:, :, 0]@e2)+np.abs(frames[:, :, 1]@e2))
    x1, x2 = states[:, :3]@e1, states[:, :3]@e2
    along, across = np.inf, np.inf
    for k in np.unique(ring):
        sel = np.nonzero((ring == k) & inner)[0]
        # Only tiles with consecutive indices are neighbours on a ring; continuous rings carry separate groups of
        # copies around each representative.
        order = sel[np.argsort(slot[sel])]
        consecutive = np.diff(slot[order]) == 1
        gaps = (half1[order[:-1]]+half1[order[1:]]-np.abs(np.diff(x1[order])))[consecutive]
        if len(gaps):
            along = min(along, float(gaps.min()))
        if k+1 in ring:
            other = np.nonzero((ring == k+1) & inner)[0]
            d1 = np.abs(x1[sel][:, None]-x1[other][None, :])
            facing = d1 < (half1[sel][:, None]+half1[other][None, :])
            o = half2[sel][:, None]+half2[other][None, :]-np.abs(x2[sel][:, None]-x2[other][None, :])
            if np.any(facing):
                across = min(across, float(o[facing].min()))
    return along, across


def ring_state(states, ring, pitch):
    """Per ring: mean osculating semi-major axis and plane normal of the given tiles, and the ring's pitch."""
    out = []
    for k in np.unique(ring):
        sel = ring == k
        h = np.cross(states[sel, :3], states[sel, 3:]).mean(axis=0)
        out.append(dict(a=float(semi_major_axis(states[sel]).mean()), h=h/np.linalg.norm(h), pitch=float(pitch[k])))
    return out


def analyse(env, t, s, ring, slot, tilt, n, v, copies=None, interior=None):
    """Clearance series, service overlaps, mid-service ring states and per-orbit corrections.

    interior(ring, index) selects the tiles, copies included, of a second
    clearance series that leaves out the segment ends or the edge rings.
    """
    clear, pairs, service = [], [], []
    angles = np.array([sun_angle(env, tt, ss) for tt, ss in zip(t, s)])
    inside, inside_pairs = [], []
    for tt, ss, ang in zip(t, s, angles):
        if copies is not None:
            full, r_, j_ = copies(ss)
        else:
            full, r_, j_ = ss, ring, slot
        d, pair = rb.own_attitude_clearance(full, tilt)
        clear.append(d); pairs.append((int(r_[pair[0]]), int(j_[pair[0]]), int(r_[pair[1]]), int(j_[pair[1]])))
        if interior is not None:
            sel = np.nonzero(interior(r_, j_))[0]
            d, pair = rb.own_attitude_clearance(full[sel], tilt)
            inside.append(d)
            inside_pairs.append(tuple(int(x) for x in (r_[sel[pair[0]]], j_[sel[pair[0]]], r_[sel[pair[1]]],
                                                       j_[sel[pair[1]]])))
        if abs(ang) <= np.radians(DESIGN['service_half_angle_deg']):
            sun = env.at(tt)['positions']['sun']-ss[:, :3].mean(axis=0)
            mask = None if copies is None else np.ones(len(full), bool)
            along, across = overlaps(full, r_, j_, sun/np.linalg.norm(sun), tilt, mask)
            service.append((tt, along, across))
    clear, service = np.array(clear), np.array(service)
    mids = [t[i]+(t[i+1]-t[i])*(-angles[i])/(angles[i+1]-angles[i]) for i in range(len(t)-1)
            if angles[i] < 0 <= angles[i+1]]
    by_orbit = []
    for tm in mids:
        j = int(np.argmin(np.abs(t-tm)))
        if copies is not None:
            full, r_, j_ = copies(s[j])
            reps = j_ == 0
            pitch = {}
            for k in np.unique(ring):
                one = (r_ == k) & (j_ == 1)
                pitch[k] = float(np.linalg.norm(full[one][0, :3]-full[(r_ == k) & reps][0, :3]))
            rs = ring_state(full[reps], r_[reps], pitch)
        else:
            pitch = {}
            m = DESIGN['interior_margin']
            for k in np.unique(ring):
                sel = (ring == k)
                q = s[j][sel][np.argsort(slot[sel]), :3]
                pitch[k] = float(np.linalg.norm(np.diff(q, axis=0), axis=1)[m:-m].mean())
            rs = ring_state(s[j], ring, pitch)
        da = np.diff([r['a'] for r in rs])
        planes = [np.degrees(np.arccos(np.clip(np.dot(rs[i]['h'], rs[i+1]['h']), -1, 1))) for i in range(len(rs)-1)]
        by_orbit.append(dict(time_h=float(t[j]/HOUR), adjacent_a_step_m=[float(da.min()), float(da.max())],
                             adjacent_plane_angle_deg=[float(min(planes)), float(max(planes))],
                             pitch_m=[float(min(pitch.values())), float(max(pitch.values()))],
                             a_steps=da.tolist(), plane_angles_deg=planes, pitches=[pitch[k] for k in sorted(pitch)]))
    corrections = []
    for a, b in zip(by_orbit[:-1], by_orbit[1:]):
        step = np.abs(np.array(b['a_steps'])-np.array(a['a_steps']))
        plane = np.radians(np.abs(np.array(b['plane_angles_deg'])-np.array(a['plane_angles_deg'])))
        pitch = np.abs(np.array(b['pitches'])-np.array(a['pitches']))
        # Restoring a radius step costs n*da/2; a plane angle v*angle; a pitch change reversed within an orbit,
        # by two along-track burns of n*(dp/3pi)/2.
        corrections.append(dict(from_h=a['time_h'], to_h=b['time_h'], a_step_change_max_m=float(step.max()),
                                plane_angle_change_max_deg=float(np.degrees(plane.max())),
                                pitch_change_max_m=float(pitch.max()), dv_radius_m_s=float(n*step.max()/2),
                                dv_plane_m_s=float(v*plane.max()), dv_pitch_m_s=float(n*pitch.max()/(3*np.pi))))
    def first_below(series, at):
        if not np.any(series < 150.):
            return dict(first_below_150_m_h=None, first_below_150_m_pair=None)
        i = int(np.argmax(series < 150.))
        return dict(first_below_150_m_h=float(t[i]/HOUR), first_below_150_m_pair=list(at[i]))
    k = int(np.argmin(clear))
    out = dict(clearance=dict(min_m=float(clear[k]), at_h=float(t[k]/HOUR), pair_ring_index=list(pairs[k]),
                              samples=len(clear), samples_below_150_m=int(np.count_nonzero(clear < 150.)),
                              **first_below(clear, pairs),
                              series_h=(t/HOUR).round(3).tolist()[::6], series_m=clear.round(1).tolist()[::6]),
               service=dict(samples=len(service), along_ring_overlap_min_m=float(service[:, 1].min()),
                            adjacent_ring_overlap_min_m=float(service[:, 2].min())),
               mid_service=by_orbit, corrections_per_orbit=corrections)
    if inside:
        inside = np.array(inside)
        k = int(np.argmin(inside))
        out['interior_clearance'] = dict(min_m=float(inside[k]), at_h=float(t[k]/HOUR),
                                         pair_ring_index=list(inside_pairs[k]),
                                         samples_below_150_m=int(np.count_nonzero(inside < 150.)),
                                         **first_below(inside, inside_pairs), series_m=inside.round(1).tolist()[::6])
    return out


class GravityOnly:
    """Representatives under the environment's gravity with finite-square quadrature and no light."""
    def __init__(self, env, tilt):
        self.env, self.tilt = env, tilt

    def acceleration(self, t, states):
        frames = rb.tile_frames(states, self.tilt)
        return finite_gravity(self.env, t, states[:, :3], frames[:, :, 0], frames[:, :, 1])


def mid_service_indices(env, t, s):
    angles = np.array([sun_angle(env, tt, ss) for tt, ss in zip(t, s)])
    return [i if -angles[i] <= angles[i+1] else i+1 for i in range(len(t)-1) if angles[i] < 0 <= angles[i+1]]


def element_drift(env, t, s, origin, a, every=6):
    """Each ring's elements about its inner neighbour at the same phase, per mid-service, and their changes (m).

    Columns follow ring_keeping.geometric_elements times a: semi-major-axis
    step, eccentricity-vector difference and plane difference. Pairs 0 and 17
    are the bundle's edges. The within-orbit ranges use every sixth sample
    (half-hourly) between mid-services; less the orbit's change they show the
    short-period part that averages out over an orbit.
    """
    keeper = ContinuousKeeping(SimpleNamespace(env=env), s[0], 0, origin, 1.)
    mids = mid_service_indices(env, t, s)
    rows = [(float(t[i]/HOUR), keeper.relative(t[i], s[i])[1:]*a) for i in mids]
    ranges = []
    for i0, i1 in zip(mids[:-1], mids[1:]):
        x = np.array([keeper.relative(t[i], s[i])[1:]*a for i in range(i0, i1+1, every)])
        ranges.append(np.stack([np.ptp(x[:, :, 0], axis=0), np.ptp(np.hypot(x[:, :, 2], x[:, :, 3]), axis=0),
                                np.ptp(np.hypot(x[:, :, 4], x[:, :, 5]), axis=0)], axis=1))
    ranges = np.array(ranges)
    changes = []
    for (h0, e0), (h1, e1) in zip(rows[:-1], rows[1:]):
        d = e1-e0
        changes.append(dict(from_h=h0, to_h=h1, a_step_change_m=np.abs(d[:, 0]).round(2).tolist(),
                            e_change_m=np.hypot(d[:, 2], d[:, 3]).round(2).tolist(),
                            plane_change_m=np.hypot(d[:, 4], d[:, 5]).round(2).tolist()))
    e = np.array([c['e_change_m'] for c in changes])
    da = np.array([c['a_step_change_m'] for c in changes])
    plane = np.array([r[:, 4:6] for _, r in rows])
    departure = np.linalg.norm(plane-plane[0], axis=2)
    return dict(mid_service_h=[h for h, _ in rows], a_step_m=[r[:, 0].round(2).tolist() for _, r in rows],
                e_separation_m=[np.hypot(r[:, 2], r[:, 3]).round(2).tolist() for _, r in rows],
                per_orbit=changes,
                interior_e_change_median_m=float(np.median(e[:, 1:-1])),
                interior_e_change_max_m=float(e[:, 1:-1].max()),
                edge_e_change_max_m=[float(e[:, 0].max()), float(e[:, -1].max())],
                interior_a_step_change_max_m=float(da[:, 1:-1].max()),
                edge_a_step_change_max_m=[float(da[:, 0].max()), float(da[:, -1].max())],
                plane_offset_m={str(k): plane[:, k].round(1).tolist() for k in (0, 8, 17)},
                plane_departure_max_m=departure.max(axis=0).round(1).tolist(),
                plane_departure_final_m=departure[-1].round(1).tolist(),
                within_orbit_range_m=dict(columns=['a_step', 'e_separation', 'plane_separation'],
                                          by_orbit_pair=ranges.round(2).tolist(),
                                          interior_median=np.median(ranges[:, 1:-1], axis=(0, 1)).round(2).tolist()))


def earth_moon_pole(env):
    """Unit pole of Earth's orbit about the Moon at the start, about which the tide turns the ring planes."""
    pole = np.cross(env.at(0.)['positions']['earth'], env.at(HOUR)['positions']['earth'])
    return pole/np.linalg.norm(pole)


def node_frame(env, t, s, a):
    """Relative tilts of adjacent rings in the regressing frame of each pair's ascending node (m).

    The tide turns all ring planes together about Earth's orbital pole. In a
    fixed frame that common regression turns every pair's tilt difference with
    it; measured from the inner ring's node on Earth's orbital plane, at the same
    phase, the difference shows only how the stack itself deforms.
    """
    pole = earth_moon_pole(env)
    mid = (len(s[0])-1)//2

    def node(state):
        h = np.cross(state[:3], state[3:])
        n = np.cross(pole, h/np.linalg.norm(h))
        return n/np.linalg.norm(n)
    rows, turn = [], []
    first = node(s[0][mid])
    for i in mid_service_indices(env, t, s):
        tilts = []
        for m in range(len(s[i])-1):
            h = np.cross(s[i][m, :3], s[i][m, 3:])
            rate = np.linalg.norm(h)/np.dot(s[i][m, :3], s[i][m, :3])
            ahead = np.arctan2(np.dot(h/np.linalg.norm(h), np.cross(s[i][m, :3], s[i][m+1, :3])),
                               np.dot(s[i][m, :3], s[i][m+1, :3]))
            matched = gravity_shift(env, t[i], s[i][m], ahead/rate)
            tilts.append(geometric_elements(matched, s[i][m+1][None], node(matched))[0][0, 4:6]*a)
        rows.append(np.array(tilts))
        n = node(s[i][mid])
        turn.append(float(np.degrees(np.arctan2(np.dot(np.cross(first, n), pole), np.dot(first, n)))))
    rows = np.array(rows)
    departure = np.linalg.norm(rows-rows[0], axis=2)
    inclination = np.degrees(np.arccos(abs(np.dot(np.cross(s[0][mid, :3], s[0][mid, 3:]), pole))/
                                       np.linalg.norm(np.cross(s[0][mid, :3], s[0][mid, 3:]))))
    return dict(inclination_to_earth_orbit_deg=float(inclination), node_turn_deg=turn,
                tilt_m={str(k): rows[:, k].round(1).tolist() for k in (0, 8, 17)},
                departure_max_m=departure.max(axis=0).round(1).tolist(),
                departure_final_m=departure[-1].round(1).tolist())


def lit_fractions(env, rings, t, s):
    """Fraction of each ring's sunlight that reaches it past the bundle's shadows, inside the service arc.

    Statistics use every sample inside the arc; the mid-service rows repeat the
    same Sun geometry each orbit, so their spread is the sliding rings' share.
    """
    def lit(i):
        _, first, bare = rings.forces.light(t[i], rings.copies(s[i])[0])
        return first[:len(s[i])]/bare[:len(s[i])]
    arc = [i for i in range(len(t)) if abs(sun_angle(env, t[i], s[i])) <= np.radians(DESIGN['service_half_angle_deg'])]
    inside = np.array([lit(i) for i in arc])
    mids = np.array([lit(i) for i in mid_service_indices(env, t, s)])
    spread = (mids.max(axis=0)-mids.min(axis=0))/mids.mean(axis=0)
    return dict(service_samples=len(arc), service_min_by_ring=inside.min(axis=0).round(4).tolist(),
                service_max_by_ring=inside.max(axis=0).round(4).tolist(),
                service_mean_by_ring=inside.mean(axis=0).round(4).tolist(),
                mid_service_by_orbit=mids.round(4).tolist(),
                mid_service_spread_over_mean_by_ring=spread.round(4).tolist(),
                interior_mid_service_spread_max=float(spread[1:-1].max()),
                interior_mid_service_spread_median=float(np.median(spread[1:-1])))


def lone_ring(env, states, ring, slot, tilt, hours=120.):
    """Seven consecutive tiles of the middle ring under gravity alone."""
    mid = (DESIGN['rings']-1)//2
    centre = (DESIGN['per_ring']-1)//2
    one = states[(ring == mid) & (slot >= centre-3) & (slot <= centre+3)]
    lag = DESIGN['pitch_m']/np.linalg.norm(one[3, 3:])
    at, _ = rb.propagate(lambda t, s: env.gravity(t, s[:, :3]), one, 0., hours*HOUR, rtol=DESIGN['rtol'],
                         atol=DESIGN['atol'], max_step=DESIGN['max_step_s'])
    rows = []
    for t in np.arange(0., hours*HOUR+1, 6*HOUR):
        s = at(t)
        f = rb.tile_frames(s, tilt)
        normal, along = [], []
        for i in range(6):
            m = f[i]+f[i+1]
            m /= np.linalg.norm(m, axis=0)
            d = (s[i+1, :3]-s[i, :3])@m
            normal.append(abs(d[2])); along.append(d[0])
        copy = [float(np.linalg.norm(rb.kepler_shift(s[3][None], np.array([j*lag]))[0][:3]-s[3+j, :3]))
                for j in (-3, -1, 1, 3)]
        rows.append(dict(hour=float(t/HOUR), shingle_normal_m=[float(min(normal)), float(max(normal))],
                         along_spacing_m=[float(min(along)), float(max(along))], copy_error_m=copy))
    return dict(hours=hours, samples=rows,
                shingle_normal_range_m=[min(r['shingle_normal_m'][0] for r in rows), max(r['shingle_normal_m'][1] for r in rows)],
                along_spacing_range_m=[min(r['along_spacing_m'][0] for r in rows), max(r['along_spacing_m'][1] for r in rows)],
                copy_error_max_m=max(max(r['copy_error_m']) for r in rows))


def main():
    wall, cpu0 = time.monotonic(), time.process_time()
    days = DESIGN['continuous_orbits']*2.0+2
    env = environment(days)
    states, ring, slot, ratio, built = initial(env)
    tilt = np.radians(DESIGN['tilt_deg'])
    a0 = semi_major_axis(states.mean(axis=0)[None])[0]
    n, v = np.sqrt(K.MOON_GM/a0**3), np.sqrt(K.MOON_GM/a0)

    segment = rb.RingBundleForces(env, tilt, np.full(len(states), ratio))
    t, s, period, cpu_a, key_a = run(env, states, segment, DESIGN['segment_orbits'], 'segment')
    m = DESIGN['interior_margin']
    part_a = analyse(env, t, s, ring, slot, tilt, n, v,
                     interior=lambda r, j: (j >= m) & (j < DESIGN['per_ring']-m))

    middle = (DESIGN['per_ring']-1)//2
    reps = states[slot == middle]
    lags = DESIGN['pitch_m']/np.linalg.norm(reps[:, 3:], axis=1)
    rings = rb.ContinuousRings(env, tilt, np.full(len(reps), ratio), lags, reach=DESIGN['copy_reach'])
    t, s, period, cpu_b, key_b = run(env, reps, rings, DESIGN['continuous_orbits'], 'continuous')
    part_b = analyse(env, t, s, np.arange(len(reps)), np.zeros(len(reps), int), tilt, n, v, copies=rings.copies,
                     interior=lambda r, j: (r > 0) & (r < len(reps)-1))
    origin = env.at(0.)['positions']['sun']
    origin = origin/np.linalg.norm(origin)
    part_b['element_drift'] = element_drift(env, t, s, origin, a0)
    part_b['lit_fraction'] = lit_fractions(env, rings, t, s)
    # A month under gravity alone needs a longer ephemeris window; the label carries the model's source, which
    # lives in this runner rather than in the dynamics files, and the window.
    month = environment(DESIGN['gravity_only_orbits']*2.0+2)
    label = ('gravity_'+hashlib.sha256(inspect.getsource(GravityOnly).encode()).hexdigest()[:8]+
             f"_{DESIGN['gravity_only_orbits']}")
    tg, sg, _, cpu_g, key_g = run(month, reps, GravityOnly(month, tilt), DESIGN['gravity_only_orbits'], label)
    part_d = dict(orbits=DESIGN['gravity_only_orbits'], run_key=key_g, label=label, propagation_cpu_s=cpu_g,
                  element_drift=element_drift(month, tg, sg, origin, a0), node_frame=node_frame(month, tg, sg, a0))
    part_c = lone_ring(env, states, ring, slot, tilt)
    shadow = rings.forces.light(0., rings.copies(reps)[0])
    loss = 1-shadow[1][:len(reps)]/shadow[2][:len(reps)]

    out = dict(schema='terluna.research.ring-bundle/1',
               producer=dict(files={f: digest(ROOT/f) for f in FILES+ANALYSIS},
                             constants=constants_used([f for f in FILES+ANALYSIS if f.endswith('.py')]),
                             inputs={'results/joint_seed.json': digest(SEED), 'de440s.bsp': digest(DEFAULT_KERNEL)}),
               evidence=('Coupled propagation of 19 nested rings from 2026-10-04 TDB: DE440s point-mass gravity with '
                         'finite-square quadrature, the filter sail force on per-tile radial-facing normals, finite-Sun '
                         '(8 points) mutual shadows in the bundle centroid frame and Moon/Earth eclipses. The finite '
                         'segment has 31 tiles per ring for two orbits; continuous rings use one representative per '
                         'ring with two-body shifted copies of its ring and the adjacent rings for six orbits; the '
                         'lone ring is seven tiles under gravity alone. Clearance and overlaps are sampled every five '
                         'minutes; the interior clearance leaves out six tiles at each segment end, or the two edge '
                         'rings of the continuous rings. The continuous representatives also fly a month under gravity '
                         'alone. Element drift compares each ring with its inner neighbour carried to the same phase '
                         'under the tide at each mid-service, on axes fixed by the starting Sun direction and, for the '
                         'tilts, in the frame of each pair\'s node on Earth\'s orbital plane. Lit fractions are each '
                         'ring\'s sunlight past the bundle\'s shadows inside the service arc and at each mid-service.'),
               reading_rule=('Continuous rings stand for complete rings; the finite segment shows what segment ends '
                             'add. A bundle passes when no sample falls below 150 m, the service overlaps stay positive '
                             'and the per-orbit corrections fit the decision gate. Rings 0 and 18 are bundle edges.'),
               design=DESIGN, built=built, mass_ratio=ratio, period_h=period/HOUR,
               starting_shadow_loss=dict(interior_rings=float(np.median(loss[:-1])), top_ring=float(loss[-1])),
               segment=dict(run_key=key_a, propagation_cpu_s=cpu_a, **part_a),
               continuous=dict(run_key=key_b, propagation_cpu_s=cpu_b, **part_b),
               continuous_gravity_only=part_d,
               lone_ring_gravity=part_c,
               resources=dict(cpu_s=time.process_time()-cpu0, wall_s=time.monotonic()-wall,
                              max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, numerical_threads=1))
    OUT.write_text(json.dumps(out, indent=1)+'\n')
    for name in ('segment', 'continuous'):
        part = out[name]
        print(name, json.dumps(dict(clearance={k: v for k, v in part['clearance'].items() if not k.startswith('series')},
                                    interior=part.get('interior_clearance'), service=part['service'],
                                    corrections_max={k: max(c[k] for c in part['corrections_per_orbit'])
                                                     for k in ('dv_radius_m_s', 'dv_plane_m_s', 'dv_pitch_m_s')}
                                    if part['corrections_per_orbit'] else None), indent=1))
    for name, drift in (('with light', part_b['element_drift']), ('gravity only', part_d['element_drift'])):
        print('element drift', name, {k: v for k, v in drift.items() if k not in ('per_orbit', 'a_step_m', 'e_separation_m',
                                                                               'mid_service_h')})
    nf = part_d['node_frame']
    print('node frame', dict(inclination=nf['inclination_to_earth_orbit_deg'], turn=nf['node_turn_deg'][-1],
                             departure_max=max(nf['departure_max_m']), departure_final=max(nf['departure_final_m'])))
    print('lit fraction', {k: v for k, v in part_b['lit_fraction'].items() if not k.endswith(('by_ring', 'by_orbit'))},
          'service mean by ring', part_b['lit_fraction']['service_mean_by_ring'])
    print('lone ring', {k: v for k, v in part_c.items() if k != 'samples'})
    print('resources', out['resources'])


if __name__ == '__main__':
    main()
