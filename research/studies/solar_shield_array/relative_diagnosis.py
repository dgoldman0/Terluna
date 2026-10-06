"""Relative-orbit diagnosis of the expanded-cycle seed formation.

Reproduces the findings of relative_orbits.md: the depth stack's spread of
orbital energy and the along-track drift it causes, an energy-matching burn
that removes the drift, the fold of the pattern through its orbital plane,
the relative orbital elements a redesign starts from, and a simplified screen
of how Earth's tide turns inclined rings against the Sun.

Every replay is gravity-only: Moon, Earth, Sun and planets as point masses
from the committed compact DE440 samples, in the Moon-centred frame the
coupled model uses (CyclingEnvironment.gravity). Sail force, mutual shadows,
finite-square contact and attitude belong to the coupled model.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.relative_diagnosis
"""
from __future__ import annotations

import hashlib
import json
import resource
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.spatial import cKDTree

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.cycling import CyclingEnvironment
from protection.dynamics.ephemeris import BODY_GM
from protection.dynamics.relative_orbit import (drift_per_revolution, energy_matching_burns, hcw_elements,
                                                minimum_rn_separation, normal_thickness, period, rtn_frame,
                                                rtn_relative, semi_major_axis)
from .cycling_search import CENTRAL

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PACKAGE = HERE/'local_continuation'
SEED = HERE/'results/joint_seed.json'
EXPANDED = HERE/'results/expanded_cycle.json'
OUT = HERE/'results/relative_diagnosis.json'
FILES = ['research/studies/solar_shield_array/relative_diagnosis.py', 'protection/dynamics/relative_orbit.py',
         'protection/dynamics/cycling.py', 'protection/dynamics/ephemeris.py',
         'research/studies/solar_shield_array/cycling_search.py']
HOUR = 3600.
RTOL, ATOL = 1e-11, 1e-4
CLOSE = 10e3


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def environment():
    """The coupled model's environment, built from the committed compact samples."""
    with np.load(PACKAGE/'ephemeris.npz', allow_pickle=False) as z:
        samples = {k: z[k] for k in z.files if not k.startswith(('p_', 'v_'))}
        samples['positions'] = {n: z['p_'+n] for n in BODY_GM}
        samples['velocities'] = {n: z['v_'+n] for n in BODY_GM}
    return CyclingEnvironment(samples, .05, CENTRAL, 'filter_only')


def propagate(env, states, t0, t1):
    n = len(states)

    def rhs(t, y):
        s = y.reshape(n, 6)
        return np.hstack([s[:, 3:], env.gravity(t, s[:, :3])]).ravel()
    sol = solve_ivp(rhs, (t0, t1), np.asarray(states, float).ravel(), method='DOP853',
                    rtol=RTOL, atol=ATOL, dense_output=True)
    if not sol.success:
        raise RuntimeError(sol.message)
    return lambda t: sol.sol(t).reshape(n, 6)


def sun_direction(env, t):
    p = env.at(t)['positions']['sun']
    return p/np.linalg.norm(p)


def sun_coordinates(env, states, t):
    """Depth along the Sun, along-track and cross-track offsets from the centroid."""
    s = sun_direction(env, t)
    centre = states.mean(axis=0)
    along = centre[3:]-np.dot(centre[3:], s)*s
    along /= np.linalg.norm(along)
    cross = np.cross(s, along)
    d = states[:, :3]-centre[:3]
    return d@s, d@along, d@cross


def energy_spread(env, states, t):
    a = semi_major_axis(states)
    depth, along, cross = sun_coordinates(env, states, t)
    design = np.c_[np.ones_like(depth), depth, along, cross]
    coef = np.linalg.lstsq(design, a, rcond=None)[0]
    fit = design@coef
    r2 = 1-np.sum((a-fit)**2)/np.sum((a-a.mean())**2)
    centre = states.mean(axis=0)
    drift = drift_per_revolution(a, a.mean())
    return dict(time_h=t/HOUR,
                semi_major_axis_mean_m=float(a.mean()), semi_major_axis_span_m=float(np.ptp(a)),
                semi_major_axis_std_m=float(a.std()), period_span_s=float(np.ptp(period(a))),
                depth_span_m=float(np.ptp(depth)),
                fit_semi_major_axis_per_depth=float(coef[1]), fit_per_along_track=float(coef[2]),
                fit_per_cross_track=float(coef[3]), fit_r_squared=float(r2),
                correlation_with_depth=float(np.corrcoef(a, depth)[0, 1]),
                centroid_angle_from_sun_deg=float(np.degrees(np.arccos(
                    np.dot(centre[:3]/np.linalg.norm(centre[:3]), sun_direction(env, t))))),
                predicted_drift_per_revolution_max_m=float(np.abs(drift).max()),
                predicted_drift_per_revolution_rms_m=float(np.sqrt(np.mean(drift**2))))


def revolution(states_at, t0, first):
    """Time at which the centroid completes one revolution in its initial orbit plane."""
    frame = rtn_frame(first.mean(axis=0))
    def angle(t):
        c = states_at(t)[:, :3].mean(axis=0)
        return np.arctan2(c@frame[1], c@frame[0])
    T = period(semi_major_axis(first.mean(axis=0)[None])[0])
    # The angle increases through pi and wraps; bracket the return to zero near one period.
    times = np.linspace(t0+.8*T, t0+1.15*T, 400)
    values = np.array([angle(t) for t in times])
    k = int(np.nonzero((values[:-1] < 0) & (values[1:] >= 0))[0][0])
    return brentq(angle, times[k], times[k+1], xtol=1e-6)


def pattern(states):
    centre = states.mean(axis=0)
    return (states[:, :3]-centre[:3])@rtn_frame(centre).T


def repeat_error(first, last):
    e = pattern(last)-pattern(first)
    norm = np.linalg.norm(e, axis=1)
    return dict(radial_max_m=float(np.abs(e[:, 0]).max()), along_track_max_m=float(np.abs(e[:, 1]).max()),
                normal_max_m=float(np.abs(e[:, 2]).max()), total_max_m=float(norm.max()),
                total_rms_m=float(np.sqrt(np.mean(norm**2))), along_track=e[:, 1])


def energy_matching(env, states, t0, refinements=2):
    """Along-track burns that equalize energy, refined from the replayed one-revolution drift."""
    target = np.full(len(states), semi_major_axis(states).mean())
    passes = []
    for k in range(refinements+1):
        burns = energy_matching_burns(states, target)
        matched = states.copy()
        matched[:, 3:] += burns
        at = propagate(env, matched, t0, t0+52*HOUR)
        t1 = revolution(at, t0, matched)
        error = repeat_error(matched, at(t1))
        dv = np.linalg.norm(burns, axis=1)
        passes.append(dict(refinement=k, burn_rms_m_s=float(np.sqrt(np.mean(dv**2))),
                           burn_max_m_s=float(dv.max()), revolution_h=t1/HOUR,
                           **{key: v for key, v in error.items() if key != 'along_track'}))
        # A tile ahead by y needs its semi-major axis raised by y/(3*pi) to fall back.
        target = target+error['along_track']/(3*np.pi)
    return passes, matched, at, t1


def fold(env, states, t0, t1, step=300.):
    at = propagate(env, states, t0, t1)
    times = np.arange(t0, t1+1, step)
    thickness, close = [], []
    for t in times:
        s = at(t)
        thickness.append(normal_thickness(s))
        close.append(len(cKDTree(s[:, :3]).query_pairs(CLOSE)))
    thickness, close = np.array(thickness), np.array(close)
    k, j = int(np.argmin(thickness)), int(np.argmax(close))
    return dict(start_h=t0/HOUR, end_h=t1/HOUR, step_s=step,
                thickness_start_m=float(thickness[0]), thickness_min_m=float(thickness[k]),
                thickness_min_h=float(times[k]/HOUR), close_pairs_start=int(close[0]),
                close_pairs_max=int(close[j]), close_pairs_max_h=float(times[j]/HOUR),
                series=dict(time_h=(times/HOUR).round(4).tolist(), thickness_m=thickness.round(1).tolist(),
                            centre_pairs_within_10_km=close.tolist()))


def element_spread(env, states, t):
    """Relative orbital elements about the centroid, with u measured from the Sun direction."""
    centre = states.mean(axis=0)
    a = semi_major_axis(centre[None])[0]
    n = np.sqrt(K.MOON_GM/a**3)
    frame = rtn_frame(centre)
    s = sun_direction(env, t)
    s_plane = s-np.dot(s, frame[2])*frame[2]
    s_plane /= np.linalg.norm(s_plane)
    # Argument of latitude of the centroid measured from the projected Sun direction.
    u = np.arctan2(np.dot(np.cross(s_plane, frame[0]), frame[2]), np.dot(s_plane, frame[0]))
    elements = hcw_elements(rtn_relative(centre, states), a, n, u)
    de, di = elements[:, 2:4], elements[:, 4:6]
    i, j = np.triu_indices(len(states), 1)
    rn = minimum_rn_separation(de[i]-de[j], di[i]-di[j], a)
    phase = np.degrees(np.abs(np.angle(np.exp(1j*(np.arctan2(de[:, 1], de[:, 0])-
                                                  np.arctan2(di[:, 1], di[:, 0]))))))
    q = [5, 50, 95]
    return dict(time_h=t/HOUR, reference_semi_major_axis_m=float(a), centroid_u_from_sun_deg=float(np.degrees(u)),
                da_m_percentiles=dict(zip(map(str, q), np.percentile(a*elements[:, 0], q).round(1).tolist())),
                de_m_percentiles=dict(zip(map(str, q), np.percentile(a*np.linalg.norm(de, axis=1), q).round(1).tolist())),
                di_m_percentiles=dict(zip(map(str, q), np.percentile(a*np.linalg.norm(di, axis=1), q).round(1).tolist())),
                de_di_phase_gap_deg_percentiles=dict(zip(map(str, q), np.percentile(phase, q).round(1).tolist())),
                pairs=int(len(rn)), pairs_rn_min_below_100_m=int(np.count_nonzero(rn < 100.)),
                pairs_rn_min_below_1_km=int(np.count_nonzero(rn < 1e3)),
                pairs_rn_min_below_10_km=int(np.count_nonzero(rn < 1e4)))


def plane_tracking(radii=(12e6, 15e6, 17e6, 19e6, 21e6), tilts=(10., 20.), days=120.):
    """Node rates of inclined retrograde circular lunar orbits, Earth and Sun on circular coplanar orbits."""
    earth_rate = np.sqrt((K.EARTH_GM+K.MOON_GM)/K.EARTH_MOON_DISTANCE**3)
    sun_rate = 2*np.pi/(365.25*K.JULIAN_DAY)

    def bodies(t):
        return [(K.EARTH_GM, K.EARTH_MOON_DISTANCE*np.array([np.cos(earth_rate*t+1.), np.sin(earth_rate*t+1.), 0.])),
                (K.SUN_GM, K.AU*np.array([np.cos(sun_rate*t), np.sin(sun_rate*t), 0.]))]
    cases, initial = [], []
    for tilt in tilts:
        for radius in radii:
            i, node = np.pi-np.radians(tilt), np.pi/2
            line = np.array([np.cos(node), np.sin(node), 0.])
            normal = np.array([np.sin(node)*np.sin(i), -np.cos(node)*np.sin(i), np.cos(i)])
            # Start at the point 90 degrees from the node, where the orbit faces the Sun.
            position = np.cross(normal, line)
            initial.append(np.r_[radius*position, np.sqrt(K.MOON_GM/radius)*np.cross(normal, position)])
            cases.append((tilt, radius))
    initial = np.array(initial)
    m = len(initial)

    def rhs(t, y):
        s = y.reshape(m, 6)
        r = s[:, :3]
        a = -K.MOON_GM*r/np.linalg.norm(r, axis=1)[:, None]**3
        for gm, p in bodies(t):
            d = p-r
            a += gm*(d/np.linalg.norm(d, axis=1)[:, None]**3-p/np.linalg.norm(p)**3)
        return np.hstack([s[:, 3:], a]).ravel()
    T = days*K.JULIAN_DAY
    sol = solve_ivp(rhs, (0., T), initial.ravel(), method='DOP853', rtol=1e-10, atol=1e-3,
                    dense_output=True, max_step=3600.)
    times = np.linspace(0., T, 2401)
    states = sol.sol(times).reshape(m, 6, -1)
    rows = []
    for k, (tilt, radius) in enumerate(cases):
        h = np.cross(states[k, :3].T, states[k, 3:].T)
        h /= np.linalg.norm(h, axis=1)[:, None]
        node = np.unwrap(np.arctan2(h[:, 0], -h[:, 1]))
        rate = np.polyfit(times/K.JULIAN_DAY, np.degrees(node), 1)[0]
        distance = np.linalg.norm(states[k, :3], axis=0)
        rows.append(dict(tilt_deg=tilt, radius_km=radius/1e3, node_rate_deg_day=float(rate),
                         node_minus_sun_deg=float(np.degrees(node[-1]-node[0]-sun_rate*T)),
                         distance_range_km=[float(distance.min()/1e3), float(distance.max()/1e3)],
                         service_fraction=float(np.arcsin(7.0e6/radius)/np.pi)))
    return dict(model='Earth and Sun on circular coplanar orbits about the Moon; point gravity; no sail force',
                days=days, sun_rate_deg_day=float(np.degrees(sun_rate)*K.JULIAN_DAY),
                service_fraction_rule='arcsin(7,000 km / radius) / pi', cases=rows)


def main():
    wall, cpu = time.monotonic(), time.process_time()
    seed = json.loads(SEED.read_text())
    arrays = {k: np.array(v) for k, v in seed['arrays'].items()}
    env = environment()
    epochs = dict(initial_epoch_state=0., first_service_end_state=6*HOUR, original_twelve_hour_end_state=12*HOUR)
    spread = {k: energy_spread(env, arrays[k], t) for k, t in epochs.items()}

    twelve = arrays['original_twelve_hour_end_state']
    at = propagate(env, twelve, 12*HOUR, 64*HOUR)
    t1 = revolution(at, 12*HOUR, twelve)
    natural = repeat_error(twelve, at(t1))
    along = natural.pop('along_track')
    slope = float(np.polyfit(semi_major_axis(twelve)-semi_major_axis(twelve).mean(), along, 1)[0])
    passes, matched, _, _ = energy_matching(env, twelve, 12*HOUR)

    expanded = json.loads(EXPANDED.read_text())['products']
    coupled_first_event_h = dict(baseline=expanded['base']['clearance_events_s'][0]/HOUR,
                                 validated_trial=expanded['actual_validated']['clearance_events_s'][0]/HOUR)
    folds = dict(from_first_service_end=fold(env, arrays['first_service_end_state'], 6*HOUR, 24*HOUR),
                 from_twelve_hour_state=fold(env, twelve, 12*HOUR, 24*HOUR, step=120.),
                 coupled_first_clearance_event_h=coupled_first_event_h)
    elements = dict(as_seeded=element_spread(env, twelve, 12*HOUR),
                    energy_matched=element_spread(env, matched, 12*HOUR))
    planes = plane_tracking()

    out = dict(schema='terluna.research.relative-diagnosis/1',
               producer=dict(files={f: digest(ROOT/f) for f in FILES}, constants=constants_used(FILES),
                             inputs={'results/joint_seed.json': digest(SEED),
                                     'local_continuation/ephemeris.npz': digest(PACKAGE/'ephemeris.npz'),
                                     'results/expanded_cycle.json': digest(EXPANDED)}),
               evidence=('Gravity-only point-mass replays of the committed expanded-cycle seed states with the '
                         'compact DE440 samples; sail force, mutual shadows, finite-square contact and attitude '
                         'are absent. The plane-tracking screen uses circular coplanar Earth and Sun orbits.'),
               reading_rule=('Diagnosis of the seed formation and of an energy-matching correction. The coupled '
                             'finite-square model with sail force verifies any design built from these values, and '
                             'the protection gates decide acceptance.'),
               units=dict(length='m unless named km', velocity='m/s', time='hours from 2026-10-04 TDB'),
               energy_spread=spread,
               natural_repeat=dict(start_h=12., revolution_h=t1/HOUR, drift_slope_per_semi_major_axis=slope,
                                   three_pi=3*np.pi, **natural),
               energy_matching=passes, fold=folds, relative_elements=elements, plane_tracking=planes,
               resources=dict(cpu_s=time.process_time()-cpu, wall_s=time.monotonic()-wall,
                              max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
                              numerical_threads=1))
    OUT.write_text(json.dumps(out, indent=1)+'\n')
    print(json.dumps(dict(spread_12h=spread['original_twelve_hour_end_state']['semi_major_axis_span_m'],
                          natural_along_max_m=natural['along_track_max_m'],
                          matched=[(p['burn_rms_m_s'], p['along_track_max_m']) for p in passes],
                          cpu_s=out['resources']['cpu_s']), indent=1))


if __name__ == '__main__':
    main()
