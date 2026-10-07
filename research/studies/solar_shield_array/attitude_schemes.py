"""Attitude schemes for the ring fleet (item 1 of the narrowing that photon_control began).

photon_control found that radial-facing shingled tiles keep station with sunlight but do not hold their attitude:
within each ring a tile lies under its neighbour along one edge, so shadows put a steady pitch torque on it. This
runner compares the ways of holding the attitude and the push each leaves on the orbit.

Clearance. Rings step outward by 1 km and slide past one another, so tiles of neighbouring rings meet at every
along-track offset. Tiles that share an attitude are parallel squares. The runner maps, against orbital phase, the
pitch (about the orbit normal) and roll (about the line of tile centres) offsets from the radial-facing attitude at
which every tile of the three nearest rings on each side, at every offset, and the two nearest tiles on the tile's own
ring stay at least 150 m away, and it follows each attitude law around the orbit.

Schemes, evaluated along the twelve kept orbits of ring_keeping_run as photon_control does:
- as built: radial-facing and shingled;
- balanced overlap: the redirecting filter covers only the middle of each tile along the ring, the 8.5 km pitch plus
  a stated overlap, and the overhanging edges are bare film that passes the light, so a tile's filter is shaded only
  across that overlap; the overlap is also evaluated with the neighbour spacing moved by the tide's swing in a free
  string (ring_bundle's lone ring);
- moving mass: a mass that carries the tile's centre of mass toward the point where radiation, gravity and the
  attitude's own torque cancel, within a stated reach;
- rolled outside service: tiles roll about the line of their centres away from the radial-facing attitude outside
  the service arc, to the side that clearance allows (one side on the day half, the other on the night half), and
  pass through radial-facing at the nodes; a one-sided roll at the angle of largest side force measures the plane
  steering that rolling gives;
- edge-on turns in pitch outside service (the earlier proposal), whose push is stated for comparison.

Torques. Radiation torque from the visible filter area's first moments (8 solar points, the bundle's shadows and
Moon and Earth eclipses) and gravity-gradient torque from a 3x3 quadrature (tile_torque), in the bundle's common
frame turned to each law's attitude; the torque that the commanded attitude itself needs, I dw/dt + w x Iw, from the
commanded frames of each tile. Trim across a tile's halves supplies up to a stated fraction of the light's pressure
on half the filter at a quarter-side lever; a momentum store takes the rest and trim unloads it.

Push. Gauss's equation for the eccentricity vector, de/dt = (2(v.f) r - (r.f) v - (r.v) f)/mu, averaged over each
orbit and taken in the frame of the Sun's direction projected into the ring plane; the forced eccentricity of
frozen_rings scales with the drive across the Sun line. The plane's turn, the radius-vector rotation f_h r/|h|, is
expressed as the height change of the sunward crossing.

Gravity and attitude torques scale with areal mass and the push inversely; they are evaluated for the kept run's
tiles (62.7 g/m2) and restated for 26 and 5 g/m2.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.attitude_schemes
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
from protection.dynamics import ring_bundle as rb
from protection.dynamics.eclipse_parallel import moments as masked_moments
from protection.dynamics.ephemeris import DEFAULT_KERNEL
from protection.dynamics.fast_parallel import moments as parallel_moments
from protection.dynamics.fleet import sun_points
from protection.dynamics.natural_pattern import retarded_sun
from protection.dynamics.optical import length, unit
from protection.dynamics.tile_torque import gravity_torque
from .frozen_rings import OUT as FROZEN_OUT
from .ring_bundle_run import DESIGN as BUNDLE, OUT as BUNDLE_OUT, SEED, environment, initial
from .ring_keeping_run import DESIGN as KEPT, FILES as KEPT_FILES, OUT as KEPT_OUT, kept

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE/'results/attitude_schemes.json'
FILES = list(dict.fromkeys(['research/studies/solar_shield_array/attitude_schemes.py',
                            'protection/dynamics/tile_torque.py', 'research/studies/solar_shield_array/frozen_rings.py']
                           + KEPT_FILES))
DESIGN = dict(every=2, suns=8, rotation=.317, clearance_m=150., neighbour_rings=3, slip_step_m=10.,
              envelope_phases_deg=[0., 25., 45., 60., 80., 90., 100., 135., 180.], service_half_angle_deg=25.,
              filter_overlaps_m=[200., 400., 800.], design_overlap_m=400., roll_deg=[30., 60., 80.],
              roll_ramp_deg=30., steer_roll_deg=35.26, ballast_reach_m=[250., 500., 1000., 1500., 2000., 3000.],
              ballast_travel_m=4500., trim_levels=[.1, .25], unload_s=3600., areal_masses_kg_m2=[.005, .026],
              turn_deg=65., turn_torques_N_m=[3e4, 1e5, 1e6])
HOUR = 3600.


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def axis_rotate(x, axis, angle):
    """Vectors x (..., 3) rotated about unit axes (..., 3) by angles (...)."""
    c, s = np.cos(angle)[..., None], np.sin(angle)[..., None]
    return x*c+np.cross(axis, x)*s+axis*np.sum(axis*x, axis=-1, keepdims=True)*(1-c)


def body_frames(states, tilt, roll=0., pitch=0.):
    """Axes (..., 3, 3), columns u, v, n: ring_bundle's radial-facing shingled frame pitched about the orbit normal
    and then rolled about the line of tile centres (the velocity), which turns a ring's tiles rigidly together."""
    states = np.asarray(states, float)
    f = rb.tile_frames(states, tilt)
    t, h = unit(states[..., 3:6]), f[..., :, 1]
    roll = np.broadcast_to(np.asarray(roll, float), f.shape[:-2])
    pitch = np.broadcast_to(np.asarray(pitch, float), f.shape[:-2])
    cols = [axis_rotate(axis_rotate(f[..., :, k], h, pitch), t, roll) for k in range(3)]
    return np.stack(cols, axis=-1)


def phase(sun, states):
    """Angle of each state from the Sun's direction projected into its orbit plane, toward the motion (rad)."""
    h = unit(np.cross(states[..., :3], states[..., 3:6]))
    s = unit(sun-states[..., :3])
    sp = unit(s-np.sum(s*h, axis=-1, keepdims=True)*h)
    w = np.cross(h, sp)
    return np.arctan2(np.sum(states[..., :3]*w, axis=-1), np.sum(states[..., :3]*sp, axis=-1))


def ramp(x, a, b):
    """Zero below a, one above b, half a cosine between."""
    z = np.clip((np.asarray(x, float)-a)/(b-a), 0., 1.)
    return .5-.5*np.cos(np.pi*z)


def rolled(u, chi, service, width, day_sign=-1.):
    """Roll (rad) at phase u: zero in the service arc and at the nodes, day_sign*chi on the day half and the opposite
    on the night half, with cosine ramps of the given width (degrees) out of service and through the nodes."""
    a = np.abs(np.angle(np.exp(1j*np.asarray(u, float))))
    s, w = np.radians(service), np.radians(width)
    day = day_sign*chi*ramp(a, s, s+w)*ramp(np.pi/2-a, 0., w)
    night = -day_sign*chi*ramp(a-np.pi/2, 0., w)
    return np.where(a <= np.pi/2, day, night)


def steering(u, chi, service, width, day_sign=-1.):
    """One-sided roll for plane steering: the rolled law on the half orbit after service only (0 < u < pi), ramped
    back to radial-facing before the anti-solar point."""
    u = np.angle(np.exp(1j*np.asarray(u, float)))
    return np.where(u > 0, rolled(u, chi, service, width, day_sign)*ramp(np.pi-u, 0., np.radians(width)), 0.)


def edge_on_pitch(u, service):
    """Pitch (rad) that turns a radial-facing tile edge-on to the Sun outside the service arc, in the ring plane."""
    u = np.angle(np.exp(1j*np.asarray(u, float)))
    s = np.radians(service)
    return np.where(np.abs(u) <= s, 0., u-np.sign(u)*np.pi/2)


def square_gap(axes, d, side):
    """Surface distance between parallel squares with shared axes (rows e1, e2, n) and centre offsets d (..., 3)."""
    x = np.abs(d@axes.T)
    return np.sqrt(np.maximum(x[..., 0]-side, 0.)**2+np.maximum(x[..., 1]-side, 0.)**2+x[..., 2]**2)


def worst_gap(u, roll, pitch, radius, step, height, pitch_m, tilt, rings, slip_step, side=rb.SIDE):
    """Smallest surface distance from a tile at phase u to its own ring's two nearest tiles and, at every along-track
    offset, to every tile of the nearest rings on each side, all sharing the attitude (local frame: radial, along,
    normal)."""
    state = np.array([radius, 0., 0., 0., 1., 0.])
    f = body_frames(state, tilt, roll, pitch)
    axes = f.T
    gaps = [square_gap(axes, np.array([[0., j*pitch_m, 0.]]), side).min() for j in (1, 2)]
    slips = np.arange(-4*side, 4*side+slip_step, slip_step)
    for m in range(1, rings+1):
        d = np.c_[np.full_like(slips, m*step), slips, np.full_like(slips, m*height*np.cos(u))]
        gaps.append(square_gap(axes, d, side).min())
    return float(min(gaps))


def envelope(radius, step, height, tilt):
    """Allowed pitch and roll offsets (degrees, contiguous spans) at the stated phases."""
    rows = []
    offsets = np.arange(-90., 91., 1.)
    for u in DESIGN['envelope_phases_deg']:
        row = dict(phase_deg=u)
        for name in ('pitch', 'roll'):
            ok = [o for o in offsets if worst_gap(np.radians(u), np.radians(o)*(name == 'roll'),
                                                  np.radians(o)*(name == 'pitch'), radius, step, height,
                                                  BUNDLE['pitch_m'], tilt, DESIGN['neighbour_rings'],
                                                  DESIGN['slip_step_m']) >= DESIGN['clearance_m']]
            spans, start = [], None
            for o in offsets:
                if o in ok and start is None:
                    start = o
                if o not in ok and start is not None:
                    spans.append([start, o-1.]); start = None
            if start is not None:
                spans.append([start, offsets[-1]])
            row[f'{name}_allowed_deg'] = spans
        rows.append(row)
    return rows


def law_gap(law, radius, step, height, tilt):
    """Smallest clearance along the orbit for a law giving (roll, pitch) at each phase, and where it falls."""
    phases = np.radians(np.arange(-180., 180., 1.))
    gaps = [worst_gap(u, *law(u), radius, step, height, BUNDLE['pitch_m'], tilt, DESIGN['neighbour_rings'],
                      DESIGN['slip_step_m']) for u in phases]
    k = int(np.argmin(gaps))
    return dict(min_gap_m=float(gaps[k]), at_phase_deg=float(np.degrees(phases[k])),
                share_of_orbit_below_clearance=float(np.mean(np.array(gaps) < DESIGN['clearance_m'])))


def shade(env, t, q, frame, sample, bodies, filter_u, filter_v):
    """Visible filter area and first moments per tile for each solar point, for filter rectangles filter_u by
    filter_v centred on each tile. Parallel shadows are exact under the in-plane scaling that turns the rectangles
    into squares; inside a body eclipse the full clear square is used."""
    out = []
    for source in sun_points(retarded_sun(sample), DESIGN['suns'], DESIGN['rotation']):
        if bodies is None:
            a, b = rb.CLEAR/filter_u, rb.CLEAR/filter_v
            area, mx, my = parallel_moments(q, frame@np.diag([a, b, 1.]), source, rb.CLEAR).T
            out.append((source, area/(a*b), mx/(a*a*b), my/(a*b*b), filter_u*filter_v))
        else:
            area, mx, my = masked_moments(q, frame, source, bodies, [K.MOON_RADIUS, K.EARTH_RADIUS]).T
            out.append((source, area, mx, my, rb.CLEAR**2))
    return out


def load(q, normals, shaded, fraction, filter_u, filter_v):
    """Sail force (N), radiation torque in the common axes (N m), normal force (N) and the torque that a full trim
    across the filter's halves would give (N m, per in-plane axis) on each tile."""
    force, torque = np.zeros_like(q), np.zeros_like(q)
    normal, cap = np.zeros(len(q)), np.zeros((len(q), 2))
    for source, area, mx, my, full in shaded:
        ray = source-q
        cosine = np.sum(unit(ray)*normals, axis=1)
        density = K.SOLAR_CONSTANT*(K.AU/length(ray))**2/DESIGN['suns']*np.abs(cosine)
        pressure = -2*fraction/K.SPEED_OF_LIGHT*density*cosine
        force += (pressure*area)[:, None]*normals
        torque += pressure[:, None]*np.c_[my, -mx, np.zeros(len(q))]
        normal += pressure*area
        lit = area/full
        cap += (np.abs(pressure)*lit)[:, None]*np.c_[filter_u*filter_v/2*filter_v/4, filter_u*filter_v/2*filter_u/4]
    return force, torque, normal, cap


def angular(frames, t):
    """Body angular velocity and acceleration (rad/s, rad/s^2) of a sequence of frames (n_t, ..., 3, 3)."""
    dR = np.gradient(frames, t, axis=0)
    W = np.einsum('t...ji,t...jk->t...ik', frames, dR)
    w = np.stack([W[..., 2, 1], W[..., 0, 2], W[..., 1, 0]], axis=-1)
    return w, np.gradient(w, t, axis=0)


def store(torque, cap, dt, trim):
    """Peak stored momentum per tile (N m s) when trim up to trim*cap takes each in-plane torque and unloads the store
    with the stated time constant; the normal axis has no trim. torque: (n_t, tiles, 3), cap: (n_t, tiles, 2)."""
    held = np.zeros(torque.shape[1:])
    peak = np.zeros(torque.shape[1])
    for j in range(len(dt)):
        want = -torque[j, :, :2]-held[:, :2]/DESIGN['unload_s']
        applied = np.clip(want, -trim*cap[j], trim*cap[j])
        held = held+(torque[j]+np.c_[applied, np.zeros(len(applied))])*dt[j]
        peak = np.maximum(peak, np.linalg.norm(held, axis=1))
    return peak


def drive(states, force, mass, period, t):
    """Orbit-averaged eccentricity-vector rate (1/s) per tile in the Sun frame (along, across the Sun line) and the
    rate of change of the sunward crossing's height (m/s) from the plane's turn."""
    r, v = states[..., :3], states[..., 3:6]
    f = force/mass[None, :, None]
    de = (2*np.sum(v*f, axis=-1)[..., None]*r-np.sum(r*f, axis=-1)[..., None]*v
          - np.sum(r*v, axis=-1)[..., None]*f)/K.MOON_GM
    h = np.cross(r, v)
    hn = np.linalg.norm(h, axis=-1)[..., None]
    turn = np.sum(f*h, axis=-1)[..., None]/hn*r/hn  # plane rotation rate vector: f_h r/|h|
    return de, turn


def main():
    wall, cpu0 = time.monotonic(), time.process_time()
    env = environment(KEPT['orbits']*2.0+2)
    states, ring, slot, ratio, built = initial(env, KEPT['radius_step_m'])
    tilt = np.radians(BUNDLE['tilt_deg'])
    reps = states[slot == (BUNDLE['per_ring']-1)//2]
    lags = BUNDLE['pitch_m']/np.linalg.norm(reps[:, 3:], axis=1)
    model = rb.ContinuousRings(env, tilt, np.full(len(reps), ratio), lags, reach=KEPT['copy_reach'])
    origin = unit(env.at(0.)['positions']['sun'])
    t, s, _, _, period, _, key, _ = kept(env, reps, model, origin)
    mass = model.forces.mass[:len(reps)]
    sigma = float(env.sigma*ratio)
    fraction = env.central_fraction
    count = len(reps)
    mid = count//2
    edges = np.isin(np.arange(count), [0, count-1])
    interior = ~edges
    radius = float(np.linalg.norm(reps[mid, :3]))
    service = DESIGN['service_half_angle_deg']
    width = DESIGN['roll_ramp_deg']

    # Clearance: the envelope and each law's smallest gap along the orbit.
    geom = (radius, KEPT['radius_step_m'], built['height_pitch_m'], tilt)
    clearance = dict(envelope=envelope(*geom))
    probe = worst_gap(0., np.radians(-30.), 0., *geom[:3], BUNDLE['pitch_m'], tilt, 1, DESIGN['slip_step_m'])
    day_sign = -1. if probe >= DESIGN['clearance_m'] else 1.
    laws = {'as_built': lambda u: (0., 0.)}
    for chi in DESIGN['roll_deg']:
        laws[f'rolled_{chi:g}_deg'] = (lambda c: lambda u: (float(rolled(u, np.radians(c), service, width, day_sign)), 0.))(chi)
    laws['steering'] = lambda u: (float(steering(u, np.radians(DESIGN['steer_roll_deg']), service, width, day_sign)), 0.)
    clearance['laws'] = {name: law_gap(law, *geom) for name, law in laws.items()}
    clearance['laws']['edge_on_pitch'] = law_gap(lambda u: (0., float(edge_on_pitch(u, service))), *geom)
    clearance['day_roll_sign'] = day_sign
    print('clearance', json.dumps(clearance['laws']), flush=True)

    # Schemes along the kept orbits.
    overlaps = sorted(set(DESIGN['filter_overlaps_m']+[DESIGN['design_overlap_m']+d for d in (-160., 160.)]))
    swing = json.loads(BUNDLE_OUT.read_text())['lone_ring_gravity']['along_spacing_range_m']
    schemes = {'as_built': dict(roll=None, filter=(rb.CLEAR, rb.CLEAR))}
    for d in overlaps:
        schemes[f'balanced_{d:g}_m'] = dict(roll=None, filter=(BUNDLE['pitch_m']+d, rb.CLEAR))
    for chi in DESIGN['roll_deg']:
        schemes[f'rolled_{chi:g}_deg'] = dict(roll=lambda u, c=chi: rolled(u, np.radians(c), service, width, day_sign),
                                              filter=(rb.CLEAR, rb.CLEAR))
    schemes['steering'] = dict(roll=lambda u: steering(u, np.radians(DESIGN['steer_roll_deg']), service, width,
                                                       day_sign), filter=(rb.CLEAR, rb.CLEAR))
    names = list(schemes)
    sun_at = np.stack([retarded_sun(env.at(x)) for x in t])
    u_all = phase(sun_at[:, None, :], s)
    # Commanded frames at every sample for the attitude's own torque.
    commanded = {}
    for name, sch in schemes.items():
        roll = 0. if sch['roll'] is None else sch['roll'](u_all)
        commanded[name] = body_frames(s, tilt, roll)
    pick = np.arange(0, len(t), DESIGN['every'])
    dt = np.gradient(t[pick])
    rec = {name: dict(force=[], torque=[], normal=[], cap=[], gravity=[]) for name in names}
    for i in pick:
        state = s[i]
        everything = model.copies(state)[0]
        q = everything[:, :3]
        sample = env.at(t[i])
        bodies = rb.eclipse_bodies(q, sample)
        u_c = phase(retarded_sun(sample), everything.mean(axis=0))
        for name, sch in schemes.items():
            roll = 0. if sch['roll'] is None else float(sch['roll'](u_c))
            common = body_frames(everything.mean(axis=0), tilt, roll)
            normals = body_frames(everything, tilt, roll)[:, :, 2]
            fu, fv = sch['filter']
            shaded = shade(env, t[i], q, common, sample, bodies, fu, fv)
            force, torque, normal, cap = load(q, normals, shaded, fraction, fu, fv)
            r = rec[name]
            r['force'].append(force[:count]); r['torque'].append(torque[:count]); r['normal'].append(normal[:count])
            r['cap'].append(cap[:count]); r['gravity'].append(gravity_torque(env, t[i], q[:count], common)*mass[:, None])
    for name in names:
        rec[name] = {k: np.array(v) for k, v in rec[name].items()}
    inertia = np.array([1., 1., 2.])*rb.SIDE**2/12
    orbit = np.minimum((t[pick]/period).astype(int), KEPT['orbits']-1)

    def stats(x, mask=interior):
        return dict(median=float(np.median(x[:, mask])), p95=float(np.percentile(x[:, mask], 95)),
                    max=float(np.max(x[:, mask])))

    def orbit_mean(x):
        return np.array([[np.average(x[orbit == k, j], weights=dt[orbit == k], axis=0) for j in range(count)]
                         for k in range(KEPT['orbits'])])

    results, required_all = {}, {}
    base_drive = None
    for name in names:
        r = rec[name]
        w, alpha = angular(commanded[name], t)
        own = (alpha*inertia+np.cross(w, w*inertia))[pick]*mass[None, :, None]
        required = own-r['gravity']-r['torque']
        required_all[name] = required
        de, turn = drive(s[pick], r['force'], mass, period, t[pick])
        sun_p = sun_at[pick][:, None, :]
        h = unit(np.cross(s[pick][..., :3], s[pick][..., 3:6]))
        sp = unit((sun_p-s[pick][..., :3])-np.sum((sun_p-s[pick][..., :3])*h, axis=-1, keepdims=True)*h)
        wdir = np.cross(h, sp)
        de_sun = np.stack([np.sum(de*sp, axis=-1), np.sum(de*wdir, axis=-1)], axis=-1)
        per_orbit = orbit_mean(de_sun)
        height_rate = -radius*np.sum(turn*wdir, axis=-1)
        strip_rate = np.sum(turn*sp, axis=-1)
        mean_drive = np.mean(np.linalg.norm(per_orbit[:, interior], axis=-1))
        if name == 'as_built':
            base_drive = per_orbit
        rad_mean, req_mean = orbit_mean(r['torque']), orbit_mean(required)
        lit = r['cap'][..., 0] > 0
        trim_need = np.max(np.abs(required[..., :2])/np.where(r['cap'] > 0, r['cap'], np.inf), axis=-1)
        cop = np.c_[(-r['torque'][..., 1]/np.where(np.abs(r['normal']) > 0, r['normal'], np.inf)).ravel(),
                    (r['torque'][..., 0]/np.where(np.abs(r['normal']) > 0, r['normal'], np.inf)).ravel()]
        strong = (np.abs(r['normal']) > .05*np.abs(r['normal']).max()).ravel()
        row = dict(
            sail_force_N=stats(np.linalg.norm(r['force'], axis=-1)),
            radiation_torque_N_m=stats(np.linalg.norm(r['torque'], axis=-1)),
            radiation_orbit_mean_N_m=dict(roll=float(np.median(np.abs(rad_mean[:, interior, 0]))),
                                          pitch=float(np.median(np.abs(rad_mean[:, interior, 1])))),
            gravity_torque_N_m=stats(np.linalg.norm(r['gravity'], axis=-1)),
            attitude_torque_N_m=stats(np.linalg.norm(own, axis=-1)),
            required_torque_N_m=stats(np.linalg.norm(required, axis=-1)),
            required_orbit_mean_N_m=dict(roll=float(np.median(np.abs(req_mean[:, interior, 0]))),
                                         pitch=float(np.median(np.abs(req_mean[:, interior, 1]))),
                                         normal=float(np.median(np.abs(req_mean[:, interior, 2])))),
            pressure_centre_offset_m=dict(median=float(np.median(np.linalg.norm(cop[strong], axis=1))),
                                          p95=float(np.percentile(np.linalg.norm(cop[strong], axis=1), 95))),
            trim_needed_lit=dict(median=float(np.median(trim_need[lit & interior[None, :]])),
                                 p95=float(np.percentile(trim_need[lit & interior[None, :]], 95))),
            store_N_m_s={f'trim_{x:g}': dict(median=float(np.median(store(required, r['cap'], dt, x)[interior])),
                                              max=float(np.max(store(required, r['cap'], dt, x)[interior])))
                         for x in DESIGN['trim_levels']},
            eccentricity_drive=dict(per_s=float(mean_drive),
                                    over_as_built=float(mean_drive/np.mean(np.linalg.norm(base_drive[:, interior],
                                                                                          axis=-1)))),
            crossing_height_rate_km_per_day=dict(
                median=float(np.median(orbit_mean(height_rate[..., None])[:, interior, 0])*86.4),
                min=float(np.min(orbit_mean(height_rate[..., None])[:, interior, 0])*86.4),
                max=float(np.max(orbit_mean(height_rate[..., None])[:, interior, 0])*86.4)),
            strip_turn_deg_per_day=dict(
                median=float(np.degrees(np.median(orbit_mean(strip_rate[..., None])[:, interior, 0]))*86400),
                min=float(np.degrees(np.min(orbit_mean(strip_rate[..., None])[:, interior, 0]))*86400),
                max=float(np.degrees(np.max(orbit_mean(strip_rate[..., None])[:, interior, 0]))*86400)))
        scaled = {}
        for areal in DESIGN['areal_masses_kg_m2']:
            k = areal/sigma
            req = own*k-r['gravity']*k-r['torque']
            scaled[f'{areal*1e3:g}_g_m2'] = dict(
                store_N_m_s={f'trim_{x:g}': float(np.median(store(req, r['cap'], dt, x)[interior]))
                             for x in DESIGN['trim_levels']},
                required_orbit_mean_N_m=float(np.median(np.linalg.norm(orbit_mean(req)[:, interior, :2], axis=-1))))
        row['by_areal_mass'] = scaled
        results[name] = row
        print(name, json.dumps({k: row[k] for k in ('radiation_orbit_mean_N_m', 'required_orbit_mean_N_m',
                                                    'store_N_m_s', 'eccentricity_drive')}), flush=True)
    # The edge-on pitch law's push: the as-built force inside the service arc only.
    u_pick = u_all[pick]
    inside = np.abs(u_pick) <= np.radians(service)
    de, _ = drive(s[pick], rec['as_built']['force']*inside[..., None], mass, period, t[pick])
    sun_p = sun_at[pick][:, None, :]
    h = unit(np.cross(s[pick][..., :3], s[pick][..., 3:6]))
    sp = unit((sun_p-s[pick][..., :3])-np.sum((sun_p-s[pick][..., :3])*h, axis=-1, keepdims=True)*h)
    de_sun = np.stack([np.sum(de*sp, axis=-1), np.sum(de*np.cross(h, sp), axis=-1)], axis=-1)
    edge_drive = np.mean(np.linalg.norm(orbit_mean(de_sun)[:, interior], axis=-1))
    edge_on = dict(eccentricity_drive_over_as_built=float(edge_drive/np.mean(np.linalg.norm(base_drive[:, interior],
                                                                                           axis=-1))))
    # Turn time for a rest-to-rest pitch turn of the stated angle under a constant torque budget (bang-bang).
    turns = {}
    for areal in [sigma]+DESIGN['areal_masses_kg_m2']:
        I = areal*rb.SIDE**2*rb.SIDE**2/12
        turns[f'{areal*1e3:.3g}_g_m2'] = {f'{tau:g}_N_m': float(2*np.sqrt(np.radians(DESIGN['turn_deg'])*I/tau)/HOUR)
                                          for tau in DESIGN['turn_torques_N_m']}
    edge_on['turn_hours_by_torque'] = turns
    edge_on['clearance'] = clearance['laws']['edge_on_pitch']
    # Moving mass on the as-built tiles: reach of the centre of mass against the store left over.
    r = rec['as_built']
    w, alpha = angular(commanded['as_built'], t)
    own = (alpha*inertia+np.cross(w, w*inertia))[pick]*mass[None, :, None]
    ballast = {}
    for areal in [sigma]+DESIGN['areal_masses_kg_m2']:
        k = areal/sigma
        other = r['torque']+r['gravity']*k-own*k
        normal = np.where(np.abs(r['normal']) > 1e-3, r['normal'], np.inf)
        target = np.stack([-other[..., 1]/normal, other[..., 0]/normal], axis=-1)
        rows = {}
        for reach in DESIGN['ballast_reach_m']:
            c = np.clip(target, -reach, reach)
            residual = own*k-r['gravity']*k-r['torque']+np.stack([c[..., 1]*r['normal'], -c[..., 0]*r['normal'],
                                                                 np.zeros_like(r['normal'])], axis=-1)
            rows[f'{reach:g}_m'] = dict(
                ballast_over_tile_mass=float(reach/(DESIGN['ballast_travel_m']-reach)),
                store_N_m_s={f'trim_{x:g}': float(np.median(store(residual, r['cap'], dt, x)[interior]))
                             for x in DESIGN['trim_levels']},
                residual_orbit_mean_N_m=float(np.median(np.linalg.norm(orbit_mean(residual)[:, interior, :2],
                                                                       axis=-1))))
        need = np.linalg.norm(target, axis=-1)[np.abs(r['normal']) > .05*np.abs(r['normal']).max()]
        ballast[f'{areal*1e3:.3g}_g_m2'] = dict(target_offset_m=dict(median=float(np.median(need)),
                                                                     p95=float(np.percentile(need, 95))),
                                                by_reach=rows)
    frozen = json.loads(FROZEN_OUT.read_text())
    last = frozen['passes'][-1]
    forced = {}
    for areal in (.0627, .026, .005):
        e = [np.hypot(*row['centre']) for row in last if abs(row['sigma_kg_m2']-areal) < 1e-6]
        forced[f'{areal*1e3:g}_g_m2'] = [float(min(e)), float(max(e))]
    out = dict(schema='terluna.research.attitude-schemes/1',
               producer=dict(files={f: digest(ROOT/f) for f in FILES},
                             constants=constants_used([f for f in FILES if f.endswith('.py')]),
                             inputs={'results/ring_keeping.json': digest(KEPT_OUT), 'results/joint_seed.json': digest(SEED),
                                     'results/ring_bundle.json': digest(BUNDLE_OUT),
                                     'results/frozen_rings.json': digest(FROZEN_OUT), 'de440s.bsp': digest(DEFAULT_KERNEL)}),
               evidence=('Clearance: parallel 10 km squares sharing each attitude, the tile\'s own ring at 8.5 km '
                         'pitch with the 1.5 degree shingle, and the three nearest rings on each side at 1 km radius '
                         'steps and the bundle\'s height step times the cosine of the phase, at every along-track '
                         'offset in 10 m steps. Schemes: the twelve kept orbits of ring_keeping_run (run key in '
                         'kept_run_key), every second five-minute sample, with the bundle\'s copies, 8 solar points '
                         'and Moon and Earth eclipses; filter rectangles by exact in-plane scaling of the parallel '
                         'shadow union (full squares inside body eclipses); gravity torque by tile_torque\'s 3x3 '
                         'quadrature in the common frame; the attitude\'s own torque from finite differences of each '
                         'tile\'s commanded frames at every sample. Push: orbit-averaged Gauss rates in the Sun frame.'),
               reading_rule=('A scheme holds attitude where the orbit-mean required torque is within trim and the '
                             'store stays small; it is safe only where its law keeps the clearance. The forced '
                             'eccentricity of a scheme is frozen_rings\' value times its drive over the as-built '
                             'drive. Gravity and the attitude\'s own torque scale with areal mass; radiation torque '
                             'and trim do not.'),
               design=DESIGN, kept_run_key=key, period_h=period/HOUR, tile_areal_mass_kg_m2=sigma,
               radius_km=radius/1e3, height_step_m=built['height_pitch_m'], tide_spacing_range_m=swing,
               clearance=clearance, schemes=results, edge_on_pitch=edge_on, moving_mass=ballast,
               forced_eccentricity_as_built=forced,
               resources=dict(cpu_s=time.process_time()-cpu0, wall_s=time.monotonic()-wall,
                              max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, numerical_threads=1))
    OUT.write_text(json.dumps(out, indent=1)+'\n')
    print(json.dumps({k: out[k] for k in ('clearance', 'edge_on_pitch', 'resources')}, indent=1)[:4000])


if __name__ == '__main__':
    main()
