"""The ring stack's radii against the tide's turning of its planes (the outer rings; the closing step of 8 October).

plane_motion flew nine rings at one radius and found the outermost strips falling behind the Sun: a ring tilted
further from the Moon's orbit plane, about the line of nodes they all share, turns more slowly under Earth's tide.
The full disk is about 1,860 nested rings, and rings that share that line and a radius cross at its two ends, so
every ring needs a radius of its own: the stack spans a couple of thousand kilometres of radius, and the tide turns a
plane faster the further out it is. The layout of radii across the stack decides whether every plane turns with the
Sun. This runner

1. calibrates the turning: rings on their forced eccentricity (26 g/m2 tiles, radial-facing, the interior shadow
   factor, frozen_rings' sail force) at tilts of -26 to +26 degrees from the Moon's orbit plane and at radii of
   18,000-22,000 km fly a year in DE440s point-mass gravity, two passes as plane_motion does; each ring's turning
   against the Sun is the yearly trend of the azimuth of its pole about the Moon's orbit pole, less the Sun's, and
   the radius at which each height's ring turns with the Sun is found from it;
2. lays out the radii of the full stack, keeping the clearance pair by pair: neighbouring strips overlap like shingles
   at apolune and perilune alike, while rings far apart in the stack meet only at the nodes. With every ring on its
   forced eccentricity ('forced'), which grows with tilt and radius toward the stack's edges, a strip further out
   swings in at perilune and the steps run away; with keeping holding the eccentricity vectors of rings that come
   close together ('held'), as the bundle keeping does, they stay near the clearance. The layouts: one radius for
   all (plane_motion's idealisation, which no stack can fly); a one-way staircase, as the bundle studies flew it; a
   two-sided layout, ranked by height from the middle with the sides alternating; and a matched layout, ranked by
   the radius at which each ring turns with the Sun, which grows toward both edges with the sides interleaving as
   their turning asks. Each ring's turning and forced eccentricity come from the calibration; the photon steering it
   would need is its strip's turn (the turning times the sine of its tilt) against the roll's reach (plane_motion's
   scaling of attitude_schemes). For the held matched layout, the eccentricities keeping would have to hold so that
   strips on their own eccentricities keep the clearance, and what holding them does to each plane's turning (the
   tide's secular node rate, as (1 + 4 e^2)/sqrt(1 - e^2) with the apse across the node line), are estimated, and
   the layout is laid out again with that change, pass after pass ('matched_coupled');
3. flies a sample of the matched layout's rings for a year on their forced eccentricities and measures them as
   plane_motion does: the turning at each ring's layout radius, the common shift and rotation, the change in overlap
   between neighbouring strips and the steering each ring would need.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.ring_layout          # about 25 CPU minutes
    OPENBLAS_NUM_THREADS=1 python -m research.studies.solar_shield_array.ring_layout --smoke  # eight days, a few rings

Each pass is checkpointed in research/runs/solar_shield_array/ring_layout under the flight code's version, so a
stopped run resumes at the next and changes to the layouts leave the flights in place.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import resource
import time
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.ephemeris import DEFAULT_KERNEL
from protection.dynamics.optical import covering_radius, unit
from . import frozen_rings as fr
from .frozen_rings import FILES as FROZEN_FILES, OUT as FROZEN_OUT, fly, geometry, measure, start
from .plane_motion import ATTITUDE, DESIGN as PLANES, OUT as PLANES_OUT, passages, split, trend
from .ring_screen import environment
from .run import SCENARIO

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE/'results/ring_layout.json'
RUN = ROOT/'research/runs/solar_shield_array/ring_layout'
FILES = ['research/studies/solar_shield_array/ring_layout.py', 'research/studies/solar_shield_array/plane_motion.py',
         'research/studies/solar_shield_array/run.py', 'research/studies/solar_shield_array/scenario.json'] + FROZEN_FILES
DESIGN = dict(plane='moon_orbit', sigma_kg_m2=.026, passes=2, middle_radius_km=20000.,
              calibration_tilts_deg=[-26., -19.5, -13., -6.5, 0., 6.5, 13., 19.5, 26.],
              calibration_radii_km=[18000., 18800., 19600., 20400., 21200., 22000.],
              height_step_m=PLANES['height_step_m'], common_shift_share=.2255,
              clearance_m=[600., 1000.], clearance_design_m=1000.,
              sample_rings=13, width_fraction=PLANES['width_fraction'], grid_days=PLANES['grid_days'],
              series_every_days=10, untilted_deg=2., matched_stretches=[.7, .85, 1.],
              runaway_margin_km=1000., sun_rate_deg_per_day=360./K.JULIAN_YEAR_DAYS, coupled_passes=6)
SMOKE = dict(days=8., calibration_tilts_deg=[-19.5, 0., 19.5], calibration_radii_km=[19400., 20600.], sample_rings=5,
             grid_days=[1., 7.])
DAY = K.JULIAN_DAY
# The version of this module's flight code (flights, turning) that the checkpoints were flown with: the module's
# digest when its flights were first run. Change it whenever that code changes; test_ring_layout pins the sources.
FLIGHT_CODE = '902004710fded91c26c61fb251db5da24bf6ec417bc43ce92254b1773aa44bf8'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def frozen_centre(frozen, radius_km, tilt_deg):
    """frozen_rings' 26 g/m2 centre in the Moon-orbit plane at the nearest radius it flew, interpolated in the size of
    the tilt between its 0 and 20 degrees (as plane_motion starts its rings)."""
    near = min({row['radius_km'] for row in frozen}, key=lambda x: abs(x-radius_km))
    ref = {row['tilt_deg']: np.array(row['centre']) for row in frozen
           if row['plane'] == DESIGN['plane'] and row['radius_km'] == near
           and abs(row['sigma_kg_m2']-DESIGN['sigma_kg_m2']) < 1e-9}
    return ref[0.]+(ref[20.]-ref[0.])*abs(tilt_deg)/20.


def turning(env, times, traj, every=6):
    """Each ring's pole azimuth about the Moon's orbit pole less the Sun's (degrees, unwrapped), sampled every
    `every` samples, and its yearly trend (degrees a day): the rate at which the tide turns the plane against the
    Sun."""
    idx = np.arange(0, len(times), every)
    out = []
    poles, suns = [], []
    for t in times[idx]:
        sample = env.at(t)
        poles.append(unit(np.cross(-sample['positions']['earth'], sample['moon_v']-sample['earth_v'])))
        suns.append(sample['positions']['sun'])
    poles, suns = np.array(poles), np.array(suns)
    s = unit(suns-np.sum(suns*poles, axis=1)[:, None]*poles)
    for j in range(traj.shape[1]):
        r, v = traj[idx, j, :3], traj[idx, j, 3:]
        h = unit(np.cross(r, v))
        hp = h-np.sum(h*poles, axis=1)[:, None]*poles
        angle = np.unwrap(np.arctan2(np.sum(np.cross(s, hp)*poles, axis=1), np.sum(s*hp, axis=1)))
        days = times[idx]/DAY
        out.append(dict(rate_deg_per_day=float(np.degrees(np.polyfit(days, angle, 1)[0])),
                        swing_deg=float(np.degrees(np.ptp(angle-np.polyval(np.polyfit(days, angle, 1), days)))),
                        lean_deg=float(np.degrees(np.mean(np.arcsin(np.clip(np.linalg.norm(hp, axis=1), 0, 1)))))))
    return out


def flights(env, sun, poles, cases, tag, frozen):
    """Fly the cases for DESIGN['passes'] passes, each starting on the centres the previous one found; checkpointed.
    Returns the last pass's rows (measure, turning) and crossings, and the function evaluations used."""
    nfev = 0
    code = ''.join(digest(ROOT/f) for f in FROZEN_FILES)+FLIGHT_CODE
    for case in cases:
        if 'e0' not in case:
            case['e0'] = [float(v) for v in frozen_centre(frozen, case['radius_km'], case['tilt_deg'])]
    for k in range(DESIGN['passes']):
        key = hashlib.sha256(json.dumps([[c['radius_km'], c['tilt_deg'], c['e0'], c['width_m']] for c in cases]).encode()
                             + json.dumps({n: fr.DESIGN[n] for n in ('days', 'sample_s', 'rtol', 'atol', 'max_step_s',
                                                                      'tile_tilt_deg', 'interior_light')}).encode()
                             + code.encode()).hexdigest()[:16]
        path = RUN/f'{tag}_pass_{k}_{key}.json'
        if path.exists():
            saved = json.loads(path.read_text())
        else:
            states = np.array([start(sun, poles[DESIGN['plane']], c['radius_km'], c['tilt_deg'], *c['e0'])[0]
                               for c in cases])
            times, traj, used = fly(env, states, [DESIGN['sigma_kg_m2']]*len(cases))
            nfev += used
            rows = [measure(env, times, traj[:, j], dict(c, pole=poles[DESIGN['plane']], ecliptic=poles['ecliptic']))
                    for j, c in enumerate(cases)]
            turns = turning(env, times, traj)
            crossings = [[{n: (v.tolist() if isinstance(v, np.ndarray) else v) for n, v in row.items()} for row in ring]
                         for ring in passages(env, times, traj, poles['ecliptic'], [c['width_m'] for c in cases])]
            saved = dict(rows=rows, turns=turns, crossings=crossings, nfev=used)
            RUN.mkdir(parents=True, exist_ok=True)
            tmp = path.with_suffix('.tmp')
            tmp.write_text(json.dumps(saved, default=float))
            tmp.rename(path)
            del traj
        print(f'{tag} pass {k}', flush=True)
        for c, r, t in zip(cases, saved['rows'], saved['turns']):
            print(f"  {c['radius_km']:8.1f} km tilt {c['tilt_deg']:+7.2f} e0 {np.hypot(*c['e0']):.4f} centre "
                  f"{np.hypot(*r['centre']):.4f} circle {r['circle_radius']:.4f} turning {t['rate_deg_per_day']:+.5f} "
                  f"deg/d", flush=True)
        if k < DESIGN['passes']-1:
            for c, r in zip(cases, saved['rows']):
                c['e0'] = [float(v) for v in r['centre']]
    return saved, nfev


class Calibration:
    """The turning (degrees a day against the Sun) and forced eccentricity by radius and tilt, interpolated."""

    def __init__(self, radii, tilts, rate, ecc):
        self.radii, self.tilts = np.asarray(radii), np.asarray(tilts)
        self.rate = RegularGridInterpolator((self.radii, self.tilts), np.asarray(rate), bounds_error=False,
                                            fill_value=None)
        self.ecc = RegularGridInterpolator((self.radii, self.tilts), np.asarray(ecc), bounds_error=False,
                                           fill_value=None)

    def at(self, radius_km, tilt_deg):
        """Turning and forced eccentricity, with radius and tilt held to the calibrated ranges."""
        radius = np.clip(np.atleast_1d(radius_km), self.radii[0], self.radii[-1])
        tilt = np.clip(np.atleast_1d(tilt_deg), self.tilts[0], self.tilts[-1])
        p = np.c_[radius, tilt]
        return self.rate(p), np.clip(self.ecc(p), 0., .6)

    def covers(self, radii, tilts):
        return bool(np.all((radii >= self.radii[0]) & (radii <= self.radii[-1]))
                    and np.all(np.abs(tilts) <= max(abs(self.tilts[0]), abs(self.tilts[-1]))))


def reach(radius_km):
    """plane_motion's photon steering at 26 g/m2: the strip turn (degrees a day) the roll gives at this radius."""
    attitude = json.loads(ATTITUDE.read_text())
    steer = attitude['schemes']['steering']
    base, ratio = attitude['tile_areal_mass_kg_m2'], np.asarray(radius_km)/attitude['radius_km']
    return abs(steer['strip_turn_deg_per_day']['median'])*base/DESIGN['sigma_kg_m2']*ratio**.5


def stack_heights(radius_km):
    aperture = covering_radius(radius_km*1e3, K.AU, SCENARIO['protected_radii']*K.MOON_RADIUS)+SCENARIO['formation_margin_m']
    top = aperture*(1+DESIGN['common_shift_share'])
    return np.arange(-top, top+1e-6, DESIGN['height_step_m'])/1e3, aperture/1e3


def evaluate(cal, radii, heights, target):
    """Each ring's tilt, turning against the target (zero: with the Sun), forced eccentricity and steering need over
    reach."""
    tilts = np.degrees(np.arcsin(np.clip(heights/radii, -1, 1)))
    rate, ecc = cal.at(radii, tilts)
    error = rate-target
    need = np.abs(error)*np.abs(np.sin(np.radians(tilts)))/reach(radii)
    return tilts, error, ecc, need


def clear(radii, heights, ecc, clearance_km, held=False):
    """Smallest clearances: neighbouring strips at apolune and perilune, radius neighbours at the nodes (km). Held:
    keeping holds the eccentricity vectors of rings that come close together, so each pair shares its mean."""
    order = np.argsort(heights)
    r, e = radii[order], ecc[order]
    if held:
        common = .5*(e[1:]+e[:-1])
        peri = np.abs(np.diff(r))*(1-common)/(1+common)
    else:
        peri = np.abs(np.diff(r*(1-e)/(1+e)))
    shingle = min(np.abs(np.diff(r)).min(), peri.min())
    by_radius = np.argsort(radii)
    rs, es = radii[by_radius], ecc[by_radius]
    nodes = np.diff(rs)*(1-.5*(es[1:]+es[:-1])) if held else np.diff(rs*(1-es))
    return dict(shingle_km=float(shingle), nodes_km=float(nodes.min()),
                holds=bool(shingle >= clearance_km-1e-6 and nodes.min() >= clearance_km-1e-6))


def outward(cal, height, prior, constraints, cap):
    """The smallest radius at or beyond `prior` that meets each (radius, eccentricity, kind) constraint of rings already
    placed inside it, for a ring at this height; None beyond the cap. Kinds: 'strip' (a neighbouring strip, clear at
    apolune and perilune), 'node' (clear at the nodes), each 'forced' (both on their forced eccentricities) or 'held'
    (keeping holds them to a common one)."""
    trial = prior
    for _ in range(60):
        e = cal.at(trial, np.degrees(np.arcsin(height/trial)))[1][0]
        cands = [prior]
        for radius, ecc, kind, mode, clearance in constraints:
            if mode == 'held':
                common = .5*(e+ecc)
                cands.append(radius+clearance*((1+common)/(1-common) if kind == 'strip' else 1/(1-common)))
            elif kind == 'strip':
                cands += [radius+clearance, (radius*(1-ecc)/(1+ecc)+clearance)*(1+e)/(1-e)]
            else:
                cands.append((radius*(1-ecc)+clearance)/(1-e))
        new = max(cands)
        if new > cap:
            return None
        if abs(new-trial) < 1e-7:
            break
        trial = new
    return trial


def build(cal, heights, rank, want, clearance_km, mode, cap):
    """Place the rings in the order of `rank`, each at its wanted radius or the first radius beyond it that keeps the
    clearance from the ring placed before (at the nodes) and from its neighbouring strips placed already; the order
    makes every ring placed later lie further out. None if the stack runs past the cap."""
    order = np.argsort(heights)
    position = np.empty(len(heights), int)
    position[order] = np.arange(len(heights))
    radii = np.full(len(heights), np.nan)
    ecc = np.full(len(heights), np.nan)
    last = None
    for k in rank:
        constraints = []
        prior = want[k]
        if last is not None:
            prior = max(prior, radii[last]+1e-3)
            constraints.append((radii[last], ecc[last], 'node', mode, clearance_km))
        for nb in (position[k]-1, position[k]+1):
            if 0 <= nb < len(heights):
                j = order[nb]
                if not np.isnan(radii[j]):
                    constraints.append((radii[j], ecc[j], 'strip', mode, clearance_km))
        radius = outward(cal, heights[k], prior, constraints, cap)
        if radius is None:
            return None
        radii[k] = radius
        ecc[k] = cal.at(radius, np.degrees(np.arcsin(heights[k]/radius)))[1][0]
        last = k
    return radii


def staircase(cal, heights, start_km, clearance_km, up=True, mode='forced', cap=np.inf):
    """Radii rising from one edge to the other, each strip clear of the one before."""
    rank = np.argsort(heights) if up else np.argsort(-heights)
    return build(cal, heights, rank, np.full(len(heights), start_km), clearance_km, mode, cap)


def two_sided(cal, heights, start_km, clearance_km, mode='forced', cap=np.inf):
    """Rings ranked by their height from the middle, the two sides alternating, so the radius grows toward both
    edges."""
    rank = np.argsort(np.abs(heights)+1e-9*np.sign(heights), kind='stable')
    return build(cal, heights, rank, np.full(len(heights), start_km), clearance_km, mode, cap)


def ideal_radii(cal, heights):
    """The radius at which a ring at each height turns with the Sun, by bisection on the calibrated turning, held to
    the calibrated radii."""
    out = []
    lo0, hi0 = cal.radii[0], cal.radii[-1]
    for h in heights:
        turn = lambda r: cal.at(r, np.degrees(np.arcsin(h/r)))[0][0]
        if turn(lo0) >= 0:
            out.append(lo0)
            continue
        if turn(hi0) <= 0:
            out.append(hi0)
            continue
        lo, hi = lo0, hi0
        for _ in range(50):
            mid = .5*(lo+hi)
            lo, hi = (mid, hi) if turn(mid) < 0 else (lo, mid)
        out.append(.5*(lo+hi))
    return np.array(out)


def matched(cal, heights, ideal, offset, stretch, clearance_km, mode='held', cap=np.inf):
    """Rings ranked by the radius at which each turns with the Sun, each placed at that radius (scaled about the
    stack's median by `stretch` and moved by `offset`) or as far beyond as the clearance demands: the radius grows
    toward both edges, the sides interleaving as their turning asks."""
    rank = np.argsort(ideal, kind='stable')
    centre = float(np.median(ideal))
    return build(cal, heights, rank, centre+(ideal-centre)*stretch+offset, clearance_km, mode, cap)


def best_layout(cal, heights, make, clearance_km, target, lo, hi, cap):
    """The layout whose worst steering need over reach is least over its starting radius (every 100 km between lo
    and hi, then every 10 km about the best); stacks that run past the cap are not taken. None if every one does."""
    best = None
    for starts in (np.arange(lo, hi+1e-6, 100.), None):
        if starts is None:
            if best is None:
                return None
            starts = np.arange(best[1]-100., best[1]+100.+1e-6, 10.)
        for s in starts:
            radii = make(cal, heights, s, clearance_km, cap)
            if radii is None:
                continue
            tilts, error, ecc, need = evaluate(cal, radii, heights, target)
            if best is None or need.max() < best[0]:
                best = (float(need.max()), float(s), radii)
    return best


def best_matched(cal, heights, ideal, clearance_km, target, mode, cap):
    best = None
    for stretch in DESIGN['matched_stretches']:
        for coarse in (True, False):
            offsets = (np.arange(-1500., 400.+1e-6, 100.) if coarse
                       else np.arange(best[1]-50., best[1]+50.+1e-6, 10.) if best and best[2] == stretch else [])
            for offset in offsets:
                radii = matched(cal, heights, ideal, offset, stretch, clearance_km, mode, cap)
                if radii is None:
                    continue
                need = evaluate(cal, radii, heights, target)[3]
                if best is None or need.max() < best[0]:
                    best = (float(need.max()), float(offset), stretch, radii)
    return best


def held_profile(radii, heights, forced, clearance_km):
    """Eccentricities that keeping could hold so that neighbouring strips on their own eccentricities keep the
    clearance at perilune: from the middle outward on each side, each ring's forced value moved no further than the
    clearance allows from its inner neighbour's held value."""
    order = np.argsort(heights)
    held = np.array(forced, float)
    middle = int(np.argmin(np.abs(heights[order])))
    for step in (1, -1):
        k = middle+step
        while 0 <= k < len(order):
            i, j = order[k-step], order[k]
            ri, rj, ei = radii[i], radii[j], held[i]
            pi = ri*(1-ei)/(1+ei)
            # perilune of j against i: rj (1-e)/(1+e) - pi >= c (j outer) or pi - rj (1-e)/(1+e) >= c (j inner)
            if rj > ri:
                q = (pi+clearance_km)/rj                      # (1-e)/(1+e) at least q
                upper = (1-q)/(1+q)
                held[j] = min(forced[j], upper)
            else:
                q = (pi-clearance_km)/rj                      # (1-e)/(1+e) at most q
                lower = (1-q)/(1+q) if q > 0 else 1.
                held[j] = max(forced[j], lower)
            k += step
    return held


def held_turning_change(ecc, held):
    """The change in each plane's turning (degrees a day) when keeping holds its eccentricity off the forced value:
    the tide's secular node rate with the apse across the node line scales as (1 + 4 e^2)/sqrt(1 - e^2)."""
    g = lambda e: (1+4*e*e)/np.sqrt(1-e*e)
    return DESIGN['sun_rate_deg_per_day']*(g(held)/g(ecc)-1)


def ideal_radii_shifted(cal, heights, shift):
    """As ideal_radii, with each height's turning moved by `shift` (degrees a day)."""
    out = []
    lo0, hi0 = cal.radii[0], cal.radii[-1]
    for h, d in zip(heights, shift):
        turn = lambda r: cal.at(r, np.degrees(np.arcsin(h/r)))[0][0]+d
        if turn(lo0) >= 0:
            out.append(lo0)
            continue
        if turn(hi0) <= 0:
            out.append(hi0)
            continue
        lo, hi = lo0, hi0
        for _ in range(50):
            mid = .5*(lo+hi)
            lo, hi = (mid, hi) if turn(mid) < 0 else (lo, mid)
        out.append(.5*(lo+hi))
    return np.array(out)


def coupled_matched(cal, heights, clearance_km, target):
    """The held matched layout with the turning change its held eccentricities bring (held_turning_change), iterated:
    each pass lays the rings out by the radius at which each turns with the Sun given the previous pass's change, then
    finds the held eccentricities and their change again, moving halfway to it. Returns the last layout, its held
    eccentricities and change, and each pass's figures."""
    shift = np.zeros(len(heights))
    passes = []
    for _ in range(DESIGN['coupled_passes']):
        ideal = ideal_radii_shifted(cal, heights, shift)
        best = None
        for stretch in DESIGN['matched_stretches']:
            for offset in np.arange(-1500., 400.+1e-6, 100.):
                radii = matched(cal, heights, ideal, offset, stretch, clearance_km, 'held')
                tilts, error, forced, _ = evaluate(cal, radii, heights, target)
                need = np.abs(error+shift)*np.abs(np.sin(np.radians(tilts)))/reach(radii)
                if best is None or need.max() < best[0]:
                    best = (float(need.max()), radii)
        radii = best[1]
        tilts, error, forced, _ = evaluate(cal, radii, heights, target)
        kept = held_profile(radii, heights, forced, clearance_km)
        change = held_turning_change(forced, kept)
        need = np.abs(error+change)*np.abs(np.sin(np.radians(tilts)))/reach(radii)
        passes.append(dict(radius_km=[float(radii.min()), float(radii.max())],
                           largest_departure=float(np.abs(kept-forced).max()),
                           turning_change_min_deg_per_day=float(change.min()),
                           steering_over_reach_max=float(need.max()), rings_over=int(np.sum(need > 1))))
        shift = .5*shift+.5*change
    return radii, kept, change, need, passes


def summary(cal, radii, heights, target, clearance_km, middle_km, held=False):
    tilts, error, ecc, need = evaluate(cal, radii, heights, target)
    return dict(within_calibration=cal.covers(radii, tilts),
                radius_km=dict(min=float(radii.min()), max=float(radii.max()), span=float(np.ptp(radii))),
                tiles_over_one_radius=float(np.mean(radii)/middle_km),
                turning_error_deg_per_day=dict(min=float(error.min()), max=float(error.max())),
                steering_over_reach=dict(max=float(need.max()), rings_over=int(np.sum(need > 1)),
                                         rings=int(len(radii)), share_over=float(np.mean(need > 1)),
                                         over_a_tenth=int(np.sum(need > .1))),
                clearance=clear(radii, heights, ecc, clearance_km, held),
                steps_km=dict(sorted_min=float(np.diff(np.sort(radii)).min()),
                              sorted_max=float(np.diff(np.sort(radii)).max())))


def sample(radii, heights, count):
    """Rings spread evenly across the stack's height, with their layout radii."""
    picks = np.unique(np.round(np.linspace(0, len(heights)-1, count)).astype(int))
    return [(float(radii[k]), float(heights[k])) for k in picks]


def measure_layout(env, sun, poles, rings, tag, frozen, cal, target):
    """Fly a layout's sample rings for a year and measure them as plane_motion does."""
    cases = []
    for radius, height in rings:
        tilt = float(np.degrees(np.arcsin(height/radius)))
        aperture = covering_radius(radius*1e3, K.AU, SCENARIO['protected_radii']*K.MOON_RADIUS)+SCENARIO['formation_margin_m']
        e = cal.at(radius, tilt)[1][0]
        base = frozen_centre(frozen, radius, tilt)
        scaled = base*e/max(np.hypot(*base), 1e-9)
        cases.append(dict(radius_km=radius, tilt_deg=tilt, design_height_m=height*1e3, aperture_m=float(aperture),
                          width_m=float(DESIGN['width_fraction']*np.sqrt(max(aperture**2-min(abs(height*1e3), aperture)**2,
                                                                             0.))),
                          e0=[float(v) for v in scaled]))
    saved, nfev = flights(env, sun, poles, cases, tag, frozen)
    grid = np.arange(DESIGN['grid_days'][0], DESIGN['grid_days'][1], 1.)
    heights = np.array([c['design_height_m'] for c in cases])
    y = np.zeros((len(grid), len(cases), 3)); r = np.zeros_like(y); xs = np.zeros((len(cases), 3))
    for a, ring in enumerate(saved['crossings']):
        for b, row in enumerate(ring):
            if len(row['t']) < 2:
                raise RuntimeError(f'{tag}: ring {a} crossed the aperture plane {b} fewer than twice')
            y[:, a, b] = np.interp(grid, np.asarray(row['t'])/DAY, row['y'])
            r[:, a, b] = np.interp(grid, np.asarray(row['t'])/DAY, row['r'])
            xs[a, b] = row['x']
    elevation = np.arcsin(np.clip(y/r, -1, 1))
    y_plane = r[0][None]*(np.sin(elevation)-np.sin(elevation[0])[None])
    c0, theta, residual, change = split(heights, xs, y_plane, DESIGN['height_step_m'])
    turns = np.array([t['rate_deg_per_day'] for t in saved['turns']])
    tilts = np.array([c['tilt_deg'] for c in cases])
    radii = np.array([c['radius_km'] for c in cases])
    need = np.abs(turns-target)*np.abs(np.sin(np.radians(tilts)))/reach(radii)
    rate = trend(grid, change)
    return dict(rings=[dict(radius_km=c['radius_km'], height_km=c['design_height_m']/1e3, tilt_deg=round(c['tilt_deg'], 3),
                            forced_eccentricity=round(float(np.hypot(*row['centre'])), 4),
                            eccentricity_circle=round(row['circle_radius'], 4),
                            turning_deg_per_day=round(t['rate_deg_per_day'], 5), steering_over_reach=round(float(n), 3))
                       for c, row, t, n in zip(cases, saved['rows'], saved['turns'], need)],
                common=dict(shift_km=dict(min=float(c0.min()/1e3), max=float(c0.max()/1e3)),
                            centred_shift_km=float(np.ptp(c0)/2e3),
                            rotation_deg=dict(min=float(np.degrees(theta.min())), max=float(np.degrees(theta.max())))),
                differential=dict(residual_km_max=float(np.abs(residual).max()/1e3),
                                  overlap_change_per_step_m=dict(opening_max=float(change.max()),
                                                                 closing_max=float(-change.min())),
                                  opening_rate_m_per_day_max=float(rate.max()),
                                  closing_rate_m_per_day_max=float(-rate.min())),
                steering_over_reach_max=float(need.max()), turning_spread_deg_per_day=float(np.ptp(turns))), nfev


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--smoke', action='store_true', help='eight days and a few rings, to check the run')
    args = parser.parse_args(argv)
    wall, cpu0 = time.monotonic(), time.process_time()
    global RUN
    if args.smoke:
        fr.DESIGN['days'] = SMOKE['days']
        DESIGN.update({k: v for k, v in SMOKE.items() if k != 'days'})
        RUN = RUN.parent/'ring_layout_smoke'
    env = environment(fr.DESIGN['days'])
    sun, poles = geometry(env)
    frozen = json.loads(FROZEN_OUT.read_text())['passes'][-1]
    # 1. Calibration.
    cases = [dict(radius_km=radius, tilt_deg=tilt, design_height_m=radius*1e3*np.sin(np.radians(tilt)),
                  width_m=0.8*7.04e6) for radius in DESIGN['calibration_radii_km'] for tilt in DESIGN['calibration_tilts_deg']]
    calibrated, nfev = flights(env, sun, poles, cases, 'calibration', frozen)
    print(f'calibration done, {time.monotonic()-wall:.0f} s', flush=True)
    nr, nt = len(DESIGN['calibration_radii_km']), len(DESIGN['calibration_tilts_deg'])
    rate = np.array([t['rate_deg_per_day'] for t in calibrated['turns']]).reshape(nr, nt)
    ecc = np.array([np.hypot(*row['centre']) for row in calibrated['rows']]).reshape(nr, nt)
    # The turning of a ring near the Moon's orbit plane is not defined; take it from its tilted neighbours.
    flat = np.abs(np.array(DESIGN['calibration_tilts_deg'])) < DESIGN['untilted_deg']
    if flat.any() and (~flat).sum() >= 2:
        for i in range(nr):
            rate[i, flat] = np.interp(np.array(DESIGN['calibration_tilts_deg'])[flat],
                                      np.array(DESIGN['calibration_tilts_deg'])[~flat], rate[i, ~flat])
    cal = Calibration(DESIGN['calibration_radii_km'], DESIGN['calibration_tilts_deg'], rate, ecc)
    middle = DESIGN['middle_radius_km']
    heights, aperture = stack_heights(middle)
    # Every plane has to turn with the Sun: a turning common to the stack would still fan the strips out, top against
    # bottom, by the sine of each tilt, year after year.
    target = 0.
    untilted = float(cal.at(middle, 0.)[0][0])
    ideal = ideal_radii(cal, heights)
    cap = float(cal.radii[-1]+DESIGN['runaway_margin_km'])
    layouts, chosen = {}, None
    for clearance in DESIGN['clearance_m']:
        c = clearance/1e3
        row = dict(one_radius=summary(cal, np.full(len(heights), middle), heights, target, c, middle))
        for mode in ('forced', 'held'):
            held = mode == 'held'
            up = best_layout(cal, heights, lambda cl, h, s, cc, cp: staircase(cl, h, s, cc, True, mode, cp), c, target,
                             middle-2000., middle+600., cap)
            down = best_layout(cal, heights, lambda cl, h, s, cc, cp: staircase(cl, h, s, cc, False, mode, cp), c,
                               target, middle-2000., middle+600., cap)
            stair = min([b for b in (up, down) if b is not None], key=lambda b: b[0], default=None)
            sided = best_layout(cal, heights, lambda cl, h, s, cc, cp: two_sided(cl, h, s, cc, mode, cp), c, target,
                                middle-2600., middle+600., cap)
            match = best_matched(cal, heights, ideal, c, target, mode, cap)
            runaway = dict(runs_past_km=cap, note='every start runs the stack past the calibrated radii and the margin')
            row[mode] = dict(
                staircase=(dict(summary(cal, stair[2], heights, target, c, middle, held), rising_upward=bool(stair is up))
                           if stair else runaway),
                two_sided=summary(cal, sided[2], heights, target, c, middle, held) if sided else runaway,
                matched=(dict(summary(cal, match[3], heights, target, c, middle, held), offset_km=match[1],
                              stretch=match[2]) if match else runaway))
            if held and match:
                tilts, error, forced, need = evaluate(cal, match[3], heights, target)
                kept = held_profile(match[3], heights, forced, c)
                change = held_turning_change(forced, kept)
                need_kept = np.abs(error+change)*np.abs(np.sin(np.radians(tilts)))/reach(match[3])
                row[mode]['matched']['held_eccentricity'] = dict(
                    largest_departure=float(np.abs(kept-forced).max()),
                    departure_at_edges=float(np.abs(kept-forced)[np.argsort(np.abs(heights))[-20:]].mean()),
                    turning_change_deg_per_day=dict(min=float(change.min()), max=float(change.max())),
                    steering_over_reach_with_held_eccentricity=dict(max=float(need_kept.max()),
                                                                    rings_over=int(np.sum(need_kept > 1))))
                if clearance == DESIGN['clearance_design_m']:
                    chosen = match[3]
                radii, kept, change, need_coupled, passes = coupled_matched(cal, heights, c, target)
                tilts, error, forced, _ = evaluate(cal, radii, heights, target)
                row[mode]['matched_coupled'] = dict(
                    summary(cal, radii, heights, target, c, middle, held),
                    held_eccentricity=dict(largest_departure=float(np.abs(kept-forced).max()),
                                           turning_change_deg_per_day=dict(min=float(change.min()),
                                                                           max=float(change.max()))),
                    steering_over_reach_with_held_eccentricity=dict(max=float(need_coupled.max()),
                                                                    rings_over=int(np.sum(need_coupled > 1)),
                                                                    over_a_tenth=int(np.sum(need_coupled > .1))),
                    passes=passes)
        layouts[f'{clearance:g}'] = row
        print(f'clearance {clearance:g} m:', json.dumps({m: {k: v.get('steering_over_reach', 'runs away')
                                                             for k, v in row[m].items()} for m in ('forced', 'held')}),
              flush=True)
    print(f'layouts done, {time.monotonic()-wall:.0f} s', flush=True)
    # 3. A year's flight of a sample of the matched layout's rings, on their forced eccentricities.
    flown = {}
    if chosen is not None:
        rings = sample(chosen, heights, DESIGN['sample_rings'])
        flown['matched'], used = measure_layout(env, sun, poles, rings, 'matched', frozen, cal, target)
        nfev += used
    one = json.loads(PLANES_OUT.read_text())
    idealised = [c for c in one['cases'] if c['plane'] == DESIGN['plane'] and c['radius_km'] == middle]
    out = dict(schema='terluna.research.ring-layout/1',
               producer=dict(files={f: digest(ROOT/f) for f in FILES},
                             constants=constants_used([f for f in FILES if f.endswith('.py')]),
                             inputs={'results/frozen_rings.json': digest(FROZEN_OUT),
                                     'results/plane_motion.json': digest(PLANES_OUT),
                                     'results/attitude_schemes.json': digest(ATTITUDE), 'de440s.bsp': digest(DEFAULT_KERNEL)}),
               evidence=('Single tiles on retrograde rings about the Moon, one year from 2026-10-04 TDB in DE440s '
                         'point-mass gravity with frozen_rings\' sail force at 26 g/m2, each started on its forced '
                         'eccentricity after one refinement pass. The turning is the yearly trend of the azimuth of each '
                         'ring\'s pole about the Moon\'s osculating orbit pole, less the Sun\'s. The layouts take the '
                         'calibration\'s turning and forced eccentricity interpolated in radius and tilt; the sampled '
                         'rings of the held matched layout are flown on their forced eccentricities and measured as '
                         'plane_motion measures its rings.'),
               reading_rule=('A ring turns with the Sun where its turning is zero; the steering it needs is its '
                             'strip\'s turn, the turning times the sine of its tilt, over the photon roll\'s reach at '
                             '26 g/m2. Rings over 1 need more than photon roll gives. layouts[clearance][mode]: forced, '
                             'every ring on its forced eccentricity; held, keeping holds the eccentricity vectors of '
                             'rings that come close together; a layout that runs_past_km found no start keeping the '
                             'stack within the calibrated radii and the margin. The held matched layout\'s '
                             'held_eccentricity estimates how far keeping must move each ring off its forced '
                             'eccentricity and what that does to its turning (secular tidal theory). Clearances: neighbouring strips at apolune and perilune, radius '
                             'neighbours at the nodes. The one-radius layout is plane_motion\'s idealisation, and its '
                             'zero clearance shows that no stack can fly it.'),
               design=DESIGN, smoke=bool(args.smoke), target_turning_deg_per_day=target,
               untilted_middle_turning_deg_per_day=untilted,
               stack=dict(rings=int(len(heights)), aperture_km=float(aperture), height_top_km=float(heights.max()),
                          ideal_radius_km=dict(min=float(ideal.min()), max=float(ideal.max()),
                                               by_height=[dict(height_km=round(float(heights[k]), 1),
                                                               radius_km=round(float(ideal[k]), 1))
                                                          for k in np.linspace(0, len(heights)-1, 17).astype(int)])),
               calibration=dict(radii_km=DESIGN['calibration_radii_km'], tilts_deg=DESIGN['calibration_tilts_deg'],
                                turning_deg_per_day=rate.round(6).tolist(), forced_eccentricity=ecc.round(5).tolist(),
                                swing_deg=[round(t['swing_deg'], 3) for t in calibrated['turns']]),
               layouts=layouts, flown=flown,
               idealised_plane_motion=dict(steering_over_reach=idealised[0]['photon_steering']['needed_over_available_by_ring']
                                           if idealised else None),
               resources=dict(cpu_s=time.process_time()-cpu0, wall_s=time.monotonic()-wall, nfev=nfev,
                              max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024, numerical_threads=1))
    target_path = RUN/'ring_layout_smoke.json' if args.smoke else OUT
    target_path.write_text(json.dumps(out, indent=1, default=float)+'\n')


if __name__ == '__main__':
    main()
