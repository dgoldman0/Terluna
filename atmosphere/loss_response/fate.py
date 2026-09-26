"""Where the escaping air goes: the cloud it forms around the Earth-Moon system, and its end.

    OPENBLAS_NUM_THREADS=1 python -m atmosphere.loss_response.fate   # writes results/escape_fate.json

Molecules leave the Moon's Hill sphere through the necks around L1 and L2, or
faster than two-body escape (tides.py). Beyond the Hill sphere they orbit Earth
near the Moon's orbit until sunlight destroys them: N2 is ionized in 8-19 days
at 1 AU by the repository's own rate (absorption.py) and dissociated as well,
and O2 is dissociated within days (Huebner and Mukherjee 2015 tabulate the
rates). The fragments of dissociation carry a few km/s, more than
Earth's escape speed at the Moon's distance (1.4 km/s), so they leave the
Earth-Moon system; the ions are swept off by the solar wind, or, while in Earth's
magnetotail, carried by its plasma.

This step launches molecules from the exobase with the outgoing Maxwellian flux
above the energy that can pass L1, in the circular restricted Earth-Moon
problem of tides.py, and follows them for 120 days, until they fall back on the
Moon, reach Earth's exobase (500 km up), or pass 1.5 million km from Earth
(Earth's Hill radius in the Sun's field, where the Sun's tide takes over).
Positions are recorded every quarter day once a molecule has left the Moon's
Hill sphere. Weighting each record by the chance the molecule is still intact,
exp(-t / tau), gives the steady cloud of escaped air: its mass per kg/s of
escape, where it lies, and how it ends. The lifetime tau is taken at 3, 10 and
20 days, bracketing O2 and N2 at quiet and active Sun.

Earth's magnetotail is taken as a cylinder of 25 Earth radii around the
anti-sunward line behind Earth; the Sun turns once a synodic month in the
rotating frame, from a random phase for each molecule. The share of the cloud's
destruction that happens inside the tail bounds the ions that Earth's
magnetosphere could keep for a while. Photon pressure, the Sun's tide inside
the region, the Moon's eccentric orbit and collisions are left out.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

from shared.constants import AVOGADRO, EARTH_MOON_DISTANCE, EARTH_RADIUS, SYNODIC_MONTH_DAYS
from atmosphere.loss_response import tides as td

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCHEMA = 'terluna.atmosphere.escape-fate/1'
EXOBASE_R = 3.0
LAMBDA = 11.0
SAMPLES = 6000
T_MAX_DAYS = 120.0
RECORD_DAYS = 0.25
LIFETIMES_DAYS = (3.0, 10.0, 20.0)
EARTH_CAPTURE = (EARTH_RADIUS + 500e3) / EARTH_MOON_DISTANCE
REGION = 1.5e9 / EARTH_MOON_DISTANCE
TAIL_RADIUS = 25.0 * EARTH_RADIUS / EARTH_MOON_DISTANCE
AIR_KG = 0.0289 / AVOGADRO
MAX_STEPS = 400000
SUN_RATE = 2 * math.pi / SYNODIC_MONTH_DAYS / td.DAY_UNITS    # the Sun turns backwards once a synodic month


def earth_distance(x, y, z):
    return np.sqrt((x + td.MU) ** 2 + y ** 2 + z ** 2)


def launch_escaping(r_c_R, lam, n, rng):
    """Launch states with the outgoing Maxwellian flux above the local L1 threshold: p(eps) ~ eps exp(-eps)."""
    state, _, v_min2, v_esc2, r_c = td.launch(r_c_R, n, rng)
    direction = state[3:] / np.linalg.norm(state[3:], axis=0)
    eps_min = lam * v_min2 / v_esc2
    exponential = rng.random(n) < eps_min / (eps_min + 1.0)
    eps = eps_min + np.where(exponential, rng.exponential(1.0, n), rng.gamma(2.0, 1.0, n))
    state[3:] = np.sqrt(eps * v_esc2 / lam) * direction
    return state, eps, eps_min, r_c


def follow(state, r_c, t_max=T_MAX_DAYS * td.DAY_UNITS, record=RECORD_DAYS * td.DAY_UNITS, eta=td.ETA):
    """Outcomes (1 Moon, 2 Earth, 3 left the region, 4 still in it), end times, and records beyond the Hill sphere."""
    s = state.copy()
    n = s.shape[1]
    outcome = np.zeros(n, dtype=int)
    t = np.zeros(n)
    left_hill = np.zeros(n, dtype=bool)
    next_record = np.full(n, record)
    records = []
    active = np.arange(n)
    for _ in range(MAX_STEPS):
        if not active.size:
            break
        st = s[:, active]
        r2, r1 = td.moon_distance(*st[:3]), earth_distance(*st[:3])
        dt = eta * np.minimum(np.minimum(np.sqrt(r2 ** 3 / td.MU), np.sqrt(r1 ** 3 / (1 - td.MU))), 1.0)
        st = td.rk4(st, dt)
        s[:, active] = st
        t[active] += dt
        r2, r1 = td.moon_distance(*st[:3]), earth_distance(*st[:3])
        left_hill[active] |= r2 > td.HILL
        due = (t[active] >= next_record[active]) & (r2 > td.HILL)
        if due.any():
            idx = active[due]
            records.append((idx, t[idx].copy(), st[:3, due].copy()))
        next_record[active] = np.where(t[active] >= next_record[active], next_record[active] + record, next_record[active])
        done = np.where(r2 <= r_c, 1, np.where(r1 <= EARTH_CAPTURE, 2, np.where(r1 >= REGION, 3,
                        np.where(t[active] >= t_max, 4, 0))))
        outcome[active] = done
        active = active[done == 0]
    outcome[active] = 4
    return outcome, t, left_hill, records


def in_tail(positions, times, phase):
    """Whether each recorded position lies in Earth's magnetotail at its time."""
    angle = phase - SUN_RATE * times                     # direction of the Sun seen from Earth, rotating frame
    anti = -np.array([np.cos(angle), np.sin(angle)])
    geo = positions[:2] - np.array([-td.MU, 0.0])[:, None]
    along = np.sum(geo * anti, axis=0)
    across2 = np.sum(geo ** 2, axis=0) - along ** 2 + positions[2] ** 2
    return (along > 0) & (across2 < TAIL_RADIUS ** 2)


def cloud(outcome, t, left_hill, records, phases, tau_days):
    """Mass per unit escape rate, where it lies and how it ends, for an intact lifetime tau."""
    tau = tau_days * td.DAY_UNITS
    escaped = left_hill
    n_esc = int(escaped.sum())
    idx = np.concatenate([r[0] for r in records])
    times = np.concatenate([r[1] for r in records])
    pos = np.concatenate([r[2] for r in records], axis=1)
    weight = np.exp(-times / tau)
    dt_record = RECORD_DAYS * 86400.0
    seconds_per_escape = float(weight.sum() * dt_record / n_esc)        # kg in the cloud per kg/s of escape
    geo = pos[:2] - np.array([-td.MU, 0.0])[:, None]
    r_geo = np.hypot(geo[0], geo[1]) * EARTH_MOON_DISTANCE / 1e3
    phase = np.degrees(np.arctan2(geo[1], geo[0]))       # angle ahead of the Moon, which sits at 0
    height = np.abs(pos[2]) * EARTH_MOON_DISTANCE / 1e3

    def quantiles(values, qs):
        order = np.argsort(values)
        cum = np.cumsum(weight[order]) / weight.sum()
        return [float(values[order][np.searchsorted(cum, q)]) for q in qs]

    # 10-30 degrees ahead of and behind the Moon, clear of its Hill sphere (9 degrees), within 50,000 km of the orbit
    near = ((np.abs(phase) > 10) & (np.abs(phase) < 30) & (np.abs(r_geo - EARTH_MOON_DISTANCE / 1e3) < 50e3)
            & (height < 50e3))
    box_m3 = 2 * (math.radians(20) * EARTH_MOON_DISTANCE) * 100e6 * 100e6
    density = float(weight[near].sum() * dt_record / n_esc / AIR_KG / box_m3 / 1e6)   # per cm^3 per kg/s
    tail = in_tail(pos, times, phases[idx])
    end = np.exp(-t / tau)
    fates = {name: float(np.sum(end[escaped & (outcome == code)]) / n_esc)
             for code, name in ((1, 'back_on_the_moon'), (2, 'into_earths_atmosphere'), (3, 'left_the_earth_moon_region'))}
    fates['destroyed_by_sunlight'] = 1.0 - sum(fates.values())
    return dict(lifetime_days=tau_days, cloud_kg_per_kg_s=seconds_per_escape,
                molecules_per_cm3_near_the_moons_orbit_per_kg_s=density,
                ahead_of_the_moon_share=float(weight[phase > 0].sum() / weight.sum()),
                angle_from_the_moon_deg_5_50_95=quantiles(phase, (0.05, 0.5, 0.95)),
                distance_from_earth_km_5_50_95=quantiles(r_geo, (0.05, 0.5, 0.95)),
                height_above_orbit_plane_km_50_95=quantiles(height, (0.5, 0.95)),
                fates=fates,
                destroyed_inside_magnetotail_share=float(np.sum((weight / tau)[tail]) / np.sum(weight / tau)))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out', type=Path, default=HERE / 'results')
    parser.add_argument('--samples', type=int, default=SAMPLES)
    args = parser.parse_args(argv)
    rng = np.random.default_rng(7)
    state, eps, eps_min, r_c = launch_escaping(EXOBASE_R, LAMBDA, args.samples, rng)
    outcome, t, left_hill, records = follow(state, r_c)
    phases = rng.random(args.samples) * 2 * math.pi
    clouds = [cloud(outcome, t, left_hill, records, phases, tau) for tau in LIFETIMES_DAYS]
    escaped_share = float(left_hill.mean())
    product = dict(
        schema=SCHEMA,
        producer=dict(domain='atmosphere', files={p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()[:16] for p in (
            'atmosphere/loss_response/fate.py', 'atmosphere/loss_response/tides.py')}),
        evidence=('Test molecules in the circular restricted Earth-Moon problem, launched with the outgoing Maxwellian '
                  'flux above the L1 energy from an exobase of 3 lunar radii at a Jeans parameter of 11, followed for '
                  '120 days; the cloud weights them by an assumed intact lifetime. No Sun, photon pressure, plasma, '
                  'collisions or orbital eccentricity.'),
        reading_rule=('clouds[i] is the steady cloud of escaped air for intact lifetime lifetimes_days[i]: its mass per '
                      'kg/s of escape, a number density near the Moon\'s orbit per kg/s of escape, where it lies '
                      '(angles ahead of the Moon, distances from Earth, heights), and the share of escaped molecules '
                      'that end each way. Shares count molecules that left the Moon\'s Hill sphere.'),
        exobase_radius_R=EXOBASE_R, jeans_lambda=LAMBDA, samples=args.samples,
        share_leaving_the_hill_sphere=escaped_share,
        outcomes_unweighted={name: int(np.sum(left_hill & (outcome == code))) for code, name in (
            (1, 'moon'), (2, 'earth'), (3, 'left_region'), (4, 'still_in_region'))},
        clouds=clouds)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / 'escape_fate.json').write_text(json.dumps(product, indent=1) + '\n')
    print(f"{args.samples} launched above the L1 energy; {escaped_share:.2f} left the Hill sphere; outcomes of those "
          f"{product['outcomes_unweighted']}")
    for c in clouds:
        f = c['fates']
        print(f"tau {c['lifetime_days']:4.0f} d: cloud {c['cloud_kg_per_kg_s']/1e3:7.1f} t per kg/s, "
              f"{c['molecules_per_cm3_near_the_moons_orbit_per_kg_s']:.2f} /cm3 near the orbit per kg/s; ahead "
              f"{c['ahead_of_the_moon_share']:.2f}; angle 5/50/95% {np.round(c['angle_from_the_moon_deg_5_50_95'], 0)}; "
              f"distance {np.round(c['distance_from_earth_km_5_50_95'], -3)} km; height 50/95% "
              f"{np.round(c['height_above_orbit_plane_km_50_95'], -3)} km; Moon {f['back_on_the_moon']:.3f}, Earth "
              f"{f['into_earths_atmosphere']:.4f}, left {f['left_the_earth_moon_region']:.3f}, destroyed "
              f"{f['destroyed_by_sunlight']:.3f}; in tail {c['destroyed_inside_magnetotail_share']:.2f}")
    return 0


if __name__ == '__main__':
    sys.exit(main())
