"""Earth's tide on molecular escape: test molecules in the Earth-Moon restricted three-body problem.

    OPENBLAS_NUM_THREADS=1 python -m atmosphere.loss_response.tides   # writes results/tidal_escape.json

Jeans escape counts the molecules that leave the exobase faster than the Moon's
two-body escape speed. Earth's tide lowers the barrier. In the frame turning
with the Moon's orbit, a molecule whose Jacobi constant is below its value at
the L1 point can leave the Moon's Hill sphere through the neck around L1 or L2
while slower than the two-body escape speed. Whether it does depends on where it
starts and where it is aimed. A uniform lowering of the barrier, the feasibility
baseline's factor K = 1 - 3/(2 xi) + 1/(2 xi^3) with xi the Hill radius over the
exobase radius, counts every such molecule as lost and gives an upper figure.

This step follows test molecules in the circular restricted three-body problem
(mu = GM_moon / (GM_earth + GM_moon), the Moon's mean distance, units of that
distance and of the orbit's angular rate). They start from random points on an
exobase sphere of radius r_c around the Moon, in directions weighted by the
cosine from the vertical (the outgoing flux of a Maxwellian), with speeds
sampled uniformly in v^2 between the lowest that can pass L1 from that point and
the two-body escape speed. The atmosphere turns with the Moon, so speeds are
taken in the rotating frame. Each molecule is followed with fourth-order
Runge-Kutta steps of a fixed fraction of its local orbital time until it falls
back inside the exobase, passes two Hill radii (escaped: it then orbits Earth
and is ionized within weeks), or 120 days pass (held).

Weighting the outcomes by the flux-weighted Maxwellian, p(eps) = eps exp(-eps)
with eps = v^2 / v_th^2 = lambda v^2 / v_esc^2, gives the escape flux relative
to Jeans:

    M(r_c, lambda) = 1 + <w P_esc> / ((1 + lambda) exp(-lambda)),
    w = eps exp(-eps) (lambda - eps_min),

where lambda is the two-body Jeans parameter at the exobase, eps_min the local
threshold, and molecules above the two-body escape speed escape as in Jeans (a
sample of them checks this). Held molecules count as returning in M and as lost
in M_held_lost. The model leaves out the Sun's tide (under 1% of Earth's at
the Hill radius), the eccentricity and inclination of the Moon's orbit,
collisions above the exobase and radiation pressure.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq

from shared.constants import EARTH_GM, EARTH_MOON_DISTANCE, MOON_GM, MOON_RADIUS

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCHEMA = 'terluna.atmosphere.tidal-escape/1'
MU = MOON_GM / (EARTH_GM + MOON_GM)
DAY_UNITS = 86400.0 * math.sqrt((EARTH_GM + MOON_GM) / EARTH_MOON_DISTANCE ** 3)   # one day in units of 1/Omega
HILL = (MU / 3.0) ** (1.0 / 3.0)
RADII_R = (1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0)
LAMBDAS = tuple(float(x) for x in np.arange(3.0, 30.5, 0.5))
SAMPLES = 16000
ABOVE_ESCAPE_SAMPLES = 1000
ETA = 0.01
T_MAX_DAYS = 120.0
ESCAPE_RADIUS_HILL = 2.0
MAX_STEPS = 400000


def moon_distance(x, y, z):
    return np.sqrt((x - 1.0 + MU) ** 2 + y ** 2 + z ** 2)


def omega(x, y, z):
    r1 = np.sqrt((x + MU) ** 2 + y ** 2 + z ** 2)
    return 0.5 * (x ** 2 + y ** 2) + (1.0 - MU) / r1 + MU / moon_distance(x, y, z)


def jacobi(s):
    return 2.0 * omega(s[0], s[1], s[2]) - (s[3] ** 2 + s[4] ** 2 + s[5] ** 2)


def l1_jacobi():
    """Jacobi constant at L1, between Earth (x = -mu) and Moon (x = 1 - mu)."""
    grad = lambda x: x - (1 - MU) / (x + MU) ** 2 + MU / (x - 1 + MU) ** 2
    x1 = brentq(grad, 0.5, 1.0 - MU - 1e-6)
    return float(2.0 * omega(x1, 0.0, 0.0)), float(x1)


def derivatives(s):
    x, y, z, vx, vy, vz = s
    r1 = np.sqrt((x + MU) ** 2 + y ** 2 + z ** 2)
    r2 = moon_distance(x, y, z)
    a1, a2 = (1.0 - MU) / r1 ** 3, MU / r2 ** 3
    return np.array([vx, vy, vz,
                     2 * vy + x - a1 * (x + MU) - a2 * (x - 1 + MU),
                     -2 * vx + y - a1 * y - a2 * y,
                     -a1 * z - a2 * z])


def rk4(s, dt):
    k1 = derivatives(s)
    k2 = derivatives(s + 0.5 * dt * k1)
    k3 = derivatives(s + 0.5 * dt * k2)
    k4 = derivatives(s + dt * k3)
    return s + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4)


def launch(r_c_R, n, rng, above_escape=False):
    """Launch states on the exobase sphere, with v^2 uniform between the L1 threshold and two-body escape."""
    r_c = r_c_R * MOON_RADIUS / EARTH_MOON_DISTANCE
    up = rng.normal(size=(3, n))
    up /= np.linalg.norm(up, axis=0)
    pos = np.array([1.0 - MU, 0.0, 0.0])[:, None] + r_c * up
    helper = np.where(np.abs(up[2]) < 0.9, np.array([0, 0, 1.0])[:, None], np.array([1.0, 0, 0])[:, None])
    e1 = np.cross(up.T, helper.T).T
    e1 /= np.linalg.norm(e1, axis=0)
    e2 = np.cross(up.T, e1.T).T
    cos_t = np.sqrt(rng.random(n))
    phi = 2 * math.pi * rng.random(n)
    direction = cos_t * up + np.sqrt(1 - cos_t ** 2) * (np.cos(phi) * e1 + np.sin(phi) * e2)
    c_l1, _ = l1_jacobi()
    v_min2 = np.maximum(2.0 * omega(*pos) - c_l1, 0.0)
    v_esc2 = 2.0 * MU / r_c
    if above_escape:
        v2 = v_esc2 * (1.0 + 0.2 * rng.random(n))
    else:
        v2 = v_min2 + rng.random(n) * (v_esc2 - v_min2)
    state = np.vstack([pos, np.sqrt(v2) * direction])
    return state, v2, v_min2, v_esc2, r_c


def follow(state, r_c, t_max=T_MAX_DAYS * DAY_UNITS, eta=ETA, r_out=ESCAPE_RADIUS_HILL * HILL):
    """Outcome per molecule: 1 returned inside the exobase, 2 escaped past r_out, 3 held at t_max."""
    s = state.copy()
    n = s.shape[1]
    outcome = np.zeros(n, dtype=int)
    t = np.zeros(n)
    c0 = jacobi(s)
    active = np.arange(n)
    for _ in range(MAX_STEPS):
        if not active.size:
            break
        st = s[:, active]
        dt = eta * np.minimum(np.sqrt(moon_distance(*st[:3]) ** 3 / MU), 1.0)
        st = rk4(st, dt)
        s[:, active] = st
        t[active] += dt
        r2 = moon_distance(*st[:3])
        done = np.where(r2 <= r_c, 1, np.where(r2 >= r_out, 2, np.where(t[active] >= t_max, 3, 0)))
        outcome[active] = done
        active = active[done == 0]
    outcome[active] = 3
    drift = np.abs(jacobi(s) - c0)
    return outcome, t, drift


def multipliers(v2, v_min2, v_esc2, outcome, lambdas=LAMBDAS):
    """Tidal escape multiplier over lambda, counting held molecules as returning and as lost."""
    out, held_lost, bound = [], [], []
    for lam in lambdas:
        eps = lam * v2 / v_esc2
        eps_min = lam * v_min2 / v_esc2
        w = eps * np.exp(-eps) * (lam - eps_min)
        jeans = (1.0 + lam) * math.exp(-lam)
        out.append(1.0 + float(np.mean(w * (outcome == 2))) / jeans)
        held_lost.append(1.0 + float(np.mean(w * (outcome >= 2))) / jeans)
        bound.append(float(np.mean((1 + eps_min) * np.exp(-eps_min))) / jeans)
    return out, held_lost, bound


def erkaev(r_c_R, lam):
    xi = HILL * EARTH_MOON_DISTANCE / (r_c_R * MOON_RADIUS)
    k = 1.0 - 1.5 / xi + 0.5 / xi ** 3
    return (1.0 + k * lam) / (1.0 + lam) * math.exp(lam * (1.0 - k))


def run_radius(r_c_R, samples=SAMPLES, seed=0):
    rng = np.random.default_rng(seed + int(round(r_c_R * 100)))
    state, v2, v_min2, v_esc2, r_c = launch(r_c_R, samples, rng)
    outcome, t, drift = follow(state, r_c)
    m, m_held, bound = multipliers(v2, v_min2, v_esc2, outcome)
    above, _, _, _, _ = launch(r_c_R, ABOVE_ESCAPE_SAMPLES, rng, above_escape=True)
    above_outcome, _, _ = follow(above, r_c)
    counts = {name: int(np.sum(outcome == code)) for code, name in ((1, 'returned'), (2, 'escaped'), (3, 'held'))}
    return dict(exobase_radius_R=r_c_R, samples=samples, counts=counts,
                barrier_lowering_share=float(np.mean(1.0 - v_min2 / v_esc2)),
                escape_fraction_above_two_body_escape=float(np.mean(above_outcome == 2)),
                max_jacobi_drift_over_barrier=float(drift.max() / np.mean(v_esc2 - v_min2)),
                median_escape_days=float(np.median(t[outcome == 2]) / DAY_UNITS) if counts['escaped'] else None,
                multiplier=m, multiplier_held_lost=m_held, energy_bound=bound,
                erkaev=[erkaev(r_c_R, lam) for lam in LAMBDAS])


def multiplier(table, r_c_R, lam, key='multiplier'):
    """Interpolate the tidal multiplier: linear in exobase radius, linear in log M over lambda."""
    radii = np.array([row['exobase_radius_R'] for row in table['radii']])
    lams = np.array(table['lambdas'])
    logm = np.log(np.array([row[key] for row in table['radii']]))
    r = float(np.clip(r_c_R, radii[0], radii[-1]))
    lam = float(np.clip(lam, lams[0], lams[-1]))
    by_radius = np.array([np.interp(lam, lams, row) for row in logm])
    return float(np.exp(np.interp(r, radii, by_radius)))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out', type=Path, default=HERE / 'results')
    parser.add_argument('--samples', type=int, default=SAMPLES)
    args = parser.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    partial = args.out / 'tidal_escape.partial.json'
    done = json.loads(partial.read_text()) if partial.exists() else {}
    for r in RADII_R:
        key = f'{r:g}'
        if key in done and done[key]['samples'] == args.samples:
            continue
        row = run_radius(r, args.samples)
        done[key] = row
        partial.write_text(json.dumps(done))
        i11 = LAMBDAS.index(11.0)
        print(f"r_c {r:3.1f} R: barrier lowered {row['barrier_lowering_share']:.3f}; outcomes {row['counts']}; "
              f"above-escape escape {row['escape_fraction_above_two_body_escape']:.3f}; M(lambda 11) "
              f"{row['multiplier'][i11]:.2f} (held lost {row['multiplier_held_lost'][i11]:.2f}, energy bound "
              f"{row['energy_bound'][i11]:.2f}, K {row['erkaev'][i11]:.2f}); Jacobi drift/barrier "
              f"{row['max_jacobi_drift_over_barrier']:.1e}", flush=True)
    c_l1, x_l1 = l1_jacobi()
    product = dict(
        schema=SCHEMA,
        producer=dict(domain='atmosphere', files={'atmosphere/loss_response/tides.py': hashlib.sha256(
            (ROOT / 'atmosphere/loss_response/tides.py').read_bytes()).hexdigest()[:16]}),
        evidence=('Test molecules from an exobase sphere in the circular restricted Earth-Moon three-body problem, '
                  'weighted by the outgoing Maxwellian flux. Collisionless above the exobase; no Sun, no orbital '
                  'eccentricity or inclination. A molecule counts as escaped past two Hill radii.'),
        reading_rule=('radii[i].multiplier[j] multiplies the two-body Jeans escape at exobase radius radii[i] '
                      '(lunar radii) and Jeans parameter lambdas[j]; interpolate linearly in radius and in log M over '
                      'lambda (tides.multiplier). multiplier_held_lost counts molecules still inside two Hill radii '
                      'after 120 days as lost; energy_bound counts every molecule above the L1 energy as lost; erkaev '
                      'is the feasibility baseline\'s uniform factor.'),
        mu=MU, hill_radius_R=HILL * EARTH_MOON_DISTANCE / MOON_RADIUS, l1_jacobi=c_l1,
        l1_distance_from_moon_km=(1.0 - MU - x_l1) * EARTH_MOON_DISTANCE / 1e3,
        lambdas=list(LAMBDAS), radii=[done[f'{r:g}'] for r in RADII_R])
    (args.out / 'tidal_escape.json').write_text(json.dumps(product, indent=1) + '\n')
    partial.unlink()
    return 0


if __name__ == '__main__':
    sys.exit(main())
