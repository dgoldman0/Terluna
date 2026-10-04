"""Independent propagation and omitted finite-extent force diagnostics."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.fleet import fleet_acceleration, square_vertices
from .fleet_run import read_raw, environment, identities, digest
from .cycling_search import HERE


def finite_extent_gravity(env, t, state, side):
    d = fleet_acceleration(env, t, state, 4*K.MOON_RADIUS, side, diagnostics=True)
    vertices = square_vertices(state[:, :3], d['normal'], d['sun'], side)
    a = (vertices[:, 1]-vertices[:, 0])/side
    b = (vertices[:, 3]-vertices[:, 0])/side
    x, w = np.polynomial.legendre.leggauss(3)
    total = np.zeros_like(state[:, :3])
    for i in range(3):
        for j in range(3):
            q = state[:, :3]+side/2*(x[i]*a+x[j]*b)
            total += env.gravity(t, q)*(w[i]*w[j]/4)
    error = np.linalg.norm(total-env.gravity(t, state[:, :3]), axis=-1)
    return dict(maximum_m_s2=float(error.max()), mean_m_s2=float(error.mean()),
                quadrature='3 by 3 Gauss-Legendre over a uniform rigid square')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--names', nargs='+', required=True)
    parser.add_argument('--replay', help='Replay eight widely spaced members from this named run')
    parser.add_argument('--output', default='fleet_validation.json')
    a = parser.parse_args(); cases = {}
    for name in a.names:
        raw, meta = read_raw(name); t = raw['t']; states = raw['state']
        side = meta['config']['side_km']*1000; env = environment(meta['config']['days'])
        force = [finite_extent_gravity(env, t[i], states[i], side) for i in (0, len(t)//2, len(t)-1)]
        rate = np.linalg.norm(np.cross(env.positions(t)[:, env.names.index('sun')],
            env.positions(t, 1)[:, env.names.index('sun')]), axis=-1)/np.linalg.norm(
            env.positions(t)[:, env.names.index('sun')], axis=-1)**2
        result = dict(raw_sha256=meta['sha256'], finite_extent_gravity=force,
            max_sun_direction_rate_rad_s=float(rate.max()), assumed_shadow_guard_rate_rad_s=3e-7,
            finite_extent_forces_in_main_propagation=False,
            reading_rule='These residual forces quantify one omitted term; they do not supply a corrected coupled fleet propagation.')
        if name == a.replay:
            ids = np.linspace(0, len(states[0])-1, 8, dtype=int)
            initial = states[0, ids]
            def f(time, flat):
                state = flat.reshape(len(ids), 6)
                acc = fleet_acceleration(env, time, state, 4*K.MOON_RADIUS, side)
                return np.column_stack([state[:, 3:], acc]).ravel()
            sol = solve_ivp(f, [t[0], t[-1]], initial.ravel(), t_eval=t, method='DOP853',
                rtol=2e-12, atol=np.tile([1e-5]*3+[1e-9]*3, len(ids)), max_step=300.)
            if not sol.success:
                raise RuntimeError(sol.message)
            reference = sol.y.T.reshape(len(t), len(ids), 6)
            error = np.linalg.norm(reference[..., :3]-states[:, ids, :3], axis=-1)
            result['independent_replay'] = dict(member_ids=ids.tolist(), days=float(t[-1]/K.JULIAN_DAY),
                max_position_difference_m=float(error.max()),
                final_position_difference_m=float(error[-1].max()),
                max_step_s=300., relative_tolerance=2e-12,
                trajectory_agreement_within_50m=bool(error.max() < 50.),
                comparison='Whole month, eight distinct phases, fourfold smaller maximum step and 100-fold tighter relative tolerance')
        cases[name] = result; raw.close()
        print(json.dumps(dict(name=name, **result)), flush=True)
    sources = {**identities(), 'research/studies/solar_shield_array/fleet_validate.py': digest(__file__)}
    product = dict(schema='terluna.research.solar-shield-fleet-validation/1',
        producer=dict(source_hashes=sources, constants=constants_used(sources)),
        evidence='Numerical centre-model replay and finite-extent gravity residual; operational tracking and material forces remain unvalidated',
        cases=cases)
    (HERE/'results'/a.output).write_text(json.dumps(product, indent=2)+'\n')


if __name__ == '__main__':
    main()
