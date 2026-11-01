"""Isolate mass-grid spreading in vapour growth with continuous characteristics.

The frozen sectional column is not changed. A fixed-air box has no collisions,
settling or fresh particle injection after t=0. Its deposition and shared vapour
cap are the existing laws. Continuous masses with constant number weights give
an independent characteristic reference. A copy of only the sectional growth
operator uses the frozen SSPRK2 stepper. Matching each grid's initial projection
and projecting its characteristics just once at the end separates repeated
mass-grid transport from initial/final pivot representation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
import scipy
from scipy.integrate import solve_ivp

from atmosphere.snow import aggregation as ag, flakes as sf, sectional
from atmosphere.snow.habits import GEOMETRIES
from shared.provenance import constants_used

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'atmosphere/snow/results/growth_transport.json'
SCHEMA = 'terluna.atmosphere.snow-growth-transport/1'
FILES = ('atmosphere/snow/growth_transport.py', 'atmosphere/snow/sectional.py',
         'atmosphere/snow/aggregation.py', 'atmosphere/snow/habits.py',
         'atmosphere/snow/flakes.py', 'atmosphere/rain/drops.py',
         'atmosphere/rain/targets.py', 'atmosphere/radiative_convective/thermodynamics.py')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _deadline(deadline):
    if deadline is not None and time.monotonic() > deadline:
        raise TimeoutError('growth comparison exceeded its stated wall-time budget')


def _validate(masses, numbers, times, production):
    m, n, times = [np.asarray(a, dtype=float) for a in (masses, numbers, times)]
    if m.ndim != 1 or n.shape != m.shape or not np.isfinite(m).all() or (m <= 0).any():
        raise ValueError('finite positive masses and matching number weights required')
    if not np.isfinite(n).all() or (n < 0).any() or n.sum() <= 0:
        raise ValueError('finite nonnegative weights with positive total required')
    if times.ndim != 1 or times.size < 2 or times[0] != 0 or not np.isfinite(times).all() or (np.diff(times) <= 0).any():
        raise ValueError('times must increase strictly from zero')
    if not math.isfinite(production) or production < 0:
        raise ValueError('production must be finite and nonnegative')
    return m, n, times


def characteristics(masses, numbers, times, growth, production, *, rtol=1e-9, atol=1e-22, deadline=None):
    """Continuous-mass trajectories; zero-weight endpoint tracers are allowed.

    dm_a/dt = G(m_a) min(1, P / sum n_a G(m_a)). The final state also integrates
    the deposited-water ledger using the same derivative, without renormalizing
    masses. Characteristic weights remain exactly the initial number weights.
    """
    masses, numbers, times = _validate(masses, numbers, times, production)
    def rhs(t, state):
        _deadline(deadline)
        g = np.asarray(growth(state[:-1]), dtype=float)
        if not np.isfinite(g).all() or (g < 0).any():
            raise ValueError('growth law must return finite nonnegative rates')
        demand = float(np.dot(numbers, g))
        rate = g * (min(1., production / demand) if demand > 0 else 0.)
        return np.r_[rate, np.dot(numbers, rate)]
    start = time.monotonic()
    initial = np.r_[masses, 0.]
    solution = solve_ivp(rhs, (0., float(times[-1])), initial, method='DOP853', t_eval=times,
                         rtol=rtol, atol=atol, max_step=float(times[-1]) / 32.)
    if not solution.success:
        raise RuntimeError('characteristic integration failed: ' + solution.message)
    result_masses = solution.y[:-1].T
    if not np.isfinite(result_masses).all() or (result_masses <= 0).any():
        raise FloatingPointError('characteristic integration produced invalid masses')
    return dict(masses=result_masses, numbers=numbers, times=times,
                deposited_mass=solution.y[-1], evaluations=solution.nfev,
                wall_seconds=time.monotonic()-start, rtol=rtol, atol=atol)


def fixed_grid(masses, numbers, times, growth_rates, production, *, max_step_s=30., cfl=.3, deadline=None):
    """The frozen column's growth equations alone, with identical ghost exit.

    No tail is removed or renormalized. The main diagnostic has negligible
    boundary loss, but number and mass crossing the absorbing ghost are retained.
    """
    masses, numbers, times = _validate(masses, numbers, times, production)
    masses = sectional._grid(masses)
    growth_rates = np.asarray(growth_rates, dtype=float)
    if growth_rates.shape != masses.shape or not np.isfinite(growth_rates).all() or (growth_rates < 0).any():
        raise ValueError('growth rates must match pivots and be finite/nonnegative')
    if not math.isfinite(max_step_s) or max_step_s <= 0 or not 0 < cfl <= 1:
        raise ValueError('positive step and 0<cfl<=1 required')
    ghost = masses[-1]**2 / masses[-2]
    spacing = np.r_[np.diff(masses), ghost-masses[-1]]
    def evaluate(n):
        _deadline(deadline)
        demand = float(np.dot(n, growth_rates))
        rate = growth_rates * (min(1., production / demand) if demand > 0 else 0.)
        coefficient = rate / spacing
        transfer = coefficient * n
        dn = -transfer
        dn[1:] += transfer[:-1]
        occupied = n > 0
        fastest = float(coefficient[occupied].max()) if occupied.any() else 0.
        limit = cfl / fastest if fastest > 0 else math.inf
        diagnostic = np.array([np.dot(n, rate), transfer[-1]*ghost, transfer[-1]])
        return dn, diagnostic, limit
    start = time.monotonic()
    state = numbers.copy()
    snapshots, ledgers = [state.copy()], [np.zeros(3)]
    ledger = np.zeros(3)
    t, steps, rejected = 0., 0, 0
    smallest, largest = math.inf, 0.
    for end in times[1:]:
        while t < end:
            state, dt, rates, _, retries = sectional._advance(state, evaluate, min(max_step_s, end-t), cfl)
            ledger += rates*dt
            t += dt; steps += 1; rejected += retries
            smallest, largest = min(smallest, dt), max(largest, dt)
        snapshots.append(state.copy()); ledgers.append(ledger.copy())
    return dict(masses=masses, numbers=np.asarray(snapshots), times=times,
                deposited_mass=np.asarray(ledgers)[:, 0], upper_exit_mass=np.asarray(ledgers)[:, 1],
                upper_exit_number=np.asarray(ledgers)[:, 2], upper_ghost_mass=ghost,
                diagnostics=dict(steps=steps, rejected_steps=rejected, minimum_step_s=smallest,
                                 maximum_step_s=largest, wall_seconds=time.monotonic()-start),
                max_step_s=max_step_s, cfl=cfl)


def project_once(pivots, masses, numbers):
    """Conservative final projection, used only as a comparison, never repeated."""
    pivots = sectional._grid(pivots)
    result = np.zeros_like(pivots)
    for mass, number in zip(masses, numbers):
        if number == 0:
            continue
        allocation, side = sectional._allocation(pivots, mass)
        if side:
            raise ValueError('characteristic support leaves the comparison grid')
        for i, weight in allocation:
            result[i] += number*weight
    return result


def moments(masses, numbers, geometry):
    m, n = np.asarray(masses), np.asarray(numbers)
    d = geometry.size(m)
    number, mass = float(n.sum()), float(np.dot(n, m))
    mean = mass/number
    variance = float(np.dot(n, (m-mean)**2)/number)
    order = np.argsort(m)
    cumulative = np.cumsum((m*n)[order])/mass
    def quantile(q):
        return float(d[order[min(np.searchsorted(cumulative, q), m.size-1)]])
    return dict(number_m3=number, mass_kg_m3=mass, mean_mass_kg=mean,
                variance_mass_kg2=variance, relative_mass_variance=variance/mean**2,
                second_mass_moment_kg2_m3=float(np.dot(n, m*m)),
                mass_weighted_log_diameter=float(np.dot(m*n, np.log(d))/mass),
                mass_d50_m=quantile(.5), mass_d90_m=quantile(.9), mass_d99_m=quantile(.99),
                mass_d999_m=quantile(.999),
                mass_share_above_diameter={str(edge):float(np.dot(m[d>edge], n[d>edge])/mass)
                                          for edge in (.001, .003, .005, .01, .15, .6)})


def characteristic_product(result, geometry, *, id, endpoint_tracers=False):
    n, m = result['numbers'], result['masses']
    initial = float(np.dot(n, m[0]))
    rows = []
    for t, mass, deposited in zip(result['times'], m, result['deposited_mass']):
        row = moments(mass, n, geometry)
        row['time_s'] = float(t)
        row['support_minimum_m'] = float(geometry.size(mass[-2] if endpoint_tracers else mass[n>0].min()))
        row['support_maximum_m'] = float(geometry.size(mass[-1] if endpoint_tracers else mass[n>0].max()))
        row['water_ledger'] = dict(initial=initial, deposited=float(deposited), remaining=row['mass_kg_m3'],
                                   residual=initial+float(deposited)-row['mass_kg_m3'])
        rows.append(row)
    return dict(id=id, rtol=result['rtol'], atol=result['atol'], evaluations=result['evaluations'],
                wall_seconds=result['wall_seconds'], snapshots=rows,
                final_masses_kg=m[-1].tolist(), number_weights_m3=n.tolist())


def physical_seeds(geometry, nodes, concentration):
    x, weights = np.polynomial.legendre.leggauss(nodes)
    lo, hi = np.log(ag.START_M)
    diameters = np.exp((lo+hi)/2 + (hi-lo)*x/2)
    # Zero-weight endpoints give the support of the continuous input, rather
    # than the slightly narrower quadrature-node range.
    return np.r_[geometry.mass(diameters), geometry.mass(np.asarray(ag.START_M))], np.r_[concentration*weights/2, 0., 0.]


def identity(inputs):
    return dict(files={p:sha(ROOT/p) for p in FILES}, inputs={p:sha(ROOT/p) for p in inputs},
                constants=constants_used(FILES), numpy=np.__version__, scipy=scipy.__version__)


def analyse(*, budget_s=180.):
    start = time.monotonic(); deadline = start+budget_s
    world, binding = sf.moon_column()
    historical_path = 'atmosphere/snow/results/repair.json'
    historical = json.loads((ROOT/historical_path).read_text())
    control = next(r for r in historical['runs'] if r['id']=='moon_heavy_m144')
    duration = float(control['diagnostics']['final_time_s'])
    times = np.array([0., 6*3600., 24*3600., 72*3600., duration])
    geometry = GEOMETRIES['bf95_dmax_dendrites_area']
    temperature_c, concentration, source_rate = -10., 1e5, 3.
    depth = world.height_of(sf.T0-15.)-world.z_freeze
    production = source_rate/3600./depth
    state = world.state(world.height_of(sf.T0+temperature_c))
    def growth(m):
        diameter = geometry.size(m)
        velocity, reynolds = ag.speeds(diameter, m, world.gravity, state, geometry.area_ratio)
        return np.maximum(ag.deposition(diameter, reynolds, state, geometry.capacitance_ratio), 0.)
    inputs = dict(binding.get('inputs', {}))
    inputs[historical_path] = sha(ROOT/historical_path)
    before = identity(inputs)
    if before['inputs'] != inputs:
        raise RuntimeError('Loaded climate/column input hashes do not match the intended bytes')
    physical = []
    for nodes, tolerance in ((64, 1e-9), (128, 1e-9), (128, 1e-11)):
        m, n = physical_seeds(geometry, nodes, concentration)
        result = characteristics(m, n, times, growth, production, rtol=tolerance,
                                 atol=1e-24 if tolerance<1e-10 else 1e-22, deadline=deadline)
        physical.append(characteristic_product(result, geometry, id=f'q{nodes}_rtol{tolerance:g}', endpoint_tracers=True))
    grid_references, runs = [], []
    for bins in (96, 144, 216):
        pivots = np.geomspace(float(geometry.mass(2e-5)), float(geometry.mass(10.)), bins)
        initial = concentration*sectional.seed_distribution(pivots, geometry)
        occupied = initial>0
        reference = characteristics(pivots[occupied], initial[occupied], times, growth, production,
                                    rtol=1e-11, atol=1e-24, deadline=deadline)
        grid_references.append(characteristic_product(reference, geometry, id=f'projected_seed_m{bins}'))
        for step in ((30., 15., 7.5) if bins==144 else (30.,)):
            result = fixed_grid(pivots, initial, times, growth(pivots), production,
                                max_step_s=step, deadline=deadline)
            initial_mass = float(np.dot(pivots, initial))
            rows = []
            for k, (t, number) in enumerate(zip(times, result['numbers'])):
                exact = reference['masses'][k]
                once = project_once(pivots, exact, reference['numbers'])
                top = np.flatnonzero(once>0)[-1]
                row = moments(pivots, number, geometry)
                once_moments = moments(pivots, once, geometry)
                char_moments = moments(exact, reference['numbers'], geometry)
                deposited, exited = float(result['deposited_mass'][k]), float(result['upper_exit_mass'][k])
                row.update(time_s=float(t), once_projected_characteristic=once_moments,
                           characteristic=char_moments,
                           once_projected_support_maximum_m=float(geometry.size(pivots[top])),
                           mass_share_beyond_once_projected_support=float(np.dot(pivots[top+1:], number[top+1:])/row['mass_kg_m3']),
                           number_share_beyond_once_projected_support=float(number[top+1:].sum()/number.sum()),
                           excess_variance_over_once_projection_kg2=row['variance_mass_kg2']-once_moments['variance_mass_kg2'],
                           water_ledger=dict(initial=initial_mass, deposited=deposited, upper_exit=exited,
                                             remaining=row['mass_kg_m3'], residual=initial_mass+deposited-exited-row['mass_kg_m3']),
                           number_residual=float(initial.sum()-number.sum()-result['upper_exit_number'][k]))
                rows.append(row)
            runs.append(dict(id=f'm{bins}_dt{step:g}', mass_bins=bins, max_step_s=step,
                             diagnostics=result['diagnostics'], snapshots=rows,
                             mass_pivots_kg=pivots.tolist(), final_numbers_m3=result['numbers'][-1].tolist(),
                             initial_numbers_m3=initial.tolist()))
            print(f'growth box m{bins} dt{step:g}: {result["diagnostics"]["wall_seconds"]:.2f}s, '
                  f'final excess-tail mass={rows[-1]["mass_share_beyond_once_projected_support"]:.6g}', flush=True)
    reference = physical[-1]['snapshots'][-1]
    keys = ('mean_mass_kg', 'second_mass_moment_kg2_m3', 'support_maximum_m')
    reference_controls = []
    for row in physical[:-1]:
        metric = {key:abs(row['snapshots'][-1][key]/reference[key]-1.) for key in keys}
        reference_controls.append(dict(case=row['id'], reference=physical[-1]['id'], relative_differences=metric,
                                       passed=max(metric.values())<1e-6))
    base, half, quarter = [next(r for r in runs if r['id']==name)['snapshots'][-1]
                           for name in ('m144_dt30', 'm144_dt15', 'm144_dt7.5')]
    temporal = {}
    for key in ('relative_mass_variance', 'mass_share_beyond_once_projected_support'):
        spatial = (base[key]-base['once_projected_characteristic'][key]
                   if key=='relative_mass_variance' else base[key])
        last_change = abs(half[key]-quarter[key])
        temporal[key] = dict(coarse_to_half=abs(base[key]-half[key]), half_to_quarter=last_change,
                             measured_spatial_difference=abs(spatial),
                             last_time_change_fraction_of_spatial=last_change/abs(spatial) if spatial else None,
                             passed=last_change<.05*abs(spatial) if spatial else last_change==0)
    after = identity(inputs)
    if after != before:
        raise RuntimeError('Producer/input identity changed during comparison')
    return dict(schema=SCHEMA, producer=before,
                settings=dict(temperature_c=temperature_c, height_m=world.height_of(sf.T0+temperature_c),
                              pressure_pa=state['p'], source_cloud_top_c=-15., source_depth_m=depth,
                              nominal_vapour_supply_mm_h=source_rate, vapour_production_kg_m3_s=production,
                              initial_number_m3=concentration, seed_diameter_range_m=list(ag.START_M),
                              duration_s=duration, duration_source=historical_path,
                              initial_demand_supply_ratio=float(np.dot(physical_seeds(geometry,128,concentration)[1],
                                  growth(physical_seeds(geometry,128,concentration)[0]))/production),
                              geometry=geometry.id, budget_s=budget_s),
                physical_characteristics=physical, grid_matched_characteristics=grid_references, runs=runs,
                controls=dict(reference=reference_controls, time=temporal), wall_seconds=time.monotonic()-start,
                evidence='Controlled condensation-only fixed-air box. Continuous-mass characteristics isolate '
                         'deposition spreading from fixed mass pivots; no collisions, settling, new seed injection '
                         'or claim about the coupled column runaway. Analytical tests and numerical controls '
                         'establish numerical behavior, not empirical lunar snow validation.',
                reading_rule='Compare each grid with its own continuously grown initial projected seeds, '
                             'projected only once for the endpoint comparison. Tail mass beyond that support '
                             'excludes the unavoidable single bracketing-pivot representation. Relative time '
                             'changes must be small compared with measured spatial spreading. The closed box '
                             'accumulates water for the stated duration and is not a realistic cloud history.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--budget-s', type=float, default=180.)
    args = parser.parse_args()
    product = analyse(budget_s=args.budget_s)
    OUT.write_text(json.dumps(product, indent=1, default=float)+'\n')
    print(f'Wrote {OUT.relative_to(ROOT)} in {product["wall_seconds"]:.2f}s', flush=True)


if __name__ == '__main__':
    main()
