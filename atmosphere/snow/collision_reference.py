"""Independent continuous-mass collision references for the inherited snow law.

Exact-event finite-volume coagulation and prescribed-bath pure birth preserve
individual mass sums. The existing fixed-pivot operator supplies comparisons;
its producers and outputs remain unchanged. No deposition, settling or chipping
is included in either side of these explicitly aggregation-only controls.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.sparse import coo_matrix, diags
from scipy.sparse.linalg import expm_multiply

from atmosphere.snow import aggregation as ag, flakes as sf, sectional
from atmosphere.snow.habits import GEOMETRIES
from shared.provenance import constants_used

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT/'atmosphere/snow/results/collision_reference.json'
SCHEMA = 'terluna.atmosphere.snow-collision-reference/1'
GEOMETRY = GEOMETRIES['bf95_dmax_dendrites_area']
FILES = ('atmosphere/snow/collision_reference.py', 'atmosphere/snow/sectional.py',
         'atmosphere/snow/aggregation.py', 'atmosphere/snow/habits.py',
         'atmosphere/snow/flakes.py', 'atmosphere/rain/drops.py',
         'atmosphere/rain/targets.py', 'atmosphere/radiative_convective/thermodynamics.py')


def identity(binding):
    sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
    inputs = {p: sha(ROOT/p) for p in binding.get('inputs', {})}
    if inputs != binding.get('inputs', {}):
        raise RuntimeError('Loaded air-column inputs changed')
    return dict(files={p: sha(ROOT/p) for p in FILES}, inputs=inputs,
                constants=constants_used(FILES))


class SnowKernel:
    """Unchanged geometry, drag and temperature-only sticking, at fixed air."""
    def __init__(self, world, temperature_c=-10.):
        self.state = world.state(world.height_of(sf.T0+temperature_c))
        self.gravity = world.gravity
        self.sticking = math.exp(.09*temperature_c)

    def __call__(self, left, right):
        left, right = np.asarray(left), np.asarray(right)
        dl, dr = GEOMETRY.size(left), GEOMETRY.size(right)
        vl = ag.speeds(dl, left, self.gravity, self.state, GEOMETRY.area_ratio)[0]
        vr = ag.speeds(dr, right, self.gravity, self.state, GEOMETRY.area_ratio)[0]
        return self.sticking*GEOMETRY.area_ratio*math.pi/4*(dl+dr)**2*np.abs(vl-vr)


def _checked_kernel(kernel, left, right):
    result = np.asarray(kernel(left, right), dtype=float)
    result = np.broadcast_to(result, np.broadcast_shapes(np.shape(left), np.shape(right))).copy()
    if not np.isfinite(result).all() or (result < 0).any():
        raise ValueError('Pair rates must be finite and nonnegative')
    return result


def coagulate(initial_masses, volume_m3, observation_times_s, kernel, *, seed=0):
    """Direct Marcus–Lushnikov events; each unordered pair has rate K/V.

    Each particle represents one particle in the finite volume. No weighted
    packets, mass grid, collision batches or time steps enter this algorithm.
    Monomer counts add with mass. Rate sums are refreshed every 32 events;
    negative roundoff sums trigger recomputation, never clipping.
    """
    mass = np.asarray(initial_masses, dtype=float).copy()
    times = np.asarray(observation_times_s, dtype=float)
    if (mass.ndim != 1 or not mass.size or not np.isfinite(mass).all() or (mass <= 0).any()
            or not math.isfinite(volume_m3) or volume_m3 <= 0
            or times.ndim != 1 or not times.size or not np.isfinite(times).all()
            or (times < 0).any() or (np.diff(times) < 0).any()):
        raise ValueError('Positive finite masses/volume and ordered nonnegative times required')
    rates = _checked_kernel(kernel, mass[:, None], mass[None, :])
    if not np.allclose(rates, rates.T, rtol=1e-13, atol=0):
        raise ValueError('Pair kernel must be symmetric')
    np.fill_diagonal(rates, 0.)  # distinct particles only, even for K(m,m)>0
    row_sum = rates.sum(axis=1)
    monomers = np.ones(mass.size, dtype=np.int64)
    rng = np.random.default_rng(seed)
    initial_mass = float(mass.sum())
    active, events, t, target = mass.size, 0, 0., 0
    snapshots = []
    maximum_rate_sum_discrepancy = 0.
    while target < times.size:
        if events % 32 == 0 or (row_sum[:active] < 0).any():
            exact = rates[:active, :active].sum(axis=1)
            denominator = float(exact.sum())
            if denominator:
                maximum_rate_sum_discrepancy = max(maximum_rate_sum_discrepancy,
                    float(np.sum(np.abs(exact-row_sum[:active]))/denominator))
            row_sum[:active] = exact
        cumulative = np.cumsum(row_sum[:active])
        total = float(cumulative[-1])/(2*volume_m3)
        event_time = t+rng.exponential(1/total) if total > 0 else math.inf
        while target < times.size and times[target] < event_time:
            snapshots.append(dict(time_s=float(times[target]), masses_kg=mass[:active].copy(),
                                  monomer_counts=monomers[:active].copy(), events=events))
            target += 1
        if target == times.size:
            break
        i = int(np.searchsorted(cumulative, rng.random()*cumulative[-1], side='right'))
        partner_cumulative = np.cumsum(rates[i, :active])
        if partner_cumulative[-1] <= 0:
            raise FloatingPointError('Roundoff selected a zero-rate row; no event admitted')
        j = int(np.searchsorted(partner_cumulative, rng.random()*partner_cumulative[-1], side='right'))
        if i == j:
            raise AssertionError('Self-collision selected')
        i, j = min(i, j), max(i, j)
        new_mass = mass[i]+mass[j]
        new_monomers = monomers[i]+monomers[j]
        row_sum[:active] -= rates[:active, i]+rates[:active, j]
        last = active-1
        if j != last:
            mass[j], monomers[j] = mass[last], monomers[last]
            rates[j, :active] = rates[last, :active]
            rates[:active, j] = rates[:active, last]
            row_sum[j] = row_sum[last]
        active -= 1
        mass[i], monomers[i] = new_mass, new_monomers
        new_rates = _checked_kernel(kernel, new_mass, mass[:active])
        new_rates[i] = 0.
        rates[i, :active] = new_rates
        rates[:active, i] = new_rates
        row_sum[:active] += new_rates
        row_sum[i] = float(new_rates.sum())
        events += 1
        t = event_time
    mass_residual = max(abs(float(s['masses_kg'].sum())-initial_mass)/initial_mass for s in snapshots)
    if mass_residual > 1e-12 or maximum_rate_sum_discrepancy > 1e-11:
        raise FloatingPointError('Exact-event mass/rate bookkeeping tolerance failed')
    return dict(snapshots=snapshots, events=events, volume_m3=volume_m3,
                maximum_relative_mass_residual=mass_residual,
                maximum_rate_sum_relative_discrepancy=maximum_rate_sum_discrepancy)


def pure_birth(initial_mass, captured_mass, bath_number_m3, times, kernel, *,
               samples=4096, seed=0, maximum_events=100000, sample_block=32, event_block=2048):
    """Exact independent capture waiting times; no diffusion approximation.

    Before the event cap is reached, rates are bounded by their computed
    maximum. A dominating Poisson count supplies a Chernoff upper bound on
    hitting that computational cap. Samples reaching it cause failure.
    """
    times = np.asarray(times, dtype=float)
    if (initial_mass <= 0 or captured_mass <= 0 or bath_number_m3 <= 0 or samples < 2
            or maximum_events < 1 or (times < 0).any() or (np.diff(times) < 0).any()):
        raise ValueError('Positive capture parameters and ordered times required')
    path = initial_mass+captured_mass*np.arange(maximum_events)
    rates = bath_number_m3*_checked_kernel(kernel, path, captured_mass)
    if (rates <= 0).any():
        raise ValueError('Capture reference requires strictly positive event rates')
    poisson_mean = float(rates.max()*times[-1])
    log_bound = (-maximum_events*math.log(maximum_events/poisson_mean)+maximum_events-poisson_mean
                 if 0 < poisson_mean < maximum_events else (float('-inf') if poisson_mean == 0 else 0.))
    if log_bound > math.log(1e-12):
        raise ValueError('Event-domain hitting probability is insufficiently bounded')
    rng = np.random.default_rng(seed)
    counts = np.zeros((times.size, samples), dtype=np.int64)
    for start in range(0, samples, sample_block):
        size = min(sample_block, samples-start)
        ages = np.zeros(size)
        for offset in range(0, maximum_events, event_block):
            stop = min(offset+event_block, maximum_events)
            waits = np.cumsum(rng.exponential(size=(size, stop-offset))/rates[offset:stop], axis=1)
            waits += ages[:, None]
            for it, end in enumerate(times):
                counts[it, start:start+size] += np.sum(waits <= end, axis=1)
            ages = waits[:, -1]
            if np.all(ages > times[-1]):
                break
        else:
            raise RuntimeError('Capture event domain exhausted; no result admitted')
    return dict(counts=counts, masses_kg=initial_mass+captured_mass*counts,
                event_cap=maximum_events, maximum_events_observed=int(counts.max()),
                dominating_poisson_mean=poisson_mean, log_event_cap_probability_bound=log_bound,
                event_cap_probability_upper_bound=math.exp(log_bound), samples=samples)


def bath_fixed_pivot(masses, initial_mass, captured_mass, bath_number_m3, times, kernel):
    """Exponentiate the frozen allocation with a prescribed replenished bath.

    Only collector/bath encounters are enabled. Their bath loss is replenished;
    the captured-water source is integrated separately from domain exits.
    """
    mass = np.asarray(masses)
    i = int(np.flatnonzero(mass == captured_mass)[0])
    j0 = int(np.flatnonzero(mass == initial_mass)[0])
    pair = np.zeros((mass.size, mass.size))
    js = np.arange(i+1, mass.size)
    pair[i, js] = pair[js, i] = kernel(mass[js], captured_mass)
    operator = sectional.CollisionOperator(mass, pair)
    columns = np.flatnonzero((operator.i == i) & (operator.j > i))
    js = operator.j[columns]
    rate = operator.rate_coefficient[columns]*bath_number_m3
    transfer = (operator.transfer[:, columns] @ diags(rate)).tocoo()
    keep = transfer.row != i
    rows = list(transfer.row[keep]); cols = list(js[transfer.col[keep]]); values = list(transfer.data[keep])
    for extra, coefficients in enumerate((operator.exit_number[1, columns]*rate,
                                         operator.exit_mass[1, columns]*rate,
                                         captured_mass*rate)):
        rows.extend([mass.size+extra]*len(js)); cols.extend(js); values.extend(coefficients)
    matrix = coo_matrix((values, (rows, cols)), shape=(mass.size+3, mass.size+3)).tocsr()
    initial = np.zeros(mass.size+3); initial[j0] = 1.
    outputs = []
    for end in times:
        state = expm_multiply(matrix*end, initial)
        split = expm_multiply(matrix*(end/2), expm_multiply(matrix*(end/2), initial))
        if (state < 0).any():
            raise FloatingPointError('Negative fixed-pivot exponential state')
        number = state[:-3]
        mass_residual = float(np.dot(number, mass)+state[-2]-initial_mass-state[-1])
        outputs.append(dict(time_s=float(end), numbers=number,
                            upper_exit_number=float(state[-3]), upper_exit_mass_kg=float(state[-2]),
                            supplied_capture_mass_kg=float(state[-1]),
                            relative_mass_residual=mass_residual/(initial_mass+state[-1]),
                            number_residual=float(number.sum()+state[-3]-1.),
                            exponential_composition_error=float(np.max(np.abs(state-split)))))
    return outputs


def _project(masses, particles, volume):
    number = np.zeros_like(masses)
    for mass in particles:
        allocation, side = sectional._allocation(masses, mass)
        if side:
            raise ValueError('Initial particles outside the fixed mass domain')
        for i, weight in allocation:
            number[i] += weight/volume
    return number


def _seed_masses(count):
    # Deterministic midpoint quadrature of the original log-uniform seed law.
    return GEOMETRY.mass(ag.START_M[0]*(ag.START_M[1]/ag.START_M[0])**((np.arange(count)+.5)/count))


def _summary(masses, number, query):
    mass, number = np.asarray(masses), np.asarray(number)
    d = GEOMETRY.size(mass)
    n, m1 = float(number.sum()), float(np.dot(mass, number))
    m2 = float(np.dot(mass*mass, number))
    order = np.argsort(d); cumulative = np.cumsum((mass*number)[order])
    quantiles = {str(p): float(d[order[min(np.searchsorted(cumulative, p*m1), mass.size-1)]])
                 for p in (.5, .9, .99)}
    locations = np.searchsorted(d[order], query, side='right')
    return dict(number_m3=n, mass_kg_m3=m1, second_mass_moment_kg2_m3=m2,
                relative_mass_variance=m2*n/(m1*m1)-1., mass_quantiles_diameter_m=quantiles,
                mass_cdf=(np.r_[0., cumulative][locations]/m1).tolist(),
                number_below_m3=np.r_[0., np.cumsum(number[order])][locations].tolist())


def _plain(value):
    if isinstance(value, np.ndarray): return value.tolist()
    if isinstance(value, np.generic): return value.item()
    if isinstance(value, dict): return {k: _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)): return [_plain(v) for v in value]
    return value


def nested_grid(factor=1):
    base = np.geomspace(float(GEOMETRY.mass(2e-5)), float(GEOMETRY.mass(10.)), 96)
    if factor == 1:
        return base
    return np.r_[np.concatenate([np.geomspace(a, b, factor+1)[:-1]
                                 for a, b in zip(base[:-1], base[1:])]), base[-1]]


def ordinary_initial():
    grid = nested_grid()
    diameter = GEOMETRY.size(grid)
    return grid[(diameter >= ag.START_M[0]) & (diameter <= ag.START_M[1])]


def ensemble_box(kernel, particle_count, repetitions, times, query, *, seed=41021):
    species = ordinary_initial()
    if particle_count % species.size:
        raise ValueError('Equal integer counts per initial mass species required')
    initial = np.repeat(species, particle_count//species.size)
    volume = particle_count/1e5
    moments = ('number_m3', 'mass_kg_m3', 'second_mass_moment_kg2_m3')
    rows = [[] for _ in times]
    max_mass_error = max_rate_error = 0.
    for repeat in range(repetitions):
        result = coagulate(initial, volume, times, kernel, seed=seed+particle_count*1000+repeat)
        max_mass_error = max(max_mass_error, result['maximum_relative_mass_residual'])
        max_rate_error = max(max_rate_error, result['maximum_rate_sum_relative_discrepancy'])
        for index, snapshot in enumerate(result['snapshots']):
            row = _summary(snapshot['masses_kg'], np.full(len(snapshot['masses_kg']), 1/volume), query)
            row['maximum_particle_mass_kg'] = float(snapshot['masses_kg'].max())
            rows[index].append(row)
    summaries = []
    for end, samples in zip(times, rows):
        mean = {key: float(np.mean([r[key] for r in samples])) for key in moments}
        stderr = {key: float(np.std([r[key] for r in samples], ddof=1)/math.sqrt(repetitions)) for key in moments}
        cdfs = np.array([r['mass_cdf'] for r in samples])
        number_below = np.array([r['number_below_m3'] for r in samples])
        summaries.append(dict(time_s=float(end), mean=mean, standard_error=stderr,
            relative_mass_variance_of_mean_population=mean[moments[2]]*mean[moments[0]]/mean[moments[1]]**2-1.,
            mass_cdf_mean=cdfs.mean(axis=0), mass_cdf_standard_error=cdfs.std(axis=0, ddof=1)/math.sqrt(repetitions),
            number_below_mean_m3=number_below.mean(axis=0),
            maximum_particle_mass_kg_per_realization=[r['maximum_particle_mass_kg'] for r in samples]))
    return dict(initial_particles=particle_count, independent_realizations=repetitions,
                volume_m3=volume, initial_masses_kg=species,
                initial_number_m3_per_species=np.full(species.size, 1e5/species.size),
                maximum_possible_particle_diameter_m=float(GEOMETRY.size(initial.sum())),
                pointwise_95_percent_hoeffding_halfwidth_mass_cdf=math.sqrt(math.log(40)/(2*repetitions)),
                maximum_relative_mass_residual=max_mass_error,
                maximum_rate_sum_relative_discrepancy=max_rate_error, snapshots=summaries)


def fixed_box(kernel, factor, times, query, *, max_step_s=30.):
    masses = nested_grid(factor)
    initial = ordinary_initial()
    # All initial species are exactly present in every nested mass grid.
    numbers = _project(masses, initial, len(initial)/1e5)
    operator = sectional.CollisionOperator(masses, kernel(masses[:, None], masses[None, :]))
    snapshots = []
    exit_mass = exit_number = 0.
    steps = retries = 0
    previous = 0.
    source_mass = float(np.dot(numbers, masses))
    for end in times:
        if end > previous:
            result = sectional.integrate_box(operator, numbers, end-previous, max_step_s=max_step_s)
            numbers = result['numbers']
            exit_mass += result['upper_exit_mass']+result['lower_exit_mass']
            exit_number += result['upper_exit_number']+result['lower_exit_number']
            steps += result['steps']; retries += result['rejected_steps']
        row = _summary(masses, numbers, query)
        row.update(time_s=float(end), numbers_m3=numbers.copy(), domain_exit_mass_kg_m3=exit_mass,
                   domain_exit_number_m3=exit_number, domain_mass_fraction=exit_mass/source_mass,
                   relative_mass_residual=(row['mass_kg_m3']+exit_mass-source_mass)/source_mass)
        snapshots.append(row)
        previous = end
    return dict(mass_bins=masses.size, refinement_factor=factor, max_step_s=max_step_s,
                mass_pivots_kg=masses, diameters_m=GEOMETRY.size(masses),
                initial_masses_kg=initial, steps=steps, rejected_steps=retries, snapshots=snapshots)


def capture_case(kernel):
    base = nested_grid()
    d = GEOMETRY.size(base)
    small = float(base[np.argmin(np.abs(np.log(d/.0002)))])
    collector = float(base[np.argmin(np.abs(np.log(d/.15)))])
    bath = 1e5
    initial_rate = float(kernel(collector, small)*bath)
    times = np.array([0., .03, .1])*collector/(small*initial_rate)
    event_cap = int(math.ceil(.4*collector/small))
    references = []
    for samples in (4096, 8192):
        birth = pure_birth(collector, small, bath, times, kernel, samples=samples,
                           seed=91021+samples, maximum_events=event_cap)
        snapshots = []
        for end, masses in zip(times, birth['masses_kg']):
            unique, counts = np.unique(masses, return_counts=True)
            snapshots.append(dict(time_s=float(end), masses_kg=unique, sample_counts=counts,
                mean_mass_kg=float(masses.mean()), mean_mass_standard_error_kg=float(masses.std(ddof=1)/math.sqrt(samples)),
                variance_mass_kg2=float(masses.var(ddof=1)),
                mass_quantiles_kg={str(p): float(np.quantile(masses, p, method='inverted_cdf')) for p in (.5, .9, .99)},
                supplied_mean_capture_mass_kg=float(masses.mean()-collector)))
        references.append(dict(samples=samples, independent_seed=91021+samples,
            number_cdf_uniform_95_percent_dkw_halfwidth=math.sqrt(math.log(40)/(2*samples)),
            event_cap=birth['event_cap'], maximum_events_observed=birth['maximum_events_observed'],
            dominating_poisson_mean=birth['dominating_poisson_mean'],
            log_event_cap_probability_bound=birth['log_event_cap_probability_bound'],
            event_cap_probability_upper_bound=birth['event_cap_probability_upper_bound'], snapshots=snapshots))
    grids = []
    for factor in (1, 2, 4):
        masses = nested_grid(factor)
        results = bath_fixed_pivot(masses, collector, small, bath, times, kernel)
        rows = []
        for row in results:
            numbers = row.pop('numbers')
            summary = _summary(masses, numbers, [])
            rows.append(dict(**row, numbers_per_initial_collector=numbers, **summary))
        grids.append(dict(mass_bins=masses.size, mass_pivots_kg=masses,
                          diameters_m=GEOMETRY.size(masses), snapshots=rows))
    return dict(initial_collector_mass_kg=collector, captured_mass_kg=small,
        initial_collector_diameter_m=float(GEOMETRY.size(collector)),
        captured_diameter_m=float(GEOMETRY.size(small)), bath_number_m3=bath,
        initial_capture_rate_s=initial_rate, times_s=times,
        duration_rule='Chosen so initial-rate capture would add 3% and 10% of collector mass. Actual state-dependent rates remain exact.',
        source_rule='Prescribed replenished small-particle bath; each collector capture adds exactly the fixed partner mass. This is not a closed column water budget.',
        continuous_mass_references=references, fixed_pivot=grids)


def analyse(*, repetitions=64):
    started = time.monotonic()
    world, binding = sf.moon_column()
    before = identity(binding)
    kernel = SnowKernel(world)
    times = np.array([0., 2., 4., 6., 12.])*3600.
    query = np.unique(np.r_[np.geomspace(5e-5, .1, 301), .001, .003, .005, .01, .15, .6, 1., 10.])
    print('Exact prescribed-bath captures', flush=True)
    capture = capture_case(kernel)
    ensemble = []
    species_count = len(ordinary_initial())
    for multiple in (32, 64, 128, 256):
        count = species_count*multiple
        ensemble.append(ensemble_box(kernel, count, repetitions, times, query))
        print(f'Exact coagulation N={count}, R={repetitions}; {time.monotonic()-started:.1f} s', flush=True)
    fixed = []
    for factor in (1, 2, 4):
        fixed.append(fixed_box(kernel, factor, times, query))
        print(f'Fixed-pivot {fixed[-1]["mass_bins"]}; {time.monotonic()-started:.1f} s', flush=True)
    time_control = fixed_box(kernel, 2, times, query, max_step_s=15.)
    if identity(binding) != before:
        raise RuntimeError('Collision-reference sources or air inputs changed; no result admitted')
    return _plain(dict(schema=SCHEMA, producer=before,
        settings=dict(temperature_c=-10., geometry=GEOMETRY.to_dict(),
                      air_state={k: float(v) for k, v in kernel.state.items()},
                      sticking=kernel.sticking, initial_number_m3=1e5,
                      aggregation=True, chipping=False, deposition=False, settling=False,
                      simulation_seconds=times[-1], wall_seconds=time.monotonic()-started),
        capture=capture, ordinary_box=dict(cdf_query_diameter_m=query, times_s=times,
                      exact_event_ensembles=ensemble, fixed_pivot=fixed, timestep_control=time_control),
        evidence='Independent finite-volume exact-event coagulation and exact prescribed-bath capture processes, compared with frozen fixed-pivot allocation. No time or mass discretization enters collision events; statistical and finite-volume errors remain distinct.',
        reading_rule='The bath samples are independent collector histories. Their DKW band applies to number CDF, not mass-weighted CDF. Closed-box standard errors use independent entire realizations, not correlated particles within one realization. A smaller Monte Carlo standard error does not remove finite-volume bias. The 12 h point is a stress control and requires its own volume-convergence assessment. Fixed-pivot domain exits and finite-volume maximum particle masses are numerical limitations, not physical snow limits.',
        initial_state_rule='All ordinary runs start from the same equally weighted physical seed masses on the base grid within 0.1–0.4 mm. Nested refinements retain those masses exactly. This is a declared discrete seed distribution, not an assertion that it exactly equals the continuous log-uniform injection law.',
        chipping_scope='Chipping is disabled in both algorithms to isolate aggregation. Reproducing the inherited nonsticking channel later requires exact smaller-parent remnants and a coherent stochastic or deterministic treatment of fractional expected chip counts; rounding them to whole particles would introduce a new model change.',
        remaining=['Finite-volume/statistical assessment at each time and for each tail observable.',
                   'Collision-grid error with the inherited smaller-parent chipping channel.',
                   'Coupled deposition, settling and climate-supported moisture histories.',
                   'Geometry, hydrodynamic interception and whole-aggregate fracture calibration.']
    ))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--repetitions', type=int, default=64)
    args = parser.parse_args()
    if args.repetitions < 2:
        parser.error('At least two independent repetitions required')
    result = analyse(repetitions=args.repetitions)
    OUT.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(f'Wrote {OUT.relative_to(ROOT)}', flush=True)


if __name__ == '__main__':
    main()
