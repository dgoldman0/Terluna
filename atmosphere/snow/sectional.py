"""Conservative sectional reference for the audited snow equations.

This independent solver retains the prescribed column, geometry, vapour supply,
settling law, sticking efficiency and smaller-particle chipping law. Number
concentrations are real densities, never integer superparticle multiplicities.
Binary products are allocated to adjacent mass pivots with weights preserving
number and mass. Sticking and breakup are separate mean-rate reaction channels;
there are no multi-encounter Bernoulli batches. The historical simulator and its
outputs remain untouched.

The allocation is the fixed-pivot method of Kumar and Ramkrishna (1996),
doi:10.1016/0009-2509(96)88489-2. Their Part II explicitly warns that coarse
geometric grids overpredict steep large-size fronts; mass-grid refinement and
upper-domain flux checks are required here. This solver verifies equations, not
large-aggregate physics. A superlinear extrapolated kernel can itself transport
mass outside every chosen finite domain.

Time integration is adaptive SSPRK2. Every forward-Euler stage must stay
nonnegative; rejection changes the step, never the population. Conservative
net transfers remove cancellation when a small collector increment mostly stays
in its original pivot. Deposition is conservative upwind transport in mass and
settling is conservative upwind transport in height. Their CFL limits and the
negative net tendency control the step. Domain exits are explicit sinks with
number AND mass ledgers, never added to snowfall. The absorbing upper ghost
pivot continues the geometric mass ratio, allowing conservative partial exits
from the last resolved pivot without an artificial gross-collision timescale. They bound numerical-domain
adequacy and are not prescribed physical size limits.
"""
from __future__ import annotations

import math

import numpy as np
from scipy.sparse import coo_matrix

from atmosphere.snow import aggregation as ag
from atmosphere.snow import flakes as sf
from atmosphere.snow.habits import GEOMETRIES

SCHEMA = 'terluna.atmosphere.snow-sectional/1'
LEDGER_KEYS = ('injected', 'deposited', 'freezing_exit', 'lower_domain_exit', 'upper_domain_exit')


def _grid(masses):
    m = np.asarray(masses, dtype=float)
    if m.ndim != 1 or m.size < 2 or not np.isfinite(m).all() or (m <= 0).any() or (np.diff(m) <= 0).any():
        raise ValueError('mass pivots must be finite, positive and strictly increasing')
    return m


def _allocation(masses, mass):
    """Number weights with exact target first moment; outside products exit."""
    if mass < masses[0]:
        return (), 'lower'
    if mass > masses[-1]:
        return (), 'upper'
    hi = int(np.searchsorted(masses, mass, side='left'))
    if masses[hi] == mass:
        return ((hi, 1.0),), None
    lo = hi - 1
    high = (mass - masses[lo]) / (masses[hi] - masses[lo])
    return ((lo, 1.0 - high), (hi, high)), None


class CollisionOperator:
    """One homogeneous cell's binary population balance on fixed mass pivots.

    kernel is symmetric m3/s. loss_per_breakup, when provided, is the kg chipped
    from the smaller parent in one nonsticking encounter; fragments have the
    supplied mass. This generic interface also supports analytic kernels.
    """

    def __init__(self, masses, kernel, *, sticking=1.0, loss_per_breakup=None, fragment_mass=None):
        self.masses = _grid(masses)
        n = self.masses.size
        self.upper_ghost_mass = self.masses[-1] ** 2 / self.masses[-2]
        kernel = np.asarray(kernel, dtype=float)
        if kernel.shape != (n, n) or not np.isfinite(kernel).all() or (kernel < 0).any():
            raise ValueError('kernel must be a finite nonnegative square matrix')
        if not np.allclose(kernel, kernel.T, rtol=1e-13, atol=0):
            raise ValueError('kernel must be symmetric')
        if not 0 <= sticking <= 1:
            raise ValueError('sticking must be in [0, 1]')
        self.i, self.j = np.triu_indices(n)
        self.rate_coefficient = kernel[self.i, self.j] * np.where(self.i == self.j, .5, 1.)
        count = self.i.size
        rows, cols, vals = [], [], []
        self.exit_number = np.zeros((2, count))
        self.exit_mass = np.zeros((2, count))
        if loss_per_breakup is not None:
            loss = np.asarray(loss_per_breakup, dtype=float)
            if loss.shape != (n, n) or not np.isfinite(loss).all() or (loss < 0).any():
                raise ValueError('breakup losses must be a finite nonnegative matrix')
            if fragment_mass is None or not math.isfinite(fragment_mass) or fragment_mass <= 0:
                raise ValueError('fragment mass must be positive')
            if (loss[self.i, self.j] > self.masses[self.i]).any():
                raise ValueError('breakup cannot remove more than the smaller parent')
        else:
            loss = None

        def put(row, col, value):
            if value:
                rows.append(row); cols.append(col); vals.append(value)

        def add_product(col, mass, number):
            if self.masses[-1] < mass < self.upper_ghost_mass:
                fraction = (mass - self.masses[-1]) / (self.upper_ghost_mass - self.masses[-1])
                put(n - 1, col, number * (1.0 - fraction))
                self.exit_number[1, col] += number * fraction
                self.exit_mass[1, col] += number * fraction * self.upper_ghost_mass
                return
            allocation, side = _allocation(self.masses, mass)
            if side:
                index = 0 if side == 'lower' else 1
                self.exit_number[index, col] += number
                self.exit_mass[index, col] += number * mass
            else:
                for row, weight in allocation:
                    put(row, col, number * weight)

        for col, (i, j) in enumerate(zip(self.i, self.j)):
            mi, mj = self.masses[i], self.masses[j]
            if sticking:
                product = mi + mj
                # For small-on-large events, compute the small transferred
                # fraction directly, avoiding subtracting nearly equal ones.
                next_mass = self.masses[j + 1] if j < n - 1 else self.upper_ghost_mass
                if product < next_mass:
                    fraction = mi / (next_mass - mj)
                    put(i, col, -sticking)
                    put(j, col, -sticking * fraction)
                    if j < n - 1:
                        put(j + 1, col, sticking * fraction)
                    else:
                        self.exit_number[1, col] += sticking * fraction
                        self.exit_mass[1, col] += sticking * fraction * next_mass
                else:
                    put(i, col, -sticking); put(j, col, -sticking)
                    add_product(col, product, sticking)
            if loss is not None and sticking < 1 and loss[i, j] > 0:
                chipped = loss[i, j]
                factor = 1.0 - sticking
                remnant = mi - chipped
                # The larger collision partner remains unchanged. For a small
                # chip, cancel the remnant's return to the same bin exactly.
                if i > 0 and remnant > self.masses[i - 1]:
                    fraction = chipped / (mi - self.masses[i - 1])
                    put(i, col, -factor * fraction)
                    put(i - 1, col, factor * fraction)
                else:
                    put(i, col, -factor)
                    if remnant > 0:
                        add_product(col, remnant, factor)
                add_product(col, fragment_mass, factor * chipped / fragment_mass)
        self.transfer = coo_matrix((vals, (rows, cols)), shape=(n, count)).tocsr()
        self.transfer.sum_duplicates()
        self.transfer.eliminate_zeros()

    def evaluate(self, numbers):
        n = np.asarray(numbers, dtype=float)
        if n.shape != self.masses.shape or not np.isfinite(n).all() or (n < 0).any():
            raise ValueError('number densities must be finite and nonnegative')
        rates = self.rate_coefficient * n[self.i] * n[self.j]
        return np.asarray(self.transfer @ rates), self.exit_number @ rates, self.exit_mass @ rates


def _negative_step(state, derivative, cfl):
    dying = derivative < 0
    return float(cfl * np.min(state[dying] / -derivative[dying])) if dying.any() else math.inf


def _advance(state, evaluate, dt_max, cfl):
    """Positive SSPRK2 step; all diagnostic rates use its same quadrature."""
    derivative, diagnostic, physical_limit = evaluate(state)
    dt = min(dt_max, physical_limit, _negative_step(state, derivative, cfl))
    rejected = 0
    while True:
        if not math.isfinite(dt) or dt <= np.finfo(float).tiny:
            raise FloatingPointError('time step underflow; numerical/physical stiffness remains unresolved')
        stage = state + dt * derivative
        if (stage < 0).any() or not np.isfinite(stage).all():
            dt *= .5; rejected += 1; continue
        derivative2, diagnostic2, limit2 = evaluate(stage)
        stage2 = stage + dt * derivative2
        if dt > limit2 * (1. + 1e-12) or (stage2 < 0).any() or not np.isfinite(stage2).all():
            dt *= .5; rejected += 1; continue
        result = .5 * state + .5 * stage2
        return result, dt, .5 * (diagnostic + diagnostic2), .5 * (state + stage), rejected


def integrate_box(operator, initial, seconds, *, max_step_s=.01, cfl=.3):
    """Closed-cell benchmark with explicit mass and number boundary exits."""
    state = np.asarray(initial, dtype=float).copy()
    if state.shape != operator.masses.shape or (state < 0).any() or not np.isfinite(state).all():
        raise ValueError('initial population must be finite, nonnegative and match pivots')
    if not math.isfinite(seconds) or seconds <= 0 or not 0 < cfl <= 1 or max_step_s <= 0:
        raise ValueError('positive duration/step and 0<cfl<=1 required')
    def evaluate(n):
        dn, exits_n, exits_m = operator.evaluate(n)
        return dn, np.concatenate((exits_n, exits_m)), math.inf
    t, steps, rejected = 0., 0, 0
    exits = np.zeros(4)
    while t < seconds:
        state, dt, rate, _, retries = _advance(state, evaluate, min(max_step_s, seconds - t), cfl)
        exits += dt * rate
        t += dt; steps += 1; rejected += retries
    return dict(numbers=state, lower_exit_number=exits[0], upper_exit_number=exits[1],
                lower_exit_mass=exits[2], upper_exit_mass=exits[3], steps=steps, rejected_steps=rejected)


def seed_distribution(masses, geometry):
    """Unit-number log-uniform injected sizes, with grid-independent mean mass.

    Gauss integration is split at every pivot, so number allocation and the
    expected physical mass are preserved across grid refinements. The injection
    normalization remains N times the speed of the geometric-mean seed, exactly
    as in the historical forcing; it is not N times each seed's own speed.
    """
    masses = _grid(masses)
    start, stop = ag.START_M
    if masses[0] > geometry.mass(start) or masses[-1] < geometry.mass(stop):
        raise ValueError('mass grid must contain the complete injected size distribution')
    pivots = geometry.size(masses)
    breaks = np.unique(np.concatenate(([start], pivots[(pivots > start) & (pivots < stop)], [stop])))
    x, w = np.polynomial.legendre.leggauss(12)
    result = np.zeros_like(masses)
    normalization = math.log(stop / start)
    for lo, hi in zip(np.log(breaks[:-1]), np.log(breaks[1:])):
        logs = (lo + hi) / 2 + (hi - lo) / 2 * x
        weights = w * (hi - lo) / (2 * normalization)
        for mass, number in zip(geometry.mass(np.exp(logs)), weights):
            allocation, side = _allocation(masses, float(mass))
            if side:
                raise ArithmeticError('seed integration unexpectedly left the mass grid')
            for i, fraction in allocation:
                result[i] += number * fraction
    return result


def _snow_operator(masses, diameter, velocity, temperature, geometry, breakup):
    relative = np.abs(velocity[:, None] - velocity[None, :])
    kernel = geometry.area_ratio * math.pi / 4.0 * (diameter[:, None] + diameter[None, :]) ** 2 * relative
    sticking = math.exp(.09 * (temperature - sf.T0))
    fragment_mass = float(geometry.mass(ag.FRAGMENT_M))
    loss = None
    if breakup:
        small = np.minimum(masses[:, None], masses[None, :])
        energy = masses[:, None] * masses[None, :] / (masses[:, None] + masses[None, :]) * relative ** 2
        count = ag.fragments(small, energy, mass_to_size=geometry.size)
        loss = np.minimum(count * fragment_mass, .5 * small)
    return CollisionOperator(masses, kernel, sticking=sticking, loss_per_breakup=loss, fragment_mass=fragment_mass)


def column(world, top_c, number_per_m3, snowfall_mm_h, *, geometry=None, layers=30, mass_bins=96,
           minimum_diameter_m=2e-5, maximum_diameter_m=10., max_step_s=30., cfl=.3,
           fill=1.5, gather=1.5, breakup=True, horizontal_area_m2=1., progress=None,
           mass_grid=None, maximum_steps=2000000, domain_mass_tolerance=1e-6):
    """Solve the same prescribed column with a conservative sectional method.

    Geometry and atmosphere are frozen inputs; there is no added fracture cap,
    cloud lifetime, turbulent kernel, sublimation, or empirical calibration.
    Densities and all per-m2 outputs are independent of horizontal_area_m2; the
    extensive ledger supplies an explicit dimensional area-scaling check.
    """
    geometry = GEOMETRIES['legacy_bf95'] if geometry is None else geometry
    values = (top_c, number_per_m3, snowfall_mm_h, minimum_diameter_m, maximum_diameter_m,
              max_step_s, cfl, fill, gather, horizontal_area_m2, domain_mass_tolerance)
    if not all(math.isfinite(float(v)) for v in values):
        raise ValueError('column inputs must be finite')
    if top_c >= 0 or number_per_m3 <= 0 or snowfall_mm_h < 0 or layers < 1 or mass_bins < 2:
        raise ValueError('subzero cloud, positive seed supply, nonnegative vapour and valid grids required')
    if not 0 < minimum_diameter_m < maximum_diameter_m or max_step_s <= 0 or not 0 < cfl <= 1:
        raise ValueError('invalid mass domain, maximum timestep or CFL')
    if fill < 0 or gather <= 0 or horizontal_area_m2 <= 0 or domain_mass_tolerance < 0:
        raise ValueError('invalid duration, area or domain tolerance')
    masses = (_grid(mass_grid) if mass_grid is not None else
              np.geomspace(float(geometry.mass(minimum_diameter_m)), float(geometry.mass(maximum_diameter_m)), mass_bins))
    diameter = geometry.size(masses)
    z_top, z_freeze = world.height_of(sf.T0 + top_c), world.z_freeze
    depth = z_top - z_freeze
    dz = depth / layers
    edges = np.linspace(z_freeze, z_top, layers + 1)
    height = .5 * (edges[:-1] + edges[1:])
    states = [world.state(z) for z in height]
    velocity, growth, operators = [], [], []
    for s in states:
        v, re = ag.speeds(diameter, masses, world.gravity, s, geometry.area_ratio)
        velocity.append(v)
        growth.append(np.maximum(ag.deposition(diameter, re, s, geometry.capacitance_ratio), 0.))
        operators.append(_snow_operator(masses, diameter, v, s['t'], geometry, breakup))
    velocity, growth = np.asarray(velocity), np.asarray(growth)
    seed_d = math.sqrt(ag.START_M[0] * ag.START_M[1])
    u_in = float(ag.speeds(np.array([seed_d]), np.array([geometry.mass(seed_d)]), world.gravity,
                         world.state(z_top), geometry.area_ratio)[0][0])
    seed_flux = number_per_m3 * u_in * seed_distribution(masses, geometry)
    crossing = depth / u_in
    gather_start, end = fill * crossing, (fill + gather) * crossing
    production = snowfall_mm_h / 3600. / depth
    virtual_top = masses[-1] * masses[-1] / masses[-2]
    spacing = np.concatenate((np.diff(masses), [virtual_top - masses[-1]]))
    sedimentation_limit = cfl * dz / float(velocity.max())
    # Diagnostic vector: 5 mass rates, 2 domain number rates, boundary number
    # flux by pivot. All rates are per square metre, not totals of an area.
    diagnostic_size = len(LEDGER_KEYS) + 2 + masses.size
    injected_rate = float(np.dot(seed_flux, masses))

    def evaluate(number):
        dn = np.zeros_like(number)
        diagnostic = np.zeros(diagnostic_size)
        diagnostic[0] = injected_rate
        limit = sedimentation_limit
        for k, operator in enumerate(operators):
            change, exit_number, exit_mass = operator.evaluate(number[k])
            dn[k] = change
            diagnostic[3:5] += exit_mass * dz
            diagnostic[5:7] += exit_number * dz
        demand = np.sum(growth * number, axis=1)
        scale = np.minimum(1., np.divide(production, demand, out=np.ones_like(demand), where=demand > 0))
        rate = growth * scale[:, None]
        transfer_coefficient = rate / spacing
        # Only occupied cells constrain this linear upwind update; the second
        # SSP stage checks newly populated cells as well, without deleting them.
        occupied = number > 0
        if occupied.any():
            largest = float(transfer_coefficient[occupied].max())
            if largest:
                limit = min(limit, cfl / largest)
        transfers = number * transfer_coefficient
        dn -= transfers
        dn[:, 1:] += transfers[:, :-1]
        diagnostic[1] = float(np.sum(number * rate) * dz)
        diagnostic[4] += float(transfers[:, -1].sum() * virtual_top * dz)
        diagnostic[6] += float(transfers[:, -1].sum() * dz)
        flux = velocity * number
        dn -= flux / dz
        dn[:-1] += flux[1:] / dz
        dn[-1] += seed_flux / dz
        diagnostic[2] = float(np.dot(flux[0], masses))
        diagnostic[7:] = flux[0]
        return dn, diagnostic, limit

    number = np.zeros((layers, masses.size))
    ledger = np.zeros(len(LEDGER_KEYS))
    domain_number = np.zeros(2)
    arrivals = np.zeros_like(masses)
    occupation = np.zeros_like(number)
    gathered_domain_mass = np.zeros(2)
    t, steps, rejected = 0., 0, 0
    smallest_dt, largest_dt = math.inf, 0.
    while t < end:
        if steps >= maximum_steps:
            raise RuntimeError('maximum_steps reached; unresolved numerical stiffness, no completed product')
        boundary = gather_start if t < gather_start else end
        dt_max = min(max_step_s, boundary - t)
        next_number, dt, diagnostic, mean_stage, retries = _advance(number, evaluate, dt_max, cfl)
        ledger += dt * diagnostic[:5]
        domain_number += dt * diagnostic[5:7]
        if t >= gather_start:
            arrivals += dt * diagnostic[7:]
            occupation += dt * mean_stage
            gathered_domain_mass += dt * diagnostic[3:5]
        number = next_number
        t += dt; steps += 1; rejected += retries
        smallest_dt, largest_dt = min(smallest_dt, dt), max(largest_dt, dt)
        if progress is not None:
            progress(t, steps)
    remaining = float(np.sum(number * masses) * dz)
    mass_ledger = dict(zip(LEDGER_KEYS, map(float, ledger)))
    mass_ledger['remaining'] = remaining
    mass_ledger['residual'] = float(ledger[0] + ledger[1] - ledger[2] - ledger[3] - ledger[4] - remaining)
    source = ledger[0] + ledger[1]
    mass_ledger['relative_residual'] = mass_ledger['residual'] / source if source else 0.
    leakage = float((ledger[3] + ledger[4]) / source) if source else 0.
    adequate = leakage <= domain_mass_tolerance
    return dict(schema=SCHEMA, masses=masses, diameters_m=diameter, arriving_numbers_m2=arrivals,
                seconds=gather * crossing, crossing_hours=crossing / 3600.,
                time_mean_number_m3=occupation / (gather * crossing), final_number_m3=number,
                layer_height_m=height, layer_edges_m=edges, velocity_m_s=velocity,
                mass_ledger_kg_m2=mass_ledger,
                extensive_mass_ledger_kg={k: v * horizontal_area_m2 for k, v in mass_ledger.items() if k != 'relative_residual'},
                domain_exit_numbers_m2=dict(lower=float(domain_number[0]), upper=float(domain_number[1])),
                gather_domain_exit_mass_kg_m2=dict(lower=float(gathered_domain_mass[0]), upper=float(gathered_domain_mass[1])),
                domain_adequacy=dict(passed=adequate, mass_fraction_exited=leakage,
                                     tolerance=domain_mass_tolerance, physical_size_cutoff=False,
                                     upper_absorbing_pivot_m=float(geometry.size(virtual_top)),
                                     interpretation='Finite mass-domain adequacy only; separate grid and timestep convergence required.'),
                diagnostics=dict(steps=steps, rejected_steps=rejected, minimum_step_s=smallest_dt,
                                 maximum_step_s=largest_dt, final_time_s=t, gather_start_s=gather_start,
                                 minimum_number_m3=float(number.min()),
                                 gathering_time_s=gather * crossing),
                settings=dict(top_c=top_c, number_per_m3=number_per_m3, snowfall_mm_h=snowfall_mm_h,
                              layers=layers, mass_bins=masses.size, minimum_diameter_m=float(diameter[0]),
                              maximum_diameter_m=float(diameter[-1]), max_step_s=max_step_s, cfl=cfl,
                              fill=fill, gather=gather, breakup=breakup, horizontal_area_m2=horizontal_area_m2,
                              geometry=geometry.id, domain_mass_tolerance=domain_mass_tolerance,
                              injection_reference_velocity_m_s=u_in,
                              injected_number_flux_m2_s=float(seed_flux.sum()),
                              injected_mass_flux_kg_m2_s=injected_rate,
                              injection_state='actual cloud-top height; independent of vertical layer count'))
