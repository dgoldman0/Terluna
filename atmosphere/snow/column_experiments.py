"""Explicit physical-closure comparisons using the frozen sectional numerics.

This experimental extension retains the checkpointed mass allocation, SSPRK2,
source quadrature and ledgers. Its column loop derives from sectional.py at
16cde54; that reference and its producer-bound artifacts remain unchanged.
Hydrodynamic interception and a declared asymptotic geometry sensitivity are
separate choices. Neither constitutes a calibrated giant-flake survival law.
"""
from __future__ import annotations
from dataclasses import dataclass
import math
import numpy as np
from atmosphere.snow import aggregation as ag, flakes as sf
from atmosphere.snow.habits import GEOMETRIES
from atmosphere.snow.sectional import CollisionOperator, _grid, _advance, seed_distribution, LEDGER_KEYS

SCHEMA = 'terluna.atmosphere.snow-closure-column/1'


@dataclass(frozen=True)
class QuadraticContinuation:
    """Continuous upper-size sensitivity; seed geometry is preserved.

    Join at an explicitly chosen 1 cm reference; this is not measured evolving
    morphology, an accepted production closure, or a maximum snowflake size.
    """
    base: object
    anchor_m: float = .01
    id: str = 'quadratic_above_1cm'

    @property
    def area_ratio(self): return self.base.area_ratio
    @property
    def capacitance_ratio(self): return self.base.capacitance_ratio
    def mass(self, d):
        d=np.asarray(d,dtype=float)
        return np.where(d<=self.anchor_m,self.base.mass(d),self.base.mass(self.anchor_m)*(d/self.anchor_m)**2)
    def size(self, m):
        m=np.asarray(m,dtype=float); anchor=self.base.mass(self.anchor_m)
        return np.where(m<=anchor,self.base.size(m),self.anchor_m*np.sqrt(m/anchor))
    def area(self,d): return self.base.area(d)
    def to_dict(self):
        return dict(id=self.id,anchor_m=self.anchor_m,base=self.base.to_dict(),
                    evidence='Proposed sensitivity only: preserve the baseline below 1 cm and continue m proportional to D squared above it. No evolving shape/contact state.')


class Kinetics:
    """Fall speeds and encounter kernel before sticking and chipping."""
    MODES=('legacy','outer_area','boehm_unity','boehm_sigma0','boehm_sigma025')
    def __init__(self,mode='legacy'):
        if mode not in self.MODES: raise ValueError('unknown kinetics mode')
        self.mode=mode
        self.efficiency_diagnostics=[]
    def description(self):
        return dict(mode=self.mode,alpha=1 if self.mode.startswith('boehm') else None,
                    sigma=.25 if self.mode=='boehm_sigma025' else 0.,
                    outer_area='projected-ice-area approximation' if self.mode=='legacy' else 'circumscribing-circle encounter area',
                    evidence='Explicit comparison; no giant-aggregate hydrodynamic or mechanical validation.')
    def __call__(self,m,d,geometry,gravity,state):
        m=np.asarray(m,dtype=float);d=np.asarray(d,dtype=float)
        if self.mode in ('legacy','outer_area'):
            v,re=ag.speeds(d,m,gravity,state,geometry.area_ratio)
            efficiency=geometry.area_ratio if self.mode=='legacy' else 1.
        else:
            from atmosphere.snow import boehm
            props=boehm.particle_properties(d,m,gravity,state,geometry.area_ratio,
                    alpha=1.,sigma=.25 if self.mode=='boehm_sigma025' else 0.)
            v,re=props['v'],props['Re']
            if self.mode=='boehm_unity': efficiency=1.
            else:
                left={key:np.broadcast_to(value,d.shape)[:,None] for key,value in props.items()}
                right={key:np.broadcast_to(value,d.shape)[None,:] for key,value in props.items()}
                efficiency=boehm.collision_efficiency(left,right,state)
                active=np.abs(v[:,None]-v[None,:])>0
                finite=np.isfinite(efficiency)
                self.efficiency_diagnostics.append(dict(temperature_k=float(state['t']),
                    minimum=float(np.min(efficiency[finite])) if finite.any() else None,
                    maximum=float(np.max(efficiency[finite])) if finite.any() else None,
                    negative_pairs=int(np.sum(active & (efficiency<0))),
                    nonfinite_active_pairs=int(np.sum(active & ~finite)),
                    above_one_pairs=int(np.sum(active & (efficiency>1)))))
        relative=np.abs(v[:,None]-v[None,:])
        geometric=math.pi/4*(d[:,None]+d[None,:])**2*relative
        if np.ndim(efficiency):
            kernel=np.zeros_like(geometric); active=geometric>0
            kernel[active]=geometric[active]*efficiency[active]
        else: kernel=geometric*efficiency
        if not np.isfinite(v).all() or (v<=0).any() or not np.isfinite(re).all():
            raise FloatingPointError('invalid velocity or Reynolds number from closure')
        if not np.isfinite(kernel).all() or (kernel<0).any():
            raise FloatingPointError('invalid active hydrodynamic kernel; no clipping applied')
        return v,re,kernel


def collision_operator(masses,velocity,kernel,temperature,geometry,breakup):
    relative=np.abs(velocity[:,None]-velocity[None,:])
    sticking=math.exp(.09*(temperature-sf.T0))
    fragment_mass=float(geometry.mass(ag.FRAGMENT_M)); loss=None
    if breakup:
        small=np.minimum(masses[:,None],masses[None,:])
        energy=masses[:,None]*masses[None,:]/(masses[:,None]+masses[None,:])*relative**2
        count=ag.fragments(small,energy,mass_to_size=geometry.size)
        loss=np.minimum(count*fragment_mass,.5*small)
    return CollisionOperator(masses,kernel,sticking=sticking,loss_per_breakup=loss,fragment_mass=fragment_mass)


def solve_column(world, top_c, number_per_m3, snowfall_mm_h, *, geometry=None, layers=30, mass_bins=96,
           minimum_diameter_m=2e-5, maximum_diameter_m=10., max_step_s=30., cfl=.3,
           fill=1.5, gather=1.5, breakup=True, horizontal_area_m2=1., progress=None,
           mass_grid=None, maximum_steps=2000000, domain_mass_tolerance=1e-6,
           kinetics=None, injection_velocity_m_s=None, growth_scheme='donor', batch_growth=True):
    """Experimental column with an explicit drag/encounter callback.

    The integration and ledger follow the checkpointed sectional reference.
    Keeping that implementation immutable allows direct baseline comparisons;
    this extension changes prescribed physics through Kinetics, without global
    monkeypatches or relabelling historical products. An optional fixed injection
    velocity holds the source flux and simulated duration across drag controls.
    """
    if growth_scheme not in ('donor', 'limited'):
        raise ValueError('growth_scheme must be donor or limited')
    kinetics = Kinetics() if kinetics is None else kinetics
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
    velocity, growth, operators, growth_operators = [], [], [], []
    if growth_scheme == 'limited':
        from atmosphere.snow.limited_growth import LimitedGrowth, prepare
        face_mass = prepare(masses)['face_masses']
        face_diameter = geometry.size(face_mass)
    for s in states:
        v, re, kernel = kinetics(masses, diameter, geometry, world.gravity, s)
        velocity.append(v)
        growth.append(np.maximum(ag.deposition(diameter, re, s, geometry.capacitance_ratio), 0.))
        operators.append(collision_operator(masses, v, kernel, s['t'], geometry, breakup))
        if growth_scheme == 'limited':
            _, face_re, _ = kinetics(face_mass, face_diameter, geometry, world.gravity, s)
            face_growth = np.maximum(ag.deposition(face_diameter, face_re, s, geometry.capacitance_ratio), 0.)
            growth_operators.append(LimitedGrowth(masses, face_growth, cfl=cfl))
    velocity, growth = np.asarray(velocity), np.asarray(growth)
    growth_batch = None
    if growth_scheme == 'limited' and batch_growth:
        from atmosphere.snow.limited_growth_batch import BatchLimitedGrowth
        growth_batch = BatchLimitedGrowth(growth_operators)
    seed_d = math.sqrt(ag.START_M[0] * ag.START_M[1])
    actual_seed_speed = float(kinetics(np.array([geometry.mass(seed_d)]), np.array([seed_d]),
                                      geometry, world.gravity, world.state(z_top))[0][0])
    u_in = actual_seed_speed if injection_velocity_m_s is None else float(injection_velocity_m_s)
    if not math.isfinite(u_in) or u_in <= 0:
        raise ValueError('injection reference velocity must be finite and positive')
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
        if growth_scheme == 'limited':
            deposited = 0.
            if growth_batch is not None:
                transfers, water_rates, step_limit = growth_batch.evaluate(number, production)
                for water_rate in water_rates:
                    deposited += water_rate
                limit = min(limit, step_limit)
            else:
                transfers = np.zeros_like(number)
                for k, operator in enumerate(growth_operators):
                    transfers[k], water_rate, step_limit = operator.evaluate(number[k], production)
                    deposited += water_rate
                    limit = min(limit, step_limit)
            diagnostic[1] = deposited * dz
        else:
            demand = np.sum(growth * number, axis=1)
            scale = np.minimum(1., np.divide(production, demand, out=np.ones_like(demand), where=demand > 0))
            rate = growth * scale[:, None]
            transfer_coefficient = rate / spacing
            occupied = number > 0
            if occupied.any():
                largest = float(transfer_coefficient[occupied].max())
                if largest:
                    limit = min(limit, cfl / largest)
            transfers = number * transfer_coefficient
            diagnostic[1] = float(np.sum(number * rate) * dz)
        dn -= transfers
        dn[:, 1:] += transfers[:, :-1]
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
                              injection_state=('actual cloud-top terminal speed' if injection_velocity_m_s is None else
                                               'fixed reference speed preserves source flux and integration duration'),
                              actual_seed_terminal_speed_m_s=actual_seed_speed,
                              implied_seed_number_per_m3=number_per_m3*u_in/actual_seed_speed,
                              growth_scheme=growth_scheme, batch_growth=batch_growth,
                              growth_quadrature=('limited face flux times adjacent pivot-mass spacing' if growth_scheme == 'limited'
                                                 else 'pivot growth rate times pivot number'),
                              kinetics=kinetics.description()))
