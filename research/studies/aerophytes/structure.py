"""Structure beyond the ideal arithmetic: joints, gusts, damage, creep and colonies.

Pure functions plus evaluate(). Pressures in Pa, lengths in m, masses in kg;
"per area" means per square metre of horizontal projected organism area. The
wind inputs are the committed GCM wind product and the CM1 equatorial ring.
Design guesses enter as arguments, with the ranges evaluate() reports.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_GM, MOON_RADIUS

ROOT = Path(__file__).resolve().parents[3]
INPUT_FILES = (
    'climate/results/gcm/global_winds_A28_dim5_moon.json',
    'climate/results/crm/ring_ring_equator.json',
    'research/studies/sky_ships/results/sky_ships.json',
    'research/studies/aerophytes/structure_sources.json',
)
HOURS_PER_YEAR = JULIAN_YEAR_DAYS * JULIAN_DAY / 3600
YEAR_S = JULIAN_YEAR_DAYS * JULIAN_DAY
# Madison curve of wood under sustained load (Wood 1951, FPL Report 1916):
# y = 108.4/x**0.04635 + 18.3, stress in % of the standard 5-minute bending
# strength that fails after x seconds. 10 years gives 62%.
MADISON_COEFFICIENT = 108.4
MADISON_EXPONENT = 0.04635
MADISON_ASYMPTOTE = 18.3
# Wallaby tail tendon creeps to rupture at 10 MPa and above against a dynamic
# yield near 144 MPa (Wang and Ker 1995): the collagen tendon's long-term level.
COLLAGEN_SUSTAINED_PA = 10e6
# Gust factor for sharper gusts than the CM1 6-km columns resolve, as the summit
# tower and sky-ship studies apply it (1.4 for a 3-second gust).
GUST_FACTOR = 1.4
CENTURY_FACTOR = 1.2


def _positive(**values):
    if any(not math.isfinite(v) or v <= 0 for v in values.values()):
        raise ValueError(f'finite positive inputs required: {tuple(values)}')


def _nonnegative(**values):
    if any(not math.isfinite(v) or v < 0 for v in values.values()):
        raise ValueError(f'finite nonnegative inputs required: {tuple(values)}')


def _fraction(name, value):
    if not math.isfinite(value) or not 0 < value <= 1:
        raise ValueError(f'{name} must lie in (0, 1]')


def _interp(x, xs, ys):
    """Linear interpolation inside a stored table; no extrapolation."""
    if not xs[0] <= x <= xs[-1]:
        raise ValueError('outside the stored table')
    for i in range(1, len(xs)):
        if x <= xs[i]:
            f = (x-xs[i-1])/(xs[i]-xs[i-1])
            return ys[i-1] + f*(ys[i]-ys[i-1])
    return ys[-1]


def gravity_at_height(height_m):
    _nonnegative(height=height_m)
    return MOON_GM / (MOON_RADIUS + height_m)**2


# --- joints, fittings and the long-term allowable -------------------------------

def duration_of_load(years, coefficient=MADISON_COEFFICIENT, exponent=MADISON_EXPONENT,
                     asymptote=MADISON_ASYMPTOTE):
    """Share of standard-test strength a material carries for `years` before failing.

    Wood's hyperbola y = 108.4/x^0.04635 + 18.3 (x in seconds), capped at the
    standard test. Wood is the measured biological analogue for a cellulose
    fibre composite; an aerophyte fibre's own curve is a design input.
    """
    _positive(years=years)
    seconds = years * YEAR_S
    return min(1., (coefficient/seconds**exponent + asymptote)/100)


def assembled_allowable(coupon_pa, joint_efficiency, renewal_years, design_factor):
    """Coupon strength carried through joints, sustained load and a design factor.

    A part renewed every `renewal_years` must hold its load that long; living
    replacement means the duration is the renewal interval, not the lifespan.
    """
    _positive(coupon=coupon_pa, renewal=renewal_years, design_factor=design_factor)
    _fraction('joint efficiency', joint_efficiency)
    if design_factor < 1:
        raise ValueError('design factor below one')
    duration = duration_of_load(renewal_years)
    allowable = coupon_pa * joint_efficiency * duration / design_factor
    return dict(coupon_pa=coupon_pa, joint_efficiency=joint_efficiency,
                renewal_years=renewal_years, duration_factor=duration,
                design_factor=design_factor, allowable_pa=allowable,
                knockdown=allowable/coupon_pa)


def breaking_length(allowable_pa, density_kg_m3, gravity_m_s2):
    """Length of fibre that hangs from itself at its allowable stress, sigma/(rho g)."""
    _positive(allowable=allowable_pa, density=density_kg_m3, gravity=gravity_m_s2)
    return allowable_pa / (density_kg_m3 * gravity_m_s2)


def equatorial_resultant(radius_m, base_pressure_pa, lift_density_kg_m3, gravity_m_s2,
                         extra_pressure_pa=0.):
    """Net vertical pressure force on the upper hemisphere of a full gas sphere.

    Internal excess p(z) = base + lift_density*g*(z - z_bottom) over the upper
    hemisphere's projection gives pi R² (base + extra + d g R) + 2/3 pi R³ d g.
    The tendons and film crossing the equator carry this force.
    """
    _positive(radius=radius_m, lift_density=lift_density_kg_m3, gravity=gravity_m_s2)
    _nonnegative(base=base_pressure_pa, extra=extra_pressure_pa)
    d, g, r = lift_density_kg_m3, gravity_m_s2, radius_m
    return math.pi*r*r*(base_pressure_pa + extra_pressure_pa + d*g*r) + 2/3*math.pi*r**3*d*g


def round_tendon_share(radius_m, lift_density_kg_m3, gravity_m_s2, tendon_allowable_pa,
                       tendon_density_kg_m3, base_pressure_pa=0.):
    """Share of gross lift taken by meridional tendons carrying the equatorial force.

    N tendons of length pi R carry the equatorial resultant along their whole
    length, so M = pi R rho_t F/sigma. The share of (4/3) pi R³ d is
    (3 pi/4) rho_t p_b/(d sigma) + (5 pi/4) rho_t g R/sigma: it grows in
    proportion to radius over the fibre's breaking length.
    """
    force = equatorial_resultant(radius_m, base_pressure_pa, lift_density_kg_m3, gravity_m_s2)
    _positive(allowable=tendon_allowable_pa, density=tendon_density_kg_m3)
    mass = math.pi*radius_m*tendon_density_kg_m3*force/tendon_allowable_pa
    lift = 4/3*math.pi*radius_m**3*lift_density_kg_m3
    return dict(radius_m=radius_m, equatorial_force_n=force, tendon_mass_kg=mass,
                tendon_kg_per_projected_m2=mass/(math.pi*radius_m**2),
                gross_lift_kg=lift, tendon_share_of_lift=mass/lift)


def round_ceiling_radius(share, lift_density_kg_m3, gravity_m_s2, tendon_allowable_pa,
                         tendon_density_kg_m3, base_pressure_pa=0.):
    """Radius at which ideal tendons take `share` of the gross lift (closed form)."""
    _fraction('share', share)
    _positive(lift=lift_density_kg_m3, gravity=gravity_m_s2,
              allowable=tendon_allowable_pa, density=tendon_density_kg_m3)
    constant = 3*math.pi/4*tendon_density_kg_m3*base_pressure_pa/(lift_density_kg_m3*tendon_allowable_pa)
    if constant >= share:
        return 0.
    return (share-constant)*tendon_allowable_pa/(5*math.pi/4*tendon_density_kg_m3*gravity_m_s2)


def end_fitting_share(fitting_radius_over_body_radius):
    """Mass of crown and base rings as a share of the tendon mass.

    Converging tendons put hoop force F/(2 pi) in a ring of radius r_f, so each
    ring weighs r_f F rho/sigma against the tendons' pi R F rho/sigma.
    """
    _positive(ratio=fitting_radius_over_body_radius)
    return 2*fitting_radius_over_body_radius/math.pi


def suspension_share(path_length_m, allowable_pa, density_kg_m3, gravity_m_s2):
    """Fibre mass per kg of hung payload carried over `path_length_m` to the structure."""
    return path_length_m / breaking_length(allowable_pa, density_kg_m3, gravity_m_s2)


def seam_mass_share(overlap_m, gore_width_m):
    """Extra film where gores overlap at their joins, as a share of film mass."""
    _positive(gore=gore_width_m)
    _nonnegative(overlap=overlap_m)
    return overlap_m/gore_width_m


# --- colonies of modules ---------------------------------------------------------

def module_optimum(barrier_areal_kg_m2, film_allowable_pa, film_density_kg_m3,
                   lift_density_kg_m3, gravity_m_s2, base_pressure_pa=0., floor_thickness_m=25e-6):
    """Module radius that minimises skin mass per gas volume, and that minimum.

    Skin per volume = 3 rho_s p0/(2 sigma) + 3 rho_s d g r/sigma + 3 mu_b/r, with
    mu_b the barrier alone and the stress-bearing film at least the floor thickness.
    Where the stress thickness exceeds the floor at the optimum, the minimum lies at
    r* = sqrt(mu_b sigma/(rho_s d g)). A colony of N modules of radius r carries the
    same skin share as one module, at any size.
    """
    _positive(barrier=barrier_areal_kg_m2, allowable=film_allowable_pa,
              density=film_density_kg_m3, lift=lift_density_kg_m3, gravity=gravity_m_s2)
    _nonnegative(base=base_pressure_pa, floor=floor_thickness_m)
    rs, d, g, s = film_density_kg_m3, lift_density_kg_m3, gravity_m_s2, film_allowable_pa
    r = math.sqrt(barrier_areal_kg_m2*s/(rs*d*g))
    stress_thickness = (base_pressure_pa+2*d*g*r)*r/(2*s)
    per_volume = 3*rs*base_pressure_pa/(2*s) + 6*math.sqrt(barrier_areal_kg_m2*rs*d*g/s)
    return dict(optimum_radius_m=r, optimum_diameter_m=2*r,
                skin_kg_per_m3=per_volume, skin_share_of_lift=per_volume/d,
                stress_thickness_at_optimum_um=stress_thickness*1e6,
                floor_binds=stress_thickness < floor_thickness_m)


def skin_per_volume(radius_m, barrier_areal_kg_m2, film_allowable_pa, film_density_kg_m3,
                    lift_density_kg_m3, gravity_m_s2, base_pressure_pa=0.):
    """Skin mass per gas volume of one spherical module of the given radius."""
    _positive(radius=radius_m)
    rs, d, g, s = film_density_kg_m3, lift_density_kg_m3, gravity_m_s2, film_allowable_pa
    return 3*rs*base_pressure_pa/(2*s) + 3*rs*d*g*radius_m/s + 3*barrier_areal_kg_m2/radius_m


def hex_raft(module_radius_m, lift_density_kg_m3):
    """A single layer of touching spheres on a hexagonal grid, per projected area."""
    _positive(radius=module_radius_m, lift=lift_density_kg_m3)
    cell = 2*math.sqrt(3)*module_radius_m**2
    volume = 4/3*math.pi*module_radius_m**3
    return dict(cell_area_m2=cell, cell_over_disc=cell/(math.pi*module_radius_m**2),
                gas_volume_per_projected_m2=volume/cell,
                gross_lift_kg_per_projected_m2=lift_density_kg_m3*volume/cell,
                skin_area_per_projected_m2=4*math.pi*module_radius_m**2/cell)


def colony_network(edge_pressure_pa, layer_thickness_m, allowable_pa, density_kg_m3):
    """Tie mass per projected area that holds a raft together against an edge gust.

    The windward edge meets q over the raft's thickness h, and the ties pass
    that force back through the raft: rho q h/sigma per projected area,
    whatever the raft's width.
    """
    _nonnegative(pressure=edge_pressure_pa)
    _positive(thickness=layer_thickness_m, allowable=allowable_pa, density=density_kg_m3)
    return density_kg_m3*edge_pressure_pa*layer_thickness_m/allowable_pa


def quilt_structure(thickness_m, cell_width_m, design_pressure_pa, film_allowable_pa,
                    tie_allowable_pa, density_kg_m3, barrier_areal_kg_m2,
                    minimum_film_thickness_m=0.):
    """A flat slab of gas cells under a shared skin, per projected area.

    Top and bottom faces bulge between cell walls (hoop t = p w/(2 sigma) on a
    cylinder of radius about w/2); the walls tie the faces together, each wall
    carrying p w per metre of its length. The rim is left out.
    """
    _positive(thickness=thickness_m, cell=cell_width_m, film=film_allowable_pa,
              tie=tie_allowable_pa, density=density_kg_m3)
    _nonnegative(pressure=design_pressure_pa, barrier=barrier_areal_kg_m2,
                 minimum=minimum_film_thickness_m)
    face_t = max(minimum_film_thickness_m, design_pressure_pa*cell_width_m/(2*film_allowable_pa))
    faces = 2*(face_t*density_kg_m3 + barrier_areal_kg_m2)
    wall_area = 2*thickness_m/cell_width_m
    ties = wall_area*design_pressure_pa*cell_width_m/tie_allowable_pa*density_kg_m3
    wall_floor = wall_area*minimum_film_thickness_m*density_kg_m3
    return dict(face_thickness_m=face_t, faces_kg_m2=faces, wall_area_per_projected_m2=wall_area,
                ties_kg_m2=ties, wall_film_kg_m2=wall_floor,
                total_kg_m2=faces+ties+wall_floor,
                gas_contact_area_per_projected_m2=2.)


# --- gusts and the body's response -----------------------------------------------

def gust_pressure(speed_m_s, air_density_kg_m3):
    _nonnegative(speed=speed_m_s)
    _positive(density=air_density_kg_m3)
    return .5*air_density_kg_m3*speed_m_s**2


def relaxation_time(radius_m, drag_coefficient, relative_speed_m_s, added_mass=.5):
    """Time for a neutrally buoyant sphere to halve its wind after a step change.

    m_eff du/dt = -(1/2) rho Cd pi R² u² with m_eff = (1 + Ca) rho (4/3) pi R³
    gives u = u0/(1 + t/tau), tau = 8 (1 + Ca) R/(3 Cd u0). A body much larger
    than about 10 m meets a gust front's full jump for tens of seconds or more.
    """
    _positive(radius=radius_m, cd=drag_coefficient, speed=relative_speed_m_s)
    _nonnegative(added_mass=added_mass)
    return 8*(1+added_mass)*radius_m/(3*drag_coefficient*relative_speed_m_s)


def design_gusts(global_winds, ring, heights_km=(1, 5, 10, 20, 30, 40)):
    """Gusts relative to a drifting body, by height, from the committed winds.

    Horizontal: the jump from the ring's median to its highest 3-hourly wind.
    Vertical: the ring's strongest updraft (a body holding its height meets it).
    Both times the 1.4 gust factor; a centuries case adds the 1.2 decades factor.
    The GCM 50-year level is the large-scale wind a drifting body rides.
    """
    rows = []
    gh = global_winds['height_km']
    for h in heights_km:
        key = f'{h}_km'
        if key not in ring['wind_aloft_m_s'] or key not in ring['vertical_wind_m_s']:
            raise ValueError(f'ring lacks height {h} km')
        aloft = ring['wind_aloft_m_s'][key]
        vertical = ring['vertical_wind_m_s'][key]
        density = _interp(h, gh, global_winds['density_mean_kg_m3'])
        jump = aloft['max'] - aloft['median']
        horizontal = GUST_FACTOR*jump
        upward = GUST_FACTOR*vertical['up_max']
        design = max(horizontal, upward)
        rows.append(dict(height_km=h, air_density_kg_m3=density,
            ring_median_m_s=aloft['median'], ring_p99_m_s=aloft['p99'], ring_max_m_s=aloft['max'],
            ring_updraft_max_m_s=vertical['up_max'], ring_vertical_abs_p99_m_s=vertical['abs_p99'],
            gcm_50yr_large_scale_m_s=_interp(h, gh, global_winds['return_50yr_m_s']),
            horizontal_jump_m_s=jump, design_horizontal_m_s=horizontal,
            design_vertical_m_s=upward, design_gust_m_s=design,
            design_pressure_pa=gust_pressure(design, density),
            century_gust_m_s=CENTURY_FACTOR*design,
            century_pressure_pa=gust_pressure(CENTURY_FACTOR*design, density)))
    return rows


def gust_stress_ratio(radius_m, lift_density_kg_m3, gravity_m_s2, gust_pressure_pa,
                      base_pressure_pa=20.):
    """Share of the design membrane stress that a storm gust adds, by size.

    Small bodies are gust-dominated and meet many load cycles; giants are
    dominated by their steady hydrostatic head, which creep governs.
    """
    head = 2*radius_m*lift_density_kg_m3*gravity_m_s2
    total = base_pressure_pa + head + gust_pressure_pa
    return dict(radius_m=radius_m, head_pa=head, gust_pa=gust_pressure_pa,
                design_pa=total, gust_share=gust_pressure_pa/total,
                steady_share=(base_pressure_pa+head)/total)


# --- damage, tear arrest and compartments ----------------------------------------

def tear_outflow(tear_length_m, opening_ratio, differential_pressure_pa, gas_density_kg_m3,
                 discharge_coefficient=.6):
    """Initial gas flow from a slit that opens to area beta L² under tension."""
    _nonnegative(length=tear_length_m, pressure=differential_pressure_pa)
    _positive(gas=gas_density_kg_m3)
    _fraction('opening ratio', opening_ratio)
    _fraction('discharge coefficient', discharge_coefficient)
    area = opening_ratio*tear_length_m**2
    speed = math.sqrt(2*differential_pressure_pa/gas_density_kg_m3)
    return dict(tear_length_m=tear_length_m, open_area_m2=area, jet_speed_m_s=speed,
                volume_m3_s=discharge_coefficient*area*speed)


def compartment_loss(supported_kg_m2, dumpable_water_kg_m2, compartments):
    """Whether dropping free water can replace the lift of one lost compartment."""
    _positive(supported=supported_kg_m2)
    _nonnegative(dumpable=dumpable_water_kg_m2)
    if compartments < 1 or int(compartments) != compartments:
        raise ValueError('whole number of compartments required')
    lost = supported_kg_m2/compartments
    return dict(compartments=compartments, lift_lost_kg_m2=lost,
                dumpable_water_kg_m2=dumpable_water_kg_m2,
                stays_aloft=dumpable_water_kg_m2 >= lost,
                minimum_compartments=math.ceil(supported_kg_m2/dumpable_water_kg_m2)
                if dumpable_water_kg_m2 > 0 else None)


def partition_walls(radius_m, compartments, wall_areal_kg_m2):
    """Orange-segment walls: N half-discs of area pi R²/2, per projected area."""
    _positive(radius=radius_m)
    _nonnegative(areal=wall_areal_kg_m2)
    if compartments < 1:
        raise ValueError('at least one compartment')
    walls = 0 if compartments == 1 else compartments
    return walls*0.5*wall_areal_kg_m2


# --- creep, fatigue and renewal ---------------------------------------------------

def renewal_burden(structure_dry_kg_m2, renewal_years, construction_c_per_kg_dry=1.39*.4):
    """Assimilate carbon per projected m² a year to renew a structure every N years."""
    _nonnegative(structure=structure_dry_kg_m2)
    _positive(renewal=renewal_years, construction=construction_c_per_kg_dry)
    return structure_dry_kg_m2*construction_c_per_kg_dry/renewal_years


def gust_cycles(storm_hours_per_year, gust_period_s, years):
    """Load cycles a body meets in storms, a count for creep-fatigue bookkeeping."""
    _nonnegative(hours=storm_hours_per_year, years=years)
    _positive(period=gust_period_s)
    per_year = storm_hours_per_year*3600/gust_period_s
    return dict(cycles_per_year=per_year, cycles=per_year*years)


def evaluate(root=ROOT, round_spare_pa=300., module_spare_pa=100.):
    """Design gusts, joint chains, ceilings, module optimum and damage cases.

    round_spare_pa and module_spare_pa are the spare superpressures the coupled round
    hulls and colony modules carry for trim (run.py passes its design guesses); the
    ceilings, module optimum and gust shares include them with the 20 Pa base.
    """
    root = Path(root)
    winds = json.loads((root/INPUT_FILES[0]).read_text())
    ring = json.loads((root/INPUT_FILES[1]).read_text())
    ships = json.loads((root/INPUT_FILES[2]).read_text())
    for data, schema in ((winds, 'terluna.climate.gcm-global-winds/1'),
                         (ring, 'terluna.climate.crm-ring/1'),
                         (ships, 'terluna.research.sky-ships/1')):
        if data.get('schema') != schema:
            raise ValueError(f'unexpected input schema {data.get("schema")}')
    air = {a['height_km']: a for a in ships['air']}
    ten = air[10]
    g = gravity_at_height(10000.)
    lift = ten['lift_hydrogen_kg_m3']
    gusts = design_gusts(winds, ring)
    gust10 = next(r for r in gusts if r['height_km'] == 10)
    film_cases = [assembled_allowable(c*1e6, j, y, f)
                  for c in (14., 20., 33.) for j in (.5, .7, .9) for y in (1., 10., 30.) for f in (1.5,)]
    reference_film = assembled_allowable(20e6, .7, 10., 1.5)
    tendon_cases = [assembled_allowable(c*1e6, j, y, 1.5)
                    for c in (100., 300., 1000.) for j in (.6, .8, .9) for y in (10., 30., 100., 300.)]
    collagen = dict(sustained_pa=COLLAGEN_SUSTAINED_PA, joint_efficiency=.8, design_factor=1.5,
                    allowable_pa=COLLAGEN_SUSTAINED_PA*.8/1.5,
                    note='Collagen tendon held below the stress at which wallaby tail tendon creeps to rupture.')
    ceilings = []
    fibres = [(f'cellulose_{c:g}_mpa_{y:g}_years', assembled_allowable(c*1e6, .8, y, 1.5)['allowable_pa'], c, y)
              for c in (100., 300., 1000.) for y in (30., 100., 300.)]
    fibres.append(('collagen_sustained', COLLAGEN_SUSTAINED_PA*.8/1.5, None, None))
    for name, allowable, c, y in fibres:
        row = dict(fibre=name, tendon_strength_mpa=c, renewal_years=y, allowable_mpa=allowable/1e6,
                   breaking_length_km=breaking_length(allowable, 1500., g)/1000)
        for share in (.1, .25, .5):
            row[f'radius_at_share_{share}_m'] = round_ceiling_radius(
                share, lift, g, allowable, 1500., 20.+round_spare_pa+gust10['design_pressure_pa'])
        row['base_term_share'] = 3*math.pi/4*1500.*(20.+round_spare_pa+gust10['design_pressure_pa'])/(lift*allowable)
        ceilings.append(row)
    barrier = 100e-6*1500  # the barrier alone; the stress film is at least its 25 um floor
    module = [dict(film_allowable_mpa=s/1e6, barrier_areal_kg_m2=b, base_pressure_pa=p,
                   **module_optimum(b, s, 1500., lift, g, p))
              for s in (reference_film['allowable_pa'], 5e6, 2.5e6)
              for b in (barrier, 25e-6*1500)
              for p in (20.+module_spare_pa, 20.+module_spare_pa+gust10['design_pressure_pa'])]
    stress = [gust_stress_ratio(r, lift, g, gust10['design_pressure_pa'], 20.+round_spare_pa)
              for r in (5., 10., 20., 50., 100., 250., 500., 1000.)]
    tears = []
    gas_density = ten['density_kg_m3'] - lift
    for radius in (20., 100., 1000.):
        crown = 20. + 2*radius*lift*g
        for length in (.1, 1., 4., 10.):
            for beta in (.05, .2):
                tears.append(dict(radius_m=radius, crown_pressure_pa=crown, opening_ratio=beta,
                                  **tear_outflow(length, beta, crown, gas_density)))
    compartments = [dict(supported_kg_m2=s, **compartment_loss(s, w, n))
                    for s in (18.7, 25., 40.) for w in (2.5, 5.) for n in (1, 4, 8, 16)]
    relax = [dict(radius_m=r, drag_coefficient=cd, relative_speed_m_s=gust10['design_gust_m_s'],
                  relaxation_s=relaxation_time(r, cd, gust10['design_gust_m_s']))
             for r in (1., 10., 20., 100., 1000.) for cd in (.47, .05)]
    return dict(
        evidence='Closed-form structural requirements on committed winds (GCM 3-day means, CM1 equatorial ring) and sky-ship air. Joint efficiencies, renewal intervals, fibre strengths and opening ratios are design guesses given as ranges; the duration-of-load curve is measured on wood and the collagen level on wallaby tendon.',
        reading_rule='Per area means per projected m² of organism. Gusts are relative to a drifting body: the horizontal jump from the ring median to its highest wind, and the strongest updraft for a body that holds its height, each times 1.4. Allowable = coupon x joint efficiency x Madison factor at the renewal interval / design factor.',
        design_gusts=gusts,
        film_allowables=film_cases, reference_film=reference_film,
        tendon_allowables=tendon_cases, collagen_tendon=collagen, round_tendon_ceilings=ceilings,
        end_fitting_share=[dict(fitting_over_radius=x, share_of_tendon_mass=end_fitting_share(x)) for x in (.02, .05, .1)],
        suspension_share=[dict(path_m=L, allowable_mpa=s/1e6,
                               fibre_kg_per_kg_payload=suspension_share(L, s, 1500., g))
                          for L in (10., 100., 1000.) for s in (28e6, 100e6)],
        seam_mass_share=[dict(overlap_m=o, gore_width_m=w, share=seam_mass_share(o, w)) for o in (.02, .05) for w in (1., 4.)],
        module_optimum=module,
        raft=[dict(module_radius_m=r, **hex_raft(r, lift)) for r in (10., 20., 30., 50.)],
        colony_network=[dict(edge_pressure_pa=q, layer_thickness_m=h, allowable_mpa=s/1e6,
                             tie_kg_per_projected_m2=colony_network(q, h, s, 1500.))
                        for q in (gust10['design_pressure_pa'], gust10['century_pressure_pa'])
                        for h in (40., 60., 100.) for s in (28e6, 85e6)],
        quilt=[dict(thickness_m=h, cell_width_m=w,
                    **quilt_structure(h, w, 20.+h*lift*g+gust10['design_pressure_pa'],
                                      reference_film['allowable_pa'], 28e6, 1500., .15, 25e-6))
               for h in (30., 40., 60.) for w in (20., 40.)],
        gust_stress=stress, relaxation=relax, tears=tears, compartments=compartments,
        gust_cycles=[dict(storm_hours_per_year=h, gust_period_s=p, years=y, **gust_cycles(h, p, y))
                     for h in (50., 500.) for p in (10., 100.) for y in (10., 300.)],
        unresolved=['Assembled biaxial wet strength, joint efficiency and creep-rupture of a candidate fibre and film.',
                    'Gusts sharper than the CM1 6-km columns; a centuries return level from longer storm records.',
                    'Fluid-structure response of soft lobes, quilts and rafts in gust fronts and updrafts.',
                    'Tear growth in a living membrane, healing rate under pressure and ripstop spacing.'])
