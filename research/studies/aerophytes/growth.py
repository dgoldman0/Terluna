"""Growth over centuries: gas costs, ages by architecture, the canopy nursery and phosphorus.

Pure functions plus evaluate(). Carbon is kg assimilate C; "per area" is per
projected m² of organism. The ages integrate a body's construction cost
against its own surplus; run.py supplies those tables from its coupled cases.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, SYNODIC_MONTH_DAYS

ROOT = Path(__file__).resolve().parents[3]
INPUT_FILES = (
    'biosphere/canopy/results/plant.json',
    'biosphere/canopy/results/fruit.json',
    'research/studies/megaforest_wind/results/reference_cases.csv',
    'climate/results/crm/ring_ring_equator.json',
    'research/studies/aerophytes/growth_sources.json',
)
YEAR_S = JULIAN_DAY * JULIAN_YEAR_DAYS
CYCLES_PER_YEAR = JULIAN_YEAR_DAYS / SYNODIC_MONTH_DAYS
H2_LHV_J_KG = 120e6
H2_MOLAR_KG = 0.002016
GLUCOSE_MOLAR_KG = 0.180156
GLUCOSE_C_FRACTION = 0.4
CONSTRUCTION_C_PER_KG_DRY = 1.39 * GLUCOSE_C_FRACTION
# Megaforest study allometry (research/studies/megaforest_wind run.py scenario()).
CROWN_RADIUS_PER_HEIGHT = 0.18
CROWN_BASE_FRACTION = 0.45
CROWN_FRONTAL_FRACTION = 0.30
CROWN_DRAG_COEFFICIENT = 0.8


def _positive(**values):
    if any(not math.isfinite(v) or v <= 0 for v in values.values()):
        raise ValueError(f'finite positive inputs required: {tuple(values)}')


def _nonnegative(**values):
    if any(not math.isfinite(v) or v < 0 for v in values.values()):
        raise ValueError(f'finite nonnegative inputs required: {tuple(values)}')


def _fraction(name, value):
    if not math.isfinite(value) or not 0 < value <= 1:
        raise ValueError(f'{name} must lie in (0, 1]')


# --- what hydrogen costs in carbon ------------------------------------------------

def glucose_per_hydrogen(mol_h2_per_mol_glucose):
    """kg glucose fermented per kg H2 at a molar yield (4 is the acetate ceiling)."""
    _positive(yield_=mol_h2_per_mol_glucose)
    return GLUCOSE_MOLAR_KG/(H2_MOLAR_KG*mol_h2_per_mol_glucose)


def fermentation_carbon_per_hydrogen(mol_h2_per_mol_glucose):
    """kg sugar carbon committed per kg H2 by dark fermentation."""
    return glucose_per_hydrogen(mol_h2_per_mol_glucose)*GLUCOSE_C_FRACTION


def photo_carbon_per_hydrogen(gpp_c_kg_m2_year, incident_w_m2, efficiency):
    """Carbon the same photons would have fixed, per kg H2 made by photolysis.

    Area-time given to hydrogen at efficiency eta of incident I fixes no carbon,
    so 1 kg H2 displaces GPP*LHV/(eta*I*year) kg C. Higher GPP raises this
    opportunity cost and the surplus together.
    """
    _positive(gpp=gpp_c_kg_m2_year, incident=incident_w_m2)
    _fraction('efficiency', efficiency)
    area_years = H2_LHV_J_KG/(efficiency*incident_w_m2*YEAR_S)
    return gpp_c_kg_m2_year*area_years


def construction_per_area(hydrogen_kg_m2, structure_dry_kg_m2, tissue_dry_kg_m2,
                          reserve_c_kg_m2, carbon_per_hydrogen):
    """Assimilate carbon to build one projected m² of body, by part."""
    _nonnegative(h2=hydrogen_kg_m2, structure=structure_dry_kg_m2, tissue=tissue_dry_kg_m2,
                 reserve=reserve_c_kg_m2, cost=carbon_per_hydrogen)
    gas = hydrogen_kg_m2*carbon_per_hydrogen
    structure = structure_dry_kg_m2*CONSTRUCTION_C_PER_KG_DRY
    tissue = tissue_dry_kg_m2*CONSTRUCTION_C_PER_KG_DRY
    return dict(gas_c_kg_m2=gas, structure_c_kg_m2=structure, tissue_c_kg_m2=tissue,
                reserve_c_kg_m2=reserve_c_kg_m2,
                total_c_kg_m2=gas+structure+tissue+reserve_c_kg_m2,
                gas_share=gas/(gas+structure+tissue+reserve_c_kg_m2) if gas+structure+tissue+reserve_c_kg_m2 else 0.)


# --- ages --------------------------------------------------------------------------

def round_body_ages(rows, start_radius_m, targets_m, growth_share=1.):
    """Years for a round body to grow from start radius to each target radius.

    rows: increasing radius with construction_c_kg_m2 (whole body built at that
    radius, per projected m²) and surplus_c_kg_m2_year. Total construction is
    C(R) pi R²; growth spends share*S(R) pi R² a year on dC/dR. Per-area values
    are interpolated linearly in radius and the time integrated on a fine grid.
    Growth stops where the surplus first falls to zero; later targets return None.
    """
    _fraction('growth share', growth_share)
    rows = sorted(rows, key=lambda r: r['radius_m'])
    if len(rows) < 2 or start_radius_m < rows[0]['radius_m']:
        raise ValueError('rows must start at or below the start radius')
    radii = [r['radius_m'] for r in rows]
    cost = [r['construction_c_kg_m2'] for r in rows]
    surplus = [r['surplus_c_kg_m2_year'] for r in rows]

    def at(x, ys):
        for i in range(1, len(radii)):
            if x <= radii[i]:
                f = (x-radii[i-1])/(radii[i]-radii[i-1])
                return ys[i-1] + f*(ys[i]-ys[i-1])
        return ys[-1]

    limit = None
    for i in range(1, len(radii)):
        if radii[i] > start_radius_m and surplus[i] <= 0:
            lo = max(radii[i-1], start_radius_m)
            s_lo = at(lo, surplus)
            limit = lo if s_lo <= 0 else lo + (radii[i]-lo)*s_lo/(s_lo-surplus[i])
            break
    end = min(limit if limit is not None else radii[-1], radii[-1])
    years, age = {}, 0.
    steps = 4000
    grid = [start_radius_m + (end-start_radius_m)*k/steps for k in range(steps+1)]
    total = [at(r, cost)*math.pi*r*r for r in grid]
    rate = [max(at(r, surplus), 0.)*math.pi*r*r*growth_share for r in grid]
    for k in range(1, len(grid)):
        mean = .5*(rate[k-1]+rate[k])
        if mean <= 0:
            break
        step = (total[k]-total[k-1])/mean
        for target in targets_m:
            if target not in years and grid[k-1] < target <= grid[k]:
                years[target] = age + step*(target-grid[k-1])/(grid[k]-grid[k-1])
        age += step
    return dict(start_radius_m=start_radius_m, growth_share=growth_share,
                years_to_radius={str(t): years.get(t) for t in targets_m},
                growth_stops_at_radius_m=limit)


def flat_body_ages(construction_c_kg_m2, surplus_c_kg_m2_year, start_area_m2, target_areas_m2,
                   growth_share=1.):
    """A flat colony whose every new m² costs the same and earns the same grows exponentially.

    dA/dt = share*S/C * A, so A doubles every C ln2/(share S) years.
    """
    _positive(construction=construction_c_kg_m2, start=start_area_m2)
    _fraction('growth share', growth_share)
    if surplus_c_kg_m2_year <= 0:
        return dict(doubling_years=None, years_to_area={str(a): None for a in target_areas_m2})
    rate = growth_share*surplus_c_kg_m2_year/construction_c_kg_m2
    return dict(doubling_years=math.log(2)/rate, relative_growth_per_year=rate,
                years_to_area={str(a): (math.log(a/start_area_m2)/rate if a >= start_area_m2 else 0.)
                               for a in target_areas_m2})


def survival(years, annual_hazard):
    """Probability of surviving `years` at a constant annual catastrophic hazard."""
    _nonnegative(years=years, hazard=annual_hazard)
    return math.exp(-annual_hazard*years)


def recruitment_need(adult_hazard_per_year, juvenile_survival):
    """Juveniles each adult must release a year to replace adult deaths."""
    _nonnegative(hazard=adult_hazard_per_year)
    _fraction('juvenile survival', juvenile_survival)
    return adult_hazard_per_year/juvenile_survival


# --- the canopy nursery -----------------------------------------------------------

def host_crown(tree_height_m):
    """Crown footprint and drag area of the megaforest study's allometry."""
    _positive(height=tree_height_m)
    radius = CROWN_RADIUS_PER_HEIGHT*tree_height_m
    depth = (1-CROWN_BASE_FRACTION)*tree_height_m
    # Elliptic silhouette sqrt(4s(1-s)) over the crown depth integrates to pi/4.
    silhouette = 2*radius*CROWN_FRONTAL_FRACTION*depth*math.pi/4
    centroid = (CROWN_BASE_FRACTION + (1-CROWN_BASE_FRACTION)/2)*tree_height_m
    return dict(tree_height_m=tree_height_m, crown_radius_m=radius,
                footprint_m2=math.pi*radius**2, footprint_ha=math.pi*radius**2/1e4,
                crown_drag_area_m2=CROWN_DRAG_COEFFICIENT*silhouette,
                crown_centroid_m=centroid)


def host_supply(tree_height_m, npp_c_kg_m2_year, diverted_share):
    """Carbon a host can pass to juveniles a year: footprint NPP times a share."""
    _nonnegative(npp=npp_c_kg_m2_year)
    _fraction('diverted share', diverted_share)
    crown = host_crown(tree_height_m)
    npp = crown['footprint_m2']*npp_c_kg_m2_year
    return dict(**crown, tree_npp_c_kg_year=npp, diverted_c_kg_year=npp*diverted_share,
                diverted_share=diverted_share)


def nursery_fill_years(hydrogen_kg, mol_h2_per_mol_glucose, diverted_c_kg_year,
                       own_hydrogen_kg_year=0., loss_kg_year=0.):
    """Years for host sugar (fermented) plus the juvenile's own photolysis to fill it.

    The juvenile loses gas through its skin while it fills: loss_kg_year at the full
    inventory, charged in full, which bounds the fill time from above.
    """
    _nonnegative(hydrogen=hydrogen_kg, diverted=diverted_c_kg_year, own=own_hydrogen_kg_year,
                 loss=loss_kg_year)
    from_host = diverted_c_kg_year/fermentation_carbon_per_hydrogen(mol_h2_per_mol_glucose)
    rate = from_host + own_hydrogen_kg_year - loss_kg_year
    return dict(hydrogen_kg=hydrogen_kg, hydrogen_from_host_kg_year=from_host,
                own_hydrogen_kg_year=own_hydrogen_kg_year, loss_kg_year=loss_kg_year,
                years=hydrogen_kg/rate if rate > 0 else None)


def host_wind_moment_ratio(juvenile_radius_m, tree_height_m, juvenile_drag_coefficient=.47,
                           height_on_tree=1.):
    """Juvenile's overturning moment over the host crown's at the same wind.

    The juvenile rides at `height_on_tree` x tree height; the crown's drag acts
    at its centroid. The host's first-failure wind falls by 1/sqrt(1 + ratio)
    when the juvenile adds to an otherwise unchanged quadratic drag.
    """
    _positive(radius=juvenile_radius_m, height=tree_height_m, cd=juvenile_drag_coefficient)
    _fraction('height on tree', height_on_tree)
    crown = host_crown(tree_height_m)
    juvenile = juvenile_drag_coefficient*math.pi*juvenile_radius_m**2*height_on_tree*tree_height_m
    tree = crown['crown_drag_area_m2']*crown['crown_centroid_m']
    ratio = juvenile/tree
    return dict(juvenile_radius_m=juvenile_radius_m, tree_height_m=tree_height_m,
                moment_ratio=ratio, threshold_factor=1/math.sqrt(1+ratio))


# --- phosphorus -------------------------------------------------------------------

def phosphorus_stock(active_dry_kg_m2, structure_dry_kg_m2, active_p_g_kg, structure_p_g_kg):
    """Phosphorus held per projected m² in active tissue and in inert structure, g."""
    _nonnegative(active=active_dry_kg_m2, structure=structure_dry_kg_m2,
                 active_p=active_p_g_kg, structure_p=structure_p_g_kg)
    active = active_dry_kg_m2*active_p_g_kg
    structure = structure_dry_kg_m2*structure_p_g_kg
    return dict(active_g_m2=active, structure_g_m2=structure, total_g_m2=active+structure)


def phosphorus_net_need(active_dry_kg_m2, active_turnover_per_year, active_p_g_kg,
                        resorption, structure_dry_kg_m2, structure_turnover_per_year,
                        structure_p_g_kg, growth_and_offspring_g_m2_year=0.):
    """Net P a body must import a year per projected m² after internal recycling."""
    _nonnegative(active=active_dry_kg_m2, turnover=active_turnover_per_year,
                 active_p=active_p_g_kg, structure=structure_dry_kg_m2,
                 structure_turnover=structure_turnover_per_year, structure_p=structure_p_g_kg,
                 growth=growth_and_offspring_g_m2_year)
    if not 0 <= resorption < 1:
        raise ValueError('resorption in [0, 1)')
    tissue = active_dry_kg_m2*active_turnover_per_year*active_p_g_kg*(1-resorption)
    structure = structure_dry_kg_m2*structure_turnover_per_year*structure_p_g_kg
    return dict(tissue_loss_g_m2_year=tissue, structure_loss_g_m2_year=structure,
                growth_and_offspring_g_m2_year=growth_and_offspring_g_m2_year,
                net_need_g_m2_year=tissue+structure+growth_and_offspring_g_m2_year)


def rain_phosphorus(rain_mm_day, p_mg_per_litre, catch_share=1.):
    """P delivered by rain on the body, g per projected m² a year."""
    _nonnegative(rain=rain_mm_day, p=p_mg_per_litre)
    _fraction('catch share', catch_share)
    return rain_mm_day*JULIAN_YEAR_DAYS*catch_share*p_mg_per_litre/1000


def capture_phosphorus(particle_ug_m3, p_mass_fraction, relative_speed_m_s, efficiency,
                       duty=1.):
    """P captured from airborne particles swept through collectors, g/m² a year."""
    _nonnegative(conc=particle_ug_m3, p=p_mass_fraction, speed=relative_speed_m_s)
    _fraction('efficiency', efficiency)
    if not 0 <= duty <= 1:
        raise ValueError('duty in [0, 1]')
    return particle_ug_m3*1e-6*p_mass_fraction*relative_speed_m_s*efficiency*duty*YEAR_S


def _plant(root):
    data = json.loads((Path(root)/INPUT_FILES[0]).read_text())
    if data.get('schema') != 'terluna.biosphere.plant-carbon-cycle/1':
        raise ValueError('unexpected plant product')
    return data


def _reference_tree(root):
    """The megaforest study's 300 m rigid reference tree in uniform wind (30 mm wet crown)."""
    with open(Path(root)/INPUT_FILES[2], newline='') as handle:
        for row in csv.DictReader(handle):
            if row['case_id'] == '300m-rigid-uniform-reference-30mm':
                return dict(case_id=row['case_id'], height_m=float(row['height_m']),
                            first_threshold_m_s=float(row['first_threshold_m_s']),
                            first_mode=row['first_mode'])
    raise ValueError('reference tree case missing')


def evaluate(root=ROOT):
    """Host canopy supply, nursery fill times, wind moments and phosphorus guesses."""
    root = Path(root)
    plant = _plant(root)
    fruit = json.loads((root/INPUT_FILES[1]).read_text())
    near = plant['results']['near']
    canopy = {}
    for band in ('0', '30', '60'):
        p = near[band]['earth_like']
        canopy[band] = dict(gpp_g_c_m2_cycle=p['gpp'], npp_g_c_m2_cycle=p['growth'],
                            carbon_use_efficiency=p['carbon_use_efficiency'],
                            gpp_c_kg_m2_year=p['gpp']*CYCLES_PER_YEAR/1000,
                            npp_c_kg_m2_year=p['growth']*CYCLES_PER_YEAR/1000)
    eq = canopy['0']
    # Rain supports 26-58% of the crop model's growth within 15 degrees (ecology people_and_land.md, finding 10).
    npp_cases = dict(clear_sky_equator=eq['npp_c_kg_m2_year'],
                     rain_limited_low=.26*eq['npp_c_kg_m2_year'],
                     rain_limited_high=.58*eq['npp_c_kg_m2_year'],
                     clear_sky_60=canopy['60']['npp_c_kg_m2_year'])
    yields = (4., 2.32, 2.1, 1.63)
    sugar = [dict(mol_h2_per_mol_glucose=y, glucose_kg_per_kg_h2=glucose_per_hydrogen(y),
                  carbon_kg_per_kg_h2=fermentation_carbon_per_hydrogen(y)) for y in yields]
    hosts = [dict(npp_case=k, **host_supply(h, v, s))
             for h in (100., 300., 500.) for k, v in npp_cases.items() for s in (.1, .3, .5)]
    moments = [dict(**host_wind_moment_ratio(r, h)) for h in (100., 300., 500.) for r in (5., 10., 20., 30., 50.)]
    reference = _reference_tree(root)
    reference_threshold = reference['first_threshold_m_s']
    for row in moments:
        if row['tree_height_m'] == 300.:
            row['reference_300m_threshold_m_s'] = reference_threshold*row['threshold_factor']
    stock = [dict(active_p_g_kg=a, structure_p_g_kg=s, structure_dry_kg_m2=d,
                  **phosphorus_stock(1.1, d, a, s))
             for a in (1., 2., 3.) for s in (.05, .3) for d in (2., 4.)]
    need = [dict(resorption=r, active_p_g_kg=a, structure_p_g_kg=s,
                 **phosphorus_net_need(1.1, .5, a, r, 3., .1, s))
            for r in (.5, .65, .8) for a in (1., 2., 3.) for s in (.05, .3)]
    rain = [dict(rain_mm_day=r, p_mg_per_litre=c, p_g_m2_year=rain_phosphorus(r, c))
            for r in (.1, 1.25, 2.24, 5.5) for c in (.005, .02, .05)]
    capture = [dict(particle_ug_m3=c, p_mass_fraction=f, relative_speed_m_s=v, efficiency=e,
                    p_g_m2_year=capture_phosphorus(c, f, v, e))
               for c in (1., 10.) for f in (.001, .005) for v in (.1, 1.) for e in (.1, .5)]
    ring = json.loads((root/INPUT_FILES[3]).read_text())
    land = ring['by_hour_angle']['land']
    dawn = [dict(hour_angle_deg=h, wind_10m_m_s=w, mixed_layer_km=m, cloud_base_km=c, rain_mm_h=r)
            for h, w, m, c, r in zip(ring['by_hour_angle']['hour_angle_deg'], land['wind_10m_m_s'],
                                     land['mixed_layer_km'], land['cloud_base_lcl_km'], land['rain_mm_h'])
            if -115 <= h <= 5]
    return dict(
        evidence='Host production from the committed canopy carbon cycle (clear sky, no water limit) and the megaforest study allometry; fermentation yields from the measured reactor range; phosphorus contents, resorption, rain P and particle loads are design guesses given as ranges.',
        reading_rule='Host carbon per tree a year; nursery fill in years. Phosphorus in g per projected m² of aerophyte. Ages are computed in run.py from the coupled cases with these functions.',
        cycles_per_year=CYCLES_PER_YEAR, canopy=canopy, host_npp_cases=npp_cases,
        fruit_share_precedent=fruit['settings']['shares'], reference_tree=reference,
        fermentation=sugar, hosts=hosts, host_wind=moments,
        phosphorus=dict(stock=stock, net_need=need, rain=rain, capture=capture,
            guesses=dict(active_p_g_kg=[1., 3.], structure_p_g_kg=[.05, .3], resorption=[.5, .8],
                         rain_p_mg_per_litre=[.005, .05], particle_ug_m3=[1., 10.],
                         particle_p_mass_fraction=[.001, .005],
                         note='Design guesses, once: tissue and membrane P contents, resorption, the P in lunar rain and the airborne particle load are unmeasured for the Open Moon.')),
        dawn_release=dawn,
        unresolved=['A water-limited canopy model with the host forests\' own allocation and xylem chemistry.',
                    'Whether host trees pass P in sap to an attached juvenile, and at what cost to the host.',
                    'Lunar rain and aerosol phosphorus; the airborne pollen, spore and insect loads at 10-20 km.'])
