"""Bounded aerial food-web requirements; python -m research.studies.aerial_ecology.run.

No population dynamics or climate is simulated. Trait values below are explicit
sensitivities, not measured lunar organisms. All masses use kg unless labelled.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

from biosphere.ecology.screen import fall_speed
from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_RADIUS, SYNODIC_MONTH_DAYS
from shared.constants import MOON_SURFACE_GRAVITY
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.aerial-ecology/1'
INPUTS = {
    'biosphere/ecology/results/people.json': 'terluna.biosphere.ecology-people/1',
    'research/studies/sky_ships/results/sky_ships.json': None,
    'shared/scenarios/population.json': 'terluna.scenario.population/1',
    'research/studies/aerial_ecology/sources.json': None,
}
OUT = HERE / 'results/aerial_ecology.json'
C_TO_P = 106 * 12.011 / 30.974  # Redfield atomic ratio, a sensitivity reference


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def growth_budget(doubling_days, favourable_fraction, removal_days, mortality_per_day=0.01):
    """Whole-cycle favourable fraction includes hydration, light and nutrients.

    The loss time includes rain and settling. Cloud fraction is not a trajectory's
    wet duty. This linear invasion criterion cannot determine a carrying capacity.
    """
    if doubling_days <= 0 or removal_days <= 0 or not 0 <= favourable_fraction <= 1:
        raise ValueError('positive times and a fraction in [0, 1] required')
    loss = 1 / removal_days + mortality_per_day
    r = favourable_fraction * math.log(2) / doubling_days - loss
    return dict(net_growth_per_day=r, minimum_favourable_fraction=loss * doubling_days / math.log(2),
                multiplier_per_lunar_cycle=math.exp(r * SYNODIC_MONTH_DAYS))


def feeding(concentration_mg_m3, speed_m_s, drag_coefficient, density=1.2,
            capture=0.5, assimilation=0.6, propulsive_efficiency=0.25, basal_w_m2=1.0):
    """Unit frontal capture area; drag coefficient includes filter and carrier.

    18 MJ/kg is assumed dry-food energy. The speed is through-air speed. The
    passive case requires an external reaction (tether or quantified wind shear).
    """
    if concentration_mg_m3 < 0 or speed_m_s < 0 or drag_coefficient < 0:
        raise ValueError('negative feeding input')
    if not (0 < capture <= 1 and 0 < assimilation <= 1 and 0 < propulsive_efficiency <= 1):
        raise ValueError('invalid efficiency')
    capture_kg_s = concentration_mg_m3 * 1e-6 * speed_m_s * capture
    gross = capture_kg_s * 18e6
    mechanical = 0.5 * density * drag_coefficient * speed_m_s**3
    propulsion = mechanical / propulsive_efficiency
    denominator = speed_m_s * capture * assimilation * 18e6
    return dict(captured_g_day_m2=capture_kg_s * JULIAN_DAY * 1e3,
                gross_food_w_m2=gross, assimilated_w_m2=gross * assimilation,
                drag_mechanical_w_m2=mechanical, propulsion_metabolic_w_m2=propulsion,
                powered_net_w_m2=gross * assimilation - propulsion - basal_w_m2,
                passive_net_w_m2=gross * assimilation - basal_w_m2,
                powered_break_even_mg_m3=(propulsion + basal_w_m2) / denominator * 1e6 if denominator else None,
                passive_break_even_mg_m3=basal_w_m2 / denominator * 1e6 if denominator else None)


def optical_depth(concentration_mg_m3, depth_m, diameter_um, density_kg_m3=1000, q_ext=2):
    """Dry monodisperse spheres, geometric extinction; not cloud forcing or Mie optics."""
    mass_extinction = 3 * q_ext / (2 * density_kg_m3 * diameter_um * 1e-6)
    return concentration_mg_m3 * 1e-6 * depth_m * mass_extinction


def settling(diameter_um):
    """Unhydrated spheres in an illustrative 10-km air parcel; rain is separate."""
    speed, rho = fall_speed(diameter_um*1e-6, 1000, MOON_SURFACE_GRAVITY,
                           100000, 287.35, 0.02897)
    return dict(diameter_um=diameter_um, speed_m_s=float(speed),
                still_air_days_per_km=float(1000 / speed / JULIAN_DAY),
                note='Stokes-Cunningham screen; hydrated droplets and rain have other sizes and loss rates.')


def aerophyte_mass(lift_kg_m3, dry_kg_m2, water_fraction=0.9, skin_kg_m2=1, usable_fraction=0.5):
    """Sphere projected area pi R²; skin covers 4 pi R².

    Only usable_fraction of gross lift is allocated to skin and wet biomass;
    the remainder is a mass reserve. Gas manufacture/leakage and storms are open.
    """
    if lift_kg_m3 <= 0 or dry_kg_m2 < 0 or not 0 <= water_fraction < 1 or not 0 < usable_fraction <= 1:
        raise ValueError('invalid mass budget')
    wet = dry_kg_m2 / (1 - water_fraction)
    radius = 3 * (wet + 4 * skin_kg_m2) / (4 * lift_kg_m3 * usable_fraction)
    return dict(wet_biomass_kg_m2=wet, minimum_radius_m=radius,
                water_kg_m2=wet - dry_kg_m2, skin_kg_per_projected_m2=4 * skin_kg_m2)


def production(coverage, npp_g_c_m2_year, dry_stock_kg_m2, harvest=0.05,
               edible=1.0, processing=1.0, c_to_p=C_TO_P):
    """Coverage is sum of projected organism areas / lunar SURFACE area, not disk area.

    Overlap is excluded in this first screen. NPP is net of autotroph respiration;
    trophic production, retained biomass and human harvest share that one budget.
    """
    if not 0 <= coverage <= 1 or not all(0 <= x <= 1 for x in (harvest, edible, processing)):
        raise ValueError('invalid fraction')
    area = 4 * math.pi * MOON_RADIUS**2 * coverage
    carbon = area * npp_g_c_m2_year / 1000
    dry = carbon / 0.45
    harvested = dry * harvest
    return dict(projected_area_m2=area, npp_kg_c_year=carbon, dry_production_kg_year=dry,
                dry_stock_kg=area * dry_stock_kg_m2,
                p_assimilation_kg_year=carbon / c_to_p,
                p_stock_kg=area * dry_stock_kg_m2 * 0.45 / c_to_p,
                harvested_dry_kg_year=harvested,
                food_energy_people=harvested * edible * processing * 4000 / (2500 * JULIAN_YEAR_DAYS),
                chemical_npp_w_per_projected_m2=npp_g_c_m2_year / 1000 / 0.45 * 18e6 / (JULIAN_YEAR_DAYS * JULIAN_DAY),
                stock_to_production_years=area * dry_stock_kg_m2 / dry if dry else None)


def run():
    data = {}
    for path, schema in INPUTS.items():
        data[path] = json.loads((ROOT / path).read_text())
        if schema and data[path]['schema'] != schema:
            raise ValueError(path)
    people = data['biosphere/ecology/results/people.json']
    air = data['research/studies/sky_ships/results/sky_ships.json']['air']
    area = 4 * math.pi * MOON_RADIUS**2
    central = production(0.01, 1000, 1)
    p_input = [area * f / 1000 for f in people['aloft']['dust_p_g_m2_yr']]
    feed = [dict(concentration_mg_m3=c, speed_m_s=v, drag_coefficient=cd, **feeding(c, v, cd))
            for c in (0.01, 0.1, 1.0) for v in (1., 3., 10.) for cd in (0.02, 0.1, 1.)]
    growth = [dict(doubling_days=t, favourable_fraction=f, removal_days=r, **growth_budget(t, f, r))
              for t in (1., 2.5, 3.6, 19.5) for f in (0.02, 0.1, 0.5) for r in (3., 30., 100.)]
    # Carbon reserve for a fully dark half-cycle: fractional carbon respiration per day.
    night = [dict(respiration_fraction_per_day=r,
                  reserve_fraction_of_initial_carbon=1 - math.exp(-r * SYNODIC_MONTH_DAYS / 2))
             for r in (0.001, 0.01, 0.03)]
    files = [Path(__file__), ROOT/'biosphere/ecology/screen.py']
    out = dict(
        schema=SCHEMA,
        producer=dict(files={str(p.relative_to(ROOT)): digest(p) for p in files},
                      inputs={p: digest(ROOT / p) for p in INPUTS}, constants=constants_used(files)),
        evidence='Analytic requirements and sensitivity scenarios. No organism, sustainable harvest, species count, '
                 'cloud feedback, atmospheric concentration or nutrient closure is validated.',
        reading_rule='All kg are dry biomass unless explicitly wet, carbon or phosphorus. Coverage uses lunar '
                     'surface area. Annual means use a Julian year. Feeding is per square metre of capture area.',
        assumptions=dict(carbon_share_dry=0.45, carbon_phosphorus_mass_ratio=C_TO_P,
                         food_kcal_kg_dry=4000, energy_kcal_person_day=2500,
                         feeding_energy_j_kg_dry=18e6, capture=0.5, assimilation=0.6,
                         propulsive_efficiency=0.25, basal_w_m2=1, mortality_per_day=0.01),
        baseline=dict(warm_air_volume_km3=people['sky']['warm_air_volume_km3'],
                      microbial_wet_t_at_1pg_cell={str(n): people['sky']['warm_air_volume_km3'] * 1e9 * n * 1e-15 / 1000
                                                 for n in (1e4, 1e5)},
                      p_fallout_analogue_kg_year=p_input,
                      single_pass_dry_kg_year=[p * C_TO_P / 0.45 for p in p_input]),
        growth=growth, feeding=feed, dark_half_cycle=night,
        settling=[settling(d) for d in (1., 10., 25.)],
        haze=[dict(concentration_mg_m3=c, depth_m=h, diameter_um=d, tau=optical_depth(c, h, d))
              for c in (0.01, 0.1, 1.) for h in (100., 1000.) for d in (1., 10., 100.)],
        aerophyte_mass=[dict(height_km=a['height_km'], dry_kg_m2=b, water_fraction=w,
                             **aerophyte_mass(a['lift_hydrogen_kg_m3'], b, w))
                        for a in air if a['height_km'] in (10, 20, 40)
                        for b in (1., 10.) for w in (0.8, 0.9, 0.95)],
        production=[dict(coverage=c, npp_g_c_m2_year=p, dry_kg_m2=b, **production(c, p, b))
                    for c in (0.0001, 0.001, 0.01) for p in (500., 1000., 2000.) for b in (1., 10.)],
        central=central,
        food_sensitivity=[dict(edible=e, processing=p, **production(0.01, 1000, 1, edible=e, processing=p))
                          for e, p in ((1., 1.), (0.5, 0.8), (0.25, 0.7))],
        trophic_production_kg_year=[dict(transfer_fraction=t,
            primary_consumer=(1 - 0.05) * central['dry_production_kg_year'] * t,
            secondary_consumer=(1 - 0.05) * central['dry_production_kg_year'] * t*t)
            for t in (0.03, 0.06, 0.15)],
        biodiversity=dict(species_count=None, reason='Energy and habitat budgets do not determine species richness.',
                          candidate_guilds=['cloud phototrophs', 'C1 heterotrophs', 'droplet grazers',
                          'photosynthetic aerophytes', 'collectors', 'grazers', 'predators', 'decomposers']),
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(out, indent=2, allow_nan=False) + '\n')
    return out


if __name__ == '__main__':
    run()
