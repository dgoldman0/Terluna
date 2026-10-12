"""A first screen of what the Open Moon's committed products already say about its resources.

    python -m resources.screen
    # -> resources/results/first_screen.json

Eight results, each read from committed products with the literature values named in sources.json:
- **Inventories.** The air, the seas, the rain-fed lakes and the crust's possible pore water (geography), the CO2 in
  the air and an Earth-like biosphere's carbon (atmospheric_co2).
- **Delivery.** The mass a 500-year build brings in, and the heat its fall through the Moon's own gravity leaves if it
  is dissipated on the Moon: from the exobase (the loss response's titania cases), the ring radius or a far terminal.
- **The oxygen's hydrogen.** Oxygen shipped as water leaves hydrogen that would burn all of it back.
- **Oxygen against rock.** Weathering rock oxidises its ferrous iron, an oxygen sink set by the weathering rate
  (lunar_cycle_ecology, rescaled to the corrected runoff) and, over the operation, by erosion (conservation).
- **The shield's film.** The ring fleet's titania and silica, the film plant's 2.3-4.6 Gt a year
  (solar_shield_array/results/array_heat.json), the fresh feed at a given recovery, and the lift from an airless Moon.
- **Flooded ore.** The share of each mare under the 28% sea (geography atlas).
- **Carbon beyond the biosphere.** What a metre of developed soil and seas with Earth-like dissolved carbon can hold.
- **Fresh seas.** How fast weathering makes the seas salty.
"""
from __future__ import annotations
import csv
import hashlib
import json
import sys
from pathlib import Path
import numpy as np

from shared.constants import MOON_GM, MOON_RADIUS
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SCHEMA = 'terluna.resources.first-screen/1'
AIR = ROOT / 'research' / 'baselines' / 'feasibility' / 'reference.json'
ATLAS = ROOT / 'geography' / 'results' / 'atlas.json'
DRAINAGE = ROOT / 'geography' / 'results' / 'drainage.json'
GROUNDWATER = ROOT / 'geography' / 'results' / 'groundwater.json'
CO2 = ROOT / 'research' / 'studies' / 'atmospheric_co2' / 'results' / 'atmospheric_co2.json'
SWEEP = ROOT / 'atmosphere' / 'loss_response' / 'results' / 'euv_sweep.csv'
HEAT = ROOT / 'research' / 'studies' / 'solar_shield_array' / 'results' / 'array_heat.json'
ECOLOGY = ROOT / 'research' / 'studies' / 'lunar_cycle_ecology' / 'results' / 'lunar_cycle_ecology.json'
GATES = ROOT / 'research' / 'studies' / 'conservation' / 'results' / 'gates.json'
LEDGERS = ROOT / 'research' / 'studies' / 'joint_synthesis' / 'results' / 'ledgers.json'
OUT = HERE / 'results' / 'first_screen.json'
INPUTS = (AIR, ATLAS, DRAINAGE, GROUNDWATER, CO2, SWEEP, HEAT, ECOLOGY, GATES, LEDGERS)

YEAR_S = 3.15576e7
BUILD_YEARS = 500.0                       # decisions register: construction within 500 years
OPERATION_YEARS = 1e9
TERMINAL_RADIUS_M = 60000e3               # a far cislunar terminal, for the largest fall
RING_RADIUS_M = 20000e3
CO2_PPM = 400.0
# Literature values (sources.json).
FEO_WT = (0.05, 0.16)                     # dry land is mostly feldspathic highland (4-6 wt%); mare soils 14-20 wt%
O2_PER_FEO = 0.5 * 31.998 / (2 * 71.844)  # 2 FeO + 1/2 O2 -> Fe2O3
TIO2_DENSITY, SIO2_DENSITY = 4000.0, 2200.0   # kg/m3: anatase to rutile, fused silica
TIO2_IN_ILMENITE = 79.866 / 151.71
RECOVERY = (0.99, 0.999)
SOIL_DENSITY = 1500.0                     # kg/m3, a developed soil's bulk density
CO2_PER_KG_SOIL = 0.19                    # kg CO2 per kg, the CaO+MgO of Apollo soils (atmospheric_co2)
SEA_DIC_MOL_KG = (0.002, 0.0023)          # Earth's surface seawater
FRESH_G_KG, EARTH_SEA_G_KG = 1.0, 35.0
SKY_FLEET_H2_T = (48000.0, 71000.0)       # liners and ferries (infrastructure review, via the joint ledger)
HCO3_PER_G = {'Ca': 2 * 61.017 / 40.078, 'Mg': 2 * 61.017 / 24.305, 'Na': 61.017 / 22.990, 'K': 61.017 / 39.098}

EVIDENCE = ('A screen on committed products (the reference air, the geography inventories, the CO2 study, the loss '
            'sweep, the array heat screen, the ecology register\'s weathering, the conservation gates) with the stated '
            'literature values. Closed-form arithmetic; no transport, geochemical or industrial model is run.')
READING_RULE = ('Masses in kg, flows in kg/s or per year, energies in J/kg, heat in W/m2 averaged over the Moon '
                '(against the design sunlight at the top of the air, 293.3 W/m2), times in years. Ranges are low-high.')


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def exobase_radii():
    with open(SWEEP) as f:
        rows = [r for r in csv.DictReader(f) if r['shield'] == 'titania_stack' and r['exobase_radius_R'] not in ('', 'nan')
                and float(r['leak_fraction']) <= 3e-4]
    radii = [float(r['exobase_radius_R']) for r in rows]
    return min(radii), max(radii)


def inventories(air, atlas, drainage, groundwater, co2):
    air_mol = air['mass_kg'] / air['mean_molar_mass']
    co2_kg = air_mol * CO2_PPM * 1e-6 * 0.04401
    crust = [c['crust_kg'] for c in groundwater['crust']]
    water = dict(seas_kg=atlas['water_mass_kg'], lakes_kg=drainage['lakes']['volume_km3'] * 1e12,
                 crust_kg=[min(crust), max(crust)])
    water['total_kg'] = [water['seas_kg'] + water['lakes_kg'], water['seas_kg'] + water['lakes_kg'] + max(crust)]
    return dict(air_kg=air['mass_kg'], n2_kg=air['N2_kg'], o2_kg=air['O2_kg'], air_mol=float(f'{air_mol:.4g}'),
                co2_kg=float(f'{co2_kg:.3g}'), co2_kg_constant_gravity_study=co2['design_air']['inventory_kg_co2'],
                biosphere_co2_kg=co2['biosphere']['earth_like_land_carbon_kg_co2'], water=water,
                water_over_air=[round(w / air['mass_kg'], 1) for w in water['total_kg']])


def delivery(inv):
    lo, hi = exobase_radii()
    gm_r = MOON_GM / MOON_RADIUS
    fall = dict(exobase=[gm_r * (1 - 1 / lo), gm_r * (1 - 1 / hi)], ring_radius=gm_r * (1 - MOON_RADIUS / RING_RADIUS_M),
                far_terminal=gm_r * (1 - MOON_RADIUS / TERMINAL_RADIUS_M))
    mass = [inv['air_kg'] + w for w in inv['water']['total_kg']]
    flow = [m / (BUILD_YEARS * YEAR_S) for m in mass]
    area = 4 * np.pi * MOON_RADIUS ** 2
    heat = [flow[0] * fall['exobase'][0] / area, flow[1] * fall['far_terminal'] / area]
    return dict(build_years=BUILD_YEARS, exobase_radius_r=[round(lo, 2), round(hi, 2)],
                fall_mj_per_kg={k: ([round(x / 1e6, 2) for x in v] if isinstance(v, list) else round(v / 1e6, 2)) for k, v in fall.items()},
                mass_kg=[float(f'{m:.3g}') for m in mass], flow_kg_s=[float(f'{f:.3g}') for f in flow],
                heat_w_m2_if_dissipated_on_moon=[round(h, 0) for h in heat],
                share_of_design_sunlight=[round(h / 293.3, 2) for h in heat])


def hydrogen(inv):
    h_kg = inv['o2_kg'] * 4 * 1.008 / 31.998          # O2 + 2 H2 -> 2 H2O: two hydrogen molecules per O2
    h2_mol = h_kg / 2.016e-3
    o2_burned = h2_mol / 2 * 31.998e-3
    return dict(water_to_carry_the_oxygen_kg=float(f'{inv["o2_kg"] * 18.015 / 15.999:.3g}'), hydrogen_kg=float(f'{h_kg:.3g}'),
                mole_fraction_if_released=round(h2_mol / (inv['air_mol'] + h2_mol), 3),
                share_of_air_oxygen_it_would_burn=round(o2_burned / inv['o2_kg'], 3),
                sky_fleet_share=[float(f'{t * 1e3 / h_kg:.2g}') for t in SKY_FLEET_H2_T])


def oxygen_sink(inv, ecology, drainage, atlas, gates):
    limb = ecology['return_limb']
    runoff_now = drainage['rivers']['discharge_to_sea_m3s'] * YEAR_S / limb['land_area_m2']
    scale = runoff_now / limb['runoff_m_per_yr_over_land']
    land = (atlas['main_land_share'] + atlas['island_share'] - drainage['lakes']['area_share']) * 4 * np.pi * MOON_RADIUS ** 2
    rock = [g * scale * land / 1e3 for g in limb['soils']['mare_soil']['rock_weathered_g_m2_yr']]   # kg/yr
    sink = [rock[0] * FEO_WT[0] * O2_PER_FEO, rock[1] * FEO_WT[1] * O2_PER_FEO]
    erosion = [v['rock_gt_per_year'] * 1e12 for v in gates['sediment_return'].values()]
    long_sink = [min(erosion) * FEO_WT[0] * O2_PER_FEO, max(erosion) * FEO_WT[1] * O2_PER_FEO]
    escape = 1.4 * YEAR_S                                                 # the warm titania case, kg/yr of all air
    return dict(runoff_rescale=round(scale, 3), rock_weathered_gt_yr=[round(r / 1e12, 2) for r in rock],
                o2_sink_mt_yr=[round(s / 1e9, 1) for s in sink],
                o2_cycle_myr=[round(inv['o2_kg'] / s / 1e6, 1) for s in sink[::-1]],
                erosion_limited_o2_mt_yr=[round(s / 1e9, 1) for s in long_sink],
                air_oxygen_per_operation=[round(s * OPERATION_YEARS / inv['o2_kg'], 1) for s in long_sink],
                times_all_escape=[round(s / escape) for s in sink])


def film(array_heat, ledgers):
    fleet, design = array_heat['fleet'], array_heat['design']
    area_m2 = [a * 1e6 for a in fleet['area_km2']]
    w = fleet['window_ring_share']
    (t_a, s_a), (t_w, s_w) = design['films']['annulus'], design['films']['window']
    tio2 = [a * ((1 - w) * t_a + w * t_w) * 1e-6 * TIO2_DENSITY for a in area_m2]
    sio2 = [a * ((1 - w) * s_a + w * s_w) * 1e-6 * SIO2_DENSITY for a in area_m2]
    upkeep = [g * 1e12 for g in fleet['film_upkeep_Gt_per_year']]
    fresh = {f'{r:g}': [u * (1 - r) / YEAR_S for u in upkeep] for r in RECOVERY}
    lift = MOON_GM / MOON_RADIUS * (1 - MOON_RADIUS / (2 * RING_RADIUS_M))     # to a circular orbit at the ring radius
    air = ledgers['air']['mass_kg']
    return dict(optical_mass_gt=fleet['optical_mass_Gt'], tio2_gt=[round(x / 1e12, 2) for x in tio2],
                sio2_gt=[round(x / 1e12, 1) for x in sio2], ilmenite_for_the_tio2_gt=[round(x / TIO2_IN_ILMENITE / 1e12, 1) for x in tio2],
                upkeep_gt_yr=fleet['film_upkeep_Gt_per_year'], upkeep_kg_s=[round(u / YEAR_S) for u in upkeep],
                fresh_feed_kg_s={k: [round(x) for x in v] for k, v in fresh.items()},
                fresh_feed_atmospheres_per_operation={k: [round(x * YEAR_S * OPERATION_YEARS / air, 1) for x in v] for k, v in fresh.items()},
                lift_to_ring_mj_kg=round(lift / 1e6, 2),
                lift_fleet_in_100_years_gw=[round(m * 1e12 * lift / (100 * YEAR_S) / 1e9) for m in fleet['optical_mass_Gt']])


def flooded_ore(atlas):
    maria = atlas['mare_flooding']
    return dict(maria=[dict(name=m['name'], flooded=m['flooded']) for m in maria],
                fully_flooded=[m['name'] for m in maria if m['flooded'] >= 0.99],
                partly_dry=[dict(name=m['name'], dry=round(1 - m['flooded'], 2)) for m in maria if m['flooded'] < 0.99])


def carbon_beyond(inv, atlas, drainage):
    land = (atlas['main_land_share'] + atlas['island_share'] - drainage['lakes']['area_share']) * 4 * np.pi * MOON_RADIUS ** 2
    metre = land * SOIL_DENSITY * CO2_PER_KG_SOIL
    seas = [(inv['water']['seas_kg'] + inv['water']['lakes_kg']) * d * 0.04401 for d in SEA_DIC_MOL_KG]
    return dict(co2_per_metre_of_soil_kg=float(f'{metre:.3g}'), per_metre_over_air_co2=round(metre / inv['co2_kg'], 1),
                seas_and_lakes_dic_co2_kg=[float(f'{s:.3g}') for s in seas],
                seas_over_air_co2=[round(s / inv['co2_kg'], 2) for s in seas])


def fresh_seas(ecology, drainage, inv):
    limb = ecology['return_limb']
    runoff_now = drainage['rivers']['discharge_to_sea_m3s'] * YEAR_S / limb['land_area_m2']
    scale = runoff_now / limb['runoff_m_per_yr_over_land']
    exports = limb['soils']['mare_soil']['export_g_m2_yr']
    cations = [sum(exports[e][i] for e in ('Ca', 'Mg', 'Na', 'K')) * scale for i in (0, 1)]      # g/m2/yr of land
    # Weathering's cations arrive with the bicarbonate that balances their charge, and with dissolved silica.
    bicarbonate = [sum(exports[e][i] * HCO3_PER_G[e] for e in HCO3_PER_G) * scale for i in (0, 1)]
    silica = [exports['Si'][i] * 60.084 / 28.086 * scale for i in (0, 1)]
    water = inv['water']['seas_kg'] + inv['water']['lakes_kg']
    rise = lambda g_m2: [g * limb['land_area_m2'] / 1e3 / water * 1e3 for g in g_m2]               # g/kg per year
    ions, dissolved = rise(cations), rise([c + b + s for c, b, s in zip(cations, bicarbonate, silica)])
    return dict(cations_g_m2_yr=[round(c, 2) for c in cations], bicarbonate_g_m2_yr=[round(b, 1) for b in bicarbonate],
                silica_g_m2_yr=[round(s, 1) for s in silica],
                rise_mg_kg_per_year=dict(cations=[round(p * 1e3, 3) for p in ions], all_dissolved=[round(p * 1e3, 3) for p in dissolved]),
                years_to_1_g_kg=dict(cations=[round(FRESH_G_KG / p, -2) for p in ions[::-1]],
                                     all_dissolved=[round(FRESH_G_KG / p, -2) for p in dissolved[::-1]]),
                years_to_earth_salinity_cations=[round(EARTH_SEA_G_KG / p, -3) for p in ions[::-1]],
                note='Calcite and silica precipitate and biology takes up silica, so the dissolved total is an upper bound '
                     'on the rise; the imported water\'s own salts come on top.')


def main(argv=None) -> int:
    air, atlas, drainage, groundwater, co2, array_heat, ecology, gates, ledgers = (
        json.loads(p.read_text()) for p in (AIR, ATLAS, DRAINAGE, GROUNDWATER, CO2, HEAT, ECOLOGY, GATES, LEDGERS))
    inv = inventories(air['inventory'], atlas, drainage, groundwater, co2)
    files = ['resources/screen.py']
    product = dict(
        schema=SCHEMA,
        producer=dict(domain='resources', files={f: digest(ROOT / f) for f in files},
                      inputs={str(p.relative_to(ROOT)): digest(p) for p in INPUTS}, constants=constants_used(files)),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        inventories=inv, delivery=delivery(inv), hydrogen=hydrogen(inv),
        oxygen_sink=oxygen_sink(inv, ecology, drainage, atlas, gates), film=film(array_heat, ledgers),
        flooded_ore=flooded_ore(atlas), carbon_beyond_the_biosphere=carbon_beyond(inv, atlas, drainage),
        fresh_seas=fresh_seas(ecology, drainage, inv))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    print(json.dumps({k: v for k, v in product.items() if k not in ('producer', 'evidence', 'reading_rule', 'flooded_ore')}, indent=1)[:7000])
    return 0


if __name__ == '__main__':
    sys.exit(main())
