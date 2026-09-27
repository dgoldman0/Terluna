#!/usr/bin/env python3
"""Ecology of the lunar cycle: the small calculations behind the ideas register (README.md).

    python -m research.studies.lunar_cycle_ecology.run      # writes research/studies/lunar_cycle_ecology/results/lunar_cycle_ecology.json

A register of design ideas with light calculations, not an ecosystem model. It reads the biosphere's plant and fruit
products, the geography atlas and drainage, the design air and the shared constants, and computes:

1. leaves: how much warmer the Moon's CO2 puts the C3-C4 crossover, and the growth the plant needs to replace its
   canopy at a range of leaf lifespans, against what the whole plant keeps at each fruit share;
2. the night: land temperatures through the night, the oxygen and CO2 that the night's respiration and the decay of a
   dusk fruit fall take from or give to the air column, and how fast dusk moves west over the ground;
3. flight: weight, air density, the induced power to hover and the energy per kilometre on the Moon against Earth, and
   the sugar a moth-sized and a bird-sized flier burn per kilometre;
4. distances from land to the nearest sea and to any standing water, from the geography atlas grid;
5. the missing return limb of the rock cycle: what weathering carries off the land, element by element, at the Moon's
   runoff and land temperature (the basalt law of the CO2 study); how fast the seas would gather each cation; and what
   the return paths can bring back: insects emerging from the seas at Lake Myvatn's highest rate (element by element),
   dimethyl sulfide from seas as active as Earth's ocean, and salt works holding the seas' composition (against the
   world's salt output); with the migrant conveyor sized for dissolved phosphorus;
6. glowing lures: their energy cost against the plant's night upkeep.
"""
from __future__ import annotations
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
if __package__ in (None, ''):                  # executed as a file path rather than with -m
    sys.path.insert(0, str(ROOT))
from shared.constants import (GAS_CONSTANT, MOON_RADIUS, MOON_SURFACE_GRAVITY, PLANCK, SPEED_OF_LIGHT,
                              STANDARD_GRAVITY, SYNODIC_MONTH_DAYS)
from atmosphere.radiative_convective import thermodynamics as th
from research.studies.atmospheric_co2.run import GLASS_FACTOR, KINETICS, R_KJ, weathering_mol_km2_yr
from biosphere.canopy.model import WORLDS

HERE = Path(__file__).resolve().parent
SCHEMA = 'terluna.research.lunar-cycle-ecology/1'
PLANT = ROOT / 'biosphere' / 'canopy' / 'results' / 'plant.json'
FRUIT = ROOT / 'biosphere' / 'canopy' / 'results' / 'fruit.json'
ATLAS = ROOT / 'geography' / 'results' / 'atlas.json'
ATLAS_GRID = ROOT / 'geography' / 'products' / 'atlas_28pct_4ppd.npz'
DRAINAGE = ROOT / 'geography' / 'results' / 'drainage.json'
EXPECTED = {PLANT: 'terluna.biosphere.plant-carbon-cycle/1', FRUIT: 'terluna.biosphere.fruit-carbon/2',
            ATLAS: 'terluna.geography.atlas/1', DRAINAGE: 'terluna.geography.drainage/1'}

YEAR_DAYS = 365.25
CYCLES_PER_YEAR = YEAR_DAYS / SYNODIC_MONTH_DAYS
DESIGN_PRESSURE_PA = 121590.0              # decisions register: 1.2 atm
DESIGN_CO2_PPM = 400.0                     # the atmosphere and climate runs
MOON_AIR_K = 294.0                         # decisions register: the Moon settles near 294 K
EARTH_PRESSURE_PA, EARTH_AIR_K = 101325.0, 288.15
GLUCOSE_KJ_PER_G_C = 2803.0 / 72.06        # energy of respiring carbohydrate, kJ per g of carbon
SUGAR_KJ_PER_G = 2803.0 / 180.16          # glucose; FAO gives 15.7 kJ/g for monosaccharides (sources.json)
EARTH_LAND_M2 = 1.489e14                   # Earth's land area, for Earth's fluxes per m2 of land
SEA_SHARE_MIN = 0.005                      # atlas water bodies counted as seas: at least 0.5% of the Moon

# Leaves on evergreen woody plants live from months to years (sources.json).
LEAF_LIFESPAN_MONTHS = (6.0, 12.0, 24.0, 36.0)

# Fliers: mass (kg), lift-to-drag ratio and flight-muscle efficiency, and Earth airspeeds (m/s) (sources.json): a
# Bogong-moth-sized flier (0.33 g; insect muscle efficiency under 10%; noctuid moths migrate at 3-5 m/s) and a small
# bird (gliding lift-to-drag up to 12.5, muscle efficiency 0.13-0.23; passerines migrate at 10-15 m/s).
FLIERS = {'moth_0.33g': (0.33e-3, 4.0, 0.10), 'bird_100g': (0.1, 10.0, 0.18)}
AIRSPEED_EARTH_M_S = {'moth_0.33g': (3.0, 5.0), 'bird_100g': (10.0, 15.0)}

# Sugar in the fruit: watermelon, raw, 6.2 g per 100 g fresh (USDA FoodData Central 167765).
FRUIT_SUGAR_PER_100G = 6.2

# Earth's natural river water (mg/L), what Earth's rivers carry off the land, as a reference scale (sources.json):
# potassium median and discharge-weighted, phosphate P, total dissolved P, and particulate P as a concentration
# (NEWS 2 and Meybeck's totals over 37,400 km3 a year).
RIVER_MG_L = dict(potassium=(1.06, 1.29), phosphate_p=(0.010, 0.0125), dissolved_p=(0.025, 0.025),
                  particulate_p=(0.18, 0.53))

# Rock: oxide weight % (S as the element) of average lunar soils (Lunar Sourcebook Table 7.15): mare, the median of six
# mare soils (Apollo 11, 12, 15 mare, 17 mare, Luna 16, 24); highland, the Apollo 16 average; KREEP-rich, the Apollo 14
# average. Earth: upper continental crust (Rudnick & Gao 2003), for the stocks only.
ROCK_WT_PCT = {'mare_soil': dict(SiO2=43.05, MgO=9.35, CaO=11.4, Na2O=0.365, K2O=0.16, P2O5=0.13, S=0.12),
               'highland_soil': dict(SiO2=45.0, MgO=5.7, CaO=15.7, Na2O=0.46, K2O=0.17, P2O5=0.11, S=0.07),
               'kreep_rich': dict(SiO2=48.1, MgO=9.4, CaO=10.7, Na2O=0.70, K2O=0.55, P2O5=0.51, S=None),
               'earth_upper_crust': dict(K2O=2.80, P2O5=0.15, CaO=3.59, MgO=2.48)}
# Element mass per unit oxide mass.
ELEMENT_OF = dict(Ca=('CaO', 40.078 / 56.077), Mg=('MgO', 24.305 / 40.304), Na=('Na2O', 2 * 22.990 / 61.979),
                  K=('K2O', 2 * 39.098 / 94.196), P=('P2O5', 2 * 30.974 / 141.94), S=('S', 1.0),
                  Si=('SiO2', 28.086 / 60.084))
SOIL_BULK_KG_M3, SOIL_DEPTH_M = 1500.0, 1.0
K_IN_K2O, P_IN_P2O5 = ELEMENT_OF['K'][1], ELEMENT_OF['P'][1]
CAO_G_MOL, MGO_G_MOL = 56.077, 40.304

# Earth's sea water at salinity 35, g/kg (DOE 1994 handbook, as tabulated for the reference composition), Earth's ocean
# area, the ocean's dimethyl sulfide emission (Lana et al. 2011, Tg S/yr), and the world's salt output (USGS 2025,
# 2024 estimate, kg/yr).
EARTH_SEAWATER_G_KG = dict(Na=10.78, Mg=1.284, Ca=0.412, K=0.399)
EARTH_OCEAN_M2 = 3.61e14
OCEAN_DMS_TG_S_YR = (17.6, 28.1, 34.4)
WORLD_SALT_KG_YR = 2.8e11

# Migrant bodies, % of dry mass (sources.json): insect P 0.62-1.3 (0.7 for a moth-sized insect), K 0.6-2.0; fish P
# 1-5%, most of it in bone. Dry mass as a share of fresh is an assumption.
BODY_P_PCT = {'insect': 0.7, 'fish': 2.5}
INSECT_PCT = dict(P=0.7, K=1.2, Mg=0.25, Ca=0.3)      # insect body, % of dry mass (Mg 0.08-0.40, Ca under 0.3)
INSECT_DRY_SHARE = 0.3

# Earth's biological and dust returns of phosphorus, for scale (sources.json): seabird colonies worldwide in the
# breeding season (Gg/yr); high-flying insect migrants over southern Britain (kg P/yr over km2); Lake Myvatn midges
# within 50 m of the shore in a high year (kg P/ha/yr) and the lake's emergence (g dry/m2/yr); African dust on the
# Amazon and the Amazon's hydrological loss (g P/ha/yr); Earth's natural dissolved and particulate river P (Tg/yr).
EARTH_P = dict(seabird_colonies_gg_yr=99.0, uk_insect_migration=(10000.0, 70000.0), myvatn_shore_kg_ha_yr=1.0,
               myvatn_emergence_g_m2_yr=(0.15, 3.7), amazon_dust_g_ha_yr=(7.0, 39.0),
               amazon_hydrological_loss_g_ha_yr=(8.0, 40.0), rivers_dissolved_tg_yr=1.0,
               rivers_particulate_tg_yr=(6.56, 20.0))
DEATHS_PER_VISIT = (0.1, 0.5)              # share of visiting migrants that die or are caught on land, a scenario

# Glowing lures (sources.json): the fungal pathway's emission (520-530 nm) and quantum yields (0.1-1.7% measured for
# fungal luciferin, 41% for firefly); the brightest engineered glowing plants (flowers, 6.47e10 photons per minute per
# cm2); a cost for making the luciferin of up to ten times the photon energy per reaction; lures that draw insects from
# 10-25 m around them.
LURE_NM = 530.0
LURE_QUANTUM_YIELD = (0.001, 0.017, 0.41)
LURE_FLUX_PER_CM2_S = 6.47e10 / 60
LURE_AREAS_CM2 = (10.0, 1e3, 1e4)
LURE_SUBSTRATE_FACTOR = (1.0, 10.0)
LURE_REACH_M = (10.0, 25.0)


def read(path):
    data = json.loads(path.read_text())
    if data.get('schema') != EXPECTED[path]:
        raise ValueError(f'{path}: expected schema {EXPECTED[path]}, found {data.get("schema")}')
    return data


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()[:16]


def c4_crossover_shift(t_k=300.0):
    """How much warmer the Moon's CO2 puts the temperature where a C4 leaf's light-limited efficiency overtakes a C3
    leaf's. C3's light-limited yield falls as Gamma*/ci rises; C4's is taken as fixed; Gamma* follows Bernacchi et al.
    (2001), so at a fixed ci/ca the crossover holds Gamma*/ca fixed."""
    moon, earth = WORLDS['moon']['co2_pa'], WORLDS['earth']['co2_pa']
    energy = KINETICS['gamma_star'][1]                                  # kJ/mol
    return dict(moon_co2_pa=moon, earth_co2_pa=earth, at_k=t_k,
                shift_k=math.log(moon / earth) * R_KJ * t_k ** 2 / energy)


def leaves(fruit):
    """Growth needed to replace the canopy at each leaf lifespan, against what the plant keeps at each fruit share."""
    canopy_c = fruit['settings']['canopy_carbon_g_m2']
    need = {f'{m:g}_months': canopy_c / (m * 30.4375 / SYNODIC_MONTH_DAYS) for m in LEAF_LIFESPAN_MONTHS}
    keeps = {}
    for name in ('earth_like', 'idle_quarter'):
        shares = fruit['results']['near']['0'][name]['day_fruit']['best']['fed']['shares']
        keeps[name] = {s: r['plant_growth'] for s, r in shares.items() if r}
    return dict(canopy_carbon_g_m2=canopy_c, leaf_replacement_g_c_per_cycle=need,
                plant_keeps_g_c_per_cycle_equator=keeps)


def air(pressure_pa, t_k, gravity):
    a = th.earthlike_air(pressure_pa, DESIGN_CO2_PPM)
    m = a.molar_mass                                      # kg/mol
    column = pressure_pa / gravity                        # kg of air per m2
    mass_share = lambda gas, mm: a.fractions[gas] * mm / m
    return dict(density_kg_m3=pressure_pa * m / (GAS_CONSTANT * t_k), column_kg_m2=column,
                o2_kg_m2=column * mass_share('O2', 0.031998), co2_kg_m2=column * mass_share('CO2', 0.04401),
                o2_mole_fraction=a.fractions['O2'])


def night(plant, fruit):
    """Night temperatures, what the night takes from the air column, and how fast dusk travels."""
    temps = {side: {band: dict(night_mean_c=c['night_mean_c'], coldest_c=c['coldest_c'])
                    for band, c in bands.items()} for side, bands in plant['climate'].items()}
    moon = air(DESIGN_PRESSURE_PA, MOON_AIR_K, MOON_SURFACE_GRAVITY)
    r = plant['results']['near']['0']['earth_like']
    fall = fruit['results']['near']['0']['earth_like']['day_fruit']['best']['fed']['shares']
    draw = {}
    for share in ('0.3', '0.6'):
        c = r['respiration_in_the_dark'] + fall[share]['fruit_carbon']        # g C per m2 over one night, at most
        draw[share] = dict(carbon_g_m2=c, o2_share_of_column=c * 32.0 / 12.011 / 1e3 / moon['o2_kg_m2'],
                           co2_share_of_column=c * 44.01 / 12.011 / 1e3 / moon['co2_kg_m2'])
    circumference = 2 * math.pi * MOON_RADIUS / 1e3
    dusk = {f'{lat}': circumference * math.cos(math.radians(lat)) / (SYNODIC_MONTH_DAYS * 24) for lat in (0, 30, 60)}
    # The dusk fruit fall as food for the night, against the plants' own night upkeep over the dark hours.
    food = {}
    for share in ('0.3', '0.6'):
        sugar = fall[share]['harvest_kg'] * 10 * FRUIT_SUGAR_PER_100G                 # g per m2 per cycle
        food[share] = dict(sugar_g_m2=sugar, mean_power_over_dark_hours_w_m2=sugar * SUGAR_KJ_PER_G * 1e3
                           / (r['dark_hours'] * 3600))
    return dict(land_temperatures_c=temps, air_column=moon, night_draw_equator=draw, dusk_speed_km_h=dusk,
                night_upkeep_w_m2=night_upkeep(plant), dusk_fruit_fall_equator=food)


def night_upkeep(plant):
    """The plant's mean respiration power in the dark hours (W per m2 of ground), equator, near side."""
    out = {}
    for name in ('earth_like', 'idle_quarter'):
        r = plant['results']['near']['0'][name]
        out[name] = r['respiration_in_the_dark'] * GLUCOSE_KJ_PER_G_C * 1e3 / (r['dark_hours'] * 3600)
    return out


def flight():
    """Moon against Earth for the same animal, and sugar per kilometre for two example fliers."""
    moon = air(DESIGN_PRESSURE_PA, MOON_AIR_K, MOON_SURFACE_GRAVITY)['density_kg_m3']
    earth = air(EARTH_PRESSURE_PA, EARTH_AIR_K, STANDARD_GRAVITY)['density_kg_m3']
    g = MOON_SURFACE_GRAVITY / STANDARD_GRAVITY
    fliers = {}
    for name, (mass, lift_drag, efficiency) in FLIERS.items():
        per_km = {w: mass * grav * 1e3 / lift_drag / efficiency / 1e3 / SUGAR_KJ_PER_G     # g of sugar per km
                  for w, grav in (('moon', MOON_SURFACE_GRAVITY), ('earth', STANDARD_GRAVITY))}
        speed = math.sqrt(g * earth / moon)
        fliers[name] = dict(mass_kg=mass, lift_to_drag=lift_drag, muscle_efficiency=efficiency,
                            sugar_g_per_km=per_km, airspeed_m_s=dict(
                                earth=AIRSPEED_EARTH_M_S[name], moon=[v * speed for v in AIRSPEED_EARTH_M_S[name]]))
    return dict(air_density_kg_m3=dict(moon=moon, earth=earth), weight_ratio=g,
                hover_induced_power_ratio=g ** 1.5 * math.sqrt(earth / moon),
                energy_per_km_ratio_at_fixed_lift_to_drag=g,
                minimum_power_speed_ratio=math.sqrt(g * earth / moon), fliers=fliers)


def distances(atlas):
    """Area-weighted distances from land cells to the nearest sea and to any standing water (great-circle km)."""
    from scipy.spatial import cKDTree
    grid = np.load(ATLAS_GRID)
    meta = json.loads(str(grid['metadata']))
    if meta['schema'] != 'terluna.geography.atlas-grid/1':
        raise ValueError(f'{ATLAS_GRID}: unexpected schema {meta["schema"]}')
    lat, lon = np.meshgrid(np.radians(grid['lat_deg']), np.radians(grid['lon_deg']), indexing='ij')
    xyz = np.stack([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)], axis=-1)
    label = grid['water_label']
    seas = {b['label'] for b in atlas['bodies'] if b['share'] >= SEA_SHARE_MIN}
    land = label == 0
    weight = np.cos(lat)[land]
    out = {}
    for name, target in (('sea', np.isin(label, list(seas))), ('any_water', label > 0)):
        chord, _ = cKDTree(xyz[target]).query(xyz[land])
        km = 2 * MOON_RADIUS / 1e3 * np.arcsin(np.minimum(chord / 2, 1.0))
        order = np.argsort(km)
        cum = np.cumsum(weight[order]) / weight.sum()
        pct = {f'p{p}': float(km[order][np.searchsorted(cum, p / 100)]) for p in (50, 75, 90, 99)}
        within = {f'{d}_km': float(weight[km <= d].sum() / weight.sum()) for d in (50, 100, 200, 500)}
        out[name] = dict(**pct, max_km=float(km.max()), land_share_within=within)
    names = sorted((b['share'], b.get('water_names', [])[:3]) for b in atlas['bodies'] if b['label'] in seas)[::-1]
    return dict(seas_counted=[dict(share=s, names=n) for s, n in names], sea_share_min=SEA_SHARE_MIN,
                grid_pixels_per_degree=int(round(1 / abs(grid['lat_deg'][1] - grid['lat_deg'][0]))), **out)


def return_limb(plant, fruit, drainage, atlas, dist, fly):
    """The rock cycle's missing return limb: what weathering carries off the land element by element, how fast the seas
    would gather it, and what the biological, air and engineered return paths can bring back."""
    moon_area = 4 * math.pi * MOON_RADIUS ** 2
    land_area = moon_area * (1 - atlas['water_share'])
    sea_area = moon_area * atlas['water_share']
    sea_kg = atlas['water_mass_kg']
    runoff_m = drainage['rivers']['discharge_to_sea_m3s'] * YEAR_DAYS * 86400 / land_area
    land_c = float(np.mean(plant['climate']['near']['0']['by_hour_angle_c']))
    co2 = weathering_mol_km2_yr(runoff_m * 1e3, land_c) / 1e6                        # mol CO2 per m2 per year
    kg = SOIL_BULK_KG_M3 * SOIL_DEPTH_M
    element = lambda ox, e: (ox.get(ELEMENT_OF[e][0]) or 0.0) / 100 * ELEMENT_OF[e][1]
    stock, soils = {}, {}
    for rock, ox in ROCK_WT_PCT.items():
        stock[rock] = {f'{e}_kg_m2': kg * element(ox, e) for e in ('Ca', 'Mg', 'Na', 'K', 'P') if ELEMENT_OF[e][0] in ox}
        if rock == 'earth_upper_crust':
            continue
        cations = ox['CaO'] / 100 / CAO_G_MOL + ox['MgO'] / 100 / MGO_G_MOL             # mol Ca + Mg per g of rock
        rock_g = [co2 / 2 / cations * f for f in GLASS_FACTOR]                          # g of rock per m2 per year
        export = {e: [r * element(ox, e) for r in rock_g] for e in ('Ca', 'Mg', 'Na', 'K', 'P', 'S', 'Si')
                  if ox.get(ELEMENT_OF[e][0]) is not None}
        # Years for the seas to reach Earth's sea-water concentration of each cation if all of it stayed dissolved.
        years = {e: [EARTH_SEAWATER_G_KG[e] * sea_kg / (x * land_area) for x in export[e]] for e in EARTH_SEAWATER_G_KG}
        soils[rock] = dict(rock_weathered_g_m2_yr=rock_g, years_to_weather_top_metre=[kg * 1e3 / r for r in rock_g],
                           export_g_m2_yr=export, sea_kg_yr={e: [x * land_area / 1e3 for x in v] for e, v in export.items()},
                           years_to_earth_seawater_if_dissolved=years)
    mare = soils['mare_soil']['export_g_m2_yr']
    rivers = {k: [runoff_m * c for c in v] for k, v in RIVER_MG_L.items()}             # g per m2 of land per year
    # Biology: what insects emerging from the seas at Mývatn's highest rate, all landing on land, would bring back of
    # each element, against that element's export (dissolved phosphorus: Earth's river load at the Moon's runoff).
    deposit = EARTH_P['myvatn_emergence_g_m2_yr'][1] * sea_area / land_area             # g dry per m2 of land per year
    target = dict(P=[rivers['phosphate_p'][0], rivers['dissolved_p'][1]], K=mare['K'], Mg=mare['Mg'], Ca=mare['Ca'])
    reach = {e: dict(returned_g_m2_yr=deposit * pct / 100, export_g_m2_yr=target[e],
                     share_of_export=[deposit * pct / 100 / x for x in target[e]][::-1])
             for e, pct in INSECT_PCT.items()}
    # Air: dimethyl sulfide from seas as active as Earth's ocean, per m2 of lunar land, against the land's sulfur export.
    dms = [x * 1e12 / EARTH_OCEAN_M2 * sea_area / land_area for x in OCEAN_DMS_TG_S_YR]     # g S per m2 of land per year
    air = dict(dms_sulfur_g_m2_yr=dms, land_sulfur_export_g_m2_yr=mare['S'],
               share_of_export=[min(dms) / mare['S'][1], max(dms) / mare['S'][0]])
    # Engineered: the soluble cations (Mg, Na, K; calcium assumed to settle as carbonate) salt works would take out
    # each year to hold the seas' composition, against the world's salt output.
    soluble = [sum(soils['mare_soil']['sea_kg_yr'][e][i] for e in ('Mg', 'Na', 'K')) for i in range(2)]
    engineered = dict(soluble_cations_kg_yr=soluble, world_salt_kg_yr=WORLD_SALT_KG_YR,
                      ratio_to_world_salt=[s / WORLD_SALT_KG_YR for s in soluble],
                      kreep_maria_flooded={m['name']: m['flooded'] for m in atlas['mare_flooding']
                                           if m['name'] in ('Oceanus Procellarum', 'Mare Imbrium')})
    per_ha = lambda x: x / 1e4
    earth = dict(
        seabird_colonies=EARTH_P['seabird_colonies_gg_yr'] * 1e9 / EARTH_LAND_M2,
        uk_high_flying_insects=EARTH_P['uk_insect_migration'][0] * 1e3 / (EARTH_P['uk_insect_migration'][1] * 1e6),
        myvatn_shore_within_50_m=per_ha(EARTH_P['myvatn_shore_kg_ha_yr'] * 1e3),
        amazon_dust=[per_ha(x) for x in EARTH_P['amazon_dust_g_ha_yr']],
        amazon_hydrological_loss=[per_ha(x) for x in EARTH_P['amazon_hydrological_loss_g_ha_yr']],
        earth_rivers_dissolved=EARTH_P['rivers_dissolved_tg_yr'] * 1e12 / EARTH_LAND_M2,
        earth_rivers_particulate=[x * 1e12 / EARTH_LAND_M2 for x in EARTH_P['rivers_particulate_tg_yr']])
    # The migrant conveyor sized for dissolved phosphorus: bodies left on land, the seas' emergence, trips and sugar.
    wild = fruit['results']['near']['0']['earth_like']['day_fruit']['best']['fed']['shares']['0.3']
    fruit_sugar = wild['harvest_kg'] * 10 * FRUIT_SUGAR_PER_100G * CYCLES_PER_YEAR     # g per m2 per year
    moth = fly['fliers']['moth_0.33g']
    moth_dry = FLIERS['moth_0.33g'][0] * 1e3 * INSECT_DRY_SHARE                          # g
    trip_km = 2 * dist['sea']['p50']
    load = target['P']
    dry = {kind: [x / (pct / 100) for x in load] for kind, pct in BODY_P_PCT.items()}
    deaths = [d / moth_dry for d in dry['insect']]
    sugar = [n / share * trip_km * moth['sugar_g_per_km']['moon'] for n in deaths for share in DEATHS_PER_VISIT]
    conveyor = dict(target_dissolved_p_g_m2_yr=load, dry_mass_left_on_land_g_m2_yr=dry,
                    sea_emergence_needed_g_dry_m2_yr=[d * land_area / sea_area for d in dry['insect']],
                    myvatn_emergence_g_m2_yr=EARTH_P['myvatn_emergence_g_m2_yr'],
                    round_trip_km=trip_km, one_way_hours={name: [dist['sea']['p50'] * 1e3 / v / 3600
                                                                 for v in f['airspeed_m_s']['moon']]
                                                          for name, f in fly['fliers'].items()},
                    moth_deaths_per_m2_yr=deaths, deaths_per_visit=DEATHS_PER_VISIT,
                    moth_flight_sugar_g_m2_yr=[min(sugar), max(sugar)], fruit_sugar_g_m2_yr_at_30pct=fruit_sugar,
                    share_of_fruit_sugar=[min(sugar) / fruit_sugar, max(sugar) / fruit_sugar])
    return dict(runoff_m_per_yr_over_land=runoff_m, land_area_m2=land_area, sea_area_m2=sea_area, sea_water_kg=sea_kg,
                land_mean_c=land_c, weathering_co2_mol_m2_yr=co2, rock_stock_top_metre=stock, soils=soils,
                earth_river_load_at_moon_runoff_g_m2_yr=rivers, biology_at_myvatn_rate=dict(
                    deposit_g_dry_m2_land_yr=deposit, insect_pct=INSECT_PCT, elements=reach),
                air=air, engineered=engineered, earth_phosphorus_returns_g_m2_yr=earth, conveyor=conveyor)


def lures(plant):
    """Chemical power of a glowing lure, per lure and per m2 of ground at one lure per reach circle, against the plant's
    night upkeep."""
    photon_j = PLANCK * SPEED_OF_LIGHT / (LURE_NM * 1e-9)
    upkeep = night_upkeep(plant)
    rows = []
    for area in LURE_AREAS_CM2:
        photons = area * LURE_FLUX_PER_CM2_S
        for qy in LURE_QUANTUM_YIELD:
            watts = [photons * photon_j / qy * f for f in LURE_SUBSTRATE_FACTOR]
            per_m2 = [w / (math.pi * r ** 2) for w in watts for r in LURE_REACH_M]
            rows.append(dict(area_cm2=area, photons_s=photons, quantum_yield=qy, lure_w=watts,
                             ground_w_m2=[min(per_m2), max(per_m2)],
                             share_of_night_upkeep_one_lure_per_m2={k: [w / v for w in watts] for k, v in upkeep.items()},
                             share_of_night_upkeep_at_reach={k: [min(per_m2) / v, max(per_m2) / v]
                                                             for k, v in upkeep.items()}))
    return dict(photon_energy_j=photon_j, night_upkeep_w_m2=upkeep, lures=rows)


EVIDENCE = ('Small calculations for an ideas register. The plant and fruit numbers are the biosphere models\' clear-sky '
            'carbon scenarios at the chosen climate; distances come from the hydrostatic atlas at 28% water; runoff from '
            'the first drainage estimate. Weathering assumes congruent dissolution of average lunar soils at the Dessert '
            'basalt rate (1 and 10 times for glassy regolith), and the sea timescales assume each cation stays dissolved; '
            'real weathering is incongruent and sea chemistry would precipitate carbonate, silica and clays. The return '
            'paths are scaled from Earth rates (Lake Myvatn\'s insect emergence, the ocean\'s dimethyl sulfide, river '
            'chemistry) and flight uses fixed lift-to-drag ratios and muscle efficiencies. The conveyor and lure figures '
            'are scenarios. None of this models an ecosystem, a population, a nutrient cycle or sea chemistry.')
READING_RULE = ('Carbon in g C per m2 of ground; elements in g per m2 of land per year unless stated. c4 gives how much '
                'warmer the Moon\'s CO2 puts the C3-C4 crossover. leaves compares the growth needed to replace the canopy '
                'at each leaf lifespan with what the plant keeps at each fruit share (equator, near side, day fruit). '
                'night gives land temperatures, the air column, the share of it the night\'s respiration plus the decay '
                'of a whole dusk fruit fall would use or add, the speed of dusk and the dusk fruit fall as food. flight '
                'gives Moon-to-Earth ratios for the same animal, sugar per kilometre and airspeeds for two fliers. '
                'distances are area-weighted percentiles over land. return_limb.soils gives, per soil, the rock weathered '
                'and each element carried off (pairs: 1 and 10 times the basalt rate), the same per year into the seas '
                '(kg), and the years for the seas to reach Earth sea water\'s concentration of each cation if it all '
                'stayed dissolved. biology_at_myvatn_rate gives what insects emerging at Myvatn\'s highest rate and '
                'landing on land would bring back of each element and its share of the export; air the sulfur that '
                'dimethyl sulfide from seas as active as Earth\'s ocean would bring; engineered the soluble cations salt '
                'works would take out to hold the seas steady, against the world\'s salt output; conveyor the migrant '
                'cycle sized for dissolved phosphorus. lures gives the chemical power of a lure (photon energy over '
                'quantum yield, times 1 or 10 for making the luciferin), alone and over the circle it draws insects from, '
                'against the plant\'s night upkeep.')


def results():
    plant, fruit, atlas, drainage = read(PLANT), read(FRUIT), read(ATLAS), read(DRAINAGE)
    dist = distances(atlas)
    fly = flight()
    return dict(
        schema=SCHEMA,
        producer=dict(files={'research/studies/lunar_cycle_ecology/run.py': digest(Path(__file__))},
                      inputs={str(p.relative_to(ROOT)): digest(p) for p in (PLANT, FRUIT, ATLAS, ATLAS_GRID, DRAINAGE)}),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        settings=dict(cycles_per_year=CYCLES_PER_YEAR, design_pressure_pa=DESIGN_PRESSURE_PA, moon_air_k=MOON_AIR_K,
                      leaf_lifespan_months=LEAF_LIFESPAN_MONTHS, fliers=FLIERS, fruit_sugar_per_100g=FRUIT_SUGAR_PER_100G,
                      river_mg_l=RIVER_MG_L, rock_wt_pct=ROCK_WT_PCT, soil=dict(bulk_kg_m3=SOIL_BULK_KG_M3,
                      depth_m=SOIL_DEPTH_M), body_p_pct=BODY_P_PCT, insect_pct=INSECT_PCT, earth_p=EARTH_P,
                      earth_seawater_g_kg=EARTH_SEAWATER_G_KG, ocean_dms_tg_s_yr=OCEAN_DMS_TG_S_YR,
                      world_salt_kg_yr=WORLD_SALT_KG_YR,
                      airspeed_earth_m_s=AIRSPEED_EARTH_M_S, insect_dry_share=INSECT_DRY_SHARE,
                      lure_nm=LURE_NM, lure_quantum_yield=LURE_QUANTUM_YIELD,
                      lure_flux_per_cm2_s=LURE_FLUX_PER_CM2_S, lure_areas_cm2=LURE_AREAS_CM2,
                      lure_substrate_factor=LURE_SUBSTRATE_FACTOR, lure_reach_m=LURE_REACH_M),
        c4=c4_crossover_shift(), leaves=leaves(fruit), night=night(plant, fruit), flight=fly, distances=dist,
        return_limb=return_limb(plant, fruit, drainage, atlas, dist, fly), lures=lures(plant))


def _round(obj, digits=4):
    if isinstance(obj, float):
        return float(f'{obj:.{digits}g}')
    if isinstance(obj, dict):
        return {str(k): _round(v, digits) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_round(v, digits) for v in obj]
    return obj


def main():
    out = HERE / 'results' / 'lunar_cycle_ecology.json'
    out.parent.mkdir(exist_ok=True)
    out.write_text(json.dumps(_round(results()), indent=1) + '\n')
    print('wrote', out)


if __name__ == '__main__':
    main()
