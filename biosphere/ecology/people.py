"""What the Open Moon's sky, land and living systems can hold for an affluent population, over the near horizons.

    python -m biosphere.ecology.people
    # -> biosphere/ecology/results/people.json

A screen on committed products with literature constants, each named beside it (people_sources.json holds the full
citations and how far each was read). It computes, in seconds:
- **The sky as habitat.** The air by height band from sea level to 90 km: volume above the ground and mass (the
  design run's air, climate/results/gcm/global_winds), temperature, cloud by day and night (the equatorial CM1 ring),
  cosmic-ray ionization against Earth's sea level (atmosphere/electricity), light (the column above, the Sun on the
  horizon, how much later the Sun sets aloft), wind and hydrogen lift.
- **Lift.** Earth's largest flyers carried to the Moon by the rule of similarity (research/studies/sky_ships), and what
  lunar gravity and the dense air do for the same animal.
- **Water and nutrients aloft.** The air's water and its turnover (the CM1 rings), and the primary production that
  phosphorus falling from the air could feed aloft, against the land's growth (biosphere/canopy/results/plant.json).
- **People in the sky.** How often a flyer passes near a point in each traffic band of the summit metropolis at the
  busy hour (research/studies/sky_fleet).
- **Land.** The Moon's land by kind (geography's atlas and drainage, the aerosol study's regions) and by latitude, with
  its rain and the share of the crop model's growth that rain supports (the design run's rain,
  climate/results/gcm/convection).
- **What a person takes.** Cropland and grazing by diet from the day fruit (biosphere/canopy/results/fruit.json),
  settlement land, nitrogen, phosphorus, methane, N2O, heat and lit sky, for Earth's practice and for a lunar design.
- **Ethylene.** What the vegetated land's ethylene settles at when soils are its only sink.
- **Capacity.** For 0.1 to 10 billion people, and the population at which each limit binds, by practice and by diet.
  Methane and N2O use this branch's first screen (results/first_screen.json): soils are methane's only sink and N2O
  has none.
- **Horizons.** The 500-year build, the first thousand years and several thousand: what accumulates and how fast.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

from shared.constants import (MOON_GM, MOON_RADIUS, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY, SYNODIC_MONTH_DAYS,
                              JULIAN_YEAR_DAYS)
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCHEMA = 'terluna.biosphere.ecology-people/1'
OUT = HERE / 'results' / 'people.json'

INPUTS = dict(
    atlas=(ROOT / 'geography' / 'results' / 'atlas.json', 'terluna.geography.atlas/1'),
    drainage=(ROOT / 'geography' / 'results' / 'drainage.json', 'terluna.geography.drainage/1'),
    air=(ROOT / 'climate' / 'results' / 'gcm' / 'global_winds_A28_dim5_moon.json', 'terluna.climate.gcm-global-winds/1'),
    rain=(ROOT / 'climate' / 'results' / 'gcm' / 'convection_A28_dim5_moon.json', 'terluna.climate.gcm-convection/1'),
    ring=(ROOT / 'climate' / 'results' / 'crm' / 'ring_ring_equator.json', 'terluna.climate.crm-ring/1'),
    ring_45n=(ROOT / 'climate' / 'results' / 'crm' / 'ring_ring_45n_lsw.json', 'terluna.climate.crm-ring/1'),
    ring_70=(ROOT / 'climate' / 'results' / 'crm' / 'ring_ring_70_45e.json', 'terluna.climate.crm-ring/1'),
    ring_80n=(ROOT / 'climate' / 'results' / 'crm' / 'ring_ring_80n_lsw.json', 'terluna.climate.crm-ring/1'),
    muons_moon=(ROOT / 'atmosphere' / 'electricity' / 'results' / 'muon_ionization_moon.json',
                'terluna.atmosphere.muon-ionization/1'),
    muons_earth=(ROOT / 'atmosphere' / 'electricity' / 'results' / 'muon_ionization_earth.json',
                 'terluna.atmosphere.muon-ionization/1'),
    canopy=(ROOT / 'biosphere' / 'canopy' / 'results' / 'canopy.json', 'terluna.biosphere.canopy-photosynthesis/1'),
    plant=(ROOT / 'biosphere' / 'canopy' / 'results' / 'plant.json', 'terluna.biosphere.plant-carbon-cycle/1'),
    fruit=(ROOT / 'biosphere' / 'canopy' / 'results' / 'fruit.json', 'terluna.biosphere.fruit-carbon/2'),
    screen=(HERE / 'results' / 'first_screen.json', 'terluna.biosphere.ecology-first-screen/1'),
    aerosol=(ROOT / 'research' / 'studies' / 'open_moon_aerosol' / 'results' / 'open_moon_aerosol.json',
             'terluna.research.open-moon-aerosol/1'),
    fleet=(ROOT / 'research' / 'studies' / 'sky_fleet' / 'results' / 'sky_fleet.json', None),
    reference_air=(ROOT / 'research' / 'baselines' / 'feasibility' / 'reference.json', None),
)

# ---- The sky -------------------------------------------------------------------------------------------------------
# Height bands above sea level, set by what lives and flies there (research/studies/sky_fleet, the CM1 ring).
BANDS_KM = (
    ('low_sky', 0.0, 10.0),        # the day's mixed layer (to 9.9 km), cloud base 4-7 km; people's wings, sky boats, ferries
    ('storm_layer', 10.0, 22.0),   # storms to their median top, 22 km
    ('cool_air', 22.0, 35.0),      # regional ships at 22-27 km; the air freezes near 26 km
    ('flight_band', 35.0, 45.0),   # long-haul liners and freighters; the port's crown at 35.6 km
    ('high_air', 45.0, 90.0),      # winged liners near 55 km, high platforms near 70 km
)
RAYLEIGH_TAU_550 = 0.71            # molecular optical depth of the design column at 550 nm (joint synthesis, haze)
EARTH_SEA_LEVEL_PA = 101325.0      # U.S. Standard Atmosphere 1976
EARTH_SEA_LEVEL_DENSITY = 1.225    # kg/m3, U.S. Standard Atmosphere 1976
EARTH_AIR_MASS_KG = 5.148e18       # Trenberth & Smith 2005, total mean mass of Earth's atmosphere
EARTH_RADIUS_KM = 6371.0           # for Earth's surface area
H2_MOLAR_MASS = 2.016e-3           # kg/mol
GAS_R = 8.314462618                # J/mol/K
NEAR_PASS_M = 20.0                 # a flyer passing within this distance of a point counts as a near pass (this screen)

# ---- Lift ----------------------------------------------------------------------------------------------------------
EARTH_FLYERS = {                   # span m (low, high), mass kg (low, high)
    'wandering albatross (Diomedea exulans)': ((3.0, 3.1), (7.8, 9.4)),         # Shaffer, Weimerskirch & Costa 2001
    'Argentavis magnificens (Miocene)': ((7.0, 7.0), (70.0, 72.0)),             # Chatterjee, Templin & Campbell 2007
    'Quetzalcoatlus northropi (Cretaceous)': ((10.0, 11.0), (200.0, 250.0)),    # Witton & Habib 2010
}

# ---- Water and nutrients aloft -------------------------------------------------------------------------------------
EARTH_WATER_VAPOUR_KG = 1.27e16          # mean mass of water vapour in Earth's air (Trenberth & Smith 2005)
EARTH_VAPOUR_RESIDENCE_DAYS = (8.5, 9.3)  # 8.9 +- 0.4 days (van der Ent & Tuinenburg 2017)
DUST_P_G_M2_YR = (0.0007, 0.028)         # phosphorus falling from the air: African dust on the Amazon, 7 g/ha a year
                                         # at least (Yu et al. 2015), to the median measured over Earth's land
                                         # (Mahowald et al. 2008, via Tipping et al. 2014)
C_TO_P_MASS = 106 * 12.011 / 30.974      # Redfield ratio, 106 C : 1 P by atoms (Redfield 1958)

# ---- Land ----------------------------------------------------------------------------------------------------------
LAT_BANDS_DEG = ((0, 15), (15, 30), (30, 45), (45, 60), (60, 75), (75, 90))
DRY_MM_DAY = 0.5                         # the GCM README's dry land
WUE_EARTH_CROPS_G_C_PER_KG = (1.05, 2.36)  # gross production per evapotranspiration at crop sites, annual means
                                         # (Tang et al. 2014); scaled by the design air's CO2 pressure against Earth's
                                         # (canopy.json), since a leaf's water-use efficiency is (Ca - Ci) / 1.6 VPD
EARTH_LAND_KM2 = 149e6                   # Earth's land (Ritchie & Roser 2019, from FAO)
EARTH_HABITABLE_KM2 = 104e6              # ice-free, non-barren land (same)
EARTH_AGRICULTURE_KM2 = 51e6             # agricultural land (same)

# ---- Food ----------------------------------------------------------------------------------------------------------
FRUIT_KCAL_PER_G = 0.30                  # the day fruit has today's watermelon's make-up (USDA, 30 kcal per 100 g)
INTAKE_KCAL_DAY = 2500.0                 # food energy a day, the EAT-Lancet basis (Willett et al. 2019)
CLOUD_SHARE = 0.89                       # the GCM's ground gets 89% of clear-sky sunlight (provisioning first screen)
LOSS_FACTOR = (1.3, 1.5)                 # production per food eaten: a third is lost or wasted on Earth (FAO 2011);
                                         # the low end a designed food chain
VARIETY_FACTOR = (1.5, 2.5)              # a varied plant diet with protein and fats against the energy-only fruit
                                         # (the provisioning screen's planning factor)
SHARES = ('0.5', '0.6')                  # fruit's share of new tissue on farms (the ecology register, D3)
WORLD_2010_PEOPLE = 7.0e9                # World Bank, SP.POP.TOTL
DIETS = {                                # cropland against an all-plant diet's, and grazing land a person (m2)
    'plant_based': dict(crop=(1.0, 1.0), grazing_m2=(0.0, 0.0)),            # vegan, 0.13 ha (Peters et al. 2016)
    'plant_rich': dict(crop=(1.08, 1.15), grazing_m2=(0.0, 1000.0)),        # lacto-ovo vegetarian, 0.14 ha, to the
                                                                            # 20% omnivore, 0.15 ha + 0.10 grazed (same)
    'earth_average': dict(crop=(1.24, 1.24), grazing_m2=(2.89e13 / WORLD_2010_PEOPLE,) * 2),
                                         # today's diet: 1,240 Mha arable and 2,890 Mha pasture against 1,000 Mha for
                                         # no animal products (Poore & Nemecek 2018)
    'earth_affluent': dict(crop=(0.34 / 0.13,) * 2, grazing_m2=(7400.0, 7400.0)),   # the US baseline: 0.34 ha of
                                         # cropland and 0.74 ha grazed against 0.13 ha vegan (Peters et al. 2016)
    'earth_affluent_fed_on_crops': dict(crop=(0.7 + 0.3 / 0.08, 0.7 + 0.3 / 0.07), grazing_m2=(0.0, 0.0)),
                                         # 30% of food energy from animals (Alexander et al. 2016) at a 7-8% feed
                                         # efficiency (Shepon et al. 2016)
}

# ---- Settlement, kept land -----------------------------------------------------------------------------------------
SETTLEMENT_M2 = {
    'metropolis': 25.0,                  # 40,000 people per km2 (decisions register, the summit metropolis)
    'town': 100.0,                       # 10,000 people per km2, a compact town (a scenario)
    'affluent_city': 1e4 / 70.0,         # 70 people per ha, developed countries' cities (Angel et al. 2010)
    'affluent_spread': 1e4 / 28.0,       # 28 people per ha, land-rich developed countries' cities (same)
}
KEPT_SHARES = {'30x30': 0.30, 'half': 0.50}  # Kunming-Montreal GBF Target 3 (CBD 2022); Half-Earth (Wilson 2016)

# ---- Nutrients -----------------------------------------------------------------------------------------------------
N_EATEN_KG_YR = 5.2                      # nitrogen eaten a person a year in the US (Leach et al. 2012)
FOOD_N_FOOTPRINT_KG_YR = (22.0, 30.0)    # N lost to the environment in feeding a person, the Netherlands to the US
                                         # (Leach et al. 2012): Earth's practice's new N
DESIGN_CHAIN_NUE = (0.40, 0.60)          # N eaten per N used for grains, starchy roots and legumes: virtual N factors
                                         # 1.4, 1.5 and 0.7 (Leach et al. 2012)
DESIGN_SEWAGE_N_RETURN = (0.5, 0.8)      # share of excreted N back on fields (a design target)
N2O_EF = (0.002, 0.005)                  # N2O-N per N applied: the low end of EF1, 0.2%, to dry climates, 0.5%
                                         # (IPCC 2019 Refinement, Table 11.1); the lunar design's fields
P_EATEN_KG_YR = (0.44, 0.60)             # excreted P, the average diet to meat eaters (Cordell et al. 2009)
P_LOST_PER_P_EATEN = {'earth_practice': (4.0, 4.0), 'lunar_design': (0.1, 0.3)}  # Earth: a fifth of mined P reaches
                                         # food (Cordell et al. 2009); lunar design: losses held to 10-30% (a target)
LAND_P_LOSS_G_M2_YR = (0.063, 0.178)     # land's natural P loss, dissolved 0.003-0.008 plus particulate 0.06-0.17
                                         # (research/studies/lunar_cycle_ecology, Earth rates at the Moon's runoff)
MO_G_PER_HA = (15.0, 120.0)              # molybdenum dressings, wheat to soybean (Kovacs et al. 2026)
MO_WORLD_T_YR = 2.6e5                    # world molybdenum mine production, 2024 (USGS 2025)
FIXATION_KWH_PER_KG_N = (36e3 / 3600 * 17.031 / 14.007,) * 2   # electrolytic ammonia, 36 GJ of electricity per
                                         # tonne of NH3 (IEA 2021), per kg of N

# ---- Gases ---------------------------------------------------------------------------------------------------------
CH4_KG_PERSON_YR = {'earth_practice': (206e9 / 7.22e9, 217e9 / 7.22e9), 'lunar_design': (0.3, 3.0)}
#   Earth: agriculture and waste, 206 (bottom-up) and 217 (top-down) Tg a year over 2008-2017 (Saunois et al. 2020),
#   over the decade's 7.22 billion people (World Bank); lunar design: no flooded paddies, no landfills, sealed
#   digesters, ruminants absent or inhibited (a design target)
N2O_EARTH_KG_N_PERSON_YR = (5.1e9 / 7.13e9, 7.3e9 / 7.13e9)
#   agriculture, direct 3.8 and indirect 1.3, to all anthropogenic N2O, 7.3 Tg N a year over 2007-2016
#   (Tian et al. 2020), over the decade's 7.13 billion people (World Bank)
ETHYLENE_EARTH_LAND_TG_YR = (11.9, 29.6)  # natural land sources: 18-45 Tg a year in all, 74% natural, 89% of that
                                          # from land (Sawada & Totsuka 1986)
ETHYLENE_SOIL_M_S = (1 / (9.82 * 3600), 10 / (9.82 * 3600))  # soil uptake through 9.82 h/m of resistance, a model
                                          # estimate (Sawada et al. 1986), and ten times it in engineered soils
                                          # (a design target)
ETHYLENE_HAZARD_PPB = 50.0               # held all season, wheat loses 36% of its yield, rice 63% (Klassen & Bugbee 2002)

# ---- Heat and light ------------------------------------------------------------------------------------------------
HEAT_KW_PERSON = {'earth_practice': (9.3, 9.3),    # US primary energy per person, 279 million Btu in 2023 (EIA)
                  'lunar_design': (4.0, 6.0)}      # 2 kW of final energy (the metropolis brief; needs stop rising near
                                                   # 60 GJ a year, Vogel et al. 2021) from fusion at 33-50% conversion
W_M2_PER_K = (1.0, 1.6)                  # the GCM's climate response (research/studies/atmospheric_co2)
K_PER_PERCENT_DIMMING = 1.9              # the GCM's response to the shield's dimming (climate/gcm)
LIT_LAND_SHARE = {'world': 0.23, 'eu': 0.88, 'us': (0.47, 0.50)}   # land under light-polluted skies: between 75 N and
                                         # 60 S, Europe, "almost half" of the US (Falchi et al. 2016)
LAND_AND_PEOPLE_2015 = {'world': (129.9e6, 7.44e9), 'eu': (3.99e6, 4.44e8), 'us': (9.147e6, 3.22e8)}  # km2, people
                                         # (World Bank, AG.LND.TOTL.K2 and SP.POP.TOTL)
LIT_DESIGN_FACTOR = 0.1                  # shielded, warm, dimmed and curfewed light: a tenth (a design target)
DIET_FOR = {'earth_practice': 'earth_affluent', 'lunar_design': 'plant_rich'}
SETTLEMENT_FOR = {'earth_practice': ('affluent_city', 'affluent_spread'), 'lunar_design': ('metropolis', 'town')}

POPULATIONS = (0.1e9, 0.3e9, 1e9, 3e9, 10e9)
HORIZON_YEARS = (500.0, 1000.0, 5000.0)
GROWTH_PER_YEAR = (0.005, 0.01)          # population growth rates for the horizon arithmetic (scenarios)
SOIL_MM_YR = (0.017, 0.036)              # natural soil formation (Montgomery 2007)

EVIDENCE = ('A screen on committed products: the design run\'s air, rain and land, the CM1 rings\' cloud and water, '
            'the muon ionization column, the canopy, plant and fruit models, the atlas and drainage, the aerosol '
            'study\'s regions, the sky fleet\'s busy hour and this branch\'s first screen of methane and N2O. The '
            'per-person needs are Earth figures from the literature or labelled design targets; the capacities follow '
            'by arithmetic. No ecosystem, soil, water, chemistry or traffic model is run, ethylene\'s soil uptake rests '
            'on one estimated soil resistance, and the greenhouse effect of the gases is not computed.')
READING_RULE = ('sky.bands give each height band above sea level its volume above the ground (km3), mass (kg), '
                'temperature (C), cloud fraction (CM1 equatorial ring, day and night), cosmic-ray ionization as a '
                'ratio to Earth\'s sea level, the direct sunlight with the Sun overhead and on the horizon, and the '
                'hours the Sun stays up longer at each end of the equatorial day. Pairs are low-high ranges. Land '
                'shares are of the whole Moon; per-person needs are per year; capacity gives each limit\'s '
                'population in people (1e15 marks no limit, a diet without grazing); land_ is the least of crops_ '
                '(food land at the model\'s yields against the rain-fed potential land), grazing_ (on the fog '
                'desert) and area_ (all of it on the dry land). Methane ppm are steady states over this branch\'s '
                'first screen\'s soil sink; N2O accumulates without a sink; ethylene settles where the land\'s '
                'source meets the soils\' uptake.')


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def sig(x, n=3):
    """Round to n significant figures, keeping floats plain."""
    if x == 0 or not np.isfinite(x):
        return float(x)
    return float(f'{x:.{n}g}')


def pair(lo, hi, n=3):
    return [sig(min(lo, hi), n), sig(max(lo, hi), n)]


def moon_area_m2():
    return 4 * math.pi * MOON_RADIUS ** 2


def gravity(height_m):
    return MOON_GM / (MOON_RADIUS + height_m) ** 2


# ---- The sky -------------------------------------------------------------------------------------------------------
def sky(air, ring, rings, muons_moon, muons_earth, molar_mass):
    h = np.array(air['height_km']) * 1e3
    cover, rho = np.array(air['coverage']), np.array(air['density_mean_kg_m3'])
    p, t = np.array(air['pressure_pa']), np.array(air['temperature_k'])
    shell = lambda z: 4 * np.pi * (MOON_RADIUS + z) ** 2

    def integral(a, b, f):
        z = np.linspace(a, b, 401)
        return float(np.trapezoid(f(z), z))

    volume = lambda a, b: integral(a, b, lambda z: np.interp(z, h, cover) * shell(z))
    mass = lambda a, b: integral(a, b, lambda z: np.interp(z, h, cover * rho) * shell(z))
    above_90 = float(p[-1] * shell(h[-1]) / gravity(h[-1]))
    total = mass(h[0], h[-1]) + above_90

    # Ionization by muons and their decay products, per gram of air (proportional to dose), against Earth's sea level.
    zm = np.array([r['height_km'] for r in muons_moon['rows']])[::-1]
    im = np.array([r['ion_pairs_per_g_s']['total'] for r in muons_moon['rows']])[::-1]
    earth_sea = muons_earth['rows'][-1]['ion_pairs_per_g_s']['total']
    ionization = lambda z_km: float(np.exp(np.interp(z_km, zm, np.log(im))) / earth_sea)

    # Light. The Rayleigh depth above a height scales with the pressure there. The Sun on an observer's horizon shines
    # through the column above along a grazing path, Chapman's sqrt(pi x / 2) times the vertical, x = (R + z) / H.
    sun_deg_per_hour = 360.0 / (SYNODIC_MONTH_DAYS * 24.0)
    p_ground = float(p[0])

    def light(z):
        tz, pz = float(np.interp(z, h, t)), float(np.interp(z, h, p))
        tau = RAYLEIGH_TAU_550 * pz / p_ground
        scale = GAS_R * tz / (molar_mass * gravity(z))
        grazing = math.sqrt(math.pi * (MOON_RADIUS + z) / (2 * scale))
        dip = math.degrees(math.acos(MOON_RADIUS / (MOON_RADIUS + z)))
        return dict(rayleigh_depth_above=sig(tau), direct_sun_overhead=sig(math.exp(-tau)),
                    direct_sun_on_horizon=sig(math.exp(-tau * grazing)), horizon_dip_deg=sig(dip),
                    sun_stays_up_longer_h_each_end_equator=sig(dip / sun_deg_per_hour))

    c = ring['circulation']
    ring_h = np.array(c['heights_km'], dtype=float)
    ha = np.array(c['hour_angle_deg'])
    cloud = np.array(c['cloud_fraction'], dtype=float)
    day = np.abs(ha) < 90.0
    cloud_day, cloud_night = cloud[:, day].mean(axis=1), cloud[:, ~day].mean(axis=1)
    speed = np.array(air['speed_p50_m_s'])
    h2_lift = rho * (1 - H2_MOLAR_MASS / molar_mass)            # kg per m3 of pure hydrogen
    freezing_km = float(np.interp(273.15, t[::-1], h[::-1]) / 1e3)
    warm_volume = volume(h[0], freezing_km * 1e3)

    bands = {}
    for name, a_km, b_km in BANDS_KM:
        a, b = a_km * 1e3, b_km * 1e3
        inside = (ring_h >= a_km) & (ring_h <= b_km)
        bands[name] = dict(
            bottom_km=a_km, top_km=b_km,
            volume_km3=sig(volume(a, b) / 1e9), mass_kg=sig(mass(a, b)), mass_share=sig(mass(a, b) / total),
            temperature_c=[round(float(np.interp(z, h, t)) - 273.15, 1) for z in (a, b)],
            density_kg_m3=[sig(float(np.interp(z, h, rho))) for z in (a, b)],
            column_above_bottom_t_m2=sig(float(np.interp(a, h, p)) / gravity(a) / 1e3),
            column_above_bottom_earths=sig(float(np.interp(a, h, p)) / gravity(a) / (EARTH_SEA_LEVEL_PA / STANDARD_GRAVITY)),
            ionization_mid_vs_earth_sea_level=sig(ionization((a_km + b_km) / 2)),
            light_at_bottom=light(a), light_at_top=light(b),
            wind_median_mid_m_s=sig(float(np.interp((a + b) / 2, h, speed))),
            hydrogen_lift_kg_m3=[sig(float(np.interp(z, h, h2_lift))) for z in (a, b)],
            cloud_fraction_ring=dict(heights_km=ring_h[inside].tolist(),
                                     day=[sig(x, 2) for x in cloud_day[inside]],
                                     night=[sig(x, 2) for x in cloud_night[inside]]))
    ring_land = ring['by_hour_angle']['land']
    lha = np.array(ring['by_hour_angle']['hour_angle_deg'])
    lday = np.abs(lha) < 90.0
    moon_area = moon_area_m2()
    earth_area = 4 * math.pi * (EARTH_RADIUS_KM * 1e3) ** 2
    return dict(
        bands=bands, freezing_height_km=sig(freezing_km), warm_air_volume_km3=sig(warm_volume / 1e9),
        warm_air_mean_depth_above_ground_km=sig(warm_volume / moon_area / 1e3),
        air_mass_kg=sig(total), above_90_km_mass_share=sig(above_90 / total),
        air_mass_vs_earth=sig(total / EARTH_AIR_MASS_KG), area_vs_earth=sig(moon_area / earth_area),
        column_ground_t_m2=sig(p_ground / MOON_SURFACE_GRAVITY / 1e3),
        column_vs_earth=sig(p_ground / MOON_SURFACE_GRAVITY / (EARTH_SEA_LEVEL_PA / STANDARD_GRAVITY)),
        ionization_ground_vs_earth_sea_level=sig(ionization(0.0)),
        equator_ring_land=dict(
            mixed_layer_km_day_night=[sig(float(np.mean(np.array(ring_land['mixed_layer_km'])[m])))
                                      for m in (lday, ~lday)],
            mixed_layer_km_max=sig(float(np.max(ring_land['mixed_layer_km']))),
            cloud_base_km_day=sig(float(np.mean(np.array(ring_land['cloud_base_lcl_km'])[lday]))),
            storm_top_km=ring['storms']['cloud_top_km'], cloud_cover=ring['cloud_cover']),
        water=dict(
            precipitable_water_mm={k: sig(np.mean(list(r['water_budget']['precipitable_water_mm'].values())))
                                   for k, r in rings.items()},
            vapour_residence_days={k: sig(r['water_budget']['residence_days']) for k, r in rings.items()},
            earth_precipitable_water_mm=sig(EARTH_WATER_VAPOUR_KG / earth_area),
            earth_vapour_residence_days=list(EARTH_VAPOUR_RESIDENCE_DAYS)))


# ---- Lift ----------------------------------------------------------------------------------------------------------
def lift(air, molar_mass):
    k = STANDARD_GRAVITY / MOON_SURFACE_GRAVITY                 # 6.04
    rho_moon = float(air['density_mean_kg_m3'][0])
    rho_ratio = rho_moon / EARTH_SEA_LEVEL_DENSITY
    g_ratio = 1 / k
    same_density_km = float(np.interp(-EARTH_SEA_LEVEL_DENSITY, -np.array(air['density_mean_kg_m3']),
                                      np.array(air['height_km'])))
    flyers = {name: dict(earth_span_m=list(span), earth_mass_kg=list(m),
                         moon_span_m=[sig(s * k) for s in span], moon_mass_t=[sig(x * k ** 3 / 1e3) for x in m])
              for name, (span, m) in EARTH_FLYERS.items()}
    return dict(
        similarity_factor=sig(k), earth_sea_level_density_height_km=sig(same_density_km), flyers=flyers,
        same_animal_at_the_ground=dict(
            weight=sig(g_ratio),
            hover_power=sig(g_ratio ** 1.5 / math.sqrt(rho_ratio)),            # induced power, (m g)^1.5 / sqrt(rho)
            energy_per_km=sig(g_ratio),                                        # at a fixed lift-to-drag ratio
            glide_speed_and_sink=sig(math.sqrt(g_ratio / rho_ratio))),         # at a fixed lift coefficient
        hydrogen_float_skin_kg_m2_per_m_radius=sig((rho_moon * (1 - H2_MOLAR_MASS / molar_mass)) / 3))


# ---- Water and nutrients aloft -------------------------------------------------------------------------------------
def aloft(plant):
    cycles_per_year = JULIAN_YEAR_DAYS / SYNODIC_MONTH_DAYS
    growth = [plant['results']['near']['0'][s]['growth'] * cycles_per_year for s in ('earth_like', 'idle_quarter')]
    production = [p * C_TO_P_MASS for p in DUST_P_G_M2_YR]
    return dict(
        dust_p_g_m2_yr=list(DUST_P_G_M2_YR), c_to_p_mass=sig(C_TO_P_MASS),
        production_on_dust_p_g_c_m2_yr=pair(*production),
        land_growth_g_c_m2_yr=pair(*growth),
        production_share_of_land_growth=pair(production[0] / growth[1], production[1] / growth[0]))


# ---- People in the sky ---------------------------------------------------------------------------------------------
def traffic(fleet):
    bands = fleet['metropolis']['bands']
    people = fleet['metropolis']['city']['population']
    out = {}
    for name, b in bands.items():
        per_m3 = b['per_km3'] / 1e9
        out[name] = dict(bottom_m=b['bottom_m'], top_m=b['top_m'], aloft=b['aloft'], speed_m_s=b['speed_m_s'],
                         per_km3=b['per_km3'],
                         near_passes_per_hour=sig(per_m3 * b['speed_m_s'] * math.pi * NEAR_PASS_M ** 2 * 3600.0),
                         aloft_per_million_people=sig(b['aloft'] / people * 1e6))
    return dict(near_pass_m=NEAR_PASS_M, busy_hour=out, city_area_km2=sig(fleet['metropolis']['city']['area_km2']),
                city_people=people)


# ---- Land ----------------------------------------------------------------------------------------------------------
def land(atlas, drainage, aerosol, rain, plant, canopy):
    area = moon_area_m2()
    seas = atlas['water_share']
    lakes = drainage['lakes']['area_share']
    land_share = atlas['main_land_share'] + atlas['island_share']
    dry = land_share - lakes
    regions = {k: sig(aerosol['regions'][k]['area_share'], 4) for k in ('wet_land', 'fog_desert', 'polar_dry_land')}

    lat, lon = np.array(rain['lat_deg']), np.array(rain['lon_deg'])
    pr, gland, w = np.array(rain['pr_mean']), np.array(rain['land_share']), np.array(rain['cell_area_share'])
    L = gland * w
    scale = dry / float(L.sum())                                # the GCM's T21 land scaled to the atlas's dry land
    runoff_ratio = drainage['rivers']['discharge_to_sea_m3s'] / (float((pr * L).sum()) / 1e3 / 86400.0 * area)
    et_share = 1.0 - runoff_ratio
    near = ((lon <= 90.0) | (lon >= 270.0))[None, :] * np.ones_like(pr, dtype=bool)

    # Earth's crop water-use efficiency, raised by the design air's CO2 pressure at the canopy model's fixed Ci/Ca.
    co2 = canopy['settings']['co2_pa']
    co2_factor = co2['moon'] / co2['earth']
    wue = [x * co2_factor for x in WUE_EARTH_CROPS_G_C_PER_KG]
    # The plant model's gross production over the cycle at 0, 30 and 60 degrees (near side, Earth-like strategy).
    gpp_lat = np.array([0.0, 30.0, 60.0])
    gpp_day = np.array([plant['results']['near'][k]['earth_like']['gpp'] for k in ('0', '30', '60')]) / SYNODIC_MONTH_DAYS
    bands, effective, tropics = [], [0.0, 0.0], [0.0, 0.0]
    for lo, hi in LAT_BANDS_DEG:
        m = ((np.abs(lat) >= lo) & (np.abs(lat) < hi))[:, None] * np.ones_like(pr, dtype=bool)
        a = float((L * m).sum())
        rain_mm = float((pr * L * m).sum() / a)
        gpp = float(np.interp((lo + hi) / 2, gpp_lat, gpp_day))
        supported = [min(1.0, rain_mm * et_share * w / gpp) for w in wue]
        rel = gpp / gpp_day[0]
        for i in (0, 1):
            effective[i] += a * scale * supported[i] * rel
            if hi <= 30:
                tropics[i] += a * scale * supported[i] * rel
        bands.append(dict(abs_lat_deg=[lo, hi], dry_land_share=sig(a * scale),
                          near_side_share_of_band=sig(float((L * m * near).sum()) / a),
                          rain_mm_day=sig(rain_mm), wet_share=sig(float((L * m * (pr >= DRY_MM_DAY)).sum()) / a),
                          gross_production_g_c_m2_day=sig(gpp), rain_supports_share_of_model=pair(*supported, n=2)))
    # The rivers, used whole, irrigate at the equator's shortfall: an upper bound on irrigated potential land.
    eq_need_mm = gpp_day[0] / np.array(wue)
    rivers_m3_yr = drainage['rivers']['discharge_to_sea_m3s'] * JULIAN_YEAR_DAYS * 86400.0
    irrigable = rivers_m3_yr / (eq_need_mm / 1e3 * JULIAN_YEAR_DAYS)            # m2 of potential land
    return dict(
        moon_km2=sig(area / 1e6, 4), seas=seas, lakes=lakes, land=sig(land_share, 4), dry_land=sig(dry, 4),
        dry_land_km2=sig(dry * area / 1e6, 4), regions=regions,
        dry_land_vs_earth=dict(land=sig(dry * area / 1e6 / EARTH_LAND_KM2),
                               habitable=sig(dry * area / 1e6 / EARTH_HABITABLE_KM2),
                               agriculture=sig(dry * area / 1e6 / EARTH_AGRICULTURE_KM2)),
        near_side_dry_land_share=sig(float((L * near).sum()) * scale), runoff_ratio=sig(runoff_ratio),
        by_latitude=bands, wue_earth_crops_g_c_per_kg=list(WUE_EARTH_CROPS_G_C_PER_KG),
        co2_pressure_factor=sig(co2_factor), wue_g_c_per_kg=pair(*wue),
        crop_water_need_mm_day_equator=pair(*eq_need_mm),
        potential_equivalent_share=pair(*effective, n=3),
        potential_within_30_deg_share=pair(*(t / e for t, e in zip(tropics, effective))),
        rivers_irrigate_potential_km2=pair(*(irrigable / 1e6)))


# ---- What a person takes -------------------------------------------------------------------------------------------
def food(fruit):
    per_day = {}
    for side, bands in fruit['results'].items():
        for band in bands:
            for strategy in ('earth_like', 'idle_quarter'):
                shares = fruit['results'][side][band][strategy]['day_fruit']['best']['fed']['shares']
                for s in SHARES:
                    per_day[f'{side}/{band}/{strategy}/{s}'] = shares[s]['harvest_kg'] * 1e3 * FRUIT_KCAL_PER_G / SYNODIC_MONTH_DAYS
    equator = [v for k, v in per_day.items() if k.split('/')[1] == '0']
    potential = [INTAKE_KCAL_DAY / max(equator), INTAKE_KCAL_DAY / min(equator)]
    plant_diet = [potential[0] / CLOUD_SHARE * LOSS_FACTOR[0] * VARIETY_FACTOR[0],
                  potential[1] / CLOUD_SHARE * LOSS_FACTOR[1] * VARIETY_FACTOR[1]]
    diets = {k: dict(crop_m2=pair(plant_diet[0] * d['crop'][0], plant_diet[1] * d['crop'][1]),
                     grazing_m2=pair(*d['grazing_m2']), crop_factor=pair(*d['crop']))
             for k, d in DIETS.items()}
    return dict(kcal_m2_day={k: sig(v) for k, v in per_day.items()}, potential_m2_person=pair(*potential),
                plant_based_m2_person=pair(*plant_diet), diets=diets)


def lit_sky_m2():
    """Land under light-polluted skies per person (Falchi et al. 2016 over the World Bank's land and people)."""
    per = {}
    for k, share in LIT_LAND_SHARE.items():
        km2, people = LAND_AND_PEOPLE_2015[k]
        lo, hi = (share, share) if not isinstance(share, tuple) else share
        per[k] = pair(lo * km2 * 1e6 / people, hi * km2 * 1e6 / people)
    affluent = pair(per['eu'][0], per['us'][1])
    return dict(by_region=per, earth_practice=affluent,
                lunar_design=pair(affluent[0] * LIT_DESIGN_FACTOR, affluent[1] * LIT_DESIGN_FACTOR))


def needs(food_table):
    lit = lit_sky_m2()
    out = {}
    for practice in ('earth_practice', 'lunar_design'):
        diet = food_table['diets'][DIET_FOR[practice]]
        if practice == 'earth_practice':
            n_new = list(FOOD_N_FOOTPRINT_KG_YR)
            n2o = list(N2O_EARTH_KG_N_PERSON_YR)
        else:
            nue, back = DESIGN_CHAIN_NUE, DESIGN_SEWAGE_N_RETURN
            n_new = [N_EATEN_KG_YR / nue[1] - back[1] * N_EATEN_KG_YR, N_EATEN_KG_YR / nue[0] - back[0] * N_EATEN_KG_YR]
            applied = [N_EATEN_KG_YR / nue[1], N_EATEN_KG_YR / nue[0]]       # all N on the fields, new and returned
            n2o = [applied[0] * N2O_EF[0], applied[1] * N2O_EF[1]]
        p_lost = [P_EATEN_KG_YR[0] * P_LOST_PER_P_EATEN[practice][0], P_EATEN_KG_YR[1] * P_LOST_PER_P_EATEN[practice][1]]
        out[practice] = dict(
            diet=DIET_FOR[practice], food_m2=diet['crop_m2'], grazing_m2=diet['grazing_m2'],
            settlement=list(SETTLEMENT_FOR[practice]),
            settlement_m2=pair(*(SETTLEMENT_M2[k] for k in SETTLEMENT_FOR[practice])),
            n_eaten_kg=N_EATEN_KG_YR, n_new_kg=pair(*n_new),
            fixation_w=pair(*(x * k * 1e3 / (JULIAN_YEAR_DAYS * 24) for x, k in zip(n_new, FIXATION_KWH_PER_KG_N))),
            p_eaten_kg=pair(*P_EATEN_KG_YR), p_lost_kg=pair(*p_lost),
            ch4_kg=pair(*CH4_KG_PERSON_YR[practice]), n2o_kg_n=pair(*n2o),
            heat_kw=list(HEAT_KW_PERSON[practice]), lit_sky_m2=lit[practice])
    return dict(practices=out, lit_sky_m2_by_region=lit['by_region'])


def ethylene(land_table, molar_mass, air, air_mol):
    """Ethylene from the land at Earth's natural rate per habitable area, with soils its only sink."""
    area = moon_area_m2()
    dry_m2 = land_table['dry_land'] * area
    source_tg = [x * dry_m2 / (EARTH_HABITABLE_KM2 * 1e6) for x in ETHYLENE_EARTH_LAND_TG_YR]
    n_air = float(air['density_mean_kg_m3'][0]) / molar_mass            # mol/m3 at the ground
    per_ppb = [dry_m2 * v * n_air * 1e-9 for v in ETHYLENE_SOIL_M_S]   # mol/s taken up per ppb in the air
    mol_s = [x * 1e12 / 28.05 / (JULIAN_YEAR_DAYS * 86400) for x in source_tg]
    steady = {'earth_soils': pair(mol_s[0] / per_ppb[0], mol_s[1] / per_ppb[0]),
              'engineered_soils': pair(mol_s[0] / per_ppb[1], mol_s[1] / per_ppb[1])}
    lifetime = [air_mol * 1e-9 / x / (JULIAN_YEAR_DAYS * 86400) for x in per_ppb[::-1]]
    return dict(source_tg_yr=pair(*source_tg), soil_uptake_m_s=pair(*ETHYLENE_SOIL_M_S), steady_ppb=steady,
                lifetime_years=dict(engineered_soils=sig(lifetime[0]), earth_soils=sig(lifetime[1])),
                hazard_ppb=ETHYLENE_HAZARD_PPB)


# ---- Capacity ------------------------------------------------------------------------------------------------------
def diet_limits(land_table, food_table):
    """The population each diet's crops and grazing allow, with 30% or half of the land kept wild."""
    area = moon_area_m2()
    fog_m2 = land_table['regions']['fog_desert'] * area
    eff_m2 = [s * area for s in land_table['potential_equivalent_share']]
    out = {}
    for name, diet in food_table['diets'].items():
        crop, graze = diet['crop_m2'], diet['grazing_m2']
        row = {}
        for kept_name, kept in KEPT_SHARES.items():
            free = 1.0 - kept
            row[f'crops_{kept_name}'] = pair(free * eff_m2[0] / crop[1], free * eff_m2[1] / crop[0])
            if graze[1]:
                row[f'grazing_{kept_name}'] = pair(free * fog_m2 / graze[1], free * fog_m2 / graze[0] if graze[0] else 1e15)
        out[name] = row
    return out


def capacity(land_table, need, screen):
    area = moon_area_m2()
    dry_m2 = land_table['dry_land'] * area
    fog_m2 = land_table['regions']['fog_desert'] * area
    wet_m2 = land_table['regions']['wet_land'] * area
    eff_m2 = [s * area for s in land_table['potential_equivalent_share']]
    methane, nitrous = screen['trace_gases']['methane'], screen['trace_gases']['nitrous_oxide']
    sink = methane['soil_sink_tg_yr_per_ppm']                  # Tg a year per ppm: low, central, high
    natural_ch4 = methane['central']['sources_tg_yr']
    natural_n2o = nitrous['sources_tg_n_yr']
    tg_n_per_ppb = nitrous['tg_n_per_ppb']
    tw_per_k = [w * area / 1e12 for w in W_M2_PER_K]            # TW of heat that warm the Moon 1 K
    land_p_loss_tg = [g * dry_m2 / 1e12 for g in LAND_P_LOSS_G_M2_YR]
    inf = float('inf')

    out = {}
    for practice, n in need['practices'].items():
        food_m2, graze, settle = n['food_m2'], n['grazing_m2'], n['settlement_m2']
        rows = []
        for people in POPULATIONS:
            ch4 = [people * x / 1e9 for x in n['ch4_kg']]          # Tg a year
            n2o = [people * x / 1e9 for x in n['n2o_kg_n']]        # Tg N a year
            heat_tw = [people * x * 1e3 / 1e12 for x in n['heat_kw']]
            warming = [heat_tw[0] / tw_per_k[1], heat_tw[1] / tw_per_k[0]]
            rows.append(dict(
                people=people,
                food_land_share_of_potential_land=pair(people * food_m2[0] / eff_m2[1], people * food_m2[1] / eff_m2[0]),
                grazing_share_of_fog_desert=pair(people * graze[0] / fog_m2, people * graze[1] / fog_m2),
                settlement_share_of_dry_land=pair(people * settle[0] / dry_m2, people * settle[1] / dry_m2),
                city_share_of_dry_land_at_metropolis_density=sig(people * SETTLEMENT_M2['metropolis'] / dry_m2),
                n_new_tg=pair(*(people * x / 1e9 for x in n['n_new_kg'])),
                p_lost_tg=pair(*(people * x / 1e9 for x in n['p_lost_kg'])),
                ch4_tg=pair(*ch4), ch4_added_steady_ppm=pair(ch4[0] / sink[2], ch4[1] / sink[0]),
                ch4_added_steady_ppm_central=pair(ch4[0] / sink[1], ch4[1] / sink[1]),
                n2o_tg_n=pair(*n2o), n2o_ppb_per_year=pair(*(x / tg_n_per_ppb for x in n2o)),
                n2o_ppm_after_years={f'{y:g}': pair(*(x / tg_n_per_ppb * y / 1e3 for x in n2o)) for y in HORIZON_YEARS},
                heat_tw=pair(*heat_tw), heat_w_m2=pair(*(x * 1e12 / area for x in heat_tw)),
                warming_k_unanswered=pair(*warming),
                dimming_percent_to_answer=pair(*(x / K_PER_PERCENT_DIMMING for x in warming)),
                lit_sky_share_of_dry_land=pair(*(people * x / dry_m2 for x in n['lit_sky_m2']))))

        limits = {}
        for kept_name, kept in KEPT_SHARES.items():
            free = 1.0 - kept
            # Crops: food land at the model's yields against the rain-fed potential land.
            by_food = [free * eff_m2[0] / food_m2[1], free * eff_m2[1] / food_m2[0]]
            # Grazing on the fog desert's plains, at Earth's grazed area a person.
            by_grazing = [free * fog_m2 / graze[1] if graze[1] else inf, free * fog_m2 / graze[0] if graze[0] else inf]
            # The land itself: crops on the wet land at its mean water-supported yield, grazing and settlement.
            by_all = [free * dry_m2 / (food_m2[1] * wet_m2 / eff_m2[0] + graze[1] + settle[1]),
                      free * dry_m2 / (food_m2[0] * wet_m2 / eff_m2[1] + graze[0] + settle[0])]
            limits[f'land_{kept_name}'] = pair(min(by_food[0], by_grazing[0], by_all[0]),
                                               min(by_food[1], by_grazing[1], by_all[1]))
            limits[f'crops_{kept_name}'] = pair(*by_food)
            limits[f'area_{kept_name}'] = pair(*by_all)
            if graze[1]:
                limits[f'grazing_{kept_name}'] = pair(by_grazing[0], min(by_grazing[1], 1e15))
            limits[f'lit_sky_{kept_name}'] = pair(free * dry_m2 / n['lit_sky_m2'][1], free * dry_m2 / n['lit_sky_m2'][0])
        limits['heat_1k_unanswered'] = pair(tw_per_k[0] * 1e12 / (n['heat_kw'][1] * 1e3),
                                            tw_per_k[1] * 1e12 / (n['heat_kw'][0] * 1e3))
        limits['methane_parity'] = pair(natural_ch4 * 1e9 / n['ch4_kg'][1], natural_ch4 * 1e9 / n['ch4_kg'][0])
        limits['n2o_parity'] = pair(natural_n2o[0] * 1e9 / n['n2o_kg_n'][1], natural_n2o[1] * 1e9 / n['n2o_kg_n'][0])
        limits['phosphorus_parity'] = pair(land_p_loss_tg[0] * 1e9 / n['p_lost_kg'][1],
                                           land_p_loss_tg[1] * 1e9 / n['p_lost_kg'][0])
        ranked = [k for k in limits if not k.startswith(('crops_', 'grazing_', 'area_'))]
        order = sorted(ranked, key=lambda k: math.sqrt(limits[k][0] * limits[k][1]))
        out[practice] = dict(rows=rows, limits_people=limits, order_by_geometric_mean=order)
    return dict(potential_land_km2=pair(*(x / 1e6 for x in eff_m2)), dry_land_km2=sig(dry_m2 / 1e6),
                fog_desert_km2=sig(fog_m2 / 1e6), tw_per_kelvin=pair(*tw_per_k), natural_ch4_tg=natural_ch4,
                natural_n2o_tg_n=natural_n2o, land_p_loss_tg=pair(*land_p_loss_tg), kept_shares=KEPT_SHARES,
                practices=out)


def mo_dressing(land_table):
    """One dressing of molybdenum over the potential cropland or over all the dry land, against the world's output."""
    dry_ha = land_table['dry_land_km2'] * 100.0
    crop_ha = [x * 100.0 for x in land_table['potential_equivalent_share_km2']]
    return dict(g_per_ha=list(MO_G_PER_HA), world_t_yr=MO_WORLD_T_YR,
                dry_land_t=pair(*(dry_ha * g / 1e6 for g in MO_G_PER_HA)),
                dry_land_years_of_world_output=pair(*(dry_ha * g / 1e6 / MO_WORLD_T_YR for g in MO_G_PER_HA)),
                potential_cropland_t=pair(crop_ha[0] * MO_G_PER_HA[0] / 1e6, crop_ha[1] * MO_G_PER_HA[1] / 1e6))


def horizons(screen):
    methane, nitrous = screen['trace_gases']['methane'], screen['trace_gases']['nitrous_oxide']
    lifetimes = methane['lifetime_years']
    out = {}
    for y in HORIZON_YEARS:
        out[f'{y:g}'] = dict(
            methane_share_of_steady_state=pair(*(1 - math.exp(-y / tau) for tau in (lifetimes[-1], lifetimes[0]))),
            natural_n2o_ppm=pair(*(x * y / 1e3 for x in nitrous['rise_ppb_per_year'])),
            natural_soil_mm=pair(*(s * y for s in SOIL_MM_YR)))
    growth = {f'{r:g}': {f'{n / 1e9:g}': sig(math.log(n / 1e8) / r) for n in POPULATIONS if n > 1e8}
              for r in GROWTH_PER_YEAR}
    return dict(by_years=out, years_from_0_1_billion=growth)


def main(argv=None) -> int:
    data = {}
    for key, (path, schema) in INPUTS.items():
        data[key] = json.loads(path.read_text())
        if schema is not None and data[key].get('schema') != schema:
            raise ValueError(f'{path}: schema {data[key].get("schema")!r}, expected {schema!r}')
    molar_mass = data['reference_air']['inventory']['mean_molar_mass']
    rings = dict(equator=data['ring'], lat_45n=data['ring_45n'], tilted_70=data['ring_70'], lat_80n=data['ring_80n'])
    land_table = land(data['atlas'], data['drainage'], data['aerosol'], data['rain'], data['plant'], data['canopy'])
    land_table['potential_equivalent_share_km2'] = pair(*(s * moon_area_m2() / 1e6
                                                          for s in land_table['potential_equivalent_share']))
    food_table = food(data['fruit'])
    need = needs(food_table)
    product = dict(
        schema=SCHEMA,
        producer=dict(domain='biosphere', files={'biosphere/ecology/people.py': digest(__file__)},
                      inputs={str(p.relative_to(ROOT)): digest(p) for p, _ in INPUTS.values()},
                      constants=constants_used(['biosphere/ecology/people.py'])),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        sky=sky(data['air'], data['ring'], rings, data['muons_moon'], data['muons_earth'], molar_mass),
        lift=lift(data['air'], molar_mass),
        aloft=aloft(data['plant']),
        traffic=traffic(data['fleet']),
        land=land_table,
        food=food_table,
        needs=need,
        ethylene=ethylene(land_table, molar_mass, data['air'], data['screen']['trace_gases']['air_mol']),
        molybdenum=mo_dressing(land_table),
        capacity=capacity(land_table, need, data['screen']),
        diet_limits_people=diet_limits(land_table, food_table),
        horizons=horizons(data['screen']))
    twilight = data['plant']['twilight']
    product['sky']['ground_twilight_hours_after_sunset_above_umol_m2_s'] = dict(
        moon={k: sig(v) for k, v in twilight['moon']['hours_after_sunset_above'].items()},
        earth={k: sig(v) for k, v in twilight['earth']['hours_after_sunset_above'].items()})
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    print(json.dumps({k: product[k] for k in ('sky', 'lift', 'aloft', 'traffic', 'land')}, indent=1)[:8000])
    return 0


if __name__ == '__main__':
    sys.exit(main())
