"""Population and living across the horizons: how many people, where, using what, and what binds first.

    python -m provisioning.population
    # -> provisioning/results/population.json

Scenario arithmetic on committed products and the literature values listed in population_sources.json:
- **Horizons.** Today's 8.2 billion, the UN's peak and 2100 figure, and the growth that reaches the author's 20-30
  billion across the Solar System by the year 3000, 5000 or 12,026, against the growth rates of history.
- **Affluence.** Energy per person in today's affluent nations and the decent-living floor; the band the scenarios
  use; living space, mobility and water at the decent-living level.
- **Heat.** What each energy per person adds to the Moon under the planned climate and to Earth; the Moon's
  population at a warming tolerance, the added dimming that answers it, and the heat over a dense city.
- **Water, food and stocks.** Runoff, cropland and material stocks per person on Earth and on the Moon.
- **The array.** Habitat shielding, spin, the whole population's energy against the fleet's light.
- **Commuting.** Transfers from the Moon's surface to the ring radius: times, energy, and what rotations cost.
- **Earth service.** Power beamed from the array to Earth, computing over the light time, and an Earth sunshade
  against the film plant's output.
- **Scenarios.** Three placements of 20, 25 and 30 billion people, each at the band's energies.
"""
from __future__ import annotations
import hashlib
import json
import math
import sys
from pathlib import Path

from shared.constants import (EARTH_GM, EARTH_MOON_DISTANCE, EARTH_RADIUS, EARTH_SIDEREAL_DAY_S, JULIAN_YEAR_DAYS,
                              MOON_GM, MOON_RADIUS, MOON_SURFACE_GRAVITY, SPEED_OF_LIGHT, STANDARD_GRAVITY,
                              SYNODIC_MONTH_DAYS)
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SCHEMA = 'terluna.provisioning.population/1'
SCREEN = HERE / 'results' / 'first_screen.json'
HEAT = ROOT / 'research' / 'studies' / 'solar_shield_array' / 'results' / 'array_heat.json'
LINKS = ROOT / 'research' / 'studies' / 'solar_shield_array' / 'results' / 'electromagnetic.json'
LEDGERS = ROOT / 'research' / 'studies' / 'joint_synthesis' / 'results' / 'ledgers.json'
FLEET = ROOT / 'research' / 'studies' / 'sky_fleet' / 'results' / 'sky_fleet.json'
ATLAS = ROOT / 'geography' / 'results' / 'atlas.json'
DRAINAGE = ROOT / 'geography' / 'results' / 'drainage.json'
SOURCES = HERE / 'population_sources.json'
OUT = HERE / 'results' / 'population.json'
INPUTS = (SCREEN, HEAT, LINKS, LEDGERS, FLEET, ATLAS, DRAINAGE, SOURCES)

HOURS_PER_YEAR = JULIAN_YEAR_DAYS * 24.0

# Population.
POP_2024 = 8.2e9                          # UN DESA 2024, WPP 2024 Summary of Results
WPP_PEAK = (10.3e9, 2085)                 # "around 10.3 billion people in the mid-2080s" (WPP 2024)
POP_2100 = 10.2e9                         # WPP 2024, the end of the century
IHME_PEAK = (9.73e9, 2064)                # Vollset et al. 2020, reference scenario
IHME_2100 = 8.79e9                        # Vollset et al. 2020
HISTORY = {0: 0.30e9, 1000: 0.31e9, 1500: 0.50e9, 1750: 0.79e9, 1800: 0.98e9}   # UN 1999, Table 1
PEAK_GROWTH_PERCENT = 2.04                # late 1960s (UN 1999)
REPLACEMENT_TFR = 2.1                     # WPP 2024
GENERATION_YEARS = 30.0                   # mean length of a generation, for the net reproduction rate (assumption)
TARGETS_BILLION = (20.0, 25.0, 30.0)      # the author's 20-30 billion across the Solar System
TARGET_YEARS = (3000, 5000, 12026)        # a thousand, three thousand and ten thousand years on
BUILD_END = 2526                          # construction within 500 years of 2026 (decisions register)
MOON_START = 1e7                          # people on the Moon when the build ends (assumption)
SETTLEMENT_GROWTH = (0.01, 0.02)          # sustained growth by births and arrivals (assumption)

# Energy per person. OWID per-capita energy, Energy Institute Statistical Review 2026, substitution method, 2025.
OWID_2025_KWH = dict(world=20257.705, eu27=32289.914, high_income=52511.64, united_states=75051.33,
                     germany=33183.234, japan=37192.17, france=37612.65, united_kingdom=25128.771,
                     switzerland=27196.69)
EU_FINAL_PJ = 36755.0                     # EU final energy consumption, 2024 (Eurostat)
EU_POPULATION = 449.2e6                   # 1 January 2024 (Eurostat)
EU_GROSS_GJ = 124.0                       # EU gross available energy per person, 2024 (Eurostat)
EU_INDUSTRY_SHARE = 0.245                 # industry's share of EU final energy, 2024 (Eurostat)
DECENT_LIVING_GJ = (13.0, 15.3, 18.4)     # final energy, decent living (Millward-Hopkins et al. 2020)
NEEDS_PLATEAU_GJ = 60.0                   # final energy at which all assessed needs are met (Vogel et al. 2021)
TWO_THOUSAND_WATT_KW = 2.0                # primary energy target of the 2000-watt society (Schulz et al. 2008)
COMPUTING_KW = (0.1, 1.0)                 # a future person's computing, run in orbit (assumption)
DATA_CENTRES_TW = 0.047                   # the world's data centres in 2024, 415 TWh (IEA 2025)
ENERGY_GROWTH = 0.02                      # a steady 2% a year, to show what growth without end reaches
EARTH_PUMPED_STORAGE_TWH = 9.0            # world pumped-hydro energy capacity, of order (as provisioning/screen.py)
THERMAL_EFFICIENCY = 0.4                  # a thermal plant's electricity per unit heat, for comparison (assumption)
DECENT_FLOOR_M2 = 15.0 + 2.48 + 1.6 + 4.7 # shelter, education, health, collective space (Velez-Henao & Pauliuk 2023)
DECENT_PKM = 8300.0                       # passenger-km a year (Velez-Henao & Pauliuk 2023)
DECENT_WATER_L_DAY = 50.0                 # hygiene water (Velez-Henao & Pauliuk 2023; Gleick 1996)
METROPOLIS_FLOOR_M2 = 70.0                # floor per person for home, work and services (habitation/summit_metropolis)
TRAVEL_HOURS_DAY = 1.1                    # travel time budget (Schafer & Victor 2000)

# Heat.
GCM_K_PER_PERCENT_SUNLIGHT = 1.9          # climate/gcm/README.md, the dimmer branches (as provisioning/screen.py)
GCM_W_M2_PER_K = (1.0, 1.6)               # research/studies/atmospheric_co2/README.md (as provisioning/screen.py)
TOLERANCE_K = (0.1, 0.5)                  # warming of the planned climate taken before any answer (assumption)
EXTRA_DIMMING_PERCENT = 1.0               # an added dimming of 1% of sunlight, for scale
METROPOLIS_DENSITY_KM2 = 40000.0          # decisions register, 2026-09-27
METROPOLIS_PEOPLE = 100e6
GCM_CELL_KM = 170.0                       # the GCM's T21 cells (habitation/summit_metropolis)
FLANNER_W_M2 = {2005: 0.028, 2040: 0.059, 2100: 0.19}   # Earth's anthropogenic heat flux (Flanner 2009)
FLANNER_CELL_W_M2 = 3.0                   # cells above it warm 0.15-0.24 K in his runs (Flanner 2009)
ECS_K = 3.0                               # equilibrium climate sensitivity, best estimate (IPCC AR6 WG1)
F2X_W_M2 = 3.93                           # forcing of doubled CO2 (IPCC AR6 WG1)
AR6_TOTAL_ERF_W_M2 = 2.72                 # total anthropogenic forcing 1750-2019 (IPCC AR6 WG1)

# Food and land.
FARMLAND_HA = 4.8e9                       # agricultural land, 2021 (FAO 2023)
CROPLAND_HA = 1.6e9                       # cropland, 2021 (FAO 2023)
PLANT_DIET_LAND_CUT = 0.76                # farmland saved by diets without animal products (Poore & Nemecek 2018)
US_DIET_HA = (0.13, 1.08)                 # least-meat and baseline diets of ten US scenarios (Peters et al. 2016)
BOUNDARIES_BILLION = 10.2                 # fed within four planetary boundaries after transformation (Gerten et al. 2020)
ICE_FREE_KM2 = 130e6                      # Earth's ice-free land (IPCC SRCCL)

# Water.
EARTH_RUNOFF_KM3 = 45500.0                # annual discharge, the measure of availability (Oki & Kanae 2006)
FALKENMARK_M3 = dict(stress=1700.0, scarcity=1000.0, absolute=500.0)   # per person a year (Falkenmark et al. 1989)
WATER_FOOTPRINT_M3 = dict(world=1385.0, united_states=2842.0)          # per person a year (Hoekstra & Mekonnen 2012)

# Materials.
STOCKS_T = dict(decent=42.0, world=115.0, developed=335.0)  # in use per person (Velez-Henao & Pauliuk 2023; Krausmann et al. 2017)
STOCKS_2010_GT = 792.0                    # global in-use stocks in 2010 (Krausmann et al. 2017)
MATTER_STREAM_KG_S = 2.7e8                # representative bulk flow (engineering/reference/industrial_architecture, s. 6.2)

# The array's habitats. NASA SP-413 (Johnson & Holbrow 1977), the Stanford torus.
TORUS_PEOPLE = 10000.0
TORUS_SHIELD_T = 9.9e6                    # table 5-2
TORUS_SHIELD_T_M2 = 4.41                  # 441 g/cm2 of lunar soil, appendix 5E
TORUS_PROJECTED_M2 = 670000.0             # table 5-1, area required
TORUS_FARM_M2 = 20.0                      # per person, agriculture and food processing
TORUS_KW = 3.0                            # electrical power per person
DENSE_SHIELD_SHARE = 0.25                 # shield per person in a design with four times the floor per shield (assumption)
ARRAY_PEOPLE = (1e7, 1e8, 1e9)
SPIN_RPM = (2.0, 4.0)                     # residents tolerate at least 4 rpm (Globus & Hall 2017)

# Commuting between the Moon's surface and the ring radius.
TRAVELLER_HEAVY_KG = 1000.0               # a traveller with 1 t of craft (assumption); the light case is the sky fleet's
ASCENT_LOSSES_M_S = 500.0                 # drag and gravity losses through the air, pending the ascent study (assumption)
ISP_S = 465.5                             # LOX/LH2, RL10B-2 in vacuum (L3Harris)
MIXTURE_RATIO = 5.88                      # oxidiser to fuel by mass, RL10B-2
ELECTROLYSIS_KWH_KG = (50.0, 52.0)        # electricity per kg of hydrogen (IEA 2023)
PERILUNE_ALTITUDE_KM = 1300.0             # the return dips below the exobase, 1,300-2,600 km up (provisioning README)
PLATFORM_KM = 70.0                        # the high platforms (sky-fleet study)
ROTATIONS = dict(weekly=JULIAN_YEAR_DAYS / 7.0, each_lunar_cycle=JULIAN_YEAR_DAYS / SYNODIC_MONTH_DAYS,
                 every_other_cycle=JULIAN_YEAR_DAYS / SYNODIC_MONTH_DAYS / 2.0, quarterly=4.0)

# Earth service.
BEAM_GHZ = (2.45, 5.8)
TAU = 2.0                                 # aperture parameter for 90-95% collection (Rodenbeck et al. 2021)
SPS_1978 = dict(power_w=5e9, transmitter_km=1.0, rectenna_km=(10.0, 13.0), peak_w_m2=230.0)  # US DOE & NASA 1978
EARTH_BEAMED_TW = (20.0, 100.0)           # Criswell's 20 TW for 10 billion (Criswell 2002), and 100 TW
SHADE_T = 20e6                            # an L1 sunshade for 1.8% of sunlight (Angel 2006)
LLCD_BPS = 622e6                          # downlink from lunar orbit, 2013 (Boroson et al. 2014)

# Where the people are, in billions: scenarios for comparison at the author's 20-30 billion.
SCENARIOS = {
    'A': dict(name='Earth-centred, 20 billion', earth=16.0, moon=2.0, array=0.1, elsewhere=1.9),
    'B': dict(name='Earth and Moon, 25 billion', earth=18.0, moon=4.0, array=0.5, elsewhere=2.5),
    'C': dict(name='spread, 30 billion', earth=20.0, moon=5.0, array=1.0, elsewhere=4.0),
}

EVIDENCE = ('Scenario arithmetic on committed products (the provisioning first screen, the shield study\'s array '
            'heat and power links, the joint ledger, the sky fleet, the atlas and drainage) with the literature '
            'values named in population_sources.json. Closed form; no demographic, climate, food or orbit model is '
            'run. Population placements are scenarios for comparison at the author\'s 20-30 billion.')
READING_RULE = ('people in billions unless named; energy per person in kW of primary energy by the substitution '
                'method unless named final; heat in TW and W/m2 averaged over the named body; warming in K; dimming '
                'in percent of the design sunlight; land in km2 or ha as named; water in km3 a year and m3 per person '
                'a year; mass in Gt; times in hours or years; ranges pair the low and high ends of their inputs.')


def digest(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def kw(kwh_per_year):
    return kwh_per_year / HOURS_PER_YEAR


def gj_to_kw(gj):
    return gj * 1e9 / (HOURS_PER_YEAR * 3600.0) / 1e3


def growth_rate(start, end, years):
    return math.log(end / start) / years


def horizons():
    history = dict(year0_to_1750=growth_rate(HISTORY[0], HISTORY[1750], 1750),
                   year1000_to_1750=growth_rate(HISTORY[1000], HISTORY[1750], 750),
                   year1500_to_1750=growth_rate(HISTORY[1500], HISTORY[1750], 250))
    paths = {}
    for target in TARGETS_BILLION:
        for year in TARGET_YEARS:
            r = growth_rate(POP_2100, target * 1e9, year - 2100)
            paths[f'{target:g}_by_{year}'] = dict(percent_per_year=float(f'{100 * r:.4g}'), doubling_years=round(math.log(2) / r),
                                                  net_reproduction_rate=round(math.exp(r * GENERATION_YEARS), 4))
    settle = {f'{100 * g:g}%': dict(years_to_1_billion=round(math.log(1e9 / MOON_START) / g),
                                    year=round(BUILD_END + math.log(1e9 / MOON_START) / g))
              for g in SETTLEMENT_GROWTH}
    return dict(today_billion=POP_2024 / 1e9, wpp_peak=dict(billion=WPP_PEAK[0] / 1e9, year=WPP_PEAK[1]),
                wpp_2100_billion=POP_2100 / 1e9, ihme_peak=dict(billion=IHME_PEAK[0] / 1e9, year=IHME_PEAK[1]),
                ihme_2100_billion=IHME_2100 / 1e9,
                history_percent_per_year={k: round(100 * v, 4) for k, v in history.items()},
                peak_growth_percent=PEAK_GROWTH_PERCENT, replacement_fertility=REPLACEMENT_TFR,
                paths_from_2100=paths, moon_settlement=dict(start_people=MOON_START, from_year=BUILD_END, **settle))


def affluence(eu_final_kw, final_share):
    primary = {k: round(kw(v), 2) for k, v in OWID_2025_KWH.items()}
    band = dict(efficient=TWO_THOUSAND_WATT_KW, european=primary['eu27'], high_income=primary['high_income'],
                american=primary['united_states'])
    return dict(primary_kw_2025=primary, eu_final_kw_2024=round(eu_final_kw, 2), eu_final_share_of_gross=round(final_share, 3),
                eu_industry_share_of_final=EU_INDUSTRY_SHARE,
                decent_living_final_kw=[round(gj_to_kw(g), 2) for g in DECENT_LIVING_GJ],
                needs_plateau_final_kw=round(gj_to_kw(NEEDS_PLATEAU_GJ), 2),
                two_thousand_watt_primary_kw=TWO_THOUSAND_WATT_KW, band_primary_kw=band,
                computing_in_orbit_kw=list(COMPUTING_KW),
                living=dict(decent_floor_m2=round(DECENT_FLOOR_M2, 1), metropolis_floor_m2=METROPOLIS_FLOOR_M2,
                            torus_projected_m2=round(TORUS_PROJECTED_M2 / TORUS_PEOPLE), decent_pkm_per_year=DECENT_PKM,
                            decent_water_l_day=DECENT_WATER_L_DAY, travel_hours_per_day=TRAVEL_HOURS_DAY)), band


def growth(array_heat):
    """What a steady growth of energy use reaches from today's, in years."""
    today_tw = POP_2024 * OWID_2025_KWH['world'] / HOURS_PER_YEAR / 1e9
    earth_w_m2_per_tw = 1e12 / (4 * math.pi * EARTH_RADIUS ** 2)
    marks = {'flanner_2100_on_earth': FLANNER_W_M2[2100] / earth_w_m2_per_tw,
             'one_w_m2_on_earth': 1.0 / earth_w_m2_per_tw,
             'fleet_intercepted_low': array_heat['tiles']['intercepted_PW'][0] * 1e3,
             'fleet_intercepted_high': array_heat['tiles']['intercepted_PW'][1] * 1e3}
    return dict(today_tw=round(today_tw, 1), percent_per_year=100 * ENERGY_GROWTH,
                years_to={k: round(math.log(v / today_tw) / ENERGY_GROWTH) for k, v in marks.items()},
                tw={k: round(v) for k, v in marks.items()})


def placements(final_share, delivered_heat):
    """Heat released on the Moon per unit of primary energy a person uses, by where the energy is made and used."""
    return {'all_on_moon_thermal': 1.0,
            'industry_in_orbit_thermal': 1.0 - EU_INDUSTRY_SHARE,
            'industry_in_orbit_power_from_orbit': (1.0 - EU_INDUSTRY_SHARE) * final_share * delivered_heat}


def heat(screen, array_heat, band, place, moon_area):
    use = screen['heat']['human_use']
    w_m2_per_tw, tw_per_k = use['w_m2_per_tw'], use['tw_for_one_kelvin']
    sunlight = array_heat['placement']['design_sunlight_global_mean_W_m2']
    # One percent of sunlight warms by 1.9 K; a TW of heat by w_m2_per_tw / (1.0-1.6 W/m2 per K).
    dimming_per_tw = [w_m2_per_tw / g / GCM_K_PER_PERCENT_SUNLIGHT for g in GCM_W_M2_PER_K[::-1]]   # percent per TW
    tw_per_percent = [1.0 / d for d in dimming_per_tw[::-1]]
    earth_area = 4 * math.pi * EARTH_RADIUS ** 2
    earth_w_m2_per_tw = 1e12 / earth_area
    moon = {}
    for name, s in band.items():
        rows = {}
        for p, f in place.items():
            h = s * f                                                     # kW per person released on the Moon
            rows[p] = dict(kw_per_person=round(h, 3),
                           billion_at_tolerance={f'{t:g}K': [round(t * x / h, 2) for x in tw_per_k] for t in TOLERANCE_K},
                           billion_answered_by_1pct_dimming=[round(EXTRA_DIMMING_PERCENT * x / h, 1) for x in tw_per_percent],
                           dimming_percent_per_billion=[round(h * d, 4) for d in dimming_per_tw],
                           city_w_m2_at_40000_km2=round(METROPOLIS_DENSITY_KM2 * h * 1e3 / 1e6, 1),
                           summit_region_w_m2=round(METROPOLIS_PEOPLE * h * 1e3 / (GCM_CELL_KM * 1e3) ** 2, 2))
        moon[name] = rows
    lam = ECS_K / F2X_W_M2
    flanner_tw = FLANNER_W_M2[2100] / earth_w_m2_per_tw
    earth = {name: dict(billion_at_flanner_2100_thermal=round(flanner_tw / s, 1),
                        k_per_billion_thermal=round(s * earth_w_m2_per_tw * lam, 5)) for name, s in band.items()}
    glow = screen['heat']['fleet_infrared_w_m2']
    return dict(moon_w_m2_per_tw=w_m2_per_tw, moon_tw_per_k=tw_per_k, design_sunlight_w_m2=sunlight,
                dimming_percent_per_tw=[round(d, 5) for d in dimming_per_tw], tw_per_percent_dimming=[round(x, 1) for x in tw_per_percent],
                fleet_glow_w_m2=glow, fleet_glow_as_tw_used=[round(g / w_m2_per_tw) for g in glow],
                earth_w_m2_per_tw=round(earth_w_m2_per_tw, 6), moon_over_earth_per_tw=round(w_m2_per_tw / earth_w_m2_per_tw, 2),
                earth_k_per_w_m2=round(lam, 3), flanner_w_m2=FLANNER_W_M2, flanner_2100_tw=round(flanner_tw, 1),
                ar6_total_forcing_w_m2=AR6_TOTAL_ERF_W_M2,
                placement_factor={k: round(v, 3) for k, v in place.items()}, moon=moon, earth=earth,
                moon_area_km2=round(moon_area / 1e6))


def night(screen, band, final_share):
    """Energy a person on the Moon uses through one night, if it must be stored, against Earth's pumped storage."""
    hours = screen['night']['night_hours']
    per_person = {name: round(s * (1 - EU_INDUSTRY_SHARE) * final_share * hours) for name, s in band.items()}
    return dict(night_hours=hours, final_kwh_per_person=per_person,
                metropolis_brief_kwh_per_person=round(screen['night']['metropolis_night_twh'] * 1e9 / METROPOLIS_PEOPLE),
                earth_pumped_storage_kwh_per_person=round(EARTH_PUMPED_STORAGE_TWH * 1e9 / POP_2024, 2))


def water_food(screen, drainage, atlas, moon_area):
    moon_runoff_km3 = drainage['rivers']['discharge_to_sea_m3s'] * HOURS_PER_YEAR * 3600 / 1e9
    land = atlas['main_land_share'] + atlas['island_share']
    dry_km2 = moon_area / 1e6 * (land - drainage['lakes']['area_share'])
    lift_kwh_m3 = 1000.0 * MOON_SURFACE_GRAVITY * screen['waters']['rivers']['mean_land_height_above_sea_m'] / 3.6e6
    plan = screen['food']['planning_m2_per_person']
    potential = screen['food']['far_side_equator_m2_per_person']
    crop_share = CROPLAND_HA / 100.0 / ICE_FREE_KM2
    world_ha = FARMLAND_HA / POP_2024
    plant_ha = [US_DIET_HA[0], world_ha * (1 - PLANT_DIET_LAND_CUT)]
    return dict(
        water=dict(earth_runoff_km3=EARTH_RUNOFF_KM3, moon_runoff_km3=round(moon_runoff_km3),
                   moon_share_of_earth=round(moon_runoff_km3 / EARTH_RUNOFF_KM3, 3), falkenmark_m3=FALKENMARK_M3,
                   earth_billion_at={k: round(EARTH_RUNOFF_KM3 * 1e9 / v / 1e9, 1) for k, v in FALKENMARK_M3.items()},
                   moon_billion_at={k: round(moon_runoff_km3 * 1e9 / v / 1e9, 2) for k, v in FALKENMARK_M3.items()},
                   earth_today_m3_per_person=round(EARTH_RUNOFF_KM3 * 1e9 / POP_2024),
                   water_footprint_m3=WATER_FOOTPRINT_M3, moon_lift_to_mean_land_kwh_m3=round(lift_kwh_m3, 2)),
        food=dict(earth_farmland_ha_per_person_now=round(world_ha, 3), plant_rich_ha_per_person=[round(x, 3) for x in plant_ha],
                  american_diet_ha_per_person=US_DIET_HA[1],
                  earth_billion_on_todays_farmland=dict(plant_rich=[round(FARMLAND_HA / x / 1e9, 1) for x in plant_ha[::-1]],
                                                        american=round(FARMLAND_HA / US_DIET_HA[1] / 1e9, 1)),
                  within_boundaries_billion=BOUNDARIES_BILLION, earth_cropland_share_of_ice_free=round(crop_share, 3),
                  moon_dry_land_km2=round(dry_km2, -3), moon_m2_per_person_planning=plan, moon_m2_per_person_potential=potential,
                  moon_billion_at_earths_cropland_share=[round(crop_share * dry_km2 * 1e6 / m / 1e9, 1) for m in plan[::-1]],
                  moon_dry_land_share_per_billion=[round(m * 1e9 / (dry_km2 * 1e6), 4) for m in plan]),
        moon_runoff_km3=moon_runoff_km3, dry_km2=dry_km2, plant_ha=plant_ha, crop_share=crop_share)


def array(array_heat, ledgers, band):
    fleet = array_heat['fleet']
    shield = TORUS_SHIELD_T / TORUS_PEOPLE
    spin = {f'{rpm:g}rpm': dict(earth_g_radius_m=round(STANDARD_GRAVITY / (rpm * 2 * math.pi / 60) ** 2),
                                moon_g_radius_m=round(MOON_SURFACE_GRAVITY / (rpm * 2 * math.pi / 60) ** 2, 1))
            for rpm in SPIN_RPM}
    stream_gt_yr = MATTER_STREAM_KG_S * HOURS_PER_YEAR * 3600 / 1e12
    habitats = {(f'{p / 1e9:g}_billion' if p >= 1e9 else f'{p / 1e6:g}_million'): dict(shield_gt=shield_gt(p / 1e9),
                                             times_fleet_mass=round(shield * p / 1e9 / max(fleet['optical_mass_Gt']), 2),
                                             years_of_film_plant_output=[round(shield * p / 1e9 / u, 1)
                                                                         for u in fleet['film_upkeep_Gt_per_year'][::-1]],
                                             days_of_matter_stream=round(shield * p / 1e9 / stream_gt_yr * 365.25, 1),
                                             farm_km2=round(TORUS_FARM_M2 * p / 1e6), power_tw=round(TORUS_KW * p / 1e9, 3))
                for p in ARRAY_PEOPLE}
    intercepted = array_heat['tiles']['intercepted_PW']
    night = ledgers['energy']['fleet_night_half_intercepts_w']
    comp = array_heat['computing']
    return dict(shield_t_per_person=shield, shield_t_m2=TORUS_SHIELD_T_M2,
                shield_area_m2_per_person=round(shield / TORUS_SHIELD_T_M2), dense_share=DENSE_SHIELD_SHARE,
                habitats=habitats, spin=spin, matter_stream_gt_per_year=round(stream_gt_yr),
                fleet_mass_gt=[round(m, 1) for m in fleet['optical_mass_Gt']], fleet_area_km2=[round(a, -6) for a in fleet['area_km2']],
                film_upkeep_gt_per_year=[round(u, 2) for u in fleet['film_upkeep_Gt_per_year']],
                intercepted_pw=[round(i) for i in intercepted], night_half_pw=round(night / 1e15),
                collector_km2_per_tw=[round(c) for c in comp['collector_km2_per_TW']],
                radiator_km2_per_tw=round(comp['radiator_panel_km2_per_TW']), data_centres_tw=DATA_CENTRES_TW,
                moon_w_m2_per_tw_in_orbit=array_heat['placement']['moon_W_m2_per_TW_released_in_orbit'])


def shield_gt(billion):
    """Habitat shielding for a number of people in the array: a denser design and the Stanford torus's."""
    per_person = TORUS_SHIELD_T / TORUS_PEOPLE
    return [round(per_person * DENSE_SHIELD_SHARE * billion, 1), round(per_person * billion, 1)]


def conic_transfer(r0, r1, apo):
    """Time and delta-v from rest on a sphere of radius r0 to a circular orbit at r1, on a conic with periapsis r0.

    apo is the apoapsis radius (an ellipse) or None for the parabola. Departure takes the periapsis speed; arrival
    turns the arriving velocity into the circular one. Vacuum, two bodies.
    """
    mu = MOON_GM
    vc = math.sqrt(mu / r1)
    if apo is None:
        v0 = math.sqrt(2 * mu / r0)
        theta = math.acos(2 * r0 / r1 - 1)
        d = math.tan(theta / 2)
        t = math.sqrt(2 * r0 ** 3 / mu) * (d + d ** 3 / 3)
        v1, gamma = math.sqrt(2 * mu / r1), theta / 2
    else:
        a, e = (r0 + apo) / 2, (apo - r0) / (apo + r0)
        v0 = math.sqrt(mu * (2 / r0 - 1 / a))
        theta = math.acos(min(1.0, (a * (1 - e * e) / r1 - 1) / e))
        ecc = 2 * math.atan(math.sqrt((1 - e) / (1 + e)) * math.tan(theta / 2))
        t = (ecc - e * math.sin(ecc)) * math.sqrt(a ** 3 / mu)
        v1 = math.sqrt(mu * (2 / r1 - 1 / a))
        gamma = math.atan2(e * math.sin(theta), 1 + e * math.cos(theta))
    dv1 = math.sqrt(v1 ** 2 + vc ** 2 - 2 * v1 * vc * math.cos(gamma))
    return t / 3600, v0, dv1


def commute(array_heat, sky_fleet):
    r0, r1 = MOON_RADIUS, array_heat['design']['orbit_radius_km'] * 1e3
    traveller = (sky_fleet['long_haul']['passenger_kg'], TRAVELLER_HEAVY_KG)   # 400 kg a long-haul passenger
    ve = ISP_S * STANDARD_GRAVITY
    h2_share = 1 / (1 + MIXTURE_RATIO)
    minimum = MOON_GM / r0 - MOON_GM / (2 * r1)                            # J/kg from rest to the circular orbit
    rp = r0 + PERILUNE_ALTITUDE_KM * 1e3
    back = math.sqrt(MOON_GM / r1) - math.sqrt(MOON_GM * (2 / r1 - 2 / (r1 + rp)))
    transfers = {}
    for name, apo in (('hohmann', r1), ('apolune_twice_ring', 2 * r1), ('parabolic', None)):
        hours, dv0, dv1 = conic_transfer(r0, r1, apo)
        up = dv0 + dv1 + ASCENT_LOSSES_M_S
        prop_per_kg = (math.exp(up / ve) - 1) + (math.exp(back / ve) - 1)  # propellant per kg arriving, both legs
        mwh = [[round(prop_per_kg * m * h2_share * k / 1e3, 2) for k in ELECTROLYSIS_KWH_KG] for m in traveller]
        transfers[name] = dict(hours=round(hours, 1), delta_v_vacuum_m_s=round(dv0 + dv1), with_losses_m_s=round(up),
                               propellant_kg_per_kg_round_trip=round(prop_per_kg, 3), round_trip_mwh_per_traveller=mwh,
                               continuous_kw_per_commuter={k: [round(n * mwh[0][0] * 1e3 / HOURS_PER_YEAR, 1),
                                                               round(n * mwh[1][1] * 1e3 / HOURS_PER_YEAR, 1)]
                                                           for k, n in ROTATIONS.items()})
    return dict(ring_radius_km=r1 / 1e3, traveller_kg=list(traveller), minimum_kwh_per_kg=round(minimum / 3.6e6, 3),
                minimum_kwh_per_traveller=[round(minimum * m / 3.6e6) for m in traveller],
                lift_to_platform_kwh=[round(m * MOON_SURFACE_GRAVITY * PLATFORM_KM * 1e3 / 3.6e6, 1) for m in traveller],
                return_delta_v_m_s=round(back), exhaust_m_s=round(ve), h2_share_of_propellant=round(h2_share, 4),
                kwh_per_kg_propellant=[round(h2_share * k, 2) for k in ELECTROLYSIS_KWH_KG],
                rotations_per_year={k: round(v, 2) for k, v in ROTATIONS.items()}, transfers=transfers,
                light_round_trip_s=round(2 * (r1 - r0) / SPEED_OF_LIGHT, 3),
                sky_fleet=dict(long_haul_mean_trip_km=sky_fleet['long_haul']['mean_trip_km'],
                               winged_liner_mean_trip_hours=sky_fleet['long_haul']['fast_winged']['mean_trip_hours']))


def earth_service(links, array_heat, film_upkeep):
    mw = next(c for c in links['link_cases'] if c['kind'] == 'microwave')
    delivered_heat = 1.0 / (mw['capture'] * mw['rx_efficiency'] * mw['bus_efficiency'])   # heat on Earth per W delivered
    geo = (EARTH_GM * EARTH_SIDEREAL_DAY_S ** 2 / (4 * math.pi ** 2)) ** (1 / 3) - EARTH_RADIUS
    rectenna_m2 = math.pi / 4 * SPS_1978['rectenna_km'][0] * SPS_1978['rectenna_km'][1] * 1e6
    lam_1978 = SPEED_OF_LIGHT / (BEAM_GHZ[0] * 1e9)
    tau_1978 = math.sqrt(math.pi / 4 * (SPS_1978['transmitter_km'] * 1e3) ** 2 * rectenna_m2) / (lam_1978 * geo)
    transmitters = {}
    for ghz in BEAM_GHZ:
        lam = SPEED_OF_LIGHT / (ghz * 1e9)
        for name, d in (('geo', geo), ('moon', EARTH_MOON_DISTANCE)):
            area = (TAU * lam * d) ** 2 / rectenna_m2
            transmitters[f'{ghz:g}GHz_{name}_km'] = round(math.sqrt(4 * area / math.pi) / 1e3, 2)
    km2_per_tw = 1e12 / (SPS_1978['power_w'] / rectenna_m2) / 1e6
    collectors = array_heat['computing']['collector_km2_per_TW']
    return dict(link_dc_efficiency=mw['dc_efficiency'], heat_on_earth_per_w_delivered=round(delivered_heat, 3),
                thermal_plant_heat_per_w=round(1 / THERMAL_EFFICIENCY, 2),
                geo_altitude_km=round(geo / 1e3), tau_1978=round(tau_1978, 2), transmitter_diameter=transmitters,
                rectenna_mean_w_m2=round(SPS_1978['power_w'] / rectenna_m2, 1), rectenna_km2_per_tw=round(km2_per_tw, -1),
                rectenna_km2={f'{tw:g}TW': round(tw * km2_per_tw, -3) for tw in EARTH_BEAMED_TW},
                rectenna_share_of_cropland={f'{tw:g}TW': round(tw * km2_per_tw / (CROPLAND_HA / 100), 3) for tw in EARTH_BEAMED_TW},
                tw_on_one_percent_of_ice_free_land=round(0.01 * ICE_FREE_KM2 / km2_per_tw, 1),
                collector_km2_per_tw_delivered=[round(c / mw['dc_efficiency'], -1) for c in collectors],
                light_time_s=round(EARTH_MOON_DISTANCE / SPEED_OF_LIGHT, 3),
                light_round_trip_s=round(2 * EARTH_MOON_DISTANCE / SPEED_OF_LIGHT, 2), llcd_bps=LLCD_BPS,
                shade_days_of_film_plant=[round(SHADE_T / 1e9 / u * 365.25, 1) for u in film_upkeep[::-1]]), delivered_heat


def scenarios(band, place, wf, heat_out, arr, delivered_heat_earth, final_share):
    earth_area = 4 * math.pi * EARTH_RADIUS ** 2
    lam = ECS_K / F2X_W_M2
    out = {}
    for key, s in SCENARIOS.items():
        total = s['earth'] + s['moon'] + s['array'] + s['elsewhere']
        rows = {}
        for name, e in band.items():
            civ_tw = total * e
            moon_tw = {p: s['moon'] * e * f for p, f in place.items()}
            d = heat_out['dimming_percent_per_tw']
            earth_thermal = s['earth'] * e * 1e12 / earth_area
            earth_beamed = s['earth'] * e * final_share * delivered_heat_earth * 1e12 / earth_area
            rows[name] = dict(
                civilization_tw=round(civ_tw), share_of_night_half=round(civ_tw * 1e12 / (arr['night_half_pw'] * 1e15), 6),
                collectors_km2=[round(civ_tw * k, -3) for k in arr['collector_km2_per_tw']],
                collector_share_of_fleet_area=[round(civ_tw * k / arr['fleet_area_km2'][0], 6) for k in arr['collector_km2_per_tw']],
                moon_tw={p: round(v, 2) for p, v in moon_tw.items()},
                moon_w_m2={p: round(v * heat_out['moon_w_m2_per_tw'], 3) for p, v in moon_tw.items()},
                moon_dimming_percent={p: [round(v * x, 3) for x in d] for p, v in moon_tw.items()},
                earth_tw=round(s['earth'] * e), earth_w_m2=dict(thermal=round(earth_thermal, 3), power_from_orbit=round(earth_beamed, 3)),
                earth_k=dict(thermal=round(earth_thermal * lam, 3), power_from_orbit=round(earth_beamed * lam, 3)))
        out[key] = dict(name=s['name'], billion=dict(earth=s['earth'], moon=s['moon'], array=s['array'],
                                                    elsewhere=s['elsewhere'], total=round(total, 2)),
                        by_energy=rows,
                        moon_water_m3_per_person=round(wf['moon_runoff_km3'] * 1e9 / (s['moon'] * 1e9)),
                        earth_water_m3_per_person=round(EARTH_RUNOFF_KM3 * 1e9 / (s['earth'] * 1e9)),
                        moon_food_share_of_dry_land=[round(s['moon'] * 1e9 * m / (wf['dry_km2'] * 1e6), 3)
                                                     for m in wf['food']['moon_m2_per_person_planning']],
                        earth_food_ha_plant_rich=[round(s['earth'] * 1e9 * x / 1e9, 2) for x in sorted(wf['plant_ha'])],
                        earth_over_boundaries=round(s['earth'] / BOUNDARIES_BILLION, 2),
                        computing_tw=[round(total * k, 1) for k in COMPUTING_KW],
                        computing_times_data_centres=[round(total * k / DATA_CENTRES_TW) for k in COMPUTING_KW],
                        stocks_gt=dict(earth=round(s['earth'] * STOCKS_T['developed']), moon=round(s['moon'] * STOCKS_T['developed']),
                                       times_2010=round((s['earth'] + s['moon']) * STOCKS_T['developed'] / STOCKS_2010_GT, 1),
                                       array_shield=shield_gt(s['array'])))
    return out


def binds_first(band, place, heat_out, wf):
    """Populations at which each limit is reached, for the European energy, ordered from the first."""
    e = band['european']
    moon_kw = e * place['industry_in_orbit_thermal']
    tw_per_k, tw_per_pct = heat_out['moon_tw_per_k'], heat_out['tw_per_percent_dimming']
    moon = [('heat, unanswered, at 0.1 K', [round(0.1 * x / moon_kw, 2) for x in tw_per_k]),
            ('water stress, 1,700 m3 a person', [wf['water']['moon_billion_at']['stress']] * 2),
            ('water scarcity, 1,000 m3 a person', [wf['water']['moon_billion_at']['scarcity']] * 2),
            ("cropland at Earth's share of land", wf['food']['moon_billion_at_earths_cropland_share']),
            ('heat, answered by 1% more dimming', [round(x / moon_kw, 1) for x in tw_per_pct])]
    flanner_tw = FLANNER_W_M2[2100] / heat_out['earth_w_m2_per_tw']
    earth = [("American diets on today's farmland", [wf['food']['earth_billion_on_todays_farmland']['american']] * 2),
             ('food within four planetary boundaries', [BOUNDARIES_BILLION] * 2),
             ("heat from thermal plants at Flanner's 2100 level", [round(flanner_tw / e, 1)] * 2),
             ('water stress, 1,700 m3 a person', [wf['water']['earth_billion_at']['stress']] * 2),
             ("plant-rich diets on today's farmland", wf['food']['earth_billion_on_todays_farmland']['plant_rich'])]
    order = lambda rows: [dict(limit=k, billion=v) for k, v in sorted(rows, key=lambda r: r[1][0])]
    return dict(energy_kw=round(e, 2), moon_kw_released=round(moon_kw, 2), moon=order(moon), earth=order(earth))


def main(argv=None) -> int:
    screen, array_heat, links, ledgers, sky_fleet, atlas, drainage = (
        json.loads(p.read_text()) for p in (SCREEN, HEAT, LINKS, LEDGERS, FLEET, ATLAS, DRAINAGE))
    moon_area = 4 * math.pi * MOON_RADIUS ** 2
    eu_final_kw = EU_FINAL_PJ * 1e15 / EU_POPULATION / (HOURS_PER_YEAR * 3600) / 1e3
    final_share = eu_final_kw / gj_to_kw(EU_GROSS_GJ)
    aff, band = affluence(eu_final_kw, final_share)
    service, delivered_heat = earth_service(links, array_heat, array_heat['fleet']['film_upkeep_Gt_per_year'])
    place = placements(final_share, delivered_heat)
    heat_out = heat(screen, array_heat, band, place, moon_area)
    wf = water_food(screen, drainage, atlas, moon_area)
    arr = array(array_heat, ledgers, band)
    files = ['provisioning/population.py']
    product = dict(
        schema=SCHEMA,
        producer=dict(domain='provisioning', files={f: digest(ROOT / f) for f in files},
                      inputs={str(p.relative_to(ROOT)): digest(p) for p in INPUTS}, constants=constants_used(files)),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        horizons=horizons(), energy_growth=growth(array_heat), affluence=aff, heat=heat_out,
        night=night(screen, band, final_share),
        water=wf['water'], food=wf['food'], materials=dict(stocks_t_per_person=STOCKS_T, stocks_2010_gt=STOCKS_2010_GT),
        array=arr, commute=commute(array_heat, sky_fleet), earth_service=service,
        scenarios=scenarios(band, place, wf, heat_out, arr, delivered_heat, final_share),
        binds_first=binds_first(band, place, heat_out, wf))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    print(json.dumps({k: product[k] for k in ('binds_first', 'commute')}, indent=1)[:6000])
    return 0


if __name__ == '__main__':
    sys.exit(main())
