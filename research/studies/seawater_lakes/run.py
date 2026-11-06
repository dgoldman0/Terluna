"""Seawater lakes on the Open Moon: where they are easiest to hold, how their salt is kept out of the fresh waters, and
the physics of a one-way channel with fish passage.

Closed-form screens on the project's products and published values; no climate, cloud or lake model runs. Seawater
here is Earth's reference composition (Millero et al. 2008, as the sea chemistry study carries it) at 35.165 g/kg.

- The fill: salt and chloride per km2 and per metre of depth, and the delivered water whose own chloride (resources'
  question 17) would salt it.
- The water: the design climate's open-water evaporation and rain (the design years and the turned run with the refit),
  and open-water evaporation in the drylands from the storm-water study's pools; the fresh inflow per km2 that holds a
  terminal lake's salinity in a dry belt, and the surplus a lake in a wet belt has to shed.
- The chemistry: what each feed water brings per km2, what calcifiers lay from it, the alkalinity they cannot lay
  without taking the lake's own calcium (the soda path's start), and the dosing, brine bleed and salt top-up that hold
  the lake within a tenth of seawater's make-up.
- Isolation: seepage through a sealed floor where the lake stands above the water table, the floor's convective
  stability, the fresh-salt interface beside a lake held below its neighbours, wind setup and spray, all with gravity
  explicit.
- The channel: the fresh discharge that arrests a salt wedge at a sill (densimetric Froude number 1), fishway jets for a
  given drop, and a fish lock's exchange and salt per cycle, on the Moon and on Earth.

    /home/kir/Documents/Projects/terluna/climate/gcm/.venv/bin/python -m research.studies.seawater_lakes.run   # < 1 s

writes results/seawater_lakes.json (schema terluna.research.seawater-lakes/1) with the sha256 of every input read.
"""
from __future__ import annotations

import hashlib
import json
import math
import re
from pathlib import Path

from research.studies.sea_chemistry.run import SEAWATER, SEAWATER_G_KG
from shared.constants import DATA, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'seawater_lakes.json'
SCHEMA = 'terluna.research.seawater-lakes/1'
SEA_CHEMISTRY = ROOT / 'research' / 'studies' / 'sea_chemistry' / 'results' / 'sea_chemistry.json'
# The turned design years at the target (A28R_301k_design_l15_snow_s924_turned_refit, years 10-19): the means of the run's
# yearly drift check, carried here so the screen regenerates from the repository; the file sits on the research drive,
# and the test checks these values against it wherever it is present.
DRIFT_CHECK = Path('/media/projectspace/terluna-research/session-tools/turned_design/drift_check_step2.out')
TURNED_REFIT_YEARS = dict(
    source='the drift check of A28R_301k_design_l15_snow_s924_turned_refit, years 10-19 '
           '(session-tools/turned_design/drift_check_step2.out on the research drive, sha256 407ea30b075c4e65)',
    years=list(range(10, 20)), land_evaporation_mm_day=1.503, land_rain_mm_day=2.1925,
    sea_evaporation_mm_day=2.7343, sea_rain_mm_day=1.6317)
TURNED_REFIT_RUN = 'A28R_301k_design_l15_snow_s924_turned_refit'

G_MOON = MOON_SURFACE_GRAVITY                  # GM / R^2 from shared/constants.json, 1.6242 m/s2
G_EARTH = STANDARD_GRAVITY
YEAR_S = DATA['time']['julian_year_days'] * 86400.0
LUNAR_DAY_D = DATA['moon']['synodic_month_days']

# Molar masses, g/mol
M_CL, M_CACO3, M_CACL2, M_NACL = 35.453, 100.087, 110.98, 58.443

# --- Transcribed product values (each with its source; the README sections were read in full) ---------------------
DESIGN_YEARS = dict(
    source='climate/gcm/README.md, "The design run (1 November)": A28R_301k_design_l15_snow_s934, years 10-19',
    sea_rain_mm_day=1.51, sea_evaporation_mm_day=2.68, land_rain_mm_day=2.30, land_rain_less_evaporation_mm_day=0.73,
    semi_arid_rain_mm_yr=243.0)
# A pool over the semi-arid land in the design years: mean depth (m), loss in the lunar day's daylight and through the
# night (mm per lunar day); its day at 0.88 of Priestley-Taylor's demand at water's albedo, 0.99 with Lake Mead's factor
POOLS = dict(source='research/studies/storm_water/README.md, phase 2b, "The night costs a pool a few millimetres"',
             pools=((0.3, 62.0, 5.6), (3.0, 47.0, 18.0)), day_factor_lake_mead=0.99 / 0.88)
DESERT_RAIN_MM_DAY = 0.17          # research/studies/coast_water/README.md: the deserts' rain at 301 K
WET_BELT_RAIN_MM_YR = (1554.0, 2548.0)   # climate/gcm README: Tranquillitatis's land (design years) to the land box at 0 E
SPRINGS = dict(source='research/studies/springs/README.md, mature Moon, central',
               band_15_30=dict(alkalinity=1.54, ca=0.33, mg=0.41), band_30_45=dict(alkalinity=1.34, ca=0.24, mg=0.41))
HELD_SEA_MG_KEPT = (0.05, 0.5)     # register, 19 October: the held seas' magnesium kept at 0.05-0.5 mmol/kg
CRUST_TOP_UP_KM3_YR = 1700.0       # geography/README.md: the top-up that holds the seas, 1,000-10,000 years, central
LAKE_GROUNDWATER_GAIN_M_YR = 0.05  # geography/README.md: lakes taking groundwater in, median 5 cm a year (mature Moon)
SHORE_TRANSMISSIVITY_M2_S = 5.2e-4  # research/studies/coast_water/README.md: the central crust at the ground by the shores
SEEPAGE_MARGINS_M = (1.0, 3.0)      # design guess: the lake's level below its fresh neighbour's
SEEPAGE_DISTANCE_M = 1000.0         # design guess: the strip of ground between the fresh water and the lake
DELIVERED_CL_KG_M3 = 0.28          # sea chemistry and resources' question 17: ice like 67P's, 0.28 g/kg of chloride

# --- Design choices and Earth calibrations ---------------------------------------------------------------------------
LAKE_T_C = 29.5            # design guess: the design years' sea surface, 302.7 K, for lake and inflow alike
EXAMPLE_LAKES = ((1.0, 10.0), (10.0, 20.0), (100.0, 30.0), (1000.0, 50.0))     # km2, mean depth m (design guesses)
REFERENCE_DEPTH_M = 20.0
TOLERANCE = 0.1            # design choice: each major ion held within a tenth of seawater's
BRINE_EVAPORATION_SHARE = 0.7   # design guess: a concentrating brine evaporates at 0.7 of fresh water's rate
HALITE_BRINE_G_KG = 300.0       # design guess: the pans take the bleed to near halite saturation
RO_KWH_M3 = (3.0, 4.0)          # recalled: seawater reverse osmosis on Earth; membrane pressure, independent of gravity
LINER_K_EARTH_M_S = (1e-9, 1e-10, 1e-11)  # design values: Earth's compacted-clay liner rule, ten and a hundred times tighter
LINER_THICKNESS_M = 1.0
FLOOR_POROSITY = 0.4            # design guess: fine lake-floor sediment
SALT_DIFFUSIVITY_M2_S = 1e-9    # recalled: NaCl in water at 25-30 C (1.5e-9) with a pore tortuosity
WATER_VISCOSITY_PA_S = 0.80e-3  # recalled: water near 29.5 C
AIR_DENSITY = 1.216e5 / (287.05 * 300.0)   # the design air, 1.216 bar at 300 K
DRAG_COEFFICIENT = 1.3e-3       # Earth calibration: neutral 10 m drag over water near 10 m/s
AIR_VISCOSITY_PA_S = 1.85e-5    # recalled: air near 300 K
LOCK = dict(length=20.0, width=4.0, depth=2.0)   # design guess: a fish lock chamber, m
LOCK_RESIDUAL = (0.01, 0.05)    # design guess: salt left after flushing through the floor, a share of the chamber's
LOCK_FLUSH_VOLUMES = 1.5        # design guess: chamber volumes of fresh water per flush
FISHWAY_SPEEDS_M_S = (1.0, 1.5) # design values for small, weak swimmers; Earth fishway practice is gravity-free in speed
LOCK_FRONT_FROUDE = 0.5         # recalled: Benjamin (1968), Shin et al. (2004): full-depth lock exchange
EEL_SPAWNING_DEPTH_M = (150.0, 200.0)   # published: Tsukamoto et al. (2011), Japanese eel


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


# --- Seawater ---------------------------------------------------------------------------------------------------------
def pure_water_density(t_c: float) -> float:
    """Standard mean ocean water's density (kg/m3) at one atmosphere (Bigg 1967, as UNESCO's EOS-80 takes it)."""
    t = t_c
    return (999.842594 + 6.793952e-2 * t - 9.095290e-3 * t ** 2 + 1.001685e-4 * t ** 3
            - 1.120083e-6 * t ** 4 + 6.536332e-9 * t ** 5)


def seawater_density(practical_salinity: float, t_c: float) -> float:
    """Seawater's density (kg/m3) at one atmosphere, UNESCO's EOS-80 (Millero and Poisson 1981)."""
    s, t = practical_salinity, t_c
    a = 8.24493e-1 - 4.0899e-3 * t + 7.6438e-5 * t ** 2 - 8.2467e-7 * t ** 3 + 5.3875e-9 * t ** 4
    b = -5.72466e-3 + 1.0227e-4 * t - 1.6546e-6 * t ** 2
    c = 4.8314e-4
    return pure_water_density(t) + a * s + b * s ** 1.5 + c * s * s


def practical_from_reference(s_r_g_kg: float) -> float:
    """Practical salinity from reference salinity, S_R = (35.16504 / 35) S_P (Millero et al. 2008)."""
    return s_r_g_kg * 35.0 / 35.16504


def reduced_gravity(g: float, rho_salt: float, rho_fresh: float) -> float:
    return g * (rho_salt - rho_fresh) / (0.5 * (rho_salt + rho_fresh))


# --- The water --------------------------------------------------------------------------------------------------------
def mm_day_to_m_yr(x: float) -> float:
    return x * 1e-3 * YEAR_S / 86400.0


def pool_evaporation_m_yr(day_mm: float, night_mm: float, day_factor: float = 1.0) -> float:
    """A pool's yearly evaporation from its loss per lunar day, its daylight loss scaled by day_factor."""
    per_day = (day_mm * day_factor + night_mm) / LUNAR_DAY_D
    return mm_day_to_m_yr(per_day)


def read_drift_check(path: Path) -> dict:
    return parse_drift_check(path.read_text())


def parse_drift_check(text: str) -> dict:
    """Means over the listed years of the land's and sea's rain and evaporation (mm/day) in a drift_check output."""
    pat = re.compile(r'year\s+(\d+): land: .*?evap (-?[\d.]+) pr ([\d.]+) .*?\| sea: .*?evap (-?[\d.]+) pr ([\d.]+)')
    rows = [tuple(float(x) for x in m.groups()) for m in pat.finditer(text)]
    if not rows:
        raise ValueError('no years read from the drift check')
    n = len(rows)
    return dict(years=[int(r[0]) for r in rows],
                land_evaporation_mm_day=sum(abs(r[1]) for r in rows) / n, land_rain_mm_day=sum(r[2] for r in rows) / n,
                sea_evaporation_mm_day=sum(abs(r[3]) for r in rows) / n, sea_rain_mm_day=sum(r[4] for r in rows) / n)


# --- The chemistry ----------------------------------------------------------------------------------------------------
def feed_budget(feed: dict, net_evaporation_m_yr: float, depth_m: float, salt_kg_m3: float) -> dict:
    """What one km2 of a terminal seawater lake receives from its feed water each year, and what holds it.

    feed: alkalinity (meq/kg), ca, mg, na_k (mmol/kg; na_k in meq/kg). The inflow equals the net evaporation (the bleed's
    small extra is left out). Calcifiers lay the feed's calcium with twice its alkalinity as calcite; alkalinity beyond
    that, the excess, they can lay only with the lake's own calcium, which then falls. Dosing calcium chloride replaces
    that calcium; the bleed exports the feed's magnesium, sodium and potassium and the dose's chloride, each held within
    TOLERANCE of seawater's; the bleed's salt is replaced from outside as seawater's salts, or, with the halite the pans
    lay returned to the lake, only the bleed's other salts are."""
    q = net_evaporation_m_yr * 1e6                     # m3 a year per km2
    ca, mg = feed['ca'], feed['mg']
    na_k = feed.get('na_k', max(feed['alkalinity'] - 2.0 * (ca + mg), 0.0))
    excess = max(feed['alkalinity'] - 2.0 * ca, 0.0)    # meq/kg the feed's calcium cannot lay
    ca_sea = SEAWATER['ca'] * 1e3                       # mmol/kg
    draw = 0.5 * excess * net_evaporation_m_yr          # mol of the lake's calcium a m2 a year (mmol/kg ~ mol/m3)
    bleed_terms = dict(mg=mg / (TOLERANCE * SEAWATER['mg'] * 1e3),
                       na_k=na_k / (TOLERANCE * SEAWATER['na'] * 1e3),
                       dose_chloride=excess / (TOLERANCE * SEAWATER['cl'] * 1e3))
    bleed = max(bleed_terms.values())
    pan_share = bleed * (1.0 - SEAWATER_G_KG / HALITE_BRINE_G_KG) / BRINE_EVAPORATION_SHARE
    nacl_kg_m3 = salt_kg_m3 * SEAWATER['na'] * M_NACL / SEAWATER_G_KG   # the halite share of the lake's salt
    return dict(
        inflow_m3_yr_km2=q, inflow_l_s_km2=q / YEAR_S * 1e3,
        calcite_from_feed_t_yr_km2=ca * q * M_CACO3 * 1e-6,
        excess_alkalinity_meq_kg=excess,
        lake_calcium_tenth_lost_yr=(TOLERANCE * ca_sea * depth_m / draw) if draw > 0 else None,
        lake_calcium_spent_yr=(ca_sea * depth_m / draw) if draw > 0 else None,
        calcium_chloride_dose_t_yr_km2=0.5 * excess * q * M_CACL2 * 1e-6,
        bleed_share_of_inflow=bleed, bleed_terms=bleed_terms, bleed_m3_yr_km2=bleed * q,
        salt_top_up_t_yr_km2=bleed * q * salt_kg_m3 * 1e-3,
        salt_top_up_with_halite_returned_t_yr_km2=bleed * q * (salt_kg_m3 - nacl_kg_m3) * 1e-3,
        dose_chloride_t_yr_km2=excess * q * M_CL * 1e-6,
        pan_area_share_of_lake=pan_share)


# --- Isolation --------------------------------------------------------------------------------------------------------
def liner_seepage_m_yr(k_earth_m_s: float, g: float, depth_m: float, thickness_m: float) -> float:
    """Seepage through a liner with free drainage below: K (h + L) / L, the conductivity scaled from Earth's by g."""
    return k_earth_m_s * (g / G_EARTH) * (depth_m + thickness_m) / thickness_m * YEAR_S


def seepage_supply_m3_yr(area_km2: float, margin_m: float) -> float:
    """Groundwater a round lake held margin_m below the fresh water around it takes in through its whole shore, by
    Darcy's law across SEEPAGE_DISTANCE_M of ground at the shores' transmissivity, T dh / L per metre of shore. A lake
    beside one fresh water on one side takes a share of this."""
    perimeter = 2.0 * math.pi * math.sqrt(area_km2 * 1e6 / math.pi)
    return SHORE_TRANSMISSIVITY_M2_S * margin_m / SEEPAGE_DISTANCE_M * perimeter * YEAR_S


def critical_permeability(g: float, drho: float, layer_m: float) -> float:
    """Permeability (m2) at which dense water over a porous floor of thickness layer_m turns convective: the porous
    Rayleigh number k drho g H / (mu phi D) at 4 pi^2 (Horton, Rogers and Lapwood)."""
    return 4.0 * math.pi ** 2 * WATER_VISCOSITY_PA_S * FLOOR_POROSITY * SALT_DIFFUSIVITY_M2_S / (drho * g * layer_m)


def wind_setup_m(u10: float, fetch_m: float, depth_m: float, g: float, rho_water: float) -> float:
    """Steady wind setup across a lake, tau L / (rho g h)."""
    return DRAG_COEFFICIENT * AIR_DENSITY * u10 ** 2 * fetch_m / (rho_water * g * depth_m)


def stokes_settling_m_s(diameter_m: float, g: float, rho_drop: float) -> float:
    return (rho_drop - AIR_DENSITY) * g * diameter_m ** 2 / (18.0 * AIR_VISCOSITY_PA_S)


# --- The channel ------------------------------------------------------------------------------------------------------
def arrest_discharge_per_width(g_reduced: float, sill_depth_m: float) -> float:
    """Fresh discharge per metre of sill width that holds the salt wedge below the sill: densimetric Froude number 1,
    q = sqrt(g' h^3)."""
    return math.sqrt(g_reduced * sill_depth_m ** 3)


def pool_drop_for_speed(speed_m_s: float, g: float) -> float:
    """Drop between fishway pools whose jet reaches speed_m_s, v = sqrt(2 g dz)."""
    return speed_m_s ** 2 / (2.0 * g)


def lock_front_speed(g_reduced: float, depth_m: float) -> float:
    return LOCK_FRONT_FROUDE * math.sqrt(g_reduced * depth_m)


def build() -> dict:
    inputs = {}
    sea_chem = json.loads(SEA_CHEMISTRY.read_text())
    inputs[str(SEA_CHEMISTRY.relative_to(ROOT))] = sha256(SEA_CHEMISTRY)
    held = [r for r in sea_chem['runs']
            if (r['corner'], r['delivery'], r['sinks'], r['sediment']) == ('central', 'fresh', 'managed', 'kept_open')]
    held_row = held[0]['rows'][-1]
    drift = {k: v for k, v in TURNED_REFIT_YEARS.items() if k != 'source'}

    # Seawater and the density contrast
    s_p = practical_from_reference(SEAWATER_G_KG)
    rho_s = seawater_density(s_p, LAKE_T_C)
    rho_f = seawater_density(held_row['salinity_g_kg'] * 35.0 / 35.16504, LAKE_T_C)
    gp = dict(moon=reduced_gravity(G_MOON, rho_s, rho_f), earth=reduced_gravity(G_EARTH, rho_s, rho_f))
    salt_kg_m3 = rho_s * SEAWATER_G_KG * 1e-3
    cl_kg_m3 = rho_s * SEAWATER['cl'] * M_CL * 1e-3
    seawater = dict(reference_salinity_g_kg=SEAWATER_G_KG, practical_salinity=s_p, temperature_c=LAKE_T_C,
                    density_kg_m3=rho_s, fresh_density_kg_m3=rho_f, density_difference_kg_m3=rho_s - rho_f,
                    reduced_gravity_m_s2=gp, salt_kg_m3=salt_kg_m3, chloride_kg_m3=cl_kg_m3,
                    g_moon_m_s2=G_MOON, g_earth_m_s2=G_EARTH)

    # The fill
    lakes = []
    for area, depth in EXAMPLE_LAKES:
        vol = area * 1e6 * depth
        lakes.append(dict(area_km2=area, mean_depth_m=depth, volume_km3=vol * 1e-9, salt_mt=vol * salt_kg_m3 * 1e-9,
                          chloride_mt=vol * cl_kg_m3 * 1e-9,
                          delivered_water_km3_for_its_chloride=vol * cl_kg_m3 / DELIVERED_CL_KG_M3 * 1e-9))
    fill = dict(salt_kt_per_km2_per_m=salt_kg_m3, chloride_kt_per_km2_per_m=cl_kg_m3,
                salt_mt_per_km3=salt_kg_m3, lakes=lakes,
                delivered_water_km3_per_km3_of_seawater=cl_kg_m3 / DELIVERED_CL_KG_M3,
                seawater_km3_yr_from_the_crust_top_up_chloride=CRUST_TOP_UP_KM3_YR * DELIVERED_CL_KG_M3 / cl_kg_m3)

    # The water
    e_pools = [pool_evaporation_m_yr(d, n) for _, d, n in POOLS['pools']]
    e_pools_mead = [pool_evaporation_m_yr(d, n, POOLS['day_factor_lake_mead']) for _, d, n in POOLS['pools']]
    e_dry = (min(e_pools), max(e_pools_mead))
    p_semi = DESIGN_YEARS['semi_arid_rain_mm_yr'] * 1e-3
    p_desert = mm_day_to_m_yr(DESERT_RAIN_MM_DAY)
    net_semi = (e_dry[0] - p_semi, e_dry[1] - p_semi)
    net_desert = (e_dry[0] - p_desert, e_dry[1] - p_desert)
    net_dry = (net_semi[0], net_desert[1])
    e_open = (mm_day_to_m_yr(DESIGN_YEARS['sea_evaporation_mm_day']), mm_day_to_m_yr(drift['sea_evaporation_mm_day']))
    e_wet = (min(e_open), max(e_open))
    p_wet = tuple(x * 1e-3 for x in WET_BELT_RAIN_MM_YR)
    surplus = (p_wet[0] - e_wet[1], p_wet[1] - e_wet[0])
    water = dict(
        design_years=DESIGN_YEARS,
        turned_refit_years=dict(run=TURNED_REFIT_RUN, source=TURNED_REFIT_YEARS['source'], **drift),
        open_water_net_evaporation_m_yr=dict(
            design_years=mm_day_to_m_yr(DESIGN_YEARS['sea_evaporation_mm_day'] - DESIGN_YEARS['sea_rain_mm_day']),
            turned_refit=mm_day_to_m_yr(drift['sea_evaporation_mm_day'] - drift['sea_rain_mm_day'])),
        dryland_lake_evaporation_m_yr=dict(pools=e_pools, pools_with_lake_mead_factor=e_pools_mead, range=e_dry,
                                           source=POOLS['source']),
        dry_belt=dict(semi_arid_rain_m_yr=p_semi, desert_rain_m_yr=p_desert, net_evaporation_semi_arid_m_yr=net_semi,
                      net_evaporation_desert_m_yr=net_desert, net_evaporation_m_yr=net_dry,
                      inflow_m3_yr_per_km2=tuple(x * 1e6 for x in net_dry),
                      inflow_l_s_per_km2=tuple(x * 1e6 / YEAR_S * 1e3 for x in net_dry)),
        wet_belt=dict(rain_m_yr=p_wet, open_water_evaporation_m_yr=e_wet, surplus_m_yr=surplus,
                      surplus_m3_yr_per_km2=tuple(x * 1e6 for x in surplus),
                      reverse_osmosis_mw_per_km2=(surplus[0] * 1e6 * RO_KWH_M3[0] / (YEAR_S / 3600.0) * 1e-3,
                                                  surplus[1] * 1e6 * RO_KWH_M3[1] / (YEAR_S / 3600.0) * 1e-3),
                      salt_lost_if_spilled_mixed_kt_yr_per_km2=tuple(x * salt_kg_m3 for x in surplus),
                      years_to_spill_a_reference_lakes_salt=tuple(REFERENCE_DEPTH_M / x for x in surplus)))

    # The chemistry, per feed, at the dry belt's two ends
    feeds = dict(
        held_fresh_sea=dict(alkalinity=held_row['alkalinity_meq_kg'], ca=held_row['ca'], mg=held_row['mg'],
                            na_k=held_row['na_k']),
        **{f'held_fresh_sea_mg_{mg}': dict(alkalinity=held_row['alkalinity_meq_kg'] + 2.0 * mg, ca=held_row['ca'],
                                          mg=mg, na_k=held_row['na_k']) for mg in HELD_SEA_MG_KEPT},
        springs_15_30=dict(SPRINGS['band_15_30']), springs_30_45=dict(SPRINGS['band_30_45']))
    chemistry = {name: {f'{q:.3f}_m_yr': feed_budget(f, q, REFERENCE_DEPTH_M, salt_kg_m3) for q in net_dry}
                 for name, f in feeds.items()}
    chemistry_feeds = feeds

    # Isolation
    perched = [dict(k_earth_m_s=k, seepage_m_yr=liner_seepage_m_yr(k, G_MOON, REFERENCE_DEPTH_M, LINER_THICKNESS_M),
                    seepage_earth_m_yr=liner_seepage_m_yr(k, G_EARTH, REFERENCE_DEPTH_M, LINER_THICKNESS_M))
               for k in LINER_K_EARTH_M_S]
    for p in perched:
        p['salt_t_yr_per_km2'] = p['seepage_m_yr'] * 1e6 * salt_kg_m3 * 1e-3
        p['share_of_reference_salt_per_yr'] = p['seepage_m_yr'] / REFERENCE_DEPTH_M
    drho = rho_s - rho_f
    isolation = dict(
        perched_lake_through_a_liner=perched,
        gaining_lake=dict(inflow_m_yr=LAKE_GROUNDWATER_GAIN_M_YR,
                          salt_reach_against_inflow_m=SALT_DIFFUSIVITY_M2_S * YEAR_S / LAKE_GROUNDWATER_GAIN_M_YR),
        seepage_from_the_fresh_neighbour=[
            dict(area_km2=area, margin_m=m, distance_m=SEEPAGE_DISTANCE_M, transmissivity_m2_s=SHORE_TRANSMISSIVITY_M2_S,
                 m3_yr=seepage_supply_m3_yr(area, m), share_of_inflow=tuple(seepage_supply_m3_yr(area, m) / (area * 1e6 * q)
                                                                         for q in net_dry))
            for area, _ in EXAMPLE_LAKES for m in SEEPAGE_MARGINS_M],
        convective_floor=dict(layer_m=10.0, critical_permeability_moon_m2=critical_permeability(G_MOON, drho, 10.0),
                              critical_permeability_earth_m2=critical_permeability(G_EARTH, drho, 10.0)),
        interface_depth_per_m_of_fresh_head=rho_f / drho,
        wind_setup_m=[dict(u10_m_s=u, fetch_km=L * 1e-3, depth_m=h,
                           moon=wind_setup_m(u, L, h, G_MOON, rho_s), earth=wind_setup_m(u, L, h, G_EARTH, rho_s))
                      for u in (10.0, 20.0) for L in (3e3, 10e3, 30e3) for h in (5.0, 20.0)],
        spray=dict(diameter_um=20.0, settling_moon_m_s=stokes_settling_m_s(20e-6, G_MOON, rho_s),
                   settling_earth_m_s=stokes_settling_m_s(20e-6, G_EARTH, rho_s),
                   drift_km_from_10_m_at_5_m_s=dict(
                       moon=10.0 / stokes_settling_m_s(20e-6, G_MOON, rho_s) * 5.0 * 1e-3,
                       earth=10.0 / stokes_settling_m_s(20e-6, G_EARTH, rho_s) * 5.0 * 1e-3)))

    # The channel
    arrest = [dict(sill_depth_m=h, moon_m2_s=arrest_discharge_per_width(gp['moon'], h),
                   earth_m2_s=arrest_discharge_per_width(gp['earth'], h),
                   speed_moon_m_s=math.sqrt(gp['moon'] * h), speed_earth_m_s=math.sqrt(gp['earth'] * h))
              for h in (0.3, 0.5, 1.0, 2.0)]
    sills = []
    for area, _ in EXAMPLE_LAKES:
        for q in net_dry:
            flow = area * 1e6 * q / YEAR_S
            sills.append(dict(area_km2=area, net_evaporation_m_yr=q, inflow_m3_s=flow,
                              width_at_0_5_m_moon=flow / arrest_discharge_per_width(gp['moon'], 0.5),
                              width_at_0_5_m_earth=flow / arrest_discharge_per_width(gp['earth'], 0.5)))
    fishway = [dict(jet_m_s=v, drop_moon_m=pool_drop_for_speed(v, G_MOON), drop_earth_m=pool_drop_for_speed(v, G_EARTH),
                    pools_for_2_m_moon=math.ceil(2.0 / pool_drop_for_speed(v, G_MOON)),
                    pools_for_2_m_earth=math.ceil(2.0 / pool_drop_for_speed(v, G_EARTH)),
                    jet_densimetric_froude_moon_1_m_pool=v / math.sqrt(gp['moon'] * 1.0))
               for v in FISHWAY_SPEEDS_M_S]
    vol = LOCK['length'] * LOCK['width'] * LOCK['depth']
    lock = dict(chamber=LOCK, volume_m3=vol,
                front_moon_m_s=lock_front_speed(gp['moon'], LOCK['depth']),
                front_earth_m_s=lock_front_speed(gp['earth'], LOCK['depth']),
                exchange_s_moon=2.0 * LOCK['length'] / lock_front_speed(gp['moon'], LOCK['depth']),
                exchange_s_earth=2.0 * LOCK['length'] / lock_front_speed(gp['earth'], LOCK['depth']),
                salt_per_cycle_unflushed_t=vol * salt_kg_m3 * 1e-3,
                salt_per_cycle_flushed_t=tuple(r * vol * salt_kg_m3 * 1e-3 for r in LOCK_RESIDUAL),
                flush_m3=LOCK_FLUSH_VOLUMES * vol,
                cycles_a_day_the_inflow_flushes_per_10_km2=tuple(10.0 * q * 1e6 / 365.25 / (LOCK_FLUSH_VOLUMES * vol)
                                                                 for q in net_dry))
    channel = dict(arrest=arrest, sills=sills, fishway=fishway, lock=lock)
    fish = dict(eel_spawning_depth_earth_m=EEL_SPAWNING_DEPTH_M,
                lunar_depth_for_the_same_pressure_m=tuple(z * G_EARTH / G_MOON for z in EEL_SPAWNING_DEPTH_M))

    return dict(
        schema=SCHEMA,
        producer=dict(file='research/studies/seawater_lakes/run.py', sha256=sha256(Path(__file__)),
                      imports={'research/studies/sea_chemistry/run.py': sha256(ROOT / 'research/studies/sea_chemistry/run.py'),
                               'shared/constants.json': sha256(ROOT / 'shared/constants.json')}),
        inputs=inputs,
        evidence=('Closed-form screens (calculation) on project products and published or recalled values: the design '
                  'years and the turned run with the refit for open water, the storm-water pools for dryland lakes, '
                  'the sea chemistry\'s held fresh sea and seawater, the springs\' water, geography\'s groundwater. '
                  'No model was run. Densities by UNESCO EOS-80 at one atmosphere. Every design guess and Earth '
                  'calibration is named in run.py; the arrest, lock-exchange and Rayleigh criteria are recalled '
                  'textbook results. Nothing here is validated against a lunar observation.'),
        reading_rule=('Per km2 values are per km2 of lake surface; rates per year. Ranges are (low, high) across the '
                      'dry belt from semi-arid land with the pools\' evaporation to desert with Lake Mead\'s factor, or '
                      'across the wet belt\'s rain. Chemistry budgets are for a lake of REFERENCE_DEPTH_M (20 m), whose '
                      'inflow equals its net evaporation; times scale with depth. Moon and Earth entries differ only in '
                      'gravity. Concentrations in mmol/kg and alkalinity in meq/kg.'),
        seawater=seawater, fill=fill, water=water, feeds=chemistry_feeds, chemistry=chemistry,
        isolation=isolation, channel=channel, fish=fish)


def main(argv=None) -> int:
    record = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(record, indent=1) + '\n')
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
