"""How much of the young Moon's terrain lasts without upkeep, and what keeping or renewing it costs.

    python -m research.studies.terrain_through_time.run      # results/terrain_through_time.json, seconds

On 13 October 2026 the author asked how much of the young terraformed Moon's small pools, lakes and landforms is kept
and how much is left to wear down, and what could reasonably rebuild the terrain: sediment-return blocks, buildings
and terraces for the built land and the gardens, and something for the wild. This runner reads the routing summary
(routing.py, which reruns geography's drainage) and works out:

- the lakes through time: how many, and how much of their area, stay open water after each year, by size and band;
- the small craters below the 1.9 km grid: their numbers at the saturation equilibrium of the Moon's small craters
  and how fast soil creep and accretion fill them, from a crater profile diffused in dimensionless time;
- the seas' fill from the rock the routing sends them and their own deposits (the sea chemistry);
- the relief: how fast each band lowers and how long its local relief takes to halve without uplift;
- upkeep: carrying the rock back uphill, in place or onto each catchment's highest third, with its lift, haulage and
  dredging, and the seas' deposits as stone for the built land;
- travertine: the calcite saturation the rivers and springs reach as they lose CO2 to the air, against the soils'
  ageing (the carbon cycle), and where it passes the saturation at which Earth's streams build tufa;
- dam builders: the channels whose slope and flows let dams stand, on the thresholds of Earth's beaver-dam capacity
  model translated to lunar gravity.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from research.studies.carbon_cycle import run as cc
from research.studies.terrain_through_time import routing as rt
from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'terrain_through_time.json'
SCHEMA = 'terluna.research.terrain-through-time/1'
FILES = ('research/studies/terrain_through_time/run.py', 'research/studies/terrain_through_time/routing.py',
         'research/studies/carbon_cycle/run.py')
ROUTING = 'research/studies/terrain_through_time/results/terrain_routing.json'
SEA_CHEM = 'research/studies/sea_chemistry/results/sea_chemistry.json'
SEA_P = 'resources/results/sea_phosphorus.json'
GATES = 'research/studies/conservation/results/gates.json'
INPUTS = {ROUTING: rt.SCHEMA, SEA_CHEM: 'terluna.research.sea-chemistry/1', SEA_P: 'terluna.resources.sea-phosphorus/1',
          GATES: None}
YEAR_S = JULIAN_DAY * JULIAN_YEAR_DAYS
G_RATIO = MOON_SURFACE_GRAVITY / STANDARD_GRAVITY
WATER_DENSITY = 1000.0

# Small craters. Pike (1977): fresh simple craters stand 0.196 D from rim to floor with a rim 0.036 D high; their
# ejecta thins as the cube of distance (McGetchin et al. 1973).
CRATER_DEPTH = 0.196
CRATER_RIM = 0.036
CRATER_DIAMETERS_M = (10.0, 30.0, 100.0, 300.0, 1000.0)
# Soil creep's diffusivity on Earth's soil-mantled hillslopes, m2 a year: humid forested Oregon Coast Range 0.003 and
# semi-arid Gabilan Mesa 0.038 (Roering et al. 2007), with the geometric middle; creep rises with moisture and with
# life (Richardson et al. 2019).
KAPPA_M2_YR = (0.003, 0.01, 0.04)
# The Moon's small craters now, by topographic diffusion under micrometeorites (Fassett & Thomson 2014), m2 a year.
KAPPA_AIRLESS_M2_YR = 5.5e-6
# Saturation equilibrium of small lunar craters: cumulative density as a share of geometric saturation, 1.54 D^-2
# per km2 with D in km (Gault 1970), 0.69-3.9% where equilibrium holds below about 1 km (Xiao & Werner 2015).
GEOMETRIC_SATURATION = 1.54
EQUILIBRIUM_SHARE = (0.0069, 0.039)

# Travertine. Earth's streams deposit calcite conspicuously above a mean annual saturation index of 0.75 (Merz-Preiss
# & Riding 1999); Plitvice's barriers grow where the water stands above three times saturation with pH above 8.
# Engineered builders, raising the saturation at their surfaces as charophytes and corals do, are taken to build from
# twice saturation (design guess).
SI_TUFA = 0.75
SI_PLITVICE = math.log10(3.0)
SI_DESIGNED = math.log10(2.0)
SPRING_FACTOR = (1.5, 2.0, 3.0)         # design guess: springs carrying 1.5-3 times the rivers' mean load
TRAVERTINE_YEARS = (0, 10, 100, 300, 1000, 3000, 1e4, 3e4, 1e5, 1e6)

# Dam builders: Macfarlane et al. (2017)'s beaver-dam capacity model, stream power in W per m of channel: dams are
# built below 150 at baseflow (Q80) and persist below 1,200 in the two-year flood; up to 190 and 2,400 they are built
# with breaches and blowouts; none on slopes over 23%. Dams a km: frequent 4-15, pervasive 15-40.
BRAT_LIMITS_W_M = dict(dams_persist=(150.0, 1200.0), with_blowouts=(190.0, 2400.0))
BRAT_SLOPE_MAX = 0.23
BRAT_DENSITY_PER_KM = dict(frequent=(4.0, 15.0), pervasive=(15.0, 40.0))
# Design guesses for the flows the 1.9 km channels carry: baseflow (Q80) and the two-year flood over the mean flow.
FLOW_REGIMES = dict(steady=(0.3, 4.0), flashy=(0.2, 10.0))
POND_HA_PER_DAM = 0.35                  # Earth's typical beaver impoundment, mostly 0.01-0.5 ha (Martin et al. 2015)

# Upkeep. Ideal lift at lunar gravity; haulage at the world's rail freight energy, 107 kJ a tonne-km (IEA 2019);
# dredging at the sea phosphorus study's hopper dredging; Earth's sand and gravel, 40-50 Gt a year (UNEP 2019); and
# loess falls of about 10 mm a century in the Holocene and 21 in the last glacial maximum on the Chinese Loess Plateau
# (Kohfeld & Harrison 2003 at 1.5 t/m3).
RAIL_KJ_PER_TONNE_KM = 107.0
EARTH_SAND_GRAVEL_GT_YR = (40.0, 50.0)
LOESS_MM_PER_CENTURY = (10.0, 21.0)

EVIDENCE = ('Calculations on committed products: the routing summary (geography\'s drainage rerun, erosion by relief '
            'and rain, units filled in log-spaced steps), the sea chemistry\'s works and the sea phosphorus study\'s '
            'dredging. Crater degradation is linear diffusion of Pike\'s fresh profile on a radial grid; travertine '
            'is the carbon cycle\'s carbonate solver at the land\'s air temperature; dam builders take Earth\'s '
            'capacity thresholds through equal Shields numbers. Earth\'s values are reference points, design guesses '
            'are named, and nothing here is validated against a lunar landscape.')
READING_RULE = ('Years from the start of the wet world. Survival shares are of the routing\'s rain-fed lakes by count '
                'and by area at the start. Rates of rock in Gt a year; power in GW; lowering in m a million years and '
                'mm a century; craters by rim-to-floor diameter in m; saturation as Omega and its log10, SI. Cases '
                'are the routing\'s: earliest (100 m/Myr, loose sediment, fast accretion), central and latest (10 '
                'm/Myr, dense, slow).')


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(path, schema=None):
    data = json.loads((ROOT / path).read_text())
    if schema and data.get('schema') != schema:
        raise ValueError(f'{path}: unexpected schema {data.get("schema")}')
    return data


rnd = rt.rnd


# --- Lakes and seas --------------------------------------------------------------------------------------------------

KEY_YEARS = ('1000', '10000', '100000', '1e+06', '1e+07')


def lakes(routing, sea_chem):
    out = {}
    for name, f in routing['fill'].items():
        s = f['survival']['all']
        out[name] = dict(
            open_share_by_count={y: s['open_count'][y] / s['lakes'] for y in KEY_YEARS},
            open_share_by_area={y: s['open_area_km2'][y] / s['area_km2'] for y in KEY_YEARS},
            median_life_yr=s['median_life_yr'], area_weighted_median_life_yr=s['area_weighted_median_life_yr'],
            median_life_by_size_yr={k: v['median_life_yr'] for k, v in f['survival'].items() if k != 'all'},
            all_but_one_percent_filled_yr=f['lakes_capacity_one_percent_left_yr'])
    out['sea_chemistry_box_lakes_full_yr'] = sea_chem['sediment']['lakes_full_yr']
    return out


def seas(routing, sea_chem):
    """When the common-level waters fill: the rock the routing sends them, settled at each case's dry density, and
    their own deposits (the sea chemistry's natural seas, 80% or 50% porous, or its held seas' compact deposits)."""
    volume = sea_chem['setting']['sea_kg'] / WATER_DENSITY
    dep = sea_chem['works_at_the_catchment_rate']['deposits_in_the_seas']
    deposits = dict(earliest=dep['natural']['deposit_km3_yr'][1], central=dep['natural']['deposit_km3_yr'][0],
                    latest=dep['managed']['deposit_km3_yr'][0])
    out = {}
    for name, d_km3 in deposits.items():
        f = routing['fill'][name]
        years = np.array([float(y) for y in f['to_common_level_cumulative_gt']])
        rock = np.array(list(f['to_common_level_cumulative_gt'].values())) * 1e12 / f['case']['dry_density_kg_m3']
        filled = rock + d_km3 * 1e9 * years
        if filled[-1] < volume:
            full = None
        else:
            i = int(np.argmax(filled >= volume))
            lo, hi = (years[i - 1], filled[i - 1]) if i else (0.0, 0.0)
            full = float(lo + (years[i] - lo) * (volume - hi) / (filled[i] - hi))
        out[name] = dict(deposits_km3_yr=d_km3, rock_share_at_start=f['to_common_level_gt_yr']['10'] /
                         f['case']['rock_gt_yr'], seas_full_yr=full)
    out['sea_chemistry_box_seas_full_yr'] = sea_chem['sediment']['seas_full_yr']
    out['sea_chemistry_box_passing_share'] = sea_chem['setting']['particle_pass']
    return out


# --- Small craters -------------------------------------------------------------------------------------------------

def crater_profile(r):
    """Pike's fresh simple crater in diameters: a bowl to the rim, ejecta thinning as the cube of distance."""
    inside = -(CRATER_DEPTH - CRATER_RIM) + CRATER_DEPTH * (r / 0.5) ** 2
    return np.where(r <= 0.5, inside, CRATER_RIM * (np.maximum(r, 0.5) / 0.5) ** -3)


def crater_degradation(fractions=(0.5, 0.25, 0.1), dr=1.0 / 400.0, r_max=4.0):
    """kappa t / D^2 at which linear diffusion leaves each fraction of a fresh crater's rim-to-floor depth."""
    r = np.arange(0.0, r_max + 0.5 * dr, dr)
    z = crater_profile(r)
    inner = r <= 0.75
    d0 = z[inner].max() - z[0]
    dt = 0.2 * dr * dr
    t, found = 0.0, {}
    rp, rm = r[1:-1] + 0.5 * dr, r[1:-1] - 0.5 * dr
    while len(found) < len(fractions):
        lap = np.zeros_like(z)
        lap[1:-1] = (rp * (z[2:] - z[1:-1]) - rm * (z[1:-1] - z[:-2])) / (r[1:-1] * dr * dr)
        lap[0] = 4.0 * (z[1] - z[0]) / (dr * dr)
        z = z + dt * lap
        t += dt
        d = z[inner].max() - z[0]
        for f in fractions:
            if f not in found and d <= f * d0:
                found[f] = t
    return {f'{f:g}': v for f, v in found.items()}


def small_craters(routing):
    nd = crater_degradation()
    half = nd['0.5']
    accretion = routing['parameters']['accretion_mm_yr']
    rows = []
    for d in CRATER_DIAMETERS_M:
        depth = CRATER_DEPTH * d
        creep = {f'{k:g}': half * d * d / k for k in KAPPA_M2_YR}
        settle = {k: 0.5 * depth / (v * 1e-3) for k, v in accretion.items()}
        both = {k: 1.0 / (1.0 / creep[f'{KAPPA_M2_YR[i]:g}'] + 1.0 / settle[k])
                for i, k in enumerate(('low', 'central', 'high'))}
        rows.append(dict(diameter_m=d, depth_m=depth,
                         per_km2_larger=[s * GEOMETRIC_SATURATION * (d / 1e3) ** -2 for s in EQUILIBRIUM_SHARE],
                         half_depth_by_creep_yr=creep, half_depth_by_accretion_yr=settle,
                         half_depth_together_yr=both,
                         half_depth_airless_yr=half * d * d / KAPPA_AIRLESS_M2_YR))
    return dict(dimensionless_kappa_t_over_d2=nd, kappa_m2_yr=list(KAPPA_M2_YR),
                equilibrium_share_of_saturation=list(EQUILIBRIUM_SHARE), craters=rows)


# --- Relief --------------------------------------------------------------------------------------------------------

def relief(routing):
    land = routing['land']
    out = {}
    for name, f in routing['fill'].items():
        if name not in ('earliest', 'central', 'latest'):
            continue
        gt = f['case']['rock_gt_yr']
        rows = {}
        for band, b in land['by_band'].items():
            e = b['lowering_m_myr_per_gt_yr'] * gt
            tau = b['mean_relief_m'] / e
            rows[band] = dict(lowering_m_myr=e, lowering_mm_per_century=0.1 * e, relief_m=b['mean_relief_m'],
                              relief_e_folding_myr=tau, relief_half_life_myr=math.log(2.0) * tau)
        out[name] = dict(rock_gt_yr=gt, land_lowering_m_myr=land['lowering_m_myr_per_gt_yr'] * gt, by_band=rows)
    return out


# --- Upkeep and built land -----------------------------------------------------------------------------------------

def upkeep(routing, sea_chem, sea_p, gates):
    c = routing['carrying_back']
    dredge = sea_p['routes']['dredging']
    fuel_mj_kg = dredge['fuel_mj_per_kg_p'][0] / (dredge['volume_m3_per_kg_p'][0] * dredge['fuel_kg_per_m3'][0])
    fuel_kg_m3 = dredge['fuel_kg_per_m3']
    rates = {k: v['rock_gt_per_year'] for k, v in gates['sediment_return'].items()}
    rates['soil loss the phosphorus particles imply'] = sea_chem['setting']['erosion_gt_yr'][1]
    kept = routing['upkeep_at_start']['kept_share']
    heads = c['heads_share_of_land']
    lowering = routing['land']['lowering_m_myr_per_gt_yr']
    rows = {}
    for name, gt in rates.items():
        kg_s = gt * 1e12 / YEAR_S
        lakes_m3 = [gt * 1e12 * kept['lakes'] / rho for rho in (1500.0, 500.0)]
        rows[name] = dict(
            rock_gt_yr=gt,
            lift_gw_in_place=kg_s * MOON_SURFACE_GRAVITY * c['mean_lift_m'] / 1e9,
            lift_gw_to_heads=kg_s * MOON_SURFACE_GRAVITY * c['mean_lift_to_heads_m'] / 1e9,
            haulage_tonne_km_yr=gt * c['tonne_km_yr_per_gt_yr'],
            haulage_gw_at_rail_energy=gt * c['tonne_km_yr_per_gt_yr'] * RAIL_KJ_PER_TONNE_KM * 1e3 / YEAR_S / 1e9,
            share_of_earth_sand_and_gravel=[gt / EARTH_SAND_GRAVEL_GT_YR[1], gt / EARTH_SAND_GRAVEL_GT_YR[0]],
            lake_dredging_km3_yr=[v / 1e9 for v in lakes_m3],
            lake_dredging_fuel_gw=[lakes_m3[0] * fuel_kg_m3[0] * fuel_mj_kg * 1e6 / YEAR_S / 1e9,
                                   lakes_m3[1] * fuel_kg_m3[1] * fuel_mj_kg * 1e6 / YEAR_S / 1e9],
            spread_in_place_mm_per_century=0.1 * lowering * gt,
            spread_on_heads_mm_per_century=0.1 * lowering * gt / heads,
            loess_falls_mm_per_century=list(LOESS_MM_PER_CENTURY))
    works = sea_chem['works_at_the_catchment_rate']
    km3 = works['minerals_km3_yr']
    stone = dict(calcite=km3['calcite'], magnesite=km3['magnesite'], opal=km3['opal'],
                 soda_and_sulphate=km3['soda'] + km3['sodium_sulphate'])
    total = sum(stone.values())
    built = dict(seas_stone_km3_yr=stone, seas_stone_total_km3_yr=total,
                 seas_stone_mt_yr={k: works['minerals_mt_yr'][k] for k in ('calcite', 'magnesite', 'opal')},
                 terraces_10_m_high_km2_per_millennium=total * 1e3 / 0.01,
                 land_rock_km3_yr={k: v['rock_gt_per_year'] * 1e12 / rt.ROCK_DENSITY / 1e9
                                   for k, v in gates['sediment_return'].items()})
    return dict(where_the_rock_settles_at_the_start=kept, dredging_fuel_mj_per_kg=fuel_mj_kg,
                lakes_holding_half_of_the_lakes_catch=routing['upkeep_at_start']['lakes_holding_half_of_it'],
                lakes_holding_nine_tenths=routing['upkeep_at_start']['lakes_holding_nine_tenths'],
                mean_lift_m=c['mean_lift_m'], mean_lift_to_heads_m=c['mean_lift_to_heads_m'],
                heads_share_of_land=heads, mean_distance_km=c['mean_distance_km'], rates=rows, built_land=built)


# --- Travertine ----------------------------------------------------------------------------------------------------

def saturation(m, st, c, law):
    ta = m * law
    sh = st['shares']
    box = dict(mass=1.0, ta=ta, ca=sh['ca'] * ta / 2.0, mg=sh['mg'] * ta / 2.0, mono=(sh['na'] + sh['k']) * ta,
               dic=0.0)
    return cc.speciate(box, st['pco2_atm'], c)


def travertine(routing, sea_chem):
    st = sea_chem['setting']
    c = cc.fresh_constants(cc.LAND_AIR_K)
    law = cc.DESSERT_A * math.exp(cc.DESSERT_B * (cc.LAND_AIR_K - 273.15)) / 1e6       # eq per kg at the catchment rate
    table = {}
    for m in (0.5, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0):
        sat = saturation(m, st, c, law)
        table[f'{m:g}'] = dict(alkalinity_meq_kg=m * law * 1e3, calcium_mmol_kg=st['shares']['ca'] * m * law * 500.0,
                               ph=sat['ph'], omega=sat['omega'])

    def needed(si):
        return math.exp(brentq(lambda x: math.log10(saturation(math.exp(x), st, c, law)['omega']) - si, -3.0, 4.0))
    thresholds = dict(tufa=SI_TUFA, plitvice=SI_PLITVICE, designed=SI_DESIGNED)
    load = {k: needed(si) for k, si in thresholds.items()}
    step = routing['channels']['log_hardness_step']
    rivers = [kv for b in routing['channels']['by_band'].values() for kv in b['by_log_hardness']]
    springs = [kv for b in routing['springs']['land_km2_by_log_hardness'].values() for kv in b]

    def share(hist, threshold):
        total = sum(v for _, v in hist)
        return sum(v for k, v in hist if 10 ** ((k + 0.5) * step) >= threshold) / total if total > 0 else 0.0

    corners = {}
    for i, corner in enumerate(('low', 'central', 'high')):
        p = dict(fresh_factor=cc.FRESH_FACTOR[i], soil_age_yr=cc.SOIL_AGE_AT_START_YR[i])
        rows = {}
        for y in TRAVERTINE_YEARS:
            f = cc.weathering_factor(float(y), p)
            rows[f'{y:g}'] = dict(
                soils_factor=f,
                rivers={k: share(rivers, v / f) for k, v in load.items()},
                spring_land_building_tufa={f'{x:g}': share(springs, load['tufa'] / (f * x)) for x in SPRING_FACTOR})
        corners[corner] = dict(fresh_factor=cc.FRESH_FACTOR[i], soil_age_at_start_yr=cc.SOIL_AGE_AT_START_YR[i],
                               by_year=rows)
    by_band = {b: {k: share(v['by_log_hardness'], t) for k, t in load.items()}
               for b, v in routing['channels']['by_band'].items()}
    return dict(law_alkalinity_meq_kg=law * 1e3, pco2_atm=st['pco2_atm'], air_k=cc.LAND_AIR_K,
                saturation_index=thresholds, spring_factor=list(SPRING_FACTOR), load_needed=load,
                at_load=table, rivers_at_the_catchment_rate_by_band=by_band, corners=corners)


# --- Dam builders --------------------------------------------------------------------------------------------------

def dam_builders(routing):
    """Channel length where dams are built at baseflow and stand through the two-year flood. Earth's thresholds in W
    per m become discharge times slope at equal Shields numbers, (QS)_moon = (QS)_earth sqrt(g_moon/g_earth)
    (conservation's runoff scaling); Earth's own discharge times slope is kept beside it."""
    ch = routing['channels']
    edges = ch['slope_bins_percent']
    step = ch['log_qs_step']
    out = {}
    for limit, (base_w, flood_w) in BRAT_LIMITS_W_M.items():
        for regime, (base_share, flood_ratio) in FLOW_REGIMES.items():
            for rule, scale in (('equal_shields_number', math.sqrt(G_RATIO)), ('earth_discharge_times_slope', 1.0)):
                base_qs = base_w / (WATER_DENSITY * STANDARD_GRAVITY) * scale
                flood_qs = flood_w / (WATER_DENSITY * STANDARD_GRAVITY) * scale
                rows = {}
                for band, b in ch['by_band'].items():
                    km = sum(length for s_bin, q_bin, length in b['by_slope_and_log_qs']
                             if edges[s_bin] < 100.0 * BRAT_SLOPE_MAX
                             and 10 ** ((q_bin + 0.5) * step) * base_share <= base_qs
                             and 10 ** ((q_bin + 0.5) * step) * flood_ratio <= flood_qs)
                    rows[band] = km
                total = sum(rows.values())
                out[f'{limit}/{regime}/{rule}'] = dict(
                    dams_hold_km=total, share_of_channels=total / ch['total_km'], by_band_km=rows,
                    dams={k: [total * v[0], total * v[1]] for k, v in BRAT_DENSITY_PER_KM.items()},
                    pond_km2_at_frequent=[total * d * POND_HA_PER_DAM / 100.0 for d in BRAT_DENSITY_PER_KM['frequent']])
    return dict(limits_w_m=BRAT_LIMITS_W_M, slope_max=BRAT_SLOPE_MAX, flow_regimes=FLOW_REGIMES,
                lunar_scale=math.sqrt(G_RATIO), density_per_km=BRAT_DENSITY_PER_KM, pond_ha_per_dam=POND_HA_PER_DAM,
                channel_km=ch['total_km'], cases=out)


# --- The product ---------------------------------------------------------------------------------------------------

def results():
    data = {p: load(p, s) for p, s in INPUTS.items()}
    routing, sea_chem, sea_p, gates = data[ROUTING], data[SEA_CHEM], data[SEA_P], data[GATES]
    return dict(
        schema=SCHEMA,
        producer=dict(files={p: sha256(ROOT / p) for p in FILES}, inputs={p: sha256(ROOT / p) for p in INPUTS}),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        lakes=lakes(routing, sea_chem),
        seas=seas(routing, sea_chem),
        small_craters=small_craters(routing),
        relief=relief(routing),
        upkeep=upkeep(routing, sea_chem, sea_p, gates),
        travertine=travertine(routing, sea_chem),
        dam_builders=dam_builders(routing))


def dumps(product) -> str:
    return json.dumps(rnd(product), indent=1, allow_nan=False) + '\n'


def main(argv=None) -> int:
    product = results()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(dumps(product))
    print(OUT.relative_to(ROOT))
    return 0


if __name__ == '__main__':
    sys.exit(main())
