"""How much oxygen the rivers carry down into the fresh seas of the design surface: where each river's water goes in
the sea it enters, the oxygen that goes below the nights' mixing with it, and that supply against the deep water's
demand.

    python -m research.studies.river_oxygen.mouths     # the seas and their large rivers (about 20 s)
    python -m research.studies.river_oxygen.waters     # the lake and sea water columns (about 35 minutes)
    python -m research.studies.river_oxygen.run        # results/river_oxygen.json (about 10 s)

A river reaches its sea as the upper water of the lake that set its temperature (mouths.py), carried down a short
course. Every six hours through eight of the design years' ten years (waters.py's columns), this module sets that
water beside the sea at the mouth:

- the river draws its lake's water from the surface down to its own channel's depth, and its temperature at the mouth
  is that water's a travel time earlier, moved towards the lower ground's warmth by the share of the difference a
  river of its depth loses on the way, and warmed by the work of its fall;
- it takes in ENTRAINED volumes of the sea's surface layer for each of its own, and the mixture settles at the first
  depth in the sea column where the sea's water is as dense as it is: within SURFACE_M of the top it spreads at the
  surface, deeper it lies in the upper layers, which the nights mix anyway;
- if it is denser than the whole sea column, 40 m down, it runs on down the floor, below the nights' mixing, and if it
  is denser than the coldest water the sea column ever holds at its foot it is colder than anything the sea makes
  there and goes to the floor of the hollow below the mouth.

The oxygen that goes below the nights' mixing is the river's and the entrained surface water's, each saturated at its
temperature under the air at the sea. It is set against the deep water's demand in the sea oxygen study's nutrient
states (research/studies/sea_oxygen), the areal demand of the floor and the deep water over the sea's area deeper
than DEEP_FROM_M.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np

from biosphere.ecology import lake as L
from research.studies.river_oxygen import waters as W
from shared.constants import MOON_SURFACE_GRAVITY

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
MOUTHS = HERE / 'results' / 'mouths.json'
WATERS = HERE / 'results' / 'waters.json'
SEA_OXYGEN = ROOT / 'research' / 'studies' / 'sea_oxygen' / 'results' / 'sea_oxygen.json'
OUT = HERE / 'results' / 'river_oxygen.json'
SCHEMA = 'terluna.research.river-oxygen/1'
YEAR_S = 31557600.0
HEAT_CAPACITY_J_M3_K = 4.18e6            # water's, by volume
WATER_CP_J_KG_K = 4180.0
# The exchange of heat between a water surface and the air above it, W/m2 per K of the water's departure from the
# temperature it would settle at: read from the lake model's own pans on the design years (a 1 m and a 3 m pan's
# swings through the lunar day at 47 N, 51 E give 24), with a bracket.
EXCHANGE_W_M2_K = (15.0, 24.0, 40.0)
# The river's channel at its mean flow: the river deltas study's mature channels (results on the research drive,
# river-deltas/mouths.csv.gz) follow depth = 0.153 Q^0.535 m and width = 3.875 Q^0.465 m at 1.69 m/s. A river below a
# lake runs near its mean, so the mean flow is taken as the flow that fills the channel; the bracket halves and
# doubles the depth.
CHANNEL_DEPTH, CHANNEL_WIDTH, VELOCITY_M_S = (0.153, 0.535), (3.875, 0.465), 1.69
DEPTH_FACTORS = (0.5, 1.0, 2.0)
# Surface water taken in by each volume of river water where it meets the sea and plunges (sources.json): none, which
# lets the most water sink; a fifth of a volume, the middle of the 0.1-0.3 found where a drowned valley confines the
# inflow; one volume and two, the range at mouths that open on wide water (1.9 at Lake Taupo, about 2 by 400 m from
# the Rhone's mouth in Lake Geneva). These seas' mouths open on wide water, so one volume is the central value.
ENTRAINED = (0.0, 0.2, 1.0, 2.0)
# Oxygen the river's own decaying organic matter would use, g for each m3 of river water: a fifth of 5 g/m3 of
# dissolved organic carbon is degradable (rivers' world mean, and the mean labile share in rivers, sources.json), at
# 2.67 g of oxygen for each g of carbon. Water that has stood years in a lake has spent most of it already, so this is
# an outer figure, reported beside the supply and not taken from it.
RIVER_DEBT_G_M3 = 2.5
# The water a river draws from its lake: from the surface down to its own channel's depth, with the lake's top metre
# alone and its whole column as the bracket.
DRAWN = ('the top metre', 'to the channel\'s depth', 'the whole column')
CENTRAL = dict(entrained=2, exchange=1, depth=1, drawn=1)
SURFACE_M = 5.0                          # water that settles within this depth spreads at the surface
SURFACE_LAYER_M = 10.0                   # the sea's layer a sinking river takes its surface water from
PLUNGE_FROUDE = 0.5                      # densimetric Froude number where an inflow plunges (the river deltas study's)
DEEP_FROM_M = 50.0                       # the sea's floor and water below this depth carry the deep demand
CLASSES = ('at the surface', 'in the upper layers', 'below the nights\' mixing', 'to the floor')
EARLIER_BASINS = dict(near_side_main='procellarum', near_side_south='nubium', humboldtianum='humboldtianum',
                      smythii_marginis='smythii', moscoviense=None, orientale=None)


def sha256(path) -> str:
    digest = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def channel(flow_m3s, depth_factor=1.0):
    depth = CHANNEL_DEPTH[0] * flow_m3s ** CHANNEL_DEPTH[1] * depth_factor
    width = CHANNEL_WIDTH[0] * flow_m3s ** CHANNEL_WIDTH[1]
    return depth, width


def kept_share(km, depth_m, small_lakes_m2_per_m3s, exchange):
    """The share of its source's departure from the local water's temperature that a river still carries at its
    mouth: what its own surface loses over the travel time, and what the small lakes on its way take."""
    travel_s = km * 1e3 / VELOCITY_M_S
    tau_s = HEAT_CAPACITY_J_M3_K * depth_m / exchange
    return float(np.exp(-travel_s / tau_s) / (1.0 + exchange * small_lakes_m2_per_m3s / HEAT_CAPACITY_J_M3_K)), travel_s


def fall_warming_k(fall_m, gravity=MOON_SURFACE_GRAVITY):
    """The warming of water by the work of its fall, g h / c."""
    return gravity * fall_m / WATER_CP_J_KG_K


def layer_mean(profile_c, top_m, bottom_m, depth_m):
    """The mean temperature of a column's water from its surface down to depth_m, at each kept time."""
    thick = np.clip(np.minimum(bottom_m, depth_m) - top_m, 0.0, None)
    thick[0] = max(thick[0], 1e-9)                                              # a draw thinner than the top layer
    return profile_c @ (thick / thick.sum())


def river_temperature(drawn_c, steps_late, kept, fall_m, lapse_k_m):
    """The river's temperature at its mouth, time by time: the water drawn from the lake steps_late earlier, plus the
    share of the warming down to the mouth's level that it has taken on, plus the work of its fall."""
    i = np.maximum(np.arange(drawn_c.size) - steps_late, 0)
    return drawn_c[i] + (1.0 - kept) * lapse_k_m * fall_m + fall_warming_k(fall_m)


def classes(river_c, sea_profile_c, top_m, bottom_m, coldest_foot_c, entrained, salinity=0.0):
    """Each time's class (an index into CLASSES), the mixture's temperature and the depth it settles at in the sea
    column (nan where it is denser than the whole column)."""
    rho = lambda t: L.density(t, salinity)                                     # noqa: E731
    surface = layer_mean(sea_profile_c, top_m, bottom_m, SURFACE_LAYER_M)
    mixed = (river_c + entrained * surface) / (1.0 + entrained)
    as_dense = rho(sea_profile_c) >= rho(mixed)[:, None]
    settles = as_dense.any(axis=1)
    depth = np.where(settles, top_m[np.argmax(as_dense, axis=1)], np.nan)
    out = np.where(settles, np.where(depth < SURFACE_M, 0, 1), 2).astype(np.int8)
    out[~settles & (rho(mixed) >= rho(coldest_foot_c))] = 3
    return out, mixed, depth, surface


def plunge_depth(flow_m3s, width_m, river_c, surface_c, gravity=MOON_SURFACE_GRAVITY, salinity=0.0):
    """Depth of water an inflow needs before it plunges, (q^2 / (g' Fp^2))^(1/3), where it is denser than the surface."""
    excess = (L.density(river_c, salinity) - L.density(surface_c, salinity)) / L.density(surface_c, salinity)
    g_reduced = gravity * np.maximum(excess, 1e-12)
    return ((flow_m3s / width_m) ** 2 / (g_reduced * PLUNGE_FROUDE ** 2)) ** (1.0 / 3.0)


def load_column(spec):
    with np.load(W.COLUMNS / spec['file']) as z:
        return {k: np.asarray(z[k], dtype=float) for k in z.files}


def one_river(m, pair, columns, sea_level, lapse_k_m, setting, entrained, exchange, depth_factor, drawn):
    lake, sea = load_column(columns[pair['lake']]), load_column(columns[pair['sea']])
    course, src = m['course'], m['course']['source']
    depth, width = channel(m['flow_m3s'], depth_factor)
    kept, travel_s = kept_share(course['km'], depth, course['small_lakes_m2_per_m3s'], exchange)
    fall = src['level_m'] - sea_level
    every = W.PROFILE_EVERY
    draw_m = {0: 1.0, 1: depth, 2: float(lake['layer_bottom_m'][-1])}[drawn]
    from_lake = layer_mean(lake['profile_t'], lake['layer_top_m'], lake['layer_bottom_m'], draw_m)
    river = river_temperature(from_lake, int(round(travel_s / (lake['dt_s'] * every))), kept, fall, lapse_k_m)
    foot = sea['profile_t'][:, -1]
    coldest = float(sea['bottom_t'].min())
    cls, mixed, settles_m, surface = classes(river, sea['profile_t'], sea['layer_top_m'], sea['layer_bottom_m'], coldest,
                                             entrained, setting['salinity_g_kg'])
    p_atm = sea['pressure_pa'][::every][:river.size] / 101325.0
    o2 = lambda t: L.oxygen_saturation(t, p_atm, setting['salinity_g_kg'], setting['o2_fraction'])   # noqa: E731
    carried = (o2(river) + entrained * o2(surface)) / 1.0                       # g O2 with each m3 of river water
    shares = [float((cls == k).mean()) for k in range(len(CLASSES))]
    down = cls >= 2
    floor = cls == 3
    sinking = L.density(river, setting['salinity_g_kg']) > L.density(surface, setting['salinity_g_kg'])
    out = dict(
        flow_m3s=m['flow_m3s'], lat=m['lat'], lon_east=m['lon_east'],
        source=dict(km=course['km'], above_the_sea_m=fall, area_km2=src.get('area_km2'), mean_depth_m=src.get('mean_depth_m'),
                    column_depth_m=columns[pair['lake']]['depth_m'], drawn_to_m=min(draw_m, float(lake['layer_bottom_m'][-1]))),
        channel=dict(depth_m=depth, width_m=width), kept_share_of_the_sources_departure=kept,
        travel_days=travel_s / 86400.0,
        river_c=dict(mean=float(river.mean()), lowest=float(river.min()), highest=float(river.max())),
        sea_surface_layer_c=dict(mean=float(surface.mean()), lowest=float(surface.min()), highest=float(surface.max())),
        sea_foot_c=dict(mean=float(foot.mean()), coldest=coldest),
        denser_than_the_seas_surface_layer_share=float(sinking.mean()),
        river_less_sea_surface_k=float((river - surface).mean()),
        river_less_sea_foot_k=dict(mean=float((river - foot).mean()), lowest=float((river - foot).min()),
                                   highest=float((river - foot).max())),
        settles_at_m=(dict(zip(('p10', 'median', 'p90'), np.nanpercentile(settles_m, [10, 50, 90]).tolist()))
                      if np.isfinite(settles_m).any() else None),
        shares=dict(zip(CLASSES, shares)),
        oxygen_below_the_nights_mixing_t_yr=float((carried * down).mean() * m['flow_m3s'] * YEAR_S / 1e6),
        oxygen_to_the_floor_t_yr=float((carried * floor).mean() * m['flow_m3s'] * YEAR_S / 1e6),
        water_below_the_nights_mixing_km3_yr=float(down.mean() * (1.0 + entrained) * m['flow_m3s'] * YEAR_S / 1e9),
        water_to_the_floor_km3_yr=float(floor.mean() * (1.0 + entrained) * m['flow_m3s'] * YEAR_S / 1e9),
        river_water_below_the_nights_mixing_km3_yr=float(down.mean() * m['flow_m3s'] * YEAR_S / 1e9),
        floor_water_colder_than_the_seas_coldest_k=float((coldest - mixed[floor]).mean()) if floor.any() else None,
        river_oxygen_g_m3=float(o2(river).mean()), sea_surface_oxygen_g_m3=float(o2(surface).mean()),
        plunge_depth_m=(float(np.median(plunge_depth(m['flow_m3s'], width, river[sinking], surface[sinking])))
                        if sinking.any() else None),
        floor_reaches_m=dict(first_lowest=m['floor']['first_lowest_depth_m'], last_hollow=m['floor']['lowest_depth_m'],
                             hollows_passed=m['floor']['hollows_passed']),
        pool_km2={str(p['height_m']): p['area_km2'] for p in m['floor']['pool']},
        columns=dict(lake=pair['lake'], sea=pair['sea'], lake_on_land_cell=pair['lake_cell_of_its_kind'],
                     sea_on_water_cell=pair['sea_cell_of_its_kind']))
    return out


def demand_g_m2_yr(setting):
    """The deep water's areal oxygen demand in each nutrient state, the floor's and the deep water's together."""
    return {k: (v['sediment'] + v['deep_water']) * 365.25 for k, v in setting['trophic'].items()}


def one_sea(name, mouths, waters, setting, lapse_k_m, entrained, exchange, depth_factor, drawn):
    sea, group = mouths['seas'][name], mouths['rivers'][name]
    rivers = [one_river(m, pair, waters['columns'], sea['level_m'], lapse_k_m, setting, entrained, exchange, depth_factor, drawn)
              for m, pair in zip(group['kept'], waters['rivers'][name])]
    kept = sum(r['flow_m3s'] for r in rivers)
    total = sea['rivers_m3s']
    scale = total / kept                                   # the small rivers take the kept rivers' shares
    shares = {c: sum(r['shares'][c] * r['flow_m3s'] for r in rivers) / kept for c in CLASSES}
    below = next(b for b in sea['below'] if b['depth_m'] == DEEP_FROM_M)
    oxygen_down = sum(r['oxygen_below_the_nights_mixing_t_yr'] for r in rivers) * scale
    oxygen_floor = sum(r['oxygen_to_the_floor_t_yr'] for r in rivers) * scale
    water_down = sum(r['water_below_the_nights_mixing_km3_yr'] for r in rivers) * scale
    debt = sum(r['river_water_below_the_nights_mixing_km3_yr'] for r in rivers) * scale * RIVER_DEBT_G_M3 * 1e3   # t a year
    demand = {k: v * below['area_km2'] for k, v in demand_g_m2_yr(setting).items()}          # t a year (g/m2 * km2)
    mean_o2 = sum(r['river_oxygen_g_m3'] * r['flow_m3s'] for r in rivers) / kept
    all_of_it = total * YEAR_S * mean_o2 / 1e6                                               # t a year
    return dict(
        rivers_m3s=total, rivers_km3_yr=total * YEAR_S / 1e9, area_km2=sea['area_km2'], volume_km3=sea['volume_km3'],
        rivers_share_of_volume_a_year=sea['rivers_share_of_volume_a_year'],
        deep_area_km2=below['area_km2'], deep_volume_km3=below['volume_km3'] + DEEP_FROM_M * below['area_km2'] / 1e3,
        shares_of_river_water=shares,
        water_below_the_nights_mixing_km3_yr=water_down,
        years_to_replace_the_deep_water=(below['volume_km3'] / water_down) if water_down > 0 else None,
        oxygen_below_the_nights_mixing_t_yr=oxygen_down, oxygen_to_the_floor_t_yr=oxygen_floor,
        the_rivers_own_organic_matter_would_use_t_yr=debt,
        demand_t_yr=demand,
        share_of_demand_met={k: oxygen_down / v for k, v in demand.items()},
        if_every_river_sank=dict(oxygen_t_yr=all_of_it, share_of_demand_met={k: all_of_it / v for k, v in demand.items()},
                                 with_two_volumes_of_surface_water={k: 3.0 * all_of_it / v for k, v in demand.items()}),
        floor_km2_the_supply_covers={k: oxygen_down / v for k, v in demand_g_m2_yr(setting).items()},
        earlier_basin=EARLIER_BASINS[name], rivers=rivers)


def build():
    mouths, waters = json.loads(MOUTHS.read_text()), json.loads(WATERS.read_text())
    if waters['columns_not_yet_run']:
        raise SystemExit(f'{len(waters["columns_not_yet_run"])} columns are not yet run: python -m research.studies.river_oxygen.waters')
    for spec in waters['columns'].values():
        if sha256(W.COLUMNS / spec['file']) != spec['sha256']:
            raise SystemExit(f'{spec["file"]} does not match results/waters.json')
    sea_oxygen = json.loads(SEA_OXYGEN.read_text())
    setting = sea_oxygen['setting']
    lapse = waters['lapse_k_per_km'] / 1e3
    names = [n for n, s in mouths['seas'].items() if not s['salt']]
    cases = {}
    for label, choice in (('central', CENTRAL),
                          ('no surface water taken in', dict(CENTRAL, entrained=0)),
                          ('a fifth of a volume of surface water taken in', dict(CENTRAL, entrained=1)),
                          ('two volumes of surface water taken in', dict(CENTRAL, entrained=3)),
                          ('weak exchange with the air', dict(CENTRAL, exchange=0)),
                          ('strong exchange with the air', dict(CENTRAL, exchange=2)),
                          ('shallow channels', dict(CENTRAL, depth=0)),
                          ('deep channels', dict(CENTRAL, depth=2)),
                          ('the lake\'s top metre drawn', dict(CENTRAL, drawn=0)),
                          ('the lake\'s whole column drawn', dict(CENTRAL, drawn=2)),
                          ('no surface water and the top metre drawn', dict(CENTRAL, entrained=0, drawn=0)),
                          ('no surface water and the whole column drawn', dict(CENTRAL, entrained=0, drawn=2))):
        args = (ENTRAINED[choice['entrained']], EXCHANGE_W_M2_K[choice['exchange']], DEPTH_FACTORS[choice['depth']],
                choice['drawn'])
        cases[label] = {n: one_sea(n, mouths, waters, setting, lapse, *args) for n in names}
        cases[label]['settings'] = dict(entrained=args[0], exchange_w_m2_k=args[1], channel_depth_factor=args[2],
                                        drawn=DRAWN[args[3]])
    central = cases['central']
    fresh = dict(rivers_m3s=sum(central[n]['rivers_m3s'] for n in names),
                 area_km2=sum(central[n]['area_km2'] for n in names),
                 volume_km3=sum(central[n]['volume_km3'] for n in names))
    for label, case in cases.items():
        down = sum(case[n]['oxygen_below_the_nights_mixing_t_yr'] for n in names)
        demand = {k: sum(case[n]['demand_t_yr'][k] for n in names) for k in setting['trophic']}
        case['all_fresh_seas'] = dict(
            oxygen_below_the_nights_mixing_t_yr=down, demand_t_yr=demand,
            share_of_demand_met={k: down / v for k, v in demand.items()},
            shares_of_river_water={c: sum(case[n]['shares_of_river_water'][c] * case[n]['rivers_m3s'] for n in names)
                                   / fresh['rivers_m3s'] for c in CLASSES},
            if_every_river_sank={k: sum(case[n]['if_every_river_sank']['oxygen_t_yr'] for n in names) / v
                                 for k, v in demand.items()})
    for label, case in cases.items():
        if label != 'central':
            for n in names:
                case[n].pop('rivers')
    below = lambda sea: sea['shares_of_river_water'][CLASSES[2]] + sea['shares_of_river_water'][CLASSES[3]]   # noqa: E731
    ranges = {}
    for n in names + ['all_fresh_seas']:
        shares = {label: below(case[n]) for label, case in cases.items()}
        met = {label: case[n]['share_of_demand_met']['oligotrophic'] for label, case in cases.items()}
        ranges[n] = dict(share_of_river_water_below_the_nights_mixing=dict(
                             least=min(shares.values()), most=max(shares.values()), most_in=max(shares, key=shares.get)),
                         share_of_oligotrophic_demand_met=dict(
                             least=min(met.values()), most=max(met.values()), most_in=max(met, key=met.get)))
    deep = {}
    for n in names:
        b = EARLIER_BASINS[n]
        if b:
            a = sea_oxygen['answers'][b]['typical']
            deep[n] = dict(basin=b, deep_age_years=a['deep_age_years'], deep_c=a['deep_t_c'],
                           deepest_tenth_oxygen_g_m3=({k: v['typical']['mean']
                                                       for k, v in sea_oxygen['bounds'][b]['oxygen'].items()}
                                                      if b in sea_oxygen['bounds'] else None))
    return dict(fresh=fresh, cases=cases, ranges=ranges, wind_mixing_in_the_sea_oxygen_study=deep, setting=setting, names=names,
                demand_g_m2_yr=demand_g_m2_yr(setting), lapse_k_per_km=waters['lapse_k_per_km'],
                design_years=waters['design_years'])


def write(result, out=OUT, seconds=None):
    record = dict(
        schema=SCHEMA,
        producer=dict(files={str(Path(f).resolve().relative_to(ROOT)): sha256(f) for f in (__file__, W.__file__, L.__file__)}),
        inputs={str(p.relative_to(ROOT)): sha256(p) for p in (MOUTHS, WATERS, SEA_OXYGEN)},
        evidence=('A calculation on products and model results; nothing is measured. The rivers\' flows and courses '
                  'are the design routing\'s on the estimated design climate; the temperatures are runs of the '
                  'project\'s lake model on the earlier design years\' weather at T21, one column for a lake and one '
                  'for the sea beside the mouth; the mixing a sinking river does is a bracket from Earth\'s lakes and '
                  'reservoirs; the demand is the sea oxygen study\'s, Earth lakes\' areal rates. The classes say where '
                  'water of a given temperature would go in a sea layered as the sea column is; no flow is computed, '
                  'and the sea column does not feel the river.'),
        reading_rule=('Flows m3/s, water km3 a year, oxygen t a year, temperatures C, depths m, areas km2. cases[label] '
                      'holds each fresh sea under one set of choices (cases[label].settings); cases.central also lists '
                      'each river. shares_of_river_water are the shares of the year\'s river water by where it goes: '
                      'at the surface (settling within 5 m of it), in the upper layers the nights mix (within the 40 m '
                      'sea column), below the nights\' mixing (denser than the whole column), and of that last the part '
                      'colder than the coldest water the sea column ever holds at its foot (to the floor). '
                      'oxygen_below_the_nights_mixing_t_yr counts the last two classes, with the surface water they took '
                      'in. demand_t_yr is the areal demand of each nutrient state over the sea\'s area deeper than 50 m, '
                      'and share_of_demand_met their ratio. if_every_river_sank puts all the river water below, alone '
                      'and with two volumes of surface water. the_rivers_own_organic_matter_would_use_t_yr is an '
                      'outer figure for what the sinking river water\'s own decay would take, not subtracted. '
                      'range_over_the_cases gives each sea\'s least and most over every case. '
                      'floor_km2_the_supply_covers is the floor whose demand '
                      'the supply equals. A river\'s pool_km2 gives the hollow below its mouth: the joined water '
                      'within each height above its lowest point.'),
        classes=list(CLASSES),
        settings=dict(entrained=list(ENTRAINED), drawn=list(DRAWN), surface_m=SURFACE_M, surface_layer_m=SURFACE_LAYER_M,
                      river_debt_g_m3=RIVER_DEBT_G_M3,
                      exchange_w_m2_k=list(EXCHANGE_W_M2_K), channel_depth=list(CHANNEL_DEPTH),
                      channel_width=list(CHANNEL_WIDTH), velocity_m_s=VELOCITY_M_S, depth_factors=list(DEPTH_FACTORS),
                      plunge_froude=PLUNGE_FROUDE, deep_from_m=DEEP_FROM_M, gravity_m_s2=MOON_SURFACE_GRAVITY,
                      fall_warming_k_per_km=fall_warming_k(1000.0), lapse_k_per_km=result['lapse_k_per_km'],
                      salinity_g_kg=result['setting']['salinity_g_kg'], o2_fraction=result['setting']['o2_fraction'],
                      design_years=result['design_years']),
        demand_g_m2_yr=result['demand_g_m2_yr'], fresh_seas=result['fresh'],
        wind_mixing_in_the_sea_oxygen_study=result['wind_mixing_in_the_sea_oxygen_study'],
        range_over_the_cases=result['ranges'], cases=result['cases'])
    if seconds is not None:
        record['seconds'] = seconds
    out.write_text(json.dumps(record, indent=1) + '\n')
    return record


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument('--out', type=Path, default=OUT)
    args = p.parse_args(argv)
    started = time.time()
    record = write(build(), args.out, round(time.time() - started))
    c = record['cases']['central']
    for n in record['cases']['central']:
        if n in ('settings', 'all_fresh_seas'):
            continue
        s = c[n]
        print(f'{n}: rivers {s["rivers_m3s"]:,.0f} m3/s; shares ' + ', '.join(f'{k} {100 * v:.1f}%' for k, v in s['shares_of_river_water'].items())
              + f'; oxygen below {s["oxygen_below_the_nights_mixing_t_yr"]:,.0f} t/yr; of the demand '
              + ', '.join(f'{k} {100 * v:.2f}%' for k, v in s['share_of_demand_met'].items()))
    print(f'wrote {args.out}, {time.time() - started:.0f} s')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
