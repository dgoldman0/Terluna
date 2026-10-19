"""What could take pressure off the shield's dimming: each lever on the greenhouse effect, the plan's added heat or
the Moon's own reflection, in W/m2 of the Moon's balance and in points of dimming against the hot world's absorbed
sunlight, with its cost and side effects from the project's products.

    python -m research.studies.heat_levers.run        # under a second

A screening on committed products: the target climates, the deep column's gas forcing and CO2 and pressure sweeps,
the domain synthesis's gas and heat ledgers, the food study's diet, the carbon cycle's import, the power screen, the
fleet's glow (read only), the CO2 study's leaf model, the GCM's clear-sky radiation check and climatology, and the
climate targets' extract for the deserts by band, the corrected targets and the extract they bind. No model is run.
"""
from __future__ import annotations
import hashlib
import io
import json
import math
import sys
from pathlib import Path

import numpy as np

from shared import constants as K

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'heat_levers.json'
README = HERE / 'README.md'
SCHEMA = 'terluna.research.heat-levers/2'
CODE = ('research/studies/heat_levers/run.py',)
TARGETS = 'research/studies/climate_targets/results/climate_targets.json'
GASES = 'atmosphere/radiative_convective/results/reduced_gases.json'
COLUMN = 'atmosphere/radiative_convective/results/summary.json'
LEDGERS = 'research/studies/domain_synthesis/results/ledgers.json'
FOOD = 'research/studies/food_systems/results/food_systems.json'
CARBON = 'research/studies/carbon_cycle/results/carbon_cycle.json'
POWER = 'provisioning/results/power.json'
GLOW = 'research/studies/solar_shield_array/results/glow.json'
CO2 = 'research/studies/atmospheric_co2/results/atmospheric_co2.json'
RADIATION = 'climate/results/gcm/radiation_check_A28_dim5_moon.json'
CLIMATOLOGY = 'climate/gcm/products/climatology_A28_corrected_970.npz'
INDUSTRY = 'resources/results/industry.json'
SCREEN = 'resources/results/first_screen.json'
EXTRACT = 'research/runs/climate_targets/extract.json'     # ignored, on the research drive; deserts by band
INPUTS = {TARGETS: 'terluna.research.climate-targets/5', GASES: 'terluna.atmosphere.reduced-gas-forcing/1',
          COLUMN: None, LEDGERS: 'terluna.research.domain-synthesis-ledgers/1', FOOD: 'terluna.research.food-systems/1',
          CARBON: 'terluna.research.carbon-cycle/1', POWER: 'terluna.provisioning.power/1',
          GLOW: 'terluna.research.fleet-glow/1', CO2: 'terluna.research.atmospheric-co2/1',
          RADIATION: 'terluna.climate.gcm-radiation-check/1', CLIMATOLOGY: 'terluna.climate.gcm-climatology/1',
          INDUSTRY: 'terluna.resources.industry/1', SCREEN: 'terluna.resources.first-screen/1',
          EXTRACT: 'terluna.research.climate-targets-extract/5'}
RESORT_K = 301.0                   # the author's resort, 18 October (the register)

# --- published values (sources.json) ---
AR6_ERF_2XCO2 = 3.93               # W/m2, AR6 WG1 ch. 7, effective forcing of doubled CO2
METHANOGENESIS_EV = 0.96           # Yvon-Durocher et al. 2014, seasonal temperature dependence of CH4 emissions
CIRRUS_EARTH_NET_W_M2 = 5.0        # AR6 WG1 ch. 4: Earth's cirrus warm the planet by about 5 W/m2 on net
CCT_EARTH_MAX_W_M2 = (1.0, 2.0)    # AR6 WG1 ch. 4: the most cirrus thinning could give, low confidence
CCT_OVERSEEDING_W_M2 = 4.9         # Tully et al. 2022: warming at the highest seeding concentration
EARTH_F_A = dict(clear=0.73, all_sky=0.48)   # Lenton & Vaughan 2009: two-way transmittance for surface albedo
EARTH_CLOUD_COVER = 0.67           # design guess near ISCCP's global mean, to carry Earth's cloud masking to the Moon
ELEVATOR_T_PER_YEAR = 584.0        # Pearson 2005 in provisioning's register: a lunar elevator of today's materials
# --- design guesses, each labelled in the README ---
HOLD_CO2_PPM = (280.0, 350.0)      # CO2 levels screened below the design's 400 ppm
CLOUD_COVER = dict(hot=0.383, lock=0.198)    # the GCM's cover, climate/gcm/README.md, years 10-19 of the two runs
CITY_SUN_WEIGHT = (0.6, 1.0)       # settled land's sunlight over the Moon's mean: highlands and poleward of 48.6 deg
CITY_ALBEDO_RISE = (0.10, 0.15)    # roofs and pavements (Menon et al. 2010; Lenton & Vaughan's settlements)
PARK_ALBEDO_RISE = (0.02, 0.05)    # silver foliage in the hybrid parks (crop albedo +0.04 in Ridgwell et al. 2009)
PARK_SUN_WEIGHT = (0.6, 1.0)
DESERT_ALBEDO_RISE = (0.10, 0.30)  # gravel, crusts or sheeting on the deserts; Lenton & Vaughan take 0.44
HIGHLAND_M = 4000.0                # land columns from this height up leave the albedo kernel (the targets' 4 km class)
DESERT_MM_DAY = 0.5                # the targets' desert: land under 0.5 mm of rain a day
INP_PER_LITRE = (1.0, 20.0)        # seeding particles in the cirrus layer (Gasparini et al. 2020: 1 and 18)
INP_DIAMETER_M = 0.5e-6
INP_DENSITY_KG_M3 = 2730.0         # plagioclase feldspar, the highlands' own mineral
INP_LIFE_DAYS = 10.0
CIRRUS_LAYER_KM = (63.0, 107.0)    # the GCM's two cloudy top layers after the longwave refit
DIET_KG_C_PER_DAY = 0.25           # carbon a person eats in a day, about 600 g of dry food
WHOLE_MODES = ('power_made', 'fully_lamp_lit')   # the blend's mixes made wholly with power; power_heavy in part
HELD_CUT = '0.9'                   # biology's cut at the methane sources, as the synthesis holds the pair

EVIDENCE = ('A screening on committed products. Forcings come from the deep column (line-by-line, clear sky, '
            'prescribed humidity; all-sky by AR6 ratios on its Earth control), the cloud effects and absorbed '
            'sunlight from two ExoPlaSim runs, costs from the domain synthesis, the food study, the carbon cycle, the '
            'power screen and the glow study (read only), surface-albedo kernels by latitude from the GCM radiation '
            'check\'s clear-sky line-by-line shares with Earth\'s cloud masking scaled by cover, and the cirrus, '
            'methane and albedo responses from published Earth studies. No lever is run in a climate model; each W/m2 '
            'is treated as absorbed sunlight at efficacy 1 except the glow, which takes the GCM\'s 0.645-0.71. The '
            'climate targets and their extract are the corrected ones.')
READING_RULE = ('W/m2 are means over the whole Moon, positive for what a lever takes off the balance. A point of '
                'dimming is 1% of the hot world\'s absorbed sunlight, 225.4 W/m2, as climate_targets counts the heat '
                'beyond the sunlight. Pairs are [low, high]. The resort is the author\'s 301 K, read from '
                'climate_targets\' resort_band. TW of heat on the Moon count 0.0264 W/m2 each, at the ring radius '
                '5e-5. Values marked estimate rest on the design guesses named in run.py; ranking pairs each lever\'s '
                'points as the README\'s table gives them.')


def raw(path) -> bytes:
    """An input's bytes, from the tree."""
    return (ROOT / path).read_bytes()


def digest(path) -> str:
    return hashlib.sha256(raw(path) if path in INPUTS else (ROOT / path).read_bytes()).hexdigest()


def load():
    data = {}
    for p, schema in INPUTS.items():
        if p.endswith('.npz'):
            z = np.load(io.BytesIO(raw(p)), allow_pickle=False)
            meta = json.loads(str(z['metadata']))
            data[p] = dict(meta=meta, **{k: np.asarray(z[k], float) for k in ('lat', 'lsm', 'clt_mean', 'pr_mm_day')})
            got = meta.get('schema')
        else:
            data[p] = json.loads(raw(p))
            got = data[p].get('schema') if isinstance(data[p], dict) else None
        if schema is not None and got != schema:
            raise ValueError(f'{p}: schema {got} is not {schema}')
    if data[TARGETS]['producer']['inputs'][EXTRACT] != digest(EXTRACT):
        raise ValueError('the extract is not the one the climate targets bind')
    return data


def pair(lo, hi):
    return sorted([float(lo), float(hi)])


class Basis:
    def __init__(self, data):
        t = data[TARGETS]
        self.absorbed = t['runs']['test']['top_of_air_w_m2']['absorbed']
        self.incoming = t['runs']['test']['top_of_air_w_m2']['incoming']
        self.area = 4.0 * math.pi * K.MOON_RADIUS ** 2
        ind = data[INDUSTRY]['places']
        self.moon_per_tw = ind['moon_surface']['heat_w_m2_per_tw_used']
        self.ring_per_tw = ind['ring_platforms']['heat_w_m2_per_tw_even']
        self.resort = next(r for r in t['resort_band']['rows'] if r['ground_k'] == RESORT_K)
        named = {x['name']: x for x in t['targets']}
        self.dimming = {'hot world': named['no further dimming']['with_the_heat_beyond_sunlight_percent'],
                        'resort, 301 K': self.resort['dimming']['with_the_heat_beyond_sunlight_percent'],
                        'lock, 295 K': named['295 K']['with_the_heat_beyond_sunlight_percent']}

    def points(self, w):
        if isinstance(w, (list, tuple)):
            return [100.0 * x / self.absorbed for x in w]
        return 100.0 * w / self.absorbed


def masking(cover):
    """Earth's cloud masking of the surface's reflection (Lenton & Vaughan's all-sky over clear-sky two-way
    transmittance at Earth's cover), carried to another cover: the share of a clear-sky change that reaches the top
    of the air."""
    cloudy = 1.0 - (1.0 - EARTH_F_A['all_sky'] / EARTH_F_A['clear']) / EARTH_CLOUD_COVER
    return 1.0 - cover * (1.0 - cloudy)


def co2_all_sky_ratio(data):
    """AR6's effective forcing of doubled CO2 over the Earth control's own clear-sky change per doubling."""
    return AR6_ERF_2XCO2 / -data[COLUMN]['earth_1.0atm']['co2_at_290k']['olr_change_per_doubling']


def gases(data, b):
    """(a) Methane and N2O below the held pair: the forcing saved, what holding them takes, and the warm Moon's
    lakes."""
    g = data[GASES]
    m = g['moon_at_design']['both']
    held_key = '1.8663+0.1'
    held = m[held_key]['allsky_w_m2']
    options = {'low pair, 0.3 and 0.03 ppm': '0.3+0.03', 'methane low, 0.3 and 0.1 ppm': '0.3+0.1',
               'N2O low, 1.87 and 0.03 ppm': '1.8663+0.03', 'lowest grid pair, 0.1 and 0.01 ppm': '0.1+0.01'}
    # each saving's own clear-sky trend with the columns' warmth, carried to the hot world's ground linearly
    col = g['columns']['moon_12atm']
    temps = sorted(float(k) for k in col)
    hot_k = data[TARGETS]['runs']['test']['ground_k']

    def to_the_hot_ground(values):
        slope = float(np.polyfit(temps, values, 1)[0])
        base = values[temps.index(295.0)]
        return (base + slope * (hot_k - 295.0)) / base
    toa = [col[f'{t:g}']['both'][held_key]['toa'] for t in temps]
    lw = data[LEDGERS]
    w_per_tw = lw['heat']['factors']['moon_w_m2_per_tw']
    out = dict(held_pair_allsky_w_m2=held, held_pair_points=b.points(held),
               forcing_at_the_hot_worlds_ground_factor=to_the_hot_ground(toa), column_temperatures_k=temps,
               held_pair_clear_sky_toa_w_m2=toa)
    opts = {}
    for name, key in options.items():
        saved = held - m[key]['allsky_w_m2']
        clear = [col[f'{t:g}']['both'][held_key]['toa'] - col[f'{t:g}']['both'][key]['toa'] for t in temps]
        factor = to_the_hot_ground(clear)
        opts[name] = dict(allsky_w_m2=m[key]['allsky_w_m2'], clear_sky_saving_by_column_w_m2=clear,
                          saving_at_the_hot_worlds_ground_factor=factor, saved_w_m2=pair(saved * factor, saved),
                          points=b.points(pair(saved * factor, saved)))
    out['options'] = opts
    pw = lw['gases']['pairs_weighed']

    def costs(p):
        c = {}
        for ch4_cut, n2o_cut in (('0.9', '0.95'), ('0.5', '0.5')):
            ch4, n2o = p['methane'][ch4_cut], p['n2o'][n2o_cut]
            tw = [f + d for f, d in zip(ch4['fan_tw'], n2o['decomposition_tw'])]
            c[f'methane cut {ch4_cut}, N2O cut {n2o_cut}'] = dict(
                methane_air_share_per_year=ch4['atmosphere_share_per_year'], methane_fan_tw=ch4['fan_tw'],
                n2o_air_share_per_year=n2o['atmosphere_share_per_year'], n2o_decomposition_tw=n2o['decomposition_tw'],
                removal_heat_w_m2=[x * w_per_tw for x in tw], removal_points=b.points([x * w_per_tw for x in tw]),
                methane_years_to_reach=ch4['years_to_reach_from_none'],
                n2o_years_to_reach=n2o['years_to_reach_from_none'])
        return c
    out['costs'] = dict(held=costs(pw['middle']), low_pair=costs(pw['low']))
    # air through the low pair's methane contactors, as a flow
    air_kg = data[SCREEN]['inventories']['air_kg']
    year = K.JULIAN_YEAR_DAYS * K.JULIAN_DAY
    share = pw['low']['methane']['0.9']['atmosphere_share_per_year']
    out['low_pair_methane_contactor_air_kg_s'] = [x * air_kg / year for x in share]
    out['warm_waters'] = warm_waters(data, b)
    return out


def warm_waters(data, b):
    """The warm Moon's waters, an indicator: the wild sources the biology leaves, scaled by methane emissions'
    seasonal temperature dependence between the lock's sea median and the resort's and the hot world's, the farms'
    and herds' methane and the soils' uptake held."""
    kb_ev = K.BOLTZMANN / K.ELEMENTARY_CHARGE
    t = data[TARGETS]
    seas = {'lock': t['runs']['design']['sea_median_c'], 'resort': b.resort['comfort']['sea_median_c'],
            'hot': t['runs']['test']['sea_median_c']}

    def arrhenius(c_from, c_to):
        return math.exp(METHANOGENESIS_EV / kb_ev * (1.0 / (c_from + 273.15) - 1.0 / (c_to + 273.15)))
    g = data[LEDGERS]['gases']
    held = g['held']['methane_by_source_cut'][HELD_CUT]
    left, removal, ppm = held['sources_left_tg_yr'], held['removal_tg_yr'], held['steady_ppm_by_soils_alone']
    # the held pair's two ends are the synthesis's central wild sources with the farms' low and high methane
    cases = g['methane']['central']
    if not all(math.isclose(x, (1.0 - float(HELD_CUT)) * c['sources_tg_yr']['total'], rel_tol=1e-3)
               for x, c in zip(left, cases)):
        raise ValueError('the held pair\'s sources are not the synthesis\'s central cases at the cut')
    wild = [x * c['sources_tg_yr']['wild'] / c['sources_tg_yr']['total'] for x, c in zip(left, cases)]
    share_per_tg = [s / r for s, r in zip(held['atmosphere_share_per_year'], removal)]
    # ecology's waters' oxygen as the targets bind it: the oligotrophic deep lakes' deepest tenth, mg/L
    deep = t['resort_band']['oxygen_source']['deep']
    oxygen = dict(lock=deep['lock'], resort=b.resort['oxygen']['deep_oligotrophic_deepest_tenth_mg_l'], hot=deep['hot'])
    out = dict(sea_median_c=seas, deep_oligotrophic_deepest_tenth_mg_l=oxygen,
               oxygen_source={k: t['resort_band']['oxygen_source'][k] for k in ('path', 'branch', 'commit', 'sha256')},
               activation_ev=METHANOGENESIS_EV, wild_sources_left_tg_yr=wild,
               sources_left_tg_yr=left,
               lock=dict(removal_tg_yr=removal, air_share_per_year=held['atmosphere_share_per_year'],
                         soils_alone_steady_ppm=ppm))
    for name in ('resort', 'hot'):
        factor = arrhenius(seas['lock'], seas[name])
        added = [w * (factor - 1.0) for w in wild]
        rem = [r + a for r, a in zip(removal, added)]
        out[name] = dict(source_factor_from_the_lock=factor, wild_added_tg_yr=added, removal_tg_yr=rem,
                         air_share_per_year=[r * s for r, s in zip(rem, share_per_tg)],
                         soils_alone_steady_ppm=[p * (x + a) / x for p, x, a in zip(ppm, left, added)])
    out['basis'] = ('an indicator: the wild lakes\' and wetlands\' methane the biology leaves at its 90% cut, scaled '
                    'by the seasonal temperature dependence of methane emissions from the lock\'s sea median to the '
                    'resort\'s at 301 K and the hot world\'s; the farms\' and herds\' methane and the soils\' uptake '
                    'held; substrate and oxygen also set the lakes\' methane')
    return out


def co2(data, b):
    """(b) CO2 below 400 ppm: the column's forcing against the C3 leaf model and water-use efficiency."""
    c = data[COLUMN]['moon_1.2atm']['co2_at_290k']
    ppm, olr, asr = c['co2_ppm'], c['olr'], c['asr']
    i280 = ppm.index(280.0)
    olr400 = olr[ppm.index(400.0)]
    # the unshielded 400-ppm row is the one whose sunlight sits with the other CO2 rows
    asr400 = min((a for p, a in zip(ppm, asr) if p == 400.0), key=lambda a: abs(a - asr[i280]))
    clear_280 = olr[i280] - olr400
    sunlight_280 = asr400 - asr[i280]        # positive: less CO2 absorbs less sunlight, which also cools
    per_doubling = -c['olr_change_per_doubling']
    ratio = co2_all_sky_ratio(data)
    photo = data[CO2]['photosynthesis']
    pa_per_ppm = 1.0 / data[CO2]['design_air']['ppm_per_pa']

    def rate(table, pa):
        xs = sorted(float(k) for k in table)
        return float(np.interp(pa, xs, [table[f'{x:g}'] for x in xs]))
    levels = {}
    for level in HOLD_CO2_PPM:
        doublings = math.log2(400.0 / level)
        infrared = clear_280 if level == 280.0 else per_doubling * doublings
        sunlight = sunlight_280 * doublings / math.log2(400.0 / 280.0)
        clear = infrared + sunlight
        saved = pair(ratio * clear, clear)
        leaf = {t: rate(v['relative_to_earth_today'], level * pa_per_ppm)
                / rate(v['relative_to_earth_today'], 400.0 * pa_per_ppm) for t, v in photo.items()}
        levels[f'{level:g} ppm'] = dict(
            partial_pressure_pa=level * pa_per_ppm, clear_sky_infrared_w_m2=infrared, clear_sky_sunlight_w_m2=sunlight,
            clear_sky_w_m2=clear, saved_w_m2=saved, points=b.points(saved),
            leaf_rate_over_400_ppm=leaf, water_use_efficiency_over_400_ppm=level / 400.0,
            growth_lost_per_point_percent=pair(100.0 * (1.0 - level / 400.0) / b.points(saved[1]),
                                               100.0 * (1.0 - min(leaf.values())) / b.points(saved[0])))
    return dict(column_ground_k=290.0, per_doubling_clear_w_m2=per_doubling, all_sky_ratio=ratio,
                sunlight_change_280_w_m2=sunlight_280, levels=levels,
                basis='clear-sky change in the 1.2 atm column at 290 K, the only ground temperature its CO2 sweep '
                      'covers: the infrared it lets out and the sunlight it no longer absorbs, both taking heat off, '
                      'the sunlight carried to 350 ppm by the log of the CO2 ratio; all-sky by the AR6 effective '
                      'forcing of doubled CO2 over the Earth control\'s own per-doubling change; growth by the CO2 '
                      'study\'s light-saturated C3 leaf model at 25-35 C and by water-use efficiency at a fixed '
                      'internal ratio')


def pressure(data, b):
    """(c) 1.0 against 1.2 atm: the clear-sky columns' net at fixed ground temperature, with the plants' CO2 kept,
    under a clear sky and with Earth's cloud masking."""
    s = data[COLUMN]

    def net(case):
        e = s[case]['earthlike_humidity']
        t = s[case]['shields']['titania_stack']['manabe_wetherald']
        return e['ts_k'], t['asr'], e['olr']
    ts, asr10, olr10 = net('moon_1.0atm')
    _, asr12, olr12 = net('moon_1.2atm')
    sunlight = [x - y for x, y in zip(asr10, asr12)]
    infrared = [x - y for x, y in zip(olr10, olr12)]
    # the 1.0 atm column holds 400 ppm, 17% less CO2 partial pressure; the plants' 48.6 Pa needs 480 ppm
    per_doubling = -s['moon_1.0atm']['co2_at_290k']['olr_change_per_doubling']
    co2_back = per_doubling * math.log2(121590.0 / 101325.0)
    ratio = co2_all_sky_ratio(data)
    t = data[TARGETS]
    hot = t['runs']['test']['ground_k']
    resort_cover = CLOUD_COVER['lock'] + b.resort['share_of_the_way'] * (CLOUD_COVER['hot'] - CLOUD_COVER['lock'])
    out = {}
    for name, temp, cover in (('lock', 295.0, CLOUD_COVER['lock']), ('resort', RESORT_K, resort_cover),
                              ('hot world', hot, CLOUD_COVER['hot']),
                              ('hot world settled', t['settled']['balance_k'][1], CLOUD_COVER['hot'])):
        a, o = float(np.interp(temp, ts, sunlight)), float(np.interp(temp, ts, infrared))
        d = a - o
        lever = -(d + co2_back)          # positive takes heat off the balance
        masked = pair(*(-(f * d) - ratio * co2_back for f in (masking(EARTH_CLOUD_COVER), masking(cover))))
        out[name] = dict(ground_k=temp, sunlight_change_w_m2=a, infrared_change_w_m2=o, column_net_change_w_m2=d,
                         with_the_plants_co2_kept_w_m2=lever, points=b.points(lever), gcm_cover=cover,
                         with_earths_cloud_masking_w_m2=masked, with_earths_cloud_masking_points=b.points(masked))
    scr = data[SCREEN]['inventories']
    o2_pa, ar_ratio, co2_pa = 21227.0, 0.0093 / 0.7808, 48.64

    def n2_pa(total):
        return (total - o2_pa - co2_pa) / (1.0 + ar_ratio)
    share = 1.0 - n2_pa(101325.0) / n2_pa(121590.0)
    return dict(by_target=out, co2_kept_w_m2=co2_back, n2_saved_share=share, n2_saved_kg=share * scr['n2_kg'],
                cloud_masking=dict(earths_cover=masking(EARTH_CLOUD_COVER), co2_all_sky_ratio=ratio),
                basis='titania-stack sunlight and Manabe-Wetherald humidity, clear sky, the 1.0 atm column less the '
                      '1.2 atm column at the same ground temperature; CO2 raised to 480 ppm at 1.0 atm to keep 48.6 '
                      'Pa. The clear sky is the upper end; with clouds, the column\'s net change takes the masking '
                      'Earth\'s clouds give the surface\'s reflection, at Earth\'s cover and at the GCM\'s (the '
                      'resort\'s interpolated by its share of the way), and the CO2 kept its all-sky ratio')


def cirrus(data, b):
    """(d) Thinning the refit's high cloud by seeding: its size in the GCM and Earth's estimates carried over."""
    t = data[TARGETS]
    cloud = t['settled']['cloud_effect_w_m2']
    hot_net = cloud['sunlight'] + cloud['infrared']
    d = t['runs']['design']['top_of_air_w_m2']
    design_net = (d['absorbed'] - d['absorbed_clear_sky']) + (d['outgoing_infrared_clear_sky'] - d['outgoing_infrared'])
    high = hot_net - design_net
    share = [x / CIRRUS_EARTH_NET_W_M2 for x in CCT_EARTH_MAX_W_M2]
    best = [0.0, share[1] * high]
    r = K.MOON_RADIUS
    lo, hi = (r + 1e3 * x for x in CIRRUS_LAYER_KM)
    volume = 4.0 / 3.0 * math.pi * (hi ** 3 - lo ** 3)
    particle = math.pi / 6.0 * INP_DIAMETER_M ** 3 * INP_DENSITY_KG_M3
    aloft = [n * 1e3 * particle * volume for n in INP_PER_LITRE]
    per_year = [m * K.JULIAN_YEAR_DAYS / INP_LIFE_DAYS for m in aloft]
    return dict(hot_world_cloud_net_w_m2=hot_net, design_case_cloud_net_w_m2=design_net,
                refit_high_cloud_net_w_m2=high, earth_share_thinning_could_remove=share,
                thinning_w_m2=best, points=b.points(best), overseeding_on_earth_w_m2=CCT_OVERSEEDING_W_M2,
                seeding=dict(layer_km=list(CIRRUS_LAYER_KM), particles_per_litre=list(INP_PER_LITRE),
                             aloft_kg=aloft, per_year_kg=per_year, basis='estimate: 0.5 um feldspar particles '
                             'through the layer, renewed every 10 days'))


def beams(data, b):
    """(e, i) Beamed power used in orbit in place of on the Moon: the food it makes, and its carbon."""
    t = data[TARGETS]['heat_beyond_sunlight']['beamed_power']
    lw = data[LEDGERS]
    arr = lw['heat']['array']
    link = data[POWER]['absorption']['link']['array_heat_with_close_computing_tw']
    scale = [w / s for w, s in zip(link['with_air'], link['synthesis'])]
    food = [f * k for f, k in zip(arr['heat_by_component_tw']['food'], scale)]
    food_w = [x * b.moon_per_tw for x in food]
    # the food the moved power makes: the blend's wholly powered mixes, and power-heavy's lamp-lit vegetables and
    # fruit and its microbial animal part, as shares of the diet's energy
    weights = lw['choices']['blend_weights']
    diet = data[FOOD]['basis']['diet']
    heavy = [diet['vegfruit_energy_share'] + diet['animal_source_protein_g_day'] * k / diet['intake_kcal_day']
             for k in diet['protein_kcal_per_g']]
    moved = [sum(weights[m] for m in WHOLE_MODES) + weights['power_heavy'] * h for h in heavy]
    people = lw['choices']['moon_billion'] * 1e9
    eaten = people * DIET_KG_C_PER_DAY * K.JULIAN_YEAR_DAYS / 1e12          # Gt C a year, the whole diet eaten
    carbon = [eaten * moved[0], eaten * moved[1] * diet['loss_factor'][1]]
    imports = [v['import_gt_c_yr']['10000'] for k, v in data[CARBON]['questions']['imports'].items()
               if len(k.split('/')) == 3 and k.split('/')[1] == 'held']
    managed = [min(imports), max(imports)]
    return dict(all_beams_w_m2=t['absorbed_sunlight_w_m2'], all_beams_points=t['dimming_percent'],
                all_beams_tw=t['heat_tw'], by_component_tw=arr['heat_by_component_tw'],
                food_tw=food, food_w_m2=food_w, food_points=b.points(food_w),
                food_in_orbit_w_m2=[x * b.ring_per_tw for x in food],
                food_carbon=dict(blend_weights=weights, power_heavy_powered_share=heavy, moved_diet_share=moved,
                                 people=people, whole_diet_eaten_gt_c_yr=eaten,
                                 loss_factor=diet['loss_factor'], arriving_gt_c_yr=carbon,
                                 managed_import_gt_c_yr_at_10000=managed,
                                 times_the_managed_import=[carbon[0] / managed[1], carbon[1] / managed[0]]),
                basis='the domain synthesis\'s array case by component, scaled to the power screen\'s heat with the '
                      'air\'s share; food made at the ring radius puts 5e-5 W/m2 a TW on the Moon. Its carbon: the '
                      'blend\'s power-made and fully lamp-lit mixes and power-heavy\'s lamp-lit vegetables and fruit '
                      'and microbial protein, as shares of the diet\'s energy, at 0.25 kg C a person a day (estimate), '
                      'from the food eaten to the food made at the food study\'s high loss factor, against the carbon '
                      'cycle\'s managed import at year 10,000, methane held, both routes and every corner')


def top_up(data, b):
    """(e, ii) The water top-up's fall: what braking above the exobase could take across the exobase's range, and the
    momentum of the inflow arriving there."""
    w = data[TARGETS]['heat_beyond_sunlight']['water_top_up']
    year = K.JULIAN_YEAR_DAYS * K.JULIAN_DAY
    flow = [x * 1e12 / year for x in w['first_millennium_central_km3_yr']]
    fall = data[SCREEN]['delivery']['fall_mj_per_kg']
    exo = [x * 1e6 for x in fall['exobase']]                 # J/kg from the exobase down, its low and high radius
    far = fall['far_terminal'] * 1e6
    if not math.isclose(w['fall_mj_per_kg'][0], fall['exobase'][0]) or \
            not math.isclose(w['fall_mj_per_kg'][1], fall['far_terminal']):
        raise ValueError('the targets\' fall is not the first screen\'s')
    floor = dict(low_flow=flow[0] * exo[0] / b.area, high_flow=[flow[1] * e / b.area for e in exo])
    above = sorted(w['first_millennium_w_m2'][1] - x for x in floor['high_flow'])
    gm_r = K.MOON_GM / K.MOON_RADIUS
    radii = [K.MOON_GM / (gm_r - e) / K.MOON_RADIUS for e in exo]
    arrival = sorted(math.sqrt(2.0 * (far - e)) for e in exo)    # the speed after the fall from the far terminal
    momentum = [flow[0] * arrival[0], flow[1] * arrival[1]]
    at_the_low_exobase = [f * arrival[1] for f in flow]
    photon = 2.0 * K.SOLAR_CONSTANT / K.SPEED_OF_LIGHT
    disc = math.pi * K.MOON_RADIUS ** 2
    elevator_kg_s = ELEVATOR_T_PER_YEAR * 1e3 / year
    shed = [flow[0] * min(far - e for e in exo), flow[1] * max(far - e for e in exo)]
    return dict(first_millennium_w_m2=w['first_millennium_w_m2'],
                first_millennium_points=b.points(w['first_millennium_w_m2']),
                flow_kg_s=flow, fall_power_tw=[flow[0] * exo[0] / 1e12, flow[1] * far / 1e12],
                exobase_fall_mj_per_kg=fall['exobase'], far_terminal_fall_mj_per_kg=fall['far_terminal'],
                exobase_radius_moon_radii=radii,
                first_screen_exobase_radius=data[SCREEN]['delivery']['exobase_radius_r'],
                exobase_floor_w_m2=floor, above_the_exobase_w_m2=[0.0, above[1]],
                above_the_exobase_at_the_high_end_w_m2=above, above_points=[0.0, b.points(above[1])],
                above_the_exobase_at_the_high_end_points=b.points(above),
                orbital_energy_to_shed_tw=[x / 1e12 for x in shed],
                arrival_speed_m_s=arrival, momentum_to_shed_n=momentum,
                momentum_at_the_low_exobase_n=at_the_low_exobase,
                reflecting_sail_lunar_discs=[m / photon / disc for m in momentum],
                elevator_shortfall=[f / elevator_kg_s for f in flow],
                basis='the top-up\'s first millennium (climate_targets) from the exobase down at the low end and '
                      'from the far cislunar terminal at the high; the exobase at the first screen\'s 1.77-2.56 lunar '
                      'radii; braking to rest at the exobase removes the speed of the fall from the far terminal')


def glow(data, b):
    """(e, iii) The glow at both of the GCM's efficacies: what the roll took, what it leaves, and the further answers
    the author set aside."""
    g = data[GLOW]
    efficacy = data[TARGETS]['heat_beyond_sunlight']['glow']['efficacy']
    albedo = g['design']['planetary_albedo']
    ans = g['results']['answers']

    def equiv(key, e):
        a = ans[key]
        ir = a['glow_W_m2']['global'] + a['infrared_return_W_m2']['global']
        return e * ir + (1.0 - albedo) * a['sunlight_return_W_m2']['global']
    rec = 'rolled_outside_service_80_rolled_layout'
    out = dict(efficacy=efficacy, efficacy_sources=data[TARGETS]['heat_beyond_sunlight']['glow']['efficacy_sources'],
               as_built_w_m2=[equiv('as_built', e) for e in efficacy],
               recommended_w_m2=[equiv(rec, e) for e in efficacy])
    out['the_roll_took_w_m2'] = [a - r for a, r in zip(out['as_built_w_m2'], out['recommended_w_m2'])]
    out['the_roll_took_points'] = b.points(out['the_roll_took_w_m2'])
    out['recommended_points'] = b.points(out['recommended_w_m2'])
    out['the_lower_efficacy_takes_points'] = out['recommended_points'][1] - out['recommended_points'][0]
    further = {}
    for name, key in (
            ('annulus absorber, 250 nm edge', 'rolled_outside_service_80_wide_gap_250_rolled_layout'),
            ('ultraviolet mirror, outer faces', 'rolled_outside_service_80_uv_mirror_outer_98_rolled_layout')):
        w = [equiv(key, e) for e in efficacy]
        saved = [r - x for r, x in zip(out['recommended_w_m2'], w)]
        further[name] = dict(w_m2=w, saved_w_m2=saved, points=b.points(saved))
    out['set_aside_on_13_october'] = further
    out['light_tw_per_point_of_design_sunlight'] = g['results']['dimmer_light_TW_per_percent']
    return out


def kernels(rc):
    """The clear-sky kernel, the change in reflected sunlight per unit of surface albedo, at each land column of the
    GCM radiation check below HIGHLAND_M: the column's reflected share less the sea's at the same latitude (the sea
    columns interpolated in latitude), over the difference of their surface albedos. Returns (|lat|, kernel, label)."""
    seas = sorted((abs(c['lat_deg']), c['line_by_line']['reflected'], c['surface_albedo']) for c in rc if not c['land'])
    if len({s[2] for s in seas}) != 1:
        raise ValueError('the sea columns do not share one albedo')
    out = []
    for c in rc:
        if c['land'] and c['height_m'] < HIGHLAND_M:
            lat = abs(c['lat_deg'])
            sea = float(np.interp(lat, [s[0] for s in seas], [s[1] for s in seas]))
            out.append((lat, (c['line_by_line']['reflected'] - sea) / (c['surface_albedo'] - seas[0][2]), c['label']))
    return sorted(out)


def kernel_by_latitude(found, lat):
    """[low, high] kernel at |lat|: linear between the columns' latitudes, the spread where two columns share one;
    poleward of the last the high end held there and the low end carried on along its last slope."""
    lats = sorted({round(x[0], 4) for x in found})
    lo = [min(k for x, k, _ in found if round(x, 4) == L) for L in lats]
    hi = [max(k for x, k, _ in found if round(x, 4) == L) for L in lats]
    a = np.abs(np.asarray(lat, float))
    low = np.interp(a, lats, lo)
    slope = (lo[-1] - lo[-2]) / (lats[-1] - lats[-2])
    low = np.where(a > lats[-1], lo[-1] + slope * (a - lats[-1]), low)
    return low, np.interp(a, lats, hi)


def albedo(data, b):
    """(f) The Moon's own reflection: roofs and pavements, silver foliage in the parks, and the deserts."""
    rc = data[RADIATION]['columns']
    found = kernels(rc)
    land = next(c for c in rc if c['land'] and 'site' in c['label'])
    sea = min((c for c in rc if not c['land']), key=lambda c: abs(c['lat_deg'] + land['lat_deg']))
    kernel = ((land['line_by_line']['reflected'] - sea['line_by_line']['reflected'])
              / (land['surface_albedo'] - sea['surface_albedo']))
    factor = {k: masking(c) for k, c in CLOUD_COVER.items()}
    k_all = pair(kernel * factor['hot'], kernel * factor['lock'])
    land_l = data[LEDGERS]['land']
    settled = [x * 1e6 / b.area for x in land_l['settled_km2']]   # share of the Moon
    parks = land_l['parks_km2'] * 1e6 / b.area
    s = b.incoming

    def forcing(share, weight, rise, k=(k_all[0], k_all[1])):
        return pair(k[0] * s * weight[0] * rise[0] * share[0], k[1] * s * weight[1] * rise[1] * share[1])
    cities = forcing(settled, CITY_SUN_WEIGHT, CITY_ALBEDO_RISE)
    park = forcing([parks, parks], PARK_SUN_WEIGHT, PARK_ALBEDO_RISE)
    # the deserts by band: the model's own land mask from the GCM, sunlight by the diurnal mean of cos(lat), the
    # kernel at each row's latitude, the deserts as the extract counts them on every cell's land and on the model's
    # own land cells
    clim = data[CLIMATOLOGY]
    lat, lsm = clim['lat'], clim['lsm']
    w = np.cos(np.radians(lat))[:, None] * np.ones_like(lsm)
    sun = (4.0 / math.pi) * np.cos(np.radians(lat))[:, None] * np.ones_like(lsm)
    k_lo, k_hi = kernel_by_latitude(found, lat)
    k_lo, k_hi = k_lo[:, None] * np.ones_like(lsm), k_hi[:, None] * np.ones_like(lsm)
    ex = data[EXTRACT]['runs']
    deserts = {}
    for name, run in (('hot world', 'A28_corrected_970_lw_ice'), ('lock', 'A28_corrected_970')):
        area, sunlit, kernelled = {}, {}, {}
        for count, v in ex[run]['land'].items():
            area[count] = sunlit[count] = 0.0
            kernelled[count] = [0.0, 0.0]
            for band, d in v['bands'].items():
                lo, hi = (float(x) for x in band.split('-'))
                m = (np.abs(lat)[:, None] >= lo) & (np.abs(lat)[:, None] < hi) & (lsm > 0.5)
                area[count] += float((w * m).sum() / w.sum()) * d['desert']
                sunlit[count] += float((w * sun * m).sum() / w.sum()) * d['desert']
                for i, k in enumerate((k_lo, k_hi)):
                    kernelled[count][i] += float((w * sun * k * m).sum() / w.sum()) * d['desert']
        ends = [min(x[0] for x in kernelled.values()), max(x[1] for x in kernelled.values())]
        flat_ends = [kernel * min(sunlit.values()), kernel * max(sunlit.values())]

        def brighter(e):
            return pair(factor['hot'] * s * DESERT_ALBEDO_RISE[0] * e[0],
                        factor['lock'] * s * DESERT_ALBEDO_RISE[1] * e[1])
        deserts[name] = dict(desert_share_of_moon=area, sunlit_desert_share=sunlit,
                             kernel_over_the_equators={c: [x / (kernel * sunlit[c]) for x in k]
                                                       for c, k in kernelled.items()},
                             w_m2=brighter(ends), points=b.points(brighter(ends)),
                             with_the_equators_kernel_w_m2=brighter(flat_ends),
                             with_the_equators_kernel_points=b.points(brighter(flat_ends)))
    # the cloud over the lock's own deserts, the design case's climatology, by area and sunlight
    dry = (clim['pr_mm_day'] < DESERT_MM_DAY) & (lsm > 0.5)
    desert_cover = float((w * sun * clim['clt_mean'] * dry).sum() / (w * sun * dry).sum())
    return dict(clear_sky_kernel=kernel, kernels_by_column=[dict(lat_deg=x, kernel=k, column=c) for x, k, c in found],
                cloud_factor=factor, all_sky_kernel=k_all, kernel_columns=[land['label'], sea['label']],
                lock_deserts_cloud_cover=desert_cover, cloud_factor_at_the_lock_deserts_cover=masking(desert_cover),
                settled_share_of_moon=settled, cities_w_m2=cities, cities_points=b.points(cities),
                parks_share_of_moon=parks, parks_w_m2=park, parks_points=b.points(park), deserts=deserts,
                basis='the change in reflected sunlight per unit surface albedo between the land and sea columns of '
                      'the GCM radiation check (line-by-line, clear sky), the equatorial pair for the cities and parks '
                      'and each desert row\'s latitude for the deserts; Earth\'s cloud masking from Lenton & Vaughan '
                      'scaled by the GCM\'s cover; the deserts by band as the extract counts them on every cell\'s '
                      'land and on the model\'s own land cells, on the model\'s land mask, the ends the least and the '
                      'most of the two counts and the kernel\'s two ends; the design guesses in run.py')


def shield(data, b):
    """The shield's own cost per point, for comparison (read only)."""
    t = data[TARGETS]
    g = data[GLOW]['results']
    return dict(light_tw_per_point=g['dimmer_light_TW_per_percent'],
                absorbing_dimmer_glow_return=dict(rolled=0.10, night_roll=0.14, as_built=0.25),
                reflecting_dimmer_window_coefficient=dict(at=0.276, after_two_points=0.32),
                growth_through_light_percent_per_point=t['costs_per_percent']['growth_through_light_percent'],
                land_rain_percent_per_point=t['costs_per_percent']['land_rain_percent'],
                basis='glow.md: a dimmer absorbing in the window films returns 25% of its dimming as glow as built, '
                      '14% with the night roll, 10% with both halves rolled; one reflecting from the window tiles adds '
                      '2dR to their light pressure, 0.276 to about 0.32 for two points; climate_targets\' costs')


def ranking(b, levers):
    """The README's ranked table in points, and what the first four save together against each target."""
    g, a, p = levers['glow']['set_aside_on_13_october'], levers['albedo'], levers['pressure']['by_target']
    answers = [x for v in g.values() for x in v['points']]
    roofs = [c + k for c, k in zip(a['cities_points'], a['parks_points'])]
    t = levers['top_up']
    rows = [
        ('the glow\'s set-aside answers', [min(answers), max(answers)]),
        ('food made in orbit', levers['beams']['food_points']),
        ('methane and N2O at the low pair', levers['gases']['options']['low pair, 0.3 and 0.03 ppm']['points']),
        ('white roofs and silver foliage', roofs),
        ('night-side cirrus thinning', levers['cirrus']['points']),
        ('the top-up\'s fall', t['above_points']),
        ('1.0 atm in the hot world', [min(p['hot world']['with_earths_cloud_masking_points']),
                                      p['hot world settled']['points']]),
        ('brighter deserts', [a['deserts']['hot world']['points'][0], a['deserts']['lock']['points'][1]]),
        ('CO2 below 400 ppm', [levers['co2']['levels']['350 ppm']['points'][0],
                               levers['co2']['levels']['280 ppm']['points'][1]])]
    four = [sum(r[1][0] for r in rows[:4]), sum(r[1][1] for r in rows[:4])]
    return dict(levers=[dict(rank=i + 1, lever=n, points=v) for i, (n, v) in enumerate(rows)],
                first_four_points=four,
                first_four_share_of_dimming_in_all={k: [four[0] / d[1], four[1] / d[0]] for k, d in b.dimming.items()},
                details=dict(top_up_high_end_across_the_exobase_points=t['above_the_exobase_at_the_high_end_points'],
                             pressure_hot_world_clear_sky_points=[p['hot world']['points'],
                                                                  p['hot world settled']['points']],
                             pressure_hot_world_masked_points=[
                                 p['hot world']['with_earths_cloud_masking_points'][0],
                                 p['hot world settled']['with_earths_cloud_masking_points'][1]],
                             pressure_lock_points=[p['lock']['points'],
                                                   p['lock']['with_earths_cloud_masking_points'][0]]))


def results():
    data = load()
    b = Basis(data)
    levers = dict(gases=gases(data, b), co2=co2(data, b), pressure=pressure(data, b), cirrus=cirrus(data, b),
                  beams=beams(data, b), top_up=top_up(data, b), glow=glow(data, b), albedo=albedo(data, b),
                  shield=shield(data, b))
    return dict(
        schema=SCHEMA,
        producer=dict(files={c: digest(c) for c in CODE}, inputs={p: digest(p) for p in INPUTS},
                      constants=dict(MOON_RADIUS=K.MOON_RADIUS, MOON_GM=K.MOON_GM, JULIAN_YEAR_DAYS=K.JULIAN_YEAR_DAYS,
                                     JULIAN_DAY=K.JULIAN_DAY, BOLTZMANN=K.BOLTZMANN,
                                     ELEMENTARY_CHARGE=K.ELEMENTARY_CHARGE, SOLAR_CONSTANT=K.SOLAR_CONSTANT,
                                     SPEED_OF_LIGHT=K.SPEED_OF_LIGHT)),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        basis=dict(absorbed_sunlight_w_m2=b.absorbed, incoming_w_m2=b.incoming, moon_w_m2_per_tw=b.moon_per_tw,
                   ring_w_m2_per_tw=b.ring_per_tw, resort_k=RESORT_K,
                   resort_share_of_the_way=b.resort['share_of_the_way'], targets_dimming_in_all_percent=b.dimming),
        **levers, ranking=ranking(b, levers))


def main(argv=None) -> int:
    product = results()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    print('wrote', OUT.relative_to(ROOT))
    return 0


if __name__ == '__main__':
    sys.exit(main())
