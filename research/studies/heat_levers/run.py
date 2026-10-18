"""What could take pressure off the shield's dimming: each lever on the greenhouse effect, the plan's added heat or
the Moon's own reflection, in W/m2 of the Moon's balance and in points of dimming against the hot world's absorbed
sunlight, with its cost and side effects from the project's products.

    python -m research.studies.heat_levers.run        # under a second

A screening on committed products: the target climates, the deep column's gas forcing and CO2 and pressure sweeps,
the domain synthesis's gas and heat ledgers, the power screen, the fleet's glow (read only), the CO2 study's leaf
model, the GCM's clear-sky radiation check and climatology, and the climate targets' extract for the deserts by band.
No model is run.
"""
from __future__ import annotations
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

from shared import constants as K

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'heat_levers.json'
SCHEMA = 'terluna.research.heat-levers/1'
CODE = ('research/studies/heat_levers/run.py',)
TARGETS = 'research/studies/climate_targets/results/climate_targets.json'
GASES = 'atmosphere/radiative_convective/results/reduced_gases.json'
COLUMN = 'atmosphere/radiative_convective/results/summary.json'
LEDGERS = 'research/studies/domain_synthesis/results/ledgers.json'
POWER = 'provisioning/results/power.json'
GLOW = 'research/studies/solar_shield_array/results/glow.json'
CO2 = 'research/studies/atmospheric_co2/results/atmospheric_co2.json'
RADIATION = 'climate/results/gcm/radiation_check_A28_dim5_moon.json'
CLIMATOLOGY = 'climate/gcm/products/climatology_A28_corrected_970.npz'
INDUSTRY = 'resources/results/industry.json'
SCREEN = 'resources/results/first_screen.json'
EXTRACT = 'research/runs/climate_targets/extract.json'     # ignored, on the research drive; deserts by band
INPUTS = {TARGETS: 'terluna.research.climate-targets/1', GASES: 'terluna.atmosphere.reduced-gas-forcing/1',
          COLUMN: None, LEDGERS: 'terluna.research.domain-synthesis-ledgers/1', POWER: 'terluna.provisioning.power/1',
          GLOW: 'terluna.research.fleet-glow/1', CO2: 'terluna.research.atmospheric-co2/1',
          RADIATION: 'terluna.climate.gcm-radiation-check/1', CLIMATOLOGY: 'terluna.climate.gcm-climatology/1',
          INDUSTRY: 'terluna.resources.industry/1', SCREEN: 'terluna.resources.first-screen/1',
          EXTRACT: 'terluna.research.climate-targets-extract/1'}

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
INP_PER_LITRE = (1.0, 20.0)        # seeding particles in the cirrus layer (Gasparini et al. 2020: 1 and 18)
INP_DIAMETER_M = 0.5e-6
INP_DENSITY_KG_M3 = 2730.0         # plagioclase feldspar, the highlands' own mineral
INP_LIFE_DAYS = 10.0
CIRRUS_LAYER_KM = (63.0, 107.0)    # the GCM's two cloudy top layers after the longwave refit
DIET_KG_C_PER_DAY = 0.25           # carbon a person eats in a day, about 600 g of dry food
POWERED_DIET_SHARE = (0.2, 0.6)    # share of the diet made with beamed power that could be made in orbit
PEOPLE_ON_THE_MOON = 8e9           # C'
BANDS = ((0, 15), (15, 30), (30, 45), (45, 60), (60, 90))

EVIDENCE = ('A screening on committed products. Forcings come from the deep column (line-by-line, clear sky, '
            'prescribed humidity; all-sky by AR6 ratios on its Earth control), the cloud effects and absorbed '
            'sunlight from two ExoPlaSim runs, costs from the domain synthesis, the power screen and the glow study '
            '(read only), surface-albedo kernels from the GCM radiation check\'s clear-sky line-by-line shares with '
            'Earth\'s cloud masking scaled by cover, and the cirrus, methane and albedo responses from published Earth '
            'studies. No lever is run in a climate model; each W/m2 is treated as absorbed sunlight at efficacy 1 '
            'except the glow, which takes the GCM\'s 0.71.')
READING_RULE = ('W/m2 are means over the whole Moon, positive for what a lever takes off the balance. A point of '
                'dimming is 1% of the hot world\'s absorbed sunlight, 225.4 W/m2, as climate_targets counts the heat '
                'beyond the sunlight. Pairs are [low, high]. TW of heat on the Moon count 0.0264 W/m2 each, at the '
                'ring radius 5e-5. Values marked estimate rest on the design guesses named in run.py.')


def digest(path) -> str:
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def load():
    data = {}
    for p, schema in INPUTS.items():
        if p.endswith('.npz'):
            z = np.load(ROOT / p, allow_pickle=False)
            meta = json.loads(str(z['metadata']))
            data[p] = dict(meta=meta, lat=np.asarray(z['lat'], float), lsm=np.asarray(z['lsm'], float))
            got = meta.get('schema')
        else:
            data[p] = json.loads((ROOT / p).read_text())
            got = data[p].get('schema') if isinstance(data[p], dict) else None
        if schema is not None and got != schema:
            raise ValueError(f'{p}: schema {got} is not {schema}')
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

    def points(self, w):
        if isinstance(w, (list, tuple)):
            return [100.0 * x / self.absorbed for x in w]
        return 100.0 * w / self.absorbed


def gases(data, b):
    """(a) Methane and N2O below the held pair: the forcing saved, what holding them takes, and the warm Moon's lakes."""
    g = data[GASES]
    m = g['moon_at_design']['both']
    held_key = '1.8663+0.1'
    held = m[held_key]['allsky_w_m2']
    options = {'low pair, 0.3 and 0.03 ppm': '0.3+0.03', 'methane low, 0.3 and 0.1 ppm': '0.3+0.1',
               'N2O low, 1.87 and 0.03 ppm': '1.8663+0.03', 'lowest grid pair, 0.1 and 0.01 ppm': '0.1+0.01'}
    # the held pair's clear-sky forcing falls as the columns warm; carried to the hot world's ground linearly
    col = g['columns']['moon_12atm']
    temps = sorted(float(k) for k in col)
    toa = [col[f'{t:g}']['both'][held_key]['toa'] for t in temps]
    slope = float(np.polyfit(temps, toa, 1)[0])
    hot_k = data[TARGETS]['runs']['test']['ground_k']
    warm_factor = (toa[temps.index(295.0)] + slope * (hot_k - 295.0)) / toa[temps.index(295.0)]
    lw = data[LEDGERS]
    w_per_tw = lw['heat']['factors']['moon_w_m2_per_tw']
    out = dict(held_pair_allsky_w_m2=held, held_pair_points=b.points(held),
               forcing_at_the_hot_worlds_ground_factor=warm_factor, column_temperatures_k=temps,
               held_pair_clear_sky_toa_w_m2=toa)
    opts = {}
    for name, key in options.items():
        saved = held - m[key]['allsky_w_m2']
        opts[name] = dict(allsky_w_m2=m[key]['allsky_w_m2'], saved_w_m2=pair(saved * warm_factor, saved),
                          points=b.points(pair(saved * warm_factor, saved)))
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
    # the warm Moon's waters: methanogenesis's temperature dependence between the targets' sea medians (indicator)
    kb_ev = K.BOLTZMANN / K.ELEMENTARY_CHARGE
    t = data[TARGETS]
    seas = {'lock': t['runs']['design']['sea_median_c'], 'hot': t['runs']['test']['sea_median_c'],
            'resort': next(x for x in t['targets'] if x['name'] == '300 K')['climate']['sea_median_c']}

    def arrhenius(c_from, c_to):
        return math.exp(METHANOGENESIS_EV / kb_ev * (1.0 / (c_from + 273.15) - 1.0 / (c_to + 273.15)))
    held_ch4 = lw['gases']['held']['methane_by_source_cut']['0.9']
    factor = arrhenius(seas['lock'], seas['hot'])
    left, removal = held_ch4['sources_left_tg_yr'], held_ch4['removal_tg_yr']
    share_per_tg = [s / r for s, r in zip(held_ch4['atmosphere_share_per_year'], removal)]
    hot_removal = [l * (factor - 1.0) + r for l, r in zip(left, removal)]
    out['warm_waters'] = dict(
        sea_median_c=seas, activation_ev=METHANOGENESIS_EV,
        source_factor_lock_to_hot=factor, source_factor_lock_to_resort=arrhenius(seas['lock'], seas['resort']),
        held_methane_removal_tg_yr=dict(lock=removal, hot=hot_removal),
        held_methane_air_share_per_year=dict(lock=held_ch4['atmosphere_share_per_year'],
                                             hot=[r * s for r, s in zip(hot_removal, share_per_tg)]),
        soils_alone_steady_ppm=dict(lock=held_ch4['steady_ppm_by_soils_alone'],
                                    hot=[x * factor for x in held_ch4['steady_ppm_by_soils_alone']]),
        basis='every wild source scaled by the seasonal temperature dependence of methane emissions between the lock\'s '
              'and the hot world\'s sea medians, with biology cutting 90% at the sources and the soils\' uptake held; '
              'an upper-end indicator, since substrate and oxygen also set the lakes\' methane')
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
    sunlight_280 = asr400 - asr[i280]
    per_doubling = -c['olr_change_per_doubling']
    earth_per_doubling = -data[COLUMN]['earth_1.0atm']['co2_at_290k']['olr_change_per_doubling']
    ratio = AR6_ERF_2XCO2 / earth_per_doubling
    photo = data[CO2]['photosynthesis']
    pa_per_ppm = 1.0 / data[CO2]['design_air']['ppm_per_pa']

    def rate(table, pa):
        xs = sorted(float(k) for k in table)
        return float(np.interp(pa, xs, [table[f'{x:g}'] for x in xs]))
    levels = {}
    for level in HOLD_CO2_PPM:
        clear = per_doubling * math.log2(400.0 / level) if level != 280.0 else clear_280
        net = clear - (sunlight_280 if level == 280.0 else 0.0)
        saved = pair(ratio * net, clear)
        leaf = {t: rate(v['relative_to_earth_today'], level * pa_per_ppm)
                / rate(v['relative_to_earth_today'], 400.0 * pa_per_ppm) for t, v in photo.items()}
        levels[f'{level:g} ppm'] = dict(
            partial_pressure_pa=level * pa_per_ppm, clear_sky_w_m2=clear, saved_w_m2=saved, points=b.points(saved),
            leaf_rate_over_400_ppm=leaf, water_use_efficiency_over_400_ppm=level / 400.0,
            growth_lost_per_point_percent=pair(100.0 * (1.0 - level / 400.0) / b.points(saved[1]),
                                               100.0 * (1.0 - min(leaf.values())) / b.points(saved[0])))
    return dict(per_doubling_clear_w_m2=per_doubling, all_sky_ratio=ratio, sunlight_change_280_w_m2=sunlight_280,
                levels=levels,
                basis='clear-sky OLR change in the 1.2 atm column at 290 K; all-sky by the AR6 effective forcing of '
                      'doubled CO2 over the Earth control\'s own per-doubling change; growth by the CO2 study\'s '
                      'light-saturated C3 leaf model at 25-35 C and by water-use efficiency at a fixed internal ratio')


def pressure(data, b):
    """(c) 1.0 against 1.2 atm: the clear-sky columns' net at fixed ground temperature, with the plants' CO2 kept."""
    s = data[COLUMN]
    out = {}

    def net(case):
        e = s[case]['earthlike_humidity']
        t = s[case]['shields']['titania_stack']['manabe_wetherald']
        return e['ts_k'], [a - o for a, o in zip(t['asr'], e['olr'])]
    ts, n10 = net('moon_1.0atm')
    _, n12 = net('moon_1.2atm')
    diff = [x - y for x, y in zip(n10, n12)]
    # the 1.0 atm column holds 400 ppm, 17% less CO2 partial pressure; the plants' 48.6 Pa needs 480 ppm
    per_doubling = -s['moon_1.0atm']['co2_at_290k']['olr_change_per_doubling']
    co2_back = per_doubling * math.log2(121590.0 / 101325.0)
    t = data[TARGETS]
    hot = t['runs']['test']['ground_k']
    for name, temp in (('lock', 295.0), ('resort', 300.0), ('hot world', hot),
                       ('hot world settled', t['settled']['balance_k'][1])):
        d = float(np.interp(temp, ts, diff))
        lever = -(d + co2_back)          # positive takes heat off the balance
        out[name] = dict(ground_k=temp, column_net_change_w_m2=d, with_the_plants_co2_kept_w_m2=lever,
                         points=b.points(lever))
    scr = data[SCREEN]['inventories']
    o2_pa, ar_ratio, co2_pa = 21227.0, 0.0093 / 0.7808, 48.64

    def n2_pa(total):
        return (total - o2_pa - co2_pa) / (1.0 + ar_ratio)
    share = 1.0 - n2_pa(101325.0) / n2_pa(121590.0)
    return dict(by_target=out, co2_kept_w_m2=co2_back, n2_saved_share=share, n2_saved_kg=share * scr['n2_kg'],
                basis='titania-stack sunlight and Manabe-Wetherald humidity, clear sky, the 1.0 atm column less the '
                      '1.2 atm column at the same ground temperature; CO2 raised to 480 ppm at 1.0 atm to keep 48.6 Pa')


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
    arr = data[LEDGERS]['heat']['array']
    link = data[POWER]['absorption']['link']['array_heat_with_close_computing_tw']
    scale = [w / s for w, s in zip(link['with_air'], link['synthesis'])]
    food = [f * k for f, k in zip(arr['heat_by_component_tw']['food'], scale)]
    food_w = [x * b.moon_per_tw for x in food]
    year_days = K.JULIAN_YEAR_DAYS
    carbon = [PEOPLE_ON_THE_MOON * DIET_KG_C_PER_DAY * year_days * s / 1e12 for s in POWERED_DIET_SHARE]
    return dict(all_beams_w_m2=t['absorbed_sunlight_w_m2'], all_beams_points=t['dimming_percent'],
                all_beams_tw=t['heat_tw'], by_component_tw=arr['heat_by_component_tw'],
                food_tw=food, food_w_m2=food_w, food_points=b.points(food_w),
                food_in_orbit_w_m2=[x * b.ring_per_tw for x in food],
                carbon_arriving_as_food_gt_c_yr=carbon,
                basis='the domain synthesis\'s array case by component, scaled to the power screen\'s heat with the '
                      'air\'s share; food made at the ring radius puts 5e-5 W/m2 a TW on the Moon; the carbon is an '
                      'estimate on the diet and the powered share named in run.py')


def top_up(data, b):
    """(e, ii) The water top-up's fall: what braking in orbit could take, and why its momentum binds."""
    w = data[TARGETS]['heat_beyond_sunlight']['water_top_up']
    year = K.JULIAN_YEAR_DAYS * K.JULIAN_DAY
    flow = [x * 1e12 / year for x in w['first_millennium_central_km3_yr']]
    low_e, high_e = (x * 1e6 for x in w['fall_mj_per_kg'])
    floor = [f * low_e / b.area for f in flow]
    above = w['first_millennium_w_m2'][1] - floor[1]
    gm_r = K.MOON_GM / K.MOON_RADIUS
    r_exo = K.MOON_GM / (gm_r - low_e)
    v = math.sqrt(K.MOON_GM / r_exo)
    momentum = [f * v for f in flow]
    photon = 2.0 * K.SOLAR_CONSTANT / K.SPEED_OF_LIGHT
    disc = math.pi * K.MOON_RADIUS ** 2
    elevator_kg_s = ELEVATOR_T_PER_YEAR * 1e3 / year
    return dict(first_millennium_w_m2=w['first_millennium_w_m2'], first_millennium_points=b.points(w['first_millennium_w_m2']),
                flow_kg_s=flow, fall_power_tw=[flow[0] * low_e / 1e12, flow[1] * high_e / 1e12],
                exobase_floor_w_m2=floor, above_the_exobase_w_m2=[0.0, above], above_points=b.points([0.0, above]),
                orbital_energy_to_shed_tw=[f * (high_e - low_e) / 1e12 for f in flow],
                exobase_radius_moon_radii=r_exo / K.MOON_RADIUS, circular_speed_m_s=v, momentum_to_shed_n=momentum,
                reflecting_sail_lunar_discs=[m / photon / disc for m in momentum],
                elevator_shortfall=[f / elevator_kg_s for f in flow])


def glow(data, b):
    """(e, iii) The glow: what the roll took, what it leaves, and the further answers the author set aside."""
    g = data[GLOW]
    eff = data[TARGETS]['heat_beyond_sunlight']['glow']['efficacy']
    albedo = g['design']['planetary_albedo']
    ans = g['results']['answers']

    def equiv(key, e=eff):
        a = ans[key]
        ir = a['glow_W_m2']['global'] + a['infrared_return_W_m2']['global']
        return e * ir + (1.0 - albedo) * a['sunlight_return_W_m2']['global']
    rec = 'rolled_outside_service_80_rolled_layout'
    out = dict(efficacy=eff, as_built_w_m2=equiv('as_built'), recommended_w_m2=equiv(rec))
    out['the_roll_took_w_m2'] = out['as_built_w_m2'] - out['recommended_w_m2']
    out['the_roll_took_points'] = b.points(out['the_roll_took_w_m2'])
    out['recommended_points'] = b.points(out['recommended_w_m2'])
    further = {}
    for name, key in (('annulus absorber, 250 nm edge', 'rolled_outside_service_80_wide_gap_250_rolled_layout'),
                      ('ultraviolet mirror, outer faces', 'rolled_outside_service_80_uv_mirror_outer_98_rolled_layout')):
        saved = out['recommended_w_m2'] - equiv(key)
        further[name] = dict(w_m2=equiv(key), saved_w_m2=saved, points=b.points(saved))
    out['set_aside_on_13_october'] = further
    corrected = equiv(rec, 0.59)
    out['at_the_corrected_radiations_efficacy_0_59'] = dict(w_m2=corrected,
                                                            change_points=b.points(corrected - out['recommended_w_m2']))
    out['light_tw_per_point_of_design_sunlight'] = g['results']['dimmer_light_TW_per_percent']
    return out


def albedo(data, b):
    """(f) The Moon's own reflection: roofs and pavements, silver foliage in the parks, and the deserts."""
    rc = data[RADIATION]['columns']
    land = next(c for c in rc if c['land'] and 'site' in c['label'])
    sea = min((c for c in rc if not c['land']), key=lambda c: abs(c['lat_deg'] + land['lat_deg']))
    kernel = ((land['line_by_line']['reflected'] - sea['line_by_line']['reflected'])
              / (land['surface_albedo'] - sea['surface_albedo']))
    cloudy = 1.0 - (1.0 - EARTH_F_A['all_sky'] / EARTH_F_A['clear']) / EARTH_CLOUD_COVER
    factor = {k: 1.0 - c * (1.0 - cloudy) for k, c in CLOUD_COVER.items()}
    k_all = pair(kernel * factor['hot'], kernel * factor['lock'])
    land_l = data[LEDGERS]['land']
    settled = [x * 1e6 / b.area for x in land_l['settled_km2']]   # share of the Moon
    parks = land_l['parks_km2'] * 1e6 / b.area
    s = b.incoming

    def forcing(share, weight, rise):
        return pair(k_all[0] * s * weight[0] * rise[0] * share[0], k_all[1] * s * weight[1] * rise[1] * share[1])
    cities = forcing(settled, CITY_SUN_WEIGHT, CITY_ALBEDO_RISE)
    park = forcing([parks, parks], PARK_SUN_WEIGHT, PARK_ALBEDO_RISE)
    # the deserts by band: land share of the Moon from the GCM's mask, sunlight by the diurnal mean of cos(lat)
    clim = data[CLIMATOLOGY]
    lat, lsm = clim['lat'], clim['lsm']
    w = np.cos(np.radians(lat))[:, None] * np.ones_like(lsm)
    sun = (4.0 / math.pi) * np.cos(np.radians(lat))[:, None] * np.ones_like(lsm)
    ex = data[EXTRACT]['runs']
    deserts = {}
    for name, run in (('hot world', 'A28_corrected_970_lw_ice'), ('lock', 'A28_corrected_970')):
        bands = ex[run]['land']['bands']
        weighted = area = 0.0
        for lo, hi in BANDS:
            m = (np.abs(lat)[:, None] >= lo) & (np.abs(lat)[:, None] < hi) & (lsm > 0.5)
            share = float((w * m).sum() / w.sum())
            frac = bands[f'{lo}-{hi}']['desert']
            area += share * frac
            weighted += float((w * sun * m).sum() / w.sum()) * frac
        rise = forcing([weighted, weighted], (1.0, 1.0), DESERT_ALBEDO_RISE)
        deserts[name] = dict(desert_share_of_moon=area, sunlit_desert_share=weighted, w_m2=rise, points=b.points(rise))
    return dict(clear_sky_kernel=kernel, cloud_factor=factor, all_sky_kernel=k_all,
                kernel_columns=[land['label'], sea['label']],
                settled_share_of_moon=settled, cities_w_m2=cities, cities_points=b.points(cities),
                parks_share_of_moon=parks, parks_w_m2=park, parks_points=b.points(park), deserts=deserts,
                basis='the change in reflected sunlight per unit surface albedo between the equatorial land and sea '
                      'columns of the GCM radiation check (line-by-line, clear sky), Earth\'s cloud masking from '
                      'Lenton & Vaughan scaled by the GCM\'s cover; the design guesses in run.py')


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


def results():
    data = load()
    b = Basis(data)
    t = data[TARGETS]
    targets = {x['name']: x['with_the_heat_beyond_sunlight_percent'] for x in t['targets']}
    return dict(
        schema=SCHEMA,
        producer=dict(files={c: digest(c) for c in CODE}, inputs={p: digest(p) for p in INPUTS},
                      constants=dict(MOON_RADIUS=K.MOON_RADIUS, MOON_GM=K.MOON_GM, JULIAN_YEAR_DAYS=K.JULIAN_YEAR_DAYS,
                                     JULIAN_DAY=K.JULIAN_DAY, BOLTZMANN=K.BOLTZMANN,
                                     ELEMENTARY_CHARGE=K.ELEMENTARY_CHARGE, SOLAR_CONSTANT=K.SOLAR_CONSTANT,
                                     SPEED_OF_LIGHT=K.SPEED_OF_LIGHT)),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        basis=dict(absorbed_sunlight_w_m2=b.absorbed, incoming_w_m2=b.incoming, moon_w_m2_per_tw=b.moon_per_tw,
                   ring_w_m2_per_tw=b.ring_per_tw, targets_dimming_in_all_percent=targets),
        gases=gases(data, b), co2=co2(data, b), pressure=pressure(data, b), cirrus=cirrus(data, b),
        beams=beams(data, b), top_up=top_up(data, b), glow=glow(data, b), albedo=albedo(data, b),
        shield=shield(data, b))


def main(argv=None) -> int:
    product = results()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    print('wrote', OUT.relative_to(ROOT))
    return 0


if __name__ == '__main__':
    sys.exit(main())
