"""Fruit designed for the lunar cycle: how fast a large fruit can grow on the whole plant, when to set it, and how much
of the plant's growth it can take while the plant keeps enough for itself.

    python -m biosphere.canopy.fruit        # writes results/fruit.json (a few minutes on one thread)

The plant is plant.py's base stand, followed step by step through the solar cycle at the chosen climate's land air
temperatures, with the Moon's twilight. The fruit's developmental program is a design variable: the degree-days it
needs from set to harvest (thermal time above a base temperature, up to an upper threshold). The design case, the day
fruit, is engineered to need exactly the degree-days of one sunlit half at its site, so a fruit set at sunrise is ripe
at sunset. The program sweep runs from 120 degree-days to today's watermelon's 560 at the equator, and today's
watermelon is kept as the baseline everywhere.

A fruit's potential dry mass follows the beta sigmoid of Yin et al. (2003) through its development. Each gram of fruit
carbon costs a quarter again in growth respiration, as the plant's own tissue does. The fruit's maintenance is 0.01 g
CH2O per g of dry matter per day at 25 C with Q10 = 2, the value crop models give roots, and every organ but the leaves
in tropical crops (WOFOST, after Spitters et al. 1989 and Penning de Vries et al. 1989).

The plant sets one fruit cohort each cycle. A fruit that takes longer than a cycle overlaps the next one, and the
steady cycle adds up every cohort's demand at each step. The fruit is fed in one of two ways:

- fed: it grows at its potential rate throughout, drawing on the plant's store in the dark;
- daylight: it grows only while the plant's photosynthesis exceeds its upkeep, taking at most that surplus. Its
  development runs on in the dark, so the growth it misses then is lost.

Size and composition are today's watermelon's: 7.5 kg at harvest (7.1-8.3 kg in Noh et al. 2013), 8.6% dry matter
(USDA FoodData Central 167765: 91.4 g of water per 100 g), 42% carbon in that dry matter (from the same composition:
7.55 g of carbohydrate at 42% carbon, 0.61 g of protein at 53%, 0.15 g of fat at 77%). Base temperature 10 C, upper
threshold 37 C (the range reported for watermelon is 10-12.7 C). Growers give 35-45 days from set to harvest with
25-35 C days and 15-20 C nights (about 24 C on average): 490-630 degree-days above 10 C, 560 as the central value.
The fastest dry-matter growth comes at half the development (0.4-0.6 tested), consistent with the third of the harvest
weight that Noh et al. found 15 days after fruit set. Harvest per m2 scales with the carbon the fruit takes, whatever
the size of each fruit.
"""
from __future__ import annotations
from dataclasses import dataclass, replace
import hashlib
import json
from pathlib import Path
import numpy as np

from shared.constants import SYNODIC_MONTH_DAYS
from biosphere.long_night import periodic_storage_requirement
from biosphere.canopy import model as canopy
from biosphere.canopy import plant

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
RESULTS = HERE / 'results'
SCHEMA = 'terluna.biosphere.fruit-carbon/2'
GRAMS = plant.GRAMS
GROWTH_RESPIRATION = plant.GROWTH_RESPIRATION
Q10 = plant.Q10
MAINTENANCE_25 = 0.01 * 12.011 / 30.026     # g C per g of dry matter per day at 25 C (0.01 g CH2O)
CANOPY_C = 5.0 * plant.LEAF_DRY_G_M2 * plant.CARBON_FRACTION     # g C in the base stand's leaves, per m2 of ground
SET_TIMES = 36                              # set times per cycle, every 10 degrees of hour angle
MAX_DAYS = 160.0
SHARES = (0.3, 0.5, 0.6, 0.7, 0.8)
PROGRAMS = (120.0, 160.0, 199.0, 230.0, 260.0, 300.0, 350.0, 420.0, 560.0)     # degree-days swept at the equator
DEGREE_DAYS_RANGE = (490.0, 630.0)
PEAK_RANGE = (0.4, 0.6)
WARMING_C = (3.0, 6.0)
CASES = (('near', '0'), ('near', '30'), ('near', '60'), ('far', '0'))
STRATEGIES = ('earth_like', 'idle_quarter')
POLICIES = ('fed', 'daylight')


@dataclass(frozen=True)
class Fruit:
    harvest_fresh_g: float = 7500.0
    dry_matter: float = 0.086
    carbon_fraction: float = 0.42
    base_c: float = 10.0
    upper_c: float = 37.0
    degree_days: float = 560.0
    peak: float = 0.5

    @property
    def dry_g(self):
        return self.harvest_fresh_g * self.dry_matter

    @property
    def carbon_g(self):
        return self.dry_g * self.carbon_fraction


def filled(s, peak):
    """Share of the harvest dry mass reached at development s (0 at set, 1 at harvest): Yin et al. (2003)."""
    s = np.clip(s, 0.0, 1.0)
    a = 1.0 / (1.0 - peak)
    return (1.0 + a * (1.0 - s)) * s ** a


def cohort(temp, dt_days, fruit, start, warming=0.0, max_days=MAX_DAYS):
    """One fruit set at step `start` of the periodic temperature trace, followed to harvest.

    Returns the cycle step of each step of its life (index), the share of its harvest mass that it can add in each
    (added) and has reached by the end of each (reached), or None if it is not ripe within max_days.
    """
    n = temp.size
    life = int(np.ceil(max_days / dt_days))
    index = (start + np.arange(life)) % n
    rate = np.clip(temp[index] + warming - fruit.base_c, 0.0, fruit.upper_c - fruit.base_c)
    s = np.cumsum(rate * dt_days) / fruit.degree_days
    if s[-1] < 1.0:
        return None
    m = int(np.argmax(s >= 1.0)) + 1
    reached = filled(s[:m], fruit.peak)
    return dict(index=index[:m], added=np.diff(np.concatenate([[0.0], reached])), reached=reached,
                temp=temp[index[:m]] + warming, days=m * dt_days)


def feed(net, dt, fruit, co, load, policy, upkeep=1.0):
    """The steady cycle with `load` fruits per m2 set once per cycle: the store, the plant's own growth and the fruit.

    net is the plant's own net carbon at each step (umol m-2 s-1 of ground), dt the step (s).
    """
    n, dt_days, c = net.size, dt / 86400.0, fruit.carbon_g
    if policy == 'fed':
        added = co['added']
    else:
        want = np.bincount(co['index'], weights=co['added'], minlength=n) * load * c * (1 + GROWTH_RESPIRATION)
        surplus = np.maximum(net, 0.0) * dt * GRAMS
        scale = np.where(surplus > 0, np.minimum(1.0, surplus / np.maximum(want, 1e-300)), 0.0)
        added = co['added'] * scale[co['index']]
    reached = np.cumsum(added)
    maintenance = upkeep * MAINTENANCE_25 * Q10 ** ((co['temp'] - 25.0) / 10.0) * fruit.dry_g * reached * dt_days
    per_step = load * (c * added * (1 + GROWTH_RESPIRATION) + maintenance)
    fruit_flux = np.bincount(co['index'], weights=per_step, minlength=n) / (GRAMS * dt)
    s = periodic_storage_requirement(net - fruit_flux, dt)
    grown = load * c * float(reached[-1])
    left = float(s['cycle_net'] * GRAMS) / (1 + GROWTH_RESPIRATION)
    return dict(load=load, size=float(reached[-1]), harvest_kg=load * fruit.harvest_fresh_g * float(reached[-1]) / 1e3,
                fruit_carbon=grown, fruit_maintenance=load * float(maintenance.sum()), plant_growth=left,
                fruit_share=grown / (grown + left) if grown + left > 0 else 0.0,
                canopy_cycles=CANOPY_C / left if left > 0 else None,
                store=float(s['worst_single_cycle_deficit'] * GRAMS), feasible=left >= 0.0)


def life(f, co):
    """When the fruit lives: harvest hour angle, nights lived through, and its potential growth while the plant is short."""
    dark = f['par'][co['index']] < plant.DARK_UMOL
    nights = int(np.sum(dark[1:] & ~dark[:-1]) + dark[0])
    short = f['net'][co['index']] <= 0.0
    return dict(days_to_harvest=co['days'], harvest_hour_angle=float(f['hour'][co['index'][-1]]), nights=nights,
                growth_while_plant_short=float(co['added'][short].sum()))


def sweep_loads(f, fruit, co, policy, top, upkeep=1.0, points=21):
    """Fruits per m2 per cycle from none to `top`, and the loads that give each fruit share of new tissue carbon."""
    rows = [feed(f['net'], f['dt'], fruit, co, x, policy, upkeep) for x in np.linspace(0.0, top, points)]
    rows = [r for r in rows if r['feasible']]
    share = np.array([r['fruit_share'] for r in rows])
    shares = {}
    for target in SHARES:
        if share.max() < target or np.any(np.diff(share) < 0):
            shares[str(target)] = None
            continue
        x = float(np.interp(target, share, [r['load'] for r in rows]))
        shares[str(target)] = feed(f['net'], f['dt'], fruit, co, x, policy, upkeep)
    return rows, shares


def scenario(f, temp, fruit, warming=0.0, upkeep=1.0, set_times=SET_TIMES):
    """Every set time for one plant flux trace, the best set times and the loads at them."""
    dt_days = f['dt'] / 86400.0
    alone = periodic_storage_requirement(f['net'], f['dt'])
    growth0 = float(alone['cycle_net'] * GRAMS) / (1 + GROWTH_RESPIRATION)
    reference = 0.5 * growth0 / fruit.carbon_g
    starts = np.arange(set_times) * (temp.size // set_times)
    rows = []
    for k in starts:
        co = cohort(temp, dt_days, fruit, int(k), warming)
        row = dict(hour_angle=float(f['hour'][k]), sun_up=bool(f['up'][k]))
        if co is not None:
            row.update(life(f, co))
            for policy in POLICIES:
                r = feed(f['net'], f['dt'], fruit, co, reference, policy, upkeep)
                row[policy] = dict(size=r['size'], store=r['store'], plant_growth=r['plant_growth'])
        rows.append(row)
    ok = [i for i, r in enumerate(rows) if r['sun_up'] and 'fed' in r]
    best = {}
    if ok:
        i_fed = min(ok, key=lambda i: (rows[i]['fed']['store'], -rows[i]['fed']['plant_growth']))
        i_day = min(ok, key=lambda i: (-round(rows[i]['daylight']['size'], 6), rows[i]['daylight']['store']))
        for policy, i in (('fed', i_fed), ('daylight', i_day)):
            co = cohort(temp, dt_days, fruit, int(starts[i]), warming)
            fed_top = growth0 / (fruit.carbon_g + feed(f['net'], f['dt'], fruit, co, 1.0, 'fed', upkeep)
                                 ['fruit_maintenance'] / (1 + GROWTH_RESPIRATION))
            loads, shares = sweep_loads(f, fruit, co, policy, fed_top * (1.0 if policy == 'fed' else 1.6), upkeep)
            peak = co['added'].max() / dt_days
            best[policy] = dict(set_hour_angle=rows[i]['hour_angle'], **life(f, co),
                                potential_peak_growth_g_dry_per_day=peak * fruit.dry_g,
                                potential_peak_growth_kg_fresh_per_day=peak * fruit.harvest_fresh_g / 1e3,
                                loads=loads, shares=shares)
    surplus = np.maximum(f['net'], 0.0)[f['up']].mean() * 86400.0 * GRAMS
    return dict(degree_days=fruit.degree_days, plant_growth_without_fruit=growth0,
                store_without_fruit=float(alone['worst_single_cycle_deficit'] * GRAMS), reference_load=reference,
                daylight_surplus_g_c_per_day=float(surplus),
                daylight_surplus_as_fruit_kg_per_day=float(surplus / (1 + GROWTH_RESPIRATION) / fruit.carbon_g
                                                           * fruit.harvest_fresh_g / 1e3),
                set_times=rows, best=best)


def calendar(f, temp, fruit):
    """Degree-days in the sunlit half, the plant's surplus and the whole cycle, and what one sunlit half can ripen."""
    dt_days = f['dt'] / 86400.0
    rate = np.clip(temp - fruit.base_c, 0.0, fruit.upper_c - fruit.base_c) * dt_days
    sunlit_days = float(f['up'].sum() * dt_days)
    surplus = f['net'] > 0.0
    sunrise = int(np.argmax(f['up'] & ~np.roll(f['up'], 1)))
    co = cohort(temp, dt_days, fruit, sunrise)
    at_sunset = float(co['reached'][min(int(f['up'].sum()), co['reached'].size) - 1]) if co else None
    return dict(sunlit_days=sunlit_days, degree_days_sunlit=float(rate[f['up']].sum()),
                surplus_days=float(surplus.sum() * dt_days), degree_days_in_surplus=float(rate[surplus].sum()),
                degree_days_cycle=float(rate.sum()),
                watermelon_fruit_zone_c_to_ripen_within_sunlit_half=fruit.base_c + fruit.degree_days / sunlit_days,
                watermelon_set_at_sunrise_share_of_harvest_mass_at_sunset=at_sunset,
                watermelon_set_at_sunrise_days_to_harvest=co['days'] if co else None)


def headline(s):
    """The compact numbers of a scenario at a fruit share of 0.6."""
    out = {}
    for policy, b in s['best'].items():
        at = b['shares'].get('0.6')
        out[policy] = dict(set_hour_angle=b['set_hour_angle'], days_to_harvest=b['days_to_harvest'],
                           harvest_hour_angle=b['harvest_hour_angle'], nights=b['nights'],
                           growth_while_plant_short=b['growth_while_plant_short'],
                           potential_peak_growth_kg_fresh_per_day=b['potential_peak_growth_kg_fresh_per_day'],
                           **({k: at[k] for k in ('load', 'harvest_kg', 'size', 'plant_growth', 'store',
                                                  'canopy_cycles')} if at else {}))
    return dict(degree_days=s['degree_days'], store_without_fruit=s['store_without_fruit'], **out)


EVIDENCE = ('A carbon budget for fruit on the whole-plant model\'s stand over the lunar cycle (plant.py: the canopy '
            'model\'s photosynthesis with twilight, leaf and stem-and-root respiration, at the chosen climate\'s land air '
            'temperatures). One fruit cohort is set per cycle and develops in degree-days above 10 C; the degree-days '
            'it needs are a design variable, from the day fruit engineered to ripen within one sunlit half to today\'s '
            'watermelon (560, from growers\' 35-45 days at about a 24 C mean). Potential dry mass follows a beta '
            'sigmoid peaking at half the development (checked against a third of the harvest weight 15 days after '
            'set; Noh et al. 2013), with today\'s watermelon\'s size and composition: 7.5 kg fresh, 8.6% dry matter, '
            '42% carbon (USDA), growth respiration of a quarter and maintenance of 0.01 g CH2O per g dry matter per day '
            'at 25 C (Q10 = 2). The day fruit is a design requirement: whether a fruit\'s program can be compressed '
            'to it is not established. Development follows temperature alone, and the fruit sets whenever it is '
            'triggered. Flowering, fruit set and abortion, seeds, water, fruit quality and the crop\'s own canopy (the '
            'stand is plant.py\'s leaf area index 5) are not modelled.')
READING_RULE = ('Carbon is g C per m2 of ground per solar cycle; a load is fruits of 7.5 kg set per m2 per cycle, one '
                'cohort per cycle, so a fruit that takes longer than a cycle overlaps the next. Hour angles are degrees '
                'from local noon (-90 sunrise, 90 sunset). results[side][band][strategy] (strategies as in plant.json) '
                'holds day_fruit, the design case whose degree-days equal the site\'s sunlit half, and watermelon, today\'s '
                'fruit. In each, set_times gives for every set hour angle the days to harvest, the harvest hour angle, the '
                'nights lived through, the share of potential growth that falls while the plant is in deficit '
                '(growth_while_plant_short), and at the reference load (fruit carbon half the plant\'s growth without '
                'fruit) the store and the plant\'s remaining growth when the fruit is fed throughout (fed), and its size '
                'as a share of full size when it grows only on the daylight surplus (daylight). best is the set time with '
                'the Sun up that needs the smallest store (fed) or grows the largest fruit (daylight), with loads swept '
                'from none upward and interpolated to fruit shares of new tissue carbon (shares). harvest_kg is fresh '
                'fruit per m2 per cycle; plant_growth the plant\'s own new tissue; canopy_cycles the cycles that growth '
                'would take to build a whole new canopy (211.5 g C). programs sweeps the degree-days at the equatorial '
                'near side (headline numbers at a fruit share of 0.6). calendar gives the degree-days of the sunlit half, '
                'of the plant\'s surplus and of the cycle. earth is the same stand with today\'s watermelon under Earth\'s '
                'sky with the same temperatures over 24 hours, harvests per 29.53 days. sensitivity varies one setting at '
                'a time at the equatorial near side, Earth-like plant.')


def run():
    hour_centres, traces, climate_sha = plant.temperatures()
    tables = plant.Tables()
    equator = plant.along_cycle(hour_centres, traces['near', '0'])
    r25 = plant.calibrate(tables, equator)
    fruit = Fruit()
    night = {name: plant.STRATEGIES[name][0] for name in STRATEGIES}
    results, calendars, flux, day_fruits = {}, {}, {}, {}
    for side, band in CASES:
        temp = plant.along_cycle(hour_centres, traces[side, band])
        for name in STRATEGIES:
            flux[side, band, name] = plant.fluxes(tables, 'moon', float(band), temp, plant.MOON_DAY_S, r25,
                                                  night=night[name]), temp
        cal = calendar(flux[side, band, 'earth_like'][0], temp, fruit)
        calendars.setdefault(side, {})[band] = cal
        day_fruits[side, band] = replace(fruit, degree_days=cal['degree_days_sunlit'])
        for name in STRATEGIES:
            f = flux[side, band, name][0]
            results.setdefault(side, {}).setdefault(band, {})[name] = dict(
                day_fruit=scenario(f, temp, day_fruits[side, band]), watermelon=scenario(f, temp, fruit))
    programs = {name: [headline(scenario(*flux['near', '0', name], replace(fruit, degree_days=dd)))
                       for dd in PROGRAMS] for name in STRATEGIES}
    f_moon, t_moon = flux['near', '0', 'earth_like']
    day = day_fruits['near', '0']
    earth_f = plant.fluxes(tables, 'earth', 0.0, equator, 86400.0, r25)
    earth = scenario(earth_f, equator, fruit, set_times=8)
    for b in earth['best'].values():
        for r in list(b['shares'].values()) + b['loads']:
            if r:
                r['harvest_kg_per_29_53_days'] = r['harvest_kg'] * SYNODIC_MONTH_DAYS
                r['plant_growth_per_29_53_days'] = r['plant_growth'] * SYNODIC_MONTH_DAYS
    r40 = plant.calibrate(tables, equator, 0.4)
    f40 = plant.fluxes(tables, 'moon', 0.0, t_moon, plant.MOON_DAY_S, r40)
    sensitivity = dict(day_fruit={}, watermelon={})
    for pk in PEAK_RANGE:
        sensitivity['day_fruit'][f'peak_{pk:g}'] = headline(scenario(f_moon, t_moon, replace(day, peak=pk)))
    sensitivity['day_fruit']['fruit_maintenance_half'] = headline(scenario(f_moon, t_moon, day, upkeep=0.5))
    sensitivity['day_fruit']['carbon_use_efficiency_0.4'] = headline(scenario(f40, t_moon, day))
    for dd in DEGREE_DAYS_RANGE:
        sensitivity['watermelon'][f'degree_days_{dd:g}'] = headline(scenario(f_moon, t_moon,
                                                                             replace(fruit, degree_days=dd)))
    for w in WARMING_C:
        sensitivity['watermelon'][f'fruit_zone_warmer_{w:g}c'] = headline(scenario(f_moon, t_moon, fruit, warming=w))
    sensitivity['watermelon']['fruit_maintenance_half'] = headline(scenario(f_moon, t_moon, fruit, upkeep=0.5))
    sensitivity['watermelon']['carbon_use_efficiency_0.4'] = headline(scenario(f40, t_moon, fruit))
    digest = lambda p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()[:16]
    return dict(
        schema=SCHEMA,
        producer=dict(domain='biosphere', files={p: digest(p) for p in (
            'biosphere/canopy/fruit.py', 'biosphere/canopy/plant.py', 'biosphere/canopy/model.py',
            'biosphere/canopy/radiation.py', 'biosphere/canopy/leaf.py', 'biosphere/canopy/leaf_optics.py',
            'biosphere/long_night.py')},
            inputs=dict(**{f'{k}_sha256': v for k, v in tables.shas.items()}, climatology_sha256=climate_sha,
                        climatology_run=plant.CLIMATE_RUN, climatology_years=plant.CLIMATE_YEARS)),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        units=dict(carbon='g C m-2 of ground per solar cycle', load='fruits of 7.5 kg m-2 per cycle',
                   harvest='kg fresh fruit m-2 per cycle', temperature='C', hour_angle='degrees from local noon'),
        settings=dict(fruit=dict(fruit.__dict__, dry_g=fruit.dry_g, carbon_g=fruit.carbon_g),
                      day_fruit_degree_days={f'{s}/{b}': v.degree_days for (s, b), v in day_fruits.items()},
                      maintenance_g_ch2o_per_g_dry_per_day_25c=0.01, q10=Q10, growth_respiration=GROWTH_RESPIRATION,
                      canopy_carbon_g_m2=CANOPY_C, stem_root_respiration_25c=r25, set_times=SET_TIMES,
                      max_days=MAX_DAYS, shares=SHARES, programs=PROGRAMS,
                      sources=dict(
                          harvest_mass='Noh et al. 2013, Physiol. Mol. Biol. Plants 19:509-514: 7.1-8.3 kg per plant at '
                                       'harvest and 2.3-3.0 kg 15 days after fruit set (Sambok-gul, Korea)',
                          composition='USDA FoodData Central 167765, watermelon, raw: 91.4 g water, 7.55 g '
                                      'carbohydrate, 0.61 g protein, 0.15 g fat per 100 g',
                          development='a seed company\'s cultivation guide (Origene Seeds): 35-45 days from fruit set '
                                      'to maturity with 25-35 C days and 15-20 C nights; ripe (red flesh) 36 days after '
                                      'pollination in the field at Lane, Oklahoma (Wechter et al. 2008, BMC Genomics '
                                      '9:275)',
                          thresholds='watermelon base temperatures of 10-12.7 C and an upper threshold of 37 C '
                                     '(review of growing-degree-day thresholds, Agricultural Water Management 2025)',
                          growth_curve='Yin et al. 2003, Annals of Botany 91:361-371 (beta sigmoid)',
                          maintenance='WOFOST system description: typical maintenance coefficients 0.03 (leaves), '
                                      '0.015 (stems), 0.01 (roots) g CH2O per g per day at 25 C (Spitters et al. '
                                      '1989), 0.01 for organs other than leaves in tropical crops (Penning de Vries '
                                      'et al. 1989)')),
        calendar=calendars, results=results, programs=programs, earth=earth, sensitivity=sensitivity)


def _clean(obj):
    if isinstance(obj, float) and not np.isfinite(obj):
        return None
    if isinstance(obj, dict):
        return {k: _clean(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_clean(v) for v in obj]
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    return obj


def main():
    product = canopy._round(_clean(run()))
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / 'fruit.json').write_text(json.dumps(product, indent=1) + '\n')
    print('wrote', RESULTS / 'fruit.json')


if __name__ == '__main__':
    main()
