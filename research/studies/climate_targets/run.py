"""What each target climate takes, across the domains: the warm Moon the GCM settles on with its longwave refitted by
height and ice optics for its cold cloud, the design temperature the plan was built for, and a target between them,
with the heat the plan puts on the Moon beyond its sunlight, which every target answers.

    python -m research.studies.climate_targets.run        # under a second

It binds the GCM extract (extract.py, run with the climate venv) and committed products of the domains: the gases,
glow and dimming study, the fleet's glow, the deep column's gas forcing and saturated limit, the light study, the power
screen, resources' first screen and geography's groundwater uptake. It is arithmetic on them; no model is run.
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
OUT = HERE / 'results' / 'climate_targets.json'
SCHEMA = 'terluna.research.climate-targets/1'
CODE = ('research/studies/climate_targets/run.py', 'research/studies/climate_targets/extract.py')
EXTRACT = 'research/runs/climate_targets/extract.json'
GASES_GLOW = 'research/studies/gases_glow_dimming/results/gases_glow_dimming.json'
GLOW = 'research/studies/solar_shield_array/results/glow.json'
GAS_FORCING = 'atmosphere/radiative_convective/results/reduced_gases.json'
COLUMN = 'atmosphere/radiative_convective/results/summary.json'
LIGHT = 'research/studies/light_limitation/results/light_limitation.json'
POWER = 'provisioning/results/power.json'
RESOURCES = 'resources/results/first_screen.json'
GROUNDWATER = 'geography/results/groundwater_uptake.json'
INPUTS = {EXTRACT: 'terluna.research.climate-targets-extract/1', GASES_GLOW: 'terluna.research.gases-glow-dimming/1',
          GLOW: 'terluna.research.fleet-glow/1', GAS_FORCING: 'terluna.atmosphere.reduced-gas-forcing/1',
          COLUMN: None, LIGHT: 'terluna.research.light-limitation/1', POWER: 'terluna.provisioning.power/1',
          RESOURCES: 'terluna.resources.first-screen/1', GROUNDWATER: 'terluna.geography.groundwater-uptake/1'}
DESIGN, TEST = 'A28_corrected_970', 'A28_corrected_970_lw_ice'
DESIGN_SCALE = 0.970                 # the design case's sunlight, as a share of what the titania stack passes
TARGETS_K = (295.0, 300.0)           # the design temperature and one between it and where the corrected GCM settles
# Windows of the test's years over which the top of the air's gain is fitted against the ground's temperature; the
# first years carry the start's quick adjustment.
FITS = ((5, 19), (0, 19))
HELD_PAIR = '1.8663+0.1'             # the gases' held pair, ppm of methane and of N2O (domain synthesis, gases.held)
KG_PER_KM3 = 1.0e12                  # water
FIRST_MILLENNIUM = ('0-10', '10-100', '100-1000')
LATER = '1000-10000'
EVIDENCE = ('Arithmetic on two ExoPlaSim runs and committed products. The two runs share the design sunlight (scale '
            '0.970); the design case flew PlaSim\'s earlier longwave and liquid optics for all cloud, the test the '
            'longwave refitted by height and ice optics for cloud colder than -38 C, and neither carries the gases, '
            'the glow or their dimming. The dimming a target needs is the test\'s top-of-air gain per kelvin times '
            'the kelvins to shed, at the test\'s albedo; the climate at a target between the two runs is interpolated '
            'linearly in the ground\'s temperature, so it mixes the two physics and carries neither the refit\'s high '
            'cloud at the cooler end nor the dimming\'s own cut to the rain. The heat beyond the sunlight takes the '
            'GCM\'s efficacies for the gases and the glow and an efficacy of 1 for heat released at the ground.')
READING_RULE = ('W/m2 are means over the whole Moon; absorbed sunlight is what the Moon takes in at the top of the air. '
                'A dimming is the share of the sunlight reaching the Moon that the fleet takes away, at the test\'s '
                'albedo, and the sunlight scale the GCM would fly for it is the design scale times what remains. '
                'Heat at the ground is in the 3-day means at the lowest model level, about 0.9 km up, pooled over '
                'land cells and means; deserts are land under 0.5 mm of rain a day. Pairs are [low, high] ends: the '
                'two fits for the feedback, the light study\'s and Earth\'s models\' ends for the costs.')


def digest(path) -> str:
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


def load():
    data = {p: json.loads((ROOT / p).read_text()) for p in INPUTS}
    for p, schema in INPUTS.items():
        if schema is not None and data[p].get('schema') != schema:
            raise ValueError(f'{p}: schema {data[p].get("schema")} is not {schema}')
    return data


def fit(progress, window):
    """The top of the air's net gain (W/m2) against the ground's temperature (K) over the window's years: the feedback
    (W/m2 per K, positive restoring) and the temperature at which the gain would vanish."""
    rows = [y for y in progress if window[0] <= y['year'] <= window[1]]
    slope, intercept = np.polyfit([y['ground_k'] for y in rows], [y['top_of_air_net_w_m2'] for y in rows], 1)
    return dict(years=list(window), feedback_w_m2_k=float(-slope), balance_k=float(-intercept / slope))


def crossing(ladder, celsius):
    """The height (km) at which the land's mean air falls through a temperature, linear between the ladder's heights."""
    for lo, hi in zip(ladder, ladder[1:]):
        if (lo['air_c'] - celsius) * (hi['air_c'] - celsius) <= 0 and lo['air_c'] != hi['air_c']:
            f = (lo['air_c'] - celsius) / (lo['air_c'] - hi['air_c'])
            return lo['height_km'] + f * (hi['height_km'] - lo['height_km'])
    return None


# Shares of time past a threshold do not interpolate: they follow the tail of the distribution.
RUNS_ONLY = ('wet_bulb_above_28', 'wet_bulb_above_31', 'wet_bulb_above_28_tropics')


def heights(ladder):
    return dict(height_of_20c_km=crossing(ladder, 20.0), freezing_level_km=crossing(ladder, 0.0))


def climate(run):
    """The figures the targets interpolate, from one run's years 10-19. A crossing height is None where the lowest
    level is already colder."""
    h = run['heat']
    return dict(ground_k=run['ground_k'], rain_mm_day=run['rain_mm_day'], land_rain_mm_day=run['land']['all']['rain_mm_day'],
                land_desert=run['land']['all']['desert'],
                wet_bulb_median_c=h['land']['wet_bulb']['median'], wet_bulb_p99_c=h['land']['wet_bulb']['p99'],
                wet_bulb_above_28=h['land']['above_28'], wet_bulb_above_31=h['land']['above_31'],
                wet_bulb_above_28_tropics=h['bands']['0-15']['above_28'],
                land_air_median_c=h['land']['air']['median'], sea_median_c=run['seas']['all']['median'],
                ladder=[dict(height_km=x['height_km'], air_c=x['air_c'], pressure_atm=x['pressure_atm'])
                        for x in run['sky_over_land']],
                **heights(run['sky_over_land']))


def interpolate(a, b, ground_k):
    """A target's climate between the two runs, linear in the ground's temperature, the sky ladder height by height."""
    f = (ground_k - a['ground_k']) / (b['ground_k'] - a['ground_k'])
    out = {k: a[k] + f * (b[k] - a[k]) for k in a
           if k not in RUNS_ONLY and k not in ('ladder', 'height_of_20c_km', 'freezing_level_km')}
    out['ladder'] = [dict(height_km=x['height_km'], air_c=x['air_c'] + f * (y['air_c'] - x['air_c']),
                          pressure_atm=x['pressure_atm'] + f * (y['pressure_atm'] - x['pressure_atm']))
                     for x, y in zip(a['ladder'], b['ladder'])]
    out.update(heights(out['ladder']))
    out['ground_k'] = ground_k
    return out, f


def heat_beyond_sunlight(data, absorbed, area):
    """The heat the plan adds on the Moon besides its sunlight, each as the absorbed sunlight that answers it."""
    gases_glow, glow = data[GASES_GLOW], data[GLOW]['results']['recommended']
    efficacy = gases_glow['glow']['held']['efficacy']
    albedo = data[GLOW]['design']['planetary_albedo']
    gas = data[GAS_FORCING]['moon_at_design']['both'][HELD_PAIR]['allsky_w_m2']
    infrared = glow['glow_W_m2']['global'] + glow['infrared_return_W_m2']['global']
    sunlight = glow['sunlight_return_W_m2']['global']
    link = data[POWER]['absorption']['link']
    power_tw = [i + a for i, a in zip(link['intake_tw'], link['absorbed_in_the_air_tw'])]

    def item(absorbed_w_m2, **kw):
        return dict(kw, absorbed_sunlight_w_m2=absorbed_w_m2,
                    dimming_percent=[100.0 * x / absorbed for x in absorbed_w_m2])

    out = dict(
        gases=item([gas, gas], held_pair_ppm=HELD_PAIR, forcing_w_m2=gas,
                   basis='all-sky forcing of the held pair in the deep column, efficacy 1 as the GCM found it '
                         '(its 1.69% over the absorbed sunlight held the design temperature)'),
        glow=item([efficacy * infrared + (1.0 - albedo) * sunlight] * 2, answer=data[GLOW]['results']['recommended']['answer'],
                  infrared_w_m2=infrared, sunlight_back_w_m2=sunlight, efficacy=efficacy,
                  basis='the recommended answer, both halves rolled 80 degrees on the re-flown layout, at the GCM\'s '
                        'held-branch efficacy for its infrared and the albedo for its returned sunlight'),
        beamed_power=item([tw * 1e12 / area for tw in power_tw], heat_tw=power_tw,
                          basis='the power screen\'s intake at the loads and the deep air\'s share of the beams, all '
                                'released at or near the ground'))
    out['answered_in_every_target'] = dict(
        absorbed_sunlight_w_m2=[sum(out[k]['absorbed_sunlight_w_m2'][e] for k in ('gases', 'glow', 'beamed_power'))
                                for e in (0, 1)])
    out['answered_in_every_target']['dimming_percent'] = [100.0 * x / absorbed for x in
                                                          out['answered_in_every_target']['absorbed_sunlight_w_m2']]
    # The water delivered to balance the thawing crust, if its fall into the Moon's well is spent on the Moon.
    fall = data[RESOURCES]['delivery']['fall_mj_per_kg']
    energy = [min(fall['exobase']) * 1e6, fall['far_terminal'] * 1e6]
    year = K.JULIAN_YEAR_DAYS * K.JULIAN_DAY
    top = {c: data[GROUNDWATER]['cases'][c]['moon']['top_up_km3_yr'] for c in ('slow', 'central', 'fast')}

    def w_m2(km3, e):
        return km3 * KG_PER_KM3 / year * e / area

    first = [min(top['central'][w] for w in FIRST_MILLENNIUM), max(top['central'][w] for w in FIRST_MILLENNIUM)]
    fast_peak = max(top['fast'][w] for w in FIRST_MILLENNIUM)
    out['water_top_up'] = dict(
        fall_mj_per_kg=[e / 1e6 for e in energy],
        first_millennium_central_km3_yr=first, to_10000_years_central_km3_yr=top['central'][LATER],
        first_millennium_fast_peak_km3_yr=fast_peak,
        first_millennium_w_m2=[w_m2(first[0], energy[0]), w_m2(first[1], energy[1])],
        to_10000_years_w_m2=[w_m2(top['central'][LATER], energy[0]), w_m2(top['central'][LATER], energy[1])],
        fast_peak_w_m2=[w_m2(fast_peak, energy[0]), w_m2(fast_peak, energy[1])],
        basis='geography\'s groundwater top-up through the first ten thousand years; the fall from the exobase to the '
              'far cislunar terminal (resources\' first screen); all of it spent on the Moon, the case without energy '
              'recovery on the way down')
    out['water_top_up']['first_millennium_dimming_percent'] = [100.0 * x / absorbed
                                                               for x in out['water_top_up']['first_millennium_w_m2']]
    return out


def results():
    data = load()
    ex = data[EXTRACT]
    design, test = ex['runs'][DESIGN], ex['runs'][TEST]
    toa = test['top_of_air_w_m2']
    absorbed, incoming = toa['absorbed'], toa['incoming']
    albedo = toa['reflected'] / incoming
    area = 4.0 * math.pi * K.MOON_RADIUS ** 2
    fits = [fit(ex['test_progress'], w) for w in FITS]
    lam = sorted(f['feedback_w_m2_k'] for f in fits)
    bal = sorted(f['balance_k'] for f in fits)
    saturated = data[COLUMN]['moon_1.2atm']['saturated']['max_olr']
    conversion = data[GLOW]['design']['conversion']
    light = data[LIGHT]['cost_through_light']['per_percent_dimming']
    rain_cut = [-data[GASES_GLOW]['runs']['gases']['rain_per_percent_dimming']['land'],
                data[LIGHT]['rain_pathway']['precipitation_cut_per_percent_dimming'][1]]
    heat = heat_beyond_sunlight(data, absorbed, area)
    every = heat['answered_in_every_target']['absorbed_sunlight_w_m2']
    a, b = climate(design), climate(test)
    targets = []
    for name, target in [('no further dimming', None)] + [(f'{t:g} K', t) for t in TARGETS_K]:
        if target is None:
            shed = [0.0, 0.0]
            state, frac = dict(b), 1.0
        else:
            # low end: the weaker feedback over the cooler balance; high end: the stronger over the warmer
            shed = [lam[0] * (bal[0] - target), lam[1] * (bal[1] - target)]
            state, frac = interpolate(a, b, target)
        dim = [100.0 * s / absorbed for s in shed]
        total = [100.0 * (s + e) / absorbed for s, e in zip(shed, every)]
        taken_tw = [d / 100.0 * incoming * area / 1e12 for d in dim]
        rain = state['land_rain_mm_day']
        targets.append(dict(
            name=name, ground_k=state['ground_k'] if target is not None else test['ground_k'],
            settles_k=bal if target is None else [target, target],
            absorbed_sunlight_to_shed_w_m2=shed, dimming_percent=dim,
            with_the_heat_beyond_sunlight_percent=total,
            sunlight_scale=[DESIGN_SCALE * (1.0 - t / 100.0) for t in total][::-1],
            sunlight_taken_tw=taken_tw, electricity_tw=[taken_tw[0] * conversion[0], taken_tw[1] * conversion[1]],
            growth_cost_through_light_percent=[light[0] * dim[0], light[1] * dim[1]],
            land_rain_mm_day=rain, land_rain_with_the_dimmings_own_cut_mm_day=[
                rain * (1.0 - rain_cut[1] / 100.0 * dim[1]), rain * (1.0 - rain_cut[0] / 100.0 * dim[0])],
            climate=state, interpolated_share_of_the_way=frac,
            margin_to_the_saturated_column_w_m2=[saturated - (absorbed - s) for s in shed],
            margin_if_the_heat_beyond_sunlight_goes_unanswered_w_m2=[saturated - (absorbed - s) - e
                                                                    for s, e in zip(shed, every[::-1])]))
    return dict(
        schema=SCHEMA,
        producer=dict(files={c: digest(c) for c in CODE}, inputs={p: digest(p) for p in INPUTS},
                      constants=dict(MOON_RADIUS=K.MOON_RADIUS, JULIAN_YEAR_DAYS=K.JULIAN_YEAR_DAYS,
                                     JULIAN_DAY=K.JULIAN_DAY)),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        runs=dict(design=dict(run=DESIGN, **climate(design), top_of_air_w_m2=design['top_of_air_w_m2']),
                  test=dict(run=TEST, **climate(test), top_of_air_w_m2=toa, planetary_albedo=albedo)),
        settled=dict(fits=fits, feedback_w_m2_k=lam, balance_k=bal,
                     saturated_column_max_outgoing_w_m2=saturated,
                     clear_sky_outgoing_over_the_saturated_column_w_m2=toa['outgoing_infrared_clear_sky'] - saturated),
        heat_beyond_sunlight=heat,
        targets=targets,
        costs_per_percent=dict(growth_through_light_percent=light, land_rain_percent=rain_cut))


def main(argv=None) -> int:
    product = results()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    for t in product['targets']:
        print(f"{t['name']:20s} shed {t['absorbed_sunlight_to_shed_w_m2'][0]:5.1f}-{t['absorbed_sunlight_to_shed_w_m2'][1]:5.1f} W/m2, "
              f"dimming {t['dimming_percent'][0]:4.1f}-{t['dimming_percent'][1]:4.1f}% "
              f"(+heat {t['with_the_heat_beyond_sunlight_percent'][0]:4.1f}-{t['with_the_heat_beyond_sunlight_percent'][1]:4.1f}%), "
              f"scale {t['sunlight_scale'][0]:.3f}-{t['sunlight_scale'][1]:.3f}, "
              f"wet-bulb {t['climate']['wet_bulb_median_c']:.1f}/{t['climate']['wet_bulb_p99_c']:.1f} C, "
              f"desert {100*t['climate']['land_desert']:.0f}%, sea {t['climate']['sea_median_c']:.1f} C, "
              f"margin {t['margin_to_the_saturated_column_w_m2'][0]:.1f}-{t['margin_to_the_saturated_column_w_m2'][1]:.1f} W/m2")
    print('wrote', OUT.relative_to(ROOT))
    return 0


if __name__ == '__main__':
    sys.exit(main())
