"""What each target climate takes, across the domains: the hot world the GCM settles on with its longwave refitted by
height and ice optics for its cold cloud, the lock at the design temperature the plan was built for, and the resort
between them, with the heat the plan puts on the Moon beyond its sunlight, which every target answers.

    python -m research.studies.climate_targets.run        # under a second

It binds the GCM extract (extract.py, run with the climate venv) and committed products of the domains: the gases,
glow and dimming study, the fleet's glow, the deep column's gas forcing and saturated limit, the light study, the power
screen, resources' first screen and geography's groundwater uptake. It is arithmetic on them; no model is run.

Beside the three targets it gives the dimming for every ground temperature from 295 to 307 K; the glow's efficacy as a
range, from the gases, glow and dimming study's held branch of 13 October and from the glow pair on the corrected
radiation of 17 October, carried through the heat beyond the sunlight; the land's heat in the day's 3-day means and
in the night's; and how much the 3-day means smooth the heat's extremes, from the design case's years written every
6.5 hours, carried to the targets as a sensitivity. For the summit tower it gives the air over the summit at every
target, the warmest target at which the tower's top keeps a winter, bare and with margins, the check against the
metropolis brief and against the test's own warming; and for the resort it sets the candidates from 298 to 303 K
against the tower's winter, comfort at the ground, the dimming, rain and deserts, ecology's waters' oxygen (read in
the tree) and the margin under a saturated column.
"""
from __future__ import annotations
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

from research.studies.climate_targets import extract as X
from shared import constants as K

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'climate_targets.json'
SCHEMA = 'terluna.research.climate-targets/5'
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
SUMMIT_TOWER = 'research/studies/summit_tower/results/summit_tower.json'
FORM = 'research/studies/summit_tower/results/summit_tower_form.json'
SITE_WINDS = 'climate/results/gcm/site_winds_A28_dim5_moon_summit.json'     # the metropolis brief's air
INPUTS = {EXTRACT: 'terluna.research.climate-targets-extract/5', GASES_GLOW: 'terluna.research.gases-glow-dimming/1',
          GLOW: 'terluna.research.fleet-glow/1', GAS_FORCING: 'terluna.atmosphere.reduced-gas-forcing/1',
          COLUMN: None, LIGHT: 'terluna.research.light-limitation/1', POWER: 'terluna.provisioning.power/1',
          RESOURCES: 'terluna.resources.first-screen/1', GROUNDWATER: 'terluna.geography.groundwater-uptake/1',
          SUMMIT_TOWER: 'terluna.research.summit-tower/1', FORM: 'terluna.research.summit-tower-form/1',
          SITE_WINDS: 'terluna.climate.gcm-site-winds/1'}
# Ecology's waters' oxygen in the three targets, on the corrected lake model with the wind's deep mixing (read in-tree
# since the fourth joint integration; the lake studies keep the climate targets they ran on, so the two do not bind
# each other).
LAKE_OXYGEN = dict(path='research/studies/lake_oxygen/results/lake_oxygen.json',
                   schema='terluna.research.lake-oxygen/2')
OLD_DESIGN = 'A28_dim5_moon'
BAND_K = tuple(298.0 + 0.5 * i for i in range(11))       # the resort's candidates, 298 to 303 K
MARGINS_M = (1500.0, 2000.0)     # winter's depth below the top: the design's winter parks span 13.5-15 km
# The metropolis brief's figures for the air over the summit (habitation/summit_metropolis/README.md)
BRIEF = dict(summit_c=12.0, lapse_k_per_km=0.9, disks_from_km=12.3, disks_c=1.0, winter_parks_km=(13.5, 15.0),
             winter_parks_c=(0.0, -1.0), freezing_above_summit_km=13.9, top_c=-9.0)
DESIGN, TEST = 'A28_corrected_970', 'A28_corrected_970_lw_ice'
# The glow pair on the corrected radiation: the gases and the glow with 13 October's 4.22% of added dimming, and with
# its 3.52%.
WW, HELD = 'A28_corrected_gases_glow', 'A28_corrected_gases_glow_held'
DESIGN_SCALE = 0.970                 # the design case's sunlight, as a share of what the titania stack passes
TARGETS_K = (295.0, 300.0)           # the design temperature and one between it and where the corrected GCM settles
CURVE_K = tuple(295.0 + 0.5 * i for i in range(25))      # every half kelvin from the design temperature to 307 K
# Windows of the test's years over which the top of the air's gain is fitted against the ground's temperature; the
# first years carry the start's quick adjustment.
FITS = ((5, 19), (0, 19))
HELD_PAIR = '1.8663+0.1'             # the gases' held pair, ppm of methane and of N2O (domain synthesis, gases.held)
KG_PER_KM3 = 1.0e12                  # water
FIRST_MILLENNIUM = ('0-10', '10-100', '100-1000')
LATER = '1000-10000'
EVIDENCE = ('Arithmetic on ExoPlaSim runs and committed products. The two runs the targets stand on share the design '
            'sunlight (scale 0.970); the design case flew PlaSim\'s earlier longwave and liquid optics for all cloud, '
            'the test the longwave refitted by height and ice optics for cloud colder than -38 C, and neither carries '
            'the gases, the glow or their dimming. The dimming a target needs is the test\'s top-of-air gain per '
            'kelvin times the kelvins to shed, at the test\'s albedo; the climate at a target between the two runs is '
            'interpolated linearly in the ground\'s temperature, so it mixes the two physics and carries neither the '
            'refit\'s high cloud at the cooler end nor the dimming\'s own cut to the rain, and a target at either end '
            'takes that run\'s own climate. The heat beyond the sunlight takes the GCM\'s efficacies for the gases and '
            'the glow and an efficacy of 1 for heat released at the ground. The glow\'s efficacy is a range, from the '
            'held branch of 13 October on the earlier radiation and from the glow pair on the corrected radiation of '
            '17 October (the design case\'s own radiation, before the longwave refit). The pair fixes only the total '
            'the gases and the glow take together, its held branch carried to the design case\'s temperature by the '
            'pair\'s own response to its dimming; the gases take their forcing, as on 13 October, and the glow the '
            'rest. Land and water are the model\'s own mask for the heat, the seas and the sky; the land\'s rain and '
            'deserts count every cell\'s land with that cell\'s own rain, the model\'s own land cells giving the dry '
            'bound, since the grid cannot resolve the coasts. The 2 m air is the lowest level\'s air carried down '
            'the dry adiabat, and the 2 m wet-bulb temperature carries the lowest level\'s humidity down with it, '
            'capped at saturation: the reading of a well-mixed lowest kilometre, as by day. At night the air near the '
            'ground cools below it in calm air, so the night\'s figures at 2 m read warm. The output\'s humidity is a '
            'snapshot at each mean\'s first step (ExoPlaSim 3.4.2, plasim.f90 line 3218), so each mean\'s humidity is '
            'that of its start and its end, and each run\'s last mean, with no next record, leaves the heat. The '
            'wet-bulb temperature is capped at the air\'s own where the air is supersaturated. The day\'s and the '
            'night\'s 3-day means are told apart by the output\'s mean cosine of the Sun\'s zenith angle against the '
            'Sun\'s synodic travel through each mean, a cut exact where the Sun\'s declination is zero. How much the '
            '3-day means smooth the extremes '
            'comes from the design case\'s two years written every 6.5 hours alone; carried to the warmer targets it '
            'is a sensitivity. The summit tower\'s air is the free air over the summit in the 3-day means, read at '
            'the T21 cells round it on the model\'s own ground, which T21 sets 4.0 km below the summit, by height above '
            'sea level with inverse-square gravity, interpolated linearly in geopotential across the model\'s layers; '
            'a target between the runs takes it interpolated by the ground\'s temperature, which mixes the two physics. '
            'The resort\'s oxygen interpolates ecology\'s three targets.')
READING_RULE = ('W/m2 are means over the whole Moon; absorbed sunlight is what the Moon takes in at the top of the air. '
                'A dimming is the share of the sunlight reaching the Moon that the fleet takes away, at the test\'s '
                'albedo, and the sunlight scale the GCM would fly for it is the design scale times what remains. '
                'Heat is in the 3-day means, each about a tenth of the lunar day, at the lowest model level, about '
                '0.9 km up, and at 2 m. Land and water are the model\'s own mask; land figures are weighted by land '
                'area and the seas\' by water area, quantiles and shares of time included. The day\'s means are those '
                'the Sun is up for more than half of; day_ and night_ figures in a climate are the land\'s in them, '
                'and a tenth\'s day and afternoon are the shares of its means in the day and in the day\'s afternoon. '
                'Each place\'s tenths are ranked by one quantity (by_wet_bulb, by_air) with the others over the same '
                'means. Deserts are land under 0.5 mm of rain a day; a run\'s land.every_cell counts every cell\'s '
                'land fraction times its area with the cell\'s own rain, the headline, and land.model_land the '
                'model\'s own land cells alone, the dry bound, which a climate carries as model_land_rain_mm_day and '
                'model_land_desert. A climate\'s freezing_level_km and '
                'height_of_20c_km are where the land\'s mean profile, a land-area mean at each height, crosses 0 and '
                '20 C; a run\'s '
                'freezing_level_over_land_km is the median of the land\'s column crossings. Pairs are [low, high] '
                'ends: the two fits for the feedback, the glow\'s two efficacies, the power screen\'s two ends and the '
                'light study\'s and Earth\'s models\' ends for the costs; figures that carry several of them pair the '
                'ends that give the least and the most. glow_efficacy holds the corrected radiation\'s efficacy and '
                'how it was found, dimming in per cent of the design case\'s sunlight; gases_and_glow_w_m2 is the '
                'total the pair fixes, and with_the_gases_at_their_points the other split of it. dimming_curve gives '
                'the dimming for each ground temperature, negative where a target lies above a fit\'s balance. '
                'smoothing compares the land\'s heat in the 6.5-hour means with the 3-day means their elevens make '
                'over the same two years; ratio_above_the_median is the extreme\'s height above the median at 6.5 '
                'hours over its height in the 3-day means, and a target\'s at_6_5_hours applies it to that target\'s '
                '99th percentiles. wet_bulb_per_kelvin_of_air is the change in the wet-bulb temperature for a kelvin '
                'of air at a fixed humidity, at a run\'s day medians at 2 m and the land\'s 5th and 95th percentile '
                'surface pressures. tower.targets gives, for every half kelvin, heights above the summit in km, the '
                'mean, warm and cold tenths of the air (C) at the top, 24 km up, and at the decks (decks_c, by height above '
                'the summit), and the tower\'s length below freezing in the mean; warmest_targets are ground '
                'temperatures (K) at which the top, or the height a margin below it, reaches 0 C. resort_band.rows are '
                'the resort\'s candidates; their oxygen comes from ecology\'s lake_oxygen product in the tree: the '
                'seas\' saturation at their median warmth, and the deep lakes\' oligotrophic deepest tenth (mg/L) and '
                'the eutrophic shallows\' share of nights below 2 mg/L, each [least, most] over its three lakes, '
                'interpolated by the ground\'s temperature.')


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


def shed(fits, target):
    """The absorbed sunlight (W/m2) a target sheds, each fit's own feedback times its own balance's distance from the
    target, as [low, high]."""
    return sorted(f['feedback_w_m2_k'] * (f['balance_k'] - target) for f in fits)


def crossing(ladder, celsius):
    """The height (km) at which the land's mean air falls through a temperature, linear between the ladder's heights."""
    for lo, hi in zip(ladder, ladder[1:]):
        if (lo['air_c'] - celsius) * (hi['air_c'] - celsius) <= 0 and lo['air_c'] != hi['air_c']:
            f = (lo['air_c'] - celsius) / (lo['air_c'] - hi['air_c'])
            return lo['height_km'] + f * (hi['height_km'] - lo['height_km'])
    return None


# Shares of time past a threshold do not interpolate: they follow the tail of the distribution.
RUNS_ONLY = ('wet_bulb_above_28', 'wet_bulb_above_31', 'wet_bulb_above_28_tropics',
             'wet_bulb_2m_above_28', 'wet_bulb_2m_above_31', 'wet_bulb_2m_above_28_tropics')
# The land's heat in the day's and the night's 3-day means: the climate's keys and where the extract holds them.
DAY_NIGHT = {f'{part}_{name}': (part, height, var, q)
             for part in ('day', 'night')
             for height, tag in (('at_2_m', '_2m'), ('lowest_level', ''))
             for name, var, q in ((f'wet_bulb{tag}_median_c', 'wet_bulb', 'median'),
                                  (f'wet_bulb{tag}_p99_c', 'wet_bulb', 'p99'),
                                  (f'land_air{tag}_median_c', 'air', 'median'),
                                  (f'land_air{tag}_p99_c', 'air', 'p99'))}
# The 99th percentiles the 6.5-hour sensitivity carries, with their medians and where the smoothing holds them.
SMOOTHED = {'wet_bulb_2m': ('wet_bulb_2m_median_c', 'wet_bulb_2m_p99_c', 'at_2_m', 'wet_bulb'),
            'land_air_2m': ('land_air_2m_median_c', 'land_air_2m_p99_c', 'at_2_m', 'air'),
            'wet_bulb': ('wet_bulb_median_c', 'wet_bulb_p99_c', 'lowest_level', 'wet_bulb'),
            'land_air': ('land_air_median_c', 'land_air_p99_c', 'lowest_level', 'air')}


def heights(ladder):
    return dict(height_of_20c_km=crossing(ladder, 20.0), freezing_level_km=crossing(ladder, 0.0))


def climate(run):
    """The figures the targets interpolate, from one run's years 10-19. A crossing height is None where the lowest
    level is already colder."""
    h = run['heat']
    g = h['land']['at_2_m']
    every, model = run['land']['every_cell']['all'], run['land']['model_land']['all']
    return dict(ground_k=run['ground_k'], rain_mm_day=run['rain_mm_day'], land_rain_mm_day=every['rain_mm_day'],
                land_desert=every['desert'], model_land_rain_mm_day=model['rain_mm_day'],
                model_land_desert=model['desert'],
                wet_bulb_median_c=h['land']['wet_bulb']['median'], wet_bulb_p99_c=h['land']['wet_bulb']['p99'],
                wet_bulb_above_28=h['land']['above_28'], wet_bulb_above_31=h['land']['above_31'],
                wet_bulb_above_28_tropics=h['bands']['0-15']['above_28'],
                land_air_median_c=h['land']['air']['median'], land_air_p99_c=h['land']['air']['p99'],
                wet_bulb_2m_median_c=g['wet_bulb']['median'], wet_bulb_2m_p99_c=g['wet_bulb']['p99'],
                wet_bulb_2m_above_28=g['above_28'], wet_bulb_2m_above_31=g['above_31'],
                wet_bulb_2m_above_28_tropics=h['bands']['0-15']['at_2_m']['above_28'],
                land_air_2m_median_c=g['air']['median'], land_air_2m_p99_c=g['air']['p99'],
                sea_median_c=run['seas']['all']['median'],
                **{k: h['land'][part][height][var][q] for k, (part, height, var, q) in DAY_NIGHT.items()},
                ladder=[dict(height_km=x['height_km'], air_c=x['air_c'], pressure_atm=x['pressure_atm'])
                        for x in run['sky_over_land']],
                **heights(run['sky_over_land']))


def comfort(run):
    """A run's own figures beside its climate: its heat place by place at both heights, its seas, its land's rain and
    deserts by band and height, the land the model runs on and the median of the land's column crossings of
    freezing."""
    return dict(heat=run['heat'], seas=run['seas'], land=run['land'], land_and_water=run['land_and_water'],
                freezing_level_over_land_km=run['freezing_level_over_land_km'],
                land_columns_that_never_freeze=run['land_columns_that_never_freeze'])


def interpolate(a, b, ground_k):
    """A target's climate between the two runs, linear in the ground's temperature, the sky ladder height by height;
    a target beyond either run takes that run's own climate."""
    f = min(max((ground_k - a['ground_k']) / (b['ground_k'] - a['ground_k']), 0.0), 1.0)
    out = {k: a[k] + f * (b[k] - a[k]) for k in a
           if k not in RUNS_ONLY and k not in ('ladder', 'height_of_20c_km', 'freezing_level_km')}
    out['ladder'] = [dict(height_km=x['height_km'], air_c=x['air_c'] + f * (y['air_c'] - x['air_c']),
                          pressure_atm=x['pressure_atm'] + f * (y['pressure_atm'] - x['pressure_atm']))
                     for x, y in zip(a['ladder'], b['ladder'])]
    out.update(heights(out['ladder']))
    out['ground_k'] = ground_k
    return out, f


def corrected_efficacy(ex, gases_glow, infrared, gas_forcing):
    """The glow's efficacy from the glow pair on the corrected radiation of 17 October. The two branches differ in
    their dimming alone, so their grounds give the model's response to it, which carries the held branch to the
    design case's temperature: the gases and the glow together hold it with that much of the design case's absorbed
    sunlight, and the pair fixes only that total. The gases take their forcing, as the heat beyond the sunlight books
    them and as 13 October's 1.69% was their forcing over that radiation's absorbed sunlight, and the glow takes the
    rest. The other split keeps the gases at their 1.69% of this radiation's absorbed sunlight, which is larger, the
    curved air's light among it; it gives the same total."""
    pair = ex['glow_pair']
    design, more, less = pair[DESIGN], pair[WW], pair[HELD]

    def dimming(r):
        return 100.0 * (1.0 - r['sunlight_scale'] / design['sunlight_scale'])

    response = (less['ground_k'] - more['ground_k']) / (dimming(more) - dimming(less))    # K for each 1% of dimming
    hold = dimming(less) - (design['ground_k'] - less['ground_k']) / response
    absorbed = ex['runs'][DESIGN]['top_of_air_w_m2']['absorbed']
    total = hold / 100.0 * absorbed
    points = gases_glow['runs']['gases']['added_dimming_percent']
    at_points = points / 100.0 * absorbed
    return dict(
        added_dimming_percent={WW: dimming(more), HELD: dimming(less)},
        ground_k={DESIGN: design['ground_k'], WW: more['ground_k'], HELD: less['ground_k']},
        ground_trend_k_per_decade={WW: more['ground_trend_k_per_decade'], HELD: less['ground_trend_k_per_decade']},
        response_k_per_percent=response, dimming_that_holds_percent=hold, design_absorbed_sunlight_w_m2=absorbed,
        gases_and_glow_w_m2=total, gases_w_m2=gas_forcing, glow_w_m2=total - gas_forcing, infrared_w_m2=infrared,
        efficacy=(total - gas_forcing) / infrared,
        with_the_gases_at_their_points=dict(gases_dimming_percent=points, gases_w_m2=at_points,
                                            glow_w_m2=total - at_points, efficacy=(total - at_points) / infrared))


def wet_bulb_per_kelvin_of_air(run):
    """How far a kelvin of air moves the wet-bulb temperature at a fixed humidity: at the land's day medians at 2 m,
    the humidity that gives that wet-bulb temperature at that air, over the land's 5th to 95th percentile surface
    pressures."""
    land = run['heat']['land']
    t = land['day']['at_2_m']['air']['median'] + 273.15
    tw = land['day']['at_2_m']['wet_bulb']['median'] + 273.15
    lv = 2.501e6 - 2370.0 * (tw - 273.15)
    pressures = [land['surface_pressure_bar']['p05'], land['surface_pressure_bar']['p95']]
    out = []
    for bar in pressures:
        p = np.array([bar * 1e5])
        # the psychrometric balance solved for the humidity, then the wet-bulb temperature a kelvin either side
        q = (X.qsat(np.array([tw]), p) - X.CP * (t - tw) / lv) / (1.0 + 1850.0 * (t - tw) / lv)
        up, down = (float(X.wet_bulb(np.array([t + d]), q, p)[0]) for d in (0.5, -0.5))
        out.append(dict(pressure_bar=bar, specific_humidity=float(q[0]), k_per_k=up - down))
    return dict(air_c=t - 273.15, wet_bulb_c=tw - 273.15, at=out,
                k_per_k=sorted(x['k_per_k'] for x in out))


def heat_beyond_sunlight(data, absorbed, area, corrected):
    """The heat the plan adds on the Moon besides its sunlight, each as the absorbed sunlight that answers it."""
    gases_glow, glow = data[GASES_GLOW], data[GLOW]['results']['recommended']
    sources = dict(held_branch_13_october=gases_glow['glow']['held']['efficacy'],
                   corrected_radiation_17_october=corrected)
    efficacy = sorted(sources.values())
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
        glow=item([e * infrared + (1.0 - albedo) * sunlight for e in efficacy], answer=glow['answer'],
                  infrared_w_m2=infrared, sunlight_back_w_m2=sunlight, efficacy=efficacy, efficacy_sources=sources,
                  basis='the recommended answer, both halves rolled 80 degrees on the re-flown layout, at the GCM\'s '
                        'efficacy for its infrared, from the glow pair on the corrected radiation of 17 October with '
                        'the gases at their forcing to the held branch of 13 October, and the albedo for its '
                        'returned sunlight'),
        beamed_power=item([tw * 1e12 / area for tw in power_tw], heat_tw=power_tw,
                          basis='the power screen\'s intake, the beams arriving at the receivers and the towns\' '
                                'crowns, and the deep air\'s share of the beams; all of it ends as heat at or near the '
                                'ground'))
    # What the air must send out besides the sunlight it absorbs, at full watts: the gases' forcing, all of the glow's
    # infrared, the returned sunlight the Moon absorbs and the beamed power's heat. The GCM's efficacy weighs the glow's
    # infrared for the ground's temperature; a column's ceiling on its outgoing infrared counts all of it.
    out['to_send_out_w_m2'] = [gas + infrared + (1.0 - albedo) * sunlight + tw * 1e12 / area for tw in power_tw]
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


def dimming_curve(fits, absorbed, every_dim):
    """The dimming for each ground temperature from the design temperature to where the hot world settles, with the
    sunlight scale of a design case for it, before and after the heat beyond the sunlight is answered too."""
    rows = []
    for target in CURVE_K:
        s = shed(fits, target)
        dim = [100.0 * x / absorbed for x in s]
        total = [d + e for d, e in zip(dim, every_dim)]
        rows.append(dict(ground_k=target, absorbed_sunlight_to_shed_w_m2=s, dimming_percent=dim,
                         with_the_heat_beyond_sunlight_percent=total,
                         sunlight_scale_of_its_design_case=[DESIGN_SCALE * (1.0 - d / 100.0) for d in dim][::-1],
                         sunlight_scale=[DESIGN_SCALE * (1.0 - t / 100.0) for t in total][::-1]))
    return dict(per_kelvin_percent=sorted(100.0 * f['feedback_w_m2_k'] / absorbed for f in fits), targets=rows)


def smoothing(eddies):
    """How much the 3-day means smooth the land's heat in the design case: each extreme in the 6.5-hour means against
    the same in the 3-day means their elevens make over the same two years, by its excess and by the ratio of its
    height above the median."""
    out = dict(run=eddies['run'], years=eddies['years'], start=eddies['start'], hours_each=eddies['hours_each'],
               means_to_a_3_day_mean=eddies['means_to_a_3_day_mean'], means=eddies['means'],
               means_kept=eddies['means_kept'], three_day_means_kept=eddies['three_day_means_kept'],
               supersaturated_lowest_level={k: eddies[k]['lowest_level']['supersaturated']
                                            for k in ('at_6_5_hours', 'in_3_day_means')})
    for height in ('at_2_m', 'lowest_level'):
        out[height] = {}
        for var in ('wet_bulb', 'air'):
            fine, coarse = eddies['at_6_5_hours'][height][var], eddies['in_3_day_means'][height][var]
            out[height][var] = dict(median=dict(at_6_5_hours=fine['median'], in_3_day_means=coarse['median']),
                                    **{q: dict(at_6_5_hours=fine[q], in_3_day_means=coarse[q],
                                               excess_k=fine[q] - coarse[q],
                                               ratio_above_the_median=(fine[q] - fine['median'])
                                               / (coarse[q] - coarse['median']))
                                       for q in ('p99', 'top')})
    return out


def at_6_5_hours(state, smooth):
    """A target's 99th percentiles in 6.5-hour means, as a sensitivity: its own median plus its 99th percentile's
    height above it, stretched by the design case's ratio."""
    return {f'{k}_p99_c': state[m] + smooth[height][var]['p99']['ratio_above_the_median'] * (state[p] - state[m])
            for k, (m, p, height, var) in SMOOTHED.items()}


def lake_oxygen():
    """Ecology's waters' oxygen, in the tree, and its sha256."""
    raw = (ROOT / LAKE_OXYGEN['path']).read_bytes()
    data = json.loads(raw)
    if data.get('schema') != LAKE_OXYGEN['schema']:
        raise ValueError('the waters\' oxygen product has another schema')
    return data, hashlib.sha256(raw).hexdigest()


def crossing_m(ladder, values, celsius):
    """The first height (m) from below at which the air falls through a temperature, linear between the ladder's
    heights; None where it stays warmer to the ladder's top."""
    v = np.asarray(values) - celsius
    for i in range(len(v) - 1):
        if v[i] >= 0.0 > v[i + 1]:
            return float(ladder[i] + v[i] / (v[i] - v[i + 1]) * (ladder[i + 1] - ladder[i]))
    return None


def tower(ex, data):
    """The summit tower's air for every target: the mean, warm and cold tenths of the 3-day means over the summit,
    interpolated between the design case and the test by the ground's temperature as the targets' climates are; the
    heights where the mean air reaches 0 and -5 C, the tower's length in sub-freezing air, the air at its decks and
    top; the warmest target at which its top keeps a winter, bare and with margins; and the checks."""
    s = ex['summit']
    L = np.asarray(s['ladder_m'])
    site = data[SUMMIT_TOWER]['site']
    if site['ground_above_sea_m'] != s['site']['ground_above_sea_m']:
        raise ValueError('the extract read another summit')
    z0 = float(site['ground_above_sea_m'])
    height = data[SUMMIT_TOWER]['port']['height_km'] * 1000.0
    top = z0 + height
    form = data[FORM]
    disks = [d for d in form['disks']]
    decks = dict(open_disk_m=next(d['z_m'] for d in disks if not d['sealed'] and not d['open_winter']),
                 open_winter_disk_m=next(d['z_m'] for d in disks if d['open_winter']),
                 first_sealed_disk_m=next(d['z_m'] for d in disks if d['sealed']),
                 highest_disk_m=max(d['z_m'] for d in disks), terminal_m=form['terminal']['z_m'], top_m=height)
    runs = s['runs']
    d, t = runs[DESIGN]['air_c']['bilinear'], runs[TEST]['air_c']['bilinear']
    dp, tp = np.asarray(runs[DESIGN]['pressure_atm']), np.asarray(runs[TEST]['pressure_atm'])
    gd, gt = ex['runs'][DESIGN]['ground_k'], ex['runs'][TEST]['ground_k']
    at = lambda values, z: float(np.interp(z, L, values))
    eddies = s['eddies']
    fine, coarse = eddies['at_6_5_hours']['air_c']['bilinear'], eddies['in_3_day_means']['air_c']['bilinear']
    excess = {q: dict(at_6_5_hours=at(fine[q], top) - at(fine['mean'], top),
                      in_3_day_means=at(coarse[q], top) - at(coarse['mean'], top))
              for q in ('warm_tenth', 'cold_tenth', 'high', 'low')}

    def column(f, q):
        return [a + f * (b - a) for a, b in zip(d[q], t[q])]

    def above(z):
        return None if z is None else (z - z0) / 1000.0

    rows = []
    for target in CURVE_K:
        f = min(max((target - gd) / (gt - gd), 0.0), 1.0)
        mean, warm, cold = (column(f, q) for q in ('mean', 'warm_tenth', 'cold_tenth'))
        z_0, z_5 = crossing_m(L, mean, 0.0), crossing_m(L, mean, -5.0)
        z_warm, z_cold = crossing_m(L, warm, 0.0), crossing_m(L, cold, 0.0)
        pressure = dp + f * (tp - dp)
        rows.append(dict(
            ground_k=target, share_of_the_way=f, freezing_above_summit_km=above(z_0),
            minus_5_above_summit_km=above(z_5), freezing_above_summit_in_the_tenths_km=dict(
                warm=above(z_warm), cold=above(z_cold)),
            tower_below_freezing_km=0.0 if z_0 is None else min(max((top - z_0) / 1000.0, 0.0), height / 1000.0),
            tower_below_freezing_in_the_warm_tenth_km=0.0 if z_warm is None else min(max((top - z_warm) / 1000.0, 0.0),
                                                                                       height / 1000.0),
            pressure_at_freezing_atm=None if z_0 is None else float(np.interp(z_0, L, pressure)),
            top_c=dict(mean=at(mean, top), warm_tenth=at(warm, top), cold_tenth=at(cold, top),
                       warm_tenth_at_6_5_hours=at(mean, top) + excess['warm_tenth']['at_6_5_hours'],
                       highest_6_5_hour_mean=at(mean, top) + excess['high']['at_6_5_hours']),
            summit_c=at(mean, z0), decks_c={k: at(mean, z0 + v) for k, v in decks.items()}))

    def warmest(z, q='mean', extra=0.0):
        """The target at which the interpolated air at a height reaches 0 C, its tenth's or with an excess added."""
        a, b = at(d[q], z) + extra, at(t[q], z) + extra
        return gd + (-a / (b - a)) * (gt - gd)

    margins = {f'freezing_{m / 1000:g}_km_below_the_top': warmest(top - m) for m in MARGINS_M}
    warmest_targets = dict(top_mean_at_0c=warmest(top), top_warm_tenth_at_0c=warmest(top, 'warm_tenth'),
                           top_warm_tenth_at_6_5_hours_at_0c=warmest(top, extra=excess['warm_tenth']['at_6_5_hours']),
                           top_highest_6_5_hour_mean_at_0c=warmest(top, extra=excess['high']['at_6_5_hours']),
                           **margins)
    lo, hi = (next(r for r in rows if r['ground_k'] == k) for k in (295.0, 303.0))
    rise = (hi['freezing_above_summit_km'] - lo['freezing_above_summit_km']) / (hi['ground_k'] - lo['ground_k'])
    # the brief's air: the earlier design case read at the site, its ground from the sea-level pressure taken as
    # isothermal (climate/gcm/site_winds.py)
    sw = data[SITE_WINDS]
    zb = np.asarray(sw['height_m']) + sw['ground_height_m']
    tb = np.asarray(sw['temperature_k']) - 273.15
    brief_air = lambda z: float(np.interp(z, zb, tb))
    zb0 = None
    for i in range(len(tb) - 1):
        if tb[i] >= 0.0 > tb[i + 1]:
            zb0 = float(zb[i] + tb[i] / (tb[i] - tb[i + 1]) * (zb[i + 1] - zb[i]))
            break

    def reading(air, freezing):
        return dict(summit_c=air(z0), disks_from_c=air(z0 + decks['open_disk_m']),
                    winter_parks_c=[air(z0 + 1000.0 * k) for k in BRIEF['winter_parks_km']], top_c=air(top),
                    freezing_above_summit_km=above(freezing), lapse_k_per_km=(air(z0) - air(top)) / (height / 1000.0))

    mine = {run: reading(lambda z, r=run: at(runs[r]['air_c']['bilinear']['mean'], z),
                         crossing_m(L, runs[run]['air_c']['bilinear']['mean'], 0.0)) for run in (DESIGN, OLD_DESIGN)}
    brief = dict(stated=BRIEF, site_product=dict(run=sw['run'], years=sw['years'], ground_height_m=sw['ground_height_m'],
                                                 **reading(brief_air, zb0)),
                 model_ground_under_the_weights_m=runs[OLD_DESIGN]['ground_under_m']['bilinear'],
                 earlier_design_case=mine[OLD_DESIGN], design_case=mine[DESIGN])
    # the test's own warming, year by year: the top's air against the ground, a line over two windows of years
    ground = {y['year']: y['ground_k'] for y in ex['test_progress']}
    years = [(y['year'], ground[y['year']], at(y['air_c'], top)) for y in s['test_years']]
    lines = {}
    for lo, hi in FITS:
        sel = [(g, a) for y, g, a in years if lo <= y <= hi]
        slope, icpt = np.polyfit([g for g, _ in sel], [a for _, a in sel], 1)
        lines[f'{lo}-{hi}'] = dict(k_per_k=float(slope), **{
            f'at_{k:g}_k': dict(line_c=float(slope * k + icpt),
                                interpolated_c=next(r['top_c']['mean'] for r in rows if r['ground_k'] == k))
            for k in (295.0, 300.0, 302.0)})
    weighting = {run: {k: dict(summit_c=at(v['mean'], z0) - at(runs[run]['air_c']['bilinear']['mean'], z0),
                               top_c=at(v['mean'], top) - at(runs[run]['air_c']['bilinear']['mean'], top))
                       for k, v in runs[run]['air_c'].items() if k != 'bilinear'} for run in (DESIGN, TEST)}
    return dict(
        site=dict(**s['site'], holding_cell=s['holding_cell'], tower_height_m=height, top_above_sea_m=top,
                  decks_above_summit_m=decks),
        orography=dict(holding_cell_land_fraction=s['block']['boundary_land_fraction'][2][2],
                       holding_cell_land_elevation_m=s['block']['boundary_land_elevation_m'][2][2],
                       model_ground_m={run: dict(holding_cell=r['ground_m'][2][2], **r['ground_under_m'],
                                                 spread_through_the_means_m=r['ground_spread_m'])
                                       for run, r in runs.items()},
                       design_levels=runs[DESIGN]['levels'], test_levels=runs[TEST]['levels']),
        runs={run: dict(ground_k=ex['runs'][run]['ground_k'],
                        **{q: dict(summit_c=at(r['air_c']['bilinear'][q], z0), top_c=at(r['air_c']['bilinear'][q], top))
                           for q in ('mean', 'warm_tenth', 'cold_tenth', 'p10', 'p90', 'low', 'high')})
              for run, r in ((DESIGN, runs[DESIGN]), (TEST, runs[TEST]))},
        targets=rows, warmest_targets=warmest_targets, freezing_rise_km_per_k_295_to_303=rise,
        top_rise_k_per_k=(at(t['mean'], top) - at(d['mean'], top)) / (gt - gd), smoothing_at_the_top_k=excess, brief=brief,
        test_years=dict(years=[dict(year=y, ground_k=g, top_c=a) for y, g, a in years], lines=lines),
        weighting_against_bilinear_k=weighting)


def resort_band(data, ex, fits, absorbed, incoming, area, clear, saturated, heat, light, rain_cut, conversion, a, b,
                smooth, towers):
    """The resort's candidates from 298 to 303 K against winter on the tower, comfort at the ground, the dimming, rain
    and deserts, the waters' oxygen and the margin under a saturated column. The oxygen interpolates ecology's three
    targets: saturation by the seas' median warmth, the lakes by the ground's temperature."""
    lakes, lakes_sha256 = lake_oxygen()
    sat = lakes['saturation']['targets']
    gd, gt = ex['runs'][DESIGN]['ground_k'], ex['runs'][TEST]['ground_k']
    points = dict(lock=gd, resort=300.0, hot=gt)

    def lake(kind, column, key, part=None):
        out = {}
        for name in points:
            vals = []
            for band in ('0-15', '15-30', '30-45'):
                x = lakes['runs'][kind][f'{name}/{band}']['oxygen'][column]
                vals.append((x[part] if part else x)[key])
            out[name] = [min(vals), max(vals)]
        return out

    deep = lake('deep', 'oligotrophic/renewal/supply', 'mean', 'deepest_tenth')
    shallow = lake('shallow', 'eutrophic/renewal/supply', 'hypoxic_share', 'nights')
    by_ground = lambda table, k: [float(np.interp(k, list(points.values()), [table[n][e] for n in points]))
                                  for e in (0, 1)]
    seas = sorted((sat[n]['sea_median_c'], sat[n]['mg_l']) for n in points)
    every = heat['answered_in_every_target']['absorbed_sunlight_w_m2']
    gas_w = heat['gases']['absorbed_sunlight_w_m2'][0]
    send = heat['to_send_out_w_m2']
    rows = []
    for target in BAND_K:
        s = shed(fits, target)
        state, frac = interpolate(a, b, target)
        dim = [100.0 * x / absorbed for x in s]
        total = [100.0 * (x + e) / absorbed for x, e in zip(s, every)]
        taken = [x / 100.0 * incoming * area / 1e12 for x in dim]
        rain = state['land_rain_mm_day']
        answered = [clear * (1.0 - (x + 100.0 * (gas_w + g + p) / absorbed) / 100.0) + out
                    for x in dim for g in heat['glow']['absorbed_sunlight_w_m2']
                    for p, out in zip(heat['beamed_power']['absorbed_sunlight_w_m2'], send)]
        t = next(r for r in towers['targets'] if r['ground_k'] == target)
        rows.append(dict(
            ground_k=target, share_of_the_way=frac,
            tower=dict(freezing_above_summit_km=t['freezing_above_summit_km'],
                       tower_below_freezing_km=t['tower_below_freezing_km'], top_c=t['top_c']['mean'],
                       top_warm_tenth_c=t['top_c']['warm_tenth']),
            comfort=dict(wet_bulb_2m_median_c=state['wet_bulb_2m_median_c'], wet_bulb_2m_p99_c=state['wet_bulb_2m_p99_c'],
                         day_wet_bulb_2m_p99_c=state['day_wet_bulb_2m_p99_c'],
                         wet_bulb_2m_p99_at_6_5_hours_c=at_6_5_hours(state, smooth)['wet_bulb_2m_p99_c'],
                         land_air_2m_median_c=state['land_air_2m_median_c'], sea_median_c=state['sea_median_c']),
            dimming=dict(own_percent=dim, with_the_heat_beyond_sunlight_percent=total,
                         sunlight_scale=[DESIGN_SCALE * (1.0 - x / 100.0) for x in total][::-1],
                         sunlight_taken_tw=taken, electricity_tw=[taken[0] * conversion[0], taken[1] * conversion[1]],
                         growth_cost_through_light_percent=[light[0] * dim[0], light[1] * dim[1]]),
            rain=dict(land_rain_mm_day=rain, with_the_dimmings_own_cut_mm_day=[
                          rain * (1.0 - rain_cut[1] / 100.0 * dim[1]), rain * (1.0 - rain_cut[0] / 100.0 * dim[0])],
                      land_desert=state['land_desert'], model_land_desert=state['model_land_desert']),
            oxygen=dict(saturation_mg_l=float(np.interp(state['sea_median_c'], [x for x, _ in seas], [y for _, y in seas])),
                        deep_oligotrophic_deepest_tenth_mg_l=by_ground(deep, target),
                        shallow_eutrophic_hypoxic_nights=by_ground(shallow, target)),
            gap_with_the_added_heat_answered_w_m2=[saturated - max(answered), saturated - min(answered)]))
    return dict(candidates_k=list(BAND_K), rows=rows,
                oxygen_source=dict(path=LAKE_OXYGEN['path'], sha256=lakes_sha256,
                                   targets_ground_k=points, deep=deep, shallow=shallow,
                                   saturation={n: dict(sea_median_c=sat[n]['sea_median_c'], mg_l=sat[n]['mg_l'])
                                               for n in points}))


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
    clear = toa['absorbed_clear_sky']
    saturated = data[COLUMN]['moon_1.2atm']['saturated']['max_olr']
    conversion = data[GLOW]['design']['conversion']
    light = data[LIGHT]['cost_through_light']['per_percent_dimming']
    rain_cut = [-data[GASES_GLOW]['runs']['gases']['rain_per_percent_dimming']['land'],
                data[LIGHT]['rain_pathway']['precipitation_cut_per_percent_dimming'][1]]
    recommended = data[GLOW]['results']['recommended']
    infrared = recommended['glow_W_m2']['global'] + recommended['infrared_return_W_m2']['global']
    gas = data[GAS_FORCING]['moon_at_design']['both'][HELD_PAIR]['allsky_w_m2']
    efficacy = corrected_efficacy(ex, data[GASES_GLOW], infrared, gas)
    heat = heat_beyond_sunlight(data, absorbed, area, efficacy['efficacy'])
    every = heat['answered_in_every_target']['absorbed_sunlight_w_m2']
    every_dim = heat['answered_in_every_target']['dimming_percent']
    send = heat['to_send_out_w_m2']
    shared = dict(dimming_percent=every_dim,
                  growth_cost_through_light_percent=[light[0] * every_dim[0], light[1] * every_dim[1]],
                  land_rain_cut_percent=[rain_cut[0] * every_dim[0], rain_cut[1] * every_dim[1]])
    smooth = smoothing(ex['eddies'])
    a, b = climate(design), climate(test)
    targets = []
    for name, target in [('no further dimming', None)] + [(f'{t:g} K', t) for t in TARGETS_K]:
        if target is None:
            s = [0.0, 0.0]
            state, frac = dict(b), 1.0
        else:
            s = shed(fits, target)
            state, frac = interpolate(a, b, target)
        dim = [100.0 * x / absorbed for x in s]
        total = [100.0 * (x + e) / absorbed for x, e in zip(s, every)]
        taken_tw = [d / 100.0 * incoming * area / 1e12 for d in dim]
        rain = state['land_rain_mm_day']
        # The clear sky's load against a saturated column's ceiling: the clear-sky sunlight left after the dimming,
        # with the added heat at full watts, answered by its own dimming or not. Its ends are the least and most load
        # over the ends of the feedback, the glow's efficacy and the beamed power, which also sets the watts to send.
        gas_w = heat['gases']['absorbed_sunlight_w_m2'][0]
        answered = [clear * (1.0 - (d + 100.0 * (gas_w + g + p) / absorbed) / 100.0) + out
                    for d in dim for g in heat['glow']['absorbed_sunlight_w_m2']
                    for p, out in zip(heat['beamed_power']['absorbed_sunlight_w_m2'], send)]
        unanswered = [clear * (1.0 - d / 100.0) + out for d in dim for out in send]
        targets.append(dict(
            name=name, ground_k=state['ground_k'] if target is not None else test['ground_k'],
            settles_k=bal if target is None else [target, target],
            absorbed_sunlight_to_shed_w_m2=s, dimming_percent=dim,
            with_the_heat_beyond_sunlight_percent=total,
            sunlight_scale_of_its_design_case=[DESIGN_SCALE * (1.0 - d / 100.0) for d in dim][::-1],
            sunlight_scale=[DESIGN_SCALE * (1.0 - t / 100.0) for t in total][::-1],
            sunlight_taken_tw=taken_tw, electricity_tw=[taken_tw[0] * conversion[0], taken_tw[1] * conversion[1]],
            growth_cost_through_light_percent=[light[0] * dim[0], light[1] * dim[1]],
            land_rain_mm_day=rain, land_rain_with_the_dimmings_own_cut_mm_day=[
                rain * (1.0 - rain_cut[1] / 100.0 * dim[1]), rain * (1.0 - rain_cut[0] / 100.0 * dim[0])],
            climate=state, interpolated_share_of_the_way=frac, at_6_5_hours=at_6_5_hours(state, smooth),
            clear_sky_absorbed_w_m2=[clear * (1.0 - d / 100.0) for d in dim][::-1],
            gap_to_the_saturated_ceiling_w_m2=sorted(saturated - clear * (1.0 - d / 100.0) for d in dim),
            gap_with_the_added_heat_answered_w_m2=[saturated - max(answered), saturated - min(answered)],
            gap_with_the_added_heat_unanswered_w_m2=[saturated - max(unanswered), saturated - min(unanswered)]))
    towers = tower(ex, data)
    band = resort_band(data, ex, fits, absorbed, incoming, area, clear, saturated, heat, light, rain_cut, conversion,
                       a, b, smooth, towers)
    return dict(
        schema=SCHEMA,
        producer=dict(files={c: digest(c) for c in CODE},
                      inputs={**{p: digest(p) for p in INPUTS}, LAKE_OXYGEN['path']: digest(LAKE_OXYGEN['path'])},
                      constants=dict(MOON_RADIUS=K.MOON_RADIUS, JULIAN_YEAR_DAYS=K.JULIAN_YEAR_DAYS,
                                     JULIAN_DAY=K.JULIAN_DAY)),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        # each run's heat place by place at both heights, its seas by band and the column freezing level, which the
        # climate README reads; they are the runs' own and do not interpolate
        runs=dict(design=dict(run=DESIGN, **climate(design), **comfort(design),
                              top_of_air_w_m2=design['top_of_air_w_m2']),
                  test=dict(run=TEST, **climate(test), **comfort(test), top_of_air_w_m2=toa, planetary_albedo=albedo)),
        settled=dict(fits=fits, feedback_w_m2_k=lam, balance_k=bal,
                     saturated_column_max_outgoing_w_m2=saturated,
                     clear_sky_outgoing_over_the_saturated_column_w_m2=toa['outgoing_infrared_clear_sky'] - saturated,
                     clear_sky_absorbed_w_m2=clear,
                     cloud_effect_w_m2=dict(sunlight=absorbed - clear,
                                            infrared=toa['outgoing_infrared_clear_sky'] - toa['outgoing_infrared'])),
        glow_efficacy=efficacy,
        shared_dimming=shared,
        heat_beyond_sunlight=heat,
        dimming_curve=dimming_curve(fits, absorbed, every_dim),
        smoothing=smooth,
        wet_bulb_per_kelvin_of_air=dict(design=wet_bulb_per_kelvin_of_air(design),
                                        test=wet_bulb_per_kelvin_of_air(test)),
        targets=targets,
        costs_per_percent=dict(growth_through_light_percent=light, land_rain_percent=rain_cut),
        tower=towers, resort_band=band)


def main(argv=None) -> int:
    product = results()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1) + '\n')
    g = product['glow_efficacy']
    print(f"glow efficacy on the corrected radiation {g['efficacy']:.3f} (response {g['response_k_per_percent']:.2f} "
          f"K a point); shared dimming {product['shared_dimming']['dimming_percent'][0]:.2f}-"
          f"{product['shared_dimming']['dimming_percent'][1]:.2f}%; dimming per kelvin "
          f"{product['dimming_curve']['per_kelvin_percent'][0]:.3f}-"
          f"{product['dimming_curve']['per_kelvin_percent'][1]:.3f}%")
    for t in product['targets']:
        print(f"{t['name']:20s} shed {t['absorbed_sunlight_to_shed_w_m2'][0]:5.1f}-{t['absorbed_sunlight_to_shed_w_m2'][1]:5.1f} W/m2, "
              f"dimming {t['dimming_percent'][0]:4.1f}-{t['dimming_percent'][1]:4.1f}% "
              f"(+heat {t['with_the_heat_beyond_sunlight_percent'][0]:4.1f}-{t['with_the_heat_beyond_sunlight_percent'][1]:4.1f}%), "
              f"scale {t['sunlight_scale'][0]:.3f}-{t['sunlight_scale'][1]:.3f}, "
              f"wet-bulb {t['climate']['wet_bulb_median_c']:.1f}/{t['climate']['wet_bulb_p99_c']:.1f} C, "
              f"desert {100*t['climate']['land_desert']:.0f}%, sea {t['climate']['sea_median_c']:.1f} C, "
              f"clear gap {t['gap_to_the_saturated_ceiling_w_m2'][0]:.1f}-{t['gap_to_the_saturated_ceiling_w_m2'][1]:.1f}, "
              f"answered {t['gap_with_the_added_heat_answered_w_m2'][0]:.1f}-{t['gap_with_the_added_heat_answered_w_m2'][1]:.1f}, "
              f"unanswered {t['gap_with_the_added_heat_unanswered_w_m2'][0]:.1f}-{t['gap_with_the_added_heat_unanswered_w_m2'][1]:.1f} W/m2")
    print('wrote', OUT.relative_to(ROOT))
    return 0


if __name__ == '__main__':
    sys.exit(main())
