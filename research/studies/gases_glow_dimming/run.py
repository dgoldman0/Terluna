"""The gases, the glow and the dimming in the GCM: whether the dimming holds the design climate, and what the rain does.

    climate/gcm/.venv/bin/python -m research.studies.gases_glow_dimming.extract   # once, from the GCM's output
    python -m research.studies.gases_glow_dimming.run                             # results/gases_glow_dimming.json

On 13 October 2026 the author asked for the GCM to run again before the gas pair and the dimming's cost are decided:
in Earth's models the rain falls when dimming answers a greenhouse warming (GeoMIP G1), and most of the Moon's land
is water-limited (the light limitation study). Three runs branch from the design case's year-29 restart (climate/gcm,
run A28_dim5_moon) and fly 20 years each on the same model:

- gases: the middle pair's warming, 3.58 W/m2 all-sky in the deep column, carried as the CO2 that warms PlaSim's air
  as much, since PlaSim's longwave radiation has water vapour, CO2 and ozone only. PlaSim's own forcing for 400 to
  1,022 ppm, from a calibration branch, is 5.27 W/m2, so the pair is 755 ppm. The sunlight is dimmed by the 1.69%
  that answers it;
- gases and glow: the same, with the fleet's glow and the Moon's own infrared its films return (the recommended answer
  on the re-flown layout, research/array-industry) entering as infrared from space at the top of the air, 10.17 W/m2
  over the day hemisphere and 0.52 over the night, and the sunlight dimmed by the 2.53% more that matches that
  infrared watt for watt in the energy balance. The sunlight the films return and the dimming that answers it cancel,
  and both are left out;
- gases and glow, held: the same infrared with the glow's dimming cut to 1.83%, the share that the second run's
  cooling and the model's response to sunlight (about 1.9 K for each 1%, climate/gcm README) put at the design
  temperature.

The control is the design case's years 20-29; each run's years 10-19 are compared with it, as climatologies
(climate/gcm/climatology.py) and as layer means and first 3-day means (extract.py). The two glow runs give the
model's response to the glow's dimming, and with it the dimming that holds the design temperature and the rain there.
The rain and the evaporation over land, band by band, carry the change into the rain-fed growth of the light
limitation study, whose water-limited plants grow with the water they transpire.
"""
from __future__ import annotations
import hashlib
import json
import math
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'gases_glow_dimming.json'
SCHEMA = 'terluna.research.gases-glow-dimming/1'
INPUTS = {
    'research/studies/light_limitation/results/light_limitation.json': 'terluna.research.light-limitation/1',
    'research/studies/domain_synthesis/results/ledgers.json': 'terluna.research.domain-synthesis-ledgers/1',
    'research/studies/gases_glow_dimming/results/glow_absorption.json': 'terluna.research.glow-absorption/1',
}
# The GCM's products, runs and the extract, generated outside Git and bound by their hashes.
CLIMATOLOGY = 'climate/gcm/products/climatology_{}.npz'
CLIMATOLOGY_SCHEMA = 'terluna.climate.gcm-climatology/1'
PROGRESS = 'climate/gcm/runs/{}/progress.json'
EXTRACT = 'research/runs/gases_glow_dimming/extract.json'
EXTRACT_SCHEMA = 'terluna.research.gases-glow-dimming-extract/1'
CONTROL = ('A28_dim5_moon', (20, 29))
RUNS = {'gases': ('A28_gases', (10, 19)), 'gases_and_glow': ('A28_gases_glow', (10, 19)),
        'gases_and_glow_held': ('A28_gases_glow_held', (10, 19))}
# Each forcing against the branch it differs from by that forcing alone, over the first 3-day means.
FAST = {'gases': ('A28_gases', 'A28_cal_control_patched'), 'gases_and_glow': ('A28_gases_glow', 'A28_gases'),
        'gases_and_glow_held': ('A28_gases_glow_held', 'A28_gases'),
        'calibration': ('A28_cal_co2_1022', 'A28_cal_control')}
FAST_MEANS = '10'
# The CO2 equivalent's calibration: PlaSim's all-sky forcing over a branch's first 3-day means, 400 to 1,022 ppm.
CALIBRATION = dict(co2_ppm=1022.0, control_co2_ppm=400.0, forcing_w_m2_first_mean=5.266,
                   forcing_w_m2_first_ten_means=5.297)
DESIGN_SCALE = 0.95                        # the design dimmer passes 0.95 of the titania stack's sunlight
DESIGN_SUNLIGHT_W_M2 = 293.312             # the design sunlight at the top of the air, over the whole Moon
DESIGN_ALBEDO = 0.280                      # the design case's planetary albedo
# The glow's recommended answer (research/array-industry a59f5a9, results/glow.json, results.recommended): both
# halves of each orbit rolled 80 degrees outside the service arc, on the layout re-flown for that law's push.
GLOW = dict(answer='rolled_outside_service_80_rolled_layout', glow_w_m2=3.59739, infrared_return_w_m2=1.74675,
            sunlight_return_w_m2=0.610045, top_infrared_w_m2=dict(day=10.1718, night=0.5165),
            source='research/array-industry a59f5a9, research/studies/solar_shield_array/results/glow.json')
GLOW_HELD_DIMMING_PERCENT = 1.83
BANDS = ((0, 15), (15, 30), (30, 45), (45, 60), (60, 75), (75, 90))
EVIDENCE = ("ExoPlaSim runs of the design case, T21 with 10 layers and a 50 m mixed-layer sea, branched from its "
            "year-29 restart: the gases as their CO2 equivalent in PlaSim's own radiation, the glow as infrared added "
            "at the top of the air and carried down by PlaSim's broadband transmissivities, and the dimming as less "
            "sunlight. The control is the design case's own years 20-29; each run's years 10-19 are its climate. "
            "PlaSim's clouds and convection are Earth-tuned (climate/gcm/README, Clouds), and the rain-fed growth "
            "takes the light limitation study's water-limited land to grow in proportion to the water it evaporates. "
            "The held point is interpolated between the two glow runs, linearly in their dimming.")
READING_RULE = ("Means are area-weighted over the stated surface; precipitation and evaporation mm a day of water, "
                "evaporation upward; changes are run less control, relative changes in per cent of the control. "
                "Bands are |latitude| ranges of land. Dimming is per cent of the design sunlight. Growth is Gt C a "
                "year; its pairs are the light limitation study's [low, high] ends. first_means are W/m2, positive "
                "downward, over the first ten 3-day means of each branch less the branch it differs from (FAST). "
                "glow.held is the plan with the temperature held at the control's; pairs scales the gases' run to "
                "the other pairs in proportion to their warming.")


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def load(path, schema=None):
    data = json.loads((ROOT / path).read_text())
    if schema and data.get('schema') != schema:
        raise ValueError(f'{path}: unexpected schema {data.get("schema")}')
    return data


def rnd(x, sig=6):
    """Six significant figures, which keep temperatures near 295 K to the thousandth of a kelvin."""
    if isinstance(x, dict):
        return {k: rnd(v, sig) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [rnd(v, sig) for v in x]
    if isinstance(x, bool) or not isinstance(x, float) or x == 0.0 or not math.isfinite(x):
        return x
    return round(x, sig - 1 - int(math.floor(math.log10(abs(x)))))


def climatology(run, years):
    """A run's climatology product over the stated years, and its hash."""
    path = ROOT / CLIMATOLOGY.format(run)
    z = np.load(path)
    meta = json.loads(str(z['metadata']))
    if meta['schema'] != CLIMATOLOGY_SCHEMA or meta['run'] != run or tuple(meta['years']) != tuple(years):
        raise ValueError(f'{path.name}: not {run} years {years}; run python -m climate.gcm.climatology {run}:'
                         f'{years[0]}-{years[1]}')
    return {k: z[k] for k in z.files if k != 'metadata'}, meta, sha256(path)


def means(c):
    """The surface means: global, sea and land, and land by band."""
    lat = c['lat'].astype(float)
    w = np.cos(np.radians(lat))[:, None] * np.ones(c['lsm'].shape)
    land = c['lsm'] > 0.5

    def avg(f, m):
        return float((f * w * m).sum() / (w * m).sum())
    every = np.ones(land.shape, bool)
    out = dict(
        surface_k=avg(c['ts_k'], every), land_surface_k=avg(c['ts_k'], land),
        rain_mm_day=dict(all=avg(c['pr_mm_day'], every), land=avg(c['pr_mm_day'], land), sea=avg(c['pr_mm_day'], ~land)),
        evaporation_mm_day=dict(all=avg(c['evap_mm_day'], every), land=avg(c['evap_mm_day'], land),
                                sea=avg(c['evap_mm_day'], ~land)),
        land_soil_water_m=avg(c['mrso_m'], land),
        cloud_cover=dict(all=avg(c['clt_mean'], every), low_deck=avg(c['low_mean'], every),
                         high_deck=avg(c['high_mean'], every)),
        water_vapour_kg_m2=avg(c['prw_kg_m2'], every))
    bands = {}
    for lo, hi in BANDS:
        m = land & ((np.abs(lat) >= lo) & (np.abs(lat) < hi))[:, None]
        bands[f'{lo}-{hi}'] = dict(rain_mm_day=avg(c['pr_mm_day'], m), evaporation_mm_day=avg(c['evap_mm_day'], m),
                                   soil_water_m=avg(c['mrso_m'], m), surface_k=avg(c['ts_k'], m))
    out['land_bands'] = bands
    return out


def change(run, control):
    """Each mean of a run against the control: the difference and, for water, the relative change in per cent."""
    if isinstance(run, dict):
        return {k: change(run[k], control[k]) for k in run}
    return dict(run=run, control=control, change=run - control,
                percent=100.0 * (run - control) / control if control else None)


def series(run, years):
    """The run's yearly global means over the compared years, and their trends a decade."""
    p = json.loads((ROOT / PROGRESS.format(run)).read_text())
    rows = [y for y in p['years'] if years[0] <= y['year'] <= years[1]]
    if len(rows) != years[1] - years[0] + 1:
        raise ValueError(f'{run}: years {years} are not all flown')
    t = np.array([r['year'] for r in rows], float)
    trend = lambda k: float(np.polyfit(t, [r[k] for r in rows], 1)[0] * 10.0)
    return dict(years=list(years), surface_k=[r['surface_k'] for r in rows],
                rain_mm_day=[r['precipitation_mm_day'] for r in rows], toa_net_w_m2=[r['toa_net_w_m2'] for r in rows],
                trend_per_decade=dict(surface_k=trend('surface_k'), rain_mm_day=trend('precipitation_mm_day'),
                                      toa_net_w_m2=trend('toa_net_w_m2')),
                configuration=dict(co2_ppm=p['configuration'].get('co2_ppm'),
                                   top_infrared_w_m2=p['configuration'].get('top_infrared_w_m2'),
                                   sunlight_scale=p['configuration']['model']['sunlight_scale'],
                                   initial_state=p['configuration'].get('initial_state'),
                                   key=p['configuration']['key'])), sha256(ROOT / PROGRESS.format(run))


def layers(extract, run, control):
    """The air's temperature, humidity and cloud by layer, and the rain's convective and large-scale parts, against
    the control."""
    a, b = extract['layers'][run], extract['layers'][control]
    large = lambda r: r['rain_m_s'] - r['convective_rain_m_s']
    return dict(sigma=a['sigma'], temperature_change_k=[x - y for x, y in zip(a['ta_k'], b['ta_k'])],
                humidity_percent=[100.0 * (x / y - 1.0) for x, y in zip(a['hus'], b['hus'])],
                cloud_change=[x - y for x, y in zip(a['cl'], b['cl'])], control_temperature_k=b['ta_k'],
                control_cloud=b['cl'], convective_share_of_rain=a['convective_rain_m_s'] / a['rain_m_s'],
                control_convective_share_of_rain=b['convective_rain_m_s'] / b['rain_m_s'],
                convective_rain_percent=100.0 * (a['convective_rain_m_s'] / b['convective_rain_m_s'] - 1.0),
                large_scale_rain_percent=100.0 * (large(a) / large(b) - 1.0))


def first_means(extract, run, base, glow_w_m2=None):
    """Where a forcing's energy goes before the climate responds: the branch less the branch it differs from, at the
    top of the air and at the ground, and what the air takes between them."""
    a, b = extract['first_means'][run][FAST_MEANS], extract['first_means'][base][FAST_MEANS]
    d = {k: a[k] - b[k] for k in a}
    out = dict(against=base, toa_net=d['ntr'], toa_sunlight=d['rst'], toa_infrared=d['rlut'],
               surface_sunlight=d['rss'], surface_infrared=d['rls'], surface_radiation=d['nbr'],
               surface_latent=d['hfls'], surface_sensible=d['hfss'], air_radiation=d['ntr'] - d['nbr'],
               surface_k=d['ts'])
    if glow_w_m2:
        out['glow_share_reaching_ground'] = d['rls'] / glow_w_m2
    return out


def growth_from(shares, light):
    """The rain-fed land's growth change, band by band, with each band's growth following its evaporation."""
    bands = light['land']['bands']
    out, total = {}, [0.0, 0.0]
    for key, b in bands.items():
        grown = [b['rain_fed_gt_c'][e] * shares[key] for e in (0, 1)]
        out[key] = dict(evaporation_percent=100.0 * shares[key], growth_change_gt_c=grown)
        total = [total[e] + grown[e] for e in (0, 1)]
    living = light['ledger']['total_gt_c']
    return dict(by_band=out, change_gt_c=total,
                share_of_living_moon_percent=[100.0 * total[e] / living[e] for e in (0, 1)])


def rain_fed_growth(diff, light):
    return growth_from({key: diff['land_bands'][key]['evaporation_mm_day']['percent'] / 100.0
                        for key in light['land']['bands']}, light)


def at(x, a, b):
    """Linear interpolation between a (at 0) and b (at 1), through nested dicts and lists."""
    if isinstance(a, dict):
        return {k: at(x, a[k], b[k]) for k in a if k in b}
    if isinstance(a, list):
        return [at(x, p, q) for p, q in zip(a, b)]
    if isinstance(a, float) and isinstance(b, float):
        return a + (b - a) * x
    return a


def glow_dimming(runs, control, absorption, light):
    """The model's response to the glow's dimming from the two glow runs, the dimming that holds the design
    temperature with the glow, and the climate there."""
    ww, held, gases = runs['gases_and_glow'], runs['gases_and_glow_held'], runs['gases']
    t = lambda r: r['means']['surface_k']['run']
    glow_ww = ww['added_dimming_percent'] - gases['added_dimming_percent']
    glow_held = held['added_dimming_percent'] - gases['added_dimming_percent']
    per_percent = (t(held) - t(ww)) / (glow_ww - glow_held)                  # K for each 1% less glow dimming
    hold = lambda target: glow_held - (target - t(held)) / per_percent       # glow dimming at a target temperature
    at_control, at_gases = hold(control['surface_k']), hold(t(gases))
    x = (glow_ww - at_control) / (glow_ww - glow_held)                       # 0 at watt for watt, 1 at the held run
    point = at(x, ww['means'], held['means'])
    growth = rain_fed_growth(point, light)
    infrared = GLOW['glow_w_m2'] + GLOW['infrared_return_w_m2']
    gcm_ground = [ww['first_means']['glow_share_reaching_ground'], held['first_means']['glow_share_reaching_ground']]
    d = absorption['at_design']
    ck_ground = [(GLOW['glow_w_m2'] * d['glow'][k]['ground'] + GLOW['infrared_return_w_m2']
                  * d['returned_infrared']['ground']) / infrared for k in ('150_k', '190_k')]
    rain = point['rain_mm_day']
    glow_rain = {k: rain[k]['percent'] - gases['means']['rain_mm_day'][k]['percent'] for k in ('all', 'land')}
    plan = gases['added_dimming_percent'] + at_control
    return dict(
        infrared_w_m2=infrared, watt_for_watt_dimming_percent=glow_ww, held_run_dimming_percent=glow_held,
        response_k_per_percent=per_percent,
        held=dict(glow_dimming_percent=at_control, plan_dimming_percent=plan, efficacy=at_control / glow_ww,
                  interpolation=x, means=point, rain_fed_growth=growth, glow_rain_percent=glow_rain,
                  glow_rain_percent_per_w_m2=glow_rain['all'] / infrared,
                  glow_growth_percent_of_living_moon=[
                      growth['share_of_living_moon_percent'][e]
                      - gases['rain_fed_growth']['share_of_living_moon_percent'][e] for e in (0, 1)],
                  through_light_percent_of_growth=[light['cost_through_light']['per_percent_dimming'][e] * plan
                                                   for e in (0, 1)]),
        at_the_gases_temperature=dict(glow_dimming_percent=at_gases, efficacy=at_gases / glow_ww),
        reaching_ground=dict(gcm_all_sky=gcm_ground, correlated_k_clear_sky_films_150_to_190_k=ck_ground))


def made_up(runs, glow, light):
    """What the plants' water design makes up: the rise in growth per litre of water the land evaporates that restores
    the rain-fed growth, band by band, for the glow's part (the held point against the gases' run) and for the whole
    plan (the held point against the design case). Rain-fed growth is the water transpired times the carbon fixed
    per litre, so a rise of r in either factor (carbon per litre, or the transpired share of the evaporation) restores
    a fall s in the evaporation when (1 + r)(1 + s) = 1."""
    held, gases = glow['held'], runs['gases']
    bands, out = light['land']['bands'], {}
    def rise(s):
        return 100.0 * (1.0 / (1.0 + s / 100.0) - 1.0)
    for key in bands:
        plan = held['means']['land_bands'][key]['evaporation_mm_day']['percent']
        glow_part = 100.0 * ((1.0 + plan / 100.0) / (1.0 + gases['means']['land_bands'][key]['evaporation_mm_day']
                                                       ['percent'] / 100.0) - 1.0)
        out[key] = dict(evaporation_percent=dict(plan=plan, glow=glow_part),
                        rise_percent=dict(plan=rise(plan), glow=rise(glow_part)))
    def overall(which):
        lost = [sum(bands[k]['rain_fed_gt_c'][e] * (1.0 - 1.0 / (1.0 + out[k]['rise_percent'][which] / 100.0))
                    for k in bands) for e in (0, 1)]
        grown = [sum(b['rain_fed_gt_c'][e] for b in bands.values()) for e in (0, 1)]
        return [100.0 * (1.0 / (1.0 - lost[e] / grown[e]) - 1.0) for e in (0, 1)]
    return dict(by_band=out, rise_percent=dict(plan=overall('plan'), glow=overall('glow')),
                reference='CAM fixes 2-6 times as much carbon per litre as C3 and 1-5 times as much as C4 '
                          '(Hartzell, Bartlett & Porporato 2018)')


def results():
    data = {p: load(p, s) for p, s in INPUTS.items()}
    light = data['research/studies/light_limitation/results/light_limitation.json']
    ledgers = data['research/studies/domain_synthesis/results/ledgers.json']
    absorption = data['research/studies/gases_glow_dimming/results/glow_absorption.json']
    extract = load(EXTRACT, EXTRACT_SCHEMA)
    ctl, _, ctl_sha = climatology(*CONTROL)
    control = means(ctl)
    hashes = {CLIMATOLOGY.format(CONTROL[0]): ctl_sha, EXTRACT: sha256(ROOT / EXTRACT)}
    held = ledgers['gases']['held']
    infrared = GLOW['glow_w_m2'] + GLOW['infrared_return_w_m2']
    light_cost = light['cost_through_light']['per_percent_dimming']
    runs = {}
    for name, (run, years) in RUNS.items():
        c, _, sha = climatology(run, years)
        hashes[CLIMATOLOGY.format(run)] = sha
        s, psha = series(run, years)
        hashes[PROGRESS.format(run)] = psha
        diff = change(means(c), control)
        added = 100.0 * (1.0 - s['configuration']['sunlight_scale'] / DESIGN_SCALE)
        rain = diff['rain_mm_day']
        runs[name] = dict(
            run=run, years=list(years), series=s, means=diff, layers=layers(extract, run, CONTROL[0]),
            first_means=first_means(extract, *FAST[name], glow_w_m2=None if name == 'gases' else infrared),
            added_dimming_percent=added,
            rain_per_percent_dimming=dict(all=rain['all']['percent'] / added, land=rain['land']['percent'] / added),
            rain_fed_growth=rain_fed_growth(diff, light),
            through_light_percent_of_growth=[light_cost[e] * added for e in (0, 1)])
    _, psha = series(*CONTROL)
    hashes[PROGRESS.format(CONTROL[0])] = psha
    glow = glow_dimming(runs, control, absorption, light)
    glow['made_up_by_plants'] = made_up(runs, glow, light)
    # The other pairs from the gases' run, in proportion to their warming.
    weighed = ledgers['gases']['pairs_weighed']
    g = runs['gases']
    pairs = {}
    for key, p in weighed.items():
        k = p['allsky_w_m2'] / weighed['middle']['allsky_w_m2']
        pairs[key] = dict(methane_ppm=p['methane_ppm'], n2o_ppm=p['n2o_ppm'], added_dimming_percent=p['added_dimming_percent'],
                          rain_percent=dict(all=k * g['means']['rain_mm_day']['all']['percent'],
                                            land=k * g['means']['rain_mm_day']['land']['percent']),
                          rain_fed_growth_percent_of_living_moon=[k * x for x in
                                                                  g['rain_fed_growth']['share_of_living_moon_percent']],
                          through_light_percent_of_growth=[light_cost[e] * p['added_dimming_percent'] for e in (0, 1)])
    calibration = first_means(extract, *FAST['calibration'])
    return dict(
        schema=SCHEMA,
        producer=dict(files={'research/studies/gases_glow_dimming/run.py': sha256(__file__)},
                      inputs={p: sha256(ROOT / p) for p in INPUTS}, ignored_inputs=hashes),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        setup=dict(control=dict(run=CONTROL[0], years=list(CONTROL[1])), calibration=CALIBRATION,
                   calibration_check_first_ten_means_w_m2=calibration['toa_net'],
                   co2_equivalent_ppm=runs['gases']['series']['configuration']['co2_ppm'],
                   gases_forcing_allsky_w_m2=held['forcing_w_m2']['allsky'],
                   gases_dimming_percent=held['added_dimming_percent'], glow=GLOW, glow_infrared_w_m2=infrared,
                   glow_watt_for_watt_dimming_percent=100.0 * infrared / ((1.0 - DESIGN_ALBEDO) * DESIGN_SUNLIGHT_W_M2),
                   glow_held_dimming_percent=GLOW_HELD_DIMMING_PERCENT,
                   design_scale=DESIGN_SCALE, design_sunlight_w_m2=DESIGN_SUNLIGHT_W_M2),
        control=rnd(control), runs=rnd(runs), glow=rnd(glow), pairs=rnd(pairs))


def dumps(product) -> str:
    return json.dumps(product, indent=2, allow_nan=False) + '\n'


def main(argv=None) -> int:
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(dumps(results()))
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
