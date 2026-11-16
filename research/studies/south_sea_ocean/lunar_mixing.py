"""The south sea's own interior mixing on the Moon, and the deep oxygen run that uses it.

    R=/media/projectspace/terluna-research/ocean/run.sh; M=research.studies.south_sea_ocean.lunar_mixing
    $R -m $M inertial-run       # one model year from year 50 with daily surface output (about 15 minutes)
    $R -m $M inputs             # the wind's, the tide's and the floor's inputs: results/lunar_mixing_inputs_1deg.json
    $R -m $M trial              # one model year with IDEMIX and the oxygen tracers (about 15 minutes)
    $R -m $M run --years 50     # fifty model years from year 50 (about 9 hours)

The 50-year circulation run mixed the sea's interior with a background diffusivity of 1e-5 m2/s, the open ocean's
value at thermocline depths on Earth. On the Moon the interior's mixing is what the power going into the sea's
internal waves can pay for, under the Moon's gravity and rotation. MITgcm's IDEMIX closure (Olbers and Eden 2013,
climate/ocean/stratified.py) turns that power into mixing by place and depth, and the background falls to the
molecular value. This module gives it the Moon's inputs:

- the wind: `inertial-run` continues the 50-year run's sea, unchanged, for one model year from its year-50 state,
  writing means every 0.99 days of the top layer's flow, the stress the model applies, the wind's work at every
  step and the mixed layer. Under the Moon's slow rotation the inertial period 2 pi / |f| is two to eight weeks over
  the sea, and the run's monthly means averaged those motions out. `inputs` takes the wind's work on the motions
  turning the inertial way within a quarter of the local inertial frequency; IDEMIX takes a fifth of it into the
  internal waves, as it does on Earth;
- the tide: the monthly tide trial's pair of runs gives the tidal current at the floor, fitted to the month's and
  the fortnight's lines; the power each line puts into internal waves follows the linear theory's form as Jayne and
  St. Laurent (2001) parameterised it on Earth, half rho kappa h^2 N U^2 with kappa = 2 pi / 10 km and h^2 the
  variance of LOLA's floor about a plane within each cell, times sqrt(1 - f^2 / omega^2), and nothing where the
  line's frequency lies below the local inertial frequency, where no internal wave can carry it;
- the floor's heat: 11 mW/m2 everywhere, the middle of the 8-14 mW/m2 predicted for the Schrodinger basin inside
  South Pole-Aitken (Santangelo et al. 2025, a search summary); the Apollo 15 and 17 sites measured 21 and 14.

The oxygen. `trial` and `run` carry four tracers from the year-50 state: `old`, the 50-year run's dye continued;
`new`, a dye started at 1 at year 50; and `used` and `used_lake`, the oxygen spent since the water last left the
surface, started at 0, held at 0 in the top layer and given the demand of the sea oxygen study's nutrient-poor
column on each level, in the ocean's form (sinking matter decaying as depth^-0.86 from 100 m) and in the lake
model's form, both as the 50-year reading's estimate placed them (circulation.oxygen_estimate). The demand does not
depend on the oxygen left, so `used` scales with the export and gives every nutrient state at once, and water
holds saturation less `used` where that is positive. Water that has stayed below the surface since year 50 started
the run counted as saturated; counting it as holding none instead lowers it by saturation times `new`, so the two
bound the oxygen until `new` has fallen.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np

from climate.ocean import lunar, mitgcm
from climate.ocean import stratified as ST
from research.studies.south_sea_ocean import circulation as C
from research.studies.south_sea_ocean import run as R

HERE = Path(__file__).resolve().parent
RUNS = C.RUNS
FROM_YEAR = 50.0                        # the 50-year run's end state starts every run here
INERTIAL_SAMPLES = 360                  # outputs a model year for the inertial band, 0.993 days apart
INERTIAL_BAND = 0.25                    # the band's half-width as a share of the local |f|
INERTIAL_BANDS_REPORTED = (0.15, 0.25, 0.40)
TIDE_LINES_DAYS = (27.32, 13.66)        # the tide's monthly and fortnightly lines (the trial's 30 days cannot part
                                        # the 27.2, 27.3 and 27.6-day months)
TIDE_SKIP_DAYS = 5                      # the trial's start, before the sea follows its tide
KAPPA_PER_M = 2.0 * np.pi / 10000.0     # Jayne and St. Laurent's (2001) topographic wavenumber, fitted on Earth
HEAT_FLOW_W_M2 = 0.011                  # the floor's heat; 0.008-0.014 predicted at Schrodinger (Santangelo et al. 2025)
HEAT_FLOW_RANGE_W_M2 = (0.006, 0.021)   # Diviner's lower bound for one region to the Apollo 15 site (search summaries)
MOLECULAR_M2_S = 1.4e-7                 # heat's molecular diffusivity in water at about 28 C: the background left
O2_PER_C = 32.0 / 12.0                  # oxygen spent per carbon, as the 50-year reading counted it
EXPORTS_G_C = (34.0, 25.0, 12.0, 6.0)   # carbon sinking below 100 m, g/m2 a year: the nutrient-poor column's, and
                                        # lower, to Earth's poorest gyres (about 12) and below
REFERENCE_STATE = 'oligotrophic'        # the sea oxygen study's nutrient-poor column, whose demand the tracers carry
TRIAL_TRACERS = ('old', 'new', 'used', 'used_lake')
RUN_TRACERS = ('new', 'used')           # the run keeps what the oxygen needs: each tracer costs machine time
INPUTS = HERE / 'results' / 'lunar_mixing_inputs_1deg.json'
MAPS = RUNS / 'lunar_mixing_inputs_1deg.npz'
CODE = (Path(__file__), HERE.parents[2] / 'climate' / 'ocean' / 'stratified.py')


def code_hash():
    h = hashlib.sha256()
    for path in CODE:
        h.update(path.read_bytes())
    return h.hexdigest()


# ---------------------------------------------------------------- restarts from the 50-year run

def restart_sea(name, step=1.0, from_year=FROM_YEAR):
    """The whole sea as the 50-year run set it, under another name, to begin from that run's pickup at the end of
    model year `from_year`. Returns the sea, its setting, the parent's folder and the iteration."""
    parent = mitgcm.RUNS / 'mitgcm' / f'south_sea_{step:g}deg'
    deltat = C.parent_deltat(step)
    iteration = int(round(from_year * C.YEAR_S / deltat))
    if not (parent / f'pickup.{iteration:010d}.data').exists():
        raise SystemExit(f'{parent} holds no pickup at the end of model year {from_year:g}')
    sea, setting = C.build(step, deltat, name=name)
    setting.update(parent=parent.name, started_from_iteration=iteration, started_from_year=from_year)
    return sea, setting, parent, iteration


def begin(sea, setting, parent, iteration):
    """Write a restarted run's folder unless it has begun already."""
    if not (sea.folder / 'restart.json').exists():
        sea.write(fresh=True)
        sea.restart_from(parent, iteration)
        (sea.folder / 'setting.json').write_text(json.dumps(setting, indent=1) + '\n')


def inertial_name(step, from_year):
    return f'south_sea_{step:g}deg_inertial_from{from_year:g}'


def inertial_run(args):
    """One model year from the year-50 state with the 50-year run's own settings, writing means for the inertial
    band."""
    sea, setting, parent, iteration = restart_sea(inertial_name(args.step, args.from_year), args.step, args.from_year)
    sea.dye = False                     # the flow does not feel the dye
    interval = C.YEAR_S / INERTIAL_SAMPLES
    sea.extra_lists = ((('UVEL', 'VVEL'), 'inertial_top', interval, 1),
                       (('oceTAUX', 'oceTAUY', 'GGL90tau', 'MXLDEPTH'), 'inertial_surface', interval, None))
    begin(sea, setting, parent, iteration)
    t0 = time.time()
    sea.run(iteration * sea.deltat + args.years * C.YEAR_S, C.YEAR_S)
    print(f'{sea.folder}: {time.time() - t0:.0f} s')
    return 0


# ---------------------------------------------------------------- the wind's input

def inertial_band(w, dt, f, half_width):
    """The part of the flows w = u + i v (time first) turning the inertial way within half_width |f| of the local
    inertial frequency. A free inertial motion turns as exp(-i f t); numpy's transform puts exp(i omega t) at
    +omega, so the band sits at omega = -f."""
    omega = 2.0 * np.pi * np.fft.fftfreq(w.shape[0], d=dt)
    spectrum = np.fft.fft(w, axis=0)
    centre = -f[None]
    keep = np.abs(omega[:, None, None] - centre) <= half_width * np.abs(centre)
    return np.fft.ifft(spectrum * keep, axis=0)


def wind_power(sea, f):
    """The wind's work on the sea from the inertial run's outputs (W/m2 on the grid): its whole work at every step
    (the model's GGL90tau), its work on the flow as the outputs resolve it, on the monthly-mean flow, and on the
    inertial band of each width in INERTIAL_BANDS_REPORTED; with the mixed layer's mean depth."""
    tops = sorted(sea.folder.glob('inertial_top.*.data'))
    surfaces = sorted(sea.folder.glob('inertial_surface.*.data'))
    if len(tops) < INERTIAL_SAMPLES or len(surfaces) != len(tops):
        raise SystemExit(f'{sea.folder} holds {len(tops)} of the {INERTIAL_SAMPLES} outputs a model year needs')
    tops, surfaces = tops[-INERTIAL_SAMPLES:], surfaces[-INERTIAL_SAMPLES:]
    u, v, tx, ty, work, mld = ([] for _ in range(6))
    for a, b in zip(tops, surfaces):
        top, surface = mitgcm.read_mds(str(a)[:-5]), mitgcm.read_mds(str(b)[:-5])
        uc, vc = C.at_centres(np.squeeze(top['UVEL']), np.squeeze(top['VVEL']))
        sx, sy = C.at_centres(np.squeeze(surface['oceTAUX']), np.squeeze(surface['oceTAUY']))
        u.append(uc), v.append(vc), tx.append(sx), ty.append(sy)
        work.append(np.squeeze(surface['GGL90tau']) * sea.rho_const)
        mld.append(np.squeeze(surface['MXLDEPTH']))
    u, v, tx, ty = (np.array(x, dtype=float) for x in (u, v, tx, ty))
    dt = C.YEAR_S / INERTIAL_SAMPLES
    out = dict(every_step=np.mean(work, axis=0), resolved=(tx * u + ty * v).mean(axis=0),
               mixed_layer_m=np.mean(mld, axis=0))
    month = lambda x: x.reshape(12, -1, *x.shape[1:]).mean(axis=1)
    out['monthly_mean_flow'] = (month(tx) * month(u) + month(ty) * month(v)).mean(axis=0)
    for width in INERTIAL_BANDS_REPORTED:
        band = inertial_band(u + 1j * v, dt, f, width)
        out[f'inertial_{width:g}'] = (tx * band.real + ty * band.imag).mean(axis=0)
    return out


# ---------------------------------------------------------------- the tide's input

def tide_floor_lines(sea_off, sea_on, skip_days=TIDE_SKIP_DAYS):
    """The mean square tidal current at the floor (m2/s2) on each of TIDE_LINES_DAYS, from the tide trial's pair:
    each cell's daily-mean difference at its deepest open level, fitted to a mean and the lines' cosines and sines."""
    volume = C.cell_volumes(sea_on)
    floor = np.maximum((volume > 0).sum(axis=0) - 1, 0)
    days, du, dv = [], [], []
    for path in sorted(sea_on.folder.glob('deep.*.data')):
        stem = path.name[:-5]
        day = int(stem.split('.')[1]) * sea_on.deltat / 86400.0
        if day <= skip_days or not (sea_off.folder / path.name).exists():
            continue
        on, off = mitgcm.read_mds(str(path)[:-5]), mitgcm.read_mds(str(sea_off.folder / stem))
        u_on, v_on = C.at_centres(on['UVEL'], on['VVEL'])
        u_off, v_off = C.at_centres(off['UVEL'], off['VVEL'])
        at_floor = lambda x: np.take_along_axis(x, floor[None], axis=0)[0]
        du.append(at_floor(u_on - u_off)), dv.append(at_floor(v_on - v_off))
        days.append(day - 0.5)                          # a daily mean stands for its day's middle
    lines = line_mean_squares(np.array(days) * 86400.0, np.array(du), np.array(dv))
    total = np.mean(np.array(du) ** 2 + np.array(dv) ** 2, axis=0)
    return dict(mean_square_by_line=lines, mean_square_all=total, days=[float(days[0]), float(days[-1])])


def line_mean_squares(t, u, v, periods_days=TIDE_LINES_DAYS):
    """Each line's mean square current (m2/s2) in series u, v (time first): a least-squares fit of a mean and each
    line's cosine and sine to both components."""
    columns = [np.ones_like(t)]
    for period in periods_days:
        omega = 2.0 * np.pi / (period * 86400.0)
        columns += [np.cos(omega * t), np.sin(omega * t)]
    design = np.stack(columns, axis=1)
    fits = {}
    for name, series in (('u', u), ('v', v)):
        coef, *_ = np.linalg.lstsq(design, series.reshape(len(t), -1), rcond=None)
        fits[name] = coef.reshape(len(columns), *series.shape[1:])
    return [0.5 * sum(fits[k][1 + 2 * n] ** 2 + fits[k][2 + 2 * n] ** 2 for k in ('u', 'v'))
            for n in range(len(periods_days))]


def floor_roughness(grid, rotation, frame):
    """Each cell's floor roughness (m2): the variance of LOLA's depths at 16 px/deg within the cell about the plane
    that fits them, where the cell holds at least 12 of the sea's pixels."""
    from geography import topography as tp
    from research.studies.south_sea_ocean import inflows as I
    from climate.gcm.turned import to_turned
    atlas = json.loads((C.ROOT / I.ATLAS).read_text())
    level = float(atlas['sea_level_m'])
    h, _, _ = tp.height_above_geoid(C.PPD, atlas['grid']['geoid_degree'])
    rows, cols = h.shape
    lat = 90.0 - (np.arange(rows) + 0.5) / C.PPD
    lon = (np.arange(cols) + 0.5) / C.PPD
    row_area = (lunar.RADIUS / 1000.0) ** 2 * np.radians(1.0 / C.PPD) ** 2 * np.cos(np.radians(lat))
    body = I.label_waters(h < level)
    named, _ = I.name_seas(body, row_area * 1e6, lat, lon, atlas['bodies'])
    pr, pc = np.nonzero(body == named[R.SEA_NAME]['label'])
    z = level - h[pr, pc]
    glat, glon = frame.to_grid(*to_turned(lat[pr], lon[pc], rotation))
    j, i, inside = C.cells(glat, glon, grid.xg, grid.yg)
    cell = (j * grid.nx + i)[inside]
    x = np.radians(glon[inside] - grid.xc[i[inside]]) * np.cos(np.radians(glat[inside])) * lunar.RADIUS
    y = np.radians(glat[inside] - grid.yc[j[inside]]) * lunar.RADIUS
    z = z[inside]
    return plane_residual_variance(cell, x, y, z, grid.nx * grid.ny).reshape(grid.ny, grid.nx)


def plane_residual_variance(cell, x, y, z, n_cells, least=12):
    """For each cell, the variance of z about the plane in x, y that fits the points it holds; 0 for a cell with
    fewer than `least` points."""
    s = lambda w: np.bincount(cell, weights=w, minlength=n_cells)
    n = s(np.ones_like(z))
    sx, sy, sz = s(x), s(y), s(z)
    sxx, syy, sxy, sxz, syz, szz = s(x * x), s(y * y), s(x * y), s(x * z), s(y * z), s(z * z)
    out = np.zeros(n_cells)
    for k in np.flatnonzero(n >= least):
        a = np.array([[n[k], sx[k], sy[k]], [sx[k], sxx[k], sxy[k]], [sy[k], sxy[k], syy[k]]])
        rhs = np.array([sz[k], sxz[k], syz[k]])
        try:
            beta = np.linalg.solve(a, rhs)
        except np.linalg.LinAlgError:
            continue
        out[k] = max(0.0, (szz[k] - beta @ rhs) / n[k])
    return out


def floor_buoyancy_frequency(sea, deep, volume):
    """N (1/s) between each column's two deepest open levels in a yearly mean, by potential density; 0 where the
    water there is not stably layered or the column has one level."""
    open_ = volume > 0
    floor = np.maximum(open_.sum(axis=0) - 1, 0)
    above = np.maximum(floor - 1, 0)
    pick = lambda x, k: np.take_along_axis(x, k[None], axis=0)[0]
    rho_f = R.seawater_density(pick(deep['SALT'], floor), pick(deep['THETA'], floor))
    rho_a = R.seawater_density(pick(deep['SALT'], above), pick(deep['THETA'], above))
    mid = np.cumsum(sea.delr) - 0.5 * sea.delr
    dz = np.maximum(mid[floor] - mid[above], 1e-9)
    n2 = lunar.GRAVITY * (rho_f - rho_a) / (sea.rho_const * dz)
    return np.where((floor > above) & (n2 > 0), np.sqrt(np.maximum(n2, 0.0)), 0.0)


def tidal_power(sea_off, sea_on, parent, n_floor, h2, f):
    """The power the tide puts into internal waves at the floor (W/m2): for each line, half rho kappa h^2 N U^2
    sqrt(1 - f^2 / omega^2) where omega exceeds |f|."""
    lines = tide_floor_lines(sea_off, sea_on)
    total = np.zeros_like(h2)
    by_line = []
    for period, u2 in zip(TIDE_LINES_DAYS, lines['mean_square_by_line']):
        e = line_power(parent.rho_const, h2, n_floor, u2, f, period)
        total += e
        by_line.append(dict(period_days=period, power=e, mean_square=u2,
                            radiating=np.abs(f) < 2.0 * np.pi / (period * 86400.0)))
    return total, by_line, lines


def line_power(rho, h2, n, u2, f, period_days, kappa=KAPPA_PER_M):
    """The power (W/m2) one tidal line puts into internal waves over a rough floor: half rho kappa h^2 N U^2 times
    sqrt(1 - f^2 / omega^2), and none where omega lies below |f|."""
    omega = 2.0 * np.pi / (period_days * 86400.0)
    return 0.5 * rho * kappa * h2 * n * u2 * np.sqrt(np.clip(1.0 - (np.asarray(f) / omega) ** 2, 0.0, None))


def inputs(args):
    """The Moon's three inputs on the sea's grid: results/lunar_mixing_inputs_1deg.json and the maps beside the
    runs (MAPS)."""
    from climate.gcm.turned import to_turned  # noqa: F401  (floor_roughness's)
    name = inertial_name(args.step, args.from_year)
    run_sea = C.sea_in(name)
    setting = json.loads((run_sea.folder / 'setting.json').read_text())
    frame = ST.Frame(*setting['frame_centre_turned'])
    f = C.coriolis(run_sea.grid, frame)['centre']
    wet, area = run_sea.grid.wet, run_sea.grid.area
    wind = wind_power(run_sea, f)
    rotation = np.array(json.loads((C.ROOT / R.FRAME).read_text())['frame']['rotation_today_to_turned'])
    h2 = floor_roughness(run_sea.grid, rotation, frame)
    parent = C.sea_in(f'south_sea_{args.step:g}deg')
    volume = C.cell_volumes(parent)
    iteration = int(setting['started_from_iteration'])
    deep = mitgcm.read_mds(str(parent.folder / f'deep.{iteration:010d}'))
    n_floor = floor_buoyancy_frequency(parent, deep, volume)
    pair = [C.sea_in(f'south_sea_{args.step:g}deg_tide_pair_{k}') for k in ('off', 'on')]
    tide, by_line, lines = tidal_power(pair[0], pair[1], parent, n_floor, h2, f)
    heat = np.where(wet, HEAT_FLOW_W_M2, 0.0)
    gw = lambda x: float((np.where(wet, x, 0.0) * area).sum() / 1e9)
    pct = lambda x, q=(5, 50, 95): [float(v) for v in np.percentile(x[wet], q)]
    idemix_wind = np.where(wet, np.maximum(wind[f'inertial_{INERTIAL_BAND:g}'], 0.0), 0.0)
    np.savez_compressed(MAPS, wind=idemix_wind, tidal=np.where(wet, tide, 0.0), geothermal=heat, h2=h2,
                        n_floor=n_floor, f=f, **{f'wind_{k}': v for k, v in wind.items()})
    turned_lat, _ = frame.to_true(*np.meshgrid(run_sea.grid.yc, run_sea.grid.xc, indexing='ij'))
    bands = {}
    for south, north in C.WARMTH_BANDS:
        inside = wet & (-turned_lat >= south) & (-turned_lat < north)
        bands[f'{south}-{north} S'] = dict(
            inertial_period_days=[float(x) for x in np.percentile(2 * np.pi / np.abs(f[inside]) / 86400.0, (0, 100))],
            inertial_mw=float((idemix_wind * area)[inside].sum() / 1e6),
            inertial_mw_m2_mean=float((idemix_wind * area)[inside].sum() / area[inside].sum() * 1e3),
            every_step_mw=float((wind['every_step'] * area)[inside].sum() / 1e6))
    out = dict(
        schema='terluna.research.south-sea-lunar-mixing-inputs/1',
        evidence=('The wind\'s work from one model year of the south sea\'s circulation run (MITgcm on one-degree '
                  'cells, lunar gravity and rotation, the design year\'s fluxes) continued from its year-50 state '
                  'with outputs every 0.99 days; the tide\'s from the 35-day pair of tide runs and LOLA\'s floor, '
                  'through a parameterisation fitted on Earth; the floor\'s heat from published measurements and a '
                  'model prediction near the sea, not measured under it. Computed inputs, not validated.'),
        reading_rule=('Powers in GW over the sea\'s wet cells unless named otherwise. wind: every_step is the '
                      'model\'s own wind work at every step (GGL90tau times the reference density); resolved is the '
                      'work on the flow as the 0.99-day means hold it; monthly_mean_flow on twelfths of the year; '
                      'inertial_<w> on the flow turning the inertial way within w |f| of the local inertial '
                      'frequency. IDEMIX takes inertial_0.25, of which its own setting passes a fifth into the '
                      'internal waves. tide: the power into internal waves by line and in all; lines below the local '
                      'inertial frequency give none. The maps are in the npz named by maps.'),
        producer=dict(module='research.studies.south_sea_ocean.lunar_mixing', code_sha256=code_hash()),
        run=name, maps=str(MAPS), inertial_samples=INERTIAL_SAMPLES, idemix_band=INERTIAL_BAND,
        wind_gw={k: gw(v) for k, v in wind.items() if k != 'mixed_layer_m'},
        wind_by_band_of_latitude=bands,
        mixed_layer_m=pct(wind['mixed_layer_m']),
        tide=dict(gw=gw(tide), lines=[dict(period_days=b['period_days'], gw=gw(b['power']),
                                           radiating_share_of_area=float(area[wet & b['radiating']].sum()
                                                                         / area[wet].sum()),
                                           floor_current_rms_m_s=pct(np.sqrt(b['mean_square'])))
                                      for b in by_line],
                  floor_current_rms_all_m_s=pct(np.sqrt(lines['mean_square_all'])), fitted_days=lines['days'],
                  roughness_rms_m=pct(np.sqrt(h2)), floor_n_per_s=pct(n_floor),
                  kappa_per_m=KAPPA_PER_M),
        floor_heat=dict(w_m2=HEAT_FLOW_W_M2, range_w_m2=list(HEAT_FLOW_RANGE_W_M2), gw=gw(heat)),
        idemix_inputs_gw=dict(wind=gw(idemix_wind), tidal=gw(tide)))
    INPUTS.write_text(json.dumps(out, indent=1) + '\n')
    print(json.dumps({k: out[k] for k in ('wind_gw', 'wind_by_band_of_latitude', 'mixed_layer_m', 'tide',
                                           'floor_heat', 'idemix_inputs_gw')}, indent=1))
    return 0


# ---------------------------------------------------------------- the runs with the Moon's mixing

def demand_by_level(parent, volume):
    """The reference column's demand on each model level (mg/L a second) in the ocean's form and the lake model's,
    as the 50-year reading's estimate placed them."""
    basin = C.hypsometry(parent, volume)
    sc = C.central_scenarios()[REFERENCE_STATE]
    out = {}
    for form, demand in (('ocean', C.demand_ocean_form(basin, sc)), ('lake', C.demand_lake_form(basin, sc))):
        per_level = np.zeros(len(parent.delr))
        per_level[:basin.n] = (demand['water'] + demand['floor']) / R.YEAR_S
        out[form] = per_level
        if form == 'ocean':
            out['export_g_o2_m2_yr'] = float(demand['export_g_m2_yr'])
    return out


def mixing_sea(name, args, months=False, carried=TRIAL_TRACERS):
    """The sea from the year-50 state with IDEMIX on the Moon's inputs, the floor's heat, a molecular background,
    and the tracers named in `carried` (of TRIAL_TRACERS)."""
    sea, setting, parent_folder, iteration = restart_sea(name, args.step, args.from_year)
    maps = np.load(MAPS)
    parent = C.sea_in(f'south_sea_{args.step:g}deg')
    volume = C.cell_volumes(parent)
    demand = demand_by_level(parent, volume)
    shape = (len(sea.delr),) + sea.grid.depth.shape
    old = mitgcm.read_mds(str(parent_folder / f'pickup_ptracers.{iteration:010d}'))['pTr01']
    sea.dye = False
    sea.diff_r = MOLECULAR_M2_S
    sea.idemix = dict(tidal=maps['tidal'], wind=maps['wind'])
    sea.geothermal = maps['geothermal']
    sea.deep_fields = sea.deep_fields + ('GGL90Kr', 'IDEMIX_K', 'IDEMIX_E')
    if months:
        sea.deep_output_s = C.MONTH_S
    every = dict(old=ST.Tracer('old', start=old), new=ST.Tracer('new', start=1.0),
                 used=ST.Tracer('used', start=0.0, source=np.broadcast_to(demand['ocean'][:, None, None], shape).copy()),
                 used_lake=ST.Tracer('used_lake', start=0.0,
                                     source=np.broadcast_to(demand['lake'][:, None, None], shape).copy()))
    sea.tracers = tuple(every[k] for k in carried)
    sea.tracers_from = iteration
    setting.update(idemix=dict(inputs=str(INPUTS), band=INERTIAL_BAND), background_diffusivity_m2_s=MOLECULAR_M2_S,
                   floor_heat_w_m2=HEAT_FLOW_W_M2, reference_state=REFERENCE_STATE,
                   reference_export_g_o2_m2_yr=demand['export_g_o2_m2_yr'], tracers=list(carried))
    return sea, setting, parent_folder, iteration


def note_year(sea, done, end):
    """After each chunk: the run's state in one line of its progress file."""
    deep = C.deep_means(sea)
    volume = C.cell_volumes(sea)
    mid = np.cumsum(sea.delr) - 0.5 * sea.delr
    eels = (mid >= 270.0) & (mid <= 600.0)
    layer = lambda x: float((x[eels] * volume[eels]).sum() / volume[eels].sum())
    line = dict(year=end * sea.deltat / C.YEAR_S, wall=time.strftime('%Y-%m-%d %H:%M'), monitor=C.monitor(sea.folder))
    if deep is not None:
        tr = {t.name: deep[f'TRAC{i:02d}'] for i, t in enumerate(sea.tracers, start=1)}
        for name, values in tr.items():
            if name in ('old', 'new'):
                line[f'{name}_water_share'] = C.old_water(sea, values, volume)
            else:
                line[f'{name}_270_600_m_mg_l'] = layer(values)
        line.update(diffusivity_by_depth_m2_s=C.by_depth(deep['GGL90Kr'], volume), theta_270_600_m_c=layer(deep['THETA']))
    with open(sea.folder / 'progress.jsonl', 'a') as f:
        f.write(json.dumps(line) + '\n')


def trial(args):
    """One model year with the Moon's mixing and the oxygen tracers, monthly means of every level: that it runs,
    what it costs, and how the internal waves' energy and the mixing settle."""
    sea, setting, parent, iteration = mixing_sea(f'south_sea_{args.step:g}deg_lunar_mixing_trial', args, months=True)
    sea.write(fresh=True)
    sea.restart_from(parent, iteration)
    (sea.folder / 'setting.json').write_text(json.dumps(setting, indent=1) + '\n')
    t0 = time.time()
    sea.run(iteration * sea.deltat + args.years * C.YEAR_S, C.MONTH_S, after_chunk=note_year)
    wall = time.time() - t0
    steps = int(round(args.years * C.YEAR_S / sea.deltat))
    print(f'{steps} steps in {wall:.0f} s: {wall / steps:.3f} s a step, '
          f'{wall / steps * C.YEAR_S / sea.deltat / 3600.0:.2f} hours a model year')
    return 0


def run(args):
    """The run itself, a model year a chunk, resumable."""
    sea, setting, parent, iteration = mixing_sea(f'south_sea_{args.step:g}deg_lunar_mixing', args, carried=RUN_TRACERS)
    begin(sea, setting, parent, iteration)
    sea.run(iteration * sea.deltat + args.years * C.YEAR_S, C.YEAR_S, after_chunk=note_year)
    return 0


def oxygen_at_exports(saturation, used, new, reference_g_c):
    """The oxygen (mg/L) where the water holds saturation less the oxygen spent, scaled from the reference export to
    each of EXPORTS_G_C, down to none: the upper bound counts the water below the surface at the run's start as
    saturated then, the lower bound as holding none (saturation times the dye less)."""
    out = {}
    for export in EXPORTS_G_C:
        high = np.clip(saturation - used * export / reference_g_c, 0.0, None)
        out[export] = dict(high=high, low=np.clip(high - saturation * new, 0.0, None))
    return out


def read(args):
    """What the oxygen run has done so far: each model year's mixing, renewal, oxygen spent and temperature by depth,
    and in its latest year the oxygen by depth at each export with its bounds, as results/lunar_mixing_1deg.json."""
    name = f'south_sea_{args.step:g}deg_lunar_mixing'
    sea = C.sea_in(name)
    setting = json.loads((sea.folder / 'setting.json').read_text())
    carried = setting['tracers']
    volume = C.cell_volumes(sea)
    mid = np.cumsum(sea.delr) - 0.5 * sea.delr
    true_lat, _ = ST.Frame(*setting['frame_centre_turned']).to_true(*np.meshgrid(sea.grid.yc, sea.grid.xc, indexing='ij'))
    estimate = json.loads((HERE / 'results' / f'circulation_{args.step:g}deg.json').read_text())['oxygen_estimate']
    saturation = estimate['saturation_mg_l']
    reference_g_c = setting['reference_export_g_o2_m2_yr'] / O2_PER_C
    eels = (mid >= 270.0) & (mid <= 600.0)
    layer = lambda x: float((x[eels] * volume[eels]).sum() / volume[eels].sum())
    years, last = [], None
    for path in sorted(sea.folder.glob('deep.*.data')):
        deep = mitgcm.read_mds(str(path)[:-5])
        tr = {k: deep[f'TRAC{i:02d}'] for i, k in enumerate(carried, start=1)}
        years.append(dict(year=int(path.name.split('.')[1]) * sea.deltat / C.YEAR_S,
                          diffusivity_by_depth_m2_s=C.by_depth(deep['GGL90Kr'], volume),
                          wave_diffusivity_by_depth_m2_s=C.by_depth(deep['IDEMIX_K'], volume),
                          new_water_by_depth=C.by_depth(tr['new'], volume),
                          new_water_share=C.old_water(sea, tr['new'], volume),
                          used_by_depth_mg_l=C.by_depth(tr['used'], volume),
                          used_270_600_m_mg_l=layer(tr['used']), theta_270_600_m_c=layer(deep['THETA']),
                          warmth=C.warmth(sea, deep['THETA'], volume, true_lat)))
        last = tr
    oxygen = None
    if last is not None:
        oxygen = {}
        for export, bounds in oxygen_at_exports(saturation, last['used'], last['new'], reference_g_c).items():
            entry = {}
            for side, o2 in bounds.items():
                by_depth = C.by_depth(o2, volume)
                low = [i for i, v in enumerate(by_depth) if v is not None and v < C.HYPOXIC_MG_L]
                entry[side] = dict(by_depth_mg_l=by_depth, eels_layer_mg_l=layer(o2),
                                   eels_layer_share_above_2_mg_l=float((volume[eels] * (o2[eels] >= C.HYPOXIC_MG_L)).sum()
                                                                       / volume[eels].sum()),
                                   under_2_mg_l_from_m=float(mid[low[0]]) if low else None)
            oxygen[f'{export:g} g C'] = entry
    out = dict(schema='terluna.research.south-sea-lunar-mixing/1',
               producer=dict(module='research.studies.south_sea_ocean.lunar_mixing', code_sha256=code_hash()),
               evidence=('MITgcm at lunar gravity and rotation on one-degree cells from the 50-year run\'s year-50 state, '
                         'with IDEMIX on the power the sea\'s own wind and tide put into internal waves, the floor\'s '
                         'heat, a molecular background, a dye and the oxygen spent at the sea oxygen study\'s demand in '
                         'an ocean\'s placement. Its constants and the demand\'s rates are Earth\'s; not validated.'),
               reading_rule=('years: each model year\'s mean at each level (level_mid_m) over the model\'s open cells: '
                             'the closure\'s diffusivity and the part the internal waves give, the dye started at 1 at '
                             'the run\'s start (new_water), the oxygen spent (mg/L, at the reference export) and the '
                             'temperature. oxygen: in the latest year, saturation at the cold end\'s surface less the '
                             'oxygen spent scaled to each carbon export below 100 m; high counts the water below the '
                             'surface at the start as saturated then, low as holding none. The eels\' layer is '
                             '270-600 m.'),
               setting=setting, level_mid_m=[float(x) for x in mid], saturation_mg_l=saturation,
               reference_export_g_c_m2_yr=reference_g_c, model_years=len(years), years=years, oxygen=oxygen)
    path = HERE / 'results' / f'lunar_mixing_{args.step:g}deg.json'
    path.write_text(json.dumps(out, indent=1) + '\n')
    print(f'{len(years)} model years read into {path}')
    return 0


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest='command', required=True)
    for name, fn in (('inertial-run', inertial_run), ('inputs', inputs), ('trial', trial), ('run', run), ('read', read)):
        s = sub.add_parser(name)
        s.add_argument('--step', type=float, default=1.0, help='the grid\'s spacing in degrees (30 km a degree)')
        s.add_argument('--from-year', type=float, default=FROM_YEAR)
        s.add_argument('--years', type=float, default=1.0)
        s.set_defaults(fn=fn)
    args = p.parse_args(argv)
    return args.fn(args)


if __name__ == '__main__':
    sys.exit(main())
