"""Carry the improved crystal component's Earth errors through the lunar snow column.

The data-sized check (data_sized_check.py) measures, at each of the ten
laboratory temperatures, how far the improved McSnow component's mass, axes
and fall speed sit from Takahashi et al.'s (1991) crystals.  This module turns
each error into a multiplier on the component's own physics:

* mass growth: the capacitance, and with it the vapour flux to the crystal;
* axis allocation: the ventilated inherent growth factor, which divides new
  ice between the polar and equatorial axes;
* fall speed: the terminal speed, which sets settling and ventilation together.

At every laboratory temperature each multiplier is sized so that the perturbed
Earth crystal removes the component's median error there.  The band's edges are
the smallest and largest of those multipliers.  The lunar column and the
archived lunar storm paths are then rerun at each edge.  With every multiplier
at one the perturbed builds reproduce the frozen component exactly; that parity
is checked.  No model parameter is fitted to laboratory data: the multipliers
measure the size of the component's errors, and the component stays frozen.
"""
from __future__ import annotations

import argparse
import contextlib
import ctypes
import difflib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import time

import numpy as np

from atmosphere.snow import data_sized_check as check
from atmosphere.snow import mcsnow_lunar_growth as lunar
from atmosphere.snow import mcsnow_weather as weather
from atmosphere.snow import takahashi_growth as lab

model, native = lunar.model, lunar.native
ROOT = lunar.ROOT
HERE = Path(__file__).resolve().parent
CACHE = ROOT / 'research/runs/snow_agent/earth_error_band'
PRODUCT = HERE / 'results/earth_error_band.json'
FACTORS = ('mass', 'axes', 'speed')
LAB_PRESSURE = 101000.
SEED_DIAMETER = 4e-6


# --- source patches ----------------------------------------------------------------

def depo_source(world):
    """The improved component's deposition source with three band multipliers."""
    original, changed = lunar.patches('plate_base', world)
    body = changed['mo_depo.f90']
    rep = model.replace_once
    body = rep(body, '  REAL(wp), PUBLIC :: supply_rate=-1._wp, requested_rate=0._wp',
               '  REAL(wp), PUBLIC :: supply_rate=-1._wp, requested_rate=0._wp\n'
               '  REAL(wp), PUBLIC :: band_factor(3)=1._wp')
    start = body.index('SUBROUTINE sp_depo( sp,')
    end = body.index('END SUBROUTINE sp_depo', start)
    part = rep(body[start:end], '    vt = vterm_snow( atmo, sp )  ! modify heat and vapour transfer by ventilation',
               '    vt = vterm_snow( atmo, sp ) * band_factor(3)  ! modify heat and vapour transfer by ventilation')
    body = body[:start] + part + body[end:]
    body = rep(body, '  C = capacitance( sp, Re )\n  growth_trace=0._wp',
               '  C = capacitance( sp, Re ) * band_factor(1)\n  growth_trace=0._wp')
    body = rep(body, '  growth_trace(14)=IGF_vent', '  IGF_vent=IGF_vent*band_factor(2)\n  growth_trace(14)=IGF_vent')
    sources = dict(changed)
    sources['mo_depo.f90'] = body
    return original, sources


def column_driver():
    text = lunar.DRIVER.read_text()
    rep = model.replace_once
    text = rep(text, '  use mo_depo, only: sp_depo, growth_trace, supply_rate, requested_rate\n',
               '  use mo_depo, only: sp_depo, growth_trace, supply_rate, requested_rate, band_factor\n')
    text = rep(text, '  real(wp) :: last_requested, last_admitted, budget\n',
               '  real(wp) :: last_requested, last_admitted, budget\n  real(wp) :: band_mass, band_axes, band_speed\n')
    text = rep(text, '    max_mass_fraction, diameter_stop_m, growth_enabled\n',
               '    max_mass_fraction, diameter_stop_m, growth_enabled, band_mass, band_axes, band_speed\n')
    text = rep(text, '  max_mass_fraction=0._wp; diameter_stop_m=1._wp; growth_enabled=.true.; native_saturation=.true.\n',
               '  max_mass_fraction=0._wp; diameter_stop_m=1._wp; growth_enabled=.true.; native_saturation=.true.\n'
               '  band_mass=1._wp; band_axes=1._wp; band_speed=1._wp\n')
    text = rep(text, '  IGF_id=growth_mode\n',
               "  if (min(band_mass,band_axes,band_speed)<=0) error stop 'Invalid band factor'\n"
               '  band_factor=[band_mass,band_axes,band_speed]\n  IGF_id=growth_mode\n')
    text = rep(text, '    speed=vterm_snow(air,particle)\n', '    speed=vterm_snow(air,particle)*band_factor(3)\n')
    text = rep(text, 'particle%p_area,vterm_snow(air,particle), &', 'particle%p_area,vterm_snow(air,particle)*band_factor(3), &')
    return text


def weather_driver():
    text = weather.DRIVER.read_text()
    rep = model.replace_once
    text = rep(text, '  use mo_depo, only: sp_depo, supply_rate\n', '  use mo_depo, only: sp_depo, supply_rate, band_factor\n')
    text = rep(text, '  subroutine weather_constants(values) bind(C)\n',
               '  subroutine band_set(values) bind(C)\n'
               '    real(c_double),intent(in) :: values(3)\n'
               "    if (minval(values)<=0) error stop 'Invalid band factor'\n"
               '    band_factor=values\n'
               '  end subroutine\n\n'
               '  subroutine weather_constants(values) bind(C)\n')
    text = rep(text, 'sp%m_i/vol,vterm_snow(air,sp), &', 'sp%m_i/vol,vterm_snow(air,sp)*band_factor(3), &')
    return text


DRIVERS = {'column': ('band_column_driver.f90', column_driver), 'weather': ('band_weather_driver.f90', weather_driver)}


def build(kind, world):
    path = CACHE / f'build_{kind}_{world}'
    if path.exists():
        return verify_build(kind, world)
    baseline = native.verify_build(model.BASELINE)
    original, sources = depo_source(world)
    path.mkdir(parents=True)
    shutil.copytree(model.BASELINE / 'source', path / 'source')
    diff = ''
    for name, body in sources.items():
        (path / 'source/src' / name).write_text(body)
        diff += ''.join(difflib.unified_diff(original[name].splitlines(True), body.splitlines(True),
                                             fromfile='native/' + name, tofile=f'band_{world}/' + name))
    driver_name, driver = DRIVERS[kind]
    (path / driver_name).write_text(driver())
    (path / 'band.patch').write_text(diff)
    for name in native.TABLE_NAMES:
        shutil.copyfile(model.BASELINE / name, path / name)
    compiler = baseline['compiler_path']
    flags = list(native.FLAGS) if kind == 'column' else [*native.FLAGS, '-fPIC', '-ffree-line-length-none']
    selected = native.selected_sources(path / 'source/src', path / driver_name, compiler)
    objects, commands = [], []
    target = path / ('crystal' if kind == 'column' else 'libband.so')
    with (path / 'build.log').open('w') as log:
        for i, src in enumerate(selected):
            obj = path / f'{i:02d}_{src.stem}.o'
            objects.append(str(obj))
            cmd = [compiler, *flags, '-I' + str(path / 'source/src'), '-I' + str(path), '-J' + str(path),
                   '-c', str(src), '-o', str(obj)]
            commands.append(cmd)
            subprocess.run(cmd, cwd=path, env=native.ENV, stdout=log, stderr=subprocess.STDOUT, check=True)
        cmd = [compiler, *flags, *([] if kind == 'column' else ['-shared']), *objects, '-o', str(target)]
        commands.append(cmd)
        subprocess.run(cmd, cwd=path, env=native.ENV, stdout=log, stderr=subprocess.STDOUT, check=True)
    record = dict(schema='terluna.atmosphere.snow-band-build/1', kind=kind, world=world,
                  baseline_record_sha256=native.sha(model.BASELINE / 'build.json'),
                  source_identity=baseline['source_identity'], table_identity=baseline['table_identity'],
                  files={str(p.relative_to(path)): native.sha(p) for p in sorted(path.rglob('*'))
                         if p.is_file() and p.suffix not in ('.o', '.mod')},
                  compiler_path=compiler, compiler_sha256=native.sha(compiler),
                  compiler_version=baseline['compiler_version'], flags=flags, commands=commands)
    model.save(path / 'build.json', record)
    return record


def verify_build(kind, world):
    path = CACHE / f'build_{kind}_{world}'
    record = json.loads((path / 'build.json').read_text())
    native.verify_build(model.BASELINE)
    if native.sha(model.BASELINE / 'build.json') != record['baseline_record_sha256']:
        raise ValueError('Baseline changed')
    for name, digest in record['files'].items():
        if native.sha(path / name) != digest:
            raise ValueError('Build artifact changed: ' + name)
    _, sources = depo_source(world)
    if any((path / 'source/src' / n).read_text() != b for n, b in sources.items()):
        raise ValueError('Build source differs from its declared patch')
    driver_name, driver = DRIVERS[kind]
    if (path / driver_name).read_text() != driver():
        raise ValueError('Driver differs from its declared patch')
    return record


# --- the storm-path library -----------------------------------------------------------

class BandLibrary:
    """The weather driver's interface over the band build; serial use only."""
    _initialized = set()

    def __init__(self, world='moon'):
        self.identity = verify_build('weather', world)
        path = CACHE / f'build_weather_{world}'
        self.lib = ctypes.CDLL(str(path / 'libband.so'))
        fp = np.ctypeslib.ndpointer(dtype=np.float64, flags='F_CONTIGUOUS')
        ip = np.ctypeslib.ndpointer(dtype=np.int32, flags='C_CONTIGUOUS')
        self.lib.weather_evolve.argtypes = [ctypes.c_int, fp, fp, fp, ctypes.c_int, ctypes.c_double,
                                            ctypes.c_double, ctypes.c_double, ctypes.c_double, fp, ip]
        self.lib.weather_evolve.restype = None
        self.lib.weather_constants.argtypes = [fp]
        self.lib.weather_constants.restype = None
        self.lib.weather_init.argtypes = []
        self.lib.weather_init.restype = None
        self.lib.band_set.argtypes = [np.ctypeslib.ndpointer(dtype=np.float64)]
        self.lib.band_set.restype = None
        old = Path.cwd()
        try:
            if str(path) not in self._initialized:
                os.chdir(path)
                self.lib.weather_init()
                self._initialized.add(str(path))
        finally:
            os.chdir(old)
        values = np.empty(4)
        self.lib.weather_constants(values)
        self.ice_density, self.Rd, self.Rv, self.gravity = values
        self.factors = (1., 1., 1.)
        self.set_band(1., 1., 1.)

    def set_band(self, mass=1., axes=1., speed=1.):
        values = np.array([mass, axes, speed], dtype=float)
        if not np.all(np.isfinite(values)) or np.any(values <= 0):
            raise ValueError('Band factors must be positive')
        self.lib.band_set(values)
        self.factors = tuple(values.tolist())

    # The two methods the storm integrator calls, unchanged in their arithmetic.
    seed = weather.Library.seed
    advance = weather.Library.advance


# --- Earth sizing of the band edges ------------------------------------------------------

def earth_history(lib, temperature, ages_min):
    """Constant laboratory air at water saturation from the 4 micrometre sphere."""
    T = 273.15 + temperature
    force = np.asfortranarray(np.array([[T], [LAB_PRESSURE], [100.]]))
    probe = lib.seed(force, SEED_DIAMETER)
    diag, _ = lib.advance(probe, force, 0.)
    force[2, 0] = diag[10, 0]  # native water saturation, as in the retained traces
    state = lib.seed(force, SEED_DIAMETER)
    out, elapsed = {}, 0.
    for age in sorted(ages_min):
        # A trial factor far from one can reach the 1 m censor; its frozen state then
        # reads as too large, which is the sign the bisection needs.  Native failures
        # (negative codes) still raise inside advance().
        diag, code = lib.advance(state, force, age * 60 - elapsed, saturated=True)
        elapsed = age * 60
        out[age] = dict(mass_kg=float(state[0, 0]), equatorial_dimension_m=float(diag[1, 0]),
                        polar_dimension_m=float(diag[2, 0]), terminal_speed_m_s=float(diag[4, 0]))
    return out


def crystals_at(temperature, items):
    return [c for c in items if c['nominal_C'] == temperature]


def median_log_error(lib, factor, value, temperature, items, quantity):
    kwargs = {f: 1. for f in FACTORS}
    kwargs[factor] = value
    lib.set_band(**kwargs)
    members = crystals_at(temperature, items)
    if quantity == 'terminal_speed_m_s':
        ages = sorted({c['age_min'] for c in members if check.speed_defined(temperature, c['age_min'])})
        hist = earth_history(lib, temperature, ages)
        logs = [math.log(hist[a]['terminal_speed_m_s'] / lab.fall_speed(temperature, a)) for a in ages]
    else:
        ages = sorted({c['age_min'] for c in members})
        hist = earth_history(lib, temperature, ages)
        logs = []
        for c in members:
            m = c['measured']
            if quantity == 'aspect_ratio':
                if 'polar_dimension_m' in m and 'equatorial_dimension_m' in m:
                    h = hist[c['age_min']]
                    logs.append(math.log((h['polar_dimension_m'] / h['equatorial_dimension_m'])
                                         / (m['polar_dimension_m'] / m['equatorial_dimension_m'])))
            elif quantity in m:
                logs.append(math.log(hist[c['age_min']][quantity] / m[quantity]))
    lib.set_band()
    return float(np.median(logs)) if logs else None


def size_factor(lib, factor, temperature, items, quantity, lo=.1, hi=10., tol=1e-4):
    """Bisection in log factor for a zero median log error; monotone by construction of each knob."""
    f_lo, f_hi = median_log_error(lib, factor, lo, temperature, items, quantity), \
        median_log_error(lib, factor, hi, temperature, items, quantity)
    if f_lo is None:
        return None
    if f_lo * f_hi > 0:
        return dict(factor=None, reason='no sign change within bounds', bounds=[lo, hi], errors=[f_lo, f_hi])
    a, b, fa = math.log(lo), math.log(hi), f_lo
    evaluations = 2
    while b - a > tol:
        m = .5 * (a + b)
        fm = median_log_error(lib, factor, math.exp(m), temperature, items, quantity)
        evaluations += 1
        if (fm > 0) == (fa > 0):
            a, fa = m, fm
        else:
            b = m
    return dict(factor=math.exp(.5 * (a + b)), evaluations=evaluations,
                unperturbed_median_log_error=median_log_error(lib, factor, 1., temperature, items, quantity))


QUANTITY_OF = {'mass': 'mass_kg', 'axes': 'aspect_ratio', 'speed': 'terminal_speed_m_s'}


def size_band():
    items = check.crystals()
    lib = BandLibrary('earth')
    sizing = {}
    for factor in FACTORS:
        rows = {}
        for t in check.NOMINAL:
            r = size_factor(lib, factor, t, items, QUANTITY_OF[factor])
            if r is not None:
                rows[format(t, 'g')] = r
        found = {t: r['factor'] for t, r in rows.items() if r.get('factor')}
        lo_t = min(found, key=found.get)
        hi_t = max(found, key=found.get)
        sizing[factor] = dict(quantity=QUANTITY_OF[factor], by_temperature=rows,
                              low=dict(factor=found[lo_t], nominal_C=float(lo_t)),
                              high=dict(factor=found[hi_t], nominal_C=float(hi_t)))
    return sizing


def earth_consistency():
    """The Earth band library at factor one against the retained component traces."""
    items = check.crystals()
    lib = BandLibrary('earth')
    manifest = {m['name']: m for m in json.loads(check.SHAPE_PRODUCT.read_text())['manifest']}
    worst = 0.
    for t in check.NOMINAL:
        ages = sorted({c['age_min'] for c in crystals_at(t, items)})
        trace = check.component_trace(t, check.COMPONENT_DT, manifest)
        hist = earth_history(lib, t, ages)
        for a in ages:
            ref = check.at_age(trace, a)
            if ref is None:
                continue
            for q in ('mass_kg', 'equatorial_dimension_m', 'polar_dimension_m', 'terminal_speed_m_s'):
                worst = max(worst, abs(hist[a][q] / ref[q] - 1))
    return dict(maximum_relative_difference=worst,
                note='Adaptive weather-driver steps against the fixed 1/128 s retained traces.')


# --- the lunar column -------------------------------------------------------------------

COLUMN_CASES = {
    'allowance_F': 'fall_main_moon_deep_plate_base_1_dt2.0',
    'allowance_0.1F': 'flux_0.1_dt2.0',
    'allowance_0.01F': 'flux_0.01_dt2.0',
    'unlimited_saturation': 'fall_main_moon_deep_plate_base_0_dt2.0',
}


def column_entries():
    entries = {e['name']: e for e in lunar.cases()}
    return {label: entries[name] for label, name in COLUMN_CASES.items()}


def run_column(label, entry, factors, dt=None):
    from shared.constants import MOON_SURFACE_GRAVITY
    verify_build('column', 'moon')
    buildpath = CACHE / 'build_column_moon'
    config = dict(entry['configuration'], band_mass=factors[0], band_axes=factors[1], band_speed=factors[2])
    if dt is not None:
        config.update(dt_s=dt, max_mass_fraction=.005)
    path = CACHE / 'column_runs' / label
    if path.exists():
        record = json.loads((path / 'run.json').read_text())
        if record['configuration'] != config:
            raise ValueError('Existing column run has different inputs: ' + label)
        return record
    path.mkdir(parents=True)
    (path / 'crystal.nml').write_text(native.namelist(config))
    for name in native.TABLE_NAMES:
        shutil.copyfile(buildpath / name, path / name)
    start = time.monotonic()
    with (path / 'run.log').open('w') as log:
        proc = subprocess.run([str(buildpath / 'crystal'), 'crystal.nml'], cwd=path, env=native.ENV,
                              stdout=log, stderr=subprocess.STDOUT)
    analysis = lunar.analyze(path, config, proc.returncode, MOON_SURFACE_GRAVITY)
    final = analysis['snapshots'][-1] if analysis['snapshots'] else {}
    keep = ('elapsed_s', 'height_above_freezing_m', 'mass_kg', 'equatorial_dimension_m', 'polar_dimension_m',
            'diameter_m', 'density_kg_m3', 'aspect_ratio', 'terminal_speed_m_s', 'integrated_deposition_kg',
            'allocated_water_kg')
    record = dict(label=label, case=entry['name'], configuration=config, returncode=proc.returncode,
                  elapsed_wall_s=time.monotonic() - start, complete=analysis['complete'],
                  validation_errors=analysis['validation_errors'], summary=analysis['summary'],
                  final={k: final[k] for k in keep if k in final},
                  file_sha256={p.name: native.sha(p) for p in sorted(path.iterdir()) if p.is_file()})
    model.save(path / 'run.json', record)
    return record


def column_parity():
    """With every factor at one, the band column must equal the retained run byte for byte."""
    rows = []
    for label in ('allowance_F', 'unlimited_saturation'):
        entry = column_entries()[label]
        record = run_column('parity_' + label, entry, (1., 1., 1.))
        original = lunar.CACHE / 'runs' / entry['name'] / 'trajectory.csv'
        same = (CACHE / 'column_runs' / ('parity_' + label) / 'trajectory.csv').read_bytes() == original.read_bytes()
        rows.append(dict(case=entry['name'], original_sha256=native.sha(original), identical=same,
                         complete=record['complete']))
    return rows


# --- the archived lunar storm -----------------------------------------------------------

def run_weather(label, factors, saturated, step_s=20., diameter=SEED_DIAMETER):
    path = CACHE / 'weather_runs' / label
    if path.exists():
        record = json.loads((path / 'run.json').read_text())
        if record['factors'] != list(factors) or record['saturated'] != saturated or record['step_s'] != step_s:
            raise ValueError('Existing storm run has different inputs: ' + label)
        return record
    identity = weather.verify_weather()
    data = np.load(weather.CACHE / 'weather/fields.npy', mmap_mode='r')
    lib = BandLibrary('moon')
    lib.set_band(*factors)
    with contextlib.redirect_stdout(open(os.devnull, 'w')):
        result = weather.integrate(lib, data, identity, diameter=diameter, step_s=step_s, saturated=saturated)
    lib.set_band()
    path.mkdir(parents=True)
    arrays = {k: v for k, v in result.items() if isinstance(v, np.ndarray)}
    np.savez_compressed(path / 'raw.npz', **arrays)
    status = result['status']
    peak = result['peak_diameter_m']
    record = dict(label=label, factors=list(factors), saturated=saturated, step_s=step_s, diameter_m=diameter,
                  wall_s=result['wall_s'], releases=int(len(peak)),
                  terminations={text: int(np.sum(status == code)) for code, text in weather.STATUSES.items() if code},
                  peak_median_m=float(np.median(peak)), peak_p90_m=float(np.percentile(peak, 90)),
                  peak_stats_m=weather.history.stats(peak), raw_sha256=native.sha(path / 'raw.npz'),
                  weather_record_sha256=native.sha(weather.CACHE / 'weather/input.json'))
    model.save(path / 'run.json', record)
    return record


def weather_parity():
    """With every factor at one, the band storm run must equal the retained archived-vapour run."""
    record = run_weather('parity_archived', (1., 1., 1.), 0)
    original = np.load(weather.CACHE / 'runs/actual_small_step20/raw.npz')
    mine = np.load(CACHE / 'weather_runs/parity_archived/raw.npz')
    names = ('peak_diameter_m', 'status', 'states', 'age_s', 'deposited_mass_kg', 'final_position_m')
    return dict(original=str((weather.CACHE / 'runs/actual_small_step20/raw.npz').relative_to(ROOT)),
                original_sha256=native.sha(weather.CACHE / 'runs/actual_small_step20/raw.npz'),
                identical={n: bool(np.array_equal(original[n], mine[n])) for n in names},
                complete=record['releases'] == 354)


# --- the experiment ---------------------------------------------------------------------

def edges(sizing):
    """Single-factor edges, each with the other two factors at one."""
    out = {}
    for i, f in enumerate(FACTORS):
        for side in ('low', 'high'):
            factors = [1., 1., 1.]
            factors[i] = sizing[f][side]['factor']
            out[f'{f}_{side}'] = tuple(factors)
    return out


def combined(single, sizing, size_of):
    """Join the three single-factor edges that most reduce and most increase size."""
    small, large = [1., 1., 1.], [1., 1., 1.]
    for i, f in enumerate(FACTORS):
        lo, hi = size_of(f + '_low'), size_of(f + '_high')
        a, b = sizing[f]['low']['factor'], sizing[f]['high']['factor']
        small[i], large[i] = (a, b) if lo <= hi else (b, a)
    return tuple(small), tuple(large)


def experiment():
    for kind, world in (('column', 'moon'), ('weather', 'moon'), ('weather', 'earth')):
        build(kind, world)
    sizing_path = CACHE / 'sizing.json'
    if sizing_path.exists():
        sizing = json.loads(sizing_path.read_text())
    else:
        sizing = dict(band=size_band(), earth_consistency=earth_consistency())
        model.save(sizing_path, sizing)
    band = sizing['band']
    single = edges(band)
    parity = dict(column=column_parity(), weather=weather_parity())
    columns = {}
    for case, entry in column_entries().items():
        runs = {'baseline': run_column(f'{case}__baseline', entry, (1., 1., 1.))}
        for name, factors in single.items():
            runs[name] = run_column(f'{case}__{name}', entry, factors)
        size = lambda name: runs[name]['final']['diameter_m']
        small, large = combined(single, band, size)
        runs['combined_small'] = run_column(f'{case}__combined_small', entry, small)
        runs['combined_large'] = run_column(f'{case}__combined_large', entry, large)
        for name in ('baseline', 'combined_small', 'combined_large'):
            factors = runs[name]['configuration']
            runs[name + '_dt0.5'] = run_column(f'{case}__{name}_dt0.5', entry,
                                               (factors['band_mass'], factors['band_axes'], factors['band_speed']), dt=.5)
        columns[case] = runs
    storms = {}
    for humidity, saturated in (('archived', 0), ('saturated', 1)):
        runs = {'baseline': run_weather(f'{humidity}__baseline', (1., 1., 1.), saturated)}
        for name, factors in single.items():
            runs[name] = run_weather(f'{humidity}__{name}', factors, saturated)
        size = lambda name: runs[name]['peak_median_m']
        small, large = combined(single, band, size)
        runs['combined_small'] = run_weather(f'{humidity}__combined_small', small, saturated)
        runs['combined_large'] = run_weather(f'{humidity}__combined_large', large, saturated)
        storms[humidity] = runs
    return sizing, parity, columns, storms


def assemble():
    sizing, parity, columns, storms = experiment()
    column_table = {}
    for case, runs in columns.items():
        column_table[case] = {name: dict(factors=[r['configuration']['band_mass'], r['configuration']['band_axes'],
                                                  r['configuration']['band_speed']],
                                         diameter_m=r['final']['diameter_m'], mass_kg=r['final']['mass_kg'],
                                         equatorial_dimension_m=r['final']['equatorial_dimension_m'],
                                         polar_dimension_m=r['final']['polar_dimension_m'],
                                         descent_h=r['final']['elapsed_s'] / 3600,
                                         status=r['summary']['status_code'] if r['summary'] else None,
                                         complete=r['complete'])
                              for name, r in runs.items()}
    storm_table = {h: {name: dict(factors=r['factors'], peak_median_m=r['peak_median_m'], peak_p90_m=r['peak_p90_m'],
                                  small_residual=r['terminations'].get('small_ice_residual'),
                                  diameter_censored=r['terminations'].get('diameter_censored'))
                       for name, r in runs.items()} for h, runs in storms.items()}

    def spread(values):
        return dict(minimum=min(values), maximum=max(values), ratio=max(values) / min(values))
    band_spread = {case: spread([v['diameter_m'] for k, v in t.items() if not k.endswith('dt0.5')])
                   for case, t in column_table.items()}
    band_spread.update({'storm_' + h: spread([v['peak_median_m'] for v in t.values()]) for h, t in storm_table.items()})
    moisture_spread = dict(
        column=spread([t['baseline']['diameter_m'] for t in column_table.values()]),
        storm=spread([t['baseline']['peak_median_m'] for t in storm_table.values()]),
        all_baselines=spread([t['baseline']['diameter_m'] for t in column_table.values()]
                             + [t['baseline']['peak_median_m'] for t in storm_table.values()]))
    refinement = max(abs(column_table[c][n + '_dt0.5']['diameter_m'] / column_table[c][n]['diameter_m'] - 1)
                     for c in column_table for n in ('baseline', 'combined_small', 'combined_large'))
    product = dict(
        schema='terluna.atmosphere.snow-earth-error-band/1',
        evidence=('Numerical sensitivity of the frozen improved McSnow component in prescribed lunar conditions, with '
                  'perturbations sized by its measured Earth errors; no population, aggregation or empirical lunar '
                  'validation.'),
        method=dict(factors={'mass': 'capacitance multiplier (vapour flux)', 'axes': 'ventilated inherent growth factor multiplier',
                             'speed': 'terminal speed multiplier (settling and ventilation)'},
                    sizing=('At each laboratory temperature, the factor at which the perturbed Earth crystal (4 micrometre '
                            'sphere, 101 kPa, water saturation) has zero median log error against the measured crystals: '
                            'mass against crystal masses, axes against crystal aspect ratios where both axes are measured, '
                            'speed against the Table 4 relation at the crystals\' ages. Edges are the smallest and largest '
                            'factors over temperature.'),
                    edges='Each edge alone; then the combination of edges that most reduces and most increases size.'),
        band=sizing['band'], earth_consistency=sizing['earth_consistency'], parity=parity,
        column=column_table, storm=storm_table, spread=dict(band=band_spread, moisture=moisture_spread),
        numerical_checks=dict(column_dt_2_vs_0p5=refinement),
        inputs=dict(data_sized_check_sha256=native.sha(check.PRODUCT), measurements_sha256=native.sha(check.MEASUREMENTS),
                    lunar_growth_product_sha256=native.sha(lunar.PRODUCT), weather_product_sha256=native.sha(weather.PRODUCT)),
        builds={f'{k}_{w}': native.sha(CACHE / f'build_{k}_{w}/build.json')
                for k, w in (('column', 'moon'), ('weather', 'moon'), ('weather', 'earth'))},
        producer_sha256=native.sha(Path(__file__)), raw_directory=str(CACHE.relative_to(ROOT)),
        limits=['The multipliers act at every temperature at once; the component\'s errors differ in sign and size '
                'between temperatures, so each edge is a bound on one kind of error rather than a corrected model.',
                'Laboratory errors are measured over at most 30 minutes at Earth gravity and pressure; lunar crystals grow '
                'for up to 54 hours at a sixth of the gravity.',
                'Single crystals in prescribed air: no population depletion, latent heating or collisions.',
                'The archived storm is the older A28_dim5 case, not the 301 K design climate.'])
    PRODUCT.write_text(json.dumps(product, separators=(',', ':'), allow_nan=False) + '\n')
    return product


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('build', 'size', 'parity', 'run'))
    args = parser.parse_args()
    if args.action == 'build':
        for kind, world in (('column', 'moon'), ('weather', 'moon'), ('weather', 'earth')):
            print(kind, world, build(kind, world)['files'].get('crystal') or 'built')
    elif args.action == 'size':
        print(json.dumps(dict(consistency=earth_consistency(), band=size_band()), indent=1))
    elif args.action == 'parity':
        print(json.dumps(dict(column=column_parity(), weather=weather_parity()), indent=1))
    else:
        product = assemble()
        print(json.dumps(dict(band={f: {s: product['band'][f][s] for s in ('low', 'high')} for f in FACTORS},
                              spread=product['spread'], parity=product['parity'],
                              refinement=product['numerical_checks']), indent=1))


if __name__ == '__main__':
    main()
