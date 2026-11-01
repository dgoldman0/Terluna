"""The storm-water study, phase 1: how much of the drylands' storm rain runs off dry ground, on four designed surfaces,
from the design years' 3-day means, a storm library and Green and Ampt's infiltration at lunar gravity.

    climate/gcm/.venv/bin/python -m research.studies.storm_water.run [--storms PRODUCT] [--extract]   # a few minutes

Each dryland cell's 3-day mean of convective rain becomes a Poisson number of point storms from a library (storms.py);
each storm infiltrates into the soil the bucket's water implies at that mean (infiltration.py), and what the surface
cannot take or hold runs off. The large-scale rain falls as one steady rain of LS_RATE_MM_H. Two libraries bracket the
storms' intensities (Walnut Gulch's measured storms above, CM1's 6 km window below) and two treatments bracket the
wetting between storms (each storm on the mean's soil, or one wetting front carried through the mean's storms). The
same chain at Earth's gravity on Walnut Gulch-like inputs is checked against Walnut Gulch's measured small-watershed
runoff. results/storm_water.json is the product; --storms takes any CM1 storm-lives product in place of the 6 km window.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np

from research.studies.storm_water import extract as EX
from research.studies.storm_water import infiltration as IN
from research.studies.storm_water import storms as ST
from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, SYNODIC_MONTH_DAYS
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OUT = HERE / 'results' / 'storm_water.json'
SCHEMA = 'terluna.research.storm-water-runoff/2'
STORMS = 'climate/results/crm/storms_box_0e_elec_corrected_storms_first.json'
CODE = ('research/studies/storm_water/run.py', 'research/studies/storm_water/storms.py',
        'research/studies/storm_water/infiltration.py', 'research/studies/storm_water/extract.py',
        'research/studies/storm_water/wetting.py', 'research/studies/storm_water/routing.py',
        'research/studies/storm_water/phase2a.py', 'geography/drainage.py', 'climate/crm/gcm_site.py')
MEAN_DAYS = SYNODIC_MONTH_DAYS / 10.0                    # each 3-day mean is a tenth of the lunar day
MM_PER_MEAN = 1000.0 * MEAN_DAYS * JULIAN_DAY            # m/s of water to mm per mean
SATURATIONS = np.array([0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.65, 0.8, 0.95])
LS_RATE_MM_H = 2.0                                       # estimate: the large-scale rain's steady rate
SURFACES = IN.MAIN + ('sealed_loam_forming',)
BRACKETS = ('walnut_gulch', 'cm1_6km')
MODES = ('independent', 'continuing')
CLASSES = ('arid', 'semi-arid', 'dry sub-humid')
MEASURED_RUNOFF = (0.03, 0.10)       # Walnut Gulch's small watersheds, 0.18-5.42 ha (Stone et al. 2008, via the plan)
WG_YEAR_MM, WG_SEASON_DAYS = 305.0, 80.0                 # Mendez et al. 2003; Osborn 1983 (July to early September)
THRESHOLDS = dict(yearly=0.05, tenth=0.10)               # the plan's leverage for moving to PlaSim's land scheme


def sha256(path) -> str:
    h = hashlib.sha256()
    with open(ROOT / path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 24), b''):
            h.update(block)
    return h.hexdigest()


def r3(v):
    """Three significant figures, for the product."""
    if isinstance(v, (list, tuple, np.ndarray)):
        return [r3(u) for u in v]
    if v is None or not np.isfinite(v):
        return None
    return float(f'{float(v):.3g}')


class Tables:
    """Each storm's runoff on fresh soil, and the carried-front sequences, for one surface, gravity and library."""

    def __init__(self, p, lib, n_max):
        self.p, self.lib = p, lib
        w, d = np.asarray(lib['weights'], float), ST.depths(lib)
        self.mean_depth = float((w * d).sum() / w.sum())
        r1, _ = IN.single_storms(p, lib, SATURATIONS)                         # (storm, saturation)
        self.ratio = (w[:, None] * r1).sum(axis=0) / (w * d).sum()
        self.ponding = (w[:, None] * (r1 > 0.05)).sum(axis=0) / w.sum()
        self.seq = IN.sequences(p, lib, SATURATIONS, n_max) if n_max else None


def interp(row, s):
    s = np.clip(s, SATURATIONS[0], SATURATIONS[-1])
    return np.interp(s, SATURATIONS, row)


def poisson_mix(table, lam, s):
    """The expectation over a Poisson number of storms (mean lam) of a table's row at saturation s."""
    n_max = table.shape[0] - 1
    top = float(lam.max()) if lam.size else 0.0
    cap = max(n_max, int(top + 10.0 * math.sqrt(top) + 10))
    p = np.exp(-lam)
    acc = p * interp(table[0], s)
    slope = table[n_max] - table[n_max - 1]
    for n in range(1, cap + 1):
        p = p * lam / n
        row = table[n] if n <= n_max else table[n_max] + (n - n_max) * slope
        acc = acc + p * interp(row, s)
    return acc


def runoff(x, tab, mode, mapping='plant_available', ls_rate=LS_RATE_MM_H, cells=slice(None), ls_cache=None):
    """Convective and large-scale runoff (mm per mean) for (mean, cell)."""
    p = tab.p
    conv_mm = x['prc'][:, cells].astype(float) * MM_PER_MEAN
    ls_mm = np.maximum(x['pr'][:, cells].astype(float) * MM_PER_MEAN - conv_mm, 0.0)
    s = IN.initial_saturation(p['soil'], x['mrso'][:, cells].astype(float), mapping)
    if mode == 'independent':
        conv = conv_mm * interp(tab.ratio, s)
        f0 = np.zeros_like(conv)
    else:
        lam = conv_mm / tab.mean_depth
        conv = poisson_mix(tab.seq[0], lam, s)
        f0 = poisson_mix(tab.seq[1], lam, s)
    key = (p['name'], p['gravity'], p['retention'], p['ks'], p['kc'], mapping, ls_rate, mode, id(tab) if mode != 'independent' else 0)
    if ls_cache is not None and key in ls_cache:
        return conv, ls_cache[key]
    ls = np.zeros_like(conv)
    sel = ls_mm > 0.02
    if ls_rate > p['k_sub'] and sel.any():
        ls[sel] = IN.steady_rain(p, ls_mm[sel], ls_rate, IN.deficit(p, s[sel]), f0[sel])
    if ls_cache is not None:
        ls_cache[key] = ls
    return conv, ls


def by_tenth(values, tenth):
    return np.stack([(values * (tenth == k)).sum(axis=0) for k in range(10)])        # (tenth, cell)


def summarise(x, off, cells=slice(None)):
    """Per cell: yearly runoff ratio and by tenth; per class: area-weighted ratios, volumes and leverage."""
    rain = x['pr'][:, cells].astype(float) * MM_PER_MEAN
    evap = -x['evap'][:, cells].astype(float) * MM_PER_MEAN
    tenth = x['tenth'][:, cells]
    ot, rt, et = by_tenth(off, tenth), by_tenth(rain, tenth), by_tenth(evap, tenth)
    with np.errstate(divide='ignore', invalid='ignore'):
        cell_year = off.sum(axis=0) / rain.sum(axis=0)
        cell_tenth = np.where(rt > 1e-3 * rain.sum(axis=0), ot / rt, np.nan)
    area, cls = x['area_km2'][cells], x['cls'][cells]
    classes = {}
    years = rain.shape[0] / 120.0
    for name in CLASSES + ('drylands',):
        m = np.ones(area.shape, bool) if name == 'drylands' else (cls == name)
        if not m.any():
            continue
        a = area[m]
        o, r, e = (off.sum(axis=0)[m] * a).sum(), (rain.sum(axis=0)[m] * a).sum(), (evap.sum(axis=0)[m] * a).sum()
        lev_t = (ot[:, m] * a).sum(axis=1) / np.maximum((et[:, m] * a).sum(axis=1), 1e-30)
        classes[name] = dict(runoff_ratio=o / r, runoff_mm_yr=o / a.sum() / years * JULIAN_YEAR_DAYS / (12 * SYNODIC_MONTH_DAYS),
                             runoff_km3_yr=o * 1e-6 / years * JULIAN_YEAR_DAYS / (12 * SYNODIC_MONTH_DAYS),
                             share_of_evaporation=o / e, largest_tenth_share_of_evaporation=float(lev_t.max()),
                             runoff_ratio_by_tenth=list((ot[:, m] * a).sum(axis=1) / np.maximum((rt[:, m] * a).sum(axis=1), 1e-30)))
    return cell_year, cell_tenth, classes


def libraries(storms_path):
    return dict(walnut_gulch=ST.walnut_gulch(), cm1_6km=ST.from_storm_lives(storms_path, ROOT))


def mean_pressure_pa(x):
    w = x['area_km2'] * (x['pr'].mean(axis=0))
    return float((x['ps'].mean(axis=0) * w).sum() / w.sum()) * 100.0


def n_needed(x, lib):
    lam = x['prc'].astype(float).max() * MM_PER_MEAN / ST.statistics(lib)['mean_depth_mm']
    return int(min(96, lam + 6.0 * math.sqrt(lam) + 6))


def walnut_gulch_check(libs, cache):
    """The chain at Earth's gravity on Walnut Gulch-like inputs: its thunderstorm rain (two thirds of 305 mm) as the
    library's storms on Rawls' sandy loam and loam at the wilting point, the winter's rain too gentle to run off."""
    out = {}
    share = ST.SEASON_MM / WG_YEAR_MM
    per_window = ST.SEASON_MM / (WG_SEASON_DAYS / MEAN_DAYS)
    for bracket in BRACKETS:
        lib = libs[bracket]
        row = {}
        for soil in ('sandy_loam', 'loam'):
            p = IN.parameters(soil, gravity=1.0)
            tab = cache(p, lib, 12)
            s0 = float(IN.initial_saturation(p['soil'], 0.0))
            lam = np.array([per_window / tab.mean_depth])
            cont = float(poisson_mix(tab.seq[0], lam, np.array([s0]))[0]) / per_window
            row[soil] = dict(independent=float(interp(tab.ratio, s0)) * share, air_dry=float(interp(tab.ratio, 0.0)) * share,
                             continuing_three_day_windows=cont * share, winter_rain_runs_off=bool(LS_RATE_MM_H > p['k_sub']))
        out[bracket] = row
    wg = out['walnut_gulch']
    lo, hi = min(v['independent'] for v in wg.values()), max(v['independent'] for v in wg.values())
    gm_chain = math.sqrt(wg['sandy_loam']['independent'] * wg['loam']['independent'])
    gm_meas = math.sqrt(MEASURED_RUNOFF[0] * MEASURED_RUNOFF[1])
    out['rule'] = ('Set before the chain was run: it passes when the Walnut Gulch storms at Earth\'s gravity on the two '
                   'Walnut Gulch-like soils, each storm on soil at the wilting point, give yearly runoff ratios whose span '
                   'overlaps the measured 3-10% of the small watersheds; the measured geometric mean over the chain\'s '
                   '(sandy loam and loam) is carried to the lunar bare-soil ratios as their Earth factor.')
    out['passes'] = bool(lo <= MEASURED_RUNOFF[1] and hi >= MEASURED_RUNOFF[0])
    out['earth_factor'] = gm_meas / gm_chain
    out['measured'] = list(MEASURED_RUNOFF)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split('\n')[0])
    ap.add_argument('--storms', default=STORMS, help='a CM1 storm-lives product for the lower bracket (repository path)')
    ap.add_argument('--extract', action='store_true', help='rerun the design-years extract first')
    a = ap.parse_args(argv)
    if a.extract or not EX.SIDECAR.is_file():
        EX.extract()
    x, meta = EX.load()
    libs = libraries(a.storms)
    pressure = mean_pressure_pa(x)
    kappa = float(IN.drop_energy_factor(pressure))
    store = {}

    def cache(p, lib, n_max):
        key = (tuple(sorted((k, v) for k, v in p.items() if not isinstance(v, dict))), lib['name'], n_max)
        if key not in store:
            store[key] = Tables(p, lib, n_max)
        return store[key]

    n_max = {b: n_needed(x, libs[b]) for b in BRACKETS}
    ls_cache = {}
    results, cell_rows = {}, {}
    for name in SURFACES:
        p = IN.parameters(name, pressure_pa=pressure)
        results[name], cell_rows[name] = {}, {}
        for b in BRACKETS:
            tab = cache(p, libs[b], n_max[b])
            results[name][b], cell_rows[name][b] = {}, {}
            for mode in MODES:
                conv, ls = runoff(x, tab, mode, ls_cache=ls_cache)
                cy, ct, cls = summarise(x, conv + ls)
                results[name][b][mode] = dict(classes=cls, convective_share_of_runoff=float(conv.sum() / max((conv + ls).sum(), 1e-30)))
                cell_rows[name][b][mode] = (cy, ct)
            print(f"{name} {b}: semi-arid {results[name][b]['independent']['classes']['semi-arid']['runoff_ratio']:.4f} "
                  f"/ {results[name][b]['continuing']['classes']['semi-arid']['runoff_ratio']:.4f}", file=sys.stderr)

    # Sensitivity of the semi-arid answer to each sub-model, the Walnut Gulch storms and independent storms at its base.
    semi = np.nonzero(x['cls'] == 'semi-arid')[0]
    xs = {k: (v[:, semi] if (isinstance(v, np.ndarray) and v.ndim == 2 and v.shape[1] == x['cls'].size) else
              (v[semi] if isinstance(v, np.ndarray) and v.ndim == 1 and v.size == x['cls'].size else v)) for k, v in x.items()}
    variants = {'base': {}, 'cm1_6km_storms': dict(bracket='cm1_6km'), 'continuing_front': dict(mode='continuing'),
                'air_dry_topsoil': dict(mapping='air_dry'), 'wetness_factor_topsoil': dict(mapping='wetness_factor'),
                'retention_0_mm': dict(retention_mm=0.0), 'retention_3_mm': dict(retention_mm=3.0),
                'large_scale_1_mm_h': dict(ls_rate=1.0), 'large_scale_4_mm_h': dict(ls_rate=4.0),
                'conductivity_half': dict(ks_factor=0.5), 'conductivity_double': dict(ks_factor=2.0),
                'suction_half': dict(suction_factor=0.5), 'suction_double': dict(suction_factor=2.0),
                'earth_gravity': dict(gravity=1.0)}
    surface_variants = {'sealed_loam': {'seal_forms_each_storm': dict(surface='sealed_loam_forming'),
                                        'seal_forms_kinetic_0.09': dict(surface='sealed_loam_forming', energy_factor=0.09),
                                        'seal_forms_kinetic_0.17': dict(surface='sealed_loam_forming', energy_factor=0.17),
                                        'seal_forms_kinetic_earth': dict(surface='sealed_loam_forming', energy_factor=1.0)},
                        'biocrust': {'crust_conductivity_0.35': dict(ratio=0.35), 'crust_conductivity_0.61': dict(ratio=0.61)}}
    sens = {}
    for name in IN.MAIN:
        sens[name] = {}
        for label, v in dict(variants, **surface_variants.get(name, {})).items():
            v = dict(v)
            b, mode, mapping, ls_rate = v.pop('bracket', 'walnut_gulch'), v.pop('mode', 'independent'), v.pop('mapping', 'plant_available'), v.pop('ls_rate', LS_RATE_MM_H)
            surf = v.pop('surface', name)
            p = IN.parameters(surf, pressure_pa=pressure, **v)
            tab = cache(p, libs[b], n_max[b] if mode == 'continuing' else 0)
            conv, ls = runoff(xs, tab, mode, mapping, ls_rate)
            sens[name][label] = float((conv + ls).sum(axis=0) @ xs['area_km2'] / ((xs['pr'].astype(float) * MM_PER_MEAN).sum(axis=0) @ xs['area_km2']))
        print(name, {k: round(v, 4) for k, v in sens[name].items()}, file=sys.stderr)

    check = walnut_gulch_check(libs, cache)
    print('walnut gulch', json.dumps(check, default=float)[:600], file=sys.stderr)

    # Storm statistics per cell.
    conv_mm = x['prc'].astype(float) * MM_PER_MEAN
    lunar_days = x['prc'].shape[0] / 10.0
    storms_cells = {}
    for b in BRACKETS:
        mu = ST.statistics(libs[b])['mean_depth_mm']
        lam = conv_mm / mu
        lt = by_tenth(lam, x['tenth'])
        ponds = {}
        for name in IN.MAIN:
            tab = cache(IN.parameters(name, pressure_pa=pressure), libs[b], n_max[b])
            s = IN.initial_saturation(tab.p['soil'], x['mrso'].astype(float))
            ponds[name] = (lam * interp(tab.ponding, s)).sum(axis=0) / np.maximum(lam.sum(axis=0), 1e-30)
        storms_cells[b] = dict(per_lunar_day=lam.sum(axis=0) / lunar_days, by_tenth=lt / lunar_days, ponding=ponds)

    from research.studies.storm_water import phase2a as P2
    phase2a, cells_2a = P2.run_phase2a(x, libs, pressure, cache, ls_cache)

    rain = x['pr'].astype(float) * MM_PER_MEAN
    by, mean = EX.composite(x)
    day = 86400e3                                                      # m/s of water to mm/day
    cells = []
    for c in range(x['cls'].size):
        row = dict(lat=r3(x['lat'][c]), lon_e=r3(x['lon'][c]), unep=str(x['cls'][c]), area_km2=r3(x['area_km2'][c]),
                   aridity=r3(x['aridity'][c]), rain_mm_yr=r3(mean['pr'][c] * day * JULIAN_YEAR_DAYS),
                   convective_share=r3(conv_mm[:, c].sum() / max(rain[:, c].sum(), 1e-30)),
                   evaporation_mm_yr=r3(-mean['evap'][c] * day * JULIAN_YEAR_DAYS), soil_water_mm=r3(1e3 * mean['mrso'][c]),
                   rain_mm_day_by_tenth=r3(by['pr'][:, c] * day), evaporation_mm_day_by_tenth=r3(-by['evap'][:, c] * day),
                   storms={b: dict(per_lunar_day=r3(storms_cells[b]['per_lunar_day'][c]),
                                   per_lunar_day_by_tenth=r3(storms_cells[b]['by_tenth'][:, c]),
                                   share_that_run_off={n: r3(storms_cells[b]['ponding'][n][c]) for n in IN.MAIN})
                           for b in BRACKETS},
                   runoff_ratio={n: {b: {m: dict(yearly=r3(cell_rows[n][b][m][0][c]), by_tenth=r3(cell_rows[n][b][m][1][:, c]))
                                         for m in MODES} for b in BRACKETS} for n in SURFACES})
        row.update(cells_2a[c])
        cells.append(row)

    surfaces_out = {}
    for name in SURFACES:
        pe, pm = IN.parameters(name, gravity=1.0, pressure_pa=pressure), IN.parameters(name, pressure_pa=pressure)
        keep = ('psi', 'ks', 'k_sub', 'zc', 'r', 'kc', 'retention')
        surfaces_out[name] = dict(design=IN.SURFACES[name], soil=IN.SOILS[IN.SURFACES[name]['soil']],
                                  earth={k: pe[k] for k in keep}, moon={k: pm[k] for k in keep},
                                  units=dict(psi='mm', ks='mm/h', k_sub='mm/h', zc='mm', r='h', kc='mm/h', retention='mm'))
    product = dict(
        schema=SCHEMA, evidence=EVIDENCE, reading_rule=READING_RULE,
        basis=dict(run=meta['run'], years=meta['years'], means=meta['means'], cells=meta['cells'],
                   land_shares=meta['land_shares'], mean_days=MEAN_DAYS, large_scale_rate_mm_h=LS_RATE_MM_H,
                   saturation_grid=list(SATURATIONS), drylands_pressure_pa=pressure, drop_energy_factor=kappa,
                   moon_gravity_in_earths=IN.MOON),
        storm_libraries={b: dict(statistics=ST.statistics(libs[b]), source=libs[b]['source'], evidence=libs[b]['evidence'],
                                 parameters=libs[b]['parameters'], fit=libs[b].get('fit')) for b in BRACKETS},
        surfaces=surfaces_out, classes=results, sensitivity_semi_arid=sens, walnut_gulch_check=check,
        thresholds=THRESHOLDS, phase2a=phase2a, cells=cells)
    inputs = {a.storms: sha256(a.storms)}
    inputs.update(phase2a['routing'].pop('drainage_sources'))
    product['producer'] = dict(study='research/studies/storm_water', files={p: sha256(p) for p in CODE},
                               inputs=inputs, ignored_inputs={str(EX.OUT.relative_to(ROOT)): meta['npz']['sha256']},
                               gcm_files=meta['files'], constants=constants_used(CODE))
    product['sizing'] = sizing(sens, check)
    OUT.parent.mkdir(exist_ok=True)
    head = json.dumps({k: v for k, v in product.items() if k != 'cells'}, indent=1, default=float)
    rows = ',\n'.join('  ' + json.dumps(c, separators=(',', ':'), default=float) for c in cells)
    OUT.write_text(head[:-2] + ',\n "cells": [\n' + rows + '\n ]\n}\n')
    return product


SUBMODELS = {
    'storm library': (('cm1_6km_storms',), 'Walnut Gulch: the family against Osborn\'s counts (fitted) and Mendez\'s annual '
                      'maxima (held out); no check of the lunar storms exists until the dry box runs'),
    'wetting between storms': (('continuing_front',), 'none: Walnut Gulch\'s storms come days apart, so its runoff '
                               'cannot separate the two treatments'),
    'soil at storm time': (('air_dry_topsoil', 'wetness_factor_topsoil'), 'Walnut Gulch chain at the wilting point and air-dry'),
    'surface retention': (('retention_0_mm', 'retention_3_mm'), 'Walnut Gulch chain (it carries the 1 mm store)'),
    'large-scale rain rate': (('large_scale_1_mm_h', 'large_scale_4_mm_h'), 'none (an estimate)'),
    'conductivity': (('conductivity_half', 'conductivity_double'), 'Walnut Gulch chain on two textures; the gravity '
                     'factor analytically'),
    'suction': (('suction_half', 'suction_double'), 'Walnut Gulch chain on two textures; the gravity factor analytically'),
    'gravity': (('earth_gravity',), 'analytic tests of each scaling'),
    'seal formation': (('seal_forms_each_storm', 'seal_forms_kinetic_0.09', 'seal_forms_kinetic_0.17', 'seal_forms_kinetic_earth'),
                       'the seal\'s resistance reproduces Rawls et al.\'s crust factor; the drops\' energy against the rain solver'),
    'biocrust conductivity': (('crust_conductivity_0.35', 'crust_conductivity_0.61'), 'none beyond the measured range of reductions'),
}


def sizing(sens, check):
    """For each sub-model, the semi-arid answer's range on each surface over its alternatives, beside the base, and
    the Earth check that bears on it with the resolution the measured runoff's spread allows."""
    spread = MEASURED_RUNOFF[1] / MEASURED_RUNOFF[0]
    out = dict(earth_resolution_factor=math.sqrt(spread), earth_factor=check['earth_factor'], submodels={})
    for name, (labels, tested) in SUBMODELS.items():
        row = {}
        for surface, v in sens.items():
            vals = [v[k] for k in labels if k in v]
            if vals:
                row[surface] = dict(base=v['base'], low=min(vals + [v['base']]), high=max(vals + [v['base']]),
                                    factor=max(vals + [v['base']]) / max(min(vals + [v['base']]), 1e-12))
        out['submodels'][name] = dict(variants=list(labels), tested_by=tested, semi_arid=row)
    return out


EVIDENCE = ('A model result of an offline chain on the design years, A28R_301k_design_l15_snow_s934 years 10-19 on '
            'PlaSim\'s T21 grid, with no feedback on the climate: each 3-day mean\'s convective rain as point storms from '
            'a library, the large-scale rain as a steady 2 mm/h (an estimate), and Green and Ampt infiltration at lunar '
            'gravity (conductivities times g, suction heads over g) into four designed surfaces whose Earth parameters are '
            'texture-class means (Rawls et al. 1983) and a seal from Rawls et al. (1990). The Walnut Gulch library is a '
            'two-parameter family fitted to Earth\'s measured point storms (Osborn 1983) and the CM1 library is CM1\'s '
            'model storms over wet equatorial ground at 6 km and ten minutes; neither is the drylands\' own storms. The '
            'soil at storm time reads PlaSim\'s 0.5 m bucket, which the arid, semi-arid and dry sub-humid cells keep '
            '0.4%, 2.7% and 8% full on average. Runoff here is what leaves '
            'the point as infiltration excess; it is not routed, and nothing of it is yet lost along slopes and channels '
            'or returned by pools. The surfaces are design choices, not measured lunar soils. Phase 2a (phase2a) '
            'replaces the two wetting treatments with a redistribution and drying model between storms at random times '
            'within each mean (wetting.py) and routes the central runoff on geography\'s 1.9 km grid into its '
            'depression tree (routing.py): run-on infiltration on 100 m slopes, sub-grid crater bowls from the terrain '
            'study\'s equilibrium counts for young and mature land, Lane\'s channel losses recast for gravity, and '
            'storms\' runoff-making cores of 5 km2, a value the Walnut Gulch basin check sets. Every routing parameter '
            'is a design choice or an Earth value; nothing in the lunar routing is measured.')
READING_RULE = ('Ratios are runoff over all rain (convective and large-scale) of the same cells and means, area-weighted '
                'over classes. walnut_gulch and cm1_6km bracket the storms\' intensities from above and below; independent '
                'and continuing bracket the wetting between a mean\'s storms (each on the mean\'s soil, or one front carried '
                'through them). sealed_loam keeps a formed seal; sealed_loam_forming forms it during each storm under the '
                'lunar drops\' kinetic energy. by_tenth lists the ten tenths of the lunar day from the one centred 162 '
                'degrees before noon to the one 162 after, noon falling between the fifth and sixth; null where a tenth '
                'brings under a thousandth of the cell\'s rain. share_of_evaporation is the runoff over the class\'s '
                'evaporation, an upper bound on the change of the slopes\' evaporation before pools return any of it. '
                'Volumes are km3 a year of 365.25 days; storms are expected point counts. In phase2a, wetting gives each '
                'surface\'s central runoff ratio with its low and high over the soil curves\' range; routing.fates give '
                'each fate of the runoff (bowls_inside, slope_soil, bowls_slope, channel_bed, bowls_channel, '
                'depressions, sea, humid_land, and the parts of the slopes\' and beds\' losses below the root zone, '
                'slope_recharge and bed_recharge) as shares of the runoff and of the rain, per hollow state (young: '
                'bowls of 10 m and more; mature: of 100 m and more); leverage gives the share of the rain leaving the '
                'soil\'s evaporating store (all fates but the soil above the root zone) over the evaporation. The cells\' '
                'runoff_ratio_2a is the central runoff; routed_walnut_gulch gives mm a year over the cell, by tenth, '
                'attributed to the cell where each loss or delivery happens. Per-depression inflow by tenth is in the '
                'ignored file named by phase2a.hollow_inflow, rebuilt by run.py.')


if __name__ == '__main__':
    main()
