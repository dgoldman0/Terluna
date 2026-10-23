"""Audit and summarize dilute native crystals in prescribed lunar weather."""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np

from atmosphere.snow import mcsnow_weather as w
from shared.constants import CM1_DRY_AIR_GAS_CONSTANT, CM1_VAPOUR_GAS_CONSTANT

stats = w.history.stats


def require(condition, message):
    if not condition:
        raise ValueError(message)


def verify_files(path, meta):
    for filename, key in (('raw.npz', 'raw_sha256'), ('producer.py', 'producer_sha256')):
        require(w.native.sha(path / filename) == meta[key], 'Run file changed: ' + filename)
    require(w.native.sha(w.__file__) == meta['producer_sha256'], 'Executed producer differs')
    for filename, key in (('build_moon/build.json', 'build_record_sha256'),
                          ('weather/input.json', 'weather_record_sha256')):
        require(w.native.sha(w.CACHE / filename) == meta[key], 'Run input changed: ' + filename)


def validate_arrays(a, meta, record, ice_density):
    n = len(record['releases'])
    require(all(np.all(np.isfinite(v)) for v in a.values()), 'Nonfinite output')
    require(a['states'].shape == (5, n), 'Wrong number of states')
    require(np.array_equal(a['release'], record['releases']), 'Release cells differ')
    require(np.all((a['status'] >= 1) & (a['status'] <= 7)), 'Unfinished or unknown status')
    state, diag, rows = a['states'], a['final_diagnostics'], a['rows']
    require(np.all(state[:3] > 0), 'Nonpositive crystal state')
    volume = math.pi / 6 * diag[1]**2 * diag[2]
    geometry_error = float(np.max(np.abs(volume / state[1] - 1)))
    require(geometry_error < 1e-11, 'Axes do not represent the retained volume')
    require(np.allclose(diag[2] / diag[1], state[2], rtol=1e-12, atol=0), 'Aspect ratio differs')
    require(np.allclose(state[0] / volume, diag[3], rtol=1e-11, atol=0), 'Density differs')
    require(np.array_equal(diag[0], np.maximum(diag[1], diag[2])), 'Maximum dimension differs')
    require(np.all(a['peak_diameter_m'] >= diag[0] * (1 - 1e-12)), 'Lost final maximum dimension')
    balance = state[0] - a['initial_mass_kg'] - a['deposited_mass_kg'] + a['sublimated_mass_kg']
    scale = a['initial_mass_kg'] + a['deposited_mass_kg'] + a['sublimated_mass_kg']
    balance_error = float(np.max(np.abs(balance) / scale))
    require(balance_error < 1e-11, 'Particle mass ledger fails')
    require(np.array_equal(balance, a['mass_balance_residual_kg']), 'Saved mass residual differs')
    require(np.all(a['moist_time_s'] <= a['age_s'] * (1 + 1e-12)), 'Moist exposure exceeds age')
    cfg = meta['configuration']
    min_mass = ice_density * math.pi / 6 * cfg['cutoff_mass_equivalent_diameter_m']**3
    require(np.all(state[0, a['status'] == 1] <= min_mass * (1 + 1e-12)), 'Residual cutoff not reached')
    require(np.all(a['peak_diameter_m'][a['status'] == 2] >= cfg['diameter_stop_m']), 'Diameter censor not reached')
    end = (record['shape'][0] - 1) * 600.
    for i in range(n):
        r = rows[rows[:, 0] == i]
        require(len(r) >= 2 and r[0, 2] == 0, 'Missing initial/terminal snapshot')
        require(np.all(np.diff(r[:, 1]) >= 0), 'Nonmonotone path')
        require(np.allclose(r[:, 1] - a['release'][i, 0], r[:, 2], rtol=0, atol=1e-7), 'Path time mismatch')
        require(r[-1, 18] == a['status'][i], 'Terminal snapshot status differs')
        require(np.array_equal(r[-1, 6:9], state[:3, i]), 'Terminal snapshot state differs')
        require(np.array_equal(r[-1, 9:14], diag[:5, i]), 'Terminal snapshot geometry differs')
        require(r[-1, 2] == a['age_s'][i] and r[-1, 1] <= end + 1e-7, 'Terminal time differs')
        require(np.max(r[:, 9]) <= a['peak_diameter_m'][i] * (1 + 1e-12), 'Lost saved peak')
    require(np.all((rows[:, 5] >= record['z_m'][0]) & (rows[:, 5] <= record['z_m'][-1])), 'Height extrapolation')
    return dict(snapshots=len(rows), mass_balance_relative=balance_error,
                geometry_relative=geometry_error,
                native_mass_balance_relative=meta['maximum_native_mass_balance_relative'])


def summarize(a, meta):
    terminal = a['status']
    return dict(configuration=meta['configuration'],
        terminations={v: int(np.sum(terminal == k)) for k, v in w.STATUSES.items() if k},
        peak_dimension_mm=stats(a['peak_diameter_m'] * 1000),
        final_dimension_mm_by_termination={v: stats(a['final_diagnostics'][0, terminal == k] * 1000)
                                          for k, v in w.STATUSES.items() if k},
        peak_exceedance_counts={str(mm): int(np.sum(a['peak_diameter_m'] >= mm / 1000))
                                for mm in (1, 5, 10, 100, 1000)},
        age_h=stats(a['age_s'] / 3600), moist_exposure_min=stats(a['moist_time_s'] / 60),
        deposited_mass_mg=stats(a['deposited_mass_kg'] * 1e6),
        sublimated_mass_mg=stats(a['sublimated_mass_kg'] * 1e6),
        final_mass_mg=stats(a['states'][0] * 1e6),
        residual_fraction=float(np.mean(terminal == 1)),
        diameter_censored_ids=np.flatnonzero(terminal == 2).tolist(),
        diameter_censored_states=[dict(id=int(i), age_h=float(a['age_s'][i]/3600),
            mass_mg=float(a['states'][0,i]*1e6), aspect_ratio=float(a['states'][2,i]),
            equatorial_mm=float(a['final_diagnostics'][1,i]*1000),
            polar_mm=float(a['final_diagnostics'][2,i]*1000)) for i in np.flatnonzero(terminal == 2)],
        substeps=meta['substeps'], wall_s=meta['wall_s'])


def compare(a, b):
    """b is the reference; retain changed outcomes and the full extreme tail."""
    pa, pb = a['peak_diameter_m'], b['peak_diameter_m']
    return dict(status_agreement_count=int(np.sum(a['status'] == b['status'])),
                paths=len(pa), peak_relative_difference=stats(np.abs(pa / pb - 1)),
                peak_mm_absolute_difference=stats(np.abs(pa - pb) * 1000),
                peak_quantile_relative_difference={str(q): float(np.quantile(pa, q) / np.quantile(pb, q) - 1)
                                                   for q in (.1, .5, .9, .99)},
                age_absolute_difference_s=stats(np.abs(a['age_s'] - b['age_s'])))


def moisture_inventory(lib, record):
    """Instantaneous stocks in the archived -20..0 C layers, not a water flux."""
    data = np.load(w.CACHE / 'weather/fields.npy', mmap_mode='r')
    interfaces = np.loadtxt(w.history.RUNS / w.history.CASE / 'input_grid_z')
    dz = np.diff(interfaces)[record['selected_levels']][:, None, None]
    columns = {key: [] for key in ('vapour', 'excess_over_ice', 'cloud_liquid')}
    humidity = {key: [] for key in ('winter', 'snow', 'moist_snow')}
    for frame in data:
        raw = np.moveaxis(frame, 0, -1).reshape(-1, len(w.history.t.FIELDS))
        air = w.weather_air(raw)
        T, p, pv = air
        mask = (T >= 253.15) & (T < 273.15)
        force = np.asfortranarray(air[:, mask])
        state = lib.seed(force, 4e-6)
        diag, _ = lib.advance(state, force, 0.)
        rhw, rhi = force[2] / diag[10], force[2] / diag[11]
        snow = raw[mask, 5] >= w.history.SNOW_THRESHOLD
        for key, select in (('winter', np.ones(len(rhi), bool)), ('snow', snow),
                            ('moist_snow', snow & (rhi >= 1))):
            humidity[key].append(np.stack((rhi[select], rhw[select]), axis=1))
        rho_d = (p - pv) / (CM1_DRY_AIR_GAS_CONSTANT * T)
        excess = np.zeros(len(T))
        excess[mask] = np.maximum(force[2] - diag[11], 0) / (CM1_VAPOUR_GAS_CONSTANT * force[0])
        for key, quantity in (('vapour', rho_d * raw[:, 2] * mask),
                              ('cloud_liquid', rho_d * raw[:, 3] * mask), ('excess_over_ice', excess)):
            columns[key].append(np.sum(quantity.reshape(frame.shape[1:]) * dz, axis=0).ravel())
    answer = dict(scope='All 211 frames and horizontal columns; whole selected cells whose centre is -20 <= T < 0 C.',
                  units='kg/m2 inventory, numerically equivalent to millimetres of liquid water; no supply per hour is inferred.',
                  column_stock={key: stats(np.concatenate(values)) for key, values in columns.items()},
                  resolved_grid_humidity={})
    for key, arrays in humidity.items():
        both = np.concatenate(arrays)
        answer['resolved_grid_humidity'][key] = dict(ice_rh=stats(both[:, 0]), water_rh=stats(both[:, 1]),
            fraction_above_water_saturation=float(np.mean(both[:, 1] > 1)),
            fraction_above_1_1_water_rh=float(np.mean(both[:, 1] > 1.1)))
    return answer


def snapshot_humidity(lib, a):
    r = a['rows']
    force = np.asfortranarray(r[:, [14, 15, 16]].T)
    state = np.asfortranarray(np.vstack((r[:, 6:9].T, r[:, 14], np.zeros(len(r)))))
    d, _ = lib.advance(state, force, 0.)
    return dict(scope='Saved ten-minute and terminal path samples, unweighted; extrema between samples can be missed.',
        ice_rh=stats(force[2] / d[11]), water_rh=stats(force[2] / d[10]),
        fraction_above_water_saturation=float(np.mean(force[2] > d[10])))


def assemble():
    record = w.verify_weather()
    lib = w.Library()
    arrays, runs, audit, identities = {}, {}, {}, {}
    matrix = w.matrix_requests()
    require(json.loads((w.CACHE / 'matrix.json').read_text()) == matrix, 'Matrix request identity differs')
    for case in matrix:
        name = case['name']; path = w.CACHE / 'runs' / name
        meta = json.loads((path / 'run.json').read_text())
        require(json.loads((path / 'request.json').read_text()) == case['options'], 'Requested options differ')
        verify_files(path, meta)
        with np.load(path / 'raw.npz', allow_pickle=False) as f:
            arrays[name] = {k: f[k] for k in f.files}
        audit[name] = validate_arrays(arrays[name], meta, record, lib.ice_density)
        runs[name] = summarize(arrays[name], meta)
        identities[name] = dict(run_record_sha256=w.native.sha(path / 'run.json'),
            raw_sha256=meta['raw_sha256'], producer_sha256=meta['producer_sha256'],
            request_sha256=w.native.sha(path / 'request.json'))
    numerical, growth, cadence, cutoff = {}, {}, {}, {}
    for size in ('small', 'large'):
        fine = arrays[f'actual_{size}_step10']
        for kind in ('actual', 'saturated', 'capped'):
            numerical[f'{kind}_{size}'] = compare(arrays[f'{kind}_{size}_step20'], arrays[f'{kind}_{size}_step10'])
        for minutes in (20, 30):
            cadence[f'{size}_{minutes}min'] = compare(arrays[f'actual_{size}_fields{minutes}'], fine)
        cutoff[size] = compare(arrays[f'actual_{size}_cutoff01'], fine)
        growth[size] = compare(arrays[f'actual_{size}_growth_refined'], fine)
    sources = {str(Path(p).relative_to(w.ROOT)): w.native.sha(p) for p in
        (w.__file__, __file__, w.DRIVER, w.history.__file__, w.history.t.__file__,
         w.lunar.__file__, w.lunar.trace.__file__, w.ROOT / 'shared/constants.json')}
    return dict(schema='terluna.atmosphere.mcsnow-weather/1', producer=sources,
        evidence='Conditional dilute single-crystal calculations in archived CM1 lunar weather; no empirical validation or coupled population forecast.',
        reading_rule='Compare paired release cells. Peak dimensions precede any later sublimation; terminal states are censor/temperature/archive exits, not ground snowfall. Retain numeric extremes as model diagnostics.',
        weather=dict(input_record_sha256=w.native.sha(w.CACHE / 'weather/input.json'),
                     fields_sha256=record['fields_sha256'], inventory=record['inventory'],
                     frames=record['frames'], z_m=record['z_m'], releases=record['releases'],
                     field_order=record['field_order'], shape=record['shape']),
        build_record_sha256=w.native.sha(w.CACHE / 'build_moon/build.json'),
        gravity_m_s2=lib.gravity, run_identity=identities, runs=runs, raw_audit=audit,
        step_refinement=numerical, growth_refinement=growth, field_cadence=cadence, residual_cutoff=cutoff,
        moisture_inventory=moisture_inventory(lib, record),
        path_humidity={size: snapshot_humidity(lib, arrays[f'actual_{size}_step10']) for size in ('small', 'large')},
        limitations=[
            'The older A28_dim5 ten-minute archive is used; newer design-climate storm fields are three-hourly.',
            'Weather is imposed. These unweighted probes do not deplete its vapour or heat the air; no number flux or per-crystal water allowance is inferred.',
            'Release cells were selected for snow and ice saturation in -20..0 C; counts are conditional diagnostics, not population fractions.',
            'The native improved-control shape law retains Earth errors. Riming, collisions, breakup, melting and crystal nucleation are not simulated.',
            'Paths end at -30 C, freezing, domain height bounds, archived time, a tiny mass residual or the declared 1 m numeric censor.',
            'Spatial interpolation and ten-minute fields cannot resolve the complete supersaturation history. Coarser-field differences are sensitivity tests, not convergence to unresolved weather.',
            'The capped sensitivity removes only supersaturation over water; the saturated control adds water everywhere. Neither is a coupled condensation calculation.'])


if __name__ == '__main__':
    w.save(w.PRODUCT, assemble())
    print(w.PRODUCT)
