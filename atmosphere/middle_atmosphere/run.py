"""Radiative-convective and photochemical equilibrium for the Earth control and the Moon under each shield.

    python -m atmosphere.middle_atmosphere.run            # every case
    python -m atmosphere.middle_atmosphere.run --only moon_1.2atm_edge_200nm
    python -m atmosphere.middle_atmosphere.run --list

Each case writes a profile table to results/profiles/ and its summary to
results/cases/; results/middle_atmosphere.csv gathers the summaries in case
order (finished cases are skipped, so a stopped run resumes, and several
runners can share the folder); results/middle_atmosphere.json carries the
schema, producer hashes and evidence statement.
"""
from __future__ import annotations
import argparse
import csv
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import sys
import time
import numpy as np

from atmosphere.radiative_convective import thermodynamics as th, fetch_inputs as rc_inputs
from atmosphere.middle_atmosphere import fetch_inputs, photolysis as ph, chemistry as ch, equilibrium as eq

HERE = Path(__file__).resolve().parent
SCHEMA = 'terluna.atmosphere.middle-atmosphere/1'
EVIDENCE = ('One-dimensional, global- and diurnal-mean, clear-sky radiative-convective equilibrium at a fixed surface '
            'temperature, coupled to steady-state O-H-N photochemistry with eddy diffusion. Moist-adiabatic troposphere '
            'with Manabe-Wetherald humidity; stratospheric water at the tropopause value unless prescribed. No clouds, '
            'circulation, waves, tides, chlorine, bromine, methane or aerosol chemistry; thermal emission in local '
            'thermodynamic equilibrium, which fails above about 1-10 Pa, except in the _nonlte cases (a two-level CO2 '
            '15-um band above 50 Pa, bounded by how the absorbed near-infrared sunlight becomes heat); no sunlight '
            'below 202 nm. Eddy diffusion is '
            'Earth-based and scaled for the Moon as stated per case. These are conditional calculations for comparing '
            'shields and compositions, not predictions of a lunar climate.')

MOON_12 = dict(planet=th.MOON, dry_pressure_pa=121590.0)
MOON_10 = dict(planet=th.MOON, dry_pressure_pa=101325.0)
MOON_KZZ = 36.0        # Earth eddy diffusion times (g_Earth/g_Moon)^2: the same mixing time per scale height
BALANCE_TEMPS = (268.0, 278.0, 298.0, 308.0)


def cases():
    out = [eq.Case('earth_control', th.EARTH, 101325.0, shield='none', stratospheric_h2o=5e-6, day_days=1.0)]
    for shield in ('titania_stack', 'edge_310nm', 'edge_200nm', 'none'):
        out.append(eq.Case(f'moon_1.2atm_{shield}', shield=shield, kzz_scale=MOON_KZZ, **MOON_12))
    out += [
        eq.Case('moon_1.0atm_edge_200nm', shield='edge_200nm', kzz_scale=MOON_KZZ, **MOON_10),
        eq.Case('moon_1.0atm_titania_stack', shield='titania_stack', kzz_scale=MOON_KZZ, **MOON_10),
        eq.Case('moon_1.2atm_edge_200nm_no_n2o', shield='edge_200nm', kzz_scale=MOON_KZZ, n2o_ppb=0.0, **MOON_12),
        eq.Case('moon_1.2atm_edge_200nm_slow_mixing', shield='edge_200nm', kzz_scale=1.0, **MOON_12),
        eq.Case('moon_1.2atm_edge_200nm_wet_stratosphere', shield='edge_200nm', kzz_scale=MOON_KZZ,
                stratospheric_h2o=5e-6, **MOON_12),
    ]
    out += [eq.Case(f'moon_1.2atm_edge_{nm}nm', shield=f'edge_{nm}nm', kzz_scale=MOON_KZZ, **MOON_12)
            for nm in (220, 230, 240)]
    out += [
        eq.Case('moon_1.0atm_edge_200nm_no_n2o', shield='edge_200nm', kzz_scale=MOON_KZZ, n2o_ppb=0.0, **MOON_10),
        eq.Case('moon_1.0atm_edge_200nm_slow_mixing', shield='edge_200nm', kzz_scale=1.0, **MOON_10),
    ]
    # The same columns at other surface temperatures, so balance temperatures are interpolated between
    # solved states rather than extrapolated from the one at 288 K.
    by_name = {c.name: c for c in out}
    # CO2 15-um cooling out of local thermodynamic equilibrium above 50 Pa, with the near-infrared sunlight
    # absorbed there either all turned to heat or only its collisionally deactivated fraction.
    for name in ('earth_control', 'moon_1.2atm_titania_stack', 'moon_1.0atm_titania_stack', 'moon_1.2atm_edge_200nm',
                 'moon_1.0atm_edge_200nm'):
        out.append(replace(by_name[name], name=f'{name}_nonlte', nonlte=True))
        out.append(replace(by_name[name], name=f'{name}_nonlte_nir_collisional', nonlte=True,
                           nir_thermalisation='collisional'))
    # CO2 is the upper air's main coolant: the design case with less of it, down to a plant-growth floor
    # (research/studies/atmospheric_co2), for its effect on the exobase.
    for co2 in (150.0, 280.0):
        out.append(replace(by_name['moon_1.2atm_titania_stack'], name=f'moon_1.2atm_titania_stack_nonlte_co2_{co2:g}',
                           nonlte=True, co2_ppm=co2))
    # The variants only need their balance with Earth-like clouds bracketed.
    variants = ('moon_1.2atm_edge_220nm', 'moon_1.2atm_edge_230nm', 'moon_1.2atm_edge_240nm', 'moon_1.2atm_edge_310nm',
                'moon_1.2atm_none', 'moon_1.2atm_edge_200nm_no_n2o', 'moon_1.0atm_edge_200nm_no_n2o',
                'moon_1.2atm_edge_200nm_slow_mixing', 'moon_1.0atm_edge_200nm_slow_mixing',
                'moon_1.2atm_edge_200nm_wet_stratosphere')
    for name, temps in (('earth_control', (278.0, 298.0)),
                        ('moon_1.2atm_titania_stack', BALANCE_TEMPS), ('moon_1.0atm_titania_stack', BALANCE_TEMPS),
                        ('moon_1.2atm_edge_200nm', BALANCE_TEMPS), ('moon_1.0atm_edge_200nm', BALANCE_TEMPS),
                        *((v, (278.0, 298.0)) for v in variants)):
        for ts in temps:
            out.append(replace(by_name[name], name=f'{name}_ts{ts:g}', surface_temperature_k=ts))
    return {c.name: c for c in out}


FIELDS = ['case', 'planet', 'dry_pressure_pa', 'shield', 'n2o_ppb', 'kzz_scale', 'stratospheric_h2o_prescribed',
          'surface_temperature_k', 'surface_pressure_pa', 'ozone_du', 'ozone_peak_ppm', 'ozone_peak_mixing_pa',
          'ozone_peak_density_pa', 'ozone_peak_density_km', 'tropopause_pa', 'tropopause_km', 'tropopause_k',
          'stratospheric_h2o', 'stratopause_k', 'stratopause_pa', 't_1000pa', 't_100pa', 't_10pa', 't_1pa',
          't_top', 'top_pressure_pa', 'top_km', 'noy_max_ppb', 'uv_index_mean', 'uv_index_subsolar',
          'uv_index_45deg', 'asr', 'olr', 'net_toa', 'asr_uv', 'asr_nir', 'uv_heating', 'day_night_k_100pa',
          'day_night_k_10pa', 'day_night_k_1pa', 'day_night_k_top', 'relaxation_days_100pa', 'relaxation_days_1pa',
          'column_heat_capacity_j_m2_k', 'night_column_cooling_k', 'outer_iterations', 'seconds']


def at_pressure(col, values, p):
    """Log-pressure interpolation of a layer quantity."""
    lp = np.log(col['layer_p_pa'])
    return float(np.interp(np.log(p), lp[::-1], values[::-1]))


def summarise(res: eq.Result, seconds):
    col, n = res.col, res.state.n
    m = ch.fixed_densities(col)['M']
    o3 = n[ch.INDEX['O3']]
    lp = col['layer_p_pa']
    k_mix, k_den = int(np.argmax(o3 / m)), int(np.argmax(o3))
    noy = sum(n[ch.INDEX[s]] for s in ('NO', 'NO2', 'NO3', 'HNO3')) + 2 * n[ch.INDEX['N2O5']]
    strat = np.arange(res.trop, col['layer_t_k'].size)
    k_sp = int(strat[np.argmax(col['layer_t_k'][strat])])
    dens = {s: n[ch.INDEX[s]] for s in ('O3', 'NO2', 'NO3', 'N2O5', 'N2O', 'HNO3', 'H2O2')}
    sun = ph.Sunlight(eq.cl.load_shield(res.case.shield))
    light = res.state.info['light']
    c = res.case
    return dict(
        case=c.name, planet=c.planet.name, dry_pressure_pa=c.dry_pressure_pa, shield=c.shield, n2o_ppb=c.n2o_ppb,
        kzz_scale=c.kzz_scale, stratospheric_h2o_prescribed=c.stratospheric_h2o or '',
        surface_temperature_k=c.surface_temperature_k, surface_pressure_pa=float(col['p_pa'][0]),
        ozone_du=ch.column_du(col, n), ozone_peak_ppm=float(o3[k_mix] / m[k_mix] * 1e6),
        ozone_peak_mixing_pa=float(lp[k_mix]), ozone_peak_density_pa=float(lp[k_den]),
        ozone_peak_density_km=float(col['layer_z_m'][k_den] / 1e3),
        tropopause_pa=float(col['p_pa'][res.trop]), tropopause_km=float(col['z_m'][res.trop] / 1e3),
        tropopause_k=float(res.t[res.trop]), stratospheric_h2o=float(col['x_h2o'][-1]),
        stratopause_k=float(col['layer_t_k'][k_sp]), stratopause_pa=float(lp[k_sp]),
        t_1000pa=at_pressure(col, col['layer_t_k'], 1000.0), t_100pa=at_pressure(col, col['layer_t_k'], 100.0),
        t_10pa=at_pressure(col, col['layer_t_k'], 10.0), t_1pa=at_pressure(col, col['layer_t_k'], 1.0),
        t_top=float(col['layer_t_k'][-1]), top_pressure_pa=float(lp[-1]), top_km=float(col['layer_z_m'][-1] / 1e3),
        noy_max_ppb=float(np.max(noy / m) * 1e9),
        uv_index_mean=ph.uv_index(light['surface_spectrum']),
        uv_index_subsolar=ph.uv_index(sun.noon_surface(col, dens, c.surface_albedo, 1.0)),
        uv_index_45deg=ph.uv_index(sun.noon_surface(col, dens, c.surface_albedo, float(np.cos(np.pi / 4)))),
        asr=res.budget['asr'], olr=res.budget['olr'], net_toa=res.budget['net_toa'], asr_uv=res.budget['asr_uv'],
        asr_nir=res.budget['asr_nir'], uv_heating=float(res.state.heating.sum()),
        **_cycle_summary(res.day_night), **_night(res), outer_iterations=len(res.history), seconds=round(seconds, 1))


def _night(res):
    """Whole-column heat capacity and the cooling that the mean outgoing longwave alone would cause over one
    night (half the solar day), with no sunlight and no heat brought in by the circulation."""
    col = res.col
    cp = eq.col_cp(col, eq.Profile(res.case).air)
    capacity = float(np.sum(col['layer_mass_kg_m2'] * cp))
    night_s = 0.5 * res.case.day_days * 86400.0
    return dict(column_heat_capacity_j_m2_k=capacity, night_column_cooling_k=res.budget['olr'] * night_s / capacity)


def _cycle_summary(cycle):
    lp = np.log(cycle['p_pa'])
    at = lambda values, p: float(np.interp(np.log(p), lp[::-1], values[::-1]))
    amp, tau = 2 * cycle['amplitude_k'], cycle['relaxation_days']
    return dict(day_night_k_100pa=at(amp, 100.0), day_night_k_10pa=at(amp, 10.0), day_night_k_1pa=at(amp, 1.0),
                day_night_k_top=float(amp[-1]), relaxation_days_100pa=at(tau, 100.0), relaxation_days_1pa=at(tau, 1.0))


def profile_rows(res: eq.Result, solver: eq.Solver):
    col, n = res.col, res.state.n
    m = ch.fixed_densities(col)['M']
    cp = eq.col_cp(col, eq.Profile(res.case).air)
    heat = lambda w: w / (col['layer_mass_kg_m2'] * cp) * eq.SECONDS_PER_DAY
    rows = []
    cycle = res.day_night
    for k in range(col['layer_p_pa'].size):
        row = dict(layer=k, p_pa=float(col['layer_p_pa'][k]), z_km=float(col['layer_z_m'][k] / 1e3),
                   t_k=float(col['layer_t_k'][k]), h2o=float(col['layer_x_h2o'][k]),
                   sw_heating_k_day=float(heat(res.sw_layers)[k]), lw_heating_k_day=float(heat(res.lw['absorbed'])[k]),
                   uv_heating_k_day=float(heat(res.state.heating)[k]))
        j = k - res.trop
        row['day_night_k'] = float(2 * cycle['amplitude_k'][j]) if j >= 0 else ''
        row['relaxation_days'] = float(cycle['relaxation_days'][j]) if j >= 0 else ''
        for s in ch.SPECIES:
            row[s] = float(n[ch.INDEX[s], k] / m[k])
        for name in ('O2 -> O + O', 'O3 -> O1D + O2', 'NO2 -> NO + O'):
            row['J(' + name.split(' ->')[0] + ('*' if 'O1D' in name else '') + ')'] = float(res.state.j_rates[name][k])
        rows.append(row)
    return rows


def producer():
    code = {}
    for folder in (HERE, HERE.parent / 'radiative_convective'):
        for p in sorted(folder.glob('*.py')):
            code[f'{folder.name}/{p.name}'] = hashlib.sha256(p.read_bytes()).hexdigest()[:16]
    inputs = {name: hashlib.sha256((folder / 'inputs.json').read_bytes()).hexdigest()[:16]
              for name, folder in (('middle_atmosphere', HERE), ('radiative_convective', HERE.parent / 'radiative_convective'))}
    shield = HERE.parents[1] / 'protection' / 'spectra' / 'shield_transmission.json'
    return dict(domain='atmosphere', model='atmosphere/middle_atmosphere', files=code, inputs_manifests=inputs,
                shield_product=hashlib.sha256(shield.read_bytes()).hexdigest()[:16])


def _write_csv(path, rows, fields=None):
    fields = fields or list(rows[0].keys())
    tmp = path.with_suffix('.tmp')
    with open(tmp, 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        for r in rows:
            writer.writerow({f: (f'{v:.6g}' if isinstance(v, float) else v) for f, v in r.items() if f in fields})
    tmp.replace(path)


def _table(out, plan):
    """Summary table in case order from the per-case records, so several runners can share a results folder."""
    rows = []
    for name in plan:
        record = out / 'cases' / f'{name}.json'
        if record.is_file():
            rows.append(json.loads(record.read_text()))
    if rows:
        _write_csv(out / 'middle_atmosphere.csv', rows, FIELDS)
    return rows


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out', type=Path, default=HERE / 'results')
    parser.add_argument('--only', nargs='*')
    parser.add_argument('--list', action='store_true')
    parser.add_argument('--verbose', action='store_true')
    parser.add_argument('--processes', type=int, default=0, help='worker processes for line-by-line sunlight')
    args = parser.parse_args(argv)
    if args.processes:
        eq.SW_CONFIG = replace(eq.optics.SHORTWAVE, processes=args.processes)
    plan = cases()
    if args.list:
        for name in plan:
            print(name)
        return 0
    for restore in (fetch_inputs.restore, rc_inputs.restore):
        if restore()['status'] != 'PASS':
            print('BLOCKED: inputs missing; run the fetch_inputs scripts with --download', file=sys.stderr)
            return 1
    for sub in ('profiles', 'cases'):
        (args.out / sub).mkdir(parents=True, exist_ok=True)
    table = args.out / 'middle_atmosphere.csv'
    solver = eq.Solver(verbose=args.verbose)
    for name, case in plan.items():
        if (args.only and name not in args.only) or (args.out / 'cases' / f'{name}.json').is_file():
            continue
        start = time.time()
        res = solver.solve(case)
        summary = summarise(res, time.time() - start)
        _write_csv(args.out / 'profiles' / f'{name}.csv', profile_rows(res, solver))
        (args.out / 'cases' / f'{name}.json').write_text(json.dumps(summary, indent=1) + '\n')
        _table(args.out, plan)
        print(f"{name:42s} O3 {summary['ozone_du']:6.1f} DU  tropopause {summary['tropopause_k']:5.1f} K  "
              f"stratopause {summary['stratopause_k']:5.1f} K  T(top) {summary['t_top']:5.1f} K  "
              f"UV index {summary['uv_index_subsolar']:4.1f}  net TOA {summary['net_toa']:6.1f} W/m2", flush=True)
    meta = dict(schema=SCHEMA, producer=producer(), evidence=EVIDENCE,
                units=dict(pressure='Pa', temperature='K', ozone='Dobson units (2.6867e16 cm^-2)', fluxes='W m-2 global mean',
                           mixing_ratios='mole fraction', photolysis='s-1 global mean', heating='K/day'),
                reading_rule=('Compare cases at the same surface temperature; net_toa is the forcing the clear column '
                              'lacks (positive: it would warm). Above ~1 Pa the cases without _nonlte overstate CO2 '
                              'cooling; each _nonlte pair bounds the upper air (all absorbed near-infrared sunlight '
                              'heats, or only its collisionally deactivated part). Heating from above (extreme '
                              'ultraviolet, conduction) is not included; use t_top only as the boundary of the '
                              'thermal column with that caveat. Keep the evidence statement with any use.'),
                table=table.name, profiles='profiles/<case>.csv')
    (args.out / 'middle_atmosphere.json').write_text(json.dumps(meta, indent=2) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
