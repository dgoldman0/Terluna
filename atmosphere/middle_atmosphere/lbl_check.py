"""How much the correlated-k CO2 cooling error moves the upper temperatures of a finished case.

    python -m atmosphere.middle_atmosphere.lbl_check moon_1.2atm_titania_stack earth_control

Above ~30 Pa (Earth) and ~5 Pa (Moon) the correlated-k scheme's CO2 15-um cooling
is 20-40% off the line-by-line value, and several times off in the top layer
(atmosphere/radiative_convective/validation.json). For each named case this
recomputes, on the case's final profile, the 500-820 cm^-1 heating line by line
at Doppler-resolving resolution and with the correlated-k scheme, adds the
difference to the absorbed sunlight as a fixed correction, and solves the
radiative equilibrium again with the chemistry held fixed. The temperature
change is the correlated-k contribution to the upper-air uncertainty; the
larger error from assuming local thermodynamic equilibrium is not addressed.
Writes results/lbl_check.json.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np

from atmosphere.radiative_convective import validate_ck as vc
from atmosphere.middle_atmosphere import equilibrium as eq, run as rn, chemistry as ch

HERE = Path(__file__).resolve().parent


def rebuild(case, results):
    rows = list(csv.DictReader(open(results / 'profiles' / f'{case.name}.csv', newline='')))
    summary = json.loads((results / 'cases' / f'{case.name}.json').read_text())
    prof = eq.Profile(case)
    layer_t = np.array([float(r['t_k']) for r in rows])
    levels = [case.surface_temperature_k]
    for tl in layer_t:
        levels.append(2 * tl - levels[-1])                   # layer temperatures are level means
    t = np.array(levels)
    trop = int(np.argmin(np.abs(np.log(prof.p / summary['tropopause_pa']))))
    o3 = np.array([float(r['O3']) for r in rows])
    col = prof.column(t, trop)
    cp = eq.col_cp(col, prof.air)
    sw = np.array([float(r['sw_heating_k_day']) for r in rows]) * col['layer_mass_kg_m2'] * cp / eq.SECONDS_PER_DAY
    return prof, t, trop, o3, col, sw


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('cases', nargs='+')
    parser.add_argument('--results', type=Path, default=HERE / 'results')
    parser.add_argument('--processes', type=int, default=6)
    args = parser.parse_args(argv)
    plan = rn.cases()
    out_path = args.results / 'lbl_check.json'
    record = json.loads(out_path.read_text()) if out_path.is_file() else dict(
        schema='terluna.atmosphere.middle-atmosphere-lbl-check/1', evidence=__doc__.split('\n\n')[1].replace('\n', ' '),
        cases={})
    solver = eq.Solver()
    for name in args.cases:
        case = plan[name]
        prof, t, trop, o3, col, sw = rebuild(case, args.results)
        col['o3'] = o3 * col['layer_column_cm2']
        net_lbl = vc.lbl_band_heating(col, 500.0, 820.0, 2e-4, ('H2O', 'CO2', 'O3'), args.processes)
        net_ck = vc.ck_band_net(solver.lw, col, (2, 3, 4))
        correction = (net_lbl[:-1] - net_lbl[1:]) - (net_ck[:-1] - net_ck[1:])      # W/m^2 per layer
        t_new, q, col_new, _ = solver.radiative_equilibrium(prof, t.copy(), trop, sw + correction, o3)
        lt_old, lt_new = col['layer_t_k'], col_new['layer_t_k']
        lp = np.log(col['layer_p_pa'])
        at = lambda v, p: float(np.interp(np.log(p), lp[::-1], v[::-1]))
        entry = dict(pressures_pa=[100.0, 30.0, 10.0, 3.0, 1.0, 0.3, 0.1])
        entry['t_ck_k'] = [round(at(lt_old, p), 2) for p in entry['pressures_pa']]
        entry['t_lbl_corrected_k'] = [round(at(lt_new, p), 2) for p in entry['pressures_pa']]
        entry['max_heating_correction_k_day'] = float(np.max(np.abs(correction / (col['layer_mass_kg_m2'] * eq.col_cp(col, prof.air)) * eq.SECONDS_PER_DAY)))
        record['cases'][name] = entry
        print(name, list(zip(entry['pressures_pa'], entry['t_ck_k'], entry['t_lbl_corrected_k'])), flush=True)
        out_path.write_text(json.dumps(record, indent=1) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
