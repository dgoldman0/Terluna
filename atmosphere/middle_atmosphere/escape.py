"""Exobase temperature and molecular escape with the middle atmosphere as the thermal column's lower boundary.

    python -m atmosphere.middle_atmosphere.escape      # after run.py; writes results/escape_coupling.csv

The thermal column (atmosphere/thermal_column.py) solves conduction, advected
energy and Jeans escape above a prescribed base temperature and pressure for a
given heat deposited above that base. Its screens so far assumed a base of
180 K at 0.1 Pa. Here the base comes from each radiative-convective case. The
main base is 0.3 Pa, where the correlated-k and line-by-line temperatures agree
within about 1 K (lbl_check.py). The model's top layer at 0.1 Pa and, where
lbl_check.json holds the case, its line-by-line corrected value are
sensitivities: that layer absorbs the line-core sunlight meant for everything
above it, and the two schemes differ there by 30-40 K. The 10 Pa base, the
highest level where CO2 emission is near local thermodynamic equilibrium, is
another; from there up the thermal column has no infrared cooling, so it gives
an upper bound. The deposited heat stays a swept input: what the shield lets
through below 200 nm, and how much of it becomes heat, is not computed here.

Carried with every row: the base temperatures come from a global-mean,
radiative-only column (no waves, tides or day-night circulation, which the
Moon's month-long day will drive strongly), LTE cooling above ~1-10 Pa, and the
thermal column's own limits (molecular N2/O2 only, no atomic oxygen, no
infrared cooling). The atomic-oxygen fraction the chemistry finds at the top is
reported because it is outside the thermal column's assumptions.
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np

from atmosphere.thermal_column import ColumnConfig, solve_column
from atmosphere.radiative_convective import thermodynamics as th

HERE = Path(__file__).resolve().parent
HEATS = (0.0, 1e-6, 3e-6, 1e-5)
BASELINE_K = 180.0


def _profile(results, case):
    with open(results / 'profiles' / f'{case}.csv', newline='') as handle:
        rows = list(csv.DictReader(handle))
    return {k: np.array([float(r[k]) if r[k] != '' else np.nan for r in rows]) for k in rows[0]}


def base_temperature(profile, p):
    lp = np.log(profile['p_pa'])
    return float(np.interp(np.log(p), lp[::-1], profile['t_k'][::-1]))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--results', type=Path, default=HERE / 'results')
    args = parser.parse_args(argv)
    table = args.results / 'middle_atmosphere.csv'
    if not table.is_file():
        print('No middle-atmosphere results; run atmosphere.middle_atmosphere.run first', file=sys.stderr)
        return 1
    cases = [r for r in csv.DictReader(open(table, newline='')) if r['planet'] == 'moon']
    check_path = args.results / 'lbl_check.json'
    checked = json.loads(check_path.read_text())['cases'] if check_path.is_file() else {}
    out = []
    for row in cases:
        prof = _profile(args.results, row['case'])
        dry = float(row['dry_pressure_pa'])
        o2 = th.earthlike_air(dry, 400.0).fractions['O2']
        o2_molecular = o2 / (o2 + th.earthlike_air(dry, 400.0).fractions['N2'])
        o_top = float(prof['O'][-1])
        bases = [('middle_atmosphere', 0.3, base_temperature(prof, 0.3)),
                 ('middle_atmosphere_top', 0.1, base_temperature(prof, 0.1)),
                 ('middle_atmosphere', 10.0, base_temperature(prof, 10.0)),
                 ('baseline_180K', 0.1, BASELINE_K)]
        if row['case'] in checked:
            c = checked[row['case']]
            bases.append(('line_by_line_corrected_top', 0.1, c['t_lbl_corrected_k'][c['pressures_pa'].index(0.1)]))
        for label, base_p, t_base in bases:
            for q in HEATS:
                    cfg = ColumnConfig(surface_pressure_pa=float(row['surface_pressure_pa']),
                                       surface_temperature_k=float(row['surface_temperature_k']),
                                       lower_temperature_k=t_base, lower_pressure_pa=base_p,
                                       oxygen_mole_fraction=o2_molecular)
                    rec = dict(case=row['case'], shield=row['shield'], dry_pressure_pa=dry, base=label,
                               base_pressure_pa=base_p, base_temperature_k=round(t_base, 2),
                               deposited_heat_w_m2=q, atomic_o_fraction_at_top=o_top)
                    try:
                        summary, _, _ = solve_column(q, cfg)
                        rec.update(status='converged', exobase_temperature_k=summary['exobase_temperature_k'],
                                   exobase_radius_R=summary['exobase_radius_R'],
                                   molecular_loss_kg_s=summary['molecular_loss_kg_s'],
                                   jeans_lambda_N2=summary['jeans_lambda_N2'],
                                   flags=' '.join(summary['domain_flags']))
                    except (ValueError, RuntimeError, FloatingPointError) as exc:
                        rec.update(status='failed', error=str(exc))
                    out.append(rec)
                    print(f"{row['case']:40s} {label:18s} p_b={base_p:5g} T_b={t_base:6.1f} q={q:7.1e} -> "
                          f"T_exo={rec.get('exobase_temperature_k', float('nan')):7.1f} K "
                          f"loss={rec.get('molecular_loss_kg_s', float('nan')):.2e} kg/s", flush=True)
    fields = list(dict.fromkeys(k for r in out for k in r))
    with open(args.results / 'escape_coupling.csv', 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(out)
    meta = dict(schema='terluna.atmosphere.escape-coupling/1', source_table=table.name,
                evidence=' '.join(part.replace('\n', ' ') for part in __doc__.split('\n\n')[1:]),
                reading_rule=('Compare exobase temperature and molecular loss between the middle-atmosphere base '
                              'and the 180 K baseline at the same deposited heat. Deposited heat is an input, not '
                              'a shield result. Flags are the thermal column\'s own domain warnings.'))
    (args.results / 'escape_coupling.json').write_text(json.dumps(meta, indent=2) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
