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
an upper bound. The deposited heat is swept, and also estimated for a filter
that passes 0.1% of sunlight below its edge (the idealised filters of the shield
product): leakage_heat() takes the WHI 2008 reference spectrum below 175 nm that
is absorbed above the base, with a heating efficiency of 0.4, for quiet-Sun and
rough solar-maximum irradiance. For the titania stack, film_heat() also takes the
film's own transmission from the protection domain's short-wave product
(protection/spectra/stack_short_wave.json): it passes only hard X-rays, all
counted as absorbed above the base (an upper bound, since the hardest reach
deeper). The 0.1% rows then stand for light that bypasses the film through gaps
in the aperture.

Atomic oxygen, which the thermal column does not hold, is carried on each solution
as a trace gas from the fraction the chemistry finds at the base, with Jeans
escape at the column's exobase. Above the base its fraction follows the steady
balance of eddy and molecular diffusion: well mixed (eddy mixing reaching the
exobase), separated from the base up (no eddy mixing; its own scale height,
28/16 of N2's), or, between them, with the chemistry's eddy diffusion profile
continued upward, so that separation sets in where molecular diffusion
overtakes it (the homopause). Oxygen made above the base, its 63-um cooling,
ions and non-thermal escape are not included, and where the separated fraction
at the exobase is large the trace-gas treatment no longer holds.

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
import re
import json
import math
from pathlib import Path
import sys
import numpy as np

from shared.constants import AVOGADRO, BOLTZMANN, MOON_GM, MOON_RADIUS
from atmosphere.thermal_column import ColumnConfig, lower_boundary, solve_column
from atmosphere.radiative_convective import thermodynamics as th
from atmosphere.middle_atmosphere import chemistry as ch

HERE = Path(__file__).resolve().parent
HEATS = (0.0, 1e-6, 3e-6, 1e-5)
SOLAR_MAXIMUM = 2.5           # rough ratio of solar-maximum to WHI 2008 irradiance below 175 nm
XRAY_SOLAR_MAXIMUM = 100.0    # the same below 10 nm, where X-rays vary far more (WHI's own active week: ~80x at 0.25 nm)
FILM_PRODUCT = Path(__file__).resolve().parents[2] / 'protection' / 'spectra' / 'stack_short_wave.json'
BASELINE_K = 180.0


def _profile(results, case):
    with open(results / 'profiles' / f'{case}.csv', newline='') as handle:
        rows = list(csv.DictReader(handle))
    return {k: np.array([float(r[k]) if r[k] != '' else np.nan for r in rows]) for k in rows[0]}


M_O = 0.016 / AVOGADRO        # kg per oxygen atom
B_O_AIR = 9.69e16             # O through N2: molecular diffusion D = b/n, b = 9.69e16 T^0.774 cm^-1 s^-1 (Banks & Kockarts 1973)
LEAK = 1e-3                   # transmission of the idealised edge filters below their edge (the shield product)
HEATING_EFFICIENCY = 0.4      # neutral heating per absorbed photon energy below 175 nm; 0.3-0.5 in the literature
O2_LYMAN_ALPHA_CM2 = 1.0e-20  # O2 absorption cross-section at 121.6 nm


def whi_quiet_sun():
    """WHI 2008 reference spectrum, quiet-Sun period: bin centres (nm) and irradiance (W m^-2 nm^-1, 0.1-nm bins)."""
    from atmosphere.middle_atmosphere import fetch_inputs
    rows = []
    for line in fetch_inputs.path('whi2008_ref_solar_irradiance_ver2.dat').read_text().splitlines():
        parts = line.split()
        if len(parts) >= 5:
            try:
                rows.append([float(v) for v in parts[:4]])
            except ValueError:
                pass
    data = np.array(rows)
    return data[:, 0], data[:, 3]


def leakage_heat(transmission=LEAK, base_pa=0.3, activity=1.0):
    """Heat (W/m^2 of lunar surface, global mean) from sunlight below 175 nm that a filter passing
    `transmission` there lets in, deposited above the base pressure.

    Irradiance from the WHI 2008 reference spectrum (Woods et al. 2009, near solar minimum; `activity`
    scales it). Everything shortward of 121 nm and the Schumann-Runge continuum (122.5-175 nm) is
    absorbed above the base; Lyman-alpha is absorbed above it in the fraction the O2 column allows along
    a mean slant path (mu = 0.5); the Schumann-Runge bands (175-200 nm) mostly reach below it and are left
    out. A heating efficiency of 0.4 converts absorbed energy to heat of the neutral gas.
    """
    w, f = whi_quiet_sun()
    band = lambda lo, hi: float(f[(w >= lo) & (w < hi)].sum() * 0.1)
    o2_column = 0.175 * base_pa / (0.0289 / AVOGADRO * MOON_GM / MOON_RADIUS ** 2) * 1e-4   # cm^-2
    lyman_fraction = -np.expm1(-O2_LYMAN_ALPHA_CM2 * o2_column / 0.5)
    absorbed = band(0.0, 121.0) + lyman_fraction * band(121.0, 122.5) + band(122.5, 175.0)
    return float(transmission * activity * 0.25 * HEATING_EFFICIENCY * absorbed)


def film_heat(activity=1.0, xray_activity=1.0):
    """Heat (W/m^2 of lunar surface, global mean) from sunlight below 175 nm that the titania stack's film
    itself transmits, from the protection domain's short-wave product. All of it is counted as absorbed
    above the base, an upper bound: the film passes only hard X-rays, and the hardest reach below it.
    `xray_activity` scales the irradiance below 10 nm and `activity` the rest."""
    product = json.loads(FILM_PRODUCT.read_text())
    w, f = whi_quiet_sun()
    sel = w < 175.0
    t = np.interp(w[sel], product['wavelength_nm'], product['transmission'])
    scale = np.where(w[sel] < 10.0, xray_activity, activity)
    return float(0.25 * HEATING_EFFICIENCY * (t * f[sel] * scale).sum() * 0.1)

def atomic_oxygen_loss(column_profile, base_fraction, air_molar_kg=0.0289, eddy_cm2_s=None):
    """Jeans escape of atomic oxygen (kg/s) carried as a trace gas on a thermal-column solution.

    Above the base the oxygen fraction f follows d ln f/dr = D/(D + K) (m_air - m_O) g / (k T), the
    zero-flux balance of eddy diffusion K and oxygen's molecular diffusion D through the air: 'mixed'
    takes K >> D throughout, 'separated' K = 0, and 'eddy' K from `eddy_cm2_s(p)` (cm^2/s at pressure
    p in Pa). Returns losses (kg/s), exobase fractions, oxygen's Jeans parameter at the exobase and,
    with an eddy profile, the homopause pressure where D = K.
    """
    r = column_profile['radius_R'] * MOON_RADIUS
    t = column_profile['temperature_K']
    p = column_profile['pressure_Pa']
    n_total = p / (BOLTZMANN * t)
    g = MOON_GM / r ** 2
    enrichment = (air_molar_kg / AVOGADRO - M_O) * g / (BOLTZMANN * t)       # d ln f/dr with no eddy mixing
    molecular = B_O_AIR * t ** 0.774 / (n_total * 1e-6)                        # cm^2/s
    weights = dict(mixed=np.zeros_like(t), separated=np.ones_like(t))
    out = {}
    if eddy_cm2_s is not None:
        eddy = eddy_cm2_s(p)
        weights['eddy'] = molecular / (molecular + eddy)
        above = np.where(molecular >= eddy)[0]
        out['homopause_pa'] = float(p[above[0]]) if above.size else float('nan')
    t_exo, r_exo = t[-1], r[-1]
    lam = MOON_GM * M_O / (BOLTZMANN * t_exo * r_exo)
    per_density = np.sqrt(BOLTZMANN * t_exo / (2 * np.pi * M_O)) * (1 + lam) * np.exp(-lam) * 4 * np.pi * r_exo ** 2 * M_O
    out['jeans_lambda'] = float(lam)
    for name, w in weights.items():
        integrand = w * enrichment
        log_ratio = float(np.sum(0.5 * (integrand[1:] + integrand[:-1]) * np.diff(r)))
        fraction = min(base_fraction * math.exp(min(log_ratio, 700.0)), 1.0)
        out[f'loss_{name}_kg_s'] = float(fraction * n_total[-1] * per_density)
        out[f'fraction_exobase_{name}'] = fraction
    return out


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
    # The *_ts<T> cases exist for balance temperatures; their upper air is not an alternative exobase base.
    cases = [r for r in csv.DictReader(open(table, newline=''))
             if r['planet'] == 'moon' and not re.search(r'_ts\d+$', r['case'])]
    check_path = args.results / 'lbl_check.json'
    checked = json.loads(check_path.read_text())['cases'] if check_path.is_file() else {}
    out = []
    for row in cases:
        prof = _profile(args.results, row['case'])
        dry = float(row['dry_pressure_pa'])
        o2 = th.earthlike_air(dry, 400.0).fractions['O2']
        o2_molecular = o2 / (o2 + th.earthlike_air(dry, 400.0).fractions['N2'])
        o_top = float(prof['O'][-1])
        o_at = lambda p_b: float(np.exp(np.interp(np.log(p_b), np.log(prof['p_pa'])[::-1],
                                                  np.log(np.maximum(prof['O'], 1e-30))[::-1])))
        bases = [('middle_atmosphere', 0.3, base_temperature(prof, 0.3)),
                 ('middle_atmosphere_top', 0.1, base_temperature(prof, 0.1)),
                 ('middle_atmosphere', 10.0, base_temperature(prof, 10.0)),
                 ('baseline_180K', 0.1, BASELINE_K)]
        if row['case'] in checked:
            c = checked[row['case']]
            bases.append(('line_by_line_corrected_top', 0.1, c['t_lbl_corrected_k'][c['pressures_pa'].index(0.1)]))
        leak_quiet, leak_max = leakage_heat(), leakage_heat(activity=SOLAR_MAXIMUM)
        film_quiet, film_max = film_heat(), film_heat(SOLAR_MAXIMUM, XRAY_SOLAR_MAXIMUM)
        mixing = ch.Mixing(scale=float(row['kzz_scale']))
        eddy = lambda p, m=mixing, pt=float(row['tropopause_pa']): m.profile(p, pt)
        for label, base_p, t_base in bases:
            leaks = (leak_quiet, leak_max) if base_p == 0.3 and row['shield'] != 'none' else ()
            films = (film_quiet, film_max) if base_p == 0.3 and row['shield'] == 'titania_stack' else ()
            for q in HEATS + leaks + films:
                    cfg = ColumnConfig(surface_pressure_pa=float(row['surface_pressure_pa']),
                                       surface_temperature_k=float(row['surface_temperature_k']),
                                       lower_temperature_k=t_base, lower_pressure_pa=base_p,
                                       oxygen_mole_fraction=o2_molecular)
                    source = ('leak 0.1%, quiet Sun' if q == leak_quiet else
                              'leak 0.1%, solar maximum' if q == leak_max else
                              'titania film, quiet Sun' if q == film_quiet else
                              'titania film, solar maximum' if q == film_max else 'swept')
                    rec = dict(case=row['case'], shield=row['shield'], dry_pressure_pa=dry, base=label,
                               base_pressure_pa=base_p, base_temperature_k=round(t_base, 2),
                               deposited_heat_w_m2=q, heat_source=source, atomic_o_fraction_at_top=o_top)
                    try:
                        summary, column_profile, _ = solve_column(q, cfg)
                        f_o = o_at(base_p)
                        o = atomic_oxygen_loss(column_profile, f_o, lower_boundary(cfg)[3], eddy)
                        rec.update(status='converged', exobase_temperature_k=summary['exobase_temperature_k'],
                                   exobase_radius_R=summary['exobase_radius_R'],
                                   molecular_loss_kg_s=summary['molecular_loss_kg_s'],
                                   jeans_lambda_N2=summary['jeans_lambda_N2'],
                                   atomic_o_fraction_at_base=f_o, jeans_lambda_O=o['jeans_lambda'],
                                   atomic_o_loss_mixed_kg_s=o['loss_mixed_kg_s'],
                                   atomic_o_loss_eddy_kg_s=o['loss_eddy_kg_s'],
                                   atomic_o_loss_separated_kg_s=o['loss_separated_kg_s'],
                                   atomic_o_fraction_at_exobase_eddy=o['fraction_exobase_eddy'],
                                   atomic_o_fraction_at_exobase_separated=o['fraction_exobase_separated'],
                                   homopause_pa=o['homopause_pa'],
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
