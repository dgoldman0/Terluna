"""Loss response of the Open Moon's air to its protection: a screening model.

    python -m atmosphere.loss_response.model      # writes atmosphere/loss_response/results/

How fast the finished 1.2 atm atmosphere loses mass as a function of how much
extreme and far ultraviolet the optical shield lets through and how much of the
solar wind reaches the air. It sizes the two protection functions against a
total loss budget (research/studies/protection_architecture/requirements.md).

Ultraviolet branch. Heat deposited above the 0.3 Pa base comes from light below
175 nm leaking past the shield (a fraction f of the WHI 2008 spectrum, through
escape.leakage_heat), the titania film's own hard X-rays (escape.film_heat), and
the sky's interplanetary Lyman-alpha glow, which no Sun-facing shield blocks.
The thermal column (atmosphere/thermal_column.py) turns that heat into an
exobase and a molecular Jeans loss, with the base temperature the middle
atmosphere finds behind each shield under each treatment of the upper air
(collisional, LTE, all near-infrared heats). Where the column leaves its domain
(no convergence, or a Jeans parameter below 3), escape is bounded above by the
energy-limited expression of the September feasibility report (section 6),
    Mdot = eta * pi * R_abs**3 * F_xuv / (G M),
with the flux below 121 nm reaching the air, eta = 0.01-0.3 and R_abs = 2-3
lunar radii.

Solar-wind branch. A screening range from scalings, not a plasma model. Pickup of
ionized exospheric gas is capped by mass loading: a flow that picks up about its
own mass flux is slowed and diverted, so pickup is 0.1-1 times the solar wind's
mass flux through the interaction cross-section. Sputtering comes from pickup
ions that fall back into the air (a fraction of the pickup, times a yield of
1-10) and from solar-wind protons that reach it (a precipitating fraction times
a yield of 0.01-0.1 N2 molecules per proton). All of it scales with the share of
each orbit the Moon spends in the solar wind rather than in Earth's magnetotail.

Not included: atomic-oxygen escape (about half again behind the 200-nm edge, per
the middle-atmosphere results), photochemical escape, hydrogen from water,
Earth's tidal lowering of the escape barrier (flagged where the exobase passes
3.5 lunar radii), the day-night circulation of the upper air, and plasma physics
beyond the scalings above.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.special import expn

from shared.constants import AVOGADRO, MOON_GM, MOON_RADIUS
from atmosphere.thermal_column import ColumnConfig, solve_column
from atmosphere.radiative_convective import thermodynamics as th
from atmosphere.middle_atmosphere import escape

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MIDDLE = ROOT / 'atmosphere' / 'middle_atmosphere' / 'results'
SCHEMA = 'terluna.atmosphere.loss-response/1'

PLANCK, LIGHT, PROTON_KG = 6.62607015e-34, 299792458.0, 1.67262192e-27
N2_KG = 0.0280134 / AVOGADRO
LYMAN_ALPHA_M = 121.567e-9

BASE_PA = 0.3
SHIELDS = {'titania_stack': 'moon_1.2atm_titania_stack', 'edge_200nm': 'moon_1.2atm_edge_200nm'}
TREATMENTS = {'collisional': '_nonlte_nir_collisional', 'lte': '', 'all_heats': '_nonlte'}
ACTIVITY = {'quiet': dict(uv=1.0, xray=1.0, glow=1.0),
            # escape.py's rough solar-maximum ratios; the glow is resonantly scattered solar
            # Lyman-alpha, whose line strengthens by roughly 1.5-2 over the cycle (1.5 assumed).
            'solar_maximum': dict(uv=escape.SOLAR_MAXIMUM, xray=escape.XRAY_SOLAR_MAXIMUM, glow=1.5)}
LEAKS = (1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 2e-3, 3e-3, 5e-3,
         1e-2, 2e-2, 3e-2, 5e-2, 0.1, 0.3, 1.0)
GLOW_RAYLEIGH = 1000.0          # interplanetary hydrogen glow at 1 AU (protection/report.md, section 9)
BUDGETS_KG_S = (1.0, 10.0, 100.0)
ETAS = (0.01, 0.1, 0.3)
ABSORPTION_RADII = (2.0, 3.0)
TIDAL_REVIEW_R = 3.5            # 0.1 of the lunar Hill radius, as the thermal column flags it

SOLAR_WIND = dict(density_m3=5e6, speed_m_s=4.0e5)    # typical at 1 AU (2-20 cm^-3, 300-800 km/s)
SW_RANGES = dict(pickup_of_cap=(0.1, 1.0),              # share of the mass-loading cap actually picked up
                 reimpact=(0.1, 0.5),                   # share of pickup ions that fall back into the air
                 ion_yield=(1.0, 10.0),                 # molecules ejected per re-impacting ion
                 proton_precipitation=(0.1, 1.0),       # share of solar-wind protons reaching the air
                 proton_yield=(0.01, 0.1),              # N2 molecules ejected per precipitating proton
                 exposure=(0.75, 0.85))                 # share of the orbit outside Earth's magnetotail


def geometric_mean(pair):
    return math.sqrt(pair[0] * pair[1])


def lyman_glow_flux(rayleigh=GLOW_RAYLEIGH):
    """W/m^2 onto a surface facing open sky: 1 R = 1e10/(4 pi) photons m^-2 s^-1 sr^-1, times pi sr."""
    return rayleigh * 1e10 / (4 * math.pi) * math.pi * PLANCK * LIGHT / LYMAN_ALPHA_M


def lyman_glow_heat(base_radius_R, rayleigh=GLOW_RAYLEIGH, base_pa=BASE_PA):
    """Heat (W/m^2 of lunar surface) the all-sky glow deposits above the base.

    The shell at the base radius intercepts the glow over its whole area. O2 absorbs it
    with escape.py's cross-section, over the column above the base taken with gravity at
    the base radius, for isotropic light from above (absorbed share 1 - 2 E3(tau)).
    escape.py's heating efficiency converts it to heat."""
    g_base = MOON_GM / (base_radius_R * MOON_RADIUS) ** 2
    o2_column_cm2 = 0.175 * base_pa / (0.0289 / AVOGADRO * g_base) * 1e-4
    tau = escape.O2_LYMAN_ALPHA_CM2 * o2_column_cm2
    absorbed = 1.0 - 2.0 * float(expn(3, tau))
    return lyman_glow_flux(rayleigh) * base_radius_R ** 2 * absorbed * escape.HEATING_EFFICIENCY


def xuv_flux(activity='quiet'):
    """Solar flux below 121 nm at 1 AU (W/m^2), from the WHI 2008 spectrum escape.py uses."""
    w, f = escape.whi_quiet_sun()
    return float(f[w < 121.0].sum() * 0.1) * ACTIVITY[activity]['uv']


def energy_limited_loss(flux_w_m2, absorption_radius_R, eta):
    """The feasibility report's energy-limited expression (kg/s)."""
    r = absorption_radius_R * MOON_RADIUS
    return eta * math.pi * r ** 3 * flux_w_m2 / MOON_GM


def base_conditions(shield, treatment):
    """Base temperature at 0.3 Pa and the column's surface state for one middle-atmosphere case."""
    case = SHIELDS[shield] + TREATMENTS[treatment]
    rows = [r for r in csv.DictReader(open(MIDDLE / 'middle_atmosphere.csv', newline='')) if r['case'] == case]
    if not rows:
        raise ValueError(f'No middle-atmosphere case {case}')
    row = rows[0]
    profile = escape._profile(MIDDLE, case)
    dry = float(row['dry_pressure_pa'])
    air = th.earthlike_air(dry, 400.0).fractions
    return dict(case=case, base_temperature_k=escape.base_temperature(profile, BASE_PA),
                surface_pressure_pa=float(row['surface_pressure_pa']),
                surface_temperature_k=float(row['surface_temperature_k']),
                oxygen_mole_fraction=air['O2'] / (air['O2'] + air['N2']))


def euv_sweep(shield, treatment, activity='quiet', glow=True, leaks=LEAKS):
    """Rows of deposited heat, exobase and molecular loss against the leak fraction f."""
    base = base_conditions(shield, treatment)
    cfg = ColumnConfig(surface_pressure_pa=base['surface_pressure_pa'],
                       surface_temperature_k=base['surface_temperature_k'],
                       lower_temperature_k=base['base_temperature_k'], lower_pressure_pa=BASE_PA,
                       oxygen_mole_fraction=base['oxygen_mole_fraction'])
    act = ACTIVITY[activity]
    film = escape.film_heat(act['uv'], act['xray']) if shield == 'titania_stack' else 0.0
    unit = escape.leakage_heat(transmission=1.0, activity=act['uv'])
    base_radius = solve_column(0.0, cfg)[0]['lower_radius_R']
    glow_heat = lyman_glow_heat(base_radius) * act['glow'] if glow else 0.0
    flux = xuv_flux(activity)
    rows, previous, in_domain = [], None, True
    for f in leaks:
        q = f * unit + film + glow_heat
        row = dict(shield=shield, treatment=treatment, activity=activity, lyman_glow=glow,
                   leak_fraction=f, base_temperature_k=round(base['base_temperature_k'], 2),
                   deposited_heat_w_m2=q, leak_heat_w_m2=f * unit, film_heat_w_m2=film,
                   glow_heat_w_m2=glow_heat)
        for eta in ETAS:
            for r_abs in ABSORPTION_RADII:
                row[f'energy_limited_eta{eta:g}_R{r_abs:g}_kg_s'] = energy_limited_loss(f * flux, r_abs, eta)
        if in_domain:
            try:
                summary, _, previous = solve_column(q, cfg, previous)
                flags = summary['domain_flags']
                if 'HYDROSTATIC_OUTFLOW_LIMIT_EXCEEDED' in flags:
                    in_domain = False
                else:
                    row.update(status='thermal_column', exobase_temperature_k=summary['exobase_temperature_k'],
                               exobase_radius_R=summary['exobase_radius_R'],
                               jeans_lambda_N2=summary['jeans_lambda_N2'],
                               molecular_loss_kg_s=summary['molecular_loss_kg_s'],
                               tidal_review=summary['exobase_radius_R'] > TIDAL_REVIEW_R)
            except (ValueError, RuntimeError, FloatingPointError):
                in_domain = False
        if not in_domain:
            row.update(status='outside_thermal_column_domain')
        rows.append(row)
    return rows


def loss_estimate(row):
    """The thermal column's loss where it holds; otherwise the central energy-limited bound (eta 0.1, 3 R)."""
    if row['status'] == 'thermal_column':
        return row['molecular_loss_kg_s']
    return row['energy_limited_eta0.1_R3_kg_s']


def solar_wind_losses(interaction_radius_R=3.0, wind=SOLAR_WIND, ranges=SW_RANGES):
    """Low, central and high solar-wind-driven loss (kg/s), with their parts."""
    area = math.pi * (interaction_radius_R * MOON_RADIUS) ** 2
    proton_flux = wind['density_m3'] * wind['speed_m_s']
    cap = proton_flux * PROTON_KG * area
    out = {}
    for label, pick in (('low', lambda p: p[0]), ('central', geometric_mean), ('high', lambda p: p[1])):
        pickup = cap * pick(ranges['pickup_of_cap'])
        ion_sputter = pickup * pick(ranges['reimpact']) * pick(ranges['ion_yield'])
        proton_sputter = (proton_flux * area * pick(ranges['proton_precipitation'])
                          * pick(ranges['proton_yield']) * N2_KG)
        exposure = pick(ranges['exposure'])
        out[label] = dict(pickup_kg_s=pickup * exposure, ion_sputtering_kg_s=ion_sputter * exposure,
                          proton_sputtering_kg_s=proton_sputter * exposure,
                          total_kg_s=(pickup + ion_sputter + proton_sputter) * exposure)
    out['mass_loading_cap_kg_s'] = cap
    out['interaction_radius_R'] = interaction_radius_R
    return out


def allowed_leak(rows, allowance_kg_s):
    """Largest leak fraction whose loss stays within the allowance (log-log interpolation); None if none."""
    if allowance_kg_s <= 0:
        return None
    f = np.array([r['leak_fraction'] for r in rows])
    loss = np.array([loss_estimate(r) for r in rows])
    if loss[0] > allowance_kg_s:
        return None
    above = np.nonzero(loss > allowance_kg_s)[0]
    if not len(above):
        return float(f[-1])
    i = above[0]
    x0, x1 = math.log(f[i - 1]), math.log(f[i])
    y0, y1 = math.log(max(loss[i - 1], 1e-30)), math.log(loss[i])
    return float(math.exp(x0 + (math.log(allowance_kg_s) - y0) * (x1 - x0) / (y1 - y0)))


def budget_table(sweeps, wind):
    out = []
    for key, rows in sweeps.items():
        shield, treatment, activity, glow = key
        if not glow:
            continue
        for budget in BUDGETS_KG_S:
            for level in ('low', 'high'):
                allowance = budget - wind[level]['total_kg_s']
                f = allowed_leak(rows, allowance)
                out.append(dict(shield=shield, treatment=treatment, activity=activity, budget_kg_s=budget,
                                solar_wind=level, solar_wind_kg_s=wind[level]['total_kg_s'],
                                allowance_for_ultraviolet_kg_s=allowance, allowed_leak_fraction=f,
                                floor_kg_s=loss_estimate(rows[0]),
                                verdict=('solar-wind loss alone exceeds the budget' if allowance <= 0 else
                                         'glow and film alone exceed the allowance' if f is None else
                                         'within budget')))
    return out


def compute():
    sweeps = {}
    for shield in SHIELDS:
        for treatment in TREATMENTS:
            for activity in ACTIVITY:
                for glow in (True, False):
                    sweeps[(shield, treatment, activity, glow)] = euv_sweep(shield, treatment, activity, glow)
    wind = solar_wind_losses()
    unshielded = {activity: {f'eta{eta:g}_R{r:g}_kg_s': energy_limited_loss(xuv_flux(activity), r, eta)
                             for eta in ETAS for r in (1.5,) + ABSORPTION_RADII}
                  for activity in ACTIVITY}
    return sweeps, wind, budget_table(sweeps, wind), unshielded


def _digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()[:16]


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out', type=Path, default=HERE / 'results')
    args = parser.parse_args(argv)
    sweeps, wind, budgets, unshielded = compute()
    args.out.mkdir(parents=True, exist_ok=True)
    rows = [r for v in sweeps.values() for r in v]
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with open(args.out / 'euv_sweep.csv', 'w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    product = dict(
        schema=SCHEMA,
        producer=dict(domain='atmosphere', files={p: _digest(ROOT / p) for p in (
            'atmosphere/loss_response/model.py', 'atmosphere/thermal_column.py',
            'atmosphere/middle_atmosphere/escape.py')}),
        evidence=('Screening model. The ultraviolet branch is the molecular thermal column (no infrared cooling, '
                  'no atomic oxygen, no tides) above middle-atmosphere base temperatures, with energy-limited '
                  'upper bounds beyond its domain; the solar-wind branch is a range from mass-loading and '
                  'sputtering scalings with assumed parameter ranges. Neither is a coupled aeronomy or plasma model.'),
        reading_rule=('Loss rates are global-mean kg/s for the finished 1.2 atm atmosphere. Use the thermal-column '
                      'loss where status is thermal_column and the energy-limited range beyond it. A budget row\'s '
                      'allowed_leak_fraction is the largest fraction of sunlight below 175 nm the optical shield may '
                      'let through while ultraviolet-driven loss plus the stated solar-wind loss stays within the budget.'),
        assumptions=dict(base_pressure_pa=BASE_PA, lyman_glow_rayleigh_quiet=GLOW_RAYLEIGH,
                         activity=ACTIVITY, solar_wind=SOLAR_WIND, solar_wind_ranges=SW_RANGES,
                         energy_limited_eta=ETAS, absorption_radii_R=ABSORPTION_RADII,
                         xuv_flux_below_121nm_w_m2={a: xuv_flux(a) for a in ACTIVITY},
                         lyman_glow_flux_w_m2=lyman_glow_flux()),
        solar_wind=wind,
        unshielded_energy_limited_kg_s=unshielded,
        budgets=budgets)
    (args.out / 'loss_response.json').write_text(json.dumps(product, indent=2) + '\n')
    for b in budgets:
        if b['treatment'] == 'lte' and b['solar_wind'] == 'low':
            f = b['allowed_leak_fraction']
            print(f"{b['shield']:14s} {b['activity']:13s} budget {b['budget_kg_s']:5g} kg/s: "
                  f"allowed leak {('%.2e' % f) if f else '-':>9s}  ({b['verdict']})")
    print('solar wind (kg/s): low %.3g, central %.3g, high %.3g' %
          tuple(wind[k]['total_kg_s'] for k in ('low', 'central', 'high')))
    return 0


if __name__ == '__main__':
    sys.exit(main())
