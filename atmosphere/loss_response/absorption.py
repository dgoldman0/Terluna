"""How far out the optical shield must reach, and what the sunlit exosphere beyond it produces.

    python -m atmosphere.loss_response.absorption   # after model.py; writes results/absorption_radius.json

The shield casts a shadow cylinder of radius R (its protected radius) along the
Sun-Moon line. Sunlight outside the cylinder reaches the upper air in two ways.

Heating of the thermosphere. A ray at impact parameter b < r_exo crosses the
collisional air below the exobase and leaves the fraction 1 - exp(-tau(b)) of
its ultraviolet there, tau being the band cross-section times the ray's column
below the exobase. Per unit area of the lunar disc this deposit equals the leak

    f_limb(R) = sum over bands of w_band * integral from R to r_exo of
                (1 - exp(-tau_band(b))) 2 b db / R_moon^2,

where w_band is the band's energy over the ultraviolet energy
escape.leakage_heat counts as absorbed above the base, so f_limb compares
directly with model.py's allowed leak fractions. Bands (WHI 2008): below 121 nm
with 2e-21 m^2 per molecule (N2, O2 and O photoabsorption is 1-3e-21 m^2 over
30-100 nm); Lyman-alpha with O2's 1e-24 m^2; and the Schumann-Runge continuum,
122.5-175 nm, with 3e-22 m^2 per O2 molecule, a flux-weighted middle of a
cross-section that falls from about 1.5e-21 m^2 at 140 nm to 1e-23 m^2 at
175 nm; the last two act on O2, 17.5% of the air. This band-mean curve is
kept as a cross-check of the middle atmosphere's limb tracing, which follows
the whole spectrum ray by ray and finds the same heating within 5% on rays
tangent between the base and the exobase. Since 2026-10-07 the protected radius
that covers the heating comes from that tracing (traced.heating_radius): the
smallest radius at which unfiltered light beyond the aperture gives at most a
tenth of the heat of the state the air settles to (half is also reported), the
state itself traced with the UV transmission over that radius's aperture. The
aperture adds the Sun's angular radius times the shield's distance, 20,000 km
for the ring fleet.

Ionization of the exosphere. Above the exobase the air is collisionless, and
its density is Chamberlain's exosphere with the column's exobase density,
temperature and radius: ballistic molecules, which rise and fall back, and
escaping ones, which the column's Jeans loss already counts; satellite orbits
are left empty. Photons below 80 nm ionize a molecule at the rate J, 2e-21 m^2
times the WHI photon flux there, scaled with activity. A ballistic molecule
ionized in sunlight, outside the shadow cylinder, becomes an ion the solar wind
can carry off in place of a molecule that would have returned. The production
of such ions is

    P_ion(R) = m J * integral over r > max(R, r_exo) of n_ballistic(r) sqrt(1 - R^2/r^2) 4 pi r^2 dr,

the square root being the share of a sphere of radius r outside the cylinder.
P_ion bounds this loss from above. The share that escapes depends on the
plasma interaction, which is not modelled here, and on the Moon's passages
through Earth's magnetotail. The two-body exosphere holds only well inside the
Moon's Hill radius (61,500 km, 35.4 lunar radii): Earth's tide, relative to the
Moon's pull, is (r / r_Hill)^3, under 4% at a third of the Hill radius and an
eighth at half. The integral runs to a third of the Hill radius, and half and
the whole are reported to show how much comes from farther out. The density
neglects depletion by ionization along each flight, a correction of order J
times the flight time, 0.05-0.13 for a flight of a day.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.special import gamma, gammainc

from shared.constants import AVOGADRO, BOLTZMANN, MOON_GM, MOON_RADIUS
from atmosphere.thermal_column import solve_column
from atmosphere.middle_atmosphere import escape
from atmosphere.loss_response import model as lr
from atmosphere.loss_response import traced

HERE = Path(__file__).resolve().parent
SCHEMA = 'terluna.atmosphere.absorption-radius/2'
AIR_KG = 0.0289 / AVOGADRO
O2_SHARE = 0.175
HILL_R = 61500e3 / MOON_RADIUS      # the Moon's Hill radius in the Earth's field, lunar radii
SHARES = (0.1, 0.5)
SUN_ANGULAR_RADIUS = 695700.0 / 149597870.7
SEPTEMBER = dict(distance_km=78000.0, protected_radius_R=3.0)
RING_FLEET = dict(distance_km=20000.0, protected_radius_R=4.0)       # decisions.md, 2026-10-07
CROSS_SECTIONS_M2 = dict(xuv=2e-21, lyman_alpha=1e-24 * O2_SHARE, schumann_runge=3e-22 * O2_SHARE)
IONIZATION_CROSS_SECTION_M2 = 2e-21
IONIZING_NM = 80.0
OUTER = dict(third_hill=HILL_R / 3, half_hill=HILL_R / 2, hill=HILL_R)
REPORT_RADII = (3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0)
EVIDENCE = ('The protected radius that covers the heating from the middle atmosphere\'s limb tracing at each budget\'s '
            'allowed UV transmission, with the band-averaged limb curve of the thermal-column profile as a '
            'cross-check; and the photoionization of the ballistic part of a two-body Chamberlain exosphere outside '
            'the shadow cylinder. '
            'The ion production bounds the pickup loss from above: the plasma interaction that decides how many ions '
            "escape is not modelled, and Earth's tides are not included.")
READING_RULE = ('heating["0.1"].radius_R is the smallest protected radius (lunar radii) at which unfiltered light beyond '
                'the aperture gives at most a tenth of the heat of the traced state at the allowed UV transmission; '
                'null where no radius up to 10 does. Its aperture adds the Sun\'s angular radius times the ring fleet\'s '
                'distance; area_over_ring_fleet compares it with the ring fleet\'s aperture (4 R at 20,000 km). The '
                'state, exobase and loss are at the ring fleet\'s 4 R. exosphere_ion_production_kg_s gives the ions '
                'produced outside shadows of several radii, '
                'counted to a third, half or all of the Hill radius, as an upper bound on pickup loss to set beside '
                'the budget; the share that escapes awaits a plasma model.')


def band_weights():
    """Energy per band divided by the energy leakage_heat counts as absorbed above the base."""
    w, f = escape.whi_quiet_sun()
    band = lambda lo, hi: float(f[(w >= lo) & (w < hi)].sum() * 0.1)
    counted = escape.leakage_heat(transmission=1.0) / (0.25 * escape.HEATING_EFFICIENCY)
    return dict(xuv=band(0.0, 121.0) / counted, lyman_alpha=band(121.0, 122.5) / counted,
                schumann_runge=band(122.5, 175.0) / counted)


def ionization_rate(activity='quiet'):
    """Photoionization rate per molecule (s^-1) from WHI photons below 80 nm."""
    w, f = escape.whi_quiet_sun()
    sel = w < IONIZING_NM
    photons = float((f[sel] * 0.1 * w[sel] * 1e-9 / (lr.PLANCK * lr.LIGHT)).sum())
    return IONIZATION_CROSS_SECTION_M2 * photons * lr.ACTIVITY[activity]['uv']


def _lower_gamma(a, x):
    return gamma(a) * gammainc(a, np.maximum(x, 0.0))


def partition(lam, lam_c):
    """Chamberlain's ballistic and escaping partition functions (no satellite particles)."""
    lam = np.asarray(lam, dtype=float)
    psi = lam ** 2 / (lam_c + lam)
    root = np.sqrt(np.maximum(lam_c ** 2 - lam ** 2, 0.0)) / lam_c
    g, gp, full = _lower_gamma(1.5, lam), _lower_gamma(1.5, lam - psi), gamma(1.5)
    ballistic = 2 / math.sqrt(math.pi) * (g - root * np.exp(-psi) * gp)
    escaping = 1 / math.sqrt(math.pi) * (full - g - root * np.exp(-psi) * (full - gp))
    return ballistic, escaping


def exosphere_density(r_R, n_c, t_c, r_c, mass_kg=AIR_KG):
    """Ballistic and escaping densities (m^-3) at radius r (lunar radii) above an exobase at r_c."""
    k = MOON_GM * mass_kg / (BOLTZMANN * t_c * MOON_RADIUS)
    lam, lam_c = k / np.asarray(r_R, dtype=float), k / r_c
    ballistic, escaping = partition(lam, lam_c)
    scale = n_c * np.exp(lam - lam_c)
    return scale * ballistic, scale * escaping


def exobase_state(profile):
    t_c = float(profile['temperature_K'][-1])
    return float(profile['pressure_Pa'][-1]) / (BOLTZMANN * t_c), t_c, float(profile['radius_R'][-1])


def ion_production(profile, activity, radii, outer_R, points=20000):
    """Ions produced (kg/s) from ballistic exospheric molecules in sunlight outside shadows of the given radii."""
    n_c, t_c, r_c = exobase_state(profile)
    r = np.geomspace(r_c, outer_R, points)
    ballistic, _ = exosphere_density(r, n_c, t_c, r_c)
    per_radius = AIR_KG * ionization_rate(activity) * ballistic * 4 * math.pi * r ** 2 * MOON_RADIUS ** 3
    return np.array([np.trapezoid(per_radius * np.sqrt(np.maximum(1.0 - (x / r) ** 2, 0.0)), r)
                     for x in np.atleast_1d(radii)])


def thermosphere_columns(profile, b_R, points=3000):
    """Column (m^-2) along each ray's path below the exobase."""
    r = np.asarray(profile['radius_R'])
    ln_n = np.log(np.asarray(profile['pressure_Pa']) / (BOLTZMANN * np.asarray(profile['temperature_K'])))
    b = np.asarray(b_R, dtype=float)[:, None]
    span = np.sqrt(np.maximum(r[-1] ** 2 - b ** 2, 0.0))
    s = span * np.linspace(0.0, 1.0, points)[None, :] ** 2
    x = np.sqrt(b ** 2 + s ** 2)
    dens = np.exp(np.interp(x, r, ln_n))
    return 2.0 * np.trapezoid(dens, s, axis=1) * MOON_RADIUS


def limb_heating_curve(profile, weights, b_points=1200):
    """Thermospheric equivalent leak f_limb against the protected radius, and each band's tau = 1 radius."""
    r_base, r_exo = float(profile['radius_R'][0]), float(profile['radius_R'][-1])
    b = np.linspace(r_base, r_exo, b_points)
    column = thermosphere_columns(profile, b)
    integrand = np.zeros_like(b)
    tau1 = {}
    for name, sigma in CROSS_SECTIONS_M2.items():
        tau = sigma * column
        integrand += weights[name] * -np.expm1(-tau) * 2 * b
        thick = np.nonzero(tau >= 1.0)[0]
        tau1[name] = float(b[thick[-1]]) if len(thick) else None
    seg = 0.5 * (integrand[1:] + integrand[:-1]) * np.diff(b)
    return b, np.r_[np.cumsum(seg[::-1])[::-1], 0.0], tau1


def smallest_radius(radii, values, target):
    """Smallest radius on the grid (interpolated) where a decreasing quantity is at most target."""
    if values[0] <= target:
        return float(radii[0])
    below = np.nonzero(values <= target)[0]
    if not len(below):
        return None
    i = int(below[0])
    return float(radii[i - 1] + (target - values[i - 1]) * (radii[i] - radii[i - 1]) / (values[i] - values[i - 1]))


def aperture(radius_R, distance_km=RING_FLEET['distance_km']):
    r = radius_R * MOON_RADIUS / 1000 + SUN_ANGULAR_RADIUS * distance_km
    design = RING_FLEET['protected_radius_R'] * MOON_RADIUS / 1000 + SUN_ANGULAR_RADIUS * distance_km
    return dict(aperture_radius_km=r, area_over_ring_fleet=(r / design) ** 2)


def column(shield, treatment, activity, leak, protected_R=None):
    """The thermal column at a UV transmission through gaps: its summary, its profile and the traced heat deposited
    above the base (W/m^2), at a protected radius (default the ring fleet's). Raises ValueError where the air runs
    away."""
    protected_R = lr.PROTECTED_R if protected_R is None else protected_R
    grid, parts = traced.parts(traced.entries(f'{shield}_{treatment}'), traced.FILM, activity,
                               lr.glow_heat(shield, treatment, activity), protected_R)
    q = traced.first_state(grid, parts, leak)
    if q is None:
        raise ValueError(f'the air runs away at a transmission of {leak:g} with {protected_R:g} lunar radii')
    summary, profile, _ = solve_column(q, lr.column_config(shield, treatment))
    return summary, profile, q


def heating_radius(shield, treatment, activity, leak, share=traced.HEATING_SHARE):
    """The smallest protected radius at which light beyond the aperture gives at most `share` of the heat."""
    return traced.heating_radius(traced.entries(f'{shield}_{treatment}'), traced.FILM, activity,
                                 lr.glow_heat(shield, treatment, activity), leak, share)


def case(shield, treatment, activity, leak, budget, weights):
    """One column at the given UV transmission, at the ring fleet's protected radius; the protected radius that
    covers the heating, for each share, from the traced state at each tabulated radius."""
    summary, profile, q = column(shield, treatment, activity, leak)
    b, heat, tau1 = limb_heating_curve(profile, weights)
    heating = {}
    for share in SHARES:
        r = heating_radius(shield, treatment, activity, leak, share)
        heating[f'{share:g}'] = dict(radius_R=r, aperture=aperture(r) if r is not None else None)
    r_heat = heating[f'{SHARES[0]:g}']['radius_R']
    grid = np.linspace(1.5, 30.0, 115)
    production, tenth = {}, {}
    for name, outer in OUTER.items():
        values = ion_production(profile, activity, REPORT_RADII + (r_heat,), outer)
        production[name] = dict(zip([f'{x:g}' for x in REPORT_RADII] + ['heating_radius'], map(float, values)))
        tenth[name] = smallest_radius(grid, ion_production(profile, activity, grid, outer), 0.1 * budget)
    return dict(shield=shield, treatment=treatment, activity=activity, budget_kg_s=budget,
                allowed_leak_fraction=leak, deposited_heat_w_m2=q,
                exobase_radius_R=summary['exobase_radius_R'], exobase_temperature_k=summary['exobase_temperature_k'],
                exobase_density_m3=exobase_state(profile)[0], jeans_lambda_N2=summary['jeans_lambda_N2'],
                molecular_jeans_loss_kg_s=summary['molecular_loss_kg_s'], tau_one_radius_R=tau1,
                limb_heating_leak_at_3R=float(np.interp(3.0, b, heat, right=0.0)),
                heating=heating, exosphere_ion_production_kg_s=production,
                radius_for_production_of_a_tenth_of_budget_R=tenth)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out', type=Path, default=HERE / 'results')
    args = parser.parse_args(argv)
    budgets = json.loads((args.out / 'loss_response.json').read_text())['budgets']
    weights = band_weights()
    rows = []
    fmt = lambda v: f'{v:5.2f}' if v is not None else '  -  '
    for bgt in budgets:
        if bgt['solar_wind'] != 'low' or bgt['allowed_leak_fraction'] is None:
            continue
        row = case(bgt['shield'], bgt['treatment'], bgt['activity'], bgt['allowed_leak_fraction'], bgt['budget_kg_s'], weights)
        rows.append(row)
        ion, tenth = row['exosphere_ion_production_kg_s'], row['radius_for_production_of_a_tenth_of_budget_R']
        heat = row['heating']['0.1']
        print(f"{row['shield']:13s} {row['treatment']:11s} {row['activity']:13s} {row['budget_kg_s']:4g} kg/s: "
              f"exobase {row['exobase_radius_R']:.2f} R, Jeans {row['molecular_jeans_loss_kg_s']:6.2f} kg/s; "
              f"heating radius {fmt(heat['radius_R'])} R; "
              f"ions there {ion['third_hill']['heating_radius']:6.1f} / {ion['half_hill']['heating_radius']:6.1f} / "
              f"{ion['hill']['heating_radius']:6.1f} kg/s, at 6 R {ion['third_hill']['6']:5.1f}; "
              f"tenth of budget at {fmt(tenth['half_hill'])} / {fmt(tenth['hill'])} R")
    product = dict(
        schema=SCHEMA,
        producer=dict(domain='atmosphere', files={p: hashlib.sha256((lr.ROOT / p).read_bytes()).hexdigest()[:16] for p in (
            'atmosphere/loss_response/absorption.py', 'atmosphere/loss_response/model.py',
            'atmosphere/loss_response/traced.py', 'atmosphere/thermal_column.py',
            'atmosphere/middle_atmosphere/results/limb_heat.json')}),
        evidence=EVIDENCE, reading_rule=READING_RULE,
        band_weights=weights, cross_sections_m2=CROSS_SECTIONS_M2,
        ionization_rate_s={a: ionization_rate(a) for a in lr.ACTIVITY},
        outer_radii_R=OUTER, cases=rows)
    (args.out / 'absorption_radius.json').write_text(json.dumps(product, indent=2) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
