"""The size-dependent fall law of lunar rain for CM1's Morrison scheme, and the lunar rate of rain's self-collection.

    python -m atmosphere.rain.fall_law        # about a second on one thread
    # -> atmosphere/rain/results/rain_fall_law.json

Morrison's scheme gives rain one fall law, V = A (RHOSU / rho)^0.54 D^B with A = 842 m^0.2/s and B = 0.8 on Earth, and
CM1's lunar builds slow it at every size by (g / 9.81)^((B + 1) / 3), 2.94. The single-drop solver slows each size by
its own factor, from 5.3 at 0.1 mm to 1.9 at 6 mm (rain_drops.json, cm1_check). The scheme reads its law only through
three moments of rain's exponential spectrum N0 exp(-lambda D): the mass-weighted and the number-weighted mean fall
speeds, A Gamma(4 + B) / 6 lambda^-B and A Gamma(1 + B) lambda^-B, which carry its water and its drops down, and the
integral that ventilates its evaporation, A^(1/2) Gamma(5/2 + B/2) lambda^-(5/2 + B/2). No power law follows the solver
through both of its regimes, the viscous fall of drizzle and the inertial fall of large drops, and refitting A and B
leaves errors in these moments as large as the present law's (power_laws in the product). This law takes the three
moments from the solver itself and tabulates them in log lambda over the scheme's range. The speeds are the solver's
terminal speeds in the Moon's air at 2 m in the resort at 301 K, carried below the smallest size solved by Stokes's D^2
and beyond the stability screen at the screen's speed, and divided by Morrison's density factor there, so that the
scheme's own factor gives them back; through the lowest 10 km the Moon's air holds its density within about a tenth.

Rain's self-collection, Seifert and Beheng's (2001) 5.78 m3/kg/s, scales with the gravitational collision kernel
pi/4 (D1 + D2)^2 |V1 - V2|. Its lunar rate is Earth's over the ratio of the kernels in Earth's air and the Moon's,
summed over pairs of the spectrum at the lunar breakup equilibrium (rain_breakup.json) from the scheme's smallest rain,
80 um, to Earth's stability screen, the sizes where both worlds hold drops. Autoconversion and accretion stay Earth's
fits.
"""
from __future__ import annotations

import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from shared.constants import (CM1_GRAVITY, CM1_RAIN_FALL_CAP, CM1_RAIN_FALL_COEFFICIENT, CM1_RAIN_FALL_EXPONENT,
                              CM1_RAIN_MEAN_DIAMETER_FLOOR, CM1_RAIN_SELF_COLLECTION, MOON_SURFACE_GRAVITY)
from shared.provenance import constants_used
from atmosphere.rain import drops, product as rain_product, targets as climate_targets

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCHEMA = 'terluna.atmosphere.rain-fall-law/1'
OUT = HERE / 'results' / 'rain_fall_law.json'
BREAKUP = HERE / 'results' / 'rain_breakup.json'
CODE = ('atmosphere/rain/fall_law.py', 'atmosphere/rain/drops.py', 'atmosphere/rain/product.py',
        'atmosphere/rain/targets.py', 'atmosphere/radiative_convective/thermodynamics.py')
TARGET = 'resort'                                   # the air the law is fitted in, at 2 m
SOLVED_POINTS, SMALLEST_SOLVED_M, SCREEN_MARGIN = 400, 1.0e-6, 0.998
QUADRATURE_POINTS, QUADRATURE_M = 8000, (1.0e-8, 1.0)   # the log grid the spectra are integrated over
TABLE_LAMBDA = (100.0, 1.0e5, 301)                 # 1/m, uniform in log: mean diameters from 10 mm to 10 um
COMPARISON_MM = (0.2, 3.0, 29)                     # mean diameters (mm) where the laws are compared
SMALLEST_RAIN_M = 80.0e-6                          # Seifert and Beheng's rain, drops above 2.6e-10 kg
KERNEL_POINTS = 600
KERNEL_MEAN_MM = (0.6, 1.0, 2.5)                   # beside the lunar equilibrium, where the kernel ratio is shown
VENTILATION = (0.78, 0.308)                        # Morrison's F1R and F2R: f = 0.78 + 0.308 Sc^(1/3) Re^(1/2)

EVIDENCE = ('A runnable calculation: the single-drop solver\'s terminal speeds, which reproduce Earth\'s measured fall '
            'speeds to 0.20 m/s, in the Moon\'s air at 2 m in the resort, integrated over Morrison\'s exponential '
            'spectra. The quadrature reproduces the closed forms of a power law to within 3e-6. The self-collection '
            'rate carries Seifert and Beheng\'s Earth constant by the ratio of geometric collision kernels with '
            'collision efficiency 1; no lunar collision has been measured, and no CM1 run has tested the law.')
READING_RULE = ('table holds, on ln lambda from ln_lambda_first in steps of ln_lambda_step (lambda in 1/m), the natural '
                'logs of the mass-weighted and number-weighted mean fall speeds (m/s) and of the ventilation integral '
                '(the integral of D^(3/2) V^(1/2) exp(-lambda D) dD, m3 s^-1/2) of rain\'s exponential spectrum in '
                'Morrison\'s reference air; the scheme multiplies the speeds by its density factor and the integral by '
                'its square root, and interpolates linearly in ln lambda, holding the ends. comparison gives, by mean '
                'diameter 1/lambda, the present law\'s moments and the best power law\'s over the solver\'s. '
                'self_collection gives the kernel ratio, Earth\'s over the Moon\'s, and the lunar rate. namelist holds '
                'the settings CM1\'s rain-law patch (climate/crm/cm1_run.py) takes.')
ASSUMPTIONS = [
    'The speeds are the solver\'s at 2 m in the resort\'s air at 301 K and 1.2 atm. Morrison\'s density factor '
    '(RHOSU / rho)^0.54 carries them to the air above, whose density through the lowest 10 km stays within about a '
    'tenth of the ground\'s.',
    'Below 1 um the speeds follow Stokes\'s D^2, which the solver meets there to 0.1%; beyond the stability screen, '
    'where Morrison\'s untruncated spectrum still places mass, they hold at the screen\'s speed, the largest drop\'s.',
    'Self-collection takes the geometric kernel with collision efficiency 1 between 80 um and Earth\'s stability '
    'screen, in each world\'s air at 2 m; its breakup offset is the rain limits\'.',
]


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def solved_speeds(column):
    """The solver's still-air terminal speeds at 2 m on a log grid from 1 um to just below the stability screen."""
    screen = drops.stable_diameter(column)
    diameter = np.geomspace(SMALLEST_SOLVED_M, SCREEN_MARGIN * screen, SOLVED_POINTS)
    return diameter, np.array([drops.terminal_speed(d, column) for d in diameter]), screen


def speed_law(diameter, speed):
    """Speed (m/s) at any diameter (m): the solved speeds interpolated in log-log, Stokes's D^2 below the smallest and
    the largest's speed beyond."""
    ld, lv = np.log(diameter), np.log(speed)

    def law(d):
        x = np.log(np.asarray(d, float))
        return np.exp(np.where(x < ld[0], lv[0] + 2.0 * (x - ld[0]), np.interp(x, ld, lv)))

    return law


def power_law(a, b):
    return lambda d: a * np.asarray(d, float)**b


GRID = np.geomspace(*QUADRATURE_M, QUADRATURE_POINTS)
LOG_STEP = math.log(QUADRATURE_M[1] / QUADRATURE_M[0]) / (QUADRATURE_POINTS - 1)


def moments(speed, lam, chunk=40):
    """Morrison's three moments of a fall law speed(D) for rain's exponential spectrum of slope lam (1/m): the
    mass-weighted and the number-weighted mean speeds (m/s) and the ventilation integral of D^(3/2) V^(1/2)
    exp(-lam D) dD (m3 s^-1/2), by the trapezoid rule in ln D."""
    lam = np.atleast_1d(np.asarray(lam, float))
    v = speed(GRID)
    parts = (GRID**3 * v, v, GRID**1.5 * np.sqrt(v))
    out = np.empty((3, lam.size))
    for start in range(0, lam.size, chunk):
        part = lam[start:start + chunk, None]
        weight = np.exp(-part * GRID) * GRID * LOG_STEP
        weight[:, (0, -1)] *= 0.5
        mass, number, vent = (weight @ f for f in parts)
        out[:, start:start + chunk] = (mass * part[:, 0]**4 / 6.0, number * part[:, 0], vent)
    return out


def closed_forms(a, b, lam):
    """The same moments of the power law a D^b, in closed form."""
    lam = np.asarray(lam, float)
    return np.array([a * math.gamma(4 + b) / 6 * lam**-b, a * math.gamma(1 + b) * lam**-b,
                     math.sqrt(a) * math.gamma(2.5 + b / 2) * lam**-(2.5 + b / 2)])


def table_lambda():
    first, last, points = TABLE_LAMBDA
    return np.geomspace(first, last, points)


def interpolate(table, lam):
    """The table's moments at lam, linear in ln lambda and held at its ends, as the patched scheme reads it."""
    x = (np.log(lam) - table['ln_lambda_first']) / table['ln_lambda_step']
    rows = np.array([table['ln_mass_weighted'], table['ln_number_weighted'], table['ln_ventilation']])
    i = np.clip(np.floor(x).astype(int), 0, rows.shape[1] - 2)
    w = np.clip(x - i, 0.0, 1.0)
    return np.exp((1.0 - w) * rows[:, i] + w * rows[:, i + 1])


def best_power_law(reference, lam):
    """The power law a D^b whose mass-weighted speed and ventilation integral come closest to the reference moments
    at lam, by least squares in their logs."""
    target = np.log(reference[(0, 2), :])

    def misfit(p):
        return float(np.sum((np.log(closed_forms(math.exp(p[0]), p[1], lam)[(0, 2), :]) - target)**2))

    start = [math.log(CM1_RAIN_FALL_COEFFICIENT * rain_product.cm1_gravity_factor(MOON_SURFACE_GRAVITY)),
             CM1_RAIN_FALL_EXPONENT]
    fit = minimize(misfit, start, method='Nelder-Mead', options=dict(xatol=1e-8, fatol=1e-12, maxiter=4000))
    return math.exp(fit.x[0]), float(fit.x[1])


def kernel_sum(speed, mean_m, smallest_m, largest_m, points=KERNEL_POINTS):
    """Collisions per unit time and volume among an exponential spectrum's drops, per N0^2: the geometric kernel
    pi/4 (D1 + D2)^2 |V1 - V2| over pairs of sizes from smallest_m to largest_m."""
    d = np.geomspace(smallest_m, largest_m, points)
    w = np.exp(-d / mean_m) * np.gradient(d)
    v = speed(d)
    kernel = np.pi / 4 * (d[:, None] + d[None, :])**2 * np.abs(v[:, None] - v[None, :])
    return 0.5 * float(w @ kernel @ w)


def ratios(found, reference):
    return {k: [float(x) for x in found[i] / reference[i]]
            for i, k in enumerate(('mass_weighted', 'number_weighted', 'ventilation'))}


def build():
    data, binding = climate_targets.load()
    air = rain_product.target_air(climate_targets.read(data))[TARGET]
    cols = rain_product.columns(air)
    moon, earth = cols['moon'], cols['earth']
    state = moon.state(0.0)
    density_factor = rain_product.cm1_density_factor(state['density_kg_m3'])
    d, v, screen = solved_speeds(moon)
    lunar = speed_law(d, v)
    reference = speed_law(d, v / density_factor)
    lam = table_lambda()
    table_moments = moments(reference, lam)
    first, last, points = TABLE_LAMBDA
    table = dict(ln_lambda_first=math.log(first), ln_lambda_step=math.log(last / first) / (points - 1), points=points,
                 ln_mass_weighted=[float(x) for x in np.log(table_moments[0])],
                 ln_number_weighted=[float(x) for x in np.log(table_moments[1])],
                 ln_ventilation=[float(x) for x in np.log(table_moments[2])])
    # The checks: the quadrature against a power law's closed forms; the interpolation between the table's points.
    present = (CM1_RAIN_FALL_COEFFICIENT * rain_product.cm1_gravity_factor(MOON_SURFACE_GRAVITY), CM1_RAIN_FALL_EXPONENT)
    quadrature = moments(power_law(*present), lam) / closed_forms(*present, lam) - 1.0
    middle = np.sqrt(lam[1:] * lam[:-1])
    between = interpolate(table, middle) / moments(reference, middle) - 1.0
    stokes = ((drops.WATER_DENSITY - state['density_kg_m3']) * MOON_SURFACE_GRAVITY / (18.0 * state['viscosity_pa_s'])
              / density_factor)
    smallest_lam = np.array([last])
    stokes_number = float(moments(reference, smallest_lam)[1][0] / (stokes * 2.0 / last**2))
    # The comparison over the mean diameters lunar rain occupies.
    mean_m = 1e-3 * np.geomspace(*COMPARISON_MM)
    solver = moments(reference, 1.0 / mean_m)
    fit = best_power_law(solver, 1.0 / mean_m)
    laws = dict(present=present, best=fit)
    comparison = dict(mean_diameter_mm=[float(x) for x in mean_m * 1e3],
                      solver=dict(mass_weighted_m_s=[float(x) for x in solver[0] * density_factor],
                                  number_weighted_m_s=[float(x) for x in solver[1] * density_factor]),
                      **{name: dict(a_si=a, b=b, over_solver=ratios(closed_forms(a, b, 1.0 / mean_m), solver))
                         for name, (a, b) in laws.items()})
    power_laws = {name: {k: [min(r), max(r)] for k, r in row['over_solver'].items()}
                  for name, row in comparison.items() if name in laws}
    # The lunar breakup equilibrium and the rain limits' cap on the mean fall speeds.
    breakup = json.loads(BREAKUP.read_text())
    equilibrium_m = breakup['lunar_limits']['moon']['equilibrium_mm'] * 1e-3
    at_equilibrium = moments(reference, np.array([1.0 / equilibrium_m]))[:, 0]
    present_equilibrium = closed_forms(*present, np.array([1.0 / equilibrium_m]))[:, 0]
    # Water lost per metre of fall, the evaporation rate over the mass-weighted speed, with Morrison's ventilation
    # in the air at 2 m.
    schmidt = state['viscosity_pa_s'] / (state['density_kg_m3'] * state['vapour_diffusivity_m2_s'])
    ventilated = (VENTILATION[1] * schmidt**(1 / 3)
                  * math.sqrt(density_factor * state['density_kg_m3'] / state['viscosity_pa_s']))
    plain = VENTILATION[0] * equilibrium_m**2
    evaporation = (plain + ventilated * present_equilibrium[2]) / (plain + ventilated * at_equilibrium[2])
    cap = CM1_RAIN_FALL_CAP * (MOON_SURFACE_GRAVITY / CM1_GRAVITY)**0.25         # the rain limits' VCAPR
    floor_lam = 1.0 / (breakup['lunar_limits']['moon']['mean_diameter_cap_mm'] * 1e-3)
    widest = moments(reference, np.array([floor_lam]))[:, 0]
    reach = np.geomspace(0.5e-3, breakup['lunar_limits']['moon']['mean_diameter_cap_mm'] * 1e-3, 200)
    capped = reach[moments(reference, 1.0 / reach)[0] >= cap]
    # Self-collection: Earth's kernel over the Moon's for the same drops.
    ed, ev, earth_screen = solved_speeds(earth)
    earth_speed = speed_law(ed, ev)
    kernel = {f'{m:.2f}': kernel_sum(earth_speed, m * 1e-3, SMALLEST_RAIN_M, earth_screen)
              / kernel_sum(lunar, m * 1e-3, SMALLEST_RAIN_M, earth_screen)
              for m in (*KERNEL_MEAN_MM, equilibrium_m * 1e3)}
    adopted = kernel[f'{equilibrium_m * 1e3:.2f}']
    lunar_rate = round(CM1_RAIN_SELF_COLLECTION / adopted, 4)
    sizes_mm = (0.01, 0.1, 0.2, 0.5, 1.0, 2.0, 3.0, 6.0, 10.0)
    product = dict(
        schema=SCHEMA, evidence=EVIDENCE, reading_rule=READING_RULE, assumptions=ASSUMPTIONS,
        air=dict(target=TARGET, height='2 m', temperature_k=state['temperature_k'], pressure_pa=state['pressure_pa'],
                 density_kg_m3=state['density_kg_m3'], density_factor=density_factor,
                 stability_screen_mm=screen * 1e3, screen_speed_m_s=float(v[-1])),
        speeds=dict(diameter_mm=list(sizes_mm),
                    solver_m_s=[float(lunar(x * 1e-3)) for x in sizes_mm],
                    present_law_m_s=[float(power_law(*present)(x * 1e-3) * density_factor) for x in sizes_mm]),
        table=table,
        checks=dict(power_law_quadrature_max_error=float(np.abs(quadrature).max()),
                    interpolation_max_error=float(np.abs(between).max()),
                    stokes_number_weighted_at_smallest_mean=stokes_number),
        comparison=comparison, power_laws=power_laws,
        lunar_equilibrium=dict(mean_diameter_mm=equilibrium_m * 1e3,
                               solver_mass_weighted_m_s=float(at_equilibrium[0] * density_factor),
                               solver_number_weighted_m_s=float(at_equilibrium[1] * density_factor),
                               present_over_solver=dict(zip(('mass_weighted', 'number_weighted', 'ventilation'),
                                                            [float(x) for x in present_equilibrium / at_equilibrium])),
                               present_over_solver_evaporation=float(evaporation),
                               present_over_solver_loss_per_metre=float(evaporation * at_equilibrium[0]
                                                                        / present_equilibrium[0]),
                               breakup_product_solver_mass_weighted_m_s=breakup['lunar_limits']['fall_speed'][
                                   'at_2m_air']['solver_mass_weighted_speed_m_s']),
        fall_speed_cap=dict(reference_air_m_s=cap, at_2m_m_s=cap * density_factor,
                            widest_mean_diameter_mm=1e3 / floor_lam,
                            mass_weighted_there_m_s=float(widest[0] * density_factor),
                            capped_from_mean_diameter_mm=float(capped[0] * 1e3) if capped.size else None),
        self_collection=dict(earth_rate_m3_kg_s=CM1_RAIN_SELF_COLLECTION, smallest_mm=SMALLEST_RAIN_M * 1e3,
                             largest_mm=earth_screen * 1e3,
                             kernel_ratio_by_mean_diameter_mm=kernel,
                             adopted_mean_diameter_mm=equilibrium_m * 1e3, kernel_ratio=adopted,
                             lunar_rate_m3_kg_s=lunar_rate),
        namelist=dict(rain_fall_law=1, rain_selfcoll=lunar_rate),
        scheme_range=dict(mean_diameter_floor_mm=CM1_RAIN_MEAN_DIAMETER_FLOOR * 1e3,
                          mean_diameter_cap_mm=breakup['lunar_limits']['moon']['mean_diameter_cap_mm']),
        producer=dict(domain='atmosphere', files={p: sha256(ROOT / p) for p in CODE},
                      **dict(binding, inputs={**binding.get('inputs', {}),
                                              str(BREAKUP.relative_to(ROOT)): sha256(BREAKUP)}),
                      constants=constants_used([ROOT / p for p in CODE])))
    return product


def main(argv=None) -> int:
    product = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(product, indent=1, allow_nan=False) + '\n')
    checks, eq, sc = product['checks'], product['lunar_equilibrium'], product['self_collection']
    print(f"checks: quadrature {checks['power_law_quadrature_max_error']:.1e}, interpolation "
          f"{checks['interpolation_max_error']:.1e}, Stokes {checks['stokes_number_weighted_at_smallest_mean']:.4f}")
    for name, row in product['power_laws'].items():
        print(f"{name} law over the solver, mean diameters {COMPARISON_MM[0]}-{COMPARISON_MM[1]} mm: "
              + ', '.join(f'{k} {lo:.3f}-{hi:.3f}' for k, (lo, hi) in row.items()))
    print(f"at the lunar equilibrium ({eq['mean_diameter_mm']:.2f} mm): mass-weighted {eq['solver_mass_weighted_m_s']:.3f} "
          f"m/s, present law over the solver {eq['present_over_solver']}")
    print(f"self-collection: kernel ratio {sc['kernel_ratio']:.3f} ({sc['kernel_ratio_by_mean_diameter_mm']}), "
          f"lunar rate {sc['lunar_rate_m3_kg_s']} m3/kg/s")
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
