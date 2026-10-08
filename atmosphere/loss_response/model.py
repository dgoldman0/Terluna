"""Loss response of the Open Moon's air to its protection: a screening model.

    python -m atmosphere.loss_response.model      # writes atmosphere/loss_response/results/

How fast the finished 1.2 atm atmosphere loses mass as a function of how much
extreme and far ultraviolet the optical shield lets through and how much of the
solar wind reaches the air. It sizes the two protection functions against a
total loss budget (research/studies/protection_architecture/requirements.md).

Ultraviolet branch. The heat deposited between the 0.3 Pa base and the
exobase comes from the middle atmosphere's limb tracing (traced.py, reading
atmosphere/middle_atmosphere/results/limb_heat.json): the light below 175 nm
that passes the shield's window stack and its 4 um annulus film, a UV
transmission f through gaps over the whole aperture of the protected radius,
and unfiltered sunlight beyond it, each followed along its slant paths through
air swollen by the heat it takes, and the sky's interplanetary Lyman-alpha
glow, which no Sun-facing shield blocks. The author adopted this traced count
on 2026-10-07 in place of a quarter of the light reaching the disk; solar
maximum takes 2.5 on the ultraviolet and FISM2's measured rise of the X-rays
below 10 nm (escape.xray_cycle). The films are the design's behind both
shields. The protected radius is the ring fleet's 4 lunar radii unless given.
The thermal column (atmosphere/thermal_column.py) turns that heat into an
exobase and a molecular Jeans loss, with the base temperature the middle
atmosphere finds behind each shield under each treatment of the upper air
(collisional, LTE, all near-infrared heats). Since 2026-10-08 the column also
radiates in the infrared (infrared.py): CO2's 15-um band at the middle
atmosphere's 400 ppm, separating above its homopause, and the base's atomic
oxygen and NO carried up, each in excess of what air at the base temperature
emits. The column takes each state's heat
where the tracing puts it in height (state_config, traced.shape): each part of
the light by the limb tables' record of where it heats, and the glow where O2
absorbs it, mostly within a few e-folds of pressure above the base. The states
themselves come from tables traced through air the column heats with its
'middle' shape, which swells more, so they and the runaway onsets lean
conservative. Earth's tide raises each species'
escape by the multiplier tides.py finds for its exobase radius and Jeans
parameter (results/tidal_escape.json); the two-body loss is kept beside it. Where the column leaves its domain
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
The exosphere step (exosphere.py) finds that the ions the sunlit exosphere
makes outside the shield's shadow, 36-470 times this cap, are carried off
rather than refused, and it replaces this branch's pickup and ion sputtering
with its own losses of the sunlit exosphere. The allowed transmissions here
leave those losses out; exosphere.py gives them with the losses included.

Not included: atomic-oxygen escape (about half again behind the 200-nm edge, per
the middle-atmosphere results; oxygen.py counts it), the infrared cooling of the
oxygen atoms made above the base, photochemical escape below the exobase (the
exosphere step counts it above), hydrogen from water,
Earth's tide in the energy-limited bound (it would raise that bound by about
1/K, 10-20%), the day-night circulation of the upper air, and plasma physics
beyond the scalings above.
"""
from __future__ import annotations
import argparse
import csv
import dataclasses
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.special import expn

from shared.constants import AVOGADRO, MOON_GM, MOON_RADIUS
from atmosphere.thermal_column import ColumnConfig, lower_boundary, solve_column
from atmosphere.radiative_convective import thermodynamics as th
from atmosphere.middle_atmosphere import escape
from atmosphere.loss_response import infrared, tides, traced

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MIDDLE = ROOT / 'atmosphere' / 'middle_atmosphere' / 'results'
SCHEMA = 'terluna.atmosphere.loss-response/3'
TIDES = HERE / 'results' / 'tidal_escape.json'
CYCLE = HERE / 'results' / 'solar_cycle.json'
TIDES_SCHEMA = 'terluna.atmosphere.tidal-escape/1'

PLANCK, LIGHT, PROTON_KG = 6.62607015e-34, 299792458.0, 1.67262192e-27
N2_KG = 0.0280134 / AVOGADRO
LYMAN_ALPHA_M = 121.567e-9

BASE_PA = 0.3
SHIELDS = {'titania_stack': 'moon_1.2atm_titania_stack', 'edge_200nm': 'moon_1.2atm_edge_200nm'}
TREATMENTS = {'collisional': '_nonlte_nir_collisional', 'lte': '', 'all_heats': '_nonlte'}
ACTIVITY = {'quiet': dict(uv=1.0, xray=1.0, glow=1.0),
            # FISM2's year around the strongest solar maximum from 1978 on, cycle 21's of 1979-80, every band and
            # the glow (Lyman-alpha) from the same year (cycle.py; the author's choice of 2026-10-07).
            'solar_maximum': json.loads(CYCLE.read_text())['solar_maximum'],
            # Two stress cases. The strongest maximum on record, cycle 19's of 1957-58, whose spectrum FISM2 builds
            # from the 10.7 cm radio flux alone; and escape.py's rough ratios, 2.5 on the ultraviolet above 10 nm and
            # FISM2's mean rise of the X-rays over the last three maxima, with 1.5 on the glow, which cycle 19's year
            # passes in the X-rays and the glow.
            'solar_maximum_cycle_19': json.loads(CYCLE.read_text())['record_maximum'],
            'solar_maximum_stress': dict(uv=escape.SOLAR_MAXIMUM, xray=escape.XRAY_SOLAR_MAXIMUM, glow=1.5)}
STRESS_CASES = ('solar_maximum_cycle_19', 'solar_maximum_stress')     # summarised apart from quiet Sun and maximum
LEAKS = (1e-6, 3e-6, 1e-5, 2e-5, 3e-5, 5e-5, 1e-4, 2e-4, 3e-4, 5e-4, 1e-3, 2e-3, 3e-3, 5e-3,
         1e-2, 2e-2, 3e-2, 5e-2, 0.1, 0.3, 1.0)
PROTECTED_R = 4.0               # the ring fleet's protected radius in lunar radii (decisions.md, 2026-10-07)
HEATING_SHAPE = 'traced'        # the thermal column's heating shape at a traced state: from the limb tables (2026-10-07)
INFRARED_COOLING = True         # the thermal column radiates CO2's, O's and NO's infrared (infrared.py, 2026-10-08)
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


def lyman_glow_tau(base_radius_R, base_pa=BASE_PA):
    """O2's vertical optical depth to Lyman-alpha above the base, with gravity at the base radius."""
    g_base = MOON_GM / (base_radius_R * MOON_RADIUS) ** 2
    o2_column_cm2 = 0.175 * base_pa / (0.0289 / AVOGADRO * g_base) * 1e-4
    return escape.O2_LYMAN_ALPHA_CM2 * o2_column_cm2


def lyman_glow_heat(base_radius_R, rayleigh=GLOW_RAYLEIGH, base_pa=BASE_PA):
    """Heat (W/m^2 of lunar surface) the all-sky glow deposits above the base.

    The shell at the base radius intercepts the glow over its whole area. O2 absorbs it
    with escape.py's cross-section, over the column above the base taken with gravity at
    the base radius, for isotropic light from above (absorbed share 1 - 2 E3(tau)).
    escape.py's heating efficiency converts it to heat."""
    tau = lyman_glow_tau(base_radius_R, base_pa)
    absorbed = 1.0 - 2.0 * float(expn(3, tau))
    return lyman_glow_flux(rayleigh) * base_radius_R ** 2 * absorbed * escape.HEATING_EFFICIENCY


def lyman_glow_shape(base_radius_R, x, base_pa=BASE_PA):
    """Cumulative share of the glow's heat below each log pressure x above the base: isotropic light from above,
    absorbed by O2 whose column falls with pressure (transmitted share 2 E3(tau))."""
    tau = lyman_glow_tau(base_radius_R, base_pa)
    floor = 2.0 * float(expn(3, tau))
    return (2.0 * expn(3, tau * np.exp(-np.asarray(x, float))) - floor) / (1.0 - floor)


def activity_of(activity):
    """An activity as a dict: quiet Sun, the measured solar maximum or a stress case by name, or a measured one as
    given (with a factor on each band of the limb tables, 'bands', a 'label', and 'glow', 'uv' and 'xray' factors;
    see cycle.py)."""
    return ACTIVITY[activity] if isinstance(activity, str) else activity


def activity_label(activity):
    return activity if isinstance(activity, str) else activity['label']


def xuv_flux(activity='quiet'):
    """Solar flux below 121 nm at 1 AU (W/m^2), from the WHI 2008 spectrum escape.py uses, at an activity (a
    measured one weights each band)."""
    w, f = escape.whi_quiet_sun()
    sel = w < 121.0
    act = activity_of(activity)
    if 'bands' in act:
        return float((f[sel] * traced.band_factors(w[sel], act)).sum() * 0.1)
    return float(f[sel].sum() * 0.1) * act['uv']


def energy_limited_loss(flux_w_m2, absorption_radius_R, eta):
    """The feasibility report's energy-limited expression (kg/s)."""
    r = absorption_radius_R * MOON_RADIUS
    return eta * math.pi * r ** 3 * flux_w_m2 / MOON_GM


_TIDAL = {}


def tidal_table():
    if 'table' not in _TIDAL:
        table = json.loads(TIDES.read_text())
        if table['schema'] != TIDES_SCHEMA:
            raise ValueError(f"unexpected tidal-escape schema {table['schema']}")
        _TIDAL['table'] = table
    return _TIDAL['table']


def with_tides(summary):
    """Molecular loss with Earth's tide: each species' two-body Jeans loss times its tidal multiplier."""
    table = tidal_table()
    m_n2 = tides.multiplier(table, summary['exobase_radius_R'], summary['jeans_lambda_N2'])
    m_o2 = tides.multiplier(table, summary['exobase_radius_R'], summary['jeans_lambda_O2'])
    return summary['N2_loss_kg_s'] * m_n2 + summary['O2_loss_kg_s'] * m_o2, m_n2


def base_conditions(shield, treatment):
    """Base temperature at 0.3 Pa and the column's surface state for one middle-atmosphere case."""
    case = SHIELDS[shield] + TREATMENTS[treatment]
    rows = [r for r in csv.DictReader(open(MIDDLE / 'middle_atmosphere.csv', newline='')) if r['case'] == case]
    if not rows:
        raise ValueError(f'No middle-atmosphere case {case}')
    row = rows[0]
    profile = escape._profile(MIDDLE, case)
    dry = float(row['dry_pressure_pa'])
    air = th.earthlike_air(dry, infrared.CO2_PPM).fractions
    return dict(case=case, base_temperature_k=escape.base_temperature(profile, BASE_PA),
                surface_pressure_pa=float(row['surface_pressure_pa']),
                surface_temperature_k=float(row['surface_temperature_k']),
                oxygen_mole_fraction=air['O2'] / (air['O2'] + air['N2']))


_COLUMN = {}


def column_config(shield, treatment, cooling=None):
    """The thermal column's configuration for one middle-atmosphere case, radiating in the infrared (infrared.py)
    unless cooling, or INFRARED_COOLING when it is None, says not."""
    cooling = INFRARED_COOLING if cooling is None else cooling
    if (shield, treatment, cooling) not in _COLUMN:
        base = base_conditions(shield, treatment)
        cfg = ColumnConfig(surface_pressure_pa=base['surface_pressure_pa'],
                           surface_temperature_k=base['surface_temperature_k'],
                           lower_temperature_k=base['base_temperature_k'], lower_pressure_pa=BASE_PA,
                           oxygen_mole_fraction=base['oxygen_mole_fraction'])
        if cooling:
            cfg = dataclasses.replace(cfg, **infrared.column_fields(base['case'], BASE_PA, cfg.lower_temperature_k,
                                                                    lower_boundary(cfg)[3]))
        _COLUMN[shield, treatment, cooling] = cfg
    return _COLUMN[shield, treatment, cooling]


_BASE_RADIUS = {}


def base_radius_R(shield, treatment):
    """The column's base radius (lunar radii) for a case, with no heat."""
    if (shield, treatment) not in _BASE_RADIUS:
        _BASE_RADIUS[shield, treatment] = solve_column(0.0, column_config(shield, treatment))[0]['lower_radius_R']
    return _BASE_RADIUS[shield, treatment]


def state_config(shield, treatment, activity, protected_R, gaps, heat, glow=True):
    """The thermal column's configuration at a traced state: the case's, with the heating shape HEATING_SHAPE names;
    'traced' takes each part of the light's heat where the limb tables put it and the glow's where O2 absorbs it."""
    cfg = column_config(shield, treatment)
    if HEATING_SHAPE != 'traced':
        return dataclasses.replace(cfg, heating_shape=HEATING_SHAPE)
    glow_w = glow_heat(shield, treatment, activity) if glow else 0.0
    below = traced.shape(traced.entries(f'{shield}_{treatment}'), traced.FILM, activity, glow_w,
                         lyman_glow_shape(base_radius_R(shield, treatment), traced.SHAPE_X), protected_R, gaps, heat)
    return dataclasses.replace(cfg, heating_shape='traced', heating_log_pressure=tuple(float(x) for x in traced.SHAPE_X),
                               heating_fraction=tuple(float(v) for v in below))


def glow_heat(shield, treatment, activity):
    """The sky's Lyman-alpha glow heat (W/m^2) for a case at an activity, at the column's base radius."""
    return lyman_glow_heat(base_radius_R(shield, treatment)) * activity_of(activity)['glow']


def euv_sweep(shield, treatment, activity='quiet', glow=True, leaks=LEAKS, protected_R=PROTECTED_R):
    """Rows of the traced state, exobase and molecular loss against the UV transmission f through gaps."""
    base = base_conditions(shield, treatment)
    act = activity_of(activity)
    unit = escape.leakage_heat(transmission=1.0, activity=act['uv'], xray_activity=act['xray'])
    glow_w = glow_heat(shield, treatment, activity) if glow else 0.0
    grid, parts = traced.parts(traced.entries(f'{shield}_{treatment}'), traced.FILM, activity, glow_w, protected_R)
    flux = xuv_flux(activity)
    # Where the air runs away inside the sweep, add the transmission at which it starts.
    stable = [f for f in leaks if traced.first_state(grid, parts, f) is not None]
    fs = list(leaks)
    if stable and len(stable) < len(leaks):
        lo, hi = max(stable), min(f for f in leaks if f not in stable)
        for _ in range(40):
            mid = math.sqrt(lo * hi)
            lo, hi = (mid, hi) if traced.first_state(grid, parts, mid) is not None else (lo, mid)
        fs = sorted(fs + [lo])
    rows, previous, in_domain = [], None, True
    for f in fs:
        row = dict(shield=shield, treatment=treatment, activity=activity_label(activity), lyman_glow=glow,
                   protected_radius_R=protected_R, leak_fraction=f,
                   base_temperature_k=round(base['base_temperature_k'], 2), disk_count_heat_w_m2=f * unit)
        for eta in ETAS:
            for r_abs in ABSORPTION_RADII:
                row[f'energy_limited_eta{eta:g}_R{r_abs:g}_kg_s'] = energy_limited_loss(f * flux, r_abs, eta)
        q = traced.first_state(grid, parts, f)
        if q is None:
            row.update(status='runaway', note='the heat the swollen air takes outgrows the limb tables')
            rows.append(row)
            continue
        share = traced.shares(grid, parts, f, q)
        row.update(deposited_heat_w_m2=q, leak_heat_w_m2=share['gaps'], film_heat_w_m2=share['window'] + share['annulus'],
                   glow_heat_w_m2=share['glow'], beyond_aperture_heat_w_m2=share['beyond_aperture'],
                   beyond_share=share['beyond_share'])
        if in_domain:
            try:
                summary, _, previous = solve_column(q, state_config(shield, treatment, activity, protected_R, f, q, glow),
                                                    previous)
                flags = summary['domain_flags']
                if 'HYDROSTATIC_OUTFLOW_LIMIT_EXCEEDED' in flags:
                    in_domain = False
                else:
                    tidal_loss, m_n2 = with_tides(summary)
                    row.update(status='thermal_column', exobase_temperature_k=summary['exobase_temperature_k'],
                               exobase_radius_R=summary['exobase_radius_R'],
                               jeans_lambda_N2=summary['jeans_lambda_N2'],
                               molecular_loss_two_body_kg_s=summary['molecular_loss_kg_s'],
                               tidal_multiplier_N2=m_n2, molecular_loss_kg_s=tidal_loss,
                               tidal_review=summary['exobase_radius_R'] > TIDAL_REVIEW_R)
            except (ValueError, RuntimeError, FloatingPointError):
                in_domain = False
        if not in_domain:
            row.update(status='outside_thermal_column_domain')
        rows.append(row)
    return rows


def loss_estimate(row):
    """The thermal column's loss where it holds; none bounded where the air runs away; otherwise the central
    energy-limited bound (eta 0.1, 3 R)."""
    if row['status'] == 'thermal_column':
        return row['molecular_loss_kg_s']
    if row['status'] == 'runaway':
        return math.inf
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
            'atmosphere/loss_response/model.py', 'atmosphere/loss_response/traced.py', 'atmosphere/thermal_column.py',
            'atmosphere/loss_response/infrared.py', 'atmosphere/radiative_convective/inputs.json',
            'atmosphere/middle_atmosphere/escape.py', 'atmosphere/loss_response/tides.py',
            'atmosphere/loss_response/results/tidal_escape.json',
            'atmosphere/loss_response/results/solar_cycle.json',
            'atmosphere/middle_atmosphere/results/limb_heat.json')}),
        evidence=('Screening model. The ultraviolet branch is the molecular thermal column (radiating CO2\'s '
                  '15-um band and the base\'s atomic oxygen and NO in excess of what air at the base temperature '
                  'emits; no atomic oxygen made above the base) above middle-atmosphere base temperatures, taking '
                  'the heat the limb tracing '
                  'finds along slant paths at the ring fleet\'s protected radius where the tracing puts it in height, '
                  'with its Jeans escape raised by '
                  'the tidal multiplier of test molecules in the Earth-Moon three-body problem, and energy-limited '
                  'upper bounds beyond its domain; the solar-wind branch is a range from mass-loading and '
                  'sputtering scalings with assumed parameter ranges. Neither is a coupled aeronomy or plasma model.'),
        reading_rule=('Loss rates are global-mean kg/s for the finished 1.2 atm atmosphere. Use the thermal-column '
                      'loss where status is thermal_column (molecular_loss_kg_s includes Earth\'s tide; '
                      'molecular_loss_two_body_kg_s leaves it out) and the energy-limited range beyond it; a runaway '
                      'row is one whose air outgrows the limb tables. A budget row\'s allowed_leak_fraction is the '
                      'largest UV transmission through gaps, the share of sunlight below 175 nm that passes over the '
                      'whole aperture of the 4 lunar-radius protected region, while ultraviolet-driven loss plus the '
                      'stated solar-wind loss stays within the budget. '
                      'It leaves out the sunlit exosphere\'s loss; exosphere_loss.json (schema '
                      'terluna.atmosphere.exosphere-loss/2) gives the allowed transmissions with it.'),
        assumptions=dict(base_pressure_pa=BASE_PA, lyman_glow_rayleigh_quiet=GLOW_RAYLEIGH,
                         activity=ACTIVITY, xray_cycle=escape.xray_cycle(), protected_radius_R=PROTECTED_R,
                         annulus_film=traced.FILM, solar_wind=SOLAR_WIND, solar_wind_ranges=SW_RANGES,
                         energy_limited_eta=ETAS, absorption_radii_R=ABSORPTION_RADII,
                         infrared_cooling=dict(radiates=INFRARED_COOLING, co2_ppm=infrared.CO2_PPM,
                                               co2_band_cm=infrared.BAND_CM),
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
