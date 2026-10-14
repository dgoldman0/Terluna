"""PlaSim's clear-sky sunlight, column by column: a copy of its arithmetic to recalibrate, offline, how it divides
the sunlight between the air and the ground against the line-by-line model.

    climate/gcm/.venv/bin/python -m climate.gcm.plasim_shortwave validate A28_dim5_sun_test --years 0 1
    climate/gcm/.venv/bin/python -m climate.gcm.plasim_shortwave fit A28_dim5 --years 15 24

PlaSim's shortwave (plasim/src/radmod.f90, swr) has two bands split at 0.75 um. Below it the clear air only
scatters, with all the Rayleigh scattering placed in the lowest layer; above it only water vapour absorbs, with
Lacis and Hansen's (1974) absorptivity of the water path scaled by pressure and temperature and magnified along
the slant; the layers combine by the adding method. `clear_sky` repeats that for one column (no ozone, no
aerosols, as the runs have them). `validate` compares it with the GCM's own clear-sky diagnostics on a run's
columns; `fit` finds the absorption setting that brings it to the line-by-line model on the columns of
radiation_check.py, for each way of raising the absorption, and writes ../results/gcm/shortwave_<run>.json.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'gcm'
SCHEMA = 'terluna.climate.gcm-shortwave/1'
BETA = 1.66                          # PlaSim's zbetta: the diffuse beam's magnification of the water path
RCOEFF_PER_RAYSCALE = 0.3856863591763577 / 1.8   # the Rayleigh coefficient the titania-filtered spectrum gives
FORMS = ('path', 'scale', 'linear')
# Each form raises the water-vapour absorptivity A(u) of the effective water path u (cm): 'path' takes A(k u),
# 'scale' k A(u) and 'linear' A(u) + k u, the last like a continuum growing with the path. 'path' keeps Lacis
# and Hansen's ceiling; the other two are capped just short of absorbing the whole band. The runner uses
# 'scale' (its SWH2OSCALE patch), which carried over best to a warmer, moister column the fit did not see.
# PlaSim's bands hold no CO2 or O2, whose strongest lines the Moon's deep column takes high up. GAS_BANDS adds them
# (the twilight copy's GAS_EDITS in exoplasim_run.py): CO2 and O2 in the near-infrared band (CO2's 2.0, 2.7 and 4.3 um,
# O2's A band and 1.27 um) and O2 in the visible (its B and gamma bands), each an absorptivity of its own
# pressure-scaled path, fitted by `gases` to the line-by-line model on the design column and combined with the
# water's (or nothing's) by multiplying transmissions. The collision-induced absorption, low in the column, stays in
# the water's calibration.
GAS_BANDS = (('co2_nir', 'CO2', 2), ('o2_nir', 'O2', 2), ('o2_vis', 'O2', 1))    # name, gas, PlaSim band
GAS_SLANTS = (1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 13.0)
GAS_P0 = 1e5
GASES = RESULTS / 'shortwave_gases.json'
DESIGN_FLUX_W_M2 = 1173.488286761312   # the design run's sunlight above the air (A28_dim5_moon's input record)


def gas_absorptivity(u, q):
    """A gas band's absorptivity of its effective path u (kg/m2 of the gas, each layer's scaled by (p / 1e5)^n), a
    share of the band's sunlight: q = (a, b, c, d, n)."""
    return np.minimum(q[0] * u / ((1.0 + q[1] * u) ** q[2] + q[3] * u), 0.999)


def path_transmissivities(path, slant, absorptivity):
    """PlaSim's layer transmissivities (top first) of an absorber whose layers hold path: down, from the direct
    beam's slant path to each layer's bottom, and up, for the light the ground returns (the whole slant path down,
    then BETA times the vertical path back up)."""
    vertical, slanted = np.cumsum(path), np.cumsum(slant * path)
    down, up = np.empty(path.size), np.empty(path.size)
    above, above_up = 1.0, 1.0 - absorptivity(slanted[-1] + BETA * vertical[-1])
    for l in range(path.size):
        down[l] = (1.0 - absorptivity(slanted[l])) / above
        above *= down[l]
        up[l] = above_up / (1.0 - absorptivity(slanted[-1] + BETA * (vertical[-1] - vertical[l])))
        above_up /= up[l]
    return down, up


def lacis_hansen(u):
    """Water vapour's absorptivity of the whole solar beam for an effective path u in cm (Lacis and Hansen 1974)."""
    return 2.9 * u / ((1.0 + 141.5 * u) ** 0.635 + 5.925 * u)


def band_absorptivity(u, solar2: float, form: str | None = None, k: float = 1.0):
    """PlaSim's absorptivity within its near-infrared band (A / zsolar2), with an optional recalibration."""
    if form in (None, 'path'):
        a = lacis_hansen(u * (k if form else 1.0))
    elif form == 'scale':
        a = k * lacis_hansen(u)
    elif form == 'linear':
        a = lacis_hansen(u) + k * u
    else:
        raise ValueError(f'form must be one of {FORMS}')
    return np.minimum(a / solar2, 0.999)


def half_levels(sigma) -> np.ndarray:
    """PlaSim's half levels from its full levels, which stand midway between them (sigma = 0 at the top)."""
    out, below = np.empty(len(sigma)), 0.0
    for l, s in enumerate(sigma):
        below = out[l] = 2.0 * s - below
    return out


def adding(r, rs, t, tu, albedo_direct: float, albedo_diffuse: float) -> np.ndarray:
    """Net downward flux at each interface, top to ground, per unit of sunlight entering the band: swr's adding
    method with layer reflectivities (direct r, diffuse rs) and transmissivities (down t, up tu)."""
    n = r.size
    down_t, down_rs = np.empty(n + 1), np.empty(n + 1)
    ta, tas, ras = 1.0, 1.0, 0.0
    for l in range(n):
        down_t[l], down_rs[l] = ta, ras
        z = 1.0 / (1.0 - ras * rs[l])
        ta, ras, tas = ta * t[l] * z, rs[l] + tu[l] * ras * t[l] * z, tu[l] * tas * z
    down_t[n], down_rs[n] = ta, ras
    up_r, up_rs = np.empty(n + 1), np.empty(n + 1)
    ra, ras = albedo_direct, albedo_diffuse
    for l in range(n - 1, -1, -1):
        up_r[l + 1], up_rs[l + 1] = ra, ras
        z = 1.0 / (1.0 - ras * rs[l])
        ra, ras = r[l] + t[l] * ra * tu[l] * z, rs[l] + tu[l] * ras * tu[l] * z
    up_r[0], up_rs[0] = ra, ras
    z = 1.0 / (1.0 - down_rs * up_rs)
    return down_t * z * (1.0 - up_r)


def clear_sky(column: dict, mu: float, form: str | None = None, k: float = 1.0) -> dict:
    """Shares of the sunlight arriving at the top of a column reflected, absorbed in the air and absorbed by the
    ground at cosine of zenith mu, as PlaSim's swr has them in clear sky. The column holds sigma (full levels, top
    first), ps (Pa), t (K), q (kg/kg), albedo (diffuse, and direct unless albedo_direct is given), rcoeff, solar1
    (the share of sunlight below 0.75 um) and gravity; with magnification (a, b, c), the direct beam's water path
    takes (mu + a) / (mu^2 + b mu + c), the Moon's (twilight.py), in place of Earth's 35 / sqrt(1 + 1224 mu^2); with
    gases (fits, fractions and the well-mixed air's magnification, as shortwave_gases.json and twilight.json hold
    them) the CO2 and O2 bands are added; with water_exponent n each layer's water path is weighted by (p / 1 bar)^n
    in place of PlaSim's p / 1 bar."""
    net = net_flux(column, mu, form, k)
    return dict(reflected=1.0 - net[0], air=net[0] - net[-1], ground=net[-1])


def net_flux(column: dict, mu: float, form: str | None = None, k: float = 1.0) -> np.ndarray:
    """clear_sky's net downward flux at each interface, top to ground, per unit of the sunlight arriving at the
    top: its differences are the layers' absorption, as PlaSim's swr heats them."""
    c = column
    sigma, n = np.asarray(c['sigma'], float), len(c['sigma'])
    dsigma = np.diff(np.concatenate([[0.0], half_levels(sigma)]))
    solar2 = 1.0 - c['solar1']
    lunar = c.get('magnification')
    m = (35.0 / np.sqrt(1.0 + 1224.0 * mu * mu) if lunar is None     # the direct beam's magnification
         else (mu + lunar[0]) / (mu * mu + lunar[1] * mu + lunar[2]))
    wv = 0.1 * dsigma * c['q'] * c['ps'] / c['gravity'] * np.sqrt(273.0 / c['t'])
    exponent = c.get('water_exponent', 1.0)
    wv = wv * sigma * c['ps'] / 1e5 if exponent == 1.0 else wv * (sigma * c['ps'] / 1e5) ** exponent
    t2, t2u = path_transmissivities(wv, m, lambda u: band_absorptivity(u, solar2, form, k))
    t1, t1u = np.ones(n), np.ones(n)
    gases = c.get('gases')
    if gases:
        air = gases.get('magnification')
        mg = m if air is None else (mu + air[0]) / (mu * mu + air[1] * mu + air[2])
        mass = dsigma * c['ps'] / c['gravity']
        for name, gas, band in GAS_BANDS:
            q = gases['fits'][name]
            path = mass * gases['fractions'][gas] * (sigma * c['ps'] / GAS_P0) ** q[4]
            down, up = path_transmissivities(path, mg, lambda u, q=q: gas_absorptivity(u, q))
            if band == 2:
                t2, t2u = t2 * down, t2u * up
            else:
                t1, t1u = t1 * down, t1u * up
    scale = c['rcoeff'] * c['ps'] / 101100.0 * 9.80665 / c['gravity']
    r1, r1s = np.zeros(n), np.zeros(n)
    r1[-1] = 1.0 - np.exp(scale * np.log(1.0 - 0.219 / (1.0 + 0.816 * mu)))
    r1s[-1] = 1.0 - np.exp(scale * np.log(1.0 - 0.144))
    direct = c.get('albedo_direct', c['albedo'])
    direct = direct(mu) if callable(direct) else direct
    visible = adding(r1, r1s, t1 - r1, t1u - r1s, direct, c['albedo'])
    infrared = adding(np.zeros(n), np.zeros(n), t2, t2u, direct, c['albedo'])
    return c['solar1'] * visible + solar2 * infrared


def ocean_direct_albedo(mu: float) -> float:
    """PlaSim's ice-free sea albedo for the direct beam (ECHAM-3)."""
    return min(0.05 / (mu + 0.15), 0.15)


def daily_split(column: dict, lat_deg: float, form: str | None = None, k: float = 1.0, nodes: int = 6) -> dict:
    """Shares of the day's clear-sky sunlight with the Sun over the equator, as radiation_check's line-by-line
    daily means take it: Gauss-Legendre nodes in hour angle over the sunlit half day, weighted by sunlight."""
    x, w = np.polynomial.legendre.leggauss(nodes)
    hours, weights = 0.25 * np.pi * (x + 1.0), 0.25 * np.pi * w
    mus = np.cos(np.radians(lat_deg)) * np.cos(hours)
    parts = [clear_sky(column, mu, form, k) for mu in mus]
    total = float(np.sum(weights * mus))
    return {key: float(np.sum(weights * mus * np.array([p[key] for p in parts])) / total) for key in ('reflected', 'air', 'ground')}


def column(f: dict, j: int, i: int, model_sunlight: dict, gravity: float, rayscale: float = 1.8,
           ocean: bool = False) -> dict:
    """A GCM column's time-mean profile as clear_sky takes it."""
    return dict(sigma=f['sigma'], ps=float(f['ps'][j, i]) * 100.0, t=f['ta'][:, j, i], q=f['hus'][:, j, i],
                albedo=float(f['alb'][j, i]) if not ocean else 0.069, rcoeff=RCOEFF_PER_RAYSCALE * rayscale,
                solar1=model_sunlight['energy_below_0.75um'], gravity=gravity,
                **(dict(albedo_direct=ocean_direct_albedo) if ocean else {}))


def _run_context(run: str, years):
    from climate.gcm.radiation_check import mean_fields
    from climate.crm.cm1_run import planet
    progress = json.loads((HERE / 'runs' / run / 'progress.json').read_text())
    return mean_fields(run, years), progress['model_sunlight'], planet()['gravity_m_s2'], progress['configuration']['model']


def validate(run: str, years) -> dict:
    """This copy of PlaSim's shortwave against the GCM's own clear-sky diagnostics, on the columns radiation_check
    picks. The GCM's are means of snapshots over the years; this takes the mean profile with the Sun over the
    equator, which suits runs with the Moon's 1.54 degree tilt."""
    from climate.gcm.radiation_check import choose_columns
    f, sunlight, gravity, model = _run_context(run, years)
    rows = []
    k = model.get('sw_water_absorption', 1.0)                       # the run's own water absorption
    for pick in choose_columns(f):
        j, i = pick['j'], pick['i']
        land = bool(f['lsm'][j, i])
        copy = daily_split(column(f, j, i, sunlight, gravity, model['rayleigh_scale'], ocean=not land), float(f['lat'][j]),
                           'scale' if k != 1.0 else None, k)
        arriving = float(f['arriving'][j, i])
        own = dict(reflected=1.0 - float(f['rstcs'][j, i]) / arriving,
                   air=float(f['rstcs'][j, i] - f['rsscs'][j, i]) / arriving, ground=float(f['rsscs'][j, i]) / arriving)
        rows.append(dict(pick, land=land, lat_deg=float(f['lat'][j]), gcm=own, copy=copy,
                         difference_points={k: 100.0 * (copy[k] - own[k]) for k in own}))
    return dict(run=run, years=list(years), columns=rows)


def gas_fit(processes: int = 2) -> dict:
    """CO2 and O2 in PlaSim's bands, fitted to the line-by-line model (lines alone, no continuum or collision-induced
    absorption) on the design column: for each band, the share of the band's sunlight the direct beam loses above each
    of the column's 100 levels at slant factors GAS_SLANTS, against (a, b, c, d) and the pressure exponent n of
    gas_absorptivity; with how the fit shares the absorption among the GCM's ten layers against the line-by-line."""
    from dataclasses import replace
    from scipy.optimize import least_squares
    from atmosphere.radiative_convective import optics, shortwave as sw
    from illumination.surface_light import model as surface
    from climate.gcm import twilight as tw
    col, air, _, _ = surface.build(surface.COLUMNS[tw.WORLD])
    shield, _ = surface.sunlight('design')
    mass, p = col['layer_mass_kg_m2'], col['layer_p_pa']
    fractions = {'CO2': air.fractions['CO2'] * 0.04401 / air.molar_mass, 'O2': air.fractions['O2'] * 0.031999 / air.molar_mass}
    energy = {1: tw.PLASIM_SOLAR1 * DESIGN_FLUX_W_M2, 2: (1.0 - tw.PLASIM_SOLAR1) * DESIGN_FLUX_W_M2}
    m = np.array(GAS_SLANTS)
    out = dict(fits={}, fractions=fractions, band_energy_w_m2=energy, checks={})
    for gas in ('CO2', 'O2'):
        cfg = replace(optics.SHORTWAVE, molecules=(gas,), continuum=False, cia=False, processes=processes)
        nu, tau = optics.optical_depth(col, cfg)
        e = sw.solar_spectrum(nu, shield=shield)[0] * np.gradient(nu)
        for name, g, band in GAS_BANDS:
            if g != gas:
                continue
            sel = nu < 1e7 / 750.0 if band == 2 else nu >= 1e7 / 750.0
            above = np.cumsum(tau[::-1, sel], 0)[::-1]                     # above each layer's bottom
            lbl = np.array([[e[sel] @ (1.0 - np.exp(-x * above[k])) for k in range(mass.size)] for x in m]) / energy[band]
            paths = lambda n: np.cumsum((mass * fractions[gas] * (p / GAS_P0) ** n)[::-1])[::-1]
            resid = lambda q: (gas_absorptivity(m[:, None] * paths(q[4])[None, :], q) - lbl).ravel() / lbl.max()
            best = min((least_squares(resid, [1.0, 1.0, 0.5, 0.0, n0], bounds=([0, 0, 0.05, 0, 0], [1e6, 1e9, 2.0, 1e6, 1.5]))
                        for n0 in (0.3, 0.6, 0.9)), key=lambda r: r.cost)
            q = [float(x) for x in best.x]
            fit = gas_absorptivity(m[:, None] * paths(q[4])[None, :], q)
            layers = lambda a: (a - np.r_[a[1:], 0.0])                      # each layer's share, bottom-up
            gcm = lambda a: tw.onto_gcm_layers(layers(a)[:, None], col['p_pa'])[:, 0] * energy[band]
            out['fits'][name] = q
            out['checks'][name] = dict(column_share_overhead=float(lbl[0, 0]), max_error=float(np.abs(fit - lbl).max()),
                                       overhead_w_m2_by_gcm_layer=dict(line_by_line=list(gcm(lbl[0])), fit=list(gcm(fit[0]))),
                                       slant_3_w_m2_by_gcm_layer=dict(line_by_line=list(gcm(lbl[3])), fit=list(gcm(fit[3]))))
    return out


def lunar_column(col: dict, gases: dict, product: dict) -> dict:
    """A column with the Moon's slant factors (twilight.json) and the CO2 and O2 bands (shortwave_gases.json)."""
    air = product['magnification_air']
    return dict(col, magnification=[product['magnification'][k] for k in 'abc'],
                gases=dict(fits=gases['fits'], fractions=gases['fractions'], magnification=[air[k] for k in 'abc']))


_COLUMNS = {}


def fit(check_path: Path, tests=(), lunar: bool = False, water_exponent: float = 1.0, forms=FORMS) -> dict:
    """For each form, the setting that brings PlaSim's clear-sky air absorption to line-by-line's on the columns of
    a radiation check (the same mean profiles, surface albedo and Sun), with the Rayleigh multiplier refitted
    to keep the reflection; the columns' remaining differences; and the differences on the columns of other
    checks (`tests`), which the fit did not see. With lunar, the columns take the Moon's slant factors and the CO2
    and O2 bands, and water_exponent weights their water paths (net_flux)."""
    from scipy.optimize import minimize

    if lunar:
        gases, product = json.loads(GASES.read_text()), json.loads((RESULTS / 'twilight.json').read_text())
        dress = lambda col: dict(lunar_column(col, gases, product), water_exponent=water_exponent)
    else:
        dress = lambda col: dict(col, water_exponent=water_exponent)

    def columns_of(check):
        key = (check['run'], tuple(check['years']), tuple((c['j'], c['i']) for c in check['columns']))
        if key not in _COLUMNS:                                   # the run's mean fields, read once
            f, sunlight, gravity, model = _run_context(check['run'], check['years'])
            _COLUMNS[key] = [(c, dict(column(f, c['j'], c['i'], sunlight, gravity), albedo=c['surface_albedo']))
                             for c in check['columns']], model
        cols, model = _COLUMNS[key]
        return [(c, dress(col)) for c, col in cols], model

    def split(cols, form, k, rayscale):
        return [daily_split(dict(col, rcoeff=RCOEFF_PER_RAYSCALE * rayscale), c['lat_deg'], form, k) for c, col in cols]

    def points(cols, got):
        return {key: [100.0 * (g[key] - c['line_by_line'][key]) for g, (c, _) in zip(got, cols)] for key in ('reflected', 'air', 'ground')}

    def summary(cols, got):
        d = points(cols, got)
        return dict(mean_points={k: float(np.mean(v)) for k, v in d.items()},
                    rms_points={k: float(np.sqrt(np.mean(np.square(v)))) for k, v in d.items()}, by_column=d)

    check = json.loads(check_path.read_text())
    cols, model = columns_of(check)
    others = {p.name: columns_of(json.loads(p.read_text()))[0] for p in tests}

    def cost(form, x):
        d = points(cols, split(cols, form, x[0], x[1]))
        return sum(a * a + r * r for a, r in zip(d['air'], d['reflected']))

    now = model['rayleigh_scale']
    out = dict(as_run=dict(rayleigh_scale=now, **summary(cols, split(cols, None, 1.0, now)),
                           tests={name: points(c, split(c, None, 1.0, now)) for name, c in others.items()}), forms={})
    starts = dict(path=5.0, scale=1.3, linear=0.002)
    for form in forms:
        best = minimize(lambda x: cost(form, x), [starts[form], now], method='Nelder-Mead',
                        options=dict(xatol=1e-4, fatol=1e-9, maxiter=400))
        out['forms'][form] = dict(k=float(best.x[0]), rayleigh_scale=float(best.x[1]), **summary(cols, split(cols, form, *best.x)),
                                  tests={name: points(c, split(c, form, *best.x)) for name, c in others.items()})
    out['columns'] = [dict(label=c['label'], lat_deg=c['lat_deg'], lon_deg=c['lon_deg'], land=c['land']) for c, _ in cols]
    out['test_columns'] = {name: [dict(run=json.loads((RESULTS / name).read_text())['run'], label=c['label'],
                                       water_kg_m2=c['water_kg_m2']) for c, _ in cs] for name, cs in others.items()}
    return dict(run=check['run'], years=check['years'], check=check_path.name, **out)


WATER_EXPONENTS = (1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3)


def design_column(product: dict, gases: dict) -> dict:
    """The twilight ledger's design column as PlaSim holds it, with the Moon's slant factors and the CO2 and O2
    bands (lunar_column)."""
    from illumination.surface_light import model as surface
    from climate.gcm import twilight as tw
    col, air, _, _ = surface.build(surface.COLUMNS[tw.WORLD])
    return lunar_column(tw.plasim_column(col, air.molar_mass, tw.SKY_ALBEDO), gases, product)


def design_layers(column: dict, product: dict, k: float, rayleigh: float, water_exponent: float = 1.0) -> np.ndarray:
    """PlaSim's global-mean absorbed sunlight in each of the GCM's layers (top first), W/m2 of the ledger's sunlight,
    on a design column at these water and Rayleigh settings, over the ledger's solar elevations."""
    import math
    from climate.gcm import twilight as tw
    c = dict(column, rcoeff=RCOEFF_PER_RAYSCALE * rayleigh, water_exponent=water_exponent)
    grid = np.radians(product['by_elevation']['elevation_deg'])
    layers = np.zeros((len(c['sigma']), grid.size))
    for i, el in enumerate(grid):
        mu = math.sin(el)
        if mu > 0.0:
            net = net_flux(c, mu, 'scale', k)
            layers[:, i] = mu * (net[:-1] - net[1:])
    return tw.sphere_share(layers.T, grid) * product['sunlight']['covered_w_m2'] / 4.0


def water_fit(check_path: Path, exponents=WATER_EXPONENTS) -> dict:
    """The water path's pressure exponent that best shares PlaSim's absorbed sunlight among the GCM's layers as the
    line-by-line column does: for each exponent, the water and Rayleigh multipliers (the 'scale' form) refitted on the
    radiation check's columns with the Moon's slant factors and the CO2 and O2 bands, then PlaSim's global-mean
    absorption by layer on the twilight ledger's design column against the line-by-line pseudo-spherical column's;
    the best exponent refined by a parabola through its neighbours and refitted."""
    product, gases = json.loads((RESULTS / 'twilight.json').read_text()), json.loads(GASES.read_text())
    column = design_column(product, gases)
    lbl = np.array(product['global_mean_w_m2']['whole_spectrum']['pseudo_spherical']['air_by_layer'])

    def trial(n):
        f = fit(check_path, lunar=True, water_exponent=n, forms=('scale',))['forms']['scale']
        layers = design_layers(column, product, f['k'], f['rayleigh_scale'], n)
        return dict(water_exponent=float(n), k=f['k'], rayleigh_scale=f['rayleigh_scale'], mean_points=f['mean_points'],
                    rms_points=f['rms_points'], air_by_layer_w_m2=list(layers),
                    difference_w_m2=list(layers - lbl), rms_layer_w_m2=float(np.sqrt(np.mean((layers - lbl) ** 2))))

    rows = [trial(n) for n in exponents]
    for r in rows:
        print(f"n {r['water_exponent']:.2f}: k {r['k']:.4f} rayleigh {r['rayleigh_scale']:.4f}  layers rms "
              f"{r['rms_layer_w_m2']:.3f} W/m2, air {r['mean_points']['air']:+.3f} points", flush=True)
    i = int(np.argmin([r['rms_layer_w_m2'] for r in rows]))
    best = rows[i]
    if 0 < i < len(rows) - 1:
        x = np.array([rows[j]['water_exponent'] for j in (i - 1, i, i + 1)])
        y = np.array([rows[j]['rms_layer_w_m2'] for j in (i - 1, i, i + 1)])
        a, b, _ = np.polyfit(x, y, 2)
        if a > 0:
            best = trial(round(float(-b / (2 * a)), 3))
            print(f"refined n {best['water_exponent']:.3f}: k {best['k']:.4f} rayleigh {best['rayleigh_scale']:.4f}  "
                  f"layers rms {best['rms_layer_w_m2']:.3f} W/m2", flush=True)
    return dict(check=check_path.name, line_by_line_w_m2=list(lbl), trials=rows, best=best)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('action', choices=('validate', 'fit', 'gases', 'water'))
    parser.add_argument('run', nargs='?')
    parser.add_argument('--years', type=int, nargs=2, default=(15, 24))
    parser.add_argument('--tests', nargs='*', default=[], help='other radiation checks in ../results/gcm to test the fit on')
    parser.add_argument('--lunar', action='store_true', help="fit with the Moon's slant factors and the CO2 and O2 bands")
    args = parser.parse_args(argv)
    if args.action == 'gases':
        result = gas_fit()
        for name, chk in result['checks'].items():
            print(f"{name}: {np.round(result['fits'][name], 4)}  overhead share {chk['column_share_overhead']:.4f}  "
                  f"max error {chk['max_error']:.4f}")
        result.update(schema='terluna.climate.gcm-shortwave-gases/1',
                      producer=dict(domain='climate', files={'gcm/plasim_shortwave.py': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]}))
        GASES.write_text(json.dumps(result, indent=1) + '\n')
        print(GASES)
        return 0
    if args.action == 'water':
        result = water_fit(RESULTS / f'radiation_check_{args.run}.json')
        result.update(schema='terluna.climate.gcm-shortwave-water/1',
                      producer=dict(domain='climate', files={'gcm/plasim_shortwave.py': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]},
                                    inputs={p: hashlib.sha256((RESULTS / p).read_bytes()).hexdigest()[:16]
                                            for p in ('twilight.json', 'shortwave_gases.json', f'radiation_check_{args.run}.json')}))
        from climate.crm.ring_crossings import rounded
        path = RESULTS / f'shortwave_water_{args.run}.json'
        path.write_text(json.dumps(rounded(result, 6), indent=1) + '\n')
        print(path)
        return 0
    if args.action == 'validate':
        result = validate(args.run, tuple(args.years))
        for r in result['columns']:
            d = r['difference_points']
            print(f"{r['label'][:40]:40s} {r['lat_deg']:6.1f}  copy - GCM: reflected {d['reflected']:+5.2f} air {d['air']:+5.2f} ground {d['ground']:+5.2f}")
    else:
        result = fit(RESULTS / f'radiation_check_{args.run}.json', [RESULTS / t for t in args.tests], args.lunar)
        a = result['as_run']
        print(f"as run: mean difference from line-by-line, points: {a['mean_points']}; tests {a['tests']}")
        for form, r in result['forms'].items():
            print(f"{form:7s} k={r['k']:.4g} rayleigh_scale={r['rayleigh_scale']:.3f}  mean {r['mean_points']}  rms {r['rms_points']}"
                  f"\n        tests {r['tests']}")
    result.update(schema=SCHEMA, action=args.action,
                  producer=dict(domain='climate', files={'gcm/plasim_shortwave.py': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]}))
    from climate.crm.ring_crossings import rounded
    RESULTS.mkdir(parents=True, exist_ok=True)
    path = RESULTS / f"shortwave_{args.action}_{args.run}{'_lunar' if args.lunar else ''}.json"
    path.write_text(json.dumps(rounded(result, 5), indent=1) + '\n')
    print(path)
    return 0


if __name__ == '__main__':
    sys.exit(main())
