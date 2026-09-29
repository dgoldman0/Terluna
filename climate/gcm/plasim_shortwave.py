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
    (the share of sunlight below 0.75 um) and gravity."""
    c = column
    sigma, n = np.asarray(c['sigma'], float), len(c['sigma'])
    dsigma = np.diff(np.concatenate([[0.0], half_levels(sigma)]))
    solar2 = 1.0 - c['solar1']
    m = 35.0 / np.sqrt(1.0 + 1224.0 * mu * mu)                      # the direct beam's magnification
    wv = 0.1 * dsigma * c['q'] * c['ps'] / c['gravity'] * np.sqrt(273.0 / c['t']) * sigma * c['ps'] / 1e5
    wvl, ywvl = np.cumsum(wv), np.cumsum(m * wv)
    a = lambda u: band_absorptivity(u, solar2, form, k)
    t2, t2u = np.empty(n), np.empty(n)
    above, above_up = 1.0, 1.0 - a(ywvl[-1] + BETA * wvl[-1])
    for l in range(n):
        t2[l] = (1.0 - a(ywvl[l])) / above
        above *= t2[l]
        t2u[l] = above_up / (1.0 - a(ywvl[-1] + BETA * (wvl[-1] - wvl[l])))
        above_up /= t2u[l]
    scale = c['rcoeff'] * c['ps'] / 101100.0 * 9.80665 / c['gravity']
    r1, r1s = np.zeros(n), np.zeros(n)
    r1[-1] = 1.0 - np.exp(scale * np.log(1.0 - 0.219 / (1.0 + 0.816 * mu)))
    r1s[-1] = 1.0 - np.exp(scale * np.log(1.0 - 0.144))
    direct = c.get('albedo_direct', c['albedo'])
    direct = direct(mu) if callable(direct) else direct
    visible = adding(r1, r1s, 1.0 - r1, 1.0 - r1s, direct, c['albedo'])
    infrared = adding(np.zeros(n), np.zeros(n), t2, t2u, direct, c['albedo'])
    top = c['solar1'] * visible[0] + solar2 * infrared[0]
    ground = c['solar1'] * visible[-1] + solar2 * infrared[-1]
    return dict(reflected=1.0 - top, air=top - ground, ground=ground)


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


def fit(check_path: Path, tests=()) -> dict:
    """For each form, the setting that brings PlaSim's clear-sky air absorption to line-by-line's on the columns of
    a radiation check (the same mean profiles, surface albedo and Sun), with the Rayleigh multiplier refitted
    to keep the reflection; the columns' remaining differences; and the differences on the columns of other
    checks (`tests`), which the fit did not see."""
    from scipy.optimize import minimize

    def columns_of(check):
        f, sunlight, gravity, model = _run_context(check['run'], check['years'])
        return [(c, dict(column(f, c['j'], c['i'], sunlight, gravity), albedo=c['surface_albedo'])) for c in check['columns']], model

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
    for form in FORMS:
        best = minimize(lambda x: cost(form, x), [starts[form], now], method='Nelder-Mead',
                        options=dict(xatol=1e-4, fatol=1e-9, maxiter=400))
        out['forms'][form] = dict(k=float(best.x[0]), rayleigh_scale=float(best.x[1]), **summary(cols, split(cols, form, *best.x)),
                                  tests={name: points(c, split(c, form, *best.x)) for name, c in others.items()})
    out['columns'] = [dict(label=c['label'], lat_deg=c['lat_deg'], lon_deg=c['lon_deg'], land=c['land']) for c, _ in cols]
    out['test_columns'] = {name: [dict(run=json.loads((RESULTS / name).read_text())['run'], label=c['label'],
                                       water_kg_m2=c['water_kg_m2']) for c, _ in cs] for name, cs in others.items()}
    return dict(run=check['run'], years=check['years'], check=check_path.name, **out)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('action', choices=('validate', 'fit'))
    parser.add_argument('run')
    parser.add_argument('--years', type=int, nargs=2, default=(15, 24))
    parser.add_argument('--tests', nargs='*', default=[], help='other radiation checks in ../results/gcm to test the fit on')
    args = parser.parse_args(argv)
    if args.action == 'validate':
        result = validate(args.run, tuple(args.years))
        for r in result['columns']:
            d = r['difference_points']
            print(f"{r['label'][:40]:40s} {r['lat_deg']:6.1f}  copy - GCM: reflected {d['reflected']:+5.2f} air {d['air']:+5.2f} ground {d['ground']:+5.2f}")
    else:
        result = fit(RESULTS / f'radiation_check_{args.run}.json', [RESULTS / t for t in args.tests])
        a = result['as_run']
        print(f"as run: mean difference from line-by-line, points: {a['mean_points']}; tests {a['tests']}")
        for form, r in result['forms'].items():
            print(f"{form:7s} k={r['k']:.4g} rayleigh_scale={r['rayleigh_scale']:.3f}  mean {r['mean_points']}  rms {r['rms_points']}"
                  f"\n        tests {r['tests']}")
    result.update(schema=SCHEMA, action=args.action,
                  producer=dict(domain='climate', files={'gcm/plasim_shortwave.py': hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]}))
    from climate.crm.ring_crossings import rounded
    RESULTS.mkdir(parents=True, exist_ok=True)
    path = RESULTS / f'shortwave_{args.action}_{args.run}.json'
    path.write_text(json.dumps(rounded(result, 5), indent=1) + '\n')
    print(path)
    return 0


if __name__ == '__main__':
    sys.exit(main())
