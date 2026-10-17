#!/usr/bin/env python3
"""PlaSim's clear-sky longwave, column by column: a copy of its arithmetic to refit, offline, how it cools the air by
height against the atmosphere domain's correlated-k model, which follows the line-by-line spectroscopy.

    climate/gcm/.venv/bin/python -m climate.gcm.plasim_longwave validate A28_corrected_970_heating --year 0
    climate/gcm/.venv/bin/python -m climate.gcm.plasim_longwave fit A28_corrected_970_heating --year 0

PlaSim's longwave (plasim/src/radmod.f90, lwr) follows Sasamori (1968): broadband absorptivities of water vapour, its
continuum (1 - exp(-TH2OC u)), carbon dioxide with the water overlap of Boer et al. (1984), and ozone, each on the
absorber path between two levels scaled by pressure, with the layers' effective Planck emission between their half
levels. The Terluna build then adds the curved air's extra escape to space (twilight_lw). `lwr` repeats that for one
cloud-free column. `validate` sets it against PlaSim's own longwave heating (the gridpoint diagnostics) on a run's
cloud-free 3-day means; `fit` finds the continuum settings that bring its heating by layer, its outgoing longwave and
the longwave reaching the ground to the correlated-k model on the same columns, and writes
../results/gcm/longwave_<run>.json.

The water continuum PlaSim has scales with the pressure-scaled water path alone. The self-broadened continuum, which
dominates the window in moist air, scales with the vapour's own pressure as well, so `lwr` also takes TH2OE, a
continuum on the path weighted by the vapour pressure over a reference of 1 kPa, which puts its absorption in the
moist air near the ground. TH2OE = 0 is PlaSim as it is.
"""
from __future__ import annotations
import argparse
import ast
import hashlib
import json
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'gcm'
SCHEMA = 'terluna.climate.gcm-longwave/1'
SBK = 5.67e-8
ADV = 1850.0 / 1004.64 - 1.0           # the vapour's heat capacity over the dry air's, less one
E_REF_PA = 1000.0                      # the vapour pressure over which TH2OE's path is weighted
ZERO = 1e-6


def half_sigma(sigma):
    """PlaSim's half levels (the top first, from 0 to 1) for full levels midway between them."""
    half = [1.0]
    for s in np.asarray(sigma)[::-1]:
        half.append(2.0 * s - half[-1])
    return np.array(half[::-1])


def lwr(sigma, ps, t, q, land, gravity, acpd, co2_ppmv=400.0, th2oc=0.004, th2oe=0.0, twlw=None, twlt=None,
        irtop=0.0, tpofmt=1.0, pexp=1.0, texp=0.0, correction=None):
    """Net longwave flux (W/m2, positive down) at the tops of PlaSim's layers and the ground, and each layer's heating
    (K/s), for one cloud-free column: sigma, t (K) and q (specific humidity) on the full levels from the top, t with the
    ground's temperature appended; land 1 or 0 sets the ground's emissivity. pexp is the power of p/p0 that scales
    each layer's absorber paths (PlaSim's own is 1) and texp the power of 295 K over the layer's temperature (PlaSim's
    own is 0). correction, the fitted escape of each layer over land and over sea (CORRECTION), is added after the
    curved air's as the GCM adds it."""
    n = len(sigma)
    sig = np.asarray(sigma, float)
    halfs = half_sigma(sig)
    sigmah = halfs[1:]                                   # sigmah(1..NLEV): the half level below each full level
    dsig = halfs[1:] - halfs[:-1]
    t = np.asarray(t, float)
    q = np.asarray(q, float)
    # 1) absorber paths per layer: water in g/cm2 and CO2 in cm-STP, scaled by pressure (PlaSim's npbroaden = 1).
    zsfac = (sig * ps / 1.0e5) ** pexp * (295.0 / t[:n]) ** texp * dsig * ps / gravity
    zq = 0.1 * zsfac * q
    zqco2 = (100.0 / 1.9635) * (0.0440098 / 0.0289644) * 1.0e-6 * zsfac * co2_ppmv
    e = q * sig * ps / (0.622 + 0.378 * q)               # the vapour's partial pressure, Pa
    zqe = zq * e / E_REF_PA
    # 2) sigma T^4 on full levels (index 1..NLEP, NLEP the ground) and half levels.
    st4 = np.zeros(n + 2)
    st4[1:] = SBK * t ** 4
    st4h = np.zeros(n + 2)
    st4h[1] = st4[1] - (st4[1] - st4[2]) * sig[0] / (sig[0] - sig[1])
    for jl in range(2, n + 1):
        jm = jl - 1
        st4h[jl] = (st4[jl] * (sig[jm - 1] - sigmah[jm - 1]) + st4[jm] * (sigmah[jm - 1] - sig[jl - 1])) \
            / (sig[jm - 1] - sig[jl - 1])
    if (st4[n] - st4h[n]) * (st4[n] - st4[n + 1]) > 0.0:
        st4h[n + 1] = st4[n + 1]
    else:
        st4h[n + 1] = st4[n] + (st4[n] - st4[n - 1]) * (1.0 - sig[n - 1]) / (sig[n - 1] - sig[n - 2])
    zh2o0a = 0.846 * 3.59e-5 ** 0.243
    zh2o0 = 0.832 * 0.0286 ** 0.26
    zco20 = 0.0676 * 0.01022 ** 0.421
    zaco2c = 0.0676 * 1.01022 ** 0.421 - zco20
    zah2oc = 0.846 * (0.01 + 3.59e-5) ** 0.243 - zh2o0a - 0.24 * np.log10(0.02)
    zth2oc = 1.0 - (0.832 * (2.0 + 0.0286) ** 0.26 - zh2o0) + 0.1196 * np.log(2.0 - 0.6931)

    def transmissivity(wv, wve, co2):
        ah2o = 0.846 * (wv + 3.59e-5) ** 0.243 - zh2o0a if wv <= 0.01 else 0.24 * np.log10(wv + 0.01) + zah2oc
        ah2o = min(ah2o + (1.0 - np.exp(-th2oc * wv)) + (1.0 - np.exp(-th2oe * wve)), 1.0)
        aco2 = 0.0676 * (co2 + 0.01022) ** 0.421 - zco20 if co2 <= 1.0 else 0.0546 * np.log10(co2) + zaco2c
        th2o = 1.0 - (0.832 * (wv + 0.0286) ** 0.26 - zh2o0) if wv <= 2.0 else max(0.0, zth2oc - 0.1196 * np.log(wv - 0.6931))
        return min(1.0 - ZERO, max(ZERO, 1.0 - ah2o - aco2 * th2o))

    zbd = np.zeros(n + 2)
    zbu = np.zeros(n + 2)
    zbd[0] = irtop
    eps = 1.0 if land > 0.5 else 0.98
    zbu[n + 1] = eps * st4[n + 1]
    dftu = np.zeros(n + 2)
    dftd = np.zeros(n + 2)
    tausf = np.zeros(n + 2)
    for jl in range(1, n + 1):
        jm = jl - 1
        tau = np.zeros(n + 2)
        wv = wve = co2 = 0.0
        for j2 in range(jl, n + 1):
            wv += zq[j2 - 1]
            wve += zqe[j2 - 1]
            co2 += zqco2[j2 - 1]
            tau[j2] = transmissivity(wv, wve, co2)
        if jl == 1:
            tau0 = 1.0 - ZERO
            for j2 in range(1, n + 1):
                tau0 = min(1.0 - ZERO, max(ZERO, tau[j2] / tau0 * tpofmt))
                if (st4[j2] - st4h[j2]) * (st4[j2] - st4h[j2 + 1]) > 0.0:
                    zbd[j2] = 0.5 * st4[j2] + 0.25 * (st4h[j2] + st4h[j2 + 1])
                    zbu[j2] = zbd[j2]
                else:
                    zbd[j2] = (st4h[j2 + 1] - tau0 * st4h[j2]) / (1.0 - tau0) - (st4h[j2] - st4h[j2 + 1]) / np.log(tau0)
                    zbu[j2] = st4h[j2] + st4h[j2 + 1] - zbd[j2]
                tau0 = tau[j2]
        dftu[jl] -= zbu[jl]
        dftd[jl] += zbd[jm]
        zbdl = zbd[jm] - zbd[jl]
        for j2 in range(jl, n + 1):
            dftu[jl] -= (zbu[j2 + 1] - zbu[j2]) * tau[j2]
            dftd[j2 + 1] += zbdl * tau[j2]
        tausf[jl] = tau[n]
    dftu[n + 1] -= zbu[n + 1]
    dftd[n + 1] += zbd[n]
    grey = (1.0 - eps) * dftd[n + 1]
    for jl in range(1, n + 1):
        dftu[jl] -= tausf[jl] * grey
    net = (dftu + dftd)[1:]                              # tops of layers 1..NLEV, then the ground
    if twlw is not None:
        total = twlw[n] * (t[n] / twlt[n]) ** 4
        net[n] -= total
        for jl in range(n - 1, -1, -1):
            total += twlw[jl] * (t[jl] / twlt[jl]) ** 4
            net[jl] -= total
    olr = float(-dftu[1])
    if correction is not None:
        coeff, tref = np.asarray(correction['land' if land > 0.5 else 'sea']), np.asarray(correction['tref'])
        total = coeff[n] * (t[n] / tref[n]) ** 4
        net[n] -= total
        for jl in range(n - 1, -1, -1):
            total += coeff[jl] * (t[jl] / tref[jl]) ** 4
            net[jl] -= total
        olr += total
    heating = -gravity * (net[1:] - net[:-1]) / (dsig * ps * acpd * (1.0 + ADV * q))
    return net, heating, dict(olr=olr, ground_down=float(dftd[n + 1]))


def twilight_terms(run):
    """The curved air's extra escape (twlw, W/m2) and the temperatures it is given at (twlt), as the run's
    twilight.dat holds them for its layers and the ground."""
    lines = (HERE / 'runs' / run / 'model' / 'twilight.dat').read_text().split('\n')
    lines = [line for line in lines if line.strip()]
    return np.array(lines[-2].split(), float), np.array(lines[-1].split(), float)


def cloud_free_columns(run, year, every=1, max_cloud=0.01):
    """The run's cloud-free 3-day mean columns at every latitude: profiles on the full levels from the top, the
    ground's temperature, land, and PlaSim's own longwave heating (K/s)."""
    import netCDF4
    path = HERE / 'runs' / run / 'model' / f'MOST.{year:05d}.nc'
    with netCDF4.Dataset(path) as d:
        sigma = np.asarray(d['lev'][:], float)
        lat = np.asarray(d['lat'][:], float)
        cl = np.asarray(d['cl'][:], float).max(axis=1)             # (time, lat, lon)
        pick = np.argwhere(cl < max_cloud)[::every]
        ta, hus, lw = (np.asarray(d[k][:], float) for k in ('ta', 'hus', 'dt_lw'))
        ps = np.asarray(d['ps'][:], float) * 100.0
        ts = np.asarray(d['ts'][:], float)
        lsm = np.asarray(d['lsm'][0], float)
    columns = [dict(lat=float(lat[j]), land=float(lsm[j, i]), ps=float(ps[k, j, i]),
                    t=np.append(ta[k, :, j, i], ts[k, j, i]), q=hus[k, :, j, i], plasim_k_s=lw[k, :, j, i])
               for k, j, i in pick]
    return sigma, columns


def run_settings(run):
    from climate.crm.cm1_run import planet
    record = json.loads((HERE / 'runs' / run / 'progress.json').read_text())
    cfg = record['configuration']
    cfg = ast.literal_eval(cfg) if isinstance(cfg, str) else cfg
    gascon = 8.314462618 / 0.0289644 if 'gascon' not in cfg else cfg['gascon']
    planet_nml = (HERE / 'runs' / run / 'model' / 'planet_namelist').read_text()
    for line in planet_nml.splitlines():
        if line.strip().upper().startswith('GASCON'):
            gascon = float(line.split('=')[1])
    return dict(gravity=planet()['gravity_m_s2'], acpd=gascon / 0.286, co2=cfg.get('co2_ppm', 400.0),
                th2oc=cfg['model']['h2o_continuum'], pressure_pa=float(cfg['pressure_pa']),
                correction=cfg.get('longwave_correction'))


def validate(run, year, every=20):
    s = run_settings(run)
    twlw, twlt = twilight_terms(run)
    sigma, columns = cloud_free_columns(run, year, every)
    ours, theirs = [], []
    for c in columns:
        _, h, _ = lwr(sigma, c['ps'], c['t'], c['q'], c['land'], s['gravity'], s['acpd'], s['co2'], s['th2oc'],
                      twlw=twlw, twlt=twlt, correction=s['correction'])
        ours.append(h)
        theirs.append(c['plasim_k_s'])
    ours, theirs = np.array(ours) * 86400.0, np.array(theirs) * 86400.0
    return dict(columns=len(columns), sigma=sigma.tolist(), copy_k_day=ours.mean(axis=0).tolist(),
                plasim_k_day=theirs.mean(axis=0).tolist(),
                rms_k_day=np.sqrt(((ours - theirs) ** 2).mean(axis=0)).tolist())


def targets(run, year, every):
    """The correlated-k model on the run's cloud-free columns (plane-parallel, the design air): heating of the GCM's
    lowest nine layers (K/day, the ground first), outgoing longwave and the longwave reaching the ground."""
    from climate.crm.gcm_heating import gcm_column_levels
    from climate.crm.ring_radiation import mole_fraction
    from atmosphere.radiative_convective import ck
    from atmosphere.radiative_convective import thermodynamics as td
    s = run_settings(run)
    sigma, columns = cloud_free_columns(run, year, every)
    model = ck.CKLongwave()
    air = td.earthlike_air(s['pressure_pa'], s['co2'])
    out = []
    for c in columns:
        q = c['q']
        p, t, x = gcm_column_levels(sigma, c['t'][:-1], q / (1.0 - q), c['ps'])
        col = td.column_from_levels(td.MOON, air, p, t, mole_fraction(x, air.molar_mass), constant_gravity=True)
        rate, f = model.heating_rate(col, 1004.64, t_surface=float(c['t'][-1]))
        out.append(dict(c, ck_k_day=rate[:sigma.size - 1], ck_olr=f['olr'], ck_ground_down=f['surface_down']))
    return sigma, out


def misfit(sigma, columns, s, th2oc, th2oe, pexp=1.0, texp=0.0):
    """Mean errors of the copy against the correlated-k model, W/m2: by layer (the lowest nine, the ground first),
    in the outgoing longwave and in the longwave reaching the ground; and the score the fit lowers."""
    halfs = half_sigma(sigma)
    heat = []
    for c in columns:
        _, h, f = lwr(sigma, c['ps'], c['t'], c['q'], c['land'], s['gravity'], s['acpd'], s['co2'], th2oc, th2oe,
                      pexp=pexp, texp=texp)
        cap = s['acpd'] * c['ps'] * (halfs[1:] - halfs[:-1]) / s['gravity'] / 86400.0       # W/m2 per K/day, top first
        layer = (h[::-1][:sigma.size - 1] * 86400.0 - c['ck_k_day']) * cap[::-1][:sigma.size - 1]
        heat.append(np.concatenate([layer, [f['olr'] - c['ck_olr'], f['ground_down'] - c['ck_ground_down']]]))
    heat = np.array(heat)
    return heat.mean(axis=0), float((heat ** 2).mean(axis=0).sum())


def layer_errors(sigma, columns, s, correction=None):
    """Each column's error against correlated-k, W/m2: heating of the lowest nine layers (the ground first, negative
    cools too much), then the outgoing longwave; and the layers' temperatures, the top first."""
    halfs = half_sigma(sigma)
    n = sigma.size
    errs, temps = [], []
    for c in columns:
        _, h, f = lwr(sigma, c['ps'], c['t'], c['q'], c['land'], s['gravity'], s['acpd'], s['co2'], s['th2oc'],
                      correction=correction)
        cap = s['acpd'] * c['ps'] * (halfs[1:] - halfs[:-1]) / s['gravity'] / 86400.0
        layer = (h[::-1][:n - 1] * 86400.0 - c['ck_k_day']) * cap[::-1][:n - 1]
        errs.append(np.append(layer, f['olr'] - c['ck_olr']))
        temps.append(c['t'])
    return np.array(errs), np.array(temps)


def fit_correction(sigma, columns, s):
    """Each layer's extra escape over land and over sea, c (T / tref)^4 as the curved air's is given, that removes
    the mean errors in the lowest nine layers; the top layer's carries the outgoing longwave to correlated-k's."""
    n = sigma.size
    errs, temps = layer_errors(sigma, columns, s)
    tref = temps.mean(axis=0)
    out = dict(tref=tref.round(3).tolist())
    for surface in ('land', 'sea'):
        sel = np.array([(c['land'] > 0.5) == (surface == 'land') for c in columns])
        coeff = np.zeros(n + 1)
        for k in range(n - 1):                               # layer k from the ground is level n - 1 - k from the top
            level = n - 1 - k
            scale = (temps[sel, level] / tref[level]) ** 4
            coeff[level] = float((errs[sel, k] * scale).sum() / (scale ** 2).sum())
        scales = (temps[sel] / tref) ** 4
        added = (coeff * scales).sum(axis=1)                  # each column's change in the outgoing longwave so far
        coeff[0] = float(-(errs[sel, -1] + added).sum() / scales[:, 0].sum())
        out[surface] = coeff.round(4).tolist()
    return out


def fit(run, year, every):
    from scipy.optimize import minimize
    s = run_settings(run)
    sigma, columns = targets(run, year, every)
    base_bias, base_score = misfit(sigma, columns, s, s['th2oc'], 0.0)
    # A coarse grid, then a simplex from its best point.
    grid = [(misfit(sigma, columns, s, c0, e0)[1], c0, e0) for c0 in (0.0, 0.002, 0.004, 0.008, 0.016, 0.032)
            for e0 in (0.0, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0)]
    _, c0, e0 = min(grid)
    best = minimize(lambda x: misfit(sigma, columns, s, max(x[0], 0.0), max(x[1], 0.0))[1], [c0, max(e0, 1e-3)],
                    method='Nelder-Mead', options=dict(xatol=1e-5, fatol=1e-3, maxiter=120))
    th2oc, th2oe = (max(float(x), 0.0) for x in best.x)
    bias, score = misfit(sigma, columns, s, th2oc, th2oe)
    pbest = minimize(lambda x: misfit(sigma, columns, s, max(x[0], 0.0), 0.0, x[1])[1], [th2oc, 0.8],
                     method='Nelder-Mead', options=dict(xatol=1e-4, fatol=0.5, maxiter=120))
    pbias, pscore = misfit(sigma, columns, s, max(float(pbest.x[0]), 0.0), 0.0, float(pbest.x[1]))
    correction = fit_correction(sigma, columns, s)
    errs0, _ = layer_errors(sigma, columns, s)
    errs1, _ = layer_errors(sigma, columns, s, correction)
    names = [f'layer_{k + 1}' for k in range(sigma.size - 1)] + ['olr', 'ground_down']
    lnames = names[:-1]
    bands = {}
    for name, sel in (('tropics', lambda c: abs(c['lat']) < 30), ('mid_latitudes', lambda c: 30 <= abs(c['lat']) < 60),
                      ('high_latitudes', lambda c: abs(c['lat']) >= 60), ('land', lambda c: c['land'] > 0.5),
                      ('sea', lambda c: c['land'] <= 0.5)):
        pick = np.array([sel(c) for c in columns])
        if pick.any():
            bands[name] = dict(columns=int(pick.sum()),
                               plasim=dict(zip(lnames, errs0[pick].mean(axis=0).round(2).tolist())),
                               corrected=dict(zip(lnames, errs1[pick].mean(axis=0).round(2).tolist())),
                               corrected_rms=dict(zip(lnames, np.sqrt((errs1[pick] ** 2).mean(axis=0)).round(2).tolist())))
    return dict(schema=SCHEMA, run=run, year=year, columns=len(columns), e_ref_pa=E_REF_PA,
                plasim=dict(th2oc=s['th2oc'], score=base_score, mean_error_w_m2=dict(zip(names, base_bias.round(2).tolist()))),
                continua=dict(th2oc=th2oc, th2oe=th2oe, score=score, mean_error_w_m2=dict(zip(names, bias.round(2).tolist()))),
                pressure_exponent=dict(th2oc=max(float(pbest.x[0]), 0.0), pexp=float(pbest.x[1]), score=pscore,
                                       mean_error_w_m2=dict(zip(names, pbias.round(2).tolist()))),
                correction=correction, bands=bands,
                evidence=('PlaSim\'s clear-sky longwave as a copy of its arithmetic, which reproduces its own heating on '
                          'the run\'s cloud-free 3-day means within 0.01 K/day on average, against the atmosphere '
                          'domain\'s correlated-k model on the same columns; both plane-parallel, without the curved '
                          'air\'s escape, which the GCM adds on top.'),
                reading_rule=('Mean errors of PlaSim against correlated-k in W/m2: by layer from the ground (negative '
                              'cools too much), in the outgoing longwave and in the longwave reaching the ground. The '
                              'score is the sum over these of the mean squared error across columns. continua and '
                              'pressure_exponent are PlaSim\'s own settings refitted; correction is each layer\'s extra '
                              'escape to space over land and over sea, W/m2 at tref, scaled by (T / tref)^4, the top '
                              'layer first and the ground last, as the GCM reads it after the curved air\'s.'))


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('action', choices=('validate', 'fit'))
    parser.add_argument('run')
    parser.add_argument('--year', type=int, default=0)
    parser.add_argument('--every', type=int, default=20)
    args = parser.parse_args(argv)
    if args.action == 'fit':
        out = fit(args.run, args.year, args.every)
        out['producer'] = dict(domain='climate', files={f'gcm/{Path(__file__).name}':
                                                        hashlib.sha256(Path(__file__).read_bytes()).hexdigest()[:16]})
        RESULTS.mkdir(parents=True, exist_ok=True)
        (RESULTS / f'longwave_{args.run}.json').write_text(json.dumps(out, indent=1) + '\n')
        for k in ('plasim', 'continua', 'pressure_exponent'):
            print(k, round(out[k]['score'], 1), out[k]['mean_error_w_m2'])
        print('correction land', out['correction']['land']); print('correction sea', out['correction']['sea'])
        for name, b in out['bands'].items():
            print(f'{name:15s} corrected', list(b['corrected'].values()), 'rms', list(b['corrected_rms'].values()))
        return 0
    out = validate(args.run, args.year, args.every)
    print(f'{out["columns"]} cloud-free columns; by level from the top, the copy against PlaSim (K/day), rms:')
    for s, a, b, r in zip(out['sigma'], out['copy_k_day'], out['plasim_k_day'], out['rms_k_day']):
        print(f'  sigma {s:.3f}: {a:+.3f} against {b:+.3f}, rms {r:.3f}')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
