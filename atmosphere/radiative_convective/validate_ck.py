"""Check the correlated-k thermal scheme against line-by-line fluxes and heating rates.

    python -m atmosphere.radiative_convective.validate_ck      # prints and writes results/ck_validation.json

Three comparisons:

1. Whole-spectrum fluxes (OLR, surface downward longwave) for Earth and lunar
   columns against the 0.01 cm^-1 line-by-line model (H2O, CO2, continuum, CIA;
   no ozone, which that configuration lacks).
2. Heating rates in the CO2 15-um region (500-820 cm^-1) and the O3 9.6-um
   region (980-1080 cm^-1), from the surface to 0.1 Pa, against line-by-line
   spectra resolved to a quarter of the Doppler width (2e-4 cm^-1 in the CO2
   bands, 4e-4 in the ozone band), which is what stratospheric and mesospheric
   cooling needs and the 0.01 cm^-1 grid does not give.
3. A saturated 330 K lunar column, which reaches the table's moist nodes.
"""
from __future__ import annotations
import argparse
import json
import math
import multiprocessing as mp
import os
from pathlib import Path
import sys
import numpy as np

from atmosphere.radiative_convective import thermodynamics as th, optics, longwave as lw, spectroscopy as sp, ck

HERE = Path(__file__).resolve().parent


def column(planet, dry_pa, ts, humidity='mw', top=0.1, levels=(70, 30), strat_k=None):
    """Moist-adiabatic column with an isothermal (or warming, ozone-like) layer above the tropopause."""
    air = th.earthlike_air(dry_pa, 400.0)
    rh = th.manabe_wetherald() if humidity == 'mw' else th.uniform_rh(1.0)
    ps = dry_pa + float(rh(np.array([1.0]), 1.0)[0]) * float(th.saturation_pressure(ts))
    p = np.exp(np.concatenate([np.linspace(math.log(ps), math.log(0.01 * ps), levels[0] + 1),
                               np.linspace(math.log(0.01 * ps), math.log(top), levels[1] + 1)[1:]]))
    ad = th.moist_adiabat(ts, ps, air, t_stop=150.0, x_max=math.log(ps / top) + 0.1)
    t = th.adiabat_temperature(ad, ps, p)
    t_trop = 200.0 if strat_k is None else strat_k
    trop = int(np.argmax(t <= t_trop))
    # Above the tropopause: warming with height to 260 K near 100 Pa, then cooling to 190 K at the top,
    # a shape that exercises both strong and weak cooling regimes.
    lp = np.log(p)
    above = np.arange(p.size) > trop
    shape = t_trop + 60.0 * np.exp(-((lp - math.log(100.0)) / 2.0) ** 2)
    shape = np.where(p < 100.0, np.maximum(shape, t_trop + 60.0 - 7.0 * (math.log(100.0) - lp)), shape)
    t = np.where(above, np.clip(shape, 150.0, 300.0), t)
    x = np.clip(rh(p, ps) * th.saturation_pressure(np.maximum(t, 100.0)) / p, 0.0, 1.0)
    x = np.where(above, x[trop], x)
    col = th.column_from_levels(planet, air, p, t, x)
    # An ozone layer: 8 ppm peak at 1000 Pa (Earth) scaled to the O2 column (Moon), 0.1 ppm floor.
    peak = 1000.0 * (planet.surface_gravity / th.EARTH.surface_gravity)
    o3 = 8e-6 * np.exp(-((np.log(col['layer_p_pa']) - math.log(peak)) / 1.5) ** 2) + 1e-7
    col['o3'] = o3 * col['layer_column_cm2']
    return col


def lbl_whole(col):
    nu, tau = optics.optical_depth(col, optics.LONGWAVE)
    up, down = lw.fluxes(nu, tau, col['t_k'], col['t_k'][0])
    f_up, f_dn = lw.integrate(nu, up), lw.integrate(nu, down)
    return f_up, f_dn


_S = {}


def _layer_sigma(args):
    lo, hi, h, molecule, p, t, x = args
    if molecule not in _S:
        _S[molecule] = sp.load_lines(molecule, ck.S_MIN[molecule])
        gids = _S[molecule].gid
        _S[molecule + '_q'] = sp.PartitionSums(gids)
        if molecule == 'H2O':
            _S['cont'] = sp.WaterContinuum()
    grid = sp.Grid.span(lo, hi, h)
    alpha = ck.resolution(p, t)[1]
    coarse = max(1, int(round(0.01 / h)))
    sig = sp.line_absorption(grid, _S[molecule], _S[molecule + '_q'], t, p, x * p if molecule == 'H2O' else 0.0,
                             pedestal=(molecule == 'H2O'), coarse_factor=coarse,
                             fine_halfwidth=min(1.0, max(40 * alpha, 0.03)))
    if molecule == 'H2O':
        s_c, f_c = _S['cont'].cross_sections(grid.nu, t, p, x)
        sig = sig + s_c + f_c
    return sig.astype(np.float64)


def lbl_band_heating(col, lo, hi, h, molecules, processes):
    """High-resolution line-by-line net-flux profile in [lo, hi] for the listed molecules (plus CIA)."""
    grid = sp.Grid.span(lo, hi, h)
    nu = grid.nu
    nl = col['layer_p_pa'].size
    tau = np.zeros((nl, nu.size), dtype=np.float32)
    columns = {'H2O': th.species_column(col, 'H2O'), 'CO2': th.species_column(col, 'CO2'), 'O3': col['o3']}
    with mp.get_context('fork').Pool(processes) as pool:
        for m in molecules:
            jobs = [(lo, hi, h, m, float(col['layer_p_pa'][k]), float(col['layer_t_k'][k]),
                     float(col['layer_x_h2o'][k])) for k in range(nl)]
            for k, sig in enumerate(pool.imap(_layer_sigma, jobs, chunksize=1)):
                tau[k] += sig * columns[m][k]
    n_tot = th.layer_number_density(col)
    for a, b, f in sp.CIA_PAIRS:
        table = sp.CollisionInduced(f)
        for k in range(nl):
            x = col['layer_x_h2o'][k]
            fa = col['dry_fractions'].get(a, 0.0) * (1 - x) if a != 'H2O' else x
            fb = col['dry_fractions'].get(b, 0.0) * (1 - x) if b != 'H2O' else x
            if fa > 0 and fb > 0:
                tau[k] += table.coefficient(nu, float(col['layer_t_k'][k])) * fa * fb * n_tot[k] ** 2 * col['layer_dz_cm'][k]
    net = np.zeros(nl + 1)
    chunk = 100000
    for a in range(0, nu.size - 1, chunk):
        sl = slice(a, min(a + chunk + 1, nu.size))            # chunks share their edge point; each has two or more
        up, down = lw.fluxes(nu[sl], tau[:, sl].astype(float), col['t_k'], col['t_k'][0])
        net += lw.integrate(nu[sl], up) - lw.integrate(nu[sl], down)
    return net


def ck_band_net(model, col, bands):
    tau = model.optical_depth(col, extra={'O3': col['o3']})
    nl = tau.shape[0]
    b = ck.band_planck(col['t_k'])
    bs = ck.band_planck(np.array([col['t_k'][0]]))[0]
    mu, wmu = lw.angles(4)
    net = np.zeros(nl + 1)
    for bi in bands:
        up = np.zeros((nl + 1, ck.NG)); dn = np.zeros((nl + 1, ck.NG))
        for m, wt in zip(mu, wmu):
            weight = 2 * math.pi * wt * m
            i_d = np.zeros(ck.NG); dl = [i_d]; terms = []
            for k in range(nl - 1, -1, -1):
                e, a, c = lw._linear_source_terms(tau[k, bi] / m)
                terms.append((e, a, c))
                i_d = i_d * e + b[k, bi] * a - (b[k, bi] - b[k + 1, bi]) * c
                dl.append(i_d)
            terms.reverse(); dl.reverse()
            i_u = np.full(ck.NG, bs[bi]); up[0] += weight * i_u
            for k in range(nl):
                e, a, c = terms[k]
                i_u = i_u * e + b[k + 1, bi] * a - (b[k + 1, bi] - b[k, bi]) * c
                up[k + 1] += weight * i_u
            for k in range(nl + 1):
                dn[k] += weight * dl[k]
        net += ((up - dn) * ck.W).sum(axis=1)
    return net


def heating(col, net):
    cp = 1005.0
    return (net[:-1] - net[1:]) / (col['layer_mass_kg_m2'] * cp) * 86400.0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--processes', type=int, default=max(1, (os.cpu_count() or 2) // 2))
    parser.add_argument('--out', type=Path, default=HERE / 'results' / 'ck_validation.json')
    args = parser.parse_args(argv)
    model = ck.CKLongwave()
    cases = dict(earth_288=column(th.EARTH, 101325.0, 288.0), moon12_288=column(th.MOON, 121590.0, 288.0),
                 moon12_330_saturated=column(th.MOON, 121590.0, 330.0, humidity='saturated'))
    record = dict(schema='terluna.validation-record/1', model='atmosphere/radiative_convective/ck.py',
                  scope=__doc__.split('\n\n')[1].replace('\n', ' '), fluxes={}, heating={})
    for name, col in cases.items():
        f_up, f_dn = lbl_whole(col)
        noo3 = dict(col); noo3['o3'] = np.zeros_like(col['o3'])
        f = model.fluxes(noo3, extra={'O3': noo3['o3']})
        record['fluxes'][name] = dict(olr_lbl=float(f_up[-1]), olr_ck=f['olr'], dlr_lbl=float(f_dn[0]),
                                      dlr_ck=f['surface_down'],
                                      max_net_flux_difference=float(np.max(np.abs((f['up'] - f['down']) - (f_up - f_dn)))))
        print(name, record['fluxes'][name], flush=True)
    for name in ('earth_288', 'moon12_288'):
        col = cases[name]
        for label, lo, hi, h, mols, bands in (('co2_15um', 500.0, 820.0, 2e-4, ('H2O', 'CO2', 'O3'), (2, 3, 4)),
                                              ('o3_9.6um', 980.0, 1080.0, 4e-4, ('H2O', 'CO2', 'O3'), (6,))):
            ref = heating(col, lbl_band_heating(col, lo, hi, h, mols, args.processes))
            got = heating(col, ck_band_net(model, col, bands))
            p = col['layer_p_pa']
            sel = p < col['p_pa'][0] * 0.2
            rows = [dict(p_pa=float(p[k]), lbl_k_day=float(ref[k]), ck_k_day=float(got[k]))
                    for k in range(p.size) if sel[k]][::3]
            record['heating'][f'{name}/{label}'] = dict(
                band_cm1=[lo, hi], lbl_step_cm1=h, rows=rows,
                max_abs_difference_k_day=float(np.max(np.abs(ref - got)[sel])),
                rms_difference_k_day=float(np.sqrt(np.mean((ref - got)[sel] ** 2))))
            print(name, label, 'max |dQ|', round(record['heating'][f'{name}/{label}']['max_abs_difference_k_day'], 3),
                  'K/day', flush=True)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(record, indent=1) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
