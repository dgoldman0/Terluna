"""PlaSim's boundary-layer mixing at lunar gravity, measured against the CM1 boxes.

PlaSim mixes the air near the ground in two steps each time step. Its vertical diffusion (fluxmod.f90, vdiff) gives
heat, vapour and momentum a diffusivity K = l^2 |dV/dz| F(Ri): Blackadar's mixing length l = lambda kappa z /
(lambda + kappa z), with lambda 160 m for momentum and 160 sqrt(7.5) = 438 m for heat and vapour and the heights z
counted at 250 K with the planet's gravity; ECHAM's stability functions with b = c = d = 5; a shear from the difference
in wind speed between levels, at least 1 m/s; and the coefficients carried to sigma with the planet's gravity. Its dry
adjustment (rainmod.f90, mkdca) then mixes potential temperature and vapour evenly through any layers where the dry
potential temperature falls with height. The CM1 boxes mix their boundary layer by resolved eddies and by CM1's simple
boundary-layer scheme (ipbl 2), whose vertical diffusivity is l^2 sqrt(max(S^2 - N^2/Pr, small)) with
l = (l_inf^-2 + (kappa z)^-2)^-1/2 and l_inf = 75 m scaled by the gravities to 453 m (crm/cm1_run.py).

By quarter of the lunar day at the two CM1 sites this module measures, in the GCM (a run's 3-day means in the site's
cells), the heat its vertical diffusion and its dry adjustment carry up through each half level, from the run's own
heating by process, and the diffusivity, Richardson number and mixing length its vertical diffusion computes on the
means, with the mixing lengths as they stand and stretched by the gravities. The run writes no vapour tendency by
process, and the diffusion recomputed on the means misstates the heat the run's own diffusion carries, so the
diffusivities on the means serve to place each half level's stability regime. In CM1 (each box's domain diagnostics
over its settled lunar day) it measures the heat and vapour carried up through each height, resolved and
parameterized, the effective diffusivity for vapour, the parameterized diffusivity, the vertical wind's spread and the
boundary layer's depth from CM1's bulk Richardson number. In both it reads the mixed layer's depth by one rule
(rise_height) and the mean profiles of virtual potential temperature and vapour.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent / 'results' / 'gcm'
CRM_RESULTS = HERE.parent / 'results' / 'crm'
SCHEMA = 'terluna.climate.gcm-pbl-mixing/1'
RUN = 'A28R_301k_cold_sbm_rh07_deep200_l15cap'   # the 15-level baseline the next convection tests compare with
YEARS = (0, 1)
KARMAN = 0.4                 # fluxmod's von Karman constant
TSCAL = 250.0                # K, the temperature vdiff counts its heights at
LAMBDA_M = 160.0             # m, vdiff_lamm, PlaSim's asymptotic mixing length for momentum
B = C = D = 5.0              # vdiff_b, vdiff_c, vdiff_d
UMIN = 1.0                   # m/s, zumin, the least wind-speed difference between levels
EARTH_G = 9.80665
LV = 2.501e6                 # J/kg, CM1's latent heat (constants.F), also used for PlaSim's vapour fluxes here
QUARTERS = (('sunrise-noon', -90.0, 0.0), ('noon-sunset', 0.0, 90.0), ('sunset-midnight', 90.0, 180.0),
            ('midnight-sunrise', -180.0, -90.0))
SITES = (('sea', 'box_near_side_sea', 'box_near_side_sea_day3'), ('land', 'box_0e_design', 'box_0e_design'))
TOP_KM = 20.0                # the profiles reported reach this high
RISE_K = 0.5                 # K, the virtual potential temperature's rise that ends a mixed layer (rise_height)
GRADIENT_FLOOR = 5.0e-8      # kg/kg per m (0.05 g/kg per km): CM1's effective diffusivity only where vapour falls faster


def quarter_name(hour_deg: float) -> str:
    h = (float(hour_deg) + 180.0) % 360.0 - 180.0
    return next(name for name, a, b in QUARTERS if a <= h < b)


def snapshot_hours(hour_means) -> np.ndarray:
    """Hour angles (degrees, -180..180) of PlaSim's gridpoint diagnostics, which are single time steps at the end of
    each mean (climate/crm/README.md, the snapshots and the exact means): the means' hour angles moved on by half the
    median step between consecutive means."""
    h = np.unwrap(np.radians(np.asarray(hour_means, float)))
    step = float(np.median(np.diff(h)))
    return (np.degrees(h + 0.5 * step) + 180.0) % 360.0 - 180.0


def vdiff_coefficients(t, q, speed, sigma, edges, g, gascon, rdbrv, akap, lamm=LAMBDA_M):
    """PlaSim's vertical diffusion coefficients for heat and vapour between each pair of full levels (top first), for
    columns t, q, speed (..., nlev): the half level's height as vdiff counts it (m), the mixing length (m), the shear
    |dV/dz| (1/s), the Richardson number, the stability factor F, the diffusivity K (m2/s) and its sigma form K_sigma
    (1/s, K (g sigma_h / R T_h)^2), as vdiff computes them (fluxmod.f90), with lamm the momentum length."""
    t, q, speed = (np.asarray(a, float) for a in (t, q, speed))
    sigma, edges = np.asarray(sigma, float), np.asarray(edges, float)
    dsig = np.diff(edges)
    lamh = lamm * math.sqrt(3.0 * D * 0.5)
    out = {k: [] for k in ('z_vdiff_m', 'l_m', 'shear_s', 'ri', 'f', 'k_m2_s', 'k_sigma_s')}
    for k in range(sigma.size - 1):
        sh = edges[k + 1]
        zz = -gascon * TSCAL * math.log(sh) / g
        dz = gascon * TSCAL * math.log(sigma[k + 1] / sigma[k]) / g
        mixh = lamh * KARMAN * zz / (lamh + KARMAN * zz)
        rdsig = 1.0 / (sigma[k + 1] - sigma[k])
        zkfac = math.sqrt(((((zz + dz) / zz) ** (1.0 / 3.0) - 1.0) / dz) ** 3 / zz)
        rthl = (dsig[k + 1] + dsig[k]) / (t[..., k] * dsig[k + 1] + t[..., k + 1] * dsig[k])
        dvds = np.maximum(UMIN, np.abs(speed[..., k] - speed[..., k + 1])) * rdsig
        qf1 = 1.0 + (1.0 / rdbrv - 1.0) * q[..., k]
        qf2 = 1.0 + (1.0 / rdbrv - 1.0) * q[..., k + 1]
        dthds = (t[..., k] * qf1 / sigma[k] ** akap - t[..., k + 1] * qf2 / sigma[k + 1] ** akap) * rdsig
        ri = sh ** akap * gascon * dthds / (sh * dvds * dvds)
        neg = np.sqrt(np.maximum(-ri, 0.0))
        denom_u = 1.0 + 3.0 * C * B * neg * zkfac * mixh * mixh
        f_unstable = 1.0 - 3.0 * B * ri / denom_u
        f_stable = 1.0 / (1.0 + 3.0 * B * np.maximum(ri, 0.0) * np.sqrt(1.0 + D * np.maximum(ri, 0.0)))
        f = np.where(ri <= 0.0, f_unstable, f_stable)
        to_z = g * sh * rthl / gascon                     # d sigma / dz at the half level (1/m), sign aside
        shear = dvds * to_z
        kk = mixh * mixh * f * shear
        out['z_vdiff_m'].append(np.full(np.shape(ri), zz))
        out['l_m'].append(np.full(np.shape(ri), mixh))
        out['shear_s'].append(shear)
        out['ri'].append(ri)
        out['f'].append(f)
        out['k_m2_s'].append(kk)
        out['k_sigma_s'].append(kk * to_z * to_z)
    return {k: np.stack(v, axis=-1) for k, v in out.items()}


def diffusion_fluxes(t, q, ps, sigma, edges, k_sigma, g, acpd, akap):
    """The heat and vapour (W/m2, upward) a sigma-form diffusivity carries through each half level between full levels
    on profiles t, q (top first), with the fluxes of potential temperature taken as vdiff takes them."""
    sigma, edges = np.asarray(sigma, float), np.asarray(edges, float)
    t, q = np.asarray(t, float), np.asarray(q, float)
    theta = t / sigma ** akap
    rdsig = 1.0 / (sigma[1:] - sigma[:-1])
    mass = np.asarray(ps, float)[..., None] / g
    heat = acpd * mass * k_sigma * edges[1:-1] ** akap * (theta[..., 1:] - theta[..., :-1]) * rdsig
    vapour = LV * mass * k_sigma * (q[..., 1:] - q[..., :-1]) * rdsig
    return heat, vapour


def half_level_heights(t, q, ps, edges, g, gascon, rdbrv):
    """Heights (m) above the ground of the half levels between full levels (top first), from virtual temperature."""
    tv = np.asarray(t, float) * (1.0 + (1.0 / rdbrv - 1.0) * np.asarray(q, float))
    edges = np.asarray(edges, float)
    n = tv.shape[-1]
    z = np.zeros(tv.shape[:-1] + (n + 1,))
    for k in range(n - 1, 0, -1):                          # the ground up; edges[0] = 0 is the top and is left out
        z[..., k] = z[..., k + 1] + gascon * tv[..., k] / g * math.log(edges[k + 1] / edges[k])
    return z[..., 1:n]


def full_level_heights(t, q, ps, sigma, edges, g, gascon, rdbrv):
    tv = np.asarray(t, float) * (1.0 + (1.0 / rdbrv - 1.0) * np.asarray(q, float))
    zh = half_level_heights(t, q, ps, edges, g, gascon, rdbrv)
    below = np.concatenate([zh, np.zeros(tv.shape[:-1] + (1,))], axis=-1)
    return below + gascon * tv / g * np.log(np.asarray(edges[1:], float) / np.asarray(sigma, float))


def gcm_site(run: str, years, box: str) -> dict:
    """The run's 3-day means in the cells of a CM1 box's site, by quarter of the lunar day: the heat its vertical
    diffusion and dry adjustment carry up through each half level (from the run's heating by process, single time steps
    at the end of each mean, binned by the time of those steps with the surface layer's heating at the same steps), the
    surface's
    sensible and latent heat, the vertical diffusion's coefficients on the means (mixing lengths as they stand and
    stretched by the gravities) with its heat flux recomputed there as a check, and the mean profiles."""
    import netCDF4
    from climate.crm import gcm_site as gs
    from climate.crm import ring_analysis as ra
    from climate.gcm import sbm_offline as so
    record = json.loads((ra.RUNS / box / 'case.json').read_text())
    cells = [tuple(cc) for cc in record['gcm_site']['cells']]
    lon_site = record['configuration']['site']['lon_deg']
    c = so.run_constants(run)
    g, gascon, rdbrv, acpd = c['ga'], c['gascon'], c['rdbrv'], c['acpd']
    akap = gascon / acpd
    adv = c['ADV']
    stretch = EARTH_G / g
    acc, acc_snap = {}, {}
    for year in range(years[0], years[1] + 1):
        with netCDF4.Dataset(HERE / 'runs' / run / 'model' / f'MOST.{year:05d}.nc') as d:
            lat, lon = np.asarray(d['lat'][:], float), np.asarray(d['lon'][:], float)
            sigma = np.asarray(d['lev'][:], float)
            czen = np.asarray(d['czen'][:], float)
            idx = [(int(np.argmin(abs(lat - a))), int(np.argmin(abs(lon - b)))) for a, b in cells]
            f = {k: np.stack([np.asarray(d[k][:, ..., j, i], float) for j, i in idx], axis=1)
                 for k in ('ta', 'hus', 'ua', 'va', 'ps', 'hfss', 'hfls', 'dt_vdif', 'dt_dca', 'dt_sfc')}
        edges = gs.plasim_edges(sigma)
        hour = (lon_site - gs.subsolar_longitude(czen, lat, lon) + 180.0) % 360.0 - 180.0
        hour_snap = snapshot_hours(hour)
        dsig = np.diff(edges)
        w = np.array([math.cos(math.radians(lat[j])) for j, _ in idx])
        for n in range(f['ta'].shape[0]):
            t, q = f['ta'][n], f['hus'][n]                                  # (cells, lev), top first
            ps = f['ps'][n] * 100.0
            speed = np.hypot(f['ua'][n], f['va'][n])
            mass = ps[:, None] * dsig[None, :] / g
            cpq = acpd * (1.0 + adv * q)
            heat_vdif = np.cumsum(cpq * mass * f['dt_vdif'][n], axis=1)[:, :-1]
            heat_dca = np.cumsum(cpq * mass * f['dt_dca'][n], axis=1)[:, :-1]
            sfc_layer = (cpq * mass * f['dt_sfc'][n]).sum(axis=1)
            co_now = vdiff_coefficients(t, q, speed, sigma, edges, g, gascon, rdbrv, akap)
            co_long = vdiff_coefficients(t, q, speed, sigma, edges, g, gascon, rdbrv, akap, LAMBDA_M * stretch)
            h_off, _ = diffusion_fluxes(t, q, ps, sigma, edges, co_now['k_sigma_s'], g, acpd, akap)
            zh = half_level_heights(t, q, ps, edges, g, gascon, rdbrv)
            zf = full_level_heights(t, q, ps, sigma, edges, g, gascon, rdbrv)
            theta = t * (1.0e5 / (sigma[None, :] * ps[:, None])) ** akap
            thv = theta * (1.0 + (1.0 / rdbrv - 1.0) * q)
            snap = acc_snap.setdefault(quarter_name(hour_snap[n]), {})
            for key, v in (('w', w), ('heat_vdiff_w_m2', heat_vdif), ('heat_dca_w_m2', heat_dca),
                           ('surface_layer_w_m2', sfc_layer)):
                snap.setdefault(key, []).append(np.asarray(v, float))
            a = acc.setdefault(quarter_name(hour[n]), {})
            rows = dict(w=w, z_half_m=zh, z_full_m=zf,
                        sensible_w_m2=-f['hfss'][n], latent_w_m2=-f['hfls'][n],
                        k_m2_s=co_now['k_m2_s'], k_stretched_m2_s=co_long['k_m2_s'], ri=co_now['ri'],
                        f=co_now['f'], shear_s=co_now['shear_s'], l_m=co_now['l_m'], l_stretched_m=co_long['l_m'],
                        heat_vdiff_on_means_w_m2=h_off,
                        theta_k=theta, theta_v_k=thv, q_g_kg=q / (1.0 - q) * 1000.0)   # mixing ratio, as CM1's
            for key, v in rows.items():
                a.setdefault(key, []).append(np.asarray(v, float))
    def weighted(a):
        w = np.concatenate(a.pop('w'))
        return w, {k: np.sum(np.concatenate(v) * (w[:, None] if np.concatenate(v).ndim > 1 else w), axis=0) / w.sum()
                   for k, v in a.items()}

    out = {}
    for name, _, _ in QUARTERS:
        a = acc[name]
        w, mean = weighted(a)
        w_snap, mean_snap = weighted(acc_snap[name])
        mean.update(mean_snap)
        mean['snapshots'] = int(w_snap.size)
        ri = np.concatenate(a['ri'])
        mean['ri_share_unstable'] = (ri < 0.0).mean(axis=0)
        mean['ri_median'] = np.median(ri, axis=0)
        mean['k_median'] = np.median(np.concatenate(a['k_m2_s']), axis=0)
        mean['outputs'] = int(w.size)
        out[name] = mean
    return dict(cells=cells, sigma=sigma.tolist(), stretch=stretch, constants=dict(g=g, gascon=gascon, acpd=acpd),
                quarters=out)


def cm1_box(case_name: str) -> dict:
    """A CM1 box's domain diagnostics over the span its onset record analyses, by quarter of the lunar day: the heat
    and vapour carried up through each w level (W/m2), resolved and parameterized; the vapour's kinematic flux and
    mean gradient, for its effective diffusivity; the parameterized diffusivity; the vertical wind's spread; the
    surface's fluxes; the boundary layer's depth (hpbl); and the mean profiles of potential temperature and vapour."""
    from climate.crm import box_budget as bb
    from climate.crm import cm1_output as co
    from climate.crm import ring_analysis as ra
    from climate.crm.cm1_run import solar_day_s
    from climate.crm.ring_crossings import local_hour_angle
    case = ra.RUNS / case_name
    record = json.loads((case / 'case.json').read_text())
    span = json.loads((CRM_RESULTS / f'box_onset_{case_name}.json').read_text())['span_days']
    start, lon = record['configuration']['start_hour_angle_deg'], record['configuration']['site']['lon_deg']
    geo = ra.case_geometry(case)
    zh, zw = geo['zh'], geo['zw']
    lay_s = ra.grads_layout(case / 'cm1out_diag_s.ctl')[1]
    lay_w = ra.grads_layout(case / 'cm1out_diag_w.ctl')[1]
    kappa = co.RD / co.CP
    acc = {}
    for path in sorted(case.glob('cm1out_diag_*_s.dat')):
        s = bb.read_diag(case, path, lay_s)
        if not span[0] * 86400.0 - 1.0 <= s['mtime'] <= span[1] * 86400.0 + 1.0 or s['mtime'] <= 0.0:
            continue
        wd = bb.read_diag(case, path.with_name(path.name.replace('_s.dat', '_w.dat')), lay_w)
        prs = wd['prsinterp'].copy()
        prs[0] = s['psfc']
        rho, exn = wd['rhointerp'], (np.maximum(prs, 1.0) / co.P00) ** kappa
        to_heat, to_vap = rho * co.CP * exn, rho * LV
        hour = local_hour_angle(start, lon, s['mtime'], solar_day_s())
        thv = s['th'] * (1.0 + (co.RV / co.RD - 1.0) * s['qv'])
        rows = dict(heat_resolved_w_m2=to_heat * wd['thfr'], heat_parameterized_w_m2=to_heat * (wd['thfs'] + wd['thfd']),
                    vapour_resolved_w_m2=to_vap * wd['qvfr'],
                    vapour_parameterized_w_m2=to_vap * (wd['qvfs'] + wd['qvfd']),
                    khv_m2_s=wd['khv'], dqvdz_per_m=wd['dqvdz'], vapour_kinematic=wd['qvfr'] + wd['qvfs'] + wd['qvfd'],
                    sigma_w_m_s=np.sqrt(np.maximum(wd['wpwp_w'], 0.0)),
                    sensible_w_m2=s['hfx'], latent_w_m2=s['qfx'] * LV,
                    hpbl_m=s['hpbl'],
                    theta_k=s['th'], theta_v_k=thv, q_g_kg=s['qv'] * 1000.0)
        a = acc.setdefault(quarter_name(hour), {})
        for key, v in rows.items():
            a.setdefault(key, []).append(np.asarray(v, float))
    out = {}
    for name, _, _ in QUARTERS:
        a = acc[name]
        mean = {k: np.mean(np.stack(v), axis=0) for k, v in a.items()}
        mean['heat_w_m2'] = mean['heat_resolved_w_m2'] + mean['heat_parameterized_w_m2']
        mean['vapour_w_m2'] = mean['vapour_resolved_w_m2'] + mean['vapour_parameterized_w_m2']
        mean['outputs'] = len(a['hpbl_m'])
        out[name] = mean
    l_inf = bb.namelist_value((case / 'namelist.input').read_text(), 'l_inf')
    return dict(case=case_name, span_days=span, zw_m=zw, zh_m=zh, l_inf_m=l_inf, quarters=out)


def rise_height(z, thv, z0, rise=RISE_K):
    """The first height above z0 (m) where the virtual potential temperature has risen RISE_K above its value at z0,
    linear between levels (z ground up), or None within the profile: one depth rule for both models' mixed layers."""
    z, thv = np.asarray(z, float), np.asarray(thv, float)
    ref = float(np.interp(z0, z, thv))
    above = z > z0
    zz, tt = np.concatenate([[z0], z[above]]), np.concatenate([[ref], thv[above]]) - ref
    k = np.nonzero(tt >= rise)[0]
    if k.size == 0:
        return None
    k = int(k[0])
    return float(zz[k - 1] + (rise - tt[k - 1]) / (tt[k] - tt[k - 1]) * (zz[k] - zz[k - 1]))


def compare(gcm: dict, cm1: dict) -> dict:
    """Side by side by quarter at the GCM's half levels below TOP_KM, the ground up: the heat each model carries up
    (W/m2) and its share of the surface's sensible heat; CM1's vapour (W/m2), its share of the surface's latent heat
    and its effective diffusivity (the flux over the mean gradient, where vapour falls with height by at least
    0.05 g/kg per km); the GCM's diffusivity and stability on its means; the mixed layer's depth by one rule in both;
    and the GCM's profiles beside CM1's at the GCM's full-level heights."""
    out = {}
    for name, _, _ in QUARTERS:
        g, m = gcm['quarters'][name], cm1['quarters'][name]
        zh = g['z_half_m']
        keep = np.where(zh <= TOP_KM * 1000.0)[0][::-1]          # the ground up
        zf = g['z_full_m']
        keepf = np.where(zf <= TOP_KM * 1000.0)[0][::-1]
        at = lambda v: np.interp(zh[keep], cm1['zw_m'], v)
        atf = lambda v: np.interp(zf[keepf], cm1['zh_m'], v)
        r = lambda v, n=1: [None if x is None or not np.isfinite(x) else round(float(x), n) for x in np.asarray(v, float)]
        g_sens, m_sens = float(g['sensible_w_m2']), float(m['sensible_w_m2'])
        g_sens_snap = float(g['surface_layer_w_m2'])
        m_lat = float(m['latent_w_m2'])
        g_heat = g['heat_vdiff_w_m2'][keep] + g['heat_dca_w_m2'][keep]
        grad = at(m['dqvdz_per_m'])
        k_eff = np.where(grad < -GRADIENT_FLOOR, -at(m['vapour_kinematic']) / np.where(grad < 0.0, grad, -1.0), np.nan)
        z0 = float(zf[keepf][0])
        out[name] = dict(
            z_half_m=r(zh[keep], 0),
            gcm=dict(outputs=g['outputs'], snapshots=g['snapshots'], sensible_w_m2=round(g_sens, 1),
                     latent_w_m2=round(float(g['latent_w_m2']), 1),
                     sensible_at_snapshots_w_m2=round(g_sens_snap, 1),
                     heat_w_m2=r(g_heat), heat_vdiff_w_m2=r(g['heat_vdiff_w_m2'][keep]),
                     heat_dca_w_m2=r(g['heat_dca_w_m2'][keep]),
                     heat_share_of_sensible=r(g_heat / g_sens_snap, 3) if abs(g_sens_snap) > 1.0 else None,
                     heat_vdiff_on_means_w_m2=r(g['heat_vdiff_on_means_w_m2'][keep]),
                     k_m2_s=r(g['k_m2_s'][keep]), k_median_m2_s=r(g['k_median'][keep]),
                     k_stretched_m2_s=r(g['k_stretched_m2_s'][keep]),
                     mixing_length_m=r(g['l_m'][keep], 0), mixing_length_stretched_m=r(g['l_stretched_m'][keep], 0),
                     shear_per_s=[float(f'{x:.3g}') for x in g['shear_s'][keep]],
                     ri_median=r(g['ri_median'][keep], 2), ri_share_unstable=r(g['ri_share_unstable'][keep], 2),
                     stability_factor=r(g['f'][keep], 2),
                     mixed_depth_m=rise_height(zf[keepf], g['theta_v_k'][keepf], z0),
                     z_full_m=r(zf[keepf], 0), theta_v_k=r(g['theta_v_k'][keepf], 2), q_g_kg=r(g['q_g_kg'][keepf], 2)),
            cm1=dict(outputs=m['outputs'], sensible_w_m2=round(m_sens, 1), latent_w_m2=round(m_lat, 1),
                     surface_heat_w_m2=round(float(m['heat_w_m2'][0]), 1),
                     surface_vapour_w_m2=round(float(m['vapour_w_m2'][0]), 1),
                     heat_w_m2=r(at(m['heat_w_m2'])), heat_resolved_w_m2=r(at(m['heat_resolved_w_m2'])),
                     heat_parameterized_w_m2=r(at(m['heat_parameterized_w_m2'])),
                     heat_share_of_sensible=r(at(m['heat_w_m2']) / m_sens, 3) if abs(m_sens) > 1.0 else None,
                     vapour_w_m2=r(at(m['vapour_w_m2'])), vapour_resolved_w_m2=r(at(m['vapour_resolved_w_m2'])),
                     vapour_parameterized_w_m2=r(at(m['vapour_parameterized_w_m2'])),
                     vapour_share_of_latent=r(at(m['vapour_w_m2']) / m_lat, 3) if abs(m_lat) > 1.0 else None,
                     k_effective_vapour_m2_s=r(k_eff, 0), khv_m2_s=r(at(m['khv_m2_s'])),
                     sigma_w_m_s=r(at(m['sigma_w_m_s']), 2), hpbl_m=round(float(m['hpbl_m']), 0),
                     mixed_depth_m=rise_height(cm1['zh_m'], m['theta_v_k'], z0),
                     theta_v_k=r(atf(m['theta_v_k']), 2), q_g_kg=r(atf(m['q_g_kg']), 2)))
    return out


def analyse(run: str = RUN, years=YEARS) -> dict:
    sites = {}
    for site, box, case in SITES:
        gcm = gcm_site(run, years, box)
        cm1 = cm1_box(case)
        sites[site] = dict(gcm_cells=gcm['cells'], cm1_case=case, cm1_span_days=cm1['span_days'],
                           cm1_l_inf_m=cm1['l_inf_m'], quarters=compare(gcm, cm1))
        stretch = gcm['stretch']
    return dict(
        schema=SCHEMA, run=run, years=list(years), stretch=stretch,
        plasim=dict(lambda_momentum_m=LAMBDA_M, lambda_heat_m=round(LAMBDA_M * math.sqrt(3.0 * D * 0.5), 1),
                    b=B, c=C, d=D, wind_floor_m_s=UMIN, height_temperature_k=TSCAL),
        sites=sites,
        evidence=('The GCM: 3-day means of one two-year run in the site cells (heating by process as PlaSim writes it; '
                  'the diffusivity recomputed with vdiff\'s formula on the means, which differ from the means of the '
                  'step-by-step values). CM1: domain means of one realization per box over its settled lunar day.'),
        reading_rule=('By site and quarter of the lunar day (hour angle from noon), at the GCM\'s half levels (z_half_m '
                      'above the ground, the ground up): heat and vapour carried upward through each (W/m2). '
                      'gcm.heat_w_m2 is vertical diffusion plus dry adjustment from the run\'s own heating, which '
                      'PlaSim writes as single time steps at the end of each 3-day mean, binned here by the time of '
                      'those steps (half a mean\'s length after its midpoint), with its share of the surface layer\'s '
                      'heating at the same steps (sensible_at_snapshots_w_m2); every other GCM value is a 3-day mean '
                      'binned by its midpoint. '
                      'heat_vdiff_on_means_w_m2 recomputes the diffusion on the 3-day means as a check, and k_m2_s '
                      '(the mean over the means), k_median_m2_s, k_stretched_m2_s (mixing lengths times stretch), '
                      'ri_median and stability_factor place each half '
                      'level\'s regime on the means. cm1 values are interpolated to the same heights, resolved and '
                      'parameterized parts beside the total; k_effective_vapour_m2_s is the vapour flux over the mean '
                      'gradient where vapour falls by at least 0.05 g/kg per km. Shares divide by the surface\'s '
                      'sensible or latent heat. mixed_depth_m: the first height above the GCM\'s lowest full level '
                      'where the mean virtual potential temperature has risen 0.5 K above its value there. Profiles '
                      'at the GCM\'s full levels (z_full_m), theta_v referred to 1000 hPa, vapour as a mixing ratio '
                      'in both models.'),
        producer=dict(domain='climate', files={
            name: hashlib.sha256((HERE.parent / name).read_bytes()).hexdigest()[:16]
            for name in ('gcm/pbl_mixing.py', 'gcm/sbm_offline.py', 'crm/gcm_site.py', 'crm/box_budget.py',
                         'crm/ring_analysis.py')}))


def jsonable(o):
    if isinstance(o, dict):
        return {k: jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return jsonable(o.tolist())
    if isinstance(o, (np.floating, np.integer)):
        return o.item()
    return o


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    p.add_argument('--run', default=RUN)
    p.add_argument('--years', type=int, nargs=2, default=list(YEARS))
    p.add_argument('--output', type=Path)
    a = p.parse_args(argv)
    res = analyse(a.run, tuple(a.years))
    text = json.dumps(jsonable(res), indent=1)
    if a.output:
        a.output.write_text(text + '\n')
    else:
        print(text)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
