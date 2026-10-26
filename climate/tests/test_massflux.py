"""The mass-flux convection scheme (massflux.py) on made-up lunar columns: its saturation is PlaSim's, its plume's
mixing solves its equation exactly, its saturated state keeps its conserved variables, the column keeps its energy
and water to round-off with the latent heat of what rains out, a moist boundary layer under a stable layer makes a
shallow cloud capped below that layer while an unstable column convects deep, a dry stable column stays quiet, the
closures scale as their equations say, PlaSim's evaporation of falling rain keeps the water it moves, and a box's
column laid on PlaSim's layers keeps its water."""
import math

import numpy as np
import pytest

from climate.gcm import massflux as mf
from climate.gcm import sbm

G = 1.6242
C = dict(mf.PLASIM, ga=G, gascon=288.38, acpd=288.38 / 0.286)
C['rdbrv'] = C['gascon'] / 461.51


def lunar_column(stable_layer=False, n=60, top_m=100000.0, ps=121500.0, rh_ground=0.75, rh=0.8, dry=False):
    """A column on n levels from the ground to top_m, top first: a dry-adiabatic mixed layer from 300 K at the ground,
    its vapour rh_ground of saturation there, up to that air's condensation level; above it the temperature falls
    1 K/km to 220 K (conditionally unstable on the Moon, whose moist adiabat falls 0.66 K/km near the ground) at rh of
    saturation, or with stable_layer the potential temperature first rises 1 K/km over 6 km; vapour a tenth of these
    when dry; pressure hydrostatic at lunar gravity."""
    kap = C['gascon'] / C['acpd']
    theta0 = 300.0 * (1.0e5 / ps) ** kap
    q0 = rh_ground * mf.qsat(300.0, ps, C)
    fine = np.linspace(0.0, top_m, 40001)
    temp, lnp = np.empty(fine.size), np.empty(fine.size)
    lnp[0], temp[0] = math.log(ps), 300.0
    lcl = None
    for i in range(1, fine.size):
        z, dz = fine[i], fine[i] - fine[i - 1]
        lnp[i] = lnp[i - 1] - dz * G / (C['gascon'] * temp[i - 1])
        p = math.exp(lnp[i])
        if lcl is None:
            temp[i] = theta0 * (p / 1.0e5) ** kap
            if mf.qsat(temp[i], p, C) <= q0:
                lcl = (z, temp[i], theta0 * 1.0)
        elif stable_layer and z < lcl[0] + 6000.0:
            temp[i] = (theta0 + 1.0e-3 * (z - lcl[0])) * (p / 1.0e5) ** kap
        else:
            temp[i] = max(temp[i - 1] - 1.0e-3 * dz, 220.0)
    zh = np.linspace(0.0, top_m, n + 1)
    zf = 0.5 * (zh[1:] + zh[:-1])
    t = np.interp(zf, fine, temp)
    pf = np.exp(np.interp(zf, fine, lnp))
    ph = np.exp(np.interp(zh, fine, lnp))
    q = np.where(zf < lcl[0], q0, [rh * mf.qsat(a, b, C) for a, b in zip(t, pf)])
    if dry:
        q = 0.1 * q
    return dict(t=t[::-1], q=q[::-1], p_full=pf[::-1], p_half=np.r_[0.0, ph[::-1][1:]], z_full=zf[::-1],
                z_half=zh[::-1], lcl_m=lcl[0])


def run(col, **kw):
    args = dict(dt=1800.0, shf=10.0, lhf=100.0, c=C, z_full=col['z_full'], z_half=col['z_half'])
    args.update(kw)
    return mf.column(col['t'], col['q'], col['p_full'], col['p_half'], G, **args)


def test_saturation_is_plasims():
    for t, p in ((300.0, 1.2e5), (273.0, 9.0e4), (230.0, 3.0e4)):
        assert mf.qsat(t, p) == pytest.approx(sbm.saturation(t, p), rel=1e-12)
        h = 1e-4
        assert mf.dqsat_dt(t, p) == pytest.approx((mf.qsat(t + h, p) - mf.qsat(t - h, p)) / (2 * h), rel=1e-5)


def test_the_plume_mixing_solves_its_equation():
    phi, ea, eb, eps, dz = 340.0, 320.0, 316.0, 2.0e-4, 3000.0
    fine = phi
    steps = 30000
    for i in range(steps):                      # fine explicit steps of d phi/dz = eps (phi_env - phi)
        z = (i + 0.5) / steps * dz
        fine += eps * (ea + (eb - ea) * z / dz - fine) * dz / steps
    assert mf.mix(phi, ea, eb, eps, dz) == pytest.approx(fine, abs=1e-4)
    assert mf.mix(phi, ea, eb, 0.0, dz) == phi
    assert mf.mix(phi, ea, ea, 1.0, dz) == pytest.approx(ea)


def test_the_saturated_state_keeps_its_conserved_variables():
    z, p = 20000.0, 8.0e4
    for sl, qt in ((C['acpd'] * 280.0 + G * z - 1.0e4, 0.012), (C['acpd'] * 262.0 + G * z, 0.006)):
        t, qv, l = mf.adjust(sl, qt, z, p, G, c=C)
        assert qv + l == pytest.approx(qt, rel=1e-12)
        assert C['acpd'] * t + G * z - mf.latent(t, c=C) * l == pytest.approx(sl, abs=1e-5)
        if l > 0.0:
            assert qv == pytest.approx(mf.qsat(t, p, C), rel=1e-10)
    t, qv, l = mf.adjust(C['acpd'] * 290.0 + G * z, 1.0e-4, z, p, G, c=C)
    assert l == 0.0 and qv == 1.0e-4 and t == pytest.approx(290.0)
    assert mf.ice_share(280.0) == 0.0 and mf.ice_share(240.0) == 1.0 and 0.0 < mf.ice_share(262.0) < 1.0


def column_budgets(col, r):
    m = np.diff(col['p_half']) / G
    heat = float(np.sum(C['acpd'] * r['dtdt'] * m))
    water = float(np.sum(r['dqdt'] * m))
    return heat, water


def test_the_column_keeps_its_energy_and_water():
    col = lunar_column()
    for melt in (True, False):
        r = run(col, melt=melt)
        assert r['type'] == 'deep' and r['m_base'] > 0.0
        heat, water = column_budgets(col, r)
        if melt:
            rain, snow = r['rain_surface'], r['snow_surface']
        else:                                   # what the scheme hands PlaSim's mkrain, by the phase it formed in
            rain, snow = float(r['net_rain'].sum()), float(r['net_snow'].sum())
        assert rain + snow > 0.0
        assert water == pytest.approx(-(rain + snow), rel=1e-9, abs=1e-15)
        assert heat == pytest.approx(C['ALV'] * rain + C['ALS'] * snow, rel=1e-9, abs=1e-9)


def test_a_moist_layer_under_a_stable_layer_makes_a_capped_shallow_cloud():
    capped = lunar_column(stable_layer=True, n=200)            # levels 0.5 km apart, so the thin cloud spans some
    r = run(capped, shf=10.0, lhf=150.0, bl_dqdt=None)
    assert r['type'] == 'shallow'
    assert r['shallow_trial']['depth'] < mf.PARAMETERS['deep_dp']
    assert capped['lcl_m'] < r['shallow_trial']['top_z'] < capped['lcl_m'] + 6000.0
    deep = run(lunar_column())
    assert deep['type'] == 'deep' and deep['trial']['depth'] > mf.PARAMETERS['deep_dp']
    assert deep['z_top'] > capped['lcl_m'] + 6000.0


def test_a_dry_stable_column_stays_quiet():
    col = lunar_column(dry=True)
    r = run(col, omega=np.zeros(60))
    assert r['type'] is None and r['m_base'] == 0.0 and not np.any(r['dtdt'])


def test_the_closures_scale_as_their_equations_say():
    col = lunar_column()
    base = run(col, prm=dict(mf.PARAMETERS, t_star=1e12))
    slow = run(col, prm=dict(mf.PARAMETERS, tau_factor=2.0 * mf.PARAMETERS['tau_factor'], t_star=1e12))
    if not slow['cfl_limited'] and not base['cfl_limited']:
        assert slow['m_base'] == pytest.approx(0.5 * base['m_base'], rel=1e-9)
    # heating the sub-cloud layer raises PCAPEbl and lowers the deep mass flux; cooling it does the opposite
    heat = np.zeros(60)
    heat[-5:] = 1.0e-5
    warmed = run(col, bl_dtdt=heat, u_bl=3.0)
    cooled = run(col, bl_dtdt=-heat, u_bl=3.0)
    plain = run(col, u_bl=3.0)
    assert warmed['pcape_bl'] > 0.0 > cooled['pcape_bl']
    assert warmed['m_base'] < plain['m_base'] < cooled['m_base']
    # shallow: the base mass flux is the sub-cloud layer's moist static energy supply over the base's excess
    capped = lunar_column(stable_layer=True, n=200)
    supply = np.zeros(200)
    supply[-1] = 2.0e-4
    r = run(capped, shf=10.0, lhf=150.0, bl_dqdt=supply)
    assert r['type'] == 'shallow' and r['supply'] > 0.0
    assert r['m_base'] == pytest.approx(min(r['supply'] / r['base_excess'], r['m_base']), rel=1e-9)


def test_rain_evaporation_moves_water_without_loss():
    col = lunar_column()
    r = run(col, melt=False)
    c = mf.Column(col['t'], col['q'], col['p_full'], col['p_half'], G, C, col['z_full'], col['z_half'])
    for gamma in (0.0, 0.03, 10.0):
        e = mf.evaporate(c, r['net_rain'], r['net_snow'], r['dtdt'], r['dqdt'], gamma)
        handed = float((r['net_rain'] + r['net_snow']).sum())
        evaporated = float(np.sum(e['dqdt'] * c.mass))
        assert e['rain'] + e['snow'] + evaporated == pytest.approx(handed, rel=1e-9)
        if gamma == 0.0:
            assert evaporated == 0.0


def test_a_box_column_laid_on_plasims_layers_keeps_its_water():
    col = lunar_column(n=111, top_m=150000.0)
    x = dict(t=col['t'], q=col['q'], p_full=col['p_full'], p_half=np.r_[col['p_full'][0] * 0.9, col['p_half'][1:]],
             bl_dtdt=np.zeros(111), bl_dqdt=np.zeros(111), conv_dtdt=np.zeros(111), conv_dqdt=np.zeros(111),
             omega=np.zeros(111), rad_dtdt=np.zeros(111), g=G, shf=0.0, lhf=0.0, u_bl=1.0)
    half = mf.grid_sets()['10']
    out = mf.to_grid(x, half)
    ph = x['p_half']
    water_box = float(np.sum(x['q'] * np.diff(ph)))
    target = np.r_[0.0, half] * ph[-1]
    covered = np.diff(np.maximum(target, ph[0]))
    assert float(np.sum(out['q'] * covered)) == pytest.approx(water_box, rel=1e-9)
    assert out['p_half'][-1] == pytest.approx(ph[-1]) and out['t'].size == half.size
