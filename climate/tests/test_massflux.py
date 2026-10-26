"""The mass-flux convection scheme (massflux.py) on made-up lunar columns: its saturation is PlaSim's, its plume's
mixing solves its equation exactly, its saturated state keeps its conserved variables, the column keeps its energy
and water to round-off with the latent heat of what rains out (PlaSim's moist heat capacity converting the heating),
a moist boundary layer under a stable layer makes a shallow cloud capped below that layer while an unstable column
convects deep, a dry stable column stays quiet, the closures scale as their equations say, the surface parcel takes
the air near the ground from surface-layer similarity in a convective surface layer only, the cloud's top pressure
lies between the plume's points, PlaSim's mkrain (its copy) keeps the water it moves and never lets the running
water fall below nil, the phases handed to it follow its own rule, and a box's column laid on PlaSim's layers keeps
its water."""
import math

import numpy as np
import pytest

from climate.gcm import massflux as mf
from climate.gcm import sbm

G = 1.6242
TRAPEZOID = getattr(np, 'trapezoid', None) or getattr(np, 'trapz')
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
    heat = float(np.sum(C['acpd'] * (1.0 + C['ADV'] * col['q']) * r['dtdt'] * m))
    water = float(np.sum(r['dqdt'] * m))
    return heat, water


def test_the_column_keeps_its_energy_and_water():
    col = lunar_column()
    for melt in (True, False, 'plasim'):
        r = run(col, melt=melt)
        assert r['type'] == 'deep' and r['m_base'] > 0.0
        heat, water = column_budgets(col, r)
        if melt is True:
            rain, snow = r['rain_surface'], r['snow_surface']
        else:                                   # what the scheme hands PlaSim's mkrain, by phase
            rain, snow = float(r['net_rain'].sum()), float(r['net_snow'].sum())
        assert rain + snow > 0.0
        assert water == pytest.approx(-(rain + snow), rel=1e-9, abs=1e-15)
        assert heat == pytest.approx(C['ALV'] * rain + C['ALS'] * snow, rel=1e-9, abs=1e-9)


def test_the_phases_handed_on_follow_mkrains_rule():
    """mkrain adds each layer's handed rain and snow to its running totals, then turns all of them to rain above
    freezing and to snow at or below it: with the scheme's 'plasim' phases the totals never go negative."""
    col = lunar_column()
    r = run(col, melt='plasim')
    t = col['t']
    rain_run, snow_run = 0.0, 0.0
    for k in range(t.size):
        rain_run += r['net_rain'][k]
        snow_run += r['net_snow'][k]
        if t[k] > C['TMELT'] and snow_run > 0.0:
            rain_run, snow_run = rain_run + snow_run, 0.0
        elif t[k] <= C['TMELT'] and rain_run > 0.0:
            rain_run, snow_run = 0.0, rain_run + snow_run
        assert rain_run >= -1e-18 and snow_run >= -1e-18, k
    assert r['evap_down'].sum() > 0.0 and np.any(r['net_rain'] < 0.0)     # the downdraught takes from the rain


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
    r = run(col, melt='plasim')
    c = mf.Column(col['t'], col['q'], col['p_full'], col['p_half'], G, C, col['z_full'], col['z_half'])
    handed = float((r['net_rain'] + r['net_snow']).sum())
    for gamma in (0.0, 0.03, 10.0):
        e = mf.evaporate(c, r['net_rain'], r['net_snow'], r['dtdt'], r['dqdt'], gamma)
        evaporated = float(np.sum(e['dqdt'] * c.mass))
        assert e['rain'] + e['snow'] + evaporated == pytest.approx(handed, rel=1e-9)
        assert e['rain'] >= 0.0 and e['snow'] >= 0.0
        if gamma == 0.0:
            assert evaporated == 0.0
        if gamma == 10.0:                       # the floor: the downdraught's water is never taken twice
            assert e['rain'] + e['snow'] == pytest.approx(0.0, abs=1e-12 * handed)
    # snow handed into a warm layer reaches the ground as rain, and the layer pays the heat of melting
    snow = np.zeros(c.n)
    snow[0] = 1.0e-4
    e = mf.evaporate(c, np.zeros(c.n), snow, np.zeros(c.n), np.zeros(c.n), 0.0)
    assert e['snow'] == 0.0 and e['rain'] == pytest.approx(1.0e-4)
    melt_heat = float(np.sum(C['acpd'] * (1.0 + C['ADV'] * c.q) * e['dtdt'] * c.mass))
    assert melt_heat == pytest.approx(-(C['ALS'] - C['ALV']) * 1.0e-4, rel=1e-12)


def test_the_surface_parcel_takes_the_air_near_the_ground():
    col = lunar_column(n=60)
    c = mf.Column(col['t'], col['q'], col['p_full'], col['p_half'], G, C, col['z_full'], col['z_half'])
    # the IFS's lowest level, 10 m on Earth, stretched by the ratio of the gravities
    assert mf.source_height(c) == pytest.approx(10.0 * 9.80665 / G)
    earth = mf.Column(col['t'], col['q'], col['p_full'], col['p_half'], 9.80665, C, col['z_full'], col['z_half'])
    assert mf.source_height(earth) == pytest.approx(10.0)
    # the integrated stability function is the integral of (1 - phi_h) / zeta for phi_h = (1 - 16 zeta)^(-1/3)
    for zeta in (-0.1, -2.0, -30.0):
        x = np.linspace(zeta, 0.0, 200001)
        f = np.r_[(1.0 - (1.0 - 16.0 * x[:-1]) ** (-1.0 / 3.0)) / x[:-1], -16.0 / 3.0]     # its limit at 0
        assert mf.psi_heat(zeta) == pytest.approx(-TRAPEZOID(f, x), rel=1e-5)
    # heated, evaporating ground: the air near it is warmer and moister than the lowest level's (z_n 833 m here);
    # cooled ground: the parcel keeps the lowest level's air
    dt_l, dq_l = mf.surface_layer_air(c, 100.0, 200.0)
    assert dt_l > 0.0 and dq_l > 0.0
    assert mf.surface_layer_air(c, -30.0, 0.0) == (0.0, 0.0)
    assert mf.surface_layer_air(c, 100.0, 200.0, dict(mf.PARAMETERS, source_height_earth=0.0)) == (0.0, 0.0)
    # free convection: the difference hardly depends on the friction velocity
    fast = mf.surface_layer_air(c, 100.0, 200.0, dict(mf.PARAMETERS, ustar=0.2))
    assert fast[0] == pytest.approx(dt_l, rel=0.05)
    # the excess: 1.5 times the flux over rho sigma_w, sigma_w taken at the source height
    dt_s, dq_s, ws = mf.surface_excess(c, 100.0, 0.0)
    rho = c.ph[-1] / (C['gascon'] * mf.virtual(c.t[-1], c.q[-1], C))
    flux = 100.0 / C['acpd'] / rho
    expect = (1.2 * 0.1 ** 3 + 1.5 * 0.4 * G / c.t[-1] * mf.source_height(c) * flux) ** (1.0 / 3.0)
    assert ws == pytest.approx(expect, rel=1e-12) and dt_s == pytest.approx(1.5 * flux / ws, rel=1e-12)


def test_the_cloud_top_pressure_lies_between_the_plumes_points():
    points = [dict(z=0.0, p=1.0e5), dict(z=1000.0, p=9.0e4), dict(z=2000.0, p=8.0e4)]
    assert mf.top_pressure(points, 1500.0) == pytest.approx(math.sqrt(9.0e4 * 8.0e4))
    assert mf.top_pressure(points, 2000.0) == pytest.approx(8.0e4)


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
