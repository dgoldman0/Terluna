"""PlaSim's vertical diffusion and dry adjustment as Python copies, and the nonlocal boundary layer beneath the diffusion
(climate/gcm/pbl_scheme.py): the copies keep what PlaSim's routines keep, the scheme's depth, velocity scales,
countergradient and entrainment fluxes follow their sources' forms in free convection, shear and stable air, and the
column keeps its heat, water and momentum."""
import math

import numpy as np
import pytest

from climate.gcm import levels as lv
from climate.gcm import pbl_scheme as ps

G_MOON = 1.6242188765608128
C = dict(ga=G_MOON, gascon=287.05, acpd=1004.6, rdbrv=287.05 / 461.51, ADV=0.8608)
AKAP = C['gascon'] / C['acpd']
HALF = np.asarray(lv.half_levels('l15'), float)
SIGMA = ps.full_sigma(HALF)
PS = 121500.0
DSIG = np.diff(np.r_[0.0, HALF])


def column(zi=4000.0, theta=300.0, lapse=1.0e-3, q_bl=0.012, q_scale=3000.0, wind=(8.0, 3.0), noise=0.0, seed=0):
    """A lunar column on the fifteen levels: a mixed layer zi deep under a stable layer (potential temperature rising
    lapse K/m), its vapour q_bl falling above with scale q_scale, the wind falling linearly from the top."""
    rng = np.random.default_rng(seed)
    z = -C['gascon'] * 285.0 / G_MOON * np.log(SIGMA)
    th = np.where(z < zi, theta, theta + lapse * (z - zi)) + rng.normal(0.0, noise, z.size)
    t = th * SIGMA ** AKAP
    q = np.where(z < zi, q_bl, q_bl * np.exp(-(z - zi) / q_scale))
    u = np.linspace(wind[0], wind[1], z.size)
    return t, q, u, np.zeros_like(u)


def kinematic(shf, lhf, rho=1.4):
    return dict(wth0=shf / (rho * C['acpd']), wq0=lhf / (rho * 2.5e6))


def test_the_lunar_length_is_earths_stretched_by_the_gravities():
    assert ps.LAMBDA_EARTH_M == 160.0 and ps.WIND_FLOOR == 1.0 and ps.STABILITY == 5.0
    assert ps.lunar_length(G_MOON) == pytest.approx(160.0 * 9.80665 / G_MOON) == pytest.approx(966.04, abs=0.01)
    assert ps.lunar_length(9.80665) == pytest.approx(160.0)


def test_in_neutral_air_vdiffs_diffusivity_is_the_mixing_length_squared_times_the_shear():
    t, q, u, v = column(zi=1.0e6, q_bl=0.0)                         # one potential temperature, no vapour
    km, kh, to_z = ps.local_coefficients(t, q, u, v, SIGMA, HALF, 3600.0, C)
    for k in range(SIGMA.size - 1):
        zz = -C['gascon'] * 250.0 / G_MOON * math.log(HALF[k])
        for lam, kk in ((160.0, km[k]), (160.0 * math.sqrt(7.5), kh[k])):
            l_mix = lam * 0.4 * zz / (lam + 0.4 * zz)
            shear = max(1.0, abs(u[k] - u[k + 1])) / (SIGMA[k + 1] - SIGMA[k]) * to_z[k]
            kin = l_mix ** 2 * shear
            assert kk == pytest.approx(3600.0 * kin * to_z[k] ** 2 / (SIGMA[k + 1] - SIGMA[k]), rel=1e-12)


@pytest.mark.parametrize('surface', [None, kinematic(8.0, 100.0), kinematic(-15.0, 0.0), kinematic(0.0, 0.0)])
def test_the_diffusion_keeps_the_columns_heat_water_and_momentum(surface):
    t, q, u, v = column(noise=0.2, seed=1)
    if surface is not None:
        surface = dict(surface, ustar=0.1)
    r = ps.vdiff(t, q, u, v, SIGMA, HALF, 3600.0, C, lamm=ps.lunar_length(G_MOON), surface=surface)
    for x, scale in ((r['dqdt'], np.abs(r['dqdt']).sum()), (r['dtdt_diffusion'], np.abs(r['dtdt_diffusion']).sum()),
                     (r['dudt'], np.abs(r['dudt']).sum())):
        assert abs(np.sum(DSIG * x)) <= 1e-12 * max(scale * DSIG.max(), 1e-300)
    # the dissipation's heating is the kinetic energy the column loses, no more
    mass = PS * DSIG / G_MOON
    ke_loss = -np.sum(mass * ((u + 3600.0 * r['dudt']) ** 2 - u ** 2) / 2.0) / 3600.0
    heat = np.sum(mass * C['acpd'] * (1.0 + C['ADV'] * q) * r['dtdt_dissipation'])
    assert heat == pytest.approx(ke_loss, rel=1e-9)


def test_in_free_convection_the_velocity_scale_and_prandtl_number_reach_their_limits():
    t, q, u, v = column()
    km, kh, to_z = ps.local_coefficients(t, q, u, v, SIGMA, HALF, 3600.0, C)
    s = kinematic(10.0, 100.0)
    p = ps.nonlocal_column(t, q, u, v, SIGMA, HALF, s['wth0'], s['wq0'], 0.0, km, kh, to_z, 3600.0, C)
    assert p['unstable'] and p['wm'] == pytest.approx(0.6 ** (1.0 / 3.0) * p['wstar'], rel=1e-12)
    assert p['pr'] == pytest.approx(0.1 * 0.4 * 7.2 / 0.6 ** (1.0 / 3.0), rel=1e-12)   # a k eps w*/w_m alone
    # w* from the buoyancy flux and the depth
    assert p['wstar'] ** 3 == pytest.approx(p['buoyancy'] * p['h'], rel=1e-12)
    # the depth lies at the mixed layer's top, between the levels bracketing 4 km, raised by the thermals' excess
    zf = p['z_full']
    assert 3500.0 < p['h'] < 6500.0 and p['excess'] > 0.0
    # the profile's diffusivity: k w z (1 - z/h)^2 inside h, nothing above
    inside = p['z_half'] < p['h']
    z = p['z_half'][inside]
    assert np.allclose(p['km_profile'][inside], 0.4 * p['wm'] * z * (1.0 - z / p['h']) ** 2, rtol=1e-12)
    assert np.allclose(p['kh_profile'][inside], p['km_profile'][inside] / p['pr'], rtol=1e-12)
    assert not np.any(p['km_profile'][~inside])
    assert zf[-1] < p['h']


def test_with_shear_the_scales_follow_u_star_and_the_prandtl_number_rises_toward_one():
    t, q, u, v = column()
    km, kh, to_z = ps.local_coefficients(t, q, u, v, SIGMA, HALF, 3600.0, C)
    s = kinematic(10.0, 100.0)
    calm = ps.nonlocal_column(t, q, u, v, SIGMA, HALF, s['wth0'], s['wq0'], 0.0, km, kh, to_z, 3600.0, C)
    windy = ps.nonlocal_column(t, q, u, v, SIGMA, HALF, s['wth0'], s['wq0'], 0.5, km, kh, to_z, 3600.0, C)
    assert windy['wm'] ** 3 == pytest.approx(0.5 ** 3 + 0.6 * windy['wstar'] ** 3, rel=1e-12)
    assert calm['pr'] < windy['pr'] < 1.0 + 0.1 * 0.4 * 7.2


def test_the_countergradient_and_entrainment_fluxes_carry_heat_up_inside_h_and_down_at_its_top():
    t, q, u, v = column()
    km, kh, to_z = ps.local_coefficients(t, q, u, v, SIGMA, HALF, 3600.0, C)
    s = kinematic(10.0, 100.0)
    p = ps.nonlocal_column(t, q, u, v, SIGMA, HALF, s['wth0'], s['wq0'], 0.1, km, kh, to_z, 3600.0, C)
    inside = p['z_half'] < p['h']
    assert np.all(p['flux_countergradient_t'][inside] > 0.0) and np.all(p['flux_countergradient_q'][inside] > 0.0)
    assert not np.any(p['flux_countergradient_t'][~inside])
    ja = p['inversion_interface']
    zf = p['z_full']
    assert zf[ja] > p['h'] >= zf[ja + 1]                          # the interface between the levels bracketing h
    ft = p['flux_entrainment_t']
    # at h the virtual flux is -0.15 times the surface buoyancy flux (and 0.75 theta_v u*^3 / (g h) more)
    th = t / SIGMA ** AKAP
    thv = th * (1.0 + (1.0 / C['rdbrv'] - 1.0) * q)
    flux_thv = -p['entrainment_velocity'] * (thv[ja] - thv[ja + 1])
    expect = -0.15 * p['wthv0'] - 0.15 * 5.0 * 0.1 ** 3 * thv[-1] / (G_MOON * p['h'])
    assert flux_thv == pytest.approx(expect, rel=1e-12)
    assert ft[ja] < 0.0 and not np.any(ft[:ja])                    # down at the inversion, none above
    below = np.arange(ja + 1, SIGMA.size - 1)
    assert np.allclose(ft[below], ft[ja] * (p['z_half'][below] / p['h']) ** 3, rtol=1e-12)
    assert p['flux_entrainment_q'][ja] > 0.0                       # moist air out of the mixed layer
    # the explicit increments: the air above h cools, the mixed layer warms, and the sums of sigma times T and q stay
    dt_inc, dq_inc = p['dt_explicit'], p['dq_explicit']
    assert dt_inc[ja] < 0.0 and abs(np.sum(DSIG * dt_inc)) < 1e-13 and abs(np.sum(DSIG * dq_inc)) < 1e-17


def test_in_stable_air_the_scheme_has_no_excess_countergradient_or_entrainment():
    t, q, u, v = column(zi=0.0, lapse=3.0e-3)
    km, kh, to_z = ps.local_coefficients(t, q, u, v, SIGMA, HALF, 3600.0, C)
    s = kinematic(-20.0, 0.0)
    p = ps.nonlocal_column(t, q, u, v, SIGMA, HALF, s['wth0'], s['wq0'], 0.3, km, kh, to_z, 3600.0, C)
    assert not p['unstable'] and p['excess'] == 0.0 and p['pr'] == 1.0
    assert not np.any(p['dt_explicit']) and not np.any(p['dq_explicit'])
    assert np.array_equal(p['km_profile'], p['kh_profile'])
    calm = ps.nonlocal_column(t, q, u, v, SIGMA, HALF, s['wth0'], s['wq0'], 0.0, km, kh, to_z, 3600.0, C)
    assert not np.any(calm['km_profile']) and np.array_equal(calm['km'], km)
    # neutral air: w_t = u*
    n = ps.nonlocal_column(t, q, u, v, SIGMA, HALF, 0.0, 0.0, 0.3, km, kh, to_z, 3600.0, C)
    inside = n['z_half'] < n['h']
    z = n['z_half'][inside]
    assert np.allclose(n['km_profile'][inside], 0.4 * 0.3 * z * (1.0 - z / n['h']) ** 2, rtol=1e-12)


def test_the_profile_raises_the_diffusivity_only_inside_h():
    t, q, u, v = column()
    km, kh, to_z = ps.local_coefficients(t, q, u, v, SIGMA, HALF, 3600.0, C)
    s = kinematic(10.0, 100.0)
    p = ps.nonlocal_column(t, q, u, v, SIGMA, HALF, s['wth0'], s['wq0'], 0.1, km, kh, to_z, 3600.0, C)
    inside = p['z_half'] < p['h']
    assert np.array_equal(p['km'][~inside], km[~inside]) and np.array_equal(p['kh'][~inside], kh[~inside])
    assert np.all(p['kh'][inside] >= kh[inside]) and np.any(p['kh'][inside] > kh[inside])


def test_a_profile_the_explicit_step_would_leave_negative_keeps_its_vapour():
    t, q, u, v = column(q_bl=1.0e-9, q_scale=1.0e9)
    km, kh, to_z = ps.local_coefficients(t, q, u, v, SIGMA, HALF, 3600.0, C)
    s = kinematic(10.0, 300.0)
    p = ps.nonlocal_column(t, q, u, v, SIGMA, HALF, s['wth0'], s['wq0'], 0.1, km, kh, to_z, 3600.0, C)
    assert p['vapour_kept'] and not np.any(p['dq_explicit']) and np.any(p['dt_explicit'])


def test_the_dry_adjustment_mixes_an_unstable_layer_and_keeps_heat_and_water():
    t, q, u, v = column(zi=0.0, lapse=3.0e-3)
    t0, q0 = t.copy(), q.copy()
    t[-1] += 3.0                                                  # a superadiabatic lowest layer
    q[-1] += 0.004
    dtd, dqd, passes = ps.dca(t, q, SIGMA, HALF, AKAP, 3600.0)
    assert passes >= 1
    tn, qn = t + 3600.0 * dtd, q + 3600.0 * dqd
    th = tn / SIGMA ** AKAP
    assert np.all(np.diff(th) <= 1e-12)                           # no layer's potential temperature below the next's
    assert np.sum(DSIG * tn) == pytest.approx(np.sum(DSIG * t), rel=1e-14)
    assert np.sum(DSIG * qn) == pytest.approx(np.sum(DSIG * q), rel=1e-14)
    assert np.isclose(qn[-1], qn[-2]) and np.array_equal(qn[:-3], q[:-3])
    dtd0, dqd0, passes0 = ps.dca(t0, q0, SIGMA, HALF, AKAP, 3600.0)
    assert passes0 == 0 and not np.any(dtd0) and not np.any(dqd0)


def test_full_levels_lie_midway_in_sigma_and_heights_follow_the_hypsometric_relation():
    assert np.allclose(SIGMA, 0.5 * (np.r_[0.0, HALF[:-1]] + HALF))
    t = np.full(SIGMA.size, 280.0)
    zf, zh = ps.heights(t, np.zeros(SIGMA.size), SIGMA, HALF, G_MOON, C['gascon'], C['rdbrv'])
    assert np.allclose(zh[:-1], C['gascon'] * 280.0 / G_MOON * np.log(1.0 / HALF[:-1])) and zh[-1] == 0.0
    assert np.allclose(zf, C['gascon'] * 280.0 / G_MOON * np.log(1.0 / SIGMA))
