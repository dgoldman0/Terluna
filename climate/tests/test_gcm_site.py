"""The GCM site forcing: the spherical harmonics it differentiates with, the continuity that gives its vertical wind,
the fits by local time, the cells it takes and its table as CM1 reads it; on the GCM's own output where it is at hand,
the budget it holds."""
import numpy as np
import pytest
from scipy.special import assoc_legendre_p_all

from climate.crm import gcm_site as gs

NLAT, NLON, RADIUS = 32, 64, 1.7374e6
SEA = dict(lat_deg=-5.0, lon_deg=-45.0, half_lat_deg=5.6, half_lon_deg=8.5, surface='sea')
BRANCH = gs.HERE.parent / 'gcm' / 'runs' / 'A28R_301k_cold_moist' / 'model' / 'MOST.00000.nc'


def grid():
    nodes, _ = np.polynomial.legendre.leggauss(NLAT)
    mu = np.sort(nodes)[::-1]                                               # rows north to south
    return np.degrees(np.arcsin(mu)), np.arange(NLON) * 360.0 / NLON, mu


def test_the_harmonics_give_back_a_harmonic_and_its_derivatives():
    sp = gs.Spectral(NLAT, NLON)
    lat, lon, mu = grid()
    lam = np.radians(lon)
    p = np.asarray(assoc_legendre_p_all(21, 21, mu, norm=True, diff_n=1))
    for n, m in ((0, 0), (3, 1), (7, 4), (15, 2), (21, 21)):
        f = p[0, n, m][:, None] * np.cos(m * lam)[None, :]
        s = sp.analyse(f)
        scale = np.abs(p[1, n, m]).max() + 1.0
        assert np.allclose(sp.synthesise(s), f, atol=1e-10)
        assert np.allclose(sp.synthesise(s, 'lambda'), -m * p[0, n, m][:, None] * np.sin(m * lam)[None, :], atol=1e-9)
        assert np.allclose(sp.synthesise(s, 'mu'), p[1, n, m][:, None] * np.cos(m * lam)[None, :], atol=1e-10 * scale)
    # A harmonic beyond the truncation does not pass.
    beyond = p[0, 21, 3][:, None] * np.cos(3 * lam)[None, :]
    assert np.allclose(gs.Spectral(NLAT, NLON, 20).project(beyond), 0.0, atol=1e-10)


def test_continuity_turns_a_rotating_flow_away_and_finds_a_spreading_one():
    sp, wind = gs.Spectral(NLAT, NLON), gs.Spectral(NLAT, NLON, gs.TRUNCATION + 1)
    lat, lon, mu = grid()
    phi = np.radians(lat)
    sigma = np.array([0.1, 0.4, 0.7, 0.9])                                 # midway between 0, 0.2, 0.6, 0.8 and 1
    ps = np.full((1, NLAT, NLON), 1.2e5)
    shape = (1, sigma.size, NLAT, NLON)
    turning = np.broadcast_to((10.0 * np.cos(phi))[None, None, :, None], shape)
    assert np.abs(gs.omega(turning, np.zeros(shape), ps, sigma, RADIUS, sp, wind)).max() < 1e-9
    # v = cos(phi) sin(phi): div = (1/a) d(v cos)/d(mu) = (1 - 3 mu^2)/a at every level, so omega = -p div.
    spreading = np.broadcast_to((np.cos(phi) * np.sin(phi))[None, None, :, None], shape)
    expected = -(sigma[None, :, None, None] * ps[:, None]) * ((1.0 - 3.0 * mu**2) / RADIUS)[None, None, :, None]
    assert np.allclose(gs.omega(np.zeros(shape), spreading, ps, sigma, RADIUS, sp, wind), expected, rtol=1e-8,
                       atol=1e-14)


def test_the_fits_undo_a_36_degree_mean_and_give_the_rate():
    centres = np.arange(10) * 36.0 - 162.0
    truth = lambda h: 2.0 + 0.5 * np.cos(np.radians(h)) - 0.3 * np.sin(np.radians(3 * h)) + 0.1 * np.cos(np.radians(4 * h))
    damp = lambda n: np.sin(np.radians(18.0 * n)) / np.radians(18.0 * n)    # a 36-degree mean of the n-th harmonic
    means = (2.0 + 0.5 * damp(1) * np.cos(np.radians(centres)) - 0.3 * damp(3) * np.sin(np.radians(3 * centres))
             + 0.1 * damp(4) * np.cos(np.radians(4 * centres)))
    coef = gs.undamped(gs.harmonic_coefficients(centres, means), 36.0)
    hours = np.linspace(-180.0, 180.0, 37)
    assert np.allclose(gs.evaluate(coef, hours), truth(hours), atol=1e-12)
    day = 2.55e6
    expected = (-0.5 * np.sin(np.radians(hours)) - 0.9 * np.cos(np.radians(3 * hours))
                - 0.4 * np.sin(np.radians(4 * hours))) * 2.0 * np.pi / day
    assert np.allclose(gs.rate(coef, hours, day), expected, atol=1e-17)
    # Snapshots need no undoing: the plain fit gives them back.
    assert np.allclose(gs.evaluate(gs.harmonic_coefficients(centres, truth(centres)), hours), truth(hours), atol=1e-12)


def test_the_cells_are_the_sites_surface():
    lat, lon = np.array([10.0, 5.0, 0.0, -5.0, -10.0]), np.arange(0.0, 360.0, 5.0)
    lsm = np.zeros((5, 72), bool)
    lsm[:, :3] = True                                                        # land from 0 to 10 E
    height = np.zeros((5, 72))
    height[2, 71] = 300.0                                                    # a raised lake at 355 E
    site = dict(lat_deg=0.0, lon_deg=0.0, half_lat_deg=5.0, half_lon_deg=5.0)
    land = gs.site_cells(lat, lon, lsm, height, dict(site, surface='land'))
    sea = gs.site_cells(lat, lon, lsm, height, dict(site, surface='sea'))
    assert land.sum() == 6 and set(np.nonzero(land)[1]) == {0, 1}
    assert sea.sum() == 2 and set(np.nonzero(sea)[1]) == {71}                # the lake stands above the sea
    with pytest.raises(ValueError):
        gs.site_cells(lat, lon, lsm, height, dict(site, surface='ice'))
    with pytest.raises(ValueError):
        gs.site_cells(lat, lon, np.ones((5, 72), bool), height, dict(site, surface='sea'))


def test_cm1_reads_the_table_from_the_ground_up():
    forcing = dict(profile=dict(z_m=[30000.0, 9000.0, 800.0]), vertical_wind_m_s=[0.0, 0.004, 0.001],
                   day_night=dict(theta_k_s=[[0.0, 1e-5, 2e-5], [0.0, 3e-5, 4e-5]],
                                  qv_kg_kg_s=[[0.0, 1e-8, 2e-8], [0.0, 3e-8, 4e-8]], hour_angle_deg=[-90.0, 90.0]))
    out = gs.to_cm1(forcing, [0.0, 400.0, 800.0, 9000.0, 60000.0])
    assert out['w_m_s'][0] == 0.0 and out['w_m_s'][-1] == 0.0
    assert out['w_m_s'][1:4] == pytest.approx([0.0005, 0.001, 0.004])
    table = out['day_night']
    assert list(table['z_m']) == [0.0, 800.0, 9000.0, 30000.0]
    assert np.allclose(table['theta_k_s'], [[2e-5, 2e-5, 1e-5, 0.0], [4e-5, 4e-5, 3e-5, 0.0]])
    assert np.allclose(table['qv_kg_kg_s'][:, 0], table['qv_kg_kg_s'][:, 1])


def test_the_site_budget_holds_on_the_gcm():
    pytest.importorskip('netCDF4')
    if not BRANCH.exists():
        pytest.skip(f'{BRANCH} is not at hand')
    f = gs.site_forcing('A28R_301k_cold_moist', (0, 0), SEA)
    assert len(f['cells']) == 6 and f['outputs'] == 120 and f['phases'] == 10
    column = f['budget']['column_vapour_mm_day']
    assert np.allclose(column['held_by_hour'], column['exact_by_hour'], rtol=1e-9, atol=1e-9)
    # Over the day the circulation's heating balances the physics' exactly, the storage's mean being nil.
    b = f['budget']
    assert np.allclose(b['circulation_heating_theta_k_day'], -np.asarray(b['physics_heating_theta_k_day']), atol=1e-9)
    z = np.asarray(f['profile']['z_m'])
    above = z > gs.TOP_M
    assert np.all(np.asarray(f['vertical_wind_m_s'])[above] == 0.0)
    assert np.allclose(np.asarray(f['vertical_wind_m_s'])[~above], np.asarray(f['gcm_vertical_wind_m_s'])[~above])
    assert np.all(np.asarray(f['day_night']['theta_k_s'])[:, above] == 0.0)
    # The sea evaporates more than it rains, and the circulation carries the difference away.
    assert b['evaporation_mm_day'] > b['rain_mm_day'] and column['exact'] < 0.0
