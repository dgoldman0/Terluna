"""The entrainment diagnosis (box_entrainment.py) on made-up fields: a cloud population mixing at a known rate gives that
rate back from qt and h, the samples and their environment, the top of each column's lowest unbroken cloud, a plume
that loses its buoyancy sooner the more it mixes, and the GCM records' rain formed, heated top and vapour above the
scheme's reference."""
import numpy as np
import pytest

from climate.crm import box_entrainment as be
from climate.crm import box_onset as bo
from climate.crm import cm1_output as co
from climate.gcm import sbm

G = 1.6242


def mixing_snapshot(eps=1.0e-3, n_core=40, columns=400):
    """Levels 200 m apart from 6 to 9 km: an environment drying and warming with height, and n_core columns of
    saturated, rising, light cloud whose qt and h relax toward the environment's at eps per metre."""
    z = np.arange(6000.0, 9001.0, 200.0)
    nz = z.size
    qt_e = 0.011 - 1.0e-6 * (z - 6000.0)
    h_e = 3.30e5 - 1.0 * (z - 6000.0)
    qt_c, h_c = np.empty(nz), np.empty(nz)
    qt_c[0], h_c[0] = 0.012, 3.32e5
    for k in range(1, nz):
        a = np.exp(-eps * (z[k] - z[k - 1]))      # the exact solution with the environment taken midway
        qt_c[k] = 0.5 * (qt_e[k] + qt_e[k - 1]) + (qt_c[k - 1] - 0.5 * (qt_e[k] + qt_e[k - 1])) * a
        h_c[k] = 0.5 * (h_e[k] + h_e[k - 1]) + (h_c[k - 1] - 0.5 * (h_e[k] + h_e[k - 1])) * a
    core = np.zeros(columns, bool)
    core[:n_core] = True
    qt = np.where(core, qt_c[:, None], qt_e[:, None])
    h = np.where(core, h_c[:, None], h_e[:, None])
    x = dict(qt=qt, h=h, qcl=np.where(core, 2.0e-4, 0.0)[None, :].repeat(nz, axis=0),
             w=np.where(core, 1.0, -0.1)[None, :].repeat(nz, axis=0),
             thrho=np.where(core, 290.5, 290.0)[None, :].repeat(nz, axis=0),
             rho=np.ones((nz, columns)), t=np.full((nz, columns), 290.0), qv=qt - 2.0e-4 * core,
             p=np.full((nz, columns), 9.0e4), th=np.full((nz, columns), 285.0))
    return z, x


def test_a_population_mixing_at_a_known_rate_gives_that_rate_back():
    z, x = mixing_snapshot(eps=1.0e-3)
    m = be.masks(x)
    assert m['core'].sum(axis=1).tolist() == [40] * z.size == m['updraft'].sum(axis=1).tolist()
    s = be.pooled([be.snapshot_sums(x, m)] * 2)
    e = be.entrainment(s, z)
    for v in ('qt', 'h'):
        assert e['core'][f'bulk_{v}'] == pytest.approx(1.0e-3, rel=0.03)
        inner = e['core'][f'eps_{v}'][1:-1]                      # centred differences; the ends take one side
        assert np.all(np.abs(inner - 1.0e-3) < 2.0e-5)
    p = be.profiles(s, 'core')
    assert p['area'][0] == pytest.approx(0.1) and p['w'][0] == pytest.approx(1.0)
    assert p['buoyancy_k'][0] == pytest.approx(0.45)          # 290.5 against the level's mean of 290.05
    # a steady mass flux detrains what it entrains
    assert np.nanmax(np.abs(e['core']['delta_qt'] - e['core']['eps_qt'])) < 1e-12


def test_a_slower_mixing_population_reads_slower():
    z, fast = mixing_snapshot(eps=1.0e-3)
    _, slow = mixing_snapshot(eps=2.0e-4)
    rate = lambda x: be.entrainment(be.pooled([be.snapshot_sums(x, be.masks(x))] * 2), z)['core']['bulk_qt']
    assert rate(slow) == pytest.approx(2.0e-4, rel=0.02) and rate(fast) > 4.0 * rate(slow)


def test_levels_with_too_few_points_are_left_out():
    z, x = mixing_snapshot(n_core=10)
    e = be.entrainment(be.pooled([be.snapshot_sums(x, be.masks(x))] * 4), z)     # 40 points, under MIN_POINTS
    assert e['core']['bulk_qt'] is None and not np.isfinite(e['core']['eps_qt']).any()


def test_the_entrainment_derivative_takes_one_side_at_the_ends():
    z = np.array([0.0, 100.0, 200.0, 300.0])
    phi_c, phi_e = np.array([4.0, 3.0, 2.0, 1.0]), np.zeros(4)
    ok = np.array([False, True, True, True])
    eps = be.fractional(z, phi_c, phi_e, ok)
    assert np.isnan(eps[0])
    assert eps[1] == pytest.approx(0.01 / 3.0) and eps[2] == pytest.approx(0.005) and eps[3] == pytest.approx(0.01)


def test_the_cloud_top_is_the_top_of_each_columns_lowest_unbroken_cloud():
    zh = np.array([5000.0, 6000.0, 7000.0, 8000.0, 30000.0, 45000.0, 50000.0])
    q = np.zeros((zh.size, 4))
    q[1:3, 0] = 1e-4                          # 6-7 km
    q[1:4, 1], q[5, 1] = 1e-4, 1e-4           # 6-8 km, and a layer above a gap
    q[5:7, 2] = 1e-4                          # ice standing alone above 40 km
    tops = be.cloud_tops(q, zh)
    assert tops.tolist() == [7000.0, 8000.0]
    summary = be.tops_summary([tops], 4)
    assert summary['cloudy_share'] == 0.5 and summary['median_m'] == 7500.0 and summary['below_16km_share'] == 1.0


def made_up_sounding():
    """test_box_onset's sounding with an unsaturated boundary layer: 15 g/kg mixed below 2 km under air cooling
    1.4 K/km and drying to 40%."""
    z = np.arange(50.0, 60000.0, 100.0)
    t = np.where(z < 2000.0, 300.0 - G / co.CP * z, 300.0 - G / co.CP * 2000.0 - 1.4e-3 * (z - 2000.0))
    t = np.maximum(t, 180.0)
    lnp = np.log(121500.0) - np.concatenate([[0.0], np.cumsum(np.diff(z) * G / (co.RD * 0.5 * (t[1:] + t[:-1])))])
    p = np.exp(lnp)
    q = np.where(z < 2000.0, 0.015, 0.4 * bo.saturation(t, p))
    return z, p, t, q


def test_a_plume_stops_lower_the_more_it_mixes():
    z, p, t, q = made_up_sounding()
    caps = []
    for eps in (0.0, 1.0e-4, 3.0e-4, 1.0e-3, 0.05):
        c = be.cap(*be.plume(z, p, t, q, t[0], q[0], eps, G, fallout=True), z)
        assert c['lfc_m'] is not None and c['lcl_m'] <= c['lfc_m'] < c['cap_m'] and c['cin_j_kg'] < 1e-6
        caps.append(c['cap_m'])
    assert all(a > b for a, b in zip(caps, caps[1:]))
    assert caps[-1] - z[np.argmin(abs(z - caps[-1]))] < 300.0 and caps[-1] < 4000.0   # buoyant at condensation only


def test_a_plume_buoyant_nowhere_above_its_condensation_has_no_cap():
    z = np.arange(0.0, 10001.0, 1000.0)
    wet = z >= 3000.0
    b = np.where(z < 3000.0, 0.01, -0.002)                    # buoyant only below its condensation level
    out = be.cap(b, wet, z)
    assert out['lcl_m'] == 3000.0 and out['lfc_m'] is None and out['cap_m'] is None
    assert out['best_buoyancy_m_s2'] == pytest.approx(-0.002) and out['cape_j_kg'] == 0.0
    b = np.where(z < 5000.0, -0.002, np.where(z < 7000.0, 0.004, -0.004))
    out = be.cap(b, wet, z)                                   # through 5 km of inhibition to a layer at 5-6 km
    assert out['lfc_m'] == 5000.0 and out['cap_m'] == pytest.approx(6500.0) and out['cin_j_kg'] == pytest.approx(10.0)


def test_an_undiluted_plume_keeps_its_moist_static_energy_and_its_condensate_weighs_it_down():
    z, p, t, q = made_up_sounding()
    b_kept, wet = be.plume(z, p, t, q, t[0], q[0], 0.0, G)
    b_dropped, _ = be.plume(z, p, t, q, t[0], q[0], 0.0, G, fallout=True)
    assert wet.any() and np.all(b_kept[wet] <= b_dropped[wet] + 1e-12)
    k = int(np.flatnonzero(wet)[10])
    h0 = co.CP * t[0] + G * z[0] + bo.LV * q[0]
    tk, qk, _ = bo.saturate(h0, q[0], z[k], p[k], G)
    assert co.CP * tk + G * z[k] + bo.LV * qk == pytest.approx(h0, rel=1e-9)


def test_the_rain_formed_follows_the_heating_and_the_top_is_the_highest_heated_layer():
    c = dict(sbm.PLASIM, acpd=1004.6, ADV=0.0, ALV=2.5e6, ga=1.62)
    dsigma = np.array([0.5, 0.3, 0.2])
    dtc = np.array([0.0, 1.0e-5, 2.0e-5])
    rain = be.formed_rain(dtc, np.zeros(3), 1.0e5, dsigma, c)
    mass = 1.0e5 * dsigma / 1.62
    assert rain == pytest.approx(1004.6 * (1.0e-5 * mass[1] + 2.0e-5 * mass[2]) / 2.5e6 * 86400.0)
    assert be.heated_top(dtc) == 1 and be.heated_top(np.zeros(3)) == -1


def test_the_vapour_above_the_reference_counts_only_the_convecting_layers():
    sigma, dsigma = np.array([0.3, 0.7, 0.95]), np.array([0.4, 0.4, 0.1])
    t, tp = np.array([240.0, 280.0, 300.0]), np.array([245.0, 279.0, 300.0])
    q = np.array([1e-4, 6e-3, 0.016])
    out = be.vapour_above_reference(t, q, sigma, dsigma, 1.2e5, tp, 1, 0.7, sbm.PLASIM)
    assert out[0] == 0.0
    expect = (q[2] - 0.7 * sbm.saturation(300.0, 1.2e5 * 0.95)) * 1.2e5 * 0.1 / sbm.PLASIM['ga']
    assert out[2] == pytest.approx(expect)
    assert be.vapour_above_reference(t, q, sigma, dsigma, 1.2e5, tp, None, 0.7, sbm.PLASIM).tolist() == [0, 0, 0]


def test_a_gcm_column_on_finer_levels_keeps_its_own_values_at_its_levels():
    sigma = np.array([0.3, 0.6, 0.9])
    t, q = np.array([250.0, 275.0, 295.0]), np.array([5e-4, 4e-3, 0.012])
    zh = np.arange(50.0, 80000.0, 500.0)
    z, p, tt, rr = be.fine_column(t, q, 1.2e5, sigma, zh, G)
    assert p[0] == pytest.approx(1.2e5 * 0.9) and tt[0] == 295.0 and rr[0] == pytest.approx(0.012 / 0.988)
    assert np.all(np.diff(z) > 0.0) and np.all(np.diff(p) < 0.0) and tt[-1] > 250.0 - 1e-9


def test_the_record_keeps_five_figures_and_writes_nan_as_none():
    out = be.compact(dict(a=[1.23456789, np.float64(np.nan)], b=np.int64(3), c='x'))
    assert out == dict(a=[1.2346, None], b=3, c='x')
