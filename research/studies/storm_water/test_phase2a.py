"""Checks of phase 2a: the redistribution model's analytic limits, its agreement with phase 1 where they must agree,
mass conservation, the routing's gravity scalings, Lane's method and its recast, the bowls, the routing's water
balance, the Walnut Gulch basin check, and the product's phase-2a fields."""
import json
import math

import numpy as np
import pytest

from research.studies.storm_water import infiltration as IN
from research.studies.storm_water import phase2a as P2
from research.studies.storm_water import routing as RT
from research.studies.storm_water import run as R
from research.studies.storm_water import storms as ST
from research.studies.storm_water import wetting as W

MOON = IN.MOON
SATS = np.array([0.1, 0.2])


@pytest.fixture(scope='module')
def wg():
    return ST.walnut_gulch()


@pytest.fixture(scope='module')
def product():
    return json.loads(R.OUT.read_text())


@pytest.mark.parametrize('name', ['sandy_loam', 'loam', 'sealed_loam', 'biocrust'])
def test_one_storm_on_fresh_soil_is_phase_one(wg, name):
    p = IN.parameters(name, pressure_pa=1.16e5)
    c = W.curve(IN.SURFACES[name]['soil'])
    tot, _, _ = W.mean_storms(p, c, wg, SATS, 1, 70.0, 0.2)
    r1, _ = IN.single_storms(p, wg, SATS)
    assert np.allclose(tot[1], r1.mean(axis=0), rtol=1e-9, atol=1e-12)


@pytest.mark.parametrize('gravity', [1.0, MOON])
def test_redistribution_conserves_water(wg, gravity):
    p = IN.parameters('loam', gravity=gravity, pressure_pa=1.16e5)
    c = W.curve('loam', gravity)
    tot, layer, evap = W.mean_storms(p, c, wg, SATS, 6, 70.0, 0.3)
    rain = np.array([n * float(np.mean(ST.depths(wg))) for n in range(7)])
    # Every storm is drawn the same number of times at each position, so each N carries N mean storms of rain.
    assert np.allclose(tot + layer + evap, rain[:, None], rtol=1e-9)


def test_redistributed_runoff_lies_between_phase_one_bounds(wg):
    p = IN.parameters('sandy_loam', pressure_pa=1.16e5)
    c = W.curve('sandy_loam')
    tot, _, _ = W.mean_storms(p, c, wg, SATS, 8, 70.0, 0.2)
    seq, _ = IN.sequences(p, wg, SATS, 8)
    r1, _ = IN.single_storms(p, wg, SATS)
    for n in (2, 4, 8):
        assert n * r1.mean(axis=0)[0] * 0.98 <= tot[n, 0] <= seq[n, 0] * 1.02


def test_constant_drive_gravity_free_redistribution_matches_its_solution(monkeypatch):
    """With K = 0 and a constant drive D, 1/(se0 - se_i) = 1/(se0(0) - se_i) + p D theta_e t / I^2."""
    c = W.curve('sandy_loam')
    d = 500.0
    monkeypatch.setattr(W, 'drive', lambda c_, a, b: np.full(np.shape(b), d))
    monkeypatch.setattr(W, 'conductivity', lambda c_, se: np.zeros(np.shape(se)))
    pr = W.Profile(np.array([0.1]))
    pr.layer, pr.se0 = np.array([8.0]), np.array([1.0])
    hours = 3.0
    W.redistribute(c, pr, np.array([hours]), 0.0, p_shape=1.7, substeps=4000)
    x = 1.0 / (1.0 / 0.9 + 1.7 * d * c['theta_e'] * hours / 8.0 ** 2)
    assert abs((pr.se0[0] - 0.1) / x - 1.0) < 0.01
    assert pr.layer[0] == 8.0


def test_gravity_free_redistribution_depends_on_capillarity_alone():
    """Far below the Moon's gravity the conductivity's gravity part vanishes: two tiny gravities with the same Ks hb
    redistribute alike, and the capillary products are gravity's invariants."""
    out = []
    for g in (1e-4, 1e-6):
        c = W.curve('loam', g)
        pr = W.Profile(np.array([0.2]))
        pr.layer, pr.se0 = np.array([12.0]), np.array([1.0])
        W.redistribute(c, pr, np.array([24.0]), 0.0)
        out.append(pr.se0[0])
    assert abs(out[0] - out[1]) < 1e-4
    for soil in ('sandy_loam', 'loam'):
        e, m = W.curve(soil, 1.0), W.curve(soil, MOON)
        assert math.isclose(e['ks'] * e['hb'], m['ks'] * m['hb'], rel_tol=1e-12)
        se = np.array([0.3, 0.6])
        assert np.allclose(W.drive(e, 0.1, se), W.drive(m, 0.1, se)) and np.allclose(W.desorptivity_sq(e, se), W.desorptivity_sq(m, se))


def test_ponded_infiltration_without_gravity_is_sorptivity():
    k, psi, dth = 10.0, 100.0, 0.35
    for t in (0.1, 1.0, 10.0):
        f = RT.ponded_infiltration(t, k * 1e-6, psi / 1e-6, dth)
        assert abs(f / math.sqrt(2.0 * k * psi * dth * t) - 1.0) < 0.01


def test_lunar_wetted_layer_stays_shallower(wg):
    """The same water redistributes with the same capillarity but a sixth of the drainage on the Moon."""
    z = {}
    for g in (1.0, MOON):
        c = W.curve('sandy_loam', g)
        pr = W.Profile(np.array([0.15]))
        pr.layer, pr.se0 = np.array([40.0]), np.array([1.0])
        W.redistribute(c, pr, np.array([200.0]), 0.0)
        z[g] = pr.depth(c)[0]
    assert z[MOON] < z[1.0]


@pytest.mark.parametrize('gravity', [1.0, MOON])
def test_redistribution_against_richards_equation(gravity):
    """The rectangular profile with p = 1.7 against a numerical solution of Richards' equation for the same curve:
    the surface saturation within 0.03 and the depth holding nine tenths of the water within 10% over a day."""
    c = W.curve('loam', gravity)
    ref = W.richards_redistribution(c, 10.0, 0.2, [1.0, 6.0, 24.0], dz_mm=1.0, depth_mm=600.0)
    pr = W.Profile(np.array([0.2]))
    pr.layer, pr.se0 = np.array([10.0]), np.array([1.0])
    last = 0.0
    for t, (se_ref, z90) in zip((1.0, 6.0, 24.0), ref):
        W.redistribute(c, pr, np.array([t - last]), 0.0, substeps=400)
        last = t
        assert abs(pr.se0[0] - se_ref) < 0.03
        assert abs(pr.depth(c)[0] / z90 - 1.0) < 0.10


def test_substeps_converge(wg):
    p = IN.parameters('loam', pressure_pa=1.16e5)
    c = W.curve('loam')
    a, _, _ = W.mean_storms(p, c, wg, SATS, 6, 70.0, 0.2)
    b, _, _ = W.mean_storms(p, c, wg, SATS, 6, 70.0, 0.2, substeps=4 * W.SUBSTEPS)
    assert np.allclose(a[6], b[6], rtol=0.03)


def test_manning_with_gravity_explicit():
    """At the same discharge q = alpha h^(5/3): h ~ alpha^(-3/5), so 1.71 times as deep and 0.58 times as fast."""
    ratio = RT.manning_alpha(MOON) / RT.manning_alpha(1.0)
    depth = ratio ** -0.6
    assert abs(depth - 1.71) < 0.01 and abs(1.0 / depth - 0.584) < 0.005
    assert math.isclose(RT.wave_speed_m_s(MOON) / RT.wave_speed_m_s(1.0), 1.0 / depth, rel_tol=1e-9)


def test_lane_recast_and_handbook_values():
    d = np.array([0.1, 1.0, 10.0])
    assert np.allclose(RT.lane_kd_inches(d, 1.0), RT.BED['k_eff_in_h'] * d)
    moon = RT.lane_kd_inches(d, MOON) / (RT.BED['k_eff_in_h'] * d)
    assert np.all(np.diff(moon) < 0) and moon[0] < 1.0 and moon[-1] > MOON
    # Below the threshold volume nothing passes; above it the reach passes a + b P with b < 1.
    assert RT.lane_outflow(10.0, 2.0, 5.0, 3.0) == 0.0
    q1, q2 = RT.lane_outflow(np.array([4e4, 8e4]), 2.0, 5.0, 3.0)
    assert 0.0 < q1 < 4e4 and q2 - q1 < 4e4
    # A Walnut Gulch-like reach (6.65 km, 15 m wide, a 2 h flow of 50 acre-ft) loses 2-5 acre-ft a mile, as measured.
    loss = (50 * RT.ACRE_FT_M3 - RT.lane_outflow(50 * RT.ACRE_FT_M3, 6.65, 15.0, 6.0)) / RT.ACRE_FT_M3 / 4.13
    assert 2.0 < loss < 5.0


def test_bowls_follow_the_terrain_study():
    young, mature = RT.bowl_classes('young'), RT.bowl_classes('mature')
    c = math.sqrt(RT.EQUILIBRIUM_SHARE[0] * RT.EQUILIBRIUM_SHARE[1]) * RT.GEOMETRIC_SATURATION
    assert math.isclose(young['number_km2'].sum(), c * (0.01 ** -2 - 1.9 ** -2), rel_tol=1e-9)
    assert RT.bowl_cover(mature) < RT.bowl_cover(young) < 1.0
    assert np.allclose(RT.bowl_survival(young, np.array([0.0]), np.array([1e3])), 1.0)
    big = RT.bowl_survival(young, np.array([1.0]), np.array([1e9]))
    small = RT.bowl_survival(young, np.array([1.0]), np.array([1.0]))
    assert small < big <= 1.0


def test_hillslope_loss_bounds_and_gravity():
    for depth in (1.0, 5.0, 20.0):
        e, _ = RT.hillslope_loss(depth, 0.3, 10.0, 1.0)
        m, _ = RT.hillslope_loss(depth, 0.3, 10.0, MOON)
        assert 0.0 <= e <= m <= 1.0
    a, _ = RT.hillslope_loss(2.0, 0.3, 10.0, MOON)
    b, _ = RT.hillslope_loss(20.0, 0.3, 10.0, MOON)
    assert b < a


def test_routing_conserves_water_on_the_synthetic_basin():
    nxt, area, up, length, mouth = P2.synthetic_basin(1.9)
    order = RT.fronts(nxt)
    seen = np.zeros(nxt.size, bool)
    for f in order:
        assert not np.any(seen[nxt[f][nxt[f] >= 0]])          # receivers come after their donors
        seen[f] = True
    local = np.stack([np.full(nxt.size, 2e4), np.full(nxt.size, 2e3)], axis=1)
    lost, caught, out, _ = RT.route_events(order, nxt, local, area, length, up, np.array([0.5, 0.5]), MOON,
                                           [RT.bowl_classes('young'), None])
    for j in range(2):
        assert math.isclose(local[:, j].sum(), lost[:, j].sum() + caught[:, j].sum() + out[mouth, j], rel_tol=1e-9)
    assert caught[:, 1].sum() == 0.0


def test_walnut_gulch_basin_check_passes(wg):
    chk = P2.walnut_gulch_basin(dict(walnut_gulch=wg), 1.9)
    assert chk['passes']
    for soil in ('sandy_loam', 'loam'):
        assert chk[soil]['basin_yield_mm_yr'] < chk[soil]['point_runoff_mm_yr']


def test_product_phase2a_fields(product):
    assert product['schema'] == R.SCHEMA == 'terluna.research.storm-water-runoff/2'
    p2 = product['phase2a']
    for name in R.SURFACES:
        for b in R.BRACKETS:
            row = p2['wetting']['surfaces'][name][b]['classes']['semi-arid']
            assert row['low'] <= row['central'] <= row['high']
            ph1 = product['classes'][name][b]
            ind = ph1['independent']['classes']['semi-arid']['runoff_ratio']
            cont = ph1['continuing']['classes']['semi-arid']['runoff_ratio']
            assert ind * 0.98 <= row['central'] <= cont * 1.02
    for b in R.BRACKETS:
        for name in IN.MAIN:
            for s in P2.STATES:
                d = p2['routing']['fates'][b][name][s]['drylands']
                fates = ('bowls_inside', 'slope_soil', 'bowls_slope', 'channel_bed', 'bowls_channel', 'depressions',
                         'sea', 'humid_land')
                assert abs(sum(d[q]['share_of_runoff'] for q in fates) - 1.0) < 0.02
                assert d['bed_recharge']['share_of_runoff'] <= d['channel_bed']['share_of_runoff'] + 1e-9
    assert p2['walnut_gulch_basin']['passes']
    c0 = product['cells'][0]
    assert set(c0['routed_walnut_gulch']) == set(IN.MAIN)
    assert len(c0['runoff_ratio_2a']['loam']['walnut_gulch']['by_tenth']) == 10
