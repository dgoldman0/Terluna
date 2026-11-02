"""Tests of the storm-water study's phase 2b: the pools' curves, energy balance, seepage and water balance, the slopes'
wear scalings, the Earth checks' recorded rules and the products."""
from __future__ import annotations

import json
import math

import numpy as np
import pytest

from research.studies.storm_water import infiltration as IN
from research.studies.storm_water import phase2b as B
from research.studies.storm_water import routing as RT
from research.studies.storm_water import run as R

PRODUCT = R.OUT


@pytest.fixture(scope='module')
def product():
    if not PRODUCT.is_file():
        pytest.skip('run python -m research.studies.storm_water.run')
    return json.loads(PRODUCT.read_text())


def test_depression_curve_against_its_cells():
    """The area rises linearly from the floor cell's height through each next cell's to all of them at the pass, and
    the volume at the pass is that curve's integral."""
    z = np.array([[10.0, 11.0, 13.0, np.nan], [5.0, np.nan, np.nan, np.nan]])
    a = np.array([[1.0, 2.0, 3.0, 0.0], [4.0, 0.0, 0.0, 0.0]])
    n_below = np.array([3, 1])
    pass_m = np.array([14.0, 7.0])
    vmax = np.array([1.0 * 1 / 2 + (1 + 3) / 2 * 2 + (3 + 6) / 2 * 1, 4.0 * 2 / 2])
    s, A, V, top = B.depression_curves(z, a, n_below, pass_m, vmax)
    assert np.allclose(s[0, :4], [0.0, 1.0, 3.0, 4.0]) and np.allclose(A[0, :4], [0.0, 1.0, 3.0, 6.0])
    assert np.allclose(V[0, 3], vmax[0]) and np.allclose(top, vmax)
    stage, area = B.on_curve(np.array([V[0, 2], 0.25]), s, A, V)
    assert stage[0] == pytest.approx(3.0) and area[0] == pytest.approx(3.0)
    assert stage[1] == pytest.approx(0.5) and area[1] == pytest.approx(1.0)        # a wedge: A = 2 stage, V = stage^2


def test_bowls_are_phase_2a_paraboloids():
    d = np.array([10.0, 100.0, 1000.0])
    s, A, V, top = B.bowl_curves(d)
    b = RT.bowl_classes('young')
    vol = math.pi / 8.0 * d ** 2 * RT.BOWL_DEPTH_RATIO * d
    assert np.allclose(top, vol)
    assert np.allclose(b['volume_m3'], math.pi / 8.0 * (b['diameter_km'] * 1000.0) ** 3 * RT.BOWL_DEPTH_RATIO)
    stage, area = B.on_curve(top / 4.0, s, A, V)                                    # a quarter of the volume: half the depth
    assert np.allclose(stage, 0.5 * RT.BOWL_DEPTH_RATIO * d) and np.allclose(area, 0.5 * math.pi * d ** 2 / 4.0)


def test_seal_passes_less_at_lunar_gravity_through_its_capillary_drive():
    """Conductivities fall to a sixth and suction heads rise sixfold: the water's weight drives a sixth as much, the
    capillary pull below the seal the same, so the Moon's seal passes between a sixth and all of Earth's."""
    d = np.array([0.1, 0.5, 1.0, 3.0])
    moon, earth = B.seal_seepage(d, IN.MOON), B.seal_seepage(d, 1.0)
    ratio = moon / earth
    assert np.all(ratio > IN.MOON) and np.all(ratio < 1.0)
    assert earth[1] * 8.64e7 == pytest.approx(5.6, rel=0.05)                         # the design's 5.5 mm a day at 0.5 m
    assert B.bank_seepage(IN.MOON) == pytest.approx(B.bank_seepage(1.0) * IN.MOON)


def test_charnock_roughness_carries_gravity():
    ce_m, z0_m = B.transfer(np.array([3.0]), IN.MOON)
    ce_e, z0_e = B.transfer(np.array([3.0]), 1.0)
    assert z0_m[0] > 2.0 * z0_e[0]
    assert ce_m[0] == pytest.approx(ce_e[0], rel=0.15)


def synthetic_forcing(cells=2):
    shape = (10, cells)
    sw = np.tile(np.array([0, 0, 60, 300, 450, 450, 300, 60, 0, 0], float)[:, None], (1, cells))
    return dict(sw_in=sw, lw_in=np.full(shape, 430.0), ta=np.full(shape, 300.0), qa=np.full(shape, 0.012),
                p=np.full(shape, 117000.0), wind=np.full(shape, 2.0), rls=np.full(shape, -15.0))


def test_energy_balance_is_penman_without_storage():
    """Under forcing held steady through the cycle the pool's balance settles where Penman's exact form has it."""
    f = synthetic_forcing()
    for k in f:
        f[k] = np.repeat(f[k][4:5], 10, axis=0)
    e, t = B.energy_balance(f, IN.MOON, depths=np.array([0.05, 1.0]), cycles=2)
    pen, tw = B.penman(f, IN.MOON)
    assert np.allclose(e[:, 0, :], pen.T, rtol=0.01) and np.allclose(e[:, 1, :], pen.T, rtol=0.02)


def test_deep_pools_shift_evaporation_into_the_night():
    f = synthetic_forcing()
    e, t = B.energy_balance(f, IN.MOON, depths=np.array([0.03, 3.0]))
    night = 1.0 - B.DAY_FRACTION
    assert (e[0, 1] * night).sum() > (e[0, 0] * night).sum()
    assert (e[0, 1] * B.DAY_FRACTION).sum() < (e[0, 0] * B.DAY_FRACTION).sum()


def synthetic_case(nm=240, inflow=2000.0, chain=True):
    """One cell, steady air: two bowls and a two-depression chain whose upper member overflows into the lower."""
    x = dict(tenth=(np.arange(nm) % 10).astype(np.int8)[:, None], pr=np.full((nm, 1), 2e-9, np.float32))
    off = np.full((nm, 1), 0.2)
    s_b, A_b, V_b, top_b = B.bowl_curves(np.array([30.0, 200.0]))
    z = np.array([[0.0, np.nan], [-50.0, -49.0]])
    a = np.array([[2000.0, 0.0], [1e6, 1e6]])
    s_d, A_d, V_d, top_d = B.depression_curves(z, a, np.array([1, 2]), np.array([0.5, 0.0]), np.array([500.0, 1.5e7]))
    s, A, V, top = B.pad_curves([(s_b, A_b, V_b, top_b), (s_d, A_d, V_d, top_d)], 3)
    pools = dict(cell=np.zeros(4, np.int64), s=s, A=A, V=V, top=top,
                 inflow=np.array([[inflow / 10] * 10, [inflow * 50 / 10] * 10, [inflow * 20 / 10] * 10, [0.0] * 10]),
                 parent_pos=np.array([-1, -1, 3 if chain else -1, -1]))
    tab = np.full((1, B.D_REF.size, 10), 0.08)                                       # mm/h at every depth and tenth
    return x, off, pools, tab


def test_pools_close_their_water_balance_and_spill_down_the_tree():
    x, off, pools, tab = synthetic_case()
    st = B.simulate(pools, x, tab, 'clay_floored', IN.MOON, off, np.array([4.0]))
    gain = st['inflow'] + st['rain']
    out = st['evap'].sum(axis=1) + st['seep'] + st['spill_out'] + (st['storage_end'] - st['storage_start'])
    spilled_in = out[3] - gain[3]
    assert spilled_in > 0                                                              # the lower depression is fed by spill
    assert (gain[:3].sum() - (out[:3].sum() + spilled_in)) == pytest.approx(0.0, abs=0.01 * gain.sum())
    assert np.all(st['stage_max'][:2] <= RT.BOWL_DEPTH_RATIO * np.array([30.0, 200.0]) + 1e-9)


def test_a_steadily_fed_pool_settles_where_its_area_balances_the_inflow():
    """dV/dt = G - L A: whatever the shape, the area settles at G / L."""
    x, off, pools, tab = synthetic_case(nm=600, inflow=200.0, chain=False)
    st = B.simulate(pools, x, tab, 'clay_floored', IN.MOON, off, np.array([4.0]))
    g = pools['inflow'][0].sum() / (B.JULIAN_YEAR_DAYS * 86400.0)                                         # m3/s
    area = st['area'][0].sum() / max(st['wet'][0], 1)
    depth = st['storage_end'][0] / area
    loss = 0.08 / 1000.0 / 3600.0 + float(B.seal_seepage(np.array([depth]), IN.MOON)[0])
    rain = 2e-9 * (1.0 - 0.2 / (2e-9 * R.MM_PER_MEAN))
    assert area == pytest.approx(g / (loss - rain), rel=0.03)


def test_persistence_classes():
    st = dict(wet=np.array([0, 5, 5, 5, 9]), dry_ever=np.array([True, True, True, True, False]),
              checks=dict(sunset=np.array([0, 1, 4, 4, 4]), sunrise=np.array([0, 0, 3, 4, 4]), before_storms=np.array([0, 0, 1, 3, 4])),
              n_at=dict(sunset=np.full(5, 4), sunrise=np.full(5, 4), before_storms=np.full(5, 4)))
    cls, share = B.persistence(st)
    assert list(cls) == [0, 1, 2, 3, 4]


def test_lunar_rills_carry_less_shear():
    """At the same discharge a lunar sheet stands 2.45^0.6 = 1.71 times as deep under a sixth of the gravity, so its
    shear is 0.28 of Earth's; the drops detach at the energy factor."""
    c = dict(runoff_mm=10.0, duration_h=0.5, capacity_mm_h=5.0)
    tex = B.TEXTURES['loam']
    im, rm = B.event_wear(c, tex, False, IN.MOON, 0.143)
    ie, re = B.event_wear(c, tex, False, 1.0, 1.0)
    assert im / ie == pytest.approx(0.143, rel=0.02)
    q = 10.0 / 1000.0 / 1800.0 * 100.0
    h_ratio = (1.0 / math.sqrt(IN.MOON)) ** 0.6
    assert h_ratio * IN.MOON == pytest.approx(0.283, abs=0.003)
    assert rm < re


def test_wepp_baselines_as_documented():
    ki, kr, tc = B.wepp(B.TEXTURES['loam'])
    assert ki == pytest.approx(2.728e6 + 1.921e5 * 10.0)
    assert kr == pytest.approx(0.00197 + 0.0003 * 10.0 + 0.03863 * math.exp(-1.84))
    assert tc == pytest.approx(2.67 + 0.065 * 20.0 - 0.058 * 10.0)


def test_lake_mead_check_is_recorded_with_its_rule():
    c = B.mead_check()
    assert c['rule']['tolerance'] == 0.25 and c['measured_w_m2'] == 164.0
    assert c['passes'] == (abs(c['ratio'] - 1.0) <= 0.25)
    assert len(B.MEAD['le']) == len(B.MEAD['ta']) == 24
    assert np.mean(B.MEAD['le']) == pytest.approx(164.0, abs=2.0)          # Table 4's monthly values give its period mean


def test_product_phase2b(product):
    p2 = product['phase2b']
    assert set(p2['pools']) == set(R.BRACKETS)
    for b in R.BRACKETS:
        for n in IN.MAIN:
            for seep in B.SEEPAGE:
                for s in B.STATES:
                    rows = p2['pools'][b][n][seep][s]
                    for kind in ('depressions', 'bowls_wall_fed', 'bowls_channel_fed'):
                        assert sum(rows[kind]['share_in_class'].values()) == pytest.approx(1.0, abs=1e-6)
                    lev = p2['leverage'][b][n][seep][s]['drylands']
                    assert lev['passes_yearly'] == (abs(lev['change_share_of_evaporation_yearly']) > R.THRESHOLDS['yearly'])
    assert p2['checks']['lake_mead']['rule']['tolerance'] == 0.25
    assert p2['checks']['night']['rule']['night_mm'] == [2.5, 30.0]
    row = product['cells'][0]['pools_2b_walnut_gulch_clay_floored']
    assert set(row) == set(IN.MAIN)


def test_wet_patch_layer_follows_the_convention(product):
    if not B.WET_PATCHES.is_file():
        pytest.skip('run python -m research.studies.storm_water.run')
    layer = json.loads(B.WET_PATCHES.read_text())
    assert layer['schema'] == B.WET_SCHEMA
    for key in ('producer', 'evidence', 'reading_rule', 'variants', 'cells', 'layers'):
        assert key in layer
    assert layer['producer']['files'] == product['producer']['files']
    assert len(layer['variants']) == len(R.BRACKETS) * len(IN.MAIN) * len(B.STATES) * len(B.SEEPAGE)
    v = layer['layers'][layer['variants'][0]]
    assert len(v['wet_ppm_by_tenth']) == len(layer['cells']) and len(v['wet_ppm_by_tenth'][0]) == 10
    if B.WET_CELLS.is_file():
        assert R.sha256(B.WET_CELLS) == layer['cells_1_9km']['sha256']
