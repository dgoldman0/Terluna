"""Tests of the storm-water study's tables for PlaSim's hollows' store: the store's curve sums phase 2b's hollows at a
common fill, the floors keep phase 2b's seepage laws, the routing shares split each cell's runoff whole, the runoff
tables are phase 1's Green and Ampt steps, and the product carries the grids and fields climate/gcm/land_store.py
reads, bound to its code and inputs."""
from __future__ import annotations

import json
import math

import numpy as np
import pytest

from research.studies.storm_water import infiltration as IN
from research.studies.storm_water import phase2b as B
from research.studies.storm_water import run as R
from research.studies.storm_water import store_tables as S


@pytest.fixture(scope='module')
def product():
    if not S.OUT.is_file():
        pytest.skip('run python -m research.studies.storm_water.store_tables')
    return json.loads(S.OUT.read_text())


def two_cell_hollows():
    """One cell with a two-cell depression and a bowl class of 50 bowls of 100 m; a second with one wedge."""
    z = np.array([[10.0, 11.0, np.nan], [5.0, np.nan, np.nan]])
    a = np.array([[1.0e6, 2.0e6, 0.0], [4.0e6, 0.0, 0.0]])
    n_below, pass_m = np.array([2, 1]), np.array([13.0, 7.0])
    vmax = np.array([1.0e6 * 1 / 2 + 3.0e6 / 2 * 2, 4.0e6 * 2 / 2])
    s, A, V, top = B.depression_curves(z, a, n_below, pass_m, vmax)
    sb, Ab, Vb, topb = B.bowl_curves(np.array([100.0]))
    s, A, V, cap = B.pad_curves([(s, A, V, top), (sb, Ab, Vb, topb)], s.shape[1])
    return dict(s=s, A=A, V=V, cap=cap, number=np.array([1.0, 1.0, 50.0]), cell=np.array([0, 1, 0]),
                kind=np.array([0, 0, 1], np.int8), wide=np.zeros(3, bool))


def test_the_store_sums_its_hollows_at_a_common_fill():
    h = two_cell_hollows()
    dry = np.array([1.0e8, 1.0e8])
    cap, area, seep = S.store_curve(h, 2, dry)
    assert cap == pytest.approx([(h['cap'][0] + 50.0 * h['cap'][2]) / dry[0], h['cap'][1] / dry[1]])
    rim = math.pi * 100.0 ** 2 / 4.0
    for k, f in enumerate(S.FILLS):
        v = f * h['cap']
        _, a = S.geometry(v, h['s'], h['A'], h['V'])
        assert a[2] == pytest.approx(rim * math.sqrt(f), rel=1e-12)              # the bowl, a paraboloid
        assert area[0, k] == pytest.approx((a[0] + 50.0 * a[2]) / dry[0], rel=1e-12, abs=1e-30)
        assert area[1, k] == pytest.approx(a[1] / dry[1], rel=1e-12, abs=1e-30)
    assert area[:, 0] == pytest.approx([0.0, 0.0]) and seep['clay_floored'][:, 0] == pytest.approx([0.0, 0.0])
    # full, the depressions cover their cells below the pass
    assert area[1, -1] == pytest.approx(4.0e6 / dry[1]) and area[0, -1] == pytest.approx((3.0e6 + 50.0 * rim) / dry[0])


def test_the_floors_keep_phase_2bs_seepage_laws():
    d = np.array([0.05, 0.5, 2.0, 5.0])
    a = np.array([1.0e4, 2.0e4, 3.0e4, 4.0e4])
    clog = np.array([2.0e4, 2.0e4, 1.0e4, 4.0e4])
    q = np.interp(d, S.DEPTH_TAB, B.seal_seepage(S.DEPTH_TAB, IN.MOON))
    assert S.seepage(d, a, clog, 'clay_floored') == pytest.approx(q * a)
    banks = S.seepage(d, a, clog, 'sandy_banked')
    assert banks[[0, 1, 3]] == pytest.approx((q * a)[[0, 1, 3]])                  # at or below the clogged stage
    assert banks[2] == pytest.approx(q[2] * 1.0e4 + B.bank_seepage(IN.MOON) * 2.0e4)
    tight = S.seepage(d, a, clog, 'tight_seal')
    assert tight == pytest.approx(np.interp(d, S.DEPTH_TAB, B.seal_seepage(S.DEPTH_TAB, IN.MOON, seal=S.TIGHT)) * a)
    assert (tight < q * a).all() and S.TIGHT['k_earth_mm_h'] == pytest.approx(B.SEAL['k_earth_mm_h'] / 10.0)


def test_the_routing_shares_split_the_runoff_whole():
    rng = np.random.default_rng(1)
    ncell = 5
    x = dict(pr=np.zeros((240, ncell)), tenth=np.tile(np.arange(10), 24)[:, None] * np.ones((1, ncell), np.int8),
             area_km2=np.ones(ncell))
    inp = {}
    for b in S.LIBRARIES:
        for n in S.SURFACES:
            off = rng.uniform(0.0, 2.0, (240, ncell))
            inp[f'off__{b}__{n}'] = off
            rt = R.by_tenth(off, x['tenth']) / 2.0 * (365.25 / (12 * 29.530588853))
            for s in S.STATES:
                w = dict(bowls_inside=0.3, bowls_slope=0.1, bowls_channel=0.2, depressions=0.1, slope_recharge=0.0,
                         bed_recharge=0.01, sea=0.002, humid_land=0.001)
                for q, v in w.items():
                    share = np.full(ncell, v)
                    if q == 'depressions':
                        share[0] = 0.6                                 # a cell receiving its neighbour's flow
                    inp[f'{q}__{b}__{n}__{s}'] = rt * share[None, :]
    shares, checks = S.routing_shares(inp, x)
    for v in shares.values():
        assert v.sum(axis=1) == pytest.approx(np.ones(ncell)) and (v >= 0.0).all()
        assert v[1] == pytest.approx([0.7, 0.003, 0.01, 0.287])
        assert v[0, 3] == 0.0 and v[0, :3] == pytest.approx(np.array([1.2, 0.003, 0.01]) / 1.213)
    assert all(c['cells_scaled'] == 1 for c in checks.values())


def test_the_runoff_tables_are_phase_1s_steps():
    t = S.runoff_tables(116000.0)
    for name, tab in t.items():
        assert tab.shape == (S.DEPTHS_MM.size, S.WETNESS.size)
        assert (tab[0] == 0.0).all() and (tab >= 0.0).all() and (tab <= 1.0).all()
        assert (np.diff(tab, axis=0) >= -1e-12).all()                          # deeper steps shed a larger share
        assert (np.diff(tab, axis=1) >= -1e-12).all()                          # so does wetter soil
        p = IN.parameters(name)
        s = IN.initial_saturation(p['soil'], 0.04 * IN.WSMAX_M)
        one = IN.steady_rain(p, np.array([20.0]), np.array([40.0]), IN.deficit(p, np.array([s])), 0.0,
                             dt_h=S.INFILTRATION_DT_H)[0] / 20.0
        assert tab[list(S.DEPTHS_MM).index(20.0), list(S.WETNESS).index(0.04)] == pytest.approx(one, rel=1e-12)
    # the sealed loam sheds even a gentle step; the sandy loam takes a step of a few millimetres
    k = list(S.DEPTHS_MM).index(5.0)
    assert t['sealed_loam'][k, 1] > 0.3 and t['sandy_loam'][k, 1] == 0.0


def test_product_store_tables(product):
    assert product['schema'] == S.SCHEMA and product['grid'] == dict(nlon=64, nlat=32)
    g = product['grids']
    assert (len(g['fill']), len(g['rain_depth_mm']), len(g['wetness'])) == (12, 19, 11)
    assert g['fill'][0] == 0.0 and g['fill'][-1] == 1.0 and g['rain_depth_mm'][0] == 0.0
    assert product['floors'] == list(S.FLOORS) and product['surfaces'] == list(S.SURFACES)
    assert set(product['runoff_tables']) == set(S.SURFACES)
    cells = product['cells']
    assert len(cells) == product['basis']['cells'] and len({(c['j'], c['i']) for c in cells}) == len(cells)
    for c in cells[::37]:
        for s in S.STATES:
            st = c['store'][s]
            assert st['capacity_m'] > 0.0 and st['area_share'][0] == 0.0
            assert (np.diff(st['area_share']) >= -1e-12).all() and st['area_share'][-1] < 1.0
            for f in S.FLOORS:
                assert len(st['seepage_m_s'][f]) == 12 and st['seepage_m_s'][f][0] == 0.0
        for b in S.LIBRARIES:
            for n in S.SURFACES:
                for s in S.STATES:
                    assert sum(c['shares'][b][n][s].values()) == pytest.approx(1.0)
    assert set(product['producer']['files']) == set(S.CODE)
    for path, digest in product['producer']['files'].items():
        assert R.sha256(path) == digest, path
    assert 'reading_rule' in product and 'evidence' in product
