"""Tests of the storm-water study's tables for PlaSim's hollows' store (schema /2): the store's curve shares a cell's
water among its hollows by their own inflow, each at most full; the floors keep phase 2b's seepage laws; the balanced
volume meets the store's mean inflow with its mean losses; the routing shares split each cell's runoff whole; the
convective runoff tables are phase 2a's Poisson mixing of its wetting tables per millimetre, exact in the bucket's
wetness between the grid's points; the large-scale tables are phase 1's steady 2 mm/h rain; and the product carries
the grids and fields climate/gcm/land_store.py reads, balanced by its own fields and bound to its code and inputs."""
from __future__ import annotations

import json
import math

import numpy as np
import pytest

from research.studies.storm_water import infiltration as IN
from research.studies.storm_water import phase2a as P2
from research.studies.storm_water import phase2b as B
from research.studies.storm_water import run as R
from research.studies.storm_water import store_tables as S


@pytest.fixture(scope='module')
def product():
    if not S.OUT.is_file():
        pytest.skip('run python -m research.studies.storm_water.store_tables')
    return json.loads(S.OUT.read_text())


def two_cell_hollows():
    """One cell with a two-cell depression and 50 bowls of 100 m; a second with one wedge; yearly inflows."""
    z = np.array([[10.0, 11.0, np.nan], [5.0, np.nan, np.nan]])
    a = np.array([[1.0e6, 2.0e6, 0.0], [4.0e6, 0.0, 0.0]])
    n_below, pass_m = np.array([2, 1]), np.array([13.0, 7.0])
    vmax = np.array([1.0e6 * 1 / 2 + 3.0e6 / 2 * 2, 4.0e6 * 2 / 2])
    s, A, V, top = B.depression_curves(z, a, n_below, pass_m, vmax)
    sb, Ab, Vb, topb = B.bowl_curves(np.array([100.0]))
    s, A, V, cap = B.pad_curves([(s, A, V, top), (sb, Ab, Vb, topb)], s.shape[1])
    return dict(s=s, A=A, V=V, top=cap, weight=np.array([1.0, 1.0, 50.0]), cell=np.array([0, 1, 0]),
                kind=np.array([0, 0, 1], np.int8), q=np.array([2.0e5, 1.0e5, 600.0]))


def test_the_store_shares_its_water_among_its_hollows_by_their_inflow():
    h = two_cell_hollows()
    dry = np.array([1.0e8, 1.0e8])
    cap, area, seep, held = S.store_curve(h, 2, dry)
    c0 = h['top'][0] + 50.0 * h['top'][2]
    assert cap == pytest.approx([c0 / dry[0], h['top'][1] / dry[1]])
    q0 = h['q'][0] + 50.0 * h['q'][2]
    rim = math.pi * 100.0 ** 2 / 4.0
    for k, f in enumerate(S.FILLS):
        v = np.minimum(f * np.array([c0 * h['q'][0] / q0, h['top'][1], c0 * h['q'][2] / q0]), h['top'])
        _, a = S.geometry(v, h['s'], h['A'], h['V'])
        assert area[0, k] == pytest.approx((a[0] + 50.0 * a[2]) / dry[0], rel=1e-12, abs=1e-30)
        assert area[1, k] == pytest.approx(a[1] / dry[1], rel=1e-12, abs=1e-30)
        assert held[0, k] == pytest.approx((v[0] + 50.0 * v[2]) / (f * c0), rel=1e-12) if f > 0 else True
    # the bowls, a paraboloid, take their inflow's share until full and then hold their capacity
    f_full = h['top'][2] / (c0 * h['q'][2] / q0)
    _, a = S.geometry(np.array([0.0, 0.0, 0.25 * f_full * c0 * h['q'][2] / q0]), h['s'], h['A'], h['V'])
    assert a[2] == pytest.approx(rim * math.sqrt(0.25 * f_full * c0 * h['q'][2] / q0 / h['top'][2]), rel=1e-12)
    assert area[:, 0] == pytest.approx([0.0, 0.0]) and seep['clay_floored'][:, 0] == pytest.approx([0.0, 0.0])
    assert area[1, -1] == pytest.approx(4.0e6 / dry[1])                     # a lone hollow fills whole


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


def test_the_balanced_volume_meets_the_mean_inflow_with_the_mean_losses():
    # a paraboloid's store (area as the root of the fill), seepage in proportion to the area, a steady pool rate:
    # inflow (1 - a) + a (rain - e - k) = 0 at a* = inflow / (inflow + e + k - rain), the fill (a* / a0)^2
    fill = S.FILLS
    a0, k, e, rain = 0.2, 3.0e-8, 2.0e-8, 1.0e-9
    area = np.array([a0 * np.sqrt(fill)] * 3)
    seep = k * area
    cap = np.array([5.0, 5.0, 5.0])
    inflow = np.array([1.0e-10, 1.0e-6, 0.0])
    e_day = np.full((3, B.D_REF.size), e)
    v = S.balanced_volume(fill, area, seep, cap, inflow, rain, e_day)
    a_star = inflow[0] / (inflow[0] + e + k - rain)
    assert v[0] == pytest.approx((a_star / a0) ** 2 * cap[0], rel=1e-9)
    assert abs(S.balance_residual(v[:1], fill, area[:1], seep[:1], cap[:1], inflow[:1], rain, e_day[:1])[0]) < 1e-9 * inflow[0]
    assert v[1] == cap[1]                                                  # gains even when full
    assert v[2] == 0.0                                                     # no runoff reaches it


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


def test_the_runoff_per_millimetre_is_phase_2as_poisson_mixing_and_exact_in_the_wetness():
    rng = np.random.default_rng(4)
    tot = np.vstack([np.zeros((1, P2.SAT_W.size)), np.cumsum(rng.uniform(0.0, 2.0, (6, P2.SAT_W.size)), axis=0)])
    tot = np.sort(tot, axis=1)                                             # wetter soil sheds more
    mu = 6.0
    s = np.array([0.12, 0.33])
    lam = np.array([2.0, 2.0])
    by_hand = sum(math.exp(-2.0) * 2.0 ** n / math.factorial(n) * np.interp(s, P2.SAT_W, tot[min(n, 6)]
                                                                           + max(n - 6, 0) * (tot[6] - tot[5]))
                  for n in range(60))
    assert S.per_mm(tot, mu, lam, s) == pytest.approx(by_hand / (2.0 * mu), rel=1e-12)
    assert S.per_mm(tot, mu, np.zeros(2), s) == pytest.approx(np.interp(s, P2.SAT_W, tot[1]) / mu, rel=1e-12)
    assert S.per_mm(tot, mu, np.array([1e-9, 1e-9]), s) == pytest.approx(np.interp(s, P2.SAT_W, tot[1]) / mu, rel=1e-6)
    # the wetness grid holds every surface's soil's breakpoints, so the tables are linear in the wetness between points
    for soil in {IN.SURFACES[n]['soil'] for n in S.SURFACES}:
        s0, s1 = IN.initial_saturation(IN.SOILS[soil], np.array([0.0, IN.WSMAX_M]))
        for sk in P2.SAT_W[(P2.SAT_W > s0) & (P2.SAT_W < s1)]:
            assert np.min(np.abs(S.WETNESS - (sk - s0) / (s1 - s0))) < 1e-11
        mids = 0.5 * (S.WETNESS[1:] + S.WETNESS[:-1])
        sat = lambda w: IN.initial_saturation(IN.SOILS[soil], w * IN.WSMAX_M)
        lam = np.full(mids.size, 1.7)
        exact = S.per_mm(tot, mu, lam, sat(mids))
        ends = S.per_mm(tot, mu, np.full(S.WETNESS.size, 1.7), sat(S.WETNESS))
        assert exact == pytest.approx(0.5 * (ends[1:] + ends[:-1]), rel=1e-10)


def test_the_large_scale_tables_are_phase_1s_steady_rain():
    t = S.large_scale_tables(116000.0)
    for name, tab in t.items():
        assert tab.shape == (S.RAIN_MM_DAY.size, S.WETNESS.size)
        assert (tab[0] == 0.0).all() and (tab >= 0.0).all() and (tab <= 1.0).all()
        assert (np.diff(tab, axis=1) >= -1e-12).all()                          # wetter soil sheds a larger share
        p = IN.parameters(name, pressure_pa=116000.0)
        s = IN.initial_saturation(p['soil'], 0.04 * IN.WSMAX_M)
        k = list(S.RAIN_MM_DAY).index(20.0)
        depth = 20.0 * R.MEAN_DAYS
        one = (IN.steady_rain(p, np.array([depth]), R.LS_RATE_MM_H, IN.deficit(p, np.array([s])), np.zeros(1))[0] / depth
               if R.LS_RATE_MM_H > p['k_sub'] else 0.0)
        assert tab[k, list(S.WETNESS).index(0.04)] == pytest.approx(one, rel=1e-12, abs=1e-15)


def test_product_store_tables(product):
    assert product['schema'] == S.SCHEMA and product['grid'] == dict(nlon=64, nlat=32)
    g = product['grids']
    assert (len(g['fill']), len(g['rain_mm_day']), len(g['wetness'])) == (12, 29, 13)
    assert g['fill'][0] == 0.0 and g['fill'][-1] == 1.0 and g['rain_mm_day'][0] == 0.0 and g['wetness'][-1] == 1.0
    assert product['floors'] == list(S.FLOORS) and product['surfaces'] == list(S.SURFACES)
    for b in S.LIBRARIES:
        for n in S.SURFACES:
            tab = np.array(product['convective_runoff'][b][n])
            assert tab.shape == (29, 13) and (tab > 0.0).all() and (tab <= 1.0).all()
            assert (np.diff(tab, axis=1) >= -1e-9).all()                       # wetter soil sheds a larger share
            r = product['checks']['reproduction_of_phase2a_runoff'][f'{b}__{n}']
            assert abs(r['direct_over_phase2a'] - 1.0) < 1e-6 and abs(r['tables_over_phase2a'] - 1.0) < 0.005
    fill = np.array(g['fill'])
    cells = product['cells']
    assert len(cells) == product['basis']['cells'] and len({(c['j'], c['i']) for c in cells}) == len(cells)
    for c in cells[::37]:
        e_day = np.array([c['pool_evaporation_m_s']])
        for b in S.LIBRARIES:
            for n in S.SURFACES:
                for s in S.STATES:
                    st = c['store'][b][n][s]
                    assert st['area_share'][0] == 0.0 and (np.diff(st['area_share']) >= -1e-12).all()
                    assert st['area_share'][-1] < 1.0 and st['capacity_m'] >= 0.0
                    assert sum(c['shares'][b][n][s].values()) == pytest.approx(1.0)
                    for f in S.FLOORS:
                        q = st['seepage_m_s'][f]
                        assert len(q) == 12 and q[0] == 0.0
                        vb = st['balanced_m'][f]
                        assert 0.0 <= vb <= st['capacity_m']
                        if 0.0 < vb < st['capacity_m']:              # the balance holds where it is not at an end
                            res = S.balance_residual(np.array([vb]), fill, np.array([st['area_share']]), np.array([q]),
                                                     np.array([st['capacity_m']]), np.array([st['inflow_m_s']]),
                                                     c['rain_m_s']['total'], e_day)[0]
                            assert abs(res) <= 1e-5 * st['inflow_m_s']
    assert set(product['producer']['files']) == set(S.CODE)
    for path, digest in product['producer']['files'].items():
        assert R.sha256(path) == digest, path
    assert 'reading_rule' in product and 'evidence' in product
