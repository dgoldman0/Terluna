"""Routing, fill-and-spill and lake balance on small synthetic terrains."""
import numpy as np

from geography import drainage as dr

ROWS, COLS = 24, 48
LAT = np.linspace(2.0, -2.0, ROWS)          # near the equator, so cells are about square
CELL_KM = 1.0


def slope_to_sea(bowl_depth=0.0):
    """Land rising away from a sea band in columns 0-5 (longitude wraps), with an optional bowl on the slope."""
    col = np.arange(COLS)
    distance = np.minimum(np.abs(col - 2.5), COLS - np.abs(col - 2.5)) - 3.0
    h = np.repeat((np.maximum(distance, 0.0) * 10.0)[None, :], ROWS, axis=0) + 0.01 * np.arange(ROWS)[:, None]
    h[:, :6] = -50.0
    if bowl_depth:
        r, c = np.meshgrid(np.arange(ROWS), col, indexing='ij')
        h -= bowl_depth * np.exp(-((r - 12) ** 2 + (c - 20) ** 2) / 8.0)
    return h, h < -1.0


def total(r, weights):
    return float((weights.ravel() * r['area']).sum())


def test_all_runoff_reaches_the_sea():
    h, outlet = slope_to_sea()
    runoff = np.ones_like(h)                       # 1 mm/day everywhere on land
    r = dr.route(h, outlet, LAT, CELL_KM, runoff, np.zeros_like(h))
    generated = total(r, np.where(outlet, 0.0, runoff) * dr.MM_DAY)
    assert not r['lake'].any()
    assert np.isclose(r['discharge'][r['outlet']].sum(), generated, rtol=1e-9)


def test_a_wet_depression_fills_and_spills():
    h, outlet = slope_to_sea(bowl_depth=80.0)
    runoff = np.ones_like(h)
    lake_net = np.full_like(h, 0.5)                # rain on the lake exceeds its evaporation
    r = dr.route(h, outlet, LAT, CELL_KM, runoff, lake_net)
    lake = r['lake'].reshape(h.shape)
    assert lake[12, 20] and lake.sum() > 5
    k = r['label'][12 * COLS + 20]
    assert not r['closed'][k]
    rim = r['level'][k]
    assert np.all(h[lake] < rim) and np.isclose(r['lake_level'][12 * COLS + 20], rim)
    land = ~outlet & ~lake
    expected = total(r, np.where(land, runoff, 0.0) * dr.MM_DAY) + total(r, np.where(lake, lake_net, 0.0) * dr.MM_DAY)
    assert np.isclose(r['discharge'][r['outlet']].sum(), expected, rtol=1e-9)


def test_a_dry_depression_keeps_a_smaller_closed_lake():
    h, outlet = slope_to_sea(bowl_depth=80.0)
    runoff = np.ones_like(h)
    lake_net = np.full_like(h, -20.0)              # strong evaporation from open water
    r = dr.route(h, outlet, LAT, CELL_KM, runoff, lake_net)
    k = r['label'][12 * COLS + 20]
    assert r['closed'][k] and r['passed_on'][k] == 0.0
    full = ~r['outlet'] & (h.ravel() < r['level'][r['label']]) & (r['label'] == k)
    kept = r['lake'] & (r['label'] == k)
    assert 0 < kept.sum() < full.sum()
    evaporation = float((20.0 * dr.MM_DAY * r['area'])[kept].sum())
    assert evaporation <= r['inflow'][k] + 1e-12
    assert evaporation + float((20.0 * dr.MM_DAY * r['area'])[full & ~kept].min()) > r['inflow'][k]
