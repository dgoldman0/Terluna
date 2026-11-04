"""Checks of the tide study: the cross-shore model against Ferris, Nielsen and the coast-water study's steady belt, its
water balance, the belts with a reach for each pixel against run.belts_for, the lakes' tides, the evaporation path,
and the product's binding to its code and inputs."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

import numpy as np
import pytest

from research.studies.coast_water import belts as B
from research.studies.coast_water import load as L
from research.studies.coast_water import tide as TD

ROOT = Path(__file__).resolve().parents[3]
PRODUCT = Path(__file__).resolve().parent / 'results' / 'coast_tide.json'
T_CENTRAL, SY = 5.17172e-4, math.sqrt(0.05 * 0.10)
DRAW = 1.92e-3 / 86400.0


def test_ferris_damping():
    k = TD.damping(TD.MONTH_S, SY, T_CENTRAL)
    w = 2.0 * math.pi / TD.MONTH_S
    assert k == pytest.approx(math.sqrt(w * SY / (2.0 * T_CENTRAL)), rel=1e-12)
    amp, lag = TD.ferris(1.0 / k, 2.0, TD.MONTH_S, SY, T_CENTRAL)
    assert amp == pytest.approx(2.0 * math.exp(-1.0), rel=1e-12)
    assert lag == pytest.approx(1.0 / w, rel=1e-12)
    # Earth's semidiurnal tide in ground six times as transmissive damps sqrt(53/6) times faster than the month here
    k_earth = TD.damping(12.42 * 3600.0, SY, 6.0 * T_CENTRAL)
    assert k_earth / k == pytest.approx(math.sqrt(TD.MONTH_S / (12.42 * 3600.0) / 6.0), rel=1e-12)


def test_model_follows_ferris_on_a_steep_face():
    k = TD.damping(TD.MONTH_S, SY, T_CENTRAL)
    x, ground = TD.planar_grid(1.0, 0.1, k, 6.0 / k)
    at = int(np.searchsorted(x, 1.0 / k))
    level = TD.sinusoid(0.1)
    r = TD.cross_shore(x, ground, SY, T_CENTRAL, level, dt=TD.MONTH_S / 480, window=TD.MONTH_S, tolerance=1e-6,
                       max_windows=8, record=np.array([at]))
    share = 0.5 * (r['high'][at] - r['low'][at]) / 0.1
    assert share == pytest.approx(math.exp(-k * x[at]), rel=0.03)
    lag = TD.lag(r['series'][:, 0], r['times'], TD.MONTH_S, level)
    assert lag == pytest.approx(k * x[at] * TD.MONTH_S / (2.0 * math.pi), rel=0.05)


def test_model_follows_nielsen_at_small_eps_and_saturates():
    small = TD.overheight_run(0.05)
    assert small['overheight_over_amplitude'] == pytest.approx(0.025, rel=0.1)
    big = TD.overheight_run(4.0)
    assert 0.5 < big['overheight_over_amplitude'] < 1.0 < big['nielsen']


def test_steady_belt_is_the_coast_water_belt():
    for slope in (0.02, 0.1):
        r = TD.belt_run(slope, T_CENTRAL, SY, DRAW, 15.0, 1.0, steps=120, tolerance=3e-3)
        assert r['without']['edge_mean_m'] == pytest.approx(r['analytic_without_m'], rel=0.01)
        assert abs(r['without']['effective_rise_m']) < 0.05
        assert r['with']['edge_mean_m'] > r['without']['edge_mean_m']


def test_water_balance_closes():
    k = TD.damping(TD.MONTH_S, SY, T_CENTRAL)
    x, ground = TD.planar_grid(0.03, 1.0, k, 900.0)
    c = DRAW / (2 * T_CENTRAL)
    w0 = float(B.width(0.03, c, 15.0))
    start = -c * (2 * w0 * np.clip(x, 0, w0) - np.clip(x, 0, w0) ** 2)
    r = TD.cross_shore(x, ground, SY, T_CENTRAL, TD.sinusoid(1.0), draw=DRAW, reach=15.0, dt=TD.MONTH_S / 240,
                       start=np.minimum(start, np.maximum(ground, 0)), tolerance=1e-4, max_windows=400)
    drawn = r['draw'].sum() * r['dx'] * TD.MONTH_S
    assert r['inflow'] - r['outflow'] == pytest.approx(drawn, rel=0.03)


def synthetic_pixels(n=40, seed=1):
    rng = np.random.default_rng(seed)
    coast = dict(cell=rng.integers(0, 2048, n), area=np.full(n, 3.0), edge=rng.uniform(0, 1.9, (n, 4)),
                 across=np.full((n, 4), 1.9), rise=rng.uniform(0, 80, (n, 4)), island=rng.random(n) < 0.3)
    coast['edge'][:, 2:] = 0.0
    m = 3 * n
    inland = dict(cell=rng.integers(0, 2048, m), area=np.full(m, 3.0), distance=rng.uniform(1.9, 12, m),
                  step_rows=rng.integers(-6, 7, m), step_cols=rng.integers(-6, 7, m), height=rng.uniform(0, 300, m),
                  water=rng.integers(1, 5, m), island=rng.random(m) < 0.3)
    return dict(coast=coast, inland=inland)


def test_belts_with_one_reach_match_run_belts_for():
    from research.studies.coast_water import run as R
    from research.studies.coast_water import tide_run as TR
    px = synthetic_pixels()
    uptake = np.full((32, 64), DRAW)
    a = R.belts_for(px, uptake, T_CENTRAL, 15.0)
    b = TR.belts_reach(px, uptake, T_CENTRAL, np.full((40, 4), 15.0), np.full(120, 15.0))
    assert np.allclose(a['belt'], b['belt']) and np.allclose(a['length'], b['length'])
    c = TR.belts_reach(px, uptake, T_CENTRAL, np.full((40, 4), 16.0), np.full(120, 16.0))
    assert c['belt'].sum() > a['belt'].sum()


def test_a_lake_at_one_point_has_no_tide():
    from research.studies.coast_water import tide_run as TR
    forcing = dict(amp=np.ones((12, 20)) * (1 + 1j), kappa=2.84)
    lat, lon = np.array([40.0]), np.array([20.0])
    means = TR.monomials(TR.unit(lat, lon))
    assert TR.lake_amplitude(forcing, means, np.array([0]), lat, lon)[0] == pytest.approx(0.0, abs=1e-12)


def test_evaporation_path_of_the_held_and_left_alone_seas():
    chemistry = json.loads((ROOT / 'research/studies/sea_chemistry/results/sea_chemistry.json').read_text())
    seas = L.seas(chemistry)
    held = L.evaporation_path(seas['held'], 0.01)
    assert held['order'] == ['calcite', 'opal']
    assert held['share_laid_by_end']['calcite'] == pytest.approx(1.0, rel=1e-3)
    left = L.evaporation_path(seas['left_alone'], L.SOIL_CO2_ATM[0])
    assert left['order'][0] == 'calcite' and 'magnesium_clays' in left['order']
    assert seas['left_alone']['ta'] > 2 * seas['left_alone']['ca']          # the soda side of the divide
    d = L.deposits(seas['held'], 0.7)
    assert d['mol_m2_yr']['calcite'] == pytest.approx(700.0 * seas['held']['ca'], rel=1e-12)
    assert d['hardpan_cemented_mm_per_century'] == pytest.approx(d['hardpan_solid_mm_per_century'] / 0.3)


def sha(path):
    return hashlib.sha256((ROOT / path).read_bytes()).hexdigest()


@pytest.mark.skipif(not PRODUCT.exists(), reason='the product has not been written')
def test_product_is_bound_to_its_code_and_inputs():
    p = json.loads(PRODUCT.read_text())
    assert p['schema'] == 'terluna.research.coast-tide/1'
    for path, h in p['producer']['files'].items():
        assert sha(path) == h, path
    for path, h in p['producer']['inputs'].items():
        assert sha(path) == h, path
    record = json.loads((ROOT / 'research/studies/coast_water/tide_inputs.json').read_text())
    assert p['producer']['terrain_slices'] == {e['name']: e.get('selected_sha256') for e in record['slices']}


@pytest.mark.skipif(not PRODUCT.exists(), reason='the product has not been written')
def test_product_figures_hold_together():
    p = json.loads(PRODUCT.read_text())
    for case in p['numerics']['rise'].values():
        for r in case['rows']:
            assert r['belt_without_m'] == pytest.approx(r['analytic_m'], rel=0.03)
    small = [r for r in p['numerics']['no_draw'] if r['eps'] <= 0.1]
    assert all(r['overheight_over_amplitude'] == pytest.approx(r['nielsen'], rel=0.15) for r in small)
    c = p['belts']['cases']['central_15m']['rain_mask.every_cell']
    assert c['belt_km2'][1] >= c['belt_km2'][0]
    cw = json.loads((ROOT / 'research/studies/coast_water/results/coast_water.json').read_text())
    assert c['belt_km2'][0] == pytest.approx(cw['belts']['deserts']['rain_mask.every_cell']['central_15m']['belt_km2'],
                                             rel=1e-4)
