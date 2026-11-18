"""Checks of the moisture map on the design climate: a map pixel is sea where most of it lies under one sea, a year's
sums move by the estimate's step as pattern scaling's kinds say, and the record holds together once computed (the run, its years and its settled frame; the classes' shares adding up; no lake
at the common level by rule; the design surface's land and water accounted for on the run's cells; the land in the
run's water cells given the nearest land cell's class and the shares over the design land adding up; the coasts' tide
from the design surface's own seas; the design surface and the tides matching their hashes)."""
import json
import math

import numpy as np
import pytest

from research.studies.atlas_round import moisture_final as MF


def test_a_years_sums_move_by_the_step_over_land_and_water():
    land = np.array([[True, False]])
    cnt = np.full((10, 1, 2), 36.0)                      # records in each tenth
    acc = {k: np.full((10, 1, 2), 36.0 * v) for k, v in
           dict(pr=2e-5, evap=-1e-5, mrso=0.1, rss=200.0, rls=-60.0, tas=300.0, ps=1215.0, ts=301.0).items()}
    year = dict(land=land, acc=acc, cnt=cnt, rh_sum=np.full((1, 2), 0.8 * 180), rh_n=np.full((1, 2), 180.0),
                **{k: np.full((1, 2), 290.0) for k in MF.EXTREMES})
    applied = {k: {'land': {'per_K': 0.02 if kind == 'ratio' else 1.1},
                   'water': {'per_K': 0.01 if kind == 'ratio' else 0.9}} for k, kind in MF.YEAR_KINDS.items()}
    out = MF.adjusted(year, applied, -0.5)
    mean = lambda a: a[0] / cnt[0]                       # noqa: E731  the first tenth's mean
    assert np.allclose(mean(out['acc']['pr']), [2e-5 * 0.99, 2e-5 * 0.995])
    assert np.allclose(mean(out['acc']['evap']), [-1e-5 * 0.99, -1e-5 * 0.995])
    assert np.allclose(mean(out['acc']['tas']), [300.0 - 0.55, 300.0 - 0.45])
    assert np.allclose(out['tas_min'], [[290.0 - 0.55, 290.0 - 0.45]])
    assert np.allclose(out['rh_sum'] / out['rh_n'], [[0.8 - 0.55, 0.8 - 0.45]])
    assert out['cnt'] is year['cnt'] and np.allclose(year['acc']['tas'][0] / cnt[0], 300.0)   # the year itself unchanged


def test_a_map_pixel_is_sea_where_most_of_it_lies_under_one_sea():
    rows, cols = 8, 16
    body = np.zeros((rows, cols), np.int8)
    body[0:4, 0:4] = 3                          # a whole pixel of one sea
    body[0:2, 4:8] = 1                          # half a pixel of sea,
    body[2, 4:8] = MF.LAKE                      # a quarter lake and a quarter land
    body[4, 0:4] = 5                            # a quarter of a pixel under a sea
    out = MF.map_grid(dict(body=body.ravel(), rows=rows, cols=cols))
    assert out.tolist() == [[3, 1, 0, 0], [0, 0, 0, 0]]


@pytest.fixture(scope='module')
def record():
    if not MF.OUT.exists():
        pytest.skip('run python -m research.studies.atlas_round.moisture_final first')
    return json.loads(MF.OUT.read_text())


def test_the_record_holds_together(record):
    p = record['producer']
    assert p['run'] == MF.RUN and p['years'] == list(MF.YEARS)
    frame_path, frame = MF.frame_record()
    assert p['frame']['orientation'] == str(frame_path.relative_to(MF.ROOT))
    assert p['frame']['tilt_deg'] == pytest.approx(frame['frame']['tilt_to_ecliptic_deg'])
    geography = MF.bm.GEOGRAPHY                                           # the settings come from the design surface,
    assert geography['atlas'] not in p['inputs'] and geography['atlas_grid'] not in p['inputs']
    assert str(MF.TF.PRODUCT) in p['inputs'] and geography['tides_grid'] not in p['inputs']   # the design seas' tide
    shares = record['shares']['by_class']
    assert sum(r['share'] for r in shares.values()) == pytest.approx(1.0, abs=5e-4)
    assert record['checks']['gap_never_wetter'] and record['checks']['robust_within_gap_and_years']
    soda = record['settings']['soda_lakes']
    assert soda['sea_level_lake_km2_in_dry_cells'] == 0
    assert (soda['closed_rainfed_lake_km2_in_robustly_dry_cells'] <= soda['closed_rainfed_lake_km2_in_dry_cells']
            <= soda['closed_rainfed_lake_km2_on_land'] <= soda['rainfed_lake_km2_on_land'])
    cover = record['land']['design_surface']
    held = cover['design_land_in_run_land_km2'] + cover['design_lake_in_run_land_km2'] + cover['design_sea_in_run_land_km2']
    assert held == pytest.approx(cover['run_land_km2'], rel=0.01)
    assert (cover['design_land_in_run_land_km2'] + cover['design_land_in_run_water_km2']
            == pytest.approx(cover['design_land_km2'], abs=2))
    summary = json.loads(MF.GS.OUT.read_text())['summary']
    moon_km2 = 4.0 * math.pi * (MF.MOON_RADIUS / 1e3) ** 2
    assert cover['design_land_km2'] / moon_km2 == pytest.approx(1.0 - summary['water_share'], abs=0.01)
    surface = p['design_surface']
    assert surface['record_sha256'] == MF.sha256(MF.GS.OUT)
    if MF.GS.PRODUCT.exists():
        assert MF.sha256(MF.GS.PRODUCT) == surface['sha256']
    assert p['tides']['record_sha256'] == MF.sha256(MF.TF.OUT)


def test_the_land_in_the_runs_water_cells_takes_the_nearest_land_class(record):
    cover = record['land']['design_surface']
    given = cover['assigned']
    assert given['km2'] == pytest.approx(cover['design_land_in_run_water_km2'], abs=2)
    assert sum(given['by_class_km2'].values()) == pytest.approx(given['km2'], abs=len(given['by_class_km2']))
    assert 0 < given['nearest_land_km']['median'] <= given['nearest_land_km']['p90'] <= given['nearest_land_km']['maximum']
    over = record['shares']['over_the_design_land']
    assert over['land_km2'] == pytest.approx(cover['design_land_km2'], abs=2)
    assert sum(over['by_class'].values()) == pytest.approx(1.0, abs=5e-4)
    assert sum(over['by_class_in_run_land_cells_only'].values()) == pytest.approx(1.0, abs=5e-4)
    assert over['assigned_share'] == pytest.approx(cover['design_land_unclassed_share'], abs=1e-3)
    keys = [c['key'] for c in record['classes']]
    assert over['drylands_share'] == pytest.approx(sum(over['by_class'][k] for k in keys[:4]), abs=5e-4)
