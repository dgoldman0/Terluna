"""The biome map binds its inputs, its classes follow from the aridity it stores, and its uncertainty flags, shares and
checks agree with its cells."""
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from biosphere.ecology import biome_map as bm
from shared.provenance import constants_changed

ROOT = Path(__file__).resolve().parents[2]
PRODUCT = json.loads(bm.OUT.read_text())
CELLS = {k: np.array(v, dtype=float if k not in ('robust', 'robust_gap', 'robust_years', 'robust_store', 'robust_all',
                                                    'upland', 'night_rain_flag') else bool)
         for k, v in PRODUCT['cells'].items() if not isinstance(v[0], list)}


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def test_product_binds_producer_inputs_and_constants():
    assert PRODUCT['schema'] == bm.SCHEMA
    for path, h in PRODUCT['producer']['files'].items():
        assert digest(ROOT / path) == h, path
    for path, h in PRODUCT['producer']['inputs'].items():
        if isinstance(h, str):                         # repository files; the GCM's year files carry their fields' hash
            assert digest(ROOT / path) == h, path
        else:
            assert len(h['fields_read_sha256']) == 32 and h['bytes'] > 0, path
    assert not constants_changed(PRODUCT['producer']['constants'])
    for path, schema in PRODUCT['producer']['input_schemas'].items():
        assert json.loads((ROOT / path).read_text())['schema'] == schema, path


def test_classes_follow_the_stored_aridity():
    ai, se = CELLS['aridity'], CELLS['aridity_se']
    scale = PRODUCT['gap']['scale']
    assert np.array_equal(bm.biome(ai), CELLS['biome'])
    assert np.array_equal(bm.unep_class(ai), CELLS['unep_class'])
    assert np.array_equal(bm.biome(scale * ai), CELLS['biome_gap'])
    lo, hi = bm.bands(ai, se, scale)
    assert np.array_equal(bm.biome(lo), CELLS['biome_band_low'])
    assert np.array_equal(bm.biome(hi), CELLS['biome_band_high'])
    stores = [np.array(s, float) for s in PRODUCT['cells']['aridity_store']]
    rob = bm.robustness(ai, se, stores, scale)
    for key in ('robust', 'robust_gap', 'robust_years', 'robust_store', 'robust_all'):
        assert np.array_equal(rob[key], CELLS[key]), key
    for stored, s in zip(PRODUCT['cells']['biome_store'], stores):
        assert np.array_equal(bm.biome(s), np.array(stored))


def test_the_years_band_is_two_standard_errors_of_the_spans():
    spans = np.array(PRODUCT['cells']['aridity_spans'], float)
    se = spans.std(axis=0, ddof=1) / math.sqrt(spans.shape[0])
    assert np.allclose(se, CELLS['aridity_se'], rtol=2e-3, atol=1e-6)
    # the ten-year record weights each span by its records in each tenth of the day, so the spans' plain average
    # differs from it slightly; the ten-year value always lies within the spans' range
    ai = CELLS['aridity']
    assert np.all((ai >= spans.min(axis=0) * 0.999) & (ai <= spans.max(axis=0) * 1.001))
    assert np.median(np.abs(spans.mean(axis=0) / ai - 1.0)) < 0.01


def test_the_gap_dries_and_the_flags_nest():
    assert np.all(CELLS['biome_gap'] <= CELLS['biome'])
    assert np.all(CELLS['biome_band_low'] <= CELLS['biome'])
    assert np.all(CELLS['biome_band_high'] >= CELLS['biome'])
    robust = CELLS['robust']
    assert np.all(~robust | (CELLS['robust_gap'] & CELLS['robust_years']))
    assert np.array_equal(CELLS['robust_all'], robust & CELLS['robust_store'])
    assert PRODUCT['gap']['scale'] == 0.83
    ratios = PRODUCT['gap']['box_ratios'].values()
    assert all(0.8 < r < 0.85 for r in ratios)


def test_shares_agree_with_the_cells():
    w = CELLS['weight'] / CELLS['weight'].sum()
    rows = PRODUCT['shares']['by_biome']
    assert abs(sum(r['share'] for r in rows.values()) - 1.0) < 1e-3
    for k, com in enumerate(bm.BIOMES):
        row = rows[com['key']]
        here = CELLS['biome'] == k
        assert abs(w[here].sum() - row['share']) < 2e-4, com['key']
        assert abs(w[here & CELLS['robust']].sum() - row['robust_share']) < 2e-4, com['key']
        assert row['robust_share'] <= row['share'] + 1e-9
        assert abs(w[CELLS['biome_gap'] == k].sum() - row['share_at_gap_scale']) < 2e-4, com['key']
    assert abs(w[CELLS['robust']].sum() - PRODUCT['shares']['robust']) < 2e-4


def test_the_design_readings_are_reproduced():
    checks = PRODUCT['checks']
    assert checks['shares_sum'] == 1.0 or abs(checks['shares_sum'] - 1.0) < 1e-3
    assert checks['gap_never_wetter'] and checks['robust_within_gap_and_years']
    if 'design_years_aridity_shares_max_difference' in checks:   # read only where the research drive was mounted
        assert checks['design_years_aridity_shares_max_difference'] < 1e-9
        assert checks['design_years_land_rain_difference_mm_day'] < 1e-9
    if 'class_sensitivity_max_difference' in checks:
        assert checks['class_sensitivity_max_difference'] <= 5e-5 + 1e-12      # its shares are rounded to 4 places


def test_classes_at_their_bounds():
    assert list(bm.biome([0.0, 0.0499, 0.05, 0.1999, 0.2, 0.4999, 0.5, 0.6499, 0.65, 0.9999, 1.0, 7.0])) == \
        [0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5]
    assert list(bm.unep_class([0.0, 0.05, 0.2, 0.5, 0.65, 1.0, 7.0])) == [0, 1, 2, 3, 4, 4, 4]
    down, up = bm.margins([0.25, 0.04, 2.0])
    assert math.isclose(down[0], 1 - 0.2 / 0.25) and math.isclose(up[0], 0.5 / 0.25 - 1)
    assert math.isinf(down[1]) and math.isclose(up[1], 0.05 / 0.04 - 1) and math.isinf(up[2])
    assert list(bm.night_regime([0.001, 0.5, 2.0, 2.0], [-30.0, 10.0, 40.0, -10.0])) == [0, 1, 2, 3]
    lo, hi = bm.bands([0.5], [0.01], 0.83)
    assert math.isclose(lo[0], 0.83 * 0.48) and math.isclose(hi[0], 0.52)


def test_the_warm_land_has_no_frost():
    warmth = PRODUCT['warmth']
    assert warmth['coldest_interval_minimum_air_c'] > 0.0
    assert warmth['coldest_3day_mean_air_c'] >= warmth['coldest_interval_minimum_air_c']
    assert np.all(CELLS['coldest_3day_air_c'] <= CELLS['air_c']) and np.all(CELLS['warmest_3day_air_c'] >= CELLS['air_c'])


def test_night_light_and_settings_are_physical():
    lux = CELLS['darkest_lux']
    assert np.all(lux > 0) and np.all(np.isfinite(lux))
    near = CELLS['earth_elevation_deg'] > 0
    assert np.array_equal(bm.night_regime(lux, CELLS['earth_elevation_deg']), CELLS['night_regime'])
    assert np.all(CELLS['night_regime'][near] != 3) and np.all(CELLS['night_regime'][~near] != 2)
    for key in ('sea_fraction', 'sea_level_lake_fraction', 'rainfed_lake_fraction', 'closed_lake_fraction'):
        assert np.all((CELLS[key] >= 0) & (CELLS[key] <= 1.0 + 1e-6)), key
    assert np.all(CELLS['closed_lake_fraction'] <= CELLS['rainfed_lake_fraction'] + 1e-9)
    assert np.all(CELLS['coast_km'] >= 0)
    assert np.all((CELLS['night_rain_share'] >= 0) & (CELLS['night_rain_share'] <= 1))
