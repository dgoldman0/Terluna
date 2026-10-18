"""The rain-rate product binds its code and inputs, agrees with the climate targets' land rain, and holds together."""
import hashlib
import json

import pytest

from atmosphere.rain import rates

PRODUCT = json.loads(rates.OUT.read_text())
TARGETS = json.loads(rates.TARGETS_PRODUCT.read_text())
STATS = ('p10', 'p50', 'p90', 'p99', 'top')


def test_header_and_bindings():
    assert PRODUCT['schema'] == rates.SCHEMA
    assert PRODUCT['evidence'] and PRODUCT['reading_rule'] and PRODUCT['assumptions']
    producer = PRODUCT['producer']
    assert set(producer['files']) == set(rates.CODE)
    for path, digest in {**producer['files'], **producer['inputs']}.items():
        assert hashlib.sha256((rates.ROOT / path).read_bytes()).hexdigest() == digest, path
    # The GCM's output lies on the research drive; its hashes are recorded, and read again only when it rebuilds.
    read = producer['read']
    for run, years in ((rates.DESIGN, rates.YEARS), (rates.TEST, rates.YEARS), (rates.SHORT, rates.SHORT_YEARS)):
        for year in range(years[0], years[1] + 1):
            assert len(read[str(rates.output(run, year).relative_to(rates.ROOT))]) == 64
    assert TARGETS['schema'] == rates.TARGETS_SCHEMA


def test_the_land_rain_is_the_climate_targets():
    runs = {r['run']: r for r in TARGETS['runs'].values() if isinstance(r, dict) and 'run' in r}
    assert set(runs) == {rates.DESIGN, rates.TEST}
    for name, run in runs.items():
        assert PRODUCT['runs'][name]['land_rain_mm_day'] == pytest.approx(run['land_rain_mm_day'], rel=1e-9)
    by_name = {t['name']: t for t in TARGETS['targets']}
    for key, name in rates.TARGETS:
        target = PRODUCT['targets'][key]
        assert target['share_of_the_way'] == by_name[name]['interpolated_share_of_the_way']
        assert target['three_day_means']['land_rain_mm_day'] == pytest.approx(by_name[name]['land_rain_mm_day'],
                                                                              rel=1e-9)


def test_the_short_years_make_the_design_cases_3_day_means():
    short = PRODUCT['short']
    assert short['means_to_a_3_day_mean'] == 11 and short['window_hours'] == pytest.approx(6.5)
    for key in ('share', 'p10', 'p50', 'p90', 'p99'):
        assert short['blocks_over_the_design_cases_own'][key] == pytest.approx(1.0, abs=0.02), key


def test_short_weather_is_wetter_while_rarer_and_the_warmer_targets_wetter():
    ratio = PRODUCT['short']['ratio']
    assert ratio['share'] < 0.7 and all(ratio[k] > 1.0 for k in ('p50', 'p90', 'p99', 'top'))
    for target in PRODUCT['targets'].values():
        for window in ('three_day_means', 'six_and_a_half_hour_means'):
            wet = target[window]['wet']
            assert all(wet[a] < wet[b] for a, b in zip(STATS, STATS[1:]))
            assert wet['p10'] >= PRODUCT['wet_threshold_mm_h']
    hot, lock = PRODUCT['targets']['hot_world'], PRODUCT['targets']['lock']
    for window in ('three_day_means', 'six_and_a_half_hour_means'):
        assert hot[window]['wet']['share'] > lock[window]['wet']['share']
        assert hot[window]['wet']['p50'] > lock[window]['wet']['p50']
