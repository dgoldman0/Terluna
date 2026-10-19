"""The rain-rate product binds its code and inputs, agrees with the climate targets' land rain, and holds together."""
import hashlib
import json

import pytest

from atmosphere.rain import rates, targets as climate_targets

PRODUCT = json.loads(rates.OUT.read_text())
DATA, BINDING = climate_targets.load()
READ = climate_targets.read(DATA)
STATS = ('p10', 'p50', 'p90', 'p99', 'top')


def test_header_and_bindings():
    assert PRODUCT['schema'] == rates.SCHEMA
    assert PRODUCT['evidence'] and PRODUCT['reading_rule'] and PRODUCT['assumptions']
    producer = PRODUCT['producer']
    assert set(producer['files']) == set(rates.CODE)
    for kind, entries in BINDING.items():
        assert producer[kind] == entries, kind
    for path, digest in {**producer['files'], **producer.get('inputs', {})}.items():
        assert hashlib.sha256((rates.ROOT / path).read_bytes()).hexdigest() == digest, path
    # The GCM's output lies on the research drive; its hashes are recorded, and read again only when it rebuilds.
    read = producer['read']
    for run, years in ((rates.DESIGN, rates.YEARS), (rates.TEST, rates.YEARS), (rates.SHORT, rates.SHORT_YEARS)):
        for year in range(years[0], years[1] + 1):
            assert len(read[str(rates.output(run, year).relative_to(rates.ROOT))]) == 64
    assert DATA['schema'] == climate_targets.SCHEMA


def test_the_land_rain_is_the_climate_targets():
    runs = {r['run']: r for r in DATA['runs'].values() if isinstance(r, dict) and 'run' in r}
    assert set(runs) == {rates.DESIGN, rates.TEST}
    for name, run in runs.items():
        assert PRODUCT['runs'][name]['land_rain_mm_day'] == pytest.approx(run['land_rain_mm_day'], rel=1e-9)
    for key in rates.TARGETS:
        target, theirs = PRODUCT['targets'][key], READ[key]
        assert target['share_of_the_way'] == theirs['share_of_the_way']
        assert target['three_day_means']['land_rain_mm_day'] == pytest.approx(theirs['land_rain_mm_day'], rel=1e-9)
        assert target['land_rain_with_the_dimmings_own_cut_mm_day'] == theirs[
            'land_rain_with_the_dimmings_own_cut_mm_day']
    # The resort is the resort band's candidate at 301 K; the lock is the design case and the hot world the test.
    shares = {key: PRODUCT['targets'][key]['share_of_the_way'] for key in rates.TARGETS}
    assert shares['lock'] == 0.0 < shares['resort'] < 1.0 == shares['hot_world']
    assert (PRODUCT['targets']['resort']['ground_k'], PRODUCT['targets']['resort']['held_in']) == (
        301.0, 'resort_band.rows')


def test_the_dimmings_own_cut_is_carried():
    """The hot world needs no dimming; the resort's and the lock's dimming cut their land's mean rain by a few to
    about sixteen per cent."""
    assert PRODUCT['targets']['hot_world']['dimmings_own_cut_percent'] == [0.0, 0.0]
    for key in ('resort', 'lock'):
        least, most = PRODUCT['targets'][key]['dimmings_own_cut_percent']
        assert 0 < least < most < 20
    assert PRODUCT['targets']['lock']['dimmings_own_cut_percent'][1] > PRODUCT['targets']['resort'][
        'dimmings_own_cut_percent'][1]


def test_the_short_years_match_the_design_cases_3_day_means():
    """The 6.5-hour years averaged eleven at a time match the design case's own 3-day means in distribution, within
    0.8% in the wet share and the quantiles to the 99th, as the READMEs state."""
    short = PRODUCT['short']
    assert short['means_to_a_3_day_mean'] == 11 and short['window_hours'] == pytest.approx(6.5)
    for key in ('share', 'p10', 'p50', 'p90', 'p99'):
        assert short['blocks_over_the_design_cases_own'][key] == pytest.approx(1.0, abs=0.008), key


def test_short_weather_is_wetter_while_rarer_and_the_warmer_targets_wetter():
    ratio = PRODUCT['short']['ratio']
    assert ratio['share'] < 0.7 and all(ratio[k] > 1.0 for k in ('p50', 'p90', 'p99', 'top'))
    for target in PRODUCT['targets'].values():
        for window in ('three_day_means', 'six_and_a_half_hour_means'):
            wet = target[window]['wet']
            assert all(wet[a] < wet[b] for a, b in zip(STATS, STATS[1:]))
            assert wet['p10'] >= PRODUCT['wet_threshold_mm_h']
    hot, resort, lock = (PRODUCT['targets'][k] for k in ('hot_world', 'resort', 'lock'))
    for window in ('three_day_means', 'six_and_a_half_hour_means'):
        assert hot[window]['wet']['share'] > resort[window]['wet']['share'] > lock[window]['wet']['share']
        assert hot[window]['wet']['p50'] > resort[window]['wet']['p50'] > lock[window]['wet']['p50']
