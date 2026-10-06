import json

import numpy as np

from shared.provenance import constants_changed
from protection.dynamics.relative_orbit import semi_major_axis
from .relative_diagnosis import EXPANDED, OUT, PACKAGE, ROOT, SEED, digest

PRODUCT = json.loads(OUT.read_text())


def test_product_is_current_with_its_code_inputs_and_constants():
    assert PRODUCT['schema'] == 'terluna.research.relative-diagnosis/1'
    producer = PRODUCT['producer']
    for name, value in producer['files'].items():
        assert digest(ROOT/name) == value, name
    inputs = {'results/joint_seed.json': SEED, 'local_continuation/ephemeris.npz': PACKAGE/'ephemeris.npz',
              'results/expanded_cycle.json': EXPANDED}
    for name, path in inputs.items():
        assert digest(path) == producer['inputs'][name], name
    assert not constants_changed(producer['constants'])


def test_energy_spread_follows_depth_and_matches_the_seed():
    seed = json.loads(SEED.read_text())['arrays']
    for key, entry in PRODUCT['energy_spread'].items():
        a = semi_major_axis(np.array(seed[key]))
        assert np.isclose(entry['semi_major_axis_span_m'], np.ptp(a))
    epoch = PRODUCT['energy_spread']['initial_epoch_state']
    assert epoch['fit_r_squared'] > 0.999
    assert 1.5 < epoch['fit_semi_major_axis_per_depth'] < 1.65
    assert 35e3 < epoch['depth_span_m'] < 37e3


def test_drift_is_secular_and_an_energy_match_removes_it():
    natural = PRODUCT['natural_repeat']
    assert natural['along_track_max_m'] > 300e3
    assert abs(natural['drift_slope_per_semi_major_axis']+natural['three_pi']) < 0.02*natural['three_pi']
    final = PRODUCT['energy_matching'][-1]
    assert final['along_track_max_m'] < 1e-3*natural['along_track_max_m']
    assert final['burn_rms_m_s'] < 0.5 and final['burn_max_m_s'] < 0.7


def test_the_pattern_folds_where_the_coupled_runs_first_meet_clearance():
    fold = PRODUCT['fold']
    first = min(fold['coupled_first_clearance_event_h'].values())
    for case in ('from_first_service_end', 'from_twelve_hour_state'):
        entry = fold[case]
        assert entry['thickness_min_m'] < 0.03*fold['from_first_service_end']['thickness_start_m']
        assert first < entry['thickness_min_h'] < first+1.
        assert first-.5 < entry['close_pairs_max_h'] < first+1.


def test_tide_turns_inclined_rings_with_the_sun_near_nineteen_thousand_km():
    cases = PRODUCT['plane_tracking']['cases']
    for tilt in (10., 20.):
        rows = sorted((c for c in cases if c['tilt_deg'] == tilt), key=lambda c: c['radius_km'])
        rates = [c['node_rate_deg_day'] for c in rows]
        assert all(np.diff(rates) > 0)
        nineteen = next(c for c in rows if c['radius_km'] == 19000.)
        assert abs(nineteen['node_minus_sun_deg']) < 10.
        fifteen = next(c for c in rows if c['radius_km'] == 15000.)
        assert fifteen['node_minus_sun_deg'] < -30.
