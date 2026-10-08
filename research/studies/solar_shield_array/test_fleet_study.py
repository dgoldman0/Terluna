import json
import numpy as np
import pytest

from shared import constants as K
from shared.provenance import constants_changed, constants_used
from protection.dynamics.fleet import initial_fleet
from .fleet_run import digest, ROOT, HERE
from .fleet_handover import hermite_state


class StaticFrame:
    def frames(self, t):
        return np.eye(3)


def test_fleet_has_distinct_simultaneous_3d_initial_states():
    state, labels = initial_fleet(StaticFrame(), [15e6, 30e6], 12, 48)
    assert len(state) == 1152
    assert len(np.unique(state[:, :3], axis=0)) == len(state)
    assert np.linalg.matrix_rank(state[:, :3]) == 3
    energy = np.sum(state[:, 3:]**2, axis=1)/2-K.MOON_GM/np.linalg.norm(state[:, :3], axis=1)
    np.testing.assert_allclose(energy[:576], -K.MOON_GM/(2*15e6))
    assert np.unique(labels[:, 0]).tolist() == [0, 1]


def test_fleet_products_bind_sources_and_do_not_claim_payload_or_continuous_coverage():
    products = list((HERE/'results').glob('fleet_*.json'))
    if not products:
        pytest.skip('Fleet product has not yet been generated')
    for path in products:
        p = json.loads(path.read_text())
        if p.get('schema') != 'terluna.research.solar-shield-fleet/1':
            continue
        for name, expected in p['producer']['source_hashes'].items():
            assert digest(ROOT/name) == expected, name
        assert not constants_changed(p['producer']['constants'])
        assert p['producer']['constants'] == constants_used(p['producer']['source_hashes'])
        assert p['areal_mass_kg_m2'] == .05
        assert p['acceptance']['continuous_full_coverage_demonstrated'] is False
        assert p['acceptance']['delivered_electrical_power_demonstrated_W'] is None
        for case in p['cases'].values():
            strict = case['selections']['collision_and_shadow_screened']
            assert strict['collision_pairs_retained'] == 0
            assert strict['possible_shadow_pairs_retained'] == 0
            assert strict['summary']['excess_intersection_area_fraction']['maximum'] < 1e-8
            assert case['constraints']['minimum_earth_beam_margin_deg'] > 0
            assert case['constraints']['minimum_moon_beam_margin_deg'] > 0
            assert case['control']['electric_impulse_m_s'] == 0
            for result in case['selections'].values():
                for key in ('ray_coverage', 'core_ray_coverage', 'all_sampled_sun_coverage'):
                    s = result['summary'][key]
                    assert 0 <= s['minimum'] <= s['mean'] <= s['maximum'] <= 1+1e-12
                assert result['optical_mass_kg'] == pytest.approx(result['physical_area_m2']*.05)


def test_handover_interpolation_recovers_a_known_cubic_trajectory():
    times = np.array([0., 2.])
    state = np.zeros((2, 1, 6))
    state[:, 0, 0] = times**3-2*times
    state[:, 0, 3] = 3*times**2-2
    for t in [0., .3, 1., 1.7, 2.]:
        value = hermite_state(t, times, state)[0]
        assert value[0] == pytest.approx(t**3-2*t)
        assert value[3] == pytest.approx(3*t*t-2)


def test_followup_products_bind_sources_and_preserve_their_evidence_states():
    for path in (HERE/'results').glob('fleet_*.json'):
        p = json.loads(path.read_text())
        if 'producer' not in p:
            continue
        for name, expected in p['producer']['source_hashes'].items():
            assert digest(ROOT/name) == expected, (path.name, name)
        assert not constants_changed(p['producer']['constants'])
        for name, expected in p['producer'].get('input_products', {}).items():
            assert digest(HERE/'results'/name) == expected
        for key in ('input_product', 'template_product'):
            if key in p['producer']:
                assert digest(HERE/'results'/p['producer'][key]) == p['producer'][key+'_sha256']
        if p['schema'].endswith('fleet-validation/1'):
            for case in p['cases'].values():
                assert case['max_sun_direction_rate_rad_s'] < case['assumed_shadow_guard_rate_rad_s']
                if 'independent_replay' in case:
                    r = case['independent_replay']
                    assert r['trajectory_agreement_within_50m'] == (r['max_position_difference_m'] < 50)
        if p['schema'].endswith('fleet-refinement/1'):
            assert p['target_radii_m']['3R'] == 3*K.MOON_RADIUS
            assert p['target_radii_m']['4R'] == 4*K.MOON_RADIUS
            for case in p['cases'].values():
                for variant in case['variants'].values():
                    for r in variant['records']:
                        assert r['all_sampled_sun_coverage'] <= r['ray_coverage']+1e-12
