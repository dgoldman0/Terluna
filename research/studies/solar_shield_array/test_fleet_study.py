import json
import numpy as np
import pytest

from shared import constants as K
from shared.provenance import constants_changed, constants_used
from protection.dynamics.fleet import initial_fleet
from .fleet_run import digest, ROOT, HERE


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
