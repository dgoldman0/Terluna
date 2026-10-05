"""Keep the local optical ledger distinct from an invented fleet power supply."""
import json
from pathlib import Path

import pytest

from shared.provenance import constants_changed
from .fleet_run import digest

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]


def product():
    return json.loads((HERE/'results/collection.json').read_text())


def test_collection_product_binds_sources_and_saved_trajectory_inputs():
    p = product()
    assert p['completed']
    for path, expected in p['producer']['source_hashes'].items():
        assert digest(ROOT/path) == expected, path
    for path, expected in p['producer']['inputs'].items():
        assert digest(HERE/'results'/path) == expected, path
    assert not constants_changed(p['producer']['constants'])
    assert p['maximum_force_ledger_disagreement_m_s2'] < 1e-14


def test_collection_ledger_preserves_energy_and_leaves_hardware_yield_open():
    p = product(); phases = p['stages'][-1]['phases']
    for row in phases.values():
        first = row['first_intercept_bolometric_equivalent_W']['energy_J']
        front = row['front_first_intercept_W']['energy_J']
        back = row['back_first_intercept_W']['energy_J']
        assert front+back == pytest.approx(first, rel=1e-12)
        assert row['redirected_band_optical_W']['energy_J'] == pytest.approx(
            first*p['redirected_spectral_fraction'], rel=1e-12)
        assert row['unshadowed_aperture_W']['energy_J'] == pytest.approx(
            first+row['mutual_shadow_equivalent_loss_W']['energy_J'], rel=1e-12)
    parts = ['service', 'departure_preparation', 'departure_turn', 'departure_coast']
    assert sum(phases[k]['redirected_band_optical_W']['energy_J'] for k in parts) == pytest.approx(
        phases['whole_saved_interval']['redirected_band_optical_W']['energy_J'], rel=1e-12)
    assert phases['departure_coast']['redirected_band_optical_W']['energy_J'] > 0
    for field in ['collection_capacity_W', 'delivered_power_W', 'fleet_propulsion_W',
                  'electrical_conversion_efficiency', 'outgoing_light_capture_fraction']:
        assert p[field] is None
    assert not p['accepted_fleet'] and not p['recurring_cycle_measured']
