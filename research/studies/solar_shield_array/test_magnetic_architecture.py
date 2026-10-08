"""Verify published comparisons and their evidence boundaries."""
import hashlib
import json
from pathlib import Path
import numpy as np
from shared import constants as K
from shared.provenance import constants_changed

ROOT = Path(__file__).resolve().parents[3]
RESULT = json.loads((Path(__file__).parent/'results/magnetic_architecture.json').read_text())
HISTORICAL = json.loads((Path(__file__).parent/'historical/inputs.json').read_text())


def test_parent_and_producer_content_is_pinned():
    # Inputs read before their file was rederived are kept byte for byte under historical/, with a reading rule.
    kept = {h['original_path']: h for h in HISTORICAL['inputs']
            if 'research/studies/solar_shield_array/results/magnetic_architecture.json' in h['read_by']}
    for key in ['source_hashes', 'input_hashes']:
        for path, sha in RESULT['producer'][key].items():
            if path in kept:
                assert kept[path]['sha256'] == sha, path
                assert hashlib.sha256((ROOT/kept[path]['path']).read_bytes()).hexdigest() == sha, path
                assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() != sha, path
                continue
            assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest() == sha, path
    assert not constants_changed(RESULT['producer']['constants'])


def test_geometry_refinement_and_momentum_conservation():
    for row in RESULT['regional']:
        assert row['angular_grid_relative_change'] < .002
        if not row['name'].startswith('historical'):
            assert max(abs(x) for x in row['route_height_range_m']) < 1e-7
        for pair in row['reference']['mutual_loads']:
            assert pair['mutual_refinement_relative'] < 1e-7
            assert pair['force_reciprocity_relative'] < 1e-9
            assert pair['angular_reciprocity_relative'] < 1e-9
        for case in row['pressure_sizing']:
            b = row['min_B_at4R_T'] if case['target_radius_R']==4 else row['min_B_at3R_T']
            pressure = (b*case['current_scale'])**2/(2*K.VACUUM_PERMEABILITY)
            assert np.isclose(pressure, case['pressure_Pa'], rtol=1e-12, atol=0)


def test_failed_held_cases_and_atmospheric_limits_are_retained():
    near = [r for r in RESULT['holding'] if r['distance_m']==15e6]
    assert near and all(not b['sufficient_scalar_closure'] for r in near for b in r['budgets'])
    assert any(r['reduced_moment_total_loss_kg_s'] > 10 for r in RESULT['atmosphere_cases'])
    assert all(r['assumed_boundary_field_gain']==2. for r in RESULT['atmosphere_cases'])
    assert not RESULT['accepted_magnetosphere']
    assert not RESULT['accepted_trajectory'] and not RESULT['accepted_return'] and not RESULT['accepted_cycle']
    assert RESULT['original_target_W']==23.835e12 and not RESULT['target_changed']
    assert not RESULT['superconducting_planetary_ring_selected']


def test_resources_and_unqualified_plasma_sensitivity():
    r = RESULT['resources']
    assert r['elapsed_s'] < 900 and r['peak_rss_MiB'] < 2048 and r['raw_output_MiB'] < 64
    assert r['numerical_threads']==1 and r['processes']==1
    p = RESULT['plasma_current']
    assert not p['accepted']
    assert len(p['field_loss_power'])==4
    assert max(r['monthly_grid_peak_change'] for r in RESULT['holding'] if 'monthly_grid_peak_change' in r) < .003


def test_extended_loop_attitude_uses_rigid_body_inertia(monkeypatch):
    from . import magnetic_architecture_run as model
    monkeypatch.setattr(model, 'relative_gravity', lambda q, *args: np.zeros_like(q))
    alpha = np.array([0., 0., 2e-9]); frame = np.eye(3)
    acceleration = np.cross(alpha, frame)
    geo = dict(frame=frame[None], frame_dd=acceleration[None])
    a = model.held_acceleration(dict(positions={}, moon_a=np.zeros((1, 3))), geo, 15e6, 2e6, 'axial')
    np.testing.assert_allclose(a['net_vector_m_s2'][0], np.cross(alpha, [15e6, 0., 0.]), atol=1e-14)
    # Thin ring lies in yz; its per-mass z moment of inertia is radius squared/2.
    np.testing.assert_allclose(a['torque_per_kg_N_m'][0], alpha*(2e6)**2/2, atol=1e-10)
