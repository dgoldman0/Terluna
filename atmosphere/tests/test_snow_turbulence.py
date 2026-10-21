"""Controls for snow conditioning, Galilean invariance and dissipation scaling."""
import hashlib
import json

import numpy as np
import pytest

from atmosphere.snow import turbulence as st
from shared.provenance import constants_changed


def test_snow_masks_do_not_call_clear_cold_air_snow_or_warm_snow_a_growth_cloud():
    t = np.array([260., 260., 280., 245., 273.15])
    q = np.array([0., 1e-4, 1e-4, 1e-4, 1e-5])
    masks = st.winter_masks(t, q)
    assert masks['all_winter_air'].tolist() == [True, True, False, False, True]
    assert masks['snow_ge_1e-05_kg_kg'].tolist() == [False, True, False, False, True]
    assert masks['snow_ge_0.0001_kg_kg'].tolist() == [False, True, False, False, False]


def test_uniform_wind_translation_cannot_become_flake_scale_deformation():
    z = np.array([100., 300., 700.])
    shape = (3, 4, 5)
    u = np.broadcast_to((.01 * z)[:, None, None], shape).copy()
    fields = dict(uinterp=u, vinterp=np.zeros(shape), winterp=np.zeros(shape))
    before = st.resolved_quantities(fields, z, 6000.)
    after = st.resolved_quantities(dict(fields, uinterp=u + 1000.), z, 6000.)
    for key in ('vertical_horizontal_wind_shear_s-1', 'horizontal_velocity_gradient_norm_s-1'):
        assert np.allclose(before[key], after[key], atol=1e-14)
    assert np.allclose(before['vertical_horizontal_wind_shear_s-1'], .01)
    assert np.all(before['horizontal_velocity_gradient_norm_s-1'] == 0)
    assert np.all(after['horizontal_wind_m_s'] > before['horizontal_wind_m_s'])


def test_increment_scaling_has_correct_viscous_and_inertial_limits_and_no_jump():
    epsilon, nu = 1e-3, 2e-5
    eta = (nu**3 / epsilon)**.25
    at = st.increment_scale(epsilon, eta, nu)
    assert st.increment_scale(epsilon, eta / 2, nu) == pytest.approx(at / 2)
    assert st.increment_scale(epsilon, eta * 8, nu) == pytest.approx(at * 2)
    assert st.increment_scale(epsilon, eta * (1 - 1e-9), nu) == pytest.approx(at)
    assert st.increment_scale(0., .6, nu) == 0.
    for args in ((-1., .6, nu), (epsilon, 0., nu), (epsilon, .6, 0.)):
        with pytest.raises(ValueError):
            st.increment_scale(*args)


def test_geometry_generalisation_preserves_the_legacy_fall_law_as_a_control():
    c = st.sf.earth_column()
    state = c.state(c.z_freeze)
    geometry = st.GEOMETRIES['legacy_bf95']
    for diameter in (.001, .01, .15, .6):
        expected = st.sf.fall_speed(diameter, st.sf.flake_mass(diameter), c.gravity, state)[0]
        assert st.terminal_speed(diameter, geometry, c.gravity, state) == pytest.approx(expected)


def test_selected_grads_reader_respects_levels_and_rejects_missing_bytes(tmp_path):
    shape = dict(x=3, y=2, z=4)
    layout = {name: (i * 24, 4) for i, name in enumerate(st.FIELDS)}
    raw = np.arange(len(st.FIELDS) * 24, dtype='<f4')
    path = tmp_path / 'data.bin'
    path.write_bytes(raw.tobytes())
    out = st.read_selected(path, shape, layout, raw.size, np.array([1, 3]), hashlib.sha256())
    assert out['th'].shape == (2, 2, 3)
    assert out['th'][0, 0, 0] == 6
    assert out['th'][1, -1, -1] == 23
    assert out['winterp'][1, -1, -1] == raw[-1]
    path.write_bytes(raw[:-1].tobytes())
    with pytest.raises(ValueError, match='size mismatch'):
        st.read_selected(path, shape, layout, raw.size, np.array([1, 3]), hashlib.sha256())


def test_empty_mask_has_no_fabricated_statistics():
    empty = st.statistics([])
    assert empty['count'] == 0
    assert empty['p99'] is None
    assert st.statistics([1., 2., 3.])['mean'] == 2.


def test_growth_environment_requires_both_moisture_and_liquid_when_named():
    masks = st.winter_masks(np.full(4, 263.), np.full(4, 1e-4), np.array([.8, 1.1, 1.1, .8]),
                            np.array([0., 0., 2e-5, 2e-5]))
    assert masks['snow_and_ice_supersaturation'].tolist() == [False, True, True, False]
    assert masks['snow_and_liquid_cloud'].tolist() == [False, False, True, True]
    assert masks['snow_and_liquid_and_ice_supersaturation'].tolist() == [False, False, True, False]


def test_contiguous_growth_depth_does_not_bridge_a_dry_layer():
    mask = np.array([[True, True], [True, False], [False, True], [True, True]])
    assert st.longest_layer_depth(mask, [1., 2., 3., 4.]).tolist() == [4., 7.]


def test_checked_product_provenance_and_inference_boundaries():
    product = json.loads(st.OUT.read_text())
    assert product['schema'] == st.SCHEMA
    assert not constants_changed(product['producer']['constants'])
    for name, expected in product['scenarios']['column_input']['inputs'].items():
        assert st.sha256(st.ROOT / name) == expected, name
    assert product['source_register_sha256'] == st.sha256(st.HERE / 'turbulence_sources.json')
    for name, expected in product['producer']['files'].items():
        assert st.sha256(st.ROOT / name) == expected, name
    assert 'Uncalibrated' in product['scenarios']['evidence']
    assert 'not dissipation estimates' in product['reading_rule']
    for case in product['cases']:
        assert case['turbulence_options']['ipbl'] == 2
        assert case['summaries']['all_winter_air']['sample_count'] > 0
        for group, summary in case['summaries'].items():
            assert 0 <= summary['share_of_winter_cells'] <= 1
            if group.startswith('snow_'):
                assert summary['sample_count'] <= case['summaries']['all_winter_air']['sample_count']
