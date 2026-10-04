"""Local-product provenance and comparison invariants, including remote use."""
import json

import numpy as np
import pytest

from shared.constants import SYNODIC_MONTH_DAYS
from shared.provenance import constants_changed
from .angular import PackedAtlas
from .local_sky import NativeAtlas
from .run import HERE, ROOT, digest
from .water_slopes import refine


def product(name):
    value = json.loads((HERE / 'results' / name).read_text())
    for path, expected in value['producer']['files'].items():
        assert digest(ROOT / path) == expected
    assert not constants_changed(value['producer'].get('constants'))
    return value


def test_native_sky_comparison_preserves_scene_and_flux_interpretation():
    p = product('native_sky.json')
    old = json.loads((HERE / 'results/directional_scenes.json').read_text())
    assert len(p['angular_views']) == 630
    assert len(p['scenes']) == 24
    assert p['numerical_settings'] == old['numerical_settings']
    for a, b in zip(p['angular_views'], old['angular_views'], strict=True):
        assert a['ground_lux'] == b['ground_lux']
        assert a['sun_direct_lux'] == b['sun_direct_lux']
        assert a['total_lux'] == pytest.approx(a['ambient_lux'] + a['sun_direct_lux'])
    for r in p['recovery_audit']:
        assert r['beam_and_flux_tables_exactly_equal']
        assert r['maximum_shared_sample_relative_difference'] < .0005
        assert r['atmosphere']['ground_albedo'] == .1
        assert r['atmosphere']['visible_filter'] == 1
    assert p['checks']['maximum_scene_total_relative_change'] < .001


@pytest.mark.skipif(not (ROOT / 'illumination/sky/data/earth_atlas.npz').exists(),
                    reason='Native sky archives are local ignored inputs; see prepare_local.py')
def test_native_grid_uses_archived_axes_and_exact_solar_frames():
    atlas = NativeAtlas()
    with pytest.raises(ValueError):
        atlas.get('moon', 75.)
    for world in ('moon', 'earth'):
        full, a = atlas.get(world, 30)
        packed, b = PackedAtlas().get(world, 30)
        assert a['stored_diffuse_lux'] == b['stored_diffuse_lux']
        assert full.horizontal > 0
        assert abs(full.horizontal / packed.horizontal - 1) < .004


def test_frequency_refinement_preserves_interpolated_wave_variance():
    f = np.array([.01, .03, .2, .5])
    e = np.array([[0, 2], [1, 3], [2, 1], [0, 0.]])
    x, y = refine(f, e)
    assert np.trapezoid(y, x, axis=0) == pytest.approx(np.trapezoid(e, f, axis=0))


def test_saved_slopes_use_second_cycle_and_bound_their_frequency_support():
    p = product('water_slopes.json')
    lo, hi = p['sampled_span_hours']
    assert lo >= SYNODIC_MONTH_DAYS * 24
    assert hi <= 2 * SYNODIC_MONTH_DAYS * 24
    assert p['checks']['maximum_dispersion_relative_residual'] < 1e-12
    assert p['checks']['maximum_slope_integral_relative_change_2_to_4'] < .005
    assert p['frequency_range_hz'][1] < .51
    for s in p['summaries'].values():
        assert s['snapshots'] == int(hi - lo) + 1
        assert 0 < s['rms_gradient']['minimum'] <= s['rms_gradient']['maximum']
        for k in ('slope_fraction_above_025hz', 'elevation_variance_fraction_above_025hz'):
            assert 0 <= s[k]['minimum'] <= s[k]['maximum'] <= 1
