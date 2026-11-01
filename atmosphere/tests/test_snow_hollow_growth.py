"""Source reproduction, coordinate-domain and measurement-error controls."""
import itertools
import shutil
import numpy as np
import pytest
from atmosphere.snow import hollow_shape as shape
from atmosphere.snow import hollow_growth_assessment as assessment


@pytest.fixture(scope='module')
def library(tmp_path_factory):
    if not shape.ARCHIVE.exists() or shutil.which('gfortran') is None:
        pytest.skip('Optional pinned 2026 source cache and gfortran required')
    build = tmp_path_factory.mktemp('hollow_shape') / 'build'
    shape.build(build)
    return shape.ShapeLibrary(build)


def test_original_shape_tables_and_wrapper(library):
    result = shape.audit(library)
    assert sum(row['numerical_entries'] for row in result['comparisons']) == 17000
    assert sum(row['wrapper_evaluations'] for row in result['comparisons']) == 1992
    assert [row['out_of_fit_hexagonal_requests'] for row in result['comparisons']] == [0, 8]


@pytest.mark.parametrize('phi', [.013, .2, 1., 6., 70.])
def test_hexagonal_coordinate_change_off_grid(library, phi):
    native = library.evaluate(phi, 'hexagonal_prism', 'native')
    wood = library.evaluate(phi / shape.NATIVE_COSINE, 'hexagonal_prism', 'wood')
    # Change of basal radius changes the normalization of prism-source terms.
    # These are six independent published fits, sampled off the demo grid.
    for field in shape.FIELDS[:6]:
        conversion = shape.NATIVE_COSINE if field.startswith('prism') else 1.
        assert native[field] == pytest.approx(conversion * wood[field], rel=2e-6)


@pytest.mark.parametrize('phi', [0., -.2, .009999, 100.001, float('inf'), float('nan')])
def test_no_silent_aspect_ratio_clipping(phi):
    with pytest.raises(ValueError):
        shape.validate_inputs(phi, 'cylinder', 'native')


def test_converted_fit_domain_checked_before_native_call():
    shape.validate_inputs(86.6, 'hexagonal_prism', 'native')
    with pytest.raises(ValueError, match='Converted Wood'):
        shape.validate_inputs(86.61, 'hexagonal_prism', 'native')
    shape.validate_inputs(100., 'hexagonal_prism', 'wood')
    shape.validate_inputs(100., 'cylinder', 'native')


def test_rejects_changed_original_archive(tmp_path):
    changed = tmp_path / 'wrong.zip'
    changed.write_bytes(b'not the published source')
    with pytest.raises(ValueError, match='hash mismatch'):
        shape.archive_members(changed)


def test_slope_interval_encloses_every_error_box_vertex():
    time = np.array([0., 1., 4., 9.])
    value = 5. + 3. * time
    error = np.array([.1, .8, .2, .7])
    result = assessment.line_fit(time, value, error)
    slopes = []
    for signs in itertools.product((-1, 1), repeat=4):
        perturbed = value + error * signs
        slopes.append(float(np.polyfit(time, perturbed, 1)[0]))
    assert result['slope'] == pytest.approx(3.)
    assert result['interval'] == pytest.approx([min(slopes), max(slopes)])


def test_slope_does_not_depend_on_clock_origin():
    time = np.array([0., 2., 5., 11.])
    value = np.array([3., 4., 6., 7.])
    first = assessment.line_fit(time, value, [.1]*4)
    shifted = assessment.line_fit(time + 1e6, value, [.1]*4)
    assert shifted['slope'] == pytest.approx(first['slope'])
    assert shifted['measurement_error_bound'] == pytest.approx(first['measurement_error_bound'])


def test_unresolved_growth_ratio_is_not_reported_as_finite_interval():
    data = np.zeros((4, 10))
    data[:, 0] = [0., 60., 120., 180.]
    data[:, 1] = [10., 10.1, 10.2, 10.3]
    data[:, 2] = [30., 32., 34., 36.]
    data[:, [5, 6, 9]] = 2.
    data[:, 7] = 3.
    result = assessment.window(data, 0., 4.)
    assert result['fitted_directional_growth_ratio'] > 0
    assert result['directional_ratio_interval'] is None


@pytest.mark.parametrize('time,value,error', [([0,0,1],[1,2,3],[0,0,0]),
    ([0,1,2],[1,2,3],[0,-1,0]), ([0,1,2],[1,float('nan'),3],[0,0,0]),
    ([0,1],[1,2],[0,0]), ([0,1,2],[1,2],[0,0,0])])
def test_invalid_measurement_series_rejected(time, value, error):
    with pytest.raises(ValueError):
        assessment.line_fit(time, value, error)


def test_archive_measurement_identity_and_duplicate_points():
    if not shape.ARCHIVE.exists():
        pytest.skip('Optional pinned source archive required')
    members = shape.archive_members()
    warm = assessment.free_fall_controls(members)
    assert warm['compared_measurements'] == 69
    assert warm['duplicate_early_column_points'] == {'a': 6, 'c': 5}
    cold = assessment.chamber_controls(members)
    assert [row['measurements'] for row in cold['cases']] == [26, 68]
    assert cold['cases'][1]['switch_elapsed_s'] == 13800
