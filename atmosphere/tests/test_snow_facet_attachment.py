"""Numerical component checks, requiring the identified ignored Fortran builds."""
import math

import pytest

from atmosphere.snow import facet_attachment as facet


@pytest.fixture(scope='module')
def libraries():
    root = facet.ROOT / 'research/runs/snow_scale/facet_attachment'
    paths = {'single': root/'single_native_final', 'double': root/'double_native_final',
             'stable': root/'double_stable_final'}
    if not all((p / 'build.json').exists() for p in paths.values()):
        pytest.skip('Pinned native facet component builds are not restored')
    return {kind: facet.AttachmentLibrary(path) for kind, path in paths.items()}


def test_published_demonstration_reproduced_in_both_precisions(libraries):
    report = facet.audit(*(libraries[k].build_dir for k in ('single', 'double')))
    assert all(row['rows'] == 991 for row in report['native_demo_comparisons'])
    assert report['maximum_relative_single_double_differences']['alpha_c'] < .001


def test_native_near_sphere_cancellation_is_reported(libraries):
    a, c = 5e-6, 5.0000005e-6
    args = (.1, 258.15, 101000., a, c, 10.)
    with pytest.raises(FloatingPointError):
        libraries['single'].evaluate(*args)
    answer = libraries['double'].evaluate(*args)
    phi = c/a
    # Stable analytical prolate capacitance, independent of the native formula.
    expected = a * math.sqrt((phi-1)*(phi+1)) / math.acosh(phi)
    assert math.isclose(answer['capacitance_m'], expected, rel_tol=1e-8)
    assert answer['boundary_capacitance_m'] > answer['capacitance_m']


def test_overflow_is_local_and_next_evaluation_recovers(libraries):
    args = (1e-6, 258.15, 101000., 5e-6, 5e-6, 15.)
    with pytest.raises(FloatingPointError):
        libraries['single'].evaluate(*args)
    low = libraries['double'].evaluate(*args)
    assert low['alpha_a'] == low['alpha_c'] == 1e-6
    normal = libraries['single'].evaluate(.1, 258.15, 101000., 5e-6, 5e-6)
    assert 0 < normal['alpha_c'] < normal['alpha_a'] <= 1


def test_temperature_endpoints_and_guard(libraries):
    library = libraries['double']
    for tc in (-70., -1.):
        result = library.evaluate(.05, 273.15 + tc, 101000., 1e-5, 1e-5)
        assert result['ice_saturation_pressure_Pa'] > 0
    with pytest.raises(ValueError):
        library.evaluate(.05, 273.15 - 70.001, 101000., 1e-5, 1e-5)


def test_native_double_rejects_valid_near_sphere_and_stable_variant_repairs_it(libraries):
    args=(.013198201116886521,258.15,101000.,5.000000062836666e-6,5.000000030477008e-6,1.)
    with pytest.raises(FloatingPointError):
        libraries['double'].evaluate(*args)
    result=libraries['stable'].evaluate(*args)
    assert min(args[3:5])<=result['capacitance_m']<=max(args[3:5])
    assert result['boundary_capacitance_m']>result['capacitance_m']


@pytest.mark.parametrize('phi',[.0001,.01,.1,.99999999,1.,1.00000001,10.,100.])
def test_stable_capacitance_radius_bounds_scaling_and_boundary(libraries,phi):
    library=libraries['stable'];a=1e-5;c=phi*a
    result=library.evaluate(.05,258.15,101000.,a,c)
    scaled=library.evaluate(.05,258.15,101000.,3*a,3*c)
    assert min(a,c)*(1-1e-14)<=result['capacitance_m']<=max(a,c)*(1+1e-14)
    assert scaled['capacitance_m']==pytest.approx(3*result['capacitance_m'],rel=1e-14)
    assert result['boundary_capacitance_m']>result['capacitance_m']


@pytest.mark.parametrize('phi',[.01,.1,.3,1.,3.,10.,100.])
def test_equivalent_capacitance_preserves_native_results_away_from_unity(libraries,phi):
    args=(.05,258.15,101000.,1e-5,phi*1e-5)
    old=libraries['double'].evaluate(*args);stable=libraries['stable'].evaluate(*args)
    for field in ('capacitance_m','boundary_capacitance_m','alpha_a','alpha_c'):
        assert stable[field]==pytest.approx(old[field],rel=2e-12)


def test_plate_and_column_coefficients_are_axis_specific(libraries):
    library = libraries['double']
    plate = library.evaluate(.065, 273.15-14.3, 101000., 50e-6, 5e-6)
    column = library.evaluate(.03, 273.15-7, 101000., 5e-6, 50e-6, 1.)
    assert plate['alpha_a'] > plate['alpha_c']
    assert column['alpha_c'] > column['alpha_a']
    assert plate['decimal_precision'] == column['decimal_precision'] == 15
