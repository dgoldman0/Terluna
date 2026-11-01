"""Check geometric mapping and rate diagnostics without treating fits as laws."""
import json
import math

import pytest

from atmosphere.snow import mcsnow_shape_diagnosis as shape
from atmosphere.snow import mcsnow_pair_controls as pair


@pytest.mark.parametrize('temperature,age', [(-5.3,4.),(-5.3,12.),(-14.4,3.5),(-14.4,12.),(-20.1,8.)])
@pytest.mark.parametrize('mapping',shape.MAPPINGS)
def test_fit_rates_agree_with_finite_differences_inside_segments(temperature,age,mapping):
    h=1e-4
    middle=shape.fit_rates(temperature,age,mapping)
    left=shape.fit_state(temperature,age-h,mapping)
    right=shape.fit_state(temperature,age+h,mapping)
    mass_rate=(right['mass_kg']-left['mass_kg'])/(120*h)
    volume_rate=(right['volume_m3']-left['volume_m3'])/(120*h)
    a_log=(math.log(right['model_a_m'])-math.log(left['model_a_m']))/(120*h)
    c_log=(math.log(right['model_c_m'])-math.log(left['model_c_m']))/(120*h)
    assert mass_rate == pytest.approx(middle['mass_rate_kg_s'],rel=1e-8)
    assert mass_rate/volume_rate == pytest.approx(middle['deposition_density_kg_m3'],rel=1e-8)
    assert c_log/a_log == pytest.approx(middle['effective_growth_ratio'],rel=1e-8)


def test_envelope_mapping_preserves_mass_and_aspect_and_identifies_density_change():
    axes=shape.fit_rates(-14.4,10.,'axes')
    volume=shape.fit_rates(-14.4,10.,'volume')
    assert axes['mass_kg'] == volume['mass_kg']
    assert volume['model_c_m']/volume['model_a_m'] == pytest.approx(axes['model_c_m']/axes['model_a_m'])
    assert volume['volume_m3']/axes['volume_m3'] == pytest.approx(shape.transport.HEX_TO_SPHEROID)
    assert volume['effective_growth_ratio'] == axes['effective_growth_ratio']
    assert volume['deposition_density_kg_m3']*shape.transport.HEX_TO_SPHEROID == pytest.approx(axes['deposition_density_kg_m3'])


@pytest.mark.parametrize('temperature,join',[(-5.3,6.),(-14.4,4.),(-12.2,7.),(-16.5,5.)])
def test_independent_fit_joins_are_not_differentiated(temperature,join):
    with pytest.raises(ValueError,match='join'):
        shape.fit_rates(temperature,join,'axes')
    assert any(r['age_min']==join for r in shape.fit_joins(temperature))


def test_incompatible_fitted_initial_state_is_retained_as_rejection():
    entries=list(shape.cases(919.))
    rejected=[r for r in entries if r['temperature_C']==-16.5 and r['initialization']!='sphere']
    assert rejected and all(not r['admissible'] for r in rejected)
    assert all(r['initial_state']['density_kg_m3']>919. for r in rejected)
    # A valid hexagonal volume repairs four other mapping exclusions without changing source axes.
    for t in (-8.6,-12.2,-18.2,-22.):
        assert shape.fit_state(t,3.,'axes')['density_kg_m3']>919.
        assert shape.fit_state(t,3.,'volume')['density_kg_m3']<919.


def test_model_axes_are_converted_only_for_reporting():
    state=dict(mass_kg=1.,equatorial_dimension_m=2.,polar_dimension_m=3.,terminal_speed_m_s=4.)
    result=shape.ratios(state,-14.4,10.,2.)
    assert result['values']==dict(mass_kg=1.,equatorial_dimension_m=1.,polar_dimension_m=1.5,terminal_speed_m_s=4.)
    assert state['equatorial_dimension_m']==2.


def test_actual_fixed_state_source_separates_mean_rate_from_shape():
    product=json.loads(shape.PRODUCT.read_text())
    rows=product['fixed_state_probes']
    for changed in [r for r in rows if r['mode']=='plate_base']:
        original=next(r for r in rows if r['mode']=='native' and
                      (r['temperature_C'],r['age_min'],r['mapping'])==
                      (changed['temperature_C'],changed['age_min'],changed['mapping']))
        for key in ('effective_growth_factor','deposition_density_kg_m3','terminal_speed_m_s','capacitance_m'):
            assert changed['source_evaluation'][key] == original['source_evaluation'][key]
        rates=changed['required_uptake_multiplier']
        assert rates['a']/rates['c'] == pytest.approx(changed['ratio']['growth_ratio'],rel=2e-13)


def test_saved_comparisons_never_extrapolate_matching_mass():
    product=json.loads(shape.PRODUCT.read_text())
    for run in product['trajectories']:
        for comparison in run['comparisons']:
            matched=comparison['mass_matched']
            if matched is None:
                continue
            lo,hi=matched['interpolation_bracket_s']
            elapsed=(matched['model_age_min']-run['initial_age_min'])*60
            assert lo-1e-10<=elapsed<=hi+1e-10
            assert matched['ratio']['mass_kg']==pytest.approx(1.,abs=2e-14)


def test_pair_input_rejects_nonphysical_values():
    for values in [(0.,1.,1.),(1.,float('nan'),1.),(1.,1.,-1.)]:
        with pytest.raises(ValueError):
            pair.particle(*values)


def test_native_pair_symmetry_equal_speeds_and_sticking():
    product=json.loads(pair.PRODUCT.read_text())
    rows=product['cases']
    for row in rows:
        result=row['result']
        assert result['collection_efficiency']==pytest.approx(result['reverse_collection_efficiency'],rel=1e-12,abs=1e-15)
        if row['left']==row['right']:
            assert result['relative_speed_m_s']==0.
            assert result['kernel_m3_s']==0.
    for row in [r for r in rows if r['group']=='shape' and r['sticking']==1]:
        unity=rows[row['case_id']-1]
        assert unity['sticking']==0 and unity['left']==row['left'] and unity['right']==row['right']
        result=row['result']; reference=unity['result']
        assert result['collection_efficiency']==pytest.approx(reference['collection_efficiency']*result['sticking'],rel=2e-13,abs=1e-16)
        assert result['kernel_m3_s']==pytest.approx(reference['kernel_m3_s']*result['sticking'],rel=2e-13,abs=1e-30)


def test_independent_collision_implementation_compares_finite_outputs():
    product=json.loads(pair.PRODUCT.read_text())
    assert len(product['independent_parity'])==30
    for row in product['independent_parity']:
        for quantity in row['differences'].values():
            # Native default-REAL constants are retained in each implementation;
            # compare with an explicit small arithmetic tolerance, never a fit tolerance.
            assert quantity['independent']==pytest.approx(quantity['native'],rel=5e-12,abs=1e-25)


def test_shape_control_preserves_unmodified_arithmetic_and_early_history():
    from atmosphere.snow import mcsnow_shape_controls as control
    product=json.loads(control.PRODUCT.read_text())
    assert len(product['unchanged_control_parity'])==12
    assert all(r['exact_raw_equality'] and r['exact_state_equality'] for r in product['unchanged_control_parity'])
    for row in product['trajectories']:
        if row['shape_control']!=2 or row['initialization']!='sphere':
            continue
        unchanged=next(r for r in product['trajectories'] if r['shape_control']==0 and
                       r['initialization']=='sphere' and r['temperature_C']==row['temperature_C'] and r['dt_s']==row['dt_s'])
        # Fitted ratios have no authority before the first three-minute observation.
        assert row['comparisons'][0]['age_min']==3.
        assert row['comparisons'][0]['age_matched']==unchanged['comparisons'][0]['age_matched']


def test_fitted_directional_control_preserves_density_and_is_labelled_as_conditional():
    from atmosphere.snow import mcsnow_shape_controls as control
    product=json.loads(control.PRODUCT.read_text())
    assert 'not independent validation' in product['evidence']
    assert len(product['source_checks'])==36
    assert all(r['maximum_density_relative_residual']<5e-13 and
               r['maximum_prescribed_ratio_relative_residual']<5e-13 for r in product['source_checks'])
    driver=control.driver_source()
    assert max(map(len,driver.splitlines()))<=132
    assert 'initial_age_s+elapsed>=ratio_start_s' in driver
