"""Conservation and controlled transfer checks; no lunar empirical validation."""
from copy import deepcopy
import json
import math

import pytest

from atmosphere.snow import mcsnow_lunar_growth as lunar
from shared.constants import MOON_SURFACE_GRAVITY


@pytest.fixture(scope='module')
def product():
    return json.loads(lunar.PRODUCT.read_text())


def output_fixture(product, name='fall_main_moon_deep_plate_base_1_dt0.5'):
    run = next(r for r in product['run_manifest'] if r['name']==name)
    case = next(r for r in product['falling_cases'] if r['name']==name)
    return deepcopy((case['points'],[case['summary']],
        [dict(native_gravity_m_s2=MOON_SURFACE_GRAVITY)],[],run['configuration'],MOON_SURFACE_GRAVITY))


def test_supply_budget_matches_independent_flux_balance():
    # A column with 12 kg/s crossing 3 m2 and 2 crystals/s/m2 allocates
    # 2 kg per crystal, regardless of residence time or path segmentation.
    water_flux, number_flux, depth = 12/3, 2., 100.
    assert lunar.supply_budget(water_flux,number_flux,depth,depth)==2.
    intervals = (0.,5.,31.,33.,100.)
    increments = [lunar.supply_budget(water_flux,number_flux,depth,b-a)
                  for a,b in zip(intervals,intervals[1:])]
    assert sum(increments)==pytest.approx(2.)
    assert lunar.supply_budget(water_flux,number_flux,depth,25.)==.5


@pytest.mark.parametrize('values',[
    (-1.,2.,100.,10.),(1.,0.,100.,10.),(1.,2.,0.,0.),
    (1.,2.,100.,101.),(1.,2.,100.,-1.),(math.nan,2.,100.,10.)])
def test_invalid_supply_account_rejected(values):
    with pytest.raises(ValueError):
        lunar.supply_budget(*values)


def test_actual_column_closes_geometry_and_deposition_accounts(product):
    lunar.validate_output(*output_fixture(product))


@pytest.mark.parametrize('damage',('volume','mass','allocation','time','gravity','exit'))
def test_inconsistent_physical_output_rejected(product,damage):
    rows,summary,air,steps,config,gravity = output_fixture(product)
    if damage=='volume':
        rows[-1]['polar_dimension_m']*=1.01
    elif damage=='mass':
        # A forged zero residual must not hide missing deposited mass.
        rows[-1]['integrated_deposition_kg']*=.9
        rows[-1]['ice_mass_balance_residual_kg']=0.
    elif damage=='allocation':
        rows[-1]['allocated_water_kg']*=1.01
    elif damage=='time':
        rows[-1]['elapsed_s']+=1.
    elif damage=='gravity':
        air[0]['native_gravity_m_s2']=9.81
    else:
        rows[-1]['height_above_freezing_m']=1.
    with pytest.raises(ValueError):
        lunar.validate_output(rows,summary,air,steps,config,gravity)


def test_balanced_mass_still_cannot_exceed_the_water_allowance(product):
    rows,summary,air,steps,config,gravity = output_fixture(product)
    last=rows[-1]
    last['integrated_deposition_kg']=last['allocated_water_kg']*1.1
    last['mass_kg']=rows[0]['mass_kg']+last['integrated_deposition_kg']
    last['ice_mass_balance_residual_kg']=0.
    last['density_kg_m3']=last['mass_kg']/last['ice_volume_m3']
    with pytest.raises(ValueError,match='integrated water allowance'):
        lunar.validate_output(rows,summary,air,steps,config,gravity)


def source_step(rows, config):
    last=rows[-1]
    dt=.02
    cap=config['water_supply_kg_m2_s']/config['cloud_depth_m']*last['terminal_speed_m_s']/config['number_flux_m2_s']
    return dict(evaluation_elapsed_s=last['elapsed_s']-dt,state_elapsed_s=last['elapsed_s'],
        step_s=dt,fall_speed_m_s=last['terminal_speed_m_s'],
        supply_rate_kg_s=cap,requested_rate_kg_s=2*cap,admitted_rate_kg_s=cap)


def test_source_gate_limits_the_admitted_rate(product):
    rows,summary,air,steps,config,gravity=output_fixture(product)
    step=source_step(rows,config)
    lunar.validate_output(rows,summary,air,[step],config,gravity)
    step['admitted_rate_kg_s']*=1.001
    with pytest.raises(ValueError,match='declared gate'):
        lunar.validate_output(rows,summary,air,[step],config,gravity)


def test_subtractive_mass_resolution_is_not_a_relative_rate_error(product):
    rows,summary,air,steps,config,gravity=output_fixture(product)
    step=source_step(rows,config)
    quantum=math.ulp(rows[-1]['mass_kg'])/step['step_s']
    # A sub-ulp mass increment rounds to one ulp, a large relative rate error
    # that remains consistent with floating-point mass storage.
    step.update(requested_rate_kg_s=.7*quantum,admitted_rate_kg_s=quantum)
    lunar.validate_output(rows,summary,air,[step],config,gravity)
    step['admitted_rate_kg_s']=10*quantum
    with pytest.raises(ValueError,match='declared gate'):
        lunar.validate_output(rows,summary,air,[step],config,gravity)


def test_source_allowance_uses_local_number_flux(product):
    rows,summary,air,steps,config,gravity=output_fixture(product)
    step=source_step(rows,config)
    step['supply_rate_kg_s']*=2
    step['admitted_rate_kg_s']=step['supply_rate_kg_s']
    with pytest.raises(ValueError,match='number-flux account'):
        lunar.validate_output(rows,summary,air,[step],config,gravity)


def test_earth_adapter_parity_and_analytic_settling(product):
    assert len(product['earth_parity'])==2
    assert all(r['exact_equality'] for r in product['earth_parity'])
    assert len(product['drift_controls'])==2
    for r in product['drift_controls']:
        assert r['actual_exit_s']==pytest.approx(r['expected_exit_s'],rel=1e-12)
        assert r['unchanged_mass']
    assert product['profile']['gravity_m_s2']==MOON_SURFACE_GRAVITY
    assert product['profile']['earth_control_gravity_m_s2']==9.81


def test_censored_runs_keep_partial_evidence_and_missing_endpoints(product):
    manifest=product['run_manifest']
    assert len(manifest)==68
    assert all(r['complete'] and r['returncode']==0 and not r['validation_errors'] for r in manifest)
    censored={r['name'] for r in manifest if r['termination']=='diameter_censored'}
    assert len(censored)==6
    for case in product['fixed_cases']:
        if case['name'] in censored:
            assert case['last_state']['diameter_m']>=1.
            assert case['points'][-1]['age_s']==1800
            assert case['points'][-1]['state'] is None
    for case in product['falling_cases']:
        if case['name'] in censored:
            assert case['final_state']['height_above_freezing_m']>0
            assert case['final_state']['diameter_m']>=1.


def test_all_supplied_paths_obey_the_integrated_flux_bound(product):
    for case in product['falling_cases']:
        if case['water_supply_kg_m2_s']<0:
            continue
        initial=case['points'][0]
        ceiling=case['water_supply_kg_m2_s']/case['number_flux_m2_s']
        for row in case['points']:
            fraction=1-row['height_above_freezing_m']/initial['height_above_freezing_m']
            assert row['mass_kg']-initial['mass_kg']<=ceiling*fraction+1e-12*ceiling


def test_transfer_effects_are_separated_and_refined(product):
    fixed={r['name']:r for r in product['fixed_cases']}
    get=lambda world,p:fixed[f'fixed_-14.4_{world}_{p}_plate_base_dt0.0078125']['points'][-1]['state']['diameter_m']
    # Independent pressure and gravity changes act in opposite directions.
    assert get('earth','lunar')>get('earth','lab')
    assert get('moon','lab')<get('earth','lab')
    assert get('moon','lunar')>get('earth','lab')
    assert len(product['refinement'])==32
    for comparison in product['refinement']:
        assert comparison['same_termination']
        assert max(abs(v) for v in comparison['relative_change'].values())<.004
