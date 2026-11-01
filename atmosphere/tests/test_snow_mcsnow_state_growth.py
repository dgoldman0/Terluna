"""Physical limits, selection boundaries and retained native controls."""
import json

import numpy as np
import pytest

from atmosphere.snow import mcsnow_state_growth as state
from atmosphere.snow import mcsnow_tuning as tuning
from atmosphere.tests.test_snow_mcsnow_tuning import synthetic_batch


def test_transition_limits_and_scale_invariance():
    assert state.stage_weight(1e-6, 0.) == 1.
    assert state.stage_weight(30e-6, 30.) == pytest.approx(.5)
    assert state.stage_weight(3e-6, 30.) < .001
    assert state.stage_weight(300e-6, 30.) > .999
    for radius in np.geomspace(1e-6,1e-3,12):
        assert state.stage_weight(radius,30.) == pytest.approx(state.stage_weight(radius*2,60.))


@pytest.mark.parametrize('r,t',[(0.,10.),(-1.,0.),(1.,-1.),(float('nan'),0.),(1.,float('inf'))])
def test_invalid_transition_geometry_is_rejected(r,t):
    with pytest.raises(ValueError): state.stage_weight(r,t)


def test_neutral_controls_preserve_existing_prescription():
    assert state.decode(state.BASE,'axis') == pytest.approx([1.,1.,1.,1.,1.,1.,0.,0.,0.])
    assert state.decode(state.BASE,'axis_density') == pytest.approx([1.,1.,1.,1.,1.,1.,0.,0.,1.])


def test_new_state_parameter_does_not_change_fit_data_selection():
    b=synthetic_batch()
    for j,(_,age) in enumerate(b.rows):
        if age not in tuning.TRAIN_AGES: b.values[j,:3] *= 1000
    b.values[:,3] *= 1000
    assert np.max(abs(tuning.residual(b,state.BASE))) == 0.


def test_native_parity_refinement_and_ledgers():
    p=json.loads(state.PRODUCT.read_text())
    assert all(c['passed'] for c in p['checks'])
    for name in ('original','constant'):
        c=next(c for c in p['checks'] if c['name']=='neutral_parity_'+name)
        assert c['maximum_relative_difference'] < 1e-10
    for a in p['assessments'].values():
        assert not a['failed_cases']
        assert a['maximum_density_fraction'] <= 1+1e-10


def test_state_search_keeps_both_starts_and_reserves_other_temperatures():
    p=json.loads(state.PRODUCT.read_text())
    for fit in p['searches']:
        spec=fit['specification']
        assert spec['training_ages'] == list(tuning.TRAIN_AGES)
        assert len(spec['cases']) == 4
        assert all('volume' not in name for name in spec['cases'])
    counts={k:v['states'] for k,v in p['assessments']['original']['scores'].items()}
    assert counts == {'training':12,'withheld_age':10,'geometry_sensitivity':10,'other_temperature':68}


def test_integrated_shape_budget_independent_of_growth_clock():
    for gamma in (.08,.2,.5,1.,3.,8.):
        ratio=1e4;phi0=.7
        phi=phi0*ratio**((gamma-1)/(gamma+2))
        assert state.equivalent_axis_ratio(phi0,phi,ratio) == pytest.approx(gamma)
    for row in state.shape_budget()['rows']:
        assert abs(row['reconstruction_relative_error']) < 1e-12


def test_temperature_variants_select_by_training_loss_and_keep_both_starts():
    p=json.loads(state.PRODUCT.read_text())
    assert len(p['temperature_searches']) == 4
    for fit in p['temperature_searches']:
        cases=fit['specification']['cases']
        assert len(cases) == 2
        assert any('sphere' in name for name in cases)
        assert any('axes' in name for name in cases)
        assert all('volume' not in name for name in cases)
        assert fit['seed_diameter_um'] in (0.,20.,40.)
        assert fit['physical_controls'][7] == fit['seed_diameter_um']
    for temperature,name in p['selected_temperature_searches'].items():
        choices=[f for f in p['temperature_searches'] if f['temperature_C']==float(temperature)]
        assert name == min(choices,key=lambda f:f['training_loss'])['name']
