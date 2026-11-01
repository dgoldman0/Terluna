"""Physical limits and reproducibility of the evolving column-density diagnostic."""
import json

import numpy as np
import pytest

from atmosphere.snow import mcsnow_column_density as column
from atmosphere.snow import mcsnow_growth_limits as limits


def test_density_transition_has_neutral_and_constant_limits():
    for radius in np.geomspace(1e-8,1e-2,20):
        assert column.density_multiplier(radius,80.,1.,.7)==pytest.approx(.7)
        assert column.density_multiplier(radius,0.,.4,.7)==pytest.approx(.7)
        value=column.density_multiplier(radius,80.,.4,.7)
        assert .28<=value<=.7
    assert column.density_multiplier(1e-8,80.,.4,.7)==pytest.approx(.28)
    assert column.density_multiplier(1e-2,80.,.4,.7)==pytest.approx(.7)


def test_density_transition_can_increase_or_decrease_without_changing_sign():
    radii=np.geomspace(1e-7,1e-3,20)
    increasing=[column.density_multiplier(r,80.,.4,.7) for r in radii]
    decreasing=[column.density_multiplier(r,80.,2.,.3) for r in radii]
    assert all(a<=b for a,b in zip(increasing,increasing[1:]))
    assert all(a>=b for a,b in zip(decreasing,decreasing[1:]))
    assert min(increasing+decreasing)>0


@pytest.mark.parametrize('bad',[[0]*9,[0]*12,[float('nan')]*10,[-.1]*10,[1.1]*11])
def test_invalid_coordinates_rejected(bad):
    with pytest.raises(ValueError):column.decode(bad)


def test_new_coordinates_leave_parent_physical_values_intact():
    x=np.array([.5,.4,.1,.3,.6,.7,.2,.1])
    parent=limits.decode(x,seed_diameter_um=40.)
    assert column.decode(np.r_[x,0.,.5])[:10]==pytest.approx(parent)
    assert column.decode(np.r_[x,.25,.5])[-2:]==pytest.approx([75.,1.])
    released=column.decode(np.r_[x,.25,.5,.75])
    assert released[7]==pytest.approx(60.)
    assert np.delete(released,7)==pytest.approx(np.delete(column.decode(np.r_[x,.25,.5]),7))


def test_saved_native_controls_refine_and_preserve_the_cold_fit():
    p=json.loads(column.PRODUCT.read_text())
    assert p['checks'] and all(c['passed'] for c in p['checks'])
    for f in p['searches']:
        if f['name'] not in p['selected_searches']:continue
        assert f['cold_control']['maximum_relative_difference']==0.
        assert not f['assessment']['failed_cases']
        assert f['assessment']['maximum_density_fraction']<=1+1e-10
        assert len(f['cold_control']['assessment']['rows'])==18
        for name,value in column.comparison_metrics(f['assessment'],f['specification']['fitting_ages_min']).items():
            assert f['metrics'][name]==value


def test_selection_retains_the_same_objective_and_separate_families():
    p=json.loads(column.PRODUCT.read_text())
    assert p['training']['ages_min']==[5.,10.,20.]
    assert p['training']['initializations']==['sphere','axes']
    for name in p['selected_searches']:
        f=next(f for f in p['searches'] if f['name']==name)
        group=[g for g in p['searches'] if len(g['coordinates'])==len(f['coordinates'])
               and g['specification']['density_early_ratio_bounds']==f['specification']['density_early_ratio_bounds']
               and g['specification']['fitting_ages_min']==f['specification']['fitting_ages_min']
               and g['specification']['fixed']==f['specification']['fixed']
               and g['specification'].get('age_weights',{})==f['specification'].get('age_weights',{})]
        assert f['training_loss']==min(g['training_loss'] for g in group)
    for f in p['searches']:
        assert set(f['specification']['fitting_ages_min'])<=set([3.,4.,5.,5.75,column.JOIN_TRAIN_AGE,10.,20.])
        assert 15. not in f['specification']['fitting_ages_min'] and 30. not in f['specification']['fitting_ages_min']
        assert all(f['coordinates'][int(i)]==v for i,v in f['specification']['fixed'].items())
        assert len(f['coordinates'])==(11 if f['specification'].get('compact_limit_fitted') else 10)


def test_sampled_early_history_meets_the_same_full_trajectory_endpoint():
    p=json.loads(column.PRODUCT.read_text())
    for f in p['searches']:
        if f['name'] not in p['selected_searches']:continue
        rows=f['early_model_history']['rows']
        assert rows[0]['values'][1:3]==pytest.approx([4e-6,4e-6])
        assert rows[0]['aspect_ratio']==pytest.approx(1.)
        endpoint=next(r for r in f['assessment']['rows'] if r['initialization']=='sphere' and r['age_min']==3.)
        assert rows[-1]['values']==endpoint['values']


def test_density_bound_change_can_preserve_the_physical_controls():
    x=np.array([.5,.4,.1,.3,.6,.7,.2,.1,.4,.1,.7])
    expected=column.decode(x)
    extended=x.copy()
    extended[9]=limits.early_coordinate(expected[-1],(.0625,4.))
    assert column.decode(extended,density_bounds=(.0625,4.))==pytest.approx(expected)


@pytest.mark.parametrize('bounds',[(0.,4.),(1.,1.),(2.,4.),(.25,float('nan'))])
def test_invalid_density_bounds_rejected(bounds):
    with pytest.raises(ValueError):column.decode(np.full(10,.5),density_bounds=bounds)


def test_seed_constraint_and_extra_age_labels_preserve_the_objective_boundary():
    p=json.loads(column.PRODUCT.read_text())
    for f in p['searches']:
        ages=f['specification']['fitting_ages_min']
        for r in f['assessment']['rows']:
            if r['initialization']=='volume':continue
            if r['initialization']=='axes' and r['age_min']==3.:
                assert r['group']=='initial_state'
            else:assert (r['group']=='training')==(r['age_min'] in ages)
        if 'additional_age_assessment' in f:
            rows=f['additional_age_assessment']['rows']
            assert len(rows)==2*len(set([4.,6.,8.,12.,25.])-set(ages))
            assert all(r['age_min'] not in ages for r in rows)
            assert f['metrics']['maximum_additional_age_relative_error']==max(abs(v-1) for r in rows for v in r['ratio'][:3])


def test_saved_objectives_match_independent_case_weighted_calculation():
    p=json.loads(column.PRODUCT.read_text())
    for f in p['searches']:
        spec=f['specification'];weights={float(age):weight for age,weight in spec.get('age_weights',{}).items()}
        loss=0.
        for start in ('sphere','axes'):
            rows=[r for r in f['assessment']['rows'] if r['initialization']==start and r['group']=='training']
            denominator=6*sum(weights.get(r['age_min'],1.) for r in rows)
            loss+=sum(weights.get(r['age_min'],1.)*sum(np.log(r['ratio'][:3])**2) for r in rows)/denominator
        assert loss==pytest.approx(f['fine_training_loss'],rel=1e-12,abs=1e-15)


def test_preferred_fit_meets_sampled_target_without_a_cold_regression():
    p=json.loads(column.PRODUCT.read_text());f=next(f for f in p['searches'] if f['name']==p['preferred_search'])
    rows=f['interval_assessment']['rows']+f['assessment']['rows']
    rows=[r for r in rows if r['initialization'] in ('sphere','axes')]
    assert max(abs(value-1) for r in rows for value in r['ratio'][:3])<.1
    assert f['cold_control']['maximum_relative_difference']==0.
    assert p['target_assessment']['attained']
    assert len(f['interval_assessment']['rows'])==222
