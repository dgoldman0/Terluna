"""The acceptance decision must not hide regressions or missing evidence."""
import copy
import math

import numpy as np
import pytest

from atmosphere.snow import mcsnow_joint_acceptance as audit
from atmosphere.snow import mcsnow_joint_drag as drag
from atmosphere.snow import mcsnow_joint_growth as growth


@pytest.fixture(scope='module')
def baseline():
    return audit.baseline_inventory()


def exact_candidate(baseline):
    result=copy.deepcopy(baseline)
    for row in result['points']:
        row['value']=row['target']
    return result


def test_frozen_baseline_is_not_an_accepted_model(baseline):
    result=audit.compare(baseline,baseline)
    assert result['point_count']==1320
    assert result['above_ten_percent']==426
    assert result['regression_count']==0
    assert not result['accepted']


def test_joint_exact_comparison_passes_without_relaxing_numerical_checks(baseline):
    result=audit.compare(exact_candidate(baseline),baseline)
    assert result['accepted']
    assert result['above_ten_percent']==result['regression_count']==0


def test_improvement_of_every_other_point_cannot_hide_one_regression(baseline):
    candidate=exact_candidate(baseline)
    index=next(i for i,r in enumerate(baseline['points']) if 1e-3<abs(r['value']/r['target']-1)<.04)
    b=baseline['points'][index]
    candidate['points'][index]['value']=b['target']*(1+abs(b['value']/b['target']-1)+.001)
    result=audit.compare(candidate,baseline)
    assert result['above_ten_percent']==0
    assert result['regression_count']==1
    assert not result['accepted']


def test_large_improvement_above_ten_percent_is_still_rejected(baseline):
    candidate=exact_candidate(baseline)
    index=next(i for i,r in enumerate(baseline['points']) if abs(r['value']/r['target']-1)>.5)
    candidate['points'][index]['value']*=1.100001
    result=audit.compare(candidate,baseline)
    assert result['regression_count']==0
    assert result['above_ten_percent']==1
    assert not result['accepted']


@pytest.mark.parametrize('change',('missing','extra','duplicate','target','identity','nan','negative'))
def test_incomplete_or_altered_comparisons_fail_closed(baseline,change):
    candidate=exact_candidate(baseline)
    if change=='missing':candidate['points'].pop()
    elif change=='extra':
        row=copy.deepcopy(candidate['points'][0]);row['age_min']+=.001;row['id']=audit.point_id(row)
        candidate['points'].append(row)
    elif change=='duplicate':candidate['points'].append(candidate['points'][0])
    elif change=='target':candidate['points'][0]['target']*=1.0001
    elif change=='identity':candidate['points'][0]['age_min']+=.001
    elif change=='nan':candidate['points'][0]['value']=float('nan')
    else:candidate['points'][0]['value']=-1
    with pytest.raises(ValueError):audit.compare(candidate,baseline)


@pytest.mark.parametrize('field,value',(('mass_balance_max',1.01e-10),('volume_residual_max',1.01e-10),
    ('maximum_density_fraction',1.000000001),('volume_residual_max',float('nan'))))
def test_numerical_failures_prevent_acceptance(baseline,field,value):
    candidate=exact_candidate(baseline);candidate['invariants'][0][field]=value
    assert not audit.compare(candidate,baseline)['accepted']


def test_missing_numerical_checks_are_rejected(baseline):
    candidate=exact_candidate(baseline);candidate['invariants'].pop()
    with pytest.raises(ValueError):audit.compare(candidate,baseline)


def test_stored_ratios_cannot_override_actual_values():
    p=audit.load_baseline();w=next(f for f in p['searches'] if f['name']==p['preferred_search'])
    for r in w['assessment']['rows']:r['ratio']=[1.]*4
    assert audit.inventory(w,p['previous_cold'])==audit.baseline_inventory()


def test_conflicting_duplicate_values_are_rejected():
    p=audit.load_baseline();w=next(f for f in p['searches'] if f['name']==p['preferred_search'])
    w['assessment']['rows'][0]['values'][0]*=1.001
    with pytest.raises(ValueError,match='duplicate'):audit.inventory(w,p['previous_cold'])


def test_original_best_is_identical_for_oblate_area_ratios():
    for phi in (.01,.3,1.):
        for q in (.1,.5,1.):
            args=(phi,1e-9,q,1.2,1.7e-5,9.81)
            assert drag.reynolds(*args)==drag.reynolds(*args,original_best=True)


def test_prism_geometry_uses_bounded_projected_width():
    assert math.sqrt(3)/2<growth.STATE['prism_width_fraction']<1
    assert 3/math.pi==growth.STATE['prism_width_fraction']


def test_joint_controls_have_neutral_extension_and_reject_bad_values():
    for temp in (-5.3,-14.4):
        x=growth.baseline_coordinates(temp);v=growth.decode(x,growth.NEUTRAL)
        assert len(v)==20
        assert np.array_equal(v[12:],np.array([0,0,4,0,1,0,0,0]))
    for x in (np.zeros(12),np.full(13,np.nan),np.full(13,-.1),np.full(13,1.1)):
        with pytest.raises(ValueError):growth.decode(x,growth.STATE)


def test_full_contract_ages_fit_new_native_driver_capacity():
    for temp in (-5.3,-14.4):
        for e in growth.entries_for(temp):
            assert len(e['sample_ages_min'])<=256
            for age in e['sample_ages_min']:
                steps=(age-e['initial_age_min'])*60*128
                assert abs(steps-round(steps))<1e-7


def test_previous_density_family_cannot_meet_all_local_targets():
    limits=audit.structural_limits()
    assert not limits['cold_constant_density']['feasible']
    assert limits['cold_constant_density']['minimum_possible_maximum_relative_error']>.26
    assert not limits['warm_monotone_density']['nondecreasing_possible']
    assert not limits['warm_monotone_density']['nonincreasing_possible']


def test_committed_candidates_have_complete_contract_and_no_selected_failure():
    import json
    product=json.loads(growth.PRODUCT.read_text())
    assert product['baseline_sha256']==audit.BASELINE_SHA256
    assert product['accepted_candidates']==[]
    assert product['selected_model'] is None
    assert all(c['decision']['point_count']==1320 for c in product['candidates'])
    for candidate in product['candidates']:
        d=candidate['decision']
        assert not d['accepted']
        assert d['above_ten_percent'] or d['regression_count']
    assert all(c['passed'] for c in product['checks'])


def test_native_drag_and_independent_translation_match_in_saved_trace():
    import json
    p=json.loads(drag.PRODUCT.read_text())
    assert p['translation_maximum_relative_error']<1e-11
    assert all(r['translation_relative_error']<1e-11 for r in p['rows'])
    for r in p['rows']:
        if r['native_state']['aspect_ratio']<=1:
            assert r['speeds']['paper_prism_mean']==r['speeds']['native']
