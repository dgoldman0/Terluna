"""Limits and selection controls for independent early/later growth rates."""
import json

import numpy as np
import pytest

from atmosphere.snow import mcsnow_growth_limits as limits
from atmosphere.snow import mcsnow_state_growth as state


def test_interpolation_preserves_constant_and_parent_limits():
    for radius in np.geomspace(1e-8,1e-2,20):
        w=state.stage_weight(radius,35.)
        assert limits.axis_multiplier(radius,35.,1.,2.1)==pytest.approx(2.1**w)
        assert limits.axis_multiplier(radius,35.,.7,.7)==pytest.approx(.7)
        assert limits.axis_multiplier(radius,0.,.2,1.8)==pytest.approx(1.8)
        assert .52 <= limits.axis_multiplier(radius,35.,.52,1.9) <= 1.9
        assert 2.1 <= limits.axis_multiplier(radius,50.,3.5,2.1) <= 3.5


def test_two_limits_allow_opposite_evolution_with_size():
    assert limits.axis_multiplier(1e-8,35.,.52,1.9)==pytest.approx(.52)
    assert limits.axis_multiplier(1e-2,35.,.52,1.9)==pytest.approx(1.9)
    assert limits.axis_multiplier(1e-8,50.,3.5,2.1)==pytest.approx(3.5)
    assert limits.axis_multiplier(1e-2,50.,3.5,2.1)==pytest.approx(2.1)


@pytest.mark.parametrize('bad',[[0]*7,[0]*9,[float('nan')]*8,[-.1]*8,[1.1]*8])
def test_invalid_coordinates_rejected(bad):
    with pytest.raises(ValueError):limits.decode(bad)


def test_neutral_early_control_and_input_order():
    x=np.r_[state.BASE,limits.early_coordinate(1.)]
    assert limits.decode(x)==pytest.approx([1.,1.,1.,1.,1.,1.,0.,0.,0.,1.])
    for value in limits.EARLY_BOUNDS:
        x[-1]=limits.early_coordinate(value)
        assert limits.decode(x)[-1]==pytest.approx(value)


def test_native_neutral_parity_refinement_and_ledgers():
    p=json.loads(limits.PRODUCT.read_text())
    assert p['checks'] and all(c['passed'] for c in p['checks'])
    assert p['parity']['original']['states']==114
    assert p['parity']['constant']['states']==114
    for fit in p['searches']:
        if fit['name'] in p['selected_searches']:
            assert not fit['assessment']['failed_cases']
            assert fit['assessment']['maximum_density_fraction']<=1+1e-10


def test_selection_uses_only_each_declared_training_group():
    p=json.loads(limits.PRODUCT.read_text())
    for name in p['selected_searches']:
        f=next(f for f in p['searches'] if f['name']==name)
        group=[r for r in p['searches'] if limits.group_key(r)==limits.group_key(f)]
        assert f['training_loss']==min(r['training_loss'] for r in group)
        assert f['training_loss']<=f.get('baseline_loss',f['previous_training_loss'])+1e-8
        for i,value in f['specification'].get('fixed',{}).items():
            assert f['coordinates'][int(i)]==value
        if source_name:=f['specification'].get('refined_from'):
            source=next(r for r in p['searches'] if r['name']==source_name)
            assert limits.fit_options(f)==limits.fit_options(source)
            assert f['training_loss']<=source['training_loss']+1e-10
            assert all(source['coordinates'][int(i)]==value for i,value in f['specification'].get('fixed',{}).items())
        assert len(f['specification']['cases'])==2
        assert all('volume' not in name for name in f['specification']['cases'])
        assert f['specification']['training_ages']==[5.,10.,20.]


def test_extended_bounds_preserve_physical_settings_under_coordinate_change():
    base=np.r_[state.BASE,limits.early_coordinate(.25)]
    extended=base.copy();extended[5]*=.5
    extended[7]=limits.early_coordinate(.25,(.0625,4.))
    assert limits.decode(extended,early_bounds=(.0625,4.),ventilation_max=4.)==pytest.approx(limits.decode(base))
