import numpy as np
import pytest
from atmosphere.snow.column_diagnosis_analysis import allocation_gap,quantile_bracket


def test_constant_demand_has_no_allocation_gap():
    demand=np.array([.1,1.,3.]);supply=np.ones(3)
    np.testing.assert_array_equal(allocation_gap(demand,np.minimum(demand,supply),supply),0)


def test_alternating_empty_and_limited_layer_has_known_gap():
    # Unit duration: empty for half, demand 1 kg/s for half; cap 0.5 kg/s.
    assert allocation_gap(.5,.25,.5)==.25


@pytest.mark.parametrize('args',[(-1,0,1),(1,2,1),(np.nan,0,1),(1,1,-1)])
def test_invalid_water_account_is_rejected(args):
    with pytest.raises(ValueError):allocation_gap(*args)


def test_mass_cdf_quantile_is_bracketed_without_interpolation():
    assert quantile_bracket([0,.9,1],[.001,.01,None],.99)==[.01,None]
    assert quantile_bracket([0,.9,.995,1],[.001,.01,.1,None],.99)==[.01,.1]
    with pytest.raises(ValueError):quantile_bracket([.8,.4,1],[.001,.01,None],.99)
