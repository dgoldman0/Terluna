import numpy as np
import pytest
from scipy.optimize._highspy._core import _Highs
from .service_search import smooth_step,cadence_account
from .service_constraints import projection,fit_free_controls


@pytest.fixture
def isolated_highs_scheduler():
    # Producers run in fresh single-thread processes. Earlier repository tests
    # may initialize HiGHS's process-global scheduler with another thread count.
    _Highs.resetGlobalScheduler(True)
    yield
    _Highs.resetGlobalScheduler(True)


def test_projection_jacobian_agrees_with_independent_finite_difference():
    q=np.array([[123000.,-80000.,15000000.],[50000.,40000.,14900000.]])
    source=np.array([695700000.,300000000.,1.5e11]);u=np.array([0.,0.,1.]);b=np.array([1.,0.,0.]);c=np.array([0.,1.,0.])
    _,jac,_=projection(q,source,u,b,c)
    for j in range(3):
        shift=np.eye(3)[j]*10
        difference=(projection(q+shift,source,u,b,c)[0]-projection(q-shift,source,u,b,c)[0])/20
        np.testing.assert_allclose(difference,jac[:,:,j],atol=1e-8)


def test_free_fit_enforces_energy_and_individual_power_without_endpoint(isolated_highs_scheduler):
    windows=np.array([[21600.,28800.]])
    rows=[(0,np.array([-1.,0.,0.]),-1.)]
    kwargs=dict(base=None,maps=None,windows=windows,mass=np.array([1e14]),caps=np.array([[.001]]),service=rows,specific=1.)
    control,info=fit_free_controls(**kwargs,total_budget_J=2e14,prefix_budget_J=2e14)
    assert info['success'],info
    np.testing.assert_allclose(control[0,0],[1.,0.,0.],atol=1e-8)
    control,info=fit_free_controls(**kwargs,total_budget_J=.5e14,prefix_budget_J=2e14)
    assert control is None and info['status']==2,info
    kwargs['caps']=np.array([[.0001]])
    control,info=fit_free_controls(**kwargs,total_budget_J=2e14,prefix_budget_J=2e14)
    assert control is None and info['status']==2,info


def test_smooth_turn_and_inventory_penalty():
    np.testing.assert_allclose(smooth_step(np.array([-1,0,.5,1,2])),[0,0,.5,1,1])
    two=cadence_account(48*3600,mass_kg=100.);ten=cadence_account(240*3600,mass_kg=100.)
    assert two['time_capacity_group_lower_bound']==8
    assert ten['time_capacity_group_lower_bound']==40
    assert ten['time_capacity_mass_lower_bound_kg']==5*two['time_capacity_mass_lower_bound_kg']
    assert not ten['spatial_coverage_and_handover_demonstrated']
