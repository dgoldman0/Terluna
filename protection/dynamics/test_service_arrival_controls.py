import numpy as np
import pytest
from scipy.optimize._highspy._core import _Highs
from .service_arrival_controls import arc_fraction,fit_service_arcs


@pytest.fixture
def isolated_highs_scheduler():
    # Match the fresh process used by each resource-bounded producer.
    _Highs.resetGlobalScheduler(True)
    yield
    _Highs.resetGlobalScheduler(True)


def test_prefix_energy_includes_burns_during_initial_service(isolated_highs_scheduler):
    windows=np.array([[0.,21600.],[21600.,43200.],[100000.,120000.]])
    np.testing.assert_allclose(arc_fraction(43200,windows),[1,1,0],atol=1e-14)
    np.testing.assert_allclose(arc_fraction(10800,windows),[.5,0,0],atol=1e-14)
    # The first service burn must satisfy the twelve-hour energy budget even
    # though the full sequence budget would allow it.
    specific=30000/(1.4*np.cos(np.pi/4));mass=np.array([1e14/specific]);gradient=np.zeros(9);gradient[0]=-1
    args=dict(base=None,maps=None,windows=windows,mass=mass,caps=np.ones((1,3))*.001,
        service=[(0,gradient,-1.)],cuts=[],total_budget_J=2e14,prefix_budget_J=.5e14)
    control,info=fit_service_arcs(**args)
    assert control is None and info['status']==2,info
    args['prefix_budget_J']=2e14;control,info=fit_service_arcs(**args)
    assert info['success'],info
    assert abs(info['prefix_translation_euclidean_J']/1e14-1)<1e-8
