"""The focused sensitivity changes only the explicitly selected capacitance."""
import numpy as np
import pytest
from atmosphere.snow import welss_column as w
from atmosphere.snow import welss_compact_capacitance as c


def test_compact_intervention_changes_only_capacitance_below_threshold():
    m=w.M_THRESHOLD*np.array([.01,.999999,1.,1.000001,100.])
    original=w.geometry(m);changed=c.compact_geometry(m,w.geometry)
    for key in original:
        if key!='C':
            np.testing.assert_array_equal(changed[key],original[key])
    np.testing.assert_array_equal(changed['C'][:2],original['D'][:2]/2)
    np.testing.assert_array_equal(changed['C'][2:],original['C'][2:])


def test_runtime_intervention_restores_original_even_after_failure(monkeypatch):
    original=w.geometry;original_solver=w.solve_ivp
    def fail(**kwargs):
        mass=w.M_THRESHOLD/100
        assert w.geometry(mass)['C']==w.geometry(mass)['D']/2
        raise RuntimeError('intentional fixture failure')
    monkeypatch.setattr(w,'integrate',fail)
    with pytest.raises(RuntimeError,match='intentional'):
        c.integrate_with_intervention()
    assert w.geometry is original
    assert w.solve_ivp is original_solver


def test_startup_step_bounds_trial_log_mass_change_without_changing_rhs():
    rhs=lambda z,y:np.array([1e6,2e6,-1.,-2.])
    state=np.array([-30.,-29.,0.,0.])
    def capture(fun,span,y,**kw):
        np.testing.assert_array_equal(fun(span[0],y),rhs(span[0],y))
        np.testing.assert_array_equal(y,state)
        assert kw['first_step']==5e-10
        return kw
    c.startup_limited_solver(capture,rhs,(3000.,0.),state,rtol=1e-7)


def test_invalid_trial_is_rejected_and_known_growth_solution_recovered():
    # A deliberately excessive initial step produces invalid internal RK stages.
    # The accepted endpoint must recover the analytic continuous-growth law.
    from scipy.integrate import solve_ivp
    a=1e-6;q=2/3;y0=np.array([-30.,0.]);counter=dict(invalid_rk_stages=0)
    rhs=lambda t,y:np.array([a*np.exp(-q*y[0]),1.])
    checked=c.domain_checked_rhs(rhs,1,counter)
    result=solve_ivp(checked,(0.,1.),y0,method='DOP853',first_step=1.,rtol=1e-10,atol=1e-11)
    assert result.success and counter['invalid_rk_stages']>0
    expected=np.log((np.exp(y0[0])**q+q*a)**(1/q))
    assert result.y[0,-1]==pytest.approx(expected,rel=1e-10)
    assert result.y[1,-1]==pytest.approx(1.)
    assert np.isfinite(result.y).all()
