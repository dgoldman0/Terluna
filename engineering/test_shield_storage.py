import numpy as np
from .shield_storage import dispatch,required_capacity


def test_separate_bus_energy_and_power_limits_and_conservation():
    t=np.arange(0.,21.);load=np.ones((len(t),2))*10;gen=np.zeros_like(load)
    out=dispatch(t,load,gen,np.array([100.,1000.]),np.array([100.,5.]),depth=1,charge=1,discharge=1)
    np.testing.assert_allclose(out['unserved_J'],[100,100])
    np.testing.assert_allclose(out['balance_residual_J'],0,atol=1e-10)
    assert out['first_unserved_s']==0


def test_recharge_loss_and_required_capacity_close():
    t=np.arange(11.);load=np.ones((11,2))*4;gen=np.zeros_like(load);gen[6:]=12
    needed=required_capacity(t,load,gen)
    out=dispatch(t,load,gen,needed,np.ones(2)*100)
    np.testing.assert_allclose(out['unserved_J'],0,atol=1e-10)
    np.testing.assert_allclose(out['balance_residual_J'],0,atol=1e-10)
    assert np.all(out['storage_loss_J']>0)
