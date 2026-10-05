"""Mass feedback preserves force physics and closes its stated scalar allocation."""
import numpy as np
import pytest

from .mass_feedback import allocated_mass_ratio,MassiveParallelPattern


def test_member_power_allocation_closes_fixed_plus_power_mass():
    prop=dict(peak_margin_factor=1.25,specific_power_W_kg=300.)
    specific=np.array([1.,10.]);ratio=allocated_mass_ratio(specific,prop)
    np.testing.assert_allclose(ratio,1.2+ratio*specific*1.25/300.)
    with pytest.raises(ValueError):allocated_mass_ratio([300.],prop)


def test_added_uniform_mass_changes_photon_acceleration_only(monkeypatch):
    gravity=np.array([[1.,2.,3.],[2.,3.,4.]])
    sail=np.array([[.1,.2,.3],[.2,.3,.4]])
    parent=MassiveParallelPattern.__mro__[1]
    monkeypatch.setattr(parent,'acceleration',lambda *a,**k:(gravity+sail,np.ones((8,2)),sail))
    model=MassiveParallelPattern(None,None,[1.2,1.5])
    force,light,new=model.acceleration(0.,np.zeros((2,6)),diagnostics=True)
    np.testing.assert_allclose(force,gravity+sail/np.array([1.2,1.5])[:,None])
    np.testing.assert_allclose(new,sail/np.array([1.2,1.5])[:,None])
    assert light.shape==(8,2)
