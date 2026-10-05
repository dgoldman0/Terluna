"""Physical boundary and partial-burn checks for return freedoms."""
import numpy as np
from scipy.spatial.transform import Rotation, RotationSpline

from .return_escape import rolled_return, applied_impulse, quintic_fraction
from .attitude_load import rigid_square_load


def test_smooth_in_plane_roll_preserves_normal_and_recovers_orientation_rate():
    start,end=43200.,166607.
    base=RotationSpline([start-1000.,end+1000.],Rotation.from_rotvec([[.3,.1,0],[.3,.1,0]]))
    command=rolled_return(base,start,end,60.,5400.)
    t=np.linspace(start,end,2001)
    np.testing.assert_allclose(command(t).as_matrix()[:,:,2],base(t).as_matrix()[:,:,2],atol=2e-13)
    np.testing.assert_allclose(command([start,end]).as_matrix(),base([start,end]).as_matrix(),atol=2e-13)
    assert np.linalg.norm(command(start,1))<1e-8
    assert np.linalg.norm(command(end,1))<1e-8
    assert np.max(np.linalg.norm(command(t,1),axis=1))<np.deg2rad(.1)
    dense=np.arange(start,start+5401.,1.)
    nominal=np.trapezoid(rigid_square_load(command,dense)['force_per_mass'],dense)
    # z-axis I/m=L^2/6; total |delta omega| = 2*1.875*angle/duration.
    expected=(10000./3)*(3.75*np.deg2rad(60)/5400.)
    np.testing.assert_allclose(nominal,expected,rtol=1e-4)


def test_partial_burn_integrals_are_additive_and_zero_outside_windows():
    u=np.array([[[3.,4.,0.],[0.,0.,2.]]]);windows=np.array([[10.,20.],[30.,40.]])
    np.testing.assert_allclose(applied_impulse(u,windows,0.,50.),[7.])
    np.testing.assert_allclose(applied_impulse(u,windows,10.,15.),[2.5])
    np.testing.assert_allclose(applied_impulse(u,windows,20.,30.),[0.])
    total=sum(applied_impulse(u,windows,a,b) for a,b in [(0.,13.),(13.,34.),(34.,50.)])
    np.testing.assert_allclose(total,[7.])
    np.testing.assert_allclose(quintic_fraction([-1.,0.,.5,1.,2.],0.,1.),[0.,0.,.5,1.,1.])
