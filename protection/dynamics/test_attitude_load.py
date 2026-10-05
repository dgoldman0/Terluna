import numpy as np
import pytest

from shared import constants as K
from .test_pattern_departure import inertial_environment
from .pattern_departure import sun_frame
from .attitude_load import rigid_square_load, disturbance_torque_bound


def test_one_axis_slew_charges_both_spinup_and_braking_force_couples():
    side=10000.; duration=3600.
    # Exact analytic command isolates the rigid-body accounting from the
    # separate 60 s rotation-spline approximation tested in pattern_departure.
    def command(t, derivative):
        f=np.clip((np.asarray(t)-1800.)/duration,0.,1.)
        result=np.zeros(f.shape+(3,))
        if derivative==1:
            result[...,1]=-(np.pi/2)/duration*30*f*f*(1-f)**2
        else:
            result[...,1]=-(np.pi/2)/duration**2*(60*f-180*f*f+120*f**3)
        return result
    times=np.arange(0.,7200.+1,2.)
    load=rigid_square_load(command,times,side)
    peak_rate=(np.pi/2)*1.875/duration
    expected=side*peak_rate/3
    assert np.trapezoid(load['force_per_mass'],times)==pytest.approx(expected,rel=3e-6)
    assert load['kinetic_energy_per_mass'].max()==pytest.approx(.5*(side**2/12)*peak_rate**2,rel=3e-5)
    assert load['torque_per_mass'][:,1].min()<0<load['torque_per_mass'][:,1].max()


def test_disturbance_envelope_keeps_gravity_and_both_sail_faces():
    env=inertial_environment();env.names=[];env.gm=[];env.central_fraction=.138
    state=np.array([[15e6,0.,0.,0.,0.,0.]])
    face=sun_frame(env,0.)
    gravity,radiation,couple=disturbance_torque_bound(env,0.,state,face)
    expected=(10000.**2/6)*2*K.MOON_GM/(15e6-10000./np.sqrt(2))**3
    assert gravity[0]==pytest.approx(expected)
    assert radiation[0]>0 and couple[0]>0
    opposite=face.copy();opposite[:,1:]*=-1
    np.testing.assert_allclose(disturbance_torque_bound(env,0.,state,opposite)[1],radiation)
    edge=face[:,[2,1,0]];edge[:,1]*=-1
    assert disturbance_torque_bound(env,0.,state,edge)[1][0]<radiation[0]*1e-4
