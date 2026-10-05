import numpy as np
import pytest
from shared import constants as K
from protection.magnetic_fields import loop_field, basis, field, loop_samples, regional_loops, mutual_load, mutual_inductance


def test_loop_matches_independent_polygon_and_axis():
    radius = 120.; current = 700.
    phi = np.arange(4096)*2*np.pi/4096
    vertices = np.c_[radius*np.cos(phi), radius*np.sin(phi), np.zeros(len(phi))]
    p = np.array([[25., 31., 90.], [210., -45., 27.], [0., 0., 230.]])
    exact = loop_field(p, radius, current)
    delta = np.roll(vertices, -1, axis=0)-vertices
    midpoint = vertices+delta/2
    r = p[:, None, :]-midpoint[None, :, :]
    numerical = K.VACUUM_PERMEABILITY*current/(4*np.pi)*np.sum(np.cross(delta, r)/np.linalg.norm(r, axis=2)[:, :, None]**3, axis=1)
    np.testing.assert_allclose(numerical, exact, rtol=2e-6, atol=1e-14)
    assert exact[-1, 2] == pytest.approx(K.VACUUM_PERMEABILITY*current*radius**2/(2*(radius**2+230**2)**1.5))


def test_translated_rotated_field_and_vacuum_maxwell():
    p = np.array([[350., 90., 60.], [110., 260., 10.]])
    q = basis([1., 2., 3.]); c = np.array([700., -120., 90.])
    np.testing.assert_allclose(loop_field(p@q.T+c, 100., centre=c, normal=q[:, 2]), loop_field(p, 100.)@q.T, rtol=1e-12, atol=1e-20)
    h = .01; jac = np.array([(loop_field(p+np.eye(3)[i]*h, 100.)-loop_field(p-np.eye(3)[i]*h, 100.))/(2*h) for i in range(3)])
    assert abs(np.einsum('ini->n', jac)).max() < 1e-18
    np.testing.assert_allclose(jac, jac.transpose(2, 1, 0), atol=1e-18)


def test_regional_circuits_on_sphere_and_total_moment():
    for count in [2, 4]:
        loops = regional_loops(500e3, 1.5e21, count)
        for loop in loops:
            p, _ = loop_samples(loop)
            np.testing.assert_allclose(np.linalg.norm(p, axis=1), K.MOON_RADIUS, atol=1e-8)
        moment = sum(l['current']*np.pi*l['radius']**2*l['normal'] for l in loops)
        np.testing.assert_allclose(moment, [0, 0, 1.5e21], atol=3e5)


def test_mutual_force_and_angular_momentum_conservation():
    a = dict(radius=100., current=1500., centre=np.array([0., 0., 0.]), normal=np.array([0., 0., 1.]))
    b = dict(radius=70., current=-800., centre=np.array([350., 20., 160.]), normal=np.array([.3, .2, 1.]))
    ab = mutual_load(a, b); ba = mutual_load(b, a)
    np.testing.assert_allclose(ab['force_N']+ba['force_N'], 0, atol=1e-12)
    total = ab['torque_N_m']+ba['torque_N_m']+np.cross(b['centre']-a['centre'], ab['force_N'])
    np.testing.assert_allclose(total, 0, atol=1e-10)
    step = .01; direction = np.array([1., 0., 0.])
    derivative = (mutual_inductance(a, {**b, 'centre':b['centre']+step*direction})-
                  mutual_inductance(a, {**b, 'centre':b['centre']-step*direction}))/(2*step)
    assert ab['force_N'][0] == pytest.approx(a['current']*b['current']*derivative, rel=1e-7)
