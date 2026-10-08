"""Momentum conservation, finite-Sun geometry and independent dynamics checks."""
from pathlib import Path
import numpy as np
import pytest
from shared import constants as K
from protection.dynamics.ephemeris import Ephemeris, DEFAULT_KERNEL, verify_kernel, InterpolatedEphemeris
from protection.dynamics.model import normalized_derivatives, geometry, target, relative_gravity, moon_acceleration_residual
from protection.dynamics.optical import disk_overlap, solar_visibility, optical_projection, covering_radius, length, light_time_allowance, arriving_ray_vectors


def test_disk_overlap_contacts_and_containment():
    assert disk_overlap(1., 1., 2.) == 0
    assert disk_overlap(1., 2., 0.) == pytest.approx(np.pi)
    assert disk_overlap(1., 1., 1.) == pytest.approx(2*np.pi/3-np.sqrt(3)/2)
    assert disk_overlap(1., 1., 0.) == pytest.approx(np.pi)


def test_finite_sun_total_annular_and_clear_eclipses():
    sun = np.array([K.AU, 0, 0])
    assert solar_visibility(sun, np.array([1e8, 0, 0])) == pytest.approx(0)
    assert solar_visibility(sun, np.array([-1e8, 0, 0])) == 1
    assert solar_visibility(sun, np.array([1e8, 2e7, 0])) == 1
    far = np.array([2e9, 0, 0])
    rs = np.arcsin(K.SUN_RADIUS/K.AU); re = np.arcsin(K.EARTH_RADIUS/length(far))
    assert solar_visibility(sun, far) == pytest.approx(1-(re/rs)**2)


@pytest.mark.parametrize("sigma", [0.001, 0.01, 0.05, 10])
def test_passive_light_cannot_supply_sunward_thrust(sigma):
    sun = np.array([K.AU, 0, 0]); required = np.array([0.001, 0, 0])
    optical, residual = optical_projection(required, sun, sigma)
    np.testing.assert_allclose(optical, 0, atol=1e-17)
    np.testing.assert_allclose(residual, required, atol=1e-17)


def test_random_optical_forces_obey_momentum_ball_and_outward_support():
    rng = np.random.default_rng(64)
    requests = rng.normal(size=(300, 3))*0.002
    sun = np.broadcast_to(np.array([K.AU, 0, 0]), requests.shape)
    optical, residual = optical_projection(requests, sun, 0.05, 0.138)
    b = 0.138*K.SOLAR_CONSTANT/(K.SPEED_OF_LIGHT*0.05)
    assert np.max(optical[:, 0]) <= 1e-16
    assert np.max(length(optical+np.array([b, 0, 0]))) <= b+1e-16
    assert np.all(length(residual) <= length(requests)+1e-16)


def test_forward_transmission_and_eclipse_carry_no_force():
    request = np.array([0.001, 0.0005, 0]); sun=np.array([K.AU,0,0])
    for q, light in [(0, 1), (1, 0)]:
        optical, residual = optical_projection(request, sun, 0.05, q, light)
        np.testing.assert_allclose(optical, 0)
        np.testing.assert_allclose(residual, request)


def test_photon_momentum_drops_with_additional_payload():
    request = np.array([-0.001, 0, 0]); sun=np.array([K.AU,0,0])
    f1=optical_projection(request,sun,0.05)[0]
    f2=optical_projection(request,sun,0.1)[0]
    np.testing.assert_allclose(f1,2*f2)


def test_shadow_radius_matches_extreme_straight_ray():
    d=78e6; L=K.AU; radius=4*K.MOON_RADIUS
    # Intersection at z=d of the ray from (x=radius,z=0) to (x=Rsun,z=L).
    extreme=radius+(K.SUN_RADIUS-radius)*d/L
    assert covering_radius(d,L,radius) == pytest.approx(extreme)
    assert covering_radius(d,L,radius,2000) == pytest.approx(extreme+2000)


def test_ray_flight_allowance_bounds_moving_target_displacement():
    d=78e6;v=np.array([400.,30000.,200.]);vs=np.array([0.,12.,0.])
    transverse=(length(v[1:])+length(vs[1:]))*d/K.SPEED_OF_LIGHT
    assert light_time_allowance(d,v,vs)>=transverse
    assert light_time_allowance(d,v,vs)>7000


def test_retarded_blocker_uses_barycentric_motion():
    sun=np.array([K.AU,0,0]);earth=np.array([300e6,0,0]);v=np.array([0,30000,0]);z=np.zeros(3)
    sv,ev=arriving_ray_vectors(sun,earth,z,v,z,z)
    np.testing.assert_allclose(sv,sun)
    assert ev[1]==pytest.approx(-300e6/K.SPEED_OF_LIGHT*30000)


def test_unit_vector_time_derivatives_against_independent_differences():
    p=np.array([3.,2.,1.]); v=np.array([0.3,-0.2,0.1]); a=np.array([0.03,0.01,-0.02])
    n,nd,ndd=normalized_derivatives(p,v,a)
    h=0.001
    def u(t):
        x=p+v*t+0.5*a*t*t
        return x/length(x)
    np.testing.assert_allclose(nd,(u(h)-u(-h))/(2*h),atol=1e-9)
    np.testing.assert_allclose(ndd,(u(h)-2*u(0)+u(-h))/h**2,atol=1e-9)
    assert np.dot(n,nd)==pytest.approx(0,abs=1e-16)
    assert np.dot(n,ndd)+np.dot(nd,nd)==pytest.approx(0,abs=1e-16)


@pytest.fixture(scope="module")
def ephem():
    pytest.importorskip("jplephem")
    if not DEFAULT_KERNEL.exists():
        pytest.skip("Pinned DE440s input not restored; python -m protection.dynamics.ephemeris --download")
    e=Ephemeris();yield e;e.close()


def test_bad_kernel_is_rejected(tmp_path):
    p=tmp_path/"bad.bsp";p.write_bytes(b"incomplete")
    with pytest.raises(ValueError,match="mismatch"):
        verify_kernel(p)


def test_true_ephemeris_inclination_eccentricity_and_moon_force_residual(ephem):
    t=np.linspace(0,366*K.JULIAN_DAY,600)
    s=ephem.sample(t)
    r=length(s["positions"]["earth"])
    assert r.min()<365e6 and r.max()>402e6
    eps=np.deg2rad(K.EARTH_OBLIQUITY_DEG)
    ecliptic_z=-s["positions"]["earth"][:,1]*np.sin(eps)+s["positions"]["earth"][:,2]*np.cos(eps)
    assert np.max(abs(ecliptic_z))>25e6
    assert length(moon_acceleration_residual(s)).max()<5e-9


def test_ephemeris_acceleration_difference_converges(ephem):
    t=np.linspace(0,366*K.JULIAN_DAY,100)
    a10=ephem.acceleration("moon",t,10);a30=ephem.acceleration("moon",t,30);a90=ephem.acceleration("moon",t,90)
    assert length(a10-a30).max()<2e-11
    assert length(a90-a30).max()<5e-11


def test_trajectory_derivatives_use_actual_phase_and_moving_frame(ephem):
    t=np.array([5.,15.,25.])*K.JULIAN_DAY
    c=np.zeros((3,7));c[0,:3]=[78000,9000,-6000];c[1,2]=1800;c[2,3]=1000
    r=target(geometry(ephem.sample(t)),c)
    h=20
    plus=target(geometry(ephem.sample(t+h)),c);minus=target(geometry(ephem.sample(t-h)),c)
    np.testing.assert_allclose(r["v"],(plus["q"]-minus["q"])/(2*h),atol=2e-6)
    np.testing.assert_allclose(r["a"],(plus["v"]-minus["v"])/(2*h),atol=2e-10)


def test_hermite_ephemeris_midpoints_match_direct_kernel(ephem):
    t=np.arange(0,3*K.JULIAN_DAY+1,600.)
    cache=InterpolatedEphemeris(ephem.sample(t))
    for ti in [300,4500,90000,201300]:
        s=ephem.sample([ti]);cached=cache.at(ti)
        for name in ["sun","earth","jupiter"]:
            assert length(cached["positions"][name]-s["positions"][name][0])<0.05


def test_retarded_source_matches_exact_kernel_evaluation(ephem):
    t=np.array([2.,50.,150.])*K.JULIAN_DAY;s=ephem.sample(t)
    q=target(geometry(s),np.array([[78000.,0,0,0,0,0,0],[0]*7,[0]*7]))["q"]
    sv,ev=arriving_ray_vectors(s["positions"]["sun"]-q,s["positions"]["earth"]-q,s["sun_v"],s["earth_v"],s["sun_a"],s["earth_a"])
    mp=ephem.state("moon",t)[0]
    exact=ephem.state("sun",t-length(sv)/K.SPEED_OF_LIGHT)[0]-mp-q
    assert length(sv-exact).max()<0.01
    te=np.maximum(np.sum((s["positions"]["earth"]-q)*sv/length(sv)[:,None],axis=-1),0)/K.SPEED_OF_LIGHT
    exact_e=ephem.state("earth",t-te)[0]-mp-q
    assert length(ev-exact_e).max()<0.01


def test_earth_and_sun_differential_gravity_against_circular_reference():
    # Independent closed expression used in the historical phase sweep.
    a=K.EARTH_MOON_DISTANCE;d=78e6;L=K.AU;year=K.JULIAN_DAY*K.JULIAN_YEAR_DAYS
    phases=np.linspace(0,2*np.pi,100,endpoint=False)
    moon=np.stack([a*np.cos(phases),a*np.sin(phases),np.zeros_like(phases)],axis=1)
    sun=np.broadcast_to(np.array([L,0,0]),moon.shape)
    q=np.broadcast_to(np.array([d,0,0]),moon.shape)
    ma=-K.EARTH_GM*moon/a**3+K.SUN_GM*(sun-moon)/length(sun-moon)[:,None]**3
    positions={"earth":-moon,"sun":sun-moon}
    # Other planets placed at negligible force, so their shared GM does not enter the identity.
    for k in K.SOLAR_SYSTEM_GM: positions[k]=np.broadcast_to(np.array([1e30,0,0]),moon.shape)
    required=np.broadcast_to(np.array([-d*(2*np.pi/year)**2,0,0]),moon.shape)-relative_gravity(q,positions,ma)
    screen=moon+q
    expected=-K.EARTH_GM*moon/a**3+K.EARTH_GM*screen/length(screen)[:,None]**3
    expected+=np.array([K.MOON_GM/d**2-d*(2*np.pi/year)**2,0,0])
    expected+=K.SUN_GM*((sun-moon)/length(sun-moon)[:,None]**3-(sun-screen)/length(sun-screen)[:,None]**3)
    np.testing.assert_allclose(required,expected,atol=3e-18)
