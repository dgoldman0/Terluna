"""Dimensional controls, limiting cases and evidence identity for activity screens."""
import json
import math

import pytest

from shared.provenance import constants_changed
from research.studies.world_economy import distinctive_activities as m


def test_complete_product_and_bound_constants_reproduce():
    result=m.results()
    assert result == json.loads((m.HERE/'results/distinctive_activities.json').read_text())
    assert not constants_changed(result['producer']['constants'])


def test_changed_input_identity_is_rejected(tmp_path):
    spec,_=m.load_inputs()
    next(iter(spec['products'].values()))['sha256']='0'*64
    path=tmp_path/'inputs.json';path.write_text(json.dumps(spec))
    with pytest.raises(ValueError,match='Input changed'):m.load_inputs(path)


def test_flight_ratios_match_dimensional_lift_drag_and_momentum_power():
    r=m.results();env=r['environment'];f=r['flight']
    mass,area,cl,cd,disk,rho0=120.,14.,.8,.04,20.,1.16
    rows=[]
    for g,rho in [(env['earth_gravity_m_s2'],rho0),
                  (env['moon_gravity_m_s2'],rho0*env['binary_density_ratio'])]:
        weight=mass*g
        speed=math.sqrt(2*weight/(rho*area*cl))
        drag=rho*speed**2*area*cd/2
        rows.append((speed,drag*speed,drag,weight*math.sqrt(weight/(2*rho*disk))))
    for i,key in enumerate(['speed_at_fixed_lift_coefficient',
            'aerodynamic_power_at_fixed_lift_and_drag_coefficients',
            'work_per_distance_at_fixed_lift_to_drag','ideal_induced_hover_power']):
        assert rows[1][i]/rows[0][i] == pytest.approx(f[key])
    assert all(v==1 for v in m.flight_scalings(1,1).values())
    dense=m.flight_scalings(1,4)
    assert dense['speed_at_fixed_lift_coefficient']==.5
    assert dense['work_per_distance_at_fixed_lift_to_drag']==1
    assert dense['parasite_drag_at_fixed_speed_and_coefficient']==4


@pytest.mark.parametrize('gravity',[m.STANDARD_GRAVITY,m.MOON_SURFACE_GRAVITY])
def test_ballistic_return_and_energy_close(gravity):
    speed=2.5
    result=m.ballistic(gravity,speed)
    t=result['return_time_s']
    assert speed*t-gravity*t*t/2 == pytest.approx(0,abs=1e-12)
    assert gravity*result['height_m'] == pytest.approx(speed**2/2)
    assert m.ballistic(gravity,0)==dict(height_m=0,return_time_s=0)


def test_rotating_requirements_meet_constraints_and_have_active_boundary():
    r=m.results();inputs=r['assumptions']['venue'];out=r['rotating_analogue'];v=out['venue']
    assert v['angular_speed_rad_s']**2*v['radius_m']==pytest.approx(m.MOON_SURFACE_GRAVITY)
    comparisons=[(v['fractional_gravity_change'],inputs['gravity_tolerance_fraction']),
                 (v['down_direction_change_degrees'],inputs['direction_tolerance_degrees']),
                 (v['coriolis_to_rim_gravity'],inputs['coriolis_tolerance_fraction'])]
    assert all(a<=b*(1+1e-12) for a,b in comparisons)
    assert any(a==pytest.approx(b) for a,b in comparisons)
    smaller=m.rotating_venue(m.MOON_SURFACE_GRAVITY,out['combined_radius_m']*.99,
        inputs['height_m'],inputs['circumferential_arc_width_m'],inputs['transverse_speed_m_s'])
    assert smaller['down_direction_change_degrees']>inputs['direction_tolerance_degrees']


def test_rotation_respects_axial_motion_and_zero_extent_control():
    a=m.rotating_venue(m.MOON_SURFACE_GRAVITY,1000,0,0,0)
    assert a['fractional_gravity_change']==0
    assert a['down_direction_change_degrees']==0
    assert a['coriolis_to_rim_gravity']==0
    # With omega on the z axis only x/y velocity contributes to 2 omega cross v.
    velocity=(3.,4.,100.)
    b=m.rotating_venue(m.MOON_SURFACE_GRAVITY,1000,2,10,math.hypot(*velocity[:2]))
    omega=b['angular_speed_rad_s']
    vector=(-2*omega*velocity[1],2*omega*velocity[0],0.)
    assert math.sqrt(sum(x*x for x in vector))/m.MOON_SURFACE_GRAVITY==pytest.approx(b['coriolis_to_rim_gravity'])


def test_natural_rotation_uses_sidereal_period():
    r=m.results()
    value=r['natural_moon_venue']['maximum_coriolis_to_gravity']
    velocity=r['assumptions']['venue']['transverse_speed_m_s']
    inferred_period=4*math.pi*velocity/(value*m.MOON_SURFACE_GRAVITY)
    assert inferred_period/m.JULIAN_DAY==pytest.approx(m.SIDEREAL_MONTH_DAYS)


def test_hydrostatic_pressure_closure_and_extra_surface_pressure_penalty():
    r=m.results();water=r['assumptions']['water'];env=r['environment']
    for name in ('earth','moon'):
        pressure=env[f'{name}_pressure_Pa'];gravity=env[f'{name}_gravity_m_s2']
        depth=r['water']['depth_m'][name]
        assert pressure+water['density_kg_m3']*gravity*depth==pytest.approx(water['target_absolute_pressure_Pa'])
        assert m.hydrostatic_depth(pressure,pressure,1000,gravity)==0
    equal_ambient=m.hydrostatic_depth(env['earth_pressure_Pa'],water['target_absolute_pressure_Pa'],1000,m.MOON_SURFACE_GRAVITY)
    assert equal_ambient/r['water']['depth_m']['earth']==pytest.approx(1/env['gravity_ratio'])
    assert r['water']['depth_m']['moon']<equal_ambient


def test_gravity_wave_shallow_deep_and_gravity_controls():
    g=m.STANDARD_GRAVITY
    shallow=m.gravity_wave(g,1e6,1)
    assert shallow['phase_speed_m_s']==pytest.approx(math.sqrt(g),rel=1e-9)
    deep=m.gravity_wave(g,10,1000)
    assert deep['phase_speed_m_s']==pytest.approx(math.sqrt(g*10/(2*math.pi)))
    a=m.gravity_wave(g,10,20);b=m.gravity_wave(g/4,10,20)
    assert b['phase_speed_m_s']/a['phase_speed_m_s']==pytest.approx(.5)
    assert b['period_s']/a['period_s']==pytest.approx(2)
    assert a['period_s']*a['phase_speed_m_s']==pytest.approx(10)


@pytest.mark.parametrize('bad',[0,-1,True,float('nan'),float('inf')])
def test_invalid_gravity_rejected(bad):
    with pytest.raises(ValueError):m.ballistic(bad,2.5)
    with pytest.raises(ValueError):m.gravity_wave(bad,10,20)


def test_invalid_geometry_and_pressure_rejected():
    with pytest.raises(ValueError):m.rotating_venue(1,100,100,10,1)
    with pytest.raises(ValueError):m.rotating_venue(1,100,1,1000,1)
    with pytest.raises(ValueError):m.hydrostatic_depth(120000,100000,1000,1)
    with pytest.raises(ValueError):m.radius_requirements(1,1,1,1,1,1,1)
