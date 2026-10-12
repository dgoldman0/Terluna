"""Independent balances and limiting cases for the aerophyte mechanics screens."""
import json
import math
from pathlib import Path
import pytest
from research.studies.aerophytes import mechanics as m

ROOT=Path(__file__).resolve().parents[3]

@pytest.fixture(scope='module')
def air():
    p=json.loads((ROOT/'research/studies/sky_ships/results/sky_ships.json').read_text())
    assert p['schema']=='terluna.research.sky-ships/1'
    return next(a for a in p['air'] if a['height_km']==10)


def test_hydrostatic_pressure_difference(air):
    rho=air['density_kg_m3'];gas=rho-air['lift_hydrogen_kg_m3'];g=m.gravity_at_height(10000)
    r=m.pressure_loads(100,rho,gas,g)
    assert math.isclose(r['hydrostatic_head_pa'],1.098*g*200)
    assert r['crown_static_pressure_pa']>0 # bottom-zero-pressure is not zero stress
    assert m.pressure_loads(200,rho,gas,g)['hydrostatic_head_pa']==2*r['hydrostatic_head_pa']
    gust=m.pressure_loads(100,rho,gas,g,gust_m_s=20)
    assert math.isclose(gust['design_differential_pressure_pa']-r['design_differential_pressure_pa'],.5*rho*400)


def test_isothermal_head_small_body_limit_and_large_body(air):
    rho=air['density_kg_m3'];gas=rho-air['lift_hydrogen_kg_m3'];p=air['pressure_atm']*101325
    for radius in (1.,100.,1000.):
        exact=m.isothermal_hydrostatic_head(radius,10000,p,rho,gas)
        local=m.pressure_loads(radius,rho,gas,m.gravity_at_height(10000))
        # Corrections are small at these scales; the exact model tends to local hydrostatics.
        assert math.isclose(exact['crown_pressure_pa'],local['hydrostatic_head_pa'],rel_tol=.006)
        assert exact['ambient_bottom_pa']>p>exact['ambient_top_pa']
    with pytest.raises(ValueError):m.isothermal_hydrostatic_head(20000,10000,p,rho,gas)


def test_sphere_hemisphere_force_balance():
    radius=30.;pressure=200.;allow=1e6
    s=m.spherical_skin(radius,pressure,allow,1500)
    # Equatorial ring tension balances the pressure force on a hemisphere.
    force_tension=2*math.pi*radius*allow*s['structural_thickness_m']
    force_pressure=pressure*math.pi*radius**2
    assert math.isclose(force_tension,force_pressure)
    assert math.isclose(s['total_structure_mass_kg'],s['area_m2']*1500*s['structural_thickness_m'])


def test_minimum_skin_and_barrier_are_not_double_counted():
    s=m.spherical_skin(1,1,1e6,1500,minimum_thickness_m=1e-4,barrier_areal_mass_kg_m2=.05)
    assert s['structural_thickness_m']==1e-4
    assert math.isclose(s['structural_areal_mass_kg_m2'],.15)
    assert math.isclose(s['total_structure_mass_kg']/s['area_m2'],.2)
    assert m.partition_mass(1,.2)==pytest.approx(math.pi*.2)


def test_hydrostatic_size_scaling(air):
    rho=air['density_kg_m3'];gas=rho-air['lift_hydrogen_kg_m3'];g=m.gravity_at_height(10000)
    values=[]
    for radius in (10,20):
        p=m.pressure_loads(radius,rho,gas,g)
        s=m.spherical_skin(radius,p['design_differential_pressure_pa'],1e6,1500)
        values.append(s)
    assert values[1]['structural_thickness_m']==pytest.approx(4*values[0]['structural_thickness_m'])
    assert values[1]['structural_mass_kg']==pytest.approx(16*values[0]['structural_mass_kg'])


def test_lobes_keep_global_tendon_cost():
    s=m.lobed_skin_tendons(100,500,1,1e6,1500,200e6,1500)
    assert s['tendon_mass_kg']>0
    assert s['tendon_count']*s['tendon_force_n']==pytest.approx(s['total_equatorial_force_n'])
    assert s['tendon_count']*s['tendon_cross_section_m2']*s['tendon_length_m']*1500==pytest.approx(s['tendon_mass_kg'])
    fine=m.lobed_skin_tendons(100,500,.5,1e6,1500,200e6,1500)
    assert fine['film_mass_kg']==pytest.approx(.5*s['film_mass_kg'])
    assert fine['tendon_mass_kg']==pytest.approx(s['tendon_mass_kg'])
    assert fine['tendon_count']>s['tendon_count']


def test_analytic_window_matches_independent_mass_budgets(air):
    lift=air['lift_hydrogen_kg_m3'];rho=air['density_kg_m3'];g=m.gravity_at_height(10000)
    win=m.spherical_payload_window(lift,g,1e6,1500,10,loaded_lift_fraction=.5,extra_pressure_pa=10,barrier_areal_mass_kg_m2=.05)
    assert win['homogeneous_sphere_screen_feasible']
    for radius in win['feasible_radius_window_m']:
        pressure=m.pressure_loads(radius,rho,rho-lift,g,bottom_overpressure_pa=10)
        skin=m.spherical_skin(radius,pressure['design_differential_pressure_pa'],1e6,1500,barrier_areal_mass_kg_m2=.05)
        budget=m.supported_payload(radius,lift,skin['total_structure_mass_kg'],loaded_lift_fraction=.5)
        assert budget['remaining_wet_payload_kg_per_projected_m2']==pytest.approx(10)
    r=win['optimum_radius_m'];a=win['polynomial_coefficients']['a'];b=win['polynomial_coefficients']['b'];c=win['polynomial_coefficients']['c']
    assert a*r-b*r*r-c> a*(1.1*r)-b*(1.1*r)**2-c
    assert a*r-b*r*r-c> a*(.9*r)-b*(.9*r)**2-c
    assert not m.spherical_payload_window(lift,g,1e6,1500,100)['homogeneous_sphere_screen_feasible']


def test_trim_conserves_displacement_and_hydrogen_component(air):
    radius=100;rho=air['density_kg_m3'];lift=air['lift_hydrogen_kg_m3'];gas=rho-lift
    volume=4*math.pi*radius**3/3
    for fraction in (0.,.25,.5,1.):
        n=m.neutral_trim(radius,rho,gas,fraction*lift*volume)
        assert n['required_lifting_mix_volume_fraction']==pytest.approx(fraction)
        assert n['required_hydrogen_volume_fraction']==pytest.approx(.98*fraction)
        assert n['actual_mixture_mass_kg']==pytest.approx(n['full_mixture_mass_kg']*fraction)
        assert n['actual_hydrogen_mass_kg']==pytest.approx(n['full_hydrogen_mass_kg']*fraction)
        assert n['total_mass_closure_kg']==pytest.approx(0,abs=1e-8)
        assert n['gas_barrier_area_spherical_minimum_m2']<=n['gas_barrier_area_full_exterior_m2']
    bad=m.neutral_trim(radius,rho,gas,1.01*lift*volume)
    assert not bad['neutral_trim_possible'] and bad['ballonet_air_mass_kg'] is None


def test_unused_lift_is_not_actual_upward_force(air):
    radius=50;lift=air['lift_hydrogen_kg_m3'];rho=air['density_kg_m3']
    support=m.supported_payload(radius,lift,1000,loaded_lift_fraction=.5)
    n=m.neutral_trim(radius,rho,rho-lift,support['remaining_wet_payload_kg']+1000)
    assert n['required_lifting_mix_volume_fraction']==pytest.approx(.5)
    assert support['unused_lift_capacity_kg']>0
    assert n['total_mass_closure_kg']==pytest.approx(0,abs=1e-8)


def test_orifice_initial_flux_limits():
    a=m.initial_hole_leak(1e-6,200,.1,1e5)
    b=m.initial_hole_leak(2e-6,800,.1,1e5)
    assert b['initial_gas_mass_kg_s']==pytest.approx(4*a['initial_gas_mass_kg_s'])
    assert m.initial_hole_leak(1e-6,0,.1,1e5)['initial_gas_mass_kg_s']==0
    with pytest.raises(ValueError):m.initial_hole_leak(1e-6,20000,.1,1e5)


@pytest.mark.parametrize('args',[(0,1e3,1e6,1500),(1,-1,1e6,1500),(1,1e3,0,1500),(1,1e3,1e6,0)])
def test_nonphysical_sphere_parameters_rejected(args):
    with pytest.raises(ValueError):m.spherical_skin(*args)


def test_residual_rain_descent_closes_force_balance(air):
    area=math.pi*100**2;g=m.gravity_at_height(10000);rho=air['density_kg_m3']
    speeds=[]
    for rain_kg_m2 in (1,5,10):
        mass=rain_kg_m2*area
        speed=m.residual_descent_speed(mass,g,rho,area)
        assert .5*rho*.47*area*speed**2==pytest.approx(mass*g)
        speeds.append(speed)
    assert speeds[1]/speeds[0]==pytest.approx(math.sqrt(5))
    assert speeds[2]/speeds[0]==pytest.approx(math.sqrt(10))
    assert m.residual_descent_speed(0,g,rho,area)==0
