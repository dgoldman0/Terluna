"""Guard the published measurement units and the controlled seed-shape change."""
import math
import pytest
from atmosphere.snow import scale_crystal_growth as growth


def test_table_full_lengths_are_converted_to_native_radii():
    # At age 3 min, the published fit is ~0.398 mm across and 0.00779 mm
    # thick, with ~0.000384 mg of ice. These bounds catch mm/m, g/kg and
    # full-axis/semiaxis mistakes without asserting the same fit formula.
    p=growth.initial('observed')
    assert .397e-3 < 2*p['equatorial_radius_m'] < .398e-3
    assert 7.78e-6 < 2*p['polar_radius_m'] < 7.80e-6
    assert 3.83e-10 < p['mass_kg'] < 3.84e-10
    assert 590 < p['density_kg_m3'] < 600
    end=growth.observed(30)
    assert 4.17e-3 < end['diameter_m'] < 4.18e-3
    assert 7.28e-8 < end['mass_kg'] < 7.30e-8


def test_spherical_control_keeps_initial_mass_and_span():
    plate=growth.initial('observed');sphere=growth.initial('sphere')
    assert plate['mass_kg']==sphere['mass_kg']
    assert plate['equatorial_radius_m']==sphere['equatorial_radius_m']
    for p in (plate,sphere):
        volume=4*math.pi/3*p['equatorial_radius_m']**2*p['polar_radius_m']
        assert volume*p['density_kg_m3']==pytest.approx(plate['mass_kg'],rel=2e-15)
    assert sphere['equatorial_radius_m']==sphere['polar_radius_m']
    assert 10 < sphere['density_kg_m3'] < 12


@pytest.mark.parametrize('age',[0,2.99,30.01,float('inf'),float('nan')])
def test_growth_fit_rejects_extrapolation(age):
    with pytest.raises(ValueError):growth.observed(age)
