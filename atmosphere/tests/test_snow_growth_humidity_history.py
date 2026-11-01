"""Distinguish unknown seed-shape sensitivities from observed initial states."""
import json
import math
import pytest
from atmosphere.snow import growth_humidity_history as history
from atmosphere.snow import growth_humidity_validation as humidity


def observation():
    return next(r for r in json.loads(humidity.DATA.read_text())['cases'] if r['id']=='G')


@pytest.mark.parametrize('phi',[.1,.3,1.])
def test_fixed_span_sensitivity_preserves_bulk_density_and_declares_mass_change(phi):
    sphere=history.configuration(observation(),'both',2,1.,919.)
    case=history.configuration(observation(),'both',2,phi,919.)
    assert case['initial_equatorial_m']==sphere['initial_equatorial_m']
    assert case['initial_polar_m']/case['initial_equatorial_m']==pytest.approx(phi)
    assert case['initial_mass_kg']/sphere['initial_mass_kg']==pytest.approx(phi)
    volume=math.pi/6*case['initial_equatorial_m']**2*case['initial_polar_m']
    assert case['initial_mass_kg']/volume==pytest.approx(919.)


def test_growth_package_and_timestep_change_preserve_experimental_forcing():
    original=history.configuration(observation(),'native',1,.3,919.)
    refined=history.configuration(observation(),'native',2,.3,919.,history.DT/4)
    assert original['growth_mode']==1 and refined['growth_mode']==2
    assert refined['dt_s']==original['dt_s']/4
    for k in ('temperature_K','vapor_pressure_Pa','pressure_Pa','initial_mass_kg','initial_equatorial_m','initial_polar_m'):
        assert original[k]==refined[k]


@pytest.mark.parametrize('growth_mode,phi',[(0,.3),(3,.3),(2,.2)])
def test_undeclared_shape_or_package_rejected(growth_mode,phi):
    with pytest.raises(ValueError):history.configuration(observation(),'native',growth_mode,phi,919.)
