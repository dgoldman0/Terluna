"""Independent physics/conservation checks for the infrastructure screen."""
import numpy as np
import pytest
from shared import constants as K
from .electromagnetic import (square, mutual, field_and_potential, link,
                             collector, storage, plasma, rigidity, coil_hardware)


def test_square_field_and_distant_dipole_limit():
    side=10000.;s=square(side)
    b,_=field_and_potential([[0.,0.,0.]],s)
    assert b[0,2]==pytest.approx(2*np.sqrt(2)*K.VACUUM_PERMEABILITY/(np.pi*side),rel=1e-12)
    d=1e6;m=mutual(s,square(side,centre=[0,0,d]),48)
    expected=-3*K.VACUUM_PERMEABILITY*side**4/(2*np.pi*d**4)
    assert m['force_N_A2'][2]==pytest.approx(expected,rel=2e-4)


def test_extended_force_energy_gradient_and_angular_momentum():
    s=square();d=np.array([7500.,4000.,6000.]);target=square(centre=d)
    a=mutual(s,target,96);b=mutual(target,s,96)
    assert np.linalg.norm(a['force_N_A2']+b['force_N_A2'])<1e-16
    assert a['mutual_H']==pytest.approx(b['mutual_H'],rel=1e-10)
    assert np.linalg.norm(a['torque_N_m_A2']+b['torque_N_m_A2']+np.cross(d,a['force_N_A2']))<1e-12
    derivative=[]
    for axis in np.eye(3):
        derivative.append((mutual(s,square(centre=d+axis),96)['mutual_H']-
                           mutual(s,square(centre=d-axis),96)['mutual_H'])/2)
    assert np.allclose(derivative,a['force_N_A2'],rtol=1e-6,atol=1e-18)


@pytest.mark.parametrize('kind',['microwave','laser'])
def test_power_partition_and_thermal_aperture(kind):
    x=link(100e6,1e5,kind)
    assert x['source_bus_W']==pytest.approx(x['delivered_W']+x['lost_beam_W']+x['tx_heat_W']+x['rx_heat_W'])
    assert x['captured_W']==pytest.approx(x['emitted_W']*.9)
    assert x['peak_flux_W_m2']<=5000.0001
    assert x['receiver_diameter_m']>100
    c=collector(x['source_bus_W'])
    assert .9*c['incident_W']==pytest.approx(c['pv_waste_heat_W']+c['pmad_heat_W']+c['bus_W'])


def test_storage_power_rating_and_plasma_units():
    a=storage(1e6,1e9);assert a['battery_mass_kg']==a['power_limited_mass_kg']
    p=plasma();assert p['ram_pressure_Pa']==pytest.approx(1.33809754076e-9)
    assert 100e3<p['ion_inertial_m']<103e3
    assert 2300<p['electron_inertial_m']<2400
    assert rigidity(1e8)==pytest.approx(1.48297,rel=1e-4)
    a=coil_hardware(1e5);b=coil_hardware(2e5)
    assert b['conductor_mass_kg']==2*a['conductor_mass_kg']
    assert b['stored_energy_J']==4*a['stored_energy_J']
