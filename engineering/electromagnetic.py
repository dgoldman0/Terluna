"""Bounded vacuum circuit, power and plasma-scale estimates, SI throughout.

No trajectory or plasma response is solved. Thin square source segments use
analytic Biot--Savart fields; target line integrals use Gauss quadrature.
The equivalent circular self-inductance/support model is an engineering
estimate, never used to replace the actual square mutual geometry.
"""
from functools import lru_cache
import numpy as np
from numpy.polynomial.legendre import leggauss
from shared import constants as K


def square(side=10000., centre=(0., 0., 0.), frame=None):
    vertices = side/2*np.array([[-1., -1., 0.], [1., -1., 0.],
                                [1., 1., 0.], [-1., 1., 0.]])
    return vertices @ (np.eye(3) if frame is None else frame).T+centre


@lru_cache(maxsize=12)
def gauss(n):
    return leggauss(n)


def quadrature(vertices, n=32):
    x, w = gauss(n)
    start = np.asarray(vertices); delta = np.roll(start, -1, axis=0)-start
    return ((start[:, None]+(x[None, :, None]+1)*delta[:, None]/2).reshape(-1, 3),
            (w[None, :, None]*delta[:, None]/2).reshape(-1, 3))


def field_and_potential(points, vertices):
    """B and A per ampere-turn of a closed polygon; excludes the conductor."""
    p = np.atleast_2d(points); b = np.zeros_like(p); a = np.zeros_like(p)
    for start, end in zip(vertices, np.roll(vertices, -1, axis=0)):
        delta = end-start; length = np.linalg.norm(delta); u = delta/length
        r = p-start; s1 = r@u; s2 = s1-length
        transverse = r-s1[:, None]*u; rho2 = np.sum(transverse**2, axis=1)
        # On a segment's extension B=0; A remains finite. No evaluation on wire.
        near = rho2 < 1e-16
        if np.any(near & (s1 >= 0) & (s2 <= 0)):
            raise ValueError('Field point lies on filament')
        safe = np.maximum(rho2, 1e-16)
        b += np.cross(u, transverse)/safe[:, None]*(s1/np.sqrt(safe+s1*s1)-s2/np.sqrt(safe+s2*s2))[:, None]
        a += (np.arcsinh(s1/np.sqrt(safe))-np.arcsinh(s2/np.sqrt(safe)))[:, None]*u
    return b*K.VACUUM_PERMEABILITY/(4*np.pi), a*K.VACUUM_PERMEABILITY/(4*np.pi)


def mutual(source, target, n=32):
    """Force/torque on target per I_source*I_target and mutual inductance."""
    points, dl = quadrature(target, n); b, a = field_and_potential(points, source)
    df = np.cross(dl, b)
    return dict(force_N_A2=df.sum(axis=0),
                torque_N_m_A2=np.cross(points-target.mean(axis=0), df).sum(axis=0),
                mutual_H=float(np.sum(a*dl)))


def radiator(heat_W, temperature_K=350., emissivity=.9, kg_m2=5.):
    """One-sided emitting area, zero view/background load; conservative extra panel."""
    area = heat_W/(emissivity*K.STEFAN_BOLTZMANN*temperature_K**4)
    return dict(heat_W=float(heat_W), area_m2=float(area), mass_kg=float(area*kg_m2))


def coil_hardware(ampere_turns, perimeter=40000., current_density=1e8,
                  density=6000., bundle_radius=.1, cryostat_radius=.25,
                  cryostat_kg_m=5., heat_leak=.05, cop=.005,
                  specific_strength=1e6):
    radius = perimeter/(2*np.pi); logarithm = np.log(8*radius/bundle_radius)
    inductance = K.VACUUM_PERMEABILITY*radius*(logarithm-2)
    energy = .5*inductance*ampere_turns**2
    tension = K.VACUUM_PERMEABILITY*ampere_turns**2*(logarithm-1)/(4*np.pi)
    conductor = perimeter*abs(ampere_turns)*density/current_density
    support = tension*perimeter/specific_strength
    cold = perimeter*2*np.pi*cryostat_radius*heat_leak
    electrical = cold/cop; rad = radiator(electrical+cold)
    return dict(ampere_turns=float(ampere_turns), conductor_mass_kg=float(conductor),
                ideal_self_support_mass_kg=float(support), cryostat_mass_kg=perimeter*cryostat_kg_m,
                inductance_H=float(inductance), stored_energy_J=float(energy),
                hoop_tension_N=float(tension), bundle_field_T=float(K.VACUUM_PERMEABILITY*abs(ampere_turns)/(2*np.pi*bundle_radius)),
                cold_load_W=float(cold), refrigerator_input_W=float(electrical),
                refrigerator_mass_kg=float(electrical/300), radiator=rad,
                installed_mass_kg=float(conductor+support+perimeter*cryostat_kg_m+rad['mass_kg']+electrical/300),
                scope='One winding; equivalent-perimeter circular self/support estimate. No corner bending, joints, launch, gimbals or fault containment mass.')


def link(delivered_W, range_m, kind='microwave'):
    if kind == 'microwave':
        wavelength=K.SPEED_OF_LIGHT/5.8e9; eta_t=.70; eta_r=.85; diameter=100.; flux=1000.
    elif kind == 'laser':
        wavelength=1.064e-6; eta_t=.50; eta_r=.50; diameter=1.; flux=5000.
    else:
        raise ValueError(kind)
    capture=.90; bus=.95
    emitted=delivered_W/(capture*eta_r*bus); source=emitted/eta_t
    # Gaussian paraxial focused spot; transmit Gaussian radius D/3.
    # A Fresnel focus, not far-field Friis, is needed when R < 2 D^2/lambda.
    diffraction=wavelength*range_m/(np.pi*diameter/3)
    thermal=np.sqrt(2*emitted/(np.pi*flux))
    w=max(diffraction,thermal); receiver_radius=w*np.sqrt(-np.log(1-capture)/2)
    tx_heat=source-emitted; rx_heat=emitted*capture-delivered_W
    txrad=radiator(tx_heat); rxrad=radiator(rx_heat)
    aperture_mass=np.pi*(diameter/2)**2*5
    receiver_mass=np.pi*receiver_radius**2*2
    # Separate converters, generation and distribution; no existing mass credit.
    converter_mass=source/1000+delivered_W/1000
    return dict(kind=kind,delivered_W=float(delivered_W),range_m=float(range_m),
        source_bus_W=float(source),emitted_W=float(emitted),captured_W=float(emitted*capture),
        lost_beam_W=float(emitted*(1-capture)),dc_efficiency=float(delivered_W/source),
        tx_efficiency=eta_t,rx_efficiency=eta_r,capture=capture,bus_efficiency=bus,
        transmitter_diameter_m=diameter,receiver_diameter_m=float(2*receiver_radius),
        receiver_area_m2=float(np.pi*receiver_radius**2),gaussian_radius_m=float(w),
        diffraction_radius_m=float(diffraction),peak_flux_W_m2=float(2*emitted/(np.pi*w*w)),
        pointing_allowance_rad=float(.1*w/range_m),fraunhofer_distance_m=float(2*diameter**2/wavelength),
        tx_heat_W=float(tx_heat),rx_heat_W=float(rx_heat),tx_radiator=txrad,rx_radiator=rxrad,
        tx_aperture_mass_kg=float(aperture_mass),receiver_mass_kg=float(receiver_mass),
        converter_mass_kg=float(converter_mass),
        link_installed_mass_kg=float(aperture_mass+receiver_mass+converter_mass+txrad['mass_kg']+rxrad['mass_kg']),
        emitter_recoil_N=float(emitted/K.SPEED_OF_LIGHT),receiver_absorption_force_N=float(emitted*capture/K.SPEED_OF_LIGHT),
        scope='Extrapolated modules, focused line of sight, intentional flux-limited spot; no blockage/availability guarantee. Aperture taper loss folded into 0.90 capture.')


def collector(bus_W, efficiency=.30, pmad=.95, specific_power=300.):
    gross=bus_W/pmad; area=gross/(K.SOLAR_CONSTANT*efficiency)
    heat=area*K.SOLAR_CONSTANT*.9-gross
    # Both faces of a separate solar array reject PV heat; do not double count a radiator.
    temperature=(heat/(2*.85*K.STEFAN_BOLTZMANN*area))**.25
    pmad_heat=gross-bus_W; rad=radiator(pmad_heat)
    return dict(bus_W=float(bus_W),gross_electric_W=float(gross),incident_W=float(area*K.SOLAR_CONSTANT),
        area_m2=float(area),array_mass_kg=float(gross/specific_power),pv_waste_heat_W=float(heat),
        equilibrium_K=float(temperature),pmad_heat_W=float(pmad_heat),pmad_radiator=rad,
        installed_mass_kg=float(gross/specific_power+bus_W/1000+rad['mass_kg']),
        absorbed_photon_force_N=float(.9*area*K.SOLAR_CONSTANT/K.SPEED_OF_LIGHT),
        scope='Unshadowed independently sun-tracking array; 300 W/kg includes array support; PMAD and its radiator added. 0.9 absorptivity, no reflected-force credit.')


def storage(energy_J, peak_W, wh_kg=200., depth=.8, discharge=.95, charge=.95, w_kg=1000.):
    energy_mass=energy_J/(wh_kg*3600*depth*discharge); power_mass=peak_W/w_kg
    return dict(delivered_energy_J=float(energy_J),peak_discharge_W=float(peak_W),
        energy_limited_mass_kg=float(energy_mass),power_limited_mass_kg=float(power_mass),
        battery_mass_kg=float(max(energy_mass,power_mass)),
        recharge_bus_energy_J=float(energy_J/(charge*discharge)),
        discharge_heat_peak_W=float(peak_W*(1/discharge-1)),
        heat_rejection=radiator(peak_W*(1/discharge-1)),
        scope='Pack-level exploratory specific energy; no cycle-life or shared radiator credit.')


def cable(power_W, voltage=100000., length=10000., loss_fraction=.01):
    resistivity=2.82e-8; density=2700.; current=power_W/voltage
    resistance=loss_fraction*power_W/current**2
    area=resistivity*2*length/resistance
    return dict(voltage_V=voltage,current_A=float(current),conductor_area_m2=float(area),
                conductor_mass_kg=float(area*2*length*density),loss_W=float(power_W*loss_fraction),
                scope='Two aluminium conductors, 10 km route each; no insulation, shielding or HV equipment mass.')


def plasma(density_cm3=5., speed_km_s=400., temperature_eV=10., imf_nT=5.):
    eps=1/(K.VACUUM_PERMEABILITY*K.SPEED_OF_LIGHT**2)
    me=K.VACUUM_PERMEABILITY*K.ELEMENTARY_CHARGE**2/(4*np.pi*K.CLASSICAL_ELECTRON_RADIUS)
    n=density_cm3*1e6; v=speed_km_s*1000; b=imf_nT*1e-9
    pressure=n*K.PROTON_MASS*v*v
    return dict(density_cm3=density_cm3,speed_km_s=speed_km_s,temperature_eV=temperature_eV,imf_nT=imf_nT,
        ram_pressure_Pa=float(pressure),balance_B_T=float(np.sqrt(2*K.VACUUM_PERMEABILITY*pressure)),
        ion_inertial_m=float(np.sqrt(K.PROTON_MASS/(K.VACUUM_PERMEABILITY*n*K.ELEMENTARY_CHARGE**2))),
        electron_inertial_m=float(np.sqrt(me/(K.VACUUM_PERMEABILITY*n*K.ELEMENTARY_CHARGE**2))),
        debye_m=float(np.sqrt(eps*temperature_eV/(n*K.ELEMENTARY_CHARGE))),
        proton_bulk_gyroradius_m=float(K.PROTON_MASS*v/(K.ELEMENTARY_CHARGE*b)),
        electron_thermal_gyroradius_m=float(np.sqrt(2*me*temperature_eV*K.ELEMENTARY_CHARGE)/(K.ELEMENTARY_CHARGE*b)),
        motional_electric_V_m=float(v*b))


def rigidity(energy_eV):
    energy=energy_eV*K.ELEMENTARY_CHARGE; rest=K.PROTON_MASS*K.SPEED_OF_LIGHT**2
    return float(np.sqrt(energy*(energy+2*rest))/(K.SPEED_OF_LIGHT*K.ELEMENTARY_CHARGE))
