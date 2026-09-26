"""Reproduce the screening calculations used in the Terluna bundle.
These are first-order engineering checks, not a climate or atmospheric escape model.
"""
import math
R_MOON=1737400.0
AREA=4*math.pi*R_MOON**2
G=1.624
R_GAS=8.314462618
MEAN_MOLAR_KG_MOL=0.02909011

def atmosphere_mass(surface_pressure_pa: float) -> float:
    return surface_pressure_pa*AREA/G

def column_mass(surface_pressure_pa: float) -> float:
    return surface_pressure_pa/G

def scale_height(temp_k: float=285.0) -> float:
    return R_GAS*temp_k/(MEAN_MOLAR_KG_MOL*G)

def pressure_at_altitude(surface_pressure_pa: float, altitude_m: float, temp_k: float=285.0) -> float:
    return surface_pressure_pa*math.exp(-altitude_m/scale_height(temp_k))

def global_equivalent_depth(water_mass_kg: float) -> float:
    return water_mass_kg/(1000.0*AREA)

def equivalent_sphere_diameter(mass_kg: float, density_kg_m3: float=1000.0) -> float:
    return 2*(3*mass_kg/(4*math.pi*density_kg_m3))**(1/3)

if __name__ == '__main__':
    print('80 kPa atmosphere mass, kg:', atmosphere_mass(80_000))
    print('80 kPa column mass, kg/m^2:', column_mass(80_000))
    print('Scale height at 285 K, km:', scale_height()/1000)
    print('3e17 kg water GED, m:', global_equivalent_depth(3e17))
