"""Conservation references for industrial use of the locked oxygen/nitrogen air.

These calculations do not determine combustion kinetics or fire/explosion limits.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from atmosphere.radiative_convective.thermodynamics import MOLAR_CP, MOLAR_MASS
from shared.constants import GAS_CONSTANT
from shared.provenance import ROOT, constants_used

HERE = Path(__file__).resolve().parent
INPUT_PATH = HERE / 'industrial_combustion_inputs.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def positive(*values):
    if any(isinstance(v, bool) or not math.isfinite(v) or v <= 0 for v in values):
        raise ValueError('Finite positive quantities required')


def fraction(value):
    positive(value)
    if value >= 1:
        raise ValueError('Mole fraction must be below one')


def air(pressure, temperature, oxygen):
    positive(pressure, temperature)
    fraction(oxygen)
    mass = oxygen*MOLAR_MASS['O2'] + (1-oxygen)*MOLAR_MASS['N2']
    cp = oxygen*MOLAR_CP['O2'] + (1-oxygen)*MOLAR_CP['N2']
    molar_density = pressure/(GAS_CONSTANT*temperature)
    return dict(pressure_Pa=pressure, temperature_K=temperature, oxygen_mole_fraction=oxygen,
                oxygen_partial_pressure_Pa=pressure*oxygen,
                nitrogen_partial_pressure_Pa=pressure*(1-oxygen),
                oxygen_mol_m3=molar_density*oxygen,
                density_kg_m3=molar_density*mass,
                oxygen_mass_fraction=oxygen*MOLAR_MASS['O2']/mass,
                heat_capacity_J_m3_K=molar_density*cp)


def methane(pressure, temperature, oxygen, excess_ratio=1.):
    """Complete CH4 + 2 O2 -> CO2 + 2 H2O; wet products stay gaseous.

    Volumes are evaluated at the specified temperature and pressure for comparison,
    not asserted to be actual furnace temperatures or stable ambient wet mixtures.
    """
    positive(pressure, temperature, excess_ratio)
    fraction(oxygen)
    if excess_ratio < 1:
        raise ValueError('Complete conversion requires excess_ratio >= 1')
    oxidizer = 2*excess_ratio/oxygen
    nitrogen = oxidizer*(1-oxygen)
    products = dict(CO2=1., H2O=2., O2=2*(excess_ratio-1), N2=nitrogen)
    wet = sum(products.values())
    dry = wet-products['H2O']
    molar_volume = GAS_CONSTANT*temperature/pressure
    oxidizer_mass = 2*excess_ratio*MOLAR_MASS['O2']+nitrogen*MOLAR_MASS['N2']
    return dict(excess_oxygen_ratio=excess_ratio, oxidizer_mol_per_mol_CH4=oxidizer,
                oxidizer_kg_per_mol_CH4=oxidizer_mass,
                oxidizer_m3_per_mol_CH4=molar_volume*oxidizer,
                products_mol_per_mol_CH4=products, wet_product_mol_per_mol_CH4=wet,
                wet_product_volume_reference_m3_per_mol_CH4=molar_volume*wet,
                dry_CO2_mole_fraction=1/dry,
                dry_oxygen_mole_fraction=products['O2']/dry)


def mixed_purge(pressure, temperature, volume, initial_oxygen, final_oxygen):
    """N2 needed for ideal continuous perfectly mixed purging at constant P,V,T."""
    positive(pressure, temperature, volume)
    fraction(initial_oxygen)
    fraction(final_oxygen)
    if final_oxygen > initial_oxygen:
        raise ValueError('Nitrogen purge cannot increase oxygen fraction')
    exchanges = math.log(initial_oxygen/final_oxygen)
    moles = pressure*volume/(GAS_CONSTANT*temperature)*exchanges
    return dict(volume_exchanges=exchanges, nitrogen_mol=moles,
                nitrogen_kg=moles*MOLAR_MASS['N2'],
                nitrogen_m3_at_vessel_conditions=volume*exchanges,
                final_oxygen_partial_pressure_Pa=pressure*final_oxygen)


def load_inputs(path=INPUT_PATH, root=ROOT):
    spec = json.loads(Path(path).read_text())
    if spec['schema'] != 'terluna.research.industrial-combustion-inputs/1':
        raise ValueError('Unexpected industrial combustion input schema')
    products = {}
    for name, identity in spec['products'].items():
        file = Path(root)/name
        if digest(file) != identity['sha256']:
            raise ValueError(f'Input changed; review and rebind: {name}')
        data = json.loads(file.read_text())
        if data['schema'] != identity['schema']:
            raise ValueError(f'Unexpected input schema: {name}')
        products[name] = data
    return spec, products


def results():
    spec, products = load_inputs()
    source = next(iter(products.values()))
    env = source['environment']; pp = source['gas_partial_pressures']
    t = env['temperature_K']; p0 = env['reference_Pa']; p1 = env['locked_Pa']
    y0 = pp['reference_oxygen_Pa']/p0; y1 = pp['locked_oxygen_Pa']/p1
    a0, a1 = air(p0,t,y0), air(p1,t,y1)
    burning = []
    for ratio in spec['methane']['excess_oxygen_ratios']:
        b0, b1 = methane(p0,t,y0,ratio), methane(p1,t,y1,ratio)
        burning.append(dict(reference=b0, locked=b1,
            oxidizer_volume_ratio=b1['oxidizer_m3_per_mol_CH4']/b0['oxidizer_m3_per_mol_CH4'],
            oxidizer_mass_ratio=b1['oxidizer_kg_per_mol_CH4']/b0['oxidizer_kg_per_mol_CH4'],
            nitrogen_per_fuel_ratio=b1['products_mol_per_mol_CH4']['N2']/b0['products_mol_per_mol_CH4']['N2'],
            wet_product_molar_ratio=b1['wet_product_mol_per_mol_CH4']/b0['wet_product_mol_per_mol_CH4'],
            wet_product_volume_ratio_at_equal_temperature=b1['wet_product_volume_reference_m3_per_mol_CH4']/b0['wet_product_volume_reference_m3_per_mol_CH4']))
    purge = spec['purge']; v = purge['vessel_m3']
    endpoints = [('equal_mole_fraction', target, target) for target in purge['target_oxygen_mole_fractions']]
    endpoints.append(('equal_partial_pressure',purge['target_oxygen_partial_pressure_Pa']/p0,
                      purge['target_oxygen_partial_pressure_Pa']/p1))
    purges = []
    for basis, target0, target1 in endpoints:
        q0,q1 = mixed_purge(p0,t,v,y0,target0), mixed_purge(p1,t,v,y1,target1)
        purges.append(dict(basis=basis, reference_target_fraction=target0, locked_target_fraction=target1,
                           reference=q0, locked=q1, nitrogen_moles_ratio=q1['nitrogen_mol']/q0['nitrogen_mol']))
    files = ['research/studies/world_economy/industrial_combustion.py',
             'research/studies/world_economy/industrial_combustion_inputs.json',
             'research/studies/world_economy/industrial_combustion_sources.json',
             'atmosphere/radiative_convective/thermodynamics.py']
    return dict(schema='terluna.research.industrial-combustion/1', evidence=spec['evidence'],
                reading_rule=spec['reading_rule'], assumptions=spec,
                producer=dict(files={f:digest(ROOT/f) for f in files}, constants=constants_used(files)),
                air=dict(reference=a0, locked=a1,
                         density_ratio=a1['density_kg_m3']/a0['density_kg_m3'],
                         heat_capacity_per_volume_ratio=a1['heat_capacity_J_m3_K']/a0['heat_capacity_J_m3_K'],
                         oxygen_per_volume_ratio=a1['oxygen_mol_m3']/a0['oxygen_mol_m3'],
                         oxygen_per_mass_ratio=a1['oxygen_mass_fraction']/a0['oxygen_mass_fraction']),
                methane=burning, purges=purges)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=ROOT/'research/runs/world_economy/industrial_combustion.json')
    args=parser.parse_args()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(results(),indent=2,allow_nan=False)+'\n')
    print(args.output)


if __name__ == '__main__':
    main()
