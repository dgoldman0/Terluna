"""Biological precedents for holding gas, carried into aerophyte numbers.

Pure functions plus evaluate(). The swim bladder's guanine-lined wall and
countercurrent rete, the Portuguese man-of-war's gas gland and goldbeater's
skin are the precedents (gas_biology.md, gas_biology_sources.json).
"""
from __future__ import annotations

import math

from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS

BARRER_SI = 3.35e-16
# Henry's-law solubility of H2 in water at 25 C, 7.7-7.8e-6 mol m^-3 Pa^-1 (Sander's compilation).
HYDROGEN_SOLUBILITY_MOL_M3_PA = 7.8e-6
H2_MOLAR_KG = 0.002016


def _positive(**values):
    if any(not math.isfinite(v) or v <= 0 for v in values.values()):
        raise ValueError(f'finite positive inputs required: {tuple(values)}')


def _unit(name, value):
    if not math.isfinite(value) or not 0 <= value <= 1:
        raise ValueError(f'{name} must lie in [0, 1]')


def platelet_barrier(volume_fraction, aspect_ratio):
    """Permeability of a film filled with aligned impermeable plates, over the plain film.

    Nielsen's tortuosity model: P/P0 = (1 - phi)/(1 + (alpha/2) phi), with alpha
    the plates' length over thickness. Guanine plates in a swim-bladder wall are
    the biological case; the model sets how far such plates could cut hydrogen loss.
    """
    _unit('volume fraction', volume_fraction)
    _positive(aspect=aspect_ratio)
    return (1-volume_fraction)/(1+aspect_ratio*volume_fraction/2)


def dissolved_gas_loss(fluid_flow_kg_m2_day, gas_partial_pressure_pa,
                       solubility_mol_m3_pa=HYDROGEN_SOLUBILITY_MOL_M3_PA, exchanger_efficiency=0.):
    """Hydrogen carried off per projected m² by fluid that equilibrates with the gas.

    Sap or water that touches the gas leaves saturated at the gas's partial
    pressure; a countercurrent exchanger returns a share of it, as a fish's rete
    mirabile returns gas to the bladder.
    """
    _positive(pressure=gas_partial_pressure_pa, solubility=solubility_mol_m3_pa)
    if not math.isfinite(fluid_flow_kg_m2_day) or fluid_flow_kg_m2_day < 0:
        raise ValueError('nonnegative fluid flow required')
    _unit('exchanger efficiency', exchanger_efficiency)
    cubic_metres = fluid_flow_kg_m2_day/1000  # m3 of water per m2 per day
    mol_day = cubic_metres*solubility_mol_m3_pa*gas_partial_pressure_pa*(1-exchanger_efficiency)
    return dict(mol_m2_s=mol_day/JULIAN_DAY, kg_m2_year=mol_day*H2_MOLAR_KG*JULIAN_YEAR_DAYS,
                exchanger_efficiency=exchanger_efficiency)


def permeation_flux(permeability_barrer, thickness_m, partial_pressure_pa):
    _positive(permeability=permeability_barrer, thickness=thickness_m, pressure=partial_pressure_pa)
    return permeability_barrer*BARRER_SI*partial_pressure_pa/thickness_m


def evaluate(hydrogen_partial_pa=.98*100108., reference_barrer=.173, reference_thickness_m=100e-6,
             gas_contact_ratio=5., transpiration_kg_m2_day=1.254):
    """Platelet cuts, and dissolved loss against permeation for the reference barrier."""
    # Swim-bladder guanine plates: about 0.02 um thick and up to 100 um broad (Lapennas and
    # Schmidt-Nielsen 1977), 2-50 um wide and under 20 nm thick in sardine (Pinsk et al. 2022):
    # aspect ratios of about 100 to 5,000 at a small volume share of the wall.
    plates = [dict(volume_fraction=f, aspect_ratio=a, permeability_ratio=platelet_barrier(f, a),
                   thickness_for_same_loss_ratio=platelet_barrier(f, a))
              for f in (.01, .02, .05, .1) for a in (100., 500., 1000., 2500., 5000.)]
    permeation = permeation_flux(reference_barrer, reference_thickness_m, hydrogen_partial_pa)*gas_contact_ratio
    dissolved = []
    for flow in (transpiration_kg_m2_day, .25, .07):
        for eff in (0., .9, .99):
            d = dissolved_gas_loss(flow, hydrogen_partial_pa, exchanger_efficiency=eff)
            dissolved.append(dict(fluid_flow_kg_m2_day=flow, **d,
                                  share_of_permeation=d['mol_m2_s']/permeation))
    return dict(
        evidence='Nielsen tortuosity arithmetic on aligned plates and Henry\'s-law solubility of hydrogen; the precedents are measured fish, cnidarian and processed-membrane systems read in the source register.',
        reading_rule='Permeability ratios are film with plates over plain film. Dissolved loss assumes the fluid that wets the gas side leaves saturated with hydrogen at the gas partial pressure, minus what a countercurrent exchanger returns; per projected m², compared with the reference barrier\'s permeation over the stated gas-contact area.',
        reference_permeation_mol_m2_s=permeation, platelets=plates, dissolved=dissolved,
        unresolved=['Whether guanine-like plates can be laid in a hydrogen barrier and kept aligned as it grows.',
                    'How much of an aerophyte\'s sap actually wets its gas cells.'])
