"""Module catalogue of the protection hardware, written as a data product.

    python -m protection.modules.catalogue        # writes protection/modules/catalogue.json

Each replaceable unit of the September 9 reference design, with its count, mass,
materials, power, consumables and lifetime scenarios. The engineering network
model reads it as demand. Values come only from the September model's stored
outputs (protection/results/), so they inherit its assumptions: a 50 g/m^2
allocation for film, coatings, framing and metrology; a 78,000 km thrust-held
position; 30 km/s exhaust at 70% efficiency and a 45-degree cant; 300 W/kg of
installed peak power; a seven-day propellant buffer; and four 100 km regional
magnet installations. The holding scheme is the reference, not a chosen design:
it fails requirement S6 of research/studies/protection_architecture/requirements.md.

The optical cell follows the report's 10 x 10 km division (protection/report.md,
section 8). Film materials use the stored layer thicknesses and the densities the
model uses for titania (3,900 kg/m^3) and silica (2,200 kg/m^3); the rest of the
allocation is framing, coatings and metrology, whose materials are unspecified.
"""
from __future__ import annotations
import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PROTECTION = HERE.parent
RESULTS = PROTECTION / 'results'
SCHEMA = 'terluna.protection.module-catalogue/1'

CELL_AREA_M2 = 1.0e8                 # 10 x 10 km (report, section 8)
AREAL_ALLOCATION_KG_M2 = 0.05        # film, coatings, framing and metrology (report, section 4)
TITANIA_KG_M2 = 1.0e-6 * 3900.0      # 1 um of titania
SILICA_KG_M2 = 10.0e-6 * 2200.0      # 10 um of silica
LIFETIME_SCENARIOS_YEARS = (20.0, 100.0, 1000.0)   # the report's maintenance cases (results/maintenance.csv)
MAGNET_PLANNING_ALLOWANCE_KG = (3.0e14, 3.0e15)    # four installations, report section 11


def _rows(name):
    with open(RESULTS / name, newline='') as handle:
        return list(csv.DictReader(handle))


def build():
    reference = json.loads((RESULTS / 'results.json').read_text())['orbits']['reference']
    area = reference['area_m2']
    cells = area / CELL_AREA_M2
    film_total = next(float(r['mass_kg']) for r in _rows('maintenance.csv') if r['component'] == 'near_filter')
    cell_mass = film_total / cells
    magnet = next(r for r in _rows('anchor_magnets.csv') if float(r['radius_km']) == 100.0)
    cryo = [float(r['electrical_W']) for r in _rows('magnet_cryogenics.csv')]
    installations = int(float(magnet['station_count']))
    modules = [
        dict(id='optical_cell', function='optical', count=cells,
             description='A 10 x 10 km cell of the titania-silica film with its anti-reflection layers, framing and metrology',
             unit_mass_kg=cell_mass,
             materials_kg=dict(titania=TITANIA_KG_M2 * CELL_AREA_M2, silica=SILICA_KG_M2 * CELL_AREA_M2,
                               framing_coatings_metrology=cell_mass - (TITANIA_KG_M2 + SILICA_KG_M2) * CELL_AREA_M2),
             lifetime_years_scenarios=list(LIFETIME_SCENARIOS_YEARS),
             failure='Local: a missing cell opens about 1e-6 of the aperture; overlaps and replacement overlays cover it'),
        dict(id='holding_pack', function='holding', count=cells,
             description="One cell's share of the power, processing, thruster and thermal hardware that holds the screen",
             unit_mass_kg=reference['power_system_mass_kg'] / cells,
             materials_kg=dict(power_and_propulsion_hardware=reference['power_system_mass_kg'] / cells),
             power_mean_w=reference['mean_power_W'] / cells, power_peak_w=reference['peak_power_W'] / cells,
             consumables_kg_s=dict(propellant=reference['propellant_kg_s'] / cells),
             exhaust_speed_m_s=30000.0, jet_efficiency=0.7, specific_power_w_kg=300.0,
             lifetime_years_scenarios=list(LIFETIME_SCENARIOS_YEARS),
             failure='Local loss of thrust; an unpowered screen falls toward the Moon within days, so packs need redundancy'),
        dict(id='propellant_buffer', function='holding', count=cells,
             description="One cell's share of the seven-day propellant buffer",
             unit_mass_kg=reference['fuel_buffer_mass_kg'] / cells,
             materials_kg=dict(propellant=reference['fuel_buffer_mass_kg'] / cells),
             lifetime_years_scenarios=[],
             failure='Consumed and refilled continuously'),
        dict(id='regional_magnet_installation', function='charged_particle', count=installations,
             description='A 100 km radius superconducting loop near a lunar pole, one of four with a combined moment of 1.5e21 A m^2',
             unit_mass_kg=float(magnet['conductor_mass_kg']) + float(magnet['ideal_support_mass_kg']),
             materials_kg=dict(superconductor_composite=float(magnet['conductor_mass_kg']),
                               tensile_support=float(magnet['ideal_support_mass_kg'])),
             stored_energy_j=float(magnet['stored_energy_J']),
             power_mean_w_range=[min(cryo) / installations, max(cryo) / installations],
             planning_allowance_kg_total=list(MAGNET_PLANNING_ALLOWANCE_KG),
             lifetime_years_scenarios=list(LIFETIME_SCENARIOS_YEARS),
             failure='Quench: the coupled stored energy must be contained or recovered'),
    ]
    return dict(
        schema=SCHEMA,
        producer=dict(domain='protection', files={p: hashlib.sha256((PROTECTION / p).read_bytes()).hexdigest()[:16]
                                                   for p in ('modules/catalogue.py', 'results/results.json',
                                                             'results/maintenance.csv', 'results/anchor_magnets.csv',
                                                             'results/magnet_cryogenics.csv')}),
        evidence=('Reference values of the September 9 protection design, taken from its stored model outputs. They are '
                  'component budgets, not a chosen or optimized design; the thrust-held screen fails requirement S6.'),
        reading_rule=('count x unit_mass_kg is a module family\'s held mass. Divide by a lifetime scenario for gross '
                      'replacement; consumables are continuous flows. The magnet mass is the ideal conductor and support '
                      'floor; the report plans for 3e14-3e15 kg in all.'),
        reference_design=dict(distance_km=reference['distance_km'], aperture_area_m2=area,
                              held_mass_kg=reference['total_mass_kg'], mean_power_w=reference['mean_power_W'],
                              peak_power_w=reference['peak_power_W'], propellant_kg_s=reference['propellant_kg_s']),
        modules=modules)


def main() -> int:
    product = build()
    (HERE / 'catalogue.json').write_text(json.dumps(product, indent=2) + '\n')
    for m in product['modules']:
        print(f"{m['id']:30s} count {m['count']:>12,.0f}  unit {m['unit_mass_kg']:.3e} kg")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
