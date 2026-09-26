"""Supply ledger for the protection hardware: what holding and replacing it asks of industry.

    python -m engineering.network.ledger      # writes engineering/network/results/protection_supply_ledger.json

Reads the protection module catalogue (protection/modules/catalogue.json) as
demand. For each lifetime scenario and material-recovery fraction it gives the
gross replacement flow, the fresh material after recovery, the irreversible
propellant and the power, by commodity, with totals over 10^9 years. It sets the
irreversible flows beside the loss of air that the optical shield prevents
(atmosphere/loss_response), the test of requirement S6. Rebuild power uses the
September report's illustrative 100 MJ/kg. The magnet installations are carried
both at their ideal conductor-and-support floor and at the report's planning
allowance of 3e14-3e15 kg.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys

from engineering.network.model import COMMODITIES

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
CATALOGUE = ROOT / 'protection' / 'modules' / 'catalogue.json'
LOSS_RESPONSE = ROOT / 'atmosphere' / 'loss_response' / 'results' / 'loss_response.json'
SCHEMA = 'terluna.engineering.protection-supply-ledger/1'

YEAR_S = 365.25 * 86400.0
SERVICE_YEARS = 1.0e9
RECOVERY = (0.99, 0.999)
REBUILD_J_KG = 1.0e8          # the report's illustrative rebuild energy (protection/report.md, section 15)

MATERIAL_COMMODITY = {
    'titania': 'structural_materials',
    'silica': 'structural_materials',
    'framing_coatings_metrology': 'structural_materials',
    'power_and_propulsion_hardware': 'precision_components',
    'superconductor_composite': 'conductors',
    'tensile_support': 'structural_materials',
    'propellant': 'propellant',
}


def replacement(module, life_years):
    """Gross replacement flow (kg/s) by commodity for one module family."""
    out = {}
    for material, kg in module['materials_kg'].items():
        c = MATERIAL_COMMODITY[material]
        out[c] = out.get(c, 0.0) + module['count'] * kg / (life_years * YEAR_S)
    return out


def ledger(catalogue, loss_response=None):
    modules = {m['id']: m for m in catalogue['modules']}
    replaceable = [m for m in catalogue['modules'] if m['lifetime_years_scenarios']]
    propellant = sum(m['count'] * m.get('consumables_kg_s', {}).get('propellant', 0.0) for m in catalogue['modules'])
    holding_power = sum(m['count'] * m.get('power_mean_w', 0.0) for m in catalogue['modules'])
    magnet = modules['regional_magnet_installation']
    scenarios = []
    for life in replaceable[0]['lifetime_years_scenarios']:
        gross = {}
        by_module = {}
        for m in replaceable:
            flows = replacement(m, life)
            by_module[m['id']] = sum(flows.values())
            for c, v in flows.items():
                gross[c] = gross.get(c, 0.0) + v
        total_gross = sum(gross.values())
        allowance = [a / (life * YEAR_S) for a in magnet['planning_allowance_kg_total']]
        for recovery in RECOVERY:
            fresh = total_gross * (1 - recovery)
            scenarios.append(dict(
                life_years=life, recovery=recovery, gross_replacement_kg_s=total_gross,
                gross_by_commodity_kg_s=gross, gross_by_module_kg_s=by_module,
                magnet_planning_allowance_gross_kg_s=allowance,
                fresh_material_kg_s=fresh, propellant_kg_s=propellant,
                irreversible_kg_s=fresh + propellant,
                rebuild_power_w=total_gross * REBUILD_J_KG, holding_power_mean_w=holding_power,
                over_service_kg=dict(fresh_material=fresh * SERVICE_YEARS * YEAR_S,
                                     propellant=propellant * SERVICE_YEARS * YEAR_S)))
    out = dict(service_years=SERVICE_YEARS, rebuild_j_kg=REBUILD_J_KG, scenarios=scenarios,
               propellant_kg_s=propellant, holding_power_mean_w=holding_power)
    if loss_response is not None:
        unshielded = loss_response['unshielded_energy_limited_kg_s']['quiet']
        low, high = min(unshielded.values()), max(unshielded.values())
        out['requirement_S6'] = dict(
            unshielded_loss_kg_s=[low, high],
            propellant_over_highest_unshielded_loss=propellant / high,
            verdict=('the holding scheme spends more material than the air it saves'
                     if propellant > high else
                     'propellant is below the highest unshielded loss; compare with the chosen budget'))
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out', type=Path, default=HERE / 'results')
    args = parser.parse_args(argv)
    catalogue = json.loads(CATALOGUE.read_text())
    loss = json.loads(LOSS_RESPONSE.read_text()) if LOSS_RESPONSE.is_file() else None
    result = ledger(catalogue, loss)
    digest = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()[:16]
    product = dict(
        schema=SCHEMA,
        producer=dict(domain='engineering', files={str(p.relative_to(ROOT)): digest(p) for p in (
            HERE / 'ledger.py', HERE / 'model.py', CATALOGUE) + ((LOSS_RESPONSE,) if loss else ())}),
        evidence=('Arithmetic on the protection module catalogue: replacement flows for the lifetime scenarios, fresh '
                  'material after recovery, irreversible propellant and power. Inherits every assumption of the '
                  'September reference design the catalogue records.'),
        reading_rule=('Flows are kg/s averaged over the service interval. irreversible_kg_s is fresh material plus '
                      'propellant, the material the protection consumes for good; compare it with the loss of air it '
                      'prevents (requirement_S6).'),
        commodities=COMMODITIES, **result)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / 'protection_supply_ledger.json').write_text(json.dumps(product, indent=2) + '\n')
    for s in result['scenarios']:
        print(f"life {s['life_years']:6g} yr, recovery {s['recovery']:.3f}: gross {s['gross_replacement_kg_s']:.3e} kg/s, "
              f"fresh {s['fresh_material_kg_s']:.3e} kg/s, propellant {s['propellant_kg_s']:.3e} kg/s")
    if 'requirement_S6' in result:
        print('S6:', result['requirement_S6']['verdict'])
    return 0


if __name__ == '__main__':
    sys.exit(main())
