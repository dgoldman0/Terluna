"""Record the resource and power records the product-route screens read, with the power cases they select.

The assays, the checkpoint flows and the recovery settings are read where the resources study keeps them; this
script writes inputs/basis.json, which holds each record's path and digest and the twelve power alternatives.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


ROOT = Path(__file__).resolve().parents[3]
SOURCES = dict(assays='resources/source_access/coproduct_assays.json',
               flows='resources/source_access/inputs/coproduct_flows.json',
               recovery='resources/source_access/coproduct_inputs.json')


def identity(root, relative):
    return dict(path=relative, sha256=digest(root / relative))


def export(root, output):
    root, output = Path(root), Path(output)
    resources = provisioning = humanities = root
    output.mkdir(parents=True, exist_ok=True)
    records = {}
    for name, relative in SOURCES.items():
        data = json.loads((resources / relative).read_text())
        records[name] = dict(**identity(resources, relative), schema=data['schema'])
    path = provisioning / 'provisioning/enclosed_power/results/power_comparison.json'
    power = json.loads(path.read_text())
    chosen = []
    for case in power['comparisons']:
        if case['technology_case'] != 'reference_assumptions':
            continue
        for row in case['costs']:
            if row['freight_case'] == 'reference' and row['imported_mass_fraction'] in (1, .1):
                chosen.append(dict(architecture=case['architecture'],
                                   imported_mass_fraction=row['imported_mass_fraction'],
                                   levelized_USD2026_kWh=row['lifecycle']['100']['levelized_USD2026_MWh'] / 1000,
                                   source_horizon=100, site=case['site'], mean_MW=case['mean_MW'],
                                   technology_case=case['technology_case'], boundary=case['boundary']))
    if len(chosen) != 12:
        raise ValueError('Expected twelve admitted power alternatives')
    records['power'] = dict(**identity(provisioning, str(path.relative_to(provisioning))),
                            schema=power['schema'])
    for name, base, relative in [
        ('protocols', resources, 'resources/source_access/protocols.md'),
        ('partial_gravity', humanities, 'research/studies/world_economy/manufacturing.md'),
        ('slow_settling', humanities, 'research/studies/world_economy/slow_settling.md'),
    ]:
        records[name] = identity(base, relative)
    data = dict(schema='terluna.research.product-location-basis/1',
                producer=dict(path='research/studies/product_allocation/export_inputs.py', sha256=digest(__file__)),
                source_identities=records, power=chosen,
                reading_rule='Exact assay, checkpoint-flow and recovery records are copied without changing their evidence state. Power entries select the reference equipment and freight cases at a 100-year costing horizon for continuous 1 GW mean service, with 100 or 10 percent imported mass. Neither these prices nor the physical flow flags qualify a new refining or manufacturing plant. Protocol records preserve separate physical and data requirements.')
    (output / 'basis.json').write_text(json.dumps(data, indent=2) + '\n')
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=ROOT)
    parser.add_argument('--output', type=Path, default=Path(__file__).resolve().parent / 'inputs')
    args = parser.parse_args()
    export(args.root, args.output)


if __name__ == '__main__':
    main()
