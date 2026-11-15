"""Admit existing resource and power records without importing their models."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def identity(root, relative, revision):
    committed = subprocess.run(['git', 'show', f'{revision}:{relative}'], cwd=root,
                               check=True, capture_output=True).stdout
    expected = hashlib.sha256(committed).hexdigest()
    if digest(root / relative) != expected:
        raise ValueError('Source differs from its recorded revision: ' + relative)
    return dict(path=relative, sha256=expected, revision=revision)


def export(resources, provisioning, humanities, output):
    resources, provisioning, humanities, output = map(Path, (resources, provisioning, humanities, output))
    output.mkdir(parents=True, exist_ok=True)
    records = {}
    for name, relative in [('assays', 'resources/source_access/coproduct_assays.json'),
                           ('flows', 'resources/source_access/inputs/coproduct_flows.json'),
                           ('recovery', 'resources/source_access/coproduct_inputs.json')]:
        path = resources / relative
        data = json.loads(path.read_text())
        records[name] = dict(**identity(resources, relative, '84886850'), schema=data['schema'])
        shutil.copyfile(path, output / (name + '.json'))
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
    records['power'] = dict(**identity(provisioning, str(path.relative_to(provisioning)), 'adeaa8ea'),
                            schema=power['schema'])
    for name, root, relative, revision in [
        ('protocols', resources, 'resources/source_access/protocols.md', '84886850'),
        ('partial_gravity', humanities, 'research/studies/world_economy/manufacturing.md', 'a63cd71b'),
        ('slow_settling', humanities, 'research/studies/world_economy/slow_settling.md', 'a63cd71b'),
    ]:
        records[name] = identity(root, relative, revision)
    data = dict(schema='terluna.research.product-location-basis/1',
                producer=dict(path='research/studies/product_allocation/export_inputs.py', sha256=digest(__file__)),
                source_identities=records, power=chosen,
                reading_rule='Exact assay, checkpoint-flow and recovery records are copied without changing their evidence state. Power entries select the reference equipment and freight cases at a 100-year costing horizon for continuous 1 GW mean service, with 100 or 10 percent imported mass. Neither these prices nor the physical flow flags qualify a new refining or manufacturing plant. Protocol records preserve separate physical and data requirements.')
    (output / 'basis.json').write_text(json.dumps(data, indent=2) + '\n')
    return data


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--resources-root', type=Path, required=True)
    parser.add_argument('--provisioning-root', type=Path, required=True)
    parser.add_argument('--humanities-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    export(args.resources_root, args.provisioning_root, args.humanities_root, args.output)


if __name__ == '__main__':
    main()
