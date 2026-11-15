"""Export dated virgin feed from a committed provisioning product, without rerunning it."""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent
FIELDS = (
    'source_feed_kg_year', 'source_retained_tails_kg_year',
    'source_metal_kg_year', 'source_virgin_metal_kg_year',
    'lunar_raw_kg_year', 'lunar_virgin_metal_kg_year',
)
CHECKPOINTS = (100, 200, 220)


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def export(product_path, revision):
    product_path = Path(product_path)
    data = json.loads(gzip.decompress(product_path.read_bytes()))
    if data['schema'] != 'terluna.provisioning.revised-industry/1':
        raise ValueError('Unexpected provisioning product schema')
    cases = []
    for case in data['cases']:
        if case['variant'] != 'combined':
            continue
        rows = []
        for year in CHECKPOINTS:
            matches = [r for r in case['physical_case']['annual'] if r['year'] == year]
            if len(matches) != 1:
                raise ValueError('Missing or duplicate physical checkpoint')
            row = matches[0]
            values = {k: row['physical'][k] for k in FIELDS}
            if any(not math.isfinite(v) or v < 0 for v in values.values()):
                raise ValueError('Invalid physical feed')
            rows.append({'year': year, **values, 'physical_feasible': row['feasible']})
        cases.append({'name': case['name'], 'portfolio': case['portfolio'],
                      'admitted_completed_path': case['admitted_completed_path'],
                      'checkpoints': rows})
    if len(cases) != 81 or len({c['name'] for c in cases}) != 81:
        raise ValueError('Expected 81 distinct combined-recovery cases')
    return {
        'schema': 'terluna.resources.coproduct-flows/1',
        'producer': {'path': 'resources/source_access/export_coproduct_flows.py', 'sha256': sha(__file__)},
        'source_product': {'path': 'provisioning/enclosed_industry/results/revised_industry.json.gz',
                           'sha256': sha(product_path), 'checkout_revision': revision,
                           'schema': data['schema'], 'producer': data['producer'],
                           'producer_source_commit': data['source_commit']},
        'evidence': 'Exact checkpoint extraction from the combined-recovery physical cases. Industrial scenarios, not observed demand or qualified coproduct capacity.',
        'units': 'Annual kg of fresh raw feed, retained tails, and generic metal, with explicit field suffixes.',
        'reading_rule': 'Case alternatives are not additive. Use raw source feed and lunar raw excavation once. Recycled generic metal is not additional fresh ore; tails are a subset of already processed feed. No bootstrap repetition, cumulative stock, linear interpolation or trace-element assay is supplied.',
        'cases': cases,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--product', type=Path, required=True)
    parser.add_argument('--revision', required=True)
    parser.add_argument('--output', type=Path, default=HERE/'inputs/coproduct_flows.json')
    args = parser.parse_args()
    args.output.write_text(json.dumps(export(args.product, args.revision), indent=2)+'\n')


if __name__ == '__main__':
    main()
