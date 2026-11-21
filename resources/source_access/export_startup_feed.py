"""Extract the original eight-year seed once, including late procurement."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def export(path):
    data = json.loads(gzip.decompress(Path(path).read_bytes()))
    if data['schema'] != 'terluna.provisioning.revised-industry/1':
        raise ValueError('Unexpected industrial product')
    out = []
    for case in data['cases']:
        if case['variant'] != 'combined':
            continue
        seed = case['physical_case']['seed']
        audit = case['seed_audit']
        opening = case['physical_case']['annual'][0]
        if seed['bootstrap_years'] != 8 or opening['year'] != 8 or opening['population'] != 12500:
            raise ValueError('Changed eight-year, 12,500-resident opening')
        row = {'name': case['name'], 'portfolio': case['portfolio'], 'opening_population': opening['population']}
        for name in ('source_feed', 'lunar_raw'):
            total = seed[name+'_kg']
            base = 8*audit['physical_annual_average'][name+'_kg_year']
            late = audit['late_preproduction_physical_total'][name+'_kg_year']
            if abs(total-base-late) > max(1e-3, total*1e-12):
                raise ValueError('Seed feed does not close')
            row[name] = {'total_kg': total, 'base_kg': base, 'late_years_6_7_kg': late}
        out.append(row)
    if len(out) != 81 or len({r['name'] for r in out}) != 81:
        raise ValueError('Expected 81 distinct combined-recovery cases')
    return {
        'schema': 'terluna.resources.startup-feed/1',
        'producer': {'path': 'resources/source_access/export_startup_feed.py', 'sha256': sha(__file__)},
        'source': {'path': 'provisioning/enclosed_industry/results/revised_industry.json.gz',
                   'sha256': sha(path), 'producer': data['producer'], 'source_commit': data['source_commit']},
        'evidence': 'Exact eight-year seed aggregates from the existing 12,500-resident opening. No dated mining or refining sequence is supplied by the source.',
        'units': 'kg cumulative over elapsed years [0,8); year-8 operating flows excluded.',
        'reading_rule': 'Base plus late procurement equals total; do not add the total again. Late procurement belongs to model years 6/7 (human-readable build years 7/8). Case alternatives are not additive.',
        'cases': out,
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--product', type=Path, required=True)
    parser.add_argument('--output', type=Path, default=HERE/'inputs/startup_feed.json')
    args = parser.parse_args()
    args.output.write_text(json.dumps(export(args.product), indent=2)+'\n')
