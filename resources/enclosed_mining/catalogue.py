"""Additional geometric screen of the exact retained JPL catalogue.

Orbital elements and assumed spherical masses do not qualify a cargo route or
an economic reserve. This extracts a reproducible comparison, not mine sites.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path

HERE = Path(__file__).resolve().parent


def summarize(catalogue, density, limits):
    rows = [dict(zip(catalogue['fields'], r)) for r in catalogue['data']]
    if int(catalogue['count']) != len(rows):
        raise ValueError('Incomplete catalogue response')
    selected = []
    for row in rows:
        if any(row.get(k) is None for k in ('diameter', 'albedo', 'a', 'e', 'i')):
            continue
        diameter, albedo, axis, eccentricity, inclination = [float(row[k]) for k in ('diameter', 'albedo', 'a', 'e', 'i')]
        if not all(math.isfinite(x) for x in (diameter, albedo, axis, eccentricity, inclination)):
            raise ValueError('Nonfinite catalogue parameter')
        if not (diameter > 0 and 0 <= albedo <= limits['maximum_albedo'] and 0 < axis <= limits['maximum_semimajor_axis_AU']
                and 0 <= eccentricity <= limits['maximum_eccentricity'] and 0 <= inclination <= limits['maximum_inclination_deg']):
            continue
        if any(name in row['full_name'] for name in limits['excluded_science_names']):
            continue
        mass = math.pi/6*(diameter*1000)**3*density
        selected.append({**{k: row[k] for k in ('spkid', 'full_name', 'diameter', 'albedo', 'a', 'e', 'i')},
                         'assumed_mass_kg': mass})
    selected.sort(key=lambda r:r['assumed_mass_kg'],reverse=True)
    return {'count':len(selected), 'mass_at_assumed_density_kg':math.fsum(r['assumed_mass_kg'] for r in selected),
            'largest_objects':selected[:10], 'selected_object_ids':[r['spkid'] for r in selected]}


def product(config, catalogue_path=None):
    source = config['catalogue_source']
    path = Path(catalogue_path or source['archive'])
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != source['sha256']:
        raise ValueError('Archived catalogue hash mismatch')
    screen = config['wider_geometry_screen']
    return {'schema':'terluna.resources.enclosed-source-screen/1',
            'producer':{'path':'resources/enclosed_mining/catalogue.py','sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
            'source':source, 'evidence':'Geometry and assumed-density screening of an incomplete observed catalogue; no route, composition, permission or reserve qualification.',
            'units':{'mass':'kg','density':'kg/m3','semimajor_axis':'AU','inclination':'degrees'},
            'density_assumption_kg_m3':screen['density_assumption_kg_m3'], 'filter':screen['limits'],
            'groups':{'wider_geometry_excluding_named_science_targets':summarize(json.loads(raw),screen['density_assumption_kg_m3'],screen['limits'])}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--catalogue',type=Path)
    parser.add_argument('--output',type=Path,default=HERE/'inputs/geometry_screen.json')
    args=parser.parse_args()
    config=json.loads((HERE/'inputs.json').read_text())
    result=product(config,args.catalogue)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(args.output)

if __name__=='__main__':
    main()
