"""Catalogue volume screen, not a reserve or a source-selection certificate."""
from __future__ import annotations
import argparse
import hashlib
import json
import math
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
ARCHIVE=Path('/media/projectspace/terluna-research/initial-industrial-system/revenue')
EXPECTED='14f20ec966039c91e884e7c4ae58508a10f0ecf87acd5e971392f912b38baaa0'


def summarize(catalogue,density=1500):
    rows=[dict(zip(catalogue['fields'],row)) for row in catalogue['data']]
    if int(catalogue['count'])!=len(rows):raise ValueError('Incomplete catalogue response')
    measured=[]
    for row in rows:
        if row.get('diameter') is None:continue
        d=float(row['diameter'])
        if d<=0:continue
        r=dict(row);r['volume_m3']=math.pi/6*(d*1000)**3
        r['assumed_mass_kg']=r['volume_m3']*density;measured.append(r)
    dark=[r for r in measured if r['albedo'] is not None and float(r['albedo'])<=.1]
    geometry=[r for r in dark if .8<=float(r['a'])<=1.3 and float(r['e'])<=.3 and float(r['i'])<=15]
    named=['Bennu','Ryugu','Itokawa','Eros']
    screened=[r for r in geometry if not any(n in r['full_name'] for n in named)]
    def aggregate(group):
        return {'count':len(group),'volume_m3':sum(r['volume_m3'] for r in group),
                'mass_at_assumed_density_kg':sum(r['assumed_mass_kg'] for r in group),
                'largest_objects':[{k:r[k] for k in ['spkid','full_name','diameter','albedo','a','e','i','assumed_mass_kg']} for r in sorted(group,key=lambda r:r['assumed_mass_kg'],reverse=True)[:10]]}
    return {'catalogue_count':len(rows),'with_diameter':len(measured),'without_diameter':len(rows)-len(measured),
            'density_assumption_kg_m3':density,'groups':{n:aggregate(g) for n,g in [('all_with_diameter',measured),('dark_with_diameter',dark),('dark_geometry_filter',geometry),('dark_geometry_excluding_named_science_targets',screened)]},
            'filter':{'maximum_albedo':.1,'semimajor_axis_AU':[.8,1.3],'maximum_eccentricity':.3,'maximum_inclination_deg':15,'named_science_targets_excluded_in_last_group':named},
            'limits':['Diameters are catalogue parameters; no uniform measurement accuracy or error propagation is claimed.','Missing diameters remain unknown; no albedo-based diameter is silently substituted.','A uniform density and spherical volume give scenario masses, not measured masses.','Darkness and this orbit filter establish neither water grade nor the assumed transport delta-v.','Excluded names illustrate science/heritage screening; all other bodies remain unapproved.','Known-body stock is incomplete and is not a million-year population forecast; dynamical replenishment must not be treated as an inexhaustible local reserve.']}


def run(path=ARCHIVE/'nea_catalogue.json'):
    raw=Path(path).read_bytes();digest=hashlib.sha256(raw).hexdigest()
    if digest!=EXPECTED:raise ValueError('JPL snapshot hash mismatch; preserve the expected source')
    meta=json.loads((ARCHIVE/'nea_retrieval.json').read_text())
    return {'schema':'terluna.research.initial-resource-inventory/1','date':'2026-11-12',
            'evidence':'Incomplete catalogue and assumed-density volume screen; no permitted or economically recoverable reserve.',
            'producer':{'path':'research/studies/initial_industrial_system/resource_inventory.py','sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},
            'external_source':{'url':meta['url'],'sha256':digest,'bytes':len(raw),'archive':str(path),'signature':meta['signature']},
            'inventory':summarize(json.loads(raw))}


def main():
    p=argparse.ArgumentParser();p.add_argument('--catalogue',type=Path,default=ARCHIVE/'nea_catalogue.json');p.add_argument('--output',type=Path,default=HERE/'results/resource_inventory.json');a=p.parse_args()
    a.output.write_text(json.dumps(run(a.catalogue),indent=2)+'\n')
if __name__=='__main__':main()
