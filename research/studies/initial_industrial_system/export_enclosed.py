"""Freeze inherited engineering inputs for the separate enclosed-industry model."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from . import evolution as old
from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS, STANDARD_GRAVITY, STEFAN_BOLTZMANN

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
COMMIT='146a558275f486c9fd9a5cf79908b7c6fdf21888'

def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def product():
    q,src,s,c=old.read_inputs()
    paths={str((HERE/n).relative_to(ROOT)) for n in ('evolution.py','economy.py','revenue.py','run.py','transport.py','economy_inputs.json')}
    paths.update(x['path'] for x in q['inputs'].values())
    identities={}
    for name in sorted(paths):
        actual=digest(ROOT/name)
        pinned=hashlib.sha256(subprocess.check_output(['git','-C',str(ROOT),'show',COMMIT+':'+name])).hexdigest()
        if actual!=pinned:raise ValueError('Changed inherited source: '+name)
        identities[name]=actual
    items=old.stocks(q,src,c)
    templates={}
    for a in items:
        kind=old.component_class(a['name']);part=q['components'][kind]
        templates[a['name']]={**a,'class':kind,
            'structural_share':max(a['inherited_local'],part['structural_mass_fraction']),
            'structural_life':max(a['life'],part['structural_life_years'])}
    city=src['settlement']['cases']['reference']
    return dict(schema='terluna.provisioning.enclosed-engineering-inputs/1',
        evidence='Inherited engineering scenario inputs, frozen separately from the former population-proportional architecture. Component partitions, source grades and trajectories are unqualified designs.',
        source_commit=COMMIT,source_files=identities,
        producer=dict(path=str(Path(__file__).relative_to(ROOT)),sha256=digest(__file__)),
        constants=dict(JULIAN_DAY=JULIAN_DAY,JULIAN_YEAR_DAYS=JULIAN_YEAR_DAYS,STANDARD_GRAVITY=STANDARD_GRAVITY,STEFAN_BOLTZMANN=STEFAN_BOLTZMANN),
        normalization=dict(population=c['demand']['population'],USD2026_to_USD2041=c['cost']['inflation_factor']),
        demand=c['demand'],templates=templates,industrial=s,economy=q,
        transport=src['transport_inputs'],revenue=src['revenue_inputs'],
        original_architecture=dict(assets=c['assets'],engineering=c['engineering'],peak_flows=c['peak_flows']),
        settlement=dict(production=src['settlement_inputs']['production'],production_demand=city['schedule']['production_demand'],
            operating_import_kg_year=city['resources']['operating_import_kg_year'],
            food_operators=city['resources']['food_operators'],industry_operators=city['resources']['industry_operators']),
        common_programme_USD2041=s['economics']['research_qualification_USD2026']*c['cost']['inflation_factor'],
        reading_rule='Prices are USD2026 unless named USD2041. Convert once using normalization. Original architecture is a coefficient/control source, never a population multiplier.')

def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--output',type=Path,default=ROOT/'provisioning/enclosed_industry/inputs/engineering.json');args=ap.parse_args()
    d=product();args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(args.output)

if __name__=='__main__':main()
