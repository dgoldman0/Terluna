"""Apply the existing finite-source assessment to bounded industry experiments."""
from __future__ import annotations

import argparse
from copy import deepcopy
import gzip
import json
from pathlib import Path

from . import run

HERE=Path(__file__).resolve().parent


def assess(physical,inventory,config):
    # An unresolved empty engineering result supplies no invented feed path.
    admitted=deepcopy(physical);admitted['cases']=[];unresolved=[]
    for case in physical['cases']:
        if not case['annual']:
            unresolved.append(dict(name=case['name'],physical_scenario=case['physical_scenario'],
                                   convergence=case['convergence'],first_capacity_shortfall=case['first_capacity_shortfall']))
        else:admitted['cases'].append(case)
    cases=run.assess_physical_cases(admitted,inventory,config)
    for case in admitted['cases']:
        key=case['physical_scenario']+':'+case['name']
        cases[key]['experiment']=deepcopy(case.get('experiment',{}))
        cases[key]['engineering_intervention']=deepcopy(case.get('engineering_intervention',{}))
    return cases,unresolved


def product():
    specification=json.loads((HERE/'bottlenecks_inputs.json').read_text())
    if specification['schema']!='terluna.resources.enclosed-bottleneck-inputs/1':
        raise ValueError('Unexpected diagnostic resource input schema')
    for path,sha in specification['reference_dependencies'].items():
        if run.digest(HERE/path)!=sha:raise ValueError('Changed reference resource method: '+path)
    config=json.loads((HERE/'inputs.json').read_text());inventory=run.inventory_inputs(config)
    source=specification['physical_product'];physical=run.load_checked(source)
    cases,unresolved=assess(physical,inventory,config)
    dependencies={**specification['reference_dependencies'],
                  'bottlenecks_inputs.json':run.digest(HERE/'bottlenecks_inputs.json'),
                  source['path']:source['sha256'],source['manifest_path']:source['manifest_sha256'],
                  config['inventory_product']['path']:config['inventory_product']['sha256'],
                  config['geometry_product']['path']:config['geometry_product']['sha256']}
    return dict(schema='terluna.resources.enclosed-mining/1',
        producer=dict(domain='resources',path='resources/enclosed_mining/bottlenecks.py',sha256=run.digest(__file__)),
        dependencies=dependencies,
        evidence='Conditional finite-source and lunar-excavation comparison for bounded commercial and component-life interventions. Existing inventory, access fractions and geometry assumptions are unchanged; production and costs come only from the pinned physical rerun.',
        units=dict(feed='kg',area='km2',time='elapsed programme year'),
        physical_products=[{**deepcopy(source), 'producer':physical['producer'],'evidence':physical['evidence'],'scenarios':physical['scenarios']}],
        seed_timing_rule='Aggregate seed feed is allocated evenly over bootstrap_years before the first operating year, once; this is the inherited average bootstrap screen, not a qualified deployment sequence.',
        reading_rule='Keys combine physical_scenario and exact diagnostic case name. Capacity-failed cases with dated requirements retain those conditional feed demands and their failure flags; empty unresolved cases receive no substitute trajectory. Hypothetical catalogue allowances establish no qualified reserve or route.',
        rules=config['rules'],access_rule=config['access_rule'],footprint_assumptions=config['lunar_footprint'],
        cases=cases,unresolved=unresolved)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=HERE/'results/bottlenecks.json.gz')
    args=parser.parse_args();result=product()
    raw=(json.dumps(result,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_bytes(gzip.compress(raw,mtime=0))
    print(args.output);print(run.digest(args.output));print('Cases:',len(result['cases']),'empty unresolved:',len(result['unresolved']))


if __name__=='__main__':main()
