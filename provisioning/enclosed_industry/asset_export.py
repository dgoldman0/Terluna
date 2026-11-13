"""Export complete physical component cohorts for a separate wealth consumer.

No factories, customer orders, cash prices or physical coefficients are changed.
Every fresh case must reproduce its admitted compact physical record before any
asset valuation is published. Quotes and remaining-service values are conditional
cost-based stock measures, not market prices, ownership titles or new cash flows.
"""
from __future__ import annotations

import argparse
from collections import defaultdict
from copy import deepcopy
import gzip
import hashlib
import json
import math
from pathlib import Path
import time

from . import bootstrap, model, run

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCHEMA = 'terluna.provisioning.enclosed-asset-stock/1'
VALUES = ('surviving_mass_kg','gross_replacement_value_USD2041',
          'net_remaining_service_value_USD2041','commissioned_replacement_value_USD2041',
          'replacement_commissioned_value_USD2041','retired_replacement_value_USD2041',
          'consumption_of_fixed_capital_proxy_USD2041','net_retirement_value_USD2041',
          'net_aging_value_USD2041','construction_WIP_value_USD2041',
          'seed_WIP_value_USD2041','manufacture_accrual_value_USD2041')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def components(data):
    """Complete functional catalogue, using the existing replacement unit quotes."""
    factor = data['normalization']['USD2026_to_USD2041']
    result = {}
    for name, a in model.catalogue(data).items():
        for part, share, life in (('structure',a['structural_share'],a['structural_life']),
                                  ('active',1-a['structural_share'],a['life'])):
            if share <= 0:
                continue
            unit = a['price'] if part == 'active' or a['civil'] else data['economy']['production']['local_fabrication_value_USD2026_kg']
            result[name+':'+part] = dict(asset=name, component=part, mass_share=share,
                functional_sector=a['owner'], destination=a['destination'], material=a['material'],
                civil=a['civil'], mean_life_years=life,
                replacement_price_USD2041_kg=unit*factor,
                price_basis='Inherited component purchase quote' if part=='active' else 'Inherited local component valuation; internal invoice is not an extra cash cost')
    return result


def export_cohorts(result, registry, shape, start=8):
    """Retain source seed ages and orders; no population replication or new renewal."""
    cohorts = []
    for key, c in registry.items():
        life=c['mean_life_years'];phases=max(1,min(8,math.ceil(life)))
        denominator=math.fsum(model.survival(j+1,life,shape) for j in range(phases))
        mass=result['seed']['stock_mass_kg'].get(c['asset'],0.)*c['mass_share']/denominator
        if mass <= 0:
            continue
        for j in range(phases):
            birth=start-j-1
            cohorts.append(dict(component=key,commission_year=birth,order_year=birth-2,
                                original_mass_kg=mass,replacement_mass_kg=0.,origin='seed'))
    for o in result['component_orders']:
        cohorts.append(dict(component=o['asset']+':'+o['component'],commission_year=o['commission_year'],
            order_year=o['order_year'],original_mass_kg=o['mass_kg'],
            replacement_mass_kg=o['replacement_kg'],origin='dated_order'))
    return sorted(cohorts,key=lambda c:(c['component'],c['commission_year'],c['origin']))


def service_curves(life, shape, maximum_age):
    """Survival and undiscounted remaining-service fractions, normalized at new."""
    end=max(maximum_age+1,math.ceil(15*life)+1)
    survival=[model.survival(age,life,shape) for age in range(end+1)]
    tails=[0.]*(end+1);running=0.
    for age in range(end,-1,-1):
        running+=survival[age];tails[age]=running
    return survival, [v/tails[0] for v in tails]


def fabrication_schedule(cohort):
    """Source two-year manufacture; seed negative dates are explicitly clipped.

    The source records seed commissioning ages but no fabrication schedule. For
    seed WIP only, its two half-orders are moved to year zero if dated earlier.
    Birth-zero equipment is commissioned immediately. This is a valuation timing
    convention; it neither revises finance's uniform seed spending nor qualifies
    the missing initial deployment sequence.
    """
    y=cohort['order_year']
    return [(max(0,y),.5),(max(0,y+1),.5)] if cohort['origin']=='seed' else [(y,.5),(y+1,.5)]


def annual_stocks(cohorts, registry, shape, last_year=220):
    grouped=defaultdict(list)
    for cohort in cohorts:
        if cohort['component'] not in registry:
            raise ValueError('Unknown component in source order')
        if cohort['original_mass_kg'] < 0 or not 0 <= cohort['replacement_mass_kg'] <= cohort['original_mass_kg']+1e-6:
            raise ValueError('Negative or inconsistent source cohort mass')
        grouped[cohort['component']].append(cohort)
    curves={key:service_curves(c['mean_life_years'],shape,last_year) for key,c in registry.items()}
    annual=[]
    for year in range(last_year+1):
        by_component={};by_sector={};totals=dict.fromkeys(VALUES,0.)
        for key,c in registry.items():
            survival,net=curves[key];unit=c['replacement_price_USD2041_kg']
            v=dict.fromkeys(VALUES,0.);age_mass=0.
            for cohort in grouped[key]:
                birth=cohort['commission_year'];mass=cohort['original_mass_kg'];cost=mass*unit;age=year-birth
                schedule=fabrication_schedule(cohort)
                v['manufacture_accrual_value_USD2041']+=cost*sum(f for y,f in schedule if y==year)
                if age<0:
                    wip=cost*sum(f for y,f in schedule if y<=year)
                    v['construction_WIP_value_USD2041']+=wip
                    if cohort['origin']=='seed':v['seed_WIP_value_USD2041']+=wip
                    continue
                alive=mass*survival[age]
                v['surviving_mass_kg']+=alive;age_mass+=age*alive
                v['gross_replacement_value_USD2041']+=cost*survival[age]
                v['net_remaining_service_value_USD2041']+=cost*net[age]
                if age==0:
                    v['commissioned_replacement_value_USD2041']+=cost
                    v['replacement_commissioned_value_USD2041']+=cohort['replacement_mass_kg']*unit
                else:
                    v['retired_replacement_value_USD2041']+=cost*(survival[age-1]-survival[age])
                    depreciation=cost*(net[age-1]-net[age])
                    retirement=cost*(survival[age-1]-survival[age])*net[age-1]/survival[age-1] if survival[age-1] else 0.
                    v['consumption_of_fixed_capital_proxy_USD2041']+=depreciation
                    v['net_retirement_value_USD2041']+=retirement
                    v['net_aging_value_USD2041']+=max(0.,depreciation-retirement)
            v['mass_weighted_age_years']=age_mass/v['surviving_mass_kg'] if v['surviving_mass_kg'] else None
            if v['net_remaining_service_value_USD2041'] > v['gross_replacement_value_USD2041']+max(.01,1e-12*v['gross_replacement_value_USD2041']):
                raise ValueError('Remaining-service value exceeds gross surviving replacement cost')
            by_component[key]={field:v[field] for field in ('surviving_mass_kg',
                'gross_replacement_value_USD2041','net_remaining_service_value_USD2041',
                'construction_WIP_value_USD2041','mass_weighted_age_years')}
            sector=by_sector.setdefault(c['functional_sector'],dict.fromkeys(VALUES,0.))
            for field in VALUES:
                totals[field]+=v[field];sector[field]+=v[field]
        annual.append(dict(year=year,totals=totals,by_functional_sector=by_sector,by_component=by_component))
    return annual


def check_stock_identities(annual, tolerance=.1):
    """Both stock definitions and construction WIP must reconcile independently."""
    previous=dict.fromkeys(VALUES,0.);largest=0.
    for row in annual:
        v=row['totals']
        residuals=[v['gross_replacement_value_USD2041']-previous['gross_replacement_value_USD2041']-v['commissioned_replacement_value_USD2041']+v['retired_replacement_value_USD2041'],
                   v['net_remaining_service_value_USD2041']-previous['net_remaining_service_value_USD2041']-v['commissioned_replacement_value_USD2041']+v['consumption_of_fixed_capital_proxy_USD2041'],
                   v['construction_WIP_value_USD2041']-previous['construction_WIP_value_USD2041']-v['manufacture_accrual_value_USD2041']+v['commissioned_replacement_value_USD2041']]
        scale=max(v['gross_replacement_value_USD2041'],1.)
        if any(abs(x)>max(tolerance,1e-12*scale) for x in residuals):
            raise ValueError('Component stock/WIP identity failed')
        largest=max(largest,*map(abs,residuals));previous=v
    return largest


def verify_physics(case,result,expected,data):
    """No approximate substitution: the full published compact case must agree."""
    compact=run.compact_case(case,result,'base')
    compact=bootstrap.apply({'cases':[compact],'rules':[]},data)['cases'][0]
    if compact != expected:
        raise ValueError('Fresh source case does not reproduce admitted compact product: '+case['name'])
    return compact


def export_case(case,result,expected,data,registry,last_year):
    verify_physics(case,result,expected,data)
    shape=data.get('lifetime_shape',3.)
    cohorts=export_cohorts(result,registry,shape,case['annual'][0]['year'])
    annual=annual_stocks(cohorts,registry,shape,last_year)
    residual=check_stock_identities(annual)
    # The physical producer defines a usable complete asset as the limiting
    # active/structure equivalent. Component masses preserve surplus separately.
    for actual in result['annual']:
        stock=annual[actual['year']]['by_component']
        for asset,mass in actual['available_asset_mass_kg'].items():
            equivalent=min(stock[k]['surviving_mass_kg']/c['mass_share'] for k,c in registry.items() if c['asset']==asset)
            if not math.isclose(equivalent,mass,rel_tol=2e-11,abs_tol=.01):
                raise ValueError('Exported component survival does not reproduce available stock')
    return dict(name=case['name'],portfolio=case['portfolio'],market=case['market'],growth=case['market_growth_name'],
                physical_scenario='base',source_case_sha256=hashlib.sha256(json.dumps(expected,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
                verification=dict(exact_compact_case_reproduced=True,available_stock_reproduced=True,maximum_stock_identity_residual_USD2041=residual),
                seed_cost_registry=deepcopy(expected['seed']['cost']),
                seed_total_resource_provision_USD2041=expected['seed']['total_resource_provision_USD2041'],
                remaining_unpriced_requirements=expected['seed']['remaining_unpriced_requirements'],
                annual_cost_registry=[{'year':r['year'],**r['cash']} for r in expected['annual']],
                cohorts=cohorts,annual=annual)


def inputs():
    cfg=json.loads((HERE/'asset_export_inputs.json').read_text())
    for name,sha in cfg['dependencies'].items():
        if digest(HERE/name)!=sha:raise ValueError('Changed admitted asset-export dependency: '+name)
    data=json.loads((HERE/'inputs/engineering.json').read_text())
    orders=json.loads((HERE/'inputs/orders.json').read_text())
    physical=json.loads(gzip.decompress((HERE/'results/industry.json.gz').read_bytes()))
    return cfg,data,orders,physical


def product(cfg,data,orders,physical,case_names=None,details=None):
    registry=components(data);expected={c['name']:c for c in physical['cases']};cases=[]
    selected=[c for c in orders['cases'] if c['market']==cfg['market'] and c['market_growth_name'] in cfg['growths'] and (case_names is None or c['name'] in case_names)]
    if not selected:raise ValueError('No selected asset cases')
    for case in selected:
        start=time.monotonic()
        result=model.solve(data,case['annual'],max_iterations=50)
        exported=export_case(case,result,expected[case['name']],data,registry,cfg['last_year'])
        if details is not None:
            directory=Path(details);directory.mkdir(parents=True,exist_ok=True)
            raw=json.dumps(exported,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
            (directory/(case['name'].replace(':','_')+'.json.gz')).write_bytes(gzip.compress(raw,mtime=0))
        compact=deepcopy(exported)
        cohorts=compact.pop('cohorts')
        compact['cohort_account']=dict(count=len(cohorts),seed_count=sum(c['origin']=='seed' for c in cohorts),
            sha256=hashlib.sha256(json.dumps(cohorts,sort_keys=True,separators=(',',':')).encode()).hexdigest(),
            detailed_output='Regenerate with --details into ignored research/runs/enclosed_industry/asset_export')
        for row in compact['annual']:
            if row['year'] not in cfg['component_checkpoint_years']:row.pop('by_component')
        cases.append(compact)
        print(case['name'],'asset export verified',round(time.monotonic()-start,2),'seconds',flush=True)
    return dict(schema=SCHEMA,producer=dict(path='provisioning/enclosed_industry/asset_export.py',sha256=digest(__file__)),
                dependencies={**{'provisioning/enclosed_industry/'+k:v for k,v in cfg['dependencies'].items()},
                              'provisioning/enclosed_industry/asset_export_inputs.json':digest(HERE/'asset_export_inputs.json')},
                source_schema=physical['schema'],source_commit=cfg['source_commit'],
                units='kg, years, constant USD2041 replacement-quote equivalents',
                evidence='Exact source-physics reproduction plus a conditional cost-based component valuation. It establishes neither market prices, legal ownership, collateral recovery nor empirical depreciation.',
                rules=cfg['rules'],component_checkpoint_years=cfg['component_checkpoint_years'],lifetime_shape=data.get('lifetime_shape',3.),components=registry,cases=cases,
                unvalued_assets=cfg['unvalued_assets'])


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--case',action='append')
    ap.add_argument('--details',type=Path,help='Ignored directory for complete per-cohort and annual component records')
    ap.add_argument('--output',type=Path,default=HERE/'results/asset_stock.json.gz');args=ap.parse_args()
    cfg,data,orders,physical=inputs();p=product(cfg,data,orders,physical,args.case,args.details)
    raw=json.dumps(p,sort_keys=True,separators=(',',':'),allow_nan=False).encode();encoded=gzip.compress(raw,mtime=0)
    args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_bytes(encoded)
    manifest=dict(schema=p['schema'],sha256=hashlib.sha256(encoded).hexdigest(),decompressed_sha256=hashlib.sha256(raw).hexdigest(),
                  bytes=len(encoded),decompressed_bytes=len(raw),case_count=len(p['cases']),dependencies=p['dependencies'],producer=p['producer'])
    args.output.with_suffix('.manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest),flush=True)


if __name__=='__main__':main()
