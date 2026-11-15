"""Produce compact, hash-pinned enclosed-industry cases from an orderbook."""
from __future__ import annotations
import argparse
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import time
from . import model

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SCHEMA='terluna.provisioning.enclosed-industry/1'


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def inputs():
    identities=json.loads((HERE/'inputs/identity.json').read_text())
    for name in ('engineering','orders'):
        if digest(HERE/'inputs'/f'{name}.json')!=identities[name+'_sha256']:raise ValueError('Unreviewed '+name+' product')
    data=json.loads((HERE/'inputs/engineering.json').read_text());orders=json.loads((HERE/'inputs/orders.json').read_text())
    if data['schema']!='terluna.provisioning.enclosed-engineering-inputs/1' or orders['schema']!='terluna.humanities.enclosed-orders/1':raise ValueError('Unexpected input schema')
    return data,orders


def compact_case(case,result,scenario):
    out={k:case[k] for k in ('name','portfolio','market','market_growth','market_growth_name')}
    if scenario!='base':out['name']+=':'+scenario
    out['physical_scenario']=scenario
    out['seed']=deepcopy(result['seed']);out['seed'].pop('stock_mass_kg');out['seed'].pop('physical_annual_average')
    out['convergence']=result['convergence'];out['first_capacity_shortfall']=result['first_capacity_shortfall'];out['annual']=[]
    for row in result['annual']:
        r={k:v for k,v in row.items() if k not in ('request','required_asset_mass_kg','available_asset_mass_kg','capacity_use')}
        required=row['required_asset_mass_kg'];stock=row['available_asset_mass_kg']
        r['capacity_summary']=dict(binding=max(row['capacity_use'],key=row['capacity_use'].get),resident_workforce=row['capacity_use']['resident_workforce'])
        r['asset_stocks_kg']={k:stock.get(k,0.) for k in ('asteroid_tugs','lunar_landers','Earth_hull_tanks','source_mine_trains','source_metal_trains','bulk_trains','metal_trains','source_solar','lunar_beam_constellation','construction_and_service_robots')}
        r['idle_asset_mass_kg']=sum(max(0.,v-required.get(k,0.)) for k,v in stock.items())
        out['annual'].append(r)
    return out


def product(data,orderbook,cases,scenario='base',details=None):
    data=deepcopy(data)
    if scenario=='poor_lunar':data['industrial']['lunar'].update(bulk_yield=.8,metal_yield=.05)
    if scenario=='poor_source':data['industrial']['source'].update(water_grade=.05,water_recovery=.5,nitrogen_grade=.0012,nitrogen_recovery=.5,metal_grade=.04,metal_recovery=.5)
    if scenario=='deterministic_lifetimes':data['lifetime_shape']=None
    out=[];cache_dir=ROOT/'research/runs/enclosed_industry/cache';cache_dir.mkdir(parents=True,exist_ok=True)
    for case in cases:
        start=time.monotonic()
        identity=json.dumps(dict(case=case,scenario=scenario,model=digest(HERE/'model.py'),runner=digest(__file__),engineering=digest(HERE/'inputs/engineering.json')),sort_keys=True,separators=(',',':')).encode()
        cache=cache_dir/(hashlib.sha256(identity).hexdigest()+'.json')
        if cache.exists() and not details:
            saved=json.loads(cache.read_text());out.append(saved);print(case['name'],scenario,'verified-cache',flush=True);continue
        try:result=model.solve(data,case['annual'],max_iterations=50)
        except (ValueError,OverflowError) as error:
            out.append(dict(name=case['name']+(':'+scenario if scenario!='base' else ''),portfolio=case['portfolio'],market=case['market'],market_growth=case['market_growth'],market_growth_name=case['market_growth_name'],physical_scenario=scenario,convergence=dict(converged=False,error=str(error)),annual=[],seed=None,first_capacity_shortfall=case['annual'][0]['year']))
            print(case['name'],scenario,'UNRESOLVED',str(error),flush=True);continue
        if details:
            details.mkdir(parents=True,exist_ok=True);(details/(case['name'].replace(':','_')+'_'+scenario+'.json')).write_text(json.dumps(result,separators=(',',':'))+'\n')
        saved=compact_case(case,result,scenario);out.append(saved);cache.write_text(json.dumps(saved,separators=(',',':'))+'\n')
        print(case['name'],scenario,'iterations',result['convergence']['iterations'],'first_shortfall',result['first_capacity_shortfall'],'seconds',round(time.monotonic()-start,2),flush=True)
    return dict(schema=SCHEMA,producer=dict(path=str(Path(__file__).relative_to(ROOT)),sha256=digest(__file__),model_path=str((HERE/'model.py').relative_to(ROOT)),model_sha256=digest(HERE/'model.py')),
        dependencies={str((HERE/'inputs'/name).relative_to(ROOT)):digest(HERE/'inputs'/name) for name in ('engineering.json','orders.json','identity.json')},
        shared_constants={k:getattr(model,k) for k in ('JULIAN_DAY','JULIAN_YEAR_DAYS','STANDARD_GRAVITY','STEFAN_BOLTZMANN')},
        evidence='A demand-sized equipment/cohort and process calculation with inherited prices, chemical freight and proposed physical coefficients. Convergence and component mass balances are numerical checks, not plant, trajectory, source-reserve or market qualification.',
        units='kg, W, persons, constant USD2041. Annual physical flows are per year; source-tug counts are expected available equivalents after cohort retirement.',
        reading_rule='Only feasible annual rows admit full requested production; failed rows retain full requirement costs but no qualified sales capacity. Seed cost and feed are separate preproduction totals, included once. Cash operating_total and capital_total each include their own resource provision; do not add resource_provision again. Service and hospitality identified costs are subsets of those totals, not additions.',
        scenarios=dict(lunar=data['industrial']['lunar'],source=data['industrial']['source'],lifetime_shape=data.get('lifetime_shape',3.)),
        rules=[
            'Settlement stock follows occupied population and dated construction. Industrial assets follow actual process, freight, power and heat requirements, with25 percent planning spare capacity.',
            'Seed industrial equipment is imported once. Gross seed purchases include losses before opening; additional procurement in years6/7 is charged once within seed total. Procurement years0–7 precede operating year8.',
            'Each order is manufactured over two years before commissioning. Expected component survival follows a Weibull distribution with shape3 and the inherited mean life; this uncalibrated distribution is an explicit scheduling sensitivity.',
            'A conservative fixed point only increases each dated asset requirement within the solver; surplus surviving equipment can be used later. Unneeded expired equipment is retired with zero sale value and no automatic replacement.',
            'Local structural retirement supplies separate recycled-material inventories. Product coproduct feed is the maximum water, nitrogen and virgin-metal requirement, never their sum.',
            'Service facilities use inherited research equipment/floor intensities as a priced proxy; hotels use the inherited per-bed assets, water/air and visitor consumption. Passenger vehicles and actual access capacity remain unqualified.',
            'Technical labour is an asset-value-scaled inherited staffing proxy plus explicit resident, commercial, hotel and service staff. No productivity learning is assumed; commercial payroll remains in the consumer account.',
            'Local component transfer values are diagnostics. Only actual imported purchases and payroll enter cash totals; industrial payroll allocated to construction is removed from operating payroll.',
            'Identified service stock and mean-renewal values are noncash credits against an unitemized facility-package quote. They are not further expenditure.',
            'Future continuation for40 years beyond the reported horizon closes near-end construction feedback. Demand and population stay at the final requested level in that numerical continuation.'
        ],cases=out)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--case',action='append');ap.add_argument('--scenario',choices=('base','poor_lunar','poor_source','deterministic_lifetimes'),default='base');ap.add_argument('--output',type=Path,default=ROOT/'research/runs/enclosed_industry/industry.json');ap.add_argument('--details',type=Path);args=ap.parse_args()
    data,book=inputs();cases=[c for c in book['cases'] if not args.case or c['name'] in args.case]
    if not cases:raise ValueError('No matching cases')
    d=product(data,book,cases,args.scenario,args.details);args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(d,indent=2,sort_keys=True)+'\n');print(args.output)

if __name__=='__main__':main()
