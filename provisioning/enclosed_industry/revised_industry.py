"""Full revised physical comparison and refurbishment-aware asset stocks.

The admitted baseline, reclamation experiment and original asset export remain
unchanged. One fresh physical solution supplies each revised cash and asset case.
"""
from __future__ import annotations
import argparse
from copy import deepcopy
import gzip
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import time
from . import asset_export, model, package, reclamation, run

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SCHEMA='terluna.provisioning.revised-industry/1'
ASSET_SCHEMA='terluna.provisioning.enclosed-asset-stock/1'
ASSET_EXTRA=('new_input_commissioned_quote_USD2041','refurbished_commissioned_value_USD2041',
             'retained_core_service_restoration_USD2041','refurbished_commissioned_mass_kg',
             'new_units_commissioned_mass_kg','repair_parts_purchase_USD2041')


def canonical(value):return (json.dumps(value,sort_keys=True,separators=(',',':'),allow_nan=False)+'\n').encode()
def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def value_digest(value):return hashlib.sha256(canonical(value)).hexdigest()


def load():
    cfg=json.loads((HERE/'revised_industry_inputs.json').read_text())
    for p,h in cfg['dependencies'].items():
        if digest(ROOT/p)!=h:raise ValueError('Changed revised-industry dependency: '+p)
    spec,data,orders,baseline=reclamation.load()
    previous=json.loads(gzip.decompress((HERE/'results/reclamation.json.gz').read_bytes()))
    rollout=json.loads((HERE/cfg['rollout_orders']).read_text())
    if rollout['schema']!='terluna.humanities.enclosed-bottleneck-orders/1':raise ValueError('Unexpected rollout schema')
    if len(orders['cases'])!=81 or len(rollout['cases'])!=6:raise ValueError('Changed admitted comparison scope')
    return cfg,spec,data,orders,baseline,previous,rollout


def captured_case(spec,data,case,variant='combined'):
    """Capture raw cohorts in an isolated copy of the unchanged producer.

    Only its Experiment constructor is specialized to retain the return value of
    the same engine.solve call. The shared imported modules are never patched.
    """
    identity=importlib.util.spec_from_file_location('provisioning.enclosed_industry._revised_reclamation',HERE/'reclamation.py')
    extension=importlib.util.module_from_spec(identity);identity.loader.exec_module(extension)
    original=extension.Experiment;captured={}
    class Capture(original):
        def __init__(self,*args,**kwargs):
            super().__init__(*args,**kwargs);captured['experiment']=self;solve=self.engine.solve
            def keep(*a,**k):
                result=solve(*a,**k);captured['raw']=result;return result
            self.engine.solve=keep
    extension.Experiment=Capture
    physical,audit=extension.solve_case(spec,data,case,variant)
    return physical,audit,captured['raw'],captured['experiment'].data


def revised_cohorts(raw,registry,shape,repair_fraction,start=8):
    cohorts=asset_export.export_cohorts(raw,registry,shape,start)
    orders={(o['asset']+':'+o['component'],o['commission_year']):o for o in raw['component_orders']}
    for c in cohorts:
        order=orders.get((c['component'],c['commission_year'])) if c['origin']=='dated_order' else None
        repaired=order.get('refurbished_mass_kg',0.) if order else 0.
        if repaired<0 or repaired>c['original_mass_kg']+1e-6:raise ValueError('Refurbished mass exceeds commissioned mass')
        c['refurbished_mass_kg']=repaired
        c['new_unit_mass_kg']=c['original_mass_kg']-repaired
        c['repair_parts_mass_kg']=repaired*repair_fraction
        c['retained_core_mass_kg']=repaired*(1-repair_fraction)
    return cohorts


def revised_annual_stocks(cohorts,registry,shape,last_year):
    """Full installed service potential, with only new input quotes in WIP."""
    annual=asset_export.annual_stocks(cohorts,registry,shape,last_year)
    asset_export.check_stock_identities(annual)
    for row in annual:
        for values in [row['totals'],*row['by_functional_sector'].values()]:
            values.update(dict.fromkeys(ASSET_EXTRA,0.))
            values['new_input_commissioned_quote_USD2041']=values['commissioned_replacement_value_USD2041']
    for cohort in cohorts:
        component=registry[cohort['component']];sector=component['functional_sector'];unit=component['replacement_price_USD2041_kg']
        restored=cohort['retained_core_mass_kg']*unit;birth=cohort['commission_year'];schedule=asset_export.fabrication_schedule(cohort)
        if 0<=birth<=last_year:
            row=annual[birth]
            for v in (row['totals'],row['by_functional_sector'][sector]):
                v['refurbished_commissioned_value_USD2041']+=cohort['refurbished_mass_kg']*unit
                v['retained_core_service_restoration_USD2041']+=restored
                v['new_input_commissioned_quote_USD2041']-=restored
                v['refurbished_commissioned_mass_kg']+=cohort['refurbished_mass_kg']
                v['new_units_commissioned_mass_kg']+=cohort['new_unit_mass_kg']
        for year,fraction in schedule:
            if 0<=year<=last_year:
                row=annual[year]
                for v in (row['totals'],row['by_functional_sector'][sector]):
                    v['manufacture_accrual_value_USD2041']-=restored*fraction
                    v['repair_parts_purchase_USD2041']+=cohort['repair_parts_mass_kg']*unit*fraction
        if restored:
            for year in range(max(0,cohort['order_year']),min(birth,last_year+1)):
                accrued=restored*sum(f for y,f in schedule if y<=year);row=annual[year]
                row['totals']['construction_WIP_value_USD2041']-=accrued
                row['by_functional_sector'][sector]['construction_WIP_value_USD2041']-=accrued
                row['by_component'][cohort['component']]['construction_WIP_value_USD2041']-=accrued
    residual=check_asset_identities(annual)
    return annual,residual


def check_asset_identities(annual):
    previous={};maximum=0.
    for row in annual:
        groups={'total':row['totals'],**row['by_functional_sector']}
        for name,v in groups.items():
            p=previous.get(name,dict.fromkeys(v,0.));scale=max(v['gross_replacement_value_USD2041'],v['construction_WIP_value_USD2041'],1.)
            residuals=[v['gross_replacement_value_USD2041']-p['gross_replacement_value_USD2041']-v['commissioned_replacement_value_USD2041']+v['retired_replacement_value_USD2041'],
                       v['net_remaining_service_value_USD2041']-p['net_remaining_service_value_USD2041']-v['commissioned_replacement_value_USD2041']+v['consumption_of_fixed_capital_proxy_USD2041'],
                       v['construction_WIP_value_USD2041']-p['construction_WIP_value_USD2041']-v['manufacture_accrual_value_USD2041']+v['new_input_commissioned_quote_USD2041'],
                       v['commissioned_replacement_value_USD2041']-v['new_input_commissioned_quote_USD2041']-v['retained_core_service_restoration_USD2041']]
            if any(abs(r)>max(.1,1e-12*scale) for r in residuals):raise ValueError('Revised asset stock identity failed')
            if v['construction_WIP_value_USD2041'] < -max(.1,1e-12*scale):raise ValueError('Negative new-input WIP')
            maximum=max(maximum,*map(abs,residuals));previous[name]=v
    return maximum


def export_assets(cfg,spec,case,physical,audit,raw,data):
    registry=asset_export.components(data);shape=data.get('lifetime_shape',3.)
    cohorts=revised_cohorts(raw,registry,shape,spec['refurbishment']['repair_part_mass_fraction'],case['annual'][0]['year'])
    annual,residual=revised_annual_stocks(cohorts,registry,shape,cfg['last_year'])
    for row in raw['annual']:
        current=annual[row['year']]['by_component']
        for name,mass in row['available_asset_mass_kg'].items():
            reproduced=min(current[k]['surviving_mass_kg']/c['mass_share'] for k,c in registry.items() if c['asset']==name)
            if not math.isclose(mass,reproduced,rel_tol=2e-11,abs_tol=.01):raise ValueError('Revised asset survival mismatch')
    result=dict(name=case['name']+':combined',baseline_case=case['name'],variant='combined',portfolio=case['portfolio'],market=case['market'],growth=case['market_growth_name'],physical_scenario='combined',
                source_case_sha256=value_digest(physical),source_loop_case_sha256=value_digest(audit),
                verification=dict(same_physical_solve=True,available_stock_reproduced=True,maximum_stock_identity_residual_USD2041=residual),
                seed_cost_registry=deepcopy(physical['seed']['cost']),seed_total_resource_provision_USD2041=physical['seed']['total_resource_provision_USD2041'],
                remaining_unpriced_requirements=physical['seed']['remaining_unpriced_requirements'],
                annual_cost_registry=[dict(year=r['year'],**r['cash']) for r in physical['annual']],annual=annual,
                cohort_account=dict(count=len(cohorts),seed_count=sum(c['origin']=='seed' for c in cohorts),sha256=value_digest(cohorts),
                                    refurbished_original_mass_kg=sum(c['refurbished_mass_kg'] for c in cohorts),reading_rule='Each complete component enters once. Repaired core and newly imported parts partition its mass; neither is an additional complete unit.'))
    for row in result['annual']:
        if row['year'] not in cfg['component_checkpoint_years']:row.pop('by_component')
    return result,registry


def baseline_seed_audit(data,case,control):
    _,seed=model.seed_stock(data,case['annual'][0]);average=seed['physical_annual_average'];source=data['industrial']['source']
    residual=control['seed']['source_feed_kg']-8*average['source_feed_kg_year']
    if residual < -max(.001,1e-10*control['seed']['source_feed_kg']):raise ValueError('Seed average exceeds frozen total feed')
    residual=max(0.,residual);products={}
    for name,field,grade,recovery in [('water','source_water_kg_year','water_grade','water_recovery'),('nitrogen','source_nitrogen_kg_year','nitrogen_grade','nitrogen_recovery'),('virgin_metal','source_virgin_metal_kg_year','metal_grade','metal_recovery')]:
        exact=8*average[field];upper=residual*source[grade]*source[recovery]
        products[name]=dict(exact_average_total_kg=exact,late_product_lower_kg=0.,late_product_upper_kg=upper,total_lower_kg=exact,total_upper_kg=exact+upper)
    return dict(name=case['name'],source_case_sha256=value_digest(control),exact_source_feed_kg=control['seed']['source_feed_kg'],exact_lunar_raw_kg=control['seed']['lunar_raw_kg'],
                average_source_feed_kg=8*average['source_feed_kg_year'],late_source_feed_kg=residual,products=products,
                evidence='Exact inherited seed average and frozen total feed. Late product amounts are physical coproduct-yield bounds, not probability intervals or an exact seed-product split.')


def plan_cases(cfg,orders,rollout):
    cases=[]
    for c in orders['cases']:cases.append(dict(case=c,baseline_case=c['name'],variant='combined',group='primary',source_override={}))
    for c in orders['cases']:
        if c['portfolio'] in cfg['sensitivity_portfolios'] and c['market']==cfg['sensitivity_market'] and c['market_growth_name'] in cfg['sensitivity_growths']:
            cases.append(dict(case=c,baseline_case=c['name'],variant='combined_poor_source',group='poor_source',source_override=cfg['poor_source']))
    for c in rollout['cases']:
        experiment=c['experiment']
        if experiment['variant']!='combined_rollout' or experiment.get('engineering_overrides',{}):raise ValueError('Unexpected rollout intervention')
        cases.append(dict(case=c,baseline_case=experiment['baseline_name'],variant='combined_rollout',group='rollout',source_override={}))
    return cases


def solve_planned(cfg,spec,data,planned,previous):
    case=planned['case'];d=deepcopy(data);d['industrial']['source'].update(planned['source_override'])
    try:
        physical,audit,raw,configured=captured_case(spec,d,case)
    except (ValueError,OverflowError) as error:
        failed=dict(name=case['name'],portfolio=case['portfolio'],market=case['market'],market_growth=case['market_growth'],market_growth_name=case['market_growth_name'],seed=None,annual=[],convergence=dict(converged=False,error=str(error)),first_capacity_shortfall=case['annual'][0]['year'])
        audit=dict(name=planned['baseline_case']+':'+planned['variant'],baseline_case=planned['baseline_case'],variant=planned['variant'],comparison_group=planned['group'],order_case_name=case['name'],source_order_case=case['name'],source_order_case_sha256=value_digest(case),source_coefficients=deepcopy(d['industrial']['source']),source_override=deepcopy(planned['source_override']),admitted_completed_path=False,physical_case=failed,annual=[],evidence='The physical solver could not complete this requested path; no completed stock or return is inferred.')
        return dict(case=audit,asset=None,registry=None)
    if planned['group']=='primary' and case['name'] in previous:
        if audit!=previous[case['name']]:raise ValueError('Previously admitted combined case changed')
    # The original nested metadata stays the unchanged orderbook identity. The
    # wrapper explicitly identifies new engineering and intervention provenance.
    audit['name']=planned['baseline_case']+':'+planned['variant'];audit['baseline_case']=planned['baseline_case'];audit['variant']=planned['variant']
    audit['comparison_group']=planned['group'];audit['order_case_name']=case['name'];audit['source_order_case']=case['name'];audit['source_order_case_sha256']=value_digest(case)
    audit['source_coefficients']=deepcopy(configured['industrial']['source']);audit['source_override']=deepcopy(planned['source_override'])
    audit['admitted_completed_path']=physical['convergence']['converged'] and physical['first_capacity_shortfall'] is None
    asset=None;registry=None
    if planned['group']=='primary' and case['market'] in cfg['asset_markets'] and case['market_growth_name'] in cfg['asset_growths']:
        asset,registry=export_assets(cfg,spec,case,physical,audit,raw,configured)
    return dict(case=audit,asset=asset,registry=registry)


def product(cfg,spec,data,orders,baseline,previous,rollout,selection=None):
    all_planned=plan_cases(cfg,orders,rollout);selected=[p for p in all_planned if selection is None or p['baseline_case']+':'+p['variant'] in selection]
    if not selected:raise ValueError('No revised cases selected')
    old_combined={c['baseline_case']:c for c in previous['cases'] if c['variant']=='combined'}
    cache=ROOT/'research/runs/enclosed_industry/revised';cache.mkdir(parents=True,exist_ok=True)
    cases=[];assets=[];registry=None
    for planned in selected:
        start=time.monotonic();identity=value_digest(dict(cfg=cfg,spec=spec,planned=planned,producer=digest(__file__)))
        path=cache/(identity+'.json.gz')
        if path.exists():
            record=json.loads(gzip.decompress(path.read_bytes()))
            if record['identity']!=identity or record['value_sha256']!=value_digest(record['value']):raise ValueError('Revised cache integrity failed')
            value=record['value'];cached=True
        else:
            value=solve_planned(cfg,spec,data,planned,old_combined);cached=False
            path.write_bytes(gzip.compress(canonical(dict(identity=identity,value_sha256=value_digest(value),value=value)),mtime=0))
        cases.append(value['case'])
        if value['asset'] is not None:
            assets.append(value['asset'])
            if registry is not None and registry!=value['registry']:raise ValueError('Revised component registry differs across primary cases')
            registry=value['registry']
        print(cases[-1]['name'],'cache' if cached else 'computed','capacity',cases[-1]['physical_case']['first_capacity_shortfall'],round(time.monotonic()-start,2),'seconds',flush=True)
    controls={c['name']:c for c in baseline['cases']}
    seed_controls=[baseline_seed_audit(data,c,controls[c['name']]) for c in orders['cases']]
    dependencies={**cfg['dependencies'],'provisioning/enclosed_industry/revised_industry_inputs.json':digest(HERE/'revised_industry_inputs.json')}
    producer=dict(path='provisioning/enclosed_industry/revised_industry.py',sha256=digest(__file__))
    physical=dict(schema=SCHEMA,producer=producer,dependencies=dependencies,source_commit=cfg['source_commit'],rules=cfg['rules'],
                  units='kg, kWh, W, persons and constant USD2041; years 0–220.',
                  evidence='Full original-market comparison under explicit recovery, refurbishment and cost hypotheses. Numerical closure and capacity screens are distinct from source access, industrial qualification and market validation.',
                  baseline=dict(path='provisioning/enclosed_industry/results/industry.json.gz',sha256=digest(HERE/'results/industry.json.gz'),case_count=len(controls),case_identities={k:value_digest(v) for k,v in controls.items()},seed_material_intervals=seed_controls),
                  scenario_parameters=dict(reclamation=spec,poor_source=cfg['poor_source']),cases=cases)
    asset_cfg=json.loads((HERE/'asset_export_inputs.json').read_text())
    stock=dict(schema=ASSET_SCHEMA,producer=producer,dependencies=dependencies,source_schema=SCHEMA,source_commit=cfg['source_commit'],
               units='kg, years, constant USD2041 replacement-quote equivalents; restored service potential is separate from purchased inputs.',
               evidence='A conditional remaining-service stock calculation from the same physical cohorts. Procurement avoided is neither a second asset nor capital income. Full-life restoration remains an engineering hypothesis.',
               rules=cfg['rules'],component_checkpoint_years=cfg['component_checkpoint_years'],lifetime_shape=data.get('lifetime_shape',3.),
               components=registry or {},cases=assets,unvalued_assets=asset_cfg['unvalued_assets']+['Retired cores awaiting inspection or repair have no separate resale or salvage valuation. Restored service value appears only at commissioning.'])
    return physical,stock


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--case',action='append');ap.add_argument('--output',type=Path,default=HERE/'results/revised_industry.json.gz');ap.add_argument('--assets',type=Path,default=HERE/'results/revised_asset_stock.json.gz');args=ap.parse_args()
    physical,assets=product(*load(),selection=args.case)
    scratch=ROOT/'research/runs/enclosed_industry/revised-product.json';scratch.parent.mkdir(parents=True,exist_ok=True);scratch.write_bytes(canonical(physical))
    manifest=package.package(scratch,args.output)
    assets['source_physical_sha256']=manifest['sha256'];assets['source_physical_path']=str(args.output.relative_to(ROOT)) if args.output.is_relative_to(ROOT) else str(args.output)
    scratch_assets=ROOT/'research/runs/enclosed_industry/revised-asset-product.json';scratch_assets.write_bytes(canonical(assets));asset_manifest=package.package(scratch_assets,args.assets)
    print(json.dumps(dict(physical=manifest,assets=asset_manifest),indent=2),flush=True)


if __name__=='__main__':main()
