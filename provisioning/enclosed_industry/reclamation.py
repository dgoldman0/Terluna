"""Costed reclamation experiments around the immutable enclosed-industry engine.

The isolated engine shares source bytes, not mutable module globals, with the
published control. Its fixed point still sizes all energy, freight and equipment.
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
from . import bootstrap, package, run

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SCHEMA='terluna.provisioning.reclamation/1'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load():
    spec=json.loads((HERE/'reclamation_inputs.json').read_text())
    if spec['refurbishment']['inspection_delay_years']<1:raise ValueError('Refurbishment requires a positive retirement-to-inspection lag')
    for group in ('goods','structural','refurbishment'):
        for key in ('collection','quality','repair_part_mass_fraction'):
            if key in spec[group] and not 0<=spec[group][key]<=1:raise ValueError('Recovery fractions must lie in [0,1]')
    for p,h in spec['dependencies'].items():
        if sha(ROOT/p)!=h:raise ValueError('Changed reclamation dependency: '+p)
    data,book=run.inputs()
    baseline=json.loads(gzip.decompress((HERE/'results/industry.json.gz').read_bytes()))
    return spec,data,book,baseline


def fresh_engine():
    spec=importlib.util.spec_from_file_location('_terluna_reclamation_engine',HERE/'model.py')
    engine=importlib.util.module_from_spec(spec);spec.loader.exec_module(engine)
    return engine


def water_account(spec,population,guest_years=0.,improved=False):
    scale=(population+guest_years)/spec['normalization_population'];streams={}
    for name,w in spec['water'].items():
        gross=w['gross_kg_year']*scale;r=w['improved_recovery' if improved else 'recovery']
        streams[name]=dict(gross_kg=gross,recovered_kg=gross*r,makeup_kg=gross*(1-r),loss_kg=gross*(1-r),recovery=r)
    extra=spec['landscape_agriculture_allowance_kg_year']*scale
    return dict(streams=streams,landscape_agriculture_makeup_and_loss_kg=extra,
                gross_processed_kg=sum(v['gross_kg'] for v in streams.values()),
                recovered_kg=sum(v['recovered_kg'] for v in streams.values()),
                operational_makeup_kg=extra+sum(v['makeup_kg'] for v in streams.values()))


def configure(data,spec,variant):
    d=deepcopy(data);improved=variant in ('closure','combined','combined_without_residual_allowance');repair=variant in ('refurbishment','combined','combined_without_residual_allowance')
    if improved:
        d['demand']['yearly_water']=water_account(spec,spec['normalization_population'],improved=True)['operational_makeup_kg']
        g=spec['goods'];factor=(1-g['collection']*g['quality'])/(1-g['baseline_recovery'])
        d['demand']['goods_virgin']*=factor
        d['settlement']['operating_import_kg_year']['everyday_carbon_feed']*=factor
        d['industrial']['lunar']['recovery']=d['economy']['production']['source_recovery']=spec['structural']['collection']*spec['structural']['quality']
    if variant=='tug_life_30':d['templates']['asteroid_tugs']['life']=30
    def template(name,dest,life,structure_life):
        a=deepcopy(d['templates']['fabrication_cells']);a.update(name=name,mass=0.,destination=dest,life=life,structural_life=structure_life)
        d['templates'][name]=a
    if improved:
        w=spec['water_equipment'];template('loop_water_equipment','surface',w['active_life'],w['structure_life'])
        for dest in ('surface','source'):template('loop_sort_'+dest,dest,15,30)
    if repair:
        r=spec['refurbishment']
        for dest in ('surface','source','depot'):template('loop_repair_'+dest,dest,r['workshop_active_life'],r['workshop_structural_life'])
    return d


def retirements(engine,data,rows,seed,events):
    """Reconstruct gross expired components without turning them into feed."""
    cats=engine.catalogue(data);start=rows[0]['year'];end=rows[-1]['year'];shape=data.get('lifetime_shape',3.)
    out={y:{} for y in range(start,end+1)};births={}
    for n,a in cats.items():
        for part,share,life in [('active',1-a['structural_share'],a['life']),('structure',a['structural_share'],a['structural_life'])]:
            if not share:continue
            phases=max(1,min(8,math.ceil(life)));den=sum(engine.survival(j+1,life,shape) for j in range(phases))
            births[n,part]=[(start-j-1,seed.get(n,0.)*share/den) for j in range(phases)]
    for ev in events:births[ev['asset'],ev['component']].append((ev['commission_year'],ev['mass_kg']))
    for (n,part),cohorts in births.items():
        life=cats[n]['life' if part=='active' else 'structural_life'];curve=[engine.survival(a,life,shape) for a in range(end+2)]
        values=[0.]*(end+1)
        for birth,mass in cohorts:
            for y in range(max(start+1,birth+1),end+1):
                age=y-birth;values[y]+=mass*(curve[age-1]-curve[age])
        for y in range(start,end+1):out[y][n,part]=values[y]
    return out


def empty_repair():
    return dict(retired_kg=0.,collected_kg=0.,uncollected_kg=0.,inspection_input_kg=0.,quality_reject_kg=0.,accepted_kg=0.,
                pending_inspection_kg=0.,qualified_inventory_kg=0.,withdrawn_for_repair_kg=0.,retained_old_parts_kg=0.,
                discarded_old_parts_kg=0.,repair_parts_imported_kg=0.,refurbished_commissioned_kg=0.,repair_WIP_kg=0.,
                imported_mass_avoided_kg=0.,import_purchase_avoided_USD2026=0.,mass_balance_residual_kg=0.)


def refurbishment(engine,data,spec,rows,flows,events,retired):
    """Use only previously retired, accepted components at the same destination."""
    cfg=spec['refurbishment'];cats=engine.catalogue(data);start=rows[0]['year'];end=rows[-1]['year'];lag=cfg['inspection_delay_years']
    audits={y:{k:empty_repair() for k in engine.DESTS} for y in range(start-2,end+1)}
    for y,f in flows.items():
        f['loop_inspected']={k:0. for k in engine.DESTS};f['loop_repaired']={k:0. for k in engine.DESTS}
        f['loop_consume_kg']={k:0. for k in engine.DESTS};f['loop_import_avoided']={}
    selected={n:[] for n in cfg['eligible_assets']}
    for ev in events:
        if ev['component']=='active' and ev['asset'] in selected:selected[ev['asset']].append(ev)
    for n,orders in selected.items():
        a=cats[n];dest=a['destination'];byyear={}
        for ev in orders:byyear.setdefault(ev['order_year'],[]).append(ev)
        bank=0.;pending=0.;wip=0.
        for y in range(start-2,end+1):
            v=audits[y][dest];opening_bank=bank;opening_pending=pending;opening_wip=wip
            retired_now=retired.get(y,{}).get((n,'active'),0.)
            inspected=retired.get(y-lag,{}).get((n,'active'),0.)*cfg['collection']
            accepted=inspected*cfg['quality'];rejected=inspected-accepted;collected=retired_now*cfg['collection']
            bank+=accepted;pending+=collected-inspected
            used=0.
            for ev in byyear.get(y,[]):
                use=min(bank,ev['mass_kg']);bank-=use;used+=use
                k=cfg['repair_part_mass_fraction'];avoided=use*(1-k);growth=(ev['mass_kg']-ev['replacement_kg'])/ev['mass_kg']
                for wy in (y,y+1):
                    if wy not in flows:continue
                    f=flows[wy];m=avoided/2;money=m*a['price']
                    f['imports'][dest]-=m;f['capital_imports'][dest]-=m*growth
                    f['external_capital']-=money*growth;f['external_replacement']-=money*(1-growth)
                    if a['owner']=='industry':f['industry_replacement_import_kg']-=m*(1-growth)
                    for owner in ('service','hospitality'):
                        if a['owner']==owner:
                            f[owner+'_external_capital']-=money*growth;f[owner+'_external_replacement']-=money*(1-growth)
                    f['loop_repaired'][dest]+=use/2
                    f['loop_import_avoided'][n]=f['loop_import_avoided'].get(n,0.)+m
                    av=audits[wy][dest];av['repair_parts_imported_kg']+=use*k/2
                    av['imported_mass_avoided_kg']+=m;av['import_purchase_avoided_USD2026']+=money
                cy=y+2
                if cy in audits:audits[cy][dest]['refurbished_commissioned_kg']+=use
                ev['refurbished_mass_kg']=use;ev['new_imported_mass_kg']=ev['mass_kg']-avoided
            # Per-asset scheduled arrivals are required here; destination totals
            # include other assets and must not be subtracted from every asset.
            kit=sum(ev.get('refurbished_mass_kg',0.)*cfg['repair_part_mass_fraction']/2 for ev in orders if ev['order_year'] in (y,y-1))
            commissioned=sum(ev.get('refurbished_mass_kg',0.) for ev in orders if ev['commission_year']==y)
            retained=used*(1-cfg['repair_part_mass_fraction']);discarded=used-retained
            wip+=retained+kit-commissioned
            values=dict(retired_kg=retired_now,collected_kg=collected,uncollected_kg=retired_now-collected,inspection_input_kg=inspected,
                        quality_reject_kg=rejected,accepted_kg=accepted,pending_inspection_kg=pending,qualified_inventory_kg=bank,
                        withdrawn_for_repair_kg=used,retained_old_parts_kg=retained,discarded_old_parts_kg=discarded,repair_WIP_kg=wip)
            for key,value in values.items():v[key]+=value
            v['mass_balance_residual_kg']+=opening_bank+opening_pending+opening_wip+retired_now+kit-(bank+pending+wip+retired_now-collected+rejected+discarded+commissioned)
            flows[y]['loop_inspected'][dest]+=inspected
    for y,f in flows.items():
        for dest in engine.DESTS:
            mass=(f['loop_inspected'][dest]+f['loop_repaired'][dest])*cfg['consumables_fraction_of_processed_mass']
            f['loop_consume_kg'][dest]=mass;f['imports'][dest]+=mass
    return audits


class Experiment:
    def __init__(self,data,spec,variant):
        self.spec=deepcopy(spec)
        if variant=='combined_without_residual_allowance':self.spec['landscape_agriculture_allowance_kg_year']=0.
        spec=self.spec;self.variant=variant;self.improved=variant in ('closure','combined','combined_without_residual_allowance');self.repair=variant in ('refurbishment','combined','combined_without_residual_allowance')
        self.data=configure(data,spec,variant);self.engine=fresh_engine();self.last_flows=None;self.last_retired=None;self.last_audit=None
        self.original_ledger=self.engine.cohort_ledger;self.original_physical=self.engine.physical;self.original_requirements=self.engine.requirements
        self.engine.cohort_ledger=self.ledger
        self.engine.physical=self.physical
        if variant!='baseline':self.engine.requirements=self.requirements

    def ledger(self,data,rows,targets,seed):
        flows,stocks,events=self.original_ledger(data,rows,targets,seed)
        ret=retirements(self.engine,data,rows,seed,events)
        audit=refurbishment(self.engine,data,self.spec,rows,flows,events,ret) if self.repair else None
        cats=self.engine.catalogue(data)
        for y,f in flows.items():
            f['loop_structural_retired']={k:0. for k in (*self.engine.DESTS,'bulk')}
            for (n,part),mass in ret.get(y,{}).items():
                if part=='structure':
                    a=cats[n];key='bulk' if a['material']=='bulk' else a['destination'];f['loop_structural_retired'][key]+=mass
        self.last_flows=flows;self.last_retired=ret;self.last_audit=audit
        return flows,stocks,events

    def physical(self,data,row,flow,bank=None):
        out=self.original_physical(data,row,flow,bank)
        if self.variant=='baseline':self.last_physical=deepcopy(out);return out
        extra={k:0. for k in self.engine.DESTS};spec=self.spec
        if self.improved:
            water=water_account(spec,row['population'],row.get('visitor_nights_year',0.)/data['constants']['JULIAN_YEAR_DAYS'],True)
            extra['surface']+=water['gross_processed_kg']*spec['water_equipment']['extra_kWh_per_kg_gross']
            retired=flow.get('loop_structural_retired',{});s=spec['structural']
            extra['surface']+=retired.get('bulk',0.)*s['bulk_extra_sort_kWh_kg']
            for dest in self.engine.DESTS:extra[dest]+=retired.get(dest,0.)*s['metal_extra_sort_kWh_kg']
            # Orbital structural metal is already shipped to the source in the
            # baseline; retain that routing for its extra sorting work.
            extra['source']+=extra['depot']+extra['leo'];extra['depot']=extra['leo']=0.
            g=spec['goods'];extra['surface']+=row['population']*g['gross_kg_person_year']*(g['collection']*g['quality']-g['baseline_recovery'])
        if self.repair:
            r=spec['refurbishment']
            for dest in self.engine.DESTS:
                extra[dest]+=flow.get('loop_inspected',{}).get(dest,0.)*r['inspection_kWh_kg']+flow.get('loop_repaired',{}).get(dest,0.)*r['remanufacture_kWh_kg']
            out['operating_material_USD2026']+=sum(flow.get('loop_consume_kg',{}).values())*data['industrial']['prices']['other_special_USD2026_kg']
        watts={k:v*1000/self.engine.HOURS for k,v in extra.items()}
        out['lunar_process_W']+=watts['surface'];out['ground_firm_W']+=watts['surface']*data['industrial']['economics']['capacity_margin']/data['industrial']['lunar']['process_availability']
        out['source_process_W']+=watts['source'];out['source_heat_W']+=watts['source'];out['depot_W']+=watts['depot']
        out['loop_extra_process_kWh_year']=extra
        out['loop_repair_input_kg_year']={k:flow.get('loop_inspected',{}).get(k,0.)+flow.get('loop_repaired',{}).get(k,0.) for k in self.engine.DESTS}
        out['loop_structural_retired_kg_year']=flow.get('loop_structural_retired',{})
        out['loop_repair_consumables_kg_year']=dict(flow.get('loop_consume_kg',{}))
        self.last_physical=deepcopy(out)
        return out

    def requirements(self,data,row,physical):
        out=self.original_requirements(data,row,physical);l=data['industrial']['lunar'];margin=data['industrial']['economics']['capacity_margin']
        def machine(name,flow,rate,mass):out[name]=math.ceil(flow*margin/(rate*self.engine.HOURS*l['process_availability']))*mass
        if self.improved:
            w=self.spec['water_equipment'];wa=water_account(self.spec,row['population'],row.get('visitor_nights_year',0.)/data['constants']['JULIAN_YEAR_DAYS'],True)
            machine('loop_water_equipment',wa['gross_processed_kg'],w['train_kg_hr'],w['hardware_kg'])
            ret=physical.get('loop_structural_retired_kg_year',{})
            # Bulk sorting is an explicit fraction of a throughput-sized bulk
            # train. Metal sorting uses the inherited fabrication-cell proxy.
            surface=ret.get('surface',0.)+ret.get('bulk',0.)*self.spec['structural']['bulk_processing_fraction_of_existing_train']*l['fabrication_cell_kg_hr']/l['bulk_train_kg_hr']
            g=self.spec['goods'];surface+=row['population']*g['gross_kg_person_year']*(g['collection']*g['quality']-g['baseline_recovery'])
            machine('loop_sort_surface',surface,l['fabrication_cell_kg_hr'],l['fabrication_cell_hardware_kg'])
            machine('loop_sort_source',sum(ret.get(k,0.) for k in ('source','depot','leo')),l['fabrication_cell_kg_hr'],l['fabrication_cell_hardware_kg'])
        if self.repair:
            r=self.spec['refurbishment']
            for dest in ('surface','source','depot'):machine('loop_repair_'+dest,physical['loop_repair_input_kg_year'][dest],r['workshop_kg_hr'],r['workshop_hardware_kg'])
            # Depot repair electricity is entirely rejected as heat, whereas
            # the inherited chemical-power law retains some chemical energy.
            kw=physical.get('loop_extra_process_kWh_year',{}).get('depot',0.)*1000/self.engine.HOURS
            p=data['industrial']['power'];z=data['industrial']['source'];flux=p['emissivity']*self.engine.STEFAN_BOLTZMANN*p['radiator_K']**4
            out['depot_radiators']+=kw*z['hydrogen_HHV_kWh_kg']/z['electrolysis_liquefaction_kWh_kg_H2']*margin/flux*p['radiator_kg_m2']
        return out


def material_account(experiment,row,previous_bank):
    e=experiment;year=row['year'];f=e.last_flows[year];phys=row['physical'];ret=e.last_retired[year];cats=e.engine.catalogue(e.data)
    structural={};endbank=phys['recycled_inventory_kg'];cfg=e.spec['structural']
    for dest in (*e.engine.DESTS,'bulk'):
        gross=sum(m for (n,part),m in ret.items() if part=='structure' and ('bulk' if cats[n]['material']=='bulk' else cats[n]['destination'])==dest)
        recovered=f['bulk_recycled'] if dest=='bulk' else f['recycled'][dest]
        collected=gross*cfg['collection'] if e.improved else recovered
        opening=previous_bank.get(dest,0.);closing=endbank[dest];used=opening+recovered-closing
        structural[dest]=dict(retired_kg=gross,collected_kg=collected,recovered_kg=recovered,
                              collection_loss_kg=gross-collected,quality_reject_kg=collected-recovered,
                              opening_inventory_kg=opening,used_kg=used,closing_inventory_kg=closing,
                              mass_balance_residual_kg=opening+gross-(closing+used+gross-recovered))
    g=e.spec['goods'];gross=g['gross_kg_person_year']*row['population'];recovery=g['collection']*g['quality'] if e.improved else g['baseline_recovery']
    collected=gross*g['collection'] if e.improved else gross*recovery
    goods=dict(gross_kg=gross,collected_kg=collected,recovered_kg=gross*recovery,collection_loss_kg=gross-collected,quality_reject_kg=collected-gross*recovery,loss_and_virgin_makeup_kg=gross*(1-recovery),
               imported_carbon_feed_kg=gross*(1-recovery)*g['carbon_share'],local_makeup_kg=gross*(1-recovery)*(1-g['carbon_share']),
               interpretation='Inherited steady annual goods-return assumption; distinct from capital scrap stocks.')
    active_total=sum(m for (n,p),m in ret.items() if p=='active')
    refurbishment_audit=e.last_audit[year] if e.repair else {k:empty_repair() for k in e.engine.DESTS}
    eligible_retired=sum(v['retired_kg'] for v in refurbishment_audit.values())
    return dict(structural=structural,goods=goods,active_retired_total_kg=active_total,
                active_outside_refurbishment_scope_kg=active_total-eligible_retired,refurbishment=refurbishment_audit,
                active_rejects_are_uncredited=True)


def solve_case(spec,data,case,variant):
    e=Experiment(data,spec,variant);result=e.engine.solve(e.data,case['annual'],max_iterations=50)
    seed_audit=dict(physical_annual_average=deepcopy(result['seed']['physical_annual_average']),late_preproduction_physical_total=deepcopy(e.last_physical),reading_rule='The inherited eight-year average plus the separately dated years 6/7 procurement total; the latter is divided by 8 only to present an average. No source deployment chronology is inferred.')
    before=run.compact_case(case,result,'base')
    before=bootstrap.apply(dict(cases=[before],rules=[]),e.data)['cases'][0]
    # Active refurbishment changes procurement, not the mass of installed
    # components. Correct the old gross-purchase diagnostics accordingly.
    if e.repair:
        I=e.data['normalization']['USD2026_to_USD2041'];r=spec['refurbishment'];cats=e.engine.catalogue(e.data)
        for row in before['annual']:
            y=row['year'];f=e.last_flows[y];extra_jobs=sum(f['loop_inspected'].values())/r['inspection_kg_worker_year']+sum(f['loop_repaired'].values())/r['remanufacture_kg_worker_year']
            extra_payroll=extra_jobs*e.data['industrial']['economics']['loaded_staff_USD2026_year']*I
            row['core_resident_jobs']+=extra_jobs;row['technical_resident_jobs']+=extra_jobs
            c=row['cash'];c['resident_payroll_USD2041']+=extra_payroll;c['operating_local_payroll_USD2041']+=extra_payroll;c['operating_total_USD2041']+=extra_payroll
            row['capacity_summary']['resident_workforce']=row['technical_resident_jobs']/max(1.,row['population']*.4)
            row['maximum_capacity_use']=max(row['maximum_capacity_use'],row['capacity_summary']['resident_workforce'])
            if row['maximum_capacity_use']>1+1e-6:
                row['feasible']=False;row['actual_sales_capacity_fraction']=0.;row['shortfalls']['resident_workforce']=row['capacity_summary']['resident_workforce']-1
            row['physical']['source_tug_component_import_kg_year']-=f['loop_import_avoided'].get('asteroid_tugs',0.)
            drivers={}
            for ev in result['component_orders']:
                if ev['component']=='active' and ev['order_year'] in (y,y-1):
                    n=ev['asset'];drivers[n]=drivers.get(n,0.)+ev.get('new_imported_mass_kg',ev['mass_kg'])/2*cats[n]['price']*I
            row['top_imported_component_purchases_USD2041']=dict(sorted(drivers.items(),key=lambda kv:kv[1],reverse=True)[:8])
            row['direct_refurbishment_jobs']=extra_jobs;row['direct_refurbishment_payroll_USD2041']=extra_payroll
        before['first_capacity_shortfall']=next((r['year'] for r in before['annual'] if not r['feasible']),None)
    rows=[];bank={}
    for row in before['annual']:
        request=case['annual'][row['year']-case['annual'][0]['year']]
        water=water_account(e.spec,row['population'],request.get('visitor_nights_year',0.)/e.data['constants']['JULIAN_YEAR_DAYS'],e.improved)
        stock=row['physical']['source_water_breakdown_kg_year']['resident_and_guest_water']-water['operational_makeup_kg']
        if stock < -0.001:raise ValueError('Water makeup exceeds the physical delivered water')
        water['rounding_residual_kg']=min(0.,stock);water['new_stock_and_construction_kg']=max(0.,stock)
        water['landscape_agriculture_allowance_basis']='Inherited 5 kg/person/day for garden, lake and food-production losses; retained beside explicit farm balances because its decomposition is unknown. Potential overlap is unquantified.'
        water['metabolic_oxygen_recovery']='Not separately resolved: inherited air demand is atmosphere leakage plus new stocks, with food/plant metabolism not represented as an ISS oxygen loop. No additional oxygen credit.'
        materials=material_account(e,row,bank);bank=row['physical']['recycled_inventory_kg']
        rows.append(dict(year=row['year'],phase='operation',water=water,materials=materials,
                         source_feed_kg=row['physical']['source_feed_kg_year'],lunar_raw_kg=row['physical']['lunar_raw_kg_year'],
                         source_water_equivalent_kg=row['physical']['source_water_kg_year'],source_nitrogen_kg=row['physical']['source_nitrogen_kg_year'],source_virgin_metal_kg=row['physical']['source_virgin_metal_kg_year'],source_water_breakdown_kg=row['physical']['source_water_breakdown_kg_year'],
                         depot_water_breakdown_kg=row['physical']['depot_water_breakdown_kg_year'],
                         chemical_exhaust_recovery_kg=0.,extra_process_kWh=row['physical'].get('loop_extra_process_kWh_year',{})))
    def seed_average(key):return seed_audit['physical_annual_average'][key]+seed_audit['late_preproduction_physical_total'][key]/8
    seedrows=[dict(year=y,phase='preproduction_average',source_water_equivalent_kg=seed_average('source_water_kg_year'),source_nitrogen_kg=seed_average('source_nitrogen_kg_year'),source_virgin_metal_kg=seed_average('source_virgin_metal_kg_year'),source_feed_kg=before['seed']['source_feed_kg']/8,lunar_raw_kg=before['seed']['lunar_raw_kg']/8,
                   water=None,materials=None,chemical_exhaust_recovery_kg=0.,reading_rule='One eighth of the inherited combined seed feed. Gross circulation, recovered products and deployment losses lack a dated seed sequence; null is not zero.') for y in range(8)]
    return before,dict(name=case['name']+':'+variant,baseline_case=case['name'],variant=variant,portfolio=case['portfolio'],market_growth_name=case['market_growth_name'],
                       legacy_landscape_allowance_kg_person_day=0. if variant=='combined_without_residual_allowance' else 5.,
                       physical_case=before,seed_audit=seed_audit,annual=seedrows+rows,
                       cumulative_source_feed_kg=sum(r['source_feed_kg'] for r in seedrows+rows),cumulative_lunar_raw_kg=sum(r['lunar_raw_kg'] for r in seedrows+rows))


def product(spec,data,book,baseline,cases=None,variants=None):
    controls={c['name']:c for c in baseline['cases']};out=[]
    cases=cases or [c for c in book['cases'] if c['portfolio'] in spec['portfolios'] and c['market']==spec['market'] and c['market_growth_name'] in spec['growth']]
    variants=variants or spec['variants'];cache=ROOT/'research/runs/enclosed_industry/reclamation';cache.mkdir(parents=True,exist_ok=True)
    for case in cases:
        for variant in variants:
            identity=hashlib.sha256(json.dumps(dict(spec=spec,producer=sha(__file__),case=case,variant=variant),sort_keys=True).encode()).hexdigest()
            target=cache/(identity+'.json.gz');start=time.monotonic()
            if target.exists():value=json.loads(gzip.decompress(target.read_bytes()))
            else:
                physical,value=solve_case(spec,data,case,variant)
                if variant=='baseline':
                    if physical!=controls[case['name']]:raise ValueError('Baseline control mismatch: '+case['name'])
                    value['baseline_exact_match']=True
                target.write_bytes(gzip.compress(json.dumps(value,separators=(',',':')).encode(),mtime=0))
            out.append(value);print(case['name'],variant,round(time.monotonic()-start,2),'seconds',flush=True)
    return dict(schema=SCHEMA,producer=dict(path=str(Path(__file__).relative_to(ROOT)),sha256=sha(__file__)),
                dependencies={**spec['dependencies'],str((HERE/'reclamation_inputs.json').relative_to(ROOT)):sha(HERE/'reclamation_inputs.json'),str((HERE/'reclamation_sources.json').relative_to(ROOT)):sha(HERE/'reclamation_sources.json')},
                evidence=spec['evidence'],units='kg, kWh, W, persons and constant USD2041; years 0..220; no U conversion.',
                reading_rule='Each physical_case contains complete priced annual accounts and capacity screens under an unchanged orderbook. New loop accounts audit subsets of those physical totals; do not add their gross or recovered masses to virgin demand. Seed annual material detail remains null where the source only supports an eight-year aggregate. A failed workforce screen retains requirement costs but admits no qualified production.',
                assumptions=spec['assumptions'],scenario_parameters={k:spec[k] for k in ('water','landscape_agriculture_allowance_kg_year','goods','structural','water_equipment','refurbishment')},cases=out)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--case',action='append');ap.add_argument('--variant',action='append');ap.add_argument('--output',type=Path,default=HERE/'results/reclamation.json.gz');args=ap.parse_args()
    spec,data,book,baseline=load();cases=[c for c in book['cases'] if c['name'] in args.case] if args.case else None
    result=product(spec,data,book,baseline,cases,args.variant)
    temp=ROOT/'research/runs/enclosed_industry/reclamation-product.json';temp.write_text(json.dumps(result,separators=(',',':'),allow_nan=False)+'\n')
    manifest=package.package(temp,args.output);print(json.dumps(manifest,indent=2))


if __name__=='__main__':main()
