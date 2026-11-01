"""Source-bound short campaigns for enduring structural infrastructure."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS
from shared.provenance import ROOT, constants_used
from .early_campaign import campaign

HERE=Path(__file__).resolve().parent
INPUT_PATH=HERE/'early_campaigns_inputs.json'


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_inputs(path=INPUT_PATH,root=ROOT):
    s=json.loads(Path(path).read_text())
    if s['schema']!='terluna.research.early-campaigns-inputs/1':raise ValueError('Unexpected campaign inputs')
    products={}
    for name,spec in s['products'].items():
        p=Path(root)/name
        if digest(p)!=spec['sha256']:raise ValueError(f'Input changed; review and rebind: {name}')
        d=json.loads(p.read_text())
        if d['schema']!=spec['schema']:raise ValueError(f'Unexpected product schema: {name}')
        products[name]=d
    return s,products


def results():
    s,p=load_inputs()
    source=p['research/studies/world_economy/industrial_markets_inputs.json']
    market=p['research/studies/world_economy/results/industrial_markets.json']
    f=next(x for x in source['families'] if x['id']=='structures')
    plant=next(x for x in source['plants'] if x['id']=='industrial')
    capital=source['capital'];ref=s['reference'];rows=[];reference=None
    shares=('variable_intermediate_share_of_benchmark_cost','variable_compensation_share_of_benchmark_cost','fixed_maintenance_share_of_benchmark_cost','capital_service_share_of_benchmark_cost')
    if tuple(capital[k] for k in shares)!=(.55,.15,.05,.25):
        raise ValueError('Campaign cost allocation changed; review helper and rebind')
    eligible=market['anchors']['array_structure_kg']*source['buyers']['array']['eligible_share_of_final_structure']
    for raw in s['raw_capacities_kg_year']:
        scale=raw/plant['gross_nameplate_kg_year']
        for name in s['cost_cases']:
            case=source['cost_cases'][name];base=f['base_production_cost_U_kg']
            freight=source['delivery_U_kg']['moon']['array']
            alternative=min(base*case[w]+source['delivery_U_kg'][w]['array'] for w in ('earth','array'))
            own=base*case['moon']+freight
            if own>=alternative:raise ValueError('Campaign requires a competitive source offer')
            price=base*case['moon']+source['price_share_of_cost_advantage_to_seller']*(alternative-own)
            for years in s['years']:
                for budget in s['buyer_budgets_U_year']:
                    for advance in s['customer_advance_fractions']:
                        result=campaign(raw_nameplate_kg_year=raw,years=years,eligible_orders_kg=eligible,
                            buyer_budget_U_year=budget,base_production_cost_U_kg=base,lunar_cost_factor=case['moon'],
                            producer_price_U_kg=price,freight_U_kg=freight,utilization=source['utilization'],accepted_yield=f['yield'],
                            power_W=plant['power_W']*scale,outbound_capacity_kg_year=plant['outbound_route_kg_year']*scale,
                            seconds_per_year=JULIAN_DAY*JULIAN_YEAR_DAYS,process_energy_J_raw_kg=f['factory_energy_MJ_per_feed_kg']*1e6,
                            local_feed_fraction=f['local_feed_fraction'],available_local_feed_kg=source['local_feed_allowance_kg'],
                            economic_rate=capital['economic_real_rate'],capital_life_years=capital['life_years'],
                            terminal_recovery_fraction=capital['closure_salvage_fraction_of_remaining_book'],
                            decommission_fraction=capital['closure_cost_fraction_of_initial_capital'],
                            customer_advance_fraction=advance,independently_funded_fraction=s['independently_funded_fraction'],ramp=s['ramp'])
                        tags=dict(raw_capacity_kg_year=raw,years=years,buyer_budget_U_year=budget,cost_case=name,advance=advance)
                        result.update(tags);result['producer_price_U_kg']=price
                        if tags==ref:reference=result
                        rows.append({k:v for k,v in result.items() if k!='annual'})
    assert reference is not None
    names=('early_campaigns.py','early_campaigns_inputs.json','early_campaign.py')
    return dict(schema='terluna.research.early-campaigns/1',
        producer=dict(files={str((HERE/n).relative_to(ROOT)):digest(HERE/n) for n in names},constants=constants_used([__file__]),input_products=s['products']),
        evidence=s['evidence'],eligible_order_kg=eligible,reference=reference,scenarios=rows,
        cost_boundary='Campaign helper uses the source cost shares55% intermediates,15% compensation,5% fixed and25% capital service; actual capital is counted once, not again as an annual capital-service expense.',
        boundaries=s['boundaries'])


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'research/runs/world_economy/early_campaigns.json')
    a=p.parse_args();a.output.parent.mkdir(parents=True,exist_ok=True);a.output.write_text(json.dumps(results(),indent=2,allow_nan=False)+'\n');print(a.output)


if __name__=='__main__':main()
