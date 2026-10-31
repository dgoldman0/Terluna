"""Durable stock and operating-asset value after natural vacuum is lost."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
from shared.provenance import ROOT
from . import industrial_stocks as stocks

HERE=Path(__file__).resolve().parent
INPUT_PATH=HERE/'industrial_legacy_inputs.json'


def digest(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_inputs(path=INPUT_PATH, root=ROOT):
    s=json.loads(Path(path).read_text())
    if s['schema']!='terluna.research.industrial-legacy-inputs/1':raise ValueError('Unexpected legacy input schema')
    products={}
    for name,spec in s['products'].items():
        p=Path(root)/name
        if digest(p)!=spec['sha256']:raise ValueError(f'Input changed; review and rebind: {name}')
        d=json.loads(p.read_text())
        if d['schema']!=spec['schema']:raise ValueError(f'Unexpected source schema: {name}')
        products[name]=d
    return s,products


def summarize(result):
    return dict(cheapest_tested_reserve_years=result['cheapest_tested_reserve_years'],
        cases=[dict(reserve_years=r['reserve_years'],initial_stock_kg=r['result']['settings']['initial_stock_kg'],
                    totals=r['result']['totals']) for r in result['cases']])


def results():
    s,products=load_inputs();industry=products['resources/results/industry.json']
    market=products['research/studies/world_economy/results/industrial_markets.json']
    reference={};scenarios=[]
    for kind,decay in s['classes'].items():
        for future in s['future_buy_costs_per_kg']:
            for decline in s['annual_real_price_declines']:
                for rate in s['real_discount_rates']:
                    r=stocks.evaluate_reserve_horizons(1,s['reserve_years'],s['analysis_years'],
                        compatibility_decay=decay,future_buy_cost_per_kg=future,
                        annual_price_decline=decline,real_discount_rate=rate,**s['normalized_costs'])
                    row=dict(stock_class=kind,future_buy_cost_per_kg=future,annual_price_decline=decline,
                             real_discount_rate=rate,**summarize(r))
                    scenarios.append(row)
                    if all(row[k]==v for k,v in s['reference'].items()):reference[kind]=row
    rates=industry['plant_traffic']['upkeep_gt_yr'];feed=[]
    physical_loss=s['normalized_costs']['physical_attrition']
    for renewal in rates:
        for fresh in s['fresh_feed_fractions']:
            annual=renewal*1e12*fresh
            initial=stocks.required_initial_stock(annual,s['analysis_years'],physical_loss)
            feed.append(dict(gross_renewal_Gt_year=renewal,fresh_fraction=fresh,
                annual_virgin_feed_kg=annual,loss_free_50_year_kg=annual*s['analysis_years'],
                attrition_adjusted_initial_kg=initial,
                equivalent_years_of_reference_finished_product_mass=initial/market['reference']['normal_saleable_capacity_kg_year'],
                rule='Mass comparison only: a separate grade-matched feed refinery, recovery loop and prepositioning budget are required; no spare capacity from the finished-product plant is inferred.'))
    a=s['installed_operating_assets'];asset=[]
    for life in a['useful_lives_years']:
        for rate in a['gross_operating_surplus_rates_after_current_expenses']:
            capital=a['complete_asset_value_U'];cfc=capital/life
            net=capital*rate-cfc
            asset.append(dict(complete_asset_value_U=capital,useful_life_years=life,
                gross_operating_surplus_rate=rate,gross_operating_surplus_U_year=capital*rate,
                cfc_and_replacement_U_year=cfc,net_operating_surplus_before_financing_U_year=net,
                lunar_owner_income_before_financing_U_year=net*a['lunar_ownership_fraction'],
                interpretation='Maintained operating asset with annual replacement funded; assumed surplus after current expenses. Array production and lunar resident income remain distinct. No guaranteed distribution or free initial asset.'))
    c=s['construction_package'];construction=[]
    for eligible in c['eligible_shares']:
        for saving in c['eligible_package_saving_fractions']:
            for extra in c['extra_enclosure_transition_cost_fractions_of_total']:
                gross=c['programme_cost_U']*eligible*saving
                penalty=c['programme_cost_U']*extra
                construction.append(dict(programme_cost_U=c['programme_cost_U'],eligible_share=eligible,
                    eligible_package_saving_fraction=saving,gross_avoided_cost_U=gross,
                    extra_enclosure_transition_cost_U=penalty,net_cost_reduction_U=gross-penalty,
                    independent_export_receipt_U=0,rule=c['rule']))
    names=('industrial_legacy.py','industrial_legacy_inputs.json','industrial_stocks.py')
    return dict(schema='terluna.research.industrial-legacy/1',
        producer=dict(files={str((HERE/n).relative_to(ROOT)):digest(HERE/n) for n in names},input_products=s['products']),
        evidence=s['evidence'],normalized_reference=reference,normalized_scenarios=scenarios,
        film_fresh_feed_reserves=feed,installed_asset_examples=asset,construction_package_sensitivities=construction,boundaries=s['boundaries'])


def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'research/runs/world_economy/industrial_legacy.json')
    a=p.parse_args();a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(results(),indent=2,allow_nan=False)+'\n');print(a.output)


if __name__=='__main__':main()
