"""Reproduce the three-world construction, production and financing comparison."""
from __future__ import annotations
import argparse
import copy
import json
from pathlib import Path
from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS
from shared.provenance import constants_used
from . import annual, construction, finance, inventory, productivity

HERE = Path(__file__).resolve().parent
INPUT_PATH = HERE / 'comparison_inputs.json'


def compact(row):
    measures = ('gross_world_product','ordinary_cfc','engineered_cfc','natural_depletion',
                'net_world_product','net_resident_income','cash_resources_after_principal',
                'essential_final_operations','resources_after_operations_and_principal')
    return dict(cost_case=row['cost_case'], water_case=row['water_case'],
        finance_case=row['finance_case'], applied_parameters=row['applied_parameters'],
        worlds={w:{k:v[k] for k in measures} for w,v in row['accounts']['worlds'].items()},
        ratios_to_earth=row['ratios_to_earth'],
        required_net_exports_share={w:v['required_net_exports_share_of_gwp'] for w,v in row['trade_requirements'].items()},
        nonnegative_cash_for_all=all(v['cash_resources_after_principal']>=0 for v in row['accounts']['worlds'].values()))


def operating_memo(inv, cost_inputs, prices):
    """Nonadditive funding examples, never additional deductions from GVA."""
    def get(path,key): return inv['sources'][path]['values'][key]
    hours=JULIAN_YEAR_DAYS*JULIAN_DAY/3600
    power=get('provisioning/results/power.json','reserves.ring_fault.delivered_total_tw')
    films=get('research/studies/solar_shield_array/results/array_heat.json','fleet.film_upkeep_Gt_per_year')
    p=get('research/studies/food_systems/results/food_blend.json','blend.8.wild_30_hybrid_20.phosphorus_mt_year')['lost']
    sky=get('research/studies/ways_of_living/results/ways_of_living.json','aerial_reference.hydrogen_upkeep')
    losses=sky['loss_bracket_per_year']
    h2_tw=[inv['population_billion']['moon_aloft'] * losses[i]*100 * sky['per_billion_per_percent_a_year']['remake_tw'][i] for i in (0,1)]
    rows={}
    for case,inputs in cost_inputs['cost_cases'].items():
        electricity=prices['delivered_electricity_U_per_MWh'][case]
        rows[case]=dict(
            delivered_electricity_U_year=[v*1e6*hours*electricity for v in power],
            lifting_hydrogen_energy_subset_U_year=[v*1e6*hours*electricity for v in h2_tw],
            shield_film_replacement_U_year=[v*1e12*inputs['shield_fabrication_U_per_kg'] for v in films],
            phosphorus_loss_replenishment_U_year=[v*1e9*prices['phosphorus_U_per_kg'][case] for v in p])
    return dict(evidence='Price and funding sensitivities on physical ranges, not necessary expenditure estimates or additive macroeconomic deductions.',
        prices=prices, cases=rows,
        overlap_rules=['Hydrogen electricity is a subset of total delivered electricity.',
            'Film replacement is the same service approximated by optical capital depreciation; do not add both.',
            'Producer electricity and nutrient purchases are already intermediate deductions in sector GVA.',
            'The collective-service allocation is a separate policy sensitivity, not the sum of these commercial bills.'])


def results():
    settings=json.loads(INPUT_PATH.read_text())
    if settings['schema']!='terluna.research.world-economy-comparison-inputs/1':
        raise ValueError('Unexpected comparison input schema')
    inv=inventory.results()
    costs=construction.calculate()
    ci=construction.load_inputs()
    pi=productivity.load_inputs(HERE/'productivity_inputs.json')
    ai=json.loads((HERE/'annual_inputs.json').read_text())
    pop=inv['population_billion']
    def production(scenario,income=100000,weight='central'):
        return productivity.evaluate_productivity(pi,earth_population=pop['earth']*1e9,
            moon_population=pop['moon']*1e9,array_population=pop['array']*1e9,
            earth_gdp_per_resident=income,scenario_id=scenario,weight_case=weight)
    ref=settings['reference']
    cost=next(c for c in costs['cases'] if c['cost_case']==ref['cost_case'] and c['water_case']==ref['water_case'])
    prod=production(ref['production_case'],ref['earth_gdp_per_resident'])
    def evaluate(c=cost,p=prod,fin='legacy_retired',**overrides):
        return annual.evaluate(c,p,ai,finance_case=fin,**overrides)
    productions={s:production(s,ref['earth_gdp_per_resident']) for s in pi['scenarios']}
    mature={s:evaluate(p=p) for s,p in productions.items()}
    financing={f:evaluate(fin=f) for f in ai['finance_cases']}
    cost_grid=[]
    for c in costs['cases']:
        r=compact(evaluate(c=c,fin='patient'))
        r['construction_in_earth_output_years']=c['combined_total_U']/prod['worlds']['earth']['gross_product']
        cost_grid.append(r)
    independent=[]
    for c in costs['cases']:
        if c['water_case']!=ref['water_case']: continue
        for s,p in productions.items():
            independent.append(dict(production_case=s,**compact(evaluate(c=c,p=p))))
    income=[]
    for x in settings['earth_gdp_per_resident_sensitivity']:
        income.append(dict(earth_gdp_per_resident=x,**compact(evaluate(p=production(ref['production_case'],x),fin='patient'))))
    sensitivities=dict(
        cost_and_water_under_patient_finance=cost_grid,
        independent_cost_and_production=independent,
        income_under_patient_finance=income,
        ordinary_depreciation=[compact(evaluate(ordinary_cfc_share=v)) for v in ai['ordinary_cfc_sensitivities']],
        incremental_asset_lives=[compact(evaluate(asset_life_scale=v)) for v in settings['asset_life_scales']],
        collective_operating_allocation=[compact(evaluate(operation_scale=v)) for v in ai['essential_final_operations']['scale_sensitivities']],
        depletion={name:compact(evaluate(natural_depletion_shares=v)) for name,v in ai['natural_depletion_sensitivities'].items()},
        foreign_creditor_share=[compact(evaluate(fin='patient',creditor_share=v)) for v in ai['creditor_share_sensitivities']],
        earth_inherited_external_debt=[compact(evaluate(earth_legacy_debt_multiple=v)) for v in ai['earth_legacy_debt']['multiple_sensitivities']],
        sector_weights={w:compact(evaluate(p=production(ref['production_case'],ref['earth_gdp_per_resident'],w))) for w in pi['weight_cases']})
    sensitivities['ordinary_portfolio_differences']={}
    for name,shares in settings['ordinary_portfolio_sensitivities'].items():
        varied=copy.deepcopy(ai)
        varied['ordinary_cfc_shares']=shares
        sensitivities['ordinary_portfolio_differences'][name]=compact(annual.evaluate(cost,prod,varied))
    sensitivities['construction_financed_fraction']=[]
    for fraction in (0,.25,.5,1):
        varied=copy.deepcopy(ai)
        varied['finance_cases']['patient']['financed_fraction']=fraction
        result=annual.evaluate(cost,prod,varied,finance_case='patient')
        sensitivities['construction_financed_fraction'].append(dict(financed_fraction=fraction,**compact(result)))
    # Show ownership income without relocating the array's production.
    array_nwp=mature[ref['production_case']]['accounts']['worlds']['array']['net_world_product']
    sensitivities['array_ownership_income_example']=compact(evaluate(additional_transfers=[
        dict(payer='array',recipient='moon',kind='profit',amount=.1*array_nwp)]))
    rate_grid=[dict(real_rate=r, payment_share_by_term={str(n):finance.annuity_payment(1,r,n)
                 for n in (100,250,500,1000)},
                 capitalized_500_year_build_multiple=finance.construction_finance(1,500,r)['debt_at_completion'])
               for r in (0,.001,.005,.01,.02)]
    names=('compare.py','comparison_inputs.json','annual.py','annual_inputs.json','accounts.py',
           'finance.py','construction.py','construction_inputs.json','productivity.py','productivity_inputs.json',
           'inventory.py','inputs.json')
    return dict(schema='terluna.research.world-economy-comparison/1',
        producer=dict(study='world_economy',files={str((HERE/n).relative_to(inventory.ROOT)):inventory.digest(HERE/n) for n in names},
                      constants=constants_used([__file__,HERE/'construction.py']),
                      input_products=json.loads((HERE/'inputs.json').read_text())['products']),
        evidence='Conditional economic scenarios using physical screens, assumed future prices, sector value added and asset lives. Numerical identities are checked; no market equilibrium, validated capacity or probability bounds are established.',
        reading_rule='U is one common transferable real unit. All products are annual. Earth=100 comparisons use the corresponding Earth measure, not a universal denominator. Negative debtor resources indicate infeasible claims; creditor receipts in those cases are not achievable consumption. GWP excludes intermediate inputs already; operations and debt are not deducted twice. Depletion is unpriced in the reference and varied explicitly.',
        population_billion=pop,reference=ref,production=productions,
        mature_accounts=mature,reference_financing=financing,sensitivities=sensitivities,
        financing_arithmetic=rate_grid,operating_funding_memo=operating_memo(inv,ci,settings['operating_prices']),
        limitations=['301K redesign not propagated through every inherited physical requirement.',
            'Construction and production scenario axes are independent, not matched probability percentiles.',
            'Construction-phase renewal and operations are funded separately; no transition production path is solved.',
            'The array sector allocation is bookkeeping; only aggregate array output is a scenario.',
            'No bilateral trade matrix, export demand, institution, hazard insurance or future exchange-rate system is solved.',
            'The optional ownership example transfers 10% of array NWP as assumed profits; no actual profit share or equity valuation is inferred.'])


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=inventory.ROOT/'research/runs/world_economy/comparison.json')
    args=parser.parse_args()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(results(),indent=2)+'\n')
    print(args.output)

if __name__=='__main__': main()
