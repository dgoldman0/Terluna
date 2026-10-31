"""Conditional lunar factory accounts with competitive buyers and finite demand.

No production-cost factor is inferred from gravity. Prices, qualification,
capacity and financing are scenarios; existing macroeconomic products stay fixed.
"""
from __future__ import annotations
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
from shared.constants import JULIAN_DAY, JULIAN_YEAR_DAYS
from shared.provenance import ROOT, constants_used
from . import finance

HERE = Path(__file__).resolve().parent
INPUT_PATH = HERE / 'industrial_markets_inputs.json'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_inputs(path=INPUT_PATH, root=ROOT):
    s = json.loads(Path(path).read_text())
    if s['schema'] != 'terluna.research.industrial-markets-inputs/1':
        raise ValueError('Unexpected industrial market inputs')
    products = {}
    for name, spec in s['products'].items():
        p = Path(root) / name
        if digest(p) != spec['sha256']:
            raise ValueError(f'Input changed; review and rebind: {name}')
        d = json.loads(p.read_text())
        if d['schema'] != spec['schema']:
            raise ValueError(f'Unexpected product schema: {name}')
        products[name] = d
    if not math.isclose(sum(f['mass_share'] for f in s['families']), 1):
        raise ValueError('Product shares must sum to one')
    c = s['capital']
    shares = ['fixed_maintenance_share_of_benchmark_cost', 'capital_service_share_of_benchmark_cost',
              'variable_intermediate_share_of_benchmark_cost', 'variable_compensation_share_of_benchmark_cost']
    if not math.isclose(sum(c[k] for k in shares), 1):
        raise ValueError('Benchmark cost shares must sum to one')
    return s, products


def anchors(products):
    construction = products['research/studies/world_economy/results/construction.json']
    case = next(c for c in construction['cases'] if c['cost_case'] == 'advanced_industrial' and c['water_case'] == 'surface_only')
    industry = products['resources/results/industry.json']
    structure = case['included_array_structure_comparison']
    return dict(array_structure_kg=structure['mass_kg'],
                array_structure_value_U=structure['separately_valued_U'],
                structure_unit_value_U_kg=structure['separately_valued_U']/structure['mass_kg'],
                gross_film_renewal_Gt_year=industry['plant_traffic']['upkeep_gt_yr'],
                ring_delivery_ideal_MJ_kg=industry['places']['moon_surface']['to_ring_orbit_mj_kg'])


def offer(s, family, case, buyer, *, after_atmosphere=False):
    """Identical-quality price below the cheaper optimized competitor's offer."""
    base = family['base_production_cost_U_kg']
    lunar_cost = base * case['moon']
    freight = s['delivery_U_kg']['moon'][buyer]
    if after_atmosphere:
        lunar_cost += s['post_atmosphere']['added_production_cost_U_kg'][family['id']]
        freight = s['post_atmosphere']['delivery_from_moon_U_kg'][buyer]
    alternatives = {w: base*case[w]+s['delivery_U_kg'][w][buyer] for w in ('earth','array')}
    competitor = min(alternatives.values())
    own = lunar_cost + freight
    advantage = competitor-own
    wins = advantage > 1e-12
    price = own + s['price_share_of_cost_advantage_to_seller']*max(0,advantage)
    return dict(delivered_lunar_cost_U_kg=own, producer_cost_U_kg=lunar_cost,
                competitor_costs_U_kg=alternatives, best_competitor_cost_U_kg=competitor,
                cost_advantage_U_kg=advantage, cost_competitive=wins,
                offered_delivered_price_U_kg=price, delivery_U_kg=freight,
                producer_receipt_U_kg=price-freight,
                remaining_allowance_for_added_production_cost_U_kg=advantage)


def run_case(s, bound, plant, case, years, *, keep_annual=False, independent_first=False):
    if not isinstance(years,int) or years < 1:
        raise ValueError('Positive integer operating years required')
    c = s['capital']; families = {f['id']:f for f in s['families']}
    normal = {k: plant['gross_nameplate_kg_year']*s['utilization']*f['mass_share']*f['yield'] for k,f in families.items()}
    energy_at_normal = sum(normal[k]/f['yield']*f['factory_energy_MJ_per_feed_kg']*1e6 for k,f in families.items())
    energy_budget = plant['power_W']*JULIAN_DAY*JULIAN_YEAR_DAYS
    power_factor = min(1,energy_budget/energy_at_normal)
    route_factor = min(1,plant['outbound_route_kg_year']/sum(normal.values()))
    crf = finance.annuity_payment(1,c['economic_real_rate'],c['life_years'])
    costs = {k:f['base_production_cost_U_kg']*case['moon'] for k,f in families.items()}
    capital = sum(normal[k]*costs[k]*c['capital_service_share_of_benchmark_cost']/crf for k in families)
    fixed = sum(normal[k]*costs[k]*c['fixed_maintenance_share_of_benchmark_cost'] for k in families)
    annual_replacement_reserve = capital/c['life_years']
    loan = capital*c['financed_fraction']
    annual_loan = finance.annuity_payment(loan,c['loan_real_rate'],c['loan_years'])
    closing_loan = finance.loan_balance_after_payments(loan,c['loan_real_rate'],c['loan_years'],years)
    age_at_close = years - ((years-1)//c['life_years'])*c['life_years']
    expected_reserve = age_at_close*annual_replacement_reserve
    remaining_book = capital*max(0,1-age_at_close/c['life_years'])
    salvage = remaining_book*c['closure_salvage_fraction_of_remaining_book']
    decommission = capital*c['closure_cost_fraction_of_initial_capital']
    closure_escrow_year = max(0,closing_loan+decommission-expected_reserve-salvage)/years
    offers = {k:{b:offer(s,f,case,b) for b in s['buyers']} for k,f in families.items()}
    stock = {b:{} for b in ('array','colonies')}
    for k,f in families.items():
        stock['array'][k] = bound['array_structure_kg']*s['buyers']['array']['eligible_share_of_final_structure']*f['mass_share']
        stock['colonies'][k] = s['buyers']['colonies']['eligible_stock_kg']*f['mass_share']
    initial_stock = copy.deepcopy(stock)
    budgets = {b:v['population']*v['gwp_per_resident_U']*v['procurement_budget_share'] for b,v in s['buyers'].items()}
    reserve = escrow = retained = 0.0
    feed_left = s['local_feed_allowance_kg']
    totals = {k:0.0 for k in ('produced_kg','local_feed_kg','imported_feed_kg','purchased_energy_J','delivered_sales_U','producer_receipts_U','freight_excluded_U','gva_U','nva_before_depletion_U','cash_assigned_U','independent_cash_assigned_U','unfunded_cash_gap_U','retained_cash_U','replacement_capex_U','loan_interest_U','loan_principal_U')}
    buyer_totals = {b:dict(delivered_kg=0.0,spending_U=0.0) for b in s['buyers']}
    annual = []
    unlevered_npv = -capital
    for year in range(1,years+1):
        capacity = {k:q*min(power_factor,route_factor) for k,q in normal.items()}
        local_at_capacity = sum(capacity[k]/f['yield']*f['local_feed_fraction'] for k,f in families.items())
        feed_factor = min(1,feed_left/local_at_capacity) if local_at_capacity else 1
        capacity = {k:q*feed_factor for k,q in capacity.items()}
        orders = {b:{} for b in s['buyers']}
        for b in orders:
            for k,f in families.items():
                requested = f['earth_order_kg_year'] if b=='earth' else min(stock[b][k],initial_stock[b][k]/s['buyers'][b]['procurement_years'])
                orders[b][k] = requested if offers[k][b]['cost_competitive'] else 0.0
            reserved = sum(q*offers[k][b]['offered_delivered_price_U_kg'] for k,q in orders[b].items())
            ratio = min(1,budgets[b]/reserved) if reserved else 1
            orders[b] = {k:q*ratio for k,q in orders[b].items()}
        allocations = {b:{k:0.0 for k in families} for b in s['buyers']}
        for k in families:
            remaining = capacity[k]
            # Product capacities are separate. Fixed buyer budget reservations are
            # not reallocated after this pass; this is a procurement screen.
            ranked = sorted(s['buyers'],key=lambda b:(s['buyers'][b]['independent_funding_fraction'] if independent_first else 0, offers[k][b]['producer_receipt_U_kg']-costs[k]),reverse=True)
            for b in ranked:
                q = min(remaining,orders[b][k]); allocations[b][k]=q; remaining -= q
                if b in stock: stock[b][k] -= q
        made = {k:sum(a[k] for a in allocations.values()) for k in families}
        local_feed = sum(made[k]/f['yield']*f['local_feed_fraction'] for k,f in families.items())
        imported_feed = sum(made[k]/f['yield']*(1-f['local_feed_fraction']) for k,f in families.items())
        energy = sum(made[k]/f['yield']*f['factory_energy_MJ_per_feed_kg']*1e6 for k,f in families.items())
        feed_left = max(0,feed_left-local_feed)
        sales = {b:sum(q*offers[k][b]['offered_delivered_price_U_kg'] for k,q in a.items()) for b,a in allocations.items()}
        exworks = {b:sum(q*offers[k][b]['producer_receipt_U_kg'] for k,q in a.items()) for b,a in allocations.items()}
        revenue = sum(exworks.values()); delivered = sum(sales.values())
        inputs = sum(made[k]*costs[k]*c['variable_intermediate_share_of_benchmark_cost'] for k in families)
        compensation = sum(made[k]*costs[k]*c['variable_compensation_share_of_benchmark_cost'] for k in families)
        payment = annual_loan if year<=c['loan_years'] else 0
        opening = finance.loan_balance_after_payments(loan,c['loan_real_rate'],c['loan_years'],year-1)
        interest = opening*c['loan_real_rate']; principal = payment-interest
        reserve += annual_replacement_reserve; escrow += closure_escrow_year
        actual_replacement = capital if year % c['life_years']==0 and year<years else 0
        reserve -= actual_replacement
        cash = revenue-inputs-compensation-fixed-annual_replacement_reserve-payment-closure_escrow_year
        assigned = max(0,cash)*c['development_share_of_positive_cash']
        retained_year = max(0,cash)-assigned; retained += retained_year
        gap = max(0,-cash)
        independent_share = sum(exworks[b]*s['buyers'][b]['independent_funding_fraction'] for b in exworks)/revenue if revenue else 0
        gva = revenue-inputs-fixed*c['fixed_intermediate_fraction']
        nva = gva-annual_replacement_reserve
        row = dict(year=year,produced_kg=sum(made.values()),local_feed_kg=local_feed,
            imported_feed_kg=imported_feed,purchased_energy_J=energy,delivered_sales_U=delivered,
            producer_receipts_U=revenue,freight_excluded_U=delivered-revenue,
            gva_U=gva,nva_before_depletion_U=nva,cash_assigned_U=assigned,
            independent_cash_assigned_U=assigned*independent_share,unfunded_cash_gap_U=gap,
            retained_cash_U=retained_year,replacement_capex_U=actual_replacement,
            loan_interest_U=interest,loan_principal_U=principal,
            replacement_reserve_U=reserve,closure_escrow_U=escrow,
            made_by_family_kg=made,buyer_allocations_kg=allocations,buyer_spending_U=sales,
            budget_unused_U={b:budgets[b]-sales[b] for b in budgets},
            cash_closure_U=revenue+gap-inputs-compensation-fixed-annual_replacement_reserve-payment-closure_escrow_year-assigned-retained_year)
        for key in totals: totals[key] += row[key]
        for b in buyer_totals:
            buyer_totals[b]['delivered_kg'] += sum(allocations[b].values())
            buyer_totals[b]['spending_U'] += sales[b]
        unlevered_npv += (revenue-inputs-compensation-fixed-actual_replacement)/(1+c['economic_real_rate'])**year
        if keep_annual: annual.append(row)
    unlevered_npv += (salvage-decommission)/(1+c['economic_real_rate'])**years
    terminal_cash = reserve+escrow+salvage+retained-closing_loan-decommission
    return dict(years=years,plant=plant['id'],normal_saleable_capacity_kg_year=sum(normal.values()),
        capacity_energy_J_year=energy_at_normal,power_capacity_factor=power_factor, route_capacity_factor=route_factor,
        assumed_outbound_route_kg_year=plant['outbound_route_kg_year'],
        actual_normal_capacity_utilization=totals['produced_kg']/years/sum(normal.values()),
        allocation_policy='independent_buyers_first' if independent_first else 'highest_producer_margin_first',
        initial_capital_U=capital,initial_borrowing_U=loan,initial_sponsor_equity_U=capital-loan,
        fixed_annual_operations_U=fixed,annual_replacement_reserve_U=annual_replacement_reserve,
        annual_closure_escrow_U=closure_escrow_year,offers=offers,buyer_budgets_U_year=budgets,
        totals=totals,annual_average={k:v/years for k,v in totals.items()},buyer_totals=buyer_totals,
        initial_eligible_buyer_stock_kg=initial_stock,remaining_eligible_buyer_stock_kg=stock,
        local_feed_remaining_kg=feed_left,unlevered_project_npv_U=unlevered_npv,
        closure=dict(remaining_book_U=remaining_book,assumed_salvage_U=salvage,
          zero_npv_terminal_recovery_fraction_of_remaining_book=(salvage-unlevered_npv*(1+c['economic_real_rate'])**years)/remaining_book if remaining_book else None,
          replacement_reserve_released_U=reserve,closure_escrow_released_U=escrow,
          retained_operating_cash_U=retained,remaining_loan_to_settle_U=closing_loan,
          decommission_cost_U=decommission,terminal_cash_after_obligations_U=terminal_cash,
          gap_U=max(0,-terminal_cash),rule='Liquidation/financial flows, not new GDP. Annual negative gaps require external funding; equity returns are not guaranteed.'),
        annual=annual)


def results():
    s,products=load_inputs();bound=anchors(products);reference=s['reference']
    scenarios=[];ref=None
    for plant in s['plants']:
        for name,case in s['cost_cases'].items():
            for years in s['operating_windows_years']:
                is_ref=plant['id']==reference['plant'] and name==reference['cost_case'] and years==reference['years']
                row=run_case(s,bound,plant,case,years,keep_annual=is_ref)
                row['cost_case']=name
                if is_ref: ref=row
                # Keep complete annual ledger once; other cases retain summaries.
                scenarios.append({k:v for k,v in row.items() if k not in ('annual','offers')})
    ref_plant=next(p for p in s['plants'] if p['id']==reference['plant'])
    ref_case=s['cost_cases'][reference['cost_case']]
    funding_priority=run_case(s,bound,ref_plant,ref_case,reference['years'],independent_first=True)
    closure_sensitivities=[]
    for fraction in (0,.3,.7,1):
        changed=copy.deepcopy(s);changed['capital']['closure_salvage_fraction_of_remaining_book']=fraction
        rr=run_case(changed,bound,ref_plant,ref_case,reference['years'])
        closure_sensitivities.append(dict(terminal_recovery_fraction=fraction,unlevered_npv_U=rr['unlevered_project_npv_U'],closure=rr['closure']))
    post={name:{f['id']:{b:offer(s,f,case,b,after_atmosphere=True) for b in s['buyers']} for f in s['families']} for name,case in s['cost_cases'].items()}
    names=('industrial_markets.py','industrial_markets_inputs.json','finance.py')
    return dict(schema='terluna.research.industrial-markets/1',
        producer=dict(files={str((HERE/n).relative_to(ROOT)):digest(HERE/n) for n in names},
                      constants=constants_used([__file__]),input_products=s['products']),
        evidence=s['evidence'],reading_rule=s['numeraire'],anchors=bound,reference=ref,
        scenarios=scenarios,independent_buyer_priority=funding_priority,
        terminal_recovery_sensitivities=closure_sensitivities,
        post_atmosphere_cost_thresholds=post,boundaries=s['boundaries'])


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,default=ROOT/'research/runs/world_economy/industrial_markets.json')
    a=parser.parse_args();a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(results(),indent=2,allow_nan=False)+'\n');print(a.output)


if __name__=='__main__':main()
