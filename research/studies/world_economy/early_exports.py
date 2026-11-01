"""Early enclosed-community exports: destination tourism and candidate factories.

Milestone capacities and customer orders are assumptions. Cash after facility
costs is distinct from GDP, net exports and a complete resident budget. New
capital is financed explicitly; it must be reconciled with the mature capital
allowance before combining this screen with the whole-world financing model.
"""
from __future__ import annotations
import argparse
import copy
import json
import math
from pathlib import Path
from shared.constants import JULIAN_YEAR_DAYS
from shared.provenance import constants_used
from . import finance, inventory

HERE = Path(__file__).resolve().parent
INPUT_PATH = HERE/'early_exports_inputs.json'


def load_inputs(path=INPUT_PATH, root=inventory.ROOT):
    settings=json.loads(Path(path).read_text())
    if settings['schema']!='terluna.research.early-exports-inputs/1':
        raise ValueError('Unexpected early export input schema')
    products={}
    for name,spec in settings['products'].items():
        p=Path(root)/name
        if inventory.digest(p)!=spec['sha256']:
            raise ValueError(f'Input changed; review and rebind: {name}')
        data=json.loads(p.read_text())
        if data['schema']!=spec['schema']:
            raise ValueError(f'Unexpected input schema: {name}')
        products[name]=data
    return settings,products


def facility_account(sales, capital, life, costs, *, operating_year=1):
    """Cash includes own loans, replacement and retained growth once each."""
    if sales<0 or capital<0 or life<=0 or operating_year<1:
        raise ValueError('Invalid facility quantity or age')
    intermediate=sales*costs['variable_intermediate_share']
    compensation=sales*costs['compensation_and_distributed_living_share']
    fixed=capital*costs['fixed_operating_share_of_capital']
    fixed_intermediate=fixed*costs['fixed_operating_intermediate_fraction']
    cfc=capital/life
    growth=sales*costs['retained_growth_share_of_sales']
    loan=capital*costs['financed_fraction']
    rate=costs['facility_real_rate']; term=costs['facility_repayment_years']
    payment=finance.annuity_payment(loan,rate,term) if operating_year<=term else 0
    opening=finance.loan_balance_after_payments(loan,rate,term,operating_year-1)
    closing=finance.loan_balance_after_payments(loan,rate,term,operating_year)
    interest=opening*rate
    principal=payment-interest
    gva=sales-intermediate-fixed_intermediate
    net_value_added=gva-cfc
    cash=sales-intermediate-compensation-fixed-cfc-growth-payment
    # Bridge from production to cash: salaries, unallocated fixed compensation,
    # new investment and all debt service use cash without being NWP deductions.
    bridge=net_value_added-compensation-(fixed-fixed_intermediate)-growth-payment
    return dict(sales_U=sales,initial_capital_U=capital,initial_borrowing_U=loan,
        initial_sponsor_or_equity_U=capital-loan,intermediate_U=intermediate,
        compensation_and_distributed_living_U=compensation,fixed_operations_U=fixed,
        direct_gva_U=gva,cfc_and_replacement_U=cfc,net_value_added_before_depletion_U=net_value_added,
        retained_growth_U=growth,facility_debt_service_U=payment,
        interest_U=interest,principal_U=principal,closing_facility_debt_U=closing,
        free_cash_U=cash,cash_bridge_residual_U=cash-bridge,
        cash_closure_residual_U=sales-(intermediate+compensation+fixed+cfc+growth+payment+cash))


def evaluate(settings,stage,case,*,operating_year=1):
    if settings['stay_nights']<=0 or stage['residents']<=0:
        raise ValueError('Positive stay and residents required')
    costs=copy.deepcopy(settings['costs'])
    costs['variable_intermediate_share']+=case['variable_operating_share_addition']
    if costs['variable_intermediate_share']+costs['compensation_and_distributed_living_share']>1:
        raise ValueError('Variable operating shares exceed receipts')
    pop=stage['residents']; nights=settings['stay_nights']
    beds=pop*settings['bed_places_per_resident']*case['bed_scale']
    occupancy=stage['occupancy_ceiling'] if case['occupancy_override'] is None else case['occupancy_override']
    if not 0<=occupancy<=1:
        raise ValueError('Occupancy must be a fraction')
    capacity_nights=min(beds*occupancy,beds*settings['spare_life_support_guest_places_per_bed'])*JULIAN_YEAR_DAYS
    price=stage['nightly_destination_spend_U']*case['price_scale']
    origin_population={'earth':settings['earth_population'],'array':stage['array_population']}
    origin_gwp={w:p*settings[w+'_gwp_per_resident_U'] for w,p in origin_population.items()}
    origin={}
    for w,p in origin_population.items():
        budget=origin_gwp[w]*stage['origin_budget_share_of_gwp'][w]
        ticket=stage['access_fare_U'][w]+nights*price
        requested=p*min(1,stage['annual_origin_participation'][w]*case['participation_scale'])
        allowed=min(requested,budget/ticket,stage['route_trip_capacity'][w]*case['route_scale'])
        origin[w]=dict(requested_trips=requested,budget_U=budget,trip_total_U=ticket,trips_before_bed_cap=allowed)
    allowed=sum(r['trips_before_bed_cap'] for r in origin.values())
    scale=min(1,capacity_nights/(allowed*nights)) if allowed else 1
    for w,r in origin.items():
        r['realized_trips']=r['trips_before_bed_cap']*scale
        r['destination_receipts_U']=r['realized_trips']*nights*price
        r['access_fares_excluded_U']=r['realized_trips']*stage['access_fare_U'][w]
        r['unspent_origin_budget_U']=r['budget_U']-r['realized_trips']*r['trip_total_U']
    trips=sum(r['realized_trips'] for r in origin.values())
    tourism_sales=sum(r['destination_receipts_U'] for r in origin.values())
    tourist_capital=beds*stage['visitor_capital_U_per_bed']*case['capital_price_scale']
    tourism=facility_account(tourism_sales,tourist_capital,costs['tourism_capital_life_years'],costs,operating_year=operating_year)
    utilization=settings['factory_utilization'] if case['factory_utilization_override'] is None else case['factory_utilization_override']
    unit_price=settings['factory_finished_price_U_per_kg']*case['factory_price_scale']
    potential_kg=pop*settings['factory_nameplate_finished_kg_per_resident_year']*settings['factory_saleable_yield']*utilization
    orders_U=sum(origin_gwp[w]*share for w,share in settings['factory_customer_budget_shares_of_gwp'].items())
    made_kg=min(potential_kg,orders_U/unit_price) if case['manufacturing_enabled'] else 0
    factory_capital=pop*settings['factory_capital_U_per_resident']*case['capital_price_scale'] if case['manufacturing_enabled'] else 0
    manufacture=facility_account(made_kg*unit_price,factory_capital,costs['manufacturing_capital_life_years'],costs,operating_year=operating_year)
    shared=pop*costs['unfunded_shared_community_allowance_U_per_resident_year']
    residual=tourism['free_cash_U']+manufacture['free_cash_U']-shared
    assigned=max(0,residual)*costs['development_share_of_positive_residual']
    gap=max(0,-residual)
    retained=max(0,residual)-assigned
    return dict(stage=stage['id'],residents=pop,array_population=stage['array_population'],operating_year=operating_year,
        visitor_bed_places=beds,capacity_visitor_nights=capacity_nights,
        realized_trips=trips,realized_visitor_nights=trips*nights,
        actual_occupancy=trips*nights/beds/JULIAN_YEAR_DAYS if beds else 0,
        average_visitors_present=trips*nights/JULIAN_YEAR_DAYS,
        origin=origin,tourism=tourism,manufacturing=manufacture,
        factory_saleable_capacity_kg=potential_kg,factory_sold_kg=made_kg,
        factory_customer_orders_budget_U=orders_U,
        extra_shared_community_allowance_U=shared,
        total_destination_and_manufacturing_sales_U=tourism_sales+made_kg*unit_price,
        combined_free_cash_after_allowance_U=residual,assigned_development_cash_U=assigned,
        unfunded_cash_gap_U=gap,retained_cash_after_assignment_U=retained,
        allocation_closure_U=residual+gap-assigned-retained,
        first_year_total_facility_capital_U=tourist_capital+factory_capital,
        first_year_total_sponsor_or_equity_U=(tourist_capital+factory_capital)*(1-costs['financed_fraction']),
        total_closing_facility_debt_U=tourism['closing_facility_debt_U']+manufacture['closing_facility_debt_U'],
        depletion_rule='No depletion value is priced in this candidate-factory screen; NVA is before depletion. Feed extraction/import cost belongs in intermediate inputs, and strategic reserves remain constraints.',
        evidence='Assumed viable product, independent customers, capacity, price and costs; no gravitational quality premium or whole-community self-sufficiency established.')


def results():
    settings,products=load_inputs()
    baseline=products['research/studies/world_economy/results/comparison.json']
    old_debt=baseline['reference_financing']['patient']['finance_history']['moon']['debt_at_full_population']
    scenarios={name:[evaluate(settings,stage,case) for stage in settings['stages']] for name,case in settings['cases'].items()}
    financial_scale=[]
    for stage in settings['stages']:
        row=evaluate(settings,stage,settings['cases']['reference'])
        cash=row['assigned_development_cash_U']
        # Only 100 operating years credited, beginning in year 1; a timing scale,
        # not an integrated construction schedule or net-debt forecast.
        future_credit=math.fsum(cash*1.001**(500-year) for year in range(1,101))
        financial_scale.append(dict(stage=stage['id'],timing_rule='Constant milestone assigned cash in operating years 1–100, valued at year 500 at 0.1% real; new facility debt remains separately reported, so this is not integrated net debt avoided.',cash_for_100_years_U=cash*100,
            debt_equivalent_at_year500_U=future_credit,share_of_old_patient_debt=future_credit/old_debt,
            outstanding_new_facility_debt_after_year1_U=row['total_closing_facility_debt_U'],
            after_facility_loan_retirement=evaluate(settings,stage,settings['cases']['reference'],operating_year=101)))
    names=('early_exports.py','early_exports_inputs.json','finance.py','inventory.py')
    return dict(schema='terluna.research.early-exports/1',producer=dict(study='world_economy',
        files={str((HERE/n).relative_to(inventory.ROOT)):inventory.digest(HERE/n) for n in names},
        constants=constants_used([__file__]),input_products=settings['products']),
        evidence=settings['evidence'],reading_rule='B=1e9U,T=1e12U. One common real U, not future dollars. Annual milestone outputs; Earth-equivalent 24 h visitor nights. Positive assigned cash follows facility costs, renewal, growth, own debt and the stated extra shared-services allowance. Negative gaps remain explicit. Capital and cash are not additions to prior macroeconomic totals.',
        scenarios=scenarios,financing_scale=financial_scale,boundaries=settings['boundaries'])


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=inventory.ROOT/'research/runs/world_economy/early_exports.json')
    args=parser.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(results(),indent=2)+'\n');print(args.output)


if __name__=='__main__':main()
