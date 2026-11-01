"""Finite structural-module campaigns, with explicit initial capital and orders.

Scenario inputs describe qualified plant and binding orders, not discovered
markets. Operations stop after the horizon or full order delivery. A customer
advance is financing with a delivery obligation; it does not alter earned sales,
unlevered project NPV, terminal full-capital cash or programme independence.
"""
from __future__ import annotations

import math


def _value(name, x, positive=False):
    if isinstance(x,bool) or not isinstance(x,(int,float)) or not math.isfinite(x) or x<0 or (positive and x==0):
        raise ValueError(f'Invalid {name}')
    return float(x)


def campaign(*, raw_nameplate_kg_year, years, eligible_orders_kg,
             buyer_budget_U_year, base_production_cost_U_kg=5,
             lunar_cost_factor=.9, producer_price_U_kg=4.65,
             freight_U_kg=.25, utilization=.7, accepted_yield=.9,
             power_W, outbound_capacity_kg_year, seconds_per_year,
             process_energy_J_raw_kg=100e6, local_feed_fraction=.9,
             available_local_feed_kg=100e12, economic_rate=.01,
             capital_life_years=40, terminal_recovery_fraction=.3,
             decommission_fraction=.02, customer_advance_fraction=0,
             independently_funded_fraction=0,
             ramp=(.25,.5,.75,1)):
    """Return annual production and capital-inclusive cash and NPV.

    Capital uses25% of normal unit cost at its40-year economic annuity; fixed
    operations use5%, variable inputs55% and compensation15%. Those disjoint
    shares reproduce the earlier structural-product cost convention. Initial
    equity/debt are not assumed: peak negative cash is a funding requirement.
    Source electricity and route capacities must be supplied; process energy is
    already priced in production cost and is not charged a second time.
    Customer advances apply only to the modeled deliverable contract, require
    separate buyer funds at t0, and reduce cash received against later invoices.
    """
    if isinstance(years,bool) or not isinstance(years,int) or years<1:
        raise ValueError('Positive integer campaign duration required')
    if isinstance(capital_life_years,bool) or not isinstance(capital_life_years,int) or capital_life_years<1:
        raise ValueError('Positive integer capital life required')
    numeric=locals().copy()
    for name,x in numeric.items():
        if name not in ('ramp','years','capital_life_years'):
            _value(name,x)
    if not ramp or any(not 0<=x<=1 for x in ramp):raise ValueError('Ramp fractions required')
    for name in ('utilization','accepted_yield','local_feed_fraction','terminal_recovery_fraction','customer_advance_fraction','independently_funded_fraction'):
        if not 0<=numeric[name]<=1:raise ValueError(f'{name} must be a fraction')
    if min(raw_nameplate_kg_year,eligible_orders_kg,base_production_cost_U_kg,lunar_cost_factor,producer_price_U_kg,accepted_yield,seconds_per_year,process_energy_J_raw_kg)<=0:
        raise ValueError('Positive scale, order, cost, yield, time and energy required')
    normal=raw_nameplate_kg_year*utilization*accepted_yield
    cost=base_production_cost_U_kg*lunar_cost_factor
    rate=economic_rate;life=capital_life_years
    crf=rate/(-math.expm1(-life*math.log1p(rate))) if rate else 1/life
    capital=normal*cost*.25/crf;fixed=normal*cost*.05
    price=producer_price_U_kg+freight_U_kg
    power_cap=power_W*seconds_per_year/process_energy_J_raw_kg*accepted_yield
    remaining=eligible_orders_kg;feed_remaining=available_local_feed_kg
    rows=[]
    for year in range(1,years+1):
        feed_cap=feed_remaining/local_feed_fraction*accepted_yield if local_feed_fraction else math.inf
        made=min(normal*ramp[min(year-1,len(ramp)-1)],power_cap,outbound_capacity_kg_year,
                 buyer_budget_U_year/price,remaining,feed_cap)
        feed=made/accepted_yield
        local_feed=feed*local_feed_fraction
        feed_remaining=max(0,feed_remaining-local_feed);remaining=max(0,remaining-made)
        revenue=made*producer_price_U_kg
        intermediate=made*cost*.55;compensation=made*cost*.15
        gva=revenue-intermediate-fixed*.5
        rows.append(dict(year=year,qualified_delivery_kg=made,raw_processed_kg=feed,
            local_feed_kg=local_feed,energy_J=feed*process_energy_J_raw_kg,
            delivered_customer_spending_U=made*price,earned_producer_revenue_U=revenue,
            intermediate_U=intermediate,compensation_U=compensation,fixed_operations_U=fixed,
            gva_U=gva,cfc_U=capital/life,nva_before_depletion_U=gva-capital/life,
            operating_cash_before_capital_U=revenue-intermediate-compensation-fixed))
        if remaining<=eligible_orders_kg*1e-14:break
    end=rows[-1]['year']
    age=end-((end-1)//life)*life
    book=capital*max(0,1-age/life)
    terminal_recovery=book*terminal_recovery_fraction
    decommission=capital*decommission_fraction
    contract_revenue=math.fsum(r['earned_producer_revenue_U'] for r in rows)
    advance=customer_advance_fraction*contract_revenue
    liability=advance
    unfunded_cash=-capital;financed_cash=advance-capital
    min_unfunded=unfunded_cash;min_financed=financed_cash
    npv=-capital
    for row in rows:
        year=row['year']
        replacement=capital if year%life==0 and year<end else 0
        terminal=terminal_recovery-decommission if year==end else 0
        capital_cash=row['operating_cash_before_capital_U']-replacement+terminal
        repayment=customer_advance_fraction*row['earned_producer_revenue_U']
        liability=max(0,liability-repayment)
        cash_receipts=row['earned_producer_revenue_U']-repayment
        unfunded_cash+=capital_cash;financed_cash+=capital_cash-repayment
        min_unfunded=min(min_unfunded,unfunded_cash);min_financed=min(min_financed,financed_cash)
        npv+=capital_cash/(1+rate)**year
        row.update(actual_replacement_U=replacement,terminal_recovery_net_U=terminal,
            capital_inclusive_cash_change_U=capital_cash,
            cash_received_after_advance_offset_U=cash_receipts,
            customer_advance_earned_release_U=repayment,
            closing_customer_contract_liability_U=liability,
            cumulative_cash_after_all_capital_U=unfunded_cash,
            cumulative_cash_with_customer_finance_U=financed_cash)
    recovery=next((r['year'] for i,r in enumerate(rows)
        if all(s['cumulative_cash_after_all_capital_U']>=0 for s in rows[i:])),None)
    total_mass=math.fsum(r['qualified_delivery_kg'] for r in rows)
    return dict(scheduled_years=years,effective_operating_years=end,
        normal_qualified_output_kg_year=normal,initial_capital_U=capital,
        initial_customer_advance_U=advance,
        additional_upfront_funding_after_advance_U=max(0,capital-advance),
        peak_funding_without_customer_advance_U=max(0,-min_unfunded),
        peak_funding_with_customer_advance_U=max(0,-min_financed),
        delivered_qualified_mass_kg=total_mass,remaining_orders_kg=remaining,
        remaining_local_feed_kg=feed_remaining,
        total_earned_producer_revenue_U=contract_revenue,
        total_customer_delivered_spending_U=total_mass*price,
        total_gva_U=math.fsum(r['gva_U'] for r in rows),
        total_nva_before_depletion_U=math.fsum(r['nva_before_depletion_U'] for r in rows),
        terminal_remaining_book_U=book,terminal_recovery_U=terminal_recovery,
        decommission_U=decommission,terminal_customer_liability_U=liability,
        final_cash_after_all_capital_U=unfunded_cash,
        final_cash_with_customer_finance_U=financed_cash,
        customer_finance_cash_closure_U=financed_cash-unfunded_cash,
        unlevered_project_npv_U=npv,sustained_capital_recovery_year=recovery,
        independently_funded_fraction=independently_funded_fraction,
        illustrative_independent_fraction_of_positive_final_surplus_U=max(0,unfunded_cash)*independently_funded_fraction,
        annual=rows,
        interpretation='Conditional qualified orders and campaign. Peak funding must be supplied by named equity/debt/advances; no financing authorization. Array programme purchases and their advances do not add independent joint-programme funding. Final full-capital surplus includes monetized terminal recovery; retained assets cannot supply that cash simultaneously.')
