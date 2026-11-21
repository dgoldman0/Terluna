"""The accounts of a financed case: costs by category and receipts by sector over reporting periods, the material a
case draws to each horizon, and whether its industry was sized without shortfall.

The staged studies' functions, taken whole from research/studies/enclosed_returns (bottlenecks.annual_costs,
revised.stratification, revised.cumulative_demand, revised.valid) as they stood at e49fca09, each written one
statement to a line with its syntax unchanged; cumulative_demand leaves out the staged studies' seed-bound option, which
the study never used.
"""
from __future__ import annotations
import math

from . import finance


# From research/studies/enclosed_returns/bottlenecks.py at e49fca09, annual_costs.
def annual_costs(request, physical, row, cfg):
    cash = physical['cash']
    costs = {key: cash[key] for key in ('operating_external_replacement_USD2041', 'operating_external_logistics_materials_USD2041', 'operating_external_services_USD2041', 'operating_local_payroll_USD2041', 'operating_local_replacement_USD2041', 'operating_resource_provision_USD2041', 'capital_external_USD2041', 'capital_local_USD2041', 'capital_resource_provision_USD2041')}
    costs.update(commercial_operating_USD2041=request['commercial_operating_USD2041'], hotel_payroll_USD2041=request['visitor_beds'] * cfg['hotel_staff_per_bed'] * cfg['hotel_staff_USD2041_year'], residual_service_renewal_USD2041=row['service_residual_renewal'], residual_service_capital_USD2041=row['service_residual_capital'], capital_allowance_USD2041=row['capital_cost_allowance'])
    if not math.isclose(sum(costs.values()), row['capital'] + row['operating_cost'], rel_tol=1e-12, abs_tol=0.01):
        raise ValueError('Cost stratification does not reconcile')
    return costs

# From research/studies/enclosed_returns/revised.py at e49fca09, stratification.
def stratification(case, physical, path, sim, cfg, periods):
    rows = []
    for request, pr, r, f in zip(case['annual'], physical['annual'], path['annual'], sim['annual']):
        costs = annual_costs(request, pr, r, cfg)
        rows.append({'year': r['year'], 'revenue': r['revenue_by_sector'], 'costs': costs, 'payout': f['owner_payout'], 'contribution': f['contribution']})
    result = []
    start = path['annual'][0]['year']
    for first, last in periods:
        selected = [r for r in rows if first <= r['year'] <= last]
        values = {}
        for rate in (0.0, 0.03):
            revenues = {k: math.fsum((r['revenue'][k] / (1 + rate) ** r['year'] for r in selected)) for k in rows[0]['revenue']}
            costs = {k: math.fsum((r['costs'][k] / (1 + rate) ** r['year'] for r in selected)) for k in rows[0]['costs']}
            costs['seed_deployment_and_allowance_USD2041'] = math.fsum((path['seed_hardware_and_carry'] / start / (1 + rate) ** y for y in range(start) if first <= y <= last))
            values[str(rate)] = {'revenue_by_sector': revenues, 'costs_by_category': costs, 'net_real_cash': sum(revenues.values()) - sum(costs.values()), 'distributions': math.fsum((r['payout'] / (1 + rate) ** r['year'] for r in selected)), 'additional_equity': math.fsum((r['contribution'] / (1 + rate) ** r['year'] for r in selected))}
        result.append({'first_year': first, 'last_year': last, 'values': values})
    return {'periods': result, 'financing_bridges': {str(h): finance.carry_bridge(path, sim, cfg, h, 0.03) for h in cfg['horizons']}, 'peak_retained_cash': max(({'year': r['year'], 'cash': r['cash']} for r in sim['annual']), key=lambda r: r['cash'])}

# From research/studies/enclosed_returns/revised.py at e49fca09, cumulative_demand.
def cumulative_demand(physical, audit, h):
    if audit is not None:
        rows = [r for r in audit['annual'] if r['year'] <= h]
        if [r['year'] for r in rows] != list(range(h + 1)):
            raise ValueError('Missing dated physical audit')
        total = lambda key: math.fsum((r[key] for r in rows))
        return {'feed_kg': total('source_feed_kg'), 'lunar_raw_kg': total('lunar_raw_kg'), 'water_kg_bounds': [total('source_water_equivalent_kg')] * 2, 'nitrogen_kg': total('source_nitrogen_kg'), 'virgin_metal_kg': total('source_virgin_metal_kg'), 'seed_rule': 'Exact inherited seed total counted once; seed years use its eight-year average.'}
    rows = [r for r in physical['annual'] if r['year'] <= h]
    seed = physical['seed']
    water = math.fsum((r['physical']['source_water_kg_year'] for r in rows))
    water_bounds = [water, water + seed['source_feed_kg'] * 0.1 * 0.7]
    seed_rule = 'Frozen compact baseline omits seed coproducts. Water lower bound excludes seed; upper bound adds seed feed times inherited 10 percent grade and 70 percent recovery.'
    return {'feed_kg': seed['source_feed_kg'] + math.fsum((r['physical']['source_feed_kg_year'] for r in rows)), 'lunar_raw_kg': seed['lunar_raw_kg'] + math.fsum((r['physical']['lunar_raw_kg_year'] for r in rows)), 'water_kg_bounds': water_bounds, 'nitrogen_kg': None, 'virgin_metal_kg': None, 'seed_rule': seed_rule}

# From research/studies/enclosed_returns/revised.py at e49fca09, valid.
def valid(physical):
    return physical['convergence']['converged'] and physical['first_capacity_shortfall'] is None
