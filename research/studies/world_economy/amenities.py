"""Bound hybrid-amenity, six-route tourism and development-finance scenarios.

The established GWP/NWP reference is preserved. Tourism is a satellite
attribution, and project cash is a separate use of realised income. Activity
curves describe independent revenue-producing activity; loan-funded project
construction output must not be recycled as independent earnings.
"""
from __future__ import annotations
import argparse
import json
import math
from pathlib import Path
from shared.constants import JULIAN_YEAR_DAYS
from shared.provenance import constants_used
from . import inventory, tourism, development, mining

HERE = Path(__file__).resolve().parent
INPUT_PATH = HERE / 'amenities_inputs.json'
WORLDS = ('earth', 'moon', 'array')


def load_inputs(path=INPUT_PATH, root=inventory.ROOT):
    settings = json.loads(Path(path).read_text())
    if settings['schema'] != 'terluna.research.world-economy-amenities-inputs/1':
        raise ValueError('Unexpected amenities input schema')
    products = {}
    for name, spec in settings['products'].items():
        p = Path(root) / name
        if inventory.digest(p) != spec['sha256']:
            raise ValueError(f'Input changed; review and rebind explicitly: {name}')
        product = json.loads(p.read_text())
        if product['schema'] != spec['schema']:
            raise ValueError(f'Unexpected input schema: {name}')
        products[name] = product
    return settings, products


def activity_curve(knots, years):
    """Piecewise linear commissioned activity; no mature Moon income at year1."""
    if len(knots) < 2 or knots[0] != [0, 0] or knots[-1] != [years, 1]:
        raise ValueError('Activity knots must span [0,0] to [years,1]')
    for (x, a), (y, b) in zip(knots, knots[1:]):
        if not (0 <= x < y <= years and 0 <= a <= b <= 1):
            raise ValueError('Activity dates and fractions must be ordered')
    values = []
    for year in range(1, years+1):
        a, b = next((a,b) for a,b in zip(knots,knots[1:]) if a[0] < year <= b[0])
        values.append(a[1] + (b[1]-a[1]) * (year-a[0]) / (b[0]-a[0]))
    return values


def tourism_case(settings, baseline, case_id, *, access_scale=1, capacity_scale=1):
    config = settings['tourism']
    case = config['cases'][case_id]
    production = baseline['production']['close_parity']['worlds']
    worlds = {}
    for w in WORLDS:
        gdp = production[w]['gross_product']
        worlds[w] = dict(outbound_budget_U=gdp*case['outbound_budget_share_of_gwp'][w],
            local_amenity_spend_U=gdp*config['local_amenity_spend_share_of_gwp'][w],
            inbound_visitor_night_capacity=config['inbound_visitor_night_capacity'][w]*capacity_scale,
            production={c:dict(direct_gva_share=s, imported_input_shares=config['imported_input_shares'][w])
                        for c,s in config['direct_gva_shares'].items()})
    routes = []
    for origin in WORLDS:
        for destination, share in case['route_budget_shares'][origin].items():
            pair = next(k for k in config['route_access_cost_U'] if set(k.split('_')) == {origin,destination})
            routes.append(dict(origin=origin,destination=destination,budget_share=share,
                stay_nights=config['destination_stay_nights'][destination],
                destination_spend_per_night_U=config['destination_spend_per_night_U'][destination],
                access_cost_per_trip_U=config['route_access_cost_U'][pair]*access_scale,
                operator_shares=config['operator_shares'][pair],
                trip_capacity=config['route_trip_capacity_per_year']*capacity_scale))
    required = baseline['reference_financing']['patient']['trade_requirements']
    answer = tourism.evaluate(worlds,routes,
        baseline_leisure_gva={w:production[w]['sector_gva']['leisure_visitor'] for w in WORLDS},
        required_net_exports={w:required[w]['required_net_exports_U'] for w in WORLDS})
    answer['applied_scenario'] = dict(case_id=case_id,access_scale=access_scale,capacity_scale=capacity_scale)
    for w, row in answer['worlds'].items():
        inbound = sum(r['realized_visits'] for r in answer['routes'] if r['destination']==w)
        row['inbound_visits'] = inbound
        row['average_inbound_visitors_present'] = row['inbound_capacity']['realized_visitor_nights']/JULIAN_YEAR_DAYS
        row['inbound_average_share_of_resident_population'] = row['average_inbound_visitors_present']/production[w]['population']
        row['outbound_trips_per_resident_year'] = row['origin_budget']['realized_visits']/production[w]['population']
        row['direct_tourism_gva_share_of_gwp'] = row['direct_tourism_gva_U']/production[w]['gross_product']
        row['local_amenity_direct_gva_share_of_gwp'] = row['local_amenity_direct_gva_U']/production[w]['gross_product']
    return answer


def development_case(settings, baseline, costs, activity_id, cash_share, rate, transfer_share, *, idle_stock_stress=False):
    config = settings['development']
    build_years = config['build_years']
    years = max(knots[-1][0] for knots in config['activity_fraction_knots'][activity_id].values())
    production = baseline['production']['close_parity']['worlds']
    curves = {w:activity_curve(knots,years) for w,knots in config['activity_fraction_knots'][activity_id].items()}
    cash = {w:[fraction*production[w]['gross_product']*cash_share for fraction in curve]
            for w,curve in curves.items()}
    result = {}
    idle_costs = {}
    for w in ('array','moon'):
        total_cost = costs[w+'_territorial_capital_U']
        draws = [total_cost/build_years if year<=build_years else 0 for year in range(1,years+1)]
        grants = [x*config['fixed_sponsor_share_of_construction'] for x in draws]
        mature_account = baseline['mature_accounts']['close_parity']['accounts']['worlds'][w]
        mature_cfc = mature_account['ordinary_cfc']+mature_account['engineered_cfc']
        idle_costs[w] = [max(0,min(year/build_years,1)-a)*mature_cfc if idle_stock_stress else 0 for year,a in enumerate(curves[w],1)]
        # Both comparisons fund the full extra carrying cost once. The zero-
        # contribution baseline cannot silently spend part of early cash.
        draws = [x+idle for x,idle in zip(draws,idle_costs[w])]
        receipts = list(cash[w])
        if w == 'moon':
            transfers = [r['cash_released']*transfer_share for r in result['array']['early']['annual']]
            receipts = [x+y for x,y in zip(receipts,transfers)]
        result[w] = development.compare_with_deferred(draws,receipts,rate,
            annual_grants=grants,release_excess_cash=True,
            mature_repayment_years=config['repayment_years'])
    summary = dict(activity_case=activity_id,cash_share=cash_share,annual_rate=rate,
                   array_transfer_share=transfer_share,full_population_year=years,idle_stock_stress=idle_stock_stress,worlds={})
    total_transfer = math.fsum(transfers)
    for w,r in result.items():
        e,d = r['early'],r['deferred']
        gross_revenue = math.fsum(curves[w])*production[w]['gross_product']
        mature = baseline['mature_accounts']['close_parity']['accounts']['worlds'][w]
        renewal = math.fsum(curves[w])*(mature['ordinary_cfc']+mature['engineered_cfc'])
        net_product = gross_revenue-renewal
        summary['worlds'][w] = dict(
            cumulative_realised_activity_gwp_U=gross_revenue,
            proportional_activity_renewal_memo_U=renewal,
            cumulative_own_project_available_cash_U=math.fsum(cash[w]),
            cumulative_extra_idle_stock_care_U=math.fsum(idle_costs[w]),
            saving_share_of_activity_net_product=math.fsum(cash[w])/net_product,
            cumulative_received_transfer_U=total_transfer if w=='moon' else 0,
            cumulative_paid_transfer_U=total_transfer if w=='array' else 0,
            cumulative_released_cash_retained_U=e['totals']['cash_released']-(total_transfer if w=='array' else 0),
            debt_at_full_population_without_early_cash_U=d['closing_debt'],
            debt_at_full_population_with_early_cash_U=e['closing_debt'],
            debt_reduction_U=r['closing_debt_reduction'],
            debt_reduction_share=r['closing_debt_reduction']/d['closing_debt'],
            annual_service_without_early_cash_U=d['final_annual_service']['payment'],
            annual_service_with_early_cash_U=e['final_annual_service']['payment'],
            first_year_remaining_debt_free=e['first_year_remaining_debt_free_through_horizon'],
            cumulative_interest_avoided_U=r['interest_avoided'],
            cumulative_cash_applied_U=e['totals']['early_cash_applied_to_project'],
            cumulative_construction_borrowing_U=e['totals']['new_construction_borrowing'],
            cumulative_sponsor_funding_U=e['totals']['construction_grant'],
            annual_cash_closure_max_U=e['max_abs_cash_closure_residual'],
            debt_avoidance_identity_residual_U=r['debt_avoidance_identity_residual'],
            checkpoints=[dict(year=year,activity_fraction=curves[w][year-1],
                             own_project_available_cash_U=cash[w][year-1],
                             **{k:e['annual'][year-1][k] for k in ('closing_debt','new_construction_borrowing','accrued_interest','principal_prepayment','cash_released')})
                         for year in sorted({1,50,100,200,300,400,500,years})])
    # Keep existing creditor residence and ordinary replacement unchanged at maturity.
    # Transfers above concern construction years only; no perpetual array tribute is assumed.
    retired = baseline['reference_financing']['legacy_retired']['accounts']['worlds']
    shares = baseline['reference_financing']['patient']['applied_parameters']
    creditor = shares['earth_creditor_share']
    resource = {w:retired[w]['net_world_product'] for w in WORLDS}
    for w in ('moon','array'):
        payment = summary['worlds'][w]['annual_service_with_early_cash_U']*creditor
        resource[w] -= payment
        resource['earth'] += payment
    summary['mature_cash_resources_U'] = resource
    summary['mature_cash_ratios_to_earth'] = {w:resource[w]/resource['earth'] for w in WORLDS}
    summary['cross_world_transfer_closure_U'] = summary['worlds']['array']['cumulative_paid_transfer_U']-summary['worlds']['moon']['cumulative_received_transfer_U']
    return summary


def mining_cases(settings, products, costs, baseline):
    """Finite raw-feed envelope; exports follow an ideal P-only separation.

    The generic export screen uses feed-equivalent kg throughout; revenue per
    feed kg is adjusted for the retained P mass. Actual exported product and
    retained nutrient masses are reported separately and conserve feed mass.
    """
    config = settings['mining']
    grades = products['resources/results/sea_phosphorus.json']['land_loss']['soil_p_grade']
    shield = next(a for a in costs['installed_assets'] if a['element']=='array_bulk_habitat_shield')
    years = settings['development']['build_years']
    cost = costs['moon_territorial_capital_U']
    grants = [cost/years*settings['development']['fixed_sponsor_share_of_construction']]*years
    rows = []
    for grade_name in config['grade_cases_for_sensitivity']:
        grade = grades[grade_name]
        for retention in config['phosphorus_retention_fractions']:
            product_fraction = 1-grade*retention
            nutrient_cap = config['phosphorus_export_allowance_kg']/(grade*(1-retention)) if retention<1 else config['feed_pool_kg']
            for independent in config['independent_buyer_fractions']:
                x = mining.finite_export_screen(
                    recoverable_inventory_kg=config['feed_pool_kg'],reserved_future_kg=config['reserved_future_kg'],
                    protected_kg=config['protected_kg'],uncertainty_buffer_kg=config['uncertainty_buffer_kg'],
                    buyer_demand_cap_kg=min(config['customer_finished_product_demand_cap_kg']/product_fraction,nutrient_cap),
                    years=config['years'],annual_export_capacity_kg=config['annual_feed_capacity_kg'],
                    price_U_per_kg=shield['unit_cost_U']*product_fraction,
                    independent_buyer_fraction=independent,**config['unit_costs_per_feed_kg'])
                series = [0.0]*years
                start = config['first_operating_year']-1
                if start < 0 or start+config['years']>years:
                    raise ValueError('Mining activity must fall within the financing horizon')
                series[start:start+config['years']] = x['annual_uniform_independent_project_cash_U']
                if min(series)<0:
                    raise ValueError('Resolve mining operating losses before assigning project cash')
                loan = development.compare_with_deferred([cost/years]*years,series,.001,annual_grants=grants)
                # Annual series are uniform and reproduced from totals, so avoid
                # storing36 copies of300 identical values in the committed product.
                for key in ('annual_uniform_export_kg','annual_uniform_independent_project_cash_U','annual_uniform_excluded_programme_funded_cash_U'):
                    x.pop(key)
                feed = x['quantity_kg']
                rows.append(dict(grade_case=grade_name,raw_feed_p_grade=grade,
                    phosphorus_retention_fraction=retention,independent_buyer_fraction=independent,
                    feed_equivalent_account=x,exported_finished_product_kg=feed*product_fraction,
                    retained_phosphorus_kg=feed*grade*retention,
                    exported_phosphorus_kg=feed*grade*(1-retention),
                    feed_mass_closure_kg=feed-feed*product_fraction-feed*grade*retention,
                    phosphorus_budget_feed_cap_kg=nutrient_cap,
                    annual_independent_cash_U=x['independent_project_cash_U']/config['years'],
                    standalone_debt_avoided_at_year500_U=loan['closing_debt_reduction'],
                    standalone_reduction_share_of_patient_moon_debt=loan['closing_debt_reduction']/loan['deferred']['closing_debt']))
    return dict(evidence=config['evidence'],geological_export_authorization_established=False,
        array_shield_finished_stock_kg=shield['quantity'],bound_finished_product_price_U_per_kg=shield['unit_cost_U'],
        rows=rows,accounting_rule=config['funding_rule'],
        reserve_rule='Bulk stock, strategic element, accessible deposit, environmental and protected-site reservations all require separate review; these synthetic mass buckets do not establish surplus lunar ore.',
        coproduct_rule='Only retained P is explicitly partitioned. K, other nutrients, water, trace metals and usable-tailings recovery impose additional unquantified restrictions; no unmodelled coproduct receives a revenue credit.')


def results():
    settings, products = load_inputs()
    baseline = products['research/studies/world_economy/results/comparison.json']
    cost_product = products['research/studies/world_economy/results/construction.json']
    dist = products['research/studies/distributions/results/distributions.json']
    ref = baseline['reference']
    costs = next(c for c in cost_product['cases'] if c['cost_case']==ref['cost_case'] and c['water_case']==ref['water_case'])
    area = dist['basis']['moon']['dry_land_km2']*dist['assumptions']['land_rules']['wild_30_hybrid_20']['hybrid']
    population = baseline['production']['close_parity']['worlds']['moon']['population']
    travel = {k:tourism_case(settings,baseline,k) for k in settings['tourism']['cases']}
    travel['connected_half_access_cost'] = tourism_case(settings,baseline,'connected',access_scale=.5)
    travel['connected_double_access_cost'] = tourism_case(settings,baseline,'connected',access_scale=2)
    travel['lunar_destination_double_capacity'] = tourism_case(settings,baseline,'lunar_destination',capacity_scale=2)
    attribution = settings['tourism']['hybrid_attribution']
    moon = travel['connected']['worlds']['moon']
    hybrid = dict(area_km2=area,m2_per_lunar_resident=area*1e6/population,
        assumed_resident_visits_year=population*attribution['resident_visits_per_year'],
        assumed_inbound_visitors_using_hybrid_year=moon['inbound_visits']*attribution['inbound_visits_with_hybrid_use_share'],
        attributed_local_direct_gva_U=moon['local_amenity_direct_gva_U']*attribution['moon_local_amenity_share'],
        attributed_visitor_direct_gva_U=moon['visitor_facing_direct_gva_U']*attribution['moon_inbound_service_share'],
        evidence='Accepted area; assumed use frequency and expenditure attribution. No new GWP, housing premium, food harvest or welfare value added.')
    hybrid['attributed_total_direct_gva_U'] = hybrid['attributed_local_direct_gva_U']+hybrid['attributed_visitor_direct_gva_U']
    hybrid['mean_resident_visits_per_km2_day'] = hybrid['assumed_resident_visits_year']/area/JULIAN_YEAR_DAYS
    park_costs = {}
    for name, case in settings['park_cost_sensitivity']['cases'].items():
        capital = area*1e6*case['initial_produced_capital_U_per_m2']
        operations = area*1e6*case['annual_operating_U_per_m2']
        renewal = capital/case['replacement_life_years']
        park_costs[name] = dict(produced_capital_U=capital,annual_operating_U=operations,
            annual_renewal_U=renewal,annual_cash_requirement_U=operations+renewal,
            share_of_reference_moon_gwp=(operations+renewal)/baseline['production']['close_parity']['worlds']['moon']['gross_product'])
    config = settings['development']
    development_rows = [development_case(settings,baseline,costs,a,s,r,t)
        for a in config['activity_fraction_knots']
        for s in config['project_available_cash_shares_of_realised_gwp']
        for r in config['annual_real_rates']
        for t in config['array_excess_transfer_to_moon_shares']]
    idle_stress = development_case(settings,baseline,costs,'phased',.05,.001,0,idle_stock_stress=True)
    reference = next(row for row in development_rows if all(row[k]==v for k,v in config['reference'].items()))
    filenames = ('amenities.py','amenities_inputs.json','tourism.py','development.py','finance.py','inventory.py','mining.py')
    return dict(schema='terluna.research.world-economy-amenities/1',
        producer=dict(study='world_economy',files={str((HERE/n).relative_to(inventory.ROOT)):inventory.digest(HERE/n) for n in filenames},
                      constants=constants_used([__file__]),input_products=settings['products']),
        evidence=settings['evidence'],baseline_commit=settings['baseline_commit'],
        reading_rule='Annual prices and days use common real U and24-hour Earth-equivalent days. Local amenity and tourism results attribute existing output; no extra GWP/NWP. Early cash is after operating/renewal/consumption and independent of the construction loan. No market feasibility is established.',
        hybrid=hybrid,park_cost_sensitivities=park_costs,tourism=travel,mining=mining_cases(settings,products,costs,baseline),
        development_reference=reference,development_scenarios=development_rows,idle_stock_funding_stress=idle_stress,
        limitations=['Hybrid water, maintenance machinery, utilities and seasonal/peak capacity remain unsized.',
            'The301K redesign changes comfort and ecological management; older amenity descriptions are not new regional comfort simulations.',
            'Tourism assumes visitor status, service-only spending and a first-tier import account; local amenities include destination services on domestic trips, with domestic passenger fares excluded.',
            'Capacity and commissioning schedules are explicit assumptions; early array real activity does not establish full lunar outdoor access.',
            'Construction remains uniformly funded over500years; no resource-feasible milestone network is solved.',
            'Early activity and saving shares presume independent demand and fully funded renewal; source-site service prices already include supplier renewal.',
            'Array surplus transfers require ownership or fiscal arrangements. Debt/equity finance alone does not eliminate claims on future income.',
            'The baseline ordinary capital allowance must accommodate park and export businesses; extra capital beyond that allowance needs additional funding.'])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,default=inventory.ROOT/'research/runs/world_economy/amenities.json')
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(results(),indent=2)+'\n')
    print(args.output)


if __name__ == '__main__':
    main()
