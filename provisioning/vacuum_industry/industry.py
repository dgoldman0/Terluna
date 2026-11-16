"""The industry for one schedule of residents and orders: the model solved, its years compacted, its accounts drawn.

build() takes the engineering coefficients, the recovery measures and a book of yearly orders, and returns two
records. The first is the industry year by year, with its physical flows, its costs, its capacity in use and the
founding deployment's cost. The second is the account of the loops: the water, the recovered material and the
refurbished machines of each year, with the material drawn from the near-Earth sources and the Moon.
"""
from __future__ import annotations
from copy import deepcopy
import math
from . import model, recovery

IDENTITY = ('name', 'portfolio', 'market', 'market_growth', 'market_growth_name')
STOCKS = ('asteroid_tugs', 'lunar_landers', 'Earth_hull_tanks', 'source_mine_trains', 'source_metal_trains', 'bulk_trains',
          'metal_trains', 'source_solar', 'lunar_beam_constellation', 'construction_and_service_robots')
DETAIL = ('request', 'required_asset_mass_kg', 'available_asset_mass_kg', 'capacity_use')


def compact(book, result):
    """The solved case without its per-asset detail: each year keeps its flows, costs and a summary of capacity."""
    out = {k: book[k] for k in IDENTITY if k in book}
    out['physical_scenario'] = 'base'
    out['seed'] = deepcopy(result['seed'])
    out['seed'].pop('stock_mass_kg')
    out['seed'].pop('physical_annual_average')
    out['convergence'] = result['convergence']
    out['first_capacity_shortfall'] = result['first_capacity_shortfall']
    out['annual'] = []
    for row in result['annual']:
        r = {k: v for k, v in row.items() if k not in DETAIL}
        required = row['required_asset_mass_kg']
        stock = row['available_asset_mass_kg']
        r['capacity_summary'] = dict(binding=max(row['capacity_use'], key=row['capacity_use'].get),
                                     resident_workforce=row['capacity_use']['resident_workforce'])
        r['asset_stocks_kg'] = {k: stock.get(k, 0.) for k in STOCKS}
        r['idle_asset_mass_kg'] = sum(max(0., v - required.get(k, 0.)) for k, v in stock.items())
        out['annual'].append(r)
    return out


def price_opening_route_fuel(case, engineering):
    """Add to the founding cost the first fuel of the Earth-route ships, bought on Earth and launched once."""
    if not case['annual'] or case.get('seed') is None:
        return case
    old = engineering['original_architecture']
    v = old['peak_flows']['Earth_vehicle']
    hull_per_ship = old['assets']['Earth_hull_tanks']['mass_kg'] / old['engineering']['counts']['Earth_route_ships']
    fraction = engineering['transport']['depots']['seed_earth_fleet_fraction']
    I = engineering['normalization']['USD2026_to_USD2041']
    t = engineering['industrial']['transport']
    seed = case['seed']
    if 'priced_initial_route_fuel' in seed:
        raise ValueError('Initial route fuel already priced')
    ships = math.ceil(case['annual'][0]['asset_stocks_kg']['Earth_hull_tanks'] / hull_per_ship - 1e-10)
    fuel = ships * fraction * v['burns']['ascent_kg']
    purchase = fuel * t['earth_propellant_USD2026_kg'] * I
    launch = fuel * model.at(t['LEO_USD2026_kg'], case['annual'][0].get('year', 0)) * I
    seed['priced_initial_route_fuel'] = dict(
        Earth_route_ship_count=ships, seed_fleet_fraction=fraction, fuel_to_LEO_kg=fuel,
        external_propellant_purchase_USD2041=purchase, Earth_launch_USD2041=launch, total_USD2041=purchase + launch,
        formula='Actual opening Earth-route hull count × inherited starting fleet fraction × inherited ascent propellant per ship. Earth propellant purchase and external Earth-to-LEO launch are paid once; no local/source feed is credited or charged for this Earth-provided stock.')
    seed['cost']['initial_route_fuel_USD2041'] = purchase + launch
    seed['cost']['total_USD2041'] += purchase + launch
    seed['remaining_unpriced_requirements'] = [
        'A qualified sequence for the first source/depot deployment, including any fuel beyond the priced initial Earth-route stock, storage filling and commissioning losses.',
        'Site and source access, passenger access and preproduction staff accommodation beyond the average material and payroll screen.',
        'Civil and service-facility details not established by the inherited per-resident, per-staff and per-bed proxies.'
    ]
    return case


# The part of the workforce that the staffing of locally made working parts may fill, and the passes it may take.
STAFFING_MARGIN = 0.97
STAFFING_PASSES = 8


def build(data, spec, book, measures='combined', residual_water_allowance=True):
    """Solve the industry for one book of orders and return (the industry by year, the account of its loops).

    Where the coefficients give a share of working parts the industry is able to make and leave its staffing open,
    the staffing is found here: the industry is solved with every working part bought, the residents' free hands in
    each year are counted, and the share is staffed as far as they reach, year by year, until every year's jobs fit.
    """
    stated = data['economy']['production'].get('local_working_parts', {})
    working = model.local_working_parts(data)
    if working['share'] == 0.0 or 'staffed' in stated:
        physical, loops, _ = solve_once(data, spec, book, measures, residual_water_allowance)
        return physical, loops
    wage = data['industrial']['economics']['loaded_staff_USD2026_year']
    years = [row['year'] for row in book['annual']]

    def with_staffing(staffed):
        changed = deepcopy(data)
        changed['economy']['production']['local_working_parts']['staffed'] = {'by_year': [[y, staffed[y]] for y in years]}
        return changed

    def free_hands(physical, e):
        """For each year, the part of the able share that the residents not otherwise employed can staff."""
        rows = {r['year']: r for r in physical['annual']}
        reach = {}
        for y in years:
            r = rows[y]
            room = STAFFING_MARGIN * model.WORKFORCE_SHARE * r['population'] - (r['technical_resident_jobs'] - r['working_parts_jobs'])
            wanted = e.flows[y]['local_working_labour_if_staffed'] / (wage * model.labour_index(data, y))
            reach[y] = 1. if wanted <= 0 else max(0., min(1., room / wanted))
        # An order is worked in its own year and the next, so both years' hands must reach.
        return {y: min(reach[y], reach.get(y + 1, reach[y])) for y in years}

    staffed = {y: 0. for y in years}
    physical, loops, e = solve_once(with_staffing(staffed), spec, book, measures, residual_water_allowance)
    staffed = free_hands(physical, e)
    for passes in range(1, STAFFING_PASSES + 1):
        physical, loops, e = solve_once(with_staffing(staffed), spec, book, measures, residual_water_allowance)
        over = [r['year'] for r in physical['annual'] if r['working_parts_jobs'] > 0 and r['capacity_summary']['resident_workforce'] > 1 + 1e-9]
        if not over:
            break
        tighter = free_hands(physical, e)
        # free_hands measures against the full able share; only ever lower a year's staffing.
        staffed = {y: min(staffed[y], tighter[y]) for y in years}
    loops['local_working_parts'] = dict(staffed_by_order_year={str(y): staffed[y] for y in years}, staffing_passes=passes,
                                        years_over_the_workforce=over)
    return physical, loops


def solve_once(data, spec, book, measures, residual_water_allowance):
    """One solution of the industry with the coefficients as given: (the industry by year, the account, the measures)."""
    closure, refurbishment = recovery.OPTIONS[measures]
    e = recovery.Recovery(spec, closure, refurbishment, residual_water_allowance)
    configured = e.configure(data)
    result = model.solve(configured, book['annual'], measures=e, max_iterations=50)
    seed_audit = dict(physical_annual_average=deepcopy(result['seed']['physical_annual_average']),
                      late_preproduction_physical_total=deepcopy(e.last_physical),
                      reading_rule='The inherited eight-year average plus the separately dated years 6/7 procurement total; the latter is divided by 8 only to present an average. No source deployment chronology is inferred.')
    before = price_opening_route_fuel(deepcopy(compact(book, result)), configured)
    # Refurbishment changes what is bought, and leaves the mass installed as it was: correct the purchases and add its labour.
    if e.refurbishment:
        I = configured['normalization']['USD2026_to_USD2041']
        r = e.spec['refurbishment']
        cats = model.catalogue(configured)
        worked = model.worked_by_year(result['component_orders'])
        for row in before['annual']:
            y = row['year']
            f = e.flows[y]
            # Workshop labour grows more productive with the rest, and its wage rises in step.
            output_per_job = model.labour_index(configured, y)
            extra_jobs = (sum(f['loop_inspected'].values()) / r['inspection_kg_worker_year'] + sum(f['loop_repaired'].values()) / r['remanufacture_kg_worker_year']) / output_per_job
            extra_payroll = extra_jobs * configured['industrial']['economics']['loaded_staff_USD2026_year'] * output_per_job * I
            row['core_resident_jobs'] += extra_jobs
            row['technical_resident_jobs'] += extra_jobs
            c = row['cash']
            c['resident_payroll_USD2041'] += extra_payroll
            c['operating_local_payroll_USD2041'] += extra_payroll
            c['operating_total_USD2041'] += extra_payroll
            row['capacity_summary']['resident_workforce'] = row['technical_resident_jobs'] / max(1., row['population'] * model.WORKFORCE_SHARE)
            row['maximum_capacity_use'] = max(row['maximum_capacity_use'], row['capacity_summary']['resident_workforce'])
            if row['maximum_capacity_use'] > 1 + 1e-6:
                row['feasible'] = False
                row['actual_sales_capacity_fraction'] = 0.
                row['shortfalls']['resident_workforce'] = row['capacity_summary']['resident_workforce'] - 1
            row['physical']['source_tug_component_import_kg_year'] -= f['loop_import_avoided'].get('asteroid_tugs', 0.)
            drivers = {}
            for ev in worked.get(y, []):
                if ev['component'] == 'active':
                    n = ev['asset']
                    drivers[n] = drivers.get(n, 0.) + ev.get('new_imported_mass_kg', ev['mass_kg'] * (1 - ev['locally_made_share'])) / 2 * cats[n]['price'] * ev.get('price_index', 1.0) * I
            row['top_imported_component_purchases_USD2041'] = dict(sorted(drivers.items(), key=lambda kv: kv[1], reverse=True)[:8])
            row['direct_refurbishment_jobs'] = extra_jobs
            row['direct_refurbishment_payroll_USD2041'] = extra_payroll
        before['first_capacity_shortfall'] = next((r['year'] for r in before['annual'] if not r['feasible']), None)
    rows = []
    bank = {}
    first_year = book['annual'][0]['year']
    for row in before['annual']:
        request = book['annual'][row['year'] - first_year]
        water = recovery.water_account(e.spec, row['population'], request.get('visitor_nights_year', 0.) / configured['constants']['JULIAN_YEAR_DAYS'], e.closure)
        stock = row['physical']['source_water_breakdown_kg_year']['resident_and_guest_water'] - water['operational_makeup_kg']
        if stock < -0.001:
            raise ValueError('Water makeup exceeds the physical delivered water')
        water['rounding_residual_kg'] = min(0., stock)
        water['new_stock_and_construction_kg'] = max(0., stock)
        water['landscape_agriculture_allowance_basis'] = 'Inherited 5 kg/person/day for garden, lake and food-production losses; retained beside explicit farm balances because its decomposition is unknown. Potential overlap is unquantified.'
        water['metabolic_oxygen_recovery'] = 'Not separately resolved: inherited air demand is atmosphere leakage plus new stocks, with food/plant metabolism not represented as an ISS oxygen loop. No additional oxygen credit.'
        materials = e.material_account(configured, row, bank)
        bank = row['physical']['recycled_inventory_kg']
        rows.append(dict(year=row['year'], phase='operation', water=water, materials=materials,
                         source_feed_kg=row['physical']['source_feed_kg_year'], lunar_raw_kg=row['physical']['lunar_raw_kg_year'],
                         source_water_equivalent_kg=row['physical']['source_water_kg_year'], source_nitrogen_kg=row['physical']['source_nitrogen_kg_year'],
                         source_virgin_metal_kg=row['physical']['source_virgin_metal_kg_year'], source_water_breakdown_kg=row['physical']['source_water_breakdown_kg_year'],
                         depot_water_breakdown_kg=row['physical']['depot_water_breakdown_kg_year'],
                         chemical_exhaust_recovery_kg=0., extra_process_kWh=row['physical'].get('loop_extra_process_kWh_year', {})))

    def seed_average(key):
        return seed_audit['physical_annual_average'][key] + seed_audit['late_preproduction_physical_total'][key] / 8
    seedrows = [dict(year=y, phase='preproduction_average', source_water_equivalent_kg=seed_average('source_water_kg_year'),
                     source_nitrogen_kg=seed_average('source_nitrogen_kg_year'), source_virgin_metal_kg=seed_average('source_virgin_metal_kg_year'),
                     source_feed_kg=before['seed']['source_feed_kg'] / 8, lunar_raw_kg=before['seed']['lunar_raw_kg'] / 8,
                     water=None, materials=None, chemical_exhaust_recovery_kg=0.,
                     reading_rule='One eighth of the inherited combined seed feed. Gross circulation, recovered products and deployment losses lack a dated seed sequence; null is not zero.')
                for y in range(8)]
    variant = measures if residual_water_allowance else measures + '_without_residual_allowance'
    loops = dict(name=book['name'] + ':' + variant, baseline_case=book['name'], variant=variant,
                 legacy_landscape_allowance_kg_person_day=5. if residual_water_allowance else 0.,
                 physical_case=before, seed_audit=seed_audit, annual=seedrows + rows,
                 cumulative_source_feed_kg=sum(r['source_feed_kg'] for r in seedrows + rows),
                 cumulative_lunar_raw_kg=sum(r['lunar_raw_kg'] for r in seedrows + rows))
    for key in ('portfolio', 'market_growth_name'):
        if key in book:
            loops[key] = book[key]
    return before, loops, e
