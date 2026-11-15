"""The industry model: its freight, its lifetimes, its cohorts and the closure of its physical flows."""
from copy import deepcopy
import json
from pathlib import Path
import pytest
from provisioning.vacuum_industry import industry, model as m

HERE = Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module')
def data():
    return json.loads((HERE / 'engineering.json').read_text())


@pytest.mark.parametrize('route,trips,burned', [
    ('Earth', 529.7385319490899, 645960632.1573431),
    ('Moon', 1388.5098881093122, 293462569.4859547),
    ('source', 12.777403176831827, 149396259.74518904)])
def test_freight_matches_the_founding_study_controls(data, route, trips, burned):
    # Controls from the founding industrial study's trip solver (commit 146a5582) with its vehicles and tanks,
    # for 100 million kg a year outward and 20 million kg a year back.
    got = m.ferry(1e8, 2e7, m.routes(data, 8)[route])
    assert got['trips'] == pytest.approx(trips, rel=1e-10)
    assert got['burned'] == pytest.approx(burned, rel=1e-10)


def test_no_cargo_means_no_trip(data):
    assert m.ferry(0, 0, m.routes(data, 8)['source'])['burned'] == 0


def test_the_lifetime_distribution_keeps_its_mean():
    life = 20.
    assert sum(m.survival(age, life) for age in range(300)) == pytest.approx(life, rel=1e-12)
    assert m.survival(20, 20, None) == 0


def test_an_order_is_made_over_the_two_years_before_service_and_retires_unreplaced_when_idle(data):
    rows = [dict(year=y, population=0.) for y in range(8, 20)]
    targets = {r['year']: ({'source_computer': 100.} if r['year'] == 10 else {}) for r in rows}
    flows, stocks, events = m.cohort_ledger(data, rows, targets, {})
    chosen = [e for e in events if e['asset'] == 'source_computer']
    assert len(chosen) == 1
    assert chosen[0]['order_year'] == 8 and chosen[0]['commission_year'] == 10
    assert flows[8]['imports']['source'] == pytest.approx(50.)
    assert flows[9]['imports']['source'] == pytest.approx(50.)
    assert stocks[10]['source_computer'] == pytest.approx(100.)
    assert 0 < stocks[19]['source_computer'] < stocks[11]['source_computer'] < 100


def test_the_founding_cohorts_are_bought_gross_of_their_losses_before_opening():
    life, shape, phases, stock = 20., 3., 8, 100.
    denominator = sum(m.survival(j + 1, life, shape) for j in range(phases))
    purchased = stock * phases / denominator
    assert purchased > stock
    assert sum(purchased / phases * m.survival(j + 1, life, shape) for j in range(phases)) == pytest.approx(stock)


def test_the_source_feed_is_the_largest_of_its_three_needs(data):
    row = dict(year=8, population=12500., source_goods=1e7, new_residential_people=0.)
    f = m.physical(data, row, m.blank_flow())
    z = data['industrial']['source']
    needs = (f['source_water_kg_year'] / (z['water_grade'] * z['water_recovery']),
             f['source_nitrogen_kg_year'] / (z['nitrogen_grade'] * z['nitrogen_recovery']),
             f['source_virgin_metal_kg_year'] / (z['metal_grade'] * z['metal_recovery']))
    assert f['source_feed_kg_year'] == pytest.approx(max(needs))
    assert sum(f['source_cargo_breakdown_kg_year'].values()) == pytest.approx(f['source_cargo_kg_year'])
    assert sum(f['source_water_breakdown_kg_year'].values()) == pytest.approx(f['source_water_kg_year'])
    assert sum(f['depot_water_breakdown_kg_year'].values()) == pytest.approx(f['depot_water_kg_year'])


def test_a_poorer_source_needs_more_feed_and_more_power_to_work_it(data):
    row = dict(year=8, population=12500., new_residential_people=0.)
    base = m.physical(data, row, m.blank_flow())
    poor = deepcopy(data)
    poor['industrial']['source']['water_grade'] /= 2
    stress = m.physical(poor, row, m.blank_flow())
    assert stress['source_feed_kg_year'] > base['source_feed_kg_year']
    assert stress['source_feed_process_W'] > base['source_feed_process_W']
    assert stress['source_product_refining_W'] == base['source_product_refining_W']


def test_the_settlement_scales_with_residents_and_the_industry_with_its_flows(data):
    assert all(data['templates'].get(n, {'owner': 'settlement'})['owner'] == 'settlement' for n in m.residential_targets(data, 250000.))
    row = dict(year=8, population=250000., new_residential_people=0.)
    need = m.requirements(data, row, m.physical(data, row, m.blank_flow()))
    assert need['asteroid_tugs'] < data['original_architecture']['assets']['asteroid_tugs']['mass_kg']


def test_recycled_material_is_used_before_fresh_and_the_rest_is_kept(data):
    flow = m.blank_flow()
    flow['local']['source'] = 100.
    flow['recycled']['source'] = 150.
    p = m.physical(data, dict(year=8, population=0.), flow)
    assert p['recycled_inventory_kg']['source'] == 50.
    assert p['source_virgin_metal_kg_year'] == 0.


def test_computing_adds_power_and_radiators_and_no_chemical_plant(data):
    row = dict(year=8, population=12500., new_residential_people=0.)
    digital_row = dict(row, compute=1e9)
    a = m.requirements(data, row, m.physical(data, row, m.blank_flow()))
    b = m.requirements(data, digital_row, m.physical(data, digital_row, m.blank_flow()))
    assert b['source_chemistry'] == pytest.approx(a['source_chemistry'])
    assert b['source_solar'] > a['source_solar']
    assert b['source_radiators'] > a['source_radiators']


def test_the_opening_route_fuel_is_priced_once_and_draws_no_source_feed(data):
    old = data['original_architecture']
    unit = old['assets']['Earth_hull_tanks']['mass_kg'] / old['engineering']['counts']['Earth_route_ships']
    case = dict(seed=dict(cost=dict(total_USD2041=10.), source_feed_kg=123.), annual=[dict(asset_stocks_kg=dict(Earth_hull_tanks=3 * unit))])
    seed = industry.price_opening_route_fuel(case, data)['seed']
    fuel = seed['priced_initial_route_fuel']
    assert fuel['Earth_route_ship_count'] == 3
    assert seed['source_feed_kg'] == 123.
    assert seed['cost']['total_USD2041'] == pytest.approx(10. + fuel['external_propellant_purchase_USD2041'] + fuel['Earth_launch_USD2041'])
    with pytest.raises(ValueError, match='already priced'):
        industry.price_opening_route_fuel(case, data)


def test_measures_are_called_at_the_three_points_and_nowhere_else(data):
    class Count:
        def __init__(self):
            self.calls = dict(ledger=0, physical=0, requirements=0)

        def after_ledger(self, data, rows, seed, flows, events):
            self.calls['ledger'] += 1

        def after_physical(self, data, row, flow, out):
            self.calls['physical'] += 1

        def after_requirements(self, data, row, physical, out):
            self.calls['requirements'] += 1
    rows = [dict(year=y, population=12500., new_residential_people=0.) for y in range(8, 14)]
    plain = m.solve(data, rows, max_iterations=50)
    counted = Count()
    watched = m.solve(data, rows, measures=counted, max_iterations=50)
    assert watched == plain
    assert counted.calls['ledger'] == plain['convergence']['iterations'] + 1
    assert counted.calls['physical'] > counted.calls['requirements'] > counted.calls['ledger']


def test_a_scheduled_coefficient_interpolates_numbers_and_steps_names():
    price = {'by_year': [[8, 200.], [58, 100.], [108, 50.]]}
    assert m.at(150., 30) == 150.
    assert m.at(price, 0) == 200. and m.at(price, 8) == 200.
    assert m.at(price, 33) == pytest.approx(150.) and m.at(price, 58) == 100.
    assert m.at(price, 500) == 50.
    leg = {'by_year': [[0, 'split'], [20, 'space']]}
    assert [m.at(leg, y) for y in (0, 19, 20, 200)] == ['split', 'split', 'space', 'space']


def with_leg(data, leg, **transport):
    changed = deepcopy(data)
    changed['industrial']['transport'].update(earth_leg=leg, **transport)
    return changed


def test_the_space_arrangement_keeps_the_founding_ferry_and_the_others_close_on_their_tanks(data):
    v = data['original_architecture']['peak_flows']['Earth_vehicle']
    space = m.earth_route(data, 8)
    assert space['leg'] == 'space' and space['dry'] == v['dry_kg'] and space['dry_scale'] == 1.
    assert space['tank_limit'] == v['burns']['peak_tank_contents_kg']
    fixed = data['transport']['chemical']['fixed_dry_to_payload'] * space['payload']
    per_stored = (v['dry_kg'] - fixed) / v['burns']['peak_tank_contents_kg']
    for leg in ('split', 'earth'):
        route = m.earth_route(with_leg(data, leg), 8)
        assert route['dry'] == pytest.approx(fixed + per_stored * route['tank_limit'], rel=1e-10)
        assert fixed < route['dry'] < space['dry']
        trip = m.burns(route['payload'], 0., route['dry'], route['dv'], route['isp'], route['margin'], route['mode'], route['retention'])
        assert trip['tank'] == pytest.approx(route['tank_limit'])
    with pytest.raises(ValueError):
        m.earth_route(with_leg(data, 'towed'), 8)


def test_a_split_trip_fuels_each_burn_where_it_begins_and_is_limited_by_the_heavier_leg(data):
    route = m.earth_route(with_leg(data, 'split'), 8)
    k = m.math.expm1(route['dv'] / (route['isp'] * m.STANDARD_GRAVITY)) * route['margin']
    trip = m.burns(100., 40., 30., route['dv'], route['isp'], route['margin'], 'split_chemical', route['retention'])
    assert trip['out'] == pytest.approx(130. * k) and trip['back'] == pytest.approx(70. * k)
    assert trip['burned'] == pytest.approx(trip['out'] + trip['back']) and trip['tank'] == max(trip['out'], trip['back'])
    assert trip['space_supply'] == pytest.approx(trip['back'] / route['retention'])
    outward = m.ferry(1e9, 2e8, route)
    homeward = m.ferry(2e8, 1e9, route)
    assert outward['trips'] == pytest.approx(1e9 / route['payload']) == homeward['trips']


def test_each_arrangement_draws_the_earth_legs_propellant_from_where_it_says(data):
    row = dict(year=8, population=250000., goods=5e7, new_residential_people=0.)
    flow = m.blank_flow()
    flow['imports']['surface'] = 2e8
    got = {leg: m.physical(with_leg(data, leg), row, flow) for leg in ('space', 'split', 'earth')}
    price = data['industrial']['transport']
    for leg, p in got.items():
        cargo = sum(p['imports_by_destination_kg_year'].values())
        launched = p['Earth_leg_propellant_launched_kg_year']
        assert p['LEO_kg_year'] == pytest.approx(cargo + launched)
        assert p['Earth_launch_USD2026'] == pytest.approx(p['LEO_kg_year'] * price['LEO_USD2026_kg'] + launched * price['earth_propellant_USD2026_kg'])
        assert sum(p['depot_water_breakdown_kg_year'].values()) == pytest.approx(p['depot_water_kg_year'])
    assert got['space']['Earth_leg_propellant_launched_kg_year'] == 0
    assert got['space']['depot_water_breakdown_kg_year']['Earth_ferry_fuel'] == got['space']['Earth']['space_supply']
    assert got['split']['Earth_leg_propellant_launched_kg_year'] == got['split']['Earth']['out']
    assert got['split']['depot_water_breakdown_kg_year']['Earth_ferry_fuel'] == got['split']['Earth']['space_supply']
    assert got['earth']['Earth_leg_propellant_launched_kg_year'] == got['earth']['Earth']['burned']
    assert got['earth']['depot_water_breakdown_kg_year']['Earth_ferry_fuel'] == 0
    # Less propellant from the depot means less water to fetch and less feed to mine.
    assert got['earth']['source_feed_kg_year'] < got['split']['source_feed_kg_year'] < got['space']['source_feed_kg_year']
    need = {leg: m.requirements(with_leg(data, leg), row, p) for leg, p in got.items()}
    assert need['earth']['depot_route_storage_controls'] == 0 < need['split']['depot_route_storage_controls'] < need['space']['depot_route_storage_controls']
    assert need['earth']['asteroid_tugs'] < need['split']['asteroid_tugs'] < need['space']['asteroid_tugs']


def test_the_arrangement_and_the_launch_price_can_change_with_the_year(data):
    staged = with_leg(data, {'by_year': [[0, 'split'], [20, 'space']]}, LEO_USD2026_kg={'by_year': [[8, 200.], [108, 100.]]})
    early = dict(year=12, population=12500., new_residential_people=0.)
    late = dict(early, year=108)
    a, b = m.physical(staged, early, m.blank_flow()), m.physical(staged, late, m.blank_flow())
    assert a['Earth_leg_propellant_launched_kg_year'] > 0 == b['Earth_leg_propellant_launched_kg_year']
    assert m.routes(staged, 12)['Earth']['leg'] == 'split' and m.routes(staged, 20)['Earth']['leg'] == 'space'
    assert b['Earth_launch_USD2026'] == pytest.approx(b['LEO_kg_year'] * 100.)


def with_working_parts(data, **stated):
    changed = deepcopy(data)
    changed['economy']['production']['local_working_parts'].update(stated)
    return changed


def test_a_locally_made_share_of_working_parts_moves_mass_and_money_from_imports_to_local_work(data):
    rows = [dict(year=y, population=0.) for y in range(8, 20)]
    targets = {r['year']: ({'source_computer': 100., 'return_capsule_buses': 40.} if r['year'] == 10 else {}) for r in rows}
    cats = m.catalogue(data)
    price = cats['source_computer']['price']
    bought, _, plain = m.cohort_ledger(data, rows, targets, {})
    half = with_working_parts(data, share=0.5, labour_fraction_of_import_price=0.2)
    made, _, events = m.cohort_ledger(half, rows, targets, {})
    for year in (8, 9):
        assert made[year]['imports']['source'] == pytest.approx(bought[year]['imports']['source'] / 2)
        assert made[year]['local']['source'] == pytest.approx(25.) and bought[year]['local']['source'] == 0
        assert made[year]['external_capital'] == pytest.approx(bought[year]['external_capital'] - 25. * price)
        assert made[year]['local_working_capital'] == pytest.approx(25. * price * 0.2)
        assert made[year]['local_working_labour_if_staffed'] == pytest.approx(50. * 0.5 * price * 0.2)
        # The capsules that return goods to Earth are bought there whatever the share.
        assert made[year]['imports']['leo'] == bought[year]['imports']['leo'] == pytest.approx(20.)
    by = {e['asset']: e for e in events}
    assert by['source_computer']['locally_made_share'] == 0.5 and by['return_capsule_buses']['locally_made_share'] == 0
    assert all(e['locally_made_share'] == 0 for e in plain)


def test_the_staffed_part_scales_the_share_made_and_leaves_the_labour_of_the_full_share_as_it_was(data):
    rows = [dict(year=y, population=0.) for y in range(8, 20)]
    targets = {r['year']: ({'source_computer': 100.} if r['year'] == 10 else {}) for r in rows}
    full, _, _ = m.cohort_ledger(with_working_parts(data, share=0.6, labour_fraction_of_import_price=0.2), rows, targets, {})
    part, _, events = m.cohort_ledger(with_working_parts(data, share=0.6, staffed=0.25, labour_fraction_of_import_price=0.2), rows, targets, {})
    assert events[0]['locally_made_share'] == pytest.approx(0.15)
    assert part[8]['local']['source'] == pytest.approx(full[8]['local']['source'] / 4)
    assert part[8]['local_working_labour_if_staffed'] == full[8]['local_working_labour_if_staffed'] > 0
    with pytest.raises(ValueError):
        m.cohort_ledger(with_working_parts(data, share=1.4), rows, targets, {})


def test_local_working_parts_are_paid_as_resident_labour_and_the_accounts_still_close(data):
    rows = [dict(year=y, population=12500. * 1.033**(y - 8), new_residential_people=500., goods=1e5, compute=1e6) for y in range(8, 30)]
    bought = m.solve(data, rows, max_iterations=50)
    made = m.solve(with_working_parts(data, share={'by_year': [[8, 0.], [18, 0.5]]}, labour_fraction_of_import_price=0.1), rows, max_iterations=50)
    wage = data['industrial']['economics']['loaded_staff_USD2026_year'] * data['normalization']['USD2026_to_USD2041']
    for before, after in zip(bought['annual'], made['annual']):
        cash = after['cash']
        assert abs(cash['accounting_residual_USD2041']) < .02
        assert before['working_parts_jobs'] == 0 and before['cash']['operating_local_replacement_USD2041'] == 0
        assert after['core_resident_jobs'] * wage == pytest.approx(cash['resident_payroll_USD2041'])
        assert cash['operating_local_replacement_USD2041'] >= 0
    late = [(b, a) for b, a in zip(bought['annual'], made['annual']) if a['year'] >= 22]
    assert all(a['working_parts_jobs'] > 0 for _, a in late)
    assert sum(a['cash']['operating_external_replacement_USD2041'] for _, a in late) < sum(b['cash']['operating_external_replacement_USD2041'] for b, _ in late)
    assert sum(a['physical']['LEO_kg_year'] for _, a in late) < sum(b['physical']['LEO_kg_year'] for b, _ in late)
