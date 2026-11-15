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
    got = m.ferry(1e8, 2e7, m.routes(data)[route])
    assert got['trips'] == pytest.approx(trips, rel=1e-10)
    assert got['burned'] == pytest.approx(burned, rel=1e-10)


def test_no_cargo_means_no_trip(data):
    assert m.ferry(0, 0, m.routes(data)['source'])['burned'] == 0


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
