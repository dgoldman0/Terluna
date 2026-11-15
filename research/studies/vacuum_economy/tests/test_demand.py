"""The customers' demand at each year's date: the world it is judged in, its lines of sale, and its answer to a price."""
from copy import deepcopy
import math
import pytest
from shared.constants import JULIAN_YEAR_DAYS
from research.studies.vacuum_economy import demand, programme as definition, run

SCENARIOS = ('low', 'central', 'high')


@pytest.fixture(scope='module')
def s():
    return run.sources()


def book(s, name='central', ramp='constant_growth', **stated):
    return demand.book(s['programme'], ramp, s['scenarios'], name, s['engineering'], **stated)


def rows(b):
    return {row['year']: row for row in b['annual']}


def world(s, name='central'):
    calendar = s['programme']['calendar']
    return demand.World(s['scenarios'], s['scenarios']['scenarios'][name], calendar['year_zero'], calendar['first_productive_year'],
                        s['engineering']['normalization']['USD2026_to_USD2041'])


def test_the_world_follows_the_united_nations_projection_and_is_held_where_it_ends(s):
    w = world(s)
    assert w.population(2050) == 9.664e9 and w.population(2084) == 10.289e9
    assert w.population(2100) == w.population(2241) == 10.180e9
    assert w.population(2062) == pytest.approx(9.664e9 + (10.250e9 - 9.664e9) * 12 / 25)


def test_income_grows_at_the_projected_rate_to_2049_and_at_the_scenarios_rate_after(s):
    for name, rate in (('low', 0.005), ('central', 0.01), ('high', 0.016)):
        w = world(s, name)
        assert w.income(2049) == pytest.approx(1.016 ** 25)
        assert w.income(2149) / w.income(2049) == pytest.approx((1 + rate) ** 100)
    w = world(s)
    assert w.output(2050) == pytest.approx(163.437e12 * 1.404 * s['engineering']['normalization']['USD2026_to_USD2041'])
    assert w.output(2150) / w.output(2050) == pytest.approx(w.population(2150) / w.population(2050) * 1.01 ** 100)


@pytest.mark.parametrize('name', SCENARIOS)
def test_what_customers_would_buy_does_not_depend_on_the_residents_or_the_founding_system(s, name):
    a, b = rows(book(s, name, 'constant_growth')), rows(book(s, name, 'front_loaded'))
    doubled = deepcopy(s['programme'])
    doubled['population']['founding_system'] *= 2
    c = rows(demand.book(doubled, 'constant_growth', s['scenarios'], name, s['engineering']))
    for year in (8, 33, 100, 200):
        for line in demand.LINES:
            assert a[year]['lines'][line] == b[year]['lines'][line] == c[year]['lines'][line]
        assert a[year]['population'] != b[year]['population'] or year in (8, 200)


@pytest.mark.parametrize('name', SCENARIOS)
def test_the_receipts_of_a_year_are_the_receipts_of_its_lines(s, name):
    for row in book(s, name, feed={y: 1e12 for y in range(8, 221)})['annual']:
        by_line = math.fsum(detail['receipts_USD2041'] for detail in row['lines'].values())
        assert row['gross_USD2041'] == pytest.approx(by_line) == pytest.approx(math.fsum(row['revenue_by_sector_USD2041'].values()))
        assert set(row['revenue_by_sector_USD2041']) == set(demand.SECTORS)
        assert row['commercial_operating_USD2041'] >= row['commercial_payroll_USD2041'] > 0


@pytest.mark.parametrize('name', SCENARIOS)
def test_cultural_receipts_grow_from_the_opening_portfolio_to_the_mature_following(s, name):
    stated = s['scenarios']['scenarios'][name]['culture']
    r, w = rows(book(s, name)), world(s, name)
    opening = s['scenarios']['lines']['culture']['opening_receipts_USD2041_year']['value']
    assert r[8]['lines']['culture']['receipts_USD2041'] == pytest.approx(opening)
    assert r[9]['lines']['culture']['receipts_USD2041'] == pytest.approx(opening * (1 + stated['growth_of_receipts']['value']))
    assert r[200]['lines']['culture']['quantity'] == pytest.approx(stated['followers_share_of_world']['value'] * w.population(2241))
    per_follower = r[200]['lines']['culture']['price_USD2041'] / r[100]['lines']['culture']['price_USD2041']
    assert per_follower == pytest.approx(w.income(2241) / w.income(2141))
    assert r[200]['revenue_by_sector_USD2041']['sports'] == pytest.approx(r[200]['revenue_by_sector_USD2041']['entertainment'])


def test_hosted_work_is_a_share_of_world_output_and_is_staffed_from_what_the_settlement_keeps(s):
    stated = s['scenarios']['scenarios']['central']['institutions']
    row, w = rows(book(s))[100], world(s)
    hosted = row['lines']['institutions']
    assert hosted['receipts_USD2041'] == pytest.approx(stated['share_of_world_output']['value'] * w.output(2141))
    kept = hosted['receipts_USD2041'] * (1 - stated['share_others_keep']['value'])
    assert hosted['quantity'] == pytest.approx(kept / stated['retained_USD2041_per_job']['value'])
    assert row['service_jobs'] > hosted['quantity']                      # the cultural staff beside them


@pytest.mark.parametrize('name', SCENARIOS)
def test_propellant_and_material_answer_the_price_and_never_pass_the_ceiling(s, name):
    elasticity = s['scenarios']['scenarios'][name]['space']['price_elasticity']['value']
    ceiling = rows(book(s, name))[50]['lines']['propellant']['ceiling_USD2041']
    at = lambda price: rows(book(s, name, prices={'propellant': {50: price}}))[50]['lines']
    assert at(ceiling / 2)['propellant']['quantity'] / at(ceiling)['propellant']['quantity'] == pytest.approx(2 ** elasticity)
    assert at(ceiling * 3)['propellant'] == at(ceiling)['propellant']
    assert at(ceiling / 2)['materials'] == at(ceiling)['materials']        # each line has its own price
    assert rows(book(s, name))[50]['propellant_sales'] == at(ceiling)['propellant']['quantity']
    assert rows(book(s, name))[50]['space_source_goods'] == at(ceiling)['materials']['quantity']


def test_the_ceiling_is_the_stated_share_of_a_kilogram_launched_from_earth_to_the_depot(s):
    transport = s['engineering']['industrial']['transport']
    to_2041 = s['engineering']['normalization']['USD2026_to_USD2041']
    delivered = demand.earth_delivered_price(s['engineering'], 8, to_2041)
    # More than the kilogram itself is launched: the propellant of the ferry's outward and return burns as well.
    assert delivered > 2 * transport['LEO_USD2026_kg'] * to_2041
    cheaper = deepcopy(s['engineering'])
    cheaper['industrial']['transport']['LEO_USD2026_kg'] = transport['LEO_USD2026_kg'] / 2
    assert demand.earth_delivered_price(cheaper, 8, to_2041) < delivered
    ceiling = rows(book(s))[8]['lines']['propellant']['ceiling_USD2041']
    assert ceiling == pytest.approx(s['scenarios']['lines']['space']['share_of_Earth_delivered_price']['value'] * delivered)


@pytest.mark.parametrize('name', SCENARIOS)
def test_visitors_answer_the_price_and_bring_their_beds_and_their_mass(s, name):
    lines = s['scenarios']['lines']['visits']
    elasticity = s['scenarios']['scenarios'][name]['visits']['price_elasticity']['value']
    first = lines['first_price_USD2041']['value']
    row = rows(book(s, name))[60]
    cheaper = rows(book(s, name, prices={'visits': {60: first / 2}}))[60]
    assert cheaper['visitors_year'] / row['visitors_year'] == pytest.approx(2 ** elasticity)
    assert row['visitor_nights_year'] == pytest.approx(row['visitors_year'] * lines['stay_nights']['value'])
    assert row['visitor_beds'] == math.ceil(row['visitor_nights_year'] / (JULIAN_YEAR_DAYS * lines['occupancy']['value']))
    assert row['passenger_launch_kg_year'] == pytest.approx(row['visitors_year'] * lines['launched_kg_per_visitor']['value'])
    assert row['passenger_cabin_kg_year'] == pytest.approx(row['visitors_year'] * lines['cabin_kg_per_visitor']['value'])
    assert rows(book(s, name))[200]['visitors_year'] > row['visitors_year']


def test_metals_come_only_from_feed_that_is_mined_and_their_price_falls_as_the_supply_grows(s):
    stated = s['scenarios']['scenarios']['central']['metals']
    assert all(row['source_goods'] == 0 for row in book(s)['annual'])
    scarce = rows(book(s, feed={100: 1e9}))[100]['lines']['metals']
    assert scarce['treated_feed_kg'] == pytest.approx(1e9)
    assert scarce['quantity'] == pytest.approx(1e9 * stated['refined_PGM_per_kg_feed']['value'] * (1 + s['scenarios']['lines']['metals']['gold_per_PGM']['value']))
    ample = rows(book(s, feed={100: 1e15}))[100]['lines']['metals']
    assert ample['treated_feed_kg'] < 1e15                                # the scenario's share of the world's supply, no more
    assert ample['price_USD2041'] < scarce['price_USD2041']
    assert all(row['source_goods'] == 0 for row in book(s, 'low', feed={y: 1e15 for y in range(8, 221)})['annual'])


def test_a_line_can_be_served_in_part_or_left_out_and_the_others_stay_as_they_were(s):
    whole = rows(book(s))[100]
    half = rows(book(s, scale={'culture': {100: 0.5}, 'visits': {100: 0.0}}))[100]
    assert half['lines']['culture']['receipts_USD2041'] == pytest.approx(whole['lines']['culture']['receipts_USD2041'] / 2)
    assert half['visitors_year'] == half['visitor_beds'] == half['passenger_cabin_kg_year'] == 0
    for line in ('institutions', 'computing', 'propellant', 'materials'):
        assert half['lines'][line] == whole['lines'][line]
    nothing = rows(book(s, scale={line: {100: 0.0} for line in demand.LINES}))[100]
    assert nothing['gross_USD2041'] == nothing['service_jobs'] == nothing['commercial_operating_USD2041'] == 0


def test_an_unknown_scenario_line_or_share_is_refused(s):
    with pytest.raises(ValueError, match='Unknown scenario'):
        book(s, 'first_pass')
    with pytest.raises(ValueError, match='No such line'):
        book(s, scale={'tourism': {8: 1.0}})
    with pytest.raises(ValueError, match='No such line'):
        book(s, prices={'culture': {8: 1.0}})
    with pytest.raises(ValueError, match='between none and all'):
        book(s, scale={'culture': {8: 1.5}})


def test_the_sellers_price_follows_the_elasticity(s):
    to_2041 = s['engineering']['normalization']['USD2026_to_USD2041']
    inputs = s['scenarios']
    assert demand.seller_price(inputs, 'high', 'propellant', 10.0, to_2041) == pytest.approx(10.0 * 1.4 / 0.4)
    assert demand.seller_price(inputs, 'high', 'visits', 1e5, to_2041) == pytest.approx(2e5)
    assert demand.seller_price(inputs, 'central', 'propellant', 10.0, to_2041) is None
    assert demand.seller_price(inputs, 'low', 'materials', 10.0, to_2041) is None
    anchor = inputs['lines']['visits']['anchor']['price_USD2022'] * inputs['money']['USD2022_to_USD2026'] * to_2041
    assert demand.seller_price(inputs, 'low', 'visits', 1e5, to_2041) == pytest.approx(anchor)


def test_every_scenario_states_a_basis_for_each_value_and_changes_only_coefficients_that_exist(s):
    for name, scenario in s['scenarios']['scenarios'].items():
        for line in ('culture', 'institutions', 'computing', 'space', 'metals', 'visits'):
            assert all(entry['basis'] for entry in scenario[line].values())
        case = dict(scenario=name, ramp='constant_growth', recovery='combined')
        changed = run.engineering(case)
        assert changed['industrial']['transport']['LEO_USD2026_kg'] == scenario['engineering']['industrial.transport.LEO_USD2026_kg']['value']
        assert s['engineering']['industrial']['transport']['LEO_USD2026_kg'] == 200


def test_a_cost_plus_price_is_the_cost_with_the_stated_margin_and_an_unknown_rule_is_refused(s):
    to_2041 = s['engineering']['normalization']['USD2026_to_USD2041']
    margin = s['scenarios']['lines']['cost_plus_margin']['value']
    for name in SCENARIOS:
        for line in demand.PRICED:
            assert demand.seller_price(s['scenarios'], name, line, 40.0, to_2041, 'cost_plus') == pytest.approx(40.0 * (1 + margin))
    with pytest.raises(ValueError, match='Unknown rule of pricing'):
        demand.seller_price(s['scenarios'], 'central', 'visits', 40.0, to_2041, 'at_cost')
