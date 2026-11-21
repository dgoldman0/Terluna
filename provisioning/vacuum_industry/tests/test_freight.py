"""Freight to the sources: the tugs' propulsion, the leg they fly, and the leg that follows the water drawn."""
from copy import deepcopy
import json
from pathlib import Path
import pytest
from provisioning.vacuum_industry import industry, model as m
from shared.constants import JULIAN_DAY, STANDARD_GRAVITY

HERE = Path(__file__).resolve().parents[1]


@pytest.fixture(scope='module')
def data():
    return json.loads((HERE / 'engineering.json').read_text())


def tugs(data, **stated):
    d = deepcopy(data)
    d['industrial']['transport'].update(stated)
    return d


def small_book(years=30):
    rows = []
    for year in range(8, 8 + years):
        people = 12500. * 1.033 ** (year - 8)
        rows.append(dict(year=year, population=people, new_residential_people=people * 1.033 ** 2 * .033,
                         goods=2e5 * (year - 7), propellant_sales=5e5 * (year - 7), compute=2e6))
    return dict(name='small', annual=rows)


def curve(step_kg=2e8):
    """A survey of water by leg: none within 600 m/s, then the same water in each 100 m/s to 20 km/s."""
    edges = [float(v) for v in range(100, 20001, 100)]
    within = [max(0.0, (v - 600.0) / 100.0) * step_kg for v in edges]
    return dict(upper_edges_m_s=edges, within_kg=within + [within[-1]])


def test_the_founding_tug_flies_its_design_leg_unchanged(data):
    route = m.routes(data, 30)['source']
    t = data['industrial']['transport']
    assert route['propulsion'] == 'chemical' and route['dv'] == t['nea_dv_oneway_m_s']
    assert route['dry'] == t['nea_payload_kg'] * t['nea_dry_to_payload']
    assert route['tank_limit'] == t['nea_payload_kg'] * data['original_architecture']['peak_flows']['nea_round_trip']['total'] * t['propellant_margin']
    # The same tug sized through the closure at its design leg is the founding tug.
    dry, tank = m.chemical_tug(t['nea_payload_kg'], route['dry'], route['tank_limit'], 2500., 2500., t['Isp_s'], t['propellant_margin'], t['tank_dry_to_contents'])
    assert dry == pytest.approx(route['dry'], rel=1e-12) and tank == pytest.approx(route['tank_limit'], rel=1e-12)


def test_a_chemical_tug_grows_its_tanks_with_the_leg_and_cannot_close_beyond_its_limit(data):
    legs = [1500., 2500., 3500., 4500., 5500.]
    routes = [m.routes(tugs(data, source_dv_oneway_m_s=v), 30)['source'] for v in legs]
    assert all(b['dry'] > a['dry'] and b['tank_limit'] > a['tank_limit'] for a, b in zip(routes, routes[1:]))
    with pytest.raises(m.Unreachable) as caught:
        m.routes(tugs(data, source_dv_oneway_m_s=7000.), 31)
    assert caught.value.year == 31


@pytest.mark.parametrize('kind', ['electric_water', 'electric_water_target', 'electric_xenon'])
def test_an_electric_tug_closes_on_its_parts_and_its_power_expels_its_propellant(data, kind):
    route = m.routes(tugs(data, source_propulsion=kind, source_dv_oneway_m_s=2000.), 30)['source']
    option = data['transport']['source_propulsion']['options'][kind]
    e = data['transport']['electric']
    assert route['dv'] == pytest.approx(2000. * data['transport']['source_propulsion']['low_thrust_dv_factor'])
    assert sum(route['parts'].values()) == pytest.approx(route['dry'], rel=1e-12)
    full = m.burns(route['payload'], 0.0, route['dry'], route['dv'], route['isp'], route['margin'], route['mode'], route['retention'])
    ve = option['Isp_s'] * STANDARD_GRAVITY
    assert route['bus_W'] * 2 * option['total_efficiency'] * e['source_powered_days'] * JULIAN_DAY == pytest.approx(full['burned'] * ve ** 2, rel=1e-9)
    chemical = m.routes(tugs(data, source_dv_oneway_m_s=2000.), 30)['source']
    assert m.ferry(1e9, 0., route)['burned'] < m.ferry(1e9, 0., chemical)['burned']


def test_propulsion_and_leg_follow_their_schedules(data):
    d = tugs(data, source_propulsion={'by_year': [[8, 'chemical'], [26, 'electric_water_target']]}, source_dv_oneway_m_s={'by_year': [[8, 1000.], [108, 3000.]]})
    assert m.routes(d, 25)['source']['propulsion'] == 'chemical'
    assert m.routes(d, 26)['source']['propulsion'] == 'electric_water_target'
    assert m.routes(d, 58)['source']['dv'] == pytest.approx(2000. * d['transport']['source_propulsion']['low_thrust_dv_factor'])
    assert m.propulsions(d) == {'chemical', 'electric_water_target'}
    assert {'electric_tug_' + p for p in m.ELECTRIC_PARTS} <= set(m.catalogue(d))
    assert not {'electric_tug_' + p for p in m.ELECTRIC_PARTS} & set(m.catalogue(data))


def test_the_leg_reads_the_water_drawn_nearest_first():
    reach = m.source_reach(dict(source_reach=curve()))
    assert reach[0][0] == 600. and reach[1][0] == 0.
    assert m.leg_for(reach, 0.0) == 600.
    assert m.leg_for(reach, 1e8) == pytest.approx(650.)
    with pytest.raises(ValueError):
        m.leg_for(reach, 1e30)


def test_water_burned_by_an_electric_tug_needs_no_electrolysis_and_xenon_comes_from_earth(data):
    row = small_book()['annual'][5]
    f = m.blank_flow()
    # With the ferries' outward burn fuelled on Earth; fuelled at the depot, each kilogram of xenon carried up costs
    # more depot water than the xenon tugs can fetch with it, and the freight does not close.
    plain = m.physical(tugs(data, source_dv_oneway_m_s=2000., earth_leg='split'), row, f)
    water = m.physical(tugs(data, source_propulsion='electric_water_target', source_dv_oneway_m_s=2000., earth_leg='split'), row, f)
    xenon = m.physical(tugs(data, source_propulsion='electric_xenon', source_dv_oneway_m_s=2000., earth_leg='split'), row, f)
    with pytest.raises(m.Unreachable):
        m.physical(tugs(data, source_propulsion='electric_xenon', source_dv_oneway_m_s=2000., earth_leg='space'), row, f)
    assert water['source_water_breakdown_kg_year']['source_ferry_fuel'] == pytest.approx(water['source_tugs']['propellant_kg_year'] / data['transport']['electric']['handling_retention'])
    assert water['source_water_breakdown_kg_year']['source_ferry_fuel'] < plain['source_water_breakdown_kg_year']['source_ferry_fuel']
    assert xenon['source_water_breakdown_kg_year']['source_ferry_fuel'] == 0.0
    assert xenon['source_tugs']['xenon_from_Earth_kg_year'] > 0.0
    assert xenon['imports_by_destination_kg_year']['depot'] == pytest.approx(plain['imports_by_destination_kg_year']['depot'] + xenon['source_tugs']['xenon_from_Earth_kg_year'])
    price = data['transport']['source_propulsion']['options']['electric_xenon']['price_USD2026_kg']
    assert xenon['operating_material_USD2026'] == pytest.approx(plain['operating_material_USD2026'] + price * xenon['source_tugs']['xenon_from_Earth_kg_year'])


def test_without_a_reach_the_tugs_fly_the_founding_leg_and_with_one_the_legs_lengthen_as_water_is_drawn(data):
    book = small_book()
    plain = m.solve(data, book['annual'], max_iterations=40)
    assert all(r['physical']['source_tugs']['leg_m_s'] == data['industrial']['transport']['nea_dv_oneway_m_s'] for r in plain['annual'])
    assert 'leg_moves_m_s' not in plain['convergence']
    reached = m.solve(dict(data, source_reach=curve(5e10)), book['annual'], max_iterations=40)
    assert reached['convergence']['converged'] and reached['convergence']['largest_leg_move_m_s'] <= m.LEG_TOLERANCE_M_S
    legs = [r['physical']['source_tugs']['leg_m_s'] for r in reached['annual']]
    assert all(b >= a for a, b in zip(legs, legs[1:])) and legs[-1] > legs[0] >= 600.
    reach = m.source_reach(dict(source_reach=curve(5e10)))
    grade = data['industrial']['source']['water_grade']
    # The founding deployment flies the leg of half the water it draws.
    assert reached['seed']['source_leg_m_s'] == pytest.approx(m.leg_for(reach, 0.5 * reached['seed']['source_feed_kg'] * grade), abs=m.LEG_TOLERANCE_M_S)
    assert 'source_leg_m_s' not in plain['seed']
    drawn = reached['seed']['source_feed_kg'] * grade
    for r in reached['annual']:
        this = r['physical']['source_feed_kg_year'] * grade
        assert r['physical']['source_tugs']['leg_m_s'] == pytest.approx(m.leg_for(reach, drawn + this / 2), abs=m.LEG_TOLERANCE_M_S)
        drawn += this


def test_the_staffing_starts_from_the_able_share_when_freight_cannot_fly_with_every_part_bought(data, monkeypatch):
    calls = []
    real = industry.solve_once

    def solve_once(d, spec, book, measures, allowance):
        staffed = d['economy']['production']['local_working_parts']['staffed']['by_year']
        calls.append(staffed[0][1])
        if not calls[1:]:
            raise m.Unreachable('not with every part bought', 20)
        return real(d, spec, book, measures, allowance)
    monkeypatch.setattr(industry, 'solve_once', solve_once)
    d = deepcopy(data)
    d['economy']['production']['local_working_parts'] = dict(share=0.5, labour_fraction_of_import_price=0.25)
    spec = json.loads((HERE / 'recovery.json').read_text())
    physical, loops = industry.build(d, spec, small_book(12), 'none')
    assert calls[0] == 0.0 and calls[1] == 1.0
    assert loops['local_working_parts']['staffing_passes'] >= 1
