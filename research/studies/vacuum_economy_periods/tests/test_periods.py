"""The vacuum economy period by period: its technology paths, the settlements from their designs, the engineering chosen
by cost, the account with Earth, residents as a result and the saved product."""
import json
import math
from pathlib import Path

import pytest

from research.studies.vacuum_economy_periods import designs, economy, engineering, money, run, sales, technology as tech

HERE = Path(__file__).resolve().parents[1]


def test_a_path_holds_its_stated_values_interpolates_and_the_static_case_holds_2041():
    t = tech.roadmap()
    years = t['years']
    for name, q in t['quantities'].items():
        for case in ('slow', 'central', 'fast'):
            for year, value in zip(years, q[case]):
                assert tech.value(name, case, year) == pytest.approx(value)
        assert all(tech.value(name, 'static', y) == pytest.approx(q['slow'][0]) for y in (2041, 2100, 2241))
    q = t['quantities']['array_W_kg']
    mid = (years[1] + years[2]) / 2
    assert tech.value('array_W_kg', 'central', mid) == pytest.approx(math.sqrt(q['central'][1] * q['central'][2]))
    for name in t['arrivals']:
        assert tech.arrival(name, 'static') is None


def test_with_2041_technology_the_founding_system_returns_its_design():
    """Four cities of 50,000 and a hundred towns of 500, at the design's own engineering, take what the design's
    production schedule says a resident takes, less the power hardware the model sizes itself."""
    city = designs.build('city', 'static', 2041, True)
    town = designs.build('town', 'static', 2041, True)
    system = json.loads(designs.CITIES.read_text())['cases']['reference']
    d = system['schedule']['production_demand']
    n = system['population']
    mix = {k: 0.8 * city[k] + 0.2 * town[k] for k in designs.STREAMS}
    assert mix['handling_kg'] == pytest.approx(d['handling'] / n, rel=1e-9)
    assert mix['bulk_kg'] == pytest.approx(d['local_bulk'] / n, rel=1e-9)
    assert mix['metal_kg'] == pytest.approx(d['local_metal'] / n, rel=1e-9)
    assert mix['pressure_kg'] == pytest.approx(d['pressure_components'] / n, rel=1e-9)
    assert mix['equipment_kg'] < d['other_imports'] / n
    need = designs.needs('city', 'static', 2041)
    assert need['end_use_W']['food'] == pytest.approx(system['utilities']['mean_load_components_W']['food'] / n, rel=1e-9)


@pytest.mark.parametrize('case', ['slow', 'central', 'fast'])
def test_improving_technology_lowers_what_a_resident_needs(case):
    years = (2041, 2070, 2100, 2150, 2241)
    food = [designs.food_W(case, y) for y in years]
    assert all(b < a for a, b in zip(food, food[1:]))
    built = [designs.build('town', case, y, False)['metal_kg'] for y in years]
    assert all(b <= a for a, b in zip(built, built[1:]))


def test_the_power_supply_is_the_cheapest_option_offered_and_the_polar_grid_waits_for_its_arrival():
    for case, year in (('central', 2060), ('central', 2120), ('fast', 2100), ('static', 2200)):
        p = engineering.Prices(case, year, 155797.0)
        s = p.firm_supply(site_km=1600)
        assert s['chosen'] == min(s['options'], key=lambda k: s['options'][k]['USD_kW_year'])
        assert ('polar_grid' in s['options']) == tech.arrived('polar_grid', case, year)
        assert ('fusion' in s['options']) == tech.arrived('fusion', case, year)


def test_a_chemical_tug_needs_more_water_the_farther_it_flies_and_stops_closing():
    p = engineering.Prices('central', 2049, 155797.0)
    drawn = [engineering.tug('chemical', leg, p)['drawn_kg'] for leg in (1000, 2000, 3000, 4000)]
    assert all(b > a for a, b in zip(drawn, drawn[1:]))
    assert engineering.tug('chemical', 9000, p) is None


def test_no_tug_is_flown_past_its_practical_limit():
    limit = json.loads((HERE / 'period_inputs.json').read_text())['transport']['tug_max_dry_kg_per_kg_year']['value']
    p = engineering.Prices('fast', 2241, 155797.0)
    flown = 0
    for kind in ('chemical', 'electric', 'nuclear_electric'):
        for leg in range(1000, 20001, 500):
            x = engineering.tug(kind, leg, p)
            if x is not None:
                flown += 1
                assert x['dry_kg_per_kg_year'] <= limit
    assert flown > 0
    assert engineering.tug('nuclear_electric', 7000, p) is None and engineering.tug('electric', 7000, p) is not None


@pytest.mark.parametrize('scenario,case,ramp', [('central', 'central', 'constant_growth'), ('low', 'static', 'front_loaded'),
                                                ('high', 'fast', 'constant_growth')])
def test_the_account_adds_up_and_the_bounds_hold(scenario, case, ramp):
    rows = economy.Run(scenario, case, ramp, u_policy='dollar').simulate()
    total = 0.0
    bound = json.loads((HERE / 'period_inputs.json').read_text())['sales']['visitor_cap_share'][case]
    for r in rows:
        assert r['net'] == pytest.approx(r['receipts'] - r['paid_to_Earth'])
        assert r['paid_to_Earth'] == pytest.approx(sum(r['to_Earth'].values()))
        total += r['net']
        assert r['cumulative'] == pytest.approx(total)
        assert r['visitor_share_of_world'] <= bound
        assert 0.0 <= r['lines_served_share'] <= 1.0
    s = economy.summary(rows)
    assert s['funding_needed_USD2041'] == pytest.approx(-min(0.0, min(r['cumulative'] for r in rows)))
    assert s['funding_needed_USD2041'] >= s['opening_funding_USD2041']


def test_better_technology_needs_less_funding_for_the_plan():
    funding = [economy.summary(economy.Run('central', c, 'constant_growth').simulate())['funding_needed_by_year_USD2041']['200']
               for c in ('static', 'slow', 'central', 'fast')]
    assert funding[0] > funding[1] > funding[3]


def test_residents_as_a_result_take_the_highest_rates_the_test_allows_to_the_plan_and_past_it():
    dollar = (('u_policy', 'dollar'),)
    rates, sched = run.solve('central', 'central', 'steady', dollar)
    assert rates['sustained'] and rates['plan_year'] is not None and rates['past_plan'] > 0
    assert sched == run.steady_schedule(rates['to_plan'], rates['past_plan'])
    assert run.sustained('central', 'central', sched, dollar)[0]
    assert not run.sustained('central', 'central', run.steady_schedule(rates['to_plan'] + 0.002, 0.0), dollar)[0]
    assert not run.sustained('central', 'central', run.steady_schedule(rates['to_plan'], rates['past_plan'] + 0.002), dollar)[0]
    rows = economy.Run('central', 'central', schedule=sched, u_policy='dollar').simulate()
    s = economy.summary(rows)
    assert s['NPV_USD2041']['0.03'][str(run.HORIZON)] >= 0 and s['deficit_decades_after_repayment'] == 0
    assert rows[run.HORIZON]['residents'] > run.planned()


def test_the_test_binds_year_by_year_where_the_present_value_alone_would_not():
    """Growth carried on at the programme's pace passes a test on the present value at 3%, since discounting shrinks
    its late deficits, and fails the test that once repaid the account pays its way decade by decade."""
    rows = economy.Run('central', 'central', 'constant_growth', u_policy='dollar').simulate()
    s = economy.summary(rows)
    assert s['NPV_USD2041']['0.03'][str(run.HORIZON)] > 0
    assert s['deficit_decades_after_repayment'] > 0
    assert not run.sustained('central', 'central', [rows[y]['residents'] for y in range(run.HORIZON + 1)], (('u_policy', 'dollar'),))[0]


def test_a_ramp_carries_on_past_the_plan_at_its_final_rate_until_the_sites_are_full():
    r = economy.Run('central', 'central', 'constant_growth', last_year=run.HORIZON + 20)
    rows = r.simulate()
    n = [x['residents'] for x in rows]
    assert n[200] == pytest.approx(run.planned())
    assert n[201] / n[200] == pytest.approx(n[200] / n[199], rel=1e-9)
    assert max(n) <= r.capacity() * (1 + 1e-12)
    assert economy.summary(rows)['sites_full_year'] is not None
    assert all(x['housing_short'] == 0 for x in rows)
    held = economy.Run('central', 'central', 'constant_growth', past_plan='hold').simulate()
    assert held[run.HORIZON]['residents'] == pytest.approx(run.planned())
    assert [x['net'] for x in held[:200]] == pytest.approx([x['net'] for x in rows[:200]])


def test_newcomers_past_the_initial_quota_live_on_the_eventual_set_of_sites():
    r = economy.Run('central', 'central', 'constant_growth')
    rows = r.simulate()
    inputs = json.loads((HERE / 'period_inputs.json').read_text())['settlements']
    assert r.capacity() == pytest.approx(4 * 50000 + inputs['eventual_crater_sites']['value'] * 750 + 50000)
    assert rows[200]['housed_beyond_initial_quota'] == 0
    assert rows[run.HORIZON]['housed_beyond_initial_quota'] > 0
    assert rows[run.HORIZON]['housed_town'] <= inputs['eventual_crater_sites']['value'] * 750 + 50000 + 1


def test_water_set_aside_lengthens_the_legs_and_brings_exhaustion_sooner():
    assert engineering.usable_water_kg('p50', 0.9) == pytest.approx(0.1 * engineering.usable_water_kg('p50'))
    for drawn in (1e12, 1e13, 1e14):
        assert engineering.leg_for(drawn, 'p50', 0.9) >= engineering.leg_for(drawn, 'p50')
    whole = run.exhaustion('central', 'central', 'ramp', (('reach', 'p50'), ('water_set_aside', 0.0)))
    kept = run.exhaustion('central', 'central', 'ramp', (('reach', 'p50'), ('water_set_aside', 0.9)))
    assert kept['water_exhausted_year'] < whole['water_exhausted_year']
    for e in (whole, kept):
        assert e['last_year'] == e['water_exhausted_year']
        if e['water_out_because'] == 'all drawn':
            assert e['contained_drawn_kg_at_end'] >= e['usable_water_kg']
        else:
            assert e['water_out_because'] == 'the rest beyond every tug' and e['contained_drawn_kg_at_end'] < e['usable_water_kg']
    assert whole['residents_when_sites_full'] == pytest.approx(whole['capacity_residents'])


def test_the_U_ledger_adds_up_and_holds_its_rules():
    """With a policy of U the reserve changes by the year's dollars and the metals kept from the zone's refining: the
    account with Earth at the year's closing rate, Earth's purchases of U net of sales, the funders' conversions, the
    reserve's interest and the dollars paid to the funders. U stays at par until the account is repaid and between its
    floor and par after; working balances are a week of the trade settled in U; the funders' measures come from their
    own flows; the conservative policy pays no yield and nothing in U to the funders."""
    rows = economy.Run('central', 'central', 'constant_growth', u_policy='liberal', u_floor=0.5, past_plan='hold').simulate()
    reserve = 0.0
    for r in rows:
        change = (r['receipts'] - r['paid_to_Earth'] + r['u_from_balances'] + r['u_from_investors'] + r['u_from_holders'] + r['u_interest']
                  - r['u_funders_converted'] - r['u_dollars_to_funders'] + r['u_metals_retained'])
        assert r['u_reserve'] - reserve == pytest.approx(change, abs=1e-6 * max(1.0, r['receipts']))
        reserve = r['u_reserve']
        assert r['u_reserve_currency'] >= 0 and r['u_reserve_metals'] >= 0
        assert 0.5 - 1e-12 <= r['u_rate'] <= 1.0 + 1e-12
        if not r['u_repaid']:
            assert r['u_rate'] == 1.0
        if r['u_rate'] == r['u_rate_opening']:
            assert r['u_balances'] * r['u_rate'] == pytest.approx(r['u_settlement_share'] * (r['receipts'] + r['paid_to_Earth']) * 7 / 365.25)
    s = economy.summary(rows)
    assert s['funding_needed_USD2041'] == pytest.approx(-min(0.0, min(r['funders_cumulative'] for r in rows)))
    plain = economy.Run('central', 'central', 'constant_growth', u_policy='conservative', past_plan='hold').simulate()
    assert all(r['u_yield_income'] == 0 and r['u_funders_U'] == 0 for r in plain)


def test_the_dollar_account_is_untouched_by_the_settings_of_U():
    """Without a policy of U the floor, the trade effect, the spending share and the credit line change nothing."""
    a = economy.Run('central', 'central', 'constant_growth', past_plan='hold', last_year=120, u_policy='dollar').simulate()
    b = economy.Run('central', 'central', 'constant_growth', past_plan='hold', last_year=120, u_policy='dollar', u_floor=0.5, u_trade_effect=-0.1,
                    u_spending_share=0.6, u_credit_line=True).simulate()
    assert [r['net'] for r in a] == [r['net'] for r in b]
    assert not any(k.startswith('u_') for k in a[50]) and 'zone_goods' not in a[50]['to_Earth']


def test_the_rate_absorbs_a_shortfall_down_to_its_floor_and_the_reserve_test_can_fail():
    """Once repaid, a lasting shortfall lowers U so that the residents' purchases from Earth close it, and the reserve
    holds its minimum while U is above its floor; held at par, the same shortfall takes the reserve below its minimum
    and fails the test; with the funders' credit line the reserve holds and the funders' flows carry the loan."""
    def ledger_run(floor, credit_line=False):
        led = money.Ledger(money.policy('liberal', dict(funders_paid_in_U=0.0)), money.holders(), floor=floor, credit_line=credit_line)
        out = []
        for y in range(0, 80):
            receipts = 0.0 if y < 8 else (2e11 if y < 40 else 1.15e11)
            out.append(led.step(y, receipts, 1e11, 2e10 if y >= 8 else 0.0, 2e11 if y >= 8 else 0.0, receipts))
        return out
    floating, par, lent = ledger_run(0.5), ledger_run(1.0), ledger_run(1.0, credit_line=True)
    late = range(45, 80)
    assert all(floating[y]['u_repaid'] for y in late)
    assert all(floating[y]['u_tests']['reserve'] for y in late)
    assert 0.5 < min(floating[y]['u_rate'] for y in late) < 0.99
    assert all(par[y]['u_rate'] == 1.0 and not par[y]['u_tests']['reserve'] for y in late)
    assert all(lent[y]['u_tests']['reserve'] for y in late) and sum(lent[y]['funders_flow'] for y in late) < 0



def test_investors_leave_no_faster_than_their_bonds_mature_and_their_stakes_unlock():
    """Investors drawn in by a high early yield, which then falls back to the zone's growth, sell for dollars only what
    matures, unlocks or is paid them in a year; what could leave within a year is the working balances, the funders'
    liquid U and the stakes unlocking, and the run test asks a quarter of it with the bonds falling due and six months of
    payments."""
    led = money.Ledger(money.policy('liberal', dict(real_yield=0.08, funders_paid_in_U=0.0, investors_staking_share=1 / 3)), money.holders())
    before, bonds = None, {}
    for y in range(0, 70):
        r = led.step(y, 2e11 if y >= 8 else 0.0, 1e11, 2e10 if y >= 8 else 0.0, 2e11 * 1.06 ** max(0, min(y, 30) - 8) if y >= 8 else 0.0,
                     1e11 if y >= 8 else 0.0)
        if before is not None and r['u_from_investors'] < 0:
            room = (before['u_bonds'] / 10 + before['u_staked'] / 2 + r['u_real_yield'] * (before['u_bonds'] + before['u_staked'])) * r['u_rate_opening']
            assert -r['u_from_investors'] <= room * (1 + 1e-9)
        if y >= 8:
            assert r['u_runnable'] == pytest.approx(r['u_balances'] + r['u_funders_U'] + r['u_holders_U']
                                                    + (r['u_staked'] + r['u_funders_staked'] + r['u_holders_staked']) / 2)
            assert r['u_falling_due'] == pytest.approx(r['u_bonds'] / 10)
        before = r
        bonds[y] = r['u_bonds']
    assert bonds[28] > 0 and bonds[38] > 0.2 * bonds[28]



def test_holders_for_use_keep_months_of_their_spending_in_U_and_stake_where_a_yield_is_paid():
    """People on Earth who spend in the zone on its culture and visits keep the set months of that spending, settled in U,
    a third staked while the liberal policy pays its yield, held within the zone's growth, whose reward they spend in the
    zone; under the conservative policy, with no yield, they stake nothing."""
    months = money.holders()['holders_for_use_months']
    for name, staked_share in (('liberal', 1 / 3), ('conservative', 0.0)):
        led = money.Ledger(money.policy(name, dict(funders_paid_in_U=0.0)), money.holders())
        for y in range(0, 40):
            r = led.step(y, 2e11 if y >= 8 else 0.0, 1e11, 1e10 if y >= 8 else 0.0, 3e11 * 1.03 ** y if y >= 8 else 0.0, 1e11 if y >= 8 else 0.0,
                         consumer_receipts=1e11 * 1.03 ** y if y >= 8 else 0.0)
        held = r['u_holders_U'] + r['u_holders_staked']
        assert held * r['u_rate'] == pytest.approx(r['u_settlement_share'] * 1e11 * 1.03 ** 39 * months / 12, rel=1e-6)
        assert r['u_holders_staked'] == pytest.approx(staked_share * held, rel=1e-6)

def test_the_dollar_account_pays_for_the_dollars_its_residents_hold():
    """Without U the residents and firms hold dollars in proportion to the zone's output, and what the account pays Earth
    for them adds up to that stock."""
    rows = economy.Run('central', 'central', 'constant_growth', past_plan='hold', last_year=150, u_policy='dollar').simulate()
    share = json.loads((HERE / 'period_inputs.json').read_text())['money']['dollars_held_share_of_output']['value']
    r = rows[-1]
    output = r['workforce'] * r['wage'] + r['receipts'] - r['to_Earth']['lines']
    assert sum(x['to_Earth']['dollars_held'] for x in rows) == pytest.approx(share * output)
    plain = economy.Run('central', 'central', 'constant_growth', past_plan='hold', last_year=150, u_policy='dollar', dollars_held_share=0.0).simulate()
    assert all('dollars_held' not in x['to_Earth'] for x in plain)

def test_the_five_tests_catch_an_unsafe_policy():
    """The funders' whole returns paid in U by a zone whose output is small beside its sales pile Earth's holdings above
    the zone's yearly output; a yield held above the zone's growth fails its test."""
    heavy = money.Ledger(money.policy('liberal', dict(funders_paid_in_U=1.0)), money.holders())
    for y in range(0, 40):
        a = heavy.step(y, 1e12 if y >= 8 else 0.0, 1e11, 0.0, 1e10 if y >= 8 else 0.0, 1e10 if y >= 8 else 0.0)
    assert a['u_holdings_USD'] > a['u_output'] and not a['u_tests']['holdings_within_output']
    steep = money.Ledger(money.policy('liberal', dict(real_yield=0.08, yield_follows_growth=False, funders_paid_in_U=0.0)), money.holders())
    for y in range(0, 60):
        b = steep.step(y, 1e11 if y >= 8 else 0.0, 1e11, 0.0, 1e12 if y >= 8 else 0.0, 1e11 if y >= 8 else 0.0)
    assert b['u_investment'] > 0 and not b['u_tests']['yield_within_growth']


def test_metals_kept_for_the_reserve_raise_the_price_of_those_sold():
    p = engineering.Prices('central', 2070, 155797.0)
    plain = sales.year('central', 'central', 29, p, 2000.0, 5000.0, 5e10)['metals']
    kept = sales.year('central', 'central', 29, p, 2000.0, 5000.0, 5e10, metals_keep_USD=0.3 * plain['receipts'])['metals']
    assert plain['receipts'] > 0 and plain['retained_value'] == 0
    assert kept['receipts'] < plain['receipts'] < kept['receipts'] + kept['retained_value']


def test_goods_bought_with_U_take_the_lines_earth_inputs_and_hands():
    rows = economy.Run('central', 'central', 'constant_growth', u_policy='liberal', past_plan='hold', last_year=150).simulate()
    r = next(r for r in rows[60:] if rows[r['year'] - 1]['u_spent_in_zone'] > 0)
    assert r['to_Earth']['zone_goods'] > 0 and r['jobs']['zone_goods'] > 0
    assert r['to_Earth']['zone_goods'] / rows[r['year'] - 1]['u_spent_in_zone'] == pytest.approx(r['to_Earth']['lines'] / r['receipts'])


def test_the_solver_takes_the_highest_rate_below_which_every_checked_rate_passes():
    assert run.highest(lambda r: r < 0.123) == pytest.approx(0.123, abs=1e-4)
    assert run.highest(lambda r: r < 0.05 or 0.06 < r < 0.2) == pytest.approx(0.05, abs=1e-4)


def test_growth_as_a_result_under_U_is_no_lower_with_a_lower_floor():
    found = []
    for floor in (1.0, 0.5):
        opts = (('u_floor', floor), ('u_policy', 'liberal'))
        rates, sched = run.solve('central', 'central', 'steady', opts)
        assert rates['sustained']
        for k in (5, 10, 15):
            assert run.sustained('central', 'central', run.steady_schedule(rates['to_plan'], rates['past_plan'] * k / 20), opts)[0]
        found.append(sched[run.HORIZON])
    assert found[1] >= found[0]


def test_the_saved_product_is_current():
    saved = json.loads((HERE / 'results/periods.json').read_text())
    assert saved['schema'] == run.SCHEMA
    assert saved['producer']['sha256'] == run.digest(HERE / 'run.py')
    assert saved['dependencies'] == run.dependencies()
    names = [c['name'] for c in saved['cases']]
    assert names[:len(run.cases())] == [c[0] for c in run.cases()]
    assert names[len(run.cases()):] in ([], ['central:central:steady:U_most_liberal_safe'])
    assert len(saved['u_policy_search']) == len(run.SEARCH_YIELDS) * len(run.SEARCH_FUNDERS_IN_U)


def test_the_tug_fleet_turns_over_gradually_and_always_carries_the_water():
    """On the programme's ramp, the fleet can always carry the year's water. When the kind the factories build first
    changes from chemical, the new kinds hold under a quarter of the fleet, and five years on the chemical tugs still
    carry most of it. Tugs are sold used only in years some are retired early or can no longer fly the leg. The year's
    receipts are the lines' receipts and those sales."""
    rows = economy.Run('central', 'central', 'constant_growth').simulate()
    for r in rows:
        if not r['exhausted'] and r['water_delivered_t'] > 0:
            assert sum(r['fleet_by_kind_t_year'].values()) >= r['water_delivered_t'] * (1 - 1e-9)
        if r['used_tugs_USD'] > 0:
            assert r['tugs_retired_for_sale_t'] > 0
        sold = sum(x['receipts'] for x in r['lines'].values() if x['served'])
        assert r['receipts'] == pytest.approx(sold + r['used_tugs_USD'])
    first = next(r for r in rows if r['date'] > 2060 and r['tug'] not in (None, 'chemical'))
    share = lambda r, k: r['fleet_by_kind_t_year'].get(k, 0.0) / sum(r['fleet_by_kind_t_year'].values())
    assert 1 - share(first, 'chemical') < 0.25
    assert share(rows[first['year'] + 5], 'chemical') > 0.5


def test_a_tug_flies_on_beyond_the_limit_on_building_until_its_leg_no_longer_closes():
    """The practical limit on a tug's mass bars building it for a leg; a tug already built still flies that leg, heavier
    per kilogram delivered, as long as the leg closes."""
    P = economy.prices('central', 2200, 150000.0)
    legs = [leg for leg in range(3000, 9001, 50) if engineering.tug('chemical', leg, P) is None]
    assert legs, 'the chemical tug reaches its limit on building somewhere below 9 km/s'
    beyond = legs[0]
    flown = engineering.tug('chemical', beyond, P, limit=False)
    built = engineering.tug('chemical', beyond - 50, P)
    assert flown is not None and built is not None
    assert flown['dry_kg_per_kg_year'] > built['dry_kg_per_kg_year']


def test_a_case_can_move_one_group_of_the_stack_while_the_rest_keep_theirs():
    """Every quantity and arrival belongs to exactly one group; a moved group follows its own case and every other
    quantity and arrival the base case; an unknown group or case is refused."""
    t = tech.roadmap()
    grouped = [name for names in tech.GROUPS.values() for name in names]
    assert sorted(grouped) == sorted(set(t['quantities']) | set(t['arrivals']))
    for group, names in tech.GROUPS.items():
        case = 'central|%s=slow' % group
        for name in t['quantities']:
            expected = 'slow' if name in names else 'central'
            assert tech.value(name, case, 2120) == tech.value(name, expected, 2120)
        for name in t['arrivals']:
            assert tech.arrival(name, case) == tech.arrival(name, 'slow' if name in names else 'central')
    assert tech.base('central|power=fast') == 'central' and tech.case_for('fusion', 'central|power=fast') == 'fast'
    with pytest.raises(ValueError):
        tech.moves('central|warp=fast')
    with pytest.raises(ValueError):
        tech.moves('central|power=quick')


def test_computing_and_visits_options_replace_only_the_scenarios_values():
    """The scenario's own computing at the opening and price response of visits leave a case as it was; a larger
    contract sells more computing and a stronger price response changes what visits earn."""
    base = economy.Run('central', 'central', 'constant_growth').simulate()
    s = sales.scenarios()['scenarios']['central']
    same = economy.Run('central', 'central', 'constant_growth', computing_W_2049=s['computing']['capacity_W_2049']['value'],
                       visits_price_elasticity=s['visits']['price_elasticity']['value']).simulate()
    assert [r['net'] for r in same] == [r['net'] for r in base]
    more = economy.Run('central', 'central', 'constant_growth', computing_W_2049=2 * s['computing']['capacity_W_2049']['value']).simulate()
    assert more[20]['lines']['computing']['receipts'] == pytest.approx(2 * base[20]['lines']['computing']['receipts'])
    other = economy.Run('central', 'central', 'constant_growth', visits_price_elasticity=2.0).simulate()
    assert other[20]['lines']['visits']['receipts'] != pytest.approx(base[20]['lines']['visits']['receipts'])


def test_the_funders_rate_of_return_is_recovered_from_known_flows():
    """A loan repaid with a tenth more a year later returns 10%; flows that lose money give a rate below zero; flows whose
    present value never changes sign give none."""
    from research.studies.vacuum_economy_periods import scoping
    assert scoping.irr([-100.0, 110.0]) == pytest.approx(0.10, abs=1e-9)
    assert scoping.irr([-100.0, 0.0, 121.0]) == pytest.approx(0.10, abs=1e-9)
    assert scoping.irr([-100.0, 95.0]) == pytest.approx(-0.05, abs=1e-9)
    assert scoping.irr([100.0, 10.0]) is None


def test_the_scoping_product_is_current():
    """The scoping study's measures were drawn from the present model and its saved product."""
    from research.studies.vacuum_economy_periods import scoping
    saved = json.loads((HERE / 'results/scoping.json').read_text())
    assert saved['schema'] == scoping.SCHEMA
    assert saved['producer']['sha256'] == scoping.digest(HERE / 'scoping.py')
    product = 'research/studies/vacuum_economy_periods/results/periods.json'
    assert saved['dependencies'] == dict(run.dependencies(), **{product: scoping.digest(HERE / 'results/periods.json')})
    assert [s['input'] for s in saved['sensitivities']] and len(saved['sensitivities']) == len(scoping.SENSITIVITIES)
    assert sorted(g['group'] for g in saved['technology_groups']) == sorted(tech.GROUPS)


def test_the_returns_by_growth_path_come_from_their_paths():
    """The funders' return for each growth path is that of the path's own flows: the base case's growth to the plan held
    there returns what a fresh run of that schedule returns, the two held paths end at the plan, and the base case's entry
    is the base case's own measure."""
    from research.studies.vacuum_economy_periods import scoping
    saved = json.loads((HERE / 'results/scoping.json').read_text())
    paths = saved['returns_by_path']
    base = saved['paths']['base']
    rows = economy.Run('central', 'central', schedule=run.steady_schedule(base['growth_rates']['to_plan'], 0.0)).simulate()
    held = scoping.irr([r.get('funders_flow', r['net']) for r in rows])
    assert paths['fastest_to_plan_held']['funders_IRR_to_2341'] == pytest.approx(held, abs=1e-9)
    for name in ('programme_pace_held', 'fastest_to_plan_held'):
        assert paths[name]['residents_2341'] == pytest.approx(run.planned())
    assert paths['base']['funders_IRR_to_2341'] == base['measures']['funders_IRR_to_2341']
    assert paths['base']['residents_2341'] > run.planned()


def test_the_funders_outcome_in_every_case_comes_from_the_case_and_the_account_adds_up_by_period():
    """Each grid cell's repayment year and rate of return are those of its own solved schedule; the base case's flows by
    period add up to its whole run; the funding to the end of the opening year is what it pays Earth less its receipts,
    with the reserve U holds by then; and the trade effects grow worse with a larger fall in Earth's purchases."""
    from research.studies.vacuum_economy_periods import scoping
    saved = json.loads((HERE / 'results/scoping.json').read_text())
    cases = {c['name']: c for c in json.loads((HERE / 'results/periods.json').read_text())['cases']}
    cell, case = saved['grid']['low:central:U'], cases['low:central:steady']
    fresh = scoping.outcome('low', 'central', scoping.saved_schedule(case), case['options'])
    assert cell['repaid'] == fresh['repaid']
    assert cell['funders_IRR_to_2341'] == pytest.approx(fresh['funders_IRR_to_2341'], abs=1e-9)
    base = saved['paths']['base']['series']
    assert sum(p['funders_flow'] for p in saved['base_periods']) == pytest.approx(base['funders_cumulative'][-1], rel=1e-9)
    f = saved['founding_capital']
    assert f['funding'] == pytest.approx(sum(f['paid_to_Earth_by_kind'].values()) - f['receipts'] + f['reserve_at_end'], rel=0.01)
    residents = [t['residents_2341'] for t in saved['trade_effect']]
    assert residents == sorted(residents, reverse=True)
    assert saved['USD2041_in_USD2026'] == pytest.approx(1 / 1.557967416600765)
