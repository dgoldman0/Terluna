"""The coupled run: its saved product, and its reproduction of the staged appraisal."""
import gzip
import json
import math
import pytest
from research.studies.enclosed_returns import markets
from research.studies.vacuum_economy import programme as definition, run

SAVED = json.loads((run.HERE / 'results/economy.json').read_text())
CASES = {c['name']: c for c in SAVED['cases']}
STAGED = [('receipts_USD2041_year', 'annual_gross'), ('operating_cash_USD2041_year', 'annual_operating_available'),
          ('cumulative_receipts_USD2041', 'cumulative_gross_receipts'), ('payouts_USD2041', 'owner_payouts'),
          ('initial_funding_USD2041', 'initial_full_funding'), ('additional_equity_USD2041', 'additional_equity_required'),
          ('all_contributions_USD2041', 'all_funder_contributions'), ('real_IRR', 'real_IRR_roundoff_screened'),
          ('population', 'population'), ('first_payout_year', 'first_payout_year')]
DISCOUNTED = [('funder_cash_NPV_USD2041', 'sponsor_distributed_cash_NPV'), ('system_cash_NPV_USD2041', 'unlevered_cash_NPV'),
              ('socioeconomic_NPV_USD2041', 'full_price_socioeconomic_NPV')]


def test_the_saved_product_names_its_producer_its_inputs_and_every_listed_case():
    assert SAVED['schema'] == run.SCHEMA
    assert SAVED['producer']['sha256'] == run.digest(run.__file__)
    assert SAVED['dependencies'] == run.dependencies()
    assert list(CASES) == [c['name'] for c in run.cases()]
    assert all(c['status'] == 'conditional_account' and c['first_capacity_shortfall'] is None for c in SAVED['cases'])


def test_the_first_pass_cases_reproduce_the_staged_appraisal_exactly():
    staged = json.loads(gzip.decompress((run.RETURNS / 'results/revised.json.gz').read_bytes()))
    combined = {c['name']: c for c in staged['primary_cases'] if c['variant'] == 'combined'}
    compared = 0
    for case in SAVED['cases']:
        assert case['scenario'] == 'first_pass' and case['recovery'] == 'combined'
        if case['ramp'] != 'first_pass':
            continue
        old = {(c['mode'], c['horizon']): c for c in combined[':'.join((case['portfolio'], case['market'], case['growth']))]['checkpoints']}
        assert len(case['checkpoints']) == len(old) == 9
        for new in case['checkpoints']:
            reference = old[new['financing'], new['horizon']]
            for ours, theirs in STAGED:
                assert new[ours] == reference[theirs]
            for ours, theirs in DISCOUNTED:
                assert new[ours] == {rate: v[theirs] for rate, v in reference['discounted'].items()}
            compared += 1
    assert compared == 81


def test_a_case_solved_afresh_equals_its_saved_record():
    case = run.cases(['first_pass:services'])[0]
    assert json.loads(json.dumps(run.solve(case))) == CASES['first_pass:services']


@pytest.mark.parametrize('name', sorted(CASES))
def test_the_periods_add_up_to_the_checkpoints(name):
    case = CASES[name]
    periods = case['periods']
    assert periods[0]['first_year'] == 0 and all(a['last_year'] + 1 == b['first_year'] for a, b in zip(periods, periods[1:]))
    for horizon in (100, 200, 220):
        at = next(c for c in case['checkpoints'] if c['financing'] == 'U_credit' and c['horizon'] == horizon)
        through = [p for p in periods if p['last_year'] <= horizon]
        assert through[-1]['last_year'] == horizon
        assert math.isclose(math.fsum(p['receipts_USD2041'] for p in through), at['cumulative_receipts_USD2041'], rel_tol=1e-12)
        assert math.isclose(math.fsum(p['distributions_USD2041'] for p in through), at['payouts_USD2041'], rel_tol=1e-9, abs_tol=1.)
    for p in periods:
        assert math.isclose(math.fsum(p['receipts_by_sector_USD2041'].values()), p['receipts_USD2041'], rel_tol=1e-12, abs_tol=1.)
        assert math.isclose(math.fsum(p['costs_by_category_USD2041'].values()), p['costs_USD2041'], rel_tol=1e-12, abs_tol=1.)
        if p['first_year'] == 0:
            assert p['paid_to_Earth_share'] is None and p['receipts_USD2041'] == 0
        else:
            assert 0 <= p['paid_to_Earth_USD2041'] <= p['costs_USD2041']
            assert math.isclose(p['paid_to_Earth_share'], p['paid_to_Earth_USD2041'] / p['costs_USD2041'], rel_tol=1e-12)


@pytest.mark.parametrize('name', sorted(CASES))
def test_material_drawn_only_grows_and_financing_changes_no_receipt(name):
    case = CASES[name]
    drawn = [case['material_drawn'][str(h)] for h in (100, 200, 220)]
    for kind in ('asteroid_feed_kg', 'lunar_raw_kg', 'water_kg', 'nitrogen_kg', 'virgin_metal_kg'):
        assert 0 < drawn[0][kind] <= drawn[1][kind] <= drawn[2][kind]
    by = {(c['financing'], c['horizon']): c for c in case['checkpoints']}
    for (financing, horizon), c in by.items():
        equity = by['equity', horizon]
        assert c['receipts_USD2041_year'] == equity['receipts_USD2041_year']
        assert c['cumulative_receipts_USD2041'] == equity['cumulative_receipts_USD2041']
        assert c['socioeconomic_NPV_USD2041'] == equity['socioeconomic_NPV_USD2041']


def test_every_portfolio_runs_on_the_staged_schedule_and_on_both_programme_ramps():
    by = {}
    for case in SAVED['cases']:
        by.setdefault(case['portfolio'], set()).add(case['ramp'])
    assert len(by) == 9 and all(ramps == {'first_pass', 'constant_growth', 'front_loaded'} for ramps in by.values())


@pytest.mark.parametrize('ramp', ['constant_growth', 'front_loaded'])
def test_the_order_book_takes_its_residents_from_the_named_ramp_and_leaves_the_staged_module_alone(ramp):
    programme = run.sources()['programme']
    case = run.cases(['first_pass:services:' + ramp])[0]
    book = run.orders(case)
    for row in book['annual']:
        assert row['population'] == definition.population(programme, ramp, row['year'])
        assert row['new_residential_people'] == definition.population(programme, ramp, row['year'] + 2) - definition.population(programme, ramp, row['year'] + 1)
    rules = markets.inputs()[0]
    assert run.staged_demand(ramp) is not markets and run.staged_demand('first_pass') is markets
    assert markets.population(100, rules) == definition.population(programme, 'first_pass', 100)


def test_an_unknown_scenario_ramp_or_case_is_refused():
    case = dict(run.cases(['first_pass:services'])[0])
    with pytest.raises(ValueError):
        run.orders(dict(case, scenario='central'))
    with pytest.raises(ValueError):
        run.orders(dict(case, ramp='stepwise'))
    with pytest.raises(ValueError):
        run.cases(['first_pass:unlisted'])


def test_a_case_changes_coefficients_in_its_own_copy_and_a_misspelt_path_is_refused():
    shared = run.sources()['engineering']
    before = json.dumps(shared, sort_keys=True)
    case = run.cases(['first_pass:diversified:constant_growth:launch_100:earth_leg_split_to_20'])[0]
    changed = run.engineering(case)
    assert changed['industrial']['transport']['LEO_USD2026_kg'] == 100
    assert changed['industrial']['transport']['earth_leg'] == {'by_year': [[0, 'split'], [20, 'space']]}
    assert json.dumps(shared, sort_keys=True) == before and shared['industrial']['transport']['earth_leg'] == 'space'
    assert run.engineering(run.cases(['first_pass:diversified'])[0]) is shared
    with pytest.raises(ValueError):
        run.engineering(dict(case, engineering={'industrial.transport.earth_legs': 'split'}))


@pytest.mark.parametrize('name', sorted(CASES))
def test_traffic_counts_the_launches_and_the_earth_legs_propellant_by_arrangement(name):
    case = CASES[name]
    leg = case.get('engineering', {}).get('industrial.transport.earth_leg', 'space')
    for horizon in ('100', '200', '220'):
        at = case['traffic'][horizon]
        assert at['launches_year'] == pytest.approx(at['launched_kg_year'] / 100000.)
        assert (at['Earth_leg_propellant_launched_kg_year'] > 0) == (leg == 'split')
        assert at['Earth_leg_propellant_launched_kg_year'] < at['launched_kg_year']


@pytest.mark.parametrize('portfolio', ['diversified', 'services'])
@pytest.mark.parametrize('price', ['', ':launch_100', ':launch_50'])
def test_fuelling_the_first_years_from_earth_lowers_the_opening_funding_and_ends_on_the_space_arrangement(portfolio, price):
    base = 'first_pass:%s:constant_growth%s' % (portfolio, price)
    space, switched, split = CASES[base], CASES[base + ':earth_leg_split_to_20'], CASES[base + ':earth_leg_split']
    funding = lambda c: next(x['initial_funding_USD2041'] for x in c['checkpoints'] if x['financing'] == 'U_credit' and x['horizon'] == 200)
    assert funding(switched) == funding(split) < 0.6 * funding(space)
    assert switched['traffic']['200']['launched_kg_year'] == pytest.approx(space['traffic']['200']['launched_kg_year'], rel=1e-3)
    assert split['traffic']['200']['launched_kg_year'] > 1.5 * space['traffic']['200']['launched_kg_year']
    assert split['material_drawn']['200']['asteroid_feed_kg'] < switched['material_drawn']['200']['asteroid_feed_kg'] <= space['material_drawn']['200']['asteroid_feed_kg']


@pytest.mark.parametrize('name', sorted(CASES))
def test_the_workforce_record_shows_working_parts_made_only_where_a_case_asks_for_them(name):
    case = CASES[name]
    asked = 'economy.production.local_working_parts.share' in case.get('engineering', {})
    for horizon in ('100', '200', '220'):
        at = case['workforce'][horizon]
        assert 0 < at['share_of_workforce_in_use'] <= 1 + 1e-9 and at['jobs'] > at['working_parts_jobs'] >= 0
        assert (at['working_parts_jobs'] > 0) == asked == (at['working_parts_made_share'] > 0)
        assert at['working_parts_made_share'] <= 0.78 + 1e-12


@pytest.mark.parametrize('portfolio', ['diversified', 'services'])
def test_making_working_parts_in_place_lowers_what_is_paid_to_earth_and_launched_and_cheaper_labour_makes_more(portfolio):
    base = CASES['first_pass:%s:constant_growth' % portfolio]
    made = [CASES['first_pass:%s:constant_growth:local_parts%s' % (portfolio, s)] for s in ('_labour_50', '', '_labour_10')]
    paid = lambda c: next(p['paid_to_Earth_share'] for p in c['periods'] if p['first_year'] == 151)
    value = lambda c: next(x['funder_cash_NPV_USD2041']['0.03'] for x in c['checkpoints'] if x['financing'] == 'U_credit' and x['horizon'] == 200)
    series = [base, *made]
    assert all(paid(a) > paid(b) for a, b in zip(series, series[1:]))
    assert all(value(a) < value(b) for a, b in zip(series, series[1:]))
    assert all(a['traffic']['200']['launches_year'] >= b['traffic']['200']['launches_year'] for a, b in zip(series, series[1:]))
    assert all(a['workforce']['100']['working_parts_made_share'] < b['workforce']['100']['working_parts_made_share'] for a, b in zip(made, made[1:]))
