"""Pricing, dropping and staffing the lines of sale, checked on accounts written by hand."""
import math
import pytest
from provisioning.vacuum_industry import model
from research.studies.vacuum_economy import demand, mix, run


@pytest.fixture(scope='module')
def s():
    return run.sources()


def decision(s, name='central'):
    return mix.Decision(s['programme'], 'constant_growth', s['scenarios'], name, s['engineering'], s['recovery'], 'combined',
                        s['periods'], hotel=62500.0)


def account(d, **lines):
    """An account in which each named line has the same receipts, cost, quantity and jobs in every year."""
    return {line: dict(years={y: dict(receipts=r, cost=c, quantity=q, jobs=j) for y in d.years}, founding=f, unit='kg')
            for line, (r, c, q, j, f) in lines.items()}


def test_a_summary_gives_each_period_its_cost_per_unit_and_charges_the_founding_once(s):
    d = decision(s)
    summary = d.summary(account(d, propellant=(30.0, 12.0, 4.0, 2.0, 1000.0)))['propellant']
    first, second = summary['periods'][0], summary['periods'][1]
    assert (first['first_year'], first['last_year']) == (8, 25) and first['quantity'] == pytest.approx(18 * 4.0)
    assert first['cost_USD2041'] == pytest.approx(18 * 12.0 + 1000.0)
    assert second['unit_cost_USD2041'] == pytest.approx(3.0) and second['mean_jobs'] == pytest.approx(2.0)
    discounted = math.fsum(18.0 * 1.03 ** -y for y in d.years) - 1000.0 * 1.03 ** -4
    assert summary['worth_USD2041'] == pytest.approx(discounted)
    assert summary['worth_per_job_USD2041'] == pytest.approx(discounted / math.fsum(2.0 * 1.03 ** -y for y in d.years))


def test_the_system_prices_a_line_from_its_cost_only_where_buyers_answer_strongly(s):
    costs = dict(propellant=(30.0, 12.0, 4.0, 2.0, 0.0), materials=(30.0, 40.0, 4.0, 2.0, 0.0), visits=(3e6, 1e6, 2.0, 1.0, 0.0))
    high = decision(s, 'high')
    prices = high.own_prices(high.summary(account(high, **costs)))
    assert prices['propellant'][100] == pytest.approx(3.0 * 1.4 / 0.4) and prices['materials'][100] == pytest.approx(10.0 * 1.4 / 0.4)
    assert prices['visits'][100] == pytest.approx(2 * 5e5)
    central = decision(s)
    prices = central.own_prices(central.summary(account(central, **costs)))
    assert 'propellant' not in prices and 'materials' not in prices          # an elasticity of one: the ceiling holds
    assert prices['visits'][100] == pytest.approx(5e5 * 1.25 / 0.25)


def test_hands_go_first_to_the_line_whose_jobs_bring_most_and_a_line_that_loses_is_dropped(s):
    d = decision(s)
    book = d.orders()
    rows = {row['year']: row for row in book['annual']}
    first = account(d, culture=(100.0, 10.0, 1.0, 50.0, 0.0), institutions=(100.0, 10.0, 1.0, 5000.0, 0.0), computing=(10.0, 30.0, 1.0, 10.0, 0.0))
    for line, record in first.items():
        for y in d.years:
            record['years'][y]['quantity'] = rows[y]['lines'][line]['quantity']
    base_jobs = {y: 0.0 for y in d.years}
    scale, order, dropped = d.allocate(d.summary(first), first, book, book, base_jobs)
    assert order == ['culture', 'institutions'] and dropped == ['computing']
    assert all(scale['computing'][y] == 0 for y in d.years)
    # At the opening 12,500 residents give 4,850 hands: culture takes its 50 and the hosted work is cut to the rest.
    room = mix.MARGIN * model.WORKFORCE_SHARE * rows[8]['population']
    assert scale['culture'][8] == 1.0 and scale['institutions'][8] == pytest.approx((room - 50.0) / 5000.0)
    assert scale['institutions'][200] == 1.0
    assert all(scale[line][y] == 0 for line in demand.LINES if line not in ('culture', 'institutions') for y in d.years)


def test_a_priced_line_is_judged_at_the_systems_price(s):
    d = decision(s, 'high')
    first = account(d, visits=(1.0, 2.0, 1.0, 1.0, 0.0))                     # a loss at the first price
    summary = d.summary(first)
    assert summary['visits']['worth_USD2041'] < 0
    priced = d.orders(d.own_prices(summary))
    for y in d.years:
        first['visits']['years'][y]['quantity'] = 1.0
    scale, order, dropped = d.allocate(summary, first, priced, priced, {y: 0.0 for y in d.years})
    assert order == ['visits'] and not dropped                              # sold at twice its cost it pays


def test_a_held_staffing_reaches_the_industry_model_and_the_coefficients_stay_as_they_were(s, monkeypatch):
    d = decision(s)
    seen = {}
    monkeypatch.setattr(mix.industry, 'build', lambda data, *rest: seen.update(data=data) or ('physical', {}))
    d.industry({'annual': []}, {'by_year': [[8, 0.5]]})
    assert seen['data']['economy']['production']['local_working_parts']['staffed'] == {'by_year': [[8, 0.5]]}
    assert 'staffed' not in s['engineering']['economy']['production']['local_working_parts']
    assert d.staffing({'local_working_parts': {'staffed_by_order_year': {'9': 0.2, '8': 0.1}}}) == {'by_year': [[8, 0.1], [9, 0.2]]}
    assert d.staffing({}) is None and d.solves == 1


def test_a_line_that_starts_late_is_served_from_its_first_year_only(s):
    d = decision(s)
    book = d.orders()
    rows = {row['year']: row for row in book['annual']}
    first = account(d, computing=(10.0, 30.0, 1.0, 10.0, 0.0))            # a loss when served from the opening
    for y in d.years:
        first['computing']['years'][y]['quantity'] = rows[y]['lines']['computing']['quantity']
    summary = d.summary(first)
    scale, order, dropped = d.allocate(summary, first, book, book, {y: 0.0 for y in d.years}, starts={'computing': (5.0, 51)})
    assert order == ['computing'] and not dropped
    assert all(scale['computing'][y] == 0 for y in d.years if y < 51) and scale['computing'][51] > 0
    scale, order, dropped = d.allocate(summary, first, book, book, {y: 0.0 for y in d.years})
    assert dropped == ['computing']


def test_a_later_start_is_judged_by_sizing_the_industry_with_the_line_left_out_until_then(s, monkeypatch):
    d = decision(s)
    sized = []

    def industry(book, staffed=None):
        # An industry that costs nothing by the year, and whose founding deployment costs 1e15 where culture is served at the opening.
        first = book['annual'][0]['lines']['culture']['quantity'] > 0
        sized.append(next((row['year'] for row in book['annual'] if row['lines']['culture']['quantity'] > 0), None))
        cash = dict(capital_total_USD2041=0.0, operating_total_USD2041=0.0, capital_external_USD2041=0.0, capital_local_USD2041=0.0)
        annual = [dict(year=row['year'], cash=cash, technical_resident_jobs=0.0) for row in book['annual']]
        return dict(annual=annual, seed=dict(cost=dict(total_USD2041=1e15 if first else 0.0), total_resource_provision_USD2041=0.0)), {}

    monkeypatch.setattr(d, 'industry', industry)
    record = dict(without=({y: 0.0 for y in d.years}, 0.0, {y: 0.0 for y in d.years}))
    scale = {line: {y: 1.0 for y in d.years} for line in demand.LINES}
    worth, start = d.late_start('culture', scale, None, None, record)
    assert sized == [26, 51, 101, 151]                                       # the first year of each later period
    assert start == 26 and worth > 0                                         # the earliest start keeps the most receipts
    rows = {row['year']: row for row in d.orders()['annual']}
    kept = math.fsum((rows[y]['lines']['culture']['receipts_USD2041'] - rows[y]['lines']['culture']['own_cost_USD2041']) * 1.03 ** -y for y in d.years if y >= 26)
    assert worth == pytest.approx(kept)


def test_the_cost_of_a_year_counts_the_financings_allowance_on_capital():
    cash = dict(capital_total_USD2041=110.0, operating_total_USD2041=40.0, capital_external_USD2041=80.0, capital_local_USD2041=20.0)
    physical = dict(annual=[dict(year=8, cash=cash)], seed=dict(cost=dict(total_USD2041=1000.0), total_resource_provision_USD2041=200.0))
    assert mix.yearly_cost(physical) == ({8: 150.0}, 1000.0)
    costs, founding = mix.yearly_cost(physical, 0.35)
    assert costs[8] == pytest.approx(150.0 + 0.35 * 100.0) and founding == pytest.approx(1000.0 + 0.35 * 800.0)


def test_a_job_taken_from_making_working_parts_costs_the_parts_it_would_have_made(s):
    case = dict(scenario='central', ramp='constant_growth', recovery='combined')
    d = mix.Decision(s['programme'], 'constant_growth', s['scenarios'], 'central', run.engineering(case), s['recovery'], 'combined',
                     s['periods'], hotel=62500.0)
    wage = 200000 * s['engineering']['normalization']['USD2026_to_USD2041']
    short = d.displaced({'by_year': [[8, 0.0], [100, 0.5], [150, 1.0]]})
    # Central labour is a quarter of the import price: a job-year makes parts worth four wages and costs one.
    assert short[100] == pytest.approx(3 * wage) and short[60] == pytest.approx(3 * wage)
    assert short[150] == short[200] == 0                                    # every able part is staffed
    assert short[8] == 0                                                     # nothing is yet able to be made
    assert all(v == 0 for v in d.displaced(None).values())
    assert all(v == 0 for v in decision(s).displaced({'by_year': [[8, 0.5]]}).values())   # every part bought
