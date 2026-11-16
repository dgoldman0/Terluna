"""The resident-led schedules: how a schedule follows from its rates, how a trial is judged, how a period's rate is
searched, and the saved product."""
import json
import math
import pytest
from research.studies.vacuum_economy import inputs, programme as definition, residents

P = definition.load()
FIRST, DONE, LAST = (P['calendar'][k] for k in ('first_productive_year', 'completion_year', 'last_year'))
OPENING, PLANNED = P['population']['opening'], P['population']['planned']


def test_a_schedule_opens_on_the_opening_population_grows_at_each_periods_rate_and_is_held_after_completion():
    rates = [0.05, 0.0, 0.02]
    rows = residents.schedule(P, rates)
    assert [y for y, _ in rows] == list(range(FIRST, LAST + 1)) and rows[0][1] == OPENING
    at = dict(rows)
    first, second, third = residents.PERIODS[:3]
    assert math.isclose(at[first[1]], OPENING * 1.05 ** (first[1] - FIRST), rel_tol=1e-12)
    assert at[second[1]] == at[first[1]]
    assert math.isclose(at[third[1]], at[first[1]] * 1.02 ** (third[1] - third[0] + 1), rel_tol=1e-12)
    assert all(at[y] == at[third[1]] for y in range(third[1], LAST + 1))      # periods without a rate have none
    limited = dict(residents.schedule(P, [0.3] * len(residents.PERIODS), PLANNED))
    reached = next(y for y in range(FIRST, LAST + 1) if limited[y] >= PLANNED)
    assert all(limited[y] == PLANNED for y in range(reached, LAST + 1))
    definition.solved(P, residents.RAMP, rows)                                  # a schedule the definition accepts
    with pytest.raises(ValueError):
        residents.schedule(P, [0.01] * (len(residents.PERIODS) + 1))
    with pytest.raises(ValueError):
        residents.schedule(P, [-0.01])


def test_the_periods_run_without_gaps_from_the_year_after_opening_to_completion():
    assert residents.PERIODS[0][0] == FIRST + 1 and residents.PERIODS[-1][1] == DONE
    assert all(a[1] + 1 == b[0] for a, b in zip(residents.PERIODS, residents.PERIODS[1:]))


def test_money_is_the_least_cash_above_the_buffer_with_the_opening_runway_and_the_equity_given():
    rules = dict(liquidity_buffer_years=1, initial_runway_years=5)
    path = dict(annual=[dict(year=y, operating_cost=10.0) for y in (8, 9, 20)])
    financed = dict(annual=[dict(year=8, pre_topup_cash=60.0, contribution=0.0), dict(year=9, pre_topup_cash=35.0, contribution=5.0),
                            dict(year=20, pre_topup_cash=30.0, contribution=0.0)])
    least, equity = residents.money(path, financed, rules)
    # Year 8 must hold five years of costs (50), year 9 four (40), year 20 one (10).
    assert least == pytest.approx(-5.0) and equity == pytest.approx(5.0)
    held = dict(annual=[dict(year=200, cash=5.0, reserves=2.0, owner_custody=1.0, debt=3.0)])
    assert residents.holding(held, 200) == pytest.approx(5.0)


def test_the_survey_pools_nest_the_wider_screen_within_the_broader_population():
    stock = inputs.pool()
    assert set(stock) == {'wider', 'broader'}
    assert all(0 < stock['wider'][k] < stock['broader'][k] for k in inputs.PERCENTILES)
    assert stock['wider']['p05'] < stock['wider']['p50'] < stock['wider']['p95']
    assert round(stock['wider']['p50'] / 1e12, 1) == 39.1 and round(stock['broader']['p50'] / 1e12) == 961


def test_the_water_drawn_is_the_feed_times_its_grade_and_passes_each_percentile_in_its_year(monkeypatch):
    monkeypatch.setattr(inputs, 'pool', lambda: dict(wider=dict(p05=50.0, p50=100.0, p95=1e9), broader=dict(p05=150.0, p50=1e9, p95=1e9)))
    loops = dict(annual=[dict(year=y, source_feed_kg=10.0) for y in range(0, LAST + 1)])
    coefficients = dict(industrial=dict(source=dict(water_grade=0.1)))
    out = inputs.water_drawn(None, loops, coefficients, (100, 200))
    assert out['contained_water_drawn_kg'] == {'100': pytest.approx(101.0), '200': pytest.approx(201.0)}
    assert out['year_drawn_past'] == dict(wider=dict(p05=50, p50=100, p95=None), broader=dict(p05=150, p50=None, p95=None))
    assert out['share_of_median_drawn_by_last_year']['wider'] == pytest.approx((LAST + 1) / 100.0)


class Fake:
    """A trial in which period j sustains any rate up to its own threshold, with the distance to it as the margin."""

    def __init__(self, thresholds, opening=True):
        self.thresholds, self.opening, self.calls = thresholds, opening, 0

    def __call__(self, rates, choice, limit=None, ramp=residents.RAMP):
        self.calls += 1
        if not rates:
            return dict(rates=[], sustained=self.opening, binding='money', margin=1.0 if self.opening else -1.0, opening_funding_USD2041=1.0)
        margin = min(self.thresholds[j] - r for j, r in enumerate(rates))
        return dict(rates=list(rates), sustained=margin >= 0, binding='money' if len(rates) % 2 else 'hands', margin=margin,
                    opening_funding_USD2041=1.0)


def solver(monkeypatch, fake):
    s = residents.Solver('central', log=lambda *a: None)
    monkeypatch.setattr(s, 'trial', fake)
    return s


def test_a_periods_search_brackets_the_highest_sustained_rate_within_the_stated_resolution(monkeypatch):
    for threshold in (0.0473, 0.0002, 0.21):
        s = solver(monkeypatch, Fake([threshold]))
        lo, hi, at_lo, at_hi, tried = s.period([], 0, None, 0.033, 2.0)
        assert lo <= threshold < hi and len(tried) <= residents.TRIALS
        assert hi - lo <= max(residents.RESOLUTION, residents.RELATIVE * lo) + 1e-12
    s = solver(monkeypatch, Fake([1.0]))
    lo, hi, _, _, tried = s.period([], 0, None, 0.033, 2.0)
    assert lo == residents.HIGHEST and hi is None                             # sustained at the highest rate tried


def test_the_march_finds_each_periods_rate_in_order_and_stops_where_the_opening_is_not_sustained(monkeypatch):
    thresholds = [0.01 * (j + 1) for j in range(len(residents.PERIODS))]
    s = solver(monkeypatch, Fake(thresholds))
    rates, record, refused = s.march(None)
    assert refused is None and len(rates) == len(residents.PERIODS)
    assert all(t * (1 - residents.RELATIVE) - residents.RESOLUTION <= r <= t for r, t in zip(rates, thresholds))
    assert [(p['first_year'], p['last_year']) for p in record] == [tuple(p) for p in residents.PERIODS]
    s = solver(monkeypatch, Fake(thresholds, opening=False))
    rates, record, refused = s.march(None)
    assert rates == [0.0] * len(residents.PERIODS) and refused['sustained'] is False


SAVED = json.loads(residents.PRODUCT.read_text()) if residents.PRODUCT.exists() else None
needs_product = pytest.mark.skipif(SAVED is None, reason='residents.json not yet produced')


@needs_product
def test_the_saved_product_names_its_producer_and_inputs_and_every_scenario():
    assert SAVED['schema'] == residents.SCHEMA
    assert SAVED['producer']['sha256'] == inputs.digest(residents.__file__)
    assert SAVED['dependencies'] == residents.dependencies()
    assert sorted(SAVED['scenarios']) == sorted(residents.SCENARIOS)
    assert SAVED['periods'] == [list(p) for p in residents.PERIODS]


@needs_product
@pytest.mark.parametrize('scenario', residents.SCENARIOS)
def test_each_saved_schedule_follows_its_rates_and_the_limited_one_holds_the_plan_from_the_year_it_is_reached(scenario):
    entry = SAVED['scenarios'][scenario]
    rates = [p['rate'] for p in entry['periods']]
    beyond = entry['schedules'][residents.RAMP]
    limited = entry['schedules'][residents.LIMITED]
    assert beyond == residents.schedule(P, rates) and limited == residents.schedule(P, rates, PLANNED)
    definition.solved(P, residents.RAMP, beyond)
    definition.solved(P, residents.LIMITED, limited)
    reached = entry['planned_population_reached_year']
    if reached is None:
        assert limited == beyond and max(v for _, v in beyond) < PLANNED
        assert entry['limited_trial'] is None
    else:
        assert dict(beyond)[reached] >= PLANNED > dict(beyond)[reached - 1]
        assert all(v == PLANNED for y, v in limited if y >= reached)
        assert entry['limited_trial']['sustained']
    for p in entry['periods']:
        assert p['lowest_unsustained_rate'] is None or p['rate'] < p['lowest_unsustained_rate']
    assert entry['final_trial']['sustained'] == entry['opening_settlement_sustained']
