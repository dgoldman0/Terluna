"""The programme definition: its ramps, and its agreement with the inputs of the staged studies."""
import json
import math
import pytest
from research.studies.enclosed_returns import markets
from research.studies.vacuum_economy import programme as definition, run

P = definition.load()
FIRST, DONE, LAST = (P['calendar'][k] for k in ('first_productive_year', 'completion_year', 'last_year'))
OPENING, FOUNDING, PLANNED = (P['population'][k] for k in ('opening', 'founding_system', 'planned'))


def people(ramp, year):
    return definition.population(P, ramp, year)


@pytest.mark.parametrize('ramp', definition.planned(P))
def test_every_planned_ramp_opens_at_the_opening_population_and_completes_on_the_plan(ramp):
    assert people(ramp, FIRST) == OPENING
    assert math.isclose(people(ramp, DONE), PLANNED, rel_tol=1e-12)
    assert math.isclose(people(ramp, LAST), PLANNED, rel_tol=1e-12)
    rows = definition.schedule(P, ramp)
    assert [r['year'] for r in rows] == list(range(FIRST, LAST + 1))
    assert all(r['added'] >= 0 for r in rows)
    assert math.isclose(OPENING + math.fsum(r['added'] for r in rows), PLANNED, rel_tol=1e-12)


def test_constant_growth_keeps_one_rate_and_passes_the_founding_system_at_about_year_100():
    rates = [definition.growth_rate(P, 'constant_growth', y) for y in range(FIRST + 1, DONE + 1)]
    assert max(rates) - min(rates) < 1e-12
    assert round(100 * rates[0], 1) == 3.3
    assert abs(people('constant_growth', 100) / FOUNDING - 1) < 0.01


def test_the_front_loaded_ramp_stands_the_founding_system_up_by_year_33_and_then_slows():
    stands = P['ramps']['front_loaded']['founding_system_year']
    assert stands == 33 and math.isclose(people('front_loaded', stands), FOUNDING, rel_tol=1e-12)
    early = definition.growth_rate(P, 'front_loaded', stands)
    late = definition.growth_rate(P, 'front_loaded', stands + 1)
    assert round(100 * early, 1) == 12.7 and round(100 * late, 1) == 2.0
    assert all(people('front_loaded', y) >= people('constant_growth', y) * (1 - 1e-12) for y in range(FIRST, LAST + 1))


def test_the_first_pass_ramp_is_the_schedule_of_the_staged_order_book():
    rules = markets.inputs()[0]
    assert all(people('first_pass', y) == markets.population(y, rules) for y in range(FIRST, LAST + 1))
    assert round(people('first_pass', 100)) == 222240


def test_the_staged_inputs_agree_with_the_definition():
    s = run.sources()                       # raises where the demand inputs or the financing rules differ
    assert s['programme'] == P and s['money']['horizons'] == P['horizons']


@pytest.mark.parametrize('change', [
    lambda p: p['calendar'].update(completion_year=5),
    lambda p: p['population'].update(founding_system=10000),
    lambda p: p['ramps']['front_loaded'].update(founding_system_year=250),
])
def test_an_inconsistent_definition_is_refused(change, tmp_path):
    altered = json.loads(json.dumps(P))
    change(altered)
    path = tmp_path / 'programme.json'
    path.write_text(json.dumps(altered))
    with pytest.raises(ValueError):
        definition.load(path)


def test_an_unknown_kind_of_ramp_is_refused():
    altered = json.loads(json.dumps(P))
    altered['ramps']['constant_growth']['kind'] = 'stepwise'
    with pytest.raises(ValueError):
        definition.population(altered, 'constant_growth', 50)


def test_the_solved_ramps_state_a_rule_and_take_their_residents_only_from_a_solution():
    solved = sorted(name for name, ramp in P['ramps'].items() if ramp['kind'] == 'solved')
    assert solved == ['resident_led', 'resident_led_beyond_plan'] and not set(solved) & set(definition.planned(P))
    assert P['ramps']['resident_led']['limit'] == 'planned' and P['ramps']['resident_led_beyond_plan']['limit'] is None
    with pytest.raises(ValueError):
        definition.population(P, 'resident_led', 50)                       # no solution given
    residents = [[y, OPENING * 1.01 ** (y - FIRST)] for y in range(FIRST, LAST + 1)]
    given = definition.solved(P, 'resident_led_beyond_plan', residents)
    assert 'residents' not in P['ramps']['resident_led_beyond_plan']        # the definition itself is left alone
    assert definition.population(given, 'resident_led_beyond_plan', FIRST) == OPENING
    assert definition.population(given, 'resident_led_beyond_plan', 120) == residents[120 - FIRST][1]
    assert definition.population(given, 'resident_led_beyond_plan', LAST + 2) == residents[-1][1]   # held after the last year


@pytest.mark.parametrize('change', [
    lambda r: r[:-1],                                                     # a year missing
    lambda r: [[r[0][0], 2 * OPENING]] + r[1:],                           # not opening on the opening population
    lambda r: r[:50] + [[r[50][0], r[49][1] / 2]] + r[51:],               # falling
])
def test_a_solution_that_misses_a_year_does_not_open_on_the_opening_population_or_falls_is_refused(change):
    residents = [[y, OPENING * 1.02 ** (y - FIRST)] for y in range(FIRST, LAST + 1)]
    with pytest.raises(ValueError):
        definition.solved(P, 'resident_led', change(residents))
    with pytest.raises(ValueError):
        definition.solved(P, 'constant_growth', residents)                 # a planned ramp takes no solution


def test_a_definition_that_lists_residents_for_a_solved_ramp_is_refused(tmp_path):
    altered = json.loads(json.dumps(P))
    altered['ramps']['resident_led']['residents'] = [[FIRST, OPENING]]
    path = tmp_path / 'programme.json'
    path.write_text(json.dumps(altered))
    with pytest.raises(ValueError):
        definition.load(path)
