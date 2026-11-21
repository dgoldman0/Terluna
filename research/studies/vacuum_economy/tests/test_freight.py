"""The tugs chosen for each scenario, and the freight a case flies."""
import json
import pytest
from research.studies.vacuum_economy import freight, inputs

SAVED = json.loads(freight.OUTPUT.read_text()) if freight.OUTPUT.exists() else None


def test_a_schedule_merges_periods_that_keep_the_same_kind():
    periods = [(0, 7), (8, 25), (26, 50), (51, 100)]
    assert freight.schedule(['chemical'] * 4, periods) == 'chemical'
    assert freight.schedule(['chemical', 'chemical', 'electric_water', 'electric_water'], periods) == {'by_year': [[0, 'chemical'], [26, 'electric_water']]}


def test_a_case_flies_the_founding_freight_in_the_first_pass_and_where_it_says_so():
    assert inputs.freight_of(dict(scenario='first_pass')) == 'founding'
    assert inputs.freight_of(dict(scenario='central')) == 'chosen'
    founding = inputs.engineering(dict(scenario='central', ramp='constant_growth', freight='founding'))
    assert 'source_reach' not in founding and 'source_propulsion' not in founding['industrial']['transport']
    chemical = inputs.engineering(dict(scenario='central', ramp='constant_growth', freight='chemical', reach='p95'))
    assert chemical['industrial']['transport']['source_propulsion'] == 'chemical'
    assert chemical['source_reach']['percentile'] == 'p95' and chemical['source_reach']['within_kg'] == inputs.reach()['within_kg']['p95']
    with pytest.raises(ValueError):
        inputs.engineering(dict(scenario='central', ramp='constant_growth', freight='sails'))


@pytest.mark.skipif(SAVED is None, reason='The choice of tugs has not been made')
def test_the_saved_choice_is_current_and_takes_the_cheapest_kind_standing_in_each_period():
    assert SAVED['dependencies'] == freight.dependencies()
    for scenario, found in SAVED['scenarios'].items():
        for key, options in (('chosen', freight.CANDIDATES), ('chosen_with_target', (*freight.CANDIDATES, freight.TARGET))):
            unreached = [(p['first_year'], p['last_year']) for p in found[key]['periods_no_kind_reaches']]
            for i, period in enumerate(found[key]['by_period']):
                standing = {k: v['periods'][i] for k, v in found['comparison'].items() if k in options and v['periods']
                            and v['periods'][i].get('available', True) and v['periods'][i]['feasible'] and v['periods'][i]['within_world_xenon_supply']}
                if not standing:
                    # No kind reaches the period's water: the schedule keeps the kind before it, and records the period.
                    assert (period['first_year'], period['last_year']) in unreached
                    assert period['kind'] == (found[key]['by_period'][i - 1]['kind'] if i else options[0])
                    continue
                assert (period['first_year'], period['last_year']) not in unreached
                assert period['kind'] in standing
                assert standing[period['kind']]['industry_cost_USD2041'] == min(p['industry_cost_USD2041'] for p in standing.values())
            assert inputs.chosen_freight()['scenarios'][scenario][key]['schedule'] == found[key]['schedule']
