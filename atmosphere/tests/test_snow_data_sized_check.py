"""The data-sized rule: repeat-crystal scatter, verdicts and the committed product."""
import json
import math

import pytest

from atmosphere.snow import data_sized_check as check


def test_nominal_assignment_uses_the_half_degree_window():
    assert check.nominal(-5.3) == -5.3
    assert check.nominal(-5.8) == -5.3
    assert check.nominal(-4.8) == -5.3
    assert check.nominal(-5.81) is None
    assert check.nominal(-2.8) is None
    assert check.nominal(-22.4) == -22.0


def crystal(row, t, age, mass):
    return dict(excel_row=row, temperature_C=t, nominal_C=check.nominal(t), age_min=age,
                measured={'mass_kg': mass, 'equatorial_dimension_m': mass, 'polar_dimension_m': mass})


def test_pooled_scatter_takes_deviations_from_each_group_mean():
    # Two groups: ln values {0, 2} and {1, 1, 4}; within-group sums of squares 2 and 6 on 1 + 2 dof.
    items = [crystal(1, -5.3, 5., 1.), crystal(2, -5.3, 5., math.e ** 2),
             crystal(3, -14.4, 10., math.e), crystal(4, -14.4, 10., math.e), crystal(5, -14.4, 10., math.e ** 4),
             crystal(6, -14.4, 12., 7.)]  # a lone crystal adds nothing
    s = check.scatter(items)
    assert s['mass_kg']['degrees_of_freedom'] == 3
    assert s['mass_kg']['pooled_sd_ln'] == pytest.approx(math.sqrt(8 / 3))
    assert s['mass_kg']['by_temperature']['-14.4']['sd_ln'] == pytest.approx(math.sqrt(3))
    assert [g['crystals'] for g in s['replicate_groups']] == [3, 2]


def test_verdict_boundary_is_two_standard_deviations():
    sd = .2
    inside = check.verdict(math.exp(2 * sd) * (1 - 1e-12), 1., sd)
    outside = check.verdict(math.exp(2 * sd) * (1 + 1e-9), 1., sd)
    assert inside['passed'] and not outside['passed']
    assert check.verdict(1., math.exp(-.5), sd)['z'] == pytest.approx(2.5)


def test_speed_verdict_is_absolute_three_centimetres_per_second():
    reference = check.lab.fall_speed(-5.3, 10.)
    assert check.speed_verdict(reference + .0299, -5.3, 10.)['passed']
    assert not check.speed_verdict(reference - .0301, -5.3, 10.)['passed']
    assert not check.speed_defined(-14.4, 2.) and check.speed_defined(-14.4, 30.)


@pytest.fixture(scope='module')
def product():
    if not check.PRODUCT.exists():
        pytest.skip('data-sized product not assembled')
    return json.loads(check.PRODUCT.read_text())


def test_committed_scatter_recomputes_from_the_measurements(product):
    s = check.scatter(check.crystals())
    for q in check.QUANTITIES:
        assert product['scatter'][q]['pooled_sd_ln'] == pytest.approx(s[q]['pooled_sd_ln'], rel=1e-12)
        assert product['scatter'][q]['degrees_of_freedom'] == s[q]['degrees_of_freedom']
    assert product['crystals']['total'] == 169


def test_committed_verdicts_agree_with_their_values(product):
    sd = {q: product['scatter'][q]['pooled_sd_ln'] for q in check.QUANTITIES}
    for name, result in product['results'].items():
        failures = 0
        for r in result['records']:
            if not r['compared']:
                continue
            for q, v in r['quantities'].items():
                assert v['ratio'] == pytest.approx(v['model'] / v['measured'], rel=1e-12)
                assert v['passed'] == (abs(math.log(v['ratio'])) <= 2 * sd[q])
                failures += not v['passed']
        for s in result['speeds']:
            assert s['passed'] == (abs(s['difference']) <= check.SPEED_TOLERANCE + 1e-12)
            failures += not s['passed']
        assert failures == len(result['outside']), name
    assert all(p['passed'] for p in product['numerical_checks']['calibration_parity'])
