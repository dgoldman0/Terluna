import math

import pytest

from engineering.heat_recovery import radiator_area, recover_heat
from shared.constants import STEFAN_BOLTZMANN


def test_recovered_work_returns_as_heat_already_counted():
    r = recover_heat(1e6, 600., 300., .5, .1)
    assert r['recovered_electric_W'] == 250000.
    assert r['processor_electric_W'] == 225000.
    assert r['recovered_overhead_W'] == 25000.
    assert r['engine_rejection_W'] == 750000.
    assert r['internal_heat_to_radiators_W'] == 1e6
    assert r['energy_balance_residual_W'] == 0.
    external = recover_heat(1e6, 600., 300., .5, .1, 10000.)
    assert external['added_internal_heat_W'] == 10000.


@pytest.mark.parametrize('hot,cold,fraction', [(300., 300., 1.), (600., 300., 0.)])
def test_recovered_work_needs_a_gradient_and_a_converter(hot, cold, fraction):
    assert recover_heat(1e6, hot, cold, fraction)['recovered_electric_W'] == 0.


def test_carnot_bound_and_linear_scaling():
    limit = recover_heat(1e6, 600., 300., 1.)
    assert limit['recovered_electric_W'] == 500000.
    for fraction in [0., .1, .5, 1.]:
        r = recover_heat(1e6, 600., 300., fraction)
        doubled = recover_heat(2e6, 600., 300., fraction)
        assert 0 <= r['recovered_electric_W'] <= limit['recovered_electric_W']
        assert doubled['processor_electric_W'] == 2*r['processor_electric_W']


def test_rejecting_heat_colder_needs_more_area():
    recovered = recover_heat(1e6, 600., 300., .5)
    current = radiator_area(recovered['internal_heat_to_radiators_W'], 300., .9)
    same_sink = radiator_area(1e6, 300., .9)
    hot_sink = radiator_area(1e6, 600., .9)
    assert current == same_sink
    assert current['emitting_area_m2']/hot_sink['emitting_area_m2'] == pytest.approx(16.)
    assert radiator_area(2e6, 300., .9)['emitting_area_m2'] == 2*current['emitting_area_m2']


def test_background_and_absorbed_sunlight_count_in_energy_balance():
    dark = radiator_area(1e6, 300., .9)
    lit = radiator_area(1e6, 300., .9, 100., 100.)
    assert lit['emitting_area_m2'] > dark['emitting_area_m2']
    assert lit['net_heat_rejected_W'] == pytest.approx(1e6)
    assert lit['gross_radiated_W'] == pytest.approx(1e6+lit['absorbed_environment_W'])
    with pytest.raises(ValueError, match='net heat rejection'):
        radiator_area(1e6, 300., .9, absorbed_external_W_m2=.9*STEFAN_BOLTZMANN*300.**4)


def test_zero_heat_requires_no_area_or_recovered_power():
    assert radiator_area(0., 300., .9)['emitting_area_m2'] == 0.
    assert recover_heat(0., 600., 300., .5)['recovered_electric_W'] == 0.


@pytest.mark.parametrize('changes', [dict(heat_W=-1), dict(heat_W=math.nan),
    dict(hot_K=math.inf), dict(cold_K=0), dict(cold_K=601),
    dict(carnot_fraction=1.01), dict(recovered_overhead_fraction=-.1),
    dict(recovered_overhead_fraction=1.01), dict(external_support_W=-1)])
def test_recovery_rejects_invalid_inputs(changes):
    args = dict(heat_W=1e6, hot_K=600., cold_K=300., carnot_fraction=.5)
    with pytest.raises(ValueError):
        recover_heat(**(args | changes))


@pytest.mark.parametrize('changes', [dict(heat_W=-1), dict(temperature_K=0),
    dict(temperature_K=math.nan), dict(emissivity=0), dict(emissivity=1.1),
    dict(background_K=-1), dict(background_K=300), dict(absorbed_external_W_m2=-1)])
def test_radiator_rejects_invalid_inputs(changes):
    args = dict(heat_W=1e6, temperature_K=300., emissivity=.9)
    with pytest.raises(ValueError):
        radiator_area(**(args | changes))
