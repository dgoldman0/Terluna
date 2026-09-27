"""Checks on the fruit model: the growth curve, development in degree-days, carbon accounting and stored results."""
import json

import numpy as np
import pytest

from biosphere.canopy import fruit as fr
from biosphere.canopy import plant


def test_growth_curve_runs_from_nothing_to_the_harvest_mass_fastest_at_its_peak():
    s = np.linspace(0.0, 1.0, 20001)
    for peak in (0.4, 0.5, 0.6):
        w = fr.filled(s, peak)
        assert w[0] == 0.0 and w[-1] == pytest.approx(1.0)
        assert np.all(np.diff(w) >= 0)
        assert s[np.argmax(np.gradient(w, s))] == pytest.approx(peak, abs=1e-3)


def test_development_takes_its_degree_days_at_a_constant_temperature():
    fruit = fr.Fruit()
    temp = np.full(plant.STEPS, 24.0)
    dt_days = 29.53 / plant.STEPS
    co = fr.cohort(temp, dt_days, fruit, start=100)
    assert co['days'] == pytest.approx(fruit.degree_days / 14.0, abs=dt_days)
    assert co['added'].sum() == pytest.approx(1.0)
    assert co['index'][0] == 100 and co['index'].max() < plant.STEPS
    assert fr.cohort(np.full(plant.STEPS, 9.0), dt_days, fruit, 0) is None      # below the base it never ripens


def flat(net_day, net_night, steps=plant.STEPS, period_s=plant.MOON_DAY_S):
    hour = (np.arange(steps) + 0.5) / steps * 360.0 - 180.0
    up = np.abs(hour) < 90.0
    net = np.where(up, net_day, net_night)
    return dict(hour=hour, up=up, dt=period_s / steps, net=net, par=np.where(up, 500.0, 0.0))


def test_a_fed_fruit_costs_its_carbon_growth_respiration_and_upkeep():
    fruit, f = fr.Fruit(), flat(10.0, -2.0)
    co = fr.cohort(np.full(plant.STEPS, 22.0), f['dt'] / 86400.0, fruit, 0)
    none = fr.feed(f['net'], f['dt'], fruit, co, 0.0, 'fed')
    some = fr.feed(f['net'], f['dt'], fruit, co, 0.3, 'fed')
    assert some['size'] == pytest.approx(1.0)
    cost = some['fruit_carbon'] + some['fruit_maintenance'] / (1 + fr.GROWTH_RESPIRATION)
    assert none['plant_growth'] - some['plant_growth'] == pytest.approx(cost, rel=1e-9)
    assert some['fruit_carbon'] == pytest.approx(0.3 * fruit.carbon_g)
    assert some['store'] > none['store']


def test_a_daylight_fruit_keeps_only_the_growth_that_falls_in_the_plants_surplus():
    fruit, f = fr.Fruit(), flat(10.0, -2.0)
    co = fr.cohort(np.full(plant.STEPS, 22.0), f['dt'] / 86400.0, fruit, 0)
    day = fr.feed(f['net'], f['dt'], fruit, co, 0.01, 'daylight')
    assert day['size'] == pytest.approx(co['added'][f['up'][co['index']]].sum(), rel=1e-9)
    always = flat(10.0, 10.0)
    assert fr.feed(always['net'], always['dt'], fruit, co, 0.01, 'daylight')['size'] == pytest.approx(1.0)
    crowded = fr.feed(f['net'], f['dt'], fruit, co, 50.0, 'daylight')
    assert crowded['size'] < day['size']                                         # many fruits share the surplus


@pytest.fixture(scope='module')
def product():
    path = fr.RESULTS / 'fruit.json'
    if not path.is_file():
        pytest.skip('results/fruit.json has not been generated')
    return json.loads(path.read_text())


def test_stored_product(product):
    assert product['schema'] == fr.SCHEMA
    for side, bands in product['results'].items():
        for band, strategies in bands.items():
            sunlit = product['calendar'][side][band]['sunlit_days']
            for name, fruits in strategies.items():
                day, melon = fruits['day_fruit'], fruits['watermelon']
                for policy in fr.POLICIES:
                    # the day fruit, set at sunrise, is ripe at sunset and never draws on the night store
                    d = day['best'][policy]
                    assert abs(d['set_hour_angle'] + 90.0) < 10.0
                    assert d['days_to_harvest'] == pytest.approx(sunlit, abs=0.1)
                    assert d['growth_while_plant_short'] < 0.01
                    for r in d['shares'].values():
                        assert r['size'] == pytest.approx(1.0, abs=0.01)
                        assert r['store'] == pytest.approx(day['store_without_fruit'], rel=0.02)
                    # today's watermelon lives through a night: shrunken by daylight alone, or a larger store
                    m = melon['best'][policy]
                    assert abs(m['set_hour_angle']) < 90.0 and m['days_to_harvest'] > 2 * sunlit
                    shares = [m['shares'][k] for k in ('0.3', '0.5', '0.7')]
                    assert np.all(np.diff([r['harvest_kg'] for r in shares]) > 0)
                    assert np.all(np.diff([r['plant_growth'] for r in shares]) < 0)
                    for r in shares:
                        if policy == 'fed':
                            assert r['size'] == pytest.approx(1.0) and r['store'] > 1.2 * melon['store_without_fruit']
                        else:
                            assert r['size'] < 0.8
    programs = product['programs']['earth_like']
    by_dd = {row['degree_days']: row for row in programs}
    assert by_dd[199.0]['daylight']['size'] > 0.99 > by_dd[300.0]['daylight']['size']
    assert by_dd[199.0]['fed']['store'] < by_dd[260.0]['fed']['store'] < by_dd[300.0]['fed']['store']
    near = product['results']['near']['0']['earth_like']['day_fruit']['best']['fed']['shares']['0.6']
    earth = product['earth']['best']['fed']['shares']['0.6']
    assert near['harvest_kg'] > earth['harvest_kg_per_29_53_days']
