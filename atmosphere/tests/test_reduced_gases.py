"""The reduced-gas forcing product: its binding, its Earth check and the arithmetic it rests on."""
import csv
import hashlib
import json
import math

import pytest

from atmosphere.radiative_convective import reduced_gases as R

pytestmark = pytest.mark.skipif(not R.OUT.exists(), reason='run python -m atmosphere.radiative_convective.reduced_gases')


@pytest.fixture(scope='module')
def product():
    return json.loads(R.OUT.read_text())


def close(a, b, rel=1e-6):
    return math.isclose(a, b, rel_tol=rel, abs_tol=1e-12)


def test_product_binds_its_code_and_inputs(product):
    assert product['schema'] == R.SCHEMA
    for path, h in product['producer']['files'].items():
        assert hashlib.sha256((R.ROOT / path).read_bytes()).hexdigest() == h, path
    manifest = R.input_records()
    for key, h in product['producer']['inputs'].items():
        if key.startswith('inputs/'):
            assert manifest[key.split('/', 1)[1]] == h, key
        else:
            assert hashlib.sha256((R.ROOT / key).read_bytes()).hexdigest() == h, key


def test_forcing_grows_with_each_gas(product):
    for column in product['columns'].values():
        for ts, cases in column.items():
            for gas in ('ch4', 'n2o'):
                toa = [c['toa'] for c in cases[gas].values()]
                assert all(b > a > 0 for a, b in zip(toa, toa[1:])), (gas, ts)


def test_the_gases_overlap(product):
    for column in product['columns'].values():
        for cases in column.values():
            for key, c in cases['both'].items():
                ch4, n2o = key.split('+')
                assert c['toa'] < cases['ch4'][ch4]['toa'] + cases['n2o'][n2o]['toa']


def test_the_earth_control_meets_ar6_within_a_factor_of_two(product):
    for gas, check in product['ar6_check'].items():
        assert 0.5 < check['ratio'] < 2.0, gas


def test_the_deep_column_takes_more_forcing_per_ppm_from_thin_amounts(product):
    earth, moon = product['columns']['earth_1atm']['295'], product['columns']['moon_12atm']['295']
    assert moon['n2o']['0.2701']['toa'] > earth['n2o']['0.2701']['toa']


def test_warming_is_forcing_over_the_responses(product):
    lam = product['responses']['heat_block_w_m2_k']
    response = product['responses']['clear_sky_w_m2_k']
    for gas in ('CH4', 'N2O'):
        ratio = product['ar6_check'][gas]['ratio']
        for c in product['moon_at_design'][gas].values():
            assert close(c['allsky_w_m2'], ratio * c['clear_sky_toa_w_m2'])
            assert close(c['warming_k_clear_sky_response'], c['clear_sky_toa_w_m2'] / response)
            assert close(c['warming_k_heat_block'][0], c['allsky_w_m2'] / lam[1])


def test_the_clear_sky_response_follows_the_inverse_climate_rows():
    rows = list(csv.DictReader(R.INVERSE.open()))
    response, used = R.clear_sky_response(rows, 'moon_1.2atm')
    assert set(used['olr']) == set(used['asr']) == {290.0, 300.0}
    assert 0.5 < response < 3.0
