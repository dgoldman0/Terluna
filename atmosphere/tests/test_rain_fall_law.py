"""The fall-law product binds its code and inputs, reproduces from them, integrates Morrison's spectra exactly enough
to give a power law's closed forms back, tabulates the solver's moments as the patched scheme reads them, and carries
rain's self-collection to the Moon by the collision kernels."""
import hashlib
import json
import math

import numpy as np
import pytest

from shared.constants import CM1_RAIN_FALL_COEFFICIENT, CM1_RAIN_FALL_EXPONENT, CM1_RAIN_SELF_COLLECTION
from shared.provenance import constants_changed
from atmosphere.rain import fall_law, product as rain, targets as climate_targets

PRODUCT = json.loads(fall_law.OUT.read_text())
DATA, BINDING = climate_targets.load()
AIR = rain.target_air(climate_targets.read(DATA))[fall_law.TARGET]


@pytest.fixture(scope='module')
def laws():
    cols = rain.columns(AIR)
    factor = rain.cm1_density_factor(cols['moon'].state(0.0)['density_kg_m3'])
    d, v, _ = fall_law.solved_speeds(cols['moon'])
    ed, ev, earth_screen = fall_law.solved_speeds(cols['earth'])
    return dict(moon=fall_law.speed_law(d, v), reference=fall_law.speed_law(d, v / factor),
                earth=fall_law.speed_law(ed, ev), earth_screen=earth_screen)


def test_header_and_bindings():
    assert PRODUCT['schema'] == fall_law.SCHEMA
    assert PRODUCT['evidence'] and PRODUCT['reading_rule'] and PRODUCT['assumptions']
    producer = PRODUCT['producer']
    assert set(producer['files']) == set(fall_law.CODE)
    for path, digest in {**producer['files'], **producer['inputs']}.items():
        assert hashlib.sha256((fall_law.ROOT / path).read_bytes()).hexdigest() == digest, path
    assert set(BINDING['inputs']) <= set(producer['inputs'])
    assert not constants_changed(producer['constants'])


def test_the_quadrature_gives_a_power_laws_closed_forms_back():
    lam = fall_law.table_lambda()
    for a, b in ((CM1_RAIN_FALL_COEFFICIENT, CM1_RAIN_FALL_EXPONENT), (150.0, 0.68), (2000.0, 1.09)):
        found = fall_law.moments(fall_law.power_law(a, b), lam)
        assert np.allclose(found, fall_law.closed_forms(a, b, lam), rtol=1e-5, atol=0.0), (a, b)
    assert PRODUCT['checks']['power_law_quadrature_max_error'] < 1e-5


def test_the_table_holds_the_solvers_moments_in_the_reference_air(laws):
    table = PRODUCT['table']
    lam = fall_law.table_lambda()
    assert table['points'] == lam.size == fall_law.TABLE_LAMBDA[2]
    assert table['ln_lambda_first'] == pytest.approx(math.log(fall_law.TABLE_LAMBDA[0]), abs=1e-12)
    nodes = fall_law.moments(laws['reference'], lam)
    assert np.allclose(fall_law.interpolate(table, lam), nodes, rtol=1e-9, atol=0.0)
    middle = np.sqrt(lam[1:] * lam[:-1])
    assert np.allclose(fall_law.interpolate(table, middle), fall_law.moments(laws['reference'], middle), rtol=1e-4,
                       atol=0.0)
    # Held at the ends beyond the table.
    ends = fall_law.interpolate(table, np.array([10.0, 1e7]))
    assert np.allclose(ends[:, 0], nodes[:, 0]) and np.allclose(ends[:, 1], nodes[:, -1])


def test_the_moments_behave_as_a_fall_law_must(laws):
    table, air = PRODUCT['table'], PRODUCT['air']
    mass, number = np.exp(table['ln_mass_weighted']), np.exp(table['ln_number_weighted'])
    assert np.all(np.diff(mass) < 0) and np.all(np.diff(number) < 0)        # wider spectra fall faster
    assert np.all(number < mass)                                              # small drops are many and slow
    assert mass.max() * air['density_factor'] < air['screen_speed_m_s']      # never faster than the largest drop
    assert 0.98 < PRODUCT['checks']['stokes_number_weighted_at_smallest_mean'] <= 1.0
    # The solver's lunar speeds at 2 m, as rain_drops.json gives them.
    drops = json.loads(rain.OUT.read_text())['fall_speed'][fall_law.TARGET]
    for mm, v in zip(PRODUCT['speeds']['diameter_mm'], PRODUCT['speeds']['solver_m_s']):
        if mm in rain.SPEED_MM:
            assert v == pytest.approx(drops['moon']['terminal_speed_m_s'][rain.SPEED_MM.index(mm)], rel=2e-3), mm


def test_no_power_law_follows_the_solver_much_better_than_the_present_one():
    present, best = PRODUCT['power_laws']['present'], PRODUCT['power_laws']['best']
    worst = lambda row: max(abs(math.log(x)) for pair in row.values() for x in pair)
    assert worst(best) > 0.5 * worst(present)
    assert PRODUCT['comparison']['present']['a_si'] == pytest.approx(
        CM1_RAIN_FALL_COEFFICIENT * rain.cm1_gravity_factor(fall_law.MOON_SURFACE_GRAVITY))


def test_self_collection_follows_the_kernels(laws):
    sc = PRODUCT['self_collection']
    lo, hi = fall_law.SMALLEST_RAIN_M, laws['earth_screen']
    same = fall_law.kernel_sum(laws['earth'], 1.0e-3, lo, hi) / fall_law.kernel_sum(laws['earth'], 1.0e-3, lo, hi)
    assert same == pytest.approx(1.0)
    ratios = [sc['kernel_ratio_by_mean_diameter_mm'][k] for k in ('0.60', '1.00', '1.68', '2.50')]
    assert all(r > 1 for r in ratios) and all(np.diff(ratios) < 0)            # lunar drops meet more slowly
    assert sc['kernel_ratio'] == sc['kernel_ratio_by_mean_diameter_mm']['1.68']
    assert sc['lunar_rate_m3_kg_s'] == round(CM1_RAIN_SELF_COLLECTION / sc['kernel_ratio'], 4)
    assert PRODUCT['namelist'] == dict(rain_fall_law=1, rain_selfcoll=sc['lunar_rate_m3_kg_s'])


def test_values_reproduce_from_the_code():
    fresh = fall_law.build()
    for key in ('ln_mass_weighted', 'ln_number_weighted', 'ln_ventilation'):
        assert np.allclose(fresh['table'][key], PRODUCT['table'][key], rtol=0.0, atol=1e-12), key
    assert fresh['namelist'] == PRODUCT['namelist']
    assert fresh['self_collection']['kernel_ratio'] == pytest.approx(PRODUCT['self_collection']['kernel_ratio'], rel=1e-12)
