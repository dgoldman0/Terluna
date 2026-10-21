"""The breakup product binds its code and inputs, reproduces from them, and carries Earth's rain limits to the Moon by
the similarity it states: the capillary length at equal air density, the air's density beyond it."""
import hashlib
import json
import math

import numpy as np
import pytest

from shared.constants import (CM1_GRAVITY, CM1_RAIN_BREAKUP_CUTOFF, CM1_RAIN_BREAKUP_RATE, CM1_RAIN_FALL_CAP,
                              CM1_RAIN_MEAN_DIAMETER_CAP, MOON_SURFACE_GRAVITY, STANDARD_GRAVITY)
from shared.provenance import constants_changed
from atmosphere.rain import breakup, product as rain, targets as climate_targets

PRODUCT = json.loads(breakup.OUT.read_text())
DATA, BINDING = climate_targets.load()
AIR = rain.target_air(climate_targets.read(DATA))


@pytest.fixture(scope='module')
def resort():
    cols = rain.columns(AIR['resort'])
    return dict(moon=breakup.Speeds(cols['moon']), earth=breakup.Speeds(cols['earth']))


def test_header_and_bindings():
    assert PRODUCT['schema'] == breakup.SCHEMA
    assert PRODUCT['evidence'] and PRODUCT['reading_rule'] and PRODUCT['assumptions']
    producer = PRODUCT['producer']
    assert set(producer['files']) == set(breakup.CODE)
    for kind, entries in BINDING.items():
        assert producer[kind] == entries, kind
    for path, digest in {**producer['files'], **producer.get('inputs', {})}.items():
        assert hashlib.sha256((breakup.ROOT / path).read_bytes()).hexdigest() == digest, path
    assert not constants_changed(producer['constants'])
    earth = PRODUCT['earth_limits']
    assert earth['breakup_cutoff_mm'] == CM1_RAIN_BREAKUP_CUTOFF * 1e3 and earth['breakup_rate_per_m'] == CM1_RAIN_BREAKUP_RATE
    assert earth['equilibrium_mm'] == pytest.approx(0.3 + 1e3 * math.log(2) / 2300)
    assert earth['mean_diameter_cap_mm'] == CM1_RAIN_MEAN_DIAMETER_CAP * 1e3 and earth['fall_speed_cap_m_s'] == CM1_RAIN_FALL_CAP


def test_earth_against_itself_is_scale_one(resort):
    earth = resort['earth']
    large, small, weights = breakup.breakup_weights(earth, breakup.EQUILIBRIUM_M)
    keep = weights['straub'] > 1e-9 * weights['straub'].max()
    s, reached = breakup.pair_scale(earth, earth, large[keep], small[keep])
    assert reached.all() and np.allclose(s, 1.0, rtol=1e-9)
    assert breakup.scale(earth, earth, breakup.EQUILIBRIUM_M)['scale'] == pytest.approx(1.0, rel=1e-9)


def test_both_coalescence_laws_coalesce_more_in_gentler_lunar_collisions(resort):
    """Over Low and List's own sizes, a pair meets more gently on the Moon than on Earth, and both laws let it
    coalesce at least as often there."""
    earth, moon = resort['earth'], resort['moon']
    large, small = np.meshgrid(np.linspace(1.2e-3, 4.6e-3, 8), np.linspace(0.4e-3, 1.0e-3, 4))
    large, small = large.ravel(), small.ravel()
    for world in (earth, moon):
        assert np.all(breakup.collision(world, large, small)[3] > 0)
    assert np.all(breakup.collision(moon, large, small)[3] < breakup.collision(earth, large, small)[3])
    assert np.all(breakup.straub(breakup.collision(moon, large, small)[3])
                  > breakup.straub(breakup.collision(earth, large, small)[3]))
    assert np.all(breakup.low_list(moon, large, small)[0] >= breakup.low_list(earth, large, small)[0] - 1e-12)


def test_the_scales_follow_the_capillary_length_and_the_air(resort):
    surface = PRODUCT['scale']['resort']['surface']
    expected = math.sqrt(STANDARD_GRAVITY * (breakup.drops.WATER_DENSITY - resort['earth'].density)
                         / (MOON_SURFACE_GRAVITY * (breakup.drops.WATER_DENSITY - resort['moon'].density)))
    assert surface['capillary'] == pytest.approx(expected, rel=1e-3)
    assert surface['stability_screen_mm']['world'] / surface['stability_screen_mm']['earth'] == pytest.approx(
        surface['capillary'], rel=2e-3)
    # At Earth's air density the equilibrium's scale is the capillary ratio, at the Moon's gravity and between.
    for row in PRODUCT['general']:
        if row['pressure_atm'] == 1.0:
            assert row['equilibrium']['over_capillary'] == pytest.approx(1.0, abs=0.01)
            assert row['cutoff']['over_capillary'] == pytest.approx(1.0, abs=0.05)
    # Denser air raises it as the density to a power near 0.7 at every gravity tried.
    alphas = [row['equilibrium']['alpha'] for row in PRODUCT['general'] if row['pressure_atm'] == 1.2]
    assert len(alphas) == 3 and all(0.6 < a < 0.8 for a in alphas) and max(alphas) - min(alphas) < 0.05
    # The targets and heights agree on the lunar scale within 1%.
    scales = [rows[h]['equilibrium']['scale'] for rows in PRODUCT['scale'].values() for h in ('surface', 'cloud_base')]
    assert max(scales) / min(scales) < 1.01
    assert all(rows[h]['cutoff']['scale'] < rows[h]['equilibrium']['scale'] for rows in PRODUCT['scale'].values()
               for h in ('surface', 'cloud_base'))


def test_the_literal_low_list_reading_breaks_lunar_drops_smaller():
    for rows in PRODUCT['scale'].values():
        for row in rows.values():
            for spectrum in ('cutoff', 'equilibrium'):
                literal = row[spectrum]['low_list_literal']
                assert 1.5 < literal['scale'] < 2.0 < row[spectrum]['scale']
                assert literal['breakup_share_covered'] > 0.6


def test_values_reproduce_from_the_code(resort):
    row = breakup.scale(resort['moon'], resort['earth'], breakup.EQUILIBRIUM_M)
    assert row == pytest.approx({k: v for k, v in PRODUCT['scale']['resort']['surface']['equilibrium'].items()
                                 if k in row}, rel=1e-9)


def test_the_lunar_limits():
    lunar = PRODUCT['lunar_limits']
    assert lunar['scale_capillary'] == math.sqrt(CM1_GRAVITY / MOON_SURFACE_GRAVITY)
    assert lunar['alpha'] == PRODUCT['scale']['resort']['surface']['equilibrium']['alpha']
    assert lunar['scale_collisional'] == pytest.approx(lunar['scale_capillary'] * 1.2**round(lunar['alpha'], 3))
    moon = lunar['moon']
    assert moon['breakup_cutoff_mm'] == pytest.approx(0.3 * lunar['scale_collisional'])
    assert moon['breakup_rate_per_m'] == pytest.approx(2300.0 / lunar['scale_collisional'])
    assert moon['equilibrium_mm'] == pytest.approx(lunar['scale_collisional'] * PRODUCT['earth_limits']['equilibrium_mm'])
    assert moon['mean_diameter_cap_mm'] == pytest.approx(2.8 * lunar['scale_capillary'])
    assert moon['fall_speed_cap_m_s'] == pytest.approx(9.1 * (MOON_SURFACE_GRAVITY / CM1_GRAVITY)**0.25)
    assert lunar['namelist']['brkup_gscale'] == 1 and lunar['namelist']['brkup_rhoexp'] == round(lunar['alpha'], 3)
    # The largest drop's speed follows the quarter power of gravity in the same air within 2%; at the lunar equilibrium
    # Morrison's mass-weighted speed passes the cap the core patch sets and stays under the one the similarity gives.
    fall = lunar['fall_speed']
    assert fall['gravity_ratio_in_the_same_air'] == pytest.approx(fall['similarity_quarter_power'], rel=0.02)
    at = fall['at_2m_air']
    assert at['terluna_cap_m_s'] < at['solver_mass_weighted_speed_m_s'] < at['similarity_cap_m_s']
    assert at['terluna_cap_m_s'] < at['morrison_mass_weighted_speed_m_s'] < at['similarity_cap_m_s']
