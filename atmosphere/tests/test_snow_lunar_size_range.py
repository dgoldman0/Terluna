"""The single-crystal size range: water per crystal sets the size, and the band brackets it."""
import json
import math

import pytest

from atmosphere.snow import lunar_size_range as sizes


@pytest.fixture(scope='module')
def product():
    if not sizes.PRODUCT.exists():
        pytest.skip('size-range product not assembled')
    return json.loads(sizes.PRODUCT.read_text())


def test_equal_water_per_crystal_gives_equal_crystals(product):
    by_share = {}
    for row in product['grid']:
        share = row['water_mm_h'] / row['flux_factor']
        by_share.setdefault(round(share, 9), []).append(row['baseline']['diameter_m'])
    repeated = [v for v in by_share.values() if len(v) > 1]
    assert repeated
    for values in repeated:
        assert max(values) == pytest.approx(min(values), rel=1e-6)


def test_band_brackets_every_crystal_and_deposition_stays_within_its_share(product):
    for row in product['grid']:
        assert row['combined_small']['diameter_m'] < row['baseline']['diameter_m'] < row['combined_large']['diameter_m']
        for edge in ('baseline', 'combined_small', 'combined_large'):
            assert row[edge]['deposited_kg'] <= row[edge]['allowed_kg'] * (1 + 1e-9)
    assert len(product['grid']) == len(sizes.WATER_MM_H) * len(sizes.FLUX_FACTORS)


def test_water_per_crystal_is_supply_over_flux(product):
    for row in product['grid']:
        assert row['water_per_crystal_kg'] == pytest.approx(
            row['water_mm_h'] / 3600. / row['number_flux_m2_s'], rel=1e-12)
        assert math.isclose(row['number_flux_m2_s'], product['column']['reference_number_flux_m2_s'] * row['flux_factor'],
                            rel_tol=1e-12)
