"""Checks of the nitrogen estimate (nitrogen.py): the published yields as the note gives them and each route's
arithmetic."""
import json

import numpy as np
import pytest

from research.studies.atmospheric_electricity import nitrogen as nt


def test_the_per_metre_yield_follows_wang_s_fit_and_earth_s_source_its_area():
    assert float(nt.per_metre(101325.0)) == pytest.approx(1.66e21, rel=0.01)       # the note's 1 atm value
    assert float(nt.per_metre(0.5 * 101325.0) / nt.per_metre(101325.0)) == pytest.approx(0.60, abs=0.01)
    # 5 (3-7) Tg N a year over Earth's 5.1e8 km2
    assert nt.EARTH_KM2 == pytest.approx(5.10e8, rel=0.005)
    assert nt.EARTH_KG_N_KM2_YR['central'] == pytest.approx(9.8, rel=0.01)
    assert nt.EARTH_KG_N_KM2_YR['range'] == pytest.approx((5.9, 13.7), rel=0.01)


def test_each_route_counts_its_own_way_and_spreads_the_nitrogen_over_the_box_and_the_year():
    energy = np.array([1.0e11, 3.0e11])
    out = nt.estimate(energy, np.array([30000.0, 31000.0]), np.array([40000.0, 41000.0]), np.array([1.0e8, 4.0e8]),
                      np.array([False, True]), np.array([0.0, 50000.0]), np.array([1.2e5, 0.4e5]), 1000.0, 365.25)
    assert out['per_flash']['total_mol'] == 500.0
    assert out['per_joule']['total_mol'] == pytest.approx(4.0e11 * 9.5e16 / 6.02214076e23)
    # channels: 10 km + 10 km for the flash in cloud; 41 km to the ground + 20 km for the strike
    assert out['channel_km']['total'] == pytest.approx(81.0)
    p_mid = np.interp([35000.0, 20500.0], [0.0, 50000.0], [1.2e5, 0.4e5])
    expect = (20000.0 * nt.per_metre(p_mid[0]) + 61000.0 * nt.per_metre(p_mid[1])) / 6.02214076e23
    assert out['per_metre']['total_mol'] == pytest.approx(expect)
    # over 1000 km2 and one year, in the design air's oxygen
    assert out['per_flash']['kg_n_per_km2_yr'] == pytest.approx(500.0 / 1000.0 * 14.007e-3 * 0.93)
    assert out['per_flash']['share_of_earth'] == pytest.approx(out['per_flash']['kg_n_per_km2_yr'] / 9.8, rel=0.01)
    # two flashes over 1000 km2 in a year, against Earth's 44 a second: 2.7 per km2 a year
    assert out['flashes_per_km2_yr']['box'] == pytest.approx(0.002)
    assert out['flashes_per_km2_yr']['earth'] == pytest.approx(2.72, rel=0.01)


def test_the_flash_sizes_come_from_the_products_and_set_the_benchmark_s_grids_side_by_side(tmp_path, monkeypatch):
    def product(name, kinds):
        stats = lambda median, mean: dict(p10=0.5 * median, median=median, p90=2.0 * median, max=5.0 * median, mean=mean)
        lightning = {kind: dict(count=count, energy_dissipated_j=stats(energy, 1.5 * energy),
                                positive_c=stats(charge, charge), negative_c=stats(charge if kind == 'in_cloud' else
                                                                                   2.0 * charge, charge))
                     for kind, (count, energy, charge) in kinds.items()}
        (tmp_path / f'elec_{name}.json').write_text(json.dumps(dict(lightning=lightning)))
    product('run', {'in_cloud': (40, 8.0e10, 100.0), 'negative_to_ground': (10, 6.0e10, 75.0)})
    product('fine', {'in_cloud': (1000, 1.0e9, 10.0)})
    product('coarse', {'in_cloud': (800, 2.0e9, 15.0)})
    monkeypatch.setattr(nt, 'CRM_RESULTS', tmp_path)
    out = nt.flash_sizes('run', ('fine', 'coarse'))
    assert out['in_cloud']['energy_j']['median'] == 8.0e10 and out['in_cloud']['charge_c_median'] == 100.0
    assert out['negative_to_ground']['charge_c_median'] == 150.0                 # a ground strike's negative charge
    assert out['median_energy_over_earth_benchmark'] == dict(in_cloud=80.0, negative_to_ground=60.0)
    ratio = out['earth_benchmark']['coarse_over_fine']
    assert ratio['flashes'] == 0.8 and ratio['median_energy_j'] == 2.0
    assert ratio['total_energy_j'] == pytest.approx(800 * 3.0e9 / (1000 * 1.5e9))
