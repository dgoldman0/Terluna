import json

import numpy as np

from shared import constants as K
from shared.provenance import constants_changed
from . import integrated_ledger as il

PRODUCT = json.loads(il.OUT.read_text())


def test_product_is_current_with_its_code_and_the_products_it_reads():
    assert PRODUCT['schema'] == 'terluna.research.integrated-ledger/1'
    producer = PRODUCT['producer']
    assert sorted(producer['files']) == sorted(il.FILES)
    for name, value in producer['files'].items():
        assert il.digest(il.ROOT/name) == value, name
    for name, path in il.PRODUCTS.items():
        assert producer['inputs'][f'results/{name}.json'] == il.digest(path), name
    assert not constants_changed(producer['constants'])


def test_every_gate_judges_both_candidates_with_evidence():
    assert len(PRODUCT['gates']) >= 8
    for gate in PRODUCT['gates']:
        for candidate in ('held_zoned', 'ring_fleet'):
            assert gate[candidate]['state'] in il.STATES and gate[candidate]['evidence']


def test_the_held_screen_figures_are_the_zoned_products():
    zoned = {c['label']: c for c in json.loads(il.PRODUCTS['zoned_aperture'].read_text())['cases']}
    case = zoned[json.loads(il.PRODUCTS['exhaust_isolation'].read_text())['design']['case']]
    held = PRODUCT['held_zoned']
    assert held['propellant_kg_s'] == case['propellant_kg_s']
    assert held['holding_power_mean_TW'] == case['mean_power_TW']
    assert np.isclose(held['propellant_per_held_kg_per_year'], held['propellant_Gt_per_year']/held['held_mass_Gt'])
    # A held kilogram burns about its own mass of propellant each year.
    assert 1. < held['propellant_per_held_kg_per_year'] < 1.2
    assert held['s6_propellant_over_unshielded_loss'][0] > 1.


def test_a_ring_screen_counts_whole_rings():
    for row in PRODUCT['ring_fleet']['by_radius']:
        assert row['tiles'] == row['rings']*row['tiles_per_ring']
        assert 0. < row['window_ring_share'] < 1.
        assert row['tiles_over_held_panels'] > 5.
        assert np.isclose(row['optical_mass_Gt'],
                          row['tiles']*il.REQUIREMENTS['tile_side_m']**2*row['mean_sigma_g_m2']*1e-3/1e12)


def test_collector_areas_follow_from_demand_and_conversion():
    electricity = PRODUCT['electricity']
    for case in electricity['cases']:
        area = case['held']['generation_TW']*1e12/(case['conversion']*K.SOLAR_CONSTANT)
        assert np.isclose(case['held']['collector_km2'], area/1e6)
    assert np.isclose(electricity['dimmer_optical_TW'], 586., rtol=.01)


def test_nesting_puts_annulus_tiles_at_the_window_tiles_sail_loading():
    annulus = PRODUCT['annulus']
    assert PRODUCT['producer']['inputs']['protection/spectra/annulus_film.json'] == il.digest(il.ANNULUS)
    for film in annulus['films'].values():
        assert np.isclose(film['matched_tile_g_m2'], 26*film['pressure_coefficient']/annulus['window_pressure_coefficient'])
        assert film['areal_mass_g_m2'] < film['matched_tile_g_m2'] < 16.


def test_the_nestable_fleet_carries_the_planes_oversizing():
    row = {r['radius_km']: r for r in PRODUCT['ring_fleet']['by_radius']}[20000.]
    motion = row['plane_motion']
    assert 1.15 < motion['ring_factor'] < 1.35 < motion['ring_factor_without_steering']
    for mass in row['matched_annulus'].values():
        assert mass['oversized_without_steering_Gt'] > mass['oversized_optical_mass_Gt'] > mass['optical_mass_Gt'] > \
            row['optical_mass_Gt']
