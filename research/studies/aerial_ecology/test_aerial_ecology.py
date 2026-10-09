"""Conservation and limiting-case checks for the aerial requirements screen."""
import importlib
import json
import math

import pytest

from shared.constants import MOON_RADIUS, SYNODIC_MONTH_DAYS
from shared.provenance import constants_changed, constants_used

model = importlib.import_module('research.studies.aerial_ecology.run')


def test_growth_compensates_losses_only_above_its_required_duty():
    threshold = model.growth_budget(3.6, 0.0, 30)['minimum_favourable_fraction']
    at_threshold = model.growth_budget(3.6, threshold, 30)
    assert at_threshold['net_growth_per_day'] == pytest.approx(0, abs=1e-14)
    assert at_threshold['multiplier_per_lunar_cycle'] == pytest.approx(1)
    assert model.growth_budget(3.6, threshold / 2, 30)['net_growth_per_day'] < 0
    # Even continuous favourable conditions cannot offset this removal rate.
    assert model.growth_budget(19.5, 1, 3)['minimum_favourable_fraction'] > 1


def test_no_growth_during_unfavourable_cycle():
    case = model.growth_budget(3.6, 0, 30, mortality_per_day=0)
    assert case['multiplier_per_lunar_cycle'] == pytest.approx(math.exp(-SYNODIC_MONTH_DAYS / 30))


def test_feeding_counts_assimilation_and_metabolic_cost():
    case = model.feeding(1, 10, 0.1)
    # One square metre sweeps 864000 m3/day; half of its 864 g is captured.
    assert case['captured_g_day_m2'] == pytest.approx(432)
    assert case['gross_food_w_m2'] == pytest.approx(90)
    assert case['assimilated_w_m2'] == pytest.approx(54)
    assert case['drag_mechanical_w_m2'] == pytest.approx(60)
    assert case['propulsion_metabolic_w_m2'] == pytest.approx(240)
    assert case['powered_net_w_m2'] == pytest.approx(-187)
    assert case['passive_net_w_m2'] == pytest.approx(53)
    balance = model.feeding(case['powered_break_even_mg_m3'], 10, 0.1)
    assert balance['powered_net_w_m2'] == pytest.approx(0, abs=1e-12)


def test_comoving_collector_has_no_sweep_energy():
    case = model.feeding(1, 0, 0.1)
    assert case['captured_g_day_m2'] == 0
    assert case['powered_break_even_mg_m3'] is None
    assert case['passive_net_w_m2'] == -1


def test_optical_depth_uses_particle_cross_section_and_compact_mass():
    # A 10 um sphere has mass rho*pi*d^3/6 and cross-section pi*d^2/4.
    diameter = 10e-6
    mass_each = 1000 * math.pi * diameter**3 / 6
    count_m3 = 0.1e-6 / mass_each
    expected = count_m3 * (2 * math.pi * diameter**2 / 4) * 1000
    assert model.optical_depth(0.1, 1000, 10) == pytest.approx(expected)
    assert expected == pytest.approx(0.03)


def test_floater_wet_mass_and_spherical_skin_exhaust_allowed_lift():
    case = model.floater_mass(1.1, 10, water_fraction=0.9, skin_kg_m2=1, usable_fraction=0.5)
    radius = case['minimum_radius_m']
    area = math.pi * radius**2
    allocated_lift = 0.5 * 1.1 * 4 * math.pi * radius**3 / 3
    assert case['wet_biomass_kg_m2'] == pytest.approx(100)
    assert case['water_kg_m2'] == pytest.approx(90)
    assert allocated_lift == pytest.approx(area * 100 + 4 * area)


def test_production_uses_surface_area_and_does_not_count_harvest_twice():
    case = model.production(0.01, 1000, 1)
    assert case['projected_area_m2'] == pytest.approx(4 * math.pi * MOON_RADIUS**2 * 0.01)
    assert 3.7e11 < case['projected_area_m2'] < 3.9e11
    assert 1.8e8 < case['food_energy_people'] < 1.9e8
    assert case['dry_stock_kg'] < case['dry_production_kg_year']
    reduced = model.production(0.01, 1000, 1, edible=0.5, processing=0.8)
    assert reduced['food_energy_people'] == pytest.approx(0.4 * case['food_energy_people'])
    assert reduced['harvested_dry_kg_year'] == case['harvested_dry_kg_year']


def test_runner_records_inputs_without_predicting_species_or_new_capacity(tmp_path, monkeypatch):
    monkeypatch.setattr(model, 'OUT', tmp_path / 'aerial.json')
    product = model.run()
    assert json.loads(model.OUT.read_text()) == product
    assert product['biodiversity']['species_count'] is None
    assert product['baseline']['p_fallout_analogue_kg_year'][0] < product['baseline']['p_fallout_analogue_kg_year'][1]
    assert product['central']['p_stock_kg'] > 0
    for path, digest in product['producer']['inputs'].items():
        assert model.digest(model.ROOT / path) == digest
    for case in product['trophic_production_kg_year']:
        assert case['secondary_consumer'] < case['primary_consumer']
        assert case['primary_consumer'] + product['central']['harvested_dry_kg_year'] < product['central']['dry_production_kg_year']


def test_saved_product_pins_current_dependencies_and_reproduces(tmp_path, monkeypatch):
    saved = json.loads((model.HERE / 'results/aerial_ecology.json').read_text())
    provenance = saved['producer']
    for group in ('files', 'inputs'):
        for path, recorded_digest in provenance[group].items():
            assert model.digest(model.ROOT / path) == recorded_digest, path
    assert not constants_changed(provenance['constants'])
    assert provenance['constants'] == constants_used([
        model.ROOT / path for path in provenance['files']
    ])
    monkeypatch.setattr(model, 'OUT', tmp_path / 'regenerated_aerial.json')
    model.run()
    assert json.loads(model.OUT.read_text()) == saved
