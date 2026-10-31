"""Check which carbon is available and keep harvest and nutrient accounts separate."""
import json

import pytest

from provisioning import aerophytes as A


def reference(cls='lesser'):
    products, _ = A.load()
    diameter = {'lesser': 100, 'sky_reef': 60, 'round_giant': 1000}[cls]
    return next(r for r in products['sky_ecology']['food_web']['production']
                if r['cls'] == cls and r['diameter_m'] == diameter
                and r['gpp_c_kg_m2_year'] == 2 and r['gas_route'] == 'photolysis')


def test_no_crop_is_taken_from_an_underfunded_host():
    for cls in ('lesser', 'sky_reef', 'round_giant'):
        row = A.allocate(reference(cls), 0, .003, .5)
        assert not row['carbon_budget_nonnegative']
        assert row['allocation_t_dry_per_projected_ha_year'] is None
        assert row['product_retained_c_kg_m2_year'] is None


def test_reference_surplus_matches_source_after_population_replacement():
    for cls in ('lesser', 'sky_reef', 'round_giant'):
        body = reference(cls)
        for old in body['replacement']:
            row = A.allocate(body, 1, old['hazard_per_year'], .25)
            assert row['unallocated_before_harvest_c_kg_m2_year'] == pytest.approx(old['surplus_left_c'])
            assert abs(row['carbon_closure_residual_kg_m2_year']) < 1e-12


def test_crop_has_construction_respiration_and_competes_with_other_uses():
    small = A.allocate(reference(), 1, .003, .25)
    large = A.allocate(reference(), 1, .003, .5)
    assert large['product_retained_c_kg_m2_year'] == pytest.approx(2 * small['product_retained_c_kg_m2_year'])
    assert large['remaining_for_other_uses_c_kg_m2_year'] < small['remaining_for_other_uses_c_kg_m2_year']
    assert large['product_construction_respiration_c_kg_m2_year'] > 0
    assert large['product_assimilate_c_kg_m2_year'] == pytest.approx(
        large['product_retained_c_kg_m2_year'] + large['product_construction_respiration_c_kg_m2_year'])
    assert large['existing_tenant_sugar_c_kg_m2_year'] == small['existing_tenant_sugar_c_kg_m2_year'] == .1


def test_dim_light_and_higher_mortality_reduce_allocatable_crop():
    original = A.allocate(reference(), 1, .001, .25)
    dim = A.allocate(reference(), .9, .001, .25)
    deaths = A.allocate(reference(), 1, .005, .25)
    key = 'allocation_t_dry_per_projected_ha_year'
    assert dim[key] < original[key]
    assert deaths[key] < original[key]


def test_phosphorus_return_is_circulation_not_zero_requirement():
    p = A.harvest_phosphorus(4, [1, 3], .99)
    assert p['exported_kg_p_ha_year'] == [4, 12]
    assert p['returned_upward_kg_p_ha_year'] == pytest.approx([3.96, 11.88])
    for e, r, u in zip(*p.values()):
        assert e == pytest.approx(r + u)
    assert A.harvest_phosphorus(4, [1, 3], 1)['unreturned_kg_p_ha_year'] == [0, 0]


def test_invalid_fractions_fail():
    with pytest.raises(ValueError):
        A.allocate(reference(), 1, .001, 1.01)
    with pytest.raises(ValueError):
        A.harvest_phosphorus(1, [1, 3], 1.01)


def test_stored_product_regenerates_from_pinned_inputs():
    assert json.loads(A.OUT.read_text()) == A.result()
    assert any(not r['carbon_budget_nonnegative'] for r in A.result()['cases'])


def test_living_area_cluster_counts_only_its_domatia_payload():
    product = A.result()
    products, settings = A.load()
    dm = next(r for r in products['sky_ecology']['tenants']['domatia']
              if r['diameter_m'] == 60 and r['food_to_tenants_c'] == .1)
    for c in product['inhabited_cluster_comparisons']:
        assert c['conditional_domatia_tenant_allowance_t'] == pytest.approx(c['domatia'] * dm['steady_tenants_t'])
        assert c['supported_floor_area_m2'] is None
        assert c['green_residual_after_routine_support_t_c_year'] == pytest.approx(
            c['green_residual_after_green_replacement_t_c_year'] - c['domatia_routine_support_t_c_year'])
        assert 0 < c['share_of_green_residual_for_routine_domatia_support'] < 1


def test_residential_mass_equivalents_scale_with_cluster_and_payload():
    a, b = A.result()['inhabited_cluster_comparisons']
    assert b['conditional_domatia_tenant_allowance_t'] == pytest.approx(10 * a['conditional_domatia_tenant_allowance_t'])
    assert b['sum_projected_area_ha'] == pytest.approx(10 * a['sum_projected_area_ha'])
    for case in (a, b):
        p = case['payload_equivalent_people_by_district_budget']
        assert p['light20'] == pytest.approx(2.5 * p['central50'])
        assert p['central50'] == pytest.approx(2 * p['heavy100'])
