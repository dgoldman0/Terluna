"""Independent geometric, pressure, selection and transport controls."""
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import minimum_spanning_tree
from scipy.spatial.distance import cdist

from shared.constants import GAS_CONSTANT, MOON_RADIUS, STANDARD_ATMOSPHERE
from shared.provenance import constants_changed
from research.studies.vacuum_crater_buildout import run as model


@pytest.fixture
def scenarios():
    return json.loads((model.HERE / 'scenarios.json').read_text())


def rows_at(coordinates, diameters=200):
    rows = np.zeros(len(coordinates), dtype=model.DTYPE)
    rows['id'] = [f'test_{i:04d}'.encode() for i in range(len(rows))]
    rows['lat'], rows['lon'] = np.asarray(coordinates).T
    rows['diameter'] = diameters
    rows['depth'] = rows['diameter'] * 55 / 750
    return rows


def test_common_denominator_downward_rounding():
    assert model.quota(8886558, .025) == 222163
    assert model.quota(8886558, .001) == 8886
    assert model.quota(40000, .025) == 1000
    assert model.quota(999, .001) == 0
    with pytest.raises(ValueError):
        model.quota(10, 1.1)


def test_longitude_wrap_and_variable_estate_sizes():
    rows = rows_at([(0, 179.99), (0, -179.99), (0, 179.98), (0, 179.8)], [200, 500, 200, 400])
    xyz = model.positions(rows)
    assert model.arc_length(np.linalg.norm(xyz[0] - xyz[1])) == pytest.approx(MOON_RADIUS * math.radians(.02))
    chosen = model.separated_sites(xyz, rows['diameter'], np.arange(4), 4, 2)
    assert chosen.tolist() == [0, 3]
    for a, b in itertools.combinations(chosen, 2):
        distance = model.arc_length(np.linalg.norm(xyz[a] - xyz[b]))
        assert distance >= rows['diameter'][a] + rows['diameter'][b]
    assert model.separated_sites(xyz[:0], rows['diameter'][:0], [], 5, 2).size == 0


def test_rank_and_selection_independent_of_input_order(scenarios):
    rows = rows_at([(0, x) for x in np.linspace(-1, 1, 60)])
    rows['depth'][::3] *= 1.1
    def selected_ids(r):
        _, _, order = model.ranks(r, scenarios['selection'])
        chosen = model.separated_sites(model.positions(r), r['diameter'], order, 20, 2)
        return r['id'][chosen].tolist()
    permutation = np.random.default_rng(19).permutation(len(rows))
    assert selected_ids(rows) == selected_ids(rows[permutation])


def test_clustered_subset_has_selected_seeds_and_reports_shortfall(scenarios):
    rows = rows_at([(0, x) for x in np.linspace(0, .5, 50)] + [(30, 20)], 300)
    xyz = model.positions(rows)
    bins, tie, _ = model.ranks(rows, scenarios['selection'])
    final = np.arange(len(rows))
    chosen, groups = model.clustered_initial(rows, xyz, final, 35, bins, tie, scenarios['selection'], (238, 336))
    assert len(chosen) == len(set(chosen)) == 35
    assert set(chosen) <= set(final)
    assert sum(g['count'] for g in groups) == len(chosen)
    assert all(g['max_seed_distance_km'] <= 25.00000001 for g in groups)
    assert all(g['seed_id'].encode() in rows['id'][chosen] for g in groups)
    chosen, _ = model.clustered_initial(rows, xyz, final, 60, bins, tie, scenarios['selection'], (238, 336))
    assert len(chosen) == 50  # The isolated 51st crater cannot form a qualifying group.


def test_spherical_tree_matches_complete_graph():
    rng = np.random.default_rng(5)
    rows = rows_at(np.column_stack((rng.uniform(-70, 70, 30), rng.uniform(-180, 180, 30))))
    xyz = model.positions(rows)
    edges, lengths = model.spanning_tree(xyz)
    exact = minimum_spanning_tree(csr_matrix(model.arc_length(cdist(xyz, xyz))))
    assert len(edges) == len(rows) - 1
    assert lengths.sum() == pytest.approx(exact.sum(), rel=1e-12)


def test_route_can_cross_future_water_with_dry_endpoints(scenarios):
    rows = rows_at([(0, -5), (0, 5)])
    xyz = model.positions(rows)
    edges, lengths = model.spanning_tree(xyz)
    water = np.zeros((180, 360), np.uint8)
    water[90, 0] = 1
    flags = model.route_water_intersections(xyz, edges, lengths, water, 1000)
    assert flags.tolist() == [True]
    result = model.network(xyz, water, scenarios['transport'])
    assert result['future_water_flagged_links'] == 1
    assert result['components']['25km']['isolated_sites'] == 2


def test_cap_and_bowl_volume_against_independent_slices(scenarios):
    c = scenarios['enclosure']
    d, depth = 300., 22.
    output = model.enclosure(d, depth, STANDARD_ATMOSPHERE, c)
    h = c['cap_height_over_diameter'] * d
    radius = ((d / 2) ** 2 + h ** 2) / (2 * h)
    n = 20000
    z = (np.arange(n) + .5) * h / n
    cap = np.sum(math.pi * (radius ** 2 - (radius - h + z) ** 2)) * h / n
    z = (np.arange(n) + .5) * depth / n
    bowl = np.sum(math.pi * (d / 2) ** 2 * z / depth) * depth / n
    assert output['cap_volume_m3'] == pytest.approx(cap, rel=1e-8)
    assert output['bowl_volume_m3'] == pytest.approx(bowl, rel=1e-12)
    inferred_moles = output['gross_gas_mass_kg'] / c['gas_molar_mass_kg_mol']
    assert inferred_moles * GAS_CONSTANT * c['reference_temperature_K'] == pytest.approx(STANDARD_ATMOSPHERE * (cap + bowl), rel=1e-8)


def test_pressure_reaction_and_hemisphere_limit(scenarios):
    c = {**scenarios['enclosure'], 'cap_height_over_diameter': .5}
    result = model.enclosure(300, 0, STANDARD_ATMOSPHERE, c)
    assert result['cap_volume_m3'] == pytest.approx(2 * math.pi * 150 ** 3 / 3)
    assert result['membrane_force_N_m'] == pytest.approx(STANDARD_ATMOSPHERE * 150 / 2)
    assert result['mean_vertical_rim_reaction_N_m'] * math.pi * 300 == pytest.approx(result['uplift_N'])
    assert result['uplift_N'] == pytest.approx(7.162242202e9, rel=1e-9)
    empty = model.enclosure(300, 0, 0, c)
    assert empty['uplift_N'] == empty['gross_gas_mass_kg'] == 0
    assert empty['shielding'][1]['weight_N'] > 0


def test_similar_smaller_caps_reduce_pressure_material_but_not_shield_mass(scenarios):
    c = scenarios['enclosure']
    big = model.enclosure(300, 22, STANDARD_ATMOSPHERE, c)
    small = model.enclosure(100, 22 / 3, STANDARD_ATMOSPHERE, c)
    assert small['equivalent_pressure_material_kg'] * 9 == pytest.approx(big['equivalent_pressure_material_kg'] / 3)
    assert small['shielding'][1]['mass_kg'] * 9 == pytest.approx(big['shielding'][1]['mass_kg'])
    higher = model.enclosure(300, 22, STANDARD_ATMOSPHERE * 1.2, c)
    assert higher['uplift_N'] == pytest.approx(big['uplift_N'] * 1.2)


def test_trip_kinematics_short_and_long():
    assert model.travel_time(.01, 150, .2) == pytest.approx(2 * math.sqrt(10 / .2))
    assert model.travel_time(250, 150, .2) == pytest.approx(250000 / (150 / 3.6) + (150 / 3.6) / .2)
    assert model.travel_time(0, 10, .2) == 0
    with pytest.raises(ValueError):
        model.travel_time(1, 0, .2)


def test_lifecycle_and_population_are_separate_quantities(scenarios):
    rows = rows_at([(0, 0), (0, 1)], [300, 300])
    result = model.engineering(rows, np.array([0, 1]), 1000, scenarios)
    assert result['replacement_counts']['200']['20'] == 9
    assert result['whole_crater_initial_scenario']['gross_air_mass_tonnes'] > 0
    assert result['assumed_water_loop_makeup_tonnes_per_day'] == pytest.approx(2)
    assert result['half_synodic_month_hours'] == pytest.approx(354.36708)
    assert result['full_architectural_occupancy_service_sensitivity'][0]['electric_MW'] == 2


def test_changed_input_rejected(tmp_path):
    path = tmp_path / 'product.json'
    path.write_text('{"schema":"example/1"}')
    digest = model.sha256(path)
    path.write_text('{"schema":"example/1","changed":true}')
    with pytest.raises(ValueError, match='hash mismatch'):
        model.read_product(path, 'example/1', digest)


def test_published_product_binds_inputs_and_constants():
    product = json.loads((model.HERE / 'results/buildout.json').read_text())
    assert product['schema'] == model.SCHEMA
    for path, expected in product['producer']['files'].items():
        assert model.sha256(model.ROOT / path) == expected
    for name, record in product['producer']['inputs'].items():
        path = Path(record['path'])
        if path.is_file():
            assert model.sha256(path) == record['sha256'], name
    assert not constants_changed(product['producer']['constants'])
    assert product['initial_clustered']['selected'] <= product['programme']['initial_sites']


def test_food_energy_normalizes_by_all_people_and_adds_services(scenarios):
    food = {'systems': {'power_made_food': {'total_kw_per_resident': [1, 2]}},
            'mixes': {'fully_lamp_lit': {'cases': {'4': {'kept_50': {
                'power_tw': {'total': [20, 40]}, 'people_billion': {'lunar': [4, 4]}}}}}}}
    result = model.food_comparison(1000000, food, scenarios['operation'])
    assert result['fully_lamp_lit_case']['additional_food_kW_per_resident'] == [5, 10]
    assert result['fully_lamp_lit_case']['services_plus_food_GW'] == [10, 15]
    assert result['power_made_development_case']['services_plus_food_GW'] == [6, 7]


def test_stronger_material_keeps_external_loads(scenarios):
    rows = rows_at([(0, 0)], 300)
    result = model.engineering(rows, np.array([0]), 795, scenarios)
    cases = result['material_improvement_cases']
    assert [x['thickness_and_pressure_material_factor'] for x in cases] == pytest.approx([1, 2/3, .5])
    assert cases[-1]['initial_equivalent_pressure_material_tonnes'] == pytest.approx(cases[0]['initial_equivalent_pressure_material_tonnes'] / 2)
    shell = model.enclosure(300, 22, STANDARD_ATMOSPHERE, scenarios['enclosure'])
    assert shell['shielding'][1]['pressure_to_balance_weight_Pa'] == pytest.approx(9043.65070469)
