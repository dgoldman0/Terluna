"""Nested crater-site scenarios, vacuum services and the transition to open communities.

Inputs are the previous census's data products, the adopted buildout plan and
explicit engineering scenarios. No full structural or routing solver is implied.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from decimal import Decimal, ROUND_FLOOR
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import time

import numpy as np
from scipy.spatial import cKDTree, ConvexHull
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import minimum_spanning_tree, connected_components

from shared.constants import (GAS_CONSTANT, MOON_RADIUS, MOON_SURFACE_GRAVITY,
                              STANDARD_ATMOSPHERE, STEFAN_BOLTZMANN,
                              SYNODIC_MONTH_DAYS, JULIAN_YEAR_DAYS)
from shared.provenance import constants_used

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DEFAULT_WORK = Path('/media/projectspace/terluna-research/vacuum-crater-buildout')
DTYPE = np.dtype([('id', 'S32'), ('lat', 'f8'), ('lon', 'f8'), ('diameter', 'f8'),
                  ('depth', 'f8'), ('rain', 'f8'), ('ground', 'f8')])
SCHEMA = 'terluna.research.vacuum-crater-buildout/1'


def sha256(path):
    with Path(path).open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def read_product(path, expected_schema, digest=None):
    if digest is not None and sha256(path) != digest:
        raise ValueError(f'Input hash mismatch: {path}')
    result = json.loads(Path(path).read_text())
    if result['schema'] != expected_schema:
        raise ValueError(f'Unexpected schema: {path}')
    return result


def quota(count, fraction):
    if count < 0 or not 0 <= fraction <= 1:
        raise ValueError('Nonnegative count and a fraction in [0,1] required')
    return int((Decimal(count) * Decimal(str(fraction))).to_integral_value(rounding=ROUND_FLOOR))


def positions(rows):
    lat, lon = np.deg2rad(rows['lat']), np.deg2rad(rows['lon'])
    return MOON_RADIUS * np.column_stack((np.cos(lat) * np.cos(lon),
                                          np.cos(lat) * np.sin(lon), np.sin(lat)))


def arc_length(chord):
    return 2 * MOON_RADIUS * np.arcsin(np.clip(np.asarray(chord) / (2 * MOON_RADIUS), 0, 1))


def chord_length(arc):
    return 2 * MOON_RADIUS * np.sin(np.asarray(arc) / (2 * MOON_RADIUS))


def read_candidates(path, census, selection):
    if sha256(path) != census['candidates']['sha256']:
        raise ValueError('Candidate export differs from its recorded bytes')
    capacity = census['candidates']['rows']
    data = np.empty(capacity, dtype=DTYPE)
    n, stats = 0, Counter()
    with gzip.open(path, 'rt', newline='') as f:
        for row in csv.DictReader(f):
            stats['input_records'] += 1
            high = row['confidence_level'] == 'high'
            inland = row['adjacent_water_cell'] == '0'
            stats['high_confidence'] += high
            stats['no_adjacent_water'] += inland
            if selection['require_high_confidence'] and not high:
                continue
            if selection['avoid_adjacent_future_water'] and not inland:
                continue
            ident = row['crater_id'].encode('ascii')
            if len(ident) > 32 or n >= capacity:
                raise ValueError('Candidate identity or row count exceeds the recorded format')
            data[n] = (ident, float(row['latitude_today_deg']), float(row['longitude_today_deg']),
                       float(row['diameter_m']), float(row['depth_m']),
                       float(row['rain_mm_day']), float(row['mean_ground_temperature_K']))
            n += 1
    if stats['input_records'] != capacity:
        raise ValueError('Candidate count differs from census')
    rows = data[:n].copy()
    if len(np.unique(rows['id'])) != n:
        raise ValueError('Repeated crater identifier')
    for field in DTYPE.names[1:]:
        if not np.all(np.isfinite(rows[field])):
            raise ValueError(f'Nonfinite candidate field: {field}')
    if np.any((rows['diameter'] < 200) | (rows['diameter'] > 1000) | (rows['depth'] <= 0)):
        raise ValueError('Candidate dimensions outside the programme')
    if np.any(np.abs(rows['lat']) > 90) or np.any(np.abs(rows['lon']) > 180):
        raise ValueError('Invalid candidate coordinates')
    stats['ranking_pool'] = n
    return rows, dict(stats)


def ranks(rows, selection):
    error = np.abs(rows['depth'] / rows['diameter'] - selection['target_depth_diameter_ratio'])
    bins = np.floor(error / selection['geometry_error_bin']).astype(int)
    tie = np.array([int.from_bytes(hashlib.sha256(x).digest()[:8], 'big')
                    for x in rows['id']], dtype=np.uint64)
    return bins, tie, np.lexsort((tie, bins))


def separated_sites(xyz, diameters, order, requested, estate_factor):
    """Greedy circular estates on the sphere; never promises a maximum packing."""
    if requested <= 0 or len(xyz) == 0:
        return np.array([], dtype=int)
    radii = np.asarray(diameters) * estate_factor / 2
    tree = cKDTree(xyz)
    blocked = np.zeros(len(xyz), bool)
    selected = []
    max_radius = float(radii.max())
    for i in order:
        if blocked[i]:
            continue
        selected.append(int(i))
        neighbours = np.asarray(tree.query_ball_point(xyz[i], chord_length(radii[i] + max_radius)), dtype=int)
        distances = arc_length(np.linalg.norm(xyz[neighbours] - xyz[i], axis=1))
        blocked[neighbours[distances < radii[i] + radii[neighbours] - 1e-7]] = True
        blocked[i] = True
        if len(selected) == requested:
            break
    return np.asarray(selected, dtype=int)


def clustered_initial(rows, xyz, final, requested, bins, tie, selection, diameter_bounds):
    """Groups have a selected seed within the stated radius of every group member."""
    eligible = final[(rows['diameter'][final] >= diameter_bounds[0]) & (rows['diameter'][final] <= diameter_bounds[1])]
    if not len(eligible) or not requested:
        return np.array([], dtype=int), []
    points = xyz[eligible]
    tree = cKDTree(points)
    radius = float(chord_length(selection['cluster_radius_km'] * 1000))
    density = tree.query_ball_point(points, radius, return_length=True, workers=1)
    order = np.lexsort((tie[eligible], bins[eligible], -density))
    taken = np.zeros(len(eligible), bool)
    result, groups = [], []
    for seed in order:
        if taken[seed]:
            continue
        near = np.asarray(tree.query_ball_point(points[seed], radius), dtype=int)
        near = near[~taken[near]]
        if len(near) < selection['cluster_min_sites']:
            continue
        neighbours = near[near != seed]
        ids = eligible[neighbours]
        sizes = np.floor(rows['diameter'][ids] / selection['initial_diameter_bin_m']).astype(int)
        neighbours = neighbours[np.lexsort((tie[ids], bins[ids], sizes))]
        group = np.concatenate(([seed], neighbours))[:min(selection['cluster_max_sites'], requested - len(result))]
        taken[group] = True
        result.extend(eligible[group].tolist())
        groups.append(dict(seed_id=rows['id'][eligible[seed]].decode(), count=len(group),
                           max_seed_distance_km=float(arc_length(np.linalg.norm(points[group] - points[seed], axis=1)).max() / 1000)))
        if len(result) == requested:
            break
    return np.array(result, dtype=int), groups


def populations(diameters, design):
    ground = math.pi / 4 * np.asarray(diameters) ** 2 * (1 - design['population_lake_to_settlement_diameter'] ** 2)
    return {name: float(np.sum(ground * p['non_lake_ground_building_share'] * p['mean_storeys'] / p['floor_m2_per_resident']))
            for name, p in design['population_cases'].items()}


def distribution(values):
    if not len(values):
        return None
    q = np.percentile(values, [0, 10, 50, 90, 100])
    return dict(zip(('min', 'p10', 'median', 'p90', 'max'), map(float, q)))


def sample_summary(rows, xyz, indices, design, requested, estate_factor):
    sub = rows[indices]
    nearest = None
    if len(indices) > 1:
        distance = cKDTree(xyz[indices]).query(xyz[indices], k=2, workers=1)[0][:, 1]
        nearest = distribution(arc_length(distance) / 1000)
    return dict(selected=len(indices), requested=requested, shortfall=max(requested - len(indices), 0),
                diameter_m=distribution(sub['diameter']), depth_m=distribution(sub['depth']),
                bands={name: int(np.sum((sub['diameter'] >= lo) & (sub['diameter'] < hi)))
                       for name, lo, hi in [('200_300', 200, 300), ('300_500', 300, 500), ('500_1000', 500, 1000.001)]},
                architectural_population=populations(sub['diameter'], design),
                circular_estate_area_km2=float(np.sum(math.pi / 4 * (estate_factor * sub['diameter']) ** 2) / 1e6),
                nearest_other_site_km=nearest,
                latitude_today_deg=distribution(sub['lat']),
                longitude_today_deg=distribution(sub['lon']))


def spanning_tree(xyz):
    """Spherical Delaunay candidate edges contain a distance-minimum spanning tree.

Chord and great-circle edge lengths have the same order. Thus an MST selected
on chord lengths is also an MST for great-circle weights on this sphere.
"""
    n = len(xyz)
    if n < 2:
        return np.empty((0, 2), int), np.array([], float)
    if n <= 16:
        a, b = np.triu_indices(n, 1)
        edges = np.column_stack((a, b))
    else:
        triangles = ConvexHull(xyz / MOON_RADIUS).simplices
        edges = np.unique(np.sort(np.concatenate((triangles[:, [0, 1]], triangles[:, [1, 2]], triangles[:, [2, 0]])), axis=1), axis=0)
    weights = arc_length(np.linalg.norm(xyz[edges[:, 0]] - xyz[edges[:, 1]], axis=1))
    graph = coo_matrix((np.r_[weights, weights], (np.r_[edges[:, 0], edges[:, 1]], np.r_[edges[:, 1], edges[:, 0]])), shape=(n, n)).tocsr()
    tree = minimum_spanning_tree(graph).tocoo()
    if tree.nnz != n - 1:
        raise ValueError('The geographic spanning tree is disconnected')
    return np.column_stack((tree.row, tree.col)), tree.data


def route_water_intersections(xyz, edges, lengths, water, spacing_m):
    flagged = np.zeros(len(edges), bool)
    for j, ((a, b), length) in enumerate(zip(edges, lengths)):
        count = max(2, math.ceil(length / spacing_m) + 1)
        theta = length / MOON_RADIUS
        t = np.linspace(0, 1, count)
        if theta < 1e-12:
            points = np.repeat(xyz[a][None], count, axis=0)
        else:
            points = (np.sin((1 - t) * theta)[:, None] * xyz[a] + np.sin(t * theta)[:, None] * xyz[b]) / math.sin(theta)
        lat = np.rad2deg(np.arctan2(points[:, 2], np.hypot(points[:, 0], points[:, 1])))
        lon = np.rad2deg(np.arctan2(points[:, 1], points[:, 0]))
        r = np.minimum(((90 - lat) / 180 * water.shape[0]).astype(int), water.shape[0] - 1)
        c = ((lon % 360) / 360 * water.shape[1]).astype(int) % water.shape[1]
        flagged[j] = np.any(water[r, c] != 0)
    return flagged


def network(xyz, water, transport):
    edges, lengths = spanning_tree(xyz)
    flags = route_water_intersections(xyz, edges, lengths, water, transport['water_sample_spacing_m'])
    components = {}
    for limit in transport['connection_radii_km']:
        for dry_only in (False, True):
            keep = (lengths <= limit * 1000) & (~flags if dry_only else True)
            pairs = edges[keep]
            graph = coo_matrix((np.ones(2 * len(pairs)), (np.r_[pairs[:, 0], pairs[:, 1]], np.r_[pairs[:, 1], pairs[:, 0]])), shape=(len(xyz), len(xyz))).tocsr()
            count, labels = connected_components(graph, directed=False)
            sizes = np.bincount(labels)
            components[f'{limit:g}km' + ('_without_water_flagged_tree_links' if dry_only else '')] = dict(components=int(count), isolated_sites=int(np.sum(sizes == 1)), largest_component_sites=int(sizes.max()) if len(sizes) else 0)
    return dict(geographic_minimum_tree_km=float(lengths.sum() / 1000), links=len(edges),
                link_length_km=distribution(lengths / 1000),
                future_water_flagged_links=int(flags.sum()),
                future_water_flagged_link_length_km=float(lengths[flags].sum() / 1000),
                components=components,
                reading_rule='Great-circle lengths, without terrain detours, grades, route capacity or redundancy. Water is sampled along links on the inherited coarse mask; unflagged links remain unverified. Removing flagged tree links gives a conservative connectivity diagnostic; alternative dry links were not searched.')


def enclosure(diameter, depth, pressure_pa, config):
    if diameter <= 0 or depth < 0 or pressure_pa < 0:
        raise ValueError('Positive diameter and nonnegative depth and pressure required')
    a, h = diameter / 2, config['cap_height_over_diameter'] * diameter
    area, roof = math.pi * a * a, math.pi * (a * a + h * h)
    curvature = (a * a + h * h) / (2 * h)
    cap = math.pi * h * (3 * a * a + h * h) / 6
    bowl = area * depth / 2
    tension = pressure_pa * curvature / 2
    thickness = tension / config['material_effective_allowable_Pa']
    gas = pressure_pa * (cap + bowl) * config['gas_molar_mass_kg_mol'] / (GAS_CONSTANT * config['reference_temperature_K'])
    return dict(diameter_m=diameter, depth_m=depth, pressure_difference_Pa=pressure_pa,
                plan_area_m2=area, roof_area_m2=roof, cap_volume_m3=cap, bowl_volume_m3=bowl,
                gross_gas_mass_kg=gas, uplift_N=pressure_pa * area,
                mean_vertical_rim_reaction_N_m=pressure_pa * diameter / 4,
                membrane_force_N_m=tension, equivalent_material_thickness_m=thickness,
                equivalent_pressure_material_kg=thickness * roof * config['material_density_kg_m3'],
                shielding=[dict(thickness_m=t, mass_kg=t * config['shield_density_kg_m3'] * roof,
                                weight_N=t * config['shield_density_kg_m3'] * roof * MOON_SURFACE_GRAVITY,
                                pressure_to_balance_weight_Pa=t * config['shield_density_kg_m3'] * roof * MOON_SURFACE_GRAVITY / area)
                           for t in config['shield_thickness_cases_m']])


def travel_time(distance_km, speed_km_h, acceleration):
    if distance_km < 0 or speed_km_h <= 0 or acceleration <= 0:
        raise ValueError('Nonnegative distance and positive speed and acceleration required')
    d, v = distance_km * 1000, speed_km_h / 3.6
    if d <= v * v / acceleration:
        return 2 * math.sqrt(d / acceleration)
    return d / v + v / acceleration


def engineering(rows, selected, population, scenarios):
    c, op = scenarios['enclosure'], scenarios['operation']
    pressure = c['reference_pressure_atm'] * STANDARD_ATMOSPHERE
    # Geometrically similar roofs have A,S proportional to D², Vcap and material to D³.
    unit = enclosure(1., 0., pressure, c)
    d2, d3 = float(np.sum(rows['diameter'][selected] ** 2)), float(np.sum(rows['diameter'][selected] ** 3))
    bowl = math.pi / 8 * float(np.sum(rows['diameter'][selected] ** 2 * rows['depth'][selected]))
    volume = unit['cap_volume_m3'] * d3 + bowl
    gas = volume * pressure * c['gas_molar_mass_kg_mol'] / (GAS_CONSTANT * c['reference_temperature_K'])
    roof = unit['roof_area_m2'] * d2
    hours = SYNODIC_MONTH_DAYS * 24 / 2
    operations = []
    for power in op['electric_load_kW_per_resident_cases']:
        mw = population * power / 1000
        operations.append(dict(electric_kW_per_resident=power, electric_MW=mw,
                               half_cycle_delivered_MWh=mw * hours,
                               installed_storage_MWh=mw * hours * (1 + op['storage_reserve_fraction']) / (op['storage_discharge_efficiency'] * op['storage_usable_fraction']),
                               ideal_radiator_area_km2={str(t): mw * 1e6 / (op['radiator_emissivity'] * STEFAN_BOLTZMANN * t ** 4) / 1e6 for t in op['radiator_temperature_K_cases']}))
    years = op['duration_cases_years']
    return dict(reference_cases=[enclosure(d, d * c['reference_depth_over_diameter'], p * STANDARD_ATMOSPHERE, c)
                                 for p in c['pressure_cases_atm'] for d in c['diameters_m']],
                whole_crater_initial_scenario=dict(gross_volume_m3=volume, gross_air_mass_tonnes=gas / 1000,
                    roof_area_km2=roof / 1e6, pressure_uplift_N=pressure * math.pi / 4 * d2,
                    equivalent_pressure_material_tonnes=unit['equivalent_pressure_material_kg'] * d3 / 1000,
                    shielding_tonnes={str(t): t * c['shield_density_kg_m3'] * roof / 1000 for t in c['shield_thickness_cases_m']},
                    annual_air_replacement_tonnes={str(f): gas * f / 1000 for f in op['annual_air_loss_fractions']},
                    cumulative_air_replacement_tonnes={str(y): {str(f): gas * f * y / 1000 for f in op['annual_air_loss_fractions']} for y in years}),
                material_improvement_cases=[dict(effective_allowable_Pa=stress,
                    thickness_and_pressure_material_factor=c['material_effective_allowable_Pa'] / stress,
                    initial_equivalent_pressure_material_tonnes=unit['equivalent_pressure_material_kg'] * d3 / 1000 * c['material_effective_allowable_Pa'] / stress) for stress in c['effective_allowable_cases_Pa']],
                larger_community_module_comparison_same_total_plan_area=[dict(span_m=d, cap_volume_relative_to_300m=d / 300,
                    ideal_pressure_material_relative_to_300m=d / 300, shielding_mass_relative_to_300m=1.) for d in c['module_span_m']],
                replacement_counts={str(y): {str(life): max(0, math.ceil(y / life) - 1) for life in op['replaceable_component_lives_years']} for y in years},
                lunar_cycles_by_duration={str(y): y * JULIAN_YEAR_DAYS / SYNODIC_MONTH_DAYS for y in years},
                full_architectural_occupancy_service_sensitivity=operations,
                service_sizing_residents=population,
                half_synodic_month_hours=hours,
                assumed_water_loop_makeup_tonnes_per_day=population * op['water_loop_throughput_kg_per_person_day'] * (1 - op['water_recovery_fraction']) / 1000,
                transition_pressure_atm=[dict(outside_atm=p, inside_atm=1., differential_atm=1 - p) for p in (0, .1, .5, 1., 1.2)],
                travel_cases=[dict(distance_km=d, speed_km_h=v, ideal_trip_hours=travel_time(d, v, scenarios['transport']['acceleration_m_s2']) / 3600,
                                   braking_distance_m=(v / 3.6) ** 2 / (2 * scenarios['transport']['acceleration_m_s2']))
                              for d in scenarios['transport']['trip_distances_km'] for v in scenarios['transport']['speed_cases_km_h']])


def food_comparison(population, food, operation):
    """Transfer production-energy accounts only; vacuum farm thermal loads remain open."""
    power_made = food['systems']['power_made_food']['total_kw_per_resident']
    lamp = food['mixes']['fully_lamp_lit']['cases']['4']['kept_50']
    lamp_per_person = np.array(lamp['power_tw']['total']) / np.array(lamp['people_billion']['lunar'])
    result = {}
    for name, per_person in [('power_made_development_case', power_made), ('fully_lamp_lit_case', lamp_per_person)]:
        added = np.asarray(per_person) * population / 1e6
        total = added + population * operation['reference_electric_load_kW_per_resident'] / 1e6
        result[name] = dict(additional_food_kW_per_resident=list(map(float, per_person)),
                           additional_food_GW=added.tolist(),
                           services_plus_food_GW=total.tolist(),
                           half_synodic_month_delivered_GWh=(total * SYNODIC_MONTH_DAYS * 12).tolist())
    result['reading_rule'] = ('Inherited process-energy scenarios, added to the 5 kW/person community-service comparison. '
                             'Production stays on the Moon. Vacuum farm structures, source/crop heat rejection and reserves are additional. '
                             'Power-made food remains a development scenario: fats and plant compounds lack measured complete production routes. '
                             'These comparisons do not select a diet or establish a complete food system.')
    return result


def export_sites(path, rows, final, initial, groups, bins):
    initial_set = set(map(int, initial))
    group_index = {}
    offset = 0
    for j, group in enumerate(groups):
        for i in initial[offset:offset + group['count']]:
            group_index[int(i)] = j
        offset += group['count']
    partial = path.with_suffix(path.suffix + '.partial')
    with partial.open('wb') as raw, gzip.GzipFile(filename='', fileobj=raw, mode='wb', mtime=0) as z, io.TextIOWrapper(z, encoding='utf-8', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['crater_id', 'latitude_today_deg', 'longitude_today_deg', 'diameter_m', 'depth_m', 'geometry_error_bin', 'initial_vacuum', 'initial_group'])
        for i in final:
            r = rows[i]
            writer.writerow([r['id'].decode(), r['lat'], r['lon'], r['diameter'], r['depth'], bins[i], int(i in initial_set), group_index.get(int(i), '')])
    partial.replace(path)


def run(work_dir=DEFAULT_WORK, plan_file=None):
    start = time.monotonic()
    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((HERE / 'inputs.json').read_text())
    paths = {name: Path(record['path']) if Path(record['path']).is_absolute() else ROOT / record['path'] for name, record in manifest['inputs'].items()}
    if plan_file is not None:
        paths['plan'] = Path(plan_file)
    for name, path in paths.items():
        if sha256(path) != manifest['inputs'][name]['sha256']:
            raise ValueError(f'Input hash mismatch: {name}')
    census = read_product(paths['census'], 'terluna.research.crater-gardens/1')
    plan = read_product(paths['plan'], 'terluna.habitation.crater-buildout-plan/1')
    design = read_product(paths['design'], 'terluna.habitation.crater-garden-design/1')
    food = read_product(paths['food'], 'terluna.research.food-systems/1')
    scenarios = read_product(HERE / 'scenarios.json', 'terluna.research.vacuum-crater-scenarios/1')
    selection = scenarios['selection']
    denominator = census['catalogue']['diameter_band_rows']
    wanted = quota(denominator, plan['eventual_fraction'])
    first = quota(denominator, plan['initial_vacuum_fraction'])
    if first > wanted:
        raise ValueError('The initial quota exceeds the eventual quota')
    rows, stages = read_candidates(paths['candidates'], census, selection)
    xyz = positions(rows)
    bins, tie, order = ranks(rows, selection)
    comparisons, final = {}, None
    for factor in selection['estate_diameter_factors']:
        selected = separated_sites(xyz, rows['diameter'], order, wanted, factor)
        comparisons[str(factor)] = sample_summary(rows, xyz, selected, design, wanted, factor)
        if factor == selection['reference_estate_diameter_factor']:
            final = selected
    if final is None:
        raise ValueError('Reference estate case missing')
    coefficient = populations(np.array([1.]), design)[plan['initial_population_calculation_case']]
    diameter_bounds = np.sqrt(np.array(plan['initial_community_population']) / coefficient)
    initial, groups = clustered_initial(rows, xyz, final, first, bins, tie, selection, diameter_bounds)
    smaller = final[(rows['diameter'][final] >= diameter_bounds[0]) & (rows['diameter'][final] <= diameter_bounds[1])]
    independent = smaller[np.lexsort((tie[smaller], bins[smaller], rows['diameter'][smaller]))][:first]
    with np.load(paths['context'], allow_pickle=False) as z:
        metadata = json.loads(str(z['metadata']))
        if metadata['schema'] != 'terluna.research.crater-gardens-context/1' or metadata['frame'] != 'today':
            raise ValueError('Regional context must use the catalogue frame')
        water = z['water']
    factor = selection['reference_estate_diameter_factor']
    initial_summary = sample_summary(rows, xyz, initial, design, first, factor)
    alternative_summary = sample_summary(rows, xyz, independent, design, first, factor)
    initial_summary['regional_groups'] = len(groups)
    initial_summary['group_size'] = distribution([g['count'] for g in groups])
    initial_summary['network'] = network(xyz[initial], water, scenarios['transport'])
    alternative_summary['network'] = network(xyz[independent], water, scenarios['transport'])
    output = work_dir / 'planned_sites.csv.gz'
    export_sites(output, rows, final, initial, groups, bins)
    (work_dir / 'initial_groups.json').write_text(json.dumps(groups, indent=2) + '\n')
    pinned = [HERE / 'run.py', HERE / 'scenarios.json', HERE / 'inputs.json']
    result = dict(schema=SCHEMA,
        evidence='Reproducible proxy-ranked and estate-separated site scenarios from a published crater census; ideal enclosure and service calculations. Site feasibility, route engineering, lifecycle reliability and supported vacuum occupancy remain unvalidated.',
        reading_rule='Both percentages use all 200–1,000 m catalogue records before exclusions, rounded down. Initial sites are nested within the eventual selection. Architectural populations are summed from actual diameters; service loads use full representative floor-space occupancy as a sizing sensitivity. All enclosure material, leakage, power and transport choices are scenarios. Domes are removed after a sustainable breathable environment is established.',
        programme=dict(denominator_craters=denominator, eventual_sites=wanted, initial_sites=first, initial_fraction_of_eventual=first / wanted,
                       initial_population_per_site=plan['initial_community_population'],
                       initial_diameter_bounds_m=diameter_bounds.tolist(), initial_enclosure=plan['initial_enclosure']),
        screening=stages, eventual_by_estate_factor=comparisons, initial_clustered=initial_summary,
        initial_smallest_span_comparison=alternative_summary,
        engineering=engineering(rows, initial, initial_summary['architectural_population']['representative'], scenarios),
        food_production_comparison=food_comparison(initial_summary['architectural_population']['representative'], food, scenarios['operation']),
        selected_sites=dict(path=str(output), rows=len(final), initial_rows=len(initial), sha256=sha256(output)),
        initial_groups=dict(path=str(work_dir / 'initial_groups.json'), sha256=sha256(work_dir / 'initial_groups.json')),
        unresolved=['Local rim/floor topography, foundations, gas barriers and radiation transport',
                    'Future shoreline elevations, groundwater and transient water trajectories',
                    'Conservation polygons, community review and wild/dark-land allocation',
                    'Terrain routes, future-water detours, redundancy, service demand and resource locations',
                    'Vacuum and transitional illumination, operational pressure/composition and refuge capacity',
                    'Heat balance, food/agriculture, imports, repair industry and multi-generation health',
                    'Enclosure dismantling engineering and material recovery'],
        producer=dict(files={str(p.relative_to(ROOT)): sha256(p) for p in pinned},
                      inputs={name: dict(path=str(path), sha256=sha256(path)) for name, path in paths.items()},
                      constants=constants_used(pinned)), seconds=time.monotonic() - start)
    (HERE / 'results/buildout.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(programme=result['programme'], screening=stages, final_selected=len(final),
                          initial_selected=len(initial), groups=len(groups), seconds=result['seconds']), indent=2))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, default=DEFAULT_WORK)
    parser.add_argument('--plan-file', type=Path, help='Exact-byte adopted plan in another checkout')
    args = parser.parse_args()
    run(args.work_dir, args.plan_file)
