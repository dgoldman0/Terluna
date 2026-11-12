"""Geometric controls, source parsing and exclusion boundaries for the crater census."""
from collections import Counter
import csv
import io
import json
import math
import os
from pathlib import Path
import zipfile

import numpy as np
import pytest

from research.studies.crater_gardens import count


@pytest.fixture
def design():
    path = Path(os.environ.get('TERLUNA_CRATER_DESIGN', count.ROOT / 'habitation/crater_gardens/design.json'))
    return json.loads(path.read_text())


def test_adopted_floor_space_and_area_scaling(design):
    assert count.population(750, design) == pytest.approx(4970.09775, rel=1e-7)
    assert count.population(1000, design) == pytest.approx(4 * count.population(500, design))
    # Independent sum of five annular building strips; terraces contribute side by side.
    inner, width = 187.5, 37.5
    floor = sum(2 * math.pi * (inner + (i + .5) * width) * width * .35 * 3 for i in range(5))
    assert floor / 70 == pytest.approx(count.population(750, design))
    for diameter, table in design['rounded_population_by_diameter_m'].items():
        for case, people in table.items():
            assert count.population(float(diameter), design, case) == pytest.approx(people, rel=.025)


def test_water_cell_edges_poles_and_antimeridian():
    r, c = count.grid_cells([90, -90, 0, 0, 0, 0], [0, 0, -180, 180, -.001, 360], (180, 360))
    assert r.tolist() == [0, 179, 90, 90, 90, 90]
    assert c.tolist() == [0, 0, 180, 180, 359, 0]
    with pytest.raises(ValueError):
        count.grid_cells([90.01], [0], (180, 360))


def test_climate_interpolation_linear_field_and_seam():
    lat, lon = np.array([60., 0., -60.]), np.array([0., 90., 180., 270.])
    field = lat[:, None] * 2 + np.array([0., 90., 180., 90.])[None, :]
    values = count.sample(field, lat, lon, [30, -30, 90, 0, 0], [45, 135, 0, -45, 315])
    assert values.tolist() == pytest.approx([105, 75, 120, 45, 45])


def catalogue(tmp_path, rows, header=None):
    path = tmp_path / 'craters.zip'
    buf = io.StringIO(newline='')
    writer = csv.writer(buf)
    writer.writerow(header or count.REQUIRED)
    writer.writerows(rows)
    with zipfile.ZipFile(path, 'w', compression=zipfile.ZIP_DEFLATED) as z:
        z.writestr('LUC-GRAS200_full.csv', buf.getvalue())
    return path


def test_stream_unsorted_inclusive_band_and_missing_depth(tmp_path):
    rows = [('A', 0, 0, 1000, 70, 'high', 0, -1), ('B', 0, 0, 199.9, 20, 'high', 0, -1),
            ('C', 0, 0, 200, '', 'medium', 0, -1), ('D', 0, 0, 500, 40, 'high', 0, -1),
            ('E', 0, 0, 1000.1, 40, 'high', 0, -1)]
    stats = Counter()
    batches = list(count.batches(catalogue(tmp_path, rows), (200, 1000), stats, batch_size=2))
    assert stats['catalogue_rows'] == 5
    assert stats['diameter_band_rows'] == 3
    assert [len(b) for b in batches] == [2, 1]
    assert np.isnan(batches[0]['depth_m'][1])
    with pytest.raises(ValueError, match='fields missing'):
        list(count.batches(catalogue(tmp_path, [], header=['diameter_m']), (200, 1000), Counter()))


def test_reference_and_small_scaled_relief_quality_and_water(design):
    rows = np.array([('A', 0, 0, 750, 40, 'high', 0, -1),
                     ('B', 0, 0, 1000, 70, 'medium', 0, 1),
                     ('C', 0, 0, 750, 39.99, 'high', 0, -1),
                     ('D', 0, 0, 750, 70.01, 'high', 0, -1),
                     ('E', 0, 0, 750, 50, 'low', 0, -1),
                     ('F', 0, 0, 750, 50, 'high', 1, -1),
                     ('G', 0, 90, 750, 50, 'high', 0, -1)], dtype=count.DTYPE)
    context = dict(water=np.array([[0, 1, 2, 0], [0, 1, 2, 0]], dtype=np.uint8),
                   lat=np.array([45., -45.]), lon=np.arange(4) * 90.,
                   pr_mm_day=np.ones((2, 4)), ts_k=np.full((2, 4), 300.))
    values = count.assess(rows, context)
    counts, selected, pop = count.tally(rows, values, design)
    assert selected.tolist() == [True, True, False, False, False, False, False]
    assert counts['reference_depth'] == 5
    assert counts['reference_quality'] == 3
    assert counts['reference_rain_1'] == counts['selected'] == 2
    assert pop['representative'] == pytest.approx(count.population(750, design) + count.population(1000, design))
    # Origin as a secondary crater is retained: physical eligibility does not depend on origin alone.
    assert selected[1]
    small = np.array([('S', 0, 0, 300, 20, 'high', 0, -1)], dtype=count.DTYPE)
    sc, mask, _ = count.tally(small, count.assess(small, context), design, scaled_selection=True)
    assert sc['reference_depth'] == 0
    assert sc['scaled_rain_1'] == sc['selected'] == 1
    assert mask[0]
    context['pr_mm_day'][:] = .99999
    assert count.tally(rows, count.assess(rows, context), design)[0]['selected'] == 0
    context['pr_mm_day'][:] = 1.
    context['ts_k'][:] = 308.150001
    assert count.tally(rows, count.assess(rows, context), design)[0]['selected'] == 0


def test_changed_input_is_an_evidence_error(tmp_path):
    path = tmp_path / 'input'
    path.write_bytes(b'original')
    record = dict(bytes=path.stat().st_size, sha256=count.sha256(path))
    assert count.verified(path, record) == path
    path.write_bytes(b'differnt')
    with pytest.raises(ValueError, match='recorded bytes'):
        count.verified(path, record)


def test_saved_count_provenance_and_candidate_floor_sums(design):
    product = count.HERE / 'results/first_count.json'
    if not product.exists():
        pytest.skip('The first catalogue count has not yet been generated')
    result = json.loads(product.read_text())
    for relative, expected in result['producer']['files'].items():
        assert count.sha256(count.ROOT / relative) == expected
    assert result['catalogue']['catalogue_rows'] == 10699651
    assert sum(g['diameter_band_rows'] for g in result['counts'].values()) == result['catalogue']['diameter_band_rows'] - result['catalogue'].get('invalid_coordinates_in_band', 0)
    for name, group in result['counts'].items():
        assert group['reference_depth'] >= group['reference_quality'] >= group['reference_dry']
        assert group['reference_dry'] >= group['reference_ground_temperature'] >= group['reference_rain_0.5']
        assert group['reference_rain_0.5'] >= group['reference_rain_1'] >= group['reference_rain_2']
        window = 'reference' if name == '500_1000' else 'scaled'
        assert group['selected'] == group[window + '_rain_1']
    candidate_file = Path(result['candidates']['path'])
    if not candidate_file.exists():
        return  # The source-derived per-crater file is intentionally kept out of Git.
    assert count.sha256(candidate_file) == result['candidates']['sha256']
    import gzip
    counts, diameter_squared = Counter(), Counter()
    with gzip.open(candidate_file, 'rt', encoding='utf-8', newline='') as stream:
        for row in csv.DictReader(stream):
            band = row['diameter_band_m']
            d, depth = float(row['diameter_m']), float(row['depth_m'])
            if band == '500_1000':
                assert 500 <= d <= 1000 and 40 <= depth <= 70
            else:
                lo, hi = (200, 300) if band == '200_300' else (300, 500)
                assert lo <= d < hi
                assert 40 * d / 750 <= depth <= 70 * d / 750
            assert row['confidence_level'] in ('high', 'medium') and row['resolution_flag'] == '0'
            assert float(row['rain_mm_day']) >= 1. - 1e-6
            counts[band] += 1
            diameter_squared[band] += d * d
    assert sum(counts.values()) == result['candidates']['rows']
    for band, number in counts.items():
        assert number == result['counts'][band]['selected']
        for name, p in design['population_cases'].items():
            # Independent floor sum from the exported list, using the adopted annular area.
            expected = math.pi / 4 * diameter_squared[band] * .75 * p['non_lake_ground_building_share'] * p['mean_storeys'] / p['floor_m2_per_resident']
            assert abs(expected - result['summed_floor_space_population'][band][name]) < .51
