"""Count published crater bowls against the adopted garden-town dimensions and regional water and climate."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import gzip
import hashlib
import io
import json
import math
from pathlib import Path
import time
import zipfile

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DEFAULT_WORK = Path('/media/projectspace/terluna-research/crater-settlements')
SCHEMA = 'terluna.research.crater-gardens/1'
DTYPE = np.dtype([('id', 'U32'), ('latitude', 'f8'), ('longitude', 'f8'),
                  ('diameter_m', 'f8'), ('depth_m', 'f8'), ('confidence', 'U12'),
                  ('resolution_flag', 'i2'), ('secondary_flag', 'i2')])
REQUIRED = ('crater_id', 'latitude', 'longitude', 'diameter_m', 'depth_m',
            'confidence_level', 'flag_resolution_boundary', 'secondary_flag')
DEPTH_WINDOWS = {'reference': (40., 70.), 'broader': (30., 100.)}
RAIN_THRESHOLDS = (0.5, 1., 2.)
GROUND_K = (283.15, 308.15)


def sha256(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def verified(path, record):
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f'Restore the recorded input: {path}')
    if path.stat().st_size != record['bytes'] or sha256(path) != record['sha256']:
        raise ValueError(f'Input differs from its recorded bytes: {path}')
    return path


def population(diameter_m, design, case='representative'):
    """Built floor for homes, work and services divided by its allowance per permanent resident."""
    p = design['population_cases'][case]
    area = math.pi / 4 * np.asarray(diameter_m) ** 2
    floor = area * (1 - design['population_lake_to_settlement_diameter'] ** 2)
    return floor * p['non_lake_ground_building_share'] * p['mean_storeys'] / p['floor_m2_per_resident']


def grid_cells(latitude, longitude, shape):
    """Cell membership on a global cylindrical pixel-centre grid; the antimeridian is periodic."""
    latitude, longitude = np.asarray(latitude), np.asarray(longitude)
    if np.any(~np.isfinite(latitude)) or np.any(~np.isfinite(longitude)) or np.any(np.abs(latitude) > 90):
        raise ValueError('Finite planetocentric coordinates with latitude within [-90, 90] are required')
    rows, cols = shape
    r = np.minimum(((90 - latitude) / 180 * rows).astype(np.int64), rows - 1)
    c = ((longitude % 360) / 360 * cols).astype(np.int64) % cols
    return r, c


def sample(field, latitude, longitude, latitude_points, longitude_points):
    """Bilinear regional climate sampling with periodic longitude and the end latitude rows held."""
    latitude, longitude = np.asarray(latitude), np.asarray(longitude)
    order = np.argsort(latitude)
    f, lat = np.asarray(field)[order], latitude[order]
    lp = np.clip(np.asarray(latitude_points), lat[0], lat[-1])
    r = np.interp(lp, lat, np.arange(lat.size))
    r0 = np.minimum(np.floor(r).astype(int), lat.size - 2)
    fr = r - r0
    c = ((np.asarray(longitude_points) - longitude[0]) / (longitude[1] - longitude[0])) % longitude.size
    c0 = np.floor(c).astype(int)
    fc = c - c0
    c1 = (c0 + 1) % longitude.size
    top = f[r0, c0] * (1 - fc) + f[r0, c1] * fc
    bottom = f[r0 + 1, c0] * (1 - fc) + f[r0 + 1, c1] * fc
    return top * (1 - fr) + bottom * fr


def batches(archive, diameter_bounds, statistics, batch_size=50000):
    """Stream the complete CSV, retaining only the diameter band in bounded batches. ZIP CRC is checked at EOF."""
    buffer = []
    with zipfile.ZipFile(archive) as z:
        members = [m for m in z.infolist() if m.filename == 'LUC-GRAS200_full.csv']
        if len(members) != 1:
            raise ValueError('The archive must contain exactly one LUC-GRAS200_full.csv')
        statistics['csv_bytes'] = members[0].file_size
        statistics['csv_crc32'] = f'{members[0].CRC:08x}'
        with z.open(members[0]) as raw, io.TextIOWrapper(raw, encoding='utf-8-sig', newline='') as stream:
            reader = csv.reader(stream)
            header = next(reader)
            if any(k not in header for k in REQUIRED):
                raise ValueError(f'Required catalogue fields missing: {header}')
            statistics['header'] = header
            ix = {name: header.index(name) for name in REQUIRED}
            for row in reader:
                statistics['catalogue_rows'] += 1
                if len(row) != len(header):
                    raise ValueError(f'Malformed catalogue row {statistics["catalogue_rows"]}')
                d = float(row[ix['diameter_m']])
                if not math.isfinite(d) or d <= 0:
                    statistics['invalid_diameter'] += 1
                    continue
                if not diameter_bounds[0] <= d <= diameter_bounds[1]:
                    continue
                statistics['diameter_band_rows'] += 1
                la, lo = float(row[ix['latitude']]), float(row[ix['longitude']])
                if not math.isfinite(la) or not math.isfinite(lo) or abs(la) > 90 or abs(lo) > 180:
                    statistics['invalid_coordinates_in_band'] += 1
                    continue
                depth_text = row[ix['depth_m']].strip()
                depth = float(depth_text) if depth_text else float('nan')
                buffer.append((row[ix['crater_id']], la, lo, d, depth,
                               row[ix['confidence_level']].lower(), int(row[ix['flag_resolution_boundary']]),
                               int(row[ix['secondary_flag']])))
                if len(buffer) == batch_size:
                    yield np.array(buffer, dtype=DTYPE)
                    buffer.clear()
            if buffer:
                yield np.array(buffer, dtype=DTYPE)


def assess(rows, context):
    """Attach broad water and climate classes. Depth is the published rim-to-floor difference."""
    r, c = grid_cells(rows['latitude'], rows['longitude'], context['water'].shape)
    water = context['water'][r, c]
    near_water = np.zeros(rows.size, bool)
    for dr in (-1, 0, 1):
        for dc in (-1, 0, 1):
            rr = np.clip(r + dr, 0, context['water'].shape[0] - 1)
            cc = (c + dc) % context['water'].shape[1]
            near_water |= context['water'][rr, cc] != 0
    args = (context['lat'], context['lon'], rows['latitude'], rows['longitude'])
    rain = sample(context['pr_mm_day'], *args)
    ground = sample(context['ts_k'], *args)
    measured = np.isfinite(rows['depth_m']) & (rows['depth_m'] > 0)
    confidence = np.isin(rows['confidence'], ('high', 'medium'))
    quality = measured & confidence & (rows['resolution_flag'] == 0)
    warm = np.isfinite(ground) & (ground >= GROUND_K[0]) & (ground <= GROUND_K[1])
    climate = warm & np.isfinite(rain) & (rain >= 1.)
    return dict(water=water, near_water=near_water, rain=rain, ground=ground, measured=measured,
                confidence=confidence, quality=quality, warm=warm, climate=climate)


def tally(rows, values, design, scaled_selection=False):
    """Return stage counts and alternative cuts for one batch; sensitivity counts share the same source rows."""
    counts = Counter()
    measured, quality, water = values['measured'], values['quality'], values['water']
    counts['positive_measured_depth'] = int(measured.sum())
    for name in ('high', 'medium', 'low'):
        counts['confidence_' + name] = int((rows['confidence'] == name).sum())
    counts['other_confidence'] = int((~np.isin(rows['confidence'], ('high', 'medium', 'low'))).sum())
    for flag in (0, 1, 2):
        counts[f'resolution_flag_{flag}'] = int((rows['resolution_flag'] == flag).sum())
    counts['dry_centres_all_depths'] = int((water == 0).sum())
    counts['diameter_band_rows'] = rows.size
    windows = dict(DEPTH_WINDOWS)
    windows['scaled'] = (40 * rows['diameter_m'] / design['representative_diameter_m'],
                         70 * rows['diameter_m'] / design['representative_diameter_m'])
    windows['all_measured'] = (0., np.inf)
    for window, (lo, hi) in windows.items():
        depth = measured & (rows['depth_m'] >= lo) & (rows['depth_m'] <= hi)
        masks = {
            'depth': depth,
            'confidence': depth & values['confidence'],
            'quality': depth & quality,
            'outside_common_water': depth & quality & ((water & 1) == 0),
            'dry': depth & quality & (water == 0),
            'ground_temperature': depth & quality & (water == 0) & values['warm'],
        }
        for stage, mask in masks.items():
            counts[f'{window}_{stage}'] = int(mask.sum())
        base = masks['ground_temperature']
        for minimum in RAIN_THRESHOLDS:
            eligible = base & np.isfinite(values['rain']) & (values['rain'] >= minimum)
            counts[f'{window}_rain_{minimum:g}'] = int(eligible.sum())
        # Compare sensitivity to source quality flags without altering the reference selection.
        counts[f'{window}_all_quality_rain_1'] = int((depth & (water == 0) & values['climate']).sum())
        counts[f'{window}_high_only_rain_1'] = int((base & values['climate'] & (rows['confidence'] == 'high')).sum())
    lo, hi = windows['scaled'] if scaled_selection else windows['reference']
    reference = quality & (rows['depth_m'] >= lo) & (rows['depth_m'] <= hi) & (water == 0) & values['climate']
    counts['selected'] = int(reference.sum())
    counts['selected_adjacent_water_cell'] = int((reference & values['near_water']).sum())
    for lo, hi in ((0, 15), (15, 30), (30, 45), (45, 60), (60, 90.01)):
        mask = reference & (np.abs(rows['latitude']) >= lo) & (np.abs(rows['latitude']) < hi)
        counts[f'selected_today_latitude_{lo:g}_{hi:g}'] = int(mask.sum())
    pop = {case: float(population(rows['diameter_m'][reference], design, case).sum())
           for case in design['population_cases']}
    return counts, reference, pop


def run(work_dir=DEFAULT_WORK, design_file=None):
    start = time.monotonic()
    manifest = json.loads((HERE / 'inputs.json').read_text())
    archive = verified(Path(work_dir) / 'inputs/LUC-GRAS200_full.zip', manifest['catalogue'])
    context_file = verified(Path(work_dir) / 'context.npz', manifest['context'])
    design_path = verified(design_file or ROOT / 'habitation/crater_gardens/design.json', manifest['design'])
    design = json.loads(design_path.read_text())
    with np.load(context_file, allow_pickle=False) as z:
        metadata = json.loads(str(z['metadata']))
        context = {k: z[k] for k in z.files if k != 'metadata'}
    if metadata['schema'] != 'terluna.research.crater-gardens-context/1' or metadata['frame'] != 'today':
        raise ValueError('The water and climate context must use the catalogue\'s present-day coordinates')
    groups = {'200_300': (design['comparison_diameter_m'][0], 300),
              '300_500': (300, design['diameter_m'][0]),
              '500_1000': tuple(design['diameter_m'])}
    statistics = Counter()
    counts = {group: Counter() for group in groups}
    populations = {group: {case: 0. for case in design['population_cases']} for group in groups}
    output = Path(work_dir) / 'candidates.csv.gz'
    partial = output.with_name(output.name + '.partial')
    with partial.open('wb') as binary, gzip.GzipFile(filename='', mode='wb', fileobj=binary, mtime=0) as zipped, io.TextIOWrapper(zipped, encoding='utf-8', newline='') as stream:
        writer = csv.writer(stream)
        writer.writerow(['diameter_band_m', 'crater_id', 'latitude_today_deg', 'longitude_today_deg', 'diameter_m', 'depth_m',
                         'confidence_level', 'resolution_flag', 'secondary_flag', 'rain_mm_day',
                         'mean_ground_temperature_K', 'adjacent_water_cell'])
        bounds = (design['comparison_diameter_m'][0], design['diameter_m'][1])
        for rows in batches(archive, bounds, statistics):
            values = assess(rows, context)
            for name, (lo, hi) in groups.items():
                inside = (rows['diameter_m'] >= lo) & ((rows['diameter_m'] <= hi) if name == '500_1000' else (rows['diameter_m'] < hi))
                sub = rows[inside]
                val = {k: v[inside] for k, v in values.items()}
                batch_counts, selected, pop = tally(sub, val, design, scaled_selection=name != '500_1000')
                counts[name].update(batch_counts)
                for key in populations[name]:
                    populations[name][key] += pop[key]
                for i in np.flatnonzero(selected):
                    row = sub[i]
                    writer.writerow([name, row['id'], row['latitude'], row['longitude'], row['diameter_m'], row['depth_m'],
                                     row['confidence'], row['resolution_flag'], row['secondary_flag'],
                                     round(float(val['rain'][i]), 6), round(float(val['ground'][i]), 6),
                                     int(val['near_water'][i])])
    if statistics['catalogue_rows'] != manifest['catalogue']['rows']:
        raise ValueError('The catalogue row count differs from its published record')
    partial.replace(output)
    population_table = {str(d): {k: round(float(population(d, design, k)), 3) for k in design['population_cases']}
                        for d in (200, 300, 500, 750, 1000)}
    result = dict(
        schema=SCHEMA,
        evidence='A streamed count of published crater detections under adopted size and floor-space assumptions, '
                 'published rim-to-floor depth estimates, source quality classes, a coarse standing-water mask and '
                 'regional climate means. Candidate counts precede terrain, water-level, foundation, access, '
                 'heritage, spacing and services assessment.',
        reading_rule='Counts retain source detections once each. The reference selection uses diameter 500–1000 m '
                     'inclusive, depth 40–70 m inclusive, high or medium detection confidence, resolution flag 0, '
                     'a dry centre cell, mean ground temperature 283.15–308.15 K and rain at least 1 mm/day. '
                     'Depth is a proxy for the representative bowl relief; lake level and bed shaping are open. '
                     'The 200–300 m and 300–500 m comparison bands use scaled relief (40/750 to 70/750 of '
                     'diameter) for their selected CSV entries; all bands also retain counts under the same '
                     'absolute depth cuts and with every positive depth. Bands are lower-inclusive and '
                     'upper-exclusive except for the 1000 m upper cap. '
                     'Candidate CSV coordinates remain in today\'s frame; the climate product was carried back '
                     'to that frame by its producer. Population sums are architectural capacity scenarios.',
        criteria=dict(diameter_groups_m=groups, depth_windows_m=DEPTH_WINDOWS,
                      scaled_depth_to_diameter=[40 / design['representative_diameter_m'], 70 / design['representative_diameter_m']],
                      rain_thresholds_mm_day=RAIN_THRESHOLDS, mean_ground_temperature_K=GROUND_K),
        catalogue=dict(statistics), counts={k: dict(v) for k, v in counts.items()}, context=metadata,
        population_by_diameter=population_table, design_source=str(design_path),
        summed_floor_space_population={g: {k: round(v) for k, v in p.items()} for g, p in populations.items()},
        candidates=dict(path=str(output), rows=sum(v['selected'] for v in counts.values()), sha256=sha256(output)),
        producer=dict(files={str(p.relative_to(ROOT)): sha256(p) for p in (Path(__file__), HERE / 'inputs.json')},
                      inputs={str(archive): sha256(archive), str(context_file): sha256(context_file),
                              'habitation/crater_gardens/design.json': sha256(design_path)}),
        seconds=round(time.monotonic() - start, 3))
    target = HERE / 'results/first_count.json'
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps(dict(catalogue_rows=statistics['catalogue_rows'], diameter_band=statistics['diameter_band_rows'],
                          selected={g: v['selected'] for g, v in counts.items()},
                          seconds=result['seconds']), indent=2))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--work-dir', type=Path, default=DEFAULT_WORK)
    parser.add_argument('--design', type=Path, help='the adopted design JSON at its recorded hash')
    args = parser.parse_args()
    run(args.work_dir, args.design)
