"""Lunar Prospector composition and formal-error screening; no reserve inference.

The original PDS table is restored by exact hash into ignored research/runs.
Oxide nonclosure, covariance, shallow sensing and systematic-error limits remain
visible. Neither recovery nor a composition/depth relationship is assumed.
"""
from __future__ import annotations

import argparse
import bisect
import gzip
import hashlib
import json
import math
from pathlib import Path
import urllib.request
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
ORDER = ('MgO', 'Al2O3', 'SiO2', 'CaO', 'TiO2', 'FeO', 'K', 'Th', 'U')
OXIDES = ORDER[:6]
NAMESPACE = {'p': 'http://pds.nasa.gov/pds4/pds/v1'}
FIELD_NAMES = ('Pixel_Index', 'Southernmost_Latitude', 'Northernmost_Latitude',
               'Westernmost_Longitude', 'Easternmost_Longitude', 'Atomic_Mass',
               'Neutron_Density', 'W_MgO', 'W_Al2O3', 'W_SiO2', 'W_CaO',
               'TiO2', 'W_FeO', 'W_K', 'W_Th', 'W_U')


def digest(raw: bytes) -> str:
    return hashlib.sha256(raw).hexdigest()


def solid_angle(bounds) -> float:
    """Spherical rectangle area in steradians, with east-positive longitude."""
    south, north, west, east = map(float, bounds)
    if not (all(math.isfinite(v) for v in (south, north, west, east))
            and -90 <= south < north <= 90 and -180 <= west < east <= 180):
        raise ValueError('Invalid pixel bounds')
    return math.radians(east-west)*(math.sin(math.radians(north))-math.sin(math.radians(south)))


def parse_pds_table(table_bytes: bytes, label_bytes: bytes) -> list[dict]:
    """Decode fixed-width PDS fields without changing source abundances.

    Covariance order is ORDER. First six coordinates use mass fractions; K,
    Th and U use ppm. Cross-covariance units are products of these units.
    """
    xml = ET.fromstring(label_bytes)
    table = xml.find('.//p:Table_Character', NAMESPACE)
    if table is None:
        raise ValueError('Missing PDS character table')
    fields = table.findall('.//p:Field_Character', NAMESPACE)
    names = [f.find('p:name', NAMESPACE).text for f in fields]
    covariance_names = [f'E[{i},{j}]' for i in range(9) for j in range(i, 9)]
    if names != list(FIELD_NAMES)+covariance_names:
        raise ValueError('Unexpected PDS abundance/covariance field order')
    layout = [(int(f.find('p:field_location', NAMESPACE).text)-1,
               int(f.find('p:field_length', NAMESPACE).text)) for f in fields]
    record_size = int(table.find('.//p:record_length', NAMESPACE).text)
    expected_rows = int(table.find('p:records', NAMESPACE).text)
    raw_lines = table_bytes.splitlines(keepends=True)
    if len(raw_lines) != expected_rows or any(len(line) != record_size for line in raw_lines):
        raise ValueError('PDS row count or fixed-width record size mismatch')
    pixels = []
    for row_number, line in enumerate(raw_lines):
        values = [float(line[start:start+size]) for start, size in layout]
        if not all(math.isfinite(v) for v in values):
            raise ValueError('Nonfinite PDS number')
        if values[0] != row_number:
            raise ValueError('Unexpected or duplicate PDS pixel index')
        abundance = dict(zip(ORDER, values[7:16]))
        if any(v < 0 for v in abundance.values()):
            raise ValueError('Negative PDS abundance')
        covariance = [[0.0]*9 for _ in ORDER]
        cursor = 16
        for i in range(9):
            for j in range(i, 9):
                covariance[i][j] = covariance[j][i] = values[cursor]
                cursor += 1
            if covariance[i][i] < 0:
                raise ValueError('Negative formal variance')
        pixels.append({'pixel_index': row_number, 'bounds_degrees': values[1:5],
                       'solid_angle_sr': solid_angle(values[1:5]),
                       'abundances': abundance, 'covariance': covariance,
                       'six_oxide_sum': math.fsum(abundance[k] for k in OXIDES)})
    return pixels


def linear_estimate(pixel: dict, coefficients: dict[str, float], sigma_multiple: float = 1.0) -> dict:
    """Formal propagation for a caller-defined linear composition expression.

    The coefficients must carry any desired unit/stoichiometric conversion.
    Intervals are untruncated formal-error screens, not reserve confidence
    bounds. No inter-pixel independence or recovery rate is assumed.
    """
    if set(coefficients)-set(ORDER):
        raise ValueError('Unknown abundance coefficient')
    a = [float(coefficients.get(k, 0.0)) for k in ORDER]
    if not all(math.isfinite(x) for x in a) or not math.isfinite(sigma_multiple) or sigma_multiple < 0:
        raise ValueError('Invalid coefficient or sigma multiple')
    mean = math.fsum(a[i]*pixel['abundances'][k] for i, k in enumerate(ORDER))
    terms = [a[i]*pixel['covariance'][i][j]*a[j] for i in range(9) for j in range(9)]
    variance = math.fsum(terms)
    tolerance = 1e-12*max(math.fsum(abs(x) for x in terms), 1e-300)
    if variance < -tolerance:
        raise ValueError('Source covariance gives a negative propagated variance')
    sigma = math.sqrt(max(variance, 0.0))
    return {'mean': mean, 'formal_sigma': sigma, 'sigma_multiple': sigma_multiple,
            'formal_interval': [mean-sigma_multiple*sigma, mean+sigma_multiple*sigma],
            'interpretation': 'Formal propagation only; systematic, depth, mineralogy, recovery and access uncertainties excluded.'}


def weighted_quantiles(values, weights, probabilities=(0.05, 0.25, 0.5, 0.75, 0.95)) -> list[float]:
    """Interpolate against sorted cumulative area; not statistical confidence."""
    values, weights, probabilities = list(values), list(weights), list(probabilities)
    if not values or len(values) != len(weights):
        raise ValueError('Values and weights must be nonempty and paired')
    if (not all(math.isfinite(v) for v in values+weights+probabilities)
            or any(w < 0 for w in weights) or any(not 0 <= p <= 1 for p in probabilities)):
        raise ValueError('Invalid quantile arguments')
    pairs = sorted((v, w) for v, w in zip(values, weights) if w > 0)
    if not pairs:
        raise ValueError('Positive total weight required')
    total = math.fsum(w for _, w in pairs)
    cumulative, running = [], 0.0
    for _, w in pairs:
        running += w
        cumulative.append(running/total)
    result = []
    for q in probabilities:
        j = bisect.bisect_left(cumulative, q)
        if j == 0:
            result.append(pairs[0][0])
        elif j >= len(pairs):
            result.append(pairs[-1][0])
        else:
            fraction = (q-cumulative[j-1])/(cumulative[j]-cumulative[j-1])
            result.append(pairs[j-1][0]+fraction*(pairs[j][0]-pairs[j-1][0]))
    return result


def summarize(pixels: list[dict]) -> dict:
    """Summarize any explicitly selected cells; do not infer site access."""
    if not pixels:
        raise ValueError('No selected pixels')
    weights = [p['solid_angle_sr'] for p in pixels]
    total = math.fsum(weights)
    summary = {}
    for k in ORDER:
        values = [p['abundances'][k] for p in pixels]
        sigmas = [linear_estimate(p, {k: 1.0})['formal_sigma'] for p in pixels]
        summary[k] = {'unit': 'mass_fraction' if k in OXIDES else 'ppm',
                      'area_weighted_mean': math.fsum(v*w for v, w in zip(values, weights))/total,
                      'minimum': min(values), 'maximum': max(values),
                      'spatial_quantiles': weighted_quantiles(values, weights),
                      'formal_sigma_spatial_quantiles': weighted_quantiles(sigmas, weights)}
    oxide_sums = [p['six_oxide_sum'] for p in pixels]
    return {'pixel_count': len(pixels), 'solid_angle_sr': total,
            'fraction_of_lunar_surface': total/(4*math.pi),
            'quantile_probabilities': [0.05, 0.25, 0.5, 0.75, 0.95],
            'quantile_rule': 'Linear interpolation against sorted cumulative spherical area; spatial distributions, not confidence intervals.',
            'composition': summary,
            'six_oxide_sum': {'minimum': min(oxide_sums), 'maximum': max(oxide_sums),
                              'spatial_quantiles': weighted_quantiles(oxide_sums, weights)},
            'aggregate_assay_uncertainty': None,
            'aggregate_assay_uncertainty_reason': 'Cross-pixel statistical covariance and common systematic errors are not supplied.'}



def latitude_band(pixels: list[dict], south: float, north: float) -> list[dict]:
    """Clip cells to a latitude band, retaining their observed composition.

    This is a geometric weight operation, not finer spatial information.
    Overlapping selections are correlated and must not be added as new data.
    """
    if not (-90 <= south < north <= 90):
        raise ValueError('Invalid latitude band')
    selected = []
    for pixel in pixels:
        lo, hi, west, east = pixel['bounds_degrees']
        clipped = [max(lo, south), min(hi, north), west, east]
        if clipped[0] < clipped[1]:
            selected.append({**pixel, 'bounds_degrees': clipped,
                             'solid_angle_sr': solid_angle(clipped)})
    return selected


def load_inputs(config=None, archive=None, download=False):
    """Restore only the expected bytes; reject mismatches, never substitute."""
    config = config or json.loads((HERE/'lunar_inputs.json').read_text())
    folder = Path(archive) if archive is not None else ROOT/config['archive_directory']
    blobs = {}
    for key, source in config['files'].items():
        path = folder/source['filename']
        if not path.exists() and download:
            raw = urllib.request.urlopen(source['url'], timeout=60).read()
            if len(raw) != source['bytes'] or digest(raw) != source['sha256']:
                raise ValueError(f'Downloaded {key} source hash or length mismatch')
            folder.mkdir(parents=True, exist_ok=True)
            path.write_bytes(raw)
        raw = path.read_bytes()
        if len(raw) != source['bytes'] or digest(raw) != source['sha256']:
            raise ValueError(f'Archived {key} source hash or length mismatch')
        blobs[key] = raw
    return config, blobs


def product(config=None, archive=None, download=False):
    config, blobs = load_inputs(config, archive, download)
    pixels = parse_pds_table(blobs['table'], blobs['label'])
    if len(pixels) != config['expected_pixel_count']:
        raise ValueError('Unexpected global pixel count')
    total = math.fsum(p['solid_angle_sr'] for p in pixels)
    if not math.isclose(total, 4*math.pi, rel_tol=0, abs_tol=1e-10):
        raise ValueError('Global pixel areas do not cover 4 pi steradians')
    return {'schema': 'terluna.resources.lunar-composition/1',
            'producer': {'path': 'resources/source_access/lunar.py', 'sha256': digest(Path(__file__).read_bytes())},
            'inputs_canonical_sha256': digest((json.dumps(config, sort_keys=True, separators=(',', ':'))+'\n').encode()),
            'inputs_hash_rule': 'SHA-256 of sorted compact JSON plus newline; hashes input values rather than formatting.',
            'source_register': {'path': 'resources/source_access/lunar_sources.json', 'sha256': digest((HERE/'lunar_sources.json').read_bytes())},
            'source': config,
            'evidence': 'Derived spatial composition and formal assay propagation from archived mission data; no mine, depth extrapolation, extraction yield or economic reserve qualification.',
            'covariance_order': list(ORDER),
            'covariance_units': 'Products of coordinate units: mass fraction for six oxides, ppm for K/Th/U.',
            'formal_interval_rule': 'mean plus/minus explicitly selected sigma multiple; no clipping, coverage probability or systematic-error inclusion.',
            'systematic_error': 'Unquantified here; original nonclosure and calibration limits retained.',
            'depth_extrapolation': None, 'mineralogy_conversion': None,
            'process_recovery': None, 'qualified_reserve': None,
            'validity': {'invalid_or_negative_pixels': 0, 'invalid_cell_rule': 'Reject invalid/negative source records; exact admitted table contains none. Preserve valid zero grades.', 'zero_grade_pixel_counts': {k: sum(p['abundances'][k] == 0 for p in pixels) for k in ORDER}},
            'global_summary': summarize(pixels),
            'latitude_bands': [{'south_degrees': south, 'north_degrees': south+30, 'summary': summarize(latitude_band(pixels, south, south+30))} for south in range(-90, 90, 30)],
            'pixels': pixels}


def write_product(result, path):
    path = Path(path)
    raw = (json.dumps(result, sort_keys=True, separators=(',', ':'), allow_nan=False)+'\n').encode()
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(gzip.compress(raw, mtime=0) if path.suffix == '.gz' else raw)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path)
    parser.add_argument('--download', action='store_true')
    parser.add_argument('--output', type=Path, default=HERE/'results/lunar.json.gz')
    args = parser.parse_args()
    result = product(archive=args.archive, download=args.download)
    write_product(result, args.output)
    print(f"{args.output}: {len(result['pixels'])} pixels; no reserve qualification")


if __name__ == '__main__':
    main()
