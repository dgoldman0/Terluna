"""Export qualified growth-repair records without changing frozen simulations.

    python -m atmosphere.snow.growth_products

The frozen assembler admits the simulations by their original source/input
identities. This second stage restores reference qualifications and adds only
explicit optical-use gates and a nominal-endpoint reporting diagnostic.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
from pathlib import Path
import tempfile

import numpy as np

from atmosphere.snow import growth_repair as frozen, repair

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'atmosphere/snow/results/growth_repair.json'
SCHEMA = 'terluna.atmosphere.snow-growth-final-export/1'
REFERENCE_IDS = ('earth_moderate_m144', 'moon_heavy_m96', 'moon_heavy_m144')
FILES = ('atmosphere/snow/growth_products.py', 'atmosphere/snow/growth_repair.py',
         'atmosphere/snow/repair.py')
ENDPOINT_RTOL = 64 * np.finfo(float).eps


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def cache_paths():
    return ([frozen.CACHE / (spec['id'] + '.json') for spec in frozen.cases()]
            + [repair.CACHE / (identity + '.json') for identity in REFERENCE_IDS])


def identity():
    """Include absent caches so a newly completed case cannot enter unnoticed."""
    return dict(files={p: sha(ROOT / p) for p in FILES},
                raw_caches={str(p.relative_to(ROOT)): sha(p) if p.exists() else None
                            for p in cache_paths()})


def _assemble_frozen():
    # The frozen assembler writes as well as returning its admitted records.
    # Redirect only that destination; source identities and calculations remain
    # untouched. The final destination is replaced only after all checks pass.
    original = frozen.OUT
    try:
        with tempfile.TemporaryDirectory(prefix='terluna-growth-assembly-') as directory:
            frozen.OUT = Path(directory) / 'intermediate.json'
            return frozen.assemble()
    finally:
        frozen.OUT = original


def endpoint_diagnostic(raw):
    """Interpret roundoff at the final pivot; never change stored distributions."""
    if 'masses_kg' in raw:
        mass_key, diameter_key = 'masses_kg', 'diameters_m'
    else:
        mass_key, diameter_key = 'arriving_mass_pivots_kg', 'arriving_diameter_pivots_m'
    masses = np.asarray(raw[mass_key], dtype=float)
    diameters = np.asarray(raw[diameter_key], dtype=float)
    numbers = np.asarray(raw['arriving_numbers_m2'], dtype=float)
    if (masses.ndim != 1 or masses.size < 2 or diameters.shape != masses.shape
            or numbers.shape != masses.shape or not np.isfinite(masses).all()
            or not np.isfinite(diameters).all() or not np.isfinite(numbers).all()
            or (masses <= 0).any() or (diameters <= 0).any() or (numbers < 0).any()
            or (np.diff(masses) <= 0).any() or (np.diff(diameters) <= 0).any()):
        raise ValueError('Invalid native arriving distribution: ' + raw['id'])
    weights = masses * numbers
    total = float(weights.sum())
    seconds = float(raw['gathering_seconds'])
    nominal = float(raw['specification']['maximum_diameter_m'])
    if not np.isfinite(weights).all() or not math.isfinite(total):
        raise ValueError('Nonfinite native arriving mass: ' + raw['id'])
    if not math.isfinite(seconds) or seconds <= 0 or not math.isfinite(nominal) or nominal <= 0:
        raise ValueError('Invalid gathering duration or nominal endpoint: ' + raw['id'])
    flux = total * 3600. / seconds
    if not math.isclose(flux, raw['summary']['snowfall_mm_h'], rel_tol=1e-12, abs_tol=0.):
        raise ValueError('Native arriving mass disagrees with saved flux: ' + raw['id'])
    matched = math.isclose(float(diameters[-1]), nominal, rel_tol=ENDPOINT_RTOL, abs_tol=0.)
    classified = diameters.copy()
    if matched:
        classified[-1] = nominal
    return dict(
        nominal_endpoint_m=nominal, stored_top_resolved_diameter_m=float(diameters[-1]),
        top_resolved_mass_kg=float(masses[-1]), endpoint_matches_within_roundoff=matched,
        classification_relative_tolerance=ENDPOINT_RTOL, classification_absolute_tolerance_m=0.,
        top_resolved_pivot_arriving_mass_share=float(weights[-1]/total) if total else None,
        arriving_mass_share_at_or_above_nominal_endpoint=
            float(weights[classified >= nominal].sum()/total) if total else None,
        arriving_mass_kg_m2=total,
        reading_rule='Only the final resolved diameter is matched to the configured endpoint, '
                     'and only within 64 binary64 machine epsilons in relative diameter. '
                     'All original threshold fields, pivots, populations and quantiles remain unchanged. '
                     'The added shares normalize surviving freezing-level arrivals over the gathering interval; '
                     'they exclude numerical-domain exits. Endpoint occupancy is a discretization diagnostic, '
                     'not evidence of physical particles at or above that size. Null shares mean no arrivals.',
    )


def _qualify(row, raw, *, reference):
    # Every field retained by the frozen assembler must match the cache whose
    # native arrays supply the new diagnostics and whose qualifications we read.
    if any(key not in raw or value != raw[key] for key, value in row.items()):
        raise RuntimeError('Assembled row differs from its raw cache: ' + row['id'])
    if reference:
        for key in ('reading_rule', 'evidence', 'numerical_status'):
            if not raw.get(key):
                raise ValueError('Missing reference qualification ' + key + ': ' + row['id'])
            row[key] = raw[key]
    row['optical_prediction_eligible'] = False
    row['optical_prediction_reason'] = (
        'This prescribed column has no calibrated storm-optics prediction. A passed finite-mass-domain '
        'check does not establish grid/time/height convergence or validated aggregate physics. '
        'Generic Earth controls are numerical transfer checks, not observed-storm or optical calibration. '
        'Any inherited extinction and 5% beam-distance fields are conditional Qext=2 homogeneous '
        'direct-beam diagnostics, not human visibility or a repaired lunar snowstorm forecast.'
    )
    if not row['domain_adequacy']['passed']:
        row['optical_prediction_reason'] += ' This case also fails its finite-mass-domain check.'
    row['nominal_endpoint_diagnostic'] = endpoint_diagnostic(raw)


def build_product():
    """Admit and qualify cached runs in memory; leave the public product intact."""
    before = identity()
    product = _assemble_frozen()
    if product.get('schema') != 'terluna.atmosphere.snow-growth-repair/1':
        raise ValueError('Unsupported frozen growth-repair product schema')
    for group, cache in (('cases', frozen.CACHE), ('references', repair.CACHE)):
        for row in product[group]:
            path = cache / (row['id'] + '.json')
            relative = str(path.relative_to(ROOT))
            if relative not in before['raw_caches'] or before['raw_caches'][relative] is None:
                raise RuntimeError('Assembler admitted an unbound cache: ' + relative)
            raw_bytes = path.read_bytes()
            if hashlib.sha256(raw_bytes).hexdigest() != before['raw_caches'][relative]:
                raise RuntimeError('Raw cache changed before qualification: ' + relative)
            raw = json.loads(raw_bytes)
            _qualify(row, raw, reference=group == 'references')
    if identity() != before:
        raise RuntimeError('Postprocessor, assembler or raw cache changed during final export')
    product['reading_rule'] += (
        ' All cases and references are explicitly ineligible for predictive optics. '
        'Reference evidence, reading rules and numerical status are restored from the verified raw caches. '
        'The original threshold summaries are preserved; nominal_endpoint_diagnostic separately reports '
        'final-pivot occupancy and classification at the nominal endpoint within floating-point tolerance.'
    )
    product['assembly'] = dict(
        schema=SCHEMA, stages=['Frozen growth_repair.assemble() source/input admission',
                               'growth_products qualification and nominal-endpoint diagnostics'],
        producer=before, identity_checked_before_and_after=True,
        reading_rule='Simulation producer records remain unchanged. This assembly identity binds the '
                     'postprocessor, frozen assemblers and exact raw caches used by the second stage. '
                     'Null raw-cache hashes identify cases absent during this export. No physical model, '
                     'distribution, original summary value or numerical coefficient is changed.',
    )
    return product


def export():
    product = build_product()
    serialized = json.dumps(product, indent=2, allow_nan=False) + '\n'
    if identity() != product['assembly']['producer']:
        raise RuntimeError('Source or cache identity changed before final product write')
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode='w', dir=OUT.parent, prefix='.growth-repair-',
                                         suffix='.json', delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(serialized)
        os.replace(temporary, OUT)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()
    return product


def main():
    result = export()
    print(json.dumps({'output': str(OUT), 'cases': [r['id'] for r in result['cases']],
                      'references': [r['id'] for r in result['references']],
                      'assembly_schema': result['assembly']['schema']}))


if __name__ == '__main__':
    main()
