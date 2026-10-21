"""One declared aggregate-geometry sensitivity on the baseline lunar mass grid.

The mass law continues quadratically above 1 cm. This is a physical-closure
sensitivity, not an adopted shape model or a calibration of giant snowflakes.
All baseline mass pivots, source forcing and numerical operators are retained.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np

from atmosphere.snow import aggregation as ag, flakes as sf
from atmosphere.snow.column_experiments import Kinetics, QuadraticContinuation, solve_column
from atmosphere.snow.habits import GEOMETRIES
from atmosphere.snow.sectional import seed_distribution
from shared.provenance import constants_used

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'atmosphere/snow/results/geometry_continuation.json'
SCHEMA = 'terluna.atmosphere.snow-geometry-continuation/1'
BASE = GEOMETRIES['bf95_dmax_dendrites_area']
GEOMETRY = QuadraticContinuation(BASE)
SPECIFICATION = dict(
    world='moon', top_c=-15., number_per_m3=1e5, snowfall_mm_h=3.,
    layers=30, mass_bins=96, baseline_minimum_diameter_m=2e-5,
    baseline_maximum_diameter_m=10., max_step_s=30., cfl=.3,
    fill=1.5, gather=1.5, breakup=True, growth_scheme='limited',
    batch_growth=True, kinetics='legacy', domain_mass_tolerance=1e-6,
    anchor_m=.01,
)
FILES = (
    'atmosphere/snow/geometry_continuation.py',
    'atmosphere/snow/column_experiments.py',
    'atmosphere/snow/limited_growth.py',
    'atmosphere/snow/limited_growth_batch.py',
    'atmosphere/snow/sectional.py',
    'atmosphere/snow/aggregation.py',
    'atmosphere/snow/habits.py',
    'atmosphere/snow/flakes.py',
    'atmosphere/rain/drops.py',
    'atmosphere/rain/targets.py',
    'atmosphere/radiative_convective/thermodynamics.py',
)


def identity(binding):
    """Bind sources, loaded column inputs and named shared constants."""
    def sha(path):
        return hashlib.sha256(Path(path).read_bytes()).hexdigest()
    inputs = {p: sha(ROOT / p) for p in binding.get('inputs', {})}
    if inputs != binding.get('inputs', {}):
        raise RuntimeError('Loaded column input identity changed')
    return dict(files={p: sha(ROOT / p) for p in FILES}, inputs=inputs,
                constants=constants_used(FILES))


def baseline_mass_grid():
    """Exactly the 96 mass pivots used by moon_limited_m96."""
    return np.geomspace(
        float(BASE.mass(SPECIFICATION['baseline_minimum_diameter_m'])),
        float(BASE.mass(SPECIFICATION['baseline_maximum_diameter_m'])),
        SPECIFICATION['mass_bins'],
    )


def column_arguments(world):
    masses = baseline_mass_grid()
    seed_d = math.sqrt(ag.START_M[0] * ag.START_M[1])
    top = world.height_of(sf.T0 + SPECIFICATION['top_c'])
    reference_speed = float(Kinetics('legacy')(
        np.array([BASE.mass(seed_d)]), np.array([seed_d]), BASE,
        world.gravity, world.state(top),
    )[0][0])
    keys = ('top_c', 'number_per_m3', 'snowfall_mm_h', 'layers', 'mass_bins',
            'max_step_s', 'cfl', 'fill', 'gather', 'breakup',
            'growth_scheme', 'batch_growth', 'domain_mass_tolerance')
    args = {k: SPECIFICATION[k] for k in keys}
    args.update(geometry=GEOMETRY, kinetics=Kinetics('legacy'), mass_grid=masses,
                minimum_diameter_m=float(GEOMETRY.size(masses[0])),
                maximum_diameter_m=float(GEOMETRY.size(masses[-1])),
                injection_velocity_m_s=reference_speed)
    return args


def export(raw, producer, wall_seconds):
    """Preserve native arrays and budgets, with a numerical-domain reading gate."""
    masses = np.asarray(raw['masses'])
    diameters = np.asarray(raw['diameters_m'])
    arrivals = np.asarray(raw['arriving_numbers_m2'])
    weights = masses * arrivals
    total = float(weights.sum())
    cdf = np.cumsum(weights) / total if total else np.zeros_like(weights)

    def quantile(p):
        return float(diameters[min(np.searchsorted(cdf, p), diameters.size-1)]) if total else None

    adequate = bool(raw['domain_adequacy']['passed'])
    native = {key: np.asarray(raw[key]).tolist() for key in (
        'masses', 'diameters_m', 'arriving_numbers_m2',
        'time_mean_number_m3', 'final_number_m3',
        'layer_height_m', 'layer_edges_m', 'velocity_m_s',
    )}
    native['masses_kg'] = native.pop('masses')
    return dict(
        schema=SCHEMA, status='completed', id='moon_quadratic_above_1cm_limited_m96',
        specification=SPECIFICATION, producer=producer,
        geometry=GEOMETRY.to_dict(), settings=raw['settings'],
        comparison=dict(
            baseline_case='moon_limited_m96',
            baseline_geometry=BASE.to_dict(),
            mass_grid_sha256=hashlib.sha256(masses.astype('<f8').tobytes()).hexdigest(),
            grid_rule='Identical mass pivots; physical diameter differs above the 1 cm join.',
            baseline_maximum_resolved_diameter_m=float(BASE.size(masses[-1])),
            actual_maximum_resolved_diameter_m=float(diameters[-1]),
            actual_upper_absorbing_diameter_m=raw['domain_adequacy']['upper_absorbing_pivot_m'],
            seed_numbers_per_injected_particle=seed_distribution(masses, BASE).tolist(),
            source_rule='The baseline cloud-top seed velocity fixes source flux and integration duration.',
        ),
        summary=dict(
            snowfall_mm_h=total * 3600. / raw['seconds'],
            arriving_number_flux_m2_s=float(arrivals.sum() / raw['seconds']),
            mass_median_m=quantile(.5), mass_d90_m=quantile(.9), mass_d99_m=quantile(.99),
            arriving_mass_share_above={str(x): float(weights[diameters >= x].sum()/total) if total else 0.
                                      for x in (.005, .01, .15, .6, 1.)},
        ),
        mass_ledger_kg_m2=raw['mass_ledger_kg_m2'],
        extensive_mass_ledger_kg=raw['extensive_mass_ledger_kg'],
        domain_exit_numbers_m2=raw['domain_exit_numbers_m2'],
        gather_domain_exit_mass_kg_m2=raw['gather_domain_exit_mass_kg_m2'],
        domain_adequacy=raw['domain_adequacy'],
        diagnostics=dict(raw['diagnostics'], wall_seconds=wall_seconds),
        gathering_seconds=raw['seconds'], crossing_hours=raw['crossing_hours'],
        native_arrays=native,
        evidence='One uncalibrated physical-closure sensitivity at fixed numerical resolution. '
                 'The 1 cm join is a chosen diagnostic parameter, not a fitted transition.',
        interpretation=(
            'Finite mass-domain check passed; mass/time/height convergence and physical validation remain unperformed.'
            if adequate else
            'Finite mass-domain check failed. Surviving fluxes and quantiles are diagnostics, not snow or optics forecasts.'
        ),
        optical_prediction_eligible=False,
        remaining=['Independent mass-grid, timestep, height and mass-domain controls.',
                   'Evolving geometry, interception applicability, whole-aggregate fracture and moist trajectories.',
                   'Physical calibration of the continuation and of any flake-size distribution.'],
    )


def run(*, progress=None):
    world, binding = sf.moon_column()
    before = identity(binding)
    args = column_arguments(world)
    start = time.monotonic()
    raw = solve_column(world, progress=progress, **args)
    if identity(binding) != before:
        raise RuntimeError('Producer changed during run; no result admitted')
    return export(raw, before, time.monotonic() - start)


def main():
    start, next_report = time.monotonic(), 0.

    def progress(t, steps):
        nonlocal next_report
        if t >= next_report:
            print(f'quadratic m96: {t/3600:.1f} model h, {steps} steps, '
                  f'{time.monotonic()-start:.1f} wall s', flush=True)
            next_report = t + 6*3600.

    result = run(progress=progress)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    print(json.dumps({k: result[k] for k in ('id', 'summary', 'domain_adequacy',
                                           'mass_ledger_kg_m2')}), flush=True)


if __name__ == '__main__':
    main()
