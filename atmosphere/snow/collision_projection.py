"""Measure fixed-pivot spreading in isolated small-on-large sticking events.

This diagnostic evaluates the frozen CollisionOperator on the actual lunar
mass grids. A unit kernel isolates mass allocation; no column is integrated.
It does not repeat the physical kernel slopes in kernel_diagnostics.py.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

from atmosphere.snow.habits import GEOMETRIES
from atmosphere.snow.sectional import CollisionOperator
from shared.provenance import constants_used

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'atmosphere/snow/results/collision_projection.json'
SCHEMA = 'terluna.atmosphere.snow-collision-projection/1'
GEOMETRY = GEOMETRIES['bf95_dmax_dendrites_area']
FILES = ('atmosphere/snow/collision_projection.py', 'atmosphere/snow/sectional.py',
         'atmosphere/snow/habits.py', 'atmosphere/snow/aggregation.py',
         'atmosphere/snow/flakes.py')


def identity():
    return dict(files={p: hashlib.sha256((ROOT/p).read_bytes()).hexdigest() for p in FILES},
                inputs={}, constants=constants_used(FILES))


def probe(masses, small_index, large_index):
    """Compare one implemented reaction with its exact unprojected moments."""
    masses = np.asarray(masses, dtype=float)
    i, j = small_index, large_index
    if not 0 <= i < j < len(masses)-1:
        raise ValueError('Distinct interior parents, ordered small to large, required')
    small, large = masses[i], masses[j]
    spacing = masses[j+1]-large
    if not small < spacing:
        raise ValueError('This diagnostic requires a product below the next collector pivot')
    kernel = np.zeros((masses.size, masses.size))
    kernel[i, j] = kernel[j, i] = 1.
    operator = CollisionOperator(masses, kernel, sticking=1.)
    number = np.zeros_like(masses)
    number[i] = number[j] = 1.
    change, exit_number, exit_mass = operator.evaluate(number)
    if np.any(exit_number) or np.any(exit_mass):
        raise AssertionError('The isolated interior event must not have domain exits')
    exact_second_moment_change = 2*small*large
    predicted_excess = small*(spacing-small)
    actual_second_moment_change = float(np.dot(change, masses*masses))
    mass_scale = float(np.dot(np.abs(change), masses))
    active = np.flatnonzero(change)
    return dict(
        small_index=int(i), large_index=int(j),
        small_mass_kg=float(small), collector_mass_kg=float(large),
        collector_next_pivot_mass_kg=float(masses[j+1]),
        product_mass_kg=float(small+large),
        number_change_per_event=float(change.sum()),
        relative_first_moment_residual=float(np.dot(change, masses)/mass_scale),
        exact_second_moment_change_kg2=float(exact_second_moment_change),
        projected_second_moment_change_kg2=actual_second_moment_change,
        excess_second_moment_change_kg2=float(actual_second_moment_change-exact_second_moment_change),
        analytical_excess_second_moment_kg2=float(predicted_excess),
        excess_formula_relative_residual=float(
            (actual_second_moment_change-exact_second_moment_change-predicted_excess)/predicted_excess),
        excess_second_moment_fraction=float(
            (actual_second_moment_change-exact_second_moment_change)/exact_second_moment_change),
        upper_pivot_allocation_number_fraction=float(small/spacing),
        collector_jump_second_moment_ratio=float(spacing/small),
        tendency=dict(indices=active.tolist(), number_change_per_event=change[active].tolist()),
    )


def analyse():
    before = identity()
    rows = []
    for bins in (96, 144, 216):
        masses = np.geomspace(float(GEOMETRY.mass(2e-5)), float(GEOMETRY.mass(10.)), bins)
        diameters = GEOMETRY.size(masses)
        i = int(np.argmin(np.abs(np.log(diameters/.0002))))
        for requested in (.01, .1, .15, .6, 1.):
            j = int(np.argmin(np.abs(np.log(diameters/requested))))
            row = probe(masses, i, j)
            row.update(mass_bins=bins, geometric_mass_ratio=float(masses[1]/masses[0]),
                       requested_small_diameter_m=.0002, requested_collector_diameter_m=requested,
                       actual_small_diameter_m=float(diameters[i]),
                       actual_collector_diameter_m=float(diameters[j]),
                       collector_next_pivot_diameter_m=float(diameters[j+1]))
            rows.append(row)
    checks = dict(
        maximum_first_moment_relative_residual=max(abs(r['relative_first_moment_residual']) for r in rows),
        maximum_number_change_error=max(abs(r['number_change_per_event']+1.) for r in rows),
        maximum_excess_formula_relative_residual=max(abs(r['excess_formula_relative_residual']) for r in rows),
        all_second_moment_excesses_positive=all(r['excess_second_moment_fraction']>0 for r in rows),
    )
    if (checks['maximum_first_moment_relative_residual'] > 1e-13
            or checks['maximum_number_change_error'] > 1e-13
            or checks['maximum_excess_formula_relative_residual'] > 1e-11
            or not checks['all_second_moment_excesses_positive']):
        raise AssertionError('Collision projection diagnostic failed its exact-moment checks')
    if identity() != before:
        raise RuntimeError('Diagnostic source changed during calculation; no result admitted')
    return dict(
        schema=SCHEMA, producer=before,
        settings=dict(geometry=GEOMETRY.to_dict(), mass_bins=[96, 144, 216],
                      minimum_diameter_m=2e-5, maximum_diameter_m=10.,
                      kernel_m3_s=1., each_parent_number_m3=1., sticking=1., breakup=False,
                      event_rate_per_m3_s=1., time_integration=False),
        derivation=dict(
            notation='A small mass s sticks to collector pivot M; the next collector pivot is M+delta, with 0<s<delta.',
            allocation='The product is represented by number fractions 1-s/delta at M and s/delta at M+delta.',
            exact_second_moment='(M+s)^2 - M^2 - s^2 = 2*M*s.',
            extra_second_moment='The interpolation adds s*(delta-s) per sticking event; number and first mass moment still close.',
            relative_excess='Extra / exact = (delta-s)/(2*M), approaching (r-1)/2 for a geometric ratio r and s/M tending to zero.',
            tagged_collector_interpretation='With a prescribed small-particle bath, the physical collector jump has second moment s^2 per encounter. '
                'The projected jump has second moment s*delta; their ratio is delta/s. This is an instantaneous generator comparison, '
                'not an integrated collector variance or a prediction of a coupled tail.',
        ),
        rows=rows, checks=checks,
        evidence='Direct evaluation of the frozen collision operator and exact single-event moment arithmetic. '
                 'The selected pairs use the nearest actual mass pivots to the requested diameters.',
        reading_rule='Positive second-moment excess is numerical redistribution, not physical fragmentation. '
                     'These local values quantify an unresolved collision-grid error that the limited deposition repair does not alter. '
                     'They do not determine how much of the full-column upper-domain loss it causes.',
        relationship_to_tests='Existing constant/additive-kernel tests verify independent coagulation solutions and time accuracy. '
            'The geometric-grid constant-kernel test already demonstrates second-moment error decreasing under refinement. '
            'No independent continuous-mass or sufficiently resolved reference has yet quantified the same error under the full snow kernel.',
        unperformed=['Time-integrated collision-only comparison under the inherited snow kernel.',
                     'Independent coupled-column mass/time/height/domain convergence.',
                     'Quantification of fragmentation projection and unresolved velocity differences within one mass cell.'],
    )


def main():
    result = analyse()
    OUT.write_text(json.dumps(result, indent=2, allow_nan=False)+'\n')
    print(json.dumps(result['checks']))


if __name__ == '__main__':
    main()
