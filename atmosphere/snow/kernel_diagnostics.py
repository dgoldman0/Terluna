"""Inspect the inherited collision kernel; this does not infer a physical snow limit.

The fixed-geometry calculation exposes how sweepout changes with collector mass
against a fixed small partner. A mass-exponent-2 scenario is matched to the
original mass at 1 cm; it is a diagnostic, not an adopted geometry or a fracture
law. See geometry_literature.md for the source and applicability of the scaling.
"""
from __future__ import annotations
from dataclasses import replace
from pathlib import Path
import hashlib
import json
import math
import numpy as np
from atmosphere.snow import aggregation as ag, flakes as sf
from atmosphere.snow.habits import GEOMETRIES
from shared.provenance import constants_used

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'atmosphere/snow/results/kernel_diagnostics.json'
SCHEMA = 'terluna.atmosphere.snow-kernel-diagnostics/1'


def encounter_kernel(diameter_m, partner_diameter_m, geometry, gravity, state):
    d = np.asarray(diameter_m, dtype=float)
    partner = np.full_like(d, partner_diameter_m)
    v = ag.speeds(d, geometry.mass(d), gravity, state, geometry.area_ratio)[0]
    vp = ag.speeds(partner, geometry.mass(partner), gravity, state, geometry.area_ratio)[0]
    kernel = geometry.area_ratio * math.pi / 4 * (d + partner)**2 * np.abs(v-vp)
    return kernel, v, vp


def local_mass_slope(diameter_m, partner_diameter_m, geometry, gravity, state, log_step=1e-4):
    d = np.asarray(diameter_m, dtype=float)
    lo = d * np.exp(-log_step / geometry.mass_exponent)
    hi = d * np.exp(log_step / geometry.mass_exponent)
    klo = encounter_kernel(lo, partner_diameter_m, geometry, gravity, state)[0]
    khi = encounter_kernel(hi, partner_diameter_m, geometry, gravity, state)[0]
    if np.any(klo <= 0) or np.any(khi <= 0):
        raise ValueError('Log slope is undefined at a zero encounter kernel')
    return np.log(khi/klo)/(2*log_step)


def analyse():
    world, binding = sf.moon_column()
    base = GEOMETRIES['bf95_dmax_dendrites_area']
    anchor = .01
    matched = replace(base, id='mass_exponent_2_matched_at_1_cm',
                      mass_coefficient=float(base.mass(anchor)/anchor**2), mass_exponent=2.,
                      evidence='Diagnostic constant exponent 2 matched to baseline mass at 1 cm; '
                               'not an observed lunar geometry or accepted physical repair.')
    rows = []
    slope_check = 0.
    for tc in (-15., -10., -2.):
        state = world.state(world.height_of(sf.T0+tc))
        for geometry in (base, matched):
            for partner in (.0002, .002):
                d = np.array([.01, .03, .1, .15, .3, .6, 1., 10., 100., 1000.])
                kernel, speed, small_speed = encounter_kernel(d, partner, geometry, world.gravity, state)
                slope = local_mass_slope(d, partner, geometry, world.gravity, state)
                refined = local_mass_slope(d, partner, geometry, world.gravity, state, 5e-5)
                slope_check = max(slope_check, float(np.max(np.abs(refined-slope))))
                rows.append(dict(temperature_c=tc, geometry=geometry.to_dict(),
                    partner_diameter_m=partner, diameter_m=d.tolist(),
                    collector_speed_m_s=speed.tolist(), partner_speed_m_s=float(small_speed[0]),
                    encounter_kernel_m3_s=kernel.tolist(), local_dlog_kernel_dlog_mass=slope.tolist()))
    assert slope_check < 1e-7
    files = ['atmosphere/snow/kernel_diagnostics.py','atmosphere/snow/aggregation.py',
             'atmosphere/snow/habits.py','atmosphere/snow/flakes.py','atmosphere/rain/drops.py',
             'atmosphere/rain/targets.py','atmosphere/radiative_convective/thermodynamics.py']
    fragment_mass=float(base.mass(ag.FRAGMENT_M))
    chip_cap=fragment_mass*ag.BREAKUP['most']
    chipping=[dict(diameter_m=d,mass_kg=float(base.mass(d)),
                   maximum_chipped_mass_kg=min(chip_cap,.5*float(base.mass(d))),
                   maximum_chipped_fraction=min(chip_cap/float(base.mass(d)),.5))
              for d in (.01,.15,.6)]
    return dict(schema=SCHEMA, rows=rows,
        chipping_cap=dict(fragment_diameter_m=ag.FRAGMENT_M, maximum_fragments=ag.BREAKUP['most'],
                          fragment_mass_kg=fragment_mass, examples=chipping,
                          interpretation='Upper bound per nonsticking encounter when the stated aggregate is '
                          'the smaller parent. The larger parent loses no mass in this closure. '
                          'This is a code-law diagnostic, not a measured mechanical survival threshold.'),
        checks=dict(log_derivative_step_halving_max_absolute_difference=slope_check,
                    anchor_mass_difference_kg=float(matched.mass(anchor)-base.mass(anchor))),
        evidence='Calculation of the inherited fixed-geometry kernel. Local slopes above one '
                 'identify superlinear sweepout against a specified partner; they neither prove '
                 'gelation in the finite, settling, fragmenting column nor isolate the cause of '
                 'the old stochastic runaway. Large dimensions are diagnostic extrapolations.',
        reading_rule='K is encounter volume per time, before temperature-only sticking efficiency. '
                     'Temperature is fixed within each row. dlogK/dlogm holds the partner fixed. '
                     'The exponent-2 scenario holds partner diameter fixed, which also changes its mass and speed. '
                     'Its finite-size encounter rates need not decrease. It is separate from numerical repair.',
        theory=dict(source='https://doi.org/10.1103/PhysRevE.70.021403',
                    scope='Westbrook et al. 2004 section IV: morphology feedback in an ideal '
                          'irreversible aggregation model; conditional guidance for kernel diagnosis.'),
        producer=dict(files={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
                      inputs=binding['inputs'], constants=constants_used(files)))


def main():
    result = analyse()
    OUT.write_text(json.dumps(result, indent=2, default=float)+'\n')
    print(f'Wrote {OUT.relative_to(ROOT)}; finite-difference check {result["checks"]}')


if __name__ == '__main__':
    main()
