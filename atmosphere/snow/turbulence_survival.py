"""How large a lunar snowflake can be and still carry no more load than Earth's largest observed flakes.

Plan step 5.  No measured strength law transfers to decimetre aggregates
(strength_literature.md), so the estimate compares loads instead: a lunar flake
that bears no more load at its weakest junction than an Earth flake that is
observed to survive should survive too, if the two are built alike from ice
crystals joined by similar contacts.  Two loads act on a falling flake.  Its
weight, m g, is carried by drag at terminal speed and passes through the
junctions; turbulent velocity differences across the flake add a fluctuating
drag, which for a small velocity increment du is 2 m g du / V.  Earth's flakes
of 4-6 cm survive a 9 mm/h sheared storm (Lawson et al. 1993, 1998), and flakes
of a centimetre or two are ordinary.  This module finds, for each such Earth
flake, the lunar span with the same peak load, over a range of turbulence.
It is a scaling estimate by analogy, not a strength calculation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from atmosphere.snow import flakes

HERE = Path(__file__).resolve().parent
PRODUCT = HERE / 'results/turbulence_survival.json'
# Aggregate mass laws, SI: McSnow's aggregates of side planes, columns and bullets (Mitchell 1996; the
# step-4 geometry), Leinonen's dendrite aggregates (as quoted in McSnow), and Brown and Francis (1995).
MASS_LAWS = {'mitchell1996_mixed_aggregates': (0.0028 * 10 ** (2 * 2.1 - 3), 2.1),
             'leinonen_dendrite_aggregates': (0.01243, 2.0),
             'brown_francis1995': (flakes.MASS_COEFFICIENT, flakes.MASS_EXPONENT)}
EARTH_FLAKES_M = (.02, .04, .06)          # ordinary large flakes; the 4-6 cm CASP II aggregates
DISSIPATION = (1e-5, 1e-4, 1e-3, 1e-2)     # m2 s-3; turbulence.md's range
DISSIPATION_RATIO = (.1, 1., 10.)          # lunar over Earth dissipation at the same flake


def velocity_increment(d, eps, nu):
    eta = (nu ** 3 / eps) ** .25
    return (eps * d) ** (1 / 3) if d >= eta else d * math.sqrt(eps / nu)


def load(column, d, eps, law):
    """Peak junction load proxy (N): weight plus the linearized turbulent drag fluctuation."""
    a, b = MASS_LAWS[law]
    s = column.state(column.z_freeze)
    m = a * d ** b
    v, _ = flakes.fall_speed(d, m, column.gravity, s)
    du = velocity_increment(d, eps, s['eta'] / s['rho'])
    w = m * column.gravity
    return dict(weight_n=w, turbulent_n=2 * w * du / v, total_n=w * (1 + 2 * du / v), fall_speed_m_s=v,
                velocity_increment_m_s=du, mass_kg=m)


def equivalent_span(moon, earth, d_earth, eps_earth, eps_moon, law, key='total_n'):
    target = load(earth, d_earth, eps_earth, law)[key]
    lo, hi = d_earth * .1, d_earth * 100.
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if load(moon, mid, eps_moon, law)[key] < target:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


def run():
    moon, binding = flakes.moon_column()
    earth = flakes.earth_column()
    rows = []
    for law in MASS_LAWS:
        for d in EARTH_FLAKES_M:
            for eps in DISSIPATION:
                for ratio in DISSIPATION_RATIO:
                    e_m = load(earth, d, eps, law)
                    dm = equivalent_span(moon, earth, d, eps, eps * ratio, law)
                    m_m = load(moon, dm, eps * ratio, law)
                    rows.append(dict(mass_law=law, earth_span_m=d, earth_dissipation=eps, moon_over_earth_dissipation=ratio,
                                     moon_equivalent_span_m=dm, ratio=dm / d,
                                     settling_only_span_m=equivalent_span(moon, earth, d, eps, eps * ratio, law, 'weight_n'),
                                     earth=e_m, moon=m_m))
    s_e, s_m = earth.state(earth.z_freeze), moon.state(moon.z_freeze)
    product = dict(
        schema='terluna.atmosphere.snow-turbulence-survival/1',
        evidence=('Scaling estimate by load equivalence with observed Earth aggregates; no measured strength of a '
                  'decimetre aggregate, no resolved lunar turbulence.'),
        air=dict(earth=dict(rho=s_e['rho'], p=s_e['p'], t=s_e['t']), moon=dict(rho=s_m['rho'], p=s_m['p'], t=s_m['t'])),
        gravity=dict(earth=earth.gravity, moon=moon.gravity), rows=rows,
        assumptions=['The weakest junction carries a load that grows with the flake, and its strength does not depend '
                     'on the flake\'s size; lunar and Earth flakes are built from similar crystals and contacts.',
                     'The largest observed Earth flakes are taken to sit at their strength limit; growth time or water '
                     'may limit them instead, which would make the lunar spans conservative.',
                     'Lunar flakes spend one to four days aloft: contacts sinter and strengthen (Demmenie et al. 2025), '
                     'while repeated sublimation can loosen them (Oraltay and Hallett 1989); neither is quantified.',
                     'The turbulent fluctuation is linearized in du/V; a lunar dissipation rate is unmeasured, so it is '
                     'varied from a tenth to ten times Earth\'s.'],
        inputs=dict(climate_binding=binding.get('inputs', {})),
        producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    PRODUCT.write_text(json.dumps(product, separators=(',', ':'), allow_nan=False) + '\n')
    return product


def main():
    argparse.ArgumentParser(description=__doc__).parse_args()
    product = run()
    for r in product['rows']:
        if r['mass_law'] == 'mitchell1996_mixed_aggregates' and r['moon_over_earth_dissipation'] == 1.:
            print(f"Earth {r['earth_span_m']*100:3.0f} cm, eps {r['earth_dissipation']:7.0e}: lunar equivalent "
                  f"{r['moon_equivalent_span_m']*100:5.1f} cm (x{r['ratio']:.2f}); settling only "
                  f"{r['settling_only_span_m']*100:5.1f} cm; du/V Earth {r['earth']['velocity_increment_m_s']/r['earth']['fall_speed_m_s']:.3f}"
                  f" Moon {r['moon']['velocity_increment_m_s']/r['moon']['fall_speed_m_s']:.3f}")


if __name__ == '__main__':
    main()
