"""Rain arriving at the ground: which drops of each bounded cloud-base spectrum reach it, at what size and speed, and
how many per unit of ground rain, in the three target climates beside Earth's air; drips released from leaves; and
centimetre drops leaving the cloud base.

    python -m atmosphere.rain.arrivals        # about two minutes on one thread
    # -> atmosphere/rain/results/rain_arrivals.json

The rain product (product.py, results/rain_drops.json) gives the share of each spectrum's water that reaches the ground
and the median arriving size and speed. A renderer of the rain's sound needs the arriving population itself: drops of
each size per square metre per second, with the speed of each at impact. This product follows the same 64 cloud-base
classes (0.1-6 mm) through the same columns, the Moon's and Earth's air at each target's 2 m state
(product.columns), and keeps the ground's side of the ledger class by class, normalized to 1 mm/h of rain at the
ground. With independent drops the arriving spectrum has the same shape at any rain rate, so a consumer multiplies by
the rate it chooses.

Drips are drops released at rest from a height, as from a leaf, falling through the 2 m air without a cloud above
them. The tail is drops of 8-12 mm leaving the cloud base: larger than the bounded spectra reach, below the Moon's
stability screen of about 16 mm and above Earth's of about 6.6 mm.
"""
from __future__ import annotations

import hashlib
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

from shared.provenance import constants_used
from atmosphere.rain import distributions, drops, product, targets as climate_targets

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
SCHEMA = 'terluna.atmosphere.rain-arrivals/1'
OUT = HERE / 'results' / 'rain_arrivals.json'
CODE = ('atmosphere/rain/arrivals.py', 'atmosphere/rain/product.py', 'atmosphere/rain/drops.py',
        'atmosphere/rain/distributions.py', 'atmosphere/rain/targets.py',
        'atmosphere/radiative_convective/thermodynamics.py')
WORLDS = ('moon', 'earth')
DRIP_DIAMETERS_MM = (1.0, 2.0, 3.0, 4.0, 5.0)
DRIP_HEIGHTS_M = (0.25, 0.5, 1.0, 2.0, 4.0)
TAIL_DIAMETERS_MM = (8.0, 9.0, 10.0, 11.0, 12.0)
MM_H_PER_M_S = distributions.MM_H_PER_M_S

EVIDENCE = ('A runnable calculation of independent single drops through the rain product\'s columns, the same solver '
            'that reproduces 34 Gunn-Kinzer fall speeds on Earth: the arriving populations are the bounded cloud-base '
            'spectra filtered by evaporation below a parcel cloud base in still air, and the drips are short falls '
            'from rest. The spectra and the drips\' sizes and heights are scenario boundary conditions.')
READING_RULE = ('Diameters are volume-equivalent, in mm; speeds are downward at the ground in still air, in m/s. For '
                'each target and world, classes gives each cloud-base size class\'s fall: whether it reaches the '
                'ground, its size and speed there and its time aloft. Each spectrum\'s ground_flux_per_mm_h gives, '
                'class by class, the drops arriving per square metre per second for each mm/h of rain at the ground: '
                'multiply by the ground rain rate. cloud_base_rain_per_ground_rain is the rain leaving the cloud base '
                'that each mm/h at the ground takes, 1 / water_survival_fraction; a spectrum delivering no water has '
                'null fluxes. drips gives, by release diameter and height, the size and speed at the ground of a drop '
                'released at rest; tail gives the falls of 8-12 mm drops from the cloud base, absent where they exceed '
                'the stability screen.')
ASSUMPTIONS = [
    'The columns, airs and spectra are the rain product\'s (results/rain_drops.json): each target\'s 2 m air, the '
    'resort\'s at 301 K, the dry adiabat at the repository\'s heat capacity and fixed vapour below a parcel cloud '
    'base, still air, and the bounded gamma spectra on 64 classes from 0.1 to 6 mm.',
    'Drops fall independently, so the arriving spectrum keeps its shape at any rain rate. The rain product\'s '
    'encounter diagnostic finds that on the Moon every surviving drop meets more than one other on the way down, so '
    'collisions and coalescence would reshape it, and at the rates of a shower the rain\'s own moistening of the air '
    'would let more water through than this column does.',
    'Drips start at rest in the 2 m air and lose at most 2% of their water on the way down; their sizes, heights and '
    'the share of the rain that drips are left to the consumer, pending the canopy\'s interception and storage '
    '(biosphere/canopy).',
    'The tail\'s drops leave the cloud base at their terminal speed, like the spectra\'s; how common they are is not '
    'computed.',
]


def sha256(path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def air_at_ground(column) -> dict:
    """The air at 2 m as the sound of the rain needs it: state, composition and vapour."""
    state = column.state(0.0)
    described = column.describe()
    return dict(gravity_m_s2=column.gravity_m_s2, temperature_k=state['temperature_k'],
                pressure_pa=state['pressure_pa'], relative_humidity=state['relative_humidity'],
                vapour_pressure_pa=state['vapour_pressure_pa'], vapour_mole_fraction=described['vapour_mole_fraction'],
                density_kg_m3=state['density_kg_m3'], dry_fractions=dict(column.air.fractions),
                dry_molar_mass_kg_mol=column.air.molar_mass, cloud_base_m=column.cloud_base_m,
                stable_diameter_mm=drops.stable_diameter(column) * 1e3)


def classes(column):
    """Each cloud-base class's fall to the ground: the grid, and the solver's results."""
    d, width = distributions.diameter_grid(*product.SPECTRUM_GRID_MM)
    runs = [drops.fall(x, column, sample_heights_m=[column.cloud_base_m, 0.0]) for x in d]
    return d, width, runs


def spectra(d, width, runs):
    """For each bounded spectrum: the ground's number flux by class per mm/h at the ground, and its ledger."""
    entry = np.array([r['initial_speed_m_s'] for r in runs])
    alive = np.array([r['survived'] for r in runs])
    surface_d = np.array([r['surface_diameter_m'] if r['survived'] else np.nan for r in runs])
    surface_v = np.array([r['surface_speed_m_s'] if r['survived'] else np.nan for r in runs])
    rows = []
    for dm, mu in product.SPECTRA:
        # One mm/h leaving the cloud base; a surviving class keeps its number flux while it shrinks.
        source = distributions.gamma_flux(d, width, entry, dm_mm=dm, shape_mu=mu, cloud_base_rain_mm_h=1.0)
        ground = distributions.surface_summary(d, source['flux_m2_s'], surface_d, surface_v, alive)
        surviving = ground['water_survival_fraction']
        row = dict(dm_mm=dm, shape_mu=mu, water_survival_fraction=surviving,
                   number_survival_fraction=ground['number_survival_fraction'])
        if surviving:
            flux = np.where(alive, source['flux_m2_s'], 0.0) / surviving
            row.update(cloud_base_rain_per_ground_rain=1.0 / surviving, ground_flux_per_mm_h=flux.tolist(),
                       arrivals_m2_s_per_mm_h=float(flux.sum()),
                       arrival_diameter_mm_by_number=ground['arrival_diameter_mm_by_number'],
                       arrival_speed_m_s_by_number=ground['arrival_speed_m_s_by_number'])
        else:
            row.update(cloud_base_rain_per_ground_rain=None, ground_flux_per_mm_h=None, arrivals_m2_s_per_mm_h=None,
                       arrival_diameter_mm_by_number=None, arrival_speed_m_s_by_number=None)
        rows.append(row)
    return rows


def drips(column):
    """Drops released at rest from each height in the 2 m air: size, speed and time at the ground."""
    rows = []
    for height in DRIP_HEIGHTS_M:
        released = replace(column, cloud_base_m=height)
        for mm in DRIP_DIAMETERS_MM:
            r = drops.fall(mm * 1e-3, released, sample_heights_m=[height, 0.0], initial_speed_m_s=0.0)
            rows.append(dict(release_diameter_mm=mm, release_height_m=height, status=r['status'],
                             surface_diameter_mm=r['surface_diameter_m'] * 1e3 if r['survived'] else None,
                             impact_speed_m_s=r['surface_speed_m_s'] if r['survived'] else None,
                             terminal_speed_m_s=r['surface_terminal_speed_m_s'] if r['survived'] else None,
                             fall_time_s=r['fall_time_s']))
    return rows


def tail(column):
    """Drops of 8-12 mm leaving the cloud base at their terminal speed, where the stability screen allows them."""
    rows = []
    screen = drops.stable_diameter(column, column.cloud_base_m)
    for mm in TAIL_DIAMETERS_MM:
        if mm * 1e-3 > screen:
            rows.append(dict(initial_diameter_mm=mm, status='above the stability screen'))
            continue
        r = drops.fall(mm * 1e-3, column, sample_heights_m=[column.cloud_base_m, 0.0])
        rows.append(dict(initial_diameter_mm=mm, status=r['status'],
                         surviving_water_fraction=1.0 - r['evaporated_mass_fraction'] if r['survived'] else 0.0,
                         surface_diameter_mm=r['surface_diameter_m'] * 1e3 if r['survived'] else None,
                         impact_speed_m_s=r['surface_speed_m_s'] if r['survived'] else None,
                         fall_time_s=r['fall_time_s']))
    return rows


def world(column) -> dict:
    d, width, runs = classes(column)
    return dict(
        air=air_at_ground(column),
        classes=dict(cloud_base_diameter_mm=(d * 1e3).tolist(), width_mm=(width * 1e3).tolist(),
                     status=[r['status'] for r in runs],
                     surface_diameter_mm=[r['surface_diameter_m'] * 1e3 if r['survived'] else None for r in runs],
                     impact_speed_m_s=[r['surface_speed_m_s'] if r['survived'] else None for r in runs],
                     leaving_speed_m_s=[r['initial_speed_m_s'] for r in runs],
                     fall_time_s=[r['fall_time_s'] for r in runs]),
        spectra=spectra(d, width, runs), drips=drips(column), tail=tail(column))


def build() -> dict:
    data, binding = climate_targets.load()
    air = product.target_air(climate_targets.read(data))
    out = dict(schema=SCHEMA, evidence=EVIDENCE, reading_rule=READING_RULE, assumptions=ASSUMPTIONS,
               grid_mm=list(product.SPECTRUM_GRID_MM),
               spectra=[dict(dm_mm=dm, shape_mu=mu) for dm, mu in product.SPECTRA],
               central_spectrum=dict(dm_mm=product.CENTRAL_SPECTRUM[0], shape_mu=product.CENTRAL_SPECTRUM[1]),
               targets={})
    for key in product.TARGETS:
        columns = product.columns(air[key])
        out['targets'][key] = dict(name_in_targets=air[key]['name_in_targets'], held_in=air[key]['held_in'],
                                   **{w: world(columns[w]) for w in WORLDS})
    out['producer'] = dict(domain='atmosphere', files={p: sha256(ROOT / p) for p in CODE}, **binding,
                           constants=constants_used([ROOT / p for p in CODE]))
    return out


def main(argv=None) -> int:
    result = build()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(result, indent=1, allow_nan=False) + '\n')
    central = (product.CENTRAL_SPECTRUM[0], product.CENTRAL_SPECTRUM[1])
    for key in product.TARGETS:
        line = []
        for w in WORLDS:
            row = next(s for s in result['targets'][key][w]['spectra'] if (s['dm_mm'], s['shape_mu']) == central)
            line.append(f"{w} {row['arrivals_m2_s_per_mm_h']:7.1f} drops/m2/s per mm/h, "
                        f"median {row['arrival_diameter_mm_by_number']['p50']:.2f} mm "
                        f"at {row['arrival_speed_m_s_by_number']['p50']:.2f} m/s")
        print(f'{key:9s} ' + '; '.join(line))
    print(f'wrote {OUT.relative_to(ROOT)}')
    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
