"""The protection's design point: the ultraviolet transmission a swarm can hold, and what it gives the air.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.protection_architecture.run   # results in studies/protection_architecture/results

Reads the swarm's ultraviolet transmission (protection/transmission, product
terluna.protection.uv-transmission/1) for three levels of design choices and
holes 1 and 4 times the particle diameter. At each transmission it solves the
loss response's thermal column (atmosphere/loss_response) behind each shield,
for each treatment of the upper air, at quiet Sun and at solar maximum, and
reports:

- the ultraviolet-driven loss, with the solar-wind screening range (low,
  central, high) added;
- the atmospheric cycle time, the atmosphere's mass over that loss: how long the
  loss and the resupply that matches it take to replace the whole atmosphere;
- the protected radius at which ultraviolet deposited in the thermosphere
  outside the shadow is a tenth of all the heat the upper air receives
  (transmitted sunlight, the film's X-rays and the sky's Lyman-alpha glow,
  expressed as an equivalent transmission), and its aperture against the
  September design's;
- the ions the sunlit exosphere makes outside that radius, counted to a third of
  the Hill radius.

The atmosphere's mass is the feasibility baseline's hydrostatic column
(research/baselines/feasibility/reference.json). As in the loss response, the
loss leaves out Earth's tidal lowering of the escape barrier and the escape of
the exosphere's ions, and atomic oxygen, photochemistry and hydrogen.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys

from atmosphere.loss_response import absorption as ab
from atmosphere.loss_response import model as lr
from atmosphere.middle_atmosphere import escape

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.protection-design-point/1'
TRANSMISSION = ROOT / 'protection' / 'transmission' / 'results' / 'uv_transmission.json'
TRANSMISSION_SCHEMA = 'terluna.protection.uv-transmission/1'
BASELINE = ROOT / 'research' / 'baselines' / 'feasibility' / 'reference.json'
YEAR_S = 365.25 * 86400.0
HOLE_FACTORS = ('1', '4')
WIND_LEVELS = ('low', 'central', 'high')


def atmosphere_mass_kg():
    return float(json.loads(BASELINE.read_text())['baseline']['mass_kg'])


def cycle_time_years(loss_kg_s, mass_kg):
    return mass_kg / loss_kg_s / YEAR_S


def point(shield, treatment, activity, transmission, weights, wind, mass_kg):
    row = lr.euv_sweep(shield, treatment, activity, glow=True, leaks=(transmission,))[0]
    uv = lr.loss_estimate(row)
    unit = escape.leakage_heat(transmission=1.0, activity=lr.ACTIVITY[activity]['uv'])
    equivalent = row['deposited_heat_w_m2'] / unit
    case = ab.case(shield, treatment, activity, transmission, uv, weights, reference=equivalent)
    heat = case['heating']['0.1']
    totals = {level: uv + wind[level]['total_kg_s'] for level in WIND_LEVELS}
    return dict(shield=shield, treatment=treatment, activity=activity, transmission=transmission,
                status=row['status'], heating_as_equivalent_transmission=equivalent,
                glow_share_of_heating=row['glow_heat_w_m2'] / row['deposited_heat_w_m2'],
                uv_driven_loss_kg_s=uv, total_loss_kg_s=totals,
                cycle_time_years={k: cycle_time_years(v, mass_kg) for k, v in totals.items()},
                exobase_radius_R=case['exobase_radius_R'], exobase_temperature_k=case['exobase_temperature_k'],
                protected_radius_R=heat['radius_R'], aperture=heat['aperture'],
                exosphere_ion_production_kg_s=case['exosphere_ion_production_kg_s']['third_hill']['heating_radius'])


def span(rows, get):
    values = [get(r) for r in rows]
    return [min(values), max(values)]


def summarise(points):
    out = {}
    for level in sorted({p['level'] for p in points}):
        out[level] = {}
        for shield in lr.SHIELDS:
            rows = [p for p in points if p['level'] == level and p['shield'] == shield]
            out[level][shield] = dict(
                transmission=span(rows, lambda r: r['transmission']),
                uv_driven_loss_kg_s=span(rows, lambda r: r['uv_driven_loss_kg_s']),
                total_loss_kg_s={w: span(rows, lambda r, w=w: r['total_loss_kg_s'][w]) for w in WIND_LEVELS},
                cycle_time_years={w: span(rows, lambda r, w=w: r['cycle_time_years'][w]) for w in WIND_LEVELS},
                glow_share_of_heating=span(rows, lambda r: r['glow_share_of_heating']),
                exobase_radius_R=span(rows, lambda r: r['exobase_radius_R']),
                protected_radius_R=span(rows, lambda r: r['protected_radius_R']),
                aperture_area_over_september=span(rows, lambda r: r['aperture']['area_over_september']),
                exosphere_ion_production_kg_s=span(rows, lambda r: r['exosphere_ion_production_kg_s']))
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('--out', type=Path, default=HERE / 'results')
    args = parser.parse_args(argv)
    product = json.loads(TRANSMISSION.read_text())
    if product['schema'] != TRANSMISSION_SCHEMA:
        raise SystemExit(f"unexpected transmission schema {product['schema']}")
    mass = atmosphere_mass_kg()
    wind = lr.solar_wind_losses()
    weights = ab.band_weights()
    points = []
    for level, lv in product['levels'].items():
        for k in HOLE_FACTORS:
            t = lv['transmission_by_hole_factor'][k]
            for shield in lr.SHIELDS:
                for treatment in lr.TREATMENTS:
                    for activity in lr.ACTIVITY:
                        points.append(dict(point(shield, treatment, activity, t, weights, wind, mass),
                                           level=level, hole_factor=float(k)))
    summary = summarise(points)
    for level, by_shield in summary.items():
        for shield, s in by_shield.items():
            lo, hi = s['total_loss_kg_s']['low'][0], s['total_loss_kg_s']['high'][1]
            c = s['cycle_time_years']
            print(f"{level:8s} {shield:13s}: transmission {s['transmission'][0]:.1e}-{s['transmission'][1]:.1e}; "
                  f"loss {lo:.3g}-{hi:.3g} kg/s (cycle {c['high'][0]:.2g}-{c['low'][1]:.2g} yr); glow share "
                  f"{s['glow_share_of_heating'][0]:.2f}-{s['glow_share_of_heating'][1]:.2f}; exobase "
                  f"{s['exobase_radius_R'][0]:.2f}-{s['exobase_radius_R'][1]:.2f} R; radius "
                  f"{s['protected_radius_R'][0]:.2f}-{s['protected_radius_R'][1]:.2f} R (area x"
                  f"{s['aperture_area_over_september'][0]:.2f}-{s['aperture_area_over_september'][1]:.2f}); ions "
                  f"{s['exosphere_ion_production_kg_s'][0]:.3g}-{s['exosphere_ion_production_kg_s'][1]:.3g} kg/s")
    digest = lambda p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()[:16]
    out = dict(
        schema=SCHEMA,
        producer=dict(study='protection_architecture', files={p: digest(p) for p in (
            'research/studies/protection_architecture/run.py', 'atmosphere/loss_response/model.py',
            'atmosphere/loss_response/absorption.py', 'atmosphere/thermal_column.py')},
            inputs={p: digest(p) for p in ('protection/transmission/results/uv_transmission.json',
                                           'research/baselines/feasibility/reference.json')}),
        evidence=('The loss response and its absorption step evaluated at the ultraviolet transmission the swarm model '
                  'gives for three levels of design choices. Screening models throughout; the loss leaves out Earth\'s '
                  'tidal barrier and the escape of exospheric ions, which can each exceed it.'),
        reading_rule=('summary[level][shield] gives ranges over the upper-air treatments, quiet Sun and solar maximum, '
                      'and holes 1 and 4 times the particle diameter. total_loss_kg_s adds the solar-wind screening '
                      'range to the ultraviolet-driven loss; cycle_time_years is the atmosphere\'s mass over that loss. '
                      'protected_radius_R is in lunar radii from the Moon\'s centre.'),
        atmosphere_mass_kg=mass, summary=summary, points=points)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / 'design_point.json').write_text(json.dumps(out, indent=2) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
