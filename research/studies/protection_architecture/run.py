"""The protection's design point: the ultraviolet transmission a swarm can hold, and what it gives the air.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.protection_architecture.run   # results in studies/protection_architecture/results

Reads the swarm's ultraviolet transmission (protection/transmission, product
terluna.protection.uv-transmission/1) for three levels of design choices and
holes 1 and 4 times the particle diameter. At each transmission it solves the
loss response's thermal column (atmosphere/loss_response) behind each shield,
for each treatment of the upper air, at quiet Sun and at solar maximum, and
reports:

- the ultraviolet-driven loss, with Earth's tide;
- the protected radius at which ultraviolet deposited in the thermosphere
  outside the shadow is a tenth of all the heat the upper air receives
  (transmitted sunlight, the film's X-rays and the sky's Lyman-alpha glow,
  expressed as an equivalent transmission), and its aperture against the
  September design's;
- what the sunlit exosphere outside that radius loses (the exosphere step:
  ions, fragments of broken molecules and the solar wind's charge exchange),
  with no magnetosphere and with the September design's magnets, low, central
  and high;
- the total loss in both cases: the ultraviolet-driven loss plus the
  exosphere's, and without a magnetosphere also the sputtering by
  solar-wind protons; and the atmospheric cycle time, the atmosphere's mass
  over that loss, the time the loss and the resupply that matches it take to
  replace the whole atmosphere;
- the central total with protected radii of 3, 4 and 6 lunar radii, where
  they cover the heating; the protected radius at which the exosphere's central
  loss falls to 0.1, 1 and 10 kg/s; and, with a protected radius of 4 lunar
  radii (or the heating radius where larger), the exosphere's loss against the
  moment of a lunar dipole and the smallest moment that holds it to 0.1, 1
  and 10 kg/s.

The atmosphere's mass is the feasibility baseline's hydrostatic column
(research/baselines/feasibility/reference.json). The loss leaves out atomic
oxygen, photochemistry below the exobase and hydrogen.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

from atmosphere.loss_response import absorption as ab
from atmosphere.loss_response import exosphere as ex
from atmosphere.loss_response import model as lr
from atmosphere.middle_atmosphere import escape

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.protection-design-point/2'
TRANSMISSION = ROOT / 'protection' / 'transmission' / 'results' / 'uv_transmission.json'
TRANSMISSION_SCHEMA = 'terluna.protection.uv-transmission/1'
BASELINE = ROOT / 'research' / 'baselines' / 'feasibility' / 'reference.json'
YEAR_S = 365.25 * 86400.0
HOLE_FACTORS = ('1', '4')
WIND_LEVELS = ('low', 'central', 'high')
SCENARIOS = ('no_magnetosphere', 'september_magnetosphere')
PROTECTED_RADII = (3.0, 4.0, 6.0)
EXOSPHERE_TARGETS_KG_S = (0.1, 1.0, 10.0)
SWEEP_RADIUS_R = 4.0


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
    _, profile, _ = ab.column(shield, treatment, activity, transmission)
    exo = ex.assess(profile, activity, heat['radius_R'], sweep_shadow_R=max(SWEEP_RADIUS_R, heat['radius_R']))
    protons = {level: wind[level]['proton_sputtering_kg_s'] for level in WIND_LEVELS}
    totals = dict(no_magnetosphere={lv: uv + protons[lv] + exo['no_magnetosphere'][lv]['total_kg_s'] for lv in WIND_LEVELS},
                  september_magnetosphere={lv: uv + exo['september_magnetosphere'][lv]['total_kg_s']
                                           for lv in WIND_LEVELS})
    curve = exo['against_radius']
    radii = curve['radius_R']
    at_radius = {}
    for x in PROTECTED_RADII:
        covers = x >= heat['radius_R']
        at_radius[f'{x:g}'] = dict(
            no_magnetosphere=uv + protons['central'] + float(np.interp(x, radii, curve['no_magnetosphere_kg_s']))
            if covers else None,
            september_magnetosphere=uv + float(np.interp(x, radii, curve['september_magnetosphere_kg_s']))
            if covers else None)
    needed = {f'{target:g}': {s: ab.smallest_radius(np.array(radii), np.array(curve[f'{s}_kg_s']), target)
                              for s in SCENARIOS} for target in EXOSPHERE_TARGETS_KG_S}
    return dict(shield=shield, treatment=treatment, activity=activity, transmission=transmission,
                status=row['status'], heating_as_equivalent_transmission=equivalent,
                glow_share_of_heating=row['glow_heat_w_m2'] / row['deposited_heat_w_m2'],
                uv_driven_loss_kg_s=uv, proton_sputtering_kg_s=protons,
                exobase_radius_R=case['exobase_radius_R'], exobase_temperature_k=case['exobase_temperature_k'],
                protected_radius_R=heat['radius_R'], aperture=heat['aperture'],
                exosphere=dict(ions_made_kg_s=exo['sunlit']['ions_made_kg_s']['central'],
                               fragments_escaping_kg_s=exo['sunlit']['fragments_escaping_kg_s']['central'],
                               charge_exchange_upper_bound_kg_s=exo['sunlit']['charge_exchange_upper_bound_kg_s']['mixed'],
                               loss_kg_s={s: {lv: exo[s][lv]['total_kg_s'] for lv in WIND_LEVELS} for s in SCENARIOS}),
                total_loss_kg_s=totals,
                cycle_time_years={s: {k: cycle_time_years(v, mass_kg) for k, v in totals[s].items()} for s in SCENARIOS},
                total_at_protected_radius_kg_s=at_radius,
                protected_radius_for_exosphere_loss_R=needed,
                moment_sweep=exo['moment_sweep'],
                moment_for_exosphere_loss_A_m2={f'{target:g}': smallest_moment(exo['moment_sweep'], target)
                                                 for target in EXOSPHERE_TARGETS_KG_S})


def smallest_moment(sweep, target):
    """The smallest swept moment at which the exosphere's central loss is within target; None if none."""
    for moment, loss in zip(sweep['moment_A_m2'], sweep['loss_kg_s']):
        if loss <= target:
            return moment
    return None


def span(rows, get):
    values = [v for v in (get(r) for r in rows) if v is not None]
    return [min(values), max(values)] if values else None


def missing(rows, get):
    return sum(get(r) is None for r in rows)


def summarise(points):
    out = {}
    for level in sorted({p['level'] for p in points}):
        out[level] = {}
        for shield in lr.SHIELDS:
            rows = [p for p in points if p['level'] == level and p['shield'] == shield]
            moments = rows[0]['moment_sweep']['moment_A_m2']
            out[level][shield] = dict(
                cases=len(rows),
                transmission=span(rows, lambda r: r['transmission']),
                uv_driven_loss_kg_s=span(rows, lambda r: r['uv_driven_loss_kg_s']),
                glow_share_of_heating=span(rows, lambda r: r['glow_share_of_heating']),
                exobase_radius_R=span(rows, lambda r: r['exobase_radius_R']),
                protected_radius_R=span(rows, lambda r: r['protected_radius_R']),
                aperture_area_over_september=span(rows, lambda r: r['aperture']['area_over_september']),
                exosphere_ions_made_kg_s=span(rows, lambda r: r['exosphere']['ions_made_kg_s']),
                exosphere_fragments_kg_s=span(rows, lambda r: r['exosphere']['fragments_escaping_kg_s']),
                exosphere_charge_exchange_upper_bound_kg_s=span(
                    rows, lambda r: r['exosphere']['charge_exchange_upper_bound_kg_s']),
                total_loss_kg_s={s: {w: span(rows, lambda r, s=s, w=w: r['total_loss_kg_s'][s][w])
                                     for w in WIND_LEVELS} for s in SCENARIOS},
                cycle_time_years={s: {w: span(rows, lambda r, s=s, w=w: r['cycle_time_years'][s][w])
                                      for w in WIND_LEVELS} for s in SCENARIOS},
                total_at_protected_radius_kg_s={
                    f'{x:g}': {s: dict(range=span(rows, lambda r, x=x, s=s: r['total_at_protected_radius_kg_s'][f'{x:g}'][s]),
                                       cases_not_covering_the_heating=missing(
                                           rows, lambda r, x=x, s=s: r['total_at_protected_radius_kg_s'][f'{x:g}'][s]))
                               for s in SCENARIOS} for x in PROTECTED_RADII},
                protected_radius_for_exosphere_loss_R={
                    f'{t:g}': {s: dict(range=span(rows, lambda r, t=t, s=s: r['protected_radius_for_exosphere_loss_R'][f'{t:g}'][s]),
                                       cases_beyond_11_5_R=missing(
                                           rows, lambda r, t=t, s=s: r['protected_radius_for_exosphere_loss_R'][f'{t:g}'][s]))
                               for s in SCENARIOS} for t in EXOSPHERE_TARGETS_KG_S},
                loss_against_moment_kg_s=dict(
                    shadow_radius_R=span(rows, lambda r: r['moment_sweep']['shadow_radius_R']), moment_A_m2=moments,
                    loss=[span(rows, lambda r, i=i: r['moment_sweep']['loss_kg_s'][i]) for i in range(len(moments))]),
                moment_for_exosphere_loss_A_m2={
                    f'{t:g}': dict(range=span(rows, lambda r, t=t: r['moment_for_exosphere_loss_A_m2'][f'{t:g}']),
                                   cases_beyond_1e23=missing(rows, lambda r, t=t: r['moment_for_exosphere_loss_A_m2'][f'{t:g}']))
                    for t in EXOSPHERE_TARGETS_KG_S})
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
            none, sept = s['total_loss_kg_s']['no_magnetosphere'], s['total_loss_kg_s']['september_magnetosphere']
            at = s['total_at_protected_radius_kg_s']
            fmt = lambda r: ('%.3g-%.3g' % tuple(r)) if r else '-'
            print(f"{level:8s} {shield:13s}: transmission {fmt(s['transmission'])}; UV-driven {fmt(s['uv_driven_loss_kg_s'])}; "
                  f"exobase {fmt(s['exobase_radius_R'])} R, radius {fmt(s['protected_radius_R'])} R; "
                  f"total none {none['low'][0]:.3g}-{none['high'][1]:.3g} (central {fmt(none['central'])}), "
                  f"September {sept['low'][0]:.3g}-{sept['high'][1]:.3g} (central {fmt(sept['central'])}); "
                  f"central at 3/4/6 R: none {fmt(at['3']['no_magnetosphere']['range'])} / {fmt(at['4']['no_magnetosphere']['range'])} / "
                  f"{fmt(at['6']['no_magnetosphere']['range'])}, September {fmt(at['3']['september_magnetosphere']['range'])} / "
                  f"{fmt(at['4']['september_magnetosphere']['range'])} / {fmt(at['6']['september_magnetosphere']['range'])}")
    digest = lambda p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()[:16]
    out = dict(
        schema=SCHEMA,
        producer=dict(study='protection_architecture', files={p: digest(p) for p in (
            'research/studies/protection_architecture/run.py', 'atmosphere/loss_response/model.py',
            'atmosphere/loss_response/absorption.py', 'atmosphere/loss_response/exosphere.py',
            'atmosphere/thermal_column.py')},
            inputs={p: digest(p) for p in ('protection/transmission/results/uv_transmission.json',
                                           'research/baselines/feasibility/reference.json',
                                           'atmosphere/loss_response/results/tidal_escape.json')}),
        evidence=('The loss response, its absorption step and its exosphere step evaluated at the ultraviolet '
                  'transmission the swarm model gives for three levels of design choices. Screening models '
                  'throughout: Earth\'s tide is included; the exosphere\'s ion fates are bounds and timescale '
                  'comparisons, not a plasma model.'),
        reading_rule=('summary[level][shield] gives ranges over the upper-air treatments, quiet Sun and solar maximum, '
                      'and holes 1 and 4 times the particle diameter. total_loss_kg_s[scenario][low|central|high] is '
                      'the ultraviolet-driven loss plus the sunlit exosphere\'s loss outside the protected radius, with '
                      'no magnetosphere (plus proton sputtering) or with the September magnets; cycle_time_years is '
                      'the atmosphere\'s mass over that loss. total_at_protected_radius_kg_s gives the central total '
                      'with a wider protected radius, null where that radius does not cover the heating. '
                      'protected_radius_R is in lunar radii from the Moon\'s centre.'),
        atmosphere_mass_kg=mass, summary=summary, points=points)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / 'design_point.json').write_text(json.dumps(out, indent=2) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
