"""The protection's design point: the ultraviolet transmission a swarm can hold, and what it gives the air.

    OPENBLAS_NUM_THREADS=1 python -m research.studies.protection_architecture.run   # results in studies/protection_architecture/results

Reads the swarm's ultraviolet transmission (protection/transmission, product
terluna.protection.uv-transmission/1) for three levels of design choices and
holes 1 and 4 times the particle diameter. At each transmission it takes the
loss response's traced state (atmosphere/loss_response, from the middle
atmosphere's limb tracing: the transmission passes gaps over the aperture of the
protected radius, the films' X-rays and unfiltered light beyond the aperture
follow their slant paths, with the sky's Lyman-alpha glow) behind each shield,
for each treatment of the upper air, at quiet Sun and at solar maximum, for the
ring fleet's protected radius of 4 lunar radii, and reports:

- the ultraviolet-driven loss, with Earth's tide, and the heat behind it;
- the smallest protected radius that covers the heating, where light beyond
  the aperture gives at most a tenth of the heat;
- what the sunlit exosphere outside the shadow loses (the exosphere step:
  ions, fragments of broken molecules and the solar wind's charge exchange),
  with no magnetosphere, with no magnetosphere but the wind held off within the
  ring fleet's unrefilled wake (a sensitivity), and with the September design's
  magnets, low, central and high;
- the total loss in each case: the ultraviolet-driven loss plus the
  exosphere's, and without a magnetosphere also the sputtering by
  solar-wind protons; and the atmospheric cycle time, the atmosphere's mass
  over that loss, the time the loss and the resupply that matches it take to
  replace the whole atmosphere;
- the central total with protected radii of 3, 4 and 6 lunar radii, each with
  its own traced state, where they cover the heating; the protected radius at
  which the exosphere's central loss falls to 0.1, 1 and 10 kg/s; and the
  exosphere's loss against the moment of a lunar dipole and the smallest moment
  that holds it to 0.1, 1 and 10 kg/s;
- the same central totals over the solar cycle (atmosphere/loss_response/cycle.py): each calendar year of FISM2's
  daily record taken as a steady state, at protected radii of 4 and 6 lunar radii, the mean over solar cycles 23 and
  24 (1997-2019) and over cycle 25 to 2025, the year of the largest loss, and the state of each span's mean spectrum.
  R1 is a long-term average, so the cycle mean is the figure it reads.

The atmosphere's mass is the feasibility baseline's hydrostatic column
(research/baselines/feasibility/reference.json). The loss leaves out atomic
oxygen, photochemistry below the exobase and hydrogen.
"""
from __future__ import annotations
import argparse
import hashlib
import json
from multiprocessing import Pool
import os
from pathlib import Path
import sys

import numpy as np

from atmosphere.loss_response import absorption as ab
from atmosphere.loss_response import cycle
from atmosphere.loss_response import exosphere as ex
from atmosphere.loss_response import model as lr
from atmosphere.loss_response import traced

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
SCHEMA = 'terluna.research.protection-design-point/3'
TRANSMISSION = ROOT / 'protection' / 'transmission' / 'results' / 'uv_transmission.json'
TRANSMISSION_SCHEMA = 'terluna.protection.uv-transmission/1'
BASELINE = ROOT / 'research' / 'baselines' / 'feasibility' / 'reference.json'
YEAR_S = 365.25 * 86400.0
HOLE_FACTORS = ('1', '4')
WIND_LEVELS = ('low', 'central', 'high')
SCENARIOS = ('no_magnetosphere', 'no_magnetosphere_wake', 'september_magnetosphere')
CURVE_SCENARIOS = ('no_magnetosphere', 'september_magnetosphere')
PROTECTED_RADII = (3.0, 4.0, 6.0)
DESIGN_RADIUS_R = lr.PROTECTED_R            # the ring fleet's (decisions.md, 2026-10-07)
EXOSPHERE_TARGETS_KG_S = (0.1, 1.0, 10.0)
CYCLE_RADII = (4.0, 6.0)
CYCLE_WORKERS = 3                            # single-threaded processes; the machine is shared


def atmosphere_mass_kg():
    return float(json.loads(BASELINE.read_text())['baseline']['mass_kg'])


def cycle_time_years(loss_kg_s, mass_kg):
    return mass_kg / loss_kg_s / YEAR_S


def at_radius(shield, treatment, activity, transmission, radius, wind, sweep=False):
    """The traced state at one protected radius, the sunlit exosphere's losses outside it and the totals; None where
    the air runs away or leaves the thermal column's domain."""
    row = lr.euv_sweep(shield, treatment, activity, glow=True, leaks=(transmission,), protected_R=radius)[0]
    if row['status'] != 'thermal_column':
        return None, row
    uv = lr.loss_estimate(row)
    _, profile, _ = ab.column(shield, treatment, activity, transmission, protected_R=radius)
    exo = ex.assess(profile, activity, radius, sweep=sweep, wake_R=ex.wake_core_R(radius))
    protons = {level: wind[level]['proton_sputtering_kg_s'] for level in WIND_LEVELS}
    totals = dict(no_magnetosphere={lv: uv + protons[lv] + exo['no_magnetosphere'][lv]['total_kg_s'] for lv in WIND_LEVELS},
                  no_magnetosphere_wake={lv: uv + protons[lv] + exo['no_magnetosphere_wake'][lv]['total_kg_s']
                                         for lv in WIND_LEVELS},
                  september_magnetosphere={lv: uv + exo['september_magnetosphere'][lv]['total_kg_s']
                                           for lv in WIND_LEVELS})
    return dict(uv=uv, exo=exo, totals=totals, covers=row['beyond_share'] <= traced.HEATING_SHARE), row



def point(shield, treatment, activity, transmission, wind, mass_kg):
    states = {x: at_radius(shield, treatment, activity, transmission, x, wind, sweep=(x == DESIGN_RADIUS_R))
              for x in PROTECTED_RADII}
    design, row = states[DESIGN_RADIUS_R]
    out = dict(shield=shield, treatment=treatment, activity=lr.activity_label(activity), transmission=transmission,
               status=row['status'],
               protected_radius_R=DESIGN_RADIUS_R,
               heating_radius_R=ab.heating_radius(shield, treatment, activity, transmission),
               wake_core_R=ex.wake_core_R(DESIGN_RADIUS_R),
               total_at_protected_radius_kg_s={
                   f'{x:g}': ({s: states[x][0]['totals'][s]['central'] for s in SCENARIOS}
                              if states[x][0] is not None and states[x][0]['covers'] else None)
                   for x in PROTECTED_RADII})
    if design is None:
        return out
    uv, exo, totals = design['uv'], design['exo'], design['totals']
    curve = exo['against_radius']
    radii = curve['radius_R']
    needed = {f'{target:g}': {s: ab.smallest_radius(np.array(radii), np.array(curve[f'{s}_kg_s']), target)
                              for s in CURVE_SCENARIOS} for target in EXOSPHERE_TARGETS_KG_S}
    heat = row['deposited_heat_w_m2']
    return dict(out, heat_w_m2=heat,
                heat_parts_w_m2=dict(gaps=row['leak_heat_w_m2'], films=row['film_heat_w_m2'], glow=row['glow_heat_w_m2'],
                                     beyond_aperture=row['beyond_aperture_heat_w_m2']),
                heating_as_disk_count_transmission=heat / (row['disk_count_heat_w_m2'] / transmission),
                glow_share_of_heating=row['glow_heat_w_m2'] / heat, beyond_share=row['beyond_share'],
                covers_heating=design['covers'], uv_driven_loss_kg_s=uv,
                proton_sputtering_kg_s={level: wind[level]['proton_sputtering_kg_s'] for level in WIND_LEVELS},
                exobase_radius_R=row['exobase_radius_R'], exobase_temperature_k=row['exobase_temperature_k'],
                exosphere=dict(ions_made_kg_s=exo['sunlit']['ions_made_kg_s']['central'],
                               fragments_escaping_kg_s=exo['sunlit']['fragments_escaping_kg_s']['central'],
                               charge_exchange_upper_bound_kg_s=exo['sunlit']['charge_exchange_upper_bound_kg_s']['mixed'],
                               charge_exchange_in_wake_upper_bound_kg_s=exo['sunlit'][
                                   'charge_exchange_in_wake_upper_bound_kg_s']['mixed'],
                               loss_kg_s={s: {lv: exo[s][lv]['total_kg_s'] for lv in WIND_LEVELS} for s in SCENARIOS}),
                total_loss_kg_s=totals,
                cycle_time_years={s: {k: cycle_time_years(v, mass_kg) for k, v in totals[s].items()} for s in SCENARIOS},
                protected_radius_for_exosphere_loss_R=needed,
                moment_sweep=exo['moment_sweep'],
                moment_for_exosphere_loss_A_m2={f'{target:g}': smallest_moment(exo['moment_sweep'], target)
                                                 for target in EXOSPHERE_TARGETS_KG_S})


def cycle_state(shield, treatment, activity, transmission, radius, wind):
    """One measured year's traced state at a protected radius and its central totals; None where the air runs away,
    leaves the thermal column or the radius does not cover the heating."""
    row = lr.euv_sweep(shield, treatment, activity, glow=True, leaks=(transmission,), protected_R=radius)[0]
    if row['status'] != 'thermal_column' or row['beyond_share'] > traced.HEATING_SHARE:
        return None
    _, profile, _ = ab.column(shield, treatment, activity, transmission, protected_R=radius)
    exo = ex.central_losses(profile, activity, [radius])
    uv, protons = row['molecular_loss_kg_s'], wind['central']['proton_sputtering_kg_s']
    return dict(heat_w_m2=row['deposited_heat_w_m2'], exobase_radius_R=row['exobase_radius_R'],
                exobase_temperature_k=row['exobase_temperature_k'], uv_driven_loss_kg_s=uv,
                total_kg_s=dict(no_magnetosphere=uv + protons + exo['none'][0],
                                no_magnetosphere_wake=uv + protons + exo['none_wake'][0],
                                september_magnetosphere=uv + exo['september'][0]))


def cycle_point(task):
    """The solar cycle of one case at one transmission: each measured year's state and central totals at each
    protected radius, the mean over each span's years (None if any year has no state), the year of the largest
    total, and the state of the span's mean spectrum."""
    shield, treatment, transmission, wind = task
    out = {}
    for name, span_ in cycle.product()['spans'].items():
        years = {a['label']: {f'{x:g}': cycle_state(shield, treatment, a, transmission, x, wind) for x in CYCLE_RADII}
                 for a in span_['years']}
        flat = {f'{x:g}': cycle_state(shield, treatment, span_['mean_spectrum'], transmission, x, wind)
                for x in CYCLE_RADII}
        summary = {}
        for x in CYCLE_RADII:
            key = f'{x:g}'
            states = {y: v[key] for y, v in years.items()}
            held = {y: v for y, v in states.items() if v is not None}
            whole = len(held) == len(states)
            summary[key] = dict(
                years_without_a_state=sorted(set(states) - set(held)),
                mean_total_kg_s={sc: float(np.mean([v['total_kg_s'][sc] for v in held.values()])) if whole else None
                                 for sc in SCENARIOS},
                mean_uv_driven_loss_kg_s=float(np.mean([v['uv_driven_loss_kg_s'] for v in held.values()]))
                if whole else None,
                largest_total_kg_s={sc: (max((v['total_kg_s'][sc], y) for y, v in held.items()) if held else None)
                                    for sc in SCENARIOS},
                mean_spectrum_total_kg_s=flat[key]['total_kg_s'] if flat[key] else None)
        out[name] = dict(summary=summary, years=years, mean_spectrum=flat)
    return out


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
            held = [r for r in rows if 'uv_driven_loss_kg_s' in r]
            moments = next(r['moment_sweep']['moment_A_m2'] for r in held) if held else []
            out[level][shield] = dict(
                cases=len(rows), cases_running_away_at_4_R=len(rows) - len(held),
                transmission=span(rows, lambda r: r['transmission']),
                heat_w_m2=span(held, lambda r: r['heat_w_m2']),
                uv_driven_loss_kg_s=span(held, lambda r: r['uv_driven_loss_kg_s']),
                glow_share_of_heating=span(held, lambda r: r['glow_share_of_heating']),
                beyond_share=span(held, lambda r: r['beyond_share']),
                exobase_radius_R=span(held, lambda r: r['exobase_radius_R']),
                exobase_temperature_k=span(held, lambda r: r['exobase_temperature_k']),
                heating_radius_R=span(rows, lambda r: r['heating_radius_R']),
                cases_without_heating_radius=missing(rows, lambda r: r['heating_radius_R']),
                exosphere_ions_made_kg_s=span(held, lambda r: r['exosphere']['ions_made_kg_s']),
                exosphere_fragments_kg_s=span(held, lambda r: r['exosphere']['fragments_escaping_kg_s']),
                exosphere_charge_exchange_upper_bound_kg_s=span(
                    held, lambda r: r['exosphere']['charge_exchange_upper_bound_kg_s']),
                exosphere_charge_exchange_in_wake_upper_bound_kg_s=span(
                    held, lambda r: r['exosphere']['charge_exchange_in_wake_upper_bound_kg_s']),
                total_loss_kg_s={s: {w: span(held, lambda r, s=s, w=w: r['total_loss_kg_s'][s][w])
                                     for w in WIND_LEVELS} for s in SCENARIOS},
                cycle_time_years={s: {w: span(held, lambda r, s=s, w=w: r['cycle_time_years'][s][w])
                                      for w in WIND_LEVELS} for s in SCENARIOS},
                total_at_protected_radius_kg_s={
                    f'{x:g}': {s: dict(range=span(rows, lambda r, x=x, s=s: (r['total_at_protected_radius_kg_s'][f'{x:g}'] or {}).get(s)),
                                       cases_not_covering_the_heating=missing(
                                           rows, lambda r, x=x, s=s: (r['total_at_protected_radius_kg_s'][f'{x:g}'] or {}).get(s)))
                               for s in SCENARIOS} for x in PROTECTED_RADII},
                protected_radius_for_exosphere_loss_R={
                    f'{t:g}': {s: dict(range=span(held, lambda r, t=t, s=s: r['protected_radius_for_exosphere_loss_R'][f'{t:g}'][s]),
                                       cases_beyond_11_5_R=missing(held, lambda r, t=t, s=s: r['protected_radius_for_exosphere_loss_R'][f'{t:g}'][s]))
                               for s in CURVE_SCENARIOS} for t in EXOSPHERE_TARGETS_KG_S},
                loss_against_moment_kg_s=dict(
                    shadow_radius_R=span(held, lambda r: r['moment_sweep']['shadow_radius_R']), moment_A_m2=moments,
                    loss=[span(held, lambda r, i=i: r['moment_sweep']['loss_kg_s'][i]) for i in range(len(moments))]),
                moment_for_exosphere_loss_A_m2={
                    f'{t:g}': dict(range=span(held, lambda r, t=t: r['moment_for_exosphere_loss_A_m2'][f'{t:g}']),
                                   cases_beyond_1e23=missing(held, lambda r, t=t: r['moment_for_exosphere_loss_A_m2'][f'{t:g}']))
                    for t in EXOSPHERE_TARGETS_KG_S})
    return out


def summarise_cycle(cycles):
    """Ranges over the upper-air treatments and hole factors of the cycle means, the largest yearly totals and the mean
    spectrum's totals, by span, level, shield and protected radius."""
    out = {}
    for name in cycle.SPANS:
        for level in sorted({c['level'] for c in cycles}):
            for shield in lr.SHIELDS:
                rows = [c[name]['summary'] for c in cycles if c['level'] == level and c['shield'] == shield]
                out.setdefault(name, {}).setdefault(level, {})[shield] = {
                    f'{x:g}': dict(
                        cases=len(rows),
                        cases_with_a_year_without_a_state=sum(bool(r[f'{x:g}']['years_without_a_state']) for r in rows),
                        mean_total_kg_s={sc: span(rows, lambda r, x=x, sc=sc: r[f'{x:g}']['mean_total_kg_s'][sc])
                                         for sc in SCENARIOS},
                        mean_uv_driven_loss_kg_s=span(rows, lambda r, x=x: r[f'{x:g}']['mean_uv_driven_loss_kg_s']),
                        largest_yearly_total_kg_s={sc: span(rows, lambda r, x=x, sc=sc: (
                            r[f'{x:g}']['largest_total_kg_s'][sc] or (None,))[0]) for sc in SCENARIOS},
                        mean_spectrum_total_kg_s={sc: span(rows, lambda r, x=x, sc=sc: (
                            r[f'{x:g}']['mean_spectrum_total_kg_s'] or {}).get(sc)) for sc in SCENARIOS})
                    for x in CYCLE_RADII}
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
    points = []
    for level, lv in product['levels'].items():
        for k in HOLE_FACTORS:
            t = lv['transmission_by_hole_factor'][k]
            for shield in lr.SHIELDS:
                for treatment in lr.TREATMENTS:
                    for activity in lr.ACTIVITY:
                        points.append(dict(point(shield, treatment, activity, t, wind, mass),
                                           level=level, hole_factor=float(k)))
    summary = summarise([p for p in points if p['activity'] not in lr.STRESS_CASES])
    stress_summary = {name: summarise([p for p in points if p['activity'] == name]) for name in lr.STRESS_CASES}
    tasks, keys = [], []
    for level, lv in product['levels'].items():
        for k in HOLE_FACTORS:
            for shield in lr.SHIELDS:
                for treatment in lr.TREATMENTS:
                    tasks.append((shield, treatment, lv['transmission_by_hole_factor'][k], wind))
                    keys.append(dict(shield=shield, treatment=treatment, level=level, hole_factor=float(k),
                                     transmission=lv['transmission_by_hole_factor'][k]))
    with Pool(CYCLE_WORKERS, initializer=os.nice, initargs=(10,)) as pool:
        cycles = [dict(key, **result) for key, result in zip(keys, pool.map(cycle_point, tasks, chunksize=1))]
    cycle_summary = summarise_cycle(cycles)
    for level, by_shield in summary.items():
        for shield, s in by_shield.items():
            none, sept = s['total_loss_kg_s']['no_magnetosphere'], s['total_loss_kg_s']['september_magnetosphere']
            wake = s['total_loss_kg_s']['no_magnetosphere_wake']
            at = s['total_at_protected_radius_kg_s']
            fmt = lambda r: ('%.3g-%.3g' % tuple(r)) if r else '-'
            print(f"{level:8s} {shield:13s}: transmission {fmt(s['transmission'])}; UV-driven {fmt(s['uv_driven_loss_kg_s'])}; "
                  f"exobase {fmt(s['exobase_radius_R'])} R, heating radius {fmt(s['heating_radius_R'])} R; "
                  f"runaway at 4 R {s['cases_running_away_at_4_R']}; central totals none {fmt(none['central'])}, "
                  f"wake {fmt(wake['central'])}, September {fmt(sept['central'])}; "
                  f"central at 3/4/6 R: none {fmt(at['3']['no_magnetosphere']['range'])} / {fmt(at['4']['no_magnetosphere']['range'])} / "
                  f"{fmt(at['6']['no_magnetosphere']['range'])}, September {fmt(at['3']['september_magnetosphere']['range'])} / "
                  f"{fmt(at['4']['september_magnetosphere']['range'])} / {fmt(at['6']['september_magnetosphere']['range'])}")
    digest = lambda p: hashlib.sha256((ROOT / p).read_bytes()).hexdigest()[:16]
    out = dict(
        schema=SCHEMA,
        producer=dict(study='protection_architecture', files={p: digest(p) for p in (
            'research/studies/protection_architecture/run.py', 'atmosphere/loss_response/model.py',
            'atmosphere/loss_response/absorption.py', 'atmosphere/loss_response/exosphere.py',
            'atmosphere/loss_response/traced.py', 'atmosphere/loss_response/cycle.py', 'atmosphere/thermal_column.py')},
            inputs={p: digest(p) for p in ('protection/transmission/results/uv_transmission.json',
                                           'research/baselines/feasibility/reference.json',
                                           'atmosphere/loss_response/results/tidal_escape.json',
                                           'atmosphere/middle_atmosphere/results/limb_heat.json',
                                           'atmosphere/loss_response/results/solar_cycle.json')}),
        evidence=('The loss response, its absorption step and its exosphere step evaluated at the ultraviolet '
                  'transmission the swarm model gives for three levels of design choices, with the heat the middle '
                  'atmosphere\'s limb tracing finds along slant paths for each protected radius, placed in height '
                  'where the tracing puts it. Screening models '
                  'throughout: Earth\'s tide is included; the exosphere\'s ion fates are bounds and timescale '
                  'comparisons, and the ring fleet\'s wake rests on an assumed refill length; no plasma is modelled.'),
        reading_rule=('summary[level][shield] gives ranges over the upper-air treatments, quiet Sun and solar maximum '
                      '(FISM2\'s year around cycle 21\'s maximum), and holes 1 and 4 times the particle diameter, at the '
                      'ring fleet\'s protected radius of 4 lunar radii; stress_summary[case] gives the same at each stress '
                      'case alone (cycle 19\'s year, the strongest on record, and the escape model\'s 2.5 on the '
                      'ultraviolet), and points every state. '
                      'cases_running_away_at_4_R counts those whose air outgrows the limb tables there. '
                      'total_loss_kg_s[scenario][low|central|high] is the ultraviolet-driven loss plus the sunlit '
                      'exosphere\'s loss outside the protected radius, with no magnetosphere (plus proton sputtering), '
                      'with no magnetosphere but the wind held off within the ring fleet\'s unrefilled wake (a '
                      'sensitivity), or with the September magnets; cycle_time_years is the atmosphere\'s mass over that '
                      'loss. total_at_protected_radius_kg_s gives the central total at each protected radius with its '
                      'own traced state, null where the air runs away or that radius does not cover the heating. '
                      'heating_radius_R is the smallest protected radius covering the heating; radii are in lunar radii '
                      'from the Moon\'s centre. cycle_summary[span][level][shield][radius] gives the same central totals '
                      'over the solar cycle: the mean of the measured years\' steady states (null if any year has none), '
                      'which leans high, the largest yearly total, and the state of the span\'s mean spectrum, which gives '
                      'the low end; R1, a long-term average, reads the cycle mean.'),
        atmosphere_mass_kg=mass, summary=summary, stress_summary=stress_summary, points=points,
        cycle_summary=cycle_summary, cycle=cycles)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / 'design_point.json').write_text(json.dumps(out, indent=2) + '\n')
    return 0


if __name__ == '__main__':
    sys.exit(main())
