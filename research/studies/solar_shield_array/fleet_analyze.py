"""Audit finite-tile coverage and exclusion constraints on actual fleet histories."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import time

import numpy as np

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.fleet import (fleet_acceleration, coverage_snapshot,
    swept_conflicts, conflict_free_subset)
from protection.dynamics.optical import unit
from .fleet_run import read_raw, environment, identities, digest, RUN
from .cycling_search import ROOT, HERE
from .cycling_analysis import compact_series_json, time_average


def summarize_coverage(records, times):
    result = {}
    for field in ['ray_coverage', 'core_ray_coverage', 'all_sampled_sun_coverage',
                  'natural_eclipse_fraction', 'excess_intersection_area_fraction']:
        values = np.array([x[field] if x[field] is not None else np.nan for x in records])
        if not np.all(np.isfinite(values)):
            raise ValueError('Fully eclipsed samples need a separate incident-flux quadrature')
        result[field] = dict(mean=time_average(values, times), minimum=float(values.min()),
                             maximum=float(values.max()))
    result['aperture_complete_at_any_sample'] = any(r['all_sampled_sun_coverage'] >= 1-1e-12 for r in records)
    result['aperture_complete_at_every_sample'] = all(r['all_sampled_sun_coverage'] >= 1-1e-12 for r in records)
    result['days_with_incomplete_aperture_on_sample_grid'] = time_average(
        [r['all_sampled_sun_coverage'] < 1-1e-12 for r in records], times)*np.ptp(times)/K.JULIAN_DAY
    return result


def analyze(name, cadence, suns, rotation, clearances=True):
    archive, meta = read_raw(name)
    t, state = archive['t'], archive['state']; config = meta['config']; count = meta['count']
    env = environment(config['days']); side = config['side_km']*1000
    selected_times = np.linspace(0, t[-1], int(np.ceil(t[-1]/cadence))+1)
    index = np.searchsorted(t, selected_times)
    if np.max(np.abs(t[index]-selected_times)) > 1e-7:
        raise ValueError('Coverage cadence must align with stored integration output')
    sun = unit(env.at(t)['positions']['sun'])
    before = time.monotonic()
    print(json.dumps(dict(stage='swept_exclusions', name=name, members=count)), flush=True)
    collision, shadow = swept_conflicts(t, state, side, sun_directions=sun)
    subsets = dict(nominal=np.arange(count),
        collision_screened=conflict_free_subset(count, collision),
        collision_and_shadow_screened=conflict_free_subset(count, set(collision) | set(shadow)))
    # Every pair retained passes the conservative guards over every stored
    # interval. The full population remains a rejected geometric probe when
    # collisions or unmodelled mutual-force shadows occur.
    selection = {}
    for label, members in subsets.items():
        chosen = set(members.tolist())
        selection[label] = dict(count=len(members), physical_area_m2=len(members)*side**2,
            optical_mass_kg=len(members)*side**2*.05,
            collision_pairs_retained=sum(i in chosen and j in chosen for i, j in collision),
            possible_shadow_pairs_retained=sum(i in chosen and j in chosen for i, j in shadow),
            member_ids=members.tolist())
    max_acc = 0.; min_margins = np.full(2, np.inf)
    max_sail = 0.; modified = deployed = 0; sail_integrals = np.zeros(count)
    last_normal = last_sail = None; max_rate = 0.; max_normal_jump = 0.
    # Controls are checked on the propagation output grid, independently of the
    # much more costly finite-Sun polygon quadrature.
    for k in range(len(t)):
        d = fleet_acceleration(env, t[k], state[k], 4*K.MOON_RADIUS, side, diagnostics=True)
        max_acc = max(max_acc, float(np.linalg.norm(d['acceleration'], axis=-1).max()))
        min_margins = np.minimum(min_margins, d['margins'].min(axis=0))
        modified += int(np.count_nonzero(d['beam_guard_changed']))
        deployed += int(np.count_nonzero(d['projected'] > 1e-7))
        sail = np.linalg.norm(d['sail'], axis=-1)
        max_sail = max(max_sail, float(sail.max()))
        if last_normal is not None:
            turn = np.arccos(np.clip(np.abs(np.sum(last_normal*d['normal'], axis=-1)), 0, 1))
            max_rate = max(max_rate, float(np.rad2deg(turn).max()/(t[k]-t[k-1])))
            max_normal_jump = max(max_normal_jump, float(np.rad2deg(turn).max()))
            sail_integrals += (last_sail+sail)*(t[k]-t[k-1])/2
        last_normal, last_sail = d['normal'], sail
    if max_acc > .15 or np.any(min_margins < -1e-10):
        raise ValueError(f'Acceleration bound or beam guard violated: {max_acc}, {min_margins}')
    results = {}
    for label, members in subsets.items():
        records = []
        for k, i in enumerate(index):
            s = state[i, members]
            d = fleet_acceleration(env, t[i], s, 4*K.MOON_RADIUS, side, diagnostics=True)
            records.append(coverage_snapshot(s, d['normal'], env.at(t[i]), side,
                4*K.MOON_RADIUS, sun_count=suns, rotation=rotation))
        results[label] = dict(**selection[label], summary=summarize_coverage(records, selected_times),
            records=records)
        print(json.dumps(dict(stage='coverage', name=name, subset=label,
            count=len(members), summary=results[label]['summary'])), flush=True)
    archive.close()
    return dict(config=config, raw_sha256=meta['sha256'], raw_metadata_sha256=digest(RUN/(name+'.json')),
        independent_initial_state_sha256=meta['initial_state_sha256'], range_km=meta['range_km'],
        collision_exclusion_pairs=len(collision), possible_mutual_shadow_pairs=len(shadow),
        minimum_pair_distance_lower_bound_m=min(collision.values(), default=None),
        example_collision_pairs=[dict(members=list(k), distance_lower_bound_m=v)
            for k, v in sorted(collision.items(), key=lambda kv: kv[1])[:12]],
        constraints=dict(physical_side_m=side, clear_optical_side_m=side-110.,
            seam_width_m=10., edge_placement_allowance_m=50., collision_clearance_m=100.,
            collision_sphere_pair_distance_m=np.sqrt(2)*side+100.,
            acceleration_bound_m_s2=.15, measured_max_acceleration_m_s2=max_acc,
            swept_guard_interval_s=float(np.max(np.diff(t))),
            minimum_moon_beam_margin_deg=float(np.rad2deg(min_margins[0])),
            minimum_earth_beam_margin_deg=float(np.rad2deg(min_margins[1])),
            normal_pointing_allowance_rad=.001, earth_guard_altitude_m=100e3),
        control=dict(electric_impulse_m_s=0., centre_ray_sail_impulse_range_m_s=[float(sail_integrals.min()),float(sail_integrals.max())],
            maximum_sail_acceleration_m_s2=max_sail,
            maximum_sampled_normal_rate_deg_s=max_rate, maximum_sampled_normal_step_deg=max_normal_jump,
            earth_or_moon_guard_modified_fraction=modified/(count*len(t)),
            deployed_fraction=deployed/(count*len(t)),
            reading_rule='Sampled command rates are lower bounds on actuator demand; no attitude torques, slew limit, flexure or hardware mass was propagated'),
        coverage_times_days=(selected_times/K.JULIAN_DAY).tolist(), selections=results,
        elapsed_analysis_s=time.monotonic()-before)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--names', nargs='+', required=True)
    p.add_argument('--cadence-hours', type=float, default=6.)
    p.add_argument('--suns', type=int, default=16)
    p.add_argument('--rotation', type=float, default=0.)
    p.add_argument('--output', default='fleet_checkpoint.json')
    a = p.parse_args()
    if Path(a.output).name != a.output:
        raise ValueError('Output must be a filename')
    sources = {**identities(), 'research/studies/solar_shield_array/fleet_analyze.py': digest(__file__),
               'research/studies/solar_shield_array/cycling_analysis.py': digest(HERE/'cycling_analysis.py')}
    cases = {}
    for name in a.names:
        cases[name] = analyze(name, a.cadence_hours*3600, a.suns, a.rotation)
    product = dict(schema='terluna.research.solar-shield-fleet/1',
        producer=dict(source_hashes=sources, constants=constants_used(sources)),
        evidence='Simultaneous independent ideal-sail IVPs in DE440; finite-square perspective shadow unions and conservative swept trajectory exclusions',
        epoch_tdb='2026-10-04', areal_mass_kg_m2=.05, protected_radius_m=4*K.MOON_RADIUS,
        parameters=vars(a), cases=cases,
        acceptance=dict(continuous_full_coverage_demonstrated=False, fleet_with_habitat_or_collector_payload_demonstrated=False,
            delivered_electrical_power_demonstrated_W=None,
            reading_rule='Nominal collision or mutual-shadow populations are rejected geometric probes. Screened subsets remove whole initial trajectories. Coverage is resolved for the stated solar quadrature and dates; no interpolation grants continuous seam closure.'),
        limitations=['Uniform solar disk; limb and spectral structure absent',
            'Central-ray SRP force; finite-Sun geometry is resolved for coverage and beam exclusion',
            'First-order tile retardation, reception-time attitude; 50 m edge margin retained',
            'Earth eclipse silhouette is a perspective-scaled circle; partial-contact geometry needs refinement before coverage acceptance',
            'Beam guard checked at output dates and continuously commanded during integration; no bounded-slew attitude actuator',
            'Conservative swept collision/shadow guards assume .15 m/s2 acceleration and 3e-7 rad/s Sun-direction bound',
            'No control/power hardware, collectors, habitats, replacements or electricity delivery included',
            'Excluding shadow pairs permits reuse of independent centre-ray trajectories; excluded full populations require mutual-radiation-force integration'],
        units=dict(position='m', velocity='m/s', time='seconds or TDB days as named', mass='kg', area='m2'))
    (HERE/'results'/a.output).write_text(compact_series_json(product))
    print('Wrote '+a.output, flush=True)


if __name__ == '__main__':
    main()
