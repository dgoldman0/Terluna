"""Inventory optical energy on the saved service and departure trajectories.

Ten-minute, one-process postprocessing budget; no trajectory search or new
hardware is propagated. Existing source-bound products remain unchanged.
"""
import json
import resource
import signal
import time

import numpy as np
from scipy.interpolate import CubicHermiteSpline
from scipy.spatial.transform import Rotation

from shared.provenance import constants_used, constants_changed
from protection.dynamics.collection import uneclipsed_collection, integrate_inventory
from protection.dynamics.parallel_shadow import ParallelPattern
from protection.dynamics.pattern_departure import departure_command, sun_frame
from protection.dynamics.optical import length
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE, SCENARIO

RUN = ROOT/'research/runs/solar_shield_array/collection'
FIELDS = ['unshadowed_aperture_W', 'first_intercept_bolometric_equivalent_W',
          'front_first_intercept_W', 'back_first_intercept_W',
          'redirected_band_optical_W', 'mutual_shadow_equivalent_loss_W']


def summarize(times, ledger, phases):
    summaries = {}
    for name, (start, end) in phases.items():
        mask = (times >= start) & (times <= end)
        row = dict(start_s=start, end_s=end, duration_s=end-start)
        for field in FIELDS:
            data = integrate_inventory(times[mask], ledger[field][mask])
            members = data.pop('per_member_energy_J')
            data['per_member_energy_range_J'] = [min(members), max(members)]
            row[field] = data
        summaries[name] = row
    return summaries


def main():
    signal.alarm(600)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    before = time.monotonic(); RUN.mkdir(parents=True, exist_ok=True)
    names = ['patterns.json', 'departure_validation.json', 'departure_refine.json']
    parents = {name: json.loads((HERE/'results'/name).read_text()) for name in names}
    sources = {}
    for parent in parents.values():
        for path, expected in parent['producer']['source_hashes'].items():
            if digest(ROOT/path) != expected: raise ValueError('Source changed: '+path)
        if constants_changed(parent['producer']['constants']): raise ValueError('Constants changed')
        sources.update(parent['producer']['source_hashes'])
    raw = []; raw_inputs = {}
    for entry in [parents['patterns.json']['replay'], parents['departure_validation.json']['runs'][-1]]:
        if digest(ROOT/entry['raw_path']) != entry['raw_sha256']: raise ValueError('Raw changed')
        raw_inputs[entry['raw_path']] = entry['raw_sha256']
        raw.append(np.load(ROOT/entry['raw_path']))
    if not np.array_equal(raw[0]['state'][-1], raw[1]['state'][0]):
        raise ValueError('Service and departure do not share the actual terminal states')
    for path in ['protection/dynamics/collection.py',
                 'research/studies/solar_shield_array/collection_run.py']:
        sources[path] = digest(ROOT/path)
    env = environment(1.); side = 10000.
    start, join, end = float(raw[0]['t'][0]), float(raw[0]['t'][-1]), float(raw[1]['t'][-1])
    r = parents['departure_refine.json']
    command = departure_command(env, join, end, r['delay_s'], r['slew_s'], r['axis'], r['sign'])
    interpolations = [CubicHermiteSpline(a['t'], a['state'][:, :, :3], a['state'][:, :, 3:]) for a in raw]
    times = np.arange(start, end+1., 60.)
    phases = dict(service=(start, join), departure_preparation=(join, join+r['delay_s']),
        departure_turn=(join+r['delay_s'], join+r['delay_s']+r['slew_s']),
        departure_coast=(join+r['delay_s']+r['slew_s'], end), whole_saved_interval=(start, end))
    out = dict(schema='terluna.research.solar-collection-inventory/1',
        producer=dict(source_hashes=sources, constants=constants_used(sources),
            inputs={name: digest(HERE/'results'/name) for name in names}, raw_inputs=raw_inputs),
        epoch_tdb=SCENARIO['epoch_tdb'], members=len(raw[0]['state'][0]),
        physical_tile_side_m=side, clear_tile_side_m=side-110.,
        base_optical_areal_mass_kg_m2=env.sigma,
        redirected_spectral_fraction=env.central_fraction,
        budget=dict(wall_s=600, processes=1, numerical_threads=1, address_space_GiB=2,
            raw_output_limit_MiB=128, source_counts=[16, 32], time_step_s=60),
        evidence='Postprocessed actual saved states; sampled finite-Sun optical energy with shared-plane mutual shadows. No collecting hardware, new trajectory or complete fleet has been simulated.',
        quantity_scope='First-intercept bolometric equivalent counts a ray at its first clear aperture; it is not filter absorption. The redirected-band ledger uses the uniform ideal spectral fraction already used for sail force. Transmitted rays receive no collection credit.',
        eclipse_scope='Conservative whole-tile Earth/Moon cone guard; all evaluated dates must be uneclipsed. Partial or total body eclipse rejects this prefix implementation and requires joint source/area visibility for force and collection.',
        geometry_scope='Same centre-ray irradiance/incidence approximation and perspective rectangle unions as the force model. Clear apertures only; supports and edges receive no collection credit. Both faces are tabulated separately.',
        assignment_scope='Phases are command stages, not an assertion that all intercepted rays serve the local receiver. Service assignments and useful coverage remain the separate parent-product measurements.',
        interpretation='Only sum simultaneously instantiated, mutually checked members. Never scale this local patch by a relaxed fleet inventory or sum incompatible trajectory columns as though they were coexisting tiles.',
        electrical_conversion_efficiency=None, outgoing_light_capture_fraction=None,
        collection_capacity_W=None, delivered_power_W=None, fleet_propulsion_W=None,
        hardware_mass_feedback_propagated=False, continuous_certificate=False,
        accepted_fleet=False, recurring_cycle_measured=False, stages=[])
    def save():
        out.update(elapsed_s=time.monotonic()-before,
            max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
        (HERE/'results/collection.json').write_text(json.dumps(out, indent=2, allow_nan=False)+'\n')
    save()
    ledgers = []; momentum_error = 0.
    for suns, rotation in [(16, .317), (32, .317+np.pi/32)]:
        ledger = {field: [] for field in FIELDS}
        for k, t in enumerate(times):
            interpolation = interpolations[0 if t <= join else 1]
            state = np.c_[interpolation(t), interpolation(t, 1)]
            frame = sun_frame(env, t) if t <= join else command(t).as_matrix()
            result = uneclipsed_collection(env, t, state, frame, suns=suns, rotation=rotation)
            for field in FIELDS: ledger[field].append(result[field])
            if suns == 16 and k % 60 == 0:
                # Same actual frame, all-source visibility and spectral fraction
                # as the saved force model; energy uses |cos| and force cos².
                constant_frame = lambda unused: Rotation.from_matrix(frame)
                model = ParallelPattern(env, constant_frame, suns=16)
                sail = model.acceleration(t, state, diagnostics=True)[2]
                error = length(result['reflected_force_N']/(env.sigma*side**2)-sail).max()
                momentum_error = max(momentum_error, float(error))
            if k % 180 == 0:
                print(json.dumps(dict(stage='collection', suns=suns, hours=t/3600,
                    elapsed_s=time.monotonic()-before)), flush=True)
        ledger = {key: np.asarray(value) for key, value in ledger.items()}
        ledgers.append(ledger)
        out['stages'].append(dict(sources=suns, disk_rotation_rad=rotation,
            sample_step_s=60., phases=summarize(times, ledger, phases)))
        save()
    coarse_time = summarize(times[::2], {k: v[::2] for k, v in ledgers[0].items()}, phases)
    validation = {}
    for phase in phases:
        field = 'redirected_band_optical_W'
        e16 = out['stages'][0]['phases'][phase][field]['energy_J']
        e32 = out['stages'][1]['phases'][phase][field]['energy_J']
        validation[phase] = dict(source_refinement_relative_energy_change=abs(e32-e16)/e32,
            time_coarsening_relative_energy_change=abs(coarse_time[phase][field]['energy_J']-e16)/e16)
    path = RUN/'service_departure.npz'
    np.savez_compressed(path, t=times, member_id=np.arange(out['members']),
        **ledgers[1], **{'source16_'+k: v for k, v in ledgers[0].items()})
    out.update(raw_path=str(path.relative_to(ROOT)), raw_sha256=digest(path),
        raw_MiB=path.stat().st_size/1024**2,
        selected_source_count=32, validation=validation,
        maximum_force_ledger_disagreement_m_s2=momentum_error,
        force_check_dates=13,
        refinement_scope='Collection quadrature only on the same saved trajectories; no 32-source dynamics replay. Time coarsening compares 60 and 120 s optical integrals. Extrema remain sampled.',
        completed=True)
    if momentum_error > 1e-14: raise ValueError('Optical energy and sail momentum disagree')
    if path.stat().st_size > 128*1024**2: raise ValueError('Raw output budget exceeded')
    if any(digest(ROOT/p) != h for p, h in sources.items()): raise ValueError('Source changed during run')
    save()
    print(json.dumps(dict(stage='complete', elapsed_s=out['elapsed_s'],
        max_rss_MiB=out['max_rss_MiB'], validation=validation,
        force_disagreement_m_s2=momentum_error)), flush=True)


if __name__ == '__main__':
    main()
