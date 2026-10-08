"""Bounded common-epoch orbit library and optimistic coverage-capacity screen."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import resource
import signal
import time

import numpy as np
from scipy.integrate import solve_ivp

from shared import constants as K
from shared.provenance import constants_used, constants_changed
from protection.dynamics.active_retime import finite_tile_acceleration
from protection.dynamics.coverage_pilot import (seed_library, receiver_cells, capacity_rows,
    cover_capacity, conditional_cost, orientation_matrices, rotation_step_angle)
from protection.dynamics.fleet import fleet_acceleration
from protection.dynamics.optical import length
from .fleet_run import environment, digest
from .cycling_search import ROOT, HERE, SCENARIO

RUN = ROOT/'research/runs/solar_shield_array/pilot'
SOURCES = ['protection/dynamics/'+x+'.py' for x in
    ['coverage_pilot','active_retime','fleet','cycling','ephemeris','optical','model']]+[
    'research/studies/solar_shield_array/'+x+'.py' for x in ['pilot_run','fleet_run','cycling_search']]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--phases', type=int, choices=[24,48], default=24)
    p.add_argument('--days', type=float, default=3.)
    p.add_argument('--max-step', type=float, default=600.)
    p.add_argument('--rtol', type=float, default=2e-11)
    p.add_argument('--reanalyze', action='store_true')
    a = p.parse_args()
    if not 0 < a.days <= 3:
        raise ValueError('Pilot horizon budget is three days')
    signal.alarm(1200)
    resource.setrlimit(resource.RLIMIT_AS, (2*1024**3, 2*1024**3))
    started = time.monotonic()
    source_hashes = {n:digest(ROOT/n) for n in SOURCES}
    env = environment(a.days)
    initial, families = seed_library(env, np.array([12,15,20])*1e6,
        np.arange(6)*np.pi/6, [-.15,0.,.15], a.phases)
    times = np.linspace(0,a.days*K.JULIAN_DAY,round(a.days*K.JULIAN_DAY/600)+1)
    shape = initial.shape
    last_log = [time.monotonic()]
    def rhs(t,y):
        if time.monotonic()-last_log[0] > 30:
            print(json.dumps(dict(stage='propagation',day=t/K.JULIAN_DAY,
                                 elapsed_s=time.monotonic()-started)),flush=True)
            last_log[0] = time.monotonic()
        state = y.reshape(shape)
        return np.c_[state[:,3:], finite_tile_acceleration(env,t,state)].ravel()
    reused = None
    if a.reanalyze:
        previous = json.loads((HERE/'results'/f'pilot_{a.phases}.json').read_text())
        for name, value in previous['producer']['source_hashes'].items():
            if name not in ['protection/dynamics/coverage_pilot.py',
                            'research/studies/solar_shield_array/pilot_run.py'] and digest(ROOT/name)!=value:
                raise ValueError('Dynamics source changed: '+name)
        if constants_changed(previous['producer']['constants']):
            raise ValueError('Dynamics constants changed')
        raw = ROOT/previous['raw_path']
        if digest(raw)!=previous['raw_sha256']:
            raise ValueError('Raw trace identity changed')
        saved = np.load(raw)
        if not np.array_equal(saved['initial'],initial) or not np.array_equal(saved['t'],times):
            raise ValueError('Initial states or output times changed')
        states = saved['state']; nfev = previous['nfev']
        reused = dict(raw_sha256=previous['raw_sha256'],
            original_producer=previous['producer'],original_elapsed_s=previous['elapsed_s'],
            previous_reuse=previous.get('trajectory_reuse'),
            reason='Refresh capacity/attitude reduction from verified trajectories; unchanged force code, constants and initial states',
            previous_max_labelled_frame_rate_deg_s=previous['maximum_secant_attitude_rate_deg_s'])
    else:
        sol = solve_ivp(rhs,[0,times[-1]],initial.ravel(),method='DOP853',t_eval=times,
            rtol=a.rtol,atol=np.tile([1e-4]*3+[1e-8]*3,len(initial)),max_step=a.max_step)
        if not sol.success:
            raise RuntimeError(sol.message)
        states = sol.y.T.reshape(-1,*shape); nfev = sol.nfev
    print(json.dumps(dict(stage='propagated',members=len(initial),nfev=nfev,
                         elapsed_s=time.monotonic()-started)),flush=True)
    normals=[]; beams=[]; orientations=[]; sails=[]; errors=[]
    for t,state in zip(times,states):
        _,d,error=finite_tile_acceleration(env,t,state,diagnostics=True)
        normals.append(d['normal']); beams.append(d['margins']); sails.append(d['sail'])
        errors.append(error)
        orientations.append(orientation_matrices(state,d['normal'],d['sun']))
    normals=np.asarray(normals); beams=np.asarray(beams); orientations=np.asarray(orientations)
    rates=np.rad2deg(rotation_step_angle(orientations[:-1],orientations[1:]))/np.diff(times)[:,None]
    radius=length(states[...,:3]); protected=4*K.MOON_RADIUS
    family_min_beam=beams.reshape(len(times),len(families),a.phases,2).min(axis=(0,2,3))
    family_max_rate=rates.reshape(len(times)-1,len(families),a.phases).max(axis=(0,2))
    # A finite-difference rate exceeding this cap proves the command needs
    # refinement. A rate below it does not certify the unsampled command.
    admissible=(family_min_beam>=0)&(family_max_rate<=.1)
    cells=receiver_cells(protected)
    matrices=[]; row_labels=[]
    sample_ids=np.flatnonzero(np.isclose(times%(3*3600),0,atol=1e-6))
    for k in sample_ids:
        matrix,labels=capacity_rows(states[k],normals[k],env.at(times[k]),cells,
            len(families),a.phases,limb_count=4)
        matrices.append(matrix)
        row_labels += [[float(times[k]/K.JULIAN_DAY),*label] for label in labels]
    matrix=np.concatenate(matrices)
    coarse=np.array([abs(day/.25-round(day/.25))<1e-8 for day,_,_ in row_labels])
    results={'six_hour':cover_capacity(matrix[coarse]),'three_hour':cover_capacity(matrix),
             'sampled_attitude_gate':cover_capacity(matrix,admissible)}
    for name,item in results.items():
        labels=np.asarray(row_labels)[coarse] if name=='six_hour' else np.asarray(row_labels)
        for key in ['bottleneck_rows','empty_rows']:
            if key in item:
                item[key+'_labels']=labels[item[key]].tolist()
    best=results['three_hour']
    cost=[]
    if best['success']:
        for multiplier,duty in [(1.,1.),(1.,.1),(1.2,.1)]:
            for dv in [1.,3.,5.,10.]:
                cost.append(conditional_cost(best['inventory_tiles'],dv,SCENARIO['propulsion'],multiplier,duty))
    # Save all initial states and returned states for refinement and witnesses.
    RUN.mkdir(exist_ok=True,parents=True)
    raw=RUN/f'phases_{a.phases}.npz'
    np.savez_compressed(raw,t=times,state=states,normal=normals,beam=beams,
        rate_deg_s=rates,matrix=matrix,row_labels=np.asarray(row_labels),initial=initial)
    if any(digest(ROOT/n)!=h for n,h in source_hashes.items()):
        raise ValueError('Producer changed during run')
    product=dict(schema='terluna.research.coverage-capacity-pilot/1',
        producer=dict(source_hashes=source_hashes,constants=constants_used(source_hashes)),
        evidence='Restricted-library summed-area capacity relaxation; independent sail-assisted IVPs; no accepted fleet',
        epoch_tdb=SCENARIO['epoch_tdb'],config=vars(a),families=families,members=len(initial),
        target_radius_m=protected,tile_side_m=10000.,areal_mass_kg_m2=.05,
        phase_copies_are_independent_ivps=True,translational_electric_acceleration=0.,
        mutually_shadowed_fleet_dynamics=False,finite_extent_gravity_quadrature='3x3 Gauss-Legendre',
        radius_range_km=[float(radius.min()/1000),float(radius.max()/1000)],
        radial_guard_passed=bool(radius.min()>protected+10000/np.sqrt(2)),
        finite_extent_acceleration_max_m_s2=float(np.max(errors)),
        minimum_beam_margin_deg=float(np.rad2deg(beams.min())),
        maximum_secant_attitude_rate_deg_s=float(rates.max()),
        family_min_beam_margin_rad=family_min_beam.tolist(),family_max_secant_rate_deg_s=family_max_rate.tolist(),
        attitude_rate_cap_deg_s=.1,passing_sampled_attitude_families=int(admissible.sum()),
        continuous_attitude_safety_certified=False,continuous_coverage_certified=False,
        finite_sun_constraints='central ray and four limb directions; no irradiance quadrature implied',
        cell_partition_area_error=1-sum(x.area for x in cells)/(np.pi*protected**2),
        covering=results,conditional_cost_scenarios=cost,
        benchmark_TW=238.35,collection_capacity_W=None,delivered_power_W=None,
        raw_path=str(raw.relative_to(ROOT)),raw_sha256=digest(raw),nfev=nfev,
        trajectory_reuse=reused,
        elapsed_s=time.monotonic()-started,max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024,
        reading_rule='Weights replicate a sampled phase measure and permit cellwise area repacking. They give no tile placement, phase-continuum certificate, collision-free fleet, recurring maneuver rate or globally valid inventory lower bound.')
    out=HERE/'results'/f'pilot_{a.phases}.json'
    out.write_text(json.dumps(product,indent=2)+'\n')
    print(json.dumps(dict(product=str(out),elapsed_s=product['elapsed_s'],
        covering={k:{key:val for key,val in v.items() if key in ['success','reason','inventory_tiles']} for k,v in results.items()},
        attitude_families=int(admissible.sum()),max_rate=float(rates.max()))),flush=True)


if __name__=='__main__':
    main()
