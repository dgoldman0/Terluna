"""Refine capacity bottlenecks and independently replay pilot witnesses."""
from __future__ import annotations
import json
import resource
import signal
import time

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.spatial import cKDTree

from shared import constants as K
from shared.provenance import constants_used, constants_changed
from protection.dynamics.active_retime import finite_tile_acceleration
from protection.dynamics.coverage_pilot import (receiver_cells, capacity_rows, cover_capacity,
    conditional_cost, orientation_matrices, rotation_step_angle, seed_library)
from protection.dynamics.cycling import sail_basis
from protection.dynamics.fleet import square_vertices
from protection.dynamics.active_global import square_pair_separation
from protection.dynamics.fleet_exclusions import swept_conflicts
from protection.dynamics.optical import length, unit
from .fleet_run import environment, digest
from .cycling_search import ROOT,HERE,SCENARIO


def load(phases):
    product=json.loads((HERE/'results'/f'pilot_{phases}.json').read_text())
    for n,h in product['producer']['source_hashes'].items():
        if digest(ROOT/n)!=h: raise ValueError('Pilot source changed: '+n)
    if constants_changed(product['producer']['constants']): raise ValueError('Constants changed')
    raw=ROOT/product['raw_path']
    if digest(raw)!=product['raw_sha256']: raise ValueError('Raw trace changed')
    return product,np.load(raw)


def main():
    signal.alarm(300)
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))
    start=time.monotonic()
    coarse,c=load(24); fine,f=load(48)
    phases=48; original_families=len(fine['families']); families=original_families
    t=f['t']; states=f['state']; normals=f['normal']; matrix=f['matrix']
    env=environment(fine['config']['days'])
    best=fine['covering']['three_hour']
    labels=f['row_labels']; generation=None
    base_states=states
    if not best['success']:
        witness=labels[best['empty_rows'][0]]
        sector=int(witness[2])%8
        # The plane normal is perpendicular to its projected track. Add
        # columns aimed at the witnessed sector, with independent phases.
        azimuth=((sector+.5)*2*np.pi/8+np.pi/2)%np.pi
        extra,descriptions=seed_library(env,np.array([15.,20.])*1e6,[azimuth],[-.1,.1],phases)
        def extra_rhs(ti,y):
            state=y.reshape(-1,6)
            return np.c_[state[:,3:],finite_tile_acceleration(env,ti,state)].ravel()
        solved=solve_ivp(extra_rhs,[0,t[-1]],extra.ravel(),method='DOP853',t_eval=t,
            max_step=600.,rtol=2e-11,atol=np.tile([1e-4]*3+[1e-8]*3,len(extra)))
        if not solved.success: raise RuntimeError(solved.message)
        extra_states=solved.y.T.reshape(len(t),len(extra),6)
        extra_normals=np.array([finite_tile_acceleration(env,ti,yi,diagnostics=True)[1]['normal']
                               for ti,yi in zip(t,extra_states)])
        extra_rows=[]
        for ki in np.flatnonzero(np.isclose(t%(3*3600),0,atol=1e-6)):
            rows,_=capacity_rows(extra_states[ki],extra_normals[ki],env.at(t[ki]),
                receiver_cells(4*K.MOON_RADIUS),len(descriptions),phases,limb_count=4)
            extra_rows.append(rows)
        matrix=np.c_[matrix,np.concatenate(extra_rows)]
        states=np.concatenate([states,extra_states],axis=1)
        normals=np.concatenate([normals,extra_normals],axis=1)
        families+=len(descriptions)
        best=cover_capacity(matrix)
        generation=dict(trigger_row=witness.tolist(),added_members=len(extra),families=descriptions,
                        result=best,nfev=solved.nfev)
        if not best['success']: raise RuntimeError('Targeted columns leave empty capacity rows; inspect saved pilot before extending')
    uniform_solution=best
    # Four independently populated initial-phase bands per plane replace the
    # uniform population constraint, within the 256-column pilot budget.
    orientations=[]; beam_min=[]
    for ti,yi in zip(t,states):
        _,d,_=finite_tile_acceleration(env,ti,yi,diagnostics=True)
        orientations.append(orientation_matrices(yi,d['normal'],d['sun']))
        beam_min.append(d['margins'].min(axis=1))
    orientations=np.array(orientations)
    rates=np.rad2deg(rotation_step_angle(orientations[:-1],orientations[1:]))/np.diff(t)[:,None]
    phases//=4;families*=4
    rows=[]
    for ki in np.flatnonzero(np.isclose(t%(3*3600),0,atol=1e-6)):
        block,_=capacity_rows(states[ki],normals[ki],env.at(t[ki]),
            receiver_cells(4*K.MOON_RADIUS),families,phases,limb_count=4)
        rows.append(block)
    matrix=np.concatenate(rows)
    free_phase_solution=cover_capacity(matrix)
    rate_limit=rates.reshape(len(t)-1,families,phases).max(axis=(0,2))
    beam_limit=np.array(beam_min).reshape(len(t),families,phases).min(axis=(0,2))
    admissible=(rate_limit<=.1)&(beam_limit>=0)
    gated_phase_solution=cover_capacity(matrix,admissible)
    best=gated_phase_solution if gated_phase_solution['success'] else free_phase_solution
    x=np.array(best['weights_million_tiles'])
    selected=np.flatnonzero(x>1e-6)
    # Refine the original optimum's most expensive date/cell; source directions
    # increase and every angular sector is bisected. A positive witness is an
    # actual violation of this necessary capacity condition.
    bottleneck=int(best['bottleneck_rows'][0]); day=float(labels[bottleneck,0])
    k=int(np.argmin(abs(t-day*K.JULIAN_DAY)))
    cells=receiver_cells(4*K.MOON_RADIUS,sectors=16)
    new,local_labels=capacity_rows(states[k],normals[k],env.at(t[k]),cells,families,phases,limb_count=8)
    capacity=new@x
    refined=cover_capacity(np.r_[matrix,new],admissible if gated_phase_solution['success'] else None)
    if refined['success']:
        x=np.array(refined['weights_million_tiles'])
        selected=np.flatnonzero(x>1e-6)
    # Global library diagnostics do not give a placement for fractional weights.
    # Examine actual sparse representatives from the selected phase families.
    member_ids=np.concatenate([np.arange(j*phases,(j+1)*phases) for j in selected])
    active=states[:,member_ids]
    directions=np.array([unit(env.at(ti)['positions']['sun']) for ti in t])
    collision,shadow=swept_conflicts(t,active,10000.,acceleration_bound=.06,
        clearance=100.,sun_directions=directions)
    # Exact static SAT witnesses complement conservative swept flags.
    intersections=[]
    for ki in range(len(t)):
        q=active[ki,:,:3]
        pairs=cKDTree(q).query_pairs(np.sqrt(2)*10000,output_type='ndarray')
        if len(pairs):
            _,d,_=finite_tile_acceleration(env,t[ki],active[ki],diagnostics=True)
            rot=orientation_matrices(active[ki],d['normal'],d['sun'])
            sep=square_pair_separation(q,rot[:,:,0],rot[:,:,1],pairs)
            for pair,s in zip(pairs[sep<=0],sep[sep<=0]):
                intersections.append(dict(day=float(t[ki]/K.JULIAN_DAY),
                    ids=member_ids[pair].tolist(),maximum_separating_projection_m=float(s)))
            if len(intersections)>=8: break
    # Highest-rate members, representatives of the selected library, and any
    # exact intersection witnesses get their own tighter IVPs.
    rate_order=member_ids[np.argsort(rates[:,member_ids].max(axis=0))[::-1]]
    ids=list(rate_order[:4])
    ids += [int(j*phases) for j in selected[:8]]
    ids += [i for witness in intersections for i in witness['ids']]
    ids=np.array(list(dict.fromkeys(ids))[:24],int)
    initial=states[0,ids]
    def rhs(ti,y):
        state=y.reshape(-1,6)
        return np.c_[state[:,3:],finite_tile_acceleration(env,ti,state)].ravel()
    replay=solve_ivp(rhs,[0,t[-1]],initial.ravel(),method='DOP853',t_eval=t,dense_output=True,
        max_step=120.,rtol=2e-12,atol=np.tile([1e-5]*3+[1e-9]*3,len(ids)))
    if not replay.success: raise RuntimeError(replay.message)
    independent=replay.y.T.reshape(len(t),len(ids),6)
    position=length(independent[:,:,:3]-states[:,ids,:3])
    velocity=length(independent[:,:,3:]-states[:,ids,3:])
    # Follow the most abrupt selected attitude at 10-second cadence. Smooth
    # interpolation of snapshots cannot replace the actual command function.
    where=np.unravel_index(np.argmax(rates[:,ids]),rates[:,ids].shape)
    ki,li=where; witness_id=int(ids[li])
    local_t=np.arange(max(0,t[ki]-600),min(t[-1],t[ki+1]+600)+1,10.)
    local_states=replay.sol(local_t).T.reshape(-1,len(ids),6)[:,li:li+1]
    rots=[]; beam=[]; local_normals=[]
    for ti,yi in zip(local_t,local_states):
        _,d,_=finite_tile_acceleration(env,ti,yi,diagnostics=True)
        rots.append(orientation_matrices(yi,d['normal'],d['sun'])[0]);beam.append(d['margins'][0])
        local_normals.append(d['normal'][0])
    local_rates=np.rad2deg(rotation_step_angle(np.asarray(rots[:-1]),np.asarray(rots[1:])))/10
    local_normals=np.array(local_normals)
    normal_rates=np.rad2deg(np.arccos(np.clip(abs(np.sum(local_normals[:-1]*local_normals[1:],axis=1)),0,1)))/10
    # Reflected beam safety along an interpolated slew is a separate gate;
    # command rate violations already suffice to reject this law at the cap.
    costs=[]
    if refined['success']:
        for multiplier,duty in [(1.,1.),(1.,.1),(1.2,.1)]:
            n=refined['inventory_tiles']
            def residual(dv):
                p=conditional_cost(n,dv,SCENARIO['propulsion'],multiplier,duty)
                return p['propulsion_TW']-238.35 if p['closed'] else 1e12
            dv=brentq(residual,0.,100.)
            costs.append(dict(fixed_mass_multiplier=multiplier,peak_duty=duty,
                break_even_dv_m_s_per_day=float(dv),
                cases=[conditional_cost(n,a,SCENARIO['propulsion'],multiplier,duty) for a in [1,3,5,10]]))
    # Coarse phases are nested in the fine library, at identical initial states.
    nested=base_states.reshape(len(t),original_families,48,6)[:,:,::2].reshape(c['state'].shape)
    repeated_error=float(length(nested[...,:3]-c['state'][...,:3]).max())
    new_labels=[[day,*label] for label in local_labels]
    worst=int(np.argmin(capacity))
    sources={**fine['producer']['source_hashes'],
        'research/studies/solar_shield_array/pilot_validate.py':digest(__file__),
        'protection/dynamics/fleet_exclusions.py':digest(ROOT/'protection/dynamics/fleet_exclusions.py'),
        'protection/dynamics/active_global.py':digest(ROOT/'protection/dynamics/active_global.py')}
    out=dict(schema='terluna.research.coverage-capacity-validation/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources),
            inputs={f'pilot_{p}.json':digest(HERE/'results'/f'pilot_{p}.json') for p in [24,48]}),
        targeted_column_generation=generation,
        phase_population_optimization=dict(bands_per_plane=4,phases_per_column=phases,
            columns=families,uniform_solution=uniform_solution,free_phase_solution=free_phase_solution,
            sampled_attitude_gated_solution=gated_phase_solution,
            passing_columns=int(admissible.sum())),
        phase_refinement=dict(coarse_inventory=coarse['covering']['three_hour'].get('inventory_tiles'),
            fine_inventory=fine['covering']['three_hour'].get('inventory_tiles'),
            selected_library_inventory=best['inventory_tiles'],nested_trajectory_max_difference_m=repeated_error),
        constraint_generation=dict(date_days=day,original_cell_sectors=8,refined_cell_sectors=16,
            original_limb_sources=4,refined_limb_sources=8,added_rows=len(new),
            minimum_old_solution_capacity=float(capacity.min()),worst_row=new_labels[worst],
            violated_rows=int(np.count_nonzero(capacity<1-1e-7)),refined_solution=refined,
            interpretation='Spatial/source adversarial refinement of the strongest dual-priced date; other dates retain the coarser partition.'),
        selected_sparse_library=dict(families=selected.tolist(),members=len(member_ids),
            inspected_inventory_tiles=float(x.sum()*1e6),
            possible_swept_collision_pairs=len(collision),possible_swept_shadow_pairs=len(shadow),
            acceleration_bound_m_s2=.06,snapshot_square_intersection_witnesses=intersections,
            reading_rule='These are the sparse representatives of the final refined allocation, not a physical realization of fractional multiplicities. Possible shadow pairs require coupled illumination before simultaneous acceptance.'),
        independent_replay=dict(members=ids.tolist(),max_step_s=120.,rtol=2e-12,
            max_position_difference_m=float(position.max()),max_velocity_difference_m_s=float(velocity.max()),
            passes_50m=bool(position.max()<50),nfev=replay.nfev),
        attitude_witness=dict(member=witness_id,selected_phase_column=int(witness_id//phases),
            selected_weight_million_tiles=float(x[witness_id//phases]),
            start_day=float(local_t[0]/K.JULIAN_DAY),
            end_day=float(local_t[-1]/K.JULIAN_DAY),cadence_s=10.,
            max_command_rate_deg_s=float(local_rates.max()),assumed_cap_deg_s=.1,
            max_normal_plane_rate_deg_s=float(normal_rates.max()),
            minimum_commanded_beam_margin_deg=float(np.rad2deg(np.min(beam))),
            slew_path_safety_certified=False),
        conditional_cost=costs,accepted_fleet=False,recurring_correction_measured=False,
        elapsed_s=time.monotonic()-start,max_rss_MiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024)
    raw=ROOT/'research/runs/solar_shield_array/pilot/validation.npz'
    np.savez_compressed(raw,t=t,initial=states[0],matrix=matrix,added_matrix=new,
        selected_weights=x,replay_ids=ids,replay=independent,
        attitude_t=local_t,attitude_state=local_states,attitude_rates=local_rates)
    out['raw_path']=str(raw.relative_to(ROOT));out['raw_sha256']=digest(raw)
    out['elapsed_s']=time.monotonic()-start
    if any(digest(ROOT/n)!=h for n,h in sources.items()): raise ValueError('Source changed during validation')
    (HERE/'results/pilot_validation.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out),flush=True)


if __name__=='__main__':main()
