"""Shoot a single finite tile to an advanced orbital phase at the same date."""
import argparse
import json
from pathlib import Path
import time

import numpy as np

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.active_formation import solar_frame
from protection.dynamics.active_retime import propagate,rotate_state,command,finite_tile_acceleration
from protection.dynamics.cycling import ray_geometry
from protection.dynamics.optical import unit,length
from .fleet_run import environment,identities,digest
from .cycling_search import HERE,ROOT,SCENARIO


def run(days,displacement,burn_hours=6,arcs=4):
    before=time.monotonic(); duration=days*K.JULIAN_DAY; burn=burn_hours*3600
    env=environment(days); basis=solar_frame(env.at(0)); u,b=basis[:,:2].T
    if np.dot(np.cross(u,b),env.frames(0)[:,2])>0:b=-b
    radius=15e6; frequency=np.sqrt(K.MOON_GM/radius**3); phase=-frequency*duration
    # Select an actual initial state that ends near the sunward passage. Each
    # replay is integrated from the common epoch in the same ephemeris.
    alignment=[]
    for _ in range(4):
        initial=np.r_[radius*(u*np.cos(phase)+b*np.sin(phase)),
            radius*frequency*(-u*np.sin(phase)+b*np.cos(phase))]
        rt,reference,_=propagate(env,initial,duration,np.zeros(3*arcs),burn)
        target_u=solar_frame(env.at(duration))[:,0]
        q=reference[-1,:3]; h=unit(np.cross(q,reference[-1,3:]))
        angle=np.arctan2(np.dot(np.cross(q,target_u),h),np.dot(q,target_u))
        alignment.append(float(angle)); phase+=angle
        if abs(angle)<1e-3:break
    target=rotate_state(reference[-1],displacement/length(reference[-1,:3]))
    scale=np.array([1e6]*3+[1e6*frequency]*3)
    parameters=np.zeros(3*arcs); history=[]
    def end(p):return propagate(env,initial,duration,p,burn)[1][-1]
    value=reference[-1]; jac=None
    for iteration in range(12):
        residual=(value-target)/scale
        history.append(dict(iteration=iteration,position_error_m=float(length(value[:3]-target[:3])),
            velocity_error_m_s=float(length(value[3:]-target[3:])),parameters_m_s=parameters.tolist()))
        print(json.dumps(dict(days=days,**history[-1])),flush=True)
        if length(value[:3]-target[:3])<2. and length(value[3:]-target[3:])<2e-4:break
        if jac is None or iteration in (4,8):
            varied=parameters+np.eye(len(parameters))*.02
            _,traces,_=propagate(env,np.broadcast_to(initial,(len(parameters),6)),duration,varied,burn)
            jac=((traces[-1]-value)/scale/.02).T
        delta=np.linalg.lstsq(jac,-residual,rcond=1e-10)[0]
        if length(delta)>20:delta*=20/length(delta)
        for backtrack in range(6):
            candidate=parameters+delta; next_value=end(candidate)
            if length((next_value-target)/scale)<length(residual):break
            delta*=.5
        response=(next_value-value)/scale
        jac+=np.outer(response-jac@delta,delta)/max(np.dot(delta,delta),1e-100)
        parameters,value=candidate,next_value
    # Independent tighter replay with smaller steps and an actual time history.
    t,state,nfev=propagate(env,initial,duration,parameters,burn,rtol=2e-12,max_step=120.,output_step=300.)
    commands=np.array([command(ti,yi,parameters,duration,burn) for ti,yi in zip(t,state)])
    diagnostics=[finite_tile_acceleration(env,ti,yi,diagnostics=True) for ti,yi in zip(t,state)]
    margins=np.array([d[1]['margins'][0] for d in diagnostics])
    dv=float(np.linalg.norm(parameters.reshape(arcs,3),axis=1).sum())
    peak=2*float(np.linalg.norm(parameters.reshape(arcs,3),axis=1).max())/burn
    prop=SCENARIO['propulsion']; cant=np.cos(np.deg2rad(prop['cant_deg'])); mass=5e6
    energy=mass*dv*prop['exhaust_velocity_m_s']/(2*prop['efficiency']*cant)
    fuel=mass*dv/(prop['exhaust_velocity_m_s']*cant)
    sample=env.at(duration); rays=ray_geometry(state[-1,:3],sample)
    sun_u=unit(rays['sun']); shift=state[-1,:3]-reference[-1,:3]
    projection=shift-sun_u*np.dot(shift,sun_u)
    folder=ROOT/'research/runs/solar_shield_array/active'; folder.mkdir(exist_ok=True,parents=True)
    raw=folder/f'retime_{days:g}d.npz'; np.savez_compressed(raw,t=t,state=state,command=commands,
        initial=initial,target=target,reference_t=rt,reference=reference)
    result=dict(days=days,requested_phase_arc_m=displacement,burn_hours=burn_hours,arcs=arcs,
        initial_state=initial.tolist(),same_date_reference_final=reference[-1].tolist(),target_state=target.tolist(),
        parameters_m_s=parameters.tolist(),shooting_history=history,initial_phase_alignment_radians=alignment,
        independent_replay_position_error_m=float(length(state[-1,:3]-target[:3])),
        independent_replay_velocity_error_m_s=float(length(state[-1,3:]-target[3:])),
        replay_agreement_within_50m=bool(length(state[-1,:3]-target[:3])<50),
        projected_position_shift_at_arrival_m=float(length(projection)),
        final_sunward_distance_m=float(rays['sunward']),
        final_shadow_centre_distance_m=float(length(rays['transverse'])),
        delta_v_m_s=dv,peak_acceleration_m_s2=peak,optical_mass_kg=mass,
        ideal_electric_energy_J=energy,mean_electric_power_over_maneuver_W=energy/duration,
        peak_electric_power_W=mass*peak*prop['exhaust_velocity_m_s']/(2*prop['efficiency']*cant),
        propellant_kg=fuel,propellant_fraction_of_base_optical_mass=fuel/mass,
        minimum_earth_beam_margin_deg=float(np.rad2deg(margins[:,1].min())),
        minimum_moon_beam_margin_deg=float(np.rad2deg(margins[:,0].min())),
        radius_range_km=[float(length(state[:,:3]).min()/1000),float(length(state[:,:3]).max()/1000)],
        max_finite_extent_gravity_correction_m_s2=float(max(d[2] for d in diagnostics)),
        raw_path=str(raw.relative_to(ROOT)),raw_sha256=digest(raw),replay_nfev=nfev,elapsed_s=time.monotonic()-before)
    print(json.dumps(dict(final=True,**result)),flush=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--days',type=float,nargs='+',default=[3.,7.])
    p.add_argument('--displacement-km',type=float,default=1400.);p.add_argument('--burn-hours',type=float,default=6.)
    p.add_argument('--arcs',type=int,default=4)
    p.add_argument('--output',default='active_retime.json');a=p.parse_args()
    sources={**identities(),
        'protection/dynamics/active_formation.py':digest(ROOT/'protection/dynamics/active_formation.py'),
        'protection/dynamics/active_retime.py':digest(ROOT/'protection/dynamics/active_retime.py'),
        'research/studies/solar_shield_array/active_retime_run.py':digest(__file__)}
    cases={str(d):run(d,a.displacement_km*1000,a.burn_hours,a.arcs) for d in a.days}
    result=dict(schema='terluna.research.active-orbital-retiming/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources)),
        evidence='Finite 10 km tile, actual ephemeris propagation and smooth electric-thrust arcs; endpoint compared with a same-date phase-rotated reference state',
        parameters=vars(a),propulsion_assumptions=SCENARIO['propulsion'],cases=cases,
        reading_rule='This measures individually reachable orbital phase changes. The target rotation specifies a boundary condition only; actual trajectories are independently propagated. No phase-shifted reference is counted as another fleet member. Multi-tile assignment, collisions, plume safety, coverage handovers, actuator dynamics and complete return schedules remain unproved. Power hardware, collector and habitat loads are excluded from the 50 g/m2 dynamics.',
        global_coverage_demonstrated=False,delivered_power_W=None)
    if any(digest(ROOT/n)!=h for n,h in sources.items()):raise ValueError('Sources changed during run')
    (HERE/'results'/a.output).write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
