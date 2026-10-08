"""Active finite-tile local choreography and separate propulsion accounting."""
import argparse
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.active_formation import (Formation, solar_frame,
    coverage_window, swept_square_conflicts)
from protection.dynamics.fleet import beam_margins
from protection.dynamics.cycling import ray_geometry, visible_sun
from protection.dynamics.optical import length
from .fleet_run import environment, identities, digest
from .cycling_search import HERE, ROOT, SCENARIO
from .cycling_analysis import compact_series_json


def source_ids():
    return {**identities(),
        'protection/dynamics/active_formation.py':digest(ROOT/'protection/dynamics/active_formation.py'),
        'research/studies/solar_shield_array/active_formation_run.py':digest(__file__),
        'research/studies/solar_shield_array/cycling_analysis.py':digest(HERE/'cycling_analysis.py')}


def run(cap,rows=25,hours=6.,suns=8,max_step=60.,rtol=2e-10,pitch=9750.,response=900.):
    start=time.monotonic(); duration=hours*3600
    env=environment(max(1.,hours/24))
    form=Formation(env,None,duration,rows=rows,cap=cap,suns=suns,pitch=pitch,response_s=response)
    basis=solar_frame(env.at(0)); u,b=basis[:,:2].T
    if np.dot(np.cross(u,b),env.frames(0)[:,2])>0:b=-b
    radius=15e6; phase=-np.pi/8
    q=radius*(u*np.cos(phase)+b*np.sin(phase))
    v=np.sqrt(K.MOON_GM/radius)*(-u*np.sin(phase)+b*np.cos(phase))
    def centre_f(t,y):
        return np.r_[y[3:],env.gravity(t,y[:3])+form.isolated_sail(t,y[None,:])[0]]
    centre=solve_ivp(centre_f,[0,duration],np.r_[q,v],dense_output=True,method='DOP853',
        rtol=2e-12,atol=[1e-5]*3+[1e-9]*3,max_step=30.)
    if not centre.success:raise RuntimeError(centre.message)
    form.centre=centre.sol
    initial,_=form.target(0)
    rng=np.random.default_rng(27041)
    perturb=rng.normal(size=(len(initial),3)); perturb*=20/np.maximum(length(perturb),1e-100)[:,None]
    initial[:,:3]+=perturb
    def stencil_guard(t,flat):
        state=flat.reshape(-1,6); target,_=form.target(t)
        return 1000.-length(state[:,:3]-target[:,:3]).max()
    stencil_guard.terminal=True; stencil_guard.direction=-1
    times=np.arange(0,duration+1,60.)
    sol=solve_ivp(form.derivative,[0,duration],initial.ravel(),t_eval=times,
        events=stencil_guard,method='DOP853',max_step=max_step,rtol=rtol,
        atol=np.tile([1e-4]*3+[1e-8]*3,len(initial)))
    if not sol.success:raise RuntimeError(sol.message)
    t=sol.t; states=sol.y.T.reshape(len(t),len(initial),6)
    control=[]; error=[]; light=[]; margins=[]; speed=[]; accel=[]
    coverage=[]; coverage_times=[]; natural_visibility=[]
    for k,(ti,state) in enumerate(zip(t,states)):
        natural,command,lit,target=form.forces(ti,state)
        control.append(length(command)); error.append(length(state[:,:3]-target[:,:3]))
        light.append(lit); accel.append(float(length(natural+command).max()))
        sample=env.at(ti); rays=ray_geometry(state[:,:3],sample)
        normal=-form.frames(ti)[:,0]
        margins.append(beam_margins(state[:,:3],rays['sun'],normal,sample,4*K.MOON_RADIUS,form.side))
        natural_visibility.append(float(visible_sun(rays['sun'],rays['earth'],rays['moon']).min()))
        if k%30==0 or k==len(t)-1:
            coverage_times.append(float(ti/3600))
            coverage.append(coverage_window(state,normal,sample,form.side,centre.sol(ti)))
        speed.append(float(length(form.frames(ti,1)[:,0])))
    control=np.asarray(control); error=np.asarray(error); light=np.asarray(light); margins=np.asarray(margins)
    elapsed=t[-1]-t[0]
    mean_acc=float(np.trapezoid(control.mean(axis=1),t)/elapsed)
    mean_force=float(mean_acc*len(initial)*form.side**2*env.sigma)
    prop=SCENARIO['propulsion']; cant=np.cos(np.deg2rad(prop['cant_deg']))
    peak_force=float(control.sum(axis=1).max()*form.side**2*env.sigma)
    power=mean_force*prop['exhaust_velocity_m_s']/(2*prop['efficiency']*cant)
    peak_power=peak_force*prop['exhaust_velocity_m_s']/(2*prop['efficiency']*cant)
    mass_flow=mean_force/(prop['exhaust_velocity_m_s']*cant)
    conflicts=swept_square_conflicts(t,states,form.frames,side=form.side)
    name=f'cap_{cap:.0e}_rows_{rows}_sun_{suns}_step_{max_step:g}'
    folder=ROOT/'research/runs/solar_shield_array/active'; folder.mkdir(parents=True,exist_ok=True)
    raw=folder/(name+'.npz'); np.savez_compressed(raw,t=t,state=states,control_m_s2=control,error_m=error,illumination=light)
    # Re-evaluate the same final state on a finer solar quadrature as a force
    # sensitivity check; this is separate from an independent trajectory replay.
    old=form.suns; _,coarse,_,_=form.forces(t[-1],states[-1]); form.suns=32
    _,fine,_,_=form.forces(t[-1],states[-1]); form.suns=old
    result=dict(config=dict(cap_m_s2=cap,rows=rows,hours=hours,solar_force_directions=suns,
        max_step_s=max_step,rtol=rtol,pitch_m=pitch,response_s=response),
        count=len(initial),side_m=form.side,clear_side_m=form.clear_side,
        optical_mass_kg=len(initial)*form.side**2*env.sigma,
        completed=sol.status==0,simulated_hours=float(t[-1]/3600),
        event_hours=[float(x/3600) for x in sol.t_events[0]],
        maximum_tracking_error_m=float(error.max()),final_tracking_error_m=float(error[-1].max()),
        formation_50m_gate=bool(error.max()<=50),
        mean_assist_acceleration_m_s2=mean_acc,peak_assist_acceleration_m_s2=float(control.max()),
        mean_delta_v_m_s=float(np.trapezoid(control.mean(axis=1),t)),
        maximum_tile_delta_v_m_s=float(np.trapezoid(control,t,axis=0).max()),
        propulsion_assumptions=prop,mean_assist_power_W=power,peak_assist_power_W=peak_power,
        propellant_rate_kg_s=mass_flow,propellant_used_kg=mass_flow*elapsed,
        added_peak_power_hardware_kg=peak_power*prop['peak_margin_factor']/prop['specific_power_W_kg'],
        mean_mutual_illumination_fraction=float(np.trapezoid(light.mean(axis=1),t)/elapsed),
        minimum_mutual_illumination_fraction=float(light.min()),
        minimum_earth_beam_margin_deg=float(np.rad2deg(margins[...,1].min())),
        minimum_moon_beam_margin_deg=float(np.rad2deg(margins[...,0].min())),
        minimum_natural_sun_visibility=float(min(natural_visibility)),
        maximum_normal_rate_rad_s=max(speed),maximum_centre_acceleration_m_s2=max(accel),
        collisions=conflicts,coverage_times_hours=coverage_times,coverage=coverage,
        solar_force_to_32_final_max_control_change_m_s2=float(length(fine-coarse).max()),
        raw_path=str(raw.relative_to(ROOT)),raw_sha256=digest(raw),
        numerical_evaluations=sol.nfev,elapsed_s=time.monotonic()-start)
    print(json.dumps({k:v for k,v in result.items() if k not in ('coverage','propulsion_assumptions')}),flush=True)
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--caps',nargs='+',type=float,default=[0,1e-5,1e-4,5e-4])
    p.add_argument('--rows',type=int,default=25); p.add_argument('--hours',type=float,default=6)
    p.add_argument('--suns',type=int,default=8); p.add_argument('--max-step',type=float,default=60.)
    p.add_argument('--rtol',type=float,default=2e-10); p.add_argument('--pitch',type=float,default=9750.)
    p.add_argument('--response',type=float,default=900.)
    p.add_argument('--output',default='active_formation.json'); a=p.parse_args()
    before=source_ids()
    cases={f'{cap:g}':run(cap,a.rows,a.hours,a.suns,a.max_step,a.rtol,a.pitch,a.response) for cap in a.caps}
    if source_ids()!=before:raise ValueError('Sources changed during propagation')
    product=dict(schema='terluna.research.active-tile-formation/1',
        producer=dict(source_hashes=before,constants=constants_used(before)),
        evidence='Independent finite-tile closed-loop ephemeris propagation; local moving 100 km diameter receiver window with explicit finite-Sun mutual photon forces',
        geometry='25 by 25 separate 10 km squares by default, 9750 m pitch, four 1 km-separated depth levels; target follows a natural 15000 km centre orbit',
        parameters=vars(a),cases=cases,
        acceptance=dict(global_four_lunar_radius_coverage_demonstrated=False,
                        delivered_electrical_power_demonstrated_W=None),
        limitations=['Local service passage only; global population, transfers, return arcs and handovers remain unresolved',
            'Bounded ideal electric acceleration; attitude hardware, structural flexure, plume interactions and thermal management omitted',
            'Power hardware and propellant are reported separately, without reinserting their mass into dynamics',
            'Finite squares share a commanded normal; collision guard includes assumed angular-rate and acceleration bounds',
            'Mutual spectral interception uses finite-Sun quadrature and a bounded nearest-neighbour stencil; integration stops at 1 km tracking error',
            'No collector or habitat loads are included'])
    (HERE/'results'/a.output).write_text(compact_series_json(product))


if __name__=='__main__':main()
