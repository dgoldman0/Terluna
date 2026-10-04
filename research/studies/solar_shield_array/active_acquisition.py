"""Close explicit local seam gaps with bounded control of separate 10 km tiles."""
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from shared import constants as K
from shared.provenance import constants_used
from protection.dynamics.active_formation import Formation,solar_frame,coverage_window,swept_square_conflicts
from protection.dynamics.fleet import beam_margins
from protection.dynamics.cycling import ray_geometry
from protection.dynamics.optical import length
from .fleet_run import environment,identities,digest
from .cycling_search import HERE,ROOT,SCENARIO
from .cycling_analysis import compact_series_json


def main():
    duration=10800.; env=environment(1)
    form=Formation(env,None,duration,rows=25,cap=5e-4,suns=8)
    basis=solar_frame(env.at(0));u,b=basis[:,:2].T
    if np.dot(np.cross(u,b),env.frames(0)[:,2])>0:b=-b
    radius=15e6; phase=-np.pi/8
    initial_c=np.r_[radius*(u*np.cos(phase)+b*np.sin(phase)),
        np.sqrt(K.MOON_GM/radius)*(-u*np.sin(phase)+b*np.cos(phase))]
    def cf(t,y):return np.r_[y[3:],env.gravity(t,y[:3])+form.isolated_sail(t,y[None,:])[0]]
    centre=solve_ivp(cf,[0,duration],initial_c,dense_output=True,method='DOP853',rtol=2e-12,
        atol=[1e-5]*3+[1e-9]*3,max_step=30.)
    if not centre.success:raise RuntimeError(centre.message)
    form.centre=centre.sol;initial,_=form.target(0)
    iy,ix=np.indices((25,25))
    # Alternate rows/columns create actual optical cracks, with every tile
    # initially 500 m from its assigned position. Depth clearance is preserved.
    offsets=np.column_stack([np.zeros(625),400*(2*(ix.ravel()%2)-1),300*(2*(iy.ravel()%2)-1)])
    initial[:,:3]+=offsets@basis.T
    def guard(t,flat):
        target,_=form.target(t)
        return 1000-length(flat.reshape(-1,6)[:,:3]-target[:,:3]).max()
    guard.terminal=True;guard.direction=-1
    sol=solve_ivp(form.derivative,[0,duration],initial.ravel(),method='DOP853',
        rtol=2e-11,atol=np.tile([1e-4]*3+[1e-8]*3,625),max_step=60.,
        t_eval=np.arange(0,duration+1,60.),events=guard)
    if not sol.success:raise RuntimeError(sol.message)
    t=sol.t;state=sol.y.T.reshape(len(t),625,6)
    controls=[];errors=[];coverage=[];margins=[]
    for k,(ti,y) in enumerate(zip(t,state)):
        _,control,_,target=form.forces(ti,y)
        controls.append(length(control));errors.append(float(length(y[:,:3]-target[:,:3]).max()))
        sample=env.at(ti);normal=-form.frames(ti)[:,0];rays=ray_geometry(y[:,:3],sample)
        margins.append(beam_margins(y[:,:3],rays['sun'],normal,sample,4*K.MOON_RADIUS,10000.))
        if k%10==0 or k==len(t)-1:
            coverage.append(dict(minutes=float(ti/60),**coverage_window(y,normal,sample,10000.,centre.sol(ti),suns=32)))
    controls=np.array(controls);margins=np.array(margins)
    prop=SCENARIO['propulsion'];cant=np.cos(np.deg2rad(prop['cant_deg']));mass=625*5e6
    dv=float(np.trapezoid(controls.mean(axis=1),t));mean=mass*dv/t[-1]
    folder=ROOT/'research/runs/solar_shield_array/active';folder.mkdir(exist_ok=True,parents=True)
    raw=folder/'acquisition.npz';np.savez_compressed(raw,t=t,state=state,control=controls,error=errors)
    sources={**identities(),'protection/dynamics/active_formation.py':digest(ROOT/'protection/dynamics/active_formation.py'),
        'research/studies/solar_shield_array/active_acquisition.py':digest(__file__),
        'research/studies/solar_shield_array/cycling_analysis.py':digest(HERE/'cycling_analysis.py')}
    result=dict(schema='terluna.research.active-seam-acquisition/1',
        producer=dict(source_hashes=sources,constants=constants_used(sources)),
        evidence='Bounded closed-loop ephemeris propagation of 625 finite squares from an explicitly gapped local formation',
        parameters=dict(rows=25,tile_side_m=10000.,pitch_m=9750.,layer_spacing_m=1000.,
            initial_alternating_offsets_m=[0,400,300],cap_m_s2=5e-4,force_solar_directions=8,coverage_solar_directions=32),
        completed=sol.status==0,simulated_hours=float(t[-1]/3600),coverage=coverage,
        errors_times_minutes=(t/60).tolist(),maximum_tile_error_m=errors,
        first_sample_with_50m_tracking_minutes=next((float(ti/60) for ti,e in zip(t,errors) if e<=50),None),
        first_sample_with_full_sun_coverage_minutes=next((c['minutes'] for c in coverage if c['all_sampled_sun_coverage']>=1-1e-10),None),
        mean_delta_v_m_s=dv,mean_assist_power_W=mean*prop['exhaust_velocity_m_s']/(2*prop['efficiency']*cant),
        propellant_used_kg=mass*dv/(prop['exhaust_velocity_m_s']*cant),propulsion_assumptions=prop,
        collisions=swept_square_conflicts(t,state,form.frames),
        minimum_earth_beam_margin_deg=float(np.rad2deg(margins[...,1].min())),
        minimum_moon_beam_margin_deg=float(np.rad2deg(margins[...,0].min())),
        raw_path=str(raw.relative_to(ROOT)),raw_sha256=digest(raw),
        reading_rule='Local 100 km diameter moving window; full-Sun fractions refer to sampled solar directions and coverage dates. The controller also pays to maintain the prescribed formation. This does not establish a global population or coverage during acquisition. Electrical loads are propulsion only; payload, hardware mass feedback, exhaust safety and actuator dynamics remain open.',
        global_coverage_demonstrated=False,delivered_power_W=None)
    (HERE/'results/active_acquisition.json').write_text(compact_series_json(result))
    print(json.dumps({k:v for k,v in result.items() if k not in ('producer','maximum_tile_error_m','errors_times_minutes')}),flush=True)


if __name__=='__main__':main()
